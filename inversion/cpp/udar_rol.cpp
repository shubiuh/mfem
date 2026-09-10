#include "udar_rol.h"
#include "ROL_StdVector.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_TypeB_QuasiNewtonAlgorithm.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace
{

/** @brief Concrete storage for model, residual, and data vectors. */
using Values = std::vector<double>;

/** @brief Throw an invalid-argument exception when a precondition fails.
    @param ok Condition to validate.
    @param message Static explanation used when @a ok is false. */
void require(bool ok, const char *message)
{
   if (!ok) { throw std::invalid_argument(message); }
}
/** @brief Access the storage of a ROL standard vector.
    @param x ROL vector that must have type @c ROL::StdVector<double>.
    @return Borrowed reference to the vector storage. */
const Values &values(const ROL::Vector<double> &x)
{ return *dynamic_cast<const ROL::StdVector<double> &>(x).getVector(); }

/** @brief Weighted least-squares objective evaluated through a forward
   callback.

    The class owns copies of all optimizer input vectors and caches residuals by
    canonical log-conductivity model. It borrows @c predict_ and @c context_
   only for the duration of @c udar_minimize. */
class DataObjective : public ROL::Objective<double>
{
   /** @brief Model and data dimensions. */
   int n_, nd_;
   /** @brief Lower bounds, upper bounds, prior, observations, and stddevs. */
   Values lo_, hi_, prior_, data_, sd_;
   /** @brief Copied finite-difference and stopping controls. */
   udar_options opt_;
   /** @brief Borrowed synchronous forward-model callback. */
   udar_predict_fn predict_;
   /** @brief Borrowed opaque callback context. */
   void *context_;
   /** @brief Residual vectors keyed by canonical log-conductivity model. */
   std::map<Values, Values> cache_;

 public:
   /** @brief Number of cache misses that invoked the forward callback. */
   int evaluations = 0;

   /** @brief Construct a callback-backed objective with copied vector inputs.
       @param n Model dimension.
       @param nd Data dimension.
       @param lo Borrowed lower-bound array of length @a n.
       @param hi Borrowed upper-bound array of length @a n.
       @param prior Borrowed prior array of length @a n.
       @param data Borrowed observation array of length @a nd.
       @param sd Borrowed positive standard-deviation array of length @a nd.
       @param opt Validated optimizer options.
       @param predict Borrowed forward callback.
       @param context Borrowed callback context. */
   DataObjective(int n, int nd, const double *lo, const double *hi,
                 const double *prior, const double *data, const double *sd,
                 udar_options opt, udar_predict_fn predict, void *context)
       : n_(n), nd_(nd), lo_(lo, lo + n), hi_(hi, hi + n),
         prior_(prior, prior + n), data_(data, data + nd), sd_(sd, sd + nd),
         opt_(opt), predict_(predict), context_(context)
   {
   }

   /** @brief Clamp a model to its bounds after a tolerance-based validation.
       @param input Candidate log-conductivity model.
       @return Canonical bounded model used as a cache key.

       Projection can exceed a bound by a few ULPs. Larger excursions are
       rejected to prevent invalid forward calls. */
   Values canonical(const Values &input) const
   {
      Values x(input);
      for (int j = 0; j < n_; ++j)
      {
         // x + (bound-x) can exceed a bound by a few ulps after projection.
         const double slack =
             64 * std::numeric_limits<double>::epsilon() *
             std::max({1.0, std::abs(lo_[j]), std::abs(hi_[j])});
         require(std::isfinite(x[j]) && x[j] >= lo_[j] - slack &&
                     x[j] <= hi_[j] + slack,
                 "ROL requested an invalid or out-of-bounds model");
         x[j] = std::clamp(x[j], lo_[j], hi_[j]);
      }
      return x;
   }

   /** @brief Compute or retrieve the weighted residual for a model.
       @param input Candidate log-conductivity model.
       @return Cached or newly computed `(prediction - data) / stddev` vector.
       @throws std::runtime_error When the forward callback reports failure. */
   Values residual(const Values &input)
   {
      const auto x = canonical(input);
      auto found = cache_.find(x);
      if (found != cache_.end()) { return found->second; }
      Values r(nd_, std::numeric_limits<double>::quiet_NaN());
      ++evaluations;
      if (predict_(n_, x.data(), nd_, r.data(), context_) != 0)
      {
         throw std::runtime_error("Forward callback failed; see Python "
                                  "exception and evaluation log");
      }
      for (int i = 0; i < nd_; ++i)
      {
         r[i] = (r[i] - data_[i]) / sd_[i];
         require(std::isfinite(r[i]),
                 "Non-finite prediction or weighted residual");
      }
      if (cache_.size() >= 256) { cache_.clear(); }
      cache_.emplace(x, r);
      return r;
   }

   /** @brief Evaluate weighted least squares plus quadratic prior penalty.
       @param x Current ROL log-conductivity model.
       @return Finite scalar objective value. */
   double value(const ROL::Vector<double> &x, double &) override
   {
      const auto m = canonical(values(x));
      const auto r = residual(m);
      double f = 0;
      for (double v : r) { f += 0.5 * v * v; }
      for (int j = 0; j < n_; ++j)
      {
         f += 0.5 * opt_.regularization * std::pow(m[j] - prior_[j], 2);
      }
      require(std::isfinite(f), "Non-finite objective");
      return f;
   }

   /** @brief Evaluate the finite-difference objective gradient.
       @param g Output gradient in log-conductivity coordinates.
       @param x Current ROL log-conductivity model.

       Central residual differences are used in the interior. Near a bound, a
       second-order inward difference preserves feasibility. */
   void gradient(ROL::Vector<double> &g, const ROL::Vector<double> &x,
                 double &) override
   {
      const auto m = canonical(values(x));
      auto &gv = *dynamic_cast<ROL::StdVector<double> &>(g).getVector();
      const auto r = residual(m);
      // Differentiate residuals, then J^T r. Central differences in the
      // interior; second-order inward differences near a parameter bound.
      for (int j = 0; j < n_; ++j)
      {
         double h = std::min(opt_.difference_step, 0.25 * (hi_[j] - lo_[j]));
         Values a(m), b(m);
         bool central = m[j] - lo_[j] >= h && hi_[j] - m[j] >= h;
         double sign = hi_[j] - m[j] >= m[j] - lo_[j] ? 1.0 : -1.0;
         if (central)
         {
            a[j] -= h;
            b[j] += h;
         }
         else
         {
            a[j] += sign * h;
            b[j] += sign * 2 * h;
         }
         require(a[j] != m[j] && b[j] != m[j],
                 "Finite-difference step is too small");
         const auto ra = residual(a), rb = residual(b);
         double grad = 0;
         for (int i = 0; i < nd_; ++i)
         {
            const double derivative =
                central ? (rb[i] - ra[i]) / (2 * h)
                        : sign * (-3 * r[i] + 4 * ra[i] - rb[i]) / (2 * h);
            grad += r[i] * derivative;
         }
         gv[j] = grad + opt_.regularization * (m[j] - prior_[j]);
         require(std::isfinite(gv[j]), "Non-finite objective gradient");
      }
   }
};
} // namespace

/** @copydoc udar_abi_version */
extern "C" int udar_abi_version(void) { return 1; }

/** @copydoc udar_minimize */
extern "C" int udar_minimize(int n, int nd, double *x, const double *lower,
                             const double *upper, const double *prior,
                             const double *observed, const double *stddev,
                             const udar_options *opt, udar_predict_fn predict,
                             void *context, const char *log_path,
                             udar_result *result, char *error, int capacity)
{
   if (error && capacity > 0) { error[0] = '\0'; }
   if (result)
   {
      *result = {};
      result->status = -1;
   }
   try
   {
      require(n > 0 && nd > 0 && x && lower && upper && prior && observed &&
                  stddev && opt && predict && result && log_path,
              "Missing inputs or invalid dimensions");
      require(
          opt->max_iterations > 0 && std::isfinite(opt->gradient_tolerance) &&
              opt->gradient_tolerance > 0 &&
              std::isfinite(opt->step_tolerance) && opt->step_tolerance > 0 &&
              std::isfinite(opt->difference_step) && opt->difference_step > 0 &&
              std::isfinite(opt->regularization) && opt->regularization >= 0,
          "Invalid optimizer options");
      for (int j = 0; j < n; ++j)
      {
         require(std::isfinite(x[j]) && std::isfinite(lower[j]) &&
                     std::isfinite(upper[j]) && std::isfinite(prior[j]) &&
                     lower[j] < upper[j] && x[j] >= lower[j] &&
                     x[j] <= upper[j] && std::isfinite(upper[j] - lower[j]),
                 "Invalid model, prior or bounds");
      }
      for (int i = 0; i < nd; ++i)
      {
         require(
             std::isfinite(observed[i]) && std::isfinite(stddev[i]) &&
                 stddev[i] > 0,
             "Observations must be finite and standard deviations positive");
      }

      std::ofstream log(std::filesystem::u8path(log_path));
      require(log.good(), "Cannot create ROL log");
      auto xv = ROL::makePtr<Values>(x, x + n);
      ROL::StdVector<double> model(xv);
      auto lv = ROL::makePtr<ROL::StdVector<double>>(
          ROL::makePtr<Values>(lower, lower + n));
      auto uv = ROL::makePtr<ROL::StdVector<double>>(
          ROL::makePtr<Values>(upper, upper + n));
      ROL::Bounds<double> bounds(lv, uv);
      DataObjective objective(n, nd, lower, upper, prior, observed, stddev,
                              *opt, predict, context);
      ROL::ParameterList settings;
      settings.sublist("Status Test")
          .set("Gradient Tolerance", opt->gradient_tolerance);
      settings.sublist("Status Test")
          .set("Step Tolerance", opt->step_tolerance);
      settings.sublist("Status Test")
          .set("Iteration Limit", opt->max_iterations);
      settings.sublist("General").set("Output Level", 1);
      settings.sublist("General").sublist("Secant").set("Type",
                                                        "Limited-Memory BFGS");
      ROL::TypeB::QuasiNewtonAlgorithm<double> algorithm(settings);
      double tol = 0;
      result->initial_objective = objective.value(model, tol);
      algorithm.run(model, objective, bounds, log);
      const auto state = algorithm.getState();
      const auto final_model = objective.canonical(*xv);
      std::copy(final_model.begin(), final_model.end(), x);
      result->status = static_cast<int>(state->statusFlag);
      result->iterations = state->iter;
      result->final_objective = objective.value(model, tol);
      result->forward_evaluations = objective.evaluations;
      result->projected_gradient_norm = state->gnorm;
      return 0;
   }
   catch (const std::exception &e)
   {
      if (error && capacity > 0)
      {
         std::snprintf(error, capacity, "%s", e.what());
      }
   }
   catch (...)
   {
      if (error && capacity > 0)
      {
         std::snprintf(error, capacity, "Unknown C++ exception");
      }
   }
   return -1;
}
