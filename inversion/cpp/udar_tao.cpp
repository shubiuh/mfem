#include "udar_tao.h"

#include <petsctao.h>

#if defined(PETSC_USE_COMPLEX) || !defined(PETSC_USE_REAL_DOUBLE) ||           \
    defined(PETSC_USE_64BIT_INDICES) ||                                        \
    defined(PETSC_HAVE_64BIT_BLAS_INDICES) || defined(PETSC_USE_DEBUG)
#error "udar_tao requires Release PETSc with real doubles and LP64 integers"
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

/** @brief Concrete storage for model, residual, and gradient vectors. */
using Values = std::vector<double>;

/** @brief Serialize PETSc initialization and TAO solves in this C ABI. */
std::mutex petsc_mutex;

/** @brief True after this library initializes PETSc for its first solve. */
bool petsc_initialized_here = false;

/** @brief Throw an invalid-argument exception when a precondition fails. */
void require(bool ok, const char *message)
{
   if (!ok) { throw std::invalid_argument(message); }
}

/** @brief Convert a PETSc error code to a C++ exception. */
void check(PetscErrorCode code, const char *operation)
{
   if (!code) { return; }
   const char *message = nullptr;
   PetscErrorMessage(code, &message, nullptr);
   throw std::runtime_error(std::string(operation) + ": " +
                            (message ? message : "PETSc error"));
}

/** @brief Initialize real, double-precision PETSc exactly once. */
void initialize_petsc()
{
   PetscBool initialized = PETSC_FALSE;
   PetscBool finalized = PETSC_FALSE;
   check(PetscInitialized(&initialized), "PetscInitialized");
   check(PetscFinalized(&finalized), "PetscFinalized");
   require(!finalized, "PETSc was already finalized in this process");
   if (!initialized)
   {
      check(PetscInitializeNoArguments(), "PetscInitializeNoArguments");
      petsc_initialized_here = true;
      check(PetscPushErrorHandler(PetscReturnErrorHandler, nullptr),
            "PetscPushErrorHandler");
   }
}

/** @brief Weighted least-squares objective evaluated through a forward
    callback. */
class DataObjective
{
   int n_, nd_;
   Values lo_, hi_, prior_, data_, sd_;
   udar_options opt_;
   udar_predict_fn predict_;
   void *context_;
   std::map<Values, Values> cache_;

 public:
   /** @brief Number of cache misses that invoked the forward callback. */
   int evaluations = 0;
   /** @brief Callback or validation error propagated out of TAO. */
   std::string callback_error;

   /** @brief Construct a callback-backed objective with copied inputs. */
   DataObjective(int n, int nd, const double *lo, const double *hi,
                 const double *prior, const double *data, const double *sd,
                 udar_options opt, udar_predict_fn predict, void *context)
       : n_(n), nd_(nd), lo_(lo, lo + n), hi_(hi, hi + n),
         prior_(prior, prior + n), data_(data, data + nd), sd_(sd, sd + nd),
         opt_(opt), predict_(predict), context_(context)
   {
   }

   /** @brief Clamp a model to its bounds after tolerance validation. */
   Values canonical(const Values &input) const
   {
      Values x(input);
      for (int j = 0; j < n_; ++j)
      {
         const double slack =
             64 * std::numeric_limits<double>::epsilon() *
             std::max({1.0, std::abs(lo_[j]), std::abs(hi_[j])});
         require(std::isfinite(x[j]) && x[j] >= lo_[j] - slack &&
                     x[j] <= hi_[j] + slack,
                 "TAO requested an invalid or out-of-bounds model");
         x[j] = std::clamp(x[j], lo_[j], hi_[j]);
      }
      return x;
   }

   /** @brief Compute or retrieve the weighted residual for a model. */
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

   /** @brief Evaluate the scalar objective at a canonical model. */
   double value(const Values &input)
   {
      const auto model = canonical(input);
      const auto r = residual(model);
      double objective = 0.0;
      for (double v : r) { objective += 0.5 * v * v; }
      for (int j = 0; j < n_; ++j)
      {
         objective +=
             0.5 * opt_.regularization * std::pow(model[j] - prior_[j], 2);
      }
      require(std::isfinite(objective), "Non-finite objective");
      return objective;
   }

   /** @brief Evaluate the objective and its feasible finite-difference
       gradient. */
   void value_gradient(const Values &input, double &objective, Values &gradient)
   {
      const auto model = canonical(input);
      const auto r = residual(model);
      objective = value(model);
      gradient.resize(n_);
      for (int j = 0; j < n_; ++j)
      {
         const double h =
             std::min(opt_.difference_step, 0.25 * (hi_[j] - lo_[j]));
         Values a(model), b(model);
         const bool central = model[j] - lo_[j] >= h && hi_[j] - model[j] >= h;
         const double sign =
             hi_[j] - model[j] >= model[j] - lo_[j] ? 1.0 : -1.0;
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
         require(a[j] != model[j] && b[j] != model[j],
                 "Finite-difference step is too small");
         const auto ra = residual(a);
         const auto rb = residual(b);
         double component = 0.0;
         for (int i = 0; i < nd_; ++i)
         {
            const double derivative =
                central ? (rb[i] - ra[i]) / (2 * h)
                        : sign * (-3 * r[i] + 4 * ra[i] - rb[i]) / (2 * h);
            component += r[i] * derivative;
         }
         gradient[j] = component + opt_.regularization * (model[j] - prior_[j]);
         require(std::isfinite(gradient[j]), "Non-finite objective gradient");
      }
   }
};

/** @brief Adapt the C++ objective to TAO's error-code callback contract. */
PetscErrorCode objective_gradient(Tao, Vec x, PetscReal *objective,
                                  Vec gradient, void *context) noexcept
{
   auto &problem = *static_cast<DataObjective *>(context);
   const PetscScalar *x_values = nullptr;
   PetscScalar *g_values = nullptr;
   PetscErrorCode code = VecGetArrayRead(x, &x_values);
   if (code) { return code; }
   code = VecGetArray(gradient, &g_values);
   if (code)
   {
      VecRestoreArrayRead(x, &x_values);
      return code;
   }
   try
   {
      PetscInt size = 0;
      check(VecGetSize(x, &size), "VecGetSize");
      Values model(x_values, x_values + size), values;
      double result = 0.0;
      problem.value_gradient(model, result, values);
      *objective = result;
      std::copy(values.begin(), values.end(), g_values);
   }
   catch (const std::exception &exception)
   {
      problem.callback_error = exception.what();
      code = PETSC_ERR_USER;
   }
   catch (...)
   {
      problem.callback_error = "Unknown exception in TAO objective callback";
      code = PETSC_ERR_USER;
   }
   const PetscErrorCode restore_gradient = VecRestoreArray(gradient, &g_values);
   const PetscErrorCode restore_model = VecRestoreArrayRead(x, &x_values);
   return code ? code : (restore_gradient ? restore_gradient : restore_model);
}

/** @brief Write one compact, machine-readable TAO monitor row. */
PetscErrorCode monitor(Tao tao, void *context) noexcept
{
   auto &log = *static_cast<std::ofstream *>(context);
   PetscInt iteration = 0;
   PetscReal objective = 0.0, gradient = 0.0, constraint = 0.0, step = 0.0;
   TaoConvergedReason reason = TAO_CONTINUE_ITERATING;
   const PetscErrorCode code = TaoGetSolutionStatus(
       tao, &iteration, &objective, &gradient, &constraint, &step, &reason);
   if (!code)
   {
      log << iteration << ',' << objective << ',' << gradient << ',' << step
          << ',' << static_cast<int>(reason) << '\n';
      log.flush();
   }
   return code;
}

/** @brief Map PETSc convergence reasons to the stable UDAR result contract. */
int status_code(TaoConvergedReason reason)
{
   switch (reason)
   {
   case TAO_CONVERGED_GATOL:
   case TAO_CONVERGED_GRTOL:
   case TAO_CONVERGED_GTTOL:
   case TAO_CONVERGED_MINF:
      return 0;
   case TAO_DIVERGED_MAXITS:
      return 1;
   case TAO_CONVERGED_STEPTOL:
      return 2;
   case TAO_DIVERGED_NAN:
      return 3;
   case TAO_CONVERGED_USER:
   case TAO_DIVERGED_USER:
      return 4;
   default:
      return 5;
   }
}

/** @brief Copy a standard vector into a sequential PETSc vector. */
void copy_to_vec(const Values &source, Vec destination)
{
   PetscScalar *values = nullptr;
   check(VecGetArray(destination, &values), "VecGetArray");
   std::copy(source.begin(), source.end(), values);
   check(VecRestoreArray(destination, &values), "VecRestoreArray");
}

/** @brief Copy a sequential PETSc vector into standard storage. */
Values copy_from_vec(Vec source)
{
   PetscInt size = 0;
   const PetscScalar *values = nullptr;
   check(VecGetSize(source, &size), "VecGetSize");
   check(VecGetArrayRead(source, &values), "VecGetArrayRead");
   Values result(values, values + size);
   check(VecRestoreArrayRead(source, &values), "VecRestoreArrayRead");
   return result;
}

} // namespace

/** @copydoc udar_abi_version */
extern "C" int udar_abi_version(void) { return 1; }

/** @copydoc udar_minimize */
extern "C" int udar_minimize(int n, int nd, double *x, const double *lower,
                             const double *upper, const double *prior,
                             const double *observed, const double *stddev,
                             const udar_options *options,
                             udar_predict_fn predict, void *context,
                             const char *log_path, udar_result *result,
                             char *error, int capacity)
{
   std::lock_guard<std::mutex> lock(petsc_mutex);
   if (error && capacity > 0) { error[0] = '\0'; }
   if (result)
   {
      *result = {};
      result->status = -1;
   }
   Tao tao = nullptr;
   Vec model = nullptr, lower_bound = nullptr, upper_bound = nullptr;
   try
   {
      require(n > 0 && nd > 0 && x && lower && upper && prior && observed &&
                  stddev && options && predict && result && log_path,
              "Missing inputs or invalid dimensions");
      require(options->max_iterations > 0 &&
                  std::isfinite(options->gradient_tolerance) &&
                  options->gradient_tolerance > 0 &&
                  std::isfinite(options->step_tolerance) &&
                  options->step_tolerance > 0 &&
                  std::isfinite(options->difference_step) &&
                  options->difference_step > 0 &&
                  std::isfinite(options->regularization) &&
                  options->regularization >= 0,
              "Invalid optimizer options");
      for (int j = 0; j < n; ++j)
      {
         require(std::isfinite(x[j]) && std::isfinite(lower[j]) &&
                     std::isfinite(upper[j]) && std::isfinite(prior[j]) &&
                     lower[j] < upper[j] && x[j] >= lower[j] &&
                     x[j] <= upper[j],
                 "Invalid model, prior or bounds");
      }
      for (int i = 0; i < nd; ++i)
      {
         require(
             std::isfinite(observed[i]) && std::isfinite(stddev[i]) &&
                 stddev[i] > 0,
             "Observations must be finite and standard deviations positive");
      }

      initialize_petsc();
      std::ofstream log(std::filesystem::u8path(log_path));
      require(log.good(), "Cannot create TAO log");
      log << "iteration,objective,projected_gradient_norm,step,tao_reason\n";
      DataObjective objective(n, nd, lower, upper, prior, observed, stddev,
                              *options, predict, context);
      const Values initial(x, x + n);
      result->initial_objective = objective.value(initial);

      check(VecCreateSeq(PETSC_COMM_SELF, n, &model), "VecCreateSeq(model)");
      check(VecDuplicate(model, &lower_bound), "VecDuplicate(lower)");
      check(VecDuplicate(model, &upper_bound), "VecDuplicate(upper)");
      copy_to_vec(initial, model);
      copy_to_vec(Values(lower, lower + n), lower_bound);
      copy_to_vec(Values(upper, upper + n), upper_bound);

      check(TaoCreate(PETSC_COMM_SELF, &tao), "TaoCreate");
      check(TaoSetType(tao, TAOBLMVM), "TaoSetType(TAOBLMVM)");
      check(TaoSetSolution(tao, model), "TaoSetSolution");
      check(TaoSetVariableBounds(tao, lower_bound, upper_bound),
            "TaoSetVariableBounds");
      check(TaoSetObjectiveAndGradient(tao, nullptr, objective_gradient,
                                       &objective),
            "TaoSetObjectiveAndGradient");
      check(TaoSetTolerances(tao, options->gradient_tolerance, 0.0, 0.0),
            "TaoSetTolerances");
      check(TaoSetMaximumIterations(tao, options->max_iterations),
            "TaoSetMaximumIterations");
      check(TaoMonitorSet(tao, monitor, &log, nullptr), "TaoMonitorSet");

      char step_tolerance[64];
      std::snprintf(step_tolerance, sizeof(step_tolerance), "%.17g",
                    options->step_tolerance);
      check(TaoSetOptionsPrefix(tao, "udar_"), "TaoSetOptionsPrefix");
      check(PetscOptionsSetValue(nullptr, "-udar_tao_steptol", step_tolerance),
            "PetscOptionsSetValue");
      const PetscErrorCode from_options = TaoSetFromOptions(tao);
      PetscOptionsClearValue(nullptr, "-udar_tao_steptol");
      check(from_options, "TaoSetFromOptions");

      const PetscErrorCode solve = TaoSolve(tao);
      if (!objective.callback_error.empty())
      {
         throw std::runtime_error(objective.callback_error);
      }
      check(solve, "TaoSolve");

      PetscInt iterations = 0;
      PetscReal final_objective = 0.0, gradient_norm = 0.0;
      PetscReal constraint_norm = 0.0, step = 0.0;
      TaoConvergedReason reason = TAO_CONTINUE_ITERATING;
      check(TaoGetSolutionStatus(tao, &iterations, &final_objective,
                                 &gradient_norm, &constraint_norm, &step,
                                 &reason),
            "TaoGetSolutionStatus");
      const auto final_model = copy_from_vec(model);
      std::copy(final_model.begin(), final_model.end(), x);
      result->status = status_code(reason);
      result->iterations = static_cast<int>(iterations);
      result->forward_evaluations = objective.evaluations;
      result->final_objective = final_objective;
      result->projected_gradient_norm = gradient_norm;

      TaoDestroy(&tao);
      VecDestroy(&upper_bound);
      VecDestroy(&lower_bound);
      VecDestroy(&model);
      return 0;
   }
   catch (const std::exception &exception)
   {
      if (error && capacity > 0)
      {
         std::snprintf(error, capacity, "%s", exception.what());
      }
   }
   catch (...)
   {
      if (error && capacity > 0)
      {
         std::snprintf(error, capacity, "Unknown C++ exception");
      }
   }
   if (tao) { TaoDestroy(&tao); }
   if (upper_bound) { VecDestroy(&upper_bound); }
   if (lower_bound) { VecDestroy(&lower_bound); }
   if (model) { VecDestroy(&model); }
   return -1;
}

/** @copydoc udar_finalize */
extern "C" int udar_finalize(void)
{
   std::lock_guard<std::mutex> lock(petsc_mutex);
   PetscBool initialized = PETSC_FALSE, finalized = PETSC_FALSE;
   if (PetscInitialized(&initialized) || PetscFinalized(&finalized))
   {
      return -1;
   }
   if (initialized && !finalized && petsc_initialized_here)
   {
      return PetscFinalize() ? -1 : 0;
   }
   return 0;
}
