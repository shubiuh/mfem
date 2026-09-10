#ifndef UDAR_ROL_H
#define UDAR_ROL_H

#ifdef _WIN32
#define UDAR_API __declspec(dllexport)
#else
#define UDAR_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C"
{
#endif

   /** @brief Synchronous forward-model callback used by the ROL objective.
       @param n Number of natural-log conductivity parameters in @a x.
       @param x Borrowed parameter vector in natural-log conductivity.
       @param ndata Number of values required in @a prediction.
       @param prediction Caller-owned output vector of @a ndata predictions.
       @param context Opaque caller-owned context forwarded unchanged.
       @return Zero on success; a nonzero value stops the optimization.

       The callback is invoked serially and must not throw across this C ABI.
       Neither @a x nor @a context is retained by the inversion core. */
   typedef int (*udar_predict_fn)(int n, const double *x, int ndata,
                                  double *prediction, void *context);

   /** @brief Controls for the bounded ROL quasi-Newton solve. */
   typedef struct
   {
      /**< Maximum ROL outer iterations; must be positive. */
      int max_iterations;
      /**< Positive projected-gradient convergence tolerance. */
      double gradient_tolerance;
      /**< Positive step convergence tolerance. */
      double step_tolerance;
      /**< Positive finite-difference step in log-conductivity space. */
      double difference_step;
      /**< Nonnegative quadratic prior weight. */
      double regularization;
   } udar_options;

   /** @brief Summary of a completed ROL optimization attempt. */
   typedef struct
   {
      /**< ROL status: 0 converged, 1 iteration limit, or 2 step limit. */
      int status;
      /**< Number of ROL outer iterations completed. */
      int iterations;
      /**< Number of distinct forward callback evaluations. */
      int forward_evaluations;
      /**< Objective at the supplied initial model. */
      double initial_objective;
      /**< Objective at the final model. */
      double final_objective;
      /**< Final ROL projected-gradient norm. */
      double projected_gradient_norm;
   } udar_result;

   /** @brief Return the version of the UDAR inversion C ABI.
       @return The ABI version, currently 1. */
   UDAR_API int udar_abi_version(void);

   /** @brief Minimize weighted data misfit plus a quadratic log-conductivity
      prior.
       @param n Number of model parameters; must be positive.
       @param ndata Number of predicted and observed data values; must be
      positive.
       @param x In/out natural-log conductivity model; caller retains ownership.
       @param lower Borrowed strict lower bounds in natural-log conductivity.
       @param upper Borrowed strict upper bounds in natural-log conductivity.
       @param prior Borrowed natural-log conductivity prior.
       @param observed Borrowed observed data vector.
       @param stddev Borrowed positive data standard deviations.
       @param options Borrowed ROL controls.
       @param predict Synchronous forward callback; caller retains ownership.
       @param context Opaque callback context; caller retains ownership.
       @param log_path UTF-8 path for the ROL iteration log.
       @param result Caller-owned result record populated on success.
       @param error Caller-owned error buffer of @a error_capacity bytes.
       @param error_capacity Capacity of @a error, including the null
      terminator.
       @return Zero when the solve completed; -1 on invalid input or failure.

       Completion does not imply convergence: inspect @c result->status. The
      core retains no pointer arguments after this function returns. */
   UDAR_API int udar_minimize(int n, int ndata, double *x, const double *lower,
                              const double *upper, const double *prior,
                              const double *observed, const double *stddev,
                              const udar_options *options,
                              udar_predict_fn predict, void *context,
                              const char *log_path, udar_result *result,
                              char *error, int error_capacity);

#ifdef __cplusplus
}
#endif
#endif
