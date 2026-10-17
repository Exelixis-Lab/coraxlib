#ifndef CORAX_MODEL_GAMMA_H_
#define CORAX_MODEL_GAMMA_H_

#include "corax/corax_core.h"

/* GAMMA discretization modes */
#define CORAX_GAMMA_RATES_MEAN 0
#define CORAX_GAMMA_RATES_MEDIAN 1
#define CORAX_GAMMA_RATES_WEIGHTS_LLOYD_MAX 2

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in gamma.c */

  /**
   * Compute incomplete gamma function.
   *
   * @param x Function value to evaluate
   * @param alpha Shape parameter alpha
   * @param ln_gamma_alpha natural logarithm of gamma function evaluated at alpha. Can be calculated using std::lgamma(alpha)
   */
  CORAX_EXPORT double corax_compute_incomplete_gamma(double x,
                                                     double alpha,
                                                     double ln_gamma_alpha);

  /**
   * Compute substitution rates for discretized gamma distribution.
   * Note that this function only works for discretization modes which use equal category weights (mean or median).
   *
   * @param alpha Shape parameter alpha
   * @param categories Number of categories
   * @param output_rates Pointer to array where category rates will be stored. Must be at least `categories` doubles long.
   * @param rates_mode Either `CORAX_GAMMA_RATES_MEAN` or `CORAX_GAMMA_RATES_MEDIAN`
   *
   * @returns `CORAX_SUCCESS` if operation succeeded and substitution rates were written, CORAX_FAILURE otherwise.
   */
  CORAX_EXPORT int corax_compute_gamma_rates(double       alpha,
                                             unsigned int categories,
                                             double      *output_rates,
                                             int          rates_mode);

  /**
   * Compute category weights and substitution rates for discretized gamma distribution.
   * 
   * @param alpha Shape parameter alpha
   * @param categories Number of categories
   * @param output_rates Pointer to array where category rates will be stored. Must be at least `categories` doubles long.
   * @param output_weights Pointer to array where category weights will be stored. Must be at least `categories` doubles long.
   * @param rates_mode Either `CORAX_GAMMA_RATES_MEAN` or `CORAX_GAMMA_RATES_MEDIAN` or `CORAX_GAMMA_RATES_WEIGHTS_LLOYD_MAX`
   *
   * @returns `CORAX_SUCCESS` if operation succeeded and substitution rates and weights were written, `CORAX_FAILURE` otherwise.
   */
  CORAX_EXPORT int corax_compute_gamma_cats(double       alpha,
                                            unsigned int categories,
                                            double      *output_rates,
                                            double      *output_weights,
                                            int          mode);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_MODEL_GAMMA_H_ */
