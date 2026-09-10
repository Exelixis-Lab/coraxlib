#ifndef CORAX_GAMMA_LAGUERRE_H_
#define CORAX_GAMMA_LAGUERRE_H_

#include "corax/corax_core.h"

#ifdef __cplusplus
extern "C"
{
#endif

CORAX_EXPORT int corax_compute_gamma_cats_laguerre(double       alpha,
                                                      unsigned int categories,
                                                      double      *output_rates,
                                                      double *output_weights);
#ifdef __cplusplus
}
#endif

#endif
