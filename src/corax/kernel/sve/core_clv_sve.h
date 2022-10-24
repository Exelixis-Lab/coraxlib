#ifndef CORAX_KERNEL_SVE_CORE_CLV_H_
#define CORAX_KERNEL_SVE_CORE_CLV_H_

#include "corax/core/common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in core_clv_sve.c */

  CORAX_EXPORT void
  corax_core_update_clv_ii_sve(unsigned int        states,
                               unsigned int        sites,
                               unsigned int        rate_cats,
                               double *            parent_clv,
                               unsigned int *      parent_scaler,
                               const double *      left_clv,
                               const double *      right_clv,
                               const double *      left_matrix,
                               const double *      right_matrix,
                               const unsigned int *left_scaler,
                               const unsigned int *right_scaler,
                               unsigned int        attrib);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_SVE_CORE_CLV_H_ */
