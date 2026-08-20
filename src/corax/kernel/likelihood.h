#ifndef CORAX_KERNEL_LIKELIHOOD_H_
#define CORAX_KERNEL_LIKELIHOOD_H_

#include "corax/core/partition.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in likelihood.c */

  /**
   * Computes the likelihood given a single CLV. This is intended to be the CLV
   * of the "root" of the tree.
   *
   * @param partition The partition to compute the likelihood for.
   *
   * @param clv_index Index of the root CLV.
   *
   * @param scaler_index Index of the scalar buffer for the root CLV.
   *
   * @param freqs_indices An array of indices which indicate the per site base
   * distribution of states.
   *
   * @param[out] persite_lnl Buffer to store the individual site likelihoods.
   * Optional. Set to `nullptr` to ignore.
   *
   * @param[out] sitecat_lh  Buffer to store the individual per-site per-category likelihoods.
   * Optional. Set to `nullptr` to ignore.
   *
   * @return The total likelihood of the partition.
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT double
  corax_compute_root_loglikelihood_sitecat(corax_partition_t * partition,
                                           unsigned int        clv_index,
                                           int                 scaler_index,
                                           const unsigned int *freqs_indices,
                                           double *            persite_lnl,
                                           double *            sitecat_lh);


  CORAX_EXPORT double
  corax_compute_root_loglikelihood(corax_partition_t * partition,
                                   unsigned int        clv_index,
                                   int                 scaler_index,
                                   const unsigned int *freqs_indices,
                                   double *            persite_lnl);

  /**
   * Computes the likelihood of an edge. It does this by "rootinng" the tree at
   * parent, and computing the likelihood from there.
   *
   * @param partition The partition to compute the likelihood for.
   *
   * @param parent_clv_index Index of the parent CLV
   *
   * @param parent_scaler_index Index of the parent CLV scaler.
   *
   * @param child_clv_index Index of the child CLV
   *
   * @param child_parent_scaler_index Index of the child CLV scaler.
   *
   * @param matrix_index Index of the probability matrix between `parent` and
   * `child`.
   *
   * @param freqs_indices An array of indices which indicate the per site base
   * distribution of states.
   *
   * @param[out] persite_lnl Buffer to store the individual site likelihoods.
   * Optional. Set to `nullptr` to ignore.
   *
   * @param[out] sitecat_lh  Buffer to store the individual per-site per-category likelihoods.
   * Optional. Set to `nullptr` to ignore.
   *
   * @return The total likelihood of the partition.
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT double
  corax_compute_edge_loglikelihood_sitecat(corax_partition_t * partition,
                                   unsigned int        parent_clv_index,
                                   int                 parent_scaler_index,
                                   unsigned int        child_clv_index,
                                   int                 child_scaler_index,
                                   unsigned int        matrix_index,
                                   const unsigned int *freqs_indices,
                                   double *            persite_lnl,
                                   double *            sitecat_lh);

  CORAX_EXPORT double
  corax_compute_edge_loglikelihood(corax_partition_t * partition,
                                   unsigned int        parent_clv_index,
                                   int                 parent_scaler_index,
                                   unsigned int        child_clv_index,
                                   int                 child_scaler_index,
                                   unsigned int        matrix_index,
                                   const unsigned int *freqs_indices,
                                   double *            persite_lnl);

  CORAX_EXPORT int
  corax_compute_node_ancestral(corax_partition_t * partition,
                               unsigned int        node_clv_index,
                               int                 node_scaler_index,
                               unsigned int        other_clv_index,
                               int                 other_scaler_index,
                               unsigned int        matrix_index,
                               const unsigned int *freqs_indices,
                               double *            ancestral);

  CORAX_EXPORT int
  corax_compute_node_ancestral_extbuf(corax_partition_t * partition,
                                      unsigned int        node_clv_index,
                                      int                 node_scaler_index,
                                      unsigned int        other_clv_index,
                                      int                 other_scaler_index,
                                      unsigned int        pmatrix_index,
                                      const unsigned int *freqs_indices,
                                      double *            ancestral,
                                      double *            temp_clv,
                                      unsigned int *      temp_scaler,
                                      double *            ident_pmat);

  /**
   * Joint (Pupko et al. 1998) maximum-likelihood ancestral state
   * reconstruction: for every internal node in `operations`, computes a
   * probability vector (softmax of the log-score the winning traceback
   * path was chosen from -- NOT a marginal/empirical-Bayes posterior)
   * such that the arg max over states at every node, taken together,
   * forms a single globally most likely joint assignment of states to
   * all internal nodes (and, per site, to the rate category).
   *
   * `operations`/`op_count` is the same postorder operation list used
   * for ordinary sum-product likelihood computation. `root_clv_index` /
   * `root_other_clv_index` / `root_matrix_index` describe the virtual
   * root edge, as for `corax_compute_edge_loglikelihood_sitecat()`.
   * `root_scaler_index` / `root_other_scaler_index` are accepted for
   * interface symmetry but unused (everything is computed in log-space
   * directly from CLV/tip values).
   *
   * `node_probs` is caller-allocated, of length `op_count`; each
   * `node_probs[i]` must point to a buffer of `sites*states` doubles
   * and is indexed by operation index (`node_probs[i]` describes
   * `operations[i].parent_clv_index`).
   */
  CORAX_EXPORT int
  corax_compute_joint_ancestral(corax_partition_t *      partition,
                                const corax_operation_t *operations,
                                unsigned int              op_count,
                                unsigned int              root_clv_index,
                                int                       root_scaler_index,
                                unsigned int              root_other_clv_index,
                                int                       root_other_scaler_index,
                                unsigned int              root_matrix_index,
                                const unsigned int *      freqs_indices,
                                double **                 node_probs);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_LIKELIHOOD_H_ */
