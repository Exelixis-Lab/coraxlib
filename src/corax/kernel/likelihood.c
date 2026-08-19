/*
    Copyright (C) 2015 Tomas Flouri, Diego Darriba

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "corax/corax.h"
#include <limits.h>
#include <float.h>

static double compute_asc_bias_correction(double       logl_base,
                                          unsigned int sum_w,
                                          unsigned int sum_w_inv,
                                          int          asc_bias_type)
{
  double logl_correction = 0.0;
  switch (asc_bias_type)
  {
  case CORAX_ATTRIB_AB_LEWIS:
    logl_correction = -(sum_w * log(1 - logl_base));
    break;
  case CORAX_ATTRIB_AB_STAMATAKIS:
    /* no need to add anything here */
    logl_correction = logl_base;
    break;
  case CORAX_ATTRIB_AB_FELSENSTEIN:
    logl_correction = sum_w_inv * log(logl_base);
    break;
  default:
    corax_set_error(CORAX_ERROR_AB_INVALIDMETHOD,
                    "Illegal ascertainment bias algorithm");
    return -INFINITY;
  }
  return logl_correction;
}

static double root_loglikelihood_asc_bias(corax_partition_t * partition,
                                          unsigned int        sites,
                                          const double *      clv,
                                          unsigned int *      scaler,
                                          const unsigned int *freqs_indices)
{
  unsigned int i, j, k;
  double       logl = 0;
  double       term_r, term;
  double       site_lk;

  const double *freqs         = NULL;
  unsigned int  states        = partition->states;
  unsigned int  states_padded = partition->states_padded;
  unsigned int  scale_factors;
  unsigned int *pattern_weights = partition->pattern_weights;
  double *      rate_weights    = partition->rate_weights;

  double       logl_correction = 0;
  unsigned int sum_w_inv       = 0;
  int          asc_bias_type   = partition->attributes & CORAX_ATTRIB_AB_MASK;

  /* point clvp to state sites */
  clv += sites * partition->rate_cats * partition->states_padded;

  /* 1. compute per-site logl for each state */
  for (i = 0; i < states; ++i)
  {
    term = 0;
    for (j = 0; j < partition->rate_cats; ++j)
    {
      freqs  = partition->frequencies[freqs_indices[j]];
      term_r = 0;
      for (k = 0; k < states; ++k) { term_r += clv[k] * freqs[k]; }
      term += term_r * rate_weights[j];
      clv += states_padded;
    }

    /* count number of scaling factors to acount for */
    scale_factors = scaler ? scaler[sites + i] : 0;

    sum_w_inv += pattern_weights[sites + i];
    if (asc_bias_type == CORAX_ATTRIB_AB_STAMATAKIS)
    {
      /* 2a. site_lk is the lnl weighted by the number of occurences */
      site_lk = log(term) * partition->pattern_weights[sites + i];
      if (scale_factors) site_lk += scale_factors * log(CORAX_SCALE_THRESHOLD);
    }
    else
    {
      /* 2b. site_lk is the actual likelihood */
      site_lk = term * pow(CORAX_SCALE_THRESHOLD, scale_factors);
    }

    /* logl_correction is the sum of weighted log likelihoods if
       asc_bias_type = Stamatakis, or the sum of likelihoods otherwise */
    logl_correction += site_lk;
  }

  /* 3. apply correction to the lnl score */
  logl += compute_asc_bias_correction(
      logl_correction, partition->pattern_weight_sum, sum_w_inv, asc_bias_type);

  return logl;
}

CORAX_EXPORT double
corax_compute_root_loglikelihood(corax_partition_t * partition,
                                 unsigned int        clv_index,
                                 int                 scaler_index,
                                 const unsigned int *freqs_indices,
                                 double *persite_lnl)
{
  return corax_compute_root_loglikelihood_sitecat(partition, clv_index,
                                                  scaler_index, freqs_indices,
                                                  persite_lnl, NULL);
}

CORAX_EXPORT double
corax_compute_root_loglikelihood_sitecat(corax_partition_t * partition,
                                         unsigned int        clv_index,
                                         int                 scaler_index,
                                         const unsigned int *freqs_indices,
                                         double *            persite_lnl,
                                         double *            sitecat_lh)
{
  double        logl = 0;
  unsigned int *scaler;
  unsigned int  identifiers;

  if (sitecat_lh)
  {
    corax_set_error(CORAX_ERROR_NOT_IMPLEMENTED,
                    "Per-category likelihoods for rooted trees not implemented!");
    return CORAX_FAILURE;
  }

  /* get scaler array if specified */
  if (scaler_index == CORAX_SCALE_BUFFER_NONE)
    scaler = NULL;
  else
    scaler = partition->scale_buffer[scaler_index];

  /* compute log-likelihood via the core function */
  if (corax_repeats_enabled(partition)
      && partition->repeats->pernode_ids[clv_index])
  {
    logl = corax_core_root_loglikelihood_repeats(
        partition->states,
        partition->sites,
        partition->rate_cats,
        partition->clv[clv_index],
        partition->repeats->pernode_site_id[clv_index],
        scaler,
        (const double **)partition->frequencies,
        partition->rate_weights,
        partition->pattern_weights,
        partition->prop_invar,
        partition->invariant,
        freqs_indices,
        persite_lnl,
        partition->attributes);
  }
  else
  {
    /* compute log-likelihood via the core function */
    logl =
        corax_core_root_loglikelihood(partition->states,
                                      partition->sites,
                                      partition->rate_cats,
                                      partition->clv[clv_index],
                                      scaler,
                                      (const double **)partition->frequencies,
                                      partition->rate_weights,
                                      partition->pattern_weights,
                                      partition->prop_invar,
                                      partition->invariant,
                                      freqs_indices,
                                      persite_lnl,
                                      partition->attributes);
  }

  /* ascertainment bias correction */
  if (partition->attributes & CORAX_ATTRIB_AB_MASK)
  {
    /* Note the assertion must be done for all rate matrices
    assert(prop_invar == 0);
    */
    identifiers =
        corax_get_sites_number(partition, clv_index) - partition->rate_cats;
    logl += root_loglikelihood_asc_bias(partition,
                                        identifiers,
                                        partition->clv[clv_index],
                                        scaler,
                                        freqs_indices);
  }

  return logl;
}

static double edge_loglikelihood_asc_bias_ti(corax_partition_t *partition,
                                             unsigned int  parent_clv_index,
                                             unsigned int *parent_scaler,
                                             unsigned int  matrix_index,
                                             const unsigned int *freqs_indices)
{
  unsigned int n, i, j;
  double       logl = 0;
  double       terma, terma_r;
  double       site_lk;

  const double *clvp          = partition->clv[parent_clv_index];
  const double *freqs         = NULL;
  const double *pmatrix       = partition->pmatrix[matrix_index];
  unsigned int  states        = partition->states;
  unsigned int  states_padded = partition->states_padded;
  unsigned int  scale_factors;
  unsigned int *pattern_weights = partition->pattern_weights;
  double *      rate_weights    = partition->rate_weights;

  double       logl_correction = 0;
  unsigned int sum_w_inv       = 0;
  int          asc_bias_type   = partition->attributes & CORAX_ATTRIB_AB_MASK;

  /* point clvp to state sites */
  clvp += partition->sites * partition->rate_cats * partition->states_padded;
  pattern_weights += partition->sites;
  if (parent_scaler) parent_scaler += partition->sites;

  /* 1. compute per-site logl for each state */
  for (n = 0; n < partition->states; ++n)
  {
    pmatrix = partition->pmatrix[matrix_index];
    terma   = 0;
    for (i = 0; i < partition->rate_cats; ++i)
    {
      freqs   = partition->frequencies[freqs_indices[i]];
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        terma_r += clvp[j] * freqs[j] * pmatrix[n];
        pmatrix += states_padded;
      }

      terma += terma_r * rate_weights[i];
      clvp += states_padded;
    }

    /* count number of scaling factors to acount for */
    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;

    sum_w_inv += pattern_weights[n];
    if (asc_bias_type == CORAX_ATTRIB_AB_STAMATAKIS)
    {
      /* 2a. site_lk is the lnl weighted by the number of occurences */
      site_lk = log(terma) * pattern_weights[n];
      if (scale_factors) site_lk += scale_factors * log(CORAX_SCALE_THRESHOLD);
    }
    else
    {
      /* 2b. site_lk is the actual likelihood */
      site_lk = terma * pow(CORAX_SCALE_THRESHOLD, scale_factors);
    }

    /* logl_correction is the sum of weighted log likelihoods if
       asc_bias_type = Stamatakis, or the sum of likelihoods otherwise */
    logl_correction += site_lk;
  }

  /* 3. apply correction to the lnl score */
  logl += compute_asc_bias_correction(
      logl_correction, partition->pattern_weight_sum, sum_w_inv, asc_bias_type);
  return logl;
}

static double edge_loglikelihood_tipinner(corax_partition_t *partition,
                                          unsigned int       parent_clv_index,
                                          int          parent_scaler_index,
                                          unsigned int child_clv_index,
                                          unsigned int matrix_index,
                                          const unsigned int *freqs_indices,
                                          double *            persite_lnl,
                                          double *            sitecat_lh)
{
  double       logl   = 0;
  unsigned int states = partition->states;

  unsigned int *parent_scaler;

  if (parent_scaler_index == CORAX_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[parent_scaler_index];

  if (states == 4)
  {
    logl = corax_core_edge_loglikelihood_ti_4x4(
        partition->sites,
        partition->rate_cats,
        partition->clv[parent_clv_index],
        parent_scaler,
        partition->tipchars[child_clv_index],
        partition->pmatrix[matrix_index],
        (const double **)partition->frequencies,
        partition->rate_weights,
        partition->pattern_weights,
        partition->prop_invar,
        partition->invariant,
        freqs_indices,
        persite_lnl,
        sitecat_lh,
        partition->attributes);
  }
  else
  {
    logl = corax_core_edge_loglikelihood_ti(
        partition->states,
        partition->sites,
        partition->rate_cats,
        partition->clv[parent_clv_index],
        parent_scaler,
        partition->tipchars[child_clv_index],
        partition->tipmap,
        partition->maxstates,
        partition->pmatrix[matrix_index],
        (const double **)partition->frequencies,
        partition->rate_weights,
        partition->pattern_weights,
        partition->prop_invar,
        partition->invariant,
        freqs_indices,
        persite_lnl,
        sitecat_lh,
        partition->attributes);
  }

  /* ascertainment bias correction */
  if (partition->attributes & CORAX_ATTRIB_AB_MASK)
  {
    /* Note the assertion must be done for all rate matrices
    assert(prop_invar == 0);
    */
    logl += edge_loglikelihood_asc_bias_ti(partition,
                                           parent_clv_index,
                                           parent_scaler,
                                           matrix_index,
                                           freqs_indices);
  }

  return logl;
}

static double edge_loglikelihood_asc_bias_ii(corax_partition_t * partition,
                                             const double *      clvp,
                                             unsigned int *      parent_scaler,
                                             unsigned int        parent_sites,
                                             const double *      clvc,
                                             unsigned int *      child_scaler,
                                             unsigned int        child_sites,
                                             unsigned int        matrix_index,
                                             const unsigned int *freqs_indices)
{
  unsigned int n, i, j, k;
  double       logl = 0;
  double       terma, terma_r, termb;
  double       site_lk;

  const double *freqs         = NULL;
  const double *pmatrix       = partition->pmatrix[matrix_index];
  unsigned int  states        = partition->states;
  unsigned int  states_padded = partition->states_padded;
  unsigned int  scale_factors;
  unsigned int *pattern_weights = partition->pattern_weights;
  double *      rate_weights    = partition->rate_weights;

  /* point clv, clvc, scalers and pattern weights to state sites */
  unsigned int parent_offset =
      parent_sites * partition->rate_cats * partition->states_padded;
  unsigned int child_offset =
      child_sites * partition->rate_cats * partition->states_padded;

  pattern_weights += partition->sites;
  clvp += parent_offset;
  clvc += child_offset;
  if (parent_scaler) parent_scaler += parent_sites;
  if (child_scaler) child_scaler += child_sites;

  double       logl_correction = 0;
  unsigned int sum_w_inv       = 0;
  int          asc_bias_type   = partition->attributes & CORAX_ATTRIB_AB_MASK;

  /* 1. compute per-site logl for each state */
  for (n = 0; n < partition->states; ++n)
  {
    pmatrix = partition->pmatrix[matrix_index];
    terma   = 0;
    for (i = 0; i < partition->rate_cats; ++i)
    {
      freqs   = partition->frequencies[freqs_indices[i]];
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        termb = 0;
        for (k = 0; k < states; ++k) { termb += pmatrix[k] * clvc[k]; }
        terma_r += clvp[j] * freqs[j] * termb;
        pmatrix += states_padded;
      }

      terma += terma_r * rate_weights[i];
      clvp += states_padded;
      clvc += states_padded;
    }

    /* count number of scaling factors to acount for */
    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

    sum_w_inv += pattern_weights[n];
    if (asc_bias_type == CORAX_ATTRIB_AB_STAMATAKIS)
    {
      /* 2a. site_lk is the lnl weighted by the number of occurences */
      site_lk = log(terma) * pattern_weights[n];
      if (scale_factors) site_lk += scale_factors * log(CORAX_SCALE_THRESHOLD);
    }
    else
    {
      /* 2b. site_lk is the actual likelihood */
      site_lk = terma * pow(CORAX_SCALE_THRESHOLD, scale_factors);
    }

    /* logl_correction is the sum of weighted log likelihoods if
       asc_bias_type = Stamatakis, or the sum of likelihoods otherwise */
    logl_correction += site_lk;
  }

  /* 3. apply correction to the lnl score */
  logl += compute_asc_bias_correction(
      logl_correction, partition->pattern_weight_sum, sum_w_inv, asc_bias_type);

  return logl;
}

static double edge_loglikelihood(corax_partition_t * partition,
                                 unsigned int        parent_clv_index,
                                 int                 parent_scaler_index,
                                 unsigned int        child_clv_index,
                                 int                 child_scaler_index,
                                 unsigned int        matrix_index,
                                 const unsigned int *freqs_indices,
                                 double *persite_lnl, double *sitecat_lh) {
  double logl = 0;

  const double *clvp = partition->clv[parent_clv_index];
  const double *clvc = partition->clv[child_clv_index];

  unsigned int *parent_scaler;
  unsigned int *child_scaler;

  if (child_scaler_index == CORAX_SCALE_BUFFER_NONE)
    child_scaler = NULL;
  else
    child_scaler = partition->scale_buffer[child_scaler_index];

  if (parent_scaler_index == CORAX_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[parent_scaler_index];

  /* compute log-likelihood via the core function */
  logl =
      corax_core_edge_loglikelihood_ii(partition->states,
                                       partition->sites,
                                       partition->rate_cats,
                                       clvp,
                                       parent_scaler,
                                       clvc,
                                       child_scaler,
                                       partition->pmatrix[matrix_index],
                                       (const double **)partition->frequencies,
                                       partition->rate_weights,
                                       partition->pattern_weights,
                                       partition->prop_invar,
                                       partition->invariant,
                                       freqs_indices,
                                       persite_lnl,
                                       sitecat_lh, partition->attributes);

  /* ascertainment bias correction */
  if (partition->attributes & CORAX_ATTRIB_AB_MASK)
  {
    /* Note the assertion must be done for all rate matrices
    assert(prop_invar == 0);
    */
    logl += edge_loglikelihood_asc_bias_ii(partition,
                                           clvp,
                                           parent_scaler,
                                           partition->sites,
                                           clvc,
                                           child_scaler,
                                           partition->sites,
                                           matrix_index,
                                           freqs_indices);
  }
  return logl;
}

static double edge_loglikelihood_repeats(corax_partition_t *partition,
                                         unsigned int       parent_clv_index,
                                         int                parent_scaler_index,
                                         unsigned int       child_clv_index,
                                         int                child_scaler_index,
                                         unsigned int       matrix_index,
                                         const unsigned int *freqs_indices,
                                         double *            persite_lnl,
                                         double *            sitecat_lh)
{
  double logl = 0;

  const double *clvp = partition->clv[parent_clv_index];
  const double *clvc = partition->clv[child_clv_index];

  const unsigned int *parent_site_id =
      corax_get_site_id(partition, parent_clv_index);
  const unsigned int *child_site_id =
      corax_get_site_id(partition, child_clv_index);
  unsigned int parent_sites =
      corax_get_sites_number(partition, parent_clv_index);
  unsigned int child_sites = corax_get_sites_number(partition, child_clv_index);
  unsigned int inv         = parent_sites > child_sites;

  unsigned int *parent_scaler;
  unsigned int *child_scaler;

  if (child_scaler_index == CORAX_SCALE_BUFFER_NONE)
    child_scaler = NULL;
  else
    child_scaler = partition->scale_buffer[child_scaler_index];

  if (parent_scaler_index == CORAX_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[parent_scaler_index];

  /* compute log-likelihood via the core function */
  logl = corax_core_edge_loglikelihood_repeats(
      partition->states,
      partition->sites,
      !inv ? parent_sites : child_sites,
      partition->rate_cats,
      inv ? clvp : clvc,
      inv ? parent_scaler : child_scaler,
      !inv ? clvp : clvc,
      !inv ? parent_scaler : child_scaler,
      partition->pmatrix[matrix_index],
      (const double **)partition->frequencies,
      partition->rate_weights,
      partition->pattern_weights,
      partition->prop_invar,
      partition->invariant,
      freqs_indices,
      persite_lnl,
      sitecat_lh,
      inv ? parent_site_id : child_site_id,
      !inv ? parent_site_id : child_site_id,
      partition->repeats->bclv_buffer,
      partition->attributes);

  /* ascertainment bias correction */
  if (partition->attributes & CORAX_ATTRIB_AB_MASK)
  {
    /* Note the assertion must be done for all rate matrices
    assert(prop_invar == 0);
    */
    logl += edge_loglikelihood_asc_bias_ii(partition,
                                           clvp,
                                           parent_scaler,
                                           parent_sites - partition->states,
                                           clvc,
                                           child_scaler,
                                           child_sites - partition->states,
                                           matrix_index,
                                           freqs_indices);
  }
  return logl;
}

CORAX_EXPORT double
corax_compute_edge_loglikelihood(corax_partition_t * partition,
                                 unsigned int        parent_clv_index,
                                 int                 parent_scaler_index,
                                 unsigned int        child_clv_index,
                                 int                 child_scaler_index,
                                 unsigned int        matrix_index,
                                 const unsigned int *freqs_indices,
                                 double             *persite_lnl)
{
  return corax_compute_edge_loglikelihood_sitecat(partition,
                                                  parent_clv_index, parent_scaler_index,
                                                  child_clv_index, child_scaler_index,
                                                  matrix_index, freqs_indices,
                                                  persite_lnl, NULL);
}

CORAX_EXPORT double
corax_compute_edge_loglikelihood_sitecat(corax_partition_t * partition,
                                 unsigned int        parent_clv_index,
                                 int                 parent_scaler_index,
                                 unsigned int        child_clv_index,
                                 int                 child_scaler_index,
                                 unsigned int        matrix_index,
                                 const unsigned int *freqs_indices,
                                 double             *persite_lnl,
                                 double             *sitecat_lh)
{
  double logl;

  if (corax_repeats_enabled(partition)
      && (partition->repeats->pernode_ids[parent_clv_index]
          || partition->repeats->pernode_ids[child_clv_index]))
  {
    return edge_loglikelihood_repeats(partition,
                                      parent_clv_index,
                                      parent_scaler_index,
                                      child_clv_index,
                                      child_scaler_index,
                                      matrix_index,
                                      freqs_indices,
                                      persite_lnl,
                                      sitecat_lh);
  }

  if (partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
    if ((parent_clv_index < partition->tips)
        || (child_clv_index < partition->tips))
      return edge_loglikelihood_tipinner(
          partition,
          (parent_clv_index < partition->tips) ? child_clv_index
                                               : parent_clv_index,
          (parent_clv_index < partition->tips) ? child_scaler_index
                                               : parent_scaler_index,
          (parent_clv_index < partition->tips) ? parent_clv_index
                                               : child_clv_index,
          matrix_index,
          freqs_indices,
          persite_lnl,
          sitecat_lh);


  logl = edge_loglikelihood(partition,
                            parent_clv_index,
                            parent_scaler_index,
                            child_clv_index,
                            child_scaler_index,
                            matrix_index,
                            freqs_indices,
                            persite_lnl,
                            sitecat_lh);

  return logl;
}

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
                                    double *            ident_pmat)
{
  if (!partition || !ancestral)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "Parameter value is NULL!");
    return CORAX_FAILURE;
  }

  if (!temp_clv || !temp_scaler || !ident_pmat)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "NULL buffer pointer");
    return CORAX_FAILURE;
  }

  if (corax_repeats_enabled(partition))
  {
    corax_set_error(
        CORAX_ERROR_EINVAL,
        "Site repeats are not compatible with ancestral state reconstruction!");
    return CORAX_FAILURE;
  }

  unsigned int n, i, j;
  unsigned int states        = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int sites         = partition->sites;
  unsigned int rate_cats     = partition->rate_cats;

  const double *pmat = partition->pmatrix[pmatrix_index];

  const double *node_clv = partition->clv[node_clv_index];

  unsigned int *node_scaler = (node_scaler_index == CORAX_SCALE_BUFFER_NONE)
                                  ? NULL
                                  : partition->scale_buffer[node_scaler_index];

  if (other_clv_index < partition->tips
      && (partition->attributes & CORAX_ATTRIB_PATTERN_TIP))
  {
    corax_core_update_clv_ti(states,
                             sites,
                             rate_cats,
                             temp_clv,
                             temp_scaler,
                             partition->tipchars[other_clv_index],
                             node_clv,
                             pmat,
                             ident_pmat,
                             node_scaler,
                             partition->tipmap,
                             partition->maxstates,
                             partition->attributes);
  }
  else
  {
    const double *other_clv = partition->clv[other_clv_index];
    unsigned int *other_scaler =
        (other_scaler_index == CORAX_SCALE_BUFFER_NONE)
            ? NULL
            : partition->scale_buffer[other_scaler_index];

    corax_core_update_clv_ii(states,
                             sites,
                             rate_cats,
                             temp_clv,
                             temp_scaler,
                             node_clv,
                             other_clv,
                             ident_pmat,
                             pmat,
                             node_scaler,
                             other_scaler,
                             partition->attributes);
  }

  double *clvp = temp_clv;
  double *ancp = ancestral;

  memset(ancestral, 0, sites * states * sizeof(double));

  /* PER-RATE scaling stuff */
  unsigned int  site_scalings;
  int           per_rate_scaling = (partition->attributes & CORAX_ATTRIB_RATE_SCALERS) ? 1 : 0;
  unsigned int *rate_scalings    = NULL;
  double        scale_minlh[CORAX_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling || partition->prop_invar)
  {
    double scale_factor = 1.0;
    for (i = 0; i < CORAX_SCALE_RATE_MAXDIFF; ++i)
    {
      scale_factor *= CORAX_SCALE_THRESHOLD;
      scale_minlh[i] = scale_factor;
    }
  }
  if (per_rate_scaling)
  {
    rate_scalings = (unsigned int *)calloc(rate_cats, sizeof(unsigned int));

    if (!rate_scalings)
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Cannot allocate space for rate scalers.");
      return CORAX_FAILURE;
    }
  }

  for (n = 0; n < sites; ++n)
  {
    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = temp_scaler[n * rate_cats + i];
        if (rate_scalings[i] < site_scalings) site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = CORAX_MIN(rate_scalings[i] - site_scalings,
                                     CORAX_SCALE_RATE_MAXDIFF);
      }
    }

    for (i = 0; i < rate_cats; ++i)
    {
      double *freqs       = partition->frequencies[freqs_indices[i]];
      double  rate_weight = partition->rate_weights[i];

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        rate_weight *= scale_minlh[rate_scalings[i] - 1];
      }

      for (j = 0; j < states; ++j)
      {
        ancp[j] += clvp[j] * freqs[j] * rate_weight;
      }

      clvp += states_padded;
    }

    // normalize probs
    double sum = 0;
    for (j = 0; j < states; ++j) sum += ancp[j];
    for (j = 0; j < states; ++j) ancp[j] /= sum;

    ancp += states;
  }

  if (rate_scalings) free(rate_scalings);

  return CORAX_SUCCESS;
}

CORAX_EXPORT int corax_compute_node_ancestral(corax_partition_t *partition,
                                              unsigned int       node_clv_index,
                                              int          node_scaler_index,
                                              unsigned int other_clv_index,
                                              int          other_scaler_index,
                                              unsigned int matrix_index,
                                              const unsigned int *freqs_indices,
                                              double *            ancestral)
{
  int          retval = CORAX_FAILURE;
  unsigned int i, j, k;

  unsigned int states        = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int sites         = partition->sites;
  unsigned int rate_cats     = partition->rate_cats;

  unsigned int clv_size =
      corax_get_clv_size(partition, node_clv_index) * sizeof(double);
  double *temp_clv =
      (double *)corax_aligned_alloc(clv_size, partition->alignment);
  unsigned int scaler_size =
      ((partition->attributes & CORAX_ATTRIB_RATE_SCALERS) ? sites * rate_cats
                                                           : sites)
      * sizeof(unsigned int);
  unsigned int *temp_scaler =
      (unsigned int *)corax_aligned_alloc(scaler_size, partition->alignment);
  unsigned int displacement =
      (states_padded - states) * (states_padded) * sizeof(double);
  unsigned int pmat_size =
      states_padded * states * rate_cats * sizeof(double) + displacement;
  double *ident_pmat =
      (double *)corax_aligned_alloc(pmat_size, partition->alignment);

  if (!temp_clv || !temp_scaler || !ident_pmat)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Cannot allocate memory");
    goto cleanup;
  }

  // init identity matrix
  memset(ident_pmat, 0, pmat_size);
  double *pmat = ident_pmat;
  for (i = 0; i < rate_cats; ++i)
  {
    for (j = 0; j < states; ++j)
    {
      for (k = 0; k < states_padded; ++k)
        pmat[j * states_padded + k] = (j == k) ? 1 : 0;
    }

    pmat += states * states_padded;
  }

  retval = corax_compute_node_ancestral_extbuf(partition,
                                               node_clv_index,
                                               node_scaler_index,
                                               other_clv_index,
                                               other_scaler_index,
                                               matrix_index,
                                               freqs_indices,
                                               ancestral,
                                               temp_clv,
                                               temp_scaler,
                                               ident_pmat);

cleanup:
  corax_aligned_free(temp_clv);
  corax_aligned_free(temp_scaler);
  corax_aligned_free(ident_pmat);

  return retval;
}

/* ===========================================================================
   NEW CODE -- append this block to the END of likelihood.c (after the
   existing corax_compute_node_ancestral() function). Nothing above this
   point in likelihood.c is modified.

   JOINT (Pupko et al. 1998) maximum-likelihood ancestral state
   reconstruction, as opposed to the EMPIRICAL-BAYES per-node posterior
   already provided by corax_compute_node_ancestral[_extbuf]().

   ---------------------------------------------------------------------
   Read this before using:

   * PROBABILITIES ONLY: this function computes, per internal node, a
     softmax-normalized joint log-score vector (see node_probs below).
     It does NOT compute or output per-site best/argmax states, and
     carries no backpointers or traceback -- that is done elsewhere in
     the codebase, against a corax_operation_t/tree structure this file
     does not have visibility into. If you need the joint MAP state
     per site again here, that logic (backpointer arrays + traceback)
     will need to be reintroduced; see version history / prior
     conversation for a previous implementation of it.

   * Supports both plain (unpacked) CLV storage AND
     CORAX_ATTRIB_PATTERN_TIP (packed tip characters + tipmap lookup),
     per corax_core_update_clv_ti()/corax_core_create_lookup() in
     core_clv.c: a tip's per-state likelihood at (site,state) is
     ((partition->tipmap[partition->tipchars[clv_index][site]] >> state)
     & 1) ? 1.0 : 0.0, independent of rate category. Site-repeats
     partitions are still rejected with an error -- they use a
     different per-node site-count/site-id indirection
     (corax_get_sites_number()/corax_get_site_id()) that this DP does
     not (yet) thread through; this mirrors the existing empirical-Bayes
     function's own unconditional rejection of repeats.

   * Rate-category heterogeneity: the rate category is treated as one
     more jointly-maximized hidden variable per site (we compare the
     best fully-reconstructed path under each rate category, weighted
     by its prior rate_weight, and keep the overall best one) purely to
     decide WHICH rate category's log-score table to report probs from
     per site. This is a modeling choice, not the per-rate
     marginalization scheme corax_compute_node_ancestral_extbuf() uses.

   * All computation is done in log-space (max-sum, i.e. Viterbi-style),
     so no separate scale-buffer bookkeeping is required for numerical
     stability.

   * Output: node_probs[op][site*states + s] -- a softmax-normalized
     version of the log-joint-score table at that node, evaluated at
     the globally-chosen best rate category for that site. This is a
     legitimate "confidence-like" vector but is NOT a marginal/
     empirical-Bayes posterior probability.
   =========================================================================== */

#include <limits.h>
#include <float.h>

static double corax_joint_safe_log(double x)
{
  return (x > 0.) ? log(x) : -INFINITY;
}

/* Read a single (site, rate, state) log-likelihood value for a clv_index
   that we do NOT compute ourselves (i.e. it is not the parent of any
   operation) -- under "plain" CLV storage this can be any node (tip or
   otherwise unvisited), under CORAX_ATTRIB_PATTERN_TIP storage this must
   be an actual tip (clv_index < partition->tips), since inner CLVs are
   always stored unpacked in partition->clv[] regardless of this
   attribute (see case_tiptip/case_tipinner in clv.c). */
static double corax_joint_leaf_logclv(const corax_partition_t *partition,
                                      unsigned int              clv_index,
                                      unsigned int              site,
                                      unsigned int              rate,
                                      unsigned int              state)
{
  if ((partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
      && clv_index < partition->tips)
  {
    unsigned char charcode = partition->tipchars[clv_index][site];
    corax_state_t smask    = partition->tipmap[charcode];

    return ((smask >> state) & 1u) ? 0.0 : -INFINITY; /* log(1) = 0 */
  }

  unsigned int  states_padded = partition->states_padded;
  unsigned int  rate_cats     = partition->rate_cats;
  const double *clv           = partition->clv[clv_index];

  double v =
      clv[(size_t)site * rate_cats * states_padded + (size_t)rate * states_padded
          + state];

  return corax_joint_safe_log(v);
}

/*
   Core per-partition joint-reconstruction probability engine.

   operations/op_count: the SAME postorder operation list used for
     ordinary sum-product likelihood computation (e.g. as built by
     corax_utree_create_operations() inside treeinfo_compute_loglh()).
     Children are guaranteed to be processed before their parents.

   root_clv_index / root_scaler_index / root_other_clv_index /
   root_other_scaler_index / root_matrix_index: the same "virtual root
     edge" quadruple used by corax_compute_edge_loglikelihood_sitecat(),
     i.e. treeinfo->root and treeinfo->root->back. (Scaler indices are
     unused here since everything is computed in log-space directly
     from CLV values, but are accepted for interface symmetry / future
     use and to make misuse obvious if someone passes CLV indices from
     an inconsistent rooting.)

   node_probs: caller-allocated array of length op_count. node_probs[i]
     must point to a buffer of size sites*states doubles. Indexed by
     OPERATION index i, i.e. node_probs[i] describes
     operations[i].parent_clv_index.
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
                              double **                 node_probs)
{
  CORAX_UNUSED(root_scaler_index);
  CORAX_UNUSED(root_other_scaler_index);

  if (!partition || !operations || !node_probs || !freqs_indices)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "Parameter value is NULL!");
    return CORAX_FAILURE;
  }

  if (corax_repeats_enabled(partition))
  {
    corax_set_error(
        CORAX_ERROR_EINVAL,
        "Site repeats are not compatible with ancestral state reconstruction!");
    return CORAX_FAILURE;
  }

  if (op_count == 0)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "op_count must be > 0");
    return CORAX_FAILURE;
  }

  unsigned int states        = partition->states;
  unsigned int states_padded = partition->states_padded;
  unsigned int sites         = partition->sites;
  unsigned int rate_cats     = partition->rate_cats;

  unsigned int i, n, r, s, k, op;
  int          retval = CORAX_FAILURE;

  /* ---- map clv_index -> operation index, for clv indices that WE
     compute (i.e. internal nodes). Anything not in this map is treated
     as a raw leaf CLV read directly from partition->clv[]/tipchars. ---- */
  unsigned int max_clv_index = root_clv_index;
  if (root_other_clv_index > max_clv_index) max_clv_index = root_other_clv_index;
  for (i = 0; i < op_count; ++i)
    if (operations[i].parent_clv_index > max_clv_index)
      max_clv_index = operations[i].parent_clv_index;

  unsigned int *clv_to_op =
      (unsigned int *)malloc(((size_t)max_clv_index + 1) * sizeof(unsigned int));

  /* per-node (per-operation) log-space joint table: op_count buffers of
     size sites * rate_cats * states. No backpointers -- states are not
     computed here. */
  double **jtable = (double **)calloc(op_count, sizeof(double *));

  double *      best_total = (double *)malloc(sites * sizeof(double));
  unsigned int *best_rate  = (unsigned int *)malloc(sites * sizeof(unsigned int));

  if (!clv_to_op || !jtable || !best_total || !best_rate)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Cannot allocate memory");
    goto cleanup;
  }

  for (i = 0; i <= max_clv_index; ++i) clv_to_op[i] = UINT_MAX;
  for (i = 0; i < op_count; ++i) clv_to_op[operations[i].parent_clv_index] = i;

  {
    size_t table_size = (size_t)sites * rate_cats * states;
    int    alloc_fail = 0;
    for (i = 0; i < op_count; ++i)
    {
      jtable[i] = (double *)malloc(table_size * sizeof(double));
      if (!jtable[i])
      {
        alloc_fail = 1;
        break;
      }
    }
    if (alloc_fail)
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC, "Cannot allocate memory");
      goto cleanup;
    }
  }

  /* default guards against a pathological, data-contradicts-model site
     (every state -INFINITY under every rate category) leaving best_rate
     read uninitialized below; rate category 0 is an arbitrary but
     well-defined fallback in that case */
  for (n = 0; n < sites; ++n)
  {
    best_total[n] = -INFINITY;
    best_rate[n]  = 0;
  }

  /* =====================================================================
     Main DP: for every rate category, run a full postorder max-sum pass,
     then combine across the virtual-root edge and keep, per site, which
     rate category gives the best total joint log-score (used below only
     to select which rate-category slice of jtable to report probs from).
     ===================================================================== */
  for (r = 0; r < rate_cats; ++r)
  {
    const double *freqs           = partition->frequencies[freqs_indices[r]];
    double        log_rate_weight = corax_joint_safe_log(partition->rate_weights[r]);

    for (op = 0; op < op_count; ++op)
    {
      const corax_operation_t *o  = &operations[op];
      unsigned int              c1 = o->child1_clv_index;
      unsigned int              c2 = o->child2_clv_index;
      unsigned int              c1_is_op = (c1 <= max_clv_index) && (clv_to_op[c1] != UINT_MAX);
      unsigned int              c2_is_op = (c2 <= max_clv_index) && (clv_to_op[c2] != UINT_MAX);
      unsigned int              op1 = c1_is_op ? clv_to_op[c1] : UINT_MAX;
      unsigned int              op2 = c2_is_op ? clv_to_op[c2] : UINT_MAX;

      const double *pmat1_base = partition->pmatrix[o->child1_matrix_index]
                                  + (size_t)r * states * states_padded;
      const double *pmat2_base = partition->pmatrix[o->child2_matrix_index]
                                  + (size_t)r * states * states_padded;

      for (n = 0; n < sites; ++n)
      {
        for (s = 0; s < states; ++s)
        {
          double best1 = -INFINITY;
          double best2 = -INFINITY;

          const double *row1 = pmat1_base + (size_t)s * states_padded;
          const double *row2 = pmat2_base + (size_t)s * states_padded;

          for (k = 0; k < states; ++k)
          {
            double child_val = c1_is_op
                ? jtable[op1][((size_t)n * rate_cats + r) * states + k]
                : corax_joint_leaf_logclv(partition, c1, n, r, k);
            double val = corax_joint_safe_log(row1[k]) + child_val;
            if (val > best1) best1 = val;
          }

          for (k = 0; k < states; ++k)
          {
            double child_val = c2_is_op
                ? jtable[op2][((size_t)n * rate_cats + r) * states + k]
                : corax_joint_leaf_logclv(partition, c2, n, r, k);
            double val = corax_joint_safe_log(row2[k]) + child_val;
            if (val > best2) best2 = val;
          }

          jtable[op][((size_t)n * rate_cats + r) * states + s] = best1 + best2;
        }
      }
    }

    /* ---- combine across the virtual-root edge for this rate category,
       purely to score how good this rate category is per site ---- */
    {
      unsigned int root_is_op  = (root_clv_index <= max_clv_index)
                                  && (clv_to_op[root_clv_index] != UINT_MAX);
      unsigned int other_is_op = (root_other_clv_index <= max_clv_index)
                                  && (clv_to_op[root_other_clv_index] != UINT_MAX);
      unsigned int root_op  = root_is_op ? clv_to_op[root_clv_index] : UINT_MAX;
      unsigned int other_op = other_is_op ? clv_to_op[root_other_clv_index] : UINT_MAX;

      const double *pmat_root_base =
          partition->pmatrix[root_matrix_index] + (size_t)r * states * states_padded;

      for (n = 0; n < sites; ++n)
      {
        double total_best = -INFINITY;

        for (s = 0; s < states; ++s)
        {
          double root_val = root_is_op
              ? jtable[root_op][((size_t)n * rate_cats + r) * states + s]
              : corax_joint_leaf_logclv(partition, root_clv_index, n, r, s);

          double        back_best = -INFINITY;
          const double *row = pmat_root_base + (size_t)s * states_padded;

          for (k = 0; k < states; ++k)
          {
            double other_val = other_is_op
                ? jtable[other_op][((size_t)n * rate_cats + r) * states + k]
                : corax_joint_leaf_logclv(partition, root_other_clv_index, n, r, k);
            double val = corax_joint_safe_log(row[k]) + other_val;
            if (val > back_best) back_best = val;
          }

          double total = corax_joint_safe_log(freqs[s]) + root_val + back_best;
          if (total > total_best) total_best = total;
        }

        double candidate = log_rate_weight + total_best;
        if (candidate > best_total[n])
        {
          best_total[n] = candidate;
          best_rate[n]  = r;
        }
      }
    }
  }

  /* =====================================================================
     Output probability vectors: softmax-normalize the log-joint-score
     table at each node, evaluated at that site's globally best rate
     category. This is a confidence-like companion score, NOT a
     marginal/empirical-Bayes posterior probability.
     ===================================================================== */
  for (op = 0; op < op_count; ++op)
  {
    for (n = 0; n < sites; ++n)
    {
      unsigned int rbest = best_rate[n];
      double       maxv  = -INFINITY;

      for (s = 0; s < states; ++s)
      {
        double v = jtable[op][((size_t)n * rate_cats + rbest) * states + s];
        if (v > maxv) maxv = v;
      }

      if (maxv == -INFINITY)
      {
        /* pathological: data contradicts the model at every state --
           fall back to a uniform vector rather than propagate NaN */
        for (s = 0; s < states; ++s)
          node_probs[op][(size_t)n * states + s] = 1.0 / states;
        continue;
      }

      double sum = 0.;
      for (s = 0; s < states; ++s)
      {
        double v = jtable[op][((size_t)n * rate_cats + rbest) * states + s];
        double e = exp(v - maxv);
        node_probs[op][(size_t)n * states + s] = e;
        sum += e;
      }
      for (s = 0; s < states; ++s) node_probs[op][(size_t)n * states + s] /= sum;
    }
  }

  retval = CORAX_SUCCESS;

cleanup:
  if (jtable)
    for (i = 0; i < op_count; ++i) free(jtable[i]);
  free(jtable);
  free(clv_to_op);
  free(best_total);
  free(best_rate);

  return retval;
}
