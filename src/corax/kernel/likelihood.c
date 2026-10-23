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

static double root_loglikelihood_asc_bias(corax_partition_t  *partition,
                                          unsigned int        sites,
                                          const double       *clv,
                                          unsigned int       *scaler,
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
  double       *rate_weights    = partition->rate_weights;

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
corax_compute_root_loglikelihood(corax_partition_t  *partition,
                                 unsigned int        clv_index,
                                 int                 scaler_index,
                                 const unsigned int *freqs_indices,
                                 double             *persite_lnl)
{
  return corax_compute_root_loglikelihood_sitecat(
      partition, clv_index, scaler_index, freqs_indices, persite_lnl, NULL);
}

CORAX_EXPORT double
corax_compute_root_loglikelihood_sitecat(corax_partition_t  *partition,
                                         unsigned int        clv_index,
                                         int                 scaler_index,
                                         const unsigned int *freqs_indices,
                                         double             *persite_lnl,
                                         double             *sitecat_lh)
{
  double        logl = 0;
  unsigned int *scaler;
  unsigned int  identifiers;

  if (sitecat_lh)
  {
    corax_set_error(
        CORAX_ERROR_NOT_IMPLEMENTED,
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
  double       *rate_weights    = partition->rate_weights;

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
                                          double             *persite_lnl,
                                          double             *sitecat_lh)
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

static double edge_loglikelihood_asc_bias_ii(corax_partition_t  *partition,
                                             const double       *clvp,
                                             unsigned int       *parent_scaler,
                                             unsigned int        parent_sites,
                                             const double       *clvc,
                                             unsigned int       *child_scaler,
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
  double       *rate_weights    = partition->rate_weights;

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

static double edge_loglikelihood(corax_partition_t  *partition,
                                 unsigned int        parent_clv_index,
                                 int                 parent_scaler_index,
                                 unsigned int        child_clv_index,
                                 int                 child_scaler_index,
                                 unsigned int        matrix_index,
                                 const unsigned int *freqs_indices,
                                 double             *persite_lnl,
                                 double             *sitecat_lh)
{
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
                                       sitecat_lh,
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
                                         double             *persite_lnl,
                                         double             *sitecat_lh)
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
corax_compute_edge_loglikelihood(corax_partition_t  *partition,
                                 unsigned int        parent_clv_index,
                                 int                 parent_scaler_index,
                                 unsigned int        child_clv_index,
                                 int                 child_scaler_index,
                                 unsigned int        matrix_index,
                                 const unsigned int *freqs_indices,
                                 double             *persite_lnl)
{
  return corax_compute_edge_loglikelihood_sitecat(partition,
                                                  parent_clv_index,
                                                  parent_scaler_index,
                                                  child_clv_index,
                                                  child_scaler_index,
                                                  matrix_index,
                                                  freqs_indices,
                                                  persite_lnl,
                                                  NULL);
}

CORAX_EXPORT double
corax_compute_edge_loglikelihood_sitecat(corax_partition_t *partition,
                                         unsigned int       parent_clv_index,
                                         int                parent_scaler_index,
                                         unsigned int       child_clv_index,
                                         int                child_scaler_index,
                                         unsigned int       matrix_index,
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
corax_compute_node_ancestral_extbuf(corax_partition_t  *partition,
                                    unsigned int        node_clv_index,
                                    int                 node_scaler_index,
                                    unsigned int        other_clv_index,
                                    int                 other_scaler_index,
                                    unsigned int        pmatrix_index,
                                    const unsigned int *freqs_indices,
                                    double             *ancestral,
                                    double             *temp_clv,
                                    unsigned int       *temp_scaler,
                                    double             *ident_pmat)
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
  unsigned int site_scalings;
  int          per_rate_scaling =
      (partition->attributes & CORAX_ATTRIB_RATE_SCALERS) ? 1 : 0;
  unsigned int *rate_scalings = NULL;
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
                                              double             *ancestral)
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

/* =============================================================================
   Joint (Pupko et al. 1998) ancestral state reconstruction.

   For every internal node in `operations`, computes a probability vector
   (softmax of the log-scores that the node's PARENT compared when it recorded
   the backpointer to this node) such that the arg max over states at every
   node, taken together, forms a single globally most likely joint assignment
   of states to all internal nodes (and, per site, to the rate category).

   Rate heterogeneity: one rate category is shared by the whole tree for a
   given site, and summing it out does not fold into a max-sum DP. Instead the
   rate is treated as one more jointly maximized hidden variable, i.e. for
   every site we report argmax_{states,r} w_r * P_r(data, states).
   This is a joint MAP over (rate category, ancestral states), NOT a joint MAP
   over ancestral states with the rate categories summed out.

   Implementation notes on performance:

     * Sites are independent, so the alignment is processed in small BLOCKS of
       sites. For each rate category the value-only DP of the whole block is
       evaluated and kept (val_all); each site's best rate is then picked from
       the root-edge totals and its traceback reads the tables of that rate.
       This needs R DP runs and no backpointer tables: the argmax is recomputed
       in the top-down traceback, where it costs O(states) per node and site.
       Working memory is O(block * ops * states).

     * The log-pmatrices are stored transposed (child state major), so the
       max-plus step is an element-wise max over the parent-state axis, which
       the compiler vectorizes.

     * A tip with a single possible state is simply one row of the log-pmatrix,
       a fully ambiguous tip (gap) uses a precomputed row maximum, and only
       partially ambiguous tips combine a few rows. A tip therefore costs
       O(states) and about half of all children in a binary tree are tips.

     * The traceback walks a precomputed flat list of (parent op, child op)
       edges, so it has no tip/inner branch, and the softmax shares its max
       with the argmax. The hot loops are instantiated for the usual
       (states, states_padded) combinations so the compiler sees constant trip
       counts.

   `root_scaler_index` / `root_other_scaler_index` are accepted for interface
   symmetry but unused: everything is computed in log-space directly from the
   tip data, so there is nothing to rescale.
    ^-- I'm unsure if this is correct --^
   =============================================================================
 */

#if defined(__GNUC__) || defined(__clang__)
#define CORAX_JOINT_FORCE_INLINE static inline __attribute__((always_inline))
#else
#define CORAX_JOINT_FORCE_INLINE static inline
#endif
 
/* possible-states bitmask of a tip at `site` */
static inline corax_state_t corax_joint_tip_mask(
    const corax_partition_t *partition, unsigned int clv_index, unsigned int site)
{
  if (!(partition->attributes & CORAX_ATTRIB_PATTERN_TIP))
  {
    /* no tip compression -> derive the mask from the (rate independent)
       0/1 tip clv */
    const double *clv = partition->clv[clv_index]
                        + (size_t)site * partition->rate_cats
                              * partition->states_padded;
    corax_state_t m = 0;
    for (unsigned int k = 0; k < partition->states; ++k)
    {
      assert(clv[k] == 0.0 || clv[k] == 1.0);
      if (clv[k] != 0.0) m |= (corax_state_t)1 << k;
    }
    return m;
  }
 
  unsigned char charcode = partition->tipchars[clv_index][site];
 
  /* For 4 states (DNA), corax's own specialized kernel
     corax_core_update_clv_ti_4x4() uses tipchars[][site] directly
     as the corax_state_t bitmask. For proteins, use it as an index */
  return (partition->states == 4) ? (corax_state_t)charcode
                                  : partition->tipmap[charcode];
}
 
/* Index of the first maximum of `v` (ties -> lowest index); *maxval is
   -INFINITY and the index 0 if every entry is -INFINITY. */
static inline unsigned int
corax_joint_argmax(const double *v, unsigned int states, double *maxval)
{
  double       best = -INFINITY;
  unsigned int arg  = 0;
  for (unsigned int k = 0; k < states; ++k)
  {
    int gt = v[k] > best;
    best   = gt ? v[k] : best;
    arg    = gt ? k : arg;
  }
  *maxval = best;
  return arg;
}
 
/* Softmax-normalize `score` (length `states`) into `out` and return its
   arg max. `out` is a uniform vector if every entry is -INFINITY (data
   contradicts the model at every state). */
static inline unsigned int
corax_joint_softmax_argmax(const double *score, unsigned int states, double *out)
{
  double       maxv;
  unsigned int arg = corax_joint_argmax(score, states, &maxv);
 
  if (maxv == -INFINITY)
  {
    for (unsigned int s = 0; s < states; ++s) out[s] = 1.0 / states;
    return arg;
  }
 
  double sum = 0.;
  for (unsigned int s = 0; s < states; ++s)
  {
    /* the maximum contributes exp(0) == 1 exactly: skip the libm call */
    double e = (score[s] == maxv) ? 1.0 : exp(score[s] - maxv);
    out[s]   = e;
    sum += e;
  }
  for (unsigned int s = 0; s < states; ++s) out[s] /= sum;
 
  return arg;
}
 
/* lt[r][k][s] = log P_r[s][k] (transposed, so that the DP max-plus step is an
   element-wise max over the parent state s) and rowmax[r][s] = max_k lt[r][k][s].
   Non-positive entries (numerical noise) count as probability 0, i.e. -inf.
   Padding lanes (s >= states) are -inf. */
static void corax_joint_build_lt(const corax_partition_t *partition,
                                 const double            *pmat,
                                 double                  *lt,
                                 double                  *rowmax)
{
  const unsigned int states = partition->states;
  const unsigned int sp     = partition->states_padded;
 
  for (unsigned int r = 0; r < partition->rate_cats; ++r)
  {
    double       *l = lt + (size_t)r * states * sp;
    double       *m = rowmax + (size_t)r * sp;
    const double *p = pmat + (size_t)r * states * sp;
 
    for (unsigned int s = 0; s < sp; ++s) m[s] = -INFINITY;
    for (size_t i = 0; i < (size_t)states * sp; ++i) l[i] = -INFINITY;
 
    for (unsigned int s = 0; s < states; ++s)
      for (unsigned int k = 0; k < states; ++k)
      {
        double pr              = p[(size_t)s * sp + k];
        double lg              = pr > 0.0 ? log(pr) : -INFINITY;
        l[(size_t)k * sp + s]  = lg;
        if (lg > m[s]) m[s] = lg;
      }
  }
}
 
/* contrib[s] = max_k ( logP[s][k] + child[k] )   for every parent state s,
   with logP stored transposed: lt[k * sp + s].
 
   child_val != NULL: inner child, `child_val[k]` is its DP score. The result
                      is written to `tmp`.
   child_val == NULL: tip child with possible-state set `mask` (child[k] is 0
                      or -inf). The result is, where possible, a pointer into
                      the precomputed tables and `tmp` is only touched for
                      partially ambiguous tips.
 
   Returns a vector of length `sp` (entries >= states are meaningless). */
static inline const double *corax_joint_contrib(const double *lt,
                                                const double *rowmax,
                                                const double *child_val,
                                                corax_state_t mask,
                                                corax_state_t full,
                                                unsigned int  states,
                                                unsigned int  sp,
                                                const double *neg_inf,
                                                double       *tmp)
{
  if (child_val)
  {
    for (unsigned int s = 0; s < sp; ++s) tmp[s] = lt[s] + child_val[0];
    for (unsigned int k = 1; k < states; ++k)
      for (unsigned int s = 0; s < sp; ++s)
      {
        double t = lt[(size_t)k * sp + s] + child_val[k];
        tmp[s]   = t > tmp[s] ? t : tmp[s];
      }
    return tmp;
  }
 
  mask &= full;
  if (mask == full) return rowmax; /* gap / fully ambiguous */
  if (!mask) return neg_inf;       /* no state possible */
  if (!(mask & (mask - 1)))        /* exactly one state: a row of lt */
  {
    unsigned int k = 0;
    while (!((mask >> k) & 1u)) ++k;
    return lt + (size_t)k * sp;
  }
 
  /* partially ambiguous: max over the possible states */
  for (unsigned int s = 0; s < sp; ++s) tmp[s] = -INFINITY;
  for (unsigned int k = 0; k < states; ++k)
    if ((mask >> k) & 1u)
      for (unsigned int s = 0; s < sp; ++s)
      {
        double t = lt[(size_t)k * sp + s];
        tmp[s]   = t > tmp[s] ? t : tmp[s];
      }
  return tmp;
}
 
/* The whole computation. `states` and `sp` (states_padded) are parameters so
   that the dispatcher below can instantiate it with compile-time constants for
   the usual sizes; this lets the compiler unroll/vectorize the inner loops. */
CORAX_JOINT_FORCE_INLINE int
corax_joint_impl(corax_partition_t       *partition,
                 const corax_operation_t *operations,
                 unsigned int             op_count,
                 unsigned int             root_clv_index,
                 unsigned int             root_other_clv_index,
                 unsigned int             root_matrix_index,
                 const unsigned int      *freqs_indices,
                 double                 **node_probs,
                 const unsigned int       states,
                 const unsigned int       sp)
{
  const unsigned int sites     = partition->sites;
  const unsigned int rate_cats = partition->rate_cats;
  const size_t       S         = states;
  const size_t       lt_block  = (size_t)rate_cats * S * sp; /* one (op, child) */
 
  const corax_state_t full = (states >= 8 * sizeof(corax_state_t))
                                 ? ~(corax_state_t)0
                                 : (((corax_state_t)1 << states) - 1);
 
  unsigned int max_clv_index = root_clv_index > root_other_clv_index
                                   ? root_clv_index
                                   : root_other_clv_index;
  for (unsigned int i = 0; i < op_count; ++i)
  {
    const corax_operation_t *o = &operations[i];
    if (o->parent_clv_index > max_clv_index) max_clv_index = o->parent_clv_index;
    if (o->child1_clv_index > max_clv_index) max_clv_index = o->child1_clv_index;
    if (o->child2_clv_index > max_clv_index) max_clv_index = o->child2_clv_index;
  }
 
  /* Block of sites processed together: at most 64, and keep the DP tables
     around 8 MiB. */
  size_t B = ((size_t)8 << 20)
             / ((size_t)rate_cats * op_count * sp * sizeof(double));
  if (B > 64) B = 64;
  if (B > sites) B = sites;
  if (B < 1) B = 1;
 
  /* All working memory is carved out of one allocation: doubles first, then
     corax_state_t, unsigned int and unsigned short (decreasing alignment). */
  const size_t n_dbl = (size_t)op_count * 2 * (lt_block + (size_t)rate_cats * sp)
                       + lt_block + 2 * (size_t)rate_cats * sp + rate_cats
                       + (size_t)rate_cats * op_count * B * sp + B * sp + B
                       + 6 * (size_t)sp;
  const size_t n_msk = ((size_t)op_count * 2 + 2) * B;
  const size_t n_u32 = (size_t)max_clv_index + 1 + 9 * (size_t)op_count + B;
 
  void *arena = malloc(n_dbl * sizeof(double) + n_msk * sizeof(corax_state_t)
                       + n_u32 * sizeof(unsigned int)
                       + (size_t)op_count * sizeof(unsigned short));
  if (!arena)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Cannot allocate memory");
    return CORAX_FAILURE;
  }
 
  /* [op][j][r][k][sp] transposed log-pmatrices of child j of op */
  double *lt = (double *)arena;
  /* [op][j][r][sp] row maxima of the above */
  double *rowmax = lt + (size_t)op_count * 2 * lt_block;
  /* same for the root edge: [r][k][sp] and [r][sp] */
  double *lt_root     = rowmax + (size_t)op_count * 2 * rate_cats * sp;
  double *rowmax_root = lt_root + lt_block;
  double *logfreq     = rowmax_root + (size_t)rate_cats * sp; /* [r][sp] */
  double *logw        = logfreq + (size_t)rate_cats * sp;     /* [r] */
  /* [r][b][op][sp] value-only DP tables of the current block */
  double *val_all = logw + rate_cats;
  /* [b][sp]: 0/-inf vector of the root if it is a tip */
  double *rootleaf   = val_all + (size_t)rate_cats * op_count * B * sp;
  double *best_total = rootleaf + B * sp; /* [b] */
  /* sp-sized scratch vectors */
  double *tmp1       = best_total + B;
  double *tmp2       = tmp1 + sp;
  double *tmpo       = tmp2 + sp;
  double *neg_inf    = tmpo + sp;
  double *root_score = neg_inf + sp;
  double *scratch    = root_score + sp;
 
  /* [op][j][b] tip masks, and the masks of the two root-edge nodes [b] */
  corax_state_t *masks      = (corax_state_t *)(scratch + sp);
  corax_state_t *root_mask  = masks + (size_t)op_count * 2 * B;
  corax_state_t *other_mask = root_mask + B;
 
  /* clv index -> op index (UINT_MAX: tip) */
  unsigned int *clv_to_op = (unsigned int *)(other_mask + B);
  /* [op*2+j]: op index of child j (UINT_MAX: tip) */
  unsigned int *child_op = clv_to_op + (size_t)max_clv_index + 1;
  /* tip children: slot (op*2+j) and clv index */
  unsigned int *tip_slot = child_op + 2 * (size_t)op_count;
  unsigned int *tip_clv  = tip_slot + 2 * (size_t)op_count;
  /* inner parent->child edges, parents in descending op order */
  unsigned int *edge_par   = tip_clv + 2 * (size_t)op_count;
  unsigned int *edge_child = edge_par + op_count;
  unsigned int *edge_slot  = edge_child + op_count;
  unsigned int *best_rate  = edge_slot + op_count; /* [b] */
  /* [op]: traceback state */
  unsigned short *op_state = (unsigned short *)(best_rate + B);
  assert((size_t)((char *)(op_state + op_count) - (char *)arena)
         == n_dbl * sizeof(double) + n_msk * sizeof(corax_state_t)
                + n_u32 * sizeof(unsigned int)
                + (size_t)op_count * sizeof(unsigned short));
 
  for (unsigned int s = 0; s < sp; ++s) neg_inf[s] = -INFINITY;
 
  /* Anything that is not the parent of an op is read as a tip */
  for (unsigned int i = 0; i <= max_clv_index; ++i) clv_to_op[i] = UINT_MAX;
  for (unsigned int i = 0; i < op_count; ++i)
    clv_to_op[operations[i].parent_clv_index] = i;
 
  unsigned int tip_count = 0;
  unsigned int edge_count = 0;
  for (unsigned int op = 0; op < op_count; ++op)
    for (unsigned int j = 0; j < 2; ++j)
    {
      const unsigned int clv =
          j ? operations[op].child2_clv_index : operations[op].child1_clv_index;
      child_op[2 * op + j] = clv_to_op[clv];
      if (clv_to_op[clv] == UINT_MAX)
      {
        tip_slot[tip_count] = 2 * op + j;
        tip_clv[tip_count]  = clv;
        ++tip_count;
      }
    }
  for (unsigned int op = op_count; op-- > 0;)
    for (unsigned int j = 0; j < 2; ++j)
      if (child_op[2 * op + j] != UINT_MAX)
      {
        edge_par[edge_count]   = op;
        edge_child[edge_count] = child_op[2 * op + j];
        edge_slot[edge_count]  = 2 * op + j;
        ++edge_count;
      }
 
  const unsigned int root_op     = clv_to_op[root_clv_index];
  const unsigned int other_op    = clv_to_op[root_other_clv_index];
  const int          root_is_op  = (root_op != UINT_MAX);
  const int          other_is_op = (other_op != UINT_MAX);
 
  for (unsigned int op = 0; op < op_count; ++op)
    for (unsigned int j = 0; j < 2; ++j)
      corax_joint_build_lt(partition,
                           partition->pmatrix[j ? operations[op].child2_matrix_index
                                                : operations[op].child1_matrix_index],
                           lt + ((size_t)op * 2 + j) * lt_block,
                           rowmax + ((size_t)op * 2 + j) * rate_cats * sp);
  corax_joint_build_lt(
      partition, partition->pmatrix[root_matrix_index], lt_root, rowmax_root);
 
  for (unsigned int r = 0; r < rate_cats; ++r)
  {
    const double *freqs = partition->frequencies[freqs_indices[r]];
    for (unsigned int s = 0; s < sp; ++s)
      logfreq[(size_t)r * sp + s] = (s < states) ? log(freqs[s]) : -INFINITY;
    logw[r] = log(partition->rate_weights[r]);
  }
 
  /* =====================================================================
     Process the alignment in blocks of B sites.
     ===================================================================== */
  for (unsigned int n0 = 0; n0 < sites; n0 += (unsigned int)B)
  {
    const unsigned int nb = (sites - n0 < B) ? (sites - n0) : (unsigned int)B;
 
    /* ---- tip masks (rate independent) + reset the rate selection ---- */
    for (unsigned int t = 0; t < tip_count; ++t)
      for (unsigned int b = 0; b < nb; ++b)
        masks[(size_t)tip_slot[t] * B + b] =
            corax_joint_tip_mask(partition, tip_clv[t], n0 + b);
 
    for (unsigned int b = 0; b < nb; ++b)
    {
      root_mask[b]  = root_is_op ? 0
                                 : corax_joint_tip_mask(
                                       partition, root_clv_index, n0 + b);
      other_mask[b] = other_is_op ? 0
                                  : corax_joint_tip_mask(
                                        partition, root_other_clv_index, n0 + b);
 
      for (unsigned int s = 0; s < sp; ++s)
        rootleaf[(size_t)b * sp + s] =
            (s < states && ((root_mask[b] >> s) & 1u)) ? 0.0 : -INFINITY;
 
      best_total[b] = -INFINITY;
      best_rate[b]  = 0;
    }
 
    /* ---- value-only DP + root-edge total, one rate category at a time ---- */
    for (unsigned int r = 0; r < rate_cats; ++r)
    {
      for (unsigned int op = 0; op < op_count; ++op)
      {
        const unsigned int c1op = child_op[2 * op];
        const unsigned int c2op = child_op[2 * op + 1];
        const size_t       o1   = (size_t)op * 2;
        const double *lt1 = lt + o1 * lt_block + (size_t)r * S * sp;
        const double *lt2 = lt + (o1 + 1) * lt_block + (size_t)r * S * sp;
        const double *rm1 = rowmax + (o1 * rate_cats + r) * sp;
        const double *rm2 = rowmax + ((o1 + 1) * rate_cats + r) * sp;
 
        for (unsigned int b = 0; b < nb; ++b)
        {
          double *vb = val_all + ((size_t)r * B + b) * op_count * sp;
 
          const double *x1 = corax_joint_contrib(
              lt1,
              rm1,
              c1op == UINT_MAX ? NULL : vb + (size_t)c1op * sp,
              masks[o1 * B + b],
              full,
              states,
              sp,
              neg_inf,
              tmp1);
          const double *x2 = corax_joint_contrib(
              lt2,
              rm2,
              c2op == UINT_MAX ? NULL : vb + (size_t)c2op * sp,
              masks[(o1 + 1) * B + b],
              full,
              states,
              sp,
              neg_inf,
              tmp2);
 
          double *out = vb + (size_t)op * sp;
          for (unsigned int s = 0; s < sp; ++s) out[s] = x1[s] + x2[s];
        }
      }
 
      /* root-edge combination for this rate category */
      const double *lf = logfreq + (size_t)r * sp;
      for (unsigned int b = 0; b < nb; ++b)
      {
        const double *vb = val_all + ((size_t)r * B + b) * op_count * sp;
 
        const double *oc = corax_joint_contrib(
            lt_root + (size_t)r * S * sp,
            rowmax_root + (size_t)r * sp,
            other_is_op ? vb + (size_t)other_op * sp : NULL,
            other_mask[b],
            full,
            states,
            sp,
            neg_inf,
            tmpo);
        const double *rv =
            root_is_op ? vb + (size_t)root_op * sp : rootleaf + (size_t)b * sp;
 
        double total_best = -INFINITY;
        for (unsigned int s = 0; s < states; ++s)
        {
          double v = lf[s] + rv[s] + oc[s];
          if (v > total_best) total_best = v;
        }
 
        double candidate = logw[r] + total_best;
        if (r == 0 || candidate > best_total[b])
        {
          best_total[b] = candidate;
          best_rate[b]  = r;
        }
      }
    }
 
    /* ---- per site: root combination at its best rate, then top-down
            traceback emitting each node's probability vector ---- */
    for (unsigned int b = 0; b < nb; ++b)
    {
      const unsigned int n   = n0 + b;
      const unsigned int rr  = best_rate[b];
      const double      *lf  = logfreq + (size_t)rr * sp;
      const double      *ltr = lt_root + (size_t)rr * S * sp;
      const double      *vb  = val_all + ((size_t)rr * B + b) * op_count * sp;
 
      const double *other_val = other_is_op ? vb + (size_t)other_op * sp : NULL;
      const double *oc        = corax_joint_contrib(ltr,
                                             rowmax_root + (size_t)rr * sp,
                                             other_val,
                                             other_mask[b],
                                             full,
                                             states,
                                             sp,
                                             neg_inf,
                                             tmpo);
      const double *rv =
          root_is_op ? vb + (size_t)root_op * sp : rootleaf + (size_t)b * sp;
 
      for (unsigned int s = 0; s < states; ++s)
        root_score[s] = lf[s] + rv[s] + oc[s];
 
      /* best root state (also needed for the other side if the root is a tip) */
      double       best_val;
      unsigned int best_s;
      if (root_is_op)
      {
        best_s = corax_joint_softmax_argmax(
            root_score, states, node_probs[root_op] + (size_t)n * states);
        best_val          = root_score[best_s];
        op_state[root_op] = (unsigned short)best_s;
      }
      else
        best_s = corax_joint_argmax(root_score, states, &best_val);
 
      if (other_is_op)
      {
        for (unsigned int k = 0; k < states; ++k)
          scratch[k] = ltr[(size_t)k * sp + best_s] + other_val[k];
 
        unsigned int arg = corax_joint_softmax_argmax(
            scratch, states, node_probs[other_op] + (size_t)n * states);
 
        /* if no root state has a finite score (data contradicts the model),
           keep state 0 for the other side, as the previous implementation did */
        op_state[other_op] = (unsigned short)((best_val == -INFINITY) ? 0 : arg);
      }
 
      /* Edges are sorted by descending parent op index and `operations` is in
         postorder (children before parents), so an op's own state is always
         set -- by the root seeding above or by the edge into it -- before an
         edge out of it is processed. */
      for (unsigned int e = 0; e < edge_count; ++e)
      {
        const unsigned int xs = op_state[edge_par[e]];
        const unsigned int cp = edge_child[e];
        const double *ltj = lt + (size_t)edge_slot[e] * lt_block
                            + (size_t)rr * S * sp;
        const double *cv = vb + (size_t)cp * sp;
 
        for (unsigned int k = 0; k < states; ++k)
          scratch[k] = ltj[(size_t)k * sp + xs] + cv[k];
 
        op_state[cp] = (unsigned short)corax_joint_softmax_argmax(
            scratch, states, node_probs[cp] + (size_t)n * states);
      }
    }
  }
 
  free(arena);
  return CORAX_SUCCESS;
}
 
CORAX_EXPORT int
corax_compute_joint_ancestral(corax_partition_t       *partition,
                              const corax_operation_t *operations,
                              unsigned int             op_count,
                              unsigned int             root_clv_index,
                              int                      root_scaler_index,
                              unsigned int             root_other_clv_index,
                              int                      root_other_scaler_index,
                              unsigned int             root_matrix_index,
                              const unsigned int      *freqs_indices,
                              double                 **node_probs)
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
 
  const unsigned int states = partition->states;
  const unsigned int sp     = partition->states_padded;
 
  if (states == 4 && sp == 4)
    return corax_joint_impl(partition, operations, op_count, root_clv_index,
                            root_other_clv_index, root_matrix_index,
                            freqs_indices, node_probs, 4, 4);
  if (states == 4 && sp == 8)
    return corax_joint_impl(partition, operations, op_count, root_clv_index,
                            root_other_clv_index, root_matrix_index,
                            freqs_indices, node_probs, 4, 8);
  if (states == 20 && sp == 20)
    return corax_joint_impl(partition, operations, op_count, root_clv_index,
                            root_other_clv_index, root_matrix_index,
                            freqs_indices, node_probs, 20, 20);
  if (states == 20 && sp == 24)
    return corax_joint_impl(partition, operations, op_count, root_clv_index,
                            root_other_clv_index, root_matrix_index,
                            freqs_indices, node_probs, 20, 24);
  return corax_joint_impl(partition, operations, op_count, root_clv_index,
                          root_other_clv_index, root_matrix_index,
                          freqs_indices, node_probs, states, sp);
}
