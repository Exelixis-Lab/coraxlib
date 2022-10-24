#include <corax/corax.h>
#include <limits.h>

CORAX_EXPORT
double corax_core_edge_loglikelihood_ii_sve(unsigned int         states,
                                            unsigned int         sites,
                                            unsigned int         rate_cats,
                                            const double *       parent_clv,
                                            const unsigned int * parent_scaler,
                                            const double *       child_clv,
                                            const unsigned int * child_scaler,
                                            const double *       pmatrix,
                                            const double *const *frequencies,
                                            const double *       rate_weights,
                                            const unsigned int *pattern_weights,
                                            const double *invar_proportion,
                                            const int *   invar_indices,
                                            const unsigned int *freqs_indices,
                                            double *            persite_lnl,
                                            unsigned int        attrib)
{
  unsigned int n, i, j, k;
  double       logl       = 0;
  double       prop_invar = 0;

  const double *clvp = parent_clv;
  const double *clvc = child_clv;
  const double *pmat;
  const double *freqs = NULL;

  double terma, terma_r, terminv;
  double site_lk, inv_site_lk;

//  unsigned int states_padded = (states + 3) & 0xFFFFFFFC;
  unsigned int states_padded = states;


  size_t displacement = (states_padded - states) * (states_padded);

  /* scaling stuff */
  unsigned int  site_scalings;
  unsigned int *rate_scalings    = NULL;
  int           per_rate_scaling = (attrib & CORAX_ATTRIB_RATE_SCALERS) ? 1 : 0;

  /* powers of scale threshold for undoing the scaling */
  double scale_minlh[CORAX_SCALE_RATE_MAXDIFF];
  if (per_rate_scaling || invar_proportion)
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
                      "Cannot allocate space for precomputation.");
      return -INFINITY;
    }
  }

  for (n = 0; n < sites; ++n)
  {
    pmat    = pmatrix;
    terma   = 0;
    terminv = 0;

    if (per_rate_scaling)
    {
      /* compute minimum per-rate scaler -> common per-site scaler */
      site_scalings = UINT_MAX;
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] =
            (parent_scaler) ? parent_scaler[n * rate_cats + i] : 0;
        rate_scalings[i] +=
            (child_scaler) ? child_scaler[n * rate_cats + i] : 0;
        if (rate_scalings[i] < site_scalings) site_scalings = rate_scalings[i];
      }

      /* compute relative capped per-rate scalers */
      for (i = 0; i < rate_cats; ++i)
      {
        rate_scalings[i] = CORAX_MIN(rate_scalings[i] - site_scalings,
                                     CORAX_SCALE_RATE_MAXDIFF);
      }
    }
    else
    {
      /* count number of scaling factors to account for */
      site_scalings = (parent_scaler) ? parent_scaler[n] : 0;
      site_scalings += (child_scaler) ? child_scaler[n] : 0;
    }

    for (i = 0; i < rate_cats; ++i)
    {
      freqs   = frequencies[freqs_indices[i]];
      terma_r = 0;

      for (j = 0; j < states; ++j)
      {
        svfloat64_t vb = svdup_f64(0.);
//        double termb0 = 0;
//        for (k = 0; k < states; k += 1)
//          termb0 += pmat[k] * clvc[k];

        for (k = 0; k < states; k += svcntd()) 
        { 
          svbool_t pg = svwhilelt_b64(k, states);
          svfloat64_t vpm = svld1(pg, &pmat[k]);
          svfloat64_t vcl = svld1(pg, &clvc[k]);

          vb = svmla_f64_x(pg, vb, vpm, vcl);
    
//          termb0 += pmat[k] * clvc[k]; 
        }

        double termb = svaddv_f64(svptrue_b64(), vb);
//        printf("termb: %.12lf     %.12lf\n", termb0, termb);

        terma_r += clvp[j] * freqs[j] * termb;
        pmat += states_padded;
      }

      /* apply per-rate scalers, if necessary */
      if (rate_scalings && rate_scalings[i] > 0)
      {
        terma_r *= scale_minlh[rate_scalings[i] - 1];
      }

      /* account for invariant sites */
      prop_invar = invar_proportion ? invar_proportion[freqs_indices[i]] : 0;
      if (prop_invar > 0)
      {
        terma += rate_weights[i] * terma_r * (1. - prop_invar);
        if (invar_indices[n] != -1)
        {
          freqs       = frequencies[freqs_indices[i]];
          inv_site_lk = freqs[invar_indices[n]];
          terminv += rate_weights[i] * inv_site_lk * prop_invar;
        }
      }
      else
      {
        terma += terma_r * rate_weights[i];
      }

      clvc += states_padded;
      clvp += states_padded;
      pmat -= displacement;
    }

    /* compute site log-likelihood and scale if necessary */
    if (site_scalings)
    {
      if (terminv > 0.)
      {
        /* IMPORTANT: undoing the scaling for non-variant likelihood term only!
         */
        unsigned int capped_scalings =
            CORAX_MIN(site_scalings, CORAX_SCALE_RATE_MAXDIFF);
        double scale_factor = scale_minlh[capped_scalings - 1];
        site_lk             = log(terma * scale_factor + terminv);
      }
      else
      {
        site_lk = log(terma);
        site_lk += site_scalings * log(CORAX_SCALE_THRESHOLD);
      }
    }
    else
    {
      site_lk = log(terma + terminv);
    }

    site_lk *= pattern_weights[n];

    /* store per-site log-likelihood */
    if (persite_lnl) persite_lnl[n] = site_lk;

    logl += site_lk;
  }

  if (rate_scalings) free(rate_scalings);

  return logl;
}


