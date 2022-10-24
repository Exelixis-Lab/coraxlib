#include <corax/corax.h>

static void fill_parent_scaler(unsigned int        scaler_size,
                               unsigned int *      parent_scaler,
                               const unsigned int *left_scaler,
                               const unsigned int *right_scaler)
{
  unsigned int i;

  if (!left_scaler && !right_scaler)
    memset(parent_scaler, 0, sizeof(unsigned int) * scaler_size);
  else if (left_scaler && right_scaler)
  {
    memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * scaler_size);
    for (i = 0; i < scaler_size; ++i) parent_scaler[i] += right_scaler[i];
  }
  else
  {
    if (left_scaler)
      memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * scaler_size);
    else
      memcpy(parent_scaler, right_scaler, sizeof(unsigned int) * scaler_size);
  }
}

CORAX_EXPORT void corax_core_update_clv_ii_sve(unsigned int states,
                                                unsigned int sites,
                                                unsigned int rate_cats,
                                                double * parent_clv,
                                                unsigned int * parent_scaler,
                                                const double * left_clv,
                                                const double * right_clv,
                                                const double * left_matrix,
                                                const double * right_matrix,
                                                const unsigned int * left_scaler,
                                                const unsigned int * right_scaler,
                                                unsigned int attrib)
{
  unsigned int i,j,k,n;

  unsigned int scale_mode;  /* 0 = none, 1 = per-site, 2 = per-rate */
  unsigned int site_scale = 0;
  unsigned int init_mask;

  const double * lmat;
  const double * rmat;

  unsigned int span = states * rate_cats;

//  unsigned int states_padded = (states+3) & 0xFFFFFFFC;
  unsigned int states_padded = states;
  unsigned int span_padded = states_padded * rate_cats;

  /* init scaling-related stuff */
  if (parent_scaler)
  {
    /* determine the scaling mode and init the vars accordingly */
    scale_mode = (attrib & CORAX_ATTRIB_RATE_SCALERS) ? 2 : 1;
    init_mask = (scale_mode == 1) ? 1 : 0;
    const unsigned int scaler_size = (scale_mode == 2) ? sites * rate_cats : sites;

    /* add up the scale vectors of the two children if available */
    fill_parent_scaler(scaler_size, parent_scaler, left_scaler, right_scaler);
  }
  else
  {
    /* scaling disabled / not required */
    scale_mode = init_mask = 0;
  }
  
  /* compute CLV */
  for (n = 0; n < sites; ++n)
  {
    lmat = left_matrix;
    rmat = right_matrix;
    site_scale = init_mask;

    for (k = 0; k < rate_cats; ++k)
    {
      unsigned int rate_scale = 1;
      for (i = 0; i < states; ++i)
      {
//        double terma = 0;
//        double termb = 0;
        svfloat64_t va = svdup_f64(0.);
        svfloat64_t vb = svdup_f64(0.);
        for (j = 0; j < states; j += svcntd())
        {
          svbool_t pg = svwhilelt_b64(j, states);
          svfloat64_t lm = svld1(pg, &lmat[j]);
          svfloat64_t rm = svld1(pg, &rmat[j]);
          svfloat64_t lc = svld1(pg, &left_clv[j]);
          svfloat64_t rc = svld1(pg, &right_clv[j]);

          va = svmla_f64_x(pg, va, lm, lc);
          vb = svmla_f64_x(pg, vb, rm, rc);

//          terma += lmat[j] * left_clv[j];
//          termb += rmat[j] * right_clv[j];
        }

        double terma = svaddv_f64(svptrue_b64(), va);
        double termb = svaddv_f64(svptrue_b64(), vb);
        parent_clv[i] = terma*termb;

        rate_scale &= (parent_clv[i] < CORAX_SCALE_THRESHOLD);

        lmat += states;
        rmat += states;
      }

      /* check if scaling is needed for the current rate category */
      if (scale_mode == 2)
      {
        /* PER-RATE SCALING: if *all* entries of the *rate* CLV were below
         * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
        if (rate_scale)
        {
          for (i = 0; i < states; ++i)
            parent_clv[i] *= CORAX_SCALE_FACTOR;
          parent_scaler[n*rate_cats + k] += 1;
        }
      }
      else
        site_scale = site_scale && rate_scale;

      parent_clv += states;
      left_clv   += states;
      right_clv  += states;
    }

    /* PER-SITE SCALING: if *all* entries of the *site* CLV were below
     * the threshold then scale (all) entries by PLL_SCALE_FACTOR */
    if (site_scale)
    {
      parent_clv -= span;
      for (i = 0; i < span; ++i)
        parent_clv[i] *= CORAX_SCALE_FACTOR;
      parent_clv += span;
      parent_scaler[n] += 1;
    }
  }
}
