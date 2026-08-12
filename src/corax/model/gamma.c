/*
    Copyright (C) 2015 Tomas Flouri, Diego Darriba, Alexandros Stamatakis

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

#include "corax/model/gamma.h"
#include "corax/corax.h"
#include "corax/core/common.h"

#define POINT_GAMMA(prob, alpha, beta)                                         \
  PointChi2(prob, 2.0 * (alpha)) / (2.0 * (beta))
#define ALPHA_MIN 0.02

static double IncompleteGamma(double x, double alpha, double ln_gamma_alpha)
{
  /* returns the incomplete gamma ratio I(x,alpha) where x is the upper
             limit of the integration and alpha is the shape parameter.
     returns (-1) if in error
     ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
     (1) series expansion     if (alpha>x || x<=1)
     (2) continued fraction   otherwise
     RATNEST FORTRAN by
     Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
     19: 285-287 (AS32)
  */
  int    i;
  double p = alpha, g = ln_gamma_alpha;
  double accurate = 1e-8, overflow = 1e30;
  double factor, gin = 0, rn = 0, a = 0, b = 0, an = 0, dif = 0, term = 0,
                 pn[6];

  if (x == 0) return (0);
  if (x < 0 || p <= 0) return (-1);

  factor = exp(p * log(x) - x - g);
  if (x > 1 && x >= p) goto l30;
  /* (1) series expansion */
  gin  = 1;
  term = 1;
  rn   = p;
l20:
  rn++;
  term *= x / rn;
  gin += term;

  if (term > accurate) goto l20;
  gin *= factor / p;
  goto l50;
l30:
  /* (2) continued fraction */
  a     = 1 - p;
  b     = a + x + 1;
  term  = 0;
  pn[0] = 1;
  pn[1] = x;
  pn[2] = x + 1;
  pn[3] = x * b;
  gin   = pn[2] / pn[3];
l32:
  a++;
  b += 2;
  term++;
  an = a * term;
  for (i = 0; i < 2; i++) pn[i + 4] = b * pn[i + 2] - an * pn[i];
  if (pn[5] == 0) goto l35;
  rn  = pn[4] / pn[5];
  dif = fabs(gin - rn);
  if (dif > accurate) goto l34;
  if (dif <= accurate * rn) goto l42;
l34:
  gin = rn;
l35:
  for (i = 0; i < 4; i++) pn[i] = pn[i + 2];
  if (fabs(pn[4]) < overflow) goto l32;

  for (i = 0; i < 4; i++) pn[i] /= overflow;

  goto l32;
l42:
  gin = 1 - factor * gin;

l50:
  return (gin);
}

static double LnGamma(double alpha)
{
  /* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
     Stirling's formula is used for the central polynomial part of the
     procedure. Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma
     function. Communications of the Association for Computing Machinery, 9:684
  */
  double x, f, z, result;

  x = alpha;
  f = 0.0;

  if (x < 7.0)
  {
    f = 1.0;
    z = alpha - 1.0;

    while ((z = z + 1.0) < 7.0) { f *= z; }
    x = z;

    assert(f != 0.0);

    f = -log(f);
  }

  z = 1 / (x * x);

  result =
      f + (x - 0.5) * log(x) - x + .918938533204673
      + (((-.000595238095238 * z + .000793650793651) * z - .002777777777778) * z
         + .083333333333333)
            / x;

  return result;
}

static double PointNormal(double prob)
{
  /* returns z so that Prob{x<z}=prob where x ~ N(0,1) and
     (1e-12)<prob<1-(1e-12) returns (-9999) if in error Odeh RE & Evans JO
     (1974) The percentage points of the normal distribution. Applied Statistics
     22: 96-97 (AS70)

     Newer methods:
       Wichura MJ (1988) Algorithm AS 241: the percentage points of the
         normal distribution.  37: 477-484.
       Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
         points of the normal distribution.  26: 118-121.

  */
  double a0 = -.322232431088, a1 = -1, a2 = -.342242088547,
         a3 = -.0204231210245;
  double a4 = -.453642210148e-4, b0 = .0993484626060, b1 = .588581570495;
  double b2 = .531103462366, b3 = .103537752850, b4 = .0038560700634;
  double y, z = 0, p = prob, p1;

  p1 = (p < 0.5 ? p : 1 - p);
  if (p1 < 1e-20) return (-9999);

  y = sqrt(log(1 / (p1 * p1)));
  z = y
      + ((((y * a4 + a3) * y + a2) * y + a1) * y + a0)
            / ((((y * b4 + b3) * y + b2) * y + b1) * y + b0);
  return (p < 0.5 ? -z : z);
}

static double PointChi2(double prob, double v)
{
  /* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
     returns -1 if in error.   0.000002<prob<0.999998
     RATNEST FORTRAN by
         Best DJ & Roberts DE (1975) The percentage points of the
         Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
     Converted into C by Ziheng Yang, Oct. 1993.
  */
  double e = .5e-6, aa = .6931471805, p = prob, g;
  double xx, c, ch, a = 0, q = 0, p1 = 0, p2 = 0, t = 0, x = 0, b = 0, s1, s2,
                    s3, s4, s5, s6;

  if (p < .000002 || p > .999998 || v <= 0) return (-1);

  g = LnGamma(v / 2);

  xx = v / 2;
  c  = xx - 1;
  if (v >= -1.24 * log(p)) goto l1;

  ch = pow((p * xx * exp(g + xx * aa)), 1 / xx);
  if (ch - e < 0) return (ch);
  goto l4;
l1:
  if (v > .32) goto l3;
  ch = 0.4;
  a  = log(1 - p);
l2:
  q  = ch;
  p1 = 1 + ch * (4.67 + ch);
  p2 = ch * (6.73 + ch * (6.66 + ch));
  t  = -0.5 + (4.67 + 2 * ch) / p1 - (6.73 + ch * (13.32 + 3 * ch)) / p2;
  ch -= (1 - exp(a + g + .5 * ch + c * aa) * p2 / p1) / t;
  if (fabs(q / ch - 1) - .01 <= 0)
    goto l4;
  else
    goto l2;

l3:
  x  = PointNormal(p);
  p1 = 0.222222 / v;
  ch = v * pow((x * sqrt(p1) + 1 - p1), 3.0);
  if (ch > 2.2 * v + 6) ch = -2 * (log(1 - p) - c * log(.5 * ch) + g);
l4:
  q  = ch;
  p1 = .5 * ch;
  if ((t = IncompleteGamma(p1, xx, g)) < 0.0)
  {
    printf("IncompleteGamma \n");
    return (-1);
  }

  p2 = p - t;
  t  = p2 * exp(xx * aa + g + p1 - c * log(ch));
  b  = t / ch;
  a  = 0.5 * t - b * c;

  s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) / 420;
  s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) / 2520;
  s3 = (210 + a * (462 + a * (707 + 932 * a))) / 2520;
  s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) / 5040;
  s5 = (84 + 264 * a + c * (175 + 606 * a)) / 2520;
  s6 = (120 + c * (346 + 127 * c)) / 5040;
  ch += t
        * (1 + 0.5 * t * s1
           - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));
  if (fabs(q / ch - 1) > e) goto l4;

  return (ch);
}

CORAX_EXPORT int corax_compute_gamma_rates(double       alpha,
                                          unsigned int categories,
                                          double      *output_rates,
                                          int          rates_mode)
{
  unsigned int i;

  double factor = alpha / alpha * categories, lnga1, alfa = alpha, beta = alpha,
         *gammaProbs;

  /* Note that ALPHA_MIN setting is somewhat critical due to   */
  /* numerical instability caused by very small rate[0] values */
  /* induced by low alpha values around 0.01 */

  if (alpha < ALPHA_MIN || categories < 1)
  {
    corax_set_error(
        CORAX_ERROR_INVALID_PARAM, "Invalid alpha value (%f)", alpha);
    return CORAX_FAILURE;
  }

  if (categories == 1) { output_rates[0] = 1.0; }
  else if (rates_mode == CORAX_GAMMA_RATES_MEDIAN)
  {
    double middle = 1.0 / (2.0 * categories), t = 0.0;

    for (i = 0; i < categories; i++)
      output_rates[i] = POINT_GAMMA((double)(i * 2 + 1) * middle, alfa, beta);

    for (i = 0; i < categories; i++) t += output_rates[i];
    for (i = 0; i < categories; i++) output_rates[i] *= factor / t;
  }
  else if (rates_mode == CORAX_GAMMA_RATES_MEAN)
  {
    gammaProbs = (double *)malloc(categories * sizeof(double));

    lnga1 = LnGamma(alfa + 1);

    for (i = 0; i < categories - 1; i++)
      gammaProbs[i] = POINT_GAMMA((i + 1.0) / categories, alfa, beta);

    for (i = 0; i < categories - 1; i++)
      gammaProbs[i] = IncompleteGamma(gammaProbs[i] * beta, alfa + 1, lnga1);

    output_rates[0] = gammaProbs[0] * factor;

    output_rates[categories - 1] = (1 - gammaProbs[categories - 2]) * factor;

    for (i = 1; i < categories - 1; i++)
      output_rates[i] = (gammaProbs[i] - gammaProbs[i - 1]) * factor;

    free(gammaProbs);
  }
  else
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM,
                    "Invalid GAMMA discretization mode (%d)",
                    rates_mode);
    return CORAX_FAILURE;
  }

  return CORAX_SUCCESS;
}

CORAX_EXPORT int
lloyd_max_init_cats(double alpha, unsigned int categories, double *output_rates)
{
  unsigned int i;

  double alfa = (alpha + 1.0) / 3.0, beta = alpha / 3.0;

  /* Note that ALPHA_MIN setting is somewhat critical due to   */
  /* numerical instability caused by very small rate[0] values */
  /* induced by low alpha values around 0.01 */

  if (alpha < ALPHA_MIN || categories < 1)
  {
    corax_set_error(
        CORAX_ERROR_INVALID_PARAM, "Invalid alpha value (%f)", alpha);
    return CORAX_FAILURE;
  }

  if (categories == 1) { output_rates[0] = 1.0; }
  double ratio = 1.0 / (categories + 1);

  for (i = 0; i < categories; i++)
    output_rates[i] = POINT_GAMMA((double)(i + 1) * ratio, alfa, beta);

  return CORAX_SUCCESS;
}

CORAX_EXPORT int corax_compute_gamma_cats_lloyd_max(double       alpha,
                                                      unsigned int categories,
                                                      double      *output_rates,
                                                      double *output_weights)
{
  /*
     Lloyd-Max discretization of a Gamma distribution G(alpha, beta).

     This function has the same argument format as PAML's original
     DiscreteGamma() function, but it produces a Lloyd-Max discretization
     instead of forcing all bins to have equal probability.

     Inputs:
        alpha         Shape parameter of the Gamma distribution.
        categories    Number of discrete categories/bins.

     Outputs:
        freqK[i]   Probability mass in Lloyd-Max bin i.
                   These probabilities are generally not equal.
        rK[i]      Representative rate for bin i.  In Lloyd-Max, this is
                   the conditional mean/centroid of the Gamma distribution
                   within that bin.

     External functions/macros assumed to be available elsewhere:
        QuantileGamma(prob, alpha, beta)
        IncompleteGamma(x, alpha, ln_gamma_alpha)
        malloc(), free(), fabs(), lgamma()

  */

  unsigned int i, iter;

  /* Maximum number of Lloyd-Max iterations */
  const unsigned int max_iter = 1000;

  /* Relative convergence tolerance for changes in rK[] */
  const double tol = 1e-4;

  double beta = alpha;

  double mean  = alpha / beta;
  double lnga  = lgamma(alpha);
  double lnga1 = lgamma(alpha + 1.0);
  double prev_rate, denom, change;
  double Fa0, Fa1, Ga0, Ga1, mass,
      moment_mass; /* Fa0 and Fa1 are Gamma CDF values with shape alpha at the
                      lower and upper bin boundaries.  Ga0 and Ga1 are
                      corresponding CDF values for a Gamma distribution with
                      shape alpha + 1.  The alpha + 1 terms are used to compute
                      the first moment inside a bin. */
  double *bound;   /* Boundary array of length K + 1.  Do not explicitly store
                      infinity in bound[K].   For the final bin, the upper CDF
                      values are set directly to 1.0. */

  if (categories <= 0 || alpha <= 0 || beta <= 0)
    return (-1); /* Basic validity check */

  /* temporary storage for bin boundaries */
  bound = (double *)malloc((categories + 1) * sizeof(double));

  if (bound == 0)
  {
    corax_set_error(
        CORAX_ERROR_MEM_ALLOC,
        "Failed to allocate temporary storage for gamma weight optimization ");
    return CORAX_FAILURE;
  }

  /*
     Initialization step.

     Lloyd-Max needs an initial ordered set of representative values.  Here we
     use the medians of K equal-probability Gamma intervals:

        p_i = (2i + 1) / (2K)
        rK[i] = Gamma quantile at p_i

     These are not the final Lloyd-Max representatives.  They are just a simple,
     stable, ordered starting point for the iteration.
  */
  lloyd_max_init_cats(alpha, categories, output_rates);
  /* Main Lloyd-Max iteration. */
  for (iter = 0; iter < max_iter; iter++)
  {

    /*
       Step 1: compute decision boundaries from current representatives.

       For squared-error Lloyd-Max quantization, the boundary between two
       adjacent representatives is the midpoint between them.

       The Gamma distribution has support [0, infinity), so the first lower
       boundary is 0.  The final upper boundary is infinity, handled later
       by using CDF value 1.0 for the last bin.
    */
    bound[0] = 0.0;
    for (i = 1; i < categories; i++)
    {
      bound[i] = (output_rates[i - 1] - output_rates[i])
                 / (log(output_rates[i - 1]) - log(output_rates[i]));
    }

    /* Track the largest relative change in any representative during this
       iteration.  This is used as the convergence criterion. */
    bool converged = true;

    /*
       Step 2: update each bin probability and representative value.

       For bin i with lower boundary a and upper boundary b:

          freqK[i] = P(a <= X < b)
                   = F_alpha(b) - F_alpha(a)

       The conditional mean is:

          rK[i] = E[X | a <= X < b]

       For Gamma(alpha, beta), the partial first moment over [a,b] is:

          E[X ; a <= X < b]
            = (alpha / beta) * [F_{alpha+1}(b) - F_{alpha+1}(a)]

       Therefore:

          rK[i] = (alpha / beta)
                  * [F_{alpha+1}(b) - F_{alpha+1}(a)]
                  / [F_alpha(b) - F_alpha(a)]
    */
    for (i = 0; i < categories; i++)
    {

      /* Save previous representative so we can measure convergence after
         updating rK[i]. */
      prev_rate = output_rates[i];

      /* Lower-bound CDF values.

         For the first bin, the lower boundary is 0.  The Gamma CDF at 0 is
         0 for positive alpha, so both CDF values are exactly 0. */
      if (i == 0)
      {
        Fa0 = 0.0;
        Ga0 = 0.0;
      }
      else
      {
        /* IncompleteGamma() expects the argument beta*x for a Gamma
           distribution with rate beta. */
        Fa0 = IncompleteGamma(bound[i] * beta, alpha, lnga);
        Ga0 = IncompleteGamma(bound[i] * beta, alpha + 1.0, lnga1);
      }

      /* Upper-bound CDF values.

         For the final bin, the upper boundary is infinity.  The CDF at
         infinity is 1, so we use 1.0 directly rather than storing infinity
         or calling IncompleteGamma() with an infinite argument. */
      if (i == categories - 1)
      {
        Fa1 = 1.0;
        Ga1 = 1.0;
      }
      else
      {
        Fa1 = IncompleteGamma(bound[i + 1] * beta, alpha, lnga);
        Ga1 = IncompleteGamma(bound[i + 1] * beta, alpha + 1.0, lnga1);
      }

      /* Probability mass in the current Lloyd-Max bin. */
      mass = Fa1 - Fa0;

      /* Corresponding alpha+1 CDF difference, used for the partial first
         moment in this bin. */
      moment_mass = Ga1 - Ga0;

      /* Store the bin probability in the output frequency array. */
      output_weights[i] = mass;

      /* Update representative rate to the conditional mean in the bin.

         If mass is effectively zero, leave rK[i] unchanged to avoid
         division by zero.  This should rarely occur for ordinary K and
         well-behaved alpha/beta values. */
      if (mass > CORAX_MISC_EPSILON)
      {
        output_rates[i] = mean * moment_mass / mass;
      }

      /* Switch to abstol if the change is close to the machine eps */
      denom  = fabs(prev_rate);
      change = fabs(output_rates[i] - prev_rate);
      if (denom >= CORAX_MISC_EPSILON) { change /= denom; }

      /* Keep the largest relative change across all bins. */
      converged &= change < tol;
    }

    /* Stop if all representative rates changed by less than the tolerance. */
    if (converged) { break; }
  }

  if (iter == max_iter)
  {
    corax_set_error(CORAX_ERROR_OPT_CONVERGE,
                    "Failed to converge for Lloyd-Max quantization of gamma "
                    "distribution after %d iterations",
                    iter);
    return CORAX_FAILURE;
  }

  free(bound); /* Release temporary storage */

  return CORAX_SUCCESS;
}

double glaguerre(long m, double b, double x)
{ /* Generalized Laguerre polynomial computed recursively.
     For use by initgammacat */
  long i;
  double gln, glnm1, glnp1; /* L_n, L_(n-1), L_(n+1) */

  if (m == 0)
    return 1.0;
  else {
    if (m == 1)
      return 1.0 + b - x;
    else {
      gln = 1.0+b-x;
      glnm1 = 1.0;
      for (i=2; i <= m; i++) {
        glnp1 = ((2*(i-1)+b+1.0-x)*gln - (i-1+b)*glnm1)/i;
        glnm1 = gln;
        gln = glnp1;
      }
      return gln;
    }
  }
} /* glaguerre */

#define EPSILON         0.00001 /* small number used in hermite root-finding */
#define maxcategs       9    /* maximum number of site types */
#define maxcategs2     11    /* maximum number of site types + 2 */
typedef double raterootarray[maxcategs2][maxcategs2];

void lgr(long m, double b, raterootarray lgroot)
{ /* For use by initgammacat.  Get roots of m-th Generalized Laguerre
     polynomial, given roots of (m-1)-th, these are to be
     stored in lgroot[m][] */
  long i;
  double upper, lower, x, y;
  bool dwn;   /* is function declining in this interval? */

  if (m == 1) {
    lgroot[1][1] = 1.0+b;
  } else {
    dwn = true;
    for (i=1; i<=m; i++) {
      if (i < m) {
        if (i == 1)
          lower = 0.0;
        else
          lower = lgroot[m-1][i-1];
        upper = lgroot[m-1][i];
      } else {                 /* i == m, must search above */
        lower = lgroot[m-1][i-1];
        x = lgroot[m-1][m-1];
        do {
          x = 2.0*x;
          y = glaguerre(m, b, x);
        } while ((dwn && (y > 0.0)) || ((!dwn) && (y < 0.0)));
        upper = x;
      }
      while (upper-lower > 0.000000001) {
        x = (upper+lower)/2.0;
        if (glaguerre(m, b, x) > 0.0) {
          if (dwn)
            lower = x;
          else
            upper = x;
        } else {
          if (dwn)
            upper = x;
          else
            lower = x;
        }        
      }
      lgroot[m][i] = (lower+upper)/2.0;
      dwn = !dwn;                /* switch for next one */
    }
  }
} /* lgr */

void initlaguerrecat(long categs, double alpha, double *rate, double *probcat)
{ /* calculate rates and probabilities to approximate Gamma distribution
     of rates with "categs" categories and shape parameter "alpha" using
     rates and weights from Generalized Laguerre quadrature */
  long i;
  raterootarray lgroot; /* roots of GLaguerre polynomials */
  double f, x, xi, y;

  alpha = alpha - 1.0;
  lgroot[1][1] = 1.0+alpha;
  for (i = 2; i <= categs; i++)
    lgr(i, alpha, lgroot);                   /* get roots for L^(a)_n */
  /* here get weights */
  /* Gamma weights are (1+a)(1+a/2) ... (1+a/n)*x_i/((n+1)^2 [L_{n+1}^a(x_i)]^2)  */
  f = 1;
  for (i = 1; i <= categs; i++)
    f *= (1.0+alpha/i);
  for (i = 1; i <= categs; i++) {
    xi = lgroot[categs][i];
    y = glaguerre(categs+1, alpha, xi);
    x = f*xi/((categs+1)*(categs+1)*y*y);
    rate[i-1] = xi/(1.0+alpha);
    probcat[i-1] = x;
  }
} /* initlaguerrecat */


double hermite(long n, double x)
{ /* calculates hermite polynomial with degree n and parameter x */
  /* seems to be unprecise for n>13 -> root finder does not converge*/
  double h1 = 1.;
  double h2 = 2. * x;
  double xx = 2. * x;
  long i;

  for (i = 1; i < n; i++) {
    xx = 2. * x * h2 - 2. * (i) * h1;
    h1 = h2;
    h2 = xx;
  }
  return xx;
} /* hermite */




double halfroot(double (*func)(long m, double x), long n, double startx,
                double delta)
{ /* searches from the bound (startx) only in one direction
     (by positive or negative delta, which results in
     other-bound=startx+delta)
     delta should be small.
     (*func) is a function with two arguments  */
  double xl;
  double xu;
  double xm = 0;
  double fu;
  double fl;
  double fm = 100000.;
  double gradient;
  bool dwn = false; 

  /* decide if we search above or below startx and escapes to trace back
     to the starting point that most often will be
     the root from the previous calculation */
  if (delta < 0) {
    xu = startx;
    xl = xu + delta;
  } else {
    xl = startx;
    xu = xl + delta;
  }
  delta = fabs(delta);
  fu = (*func)(n, xu);
  fl = (*func)(n, xl);
  gradient = (fl-fu)/(xl-xu);
  while(fabs(fm) > EPSILON) {        /* is root outside of our bracket?*/
    if ((fu<0.0 && fl<0.0) || (fu>0.0 && fl > 0.0)) {
      xu += delta;
      fu = (*func)(n, xu);
      fl = (*func)(n, xl);
      gradient = (fl-fu)/(xl-xu);
      dwn = (gradient < 0.0) ? true : false;
    } else {
      xm = xl - fl / gradient;
      fm = (*func)(n, xm);
      if (dwn) {
        if (fm > 0.) {
          xl = xm;
          fl = fm;
        } else {
          xu = xm;
          fu = fm;
        }
      } else {
        if (fm > 0.) {
          xu = xm;
          fu = fm;
        } else {
          xl = xm;
          fl = fm;
        }
      }
      gradient = (fl-fu)/(xl-xu);
    }
  }
  return xm;
} /* halfroot */

void root_hermite(long n, double *hroot)
{ /* find roots of Hermite polynmials */
  long z;
  long ii;
  long start;

  if (n % 2 == 0) {
    start = n/2;
    z = 1;
  } else {
    start = n/2 + 1;
    z=2;
    hroot[start-1] = 0.0;
  }
  for (ii = start; ii < n; ii++) {         /* search only upwards*/
    hroot[ii] = halfroot(hermite, n, hroot[ii-1]+EPSILON, 1./n);
    hroot[start - z] = -hroot[ii];
    z++;
  }
} /* root_hermite */

double logfac (long n)
{ /* log(n!) values were calculated with Mathematica
     with a precision of 30 digits */
  long i;
  double x;

  switch (n)
    {
    case 0:
      return 0.;
    case 1:
      return 0.;
    case 2:
      return 0.693147180559945309417232121458;
    case 3:
      return 1.791759469228055000812477358381;
    case 4:
      return 3.1780538303479456196469416013;
    case 5:
      return 4.78749174278204599424770093452;
    case 6:
      return 6.5792512120101009950601782929;
    case 7:
      return 8.52516136106541430016553103635;
    case 8:
      return 10.60460290274525022841722740072;
    case 9:
      return 12.80182748008146961120771787457;
    case 10:
      return 15.10441257307551529522570932925;
    case 11:
      return 17.50230784587388583928765290722;
    case 12:
      return 19.98721449566188614951736238706;
    default:
      x = 19.98721449566188614951736238706;
      for (i = 13; i <= n; i++)
        x += log(i);
      return x;
    }
} /* logfac */

void hermite_weight(long n, double * hroot, double * weights)
{
  /* calculate the weights for the hermite polynomial at the roots
     using formula from Abramowitz and Stegun chapter 25.4.46 p.890 */
  long i;
  double hr2;
  double numerator;

  numerator = exp(0.6931471805599 * ( n-1.) + logfac(n)) / (n*n);
  for (i = 0; i < n; i++) {
    hr2 = hermite(n-1, hroot[i]);
    weights[i] = numerator / (hr2*hr2);
  }
} /* hermiteweight */


#define SQRT2  1.4142135623730950488
void inithermitcat(long categs, double alpha, double *rate, double *probcat)
{ /* calculates rates and probabilities */
  long i;
  double *hroot;
  double std;

  std = SQRT2 /sqrt(alpha);
  hroot = (double *) malloc((categs+1) * sizeof(double));
  root_hermite(categs, hroot);         /* calculate roots */
  hermite_weight(categs, hroot, probcat);  /* set weights */
  for (i=0; i<categs; i++) {           /* set rates */
    rate[i] = 1.0 + std * hroot[i];
    probcat[i] = probcat[i];
  }
  free(hroot);
} /* inithermitcat */

CORAX_EXPORT int corax_compute_gamma_cats_laguerre(double       alpha,
                                                      unsigned int categories,
                                                      double      *output_rates,
                                                      double *output_weights)
{
  if ( categories > maxcategs) {
      corax_set_error(CORAX_ERROR_INVALID_PARAM, "gamma mode with laguerre quadrature supports only up to %u categories, requested %u", maxcategs, categories);
      return CORAX_FAILURE;
  }

  if (alpha >= 100.0)
    inithermitcat(categories, alpha, output_rates, output_weights);
  else
    initlaguerrecat(categories, alpha, output_rates, output_weights);

  return CORAX_SUCCESS;
}

CORAX_EXPORT int corax_compute_gamma_cats(double       alpha,
                                        unsigned int categories,
                                        double      *output_rates,
                                        double      *output_weights,
                                        int          mode)
{
    switch (mode) {
        case CORAX_GAMMA_RATES_MEAN:
        case CORAX_GAMMA_RATES_MEDIAN:
            return corax_compute_gamma_rates(alpha, categories, output_rates, mode);
        case CORAX_GAMMA_RATES_WEIGHTS_LLOYD_MAX:
            return corax_compute_gamma_cats_lloyd_max(alpha, categories, output_rates, output_weights);
        case CORAX_GAMMA_RATES_WEIGHTS_LAGUERRE:
            return corax_compute_gamma_cats_laguerre(alpha, categories, output_rates, output_weights);
        default:
            corax_set_error(CORAX_ERROR_INVALID_PARAM, "Invalid gamma mode %i", mode);
            return CORAX_FAILURE;
    }

    return CORAX_SUCCESS;
}
