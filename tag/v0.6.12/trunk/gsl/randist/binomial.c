/* randist/binomial.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

/* The binomial distribution has the form,

   prob(k) =  n!/(k!(n-k)!) *  p^k (1-p)^(n-k) for k = 0, 1, ..., n

   This is the algorithm from Knuth */

unsigned int
gsl_ran_binomial (const gsl_rng * r, double p, unsigned int n)
{
  unsigned int i, a, b, k = 0;

  while (n > 10)        /* This parameter is tunable */
    {
      double X;
      a = 1 + (n / 2);
      b = 1 + n - a;

      X = gsl_ran_beta (r, (double) a, (double) b);

      if (X >= p)
        {
          n = a - 1;
          p /= X;
        }
      else
        {
          k += a;
          n = b - 1;
          p = (p - X) / (1 - X);
        }
    }

  for (i = 0; i < n; i++)
    {
      double u = gsl_rng_uniform (r);
      if (u < p)
        k++;
    }

  return k;
}

double
gsl_ran_binomial_pdf (const unsigned int k, const double p,
                      const unsigned int n)
{
  if (k > n)
    {
      return 0;
    }
  else
    {
      double P;

      double ln_Cnk = gsl_sf_lnchoose (n, k);
      P = ln_Cnk + k * log (p) + (n - k) * log (1 - p);
      P = exp (P);

      return P;
    }
}
