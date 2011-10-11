//
// $File: gsl.i $
// $LastChangedDate: 2009-11-09 17:36:19 -0800 (Mon, 09 Nov 2009) $
// $Rev: 3106 $
//
// This file is part of simuPOP, a forward-time population genetics
// simulation environment. Please visit http://simupop.sourceforge.net
// for details.
//
// Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//

%define DOCSTRING
"This module exposes the following GSL (GUN Scientific Library) functions used
by simuPOP to the user interface. Although more functions may be added from time
to time, this module is not intended to become a complete wrapper for GSL. Please
refer to the GSL reference manual (http://www.gnu.org/software/gsl/manual/html_node/)
for details about these functions. Note that random number generation functions
are wrapped into the simuPOP.RNG class.

- ``gsl_cdf_gaussian_P(x, sigma)``
- ``gsl_cdf_gaussian_Q(x, sigma)``
- ``gsl_cdf_gaussian_Pinv(P, sigma)``
- ``gsl_cdf_gaussian_Qinv(Q, sigma)``
- ``gsl_cdf_ugaussian_P(x)``
- ``gsl_cdf_ugaussian_Q(x)``
- ``gsl_cdf_ugaussian_Pinv(P)``
- ``gsl_cdf_ugaussian_Qinv(Q)``
- ``gsl_cdf_exponential_P(x, mu)``
- ``gsl_cdf_exponential_Q(x, mu)``
- ``gsl_cdf_exponential_Pinv(P, mu)``
- ``gsl_cdf_exponential_Qinv(Q, mu)``
- ``gsl_cdf_chisq_P(x, nu)``
- ``gsl_cdf_chisq_Q(x, nu)``
- ``gsl_cdf_chisq_Pinv(P, nu)``
- ``gsl_cdf_chisq_Qinv(Q, nu)``
- ``gsl_cdf_gamma_P(x, a, b)``
- ``gsl_cdf_gamma_Q(x, a, b)``
- ``gsl_cdf_gamma_Pinv(P, a, b)``
- ``gsl_cdf_gamma_Qinv(Q, a, b)``
- ``gsl_ran_gamma_pdf(x, a, b)``
- ``gsl_cdf_beta_P(x, a, b)``
- ``gsl_cdf_beta_Q(x, a, b)``
- ``gsl_cdf_beta_Pinv(P, a, b)``
- ``gsl_cdf_beta_Qinv(Q, a, b)``
- ``gsl_ran_beta_pdf(x, a, b)``
- ``gsl_cdf_binomial_P(k, p, n)``
- ``gsl_cdf_binomial_Q(k, p, n)``
- ``gsl_ran_binomial_pdf(k, p, n)``
- ``gsl_cdf_poisson_P(k, mu)``
- ``gsl_cdf_poisson_Q(k, mu)``
- ``gsl_ran_poisson_pdf(k, mu)``
"
%enddef

%module(docstring=DOCSTRING) gsl

%{

#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_cdf.h"

%}

%inline %{

void my_error_handler(const char * reason, const char * file,
                       int line, int gsl_errno)
{
    fprintf(stderr, "GSL Error %d:\t%s", gsl_errno, reason);
}


int initialize()
{
    gsl_set_error_handler(&my_error_handler);
    return 1;
}

// this function is copied here to avoid the inclusion of
// rng related functions from the inclusion of gsl/randist/gamma.c.
double gsl_ran_gamma_pdf (const double x, const double a, const double b)
{
  if (x < 0)
    {
      return 0 ;
    }
  else if (x == 0)
    {
      if (a == 1)
        return 1/b ;
      else
        return 0 ;
    }
  else if (a == 1)
    {
      return exp(-x/b)/b ;
    }
  else
    {
      double p;
      double lngamma = gsl_sf_lngamma (a);
      p = exp ((a - 1) * log (x/b) - x/b - lngamma)/b;
      return p;
    }
}

/* this function is copied from randist/binomial.c to avoid including
 * more files. log1p is replaced by gsl_log1p to avoid a compiling program
 *  caused by the order of inclusion file.
 */
extern double gsl_log1p (const double x);

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

      if (p == 0) 
        {
          P = (k == 0) ? 1 : 0;
        }
      else if (p == 1)
        {
          P = (k == n) ? 1 : 0;
        }
      else
        {
          double ln_Cnk = gsl_sf_lnchoose (n, k);
          P = ln_Cnk + k * log (p) + (n - k) * gsl_log1p (-p);
          P = exp (P);
        }

      return P;
    }
}

double
gsl_ran_beta_pdf (const double x, const double a, const double b)
{
  if (x < 0 || x > 1)
    {
      return 0 ;
    }
  else 
    {
      double p;

      double gab = gsl_sf_lngamma (a + b);
      double ga = gsl_sf_lngamma (a);
      double gb = gsl_sf_lngamma (b);
      
      if (x == 0.0 || x == 1.0) 
        {
          p = exp (gab - ga - gb) * pow (x, a - 1) * pow (1 - x, b - 1);
        }
      else
        {
          p = exp (gab - ga - gb + log(x) * (a - 1)  + gsl_log1p(-x) * (b - 1));
        }

      return p;
    }
}


double
gsl_ran_poisson_pdf (const unsigned int k, const double mu)
{
  double p;
  double lf = gsl_sf_lnfact (k); 

  p = exp (log (mu) * k - lf - mu);
  return p;
}
%}

%init
%{
    initialize();
%}

extern double gsl_cdf_gaussian_P(double x, double sigma); 
extern double gsl_cdf_gaussian_Q(double x, double sigma); 
extern double gsl_cdf_gaussian_Pinv(double P, double sigma); 
extern double gsl_cdf_gaussian_Qinv(double Q, double sigma); 
extern double gsl_cdf_ugaussian_P(double x); 
extern double gsl_cdf_ugaussian_Q(double x); 
extern double gsl_cdf_ugaussian_Pinv(double P); 
extern double gsl_cdf_ugaussian_Qinv(double Q); 

extern double gsl_cdf_exponential_P(double x, double mu); 
extern double gsl_cdf_exponential_Q(double x, double mu); 
extern double gsl_cdf_exponential_Pinv(double P, double mu); 
extern double gsl_cdf_exponential_Qinv(double Q, double mu); 

extern double gsl_cdf_chisq_P(double x, double nu); 
extern double gsl_cdf_chisq_Q(double x, double nu); 
extern double gsl_cdf_chisq_Pinv(double P, double nu); 
extern double gsl_cdf_chisq_Qinv(double Q, double nu); 

extern double gsl_cdf_gamma_P(double x, double a, double b); 
extern double gsl_cdf_gamma_Q(double x, double a, double b); 
extern double gsl_cdf_gamma_Pinv(double P, double a, double b); 
extern double gsl_cdf_gamma_Qinv(double Q, double a, double b); 

extern double gsl_cdf_binomial_P(unsigned int k, double p, unsigned int n); 
extern double gsl_cdf_binomial_Q(unsigned int k, double p, unsigned int n); 
/* This function is copied into this interface file to avoid inclusion of other
   cdf functions. */
/* extern double gsl_ran_binomial_pdf(unsigned int k, double p, unsigned int n); */

extern double gsl_cdf_beta_P(double x, double a, double b); 
extern double gsl_cdf_beta_Q(double x, double a, double b); 
extern double gsl_cdf_beta_Pinv(double P, double a, double b); 
extern double gsl_cdf_beta_Qinv(double Q, double a, double b); 

extern double gsl_cdf_poisson_P(const unsigned int k, const double mu);
extern double gsl_cdf_poisson_Q(const unsigned int k, const double mu);
extern double gsl_ran_poisson_pdf(const unsigned int k, const double mu);
