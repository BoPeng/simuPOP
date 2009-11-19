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

%module gsl

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
