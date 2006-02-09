/* -*- Mode: C++; c-basic-offset: 4; -*- 
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "dist_funcs.hh"

#ifndef CMATH_INCLUDED
# include <cmath>
# define CMATH_INCLUDED
#endif
#ifndef CSTDLIB_INCLUDED
# include <cstdlib>
# define CSTDLIB_INCLUDED
#endif
#ifndef ALGORITHM_INCLUDED
# include <algorithm>
# define ALGORITHM_INCLUDED
#endif

#include <iostream>

namespace Distribution_functions
{
    double uniform()
    {
	return double(std::rand())/double(RAND_MAX);
    }

    double expdev(const double param)
    {
	return -log(double(std::rand())/double(RAND_MAX))/param;
    }

    double expdev(const int fac, const double param)
    {
	double u = uniform();
	return -log(u)/param/fac;
    }

  
    double expdist(const double param, const double x)
    {
	return 1.0-exp(-param*x);
    }

    int random_sign()
    {
	int r = -1;
	if (uniform()>0.5) r = 1;
	return r;
    }
  
    int uniform(double part_1, double part_2)
    {
	double r = uniform()*(part_1+part_2);
	if (r<part_1) return 0;
	else          return 1;
    }

    int uniform(double part_1, double part_2, double part_3)
    {
	double r = uniform()*(part_1+part_2+part_3);
	if (r<part_1) return 0;
	else if (r<part_1+part_2) return 1;
	return 2;
    }
  
    int irand(int n)
    {
	int r = int((double(std::rand()) / RAND_MAX) * n);
	return std::min(n-1, r);
    }
  
    void two_int_rand(int& n1, int& n2, int n)
    {
	n1 = irand(n);
	n2 = irand(n-1);
	if (n1==n2) n2 = n-1;
    }
}
