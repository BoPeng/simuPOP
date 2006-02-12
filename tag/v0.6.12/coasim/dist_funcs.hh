/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__DIST_FUNCTIONS_HH_INCLUDED
#define CORE__DIST_FUNCTIONS_HH_INCLUDED

namespace Distribution_functions
{
  double expdev(const double param);

  double expdev(const int fac, const double param);

  double expdist(const double param, const double x);

  double uniform();

  int uniform(double part_1, double part_2);

  int uniform(double part_1, double part_2, double part_3);

  int random_sign();

  int irand(int n);

  void two_int_rand(int& n1, int& n2, int n);
}
#endif
