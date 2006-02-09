/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "marker.hh"

#ifndef IOSTREAM_INCLUDED
# include <iostream>
# define IOSTREAM_INCLUDED
#endif

namespace core
{
  Marker::Marker(const Marker &other)
    : m_position(other.m_position), m_values(other.m_values)
  {
  }

  Marker::~Marker()
  {
  }

  void Marker::to_text(std::ostream &os) const
  {
    os << this->type();
  }
}
