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

core::Marker::Marker(const Marker &other)
    : i_position(other.i_position), i_values(other.i_values)
{
}

core::Marker::~Marker() 
{
}

void 
core::Marker::to_text(std::ostream &os) const
{
    os << this->type();
}

