/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "retired_interval.hh"

#ifndef CORE__NODE_HH_INCLUDED
# include "node.hh"
#endif
#ifndef CORE__CONFIGURATION_HH_INCLUDED
# include "configuration.hh"
#endif
#ifndef CORE__MARKER_HH_INCLUDED
# include "marker.hh"
#endif

using namespace core;

void
core::RetiredInterval::calc_surface() const
{
  i_surface = i_top_node->surface_at_point(start());
}


void
core::RetiredInterval::mutate(const Configuration &conf,
unsigned int marker_index)  const
{
  const Marker &marker = conf.marker(marker_index);
  i_top_node->initialize_marker(marker_index, marker);
  std::auto_ptr<Mutator> mutator(marker.create_mutator(conf,*this));
  i_top_node->mutate_marker(marker_index,*mutator);
}


void
core::RetiredInterval::to_xml(std::ostream &os) const
{
  os << "  <interval_node id=\"i_" << this << "\">" << std::endl
    << "    <child ref=\"i_" << top_node() << "\"/>" << std::endl
    << "    <interval start=\"" << start() << "\" end=\"" << end() << "\"/>\n"
    << "  </interval_node>" << std::endl;
}
