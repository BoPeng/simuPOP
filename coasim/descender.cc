/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "descender.hh"

#ifndef CORE__NODE_HH_INCLUDED
# include "node.hh"
#endif
#ifndef CORE__MARKER_HH_INCLUDED
# include "marker.hh"
#endif

using namespace core;

// FIXME: the priority system right now has two possible priorities --
// this could of course be generalized if at some point it is needed.
// It would require a re-write of this, though.

void core::Descender::evolve(ARG &arg) const
{
  std::vector<RetiredInterval>::const_iterator rm_itr, rm_begin, rm_end;
  rm_begin = arg.retired_intervals().begin();
  rm_end   = arg.retired_intervals().end();

  // Both markers and retired intervals are sorted, so we can perform
  // a merge of the two to do the mutating

  // Markers to be mutated first
  int m = 0;
  for (rm_itr = rm_begin; rm_itr != rm_end; ++rm_itr)
  {
    for (; m < m_conf.no_markers(); ++m)
    {
      if (!m_conf.is_first_marker(m))           continue;
      if (m_conf.position(m) < rm_itr->start()) continue;
      if (m_conf.position(m) > rm_itr->end())   break;

      first_retry:                                // handle retries when wrong freqs
      try { rm_itr->mutate(m_conf,m); }
      catch (Mutator::retry_mutation&)
      {
        goto first_retry;
      }
    }
  }

  // Remaining markers
  m = 0;
  for (rm_itr = rm_begin; rm_itr != rm_end; ++rm_itr)
  {
    for ( ; m < m_conf.no_markers(); ++m)
    {
      if (!m_conf.is_plain_marker(m))           continue;
      if (m_conf.position(m) < rm_itr->start()) continue;
      if (m_conf.position(m) > rm_itr->end())   break;

      plain_retry:                                // handle retries when wrong freqs
      try { rm_itr->mutate(m_conf,m); }
      catch (Mutator::retry_mutation&)
      {
        goto plain_retry;
      }
    }
  }
}
