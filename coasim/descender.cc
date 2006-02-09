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
  std::vector<RetiredInterval>::const_iterator ri_itr, ri_begin, ri_end;
  ri_begin = arg.retired_intervals().begin();
  ri_end   = arg.retired_intervals().end();

  // Both markers and retired intervals are sorted, so we can perform
  // a merge of the two to do the mutating

  // Markers to be mutated first
  int m = 0;
  for (ri_itr = ri_begin; ri_itr != ri_end; ++ri_itr)
  {
    for (; m < i_conf.no_markers(); ++m)
    {
      if (!i_conf.is_first_marker(m))           continue;
      if (i_conf.position(m) < ri_itr->start()) continue;
      if (i_conf.position(m) > ri_itr->end())   break;

      first_retry:                                // handle retries when wrong freqs
      try { ri_itr->mutate(i_conf,m); }
      catch (Mutator::retry_mutation&)
      {
        goto first_retry;
      }
    }
  }

  // Remaining markers
  m = 0;
  for (ri_itr = ri_begin; ri_itr != ri_end; ++ri_itr)
  {
    for ( ; m < i_conf.no_markers(); ++m)
    {
      if (!i_conf.is_plain_marker(m))           continue;
      if (i_conf.position(m) < ri_itr->start()) continue;
      if (i_conf.position(m) > ri_itr->end())   break;

      plain_retry:                                // handle retries when wrong freqs
      try { ri_itr->mutate(i_conf,m); }
      catch (Mutator::retry_mutation&)
      {
        goto plain_retry;
      }
    }
  }
}
