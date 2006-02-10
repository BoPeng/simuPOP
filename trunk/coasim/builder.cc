/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "builder.hh"

#ifndef CORE__DIST_FUNCTIONS_HH_INCLUDED
# include "dist_funcs.hh"
#endif
#ifndef CORE__NODE_HH_INCLUDED
# include "node.hh"
#endif
#ifndef CORE__BUILDER_EVENTS_HH_INCLUDED
# include "builder_events.hh"
#endif
#ifndef CORE__EPOCHS_HH_INCLUDED
# include "epochs.hh"
#endif

using namespace core;

ARG * Builder::build(BuilderMonitor *callbacks,
bool keep_empty_intervals) const
{
  using namespace Distribution_functions;

  std::auto_ptr<ARG> arg(new ARG(m_conf, keep_empty_intervals));

  double time = 0.0;
  State state(*arg, callbacks,
    m_conf.pop_sizes_begin(), m_conf.pop_sizes_end());

  Scheduler scheduler;
  unsigned int pop_no;
  std::vector<Population>::iterator j;
  for (j = state.populations().begin(), pop_no = 0;
    j != state.populations().end(); ++j, ++pop_no)
  {
    scheduler.add_event(j->coalescence_event());
    if (m_conf.growth() > 0)
      scheduler.add_event(new Growth(pop_no,m_conf.growth(), 0));
  }

  if (m_conf.rho() > 0)
    scheduler.add_event(new RecombinationEvent(m_conf.rho()));
  if (m_conf.gamma() > 0)
    scheduler.add_event(new GeneConversionEvent(m_conf.gamma(),
      m_conf.Q()));

  std::vector<Event*>::const_iterator i;
  m_conf.sort_events();
  for (i = m_conf.epochs_begin(); i != m_conf.epochs_end(); ++i)
  {
    Event *e = (*i)->copy();
    scheduler.add_event(e);
  }

  while (state.total_population_size() > 1)
  {
    Scheduler::time_event_t e = scheduler.next_event(state, time);
    assert(e.second);
    e.second->update_state(scheduler, state, e.first);
    time = e.first;
  }

  arg->sort_retired_intervals();                  // NB! important, since the
  // remaining functions rely on the
  // retired intervals being sorted!
  return arg.release();
}
