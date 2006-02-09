/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "simulator.hh"

#ifndef CORE__BUILDER_HH_INCLUDED
# include "builder.hh"
#endif
#ifndef CORE__DESCENDER_HH_INCLUDED
# include "descender.hh"
#endif
#ifndef CORE__MARKER_HH_INCLUDED
# include "marker.hh"
#endif
#ifndef CORE__MONITOR_HH_INCLUDED
# include "monitor.hh"
#endif
#ifndef CORE__NODE_HH_INCLUDED
# include "node.hh"
#endif

#ifndef SYS_TIME_H_INCLUDED
# include <sys/time.h>
# define SYS_TIME_H_INCLUDED
#endif

using namespace core;

ARG *
core::Simulator::simulate(const Configuration &conf,
BuilderMonitor *build_callbacks,
bool keep_empty_intervals,
unsigned int random_seed)
{
  Builder builder(conf);
  Descender descender(conf);
  ARG *arg = 0;

  // set rand seed
  if (!random_seed)
  {
    // use time if no seed is explicitly given
    struct timeval tv; struct timezone tz;
    gettimeofday(&tv,&tz);
    random_seed = tv.tv_usec;
  }
  std::srand(random_seed);

  try
  {

    retry:
    try
    {
      arg = builder.build(build_callbacks, keep_empty_intervals);
      descender.evolve(*arg);
    }
    catch (Mutator::retry_arg&)
    {
      delete arg; arg = 0;
      goto retry;
    }

  }
  catch(SimulationMonitor::AbortSimulation&)
  {
    if (arg) delete arg; arg = 0;
  }

  return arg;
}
