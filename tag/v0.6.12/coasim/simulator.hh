/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__SIMULATOR_HH_INCLUDED
#define CORE__SIMULATOR_HH_INCLUDED

namespace core
{

  class ARG;
  class Configuration;
  class BuilderMonitor;

  namespace Simulator
  {
    // Returns the resulting ARG, or 0 if the simulation was aborted
    ARG *simulate(const Configuration &conf,
      BuilderMonitor *build_callbacks = 0,
      bool keep_empty_intervals = false,
      unsigned int random_seed = 0);
  }

}
#endif
