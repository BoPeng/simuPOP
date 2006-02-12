/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__MONITOR_HH_INCLUDED
#define CORE__MONITOR_HH_INCLUDED

#ifndef STDEXCEPT_INCLUDED
# include <stdexcept>
# define STDEXCEPT_INCLUDED
#endif

namespace core
{
  class SimulationMonitor
  {

    public:

      // thrown this if the monitor wants the simulation aborted
      class AbortSimulation : public std::exception {};

      virtual void start_arg_building(unsigned int no_leaves) = 0;
      virtual void builder_update(unsigned int no_nodes,
        unsigned int no_top_nodes,
        unsigned long int no_iterations,
        double cur_time,
        unsigned int no_coal_events,
        unsigned int no_gene_conv_events,
        unsigned int no_recomb_events) = 0;

      virtual void builder_termination(unsigned int no_nodes,
        unsigned int no_top_nodes,
        unsigned long int no_iterations,
        double cur_time,
        unsigned int no_coal_events,
        unsigned int no_gene_conv_events,
        unsigned int no_recomb_events) = 0;

      virtual void start_mutating() = 0;
      virtual void mutator_update(unsigned int marker_no) = 0;
      virtual void retry_mutation() = 0;
      virtual void retry_arg_building() = 0;

      virtual void simulation_terminated() = 0;
  };

}
#endif
