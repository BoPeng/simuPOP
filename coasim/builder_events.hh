/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__BUILDER_EVENTS_HH_INCLUDED
#define CORE__BUILDER_EVENTS_HH_INCLUDED

#ifndef CORE__CONFIGURATION_HH_INCLUDED
# include "configuration.hh"
#endif

#ifndef VECTOR_INCLUDED
# include <vector>
# define VECTOR_INCLUDED
#endif
#ifndef CASSERT_INCLUDED
# include <cassert>
# define CASSERT_INCLUDED
#endif
#ifndef STDEXCEPT_INCLUDED
# include <stdexcept>
# define STDEXCEPT_INCLUDED
#endif

namespace core
{

  class ARG;
  class Node;
  class BuilderMonitor;
  class Scheduler;
  class CoalescenceEvent;
  class State;

  class Population
  {
    public:
      Population(ARG &arg,
        int initial_population_size,
        CoalescenceEvent *coal_event,
        double scale_fraction = 1);

      void push(Node *n)
      {
        m_nodes.push_back(n);
      }

      Node *pop_random();

      Node *pop_last();

      size_t size() const
      {
        return m_nodes.size();
      }

      double scale_fraction() const
      {
        return m_scale_fraction;
      }

      CoalescenceEvent *coalescence_event()
      {
        return m_coal_event;
      }

      void coalescence_event(CoalescenceEvent *e)
      {
        m_coal_event = e;
      }

    private:

      CoalescenceEvent *m_coal_event;

      double m_scale_fraction;

      std::vector<Node*>m_nodes;
  };

  class State
  {
    public:
      template <typename Itr>
        State(ARG &arg, BuilderMonitor *callbacks,
        Itr sizes_begin, Itr sizes_end);

      std::vector<Population> &populations()
      {
        return m_populations;
      }

      ARG                     &arg()
      {
        return m_arg;
      }

      BuilderMonitor          *callbacks()
      {
        return m_callbacks;
      }

      size_t total_population_size() const;

      // returns a random population, where each population is
      // weighted with its size.
      Population &random_population();

    private:

      ARG &m_arg;

      std::vector<Population> m_populations;

      BuilderMonitor *m_callbacks;

  };

  // thrown if a population index is given outside the allowed
  // values
  struct illegal_population : std::logic_error
  {
    illegal_population() : std::logic_error("Illegal populaiton id") {}
  };

  class CoalescenceEvent : public Event
  {
    protected:
      int m_population;

    public:
      CoalescenceEvent(int population) : m_population(population)
      {
        if (population < 0) throw illegal_population();
      }

      virtual Event *copy() const;

      int population() const
      {
        return m_population;
      }

      virtual double waiting_time(State &s, double current_time);

      virtual double event_time  (State &s, double current_time);

      virtual void   update_state(Scheduler &scheduler, State &s,
        double event_time);

      virtual void print_(std::ostream &os) const;
  };

  class RecombinationEvent : public Event
  {
    public:
      RecombinationEvent(double rho) : m_rho(rho)
        {}

      virtual Event *copy() const;

      virtual double event_time  (State &s, double current_time);

      virtual void   update_state(Scheduler &scheduler, State &s,
        double event_time);

      virtual void print_(std::ostream &os) const;

    private:
      double m_rho;
  };

  class GeneConversionEvent : public Event
  {
    double m_gamma;
    double m_Q;
    public:
      GeneConversionEvent(double gamma, double Q) : m_gamma(gamma), m_Q(Q)
        {}
      virtual Event *copy() const;

      virtual double event_time  (State &s, double current_time);

      virtual void   update_state(Scheduler &scheduler, State &s,
        double event_time);

      virtual void print_(std::ostream &os) const;
  };

  // a merge moves a set of populations to the first and disables
  // all events in the remaining.
  class PopulationMerge : public Event
  {
    public:
      template <class Itr>
        PopulationMerge(Itr pop_begin, Itr pop_end, double merge_time)
        : m_populations(pop_begin, pop_end), m_merge_time(merge_time)
      {
        assert(m_populations.size() > 1);
        std::vector<int>::const_iterator i;
        for (i = m_populations.begin(); i != m_populations.end(); ++i, ++pop_begin)
          assert(*i >= 0);
        assert(merge_time >= 0);
      }

      virtual Event *copy() const;

      const std::vector<int> &populations() const
      {
        return m_populations;
      }

      double merge_time()  const
      {
        return m_merge_time;
      }

      virtual double event_time  (State &s, double current_time);

      virtual void   update_state(Scheduler &scheduler, State &s,
        double event_time);

      virtual double earliest_event() const;

      virtual void print_(std::ostream &os) const;

    private:

      std::vector<int> m_populations;

      double m_merge_time;

  };

  // FIXME: This initialization is not optimal
  template <typename Itr>
    State::State(ARG &arg, BuilderMonitor *callbacks,
    Itr sizes_begin, Itr sizes_end)
    : m_arg(arg), m_callbacks(callbacks)
  {
    for (int p_no = 0; sizes_begin != sizes_end; ++p_no, ++sizes_begin)
    {
      Population p(arg, *sizes_begin, new CoalescenceEvent(p_no));
      m_populations.push_back(p);
    }
  }

}
#endif
