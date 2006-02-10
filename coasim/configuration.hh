/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__CONFIGURATION_HH_INCLUDED
#define CORE__CONFIGURATION_HH_INCLUDED

#ifndef CORE__MARKER_HH_INCLUDED
# include "marker.hh"
#endif

#ifndef STDEXCEPT_INCLUDED
# include <stdexcept>
# define STDEXCEPT_INCLUDED
#endif
#ifndef VECTOR_INCLUDED
# include <vector>
# define VECTOR_INCLUDED
#endif
#ifndef LIST_INCLUDED
# include <list>
# define LIST_INCLUDED
#endif
#ifndef CASSERT_INCLUDED
# include <cassert>
# define VECTOR_INCLUDED
#endif
#ifndef LIMITS_INCLUDED
# include <limits>
# define LIMITS_INCLUDED
#endif

namespace core
{
  class SimulationMonitor;

  class Scheduler;
  class State;

  struct Event
  {
    Event()
    {
    }

    virtual ~Event();

    // polymorphic copying
    virtual Event *copy() const = 0;

    virtual double earliest_event() const;

    virtual double event_time(State &s, double current_time)
      = 0;

    virtual void update_state(Scheduler &scheduler, State &s,
      double event_time)
      = 0;

    virtual void print_(std::ostream &os) const = 0;
  };

  inline std::ostream &
    operator<<(std::ostream &os, const Event &e)
  {
    e.print_(os);
    return os;
  }

  class Scheduler
  {
    public:
      ~Scheduler();

      void add_event(Event *event);

      void remove_event(Event *event);

      typedef std::pair<double,Event*> time_event_t;

      time_event_t next_event(State &s, double current_time);

    private:

      std::list<Event*> m_events;

  };

  // Exception thrown if the configuration is initialized with
  // population sizes <= 0
  struct non_pos_pop_size : public std::logic_error
  {
    non_pos_pop_size() :
    std::logic_error("Population sizes must be greater than 0."){}
  };

  // Exception thrown if the configuration is initialized with
  // un-sorted positions
  struct out_of_sequence : public std::logic_error
  {
    out_of_sequence() : std::logic_error("Marker positions not sorted."){}
  };

  // Exception thrown if we try to add to a value set before the type
  // of the marker has been initialized
  struct uninitialized_marker : public std::logic_error
  {
    uninitialized_marker() : std::logic_error("uninitialized marker."){}
  };

  class Configuration
  {
    public:

      // Need more work to make the python binding work,
      // but this has already broken other (scheme) bindings...
      // Maybe I should define two constructors?
      typedef std::vector< unsigned int > PopSizeVec;
      typedef std::vector< unsigned int >::const_iterator pop_size_itr_t;
      typedef std::vector< Marker* > MarkerVec;
      typedef std::vector< Event* > EpochVec;

      // template <typename PopSizeVec, typename MarkerVec,  typename EpochVec>
      Configuration(
        const PopSizeVec& popSizes,
        const MarkerVec& markers,
        const EpochVec&  epochs,
        double rho,
        double Q,
        double gamma,
        double growth)
        : m_rho(rho), m_Q(Q), m_gamma(gamma), m_growth(growth)
      {
        // quick variable conversion
        PopSizeVec::const_iterator ps_begin = popSizes.begin();
        PopSizeVec::const_iterator ps_end = popSizes.end();
        MarkerVec::const_iterator m_begin = markers.begin();
        MarkerVec::const_iterator m_end = markers.end();
        EpochVec::const_iterator e_begin = epochs.begin();
        EpochVec::const_iterator e_end = epochs.end();

        m_no_markers = m_end - m_begin;
        m_no_leaves = 0;

        for ( ; ps_begin!=ps_end ; ++ps_begin)
        {
          if (*ps_begin <= 0) throw non_pos_pop_size();

          m_pop_sizes.push_back(*ps_begin);
          m_no_leaves += *ps_begin;
        }

        m_first_markers = new const Marker* [m_no_markers];
        for (int m = 0; m < m_no_markers; ++m) m_first_markers[m] = 0;

        m_plain_markers = new const Marker* [m_no_markers];
        for (int m = 0; m < m_no_markers; ++m) m_plain_markers[m] = 0;

        for (int m = 0; m < m_no_markers; ++m)
          set_marker(m, *(m_begin++));

        for (int m = 1; m < m_no_markers; ++m)
          if (position(m-1) >= position(m)) throw out_of_sequence();

        for ( ; e_begin != e_end; ++e_begin)
          m_events.push_back((*e_begin)->copy());
      }

      // the original constructor
      // Ingored in python binding
      template <typename PopSizeItr, typename MarkerItr,  typename EpochItr>
        Configuration::Configuration(PopSizeItr ps_begin, PopSizeItr ps_end,
        MarkerItr  m_begin,  MarkerItr  m_end,
        EpochItr   e_begin,  EpochItr   e_end,
        double rho, double Q, double gamma,
        double growth)
        : m_rho(rho), m_Q(Q), m_gamma(gamma), m_growth(growth)
      {
        m_no_markers = m_end - m_begin;
        m_no_leaves = 0;

        for ( ; ps_begin!=ps_end ; ++ps_begin)
        {
          if (*ps_begin <= 0) throw non_pos_pop_size();

          m_pop_sizes.push_back(*ps_begin);
          m_no_leaves += *ps_begin;
        }

        m_first_markers = new const Marker* [m_no_markers];
        for (int m = 0; m < m_no_markers; ++m) m_first_markers[m] = 0;

        m_plain_markers = new const Marker* [m_no_markers];
        for (int m = 0; m < m_no_markers; ++m) m_plain_markers[m] = 0;

        for (int m = 0; m < m_no_markers; ++m)
          set_marker(m, *(m_begin++));

        for (int m = 1; m < m_no_markers; ++m)
          if (position(m-1) >= position(m)) throw out_of_sequence();

        for ( ; e_begin != e_end; ++e_begin)
          m_events.push_back((*e_begin)->copy());
      }

      ~Configuration();

      const std::vector<unsigned int> pop_sizes() const
      {
        return m_pop_sizes;
      }

      pop_size_itr_t pop_sizes_begin() const
      {
        return m_pop_sizes.begin();
      }

      pop_size_itr_t pop_sizes_end()   const
      {
        return m_pop_sizes.end();
      }

      const unsigned int no_leaves() const
      {
        return m_no_leaves;
      }

      // Number of markers for the configuration
      int no_markers() const
      {
        return m_no_markers;
      }

      // The positions of the markers
      double position(int index) const
      {
        return marker(index).position();
      }

      // Accessors to markers
      const Marker& first_marker(int index) const
      {
        if (index >= no_markers())
          throw std::out_of_range("No marker at index");
        if (!m_first_markers[index])
          throw uninitialized_marker();
        else
          return *m_first_markers[index];
      }

      const Marker& plain_marker(int index) const
      {
        if (index >= no_markers())
          throw std::out_of_range("No marker at index");
        if (!m_plain_markers[index])
          throw uninitialized_marker();
        else
          return *m_plain_markers[index];
      }

      const Marker& marker(int index) const
      {
        if (is_first_marker(index))
          return first_marker(index);
        else
          return plain_marker(index);
      }

      bool is_first_marker(int index) const
      {
        return m_first_markers[index] != 0;
      }

      bool is_plain_marker(int index) const
      {
        return m_plain_markers[index] != 0;
      }

      // Insert a marker at the position index -- this method only borrows
      // the reference, so don't free or change the marker after setting
      // it and before deleting the configuration -- and remember to free
      // it yourself after use of the configuration.  The run_first flag
      // is used to prioritize the order of marker-mutations; all markers
      // set with run_first are mutated before all markers without
      // run_first.
      void set_marker(int index, const Marker *marker)
      {
        assert(marker != 0);

        if (index >= no_markers())
          throw std::out_of_range("No marker at index");

        if (marker->run_first())
        {
          m_first_markers[index] = marker->copy();
          m_plain_markers[index] = 0;
        }
        else
        {
          m_first_markers[index] = 0;
          m_plain_markers[index] = marker->copy();
        }
      }

      // Parameters for building the ARG and assigning mutations
      double rho()    const
      {
        return m_rho;
      }

      double Q()      const
      {
        return m_Q;
      }

      double gamma()  const
      {
        return m_gamma;
      }

      double growth() const
      {
        return m_growth;
      }

      typedef std::vector<Event*>::const_iterator epoch_iterator;

      epoch_iterator epochs_begin() const
      {
        return m_events.begin();
      }

      epoch_iterator epochs_end()   const
      {
        return m_events.end();
      }

      void sort_events() const;                   // sort events wrt. their start time
      // (if it is fixed).  Prevents
      // incorrect nesting of non-nested
      // epochs...

    private:
      // Disable these
      Configuration(const Configuration&);

      Configuration& operator = (const Configuration&);

      int m_no_markers;
      unsigned int m_no_leaves;

      std::vector<unsigned int> m_pop_sizes;

      const Marker** m_first_markers;
      const Marker** m_plain_markers;

      double m_rho;
      double m_Q;
      double m_gamma;
      double m_growth;

      mutable std::vector<Event*> m_events;
  };

}
#endif
