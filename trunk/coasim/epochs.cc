/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "epochs.hh"

#ifndef CORE__BUILDER_HH_INCLUDED
#  include "builder.hh"
#endif
#ifndef CORE__BUILDER_EVENTS_HH_INCLUDED
# include "builder_events.hh"
#endif
#ifndef CORE__DIST_FUNCTIONS_HH_INCLUDED
# include "dist_funcs.hh"
#endif

#ifndef IOSTREAM_INCLUDED
# include <iostream>
# define IOSTREAM_INCLUDED
#endif
#ifndef LIMITS_INCLUDED
# include <limits>
# define LIMITS_INCLUDED
#endif
#ifndef CMATH_INCLUDED
# include <cmath>
# define CMATH_INCLUDED
#endif

namespace core
{

  Epoch::~Epoch() {}

  double Epoch::event_time(State &s, double current_time)
  {
    if (current_time < start_time()) start_time();
    return std::min(nested_event_time(s, current_time), end_time());
  }

  void Epoch::update_state(Scheduler &scheduler, State &s, double event_time)
  {
    if (event_time >= end_time())
    {
      // remove epoch
      scheduler.remove_event(this);
      delete this;
    }
    else
      nested_update_state(scheduler, s, event_time);
  }

  double Epoch::earliest_event() const
  {
    return start_time();
  }

  CoalescenceEpoch::~CoalescenceEpoch()
  {
    delete m_underlying;
  }

  double CoalescenceEpoch::nested_event_time(State &s, double current_time)
  {
    // use start time until we have executed the first event that
    // pushes this to the coalescence event stack.
    if (!m_underlying) return start_time();
    return CoalescenceEvent::event_time(s, current_time);
  }

  void
    CoalescenceEpoch::nested_update_state(Scheduler &scheduler, State &s,
    double event_time)
  {
    CoalescenceEvent::update_state(scheduler, s, event_time);
  }

  double
    CoalescenceEpoch::event_time(State &s, double current_time)
  {
    if (current_time < start_time()) return start_time();
    return std::min(nested_event_time(s, current_time), end_time());
  }

  void
    CoalescenceEpoch::update_state(Scheduler &scheduler, State &s,
    double event_time)
  {
    if (start_time() <= event_time and !m_underlying)
      // we are entering the epoch -- push this epoch on top of the stack
    {
      Population &p = s.populations().at(population());
      m_underlying = p.coalescence_event();
      p.coalescence_event(this);
      scheduler.remove_event(m_underlying);

      enter_callback(s, event_time);
    }
    else if (end_time() > 0 and end_time() <= event_time)
      // we are leaving, now pop the epoch
    {
      leave_callback(s, event_time);

      Population &p = s.populations().at(population());
      p.coalescence_event(m_underlying);
      scheduler.remove_event(this);
      scheduler.add_event(m_underlying);
      m_underlying = 0;
      delete this;
    }
    else
      // inside the epoch we just propagate
      nested_update_state(scheduler, s, event_time);
  }

  void
    CoalescenceEpoch::enter_callback(State &s, double current_time)
  {
  }

  void
    CoalescenceEpoch::leave_callback(State &s, double current_time)
  {
  }

  double
    CoalescenceEpoch::earliest_event() const
  {
    return start_time();
  }

  double
    BottleNeckEpoch::waiting_time(State &s, double c_time)
  {
    return m_scale_fraction * basic_waiting_time(s, c_time);
  }

  Event *
    BottleNeckEpoch::copy() const
  {
    return new BottleNeckEpoch(*this);
  }

  void
    BottleNeckEpoch::enter_callback(State &s, double current_time)
  {
    BuilderMonitor *callbacks = s.callbacks();
    if (!callbacks) return;
    callbacks->bottleneck_callback(population(), true,
      current_time,
      s.total_population_size());
  }

  void
    BottleNeckEpoch::leave_callback(State &s, double current_time)
  {
    BuilderMonitor *callbacks = s.callbacks();
    if (!callbacks) return;
    callbacks->bottleneck_callback(population(), false,
      current_time,
      s.total_population_size());
  }

  void
    BottleNeckEpoch::print_(std::ostream &os) const
  {
    os << "BottleNeckEpoch(" << population() << ", "
      << start_time() << ", " << end_time() << ", "
      << m_scale_fraction << ')';
  }

  double
    GrowthEpoch::waiting_time(State &s, double current_time)
  {
    using namespace Distribution_functions;
    double t_k_star = basic_waiting_time(s, current_time);
    double v_k_1 = current_time - start_time();
    return log(1.0+m_beta*t_k_star*std::exp(-m_beta*v_k_1))/m_beta;
  }

  Event *
    GrowthEpoch::copy() const
  {
    return new GrowthEpoch(*this);
  }

  void
    GrowthEpoch::enter_callback(State &s, double current_time)
  {
    BuilderMonitor *callbacks = s.callbacks();
    if (!callbacks) return;
    callbacks->growth_callback(population(), true,
      current_time, s.total_population_size());
  }

  void
    GrowthEpoch::leave_callback(State &s, double current_time)
  {
    BuilderMonitor *callbacks = s.callbacks();
    if (!callbacks) return;
    callbacks->growth_callback(population(), false,
      current_time, s.total_population_size());
  }

  void
    GrowthEpoch::print_(std::ostream &os) const
  {
    os << "GrowthEpoch(" << population() << ", "
      << start_time() << ", " << end_time() << ", "
      << m_beta << ')';
  }

  Migration::Migration(int source, int destination,
    double migration_rate,
    double start_time, double end_time)
    : Epoch(start_time, end_time),
    m_source(source), m_destination(destination),
    m_migration_rate(migration_rate)
  {
    assert(source >= 0);
    assert(destination >= 0);
    assert(migration_rate >= 0);
    assert(end_time >= 0);
  }

  Event *
    Migration::copy() const
  {
    return new Migration(*this);
  }

  double
    Migration::nested_event_time(State &s, double current_time)
  {
    Population &src = s.populations()[m_source];
    unsigned int k = src.size();
    if (k < 1) return std::numeric_limits<double>::max();
    double rate = m_migration_rate*k/2;
    return Distribution_functions::expdev(rate);
  }

  void
    Migration::nested_update_state(Scheduler &scheduler, State &s,
    double event_time)
  {
    BuilderMonitor *callbacks = s.callbacks();
    if (callbacks)
      callbacks->migration_callback(m_source, m_destination, event_time,
        s.total_population_size());

    Population &src = s.populations()[m_source];
    Population &dst = s.populations()[m_destination];
    dst.push(src.pop_random());

  }

  void Migration::print_(std::ostream &os) const
  {
    os << "Migration(" << source() << ", " << destination() << ", "
      << migration_rate() << ", "
      << start_time() << ", " << end_time() << ')';
  }
}
