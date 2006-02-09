/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__EPOCHS_HH_INCLUDED
#define CORE__EPOCHS_HH_INCLUDED

// the abstract epoch class is defined in the configuration, since it
// is a configuration thingy.
#ifndef CORE__CONFIGURATION_HH_INCLUDED
# include "configuration.hh"
#endif
#ifndef CORE__BUILDER_EVENTS_HH_INCLUDED
# include "builder_events.hh"
#endif

namespace core
{

  class Epoch : public Event
  {
    double i_start_time, i_end_time;

    // override this one to provide the epoch
    virtual double nested_event_time  (State &s, double current_time) = 0;
    virtual void   nested_update_state(Scheduler &scheduler, State &s,
      double event_time) = 0;

    public:
      // Steals the epoch -- will delete it when the epoch ends!
      Epoch(double start_time, double end_time)
        : i_start_time(start_time),
        i_end_time(end_time > 0 ? end_time
        : std::numeric_limits<double>::max())
      {
        assert(i_start_time >= 0);
        assert(i_start_time < i_end_time);
      }
      virtual ~Epoch();

      double start_time() const { return i_start_time; }
      double end_time()   const { return i_end_time; }
      virtual double earliest_event() const;

      virtual double event_time  (State &s, double current_time);
      virtual void   update_state(Scheduler &scheduler, State &s,
        double event_time);
  };

  class Migration : public Epoch
  {
    int i_source, i_destination;
    double i_migration_rate;

    virtual double nested_event_time  (State &s, double current_time);
    virtual void   nested_update_state(Scheduler &scheduler, State &s,
      double event_time);

    public:
      Migration(int source, int destination,
        double migration_rate,
        double start_time, double end_time);

      int    source()         const { return i_source; }
      int    destination()    const { return i_destination; }
      double migration_rate() const { return i_migration_rate; }

      virtual Event *copy() const;
      virtual void print_(std::ostream &os) const;
  };

  class CoalescenceEpoch : public CoalescenceEvent
  {
    double i_start_time, i_end_time;
    CoalescenceEvent *i_underlying;

    protected:
      double basic_waiting_time(State &s, double current_time)
      {
        assert(i_underlying);
        return i_underlying->waiting_time(s, current_time);
      }

      // override this to implement the change in waiting time for
      // this epoch.
      virtual double waiting_time(State &s, double current_time) = 0;

      // these are overridden to compose the waiting times and
      // perform a basic coalescence event
      virtual double nested_event_time(State &s, double current_time);
      virtual void   nested_update_state(Scheduler &scheduler, State &s,
        double event_time);

      // handles the push/pop nature of coalescence epochs
      virtual double event_time  (State &s, double current_time);
      virtual void   update_state(Scheduler &scheduler, State &s,
        double event_time);

    public:
      CoalescenceEpoch(int population, double start_time, double end_time)
        : CoalescenceEvent(population),
        i_start_time(start_time),
        i_end_time((end_time > 0) ? end_time
        : std::numeric_limits<double>::max()),
        i_underlying(0)
      {
        assert(i_start_time >= 0);
        assert(i_start_time < i_end_time);
      }
      virtual ~CoalescenceEpoch();

      double start_time() const { return i_start_time; }
      double end_time()   const { return i_end_time; }
      virtual double earliest_event() const;

      virtual void enter_callback(State &s, double current_time);
      virtual void leave_callback(State &s, double current_time);
  };

  class BottleNeck : public CoalescenceEpoch
  {
    double i_scale_fraction;
    virtual double waiting_time(State &s, double current_time);

    public:
      BottleNeck(int population, double scale_fraction,
        double start_time, double end_time = -1)
        : CoalescenceEpoch(population, start_time, end_time),
        i_scale_fraction(scale_fraction)
      {
        assert(scale_fraction > 0);
      }

      double scale_fraction() const { return i_scale_fraction; }

      virtual Event *copy() const;

      virtual void enter_callback(State &s, double current_time);
      virtual void leave_callback(State &s, double current_time);

      virtual void print_(std::ostream &os) const;
  };

  class Growth : public CoalescenceEpoch
  {
    double i_beta;

    virtual double waiting_time(State &s, double current_time);

    public:
      Growth(int population, double beta,
        double start_time, double end_time = -1)
        : CoalescenceEpoch(population, start_time, end_time),
        i_beta(beta)
      {
        assert(beta > 0);
      }

      double beta()        const { return i_beta; }

      virtual Event *copy() const;

      virtual void enter_callback(State &s, double current_time);
      virtual void leave_callback(State &s, double current_time);
      virtual void print_(std::ostream &os) const;
  };

}
#endif                                            // CORE__EPOCHS_HH_INCLUDED
