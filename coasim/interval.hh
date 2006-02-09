/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__INTERVAL_HH_INCLUDED
#define CORE__INTERVAL_HH_INCLUDED

#ifndef STDEXCEPT_INCLUDED
# include <stdexcept>
# define STDEXCEPT_INCLUDED
#endif
#ifndef VECTOR_INCLUDED
# include <vector>
# define VECTOR_INCLUDED
#endif
#ifndef IOSTREAM_INCLUDED
# include <iostream>
# define IOSTREAM_INCLUDED
#endif

#include "configuration.hh"

namespace core
{

  // exception thrown if the interval is empty (or negative) -- since
  // it must be closed in one end and open in the other, it cannot be
  // a point.
  struct empty_interval : public std::logic_error
  {
    empty_interval() : std::logic_error("empty interval.") {}
  };

  // exception thrown in the interval is not in the range [0,1)
  struct interval_out_of_range : public std::logic_error
  {
    interval_out_of_range() : std::logic_error("interval out of range.") {}
  };

  /* A sub-interval on the real interval [0,1), closed on the left and
   * open on the right. */
  class Interval
  {
    public:

      Interval(double start, double end, unsigned int leaf_contacts = 0)
        throw(empty_interval,interval_out_of_range);
      ~Interval() {}

      double start()  const { return m_start; }
      double end()    const { return m_end; }
      double length() const { return m_end - m_start; }

      unsigned int leaf_contacts() const { return m_leaf_contacts; }

      // returns whether point is in the interval
      bool contains_point(double point) const;
      // returns whether i overlaps this interval
      bool overlaps(const Interval &i) const;

      void print_(std::ostream &os) const
        { os << '[' << m_start << ':' << m_end << ")<" << m_leaf_contacts << '>'; }

      bool operator == (const Interval &i) const;
      bool operator != (const Interval &i) const;

    private:
      void check_empty() const throw(empty_interval);
      void check_range() const throw(interval_out_of_range);

      double m_start;
      double m_end;
      unsigned int m_leaf_contacts;               // number of leaf nodes that this
      // interval connects to
  };

  inline bool Interval::contains_point(double point) const
    { return (start() <= point) and (point < end()); }

  inline bool Interval::overlaps(const Interval &i) const
    { return (start() <= i.start()) ? (i.start() < end()) : (start() < i.end()); }

  inline bool Interval::operator == (const Interval &i) const
  {
    return (m_start == i.m_start) and (m_end == i.m_end)
      and (m_leaf_contacts == i.m_leaf_contacts);
  }

  inline bool Interval::operator != (const Interval &i) const
    { return !(*this == i); }

  inline std::ostream & operator << (std::ostream &os, const Interval &i)
    { i.print_(os); return os; }

  // exception thrown if we try to copy an empty or inverted
  // (stop<=start) sub-interval.
  struct illegal_interval : public std::logic_error
  {
    illegal_interval() : std::logic_error("illegal interval.") {}
  };

  /* A set of non-overlapping intervals. */
  class Intervals
  {
    public:

      // add an interval to the Intervals -- the added interval must start
      // later than the previously added intervals.
      void add(const Interval &i) throw(out_of_sequence);
      void add(double start, double end, int contacts = 0)
        throw(out_of_sequence, empty_interval,
        interval_out_of_range);

      // these methods looks up the intervals NOT CHECKING if the index is
      // within range!
      const Interval& interval(int index)     const { return m_intervals[index]; }
      const Interval& operator [] (int index) const { return interval(index); }

      int size() const { return m_intervals.size(); }

      // checking the predicates on the relevant intervals.  Throws an
      // exception if point is outside [0,1).
      bool contains_point(double pos)  const throw(std::out_of_range);
      bool overlaps(const Interval &i) const throw(std::out_of_range);

      // the first point in the first interval (the left most point)
      double first_point() const throw(std::out_of_range);
      // the last point in the last interval (the right most point)
      double last_point()  const throw(std::out_of_range);

      // number of leaves reached at point -- just a propagation to the
      // corresponding interval (or 0 if intervals does not contain
      // point).
      unsigned int leaf_contacts(double point) const;

      void reset() { m_intervals.clear(); }

      // copy the intervals between start and stop, trunkating intervals
      // that overlap start and stop.
      Intervals copy(double start, double stop) const throw(illegal_interval);

      // merge this and i, splitting overlapping intervals
      Intervals merge(const Intervals& i) const;
      Intervals operator | (const Intervals &i) const;

      // This method adds two intervals where all Interval on the one
      // Intervals comes before all Interval on the second Intervals;
      // throws an exception if one interval does not come before the
      // other.
      Intervals add_intervals(const Intervals &i) const throw(out_of_sequence);
      Intervals operator + (const Intervals &i)   const throw(out_of_sequence);

      void print_(std::ostream &os) const;

    private:

      // INVARIANT: The _intervals vector contains the non-overlapping
      // intervals in sorted order, wrt to < on intervals.
      std::vector<Interval> m_intervals;

      std::vector<Interval>::const_iterator interval_starting_before(double point) const;
      std::vector<Interval>::const_iterator interval_starting_after(double point) const;
      typedef bool (Interval::*interval_predicate_t)(double point) const;
      bool check_predicate(double point, interval_predicate_t predicate) const;

      // Adds the intervals where first comes before second. If first
      // overlaps second, an exception is thrown
      static Intervals add_ordered_intervals(Intervals const &first,
        Intervals const &second)
        throw(out_of_sequence);

  };

  inline double Intervals::first_point() const throw(std::out_of_range)
  {
    if (m_intervals.size() == 0) throw std::out_of_range("no intervals!");
    return m_intervals.front().start();
  }
  inline double Intervals::last_point() const throw(std::out_of_range)
  {
    if (m_intervals.size() == 0) throw std::out_of_range("no intervals!");
    return m_intervals.back().end();
  }

  inline bool Intervals::contains_point(double point) const throw(std::out_of_range)
  {
    if (point < 0.0 or 1.0 <= point)
      throw std::out_of_range("checking point out of the [0,1) range.");
    if (size() == 0)           return false;

    if (point == 1.0)                             // special case, needed to check for endpoint in 1.0
      return (m_intervals.back().contains_point)(point);

    if (point < first_point()) return false;
    if (point > last_point())  return false;

    return interval_starting_before(point)->contains_point(point);
  }

  inline void Interval::check_empty() const throw(empty_interval)
  {
    if (length() <= 0.0) throw empty_interval();
  }

  inline void Interval::check_range() const throw(interval_out_of_range)
  {
    if ((start() < 0 or 1 <= start()) or (end() <= 0 or 1 < end()))
      throw interval_out_of_range();
  }

  inline Intervals Intervals::operator | (const Intervals &i) const
    { return merge(i); }

  inline Intervals Intervals::operator + (const Intervals &in) const
    throw(out_of_sequence)
    { return add_intervals(in); }

  inline std::ostream & operator << (std::ostream &os, const Intervals &is)
    { is.print_(os); return os; }

}                                                 // namespace core
#endif                                            // CORE__INTERVAL_HH_INCLUDED
