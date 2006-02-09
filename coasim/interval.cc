/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "interval.hh"

#ifndef VALARRAY_INCLUDED
# include <valarray>
# define VALARRAY_INCLUDED
#endif
#ifndef ALGORITHM_INCLUDED
# include <algorithm>
# define ALGORITHM_INCLUDED
#endif

#ifndef CASSERT_INCLUDED
# include <cassert>
# define CASSERT_INCLUDED
#endif

using std::min;
using std::max;
using namespace core;

namespace
{
  using std::binary_function;

  struct interval_start_less :
  public binary_function<const Interval&,const Interval&,bool>
  {
    bool operator () (const Interval &i1, const Interval &i2) const
      { return i1.start() <= i2.start(); }
  };
  struct interval_less :
  public binary_function<const Interval&,const Interval&,bool>
  {
    bool operator () (const Interval &i1, const Interval &i2) const
    {
      return (i1.start() < i2.start())
        or ((i1.start() == i2.start()) and (i1.end() < i2.end()));
    }
  };
}


core::Interval::Interval(double start, double end, unsigned int leaf_contacts)
throw(empty_interval,interval_out_of_range)
: m_start(start), m_end(end), m_leaf_contacts(leaf_contacts)
{
  check_empty();
  check_range();
}


void core::Intervals::add(const Interval &i) throw(out_of_sequence)
{
  if (m_intervals.size() == 0 or m_intervals.back().end() <= i.start())
    m_intervals.push_back(i);
  else
    throw out_of_sequence();
}


void core::Intervals::add(double start, double end, int contacts)
throw(out_of_sequence,
empty_interval, interval_out_of_range)
{
  Interval tmp(start,end,contacts);
  add(tmp);
}


// find the first interval that starts no later than point
std::vector<Interval>::const_iterator
Intervals::interval_starting_before(double point) const
{
  Interval dummy_interval(point,1.0);
  std::vector<Interval>::const_iterator itr;
  itr = lower_bound(m_intervals.begin(),m_intervals.end(), dummy_interval,
    interval_start_less());

  // itr now points to the first interval that starts at or after point
  if (itr != m_intervals.begin()) return --itr;
  else return itr;
}


// find the first interval that starts no earlier than point
std::vector<Interval>::const_iterator
Intervals::interval_starting_after(double point) const
{
  // special case to be able to handle the endpoint
  if (point == 1.0) return m_intervals.end();

  Interval dummy_interval(point,1.0);
  std::vector<Interval>::const_iterator itr;
  itr = upper_bound(m_intervals.begin(),m_intervals.end(), dummy_interval,
    interval_start_less());

  // itr now points to the first interval that starts at or after point
  if (itr != m_intervals.end() and itr->start() == point) return ++itr;
  else return itr;
}


// check the (point)-predicate around point
bool Intervals::check_predicate(double point,
interval_predicate_t predicate) const
{
  using std::bind2nd;   using std::mem_fun_ref;

  if (point == 1.0)                               // special case, needed to check for endpoint in 1.0
    return (m_intervals.back().*predicate)(point);

  if (point < 0.0 or 1.0 <= point)
    throw std::out_of_range("checking point out of the [0,1) range.");

  if (point < first_point()) return false;
  if (point > last_point())  return false;

  std::vector<Interval>::const_iterator start, stop, res;
  start = interval_starting_before(point);
  stop = interval_starting_after(point);
  res = find_if(start,stop, bind2nd(mem_fun_ref(predicate),point));
  return res != stop;
}


// FIXME: can this be handled without the two searches
// (starting_before and starting_after)?
bool Intervals::overlaps(const Interval &i) const throw(std::out_of_range)
{
  using std::bind2nd;   using std::mem_fun_ref;

  std::vector<Interval>::const_iterator start, stop, res;
  start = interval_starting_before(i.start());
  stop = interval_starting_after(i.end());

  for (res = start; res != stop; ++res)
  {
    if (res->overlaps(i)) return true;
  }
  return false;

#if 0
  res = find_if(start,stop, bind2nd(mem_fun_ref(&Interval::overlaps),i));
  return res != stop;
#endif
}


// Adds the intervals where first comes before second
Intervals Intervals::add_ordered_intervals(Intervals const &first,
Intervals const &second)
throw(out_of_sequence)
{
  if (first.m_intervals.back().overlaps(second.m_intervals.front()))
    throw out_of_sequence();

  Intervals result(first);
  std::copy(second.m_intervals.begin(), second.m_intervals.end(),
    std::back_inserter(result.m_intervals));

  return result;
}


// this operator adds two intervals where all Interval on the one
// Intervals comes before all Interval on the second Intervals
Intervals Intervals::add_intervals(const Intervals &i) const
throw(out_of_sequence)
{
  if (size() == 0)   return i;
  if (i.size() == 0) return *this;

  // At this point we know that both intervals are non-empty.  The
  // invariant of intervals gives us that they are ordered and
  // non-overlapping, if we concatenate them in the right order, the
  // invariant is still true.  If one does not not come before the
  // other, we must throw an exception

  Intervals result;

  if (m_intervals.back().start() <= i.m_intervals.front().start())
    return add_ordered_intervals(*this,i);
  else if (i.m_intervals.back().start() <= m_intervals.front().start())
    return add_ordered_intervals(i,*this);
  else
    throw out_of_sequence();

  return result;
}


unsigned int Intervals::leaf_contacts(double point) const
{
  if (point == 1.0)                               // special case for endpoint in 1.0
  {
    const Interval &i = m_intervals.back();
    if (i.contains_point(point)) return i.leaf_contacts();
    else return 0;
  }

  std::vector<Interval>::const_iterator i = interval_starting_before(point);
  if (i->contains_point(point)) return i->leaf_contacts();
  else return 0;
}


// Copy the intervals between start and stop, trunkating the
// end-intervals to start and stop.
Intervals Intervals::copy(double start, double stop) const
throw(illegal_interval)
{
  std::vector<Interval>::const_iterator first, last, itr;
  Intervals result;

  if (stop < start)  throw illegal_interval();
  if (start < 0.0)   throw illegal_interval();
  if (1.0 < stop)    throw illegal_interval();

  // empty copies...
  if (start == stop or start == 1.0 or stop == 0.0) return result;

  first = interval_starting_before(start);
  last = interval_starting_after(stop);

  // first points to the right-most interval that starts *before*
  // start, or the very first interval if no interval starts before
  // start.  last points to the left-most interval that starts *after*
  // stop or m_intervals.end()

  /* -- handle first ------------------------------------- */
  // if first contains start, we cut [start,first->end) -- or
  // [start,stop) if stop < first->end -- otherwise we just skip the
  // first interval; the next must be completely included as it must
  // start *after* start, according to the specification of
  // interval_starting_before: if it did not it would be to the right
  // of the right-most interval that starts before start

  if (start < first->end())
  {
    double s = max(first->start(),start);
    double e = min(first->end(),  stop);
    if (s < e) result.add(s, e, first->leaf_contacts());
  }

  if (first == last) return result;               // no more...

  /* -- handle the rest ---------------------------------- */
  // the only special case is the last interval, where we must make
  // sure to cut at point stop

  for (itr = first + 1; itr != last; ++itr)
  {
    if (itr->end() <= stop)
      result.add(*itr);
    else
    {
      // cut [itr->start(),stop) and add it, then terminate
      if (stop < itr->end() and itr->start() < stop)
        result.add(itr->start(),stop,itr->leaf_contacts());
      break;
    }
  }

#if 0
  std::cout << "copy of [" << start << ',' << stop << "):\n";
  for (itr = result._intervals.begin(); itr != result._intervals.end(); ++itr)
    std::cout << '[' << itr->start() << ',' << itr->end() << ")\n";
#endif

  return result;
}


Intervals Intervals::merge(const Intervals& i) const
{
  std::vector<Interval> tmp_merge;

  std::merge(m_intervals.begin(), m_intervals.end(),
    i.m_intervals.begin(), i.m_intervals.end(),
    std::back_inserter(tmp_merge),
    interval_less());

  if (tmp_merge.size() == 0) return Intervals();  // the empty merge

  std::vector<Interval> res_intervals;
  std::vector<Interval>::const_iterator itr = tmp_merge.begin();

  res_intervals.push_back(*itr);
  for (++itr; itr != tmp_merge.end(); ++itr)
  {
    if (! res_intervals.back().overlaps(*itr))
      res_intervals.push_back(*itr);              // just move on to the next
    else
    {
      // the last interval shouldn't really have been added -- it
      // overlaps the next -- so handle that
      Interval interval = res_intervals.back(); res_intervals.pop_back();

      double point0 = interval.start();
      double point1 = interval.end();
      double point2 = itr->start();
      double point3 = itr->end();

      // INVARIANT:
      assert(point0 < point1);  assert(point2 < point3);
      assert(point0 <= point2); assert(point2 < point1);
      // UNKNOWN:   point1 <= point3 or point3 < point1 ?

      if (point3 < point1)
      {
        // INVARIANT: point0 <= point2 < point3 < point1
        if (point0 < point2)
          res_intervals.push_back(Interval(point0,
            point2,
            interval.leaf_contacts()));
        res_intervals.push_back(Interval(point2,
          point3,
          interval.leaf_contacts()
          +itr->leaf_contacts()));
        res_intervals.push_back(Interval(point3,
          point1,
          interval.leaf_contacts()));
      }
      else
      {
        // INVARIANT: point0 <= point2 < point1 <= point3
        if (point0 < point2)
          res_intervals.push_back(Interval(point0,
            point2,
            interval.leaf_contacts()));
        res_intervals.push_back(Interval(point2,
          point1,
          interval.leaf_contacts()
          +itr->leaf_contacts()));
        if (point1 < point3)
          res_intervals.push_back(Interval(point1,
            point3,
            itr->leaf_contacts()));
      }
    }
  }

#if 0
  std::cout << "Merge:\n";
  for (itr = res_intervals.begin(); itr != res_intervals.end(); ++itr)
    std::cout << '[' << itr->start() << ',' << itr->end() << ")\n";
#endif

  Intervals res; res.m_intervals = res_intervals;
  return res;
}


void Intervals::print_(std::ostream &os) const
{
  std::vector<Interval>::const_iterator i;
  for (i = m_intervals.begin(); i != m_intervals.end(); ++i)
    os << *i << ' ';
}
