/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__MARKER_HH_INCLUDED
#define CORE__MARKER_HH_INCLUDED

#ifndef STDEXCEPT_INCLUDED
# include <stdexcept>
# define STDEXCEPT_INCLUDED
#endif
#ifndef VECTOR_INCLUDED
# include <vector>
# define VECTOR_INCLUDED
#endif
#ifndef CASSERT_INCLUDED
# include <cassert>
# define CASSERT_INCLUDED
#endif

#include <sstream>
using std::ostringstream;

namespace core
{

  class Configuration;
  class RetiredInterval;
  class Node;

  // Abstract class for mutating the ARG
  class Mutator
  {
    public:
      // Exception thrown if a mutation should be re-tried (for example,
      // because the mutation frequency is outside the desired range
      struct retry_mutation : public std::exception
      {
        retry_mutation()
        {
        }
      };

      // Exception thrown if the entire ARG should be rebuilt for another
      // try (for example, because the trait mutation frequency is outside
      // the desired range
      struct retry_arg : public std::exception
      {
        retry_arg()
        {
        }
      };

      virtual int mutate( const Node &parent, const Node &child, int parent_allele) = 0;
  };

  // Exception thrown if we try to add a value to a value set that
  // doesn't fit the type of the value set
  struct illegal_value : public std::logic_error
  {
    illegal_value() : std::logic_error("illegal marker value.") {}
  };

  // Exception thrown if a marker is given a position outside [0,1).
  struct illegal_position : public std::logic_error
  {
    illegal_position() : std::logic_error("illeval position") {}
  };

  // Abstract class for the different possible marker types
  class Marker
  {
    public:

      // polymorphic copying
      virtual Marker *copy() const = 0;

virtual ~Marker()
      {
      }

      double position() const
      {
        return m_position;
      }

      void position(double position)
      {
        if (position < 0 or 1 <= position)
          throw illegal_position();

        m_position = position;
      }

      virtual bool run_first() const = 0;

      virtual int default_value() const = 0;

      // creates a new mutator -- the mutator must be deleted after use.
      virtual Mutator *create_mutator(const Configuration &conf,
        const RetiredInterval &ri) const = 0;

      virtual std::string __repr__() const =0;
 
      virtual const char * type() const = 0;

    protected:

      Marker(double position)
        : m_position(position)
      {
        if (position < 0 or 1 <= position)
          throw illegal_position();
      };

      double m_position;

      std::vector<int> m_values;

      Marker(const Marker &other)
    : m_position(other.m_position), m_values(other.m_values)
  {
  }

    private:

      // Disable
      Marker &operator = (const Marker&);
  };

}
#endif
