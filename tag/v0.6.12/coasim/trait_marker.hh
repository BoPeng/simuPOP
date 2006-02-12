/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__TRAIT_MARKER_HH_INCLUDED
#define CORE__TRAIT_MARKER_HH_INCLUDED

#ifndef CORE__MARKER_HH_INCLUDED
# include "marker.hh"
#endif

namespace core
{

  class TraitMarker : public Marker
  {
    public:
      TraitMarker(double position, double low_freq, double high_freq)
        : Marker(position), m_low_freq(low_freq), m_high_freq(high_freq)
        { m_values.push_back(0); m_values.push_back(1); }

      virtual Marker *copy() const;

      virtual bool run_first() const;

      virtual int default_value() const;

      virtual void add_value(int value)
      {                                           // don't add to trait markers
        throw illegal_value();
      }

      virtual Mutator *create_mutator(const Configuration &conf,
        const RetiredInterval &ri) const;

      double low_freq()  const { return m_low_freq; }
      double high_freq() const { return m_high_freq; }

      virtual const char * type() const;

  virtual std::string __repr__() const
      {
        ostringstream os;
        os << this->type() << "/" << position();
        return os.str();
      }

    private:
      double m_low_freq, m_high_freq;             // allowed range of mutation frequencies
  };

}
#endif
