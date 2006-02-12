/* -*- Mode: C++; c-basic-offset: 4; -*-
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__MICRO_SATELLITE_MARKER_HH_INCLUDED
#define CORE__MICRO_SATELLITE_MARKER_HH_INCLUDED

#ifndef CORE__MARKER_HH_INCLUDED
# include "marker.hh"
#endif

namespace core
{

  class MicroSatelliteMarker : public Marker
  {
    public:
      MicroSatelliteMarker(double position, double theta, int K)
        : Marker(position), m_theta(theta), m_K(K)
      {
        assert(m_K > 0);                          // FIXME: exception?
      }

      virtual Marker *copy() const;

      virtual bool run_first() const;

      virtual int default_value() const;

      virtual Mutator *create_mutator(const Configuration &conf,
        const RetiredInterval &ri) const;

      double theta() const
      {
        return m_theta;
      }

      int K() const
      {
        return m_K;
      }

      virtual const char * type() const;
      
  virtual std::string __repr__() const
      {
        ostringstream os;
        os << this->type() << "/" << position();
        return os.str();
      }

    private:
      // mutation parameter
      double m_theta;

      int m_K;
  };

}
#endif
