/* -*- Mode: C++; c-basic-offset: 4; -*- 
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__SNP_MARKER_HH_INCLUDED
#define CORE__SNP_MARKER_HH_INCLUDED

#ifndef CORE__MARKER_HH_INCLUDED
# include "marker.hh"
#endif

namespace core {

    class SNPMarker : public Marker
    {
    public:
	SNPMarker(double position, double low_freq, double high_freq) 
	    : Marker(position), i_low_freq(low_freq), i_high_freq(high_freq)
	{ i_values.push_back(0); i_values.push_back(1); }

	virtual Marker *copy() const;

	virtual bool run_first() const;

	virtual int default_value() const;

	virtual void add_value(int value) throw(illegal_value)
	{ throw illegal_value(); } // don't add to SNP markers

	virtual Mutator *create_mutator(const Configuration &conf,
					const RetiredInterval &ri) const;

	double low_freq()  const { return i_low_freq; }
	double high_freq() const { return i_high_freq; }

	virtual const char * type() const;

    private:
	double i_low_freq, i_high_freq; // allowed range of mutation frequencies
    };

}

#endif
