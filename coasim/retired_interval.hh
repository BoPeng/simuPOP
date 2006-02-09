/* -*- Mode: C++; c-basic-offset: 4; -*- 
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__RETIRED_INTERVAL_HH_INCLUDED
#define CORE__RETIRED_INTERVAL_HH_INCLUDED

#ifndef CORE__INTERVAL_HH_INCLUDED
# include "interval.hh"
#endif

namespace core {

    class Node;
    class Configuration;

    // -- Intervals that are retired because they connect to all leaves -----
    class RetiredInterval : public Interval
    {
    public:
	struct null_top_node : public std::logic_error {
	    null_top_node() : std::logic_error("null top node") {}
	};
  
	RetiredInterval(const Interval &interval, Node *const top_node)
	    throw(null_top_node)
	    : Interval(interval), i_surface(-1.0), i_top_node(top_node)
	{ if (top_node == 0) throw null_top_node(); }
  
	Node *top_node() const { return i_top_node; }
	double surface() const
	{ if (i_surface < 0.0) calc_surface(); return i_surface; }

	void mutate(const Configuration &conf, unsigned int marker_index) const;
  
	void to_xml(std::ostream &os) const;
  
    private:
	void calc_surface() const;

	mutable double i_surface;
	Node *i_top_node;
    };


    inline std::ostream & operator << (std::ostream &os, const RetiredInterval &i)
    { i.print(os); return os; }

}

#endif // RETIRED_INTERVAL_HH_INCLUDED
