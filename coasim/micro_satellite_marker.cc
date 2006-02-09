/* -*- Mode: C++; c-basic-offset: 4; -*- 
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "micro_satellite_marker.hh"

#ifndef CORE__DIST_FUNCTIONS_HH_INCLUDED
# include "dist_funcs.hh"
#endif
#ifndef CORE__NODE_HH_INCLUDED
# include "node.hh"
#endif


using namespace core;

Marker *
core::MicroSatelliteMarker::copy() const
{
    return new MicroSatelliteMarker(*this);
}

bool
core::MicroSatelliteMarker::run_first() const
{
    return false;
}

int
core::MicroSatelliteMarker::default_value() const 
{
    return 0;
}


namespace
{
    class MicroSatelliteMutator : public Mutator {
    public:
	MicroSatelliteMutator(const MicroSatelliteMarker &marker);
	bool edge_has_mutation(double parent_time, double child_time);
	int  mutate_to(const Node &n, int parent_allele)
	    throw (retry_mutation, retry_arg);
	virtual int  mutate(const Node &parent, 
			    const Node &child, 
			    int parent_allele);

    private:
	const MicroSatelliteMarker &i_marker;
    };

    MicroSatelliteMutator::MicroSatelliteMutator(const MicroSatelliteMarker &marker)
	: i_marker(marker)
    {
    }

    bool MicroSatelliteMutator::edge_has_mutation(double parent_time, 
						  double child_time) 
    {
	using namespace Distribution_functions;
	double time = parent_time - child_time;
	return uniform() < expdist(i_marker.theta()/2,time);
    }
    
    int MicroSatelliteMutator::mutate_to(const Node &n, 
					 int parent_allele)
	throw (retry_mutation, retry_arg) 
    {
	using namespace Distribution_functions;
	int new_value = irand(i_marker.K()-1);
	if (new_value == parent_allele) new_value = i_marker.K()-1;
	return new_value;
    }

    int MicroSatelliteMutator::mutate(const Node &parent, const Node &child,
				      int parent_allele)
    {
	if (edge_has_mutation(parent.time(), child.time()))
	    return mutate_to(parent, parent_allele);
	else
	    return parent_allele;
    }
}

Mutator *
core::MicroSatelliteMarker::create_mutator(const Configuration   &conf,
					   const RetiredInterval &ri) const
{

    return new MicroSatelliteMutator(*this);
}


const char *
core::MicroSatelliteMarker::type() const
{
    return "ms";
}
