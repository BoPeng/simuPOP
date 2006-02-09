/* -*- Mode: C++; c-basic-offset: 4; -*- 
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#include "trait_marker.hh"

#ifndef CORE__DIST_FUNCTIONS_HH_INCLUDED
# include "dist_funcs.hh"
#endif
#ifndef CORE__RETIRED_INTERVAL_HH_INCLUDED
# include "retired_interval.hh"
#endif
#ifndef CORE__NODE_HH_INCLUDED
# include "node.hh"
#endif

using namespace core;

Marker *
core::TraitMarker::copy() const
{
    return new TraitMarker(*this);
}

bool
core::TraitMarker::run_first() const
{
    return true;
}

int
core::TraitMarker::default_value() const
{
    return 0;
}

namespace
{
    class TraitMutator : public Mutator {
	const TraitMarker &i_marker;

	// allowed frequencies, translated into leaf node counts
	unsigned int i_low_leaf_count;
	unsigned int i_high_leaf_count;

	double i_mutation_point;	// the point on the "surface" where
				// the mutation sits
	double i_surface_so_far;	// the surface seen so far
    public:
	TraitMutator(const TraitMarker &marker,
		     unsigned int low_leaf_count,
		     unsigned int high_leaf_count,
		     double mutation_point);
	bool edge_has_mutation(double parent_time, double child_time);
	int  mutate_to(const Node &n, int parent_allele);
	virtual int  mutate(const Node &parent, 
			    const Node &child, 
			    int parent_allele);
    };

    TraitMutator::TraitMutator(const TraitMarker &marker,
			       unsigned int low_leaf_count, 
			       unsigned int high_leaf_count,
			       double mutation_point)
	: i_marker(marker),
	  i_low_leaf_count(low_leaf_count), i_high_leaf_count(high_leaf_count),
	  i_mutation_point(mutation_point), i_surface_so_far(0.0)
    {
    }

    bool TraitMutator::edge_has_mutation(double parent_time, double child_time) 
    {
	double edge_length = parent_time - child_time;
	bool mutate = (i_surface_so_far <= i_mutation_point
		       and
		       i_mutation_point < i_surface_so_far+edge_length);
	i_surface_so_far += edge_length;
	return mutate;
    }

    int TraitMutator::mutate_to(const Node &child, int parent_allele)
    {
	// check frequency
	unsigned int leaf_count = child.leaves_at_point(i_marker.position());
	if ((leaf_count < i_low_leaf_count) or (i_high_leaf_count < leaf_count))
	    throw Mutator::retry_arg();
	return !parent_allele;
    }

    int TraitMutator::mutate(const Node &parent, const Node &child,
			     int parent_allele)
    {
	if (edge_has_mutation(parent.time(), child.time()))
	    return mutate_to(child, parent_allele);
	else
	    return parent_allele;
    }
}

static inline void swap(unsigned int &i, unsigned int &j)
{ unsigned int tmp = i; i = j; j = tmp; }

Mutator *TraitMarker::create_mutator(const Configuration   &conf,
				     const RetiredInterval &ri) const
{
    unsigned int low_leaf_count
	= static_cast<unsigned int>(ceil(i_low_freq*conf.no_leaves()));
    unsigned int high_leaf_count
	= static_cast<unsigned int>(floor(i_high_freq*conf.no_leaves()));

    // possible due to ceil/floor
    if (high_leaf_count < low_leaf_count) swap(low_leaf_count,high_leaf_count);

    double mutation_point = ri.surface() * Distribution_functions::uniform();

    return new TraitMutator(*this,low_leaf_count,high_leaf_count,mutation_point);
}

const char *
core::TraitMarker::type() const
{
    return "trait";
}
