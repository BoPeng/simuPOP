/* -*- Mode: C++; c-basic-offset: 4; -*- 
 *
 *  CoaSim -- A coalescence process simulator
 *
 *  Copyright (C) 2004, 2005, 2006 by Bioinformatics ApS
 *                                    and Thomas Mailund <mailund@mailund.dk>
 */

#ifndef CORE__BUILDER_HH_INCLUDED
#define CORE__BUILDER_HH_INCLUDED

#ifndef VECTOR_INCLUDED
# include <vector>
# define VECTOR_INCLUDED
#endif


namespace core {

    class CoalescentNode;
    class RecombinationNode;
    class GeneConversionNode;
    class ARG;
    class Configuration;

    struct BuilderMonitor
    {
	// FIXME: should these be k for each population as well?
	// Should the population be reported?

	virtual void coalescence_callback(CoalescentNode *n,
					  int k) = 0;
	virtual void recombination_callback(RecombinationNode *n1,
					    RecombinationNode *n2,
					    int k) = 0;
	virtual void gene_conversion_callback(GeneConversionNode *n1,
					      GeneConversionNode *n2,
					      int k) = 0;

	virtual void bottleneck_callback(int pop, bool entering,
					 double time, int k) = 0;
	virtual void growth_callback(int pop, bool entering,
				     double time, int k) = 0;

	virtual void migration_callback(int pop1, int pop2,
					double time, int k) = 0;

	// FIXME: here it should *really* be ks for the two
	// populations...
	virtual void population_merge_callback(const std::vector<int> &pops,
					       double time, int k) = 0;
    };

    class Builder
    {
    public:
	Builder(const Configuration &conf) : i_conf(conf) {};
	~Builder() {};

	// Builds an ARG.  The ARG is dynamically allocated and must be
	// deleted after use.
	ARG *build(BuilderMonitor    *callbacks = 0,
		   bool keep_empty_intervals = false) const;

    private:
	const Configuration &i_conf;
    };


}

#endif
