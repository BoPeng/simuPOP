/**
 *  $File: mutator.h $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MUTATOR_H
#define _MUTATOR_H
/**
   \file
   \brief head file of class mutator:public baseOperator
 */
/// for hybrid mutator
#include "operator.h"

namespace simuPOP {

/** Class \c mutator is the base class of all mutators. It handles all the work
 *  of picking an allele at specified loci from certain (virtual) subpopulation
 *  with certain probability, and calling a derived mutator to mutate the
 *  allele. Alleles can be changed before and after mutation if existing allele
 *  numbers do not match those of a mutation model.
 */
class mutator : public baseOperator
{
public:
	/** A mutator mutates alleles from one state to another with given
	 *  probability. This base mutator does not perform any mutation but it
	 *  defines common behaviors of all mutators.
	 *
	 *  By default, a mutator mutates all alleles in all populations of a
	 *  simulator at all generations. A number of parameters can be used to
	 *  restrict mutations to certain generations (parameters \e begin, \e end,
	 *  \e step and \e at), replicate populations (parameter \e rep), (virtual)
	 *  subpopulations (parameter \e subPops) and loci (parameter \e loci).
	 *  Please refer to class \c baseOperator for a detailed explanation of
	 *  these parameters.
	 *
	 *  Parameter \e rate or its equivalence specifies the probability that a
	 *  a mutation event happens. The exact form and meaning of \e rate is
	 *  mutator-specific. If a single rate is specified, it will be applied to
	 *  all \e loci. If a list of mutation rates are given, they will be applied
	 *  to each locus specified in parameter \e loci. Note that not all mutators
	 *  allow specification of multiple mutation rate, especially when the
	 *  mutation rate itself is a list or matrix.
	 *
	 *  Alleles at a locus are non-negative numbers 0, 1, ... up to the maximum
	 *  allowed allele for the loaded module (1 for binary, 255 for short and
	 *  65535 for long modules). Whereas some general mutation models treat
	 *  alleles as numbers, other models assume specific interpretation of
	 *  alleles. For example, an \c acgtMutator assumes alleles \c 0, \c 1,
	 *  \c 2 and \c 3 as nucleotides \c A, \c C, \c G and \c T. Using a mutator
	 *  that is incompatible with your simulation will certainly yield erroneous
	 *  results.
	 *
	 *  If your simulation assumes different alleles with a mutation model, you
	 *  can map an allele to the allele used in the model and map the mutated
	 *  allele back. This is achieved using a \e mapIn list with its \c i-th
	 *  item being the corresponding allele of real allele \c i, and a
	 *  \e mapOut list with its \e i-th item being the real allele of allele
	 *  \c i assumed in the model. For example <tt>mapIn=[0, 0, 1]</tt> and
	 *  <tt>mapOut=[1, 2]</tt> would allow the use of a \c snpMutator to mutate
	 *  between alleles 1 and 2, instead of 0 and 1. Parameters \e mapIn and
	 *  \e mapOut also accept a user-defined Python function that returns
	 *  a corresponding allele for a given allele. This allows easier mapping
	 *  between a large number of alleles and advanced models such as random
	 *  emission of alleles.
	 *
	 *  A mutator keeps track of number of mutation events happens at each
	 *  locus. This count does not have to correspond to number of new mutants
	 *  because some mutation events do not create new mutants.
	 */
	mutator(const floatList & rates = floatList(), const uintList & loci = uintList(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		const stringFunc & output = ">", int stage = PostMating,
		int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList())
		: baseOperator(output, stage, begin, end, step, at, rep, subPops, infoFields),
		m_rates(rates.elems()), m_loci(loci.elems()), m_mapIn(mapIn), m_mapOut(mapOut),
		m_bt(rng()), m_initialized(false), m_mutCount(0)
	{
		if (m_rates.empty() )
			throw ValueError("You should specify a rate, or a sequence of rate.");

		if (m_rates.size() > 1 && m_loci.empty())
			throw ValueError("If you use variable rates, you should specify loci for each of the rate.");

		if (m_rates.size() > 1 && !m_loci.empty() && m_rates.size() != m_loci.size() )
			throw ValueError("If both rates and loci are specified, they should have the same length.");
	}


	/// destructor
	virtual ~mutator()
	{
	}


	/// deep copy of a \c mutator
	virtual baseOperator * clone() const
	{
		return new mutator(*this);
	}


	/// CPPONLY set an array of mutation rates
	void setRate(const vectorf & rates, const vectorlu & loci = vectorlu())
	{
		if (rates.size() != 1 && rates.size() != loci.size() )
			throw ValueError("If you specify more than one rate values, you should also specify corresponding applicable loci");

		m_rates = rates;
		if (!loci.empty())
			m_loci = loci;

		m_initialized = false;
	}


	/** Return number of mutation events at all loci, including loci that are
	 *  listed in parameter \e loci. Note that not all mutation events leads
	 *  to a new mutant.
	 */
	vectoru mutationCounts()
	{
		return m_mutCount;
	}


	/// CPPONLY
	virtual void mutate(AlleleRef allele)
	{
		throw SystemError("You are not supposed to call this base mutator funciton.");
	};

	/// Apply a mutator
	virtual bool apply(population & pop);

private:
	/// initialize bernulli trial according to pop size etc
	virtual void initialize(population & pop);

protected:
	/// mutation rates
	vectorf m_rates;

	/// applicable loci.
	vectorlu m_loci;

	uintListFunc m_mapIn;

	uintListFunc m_mapOut;

	/// bernulli trials. bitSet mutation results.
	BernulliTrials m_bt;

	/// initialized? the first apply() call will trigger an initialization process.
	bool m_initialized;

	/// report the number of mutation events
	vectoru m_mutCount;
};

/** A matrix mutator mutates alleles \c 0, \c 1, ..., \c n-1 using a \c n by
 *  \c n matrix, which specifies the probability at which each allele mutates
 *  to another. Conceptually speaking, this mutator goes through all mutable
 *  allele and mutate it to another state according to probabilities
 *  \f$p_{i0}$, \f$p_{i1}$, ... and \f$p_{i,n-1}$. Most alleles will not mutate
 *  because \f$p_{ii}$ is usually close to 1. Only one mutation rate matrix
 *  can be specified which will be used for all specified loci.
 #
 *  <funcForm>MatrixMutate</funcForm>
 */
class matrixMutator : public mutator
{
public:
	/** Create a mutator that mutates alleles \c 0, \c 1, ..., \c n-1 using a
	 *  \c n by \c n matrix \c rate. Item \f$p_{ij}$ of this matrix specifies
	 *  the probability at which allele \e i mutates to allele \e j. \f$p_{ii}$
	 *  are ignored because they are automatically determined by \f$p_{ij}$,
	 *  \f$j=0,...,n-1$. Only one mutation rate matrix can be specified which
	 *  will be used for all loci in the applied population, or loci specified
	 *  by parameter \e loci. Please refer to classes \c mutator and
	 *  \c baseOperator for detailed explanation of other parameters.
	 */
	matrixMutator(const matrix & rate, const uintList & loci = uintList(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		const stringFunc & output = ">", int stage = PostMating,
		int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList());

	/// destructor.
	~matrixMutator()
	{
	}


	/// CPPONLY
	virtual void mutate(AlleleRef allele);

	/// deep copy of a \c matrixMutator
	virtual baseOperator * clone() const
	{
		return new matrixMutator(*this);
	}


	///
	virtual string __repr__()
	{
		return "<simuPOP::matrix mutator";
	}


private:
	vector<weightedSampler> m_sampler;
};

/** This mutator implements a \e k-allele mutation model that assumes \e k
 *  allelic states (alleles 0, 1, 2, ..., \e k-1) at each locus. When a
 *  mutation event happens, it mutates an allele to any other states with equal
 *  probability.
 *  <funcForm>KamMutate</funcForm>
 */
class kamMutator : public mutator
{
public:
	/** Create a k-allele mutator that mutates alleles to one of the other
	 *  <tt>k-1</tt> alleles with equal probability. If parameter \e k is
	 *  unspecified, it will be assumed to be the number of allowed allelic
	 *  states of the current simuPOP module. This mutator by default applies
	 *  to all loci unless parameter \e loci is specified. A single mutation
	 *  rate will be used for all loci if a single value of parameter \e rates
	 *  is given. Otherwise, a list of mutation rates can be specified for each
	 *  locus in parameter \e loci. Please refer to classes \c mutator and
	 *  \c baseOperator for descriptions of other parameters.
	 */
	kamMutator(UINT k = 0, const floatList & rates = floatList(), const uintList & loci = uintList(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const stringList & infoFields = stringList())
		: mutator(rates, loci, mapIn, mapOut, output, stage, begin, end, step, at,
		          rep, subPops, infoFields), m_k(k)
	{
		if (m_k == 0)
			m_k = MaxAllele() + 1;
#ifndef BINARYALLELE
		if (m_k > MaxAllele() )
			throw ValueError("maxAllele exceeds population max allele.");
#endif
	}


	~kamMutator()
	{
	}


	/// CPPONLY
	virtual void mutate(AlleleRef allele);

	/// deep copy of a \c kamMutator
	virtual baseOperator * clone() const
	{
		return new kamMutator(*this);
	}


	/// used by Python print function to print out the general information of the \c kamMutator
	virtual string __repr__()
	{
		return "<simuPOP::k-allele model mutator K=" +
		       toStr(m_k) + ">" ;
	}


private:
	UINT m_k;

};


/** A stepwise mutation model treats alleles at a locus as the number of
 *  tandem repeats of microsatellite or minisatellite markers. When a mutation
 *  event happens, the number of repeats (allele) either increase or decrease.
 *  A standard stepwise mutation model increases of decreases an allele by 1
 *  with equal probability. More complex models (generalized stepwise mutation
 *  model) are also allowed. Note that an allele cannot be mutated beyond
 *  boundaries (0 and maximum allowed allele).
 *  <funcForm>SmmMutate</funcForm>
 */
class smmMutator : public mutator
{
public:
	/** Create a stepwise mutation mutator that mutates an allele by increasing
	 *  or decreasing it. This mutator by default applies to all loci unless
	 *  parameter \e loci is specified. A single mutation rate will be used for
	 *  all loci if a single value of parameter \e rates is given. Otherwise, a
	 *  list of mutation rates can be specified for each locus in parameter
	 *  \e loci.
	 *
	 *  When a mutation event happens, this operator increases or decreases an
	 *  allele by \e mutStep steps. Acceptable input of parameter \e mutStep
	 *  include
	 *  \li A number: This is the default mode with default value 1.
	 *  \li <tt>(GeometricDistribution, p)</tt>: The number of steps follows a
	 *		a geometric distribution with parameter \e p. The mean and variance
	 *		of steps are \f$\frac{p}{1-p}\f$ and 
	 *		\f$\frac{p}{\left(1-p\right)^{2}}\f$ respectively.
	 *  \li A Python function: This user defined function accepts the allele
	 *		being mutated and return the steps to mutate.
	 *
	 *	The mutation process is usually neutral in the sense that mutating up
	 *  and down is equally likely. You can adjust parameter \e incProb to
	 *  change this behavior. 
	 *
	 *  If you need to use other generalized stepwise mutation models, you can
	 *  implement them using a \c pyMutator. If performance becomes a concern,
	 *  I may add them to this operator if provided with a reliable reference.
	 */
	smmMutator(const floatList & rates = floatList(), const uintList & loci = uintList(),
		double incProb=0.5, UINT maxAllele=0, const floatListFunc & mutStep = floatListFunc(1),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(), const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList());

	~smmMutator()
	{
	}


	/// mutate according to the SMM model
	/// CPPONLY
	virtual void mutate(AlleleRef allele);

	/// deep copy of a \c smmMutator
	virtual baseOperator * clone() const
	{
		return new smmMutator(*this);
	}


	/// used by Python print function to print out the general information of the \c smmMutator
	virtual string __repr__()
	{
		return "<simuPOP::step-wise mutation model mutator>" ;
	}


private:
	double m_incProb;

	UINT m_maxAllele;

	floatListFunc m_mutStep;
};


/// A hybrid mutator
/**
   Parameters such as mutation rate of this operator are set just like others and you are supposed to
   provide a Python function to return a new allele state given an old state. \c pyMutator
   will choose an allele as usual and call your function to mutate it to another allele.
   <funcForm>PyMutate</funcForm>
 */
class pyMutator : public mutator
{
public:
	/// create a \c pyMutator
	/**
	 */
	pyMutator(const floatList & rate = floatList(), const uintList & loci = uintList(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		PyObject * func = NULL,
		const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const stringList & infoFields = stringList())
		: mutator(rate, loci, mapIn, mapOut,
		          output, stage, begin, end, step, at, rep, subPops, infoFields),
		m_func(func)
	{
		DBG_ASSERT(m_func.isValid(), ValueError,
			"Passed variable is not a callable python function.");
	}


	/// deep copy of a \c pyMutator
	virtual baseOperator * clone() const
	{
		return new pyMutator(*this);
	}


	/// mutate according to the mixed model
	virtual void mutate(AlleleRef allele);

	/// used by Python print function to print out the general information of the \c pyMutator
	virtual string __repr__()
	{
		return "<simuPOP::python mutator>" ;
	}


private:
	pyFunc m_func;
};

/// point mutator
/**
   Mutate specified individuals at specified loci to a spcified allele.
   I.e., this is a non-random mutator used to introduce diseases etc.
   \c pointMutator, as its name suggest, does point mutation. This mutator will turn
   alleles at \c loci on the first chromosome copy to \c toAllele for individual \c inds.
   You can specify \c atPloidy to mutate other, or all ploidy copies.

   <funcForm>PointMutate</funcForm>
 */
class pointMutator : public baseOperator
{
public:
	/// create a \c pointMutator
	/**
	   \param inds individuals who will mutate
	   \param toAllele allele that will be mutate to

	   Please see class \c mutator for the descriptions of other parameters.
	 */
	pointMutator(const uintList & loci, Allele toAllele,
		vectoru atPloidy = vectoru(),
		vectorlu inds = vectorlu(),
		const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const stringList & infoFields = stringList())
		: baseOperator(output, stage, begin, end, step, at, rep, subPops, infoFields),
		m_loci(loci.elems()), m_toAllele(toAllele),
		m_atPloidy(atPloidy), m_inds(inds), m_mutCount(0)
	{
		if (m_atPloidy.empty())
			m_atPloidy.push_back(0);
	}


	/// destructor
	virtual ~pointMutator()
	{
	}


	/// deep copy of a \c pointMutator
	virtual baseOperator * clone() const
	{
		return new pointMutator(*this);
	}


	/// apply a \c pointMutator
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c pointMutator
	virtual string __repr__()
	{
		return "<simuPOP::point mutator>" ;
	}


	/// return mutation count at \c locus
	ULONG mutationCount(size_t locus)
	{
		DBG_ASSERT(locus < m_mutCount.size(), IndexError,
			"locus index " + toStr(locus) + " is out of range");
		return m_mutCount[locus];
	}


	/// return mutation counts
	vectoru mutationCounts()
	{
		return m_mutCount;
	}


private:
	/// applicable loci.
	vectorlu m_loci;
	Allele m_toAllele;
	vectoru m_atPloidy;
	vectorlu m_inds;
	/// report the number of mutation events
	vectoru m_mutCount;
};

}
#endif
