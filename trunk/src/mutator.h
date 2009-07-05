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
	 *  Some mutation models are context dependent. Namely, how an allele
	 *  mutates will depend on its adjecent alleles. Whereas most simuPOP
	 *  mutators are context independent, some of them accept a parameter
	 *  \e context which is the number of alleles to the left and right of
	 *  the mutated allele. For example \e context=1 will make the alleles to
	 *  the immediate left and right to a mutated allele available to a
	 *  mutator. These alleles will be mapped in if parameter \e mapIn is
	 *  defined. How exactly a mutator makes use of these information is
	 *  mutator dependent.
	 */
	mutator(const floatList & rates = floatList(), const uintList & loci = uintList(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		int context = 0, const stringFunc & output = ">", int stage = PostMating,
		int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops,
		const stringList & infoFields = stringList())
		: baseOperator(output, stage, begin, end, step, at, reps, subPops, infoFields),
		m_rates(rates.elems()), m_loci(loci.elems()), m_mapIn(mapIn), m_mapOut(mapOut),
		m_context(context * 2), m_bt(GetRNG()), m_initialized(false)
	{
		// NOTE: empty rates is allowed because a mutator might be
		// used in a mixed mutator.
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


	/// CPPONLY
	double mutRate(UINT loc)
	{
		vectorlu::iterator it = find(m_loci.begin(), m_loci.end(), loc);

		DBG_ASSERT(it != m_loci.end(), RuntimeError,
			"Failed to find mutation rate for locus " + toStr(loc));
		DBG_ASSERT(m_rates.size() == m_loci.size(), SystemError,
			"Incorrect rate and loci size");
		return m_rates[it - m_loci.begin()];
	}


	/// CPPONLY
	virtual void mutate(AlleleRef allele, UINT locus)
	{
		throw SystemError("You are not supposed to call this base mutator funciton.");
	};

	/// CPPONLY
	/// Get the context of mutated allele. If an allele is invalid, -1 will
	/// be used (this is the case for the first and last loci on a chromosome).
	/// These is certainly more efficient if it is squeezed in the apply function,
	/// with a number of flags defined in the initialization stage. However, for
	/// a rarely used feature, performance should be a secondary consideration.
	void fillContext(const population & pop, IndAlleleIterator ptr, UINT locus);

	/// CPPONLY
	void setContext(int context)
	{
		m_context.resize(context * 2);
	}


	/// CPPONLY
	vectori & context()
	{
		return m_context;
	}


	/// Apply a mutator
	virtual bool apply(population & pop);

	/// CPPONLY initialize bernulli trial according to pop size etc
	virtual void initialize(population & pop);

protected:
	/// mutation rates
	vectorf m_rates;

	/// applicable loci.
	vectorlu m_loci;

	uintListFunc m_mapIn;

	uintListFunc m_mapOut;

	vectori m_context;

	/// bernulli trials. bitSet mutation results.
	BernulliTrials m_bt;

	/// initialized? the first apply() call will trigger an initialization process.
	bool m_initialized;
};

/** A matrix mutator mutates alleles \c 0, \c 1, ..., \c n-1 using a \c n by
 *  \c n matrix, which specifies the probability at which each allele mutates
 *  to another. Conceptually speaking, this mutator goes through all mutable
 *  allele and mutate it to another state according to probabilities in the
 *  corresponding row of the rate matrix. Only one mutation rate matrix can
 *  be specified which will be used for all specified loci.
 #
 *  <funcForm>MatrixMutate</funcForm>
 */
class matrixMutator : public mutator
{
public:
	/** Create a mutator that mutates alleles \c 0, \c 1, ..., \c n-1 using a
	 *  \c n by \c n matrix \c rate. Item <tt>(i,j)</tt> of this matrix
	 *  specifies the probability at which allele \e i mutates to allele \e j.
	 *  Diagnal items <tt>(i, i)</tt> are ignored because they are
	 *  automatically determined by other probabilities. Only one mutation rate
	 *  matrix can be specified which will be used for all loci in the applied
	 *  population, or loci specified by parameter \e loci. Please refer to
	 *  classes \c mutator and \c baseOperator for detailed explanation of
	 *  other parameters.
	 */
	matrixMutator(const matrix & rate, const uintList & loci = uintList(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		const stringFunc & output = ">", int stage = PostMating,
		int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops,
		const stringList & infoFields = stringList());

	/// destructor.
	~matrixMutator()
	{
	}


	/// CPPONLY
	virtual void mutate(AlleleRef allele, UINT locus);

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
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops, const stringList & infoFields = stringList())
		: mutator(rates, loci, mapIn, mapOut, 0, output, stage, begin, end, step, at,
		          reps, subPops, infoFields), m_k(k)
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
	virtual void mutate(AlleleRef allele, UINT locus);

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
	 *		a geometric distribution with parameter \e p.
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
		double incProb = 0.5, UINT maxAllele = 0, const floatListFunc & mutStep = floatListFunc(1),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(), const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops,
		const stringList & infoFields = stringList());

	~smmMutator()
	{
	}


	/// mutate according to the SMM model
	/// CPPONLY
	virtual void mutate(AlleleRef allele, UINT locus);

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


/** This hybrid mutator accepts a Python function that determines how to mutate
 *  an allele when an mutation event happens.
 *  <funcForm>PyMutate</funcForm>
 */
class pyMutator : public mutator
{
public:
	/** Create a hybrid mutator that uses a user-provided function to mutate an
	 *  allele when a mutation event happens. This function (parameter \e func)
	 *  accepts the allele to be mutated and return a mutated allele. If
	 *  \e context is specified, the \e context alleles to the left and to the
	 *  right of the mutated alleles will be passed to this function as the
	 *  second parameter. Invalid context alleles (e.g. left allele to the
	 *  first locus of a chromosome) will be marked by -1. The passed, returned
	 *  and context alleles might be changed if parameters \e mapIn and
	 *  \e mapOut are used although allele mappings, if needed, are usually
	 *  handled in \e func as well. This mutator by default applies to all loci
	 *  unless parameter \e loci is specified. A single mutation rate will be
	 *  used for all loci if a single value of parameter \e rates is given.
	 *  Otherwise, a list of mutation rates can be specified for each locus in
	 *  parameter \e loci. Please refer to classes \c mutator and
	 *  \c baseOperator for descriptions of other parameters.
	 */
	pyMutator(const floatList & rates = floatList(), const uintList & loci = uintList(),
		PyObject * func = NULL, int context = 0, const uintListFunc & mapIn = uintListFunc(),
		const uintListFunc & mapOut = uintListFunc(), const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops, const stringList & infoFields = stringList())
		: mutator(rates, loci, mapIn, mapOut, context, output, stage, begin, end, step, at, reps, subPops, infoFields),
		m_func(func), m_contextObj(NULL)
	{
		DBG_ASSERT(m_func.isValid(), ValueError,
			"Passed variable is not a callable python function.");
	}


	/// deep copy of a \c pyMutator
	virtual baseOperator * clone() const
	{
		return new pyMutator(*this);
	}


	/// CPPONLY
	virtual void mutate(AlleleRef allele, UINT locus);

	/// used by Python print function to print out the general information of the \c pyMutator
	virtual string __repr__()
	{
		return "<simuPOP::python mutator>" ;
	}


private:
	pyFunc m_func;

	PyObject * m_contextObj;
};


/** This mixed mutator accepts a list of mutators and use one of them to mutate
 *  an allele when an mutation event happens.
 *  <funcForm>MixedMutate</funcForm>
 */
class mixedMutator : public mutator
{
public:
	/** Create a mutator that randomly chooses one of the specified \e mutators
	 *  to mutate an allele when a mutation event happens. The mutators are
	 *  choosen according to a list of probabilities (\parameter \e prob) that
	 *  should add up to \c 1. The passed and returned alleles might be changed
	 *  if parameters \e mapIn and \e mapOut are used. Most parameters,
	 *  including \e loci, \e mapIn, \e mapOut, \e rep, and \e subPops of
	 *  mutators specified in parameter \e mutators are ignored. This mutator
	 *  by default applies to all loci unless parameter \e loci is specified.
	 *  Please refer to classes \c mutator and \c baseOperator for descriptions
	 *  of other parameters.
	 */
	mixedMutator(const floatList & rates = floatList(), const uintList & loci = uintList(),
		const opList & mutators = opList(), const vectorf & prob = vectorf(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		int context = 0, const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops, const stringList & infoFields = stringList())
		: mutator(rates, loci, mapIn, mapOut, context, output, stage, begin, end, step, at, reps, subPops, infoFields),
		m_mutators(mutators), m_sampler(GetRNG())
	{
		DBG_FAILIF(m_mutators.size() != prob.size(), ValueError,
			"Please specify a probability for each passed mutator.");

		DBG_ASSERT(fcmp_eq(std::accumulate(prob.begin(), prob.end(), 0.), 1.), ValueError,
			"Passed probabilities should add up to 1.");

		m_sampler.set(prob);
	}


	/// deep copy of a \c mixedMutator
	virtual baseOperator * clone() const
	{
		return new mixedMutator(*this);
	}


	/// CPPONLY: initialize all passed mutators
	void initialize(population & pop);


	/// CPPONLY
	virtual void mutate(AlleleRef allele, UINT locus);

	/// used by Python print function to print out the general information of the \c mixedMutator
	virtual string __repr__()
	{
		return "<simuPOP::mixed mutator>" ;
	}


private:
	opList m_mutators;

	weightedSampler m_sampler;
};


/** This context-dependent mutator accepts a list of mutators and use one of
 *  them to mutate an allele depending on the context of the mutated allele.
 *  <funcForm>ContextMutate</funcForm>
 */
class contextMutator : public mutator
{
public:
	/** Create a mutator that choose one of the specified \e mutators to mutate
	 *  an allele when a mutation event happens. The mutators are choosen
	 *  according to the context of the mutated allele, which is specified as
	 *  a list of alleles to the left and right of an allele (\parameter
	 *  \e contexts). For example, <tt>contexts=[(0,0), (0,1), (1,1)]</tt>
	 *  indicates which mutators should be used to mutate allele \c X in the
	 *  context of \c 0X0, \c 0X1, and \c 1X1. A context can include more than
	 *  one alleles at both left and right sides of a mutated allele but all
	 *  contexts should have the same (even) number of alleles. If an allele
	 *  does not have full context (e.g. when a locus is the first locus on a
	 *  chromosome), unavailable alleles will be marked as -1. There should be
	 *  a mutator for each context but an additional mutator can be specified
	 *  as the default mutator for unmatched contexts. If parameters \e mapIn
	 *  is specified, both mutated allele and its context alleles will be
	 *  mapped. Most parameters, including \e loci, \e mapIn, \e mapOut,
	 *  \e rep, and \e subPops of mutators specified in parameter \e mutators
	 *  are ignored. This mutator by default applies to all loci unless
	 *  parameter \e loci is specified. Please refer to classes \c mutator and
	 *  \c baseOperator for descriptions of other parameters.
	 */
	contextMutator(const floatList & rates = floatList(), const uintList & loci = uintList(),
		const opList & mutators = opList(), const intMatrix & contexts = intMatrix(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops, const stringList & infoFields = stringList())
		: mutator(rates, loci, mapIn, mapOut, 0, output, stage, begin, end, step, at, reps, subPops, infoFields),
		m_mutators(mutators), m_contexts(contexts)
	{
		if (m_contexts.size() != 0) {
			DBG_FAILIF(m_contexts[0].size() / 2 * 2 != m_contexts[0].size(), ValueError,
				"A context should be balanced, namely having the same number of alleles from the left and right of the mutated allele.");
			setContext(m_contexts[0].size() / 2);
		}
		for (size_t i = 1; i < m_contexts.size(); ++i) {
			DBG_FAILIF(m_contexts[i].size() != m_contexts.size(), ValueError,
				"All contexts should have the same length");
		}
		DBG_FAILIF(m_mutators.size() != m_contexts.size() && m_mutators.size() != m_contexts.size() + 1,
			ValueError,
			"Please specify a context for each passed mutator (a default mutator is allowed at the end).");
	}


	/// deep copy of a \c context-dependentMutator
	virtual baseOperator * clone() const
	{
		return new contextMutator(*this);
	}


	/// CPPONLY: initialize all passed mutators
	void initialize(population & pop);


	/// CPPONLY
	virtual void mutate(AlleleRef allele, UINT locus);

	/// used by Python print function to print out the general information of the \c context-dependentMutator
	virtual string __repr__()
	{
		return "<simuPOP::context-dependent mutator>" ;
	}


private:
	opList m_mutators;

	intMatrix m_contexts;
};


/** A point mutator is different from all other mutators because mutations in
 *  this mutator do not happen randomly. Instead, it happens to specific loci
 *  and mutate an allele to a specific state, regardless of its original state.
 *  This mutator is usually used to introduce a mutant to a population.
 *  <funcForm>PointMutate</funcForm>
 */
class pointMutator : public baseOperator
{
public:
	/** Create a point mutator that mutates alleles at specified \e loci to
	 *  a given \e allele of individuals \e inds. If there are multiple alleles
	 *  at a locus (e.g. individuals in a diploid population), only the first
	 *  allele is mutated unless indexes of alleles are listed in parameter
	 *  \e ploidy. Please refer to class \c baseOperator for detailed
	 *  descriptions of other parameters.
	 */
	pointMutator(const uintList & loci, Allele allele, const uintList & ploidy = uintList(),
		const uintList & inds = uintList(), const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops, const stringList & infoFields = stringList())
		: baseOperator(output, stage, begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci.elems()), m_allele(allele), m_ploidy(ploidy.elems()), m_inds(inds.elems())
	{
		if (m_ploidy.empty())
			m_ploidy.push_back(0);
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


private:
	/// applicable loci.
	vectorlu m_loci;
	Allele m_allele;
	vectorlu m_ploidy;
	vectorlu m_inds;
};

}
#endif
