/**
 *  $File: BaseMutator.h $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
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
   \brief head file of class mutator:public BaseOperator
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
class BaseMutator : public BaseOperator
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
	 *  Parameter \e loci can be a list of loci indexes, names or \c ALL_AVAIL.
	 *  Please refer to class \c BaseOperator for a detailed explanation of
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
	 *  alleles. For example, an \c AcgtMutator assumes alleles \c 0, \c 1,
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
	 *  <tt>mapOut=[1, 2]</tt> would allow the use of a \c SNPMutator to mutate
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
	BaseMutator(const floatList & rates = vectorf(), const lociList & loci = lociList(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		int context = 0, const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr())
		: BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_rates(rates.elems()), m_loci(loci), m_mapIn(mapIn), m_mapOut(mapOut),
		m_context(context * 2)
	{
		// NOTE: empty rates is allowed because a mutator might be
		// used in a mixed mutator.
		if (m_rates.size() > 1 && m_loci.size() == 0)
			throw ValueError("If you use variable rates, you should specify loci for each of the rate.");

		if (m_rates.size() > 1 && !m_loci.empty() && m_rates.size() != m_loci.size())
			throw ValueError("If both rates and loci are specified, they should have the same length.");
	}


	/// destructor
	virtual ~BaseMutator()
	{
	}


	/// HIDDEN Deep copy of a \c mutator
	virtual BaseOperator * clone() const
	{
		return new BaseMutator(*this);
	}


	/// CPPONLY set an array of mutation rates
	void setRate(const vectorf & rates, const lociList & loci)
	{
		if (rates.size() != 1 && rates.size() != loci.size())
			throw ValueError("If you specify more than one rate values, you should also specify corresponding applicable loci");

		m_rates = rates;
		m_loci = loci;
	}


	/// CPPONLY
	double mutRate(size_t loc) const;

	/// CPPONLY
	virtual void mutate(AlleleRef /* allele */, size_t /* locus */) const
	{
		throw SystemError("You are not supposed to call this base mutator funciton.");
	};

	/// CPPONLY
	/// Get the context of mutated allele. If an allele is invalid, -1 will
	/// be used (this is the case for the first and last loci on a chromosome).
	/// These is certainly more efficient if it is squeezed in the apply function,
	/// with a number of flags defined in the initialization stage. However, for
	/// a rarely used feature, performance should be a secondary consideration.
	void fillContext(const Population & pop, IndAlleleIterator ptr, size_t locus) const;

	/// CPPONLY
	void setContext(size_t context)
	{
		m_context.resize(context * 2);
	}


	/// CPPONLY
	vectoru & context() const
	{
		return m_context;
	}


	/// HIDDEN Apply a mutator
	virtual bool apply(Population & pop) const;

protected:
	/// This cannot be const because some mutators
	/// needs to determine these things later.

	/// mutation rates
	vectorf m_rates;

	/// applicable loci.
	lociList m_loci;

	const uintListFunc m_mapIn;

	const uintListFunc m_mapOut;

	// Be careful about this variable, which is not constant.
	mutable vectoru m_context;
};

/** A matrix mutator mutates alleles \c 0, \c 1, ..., \c n-1 using a \c n by
 *  \c n matrix, which specifies the probability at which each allele mutates
 *  to another. Conceptually speaking, this mutator goes through all mutable
 *  allele and mutate it to another state according to probabilities in the
 *  corresponding row of the rate matrix. Only one mutation rate matrix can
 *  be specified which will be used for all specified loci.
 #
 */
class MatrixMutator : public BaseMutator
{
public:
	/** Create a mutator that mutates alleles \c 0, \c 1, ..., \c n-1 using a
	 *  \c n by \c n matrix \c rate. Item <tt>(i,j)</tt> of this matrix
	 *  specifies the probability at which allele \e i mutates to allele \e j.
	 *  Diagnal items <tt>(i, i)</tt> are ignored because they are
	 *  automatically determined by other probabilities. Only one mutation rate
	 *  matrix can be specified which will be used for all loci in the applied
	 *  population, or loci specified by parameter \e loci. If alleles other
	 *  than \c 0, \c 1, ..., \c n-1 exist in the population, they will not be
	 *  mutated although a warning message will be given if debugging code
	 *  \c DBG_WARNING is turned on. Please refer to classes \c mutator and
	 *  \c BaseOperator for detailed explanation of other parameters.
	 */
	MatrixMutator(const floatMatrix & rate, const lociList & loci = lociList(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr());

	/// destructor.
	~MatrixMutator()
	{
	}


	/// CPPONLY
	virtual void mutate(AlleleRef allele, size_t locus) const;

	/// HIDDEN Deep copy of a \c MatrixMutator
	virtual BaseOperator * clone() const
	{
		return new MatrixMutator(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void) format; /* avoid warning about unused parameter */
		return "<simuPOP.MatrixMutator>";
	}


private:
	mutable vector<WeightedSampler> m_sampler;
};

/** This mutator implements a \e k-allele mutation model that assumes \e k
 *  allelic states (alleles 0, 1, 2, ..., \e k-1) at each locus. When a
 *  mutation event happens, it mutates an allele to any other states with equal
 *  probability.
 */
class KAlleleMutator : public BaseMutator
{
public:
	/** Create a k-allele mutator that mutates alleles to one of the other
	 *  <tt>k-1</tt> alleles with equal probability. This mutator by default
	 *  applies to all loci unless parameter \e loci is specified. A single
	 *  mutation rate will be used for all loci if a single value of parameter
	 *  \e rates is given. Otherwise, a list of mutation rates can be specified
	 *  for each locus in parameter \e loci. If the mutated allele is larger
	 *  than or equal to \c k, it will not be mutated. A warning message will
	 *  be displayed if debugging code \c DBG_WARNING is turned on. Please
	 *  refer to classes \c mutator and \c BaseOperator for descriptions of
	 *  other parameters.
	 */
	KAlleleMutator(UINT k, const floatList & rates = vectorf(), const lociList & loci = lociList(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr())
		: BaseMutator(rates, loci, mapIn, mapOut, 0, output, begin, end, step, at,
		              reps, subPops, infoFields), m_k(k)
	{
#ifndef BINARYALLELE
		if (m_k > 1 && static_cast<ULONG>(m_k - 1) > ModuleMaxAllele)
			throw ValueError("maxAllele exceeds Population max allele.");
#endif
	}


	~KAlleleMutator()
	{
	}


	/// CPPONLY
	virtual void mutate(AlleleRef allele, size_t locus) const;

	/// HIDDEN Deep copy of a \c KAlleleMutator
	virtual BaseOperator * clone() const
	{
		return new KAlleleMutator(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void) format; /* avoid warning about unused parameter */
		return "<simuPOP.KAlleleMutator> A k-allele mutation model with K=" +
		       toStr(m_k);
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
 */
class StepwiseMutator : public BaseMutator
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
	 *  \li <tt>(GEOMETRIC_DISTRIBUTION, p)</tt>: The number of steps follows a
	 *		a geometric distribution with parameter \e p.
	 *  \li A Python function: This user defined function accepts the allele
	 *		being mutated and return the steps to mutate.
	 *
	 *	The mutation process is usually neutral in the sense that mutating up
	 *  and down is equally likely. You can adjust parameter \e incProb to
	 *  change this behavior.
	 *
	 *  If you need to use other generalized stepwise mutation models, you can
	 *  implement them using a \c PyMutator. If performance becomes a concern,
	 *  I may add them to this operator if provided with a reliable reference.
	 */
	StepwiseMutator(const floatList & rates = vectorf(), const lociList & loci = lociList(),
		double incProb = 0.5, UINT maxAllele = 0, const floatListFunc & mutStep = floatListFunc(1.0),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(), const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr());

	~StepwiseMutator()
	{
	}


	/// mutate according to the SMM model
	/// CPPONLY
	virtual void mutate(AlleleRef allele, size_t locus) const;

	/// HIDDEN Deep copy of a \c StepwiseMutator
	virtual BaseOperator * clone() const
	{
		return new StepwiseMutator(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void) format; /* avoid warning about unused parameter */
		return "<simuPOP.StepwiseMutator> a step-wise mutation model mutator";
	}


private:
	double m_incProb;

	UINT m_maxAllele;

	floatListFunc m_mutStep;
};


/** This hybrid mutator accepts a Python function that determines how to mutate
 *  an allele when an mutation event happens.
 */
class PyMutator : public BaseMutator
{
public:
	/** Create a hybrid mutator that uses a user-provided function to mutate an
	 *  allele when a mutation event happens. This function (parameter \e func)
	 *  accepts the allele to be mutated as parameter \c allele and optional
	 *  array of alleles as parameter \c context, which are \e context alleles
	 *  the left and right of the mutated allele. Invalid context alleles (e.g.
	 *  left allele to the first locus of a chromosome) will be marked by -1.
	 *  The return value of this function will be used to mutate the passed
	 *  allele. The passed, returned and context alleles might be altered if
	 *  parameter \e mapIn and \e mapOut are used. This mutator by default
	 *  applies to all loci unless parameter \e loci is specified. A single
	 *  mutation rate will be used for all loci if a single value of
	 *  parameter \e rates is given. Otherwise, a list of mutation rates can
	 *  be specified for each locus in parameter \e loci. Please refer to
	 *  classes \c mutator and \c BaseOperator for descriptions of other
	 *  parameters.
	 */
	PyMutator(const floatList & rates = vectorf(), const lociList & loci = lociList(),
		PyObject * func = NULL, int context = 0, const uintListFunc & mapIn = uintListFunc(),
		const uintListFunc & mapOut = uintListFunc(), const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr())
		: BaseMutator(rates, loci, mapIn, mapOut, context, output, begin, end, step, at, reps, subPops, infoFields),
		m_func(func)
	{
		DBG_ASSERT(m_func.isValid(), ValueError,
			"Passed variable is not a callable python function.");
	}


	/// HIDDEN Deep copy of a \c PyMutator
	virtual BaseOperator * clone() const
	{
		return new PyMutator(*this);
	}


	/// CPPONLY
	virtual void mutate(AlleleRef allele, size_t locus) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void) format; /* avoid warning about unused parameter */
		return "<simuPOP.PyMutator>" ;
	}


private:
	pyFunc m_func;
};


/** This mixed mutator accepts a list of mutators and use one of them to mutate
 *  an allele when an mutation event happens.
 */
class MixedMutator : public BaseMutator
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
	 *  Please refer to classes \c mutator and \c BaseOperator for descriptions
	 *  of other parameters.
	 */
	MixedMutator(const floatList & rates = vectorf(), const lociList & loci = lociList(),
		const opList & mutators = opList(), const vectorf & prob = vectorf(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		int context = 0, const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr())
		: BaseMutator(rates, loci, mapIn, mapOut, context, output, begin, end, step, at, reps, subPops, infoFields),
		m_mutators(mutators), m_sampler()
	{
		DBG_FAILIF(m_mutators.size() != prob.size(), ValueError,
			"Please specify a probability for each passed mutator.");

		DBG_ASSERT(fcmp_eq(std::accumulate(prob.begin(), prob.end(), 0.), 1.), ValueError,
			"Passed probabilities should add up to 1.");

		m_sampler.set(prob.begin(), prob.end());
	}


	/// HIDDEN Deep copy of a \c MixedMutator
	virtual BaseOperator * clone() const
	{
		return new MixedMutator(*this);
	}


	/// CPPONLY
	virtual void mutate(AlleleRef allele, size_t locus) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void) format; /* avoid warning about unused parameter */
		return "<simuPOP.MixedMutator>" ;
	}


private:
	const opList m_mutators;

	mutable WeightedSampler m_sampler;
};


/** This context-dependent mutator accepts a list of mutators and use one of
 *  them to mutate an allele depending on the context of the mutated allele.
 */
class ContextMutator : public BaseMutator
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
	 *  \c BaseOperator for descriptions of other parameters.
	 */
	ContextMutator(const floatList & rates = vectorf(), const lociList & loci = lociList(),
		const opList & mutators = opList(), const intMatrix & contexts = intMatrix(),
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr())
		: BaseMutator(rates, loci, mapIn, mapOut, 0, output, begin, end, step, at, reps, subPops, infoFields),
		m_mutators(mutators), m_contexts(contexts.elems())
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


	/// HIDDEN Deep copy of a \c context-dependentMutator
	virtual BaseOperator * clone() const
	{
		return new ContextMutator(*this);
	}


	/// CPPONLY
	virtual void mutate(AlleleRef allele, size_t locus) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void) format; /* avoid warning about unused parameter */
		return "<simuPOP.ContextMutator> context-dependent mutator>" ;
	}


private:
	opList m_mutators;

	matrixi m_contexts;
};


/** A point mutator is different from all other mutators because mutations in
 *  this mutator do not happen randomly. Instead, it happens to specific loci
 *  and mutate an allele to a specific state, regardless of its original state.
 *  This mutator is usually used to introduce a mutant to a population.
 */
class PointMutator : public BaseOperator
{
public:
	/** Create a point mutator that mutates alleles at specified \e loci to
	 *  a given \e allele of individuals \e inds. If there are multiple alleles
	 *  at a locus (e.g. individuals in a diploid population), only the first
	 *  allele is mutated unless indexes of alleles are listed in parameter
	 *  \e ploidy. This operator is by default applied to individuals in the
	 *  first subpopulation but you can apply it to a different or more than
	 *  one (virtual) subpopulations using parameter *subPops* (``AllAvail`` is
	 *  also accepted). Please refer to class \c BaseOperator for detailed
	 *  descriptions of other parameters.
	 */
	PointMutator(const lociList & loci, Allele allele, const uintList & ploidy = vectoru(1, 0),
		const uintList & inds = vectoru(), const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = 0,
		const stringList & infoFields = vectorstr())
		: BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci), m_allele(allele), m_ploidy(ploidy.elems()), m_inds(inds.elems())
	{
	}


	/// destructor
	virtual ~PointMutator()
	{
	}


	/// HIDDEN Deep copy of a \c PointMutator
	virtual BaseOperator * clone() const
	{
		return new PointMutator(*this);
	}


	/// HIDDEN apply a \c PointMutator
	virtual bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void) format; /* avoid warning about unused parameter */
		return "<simuPOP.PointMutator>" ;
	}


private:
	/// applicable loci.
	lociList m_loci;
	Allele m_allele;
	vectoru m_ploidy;
	vectoru m_inds;
};

}
#endif
