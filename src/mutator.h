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
	 *  Parameter \e loci can be a list of loci indexes, names, list
	 *  of chromosome position pairs, \c ALL_AVAIL, or a function with optional
	 *  parameter \c pop that will be called at each ganeeration to determine
	 *  indexes of loci. Please refer to class \c BaseOperator for a detailed
	 *  explanation of these parameters.
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
	 *  If a valid information field is specified for parameter \e infoFields
	 *  (default to \c ind_id) for modules with lineage allele type, the lineage
	 *  of the mutated alleles will be the ID (stored in the first field of
	 *  \e infoFields) of individuals that harbor the mutated alleles if \e
	 *  lineageMode is set to \c FROM_INFO (default). If \e lineageMode is
	 *  set to \c FROM_INFO_SIGNED, the IDs will be assigned a sign depending
	 *  on the ploidy the mutation happens (1 for ploidy 0, -1 for ploidy 1,
	 *  etc). The lineage information will be transmitted along with the
	 *  alleles so this feature allows you to track the source of mutants
	 *  during evolution.A
	 *
	 *  A mutator by default does not produce any output. However, if
	 *  an non-empty output is specified, the operator will output generation
	 *  number, locus, ploidy, original allele, mutant, and values of all
	 *  information field specified by parameter \c infoFields (e.g. individual
	 *  ID if \c ind_id is specified).
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
		int context = 0, const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "ind_id"), int lineageMode = FROM_INFO)
		: BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_rates(rates.elems()), m_loci(loci), m_mapIn(mapIn), m_mapOut(mapOut),
		m_lineageMode(lineageMode), m_context(context * 2)
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
	virtual Allele mutate(Allele /* allele */, size_t /* locus */) const
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

	int m_lineageMode;

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
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "ind_id"), int lineageMode = FROM_INFO);

	/// destructor.
	~MatrixMutator()
	{
	}


	/// CPPONLY
	virtual Allele mutate(Allele allele, size_t locus) const;

	/// HIDDEN Deep copy of a \c MatrixMutator
	virtual BaseOperator * clone() const
	{
		return new MatrixMutator(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
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
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "ind_id"), int lineageMode = FROM_INFO)
		: BaseMutator(rates, loci, mapIn, mapOut, 0, output, begin, end, step, at,
		              reps, subPops, infoFields, lineageMode), m_k(k)
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
	virtual Allele mutate(Allele allele, size_t locus) const;

	/// HIDDEN Deep copy of a \c KAlleleMutator
	virtual BaseOperator * clone() const
	{
		return new KAlleleMutator(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return (boost::format("<simuPOP.KAlleleMutator> A k-allele mutation model with K=%1%") % m_k).str();
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
		const uintListFunc & mapIn = uintListFunc(), const uintListFunc & mapOut = uintListFunc(),
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "ind_id"), int lineageMode = FROM_INFO);

	~StepwiseMutator()
	{
	}


	/// mutate according to the SMM model
	/// CPPONLY
	virtual Allele mutate(Allele allele, size_t locus) const;

	/// HIDDEN Deep copy of a \c StepwiseMutator
	virtual BaseOperator * clone() const
	{
		return new StepwiseMutator(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
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
	 *  accepts the allele to be mutated as parameter \c allele, locus index
	 *  \c locus, and optional array of alleles as parameter \c context, which
	 *  are \e context alleles the left and right of the mutated allele. Invalid
	 *  context alleles (e.g. left allele to the first locus of a chromosome)
	 *  will be marked by -1. The return value of this function will be used
	 *  to mutate the passed allele. The passed, returned and context alleles
	 *  might be altered if parameter \e mapIn and \e mapOut are used. This
	 *  mutator by default applies to all loci unless parameter \e loci is
	 *  specified. A single mutation rate will be used for all loci if a
	 *  single value of parameter \e rates is given. Otherwise, a list of
	 *  mutation rates can be specified for each locus in parameter \e loci.
	 *  Please refer to classes \c mutator and \c BaseOperator for descriptions
	 *  of other parameters.
	 */
	PyMutator(const floatList & rates = vectorf(), const lociList & loci = lociList(),
		PyObject * func = NULL, int context = 0, const uintListFunc & mapIn = uintListFunc(),
		const uintListFunc & mapOut = uintListFunc(), const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "ind_id"), int lineageMode = FROM_INFO)
		: BaseMutator(rates, loci, mapIn, mapOut, context, output, begin, end,
		              step, at, reps, subPops, infoFields, lineageMode),
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
	virtual Allele mutate(Allele allele, size_t locus) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
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
		int context = 0, const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "ind_id"), int lineageMode = FROM_INFO)
		: BaseMutator(rates, loci, mapIn, mapOut, context, output, begin, end,
		              step, at, reps, subPops, infoFields, lineageMode),
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
	virtual Allele mutate(Allele allele, size_t locus) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
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
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "ind_id"), int lineageMode = FROM_INFO)
		: BaseMutator(rates, loci, mapIn, mapOut, 0, output, begin, end,
		              step, at, reps, subPops, infoFields, lineageMode),
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
	virtual Allele mutate(Allele allele, size_t locus) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
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
		const uintList & inds = vectoru(), const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = 0,
		const stringList & infoFields = vectorstr(1, "ind_id"), int lineageMode = FROM_INFO)
		: BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_lineageMode(lineageMode), m_loci(loci), m_allele(allele),
		m_ploidy(ploidy.elems()), m_inds(inds.elems())
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
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.PointMutator>" ;
	}


private:
	/// applicable loci.
	int m_lineageMode;

	lociList m_loci;
	Allele m_allele;
	vectoru m_ploidy;
	vectoru m_inds;
};


/** This operator checks all or specified loci of a population and revert all
 *  mutants at loci to wildtype alleles if they are fixed in the population.
 *  If a list of (virtual) subpopulations are specified, alleles are reverted
 *  if they are fixed in each subpopulation, regardless if the alleles are
 *  fixed in other subpopulations.
 */
class RevertFixedSites : public BaseOperator
{
public:
	/** Create an operator to set all alleles to zero at specified (parameter
	 *  \e loci) or all loci if they are fixed (having one non-zero allele) at these
	 *  loci. If parameter \e subPops are specified, only individuals in these
	 *  subpopulations are considered.
	 */
	RevertFixedSites(const lociList & loci = lociList(),
		const stringFunc & output = "", int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "ind_id"))
		: BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci)
	{
	}


	/// destructor
	virtual ~RevertFixedSites()
	{
	}


	/// HIDDEN Deep copy of a Migrator
	virtual BaseOperator * clone() const
	{
		return new RevertFixedSites(*this);
	}


	/// HIDDEN apply the Migrator to populaiton \e pop.
	virtual bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "Revert fixed alleles to wildtype allele 0 if they are fixed in the population.";
	}


private:
	lociList m_loci;

};


/** This is an infite site mutation model in mutational space. The alleles
 *  in the population is assumed to be locations of mutants. A mutation
 *  rate is given that mutate alleles in 'regions'. If number of mutants
 *  for an individual exceed the number of loci, 10 loci will be added
 *  to everyone in the population.
 */
class FiniteSitesMutator : public BaseOperator
{
public:
	/** This operator accepts a list of ranges which is the 'real range' of
	 *  each chromosome. Mutation happens with muation rate \e rate and mutants
	 *  will be recorded to the population (instead of alleles). By default,
	 *  this mutator assumes an finite-allele model where all mutations are
	 *  allowed and if a mutant (allele 1) is mutated, it will be mutated to
	 *  allele 0 (back mutation). Alternatively (\e model = 2), an
	 *  infinite-sites mutation model can be used where mutations can happen
	 *  only at a new locus. Mutations happen at a locus with existing mutants
	 *  will be moved to a random locus without existing mutant. A warning
	 *  message will be printed if there is no vacant locus available. If a
	 *  valid \e output is given, mutants will be outputted in the format of
	 *  "gen mutant ind type" where type is 0 for forward (0->1), 1 for
	 *  backward (1->0), 2 for relocated mutations, and 3 for ignored mutation
	 *  because no vacent locus is available. The second mode  has the
	 *  advantage that all mutants in the simulated population can be traced
	 *  to a single mutation event. If the regions are reasonably wide and
	 *  mutation rates are low, these two mutation models should yield
	 *  similar results.
	 */
	FiniteSitesMutator(double rate,
		// FIXME: we should not have ranges, because ranges are just
		// chromosomes, so, ranges=[1, 63000], [5000, 50000] shouldbe
		// loci=[63000, 450000], lociPos=range(63000) + range(5000, 500000)
		//
		// The problem is that you need to optimize the storage of loci positions
		// by saving ranges, not every locations in genoStru.h/cpp.
		//
		const intMatrix & ranges,
		int model = 1,
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "ind_id"), int lineageMode = FROM_INFO) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_rate(rate), m_ranges(ranges), m_model(model)
	{
		(void)lineageMode;  // avoid a compiler warning
		const matrixi & rngs = m_ranges.elems();

		for (size_t i = 0; i < rngs.size(); ++i) {
			DBG_FAILIF(rngs[i].size() != 2, ValueError, "Ranges should have two elements");
			for (size_t j = i + 1; j < rngs.size(); ++j) {
				DBG_FAILIF(rngs[j].size() != 2, ValueError, "Ranges should have two elements");
				if (i == j)
					continue;
				if (rngs[i][0] >= rngs[j][0] && rngs[i][0] <= rngs[j][1])
					throw ValueError("Overlapping ranges are currently not supported because of potential conflict of mutant locations on different chromosomes.");
				if (rngs[i][1] >= rngs[j][0] && rngs[i][1] <= rngs[j][1])
					throw ValueError("Overlapping ranges are currently not supported because of potential conflict of mutant locations on different chromosomes.");
			}
		}
#ifdef BINARYALLELE
		DBG_FAILIF(true, ValueError, "This operator does not work in binary allele type.");
#endif
	}


	/// destructor.
	~FiniteSitesMutator()
	{
	}


	virtual bool apply(Population & pop) const;

	/// HIDDEN Deep copy of a \c FiniteSitesMutator
	virtual BaseOperator * clone() const
	{
		return new FiniteSitesMutator(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.FiniteSitesMutator>";
	}


private:
	size_t locateVacantLocus(Population & pop, size_t beg, size_t end, std::set<size_t> & mutants) const;

private:
	const double m_rate;

	const intMatrix m_ranges;

	const int m_model;
};


#ifdef LONGALLELE


/** This operator looks into a population in mutational space and revert a mutant
 *  to wildtype allele if it is fixed in the population. If a valid output is
 *  specifieid, fixed alleles will be outputed with a leading generation number.
 */
class MutSpaceRevertFixedSites : public BaseOperator
{
public:
	/** Create an operator to revert alleles at fixed loci from value 1 to 0.
	 *  Parameter \e subPops is ignored.
	 */
	MutSpaceRevertFixedSites(const stringFunc & output = "", int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr())
		: BaseOperator(output, begin, end, step, at, reps, subPops, infoFields)
	{
	}


	/// destructor
	virtual ~MutSpaceRevertFixedSites()
	{
	}


	/// HIDDEN Deep copy of a Migrator
	virtual BaseOperator * clone() const
	{
		return new MutSpaceRevertFixedSites(*this);
	}


	/// HIDDEN apply the Migrator to populaiton \e pop.
	virtual bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "Revert fixed alleles to wildtype allele if it is fixed in the population.";
	}
};



/** This is an infite site mutation model in mutational space. The alleles
 *  in the population is assumed to be locations of mutants. A mutation
 *  rate is given that mutate alleles in 'regions'. If number of mutants
 *  for an individual exceed the number of loci, 10 loci will be added
 *  to everyone in the population.
 */
class MutSpaceMutator : public BaseOperator
{
public:
	/** This operator accepts a list of ranges which is the 'real range' of
	 *  each chromosome. Mutation happens with muation rate \e rate and mutants
	 *  will be recorded to the population (instead of alleles). By default,
	 *  this mutator assumes an finite-allele model where all mutations are
	 *  allowed and if a mutant (allele 1) is mutated, it will be mutated to
	 *  allele 0 (back mutation). Alternatively (\e model = 2), an
	 *  infinite-sites mutation model can be used where mutations can happen
	 *  only at a new locus. Mutations happen at a locus with existing mutants
	 *  will be moved to a random locus without existing mutant. A warning
	 *  message will be printed if there is no vacant locus available. If a
	 *  valid \e output is given, mutants will be outputted in the format of
	 *  "gen mutant ind type" where type is 0 for forward (0->1), 1 for
	 *  backward (1->0), 2 for relocated mutations, and 3 for ignored mutation
	 *  because no vacent locus is available. The second mode  has the
	 *  advantage that all mutants in the simulated population can be traced
	 *  to a single mutation event. If the regions are reasonably wide and
	 *  mutation rates are low, these two mutation models should yield
	 *  similar results.
	 */
	MutSpaceMutator(double rate, const intMatrix & ranges, int model = 1,
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_rate(rate), m_ranges(ranges), m_model(model)
	{
		const matrixi & rngs = m_ranges.elems();

		for (size_t i = 0; i < rngs.size(); ++i) {
			DBG_FAILIF(rngs[i].size() != 2, ValueError, "Ranges should have two elements");
			for (size_t j = i + 1; j < rngs.size(); ++j) {
				DBG_FAILIF(rngs[j].size() != 2, ValueError, "Ranges should have two elements");
				if (i == j)
					continue;
				if (rngs[i][0] >= rngs[j][0] && rngs[i][0] <= rngs[j][1])
					throw ValueError("Overlapping ranges are currently not supported because of potential conflict of mutant locations on different chromosomes.");
				if (rngs[i][1] >= rngs[j][0] && rngs[i][1] <= rngs[j][1])
					throw ValueError("Overlapping ranges are currently not supported because of potential conflict of mutant locations on different chromosomes.");
			}
		}
#  ifdef BINARYALLELE
		DBG_FAILIF(true, ValueError, "This operator does not work in binary allele type.");
#  endif
	}


	/// destructor.
	~MutSpaceMutator()
	{
	}


	virtual bool apply(Population & pop) const;

	/// HIDDEN Deep copy of a \c MutSpaceMutator
	virtual BaseOperator * clone() const
	{
		return new MutSpaceMutator(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.MutSpaceMutator>";
	}


private:
	size_t locateVacantLocus(Population & pop, size_t beg, size_t end, std::set<size_t> & mutants) const;

private:
	const double m_rate;

	const intMatrix m_ranges;

	const int m_model;
};


#endif

}
#endif
