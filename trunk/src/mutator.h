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
/// Base class of all mutators
/**
   The base class of all functional mutators. It is not supposed to be called directly.
   \n
   Every mutator can specify \c rate (equal rate or different rates for different
   loci) and a vector of applicable loci (default to all but should have the same
   length as \c rate if \c rate has length greater than one).
   \n
   Maximum allele can be specified as well but more parameters, if needed, should
   be implemented by individual mutator classes.
   \n
   There are numbers of possible allelic states. Most theoretical studies assume an infinite
   number of allelic states to avoid any homoplasy. If it facilitates any analysis,
   this is however extremely unrealistic.
 */
class mutator : public baseOperator
{
public:
	/// create a mutator, do not call this constructor directly
	/**
	   All mutators have the following common parameters. However, the actual meaning
	   of these parameters may vary according to different models. The only differences
	   between the following mutators are the way they actually mutate an allele, and
	   corresponding input parameters. The number of mutation events at each locus is
	            recorded and can be accessed from the \c mutationCount or \c mutationCounts
	            functions.

	   \param rate can be a number (uniform rate) or an array of mutation rates (the same length as \c loci)
	   \param loci a vector of locus indexes. Will be ignored only when single rate is specified.
	    Default to all loci.
	   \param maxAllele maximum allowed allele. Interpreted by each sub mutator class. Default to \c pop.maxAllele().
	 */
	mutator(const vectorf & rate = vectorf(),
		const vectoru & loci = vectoru(),
		UINT maxAllele = 0,
		const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr())
		: baseOperator(output, stage, begin, end, step, at, rep, subPops, infoFields),
		m_rate(rate), m_maxAllele(maxAllele), m_loci(loci),
		m_bt(rng()), m_initialized(false), m_mutCount(0)
	{
		if (m_rate.empty() )
			throw ValueError("You should specify a rate, or a sequence of rate.");

		if (rate.size() > 1 && loci.empty())
			throw ValueError("If you use variable rates, you should specify loci for each of the rate.");

		if (rate.size() > 1 && !loci.empty() && rate.size() != loci.size() )
			throw ValueError("If both rates and loci are specified, they should have the same length.");

#ifdef BINARYALLELE
		DBG_WARNING(maxAllele > 1, "MaxAllele for binary libraries must be 1");
		m_maxAllele = 1;
#else
		DBG_ASSERT(maxAllele <= ModuleMaxAllele, ValueError,
			"The maximum allele number exceeds " + toStr(ModuleMaxAllele)
			+ ". \nIf you need longer allele size, please use simuPOP_la libraries.");
#endif
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


	/// return the mutation rate
	vectorf rate()
	{
		return m_rate;
	}


	/// set an array of mutation rates
	void setRate(const vectorf & rate, const vectoru & loci = vectoru())
	{
		if (rate.size() != 1 && rate.size() != loci.size() )
			throw ValueError("If you specify more than one rate values, you should also specify corresponding applicable loci");

		m_rate = rate;
		if (!loci.empty())
			m_loci = loci;

		m_initialized = false;
	}


	/// return maximum allowable allele number
	UINT maxAllele()
	{
		return m_maxAllele;
	}


	/// set maximum allowable allele
	void setMaxAllele(UINT maxAllele)
	{
#ifndef BINARYALLELE
		m_maxAllele = maxAllele;
#endif
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


	/// describe how to mutate a single allele
	virtual void mutate(AlleleRef allele)
	{
		throw SystemError("You are not supposed to call this base mutator funciton.");
	};

	/// apply a mutator
	virtual bool apply(population & pop);

private:
	/// initialize bernulli trial according to pop size etc
	virtual void initialize(population & pop);

private:
	/// mutation rates
	vectorf m_rate;

	/// maxAllele
	UINT m_maxAllele;

	/// applicable loci.
	vectoru m_loci;

	/// bernulli trials. bitSet mutation results.
	BernulliTrials m_bt;

	/// initialized? the first apply() call will trigger an initialization process.
	bool m_initialized;

	/// report the number of mutation events
	vectoru m_mutCount;
};

/** 
 *
 */
class matrixMutator : public mutator
{
public:
	matrixMutator(const matrix & rate, const vectoru & loci = vectoru(),
		const stringFunc & output = ">", int stage = PostMating,
		int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(),
		const vectorstr & infoFields = vectorstr());

	~matrixMutator()
	{
	}


	/// mutate to a state other than current state with equal probability
	virtual void mutate(AlleleRef allele);

	/// deep copy of a \c matrixMutator
	virtual baseOperator * clone() const
	{
		return new matrixMutator(*this);
	}


	/// used by Python print function to print out the general information of the \c kamMutator
	virtual string __repr__()
	{
		return "<simuPOP::matrix mutator";
	}


private:
	vector<weightedSampler> m_sampler;
};

/**
   This mutator mutate an allele to another allelic state with equal probability.
   The specified mutation rate is actually the 'probability to mutate'. So the
   mutation rate to any other allelic state is actually \f$ \frac{rate}{K-1} \f$, where
   \f$ K \f$ is specified by parameter \c maxAllele.
   <funcForm>KamMutate</funcForm>
   \sa Crow & Kimura 1970
 */
class kamMutator : public mutator
{
public:
	/// create a K-Allele Model mutator
	/**
	   \param rate mutation rate. It is the 'probability to mutate'. The actual
	    mutation rate to any of the other \c K-1 allelic states are <tt>rate/(K-1)</tt>.
	   \param maxAllele maximum allele that can be mutated to. For binary libraries,
	   allelic states will be <tt>[0, maxAllele]</tt>. Otherwise, they are <tt>[1, maxAllele]</tt>.

	   Please see class \c mutator for the descriptions of other parameters.
	 */
	kamMutator(const vectorf & rate = vectorf(),
		const vectoru & loci = vectoru(),
		UINT maxAllele = 0,
		const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr())
		: mutator(rate, loci, maxAllele,
		          output, stage, begin, end, step, at, rep, subPops, infoFields)
	{
	}


	~kamMutator()
	{
	}


	/// mutate to a state other than current state with equal probability
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
		       toStr(this->maxAllele()) + ">" ;
	}


};

/// The stepwise mutation model
/**
   The <em>Stepwise Mutation Model</em> (SMM) assumes that alleles are represented by integer values
   and that a mutation either increases or decreases the allele value by one.
   For variable number tandem repeats(VNTR) loci, the allele value is generally
   taken as the number of tandem repeats in the DNA sequence.
   <funcForm>SmmMutate</funcForm>
   \sa Kimura & Ohta 1978
 */
class smmMutator : public mutator
{
public:
	/// create a SMM mutator
	/**
	   The SMM is developed for allozymes. It  provides better description
	   for these kinds of evolutionary processes.

	   \param incProb probability to increase allele state. Default to \c 0.5.

	   Please see class \c mutator for the descriptions of other parameters.

	 */
	smmMutator(const vectorf & rate = vectorf(),
		const vectoru & loci = vectoru(),
		UINT maxAllele = 0, double incProb = 0.5,
		const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr())
		: mutator(rate, loci, maxAllele,
		          output, stage, begin, end, step, at, rep, subPops, infoFields),
		m_incProb(incProb)
	{
#ifdef BINARYALLELE
		DBG_WARNING(true, "Symetric stepwise mutation does not work well on two state alleles.");
#endif
		DBG_ASSERT(fcmp_ge(incProb, 0.) && fcmp_le(incProb, 1.),
			ValueError, "Inc probability should be between [0,1], given " + toStr(incProb));
	}


	~smmMutator()
	{
	}


	/// mutate according to the SMM model
	/// CPPONLY
	virtual void mutate(AlleleRef allele)
	{
		// inc
		if (rng().randUniform01() < m_incProb) {
			if (AlleleUnsigned(allele) < this->maxAllele() )
				AlleleInc(allele);
		}
		// dec (use !=0 instead of > 0 to avoid warning inbinary mode
		else if (allele != 0)
			AlleleDec(allele);
	}


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
	/// probability to increase allele state
	double m_incProb;
};

/// generalized stepwise mutation model
/**
   The <em>Generalized Stepwise Mutation model</em> (GSM) is an extension to the stepwise
   mutation model. This model assumes that alleles are represented by integer values
   and that a mutation either increases or decreases the allele value by a random value.
   In other words, in this model the change in the allelic state is drawn from a random
   distribution. A <em>geometric generalized stepwise model</em> uses a geometric distribution with
   parameter \f$ p \f$, which has mean \f$ \frac{p}{1-p} \f$ and variance \f$ \frac{p}{\left(1-p\right)^{2}} \f$. \n

   \c gsmMutator implements both models. If you specify a Python function without a
   parameter, this mutator will use its return value each time a mutation occur; otherwise,
   a parameter \f$ p \f$ should be provided and the mutator will act as a geometric generalized stepwise model.

   <funcForm>GsmMutate</funcForm>
   \sa Kimura & Ohta 1978
 */
class gsmMutator : public mutator
{
public:
	/// create a \c gsmMutator
	/**
	   The GSM model is developed for allozymes.
	   It  provides better description for these kinds of evolutionary processes.

	   \param incProb probability to increase allele state. Default to \c 0.5.
	   \param func a function that returns the number of steps. This function
	    does not accept any parameter.

	   Please see class \c mutator for the descriptions of other parameters.
	 */
	gsmMutator(const vectorf & rate = vectorf(),
		const vectoru & loci = vectoru(),
		UINT maxAllele = 0, double incProb = 0.5, double p = 0, PyObject * func = NULL,
		const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr())
		: mutator(rate, loci, maxAllele,
		          output, stage, begin, end, step, at, rep, subPops, infoFields),
		m_incProb(incProb), m_p(p), m_func(func)
	{
		DBG_ASSERT(fcmp_ge(incProb, 0.) && fcmp_le(incProb, 1.),
			ValueError, "Inc probability should be between [0,1], given " + toStr(incProb));

#ifdef BINARYALLELE
		DBG_WARNING(true, "Generalized stepwise mutation does not work well on two state alleles.");
#endif

		if (!m_func.isValid()) {
			DBG_ASSERT(fcmp_ge(p, 0.) && fcmp_le(p, 1.),
				ValueError, "Parameter p of a geometric distribution should be between [0,1], given " + toStr(m_p));
		}
	}


	~gsmMutator()
	{
	}


	/// deep copy of a \c gsmMutator
	virtual baseOperator * clone() const
	{
		return new gsmMutator(*this);
	}


	/// mutate according to the GSM model
	virtual void mutate(AlleleRef allele);

	/// used by Python print function to print out the general information of the \c gsmMutator
	virtual string __repr__()
	{
		return "<simuPOP::generalized step-wise mutator>" ;
	}


private:
	/// probability to increase allele state
	double m_incProb;

	/// parameter for geometric gsm
	double m_p;

	/// the function to return random number
	pyFunc m_func;
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
	pyMutator(const vectorf & rate = vectorf(),
		const vectoru & loci = vectoru(), UINT maxAllele = 0,
		PyObject * func = NULL,
		const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr())
		: mutator(rate, loci, maxAllele,
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
	pointMutator(
		const vectoru & loci,
		Allele toAllele,
		vectoru atPloidy = vectoru(),
		vectorlu inds = vectorlu(),
		const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr())
		: baseOperator(output, stage, begin, end, step, at, rep, subPops, infoFields),
		m_loci(loci), m_toAllele(toAllele),
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
	vectoru m_loci;
	Allele m_toAllele;
	vectoru m_atPloidy;
	vectorlu m_inds;
	/// report the number of mutation events
	vectoru m_mutCount;
};

}
#endif
