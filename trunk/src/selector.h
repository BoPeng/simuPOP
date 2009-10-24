/**
 *  $File: selector.h $
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

#ifndef _SELECTOR_H
#define _SELECTOR_H
/**
   \file
   \brief head file of class selector
 */
#include "utility.h"
#include "operator.h"

#include "boost/tuple/tuple.hpp"
#include <numeric>
using std::min;

namespace simuPOP {

/** This class is the base class to all selectors, namely operators that
 *  perform natural selection. It defines a common interface for all selectors.
 *
 *  A selector can be applied before mating or during mating. If a selector is
 *  applied to one or more (virtual) subpopulations of a parental population
 *  before mating, it sets individual fitness values to all involved parents to
 *  an information field (default to \e fitness). When a mating scheme that
 *  supports natural selection is applied to the parental population, it will
 *  select parents with probabilities that are proportional to individual
 *  fitness stored in an information field (default to \e fitness). Individual
 *  fitness is considered \b relative fitness and can be any non-negative
 *  number. This simple process has some implications that can lead to advanced
 *  usages of natural selection in simuPOP:
 *  \li It is up to the mating scheme how to handle individual fitness. Some
 *      mating schemes do not support natural selection at all.
 *  \li A mating scheme performs natural selection according to fitness values
 *       stored in an information field. It does not care how these values are
 *       set. For example, fitness values can be inherited from a parent using
 *       a tagging operator, or set directly using a Python operator.
 *  \li A mating scheme can treat any information field as fitness field. If an
 *       specified information field does not exist, or if all individuals have
 *       the same fitness values (e.g. 0), the mating scheme selects parents
 *       randomly.
 *  \li Multiple selectors can be applied to the same parental generation.
 *       Individual fitness is determined by the last fitness value it is
 *       assigned.
 *  \li A selection operator can be applied to virtual subpopulations and set
 *       fitness values only to part of the individuals.
 *  \li Individuals with zero fitness in a subpopulation with anyone having a
 *       positive fitness value will not be selected to produce offspring. This
 *       can sometimes lead to unexpected behaviors. For example, if you only
 *       assign fitness value to part of the individuals in a subpopulation,
 *       the rest of them will be effectively discarded. If you migrate
 *       individuals with valid fitness values to a subpopulation with all
 *       individuals having zero fitness, the migrants will be the only mating
 *       parents.
 *  \li It is possible to assign multiple fitness values to different
 *       information fields so that different homogeneous mating schemes can
 *       react to different fitness schemes when they are used in a
 *       heterogeneous mating scheme.
 *  \li You can apply a selector to the offspring generation using the
 *       \e postOps parameter of \c simulator.evolve, these fitness values will
 *       be used when the offspring generation becomes parental generation in
 *       the next generation.
 *
 *  Alternatively, a selector can be used as a during mating operator. In this
 *  case, it caculates fitness value for each offspring which will be treated
 *  as \b absolute fitness, namely the probability for each offspring to
 *  survive. This process uses the fact that an individual will be discarded
 *  when any of the during mating operators returns \e False. It is important
 *  to remember that:
 *  \li Individual fitness needs to be between 0 and 1 in this case.
 *  \li This method applies natural selection to offspring instead of parents.
 *       These two implementation can be identical or different depending on
 *       the mating scheme used.
 *  \li Seleting offspring is less efficient than the selecting parents,
 *       especially when fitness values are low.
 *
 *  It is worth noting that a selector used as a during-mating operator does
 *  not support parameter \e subPops. If you need to apply different selection
 *  scheme to different virtual subpopulations, you can use different selectors
 *  in a heterogeneous mating scheme.
 */
class selector : public baseOperator
{
public:
	/** Create a base selector object. This operator should not be created
	 *  directly.
	 */
	selector(int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness"))
		: baseOperator("", begin, end, step, at, reps, subPops, infoFields)
	{
	}


	/// destructor
	virtual ~selector()
	{
	}


	/// deep copy of a selector
	virtual baseOperator * clone() const
	{
		return new selector(*this);
	}


	/// CPPONLY
	virtual double indFitness(individual *, ULONG gen)
	{
		///
		throw ValueError("This selector is not supposed to be called directly");
		return 1.;
	}


	/// set fitness to all individuals. No selection will happen!
	bool apply(population & pop);

	/// CPPONLY
	bool applyDuringMating(population & pop, RawIndIterator offspring,
	                       individual * dad = NULL, individual * mom = NULL)
	{
		double fitness = indFitness(& * offspring, pop.gen());

		DBG_FAILIF(fcmp_lt(fitness, 0) || fcmp_gt(fitness, 1), ValueError,
			"Fitness (probability for an offspring to survive) must be between 0 and 1 if a selector is used as a during-mating operator.");
		// accept an individual according to its fitness.
		return GetRNG().randUniform() < fitness;
	}


	/// HIDDEN
	string description()
	{
		return "<simuPOP.selector>" ;
	}


};


/** This selector assigns individual fitness values using a user-specified
 *  dictionary. This operator can be applied to populations with arbitrary
 *  number of homologous chromosomes.
 *  <applicability>all ploidy</applicability>
 */
class mapSelector : public selector
{
public:
	/** Create a selector that assigns individual fitness values using a
	 *  dictionary \e fitness with genotype at \e loci as keys, and fitness
	 *  as values. For each individual (parents if this operator is applied
	 *  before mating, and offspring if this operator is applied during
	 *  mating), genotypes at \e loci are collected one by one (e.g.
	 *  p0_loc0, p1_loc0, p0_loc1, p1_loc1... for a diploid individual) and
	 *  are looked up in the dictionary. If a genotype cannot be found, it
	 *  will be looked up again without phase information (e.g.
	 *  <tt>(1,0)</tt> will match key <tt>(0,1)</tt>). If the genotype
	 *  still can not be found, a \c ValueError will be raised. This
	 *  operator supports sex chromosomes and haplodiploid populations. In
	 *  these cases, only valid genotypes should be used to generator the
	 *  dictionary keys.
	 */
	mapSelector(const uintList & loci, const tupleDict & fitness,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		selector(begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci.elems()), m_dict(fitness)
	{
	};

	virtual ~mapSelector()
	{
	}


	/// deep copy of a map selector
	virtual baseOperator * clone() const
	{
		return new mapSelector(*this);
	}


	/** CPPONLY
	 *  calculate/return the fitness value, currently assuming diploid
	 */
	virtual double indFitness(individual * ind, ULONG gen);

	/// HIDDEN
	string description()
	{
		return "<simuPOP.selector::map selector>" ;
	}


private:
	///
	vectoru m_loci;

	/// fitness for each genotype
	tupleDict m_dict;
};

/** This operator is called a 'multi-allele' selector because it groups
 *  multiple alleles into two groups: wildtype and non-wildtype alleles.
 *  Alleles in each allele group are assumed to have the same effect on
 *  individual fitness. If we denote all wildtype alleles as \c A, and all
 *  non-wildtype alleles \c a, this operator assign individual fitness
 *  according to genotype \c AA, \c Aa, \c aa in the diploid case, and
 *  \c A and \c a in the haploid case.
 */
class maSelector : public selector
{
public:
	/** Creates a multi-allele selector that groups multiple alleles into a
	 *  wildtype group (with alleles \e wildtype, default to \c [0]), and a
	 *  non-wildtype group. A list of fitness values is specified through
	 *  parameter \e fitness, for genotypes at one or more \e loci. If we
	 *  denote wildtype alleles using capital letters \c A, \c B ... and
	 *  non-wildtype alleles using small letters \c a, \c b ..., the fitness
	 *  values should be for
	 *  \li genotypes A and a for the haploid single-locus case,
	 *  \li genotypes AB, Ab, aB and bb for haploid two=locus cases,
	 *  \li genotypes AA, Aa and aa for diploid single-locus cases,
	 *  \li genotypes AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, and aabb
	 *       for diploid two-locus cases,
	 *  \li and in general 2**n for diploid and 3**n for haploid cases if there
	 *       are \c n loci.
	 *
	 *  This operator does not support haplodiploid populations and sex
	 *  chromosomes.
	 */
	maSelector(const uintList & loci, const vectorf & fitness, const uintList & wildtype = vectoru(1, 0),
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		selector(begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci.elems()), m_fitness(fitness), m_wildtype(wildtype.elems())
	{
		DBG_WARNING(m_wildtype.empty(), "No wild type allele is defined.");
	};

	virtual ~maSelector()
	{
	}


	/// deep copy of a \c maSelector
	virtual baseOperator * clone() const
	{
		return new maSelector(*this);
	}


	/// CPPONLY
	/// calculate/return the fitness value, currently assuming diploid
	virtual double indFitness(individual * ind, ULONG gen);

	/// HIDDEN
	string description()
	{
		return "<simuPOP.selector::multiple-alleles selector>" ;
	}


private:
	/// one locus
	vectoru m_loci;

	/// fitness for each genotype
	vectorf m_fitness;

	///
	vectoru m_wildtype;
};


/// selection according to genotypes at multiple loci in a multiplicative model
/**
   This selector is a 'multiple-locus model' selector. The selector takes a vector of
   selectors (can not be another \c mlSelector) and evaluate the fitness of an
   individual as the product or sum of individual fitness values. The mode is
   determined by parameter \c mode, which takes one of the following values
   \li \c Multiplicative: the fitness is calculated as \f$ f=\prod_{i}f_{i} \f$, where \f$ f_{i} \f$
   is the single-locus fitness value.
   \li \c Additive: the fitness is calculated as
   \f$ f=\max\left(0,1-\sum_{i}(1-f_{i})\right) \f$.
   \f$ f \f$ will be set to \c 0 when \f$ f<0 \f$.

   <funcForm>MlSelect</funcForm>
 */
class mlSelector : public selector
{
public:
	/// create a multiple-locus selector
	/**
	   \param selectors a list of selectors

	   Please refer to \c mapSelector for other parameter descriptions.

	 */
	mlSelector(const opList & selectors, int mode = Multiplicative,
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		selector(begin, end, step, at, reps, subPops, infoFields),
		m_selectors(selectors), m_mode(mode)
	{
		DBG_FAILIF(selectors.empty(), ValueError, "Please specify at least one selector.");
	};

	virtual ~mlSelector()
	{
	}


	/// deep copy of a \c mlSelector
	virtual baseOperator * clone() const
	{
		return new mlSelector(*this);
	}


	/** CPPONLY
	 *  calculate/return the fitness value, currently assuming diploid
	 */
	virtual double indFitness(individual * ind, ULONG gen);

	/// HIDDEN
	string description()
	{
		return "<simuPOP.selector::multiple-loci selector>" ;
	}


private:
	/// a list of selectors
	opList m_selectors;

	/// mode
	int m_mode;
};

/// selection using user provided function
/**
   This selector assigns fitness values by calling a user provided function.
   It accepts a list of loci and a Python function \c func. For
   each individual, this operator will pass the genotypes at these loci,
   generation number, and optionally values at some information fields
   to this function. The return value is treated as the fitness value.
   The genotypes are arranged in the order
   of <tt>0-0,0-1,1-0,1-1</tt> etc. where X-Y represents locus X - ploidy Y.
   More specifically, \c func can be
   \li <tt>func(geno, gen)</tt> if \c infoFields has length 0 or 1.
   \li <tt>func(geno, gen, fields)</tt> when \c infoFields has more than 1 fields.
    Values of fields 1, 2, ... will be passed.
   Both \c geno and \c fields should be a list.

   <funcForm>PySelect</funcForm>
 */
class pySelector : public selector
{
public:
	/// create a Python hybrid selector
	/**
	   \param loci susceptibility loci. The genotype at these loci will be
	    passed to \c func.
	   \param func a Python function that accepts genotypes at specified loci,
	    generation number, and optionally information fields. It returns the fitness value.
	   \param output and other parameters please refer to help (<tt>baseOperator.__init__</tt>)
	   \param infoFields if specified, the first field should be the information
	    field to save calculated fitness value (should be 'fitness' in most cases).
	    The values of the rest of the information fields (if available) will also
	    be passed to the user defined penetrance function.
	 */
	// provide locus and fitness for 11, 12, 13 (in the form of dictionary)
	pySelector(uintList loci, PyObject * func,
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		selector(begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci.elems()), m_func(func), m_alleles(0), m_len(0), m_numArray(NULL)
	{
		if (!m_func.isValid())
			throw ValueError("Passed variable is not a callable python function.");

		DBG_FAILIF(m_loci.empty(), ValueError,
			"Please specify susceptibility loci");
	};

	/// CPPONLY
	pySelector(const pySelector & rhs) :
		selector(rhs),
		m_loci(rhs.m_loci),
		m_func(rhs.m_func),
		m_alleles(rhs.m_alleles),
		m_info(rhs.m_info),
		m_len(rhs.m_len),
		m_numArray(NULL),
		m_infoArray(NULL)
	{
	}


	/// deep copy of a \c pySelector
	virtual baseOperator * clone() const
	{
		return new pySelector(*this);
	}


	/** CPPONLY
	 *  calculate/return the fitness value, currently assuming diploid
	 */
	virtual double indFitness(individual * ind, ULONG gen);

	/// HIDDEN
	string description()
	{
		return "<simuPOP.selector::python selector>" ;
	}


private:
	/// susceptibility loci
	vectoru m_loci;

	/// user supplied python function
	pyFunc m_func;

	/// copy of alleles of each individual a time.
	vectora m_alleles;

	/// copy of information fields
	vectorinfo m_info;

	/// length of m_alleles
	int m_len;

	/// the object that passed to func
	PyObject * m_numArray;

	// the object that passed to func
	PyObject * m_infoArray;

};
}
#endif
