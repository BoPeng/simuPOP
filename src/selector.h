/**
 *  $File: selector.h $
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

#ifndef _SELECTOR_H
#define _SELECTOR_H
/**
   \file
   \brief head file of class selector
 */
#include "utility.h"
#include "operator.h"

#include "boost_pch.hpp"
#include <numeric>
using std::min;

#if TR1_SUPPORT == 0
#  include <map>
#elif TR1_SUPPORT == 1
#  include <unordered_map>
#else
#  include <tr1/unordered_map>
#endif

#include <set>

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
 *       individual fitness is determined by the last fitness value it is
 *       assigned.
 *  \li A selection operator can be applied to virtual subpopulations and set
 *       fitness values only to part of the individuals.
 *  \li individuals with zero fitness in a subpopulation with anyone having a
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
 *       \e postOps parameter of \c Simulator.evolve, these fitness values will
 *       be used when the offspring generation becomes parental generation in
 *       the next generation.
 *
 *  Alternatively, a selector can be used as a during mating operator. In this
 *  case, it caculates fitness value for each offspring which will be treated
 *  as \b absolute fitness, namely the probability for each offspring to
 *  survive. This process uses the fact that an individual will be discarded
 *  when any of the during mating operators returns \e False. It is important
 *  to remember that:
 *  \li individual fitness needs to be between 0 and 1 in this case.
 *  \li Fitness values are not stored so the population does not need an
 *       information field \e fitness.
 *  \li This method applies natural selection to offspring instead of parents.
 *       These two implementation can be identical or different depending on
 *       the mating scheme used.
 *  \li Seleting offspring is less efficient than the selecting parents,
 *       especially when fitness values are low.
 *  \li Parameter \e subPops are applied to the offspring population and is
 *       used to judge if an operator should be applied. It thus does not make
 *       sense to apply a selector to a virtual subpopulation with affected
 *       individuals.
 */
class BaseSelector : public BaseOperator
{
public:
	/** Create a base selector object. This operator should not be created
	 *  directly.
	 */
	BaseSelector(const stringFunc & output = "", int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness"))
		: BaseOperator(output, begin, end, step, at, reps, subPops, infoFields)
	{
	}


	/// destructor
	virtual ~BaseSelector()
	{
	}


	/// HIDDEN Deep copy of a selector
	virtual BaseOperator * clone() const
	{
		return new BaseSelector(*this);
	}


	/// CPPONLY
	virtual double indFitness(Population & /* pop */, RawIndIterator) const
	{
		///
		throw ValueError("This selector is not supposed to be called directly");
		return 1.;
	}


	/// HIDDEN set fitness to all individuals. No selection will happen!
	bool apply(Population & pop) const;

	/// CPPONLY
	bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
	                       Individual * dad = NULL, Individual * mom = NULL) const
	{
		(void)pop;  // avoid warning about unused parameter
		(void)dad;  // avoid warning about unused parameter
		(void)mom;  // avoid warning about unused parameter
		// if offspring does not belong to subPops, do nothing, but does not fail.
		if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
			return true;

		double fitness = indFitness(offPop, offspring);

		DBG_FAILIF(fcmp_lt(fitness, 0) || fcmp_gt(fitness, 1), ValueError,
			"Fitness (probability for an offspring to survive) must be between 0 and 1 if a selector is used as a during-mating operator.");
		// accept an individual according to its fitness.
		return getRNG().randUniform() < fitness;
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.selector>" ;
	}


};


/** This selector assigns individual fitness values using a user-specified
 *  dictionary. This operator can be applied to populations with arbitrary
 *  number of homologous chromosomes.
 */
class MapSelector : public BaseSelector
{
public:
	/** Create a selector that assigns individual fitness values using a
	 *  dictionary \e fitness with genotype at \e loci as keys, and fitness
	 *  as values. Parameter \e loci can be a list of indexes, loci names, list
	 *  of chromosome position pairs, \c ALL_AVAIL, or a function with optional
	 *  parameter \c pop that will be called at each ganeeration to determine
	 *  indexes of loci. For each individual (parents if this operator is applied
	 *  before mating, and offspring if this operator is applied during
	 *  mating), genotypes at \e loci are collected one by one (e.g.
	 *  p0_loc0, p1_loc0, p0_loc1, p1_loc1... for a diploid individual, with
	 *  number of alleles varying for sex and mitochondrial DNAs) and
	 *  are looked up in the dictionary. If a genotype cannot be found, it
	 *  will be looked up again without phase information (e.g.
	 *  <tt>(1,0)</tt> will match key <tt>(0,1)</tt>). If the genotype
	 *  still can not be found, a \c ValueError will be raised. This
	 *  operator supports sex chromosomes and haplodiploid populations. In
	 *  these cases, only valid genotypes should be used to generator the
	 *  dictionary keys.
	 */
	MapSelector(const lociList & loci, const tupleDict & fitness,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		BaseSelector("", begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci), m_dict(fitness)
	{
	};

	virtual ~MapSelector()
	{
	}


	/// HIDDEN Deep copy of a map selector
	virtual BaseOperator * clone() const
	{
		return new MapSelector(*this);
	}


	/** CPPONLY
	 *  calculate/return the fitness value, currently assuming diploid
	 */
	virtual double indFitness(Population & pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.MapSelector>" ;
	}


	/// CPPONLY
	bool parallelizable() const
	{
		return true;
	}


private:
	///
	const lociList m_loci;

	/// fitness for each genotype
	const tupleDict m_dict;
};

/** This operator is called a 'multi-allele' selector because it groups
 *  multiple alleles into two groups: wildtype and non-wildtype alleles.
 *  Alleles in each allele group are assumed to have the same effect on
 *  individual fitness. If we denote all wildtype alleles as \c A, and all
 *  non-wildtype alleles \c a, this operator assign individual fitness
 *  according to genotype \c AA, \c Aa, \c aa in the diploid case, and
 *  \c A and \c a in the haploid case.
 */
class MaSelector : public BaseSelector
{
public:
	/** Creates a multi-allele selector that groups multiple alleles into a
	 *  wildtype group (with alleles \e wildtype, default to <tt>[0]</tt>), and
	 *  a non-wildtype group. A list of fitness values is specified through
	 *  parameter \e fitness, for genotypes at one or more \e loci. Parameter
	 *  \e loci can be a list of indexes, loci names , list of chromosome position
	 *  pairs, \c ALL_AVAIL, or a function with optional parameter \c pop that
	 *  will be called at each ganeeration to determine indexes of loci. If we
	 *  denote wildtype alleles using capital letters \c A, \c B ... and
	 *  non-wildtype alleles using small letters \c a, \c b ..., the fitness
	 *  values should be for
	 *  \li genotypes \c A and \c a for the haploid single-locus case,
	 *  \li genotypes \c AB, \c Ab, \c aB and \c bb for haploid two=locus cases,
	 *  \li genotypes \c AA, \c Aa and \c aa for diploid single-locus cases,
	 *  \li genotypes \c AABB, \c AABb, \c AAbb, \c AaBB, \c AaBb, \c Aabb,
	 *       \c aaBB, \c aaBb, and \c aabb for diploid two-locus cases,
	 *  \li and in general 2**n for diploid and 3**n for haploid cases if there
	 *       are \c n loci.
	 *
	 *  This operator does not support haplodiploid populations, sex and
	 *  mitochondrial chromosomes.
	 */
	MaSelector(const lociList & loci, const vectorf & fitness, const uintList & wildtype = vectoru(1, 0),
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		BaseSelector("", begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci), m_fitness(fitness), m_wildtype(wildtype.elems())
	{
		DBG_WARNIF(m_wildtype.empty(), "No wild type allele is defined.");
	};

	virtual ~MaSelector()
	{
	}


	/// HIDDEN Deep copy of a \c MaSelector
	virtual BaseOperator * clone() const
	{
		return new MaSelector(*this);
	}


	/// CPPONLY
	/// calculate/return the fitness value, currently assuming diploid
	virtual double indFitness(Population & pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.MaSelector>" ;
	}


	/// CPPONLY
	bool parallelizable() const
	{
		return true;
	}


private:
	/// one locus
	const lociList m_loci;

	/// fitness for each genotype
	const vectorf m_fitness;

	///
	const vectoru m_wildtype;
};


/** This selector is created by a list of selectors. When it is applied to an
 *  individual, it applies these selectors to the individual, obtain a list of
 *  fitness values, and compute a combined fitness value from them. ADDITIVE,
 *  multiplicative, and a heterogeneour multi-locus model are supported.
 */
class MlSelector : public BaseSelector
{
public:
	/** Create a multiple-locus selector from a list selection operator
	 *  \e selectors. When this operator is applied to an individual (parents
	 *  when used before mating and offspring when used during mating), it
	 *  applies these operators to the individual and obtain a list of (usually
	 *  single-locus) fitness values. These fitness values are combined to a
	 *  single fitness value using
	 *  \li <em>Prod(f_i)</em>, namely the product of individual fitness if
	 *       \e mode = \c MULTIPLICATIVE,
	 *  \li <em>1-sum(1 - f_i)</em> if \e mode = \c ADDITIVE,
	 *  \li <em>1-Prod(1 - f_i)</em> if \e mode = \c HETEROGENEITY, and
	 *  \li <em>exp(- sum(1 - f_i))</em> if \e mode = \c EXPONENTIAL,
	 *
	 *  zero will be returned if the combined fitness value is less than zero.
	 *
	 *  Applicability parameters (begin, end, step, at, reps, subPops) could be
	 *  used in both \c MlSelector and selectors in parameter \e ops, but
	 *  parameters in \c MlSelector will be interpreted first.
	 */
	MlSelector(const opList & ops, int mode = MULTIPLICATIVE,
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		BaseSelector("", begin, end, step, at, reps, subPops, infoFields),
		m_selectors(ops), m_mode(mode)
	{
		DBG_FAILIF(ops.empty(), ValueError, "Please specify at least one selector.");
	};

	virtual ~MlSelector()
	{
	}


	/// HIDDEN Deep copy of a \c MlSelector
	virtual BaseOperator * clone() const
	{
		return new MlSelector(*this);
	}


	/** CPPONLY
	 *  calculate/return the fitness value, currently assuming diploid
	 */
	virtual double indFitness(Population & pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.MlSelector>" ;
	}


	/// CPPONLY
	bool parallelizable() const
	{
		return true;
	}


private:
	/// a list of selectors
	const opList m_selectors;

	/// mode
	const int m_mode;
};


/** This selector assigns fitness values by calling a user provided function.
 *  It accepts a list of loci (parameter \e loci) and a Python function \c func
 *  which should be defined with one or more of parameters \c geno, \c mut,
 *  \c gen, \c ind, \c pop or names of information fields. Parameter \e loci can
 *  be a list of loci indexes, names, list of chromosome position pairs,
 *  \c ALL_AVAIL, or a function with optional parameter \c pop that will be
 *  called at each ganeeration to determine indexes of loci. When this operator is applied
 *  to a population, it passes genotypes or mutants at specified loci, generation
 *  number, a reference to an individual, a reference to the current population
 *  (usually used to retrieve population variable), and values at specified
 *  information fields to respective parameters of this function. Genotypes are
 *  passed as a tuple of alleles arranged locus by locus (in the order of
 *  A1,A2,B1,B2 for loci A and B). Mutants are passed as a default dictionary of loci
 *  index (with respect to all genotype of individuals, not just the first ploidy)
 *  and alleles. The returned value will be used to determine the fitness of each
 *  individual.
 */
class PySelector : public BaseSelector
{
public:
	/** Create a Python hybrid selector that passes genotype at specified
	 *  \e loci, values at specified information fields (if requested) and
	 *  a generation number to a user-defined function \e func. The return
	 *  value will be treated as individual fitness.
	 */
	PySelector(PyObject * func, lociList loci = vectoru(),
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const stringFunc & output = "",
		const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		BaseSelector(output, begin, end, step, at, reps, subPops, infoFields),
		m_func(func), m_loci(loci)
	{
		DBG_ASSERT(m_func.isValid(), ValueError, "Passed variable is not a callable python function.");
	}


	/// HIDDEN Deep copy of a \c PySelector
	virtual BaseOperator * clone() const
	{
		return new PySelector(*this);
	}


	/** CPPONLY
	 *  calculate/return the fitness value, currently assuming diploid
	 */
	virtual double indFitness(Population & pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.PySelector>" ;
	}


	/// CPPONLY
	bool parallelizable() const
	{
		return false;
	}


private:
	/// user supplied python function
	const pyFunc m_func;

	/// susceptibility loci
	const lociList m_loci;

};


/** This selector is a multi-locus Python selector that assigns fitness to
 *  individuals by combining locus and genotype specific fitness values.
 *  It differs from a \c PySelector in that the python function is responsible
 *  for assigning fitness values for each gentoype type at each locus, which
 *  can potentially be random, and locus or gentoype-specific.
 */
class PyMlSelector : public BaseSelector
{
public:
	/** Create a selector that assigns individual fitness values by combining
	 *  locus-specific fitness values that are determined by a Python call-back
	 *  function. The callback function accepts parameter \e loc, \e alleles
	 *  (both optional) and returns location- or genotype-specific fitness values
	 *  that can be constant or random. The fitness values for each genotype will
	 *  be cached so the same fitness values will be assigned to genotypes with
	 *  previously assigned values. Note that a function that does not examine the
	 *  genotype naturally assumes a dominant model where genotypes with one or
	 *  two mutants have the same fitness effect. Because genotypes at a locus
	 *  are passed separately and in no particular order, this function is also
	 *  responsible for assigning consistent fitness values for genotypes at the
	 *  same locus (a class is usually used). This operator currently ignores
	 *  chromosome types so unused alleles will be passed for loci on sex or
	 *  mitochondrial chromosomes. It also ignores phase of genotype so it will
	 *  use the same fitness value for genotype (a,b) and (b,a).
	 *
	 *  Individual fitness will be combined in \c ADDITIVE, \c MULTIPLICATIVE,
	 *  \c HETEROGENEITY, or \c EXPONENTIAL mode from fitness values of loci with
	 *  at least one non-zero allele (See \c MlSelector for details). If an output
	 *  is given, location, genotype, fitness and generation at which the new
	 *  genotype is assgined the value will be written to the output, in the
	 *  format of 'loc a1 a2 fitness gen' for loci on autosomes of diploid
	 *  populations.
	 */
	PyMlSelector(PyObject * func, int mode = EXPONENTIAL,
		const lociList & loci = lociList(), const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness"));

	virtual ~PyMlSelector()
	{
	}


	/// HIDDEN Deep copy of a map selector
	virtual BaseOperator * clone() const
	{
		return new PyMlSelector(*this);
	}


	/** CPPONLY
	 *  calculate/return the fitness value, currently assuming diploid
	 */
	virtual double indFitness(Population & pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.PyMlSelector>" ;
	}


	/// CPPONLY
	bool apply(Population & pop) const;

	typedef std::pair<size_t, vectora> LocGenotype;

private:
	double getFitnessValue(const LocGenotype &) const;

	///
	pyFunc m_func;

	int m_mode;

	lociList m_loci;

	int m_searchMode;

	/// tr1 map cannot be used because LocMutant etc are not hashable
	typedef std::map<LocGenotype, double> GenoSelMap;

	mutable std::vector<LocGenotype> m_newGenotypes;
	mutable GenoSelMap m_fitnessFactory;
};



#ifdef LONGALLELE

/** This selector assumes that alleles are mutant locations in the mutational
 *  space and assign fitness values to them according to a random distribution.
 *  The overall individual fitness is determined by either an additive, an
 *  multiplicative or an exponential model.
 */
class MutSpaceSelector : public BaseSelector
{
public:
	/** Create a selector that assigns individual fitness values according to
	 *  random fitness effects. \e selDist can be
	 *  \li <tt>(CONSTANT, s, h)</tt> where s will be used for all mutants. The
	 *      fitness value for genotypes AA, Aa and aa will be (1, 1-hs, 1-s).
	 *      If h is unspecified, a default value h=0.5 (additive model) will
	 *      be used.
	 *  \li <tt>(GAMMA_DISTRIBUTION, theta, k, h</tt> where s follows a gamma
	 *      distribution with scale parameter theta and shape parameter k.
	 *      Fitness values for genotypes AA, Aa and aa will be 1, 1-hs and 1-s.
	 *      A default value h=0.5 will be used if h is unspecified.
	 *  \li a Python function, which will be called when selection coefficient
	 *      of a new mutant is needed. This function should return a single
	 *      value s (with default value h=0.5) or a sequence of (h, s). Mutant
	 *      location will be passed to this function if it accepts a parameter
	 *      \c loc. This allows the definition of site-specific selection
	 *      coefficients.
	 *  Individual fitness will be combined in \c ADDITIVE,
	 *     \c MULTIPLICATIVE or \c EXPONENTIAL mode. (See \c MlSelector for
	 *     details).
	 *  If an output is given, mutants and their fitness values will be written
	 *  to the output, in the form of 'mutant s h'.
	 */
	MutSpaceSelector(const floatListFunc & selDist, int mode = EXPONENTIAL,
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		BaseSelector(output, begin, end, step, at, reps, subPops, infoFields),
		m_selDist(selDist), m_mode(mode), m_selFactory(), m_additive(true)
	{
		if (m_selDist.size() == 0) {
			DBG_FAILIF(!m_selDist.func().isValid(), ValueError,
				"Please specify either a distribution with parameter or a function.");
		} else if (static_cast<int>(m_selDist[0]) == CONSTANT) {
			DBG_FAILIF(m_selDist.size() < 2, ValueError, "At least one parameter is needed for constant selection coefficient.");
		} else if (static_cast<int>(m_selDist[0]) == GAMMA_DISTRIBUTION) {
			DBG_FAILIF(m_selDist.size() < 3, ValueError, "At least two parameters are needed for gamma distribution.");
		}
	}


	virtual ~MutSpaceSelector()
	{
	}


	/// HIDDEN Deep copy of a map selector
	virtual BaseOperator * clone() const
	{
		return new MutSpaceSelector(*this);
	}


	/** CPPONLY
	 *  calculate/return the fitness value, currently assuming diploid
	 */
	virtual double indFitness(Population & pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.MutSpaceSelector>" ;
	}


	/// CPPONLY
	bool apply(Population & pop) const;

	typedef std::pair<double, double> SelCoef;

private:
	SelCoef getFitnessValue(size_t mutant) const;


	double randomSelAddFitness(GenoIterator it, GenoIterator it_end, bool maleChrX) const;

	double randomSelExpFitness(GenoIterator it, GenoIterator it_end, bool maleChrX) const;

	// extended models does not assume additivity (h != 0.5)
	double randomSelMulFitnessExt(GenoIterator it, GenoIterator it_end, bool maleChrX) const;

	double randomSelAddFitnessExt(GenoIterator it, GenoIterator it_end, bool maleChrX) const;

	double randomSelExpFitnessExt(GenoIterator it, GenoIterator it_end, bool maleChrX) const;

private:
	///
	floatListFunc m_selDist;

	int m_mode;
	///
#  if TR1_SUPPORT == 0
	typedef std::map<unsigned int, SelCoef> SelMap;
	typedef std::map<unsigned int, int> MutCounter;
#  elif TR1_SUPPORT == 1
	// this is faster than std::map
	typedef std::unordered_map<size_t, SelCoef> SelMap;
	typedef std::unordered_map<size_t, size_t> MutCounter;
#  else
	// this is faster than std::map
	typedef std::tr1::unordered_map<size_t, SelCoef> SelMap;
	typedef std::tr1::unordered_map<size_t, size_t> MutCounter;
#  endif
	mutable SelMap m_selFactory;
	mutable vectoru m_newMutants;
	// whether or not all markers are additive.
	mutable bool m_additive;
};



#endif

}
#endif
