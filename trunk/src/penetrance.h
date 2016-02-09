/**
 *  $File: penetrance.h $
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

#ifndef _PENETRANCE_H
#define _PENETRANCE_H
/**
   \file
   \brief head file of class penetrance operator:public BaseOperator
 */
#include "utility.h"
#include "operator.h"

#include "boost_pch.hpp"
#include <numeric>
using std::min;

namespace simuPOP {

/** A penetrance model models the probability that an individual has a certain
 *  disease provided that he or she has certain genetic (genotype) and
 *  environmental (information field) riske factors. A penetrance operator
 *  calculates this probability according to provided information and set his
 *  or her affection status randomly. For example, an individual will have
 *  probability 0.8 to be affected if the penetrance is 0.8. This class is the
 *  base class to all penetrance operators and defines a common interface for
 *  all penetrance operators.
 *
 *  A penetrance operator can be applied at any stage of an evolutionary cycle.
 *  If it is applied before or after mating, it will set affection status of
 *  all parents and offspring, respectively. If it is applied during mating, it
 *  will set the affection status of each offspring. You can also apply a
 *  penetrance operator to an individual using its \c applyToIndividual
 *  member function.
 *
 *  By default, a penetrance operator assigns affection status of individuals
 *  but does not save the actual penetrance value. However, if an information
 *  field is specified, penetrance values will be saved to this field for
 *  future analysis.
 *
 *  When a penetrance operator is applied to a population, it is only applied
 *  to the current generation. You can, however, use parameter \e ancGens to
 *  set affection status for all ancestral generations (\c ALL_AVAIL), or
 *  individuals in specified generations if a list of ancestral generations
 *  is specified. Note that this parameter is ignored if the operator is
 *  applied during mating.
 */
class BasePenetrance : public BaseOperator
{
public:
	/** Create a base penetrance operator. This operator assign individual
	 *  affection status in the present generation (default). If \c ALL_AVAIL
	 *  or a list of ancestral generations are spcified in parameter \e ancGens,
	 *  individuals in specified ancestral generations will be processed. A
	 *  penetrance operator can be applied to specified (virtual)
	 *  subpopulations (parameter \e subPops) and replicates (parameter
	 *  \e reps). If an informatio field is given, penetrance value will be
	 *  stored in this information field of each individual.
	 */
	BasePenetrance(const uintList & ancGens = uintList(NULL),
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr())
		: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
		m_ancGens(ancGens)
	{
	}


	/// destructor
	virtual ~BasePenetrance()
	{
	}


	/// HIDDEN Deep copy of a penetrance operator
	virtual BaseOperator * clone() const
	{
		return new BasePenetrance(*this);
	}


	/** CPPONLY
	 *  calculate/return penetrance etc.
	 */
	virtual double penet(Population *, RawIndIterator) const
	{
		throw ValueError("This penetrance calculator is not supposed to be called directly");
		return 1.;
	}


	/// set penetrance to all individuals and record penetrance if requested
	virtual bool apply(Population & pop) const;

	/** Apply the penetrance operator to a single individual \e ind and set his
	 *  or her affection status. A population reference can be passed if the
	 *  penetrance model depends on population properties such as generation
	 *  number. This function returns the affection status.
	 */
	virtual bool applyToIndividual(Individual * ind, Population * pop = NULL);

	/// set penetrance to all individuals
	/// CPPONLY
	virtual bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.penetrance>" ;
	}


private:
	/// how to handle ancestral gen
	const uintList m_ancGens;
};

/** This penetrance operator assigns individual affection status using a
 *  user-specified penetrance dictionary.
 */
class MapPenetrance : public BasePenetrance
{
public:
	/** Create a penetrance operator that get penetrance value from a
	 *  dictionary \e penetrance with genotype at \e loci as keys, and \e
	 *  penetrance as values. For each individual, genotypes at \e loci are
	 *  collected one by one (e.g. p0_loc0, p1_loc0, p0_loc1, p1_loc1... for
	 *  a diploid individual) and are looked up in the dictionary. Parameter
	 *  \e loci can be a list of loci indexes, names, list of chromosome
	 *  position pairs, \c ALL_AVAIL, or a function with optional parameter
	 *  \c pop that will be called at each ganeeration to determine indexes of loci.
	 *  If a genotype cannot be found, it will be looked up again without phase
	 *  information (e.g. <tt>(1,0)</tt> will match key <tt>(0,1)</tt>). If the
	 *  genotype still can not be found, a \c ValueError will be raised. This
	 *  operator supports sex chromosomes and haplodiploid populations. In
	 *  these cases, only valid genotypes should be used to generator the
	 *  dictionary keys.
	 */
	MapPenetrance(const lociList & loci, const tupleDict & penetrance,
		const uintList & ancGens = uintList(NULL), int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		BasePenetrance(ancGens, begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci), m_dict(penetrance)
	{
	};

	virtual ~MapPenetrance()
	{
	}


	/// HIDDEN Deep copy of a map penetrance operator
	virtual BaseOperator * clone() const
	{
		return new MapPenetrance(*this);
	}


	/// currently assuming diploid
	/// CPPONLY
	virtual double penet(Population * pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.MapPenetrance> map penetrance";
	}


	/// CPPONLY
	bool parallelizable() const
	{
		return true;
	}


private:
	/// one locus
	const lociList m_loci;

	/// penetrance for each genotype
	const tupleDict m_dict;
};

/** This operator is called a 'multi-allele' penetrance operator because it
 *  groups multiple alleles into two groups: wildtype and non-wildtype alleles.
 *  Alleles in each allele group are assumed to have the same effect on
 *  individual penetrance. If we denote all wildtype alleles as \c A, and all
 *  non-wildtype alleles \c a, this operator assign Individual penetrance
 *  according to genotype \c AA, \c Aa, \c aa in the diploid case, and
 *  \c A and \c a in the haploid case.
 */
class MaPenetrance : public BasePenetrance
{
public:
	/** Creates a multi-allele penetrance operator that groups multiple alleles
	 *  into a wildtype group (with alleles \e wildtype, default to <tt>[0]</tt>),
	 *  and a non-wildtype group. A list of penetrance values is specified
	 *  through parameter \e penetrance, for genotypes at one or more \e loci.
	 *  Parameter \e loci can be a list of loci indexes, names, list of chromosome
	 *  position pairs, \c ALL_AVAIL, or a function with optional parameter \c pop
	 *  that will be called at each ganeeration to determine indexes of loci.
	 *  If we denote wildtype alleles using capital letters \c A, \c B ... and
	 *  non-wildtype alleles using small letters \c a, \c b ..., the penetrance
	 *  values should be for
	 *  \li genotypes \c A and \c a for the haploid single-locus case,
	 *  \li genotypes \c AB, \c Ab, \c aB and \c bb for haploid two=locus cases,
	 *  \li genotypes \c AA, \c Aa and \c aa for diploid single-locus cases,
	 *  \li genotypes \c AABB, \c AABb, \c AAbb, \c AaBB, \c AaBb, \c Aabb,
	 *       \c aaBB, \c aaBb, and \c aabb for diploid two-locus cases,
	 *  \li and in general 2**n for diploid and 3**n for haploid cases if there
	 *       are \c n loci.
	 *
	 *  This operator does not support haplodiploid populations and sex
	 *  chromosomes.
	 */
	MaPenetrance(const lociList & loci, const vectorf & penetrance, const uintList & wildtype = vectoru(1, 0),
		const uintList & ancGens = uintList(NULL), int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		BasePenetrance(ancGens, begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci), m_penetrance(penetrance), m_wildtype(wildtype.elems())
	{
		DBG_ASSERT(m_penetrance.size() == static_cast<UINT>(pow(static_cast<double>(3),
																static_cast<double>(m_loci.size()))),
			ValueError, "Please specify penetrance for each combination of genotype.");
	};

	virtual ~MaPenetrance()
	{
	}


	/// HIDDEN Deep copy of a multi-allele penetrance operator
	virtual BaseOperator * clone() const
	{
		return new MaPenetrance(*this);
	}


	/** CPPONLY
	 * currently assuming diploid
	 */
	virtual double penet(Population * pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const;

	/// CPPONLY
	bool parallelizable() const
	{
		return true;
	}


private:
	/// one locus
	const lociList m_loci;

	/// penetrance for each genotype
	const vectorf m_penetrance;

	///
	const vectoru m_wildtype;
};

/** This penetrance operator is created by a list of penetrance operators. When
 *  it is applied to an individual, it applies these penetrance operators to
 *  the individual, obtain a list of penetrance values, and compute a combined
 *  penetrance value from them and assign affection status accordingly.
 *  ADDITIVE, multiplicative, and a heterogeneour multi-locus model are
 *  supported. Please refer to Neil Rish (1989) "Linkage Strategies for
 *  Genetically Complex Traits" for some analysis of these models.
 *
 */
class MlPenetrance : public BasePenetrance
{
public:
	/** Create a multiple-locus penetrance operator from a list penetrance
	 *  operator \e ops. When this operator is applied to an individual (parents
	 *  when used before mating and offspring when used during mating), it
	 *  applies these operators to the individual and obtain a list of (usually
	 *  single-locus) penetrance values. These penetrance values are combined to a
	 *  single penetrance value using
	 *  \li <em>Prod(f_i)</em>, namely the product of individual penetrance if
	 *       \e mode = \c MULTIPLICATIVE,
	 *  \li <em>sum(f_i)</em> if \e mode = \c ADDITIVE, and
	 *  \li <em>1-Prod(1 - f_i)</em> if \e mode = \c HETEROGENEITY
	 *
	 *  0 or 1 will be returned if the combined penetrance value is less than
	 *  zero or greater than 1.
	 *
	 *  Applicability parameters (begin, end, step, at, reps, subPops) could be
	 *  used in both \c MlSelector and selectors in parameter \e ops, but
	 *  parameters in \c MlSelector will be interpreted first.
	 */
	MlPenetrance(const opList & ops, int mode = MULTIPLICATIVE,
		const uintList & ancGens = uintList(NULL), int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		BasePenetrance(ancGens, begin, end, step, at, reps, subPops, infoFields),
		m_peneOps(ops), m_mode(mode)
	{
		DBG_FAILIF(ops.empty(), ValueError, "Please specify at least one penetrance operator.");
	};

	virtual ~MlPenetrance()
	{
	}


	/// HIDDEN Deep copy of a multi-loci penetrance operator
	virtual BaseOperator * clone() const
	{
		return new MlPenetrance(*this);
	}


	/** CPPONLY
	 *  currently assuming diploid
	 */
	virtual double penet(Population * pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.MlPenetrance> multiple-loci penetrance>" ;
	}


	/// CPPONLY
	bool parallelizable() const
	{
		return true;
	}


private:
	/// a list of peneOps
	const opList m_peneOps;

	/// mode
	const int m_mode;
};

/** This penetrance operator assigns penetrance values by calling a user
 *  provided function. It accepts a list of loci (parameter \c loci),
 *  and a Python function \c func which should be defined with one or more of
 *  parameters \c geno, \c mut, \c gen, \c ind, \c pop, or names of
 *  information fields. When this operator is applied to a population, it
 *  passes genotypes or mutants (non-zero alleles) at specified loci at
 *  specified loci, generation number, a reference to an individual, a
 *  reference to the current population (usually used to retrieve population
 *  variables) and values at specified information fields to respective
 *  parameters of this function. Genotypes of each individual are
 *  passed as a tuple of alleles arranged locus by locus (in the order of
 *  A1,A2,B1,B2 for loci A and B). Mutants are passed as a default dictionary
 *  of loci index (with respect to all genotype of individuals, not just
 *  the first ploidy) and alleles. The returned penetrance values will be
 *  used to determine the affection status of each individual.
 */
class PyPenetrance : public BasePenetrance
{
public:
	/** Create a Python hybrid penetrance operator that passes genotype at
	 *  specified \e loci, values at specified information fields (if
	 *  requested), and a generation number to a user-defined function \e func.
	 *  Parameter \e loci can be a list of loci indexes, names, list
	 *  of chromosome position pairs, \c ALL_AVAIL, or a function with optional
	 *  parameter \c pop that will be called at each ganeeration to determine
	 *  indexes of loci. The return value will be treated as Individual penetrance.
	 */
	PyPenetrance(PyObject * func,
		const lociList & loci = vectoru(),
		const uintList & ancGens = uintList(NULL),
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		BasePenetrance(ancGens, begin, end, step, at, reps, subPops, infoFields),
		m_func(func), m_loci(loci)
	{
		DBG_ASSERT(m_func.isValid(), ValueError, "Passed variable is not a callable python function.");
	};


	/// HIDDEN Deep copy of a Python penetrance operator
	virtual BaseOperator * clone() const
	{
		return new PyPenetrance(*this);
	}


	/** CPPONLY
	 *  currently assuming diploid
	 */
	virtual double penet(Population * pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.PyPenetrance> python penetrance>" ;
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


/** This penetrance operator is a multi-locus Python penetrance operator that
 *  assigns penetrance values by combining locus and genotype specific penetrance
 *  values. It differs from a \c PyPenetrance in that the python function is
 *  responsible for penetrance values values for each gentoype type at each locus,
 *  which can potentially be random, and locus or gentoype-specific.
 */
class PyMlPenetrance : public BasePenetrance
{
public:
	/** Create a penetrance operator that assigns individual affection status
	 *  according to penetrance values combined from locus-specific penetrance
	 *  values that are determined by a Python call-back function. The callback
	 *  function accepts parameter \e loc, \e alleles (both optional) and returns
	 *  location- or genotype-specific penetrance values that can be constant or
	 *  random. The penetrance values for each genotype will be cached so the
	 *  same penetrance values will be assigned to genotypes with previously
	 *  assigned values. Note that a function that does not examine the genotype
	 *  naturally assumes a dominant model where genotypes with one or two mutants
	 *  have the same penetrance value. Because genotypes at a locus are passed
	 *  separately and in no particular order, this function is also responsible
	 *  for assigning consistent fitness values for genotypes at the same locus
	 *  (a class is usually used). This operator currently ignores chromosome
	 *  types so unused alleles will be passed for loci on sex or mitochondrial
	 *  chromosomes. This operator also ignores the phase of genotype so genotypes
	 *  (a,b) and (b,a) are assumed to have the same fitness effect.
	 *
	 *  Individual penetrance will be combined in \c ADDITIVE, \c MULTIPLICATIVE,
	 *  or \c HETEROGENEITY mode from penetrance values of loci with
	 *  at least one non-zero allele (See \c MlPenetrance for details).
	 */
	PyMlPenetrance(PyObject * func, int mode = MULTIPLICATIVE,
		const lociList & loci = lociList(),
		const uintList & ancGens = uintList(NULL),
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr());

	virtual ~PyMlPenetrance()
	{
	}


	/// HIDDEN Deep copy of a map penetrance operator
	virtual BaseOperator * clone() const
	{
		return new PyMlPenetrance(*this);
	}


	/** CPPONLY
	 *  calculate/return the penetrance value, currently assuming diploid
	 */
	double penet(Population * pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.PyMlPenetrance>" ;
	}


	typedef std::pair<size_t, vectora> LocGenotype;

private:
	double getPenetranceValue(const LocGenotype &) const;

	///
	pyFunc m_func;

	int m_mode;

	lociList m_loci;

	int m_searchMode;

	/// tr1 map cannot be used because LocMutant etc are not hashable
	typedef std::map<LocGenotype, double> GenoPenetranceMap;

	mutable GenoPenetranceMap m_penetFactory;
};


}
#endif
