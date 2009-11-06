/**
 *  $File: penetrance.h $
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

#ifndef _PENETRANCE_H
#define _PENETRANCE_H
/**
   \file
   \brief head file of class penetrance operator:public baseOperator
 */
#include "utility.h"
#include "operator.h"

#include "boost/tuple/tuple.hpp"
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
 *  will set the affection status of each offspring.
 *
 *  By default, a penetrance operator assigns affection status of individuals
 *  but does not save the actual penetrance value. However, if an information
 *  field is specified, penetrance values will be saved to this field for
 *  future analysis.
 *
 *  When a penetrance operator is applied to a population, it is only applied
 *  to the current generation. You can, however, use parameter \e ancGen=-1 to
 *  set affection status for all ancestral generations, or a generation index
 *  to apply to only ancestral generation younger than \e ancGen. Note that
 *  this parameter is ignored if the operator is applied during mating.
 */
class basePenetrance : public baseOperator
{
public:
	/** Create a base penetrance operator. If \e ancGen=0 (default), only the
	 *  current generation will be applied. If \e ancGen=-1, affection status
	 *  of all ancestral generations will be set. If a positive number is
	 *  given, ancestral generations with index <= ancGen will be applied. A
	 *  penetrance operator can be applied to specified (virtual)
	 *  subpopulations (parameter \e subPops) and replicates (parameter
	 *  \e reps). If an informatio field is given, penetrance value will be
	 *  stored in this information field of each individual.
	 */
	basePenetrance(int ancGen = 0,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr())
		: baseOperator("", begin, end, step, at, reps, subPops, infoFields),
		m_ancGen(ancGen)
	{
	}


	/// destructor
	virtual ~basePenetrance()
	{
	}


	/// deep copy of a penetrance operator
	virtual baseOperator * clone() const
	{
		return new basePenetrance(*this);
	}


	/** CPPONLY
	 *  calculate/return penetrance etc.
	 */
	virtual double penet(individual *, ULONG gen)
	{
		throw ValueError("This penetrance calculator is not supposed to be called directly");
		return 1.;
	}


	/// set penetrance to all individuals and record penetrance if requested
	virtual bool apply(population & pop);

	/// set penetrance to all individuals
	/// CPPONLY
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// HIDDEN
	string describe()
	{
		return "<simuPOP.penetrance>" ;
	}


private:
	/// how to handle ancestral gen
	int m_ancGen;
};

/** This penetrance operator assigns individual affection status using a
 *  user-specified penetrance dictionary.
 *  <applicability>all ploidy</applicability>
 *  <funcForm>MapPenetrance</funcForm>
 */
class mapPenetrance : public basePenetrance
{
public:
	/** Create a penetrance operator that get penetrance value from a
	 *  dictionary \e penetrance with genotype at \e loci as keys, and \e
	 *  penetrance as values. For each individual, genotypes at \e loci are
	 *  collected one by one (e.g. p0_loc0, p1_loc0, p0_loc1, p1_loc1... for
	 *  a diploid individual) and are looked up in the dictionary. If a
	 *  genotype cannot be found, it will be looked up again without phase
	 *  information (e.g. <tt>(1,0)</tt> will match key <tt>(0,1)</tt>). If the
	 *  genotype still can not be found, a \c ValueError will be raised. This
	 *  operator supports sex chromosomes and haplodiploid populations. In
	 *  these cases, only valid genotypes should be used to generator the
	 *  dictionary keys.
	 */
	mapPenetrance(const uintList & loci, const tupleDict & penetrance,
		int ancGen = 0, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		basePenetrance(ancGen, begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci.elems()), m_dict(penetrance)
	{
	};

	virtual ~mapPenetrance()
	{
	}


	/// deep copy of a map penetrance operator
	virtual baseOperator * clone() const
	{
		return new mapPenetrance(*this);
	}


	/// currently assuming diploid
	/// CPPONLY
	virtual double penet(individual * ind, ULONG gen);

	/// HIDDEN
	string describe()
	{
		return "<simuPOP.penetrance::map penetrance>" ;
	}


private:
	/// one locus
	vectoru m_loci;

	/// penetrance for each genotype
	tupleDict m_dict;
};

/** This operator is called a 'multi-allele' penetrance operator because it
 *  groups multiple alleles into two groups: wildtype and non-wildtype alleles.
 *  Alleles in each allele group are assumed to have the same effect on
 *  individual penetrance. If we denote all wildtype alleles as \c A, and all
 *  non-wildtype alleles \c a, this operator assign individual penetrance
 *  according to genotype \c AA, \c Aa, \c aa in the diploid case, and
 *  \c A and \c a in the haploid case.
 *  <funcForm>MaPenetrance</funcForm>
 */
class maPenetrance : public basePenetrance
{
public:
	/** Creates a multi-allele penetrance operator that groups multiple alleles
	 *  into a wildtype group (with alleles \e wildtype, default to <tt>[0]</tt>),
	 *  and a non-wildtype group. A list of penetrance values is specified
	 *  through parameter \e penetrance, for genotypes at one or more \e loci.
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
	maPenetrance(const uintList & loci, const vectorf & penetrance, const uintList & wildtype = vectoru(1, 0),
		int ancGen = 0, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		basePenetrance(ancGen, begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci.elems()), m_penetrance(penetrance), m_wildtype(wildtype.elems())
	{
		DBG_ASSERT(m_penetrance.size() == static_cast<UINT>(pow(static_cast<double>(3),
																static_cast<double>(m_loci.size()))),
			ValueError, "Please specify penetrance for each combination of genotype.");
	};

	virtual ~maPenetrance()
	{
	}


	/// deep copy of a multi-allele penetrance operator
	virtual baseOperator * clone() const
	{
		return new maPenetrance(*this);
	}


	/** CPPONLY
	 * currently assuming diploid
	 */
	virtual double penet(individual * ind, ULONG gen);

	/// HIDDEN
	string describe();

private:
	/// one locus
	vectoru m_loci;

	/// penetrance for each genotype
	vectorf m_penetrance;

	///
	vectoru m_wildtype;
};

/** This penetrance operator is created by a list of penetrance operators. When
 *  it is applied to an individual, it applies these penetrance operators to
 *  the individual, obtain a list of penetrance values, and compute a combined
 *  penetrance value from them and assign affection status accordingly.
 *  Additive, multiplicative, and a heterogeneour multi-locus model are
 *  supported. Please refer to Neil Rish (1989) "Linkage Strategies for
 *  Genetically Complex Traits" for some analysis of these models.
 *
 *  <funcForm>MlPenetrance</funcForm>
 */
class mlPenetrance : public basePenetrance
{
public:
	/** Create a multiple-locus penetrance operator from a list penetrance
	 *  operator \e ops. When this operator is applied to an individual (parents
	 *  when used before mating and offspring when used during mating), it
	 *  applies these operators to the individual and obtain a list of (usually
	 *  single-locus) penetrance values. These penetrance values are combined to a
	 *  single penetrance value using
	 *  \li <em>Prod(f_i)</em>, namely the product of individual penetrance if
	 *       \e mode = \c Multiplicative,
	 *  \li <em>sum(f_i)</em> if \e mode = \c Additive, and
	 *  \li <em>1-Prod(1 - f_i)</em> if \e mode = \c Heterogeneity
	 *
	 *  0 or 1 will be returned if the combined penetrance value is less than
	 *  zero or greater than 1.
	 */
	mlPenetrance(const opList & ops, int mode = Multiplicative,
		int ancGen = 0, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		basePenetrance(ancGen, begin, end, step, at, reps, subPops, infoFields),
		m_peneOps(ops), m_mode(mode)
	{
		DBG_FAILIF(ops.empty(), ValueError, "Please specify at least one penetrance operator.");
	};

	virtual ~mlPenetrance()
	{
	}


	/// deep copy of a multi-loci penetrance operator
	virtual baseOperator * clone() const
	{
		throw ValueError("Multi-loci penetrance operator can not be nested.");
	}


	/** CPPONLY
	 *  currently assuming diploid
	 */
	virtual double penet(individual * ind, ULONG gen);

	/// HIDDEN
	string describe()
	{
		return "<simuPOP.penetrance::multiple-loci penetrance>" ;
	}


private:
	/// a list of peneOps
	opList m_peneOps;

	/// mode
	int m_mode;
};

/** This penetrance operator assigns penetrance values by calling a user
 *  provided function. It accepts a list of loci and a Python function \c func.
 *  For each individual, this operator passes the genotypes at these loci, and
 *  a generation number to this function. The return value is treated as the
 *  penetrance value. Optionally, several information fields can be given to
 *  parameter \e paramFields. In this case, the user-defined Python function
 *  should accept a second parameter that is a list of values at these
 *  information fields. In another word, a user-defined function in the form
 *  of
 *  \li <tt>func(geno, gen)</tt> is needed if \c paramFields is empty, or
 *  \li <tt>func(geno, fields, gen)</tt> is needed if \c paramFields has some
 *      information fields.
 *
 *  If you need to pass sex or affection status to this function, you should
 *  define an information field (e.g. sex) and sync individual property with
 *  this field using operator \e infoExec (e.g.
 *  <tt>infoExec('sex=ind.sex', exposeInd='ind')</tt>.
 *
 *  <funcForm>PyPenetrance</funcForm>
 */
class pyPenetrance : public basePenetrance
{
public:
	/** Create a Python hybrid penetrance operator that passes genotype at specified
	 *  \e loci, optional values at specified information fields (parameter
	 *  \e paramFields), and a generation number to a user-defined function
	 *  \e func. The return value will be treated as individual penetrance.
	 */
	pyPenetrance(const uintList & loci, PyObject * func, int ancGen = 0,
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & paramFields = vectorstr(), const stringList & infoFields = vectorstr()) :
		basePenetrance(ancGen, begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci.elems()), m_func(func),  m_paramFields(paramFields.elems()),
		m_genotype(NULL), m_info(NULL)
	{
		if (!m_func.isValid())
			throw ValueError("Passed variable is not a callable python function.");

		DBG_FAILIF(m_loci.empty(), ValueError,
			"Please specify susceptibility loci");
	};

	/// destructor
	virtual ~pyPenetrance()
	{
		Py_XDECREF(m_genotype);
		Py_XDECREF(m_info);
	}


	/// CPPONLY
	pyPenetrance(const pyPenetrance & rhs) :
		basePenetrance(rhs),
		m_loci(rhs.m_loci),
		m_func(rhs.m_func),
		m_paramFields(rhs.m_paramFields),
		m_genotype(NULL),
		m_info(NULL)
	{
	}


	/// deep copy of a Python penetrance operator
	virtual baseOperator * clone() const
	{
		return new pyPenetrance(*this);
	}


	/** CPPONLY
	 *  currently assuming diploid
	 */
	virtual double penet(individual * ind, ULONG gen);

	/// HIDDEN
	string describe()
	{
		return "<simuPOP.penetrance::python penetrance>" ;
	}


private:
	/// susceptibility loci
	vectoru m_loci;

	/// user supplied python function
	pyFunc m_func;

	/// copy of information fields
	vectorstr m_paramFields;

	/// the object that passed to func
	PyObject * m_genotype;

	// the object that passed to func
	PyObject * m_info;
};

}
#endif
