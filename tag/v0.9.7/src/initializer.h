/**
 *  $File: initializer.h $
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

#ifndef _INITIALIZER_H
#define _INITIALIZER_H
/**
   \file
   \brief head file of class initializer:public baseOperator
 */
#include "utility.h"
#include "operator.h"

#include <numeric>
using std::string;
using std::accumulate;

namespace simuPOP {

/** This operator initialize sex of individuals, either randomly or use a list
 *  of sexes. For convenience, the function of this operator is included in
 *  other \e initializers such as \c initByFreq and \c initByValue so that you
 *  do not have to intiailize sexes separately from genotype.
 *  <funcForm>InitSex</funcForm>
 */
class initSex : public baseOperator
{
public:
	/** Create an operator that initialize individual sex to \c Male or \c Female.
	 *  By default, it assign sex to individuals randomly, with equal probability
	 *  of having a male or a female. This probabability can be adjusted through
	 *  parameter \e maleFreq. Alternatively, a fixed sequence of sexes can be
	 *  assigned. For example, if <tt>sex=[Male, Female]</tt>, individuals will
	 *  be assigned \c Male and \c Female successively. Parameter \e maleFreq
	 *  is ignored if \e sex is given. If a list of (virtual) subpopulation is
	 *  specified in parameter \e subPop, only individuals in these
	 *  subpopulations will be initialized. Note that the \e sex sequence, if
	 *  used, is assigned repeatedly regardless of subpopulation boundaries.
	 */
	initSex(double maleFreq = 0.5, const intList & sex = vectori(),
		int stage = PreMating, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr())
		: baseOperator("", stage, begin, end, step, at, reps, subPops, infoFields),
		m_maleFreq(maleFreq), m_sex(sex.elems())
	{
		if (!m_sex.empty()) {
			for (vectori::iterator it = m_sex.begin(); it != m_sex.end(); ++it) {
				DBG_ASSERT(*it == int(Male) || *it == int(Female),
					ValueError, "Parameter sex must be an array of Male or Female. ");
			}
		}
	}


	/// destructor
	virtual ~initSex()
	{
	}


	/// deep copy of an \c initSex operator.
	virtual baseOperator * clone() const
	{
		return new initSex(*this);
	}


	/// HIDDEN
	string description()
	{
		return "<simuPOP.initSex>";
	}


	/// apply this operator to population \e pop
	bool apply(population & pop);

protected:
	/// sex frequency
	double m_maleFreq;

	/// specify sex
	vectori m_sex;
};


/** This operator initializes given information fields with a sequence of
 *  values, or a user-provided function such as \c random.random.
 *  <funcForm>InitInfo</funcForm>
 */
class initInfo : public baseOperator
{
public:
	/** Create an operator that initialize individual information fields
	 *  \e infoFields using a sequence of values or a user-defined function.
	 *  If a list of values are given, it will be used sequentially for all
	 *  individuals. The values will be reused if its length is less than the
	 *  number of individuals. The values will be assigned repeatedly
	 *  regardless of subpopulation boundaries. If a Python function is given,
	 *  it will be called, without any argument, whenever a value is needed. If
	 *  a list of (virtual) subpopulation is specified in parameter \e subPop,
	 *  only individuals in these subpopulations will be initialized.
	 */
	initInfo(const floatListFunc & values,
		int stage = PreMating, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr())
		: baseOperator("", stage, begin, end, step, at, reps, subPops, infoFields),
		m_values(values)
	{
		DBG_FAILIF(m_values.empty() && !m_values.func().isValid(), ValueError, "Please specify a list of values or a Python function.");
	}


	/// destructor
	virtual ~initInfo()
	{
	}


	/// deep copy of an \c initInfo operator.
	virtual baseOperator * clone() const
	{
		return new initInfo(*this);
	}


	/// HIDDEN
	string description()
	{
		return "<simuPOP.initInfo>";
	}


	/// apply this operator to population \e pop
	bool apply(population & pop);

protected:
	/// specify sex
	floatListFunc m_values;
};


/** This operator assigns alleles at all or part of loci with given allele
 *  frequencies. Alternatively, an individual can be initialized and be copied
 *  to all individuals in the same (virtual) subpopulations.
 *  <funcForm>InitByFreq</funcForm>
 */
class initByFreq : public baseOperator
{
public:
	/** This function creates an initializer that initializes individual
	 *  genotypes randomly, using allele frequencies specified in parameter
	 *  \e alleleFreq. Elements in \e alleleFreq specifies the allele
	 *  frequencies of allele \c 0, \c 1, ... respectively. These frequencies
	 *  should add up to \c 1. If \e loci, \e ploidy and/or \e subPop are
	 *  specified, only specified loci, ploidy, and individuals in these
	 *  (virtual) subpopulations will be initialized. If \e identicalInds is
	 *  \c True, the first individual in each (virtual) subpopulation will be
	 *  initialized randomly, and be copied to all other individuals in this
	 *  (virtual) subpopulation. If a list of frequencies are given, they will
	 *  be used for each (virtual) subpopulation. This operator initializes all
	 *  chromosomes, including unused genotype locations and customized
	 *  chromosomes.
	 */
	initByFreq(const matrix & alleleFreq = matrix(), const uintList & loci = uintList(),
		const uintList & ploidy = uintList(), bool identicalInds = false,
		int stage = PreMating, int begin = 0, int end = 1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr());


	~initByFreq()
	{
	}


	/// deep copy of the operator \c initByFreq
	virtual baseOperator * clone() const
	{
		return new initByFreq(*this);
	}


	/// HIDDEN
	string description()
	{
		return "<simuPOP.initByFreq>";
	}


	/// apply this operator to population \e pop
	bool apply(population & pop);

private:
	/// allele frequencies (assume all loci are the same for a subPop
	matrix m_alleleFreq;

	///
	bool m_identicalInds;

	//
	uintList m_loci;

	//
	uintList m_ploidy;
};

/** This operator initialize individuals by given values.
 *  <funcForm>InitByValue</funcForm>
 */
class initByValue : public baseOperator
{
public:
	/** This function creates an initializer that initializes individual
	 *  genotypes with given genotype \e value. If \e loci, \e ploidy
	 *  and/or \e subPop are specified, only specified loci, ploidy, and
	 *  individuals in these (virtual) subpopulations will be initialized.
	 *  \e value can be used to initialize given \e loci, all loci, and all
	 *  homologous copies of these loci. If \e proportions (a list of positive
	 *  numbers that add up to \c 1) is given, \e value should be a list of
	 *  values that will be assigned randomly according to their respective
	 *  proportion. If a list of values are given without \e proportions,
	 *  they will be used for each (virtual) subpopulations. This operator
	 *  initializes all chromosomes, including unused genotype locations and
	 *  customized chromosomes.
	 */
	initByValue(intMatrix value = intMatrix(),
		const uintList & loci = uintList(), const uintList & ploidy = uintList(),
		const floatList & proportions = vectorf(),
		int stage = PreMating, int begin = 0, int end = 1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr());

	~initByValue()
	{
	}


	/// deep copy of the operator \c initByValue
	virtual baseOperator * clone() const
	{
		return new initByValue(*this);
	}


	/// HIDDEN
	string description()
	{
		return "<simuPOP.initByValue>";
	}


	/// apply this operator to population \e pop
	bool apply(population & pop);

private:
	/// allele frequencies (assume all loci are the same
	intMatrix m_value;

	/// if assign randomly
	vectorf m_proportion;

	//
	uintList m_loci;

	//
	uintList m_ploidy;
};

}
#endif
