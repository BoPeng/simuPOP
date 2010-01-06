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
   \brief head file of class initializer:public BaseOperator
 */
#include "utility.h"
#include "operator.h"

#include <numeric>
using std::string;
using std::accumulate;

namespace simuPOP {

/** This operator initializes sex of individuals, either randomly or use a list
 *  of sexes.
 *  <funcForm>initSex</funcForm>
 */
class InitSex : public BaseOperator
{
public:
	/** Create an operator that initializes individual sex to \c MALE or
	 *  \c FEMALE. By default, it assign sex to individuals randomly, with
	 *  equal probability of having a male or a female. This probabability can
	 *  be adjusted through parameter \e maleFreq or be made to exact
	 *  proportions by specifying parameter \e maleProp. Alternatively, a fixed
	 *  sequence of sexes can be assigned. For example, if
	 *  <tt>sex=[MALE, FEMALE]</tt>, individuals will be assigned \c MALE and
	 *  \c FEMALE successively. Parameter \e maleFreq or \e maleProp are
	 *  ignored if \e sex is given. If a list of (virtual) subpopulation is
	 *  specified in parameter \e subPop, only individuals in these
	 *  subpopulations will be initialized. Note that the \e sex sequence, if
	 *  used, is assigned repeatedly regardless of (virtual) subpopulation
	 *  boundaries.
	 */
	InitSex(double maleFreq = 0.5, double maleProp = -1, const intList & sex = vectori(),
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr())
		: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
		m_maleFreq(maleFreq), m_maleProp(maleProp), m_sex(sex.elems())
	{
		if (!m_sex.empty()) {
			vectori::const_iterator it = m_sex.begin();
			for (; it != m_sex.end(); ++it) {
				DBG_ASSERT(*it == int(MALE) || *it == int(FEMALE),
					ValueError, "Parameter sex must be an array of MALE or FEMALE. ");
			}
		}
	}


	/// destructor
	virtual ~InitSex()
	{
	}


	/// HIDDEN Deep copy of an \c InitSex operator.
	virtual BaseOperator * clone() const
	{
		return new InitSex(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;


	/// HIDDEN apply this operator to population \e pop
	bool apply(Population & pop) const;

protected:
	/// sex frequency
	const double m_maleFreq;

	const double m_maleProp;

	/// specify sex
	const vectori m_sex;
};


/** This operator initializes given information fields with a sequence of
 *  values, or a user-provided function such as \c random.random.
 *  <funcForm>initInfo</funcForm>
 */
class InitInfo : public BaseOperator
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
	InitInfo(const floatListFunc & values,
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr())
		: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
		m_values(values)
	{
		DBG_FAILIF(m_values.empty() && !m_values.func().isValid(), ValueError, "Please specify a list of values or a Python function.");
	}


	/// destructor
	virtual ~InitInfo()
	{
	}


	/// HIDDEN Deep copy of an \c InitInfo operator.
	virtual BaseOperator * clone() const
	{
		return new InitInfo(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;

	/// HIDDEN apply this operator to population \e pop
	bool apply(Population & pop) const;

protected:
	/// specify sex
	const floatListFunc m_values;
};


/** This operator assigns alleles at all or part of loci with given allele
 *  frequencies, proportions or values.
 *  <funcForm>initGenotype</funcForm>
 */
class InitGenotype : public BaseOperator
{
public:
	/** This function creates an initializer that initializes individual
	 *  genotypes either randomly with specified allele frequencies
	 *  (parameter \e freq) or proportions (parameter \e prop), or with fixed
	 *  genotypes (\e genotype). Elements in \e freq specifies the allele
	 *  frequencies of allele \c 0, \c 1, ... respectively. These frequencies
	 *  should add up to \c 1. Elements in \e prop specified the proportions of
	 *  alleles. Parameter \e prop is similar to \e freq except that \e prop
	 *  guarantees exact proportions of alleles, although alleles with small
	 *  proportions might not be allocated at all. Parameter \e genotype
	 *  specifies a list of genotype that will be assigned repeatedly to all
	 *  individuals (similar to \c Population.setGenotype() except that this
	 *  operator supports parameter \e ploidy). If \e loci, \e ploidy and/or
	 *  \e subPop are specified, only specified loci, ploidy, and individuals
	 *  in these (virtual) subpopulations will be initialized. This operator
	 *  initializes all chromosomes, including unused genotype locations and
	 *  customized chromosomes.
	 */
	InitGenotype(const vectorf & freq = vectorf(),
		const uintList & genotype = uintList(), const vectorf & prop = vectorf(),
		const uintList & loci = uintList(),
		const uintList & ploidy = uintList(),
		int begin = 0, int end = 1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr());


	~InitGenotype()
	{
	}


	/// HIDDEN Deep copy of the operator \c InitGenotype
	virtual BaseOperator * clone() const
	{
		return new InitGenotype(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;


	/// HIDDEN apply this operator to population \e pop
	bool apply(Population & pop) const;

private:
	/// allele frequencies (assume all loci are the same for a subPop
	const vectorf m_freq;
	const vectoru m_genotype;
	const vectorf m_prop;
	//
	const uintList m_loci;

	//
	const uintList m_ploidy;
};


}
#endif
