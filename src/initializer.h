/**
 *  $File: initializer.h $
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
 */
class InitSex : public BaseOperator
{
public:
	/** Create an operator that initializes individual sex to \c MALE or
	 *  \c FEMALE. By default, it assigns sex to individuals randomly, with
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
	 *  boundaries so that you can assign \e sex to all individuals in a
	 *  population.
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
				PARAM_ASSERT(*it == int(MALE) || *it == int(FEMALE),
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
		PARAM_FAILIF(m_values.empty() && !m_values.func().isValid(), ValueError, "Please specify a list of values or a Python function.");
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
 *  frequencies, proportions or values. This operator initializes all
 *  chromosomes, including unused genotype locations and customized
 *  chromosomes.
 */
class InitGenotype : public BaseOperator
{
public:
	/** This function creates an initializer that initializes individual
	 *  genotypes with random alleles, genotypes, or haplotypes with specified
	 *  frequencies (parameter \e freq) or proportions (parameter \e prop).
	 *  If parameter \e genotypes or \e haplotypes is not specified, \e freq
	 *  specifies the allele frequencies of alleles \c 0, \c 1, \c 2...
	 *  respectively. Alternatively, you can use parameter \e prop to
	 *  specified the exact proportions of alleles \c 0, \c 1, ..., although
	 *  alleles with small proportions might not be assigned at all.

	 *  Values of parameter \e prob or \e prop should add up to 1. In addition
	 *  to a vector, parameter \e prob and \e prop can also be a function that
	 *  accepts optional parameters \e loc, \e subPop or \e vsp and returns a
	 *  list of requencies for alleles  \c 0, \c 1, etc, or a number for frequency
	 *  of allele \c 0 as a speciail case for each locus, subpopulation (parameter
	 *  \e subPop), or virtual subpopulations (parameter \e vsp, pass as a tuple).

	 *  If parameter \e genotypes is specified, it should contain a list of genotypes
	 *  (alleles on different strand of chromosomes) with length equal to population
	 *  ploidy. Parameter \e prob and \e prop then specifies frequencies or proportions
	 *  of each genotype, which can vary for each subpopulation but not each locus if
	 *  the function form of parameters is used.

	 *  If parameter \e haplotypes is specified, it should contain a list of 
	 *  haplotypes (alleles on the same strand of chromosome) and parameter \e prob
	 *  or \e prop specifies frequencies or proportions of each haplotype.

	 *  If \e loci, \e ploidy and/or \e subPop are specified, only specified loci,
	 *  ploidy, and individuals in these (virtual) subpopulations will be initialized.
	 *  Parameter \e loci can be a list of loci indexes, names or \c ALL_AVAIL.
	 *  If the length of a haplotype is not enough to fill all loci, the
	 *  haplotype will be reused. If a list (or a single) haplotypes are
	 *  specified without \e freq or \e prop, they are used with equal
	 *  probability.
	 *
	 *  In the last case, if a sequence of genotype is specified through parameter
	 *  \e genotype (not \e genotypes), it will be used repeatedly to initialize all
	 *  alleles sequentially. This works similar to function \c Population.setGenotype()
	 *  except that you can limit the initialization to certain \e loci and \e ploidy.
	 */
	InitGenotype(const floatListFunc & freq = vectorf(),
		const uintList & genotype = vectoru(),
		const floatListFunc & prop = vectorf(),
		const intMatrix & haplotypes = intMatrix(),
		const intMatrix & genotypes = intMatrix(),
		const lociList & loci = lociList(),
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
	vectorf getFreqOrProp(size_t loc, const vspID & vsp) const;

private:
	/// allele frequencies (assume all loci are the same for a subPop
	const floatListFunc m_freq;
	const vectoru m_genotype;
	const floatListFunc m_prop;
	const matrixi m_haplotypes;
	const matrixi m_genotypes;
	//
	const lociList m_loci;

	//
	const uintList m_ploidy;
};

/** This operator assigns lineages at all or part of loci with given
 *  values. This operator initializes all chromosomes, including unused
 *  lineage locations and customized chromosomes.
 */
class InitLineage : public BaseOperator
{
public:
	/** This function creates an initializer that initializes lineages
	 *  with either a specified set of values or from the field \e infoFields
	 *  (default to \c ind_id), whose value will be saved as the lineage of
	 *  modified alleles. If a list of values is specified in parameter \e
	 *  lineage, each value in this list is applied to one or more alleles so
	 *  that each allele (\c PER_ALLELE, default mode), alleles on each
	 *  chromosome (\c PER_CHROMOSOME), on chromosomes of each ploidy (\c
	 *  PER_PLOIDY), or for each individual (\c PER_INDIVIDUAL) have the same
	 *  lineage. A single value is allowed and values in \e lineage will be
	 *  re-used if not enough values are provided. If an empty list is provided,
	 *  values 1, 2, 3, .. will be used to provide an unique identify for each
	 *  allele, genotype, chromosome, etc. If a valid field is specified (default
	 *  to \c ind_id), the value of this field will be used for all alleles of
	 *  each individual if \e mode is set to  \c FROM_INFO, or be adjusted to
	 *  produce positive values for alleles on the frist ploidy, and negative
	 *  values for the second ploidy (and so on) if \e mode equals to
	 *  \c FROM_INFO_SIGNED. If \e loci, \e ploidy and/or \e subPops are
	 *  specified, only specified loci, ploidy, and individuals in these
	 *  (virtual) subpopulations will be initialized.
	 */
	InitLineage(const intList & lineage = vectori(), int mode = PER_ALLELE,
		const lociList & loci = lociList(), const uintList & ploidy = uintList(),
		int begin = 0, int end = 1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr(1, "ind_id"));

	~InitLineage()
	{
	}


	/// HIDDEN Deep copy of the operator \c InitLineage
	virtual BaseOperator * clone() const
	{
		return new InitLineage(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;


	/// HIDDEN apply this operator to population \e pop
	bool apply(Population & pop) const;

private:
	const vectori m_lineage;

	//
	const lociList m_loci;

	//
	const uintList m_ploidy;

	//
	const int m_mode;
};


}
#endif
