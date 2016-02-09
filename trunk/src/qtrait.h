/**
 *  $File: qtrait.h $
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

#ifndef _QTRAIT_H
#define _QTRAIT_H
/**
   \file
   \brief head file of class selector:public BaseOperator
 */
#include "utility.h"
#include "operator.h"

#include "boost_pch.hpp"
#include <numeric>
using std::min;

namespace simuPOP {

/** A quantitative trait in simuPOP is simply an information field. A
 *  quantitative trait model simply assigns values to one or more information
 *  fields (called trait fields) of each individual according to its genetic
 *  (genotype) and environmental (information field) factors. It can be applied
 *  at any stage of an evolutionary cycle. If a quantitative trait operator is
 *  applied before or after mating, it will set the trait fields of all parents
 *  and offspring. If it is applied during mating, it will set the trait fields
 *  of each offspring.
 *
 *  When a quantitative trait operator is applied to a population, it is only
 *  applied to the current generation. You can, however, use parameter
 *  \e ancGen=-1 to set the trait field of all ancestral generations, or a
 *  generation index to apply to only ancestral generation younger than
 *  \e ancGen. Note that this parameter is ignored if the operator is applied
 *  during mating.
 */
class BaseQuanTrait : public BaseOperator
{
public:
	/** Create a base quantitative trait operator. This operator assigns one
	 *  or more quantitative traits to trait fields in the present generation
	 *  (default). If \c ALL_AVAIL or a list of ancestral generations are
	 *  specified, this operator will be applied to individuals in these
	 *  generations as well. A quantitative trait operator can be applied to
	 *  specified (virtual) subpopulations (parameter \e subPops) and
	 *  replicates (parameter \e reps).
	 */
	BaseQuanTrait(const uintList & ancGens = uintList(NULL),  int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr())
		: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
		m_ancGens(ancGens)
	{
		DBG_ASSERT(infoSize() >= 1, ValueError,
			"Please specify at least one quantitative trait field");
	}


	/// destructor
	virtual ~BaseQuanTrait()
	{
	}


	/// HIDDEN Deep copy of a quantitative trait operator
	virtual BaseOperator * clone() const
	{
		return new BaseQuanTrait(*this);
	}


	/** CPPONLY
	 *  calculate/return quantitative trait etc.
	 */
	virtual void qtrait(Individual *, size_t /* gen */, vectorf & /* traits */) const
	{
		///
		throw ValueError("This quantitative trait calculator is not supposed to be called directly");
	}


	/// set \c qtrait to all individual
	bool apply(Population & pop) const;

	/// CPPONLY
	bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.BaseQuanTrait> quantitative trait>" ;
	}


private:
	/// how to handle ancestral gen
	const uintList m_ancGens;

};


/** This quantitative trait operator assigns a trait field by calling a user
 *  provided function. It accepts a list of loci (parameter \e loci), and
 *  a Python function \c func which should be defined with one or more of
 *  parameters \c geno, \c mut, \c gen, \c ind, or names of information fields.
 *  When this operator is applied to a population, it passes genotypes or
 *  mutants (non-zero alleles) of each individual at specified loci,
 *  generation number, a reference to an individual, and values at
 *  specified information fields to respective parameters of this function.
 *  Genotypes of each individual are passed as a tuple of alleles arranged
 *  locus by locus (in the order of A1,A2,B1,B2 for loci A and B). Mutants are
 *  passed as a default dictionary of loci index (with respect to all genotype of
 *  individuals, not just the first ploidy) and alleles. The return values
 *  will be assigned to specified trait fields.
 */
class PyQuanTrait : public BaseQuanTrait
{
public:
	/** Create a Python hybrid quantitative trait operator that passes genotype
	 *  at specified \e loci, optional values at specified information fields
	 *  (if requested), and an optional generation number to a user-defined
	 *  function \e func. Parameter \e loci can be a list of loci indexes,
	 *  names, or \c ALL_AVAIL. The return value will be assigned to specified
	 *  trait fields (\e infoField). If only one trait field is specified, a
	 *  number or a sequence of one element is acceptable. Otherwise, a
	 *  sequence of values will be accepted and be assigned to each trait
	 *  field.
	 */
	PyQuanTrait(PyObject * func, const lociList & loci = vectoru(),
		const uintList ancGens = uintList(NULL), int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		BaseQuanTrait(ancGens, begin, end, step, at, reps, subPops, infoFields),
		m_func(func), m_loci(loci)
	{
		DBG_ASSERT(m_func.isValid(), ValueError, "Passed variable is not a callable python function.");

	};

	/// HIDDEN Deep copy of a Python quantitative trait operator
	virtual BaseOperator * clone() const
	{
		return new PyQuanTrait(*this);
	}


	/** CPPONLY
	 *  currently assuming diploid
	 */
	virtual void qtrait(Individual * ind, size_t gen, vectorf & traits) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.PyQuanTrait> a hybrid quantitative trait model";
	}


private:
	/// user supplied python function
	const pyFunc m_func;

	/// susceptibility loci
	const lociList m_loci;
};

}
#endif
