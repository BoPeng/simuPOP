/**
 *  $File: qtrait.h $
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

#ifndef _QTRAIT_H
#define _QTRAIT_H
/**
   \file
   \brief head file of class selector:public baseOperator
 */
#include "utility.h"
#include "operator.h"

#include "boost/tuple/tuple.hpp"
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
class baseQuanTrait : public baseOperator
{
public:
	/** Create a base quantitative trait operator. If \e ancGen=0 (default),
	 *  only the current generation will be applied. If \e ancGen=-1, the trait
	 *  fields (\e infoFields) of all ancestral generations will be set. If a
	 *  positive number is given, ancestral generations with index <= ancGen
	 *  will be applied. A quantitative trait operator can be applied to
	 *  specified (virtual) subpopulations (parameter \e subPops) and
	 *  replicates (parameter \e reps).
	 */
	baseQuanTrait(int ancGen = -1,  int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr())
		: baseOperator("", begin, end, step, at, reps, subPops, infoFields),
		m_ancGen(ancGen)
	{
		DBG_ASSERT(infoSize() >= 1, ValueError,
			"Please specify at least one quantitative trait field");
	}


	/// destructor
	virtual ~baseQuanTrait()
	{
	}


	/// deep copy of a quantitative trait operator
	virtual baseOperator * clone() const
	{
		return new baseQuanTrait(*this);
	}


	/** CPPONLY
	 *  calculate/return quantitative trait etc.
	 */
	virtual void qtrait(individual *, ULONG gen, vectorf & traits)
	{
		///
		throw ValueError("This quantitative trait calculator is not supposed to be called directly");
	}


	/// set \c qtrait to all individual
	bool apply(population & pop);

	/// CPPONLY
	bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// HIDDEN
	string describe(bool format = true)
	{
		return "<simuPOP.baseQuanTrait> quantitative trait>" ;
	}


private:
	/// how to handle ancestral gen
	int m_ancGen;

};


/** This quantitative trait operator assigns a trait field by calling a user
 *  provided function. It accepts a list of loci (parameter \e loci), and
 *  a Python function \c func which should be defined with one or more of
 *  parameters \c geno, \c gen, \c ind, or names of information fields. When
 *  this operator is applied to a population, it passes genotypes at specified
 *  loci, generation number, a reference to an individual, and values at
 *  specified information fields to respective parameters of this function.
 *  The return values will be assigned to specified trait fields.
 *
 *  <funcForm>pyQuanTrait</funcForm>
 */
class PyQuanTrait : public baseQuanTrait
{
public:
	/** Create a Python hybrid quantitative trait operator that passes genotype
	 *  at specified \e loci, optional values at specified information fields
	 *  (if requested), and an optional generation number to a user-defined
	 *  function \e func. The return value will be assigned to specified trait
	 *  fields (\e infoField). If only one trait field is specified, a number
	 *  or a sequence of one element is acceptable. Otherwise, a sequence of
	 *  values will be accepted and be assigned to each trait field.
	 */
	PyQuanTrait(PyObject * func, const uintList & loci = vectoru(),
		int ancGen = 0, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		baseQuanTrait(ancGen, begin, end, step, at, reps, subPops, infoFields),
		m_func(func), m_loci(loci.elems())
	{
		DBG_ASSERT(m_func.isValid(), ValueError, "Passed variable is not a callable python function.");

		DBG_FAILIF(m_loci.empty(), ValueError,
			"Please specify susceptibility loci");
	};

	/// deep copy of a Python quantitative trait operator
	virtual baseOperator * clone() const
	{
		return new PyQuanTrait(*this);
	}


	/** CPPONLY
	 *  currently assuming diploid
	 */
	virtual void qtrait(individual * ind, ULONG gen, vectorf & traits);

	/// HIDDEN
	string describe(bool format = true)
	{
		return "<simuPOP.PyQuanTrait> a hybrid quantitative trait model";
	}


private:
	/// user supplied python function
	pyFunc m_func;

	/// susceptibility loci
	vectoru m_loci;
};

}
#endif
