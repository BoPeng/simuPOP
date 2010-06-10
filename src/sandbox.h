/**
 *  $File: sandbox.h $
 *  $LastChangedDate: 2010-06-04 13:29:09 -0700 (Fri, 04 Jun 2010) $
 *  $Rev: 3579 $
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

#ifndef _SANDBOX_H
#define _SANDBOX_H
/**
   \file
   \brief head file of module sandbox
 */
#include "selector.h"


namespace simuPOP {

/** This selector assigns individual fitness values using a user-specified
 *  dictionary. This operator can be applied to populations with arbitrary
 *  number of homologous chromosomes.
 */
class InfSitesSelector : public BaseSelector
{
public:
	/** Create a selector that assigns individual fitness values using a
	 *  dictionary \e fitness with genotype at \e loci as keys, and fitness
	 *  as values. Parameter \e loci can be a list of indexes, loci names or
	 *  \c ALL_AVAIL. For each individual (parents if this operator is applied
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
	InfSitesSelector(const lociList & loci, const tupleDict & fitness,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		BaseSelector(begin, end, step, at, reps, subPops, infoFields),
		m_loci(loci), m_dict(fitness)
	{
	};

	virtual ~InfSitesSelector()
	{
	}


	/// HIDDEN Deep copy of a map selector
	virtual BaseOperator * clone() const
	{
		return new InfSitesSelector(*this);
	}


	/** CPPONLY
	 *  calculate/return the fitness value, currently assuming diploid
	 */
	virtual double indFitness(Individual * ind, ULONG gen) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		return "<simuPOP.InfSitesSelector>" ;
	}


private:
	///
	const lociList m_loci;

	/// fitness for each genotype
	const tupleDict m_dict;
};

}
#endif
