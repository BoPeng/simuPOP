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

#if TR1_SUPPORT == 0
#  include <map>
#elif TR1_SUPPORT == 1
#  include <unordered_map>
#else
#  include <tr1/unordered_map>
#endif

namespace simuPOP {

/** This selector assumes that alleles are mutant locations in the mutational
 *  space and assign fitness values to them according to a random distribution.
 *  The overall individual fitness is determined by either an additive, an
 *  multiplicative or an exponential model.
 */
class InfSitesSelector : public BaseSelector
{
public:
	/** Create a selector that assigns individual fitness values according to
	 *  random fitness effects. \e selDist can be
	 *  \li <tt>(CONSTANT, s)</tt> where s will be used for all mutants.
	 *  \li <tt>(GAMMA_DISTRIBUTION, theta, k</tt> where -s (note the - sign)
	 *      follows a gamma distribution with scale parameter theta and shape
	 *      parameter k.
	 *  \li a Python function, which will be called when fitness value of a new
	 *      mutant is needed.
	 *  Individual fitness (1+s_i) will be combined in \c ADDITIVE,
	 *     \c MULTIPLICATIVE or \c EXPONENTIAL mode. (See \c MlSelector for
	 *     details).
	 */
	InfSitesSelector(const floatListFunc & selDist, int mode = EXPONENTIAL,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		BaseSelector(begin, end, step, at, reps, subPops, infoFields),
		m_selDist(selDist), m_mode(mode), m_selFactory()
	{
		if (m_selDist.size() == 0) {
			DBG_FAILIF(!m_selDist.func().isValid(), ValueError,
				"Please specify either a distribution with parameter or a function.");
		} else if (static_cast<int>(m_selDist[0]) == CONSTANT) {
			DBG_FAILIF(m_selDist.size() != 2, ValueError, "One parameters are needed for gamma distribution.");
		} else if (static_cast<int>(m_selDist[0]) == GAMMA_DISTRIBUTION) {
			DBG_FAILIF(m_selDist.size() != 3, ValueError, "Two parameters are needed for gamma distribution.");
		}
	}


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


	/** Return a dictionary of selection coefficient for each mutant.
	 */
	intDict selCoef() const;

private:
	double getFitnessValue(int mutant) const;

	double randomSelMulFitness(GenoIterator it, GenoIterator it_end) const;

	double randomSelAddFitness(GenoIterator it, GenoIterator it_end) const;

	double randomSelExpFitness(GenoIterator it, GenoIterator it_end) const;

private:
	///
	floatListFunc m_selDist;

	///
	int m_mode;

#if TR1_SUPPORT == 0
	typedef std::map<unsigned int, double> SelMap;
#else
	// this is faster than std::map
	typedef std::tr1::unordered_map<unsigned int, double> SelMap;
#endif
	mutable SelMap m_selFactory;
};


/** This is an infite site mutation model in mutational space. The alleles
 *  in the population is assumed to be locations of mutants. A mutation
 *  rate is given that mutate alleles in 'regions'. If number of mutants
 *  for an individual exceed the number of loci, 10 loci will be added
 *  to everyone in the population.
 */
class InfSitesMutator : public BaseOperator
{
public:
	/** This operator accepts a list of ranges which is the 'real range' of
	 *  each chromosome. Mutation happens with muation rate \e rate and mutants
	 *  will be recorded to the population (instead of alleles).
	 */
	InfSitesMutator(double rate, const intMatrix & ranges,
		const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_rate(rate), m_ranges(ranges)
	{
#ifdef BINARYALLELE
		DBG_FAILIF(true, ValueError, "This operator does not work in binary allele type.");
#endif
	}


	/// destructor.
	~InfSitesMutator()
	{
	}


	virtual bool apply(Population & pop) const;

	/// HIDDEN Deep copy of a \c InfSitesMutator
	virtual BaseOperator * clone() const
	{
		return new InfSitesMutator(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		return "<simuPOP.InfSitesMutator>";
	}


private:
	double m_rate;

	intMatrix m_ranges;

};


}
#endif

