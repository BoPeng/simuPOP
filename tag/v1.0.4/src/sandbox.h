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

#include <set>

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
	 *  \li a Python function, which will be called when selection coefficient
	 *      of a new mutant is needed. Mutant location will be passed to this
	 *      function if it accepts a parameter \c loc. This allows the
	 *      definition of site-specific selection coefficients.
	 *  Individual fitness (1+s_i) will be combined in \c ADDITIVE,
	 *     \c MULTIPLICATIVE or \c EXPONENTIAL mode. (See \c MlSelector for
	 *     details).
	 *  If an output is given, mutants and their fitness values will be written
	 *  to the output, in the form of 'mutant fitness'.
	 */
	InfSitesSelector(const floatListFunc & selDist, int mode = EXPONENTIAL,
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		BaseSelector(output, begin, end, step, at, reps, subPops, infoFields),
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
	virtual double indFitness(Population & pop, Individual * ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		return "<simuPOP.InfSitesSelector>" ;
	}


	/// CPPONLY
	bool apply(Population & pop) const;

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
	mutable vectoru m_newMutants;
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
	 *  will be recorded to the population (instead of alleles). By default,
	 *  this mutator assumes an infinite-sites mutation model so that mutations
	 *  can happen only at a new locus. Mutations happen at a locus with
	 *  existing mutant will be ignored. Alternatively, if \e model=2,
	 *  all mutations are allowed and if a mutant (allele 1) is mutated, it
	 *  will be mutated to allele 0 (back mutation). If an \e output is given,
	 *  mutants will be outputted in the format of "gen mutant ind type" where
	 *  type is 0 for forward (0->1), 1 for backward (1->0) and 2 for invalid
	 *  (ignored) mutations. The first mode has the advantage that all mutants
	 *  in the simulated population can be traced to a single mutation event.
	 *  If the regions are reasonably wide and mutation rates are low, these
	 *  two mutation models should yield similar results.
	 */
	InfSitesMutator(double rate, const intMatrix & ranges, int model = 1,
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_rate(rate), m_ranges(ranges), m_model(model), m_mutants()
	{
		const matrixi & rngs = m_ranges.elems();

		for (size_t i = 0; i < rngs.size(); ++i) {
			DBG_FAILIF(rngs[i].size() != 2, ValueError, "Ranges should have two elements");
			for (size_t j = i + 1; j < rngs.size(); ++j) {
				DBG_FAILIF(rngs[j].size() != 2, ValueError, "Ranges should have two elements");
				if (i == j)
					continue;
				if (rngs[i][0] >= rngs[j][0] && rngs[i][0] <= rngs[j][1])
					throw ValueError("Overlapping ranges are currently not supported because of potential conflict of mutant locations on different chromosomes.");
				if (rngs[i][1] >= rngs[j][0] && rngs[i][1] <= rngs[j][1])
					throw ValueError("Overlapping ranges are currently not supported because of potential conflict of mutant locations on different chromosomes.");
			}
		}
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
	const double m_rate;

	const intMatrix m_ranges;

	const int m_model;

	/// used to record mutated loci
	mutable std::set<ULONG> m_mutants;
};


}
#endif

