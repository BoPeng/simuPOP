/**
 *  $File: sandbox.h $
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

#ifndef _SANDBOX_H
#define _SANDBOX_H
/**
   \file
   \brief head file of module sandbox
 */
#include "selector.h"
#include "transmitter.h"

#if TR1_SUPPORT == 0
#  include <map>
#elif TR1_SUPPORT == 1
#  include <unordered_map>
#else
#  include <tr1/unordered_map>
#endif

#include <set>

namespace simuPOP {
namespace sandbox {

#ifdef LONGALLELE

/** This operator looks into a population in mutational space and revert a mutant
 *  to wildtype allele if it is fixed in the population. If a valid output is
 *  specifieid, fixed alleles will be outputed with a leading generation number.
 */
class RevertFixedSites : public BaseOperator
{
public:
	/** Create an operator to revert alleles at fixed loci from value 1 to 0.
	 *  Parameter \e subPops is ignored.
	 */
	RevertFixedSites(const stringFunc & output = "", int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr())
		: BaseOperator(output, begin, end, step, at, reps, subPops, infoFields)
	{
	}


	/// destructor
	virtual ~RevertFixedSites()
	{
	}


	/// HIDDEN Deep copy of a Migrator
	virtual BaseOperator * clone() const
	{
		return new RevertFixedSites(*this);
	}


	/// HIDDEN apply the Migrator to populaiton \e pop.
	virtual bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "Revert fixed alleles to wildtype allele if it is fixed in the population.";
	}


};

/** This selector assumes that alleles are mutant locations in the mutational
 *  space and assign fitness values to them according to a random distribution.
 *  The overall individual fitness is determined by either an additive, an
 *  multiplicative or an exponential model.
 */
class MutSpaceSelector : public BaseSelector
{
public:
	/** Create a selector that assigns individual fitness values according to
	 *  random fitness effects. \e selDist can be
	 *  \li <tt>(CONSTANT, s, h)</tt> where s will be used for all mutants. The
	 *      fitness value for genotypes AA, Aa and aa will be (1, 1-hs, 1-s).
	 *      If h is unspecified, a default value h=0.5 (additive model) will
	 *      be used.
	 *  \li <tt>(GAMMA_DISTRIBUTION, theta, k, h</tt> where s follows a gamma
	 *      distribution with scale parameter theta and shape parameter k.
	 *      Fitness values for genotypes AA, Aa and aa will be 1, 1-hs and 1-s.
	 *      A default value h=0.5 will be used if h is unspecified.
	 *  \li a Python function, which will be called when selection coefficient
	 *      of a new mutant is needed. This function should return a single
	 *      value s (with default value h=0.5) or a sequence of (h, s). Mutant
	 *      location will be passed to this function if it accepts a parameter
	 *      \c loc. This allows the definition of site-specific selection
	 *      coefficients.
	 *  Individual fitness will be combined in \c ADDITIVE,
	 *     \c MULTIPLICATIVE or \c EXPONENTIAL mode. (See \c MlSelector for
	 *     details).
	 *  If an output is given, mutants and their fitness values will be written
	 *  to the output, in the form of 'mutant s h'.
	 */
	MutSpaceSelector(const floatListFunc & selDist, int mode = EXPONENTIAL,
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("fitness")) :
		BaseSelector(output, begin, end, step, at, reps, subPops, infoFields),
		m_selDist(selDist), m_mode(mode), m_selFactory(), m_additive(true)
	{
		if (m_selDist.size() == 0) {
			DBG_FAILIF(!m_selDist.func().isValid(), ValueError,
				"Please specify either a distribution with parameter or a function.");
		} else if (static_cast<int>(m_selDist[0]) == CONSTANT) {
			DBG_FAILIF(m_selDist.size() < 2, ValueError, "At least one parameter is needed for constant selection coefficient.");
		} else if (static_cast<int>(m_selDist[0]) == GAMMA_DISTRIBUTION) {
			DBG_FAILIF(m_selDist.size() < 3, ValueError, "At least two parameters are needed for gamma distribution.");
		}
	}


	virtual ~MutSpaceSelector()
	{
	}


	/// HIDDEN Deep copy of a map selector
	virtual BaseOperator * clone() const
	{
		return new MutSpaceSelector(*this);
	}


	/** CPPONLY
	 *  calculate/return the fitness value, currently assuming diploid
	 */
	virtual double indFitness(Population & pop, RawIndIterator ind) const;

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.MutSpaceSelector>" ;
	}


	/// CPPONLY
	bool apply(Population & pop) const;

	typedef std::pair<double, double> SelCoef;

private:
	SelCoef getFitnessValue(size_t mutant) const;


	double randomSelAddFitness(GenoIterator it, GenoIterator it_end, bool maleChrX) const;

	double randomSelExpFitness(GenoIterator it, GenoIterator it_end, bool maleChrX) const;

	// extended models does not assume additivity (h != 0.5)
	double randomSelMulFitnessExt(GenoIterator it, GenoIterator it_end, bool maleChrX) const;

	double randomSelAddFitnessExt(GenoIterator it, GenoIterator it_end, bool maleChrX) const;

	double randomSelExpFitnessExt(GenoIterator it, GenoIterator it_end, bool maleChrX) const;

private:
	///
	floatListFunc m_selDist;

	int m_mode;
	///
#  if TR1_SUPPORT == 0
	typedef std::map<unsigned int, SelCoef> SelMap;
	typedef std::map<unsigned int, int> MutCounter;
#  else
	// this is faster than std::map
	typedef std::tr1::unordered_map<size_t, SelCoef> SelMap;
	typedef std::tr1::unordered_map<size_t, size_t> MutCounter;
#  endif
	mutable SelMap m_selFactory;
	mutable vectoru m_newMutants;
	// whether or not all markers are additive.
	mutable bool m_additive;
};


/** This is an infite site mutation model in mutational space. The alleles
 *  in the population is assumed to be locations of mutants. A mutation
 *  rate is given that mutate alleles in 'regions'. If number of mutants
 *  for an individual exceed the number of loci, 10 loci will be added
 *  to everyone in the population.
 */
class MutSpaceMutator : public BaseOperator
{
public:
	/** This operator accepts a list of ranges which is the 'real range' of
	 *  each chromosome. Mutation happens with muation rate \e rate and mutants
	 *  will be recorded to the population (instead of alleles). By default,
	 *  this mutator assumes an finite-allele model where all mutations are
	 *  allowed and if a mutant (allele 1) is mutated, it will be mutated to
	 *  allele 0 (back mutation). Alternatively (\e model = 2), an
	 *  infinite-sites mutation model can be used where mutations can happen
	 *  only at a new locus. Mutations happen at a locus with existing mutants
	 *  will be moved to a random locus without existing mutant. A warning
	 *  message will be printed if there is no vacant locus available. If a
	 *  valid \e output is given, mutants will be outputted in the format of
	 *  "gen mutant ind type" where type is 0 for forward (0->1), 1 for
	 *  backward (1->0), 2 for relocated mutations, and 3 for ignored mutation
	 *  because no vacent locus is available. The second mode  has the
	 *  advantage that all mutants in the simulated population can be traced
	 *  to a single mutation event. If the regions are reasonably wide and
	 *  mutation rates are low, these two mutation models should yield
	 *  similar results.
	 */
	MutSpaceMutator(double rate, const intMatrix & ranges, int model = 1,
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_rate(rate), m_ranges(ranges), m_model(model)
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
#  ifdef BINARYALLELE
		DBG_FAILIF(true, ValueError, "This operator does not work in binary allele type.");
#  endif
	}


	/// destructor.
	~MutSpaceMutator()
	{
	}


	virtual bool apply(Population & pop) const;

	/// HIDDEN Deep copy of a \c MutSpaceMutator
	virtual BaseOperator * clone() const
	{
		return new MutSpaceMutator(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.MutSpaceMutator>";
	}


private:
	size_t locateVacantLocus(Population & pop, size_t beg, size_t end, std::set<size_t> & mutants) const;

private:
	const double m_rate;

	const intMatrix m_ranges;

	const int m_model;
};


/** This during mating operator recombine chromosomes, which records mutant
 *  locations, using a fixed recombination rate (per base pair).
 */
class MutSpaceRecombinator : public GenoTransmitter
{
public:
	/** Create a Recombinator (a mendelian genotype transmitter with
	 *  recombination and gene conversion) that passes genotypes from parents
	 *  (or a parent in case of self-fertilization) to offspring. A
	 *  recombination \e rate in the unit of base pair is needed.
	 */
	MutSpaceRecombinator(double rate, const intMatrix & ranges,
		const stringFunc & output = "", int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr())
		: GenoTransmitter(output, begin, end, step, at, reps, subPops, infoFields),
		m_rate(rate), m_ranges(ranges)
	{
		DBG_FAILIF(rate > 0.5 || rate < 0, ValueError, "Recombination rate should be between 0 and 0.5");
#  ifdef BINARYALLELE
		DBG_FAILIF(true, ValueError, "This operator does not work in binary allele type.");
#  endif
	}


	/// HIDDEN Deep copy of a Recombinator
	virtual BaseOperator * clone() const
	{
		return new MutSpaceRecombinator(*this);
	}


	virtual ~MutSpaceRecombinator()
	{
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.MutSpaceRecombinator>";
	}


	/** CPPONLY
	 *  Apply the Recombinator during mating
	 */
	virtual bool applyDuringMating(Population & pop, Population & offPop,
		RawIndIterator offspring,
		Individual * dad, Individual * mom) const;

private:
#  if TR1_SUPPORT == 0
	typedef std::map<unsigned int, int> MutCounter;
#  else
	// this is faster than std::map
	typedef std::tr1::unordered_map<unsigned int, int> MutCounter;
#  endif
	// use when m_rate = 0.5
	void transmitGenotype0(Population & pop, Population & offPop, const Individual & parent,
		size_t offIndex, int ploidy) const;

	// use when m_rate < 1e-4
	void transmitGenotype1(Population & pop, Population & offPop, const Individual & parent,
		size_t offIndex, int ploidy) const;

private:
	/// recombination rate
	const double m_rate;
	const intMatrix m_ranges;
};

#endif
}
}
#endif

