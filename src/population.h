/**
 *  $File: population.h $
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
#ifndef _POPULATION_H
#define _POPULATION_H

/**
   \file
   \brief head file of class Population
 */

#include "utility.h"
#include <vector>
#include <iterator>
#include <numeric>
using std::vector;
using std::accumulate;
using std::pair;

#include <functional>
using std::equal_to;

// used to save history Population
// 0 (first parental) 1, ...., n
#include <deque>
using std::deque;

#include "boost_pch.hpp"
#include "individual.h"
#include "virtualSubPop.h"


namespace simuPOP {


/** this class implements a Python itertor class that can be used to iterate
 *  through individuals in a (sub)population. If allInds are true,
 *  visiblility of individuals will not be checked. Otherwise, a functor
 *  will be used to check if indiviudals belong to a specified virtual
 *  subpopulation.
 *
 *  An instance of this class is returned by
 *  population::Individuals() and Population::Individuals(subPop)
 */
class pyIndIterator
{
public:
	pyIndIterator(vector<Individual>::iterator const begin,
		vector<Individual>::iterator const end,
		bool allInds, vspFunctor func) :
		m_begin(begin),
		m_index(begin),
		m_end(end),
		m_allInds(allInds),
		m_isIteratable(func)
	{
		// m_index does not have to be pointed to the first
		// valid individual. the next() function will return
		// so.
	}


	~pyIndIterator()
	{
	}


	pyIndIterator __iter__()
	{
		return *this;
	}


	// python 2.x uses next()
	Individual & next()
	{
		// this is the easy (and faster) case
		if (m_allInds) {
			if (m_index == m_end)
				throw StopIteration("");
			else
				return *m_index++;
		}
		// check the visibility of individuals
		do {
			if (m_index == m_end)
				throw StopIteration("");
			else if (m_isIteratable(m_index - m_begin)) {
				return *m_index++;
			} else
				++m_index;
		} while (true);
	}


	// python 3.x uses __next__ instead of next.
	Individual & __next__()
	{
		return next();
	}


private:
	// current (initial individual)
	vector<Individual>::iterator m_begin;

	// current (initial individual)
	vector<Individual>::iterator m_index;

	// ending idx
	vector<Individual>::iterator m_end;

	//
	bool m_allInds;
	//
	vspFunctor m_isIteratable;
};


class Pedigree;


/**
 *  A simuPOP population consists of individuals of the same genotypic
 *  structure, organized by generations, subpopulations and virtual
 *  subpopulations. It also contains a Python dictionary that is used to
 *  store arbitrary population variables.
 *
 *  In addition to genotypic structured related functions provided by the
 *  \c GenoStruTrait class, the population class provides a large number
 *  of member functions that can be used to
 *  \li Create, copy and compare populations.
 *  \li Manipulate subpopulations. A population can be divided into several
 *    subpopulations. Because individuals only mate with individuals within
 *    the same subpopulation, exchange of genetic information across
 *    subpopulations can only be done through migration. A number of functions
 *    are provided to access subpopulation structure information, and to merge
 *    and split subpopulations.
 *  \li Define and access virtual subpopulations. A <em>virtual subpopulation
 *    splitter</em> can be assigned to a population, which defines groups of
 *    individuals called <em>virtual subpopulations</em> (VSP) within each
 *    subpopulation.
 *  \li Access individuals individually, or through iterators that iterate
 *    through individuals in (virtual) subpopulations.
 *  \li Access genotype and information fields of individuals at the population
 *    level. From a population point of view, all genotypes are arranged
 *    sequentially individual by individual. Please refer to class \c Individual
 *    for an introduction to genotype arragement of each individual.
 *  \li Store and access <em>ancestral generations</em>. A population can save
 *    arbitrary number of ancestral generations. It is possible to directly
 *    access an ancestor, or make an ancestral generation the current generation
 *    for more efficient access.
 *  \li Insert or remove loci, resize (shrink or expand) a population, sample
 *    from a population, or merge with other populations.
 *  \li Manipulate population variables and evaluate expressions in this
 *    <em>local namespace</em>.
 *  \li Save and load a population.
 */
class Population : public GenoStruTrait
{
public:
#define HAPLODIPLOID 2.5

	/** @name  constructors and destructor */
	//@{

	/** The following parameters are used to create a population object:
	 *
	 *  \param size A list of subpopulation sizes. The length of this list
	 *    determines the number of subpopulations of this population. If
	 *    there is no subpopulation, <tt>size=[popSize]</tt> can be
	 *    written as <tt>size=popSize</tt>.
	 *  \param ploidy Number of homologous sets of chromosomes. Default to
	 *    \c 2 (diploid). For efficiency considerations, all chromosomes have
	 *    the same number of homologous sets, even if some customized
	 *    chromosomes or some individuals (e.g. males in a haplodiploid
	 *    population) have different numbers of homologous sets. The first
	 *    case is handled by setting \e chromTypes of each chromosome. Only
	 *    the haplodiploid populations are handled for the second case, for
	 *    which <tt>ploidy=HAPLODIPLOID</tt> should be used.
	 *  \param loci A list of numbers of loci on each chromosome. The length of
	 *    this parameter determines the number of chromosomes. If there is only
	 *    one chromosome, \c numLoci instead of <tt>[numLoci]</tt> can be used.
	 *  \param chromTypes A list that specifies the type of each chromosome,
	 *    which can be \c AUTOSOME, \c CHROMOSOME_X, \c CHROMOSOME_Y, or
	 *    \c CUSTOMIZED. All chromosomes are assumed to be autosomes if
	 *    this parameter is ignored. Sex chromosome can only be specified in a
	 *    diploid population where the sex of an individual is determined by
	 *    the existence of these chromosomes using the \c XX (\c FEMALE) and
	 *    \c XY (\c MALE) convention. Both sex chromosomes have to be available
	 *    and be specified only once. Because chromosomes \c X and \c Y are
	 *    treated as two chromosomes, recombination on the pseudo-autosomal
	 *    regions of the sex chromsomes is not supported. \c CUSTOMIZED
	 *    chromosomes are special chromosomes whose inheritance patterns are
	 *    undefined. They rely on user-defined functions and operators to be
	 *    passed from parents to offspring. Multiple customized chromosomes
	 *    have to be arranged consecutively.
	 *  \param lociPos Positions of all loci on all chromosome, as a list of
	 *    float numbers. Default to \c 1, \c 2, ... etc on each chromosome.
	 *    \e lociPos should be arranged chromosome by chromosome. If \c lociPos
	 *    are not in order within a chromosome, they will be re-arranged along
	 *    with corresponding \e lociNames (if specified).
	 *  \param ancGen Number of the most recent ancestral generations to keep
	 *    during evolution. Default to \c 0, which means only the current
	 *    generation will be kept. If it is set to \c -1, all ancestral
	 *    generations will be kept in this population (and exhaust your
	 *    computer RAM quickly).
	 *  \param chromNames A list of chromosome names. Default to \c '' for all
	 *    chromosomes.
	 *  \param alleleNames A list or a nested list of allele names. If a list
	 *    of alleles is given, it will be used for all loci in this population.
	 *    For example, <tt>alleleNames=('A','C','T','G')</tt> gives names \c A,
	 *    \c C, \c T, and \c G to alleles \c 0, \c 1, \c 2, and \c 3
	 *    respectively. If a nested list of names is given, it should specify
	 *    alleles names for all loci.
	 *  \param lociNames A list of names for each locus. It can be empty or a
	 *    list of unique names for each locus. If loci are not specified in
	 *    order, loci names will be rearranged according to their position on
	 *    the chromosome.
	 *  \param subPopNames A list of subpopulation names. All subpopulations
	 *    will have name \c '' if this parameter is not specified.
	 *  \param infoFields Names of information fields (named float number) that
	 *    will be attached to each individual.
	 */
	Population(const uintList & size = vectoru(),
		float ploidy = 2,
		const uintList & loci = vectoru(),
		const uintList & chromTypes = vectoru(),
		const floatList & lociPos = vectorf(),
		int ancGen = 0,
		const stringList & chromNames = vectorstr(),
		const stringMatrix & alleleNames = stringMatrix(),
		const stringList & lociNames = vectorstr(),
		const stringList & subPopNames = vectorstr(),
		const stringList & infoFields = vectorstr());

	/// CPPONLY copy constructor
	Population(const Population & rhs);

	/** Create a cloned copy of a population. Note that Python statement
	 *  <tt>pop1 = pop</tt> only creates a reference to an existing population
	 *  \c pop.
	 *  <group>8-pop</group>
	 */
	Population * clone() const;

	/** Swap the content of two population objects, which can be handy in some
	 *  particular circumstances. For example, you could swap out a population
	 *  in a simulator.
	 *  <group>1-pop</group>
	 */
	void swap(Population & rhs)
	{
		GenoStruTrait::swap(rhs);
		std::swap(m_popSize, rhs.m_popSize);

		m_subPopSize.swap(rhs.m_subPopSize);
		m_subPopNames.swap(rhs.m_subPopNames);
		m_subPopIndex.swap(rhs.m_subPopIndex);
		m_genotype.swap(rhs.m_genotype);
#ifdef LINEAGE
		m_lineage.swap(rhs.m_lineage);
#endif
		m_info.swap(rhs.m_info);
		m_inds.swap(rhs.m_inds);
		std::swap(m_ancestralGens, rhs.m_ancestralGens);
		m_vars.swap(rhs.m_vars);
		m_ancestralPops.swap(rhs.m_ancestralPops);
		std::swap(m_curAncestralGen, rhs.m_curAncestralGen);
		std::swap(m_indOrdered, rhs.m_indOrdered);
		std::swap(m_vspSplitter, rhs.m_vspSplitter);
		std::swap(rhs.m_gen, m_gen);
		std::swap(rhs.m_rep, m_rep);
#ifdef MUTANTALLELE
		// vectorm must be setGenoPtr after swap
		GenoIterator ptr = m_genotype.begin();
		for (size_t i = 0; i < m_inds.size(); ++i, ptr += genoSize())
			m_inds[i].setGenoPtr(ptr);
		ptr = rhs.m_genotype.begin();
		for (size_t i = 0; i < rhs.m_inds.size(); ++i, ptr += rhs.genoSize())
			rhs.m_inds[i].setGenoPtr(ptr);
#endif
	}


	/// destroy a population
	~Population();

	/** CPPONLY
	 * Validate if a population is in good shape. This is mostly used
	 * to detect if scratch population is prepared properly during
	 * evolution
	 */
	void validate(const string & msg) const;

	/** CPPONLY
	 * Fix a population, resize it if necessary. The content
	 * of the population will be cleared.
	 */
	void fitSubPopStru(const vectoru & newSubPopSizes,
		const vectorstr & newSubPopNames);

	/** if a population has any activated virtual subpopulations
	 *  CPPONLY
	 */
	bool hasActivatedVirtualSubPop() const;

	/** if a subpopulation has any activated virtual subpopulation
	 *  CPPONLY
	 */
	bool hasActivatedVirtualSubPop(size_t subPop) const;

	/** CPPONLY because this is simply numVirtualSubPop() != 0.
	 *  Return True if virtual subpopulations are defined for this population.
	 *  <group>3-VSP</group>
	 */
	bool hasVirtualSubPop() const;

	/** Return the virtual splitter associated with the population,
	 *  \c None will be returned if there is no splitter.
	 */
	BaseVspSplitter * virtualSplitter() const { return m_vspSplitter; }

	/** Set a VSP \e splitter to the population, which defines the same VSPs
	 *  for all subpopulations. If different VSPs are needed for different
	 *  subpopulations, a \c CombinedSplitter can be used to make these VSPs
	 *  available to all subpopulations.
	 *  <group>3-VSP</group>
	 */
	void setVirtualSplitter(BaseVspSplitter * splitter);

	/** Return the number of virtual subpopulations (VSP) defined by a VSP
	 *  splitter. Return \c 0 if no VSP is defined.
	 *  <group>3-VSP</group>
	 */
	size_t numVirtualSubPop() const;

	/// HIDDEN activate a virtual subpopulation.
	void activateVirtualSubPop(vspID subPop) const;

	/** HIDDEN
	 *  deactivate virtual subpopulations in a given
	 *  subpopulation. In another word, all individuals
	 *  will become visible.
	 */
	void deactivateVirtualSubPop(size_t subPop) const;

	// allow compaison of populations in python
	// only equal or unequal, no greater or less than
	/// a python function used to compare the population objects
	int __cmp__(const Population & rhs) const;

	/** HIDDEN
	 *  adapt the current population to anther population structure.
	 *  population size might or might not be changed.
	 */
	void fitGenoStru(size_t stru);

	/** HIDDEN
	 *  set population/subpopulation structure given subpopulation sizes
	 *  \param newSubPopSizes an array of new subpopulation sizes. The overall
	 *    population size should not changed.
	 *  <group>2-subpop</group>
	 */
	void setSubPopStru(const vectoru & newSubPopSizes, const vectorstr & newSubPopNames);

	/** Return the number of subpopulations in a population. Return 1 if there
	 *  is no subpopulation structure.
	 *  <group>2-subpop</group>
	 */
	size_t numSubPop() const
	{
		return static_cast<size_t>(m_subPopSize.size());
	}


	/** Return the size of a subpopulation (<tt>subPopSize(sp)</tt>) or a
	 *  virtual subpopulation (<tt>subPopSize([sp, vsp])</tt>) in the current
	 *  generation (default) or a specified ancestral generation \e ancGen. If
	 *  no \e subpop is given, it is the same as <tt>popSize(ancGen, sex)</tt>. 
	 *  Population and virtual subpopulation names can be used. This function
	 *  by default returns number of all individuals (<tt>sex=ANY_SEX</tt>), but it
	 *  will return number of males (if <tt>sex=MALE_ONLY</tt>), number of females
	 *  (if <tt>sex=MALE_ONLY</tt>), and number of male/female pairs (if
	 *  <tt>sex=PAIR_ONLY</tt>) which is essentially less of the number of males
	 *  and females.
	 *  <group>2-subpopsize</grouplociList()>
	 */
	size_t subPopSize(vspID subPop = vspID(), int ancGen = -1, SexChoice sex = ANY_SEX) const;


	/** Return the index of the first subpopulation with name \e name. An
	 *  \c IndexError will be raised if subpopulations are not named, or
	 *  if no subpopulation with name \e name is found. Virtual subpopulation
	 *  name is not supported.
	 *  <group>2-subpopname</group>
	 */
	size_t subPopByName(const string & name) const;

	/** Return the "spName - vspName" (virtual named subpopulation), "" (unnamed
	 *  non-virtual subpopulation), "spName" (named subpopulation) or "vspName"
	 *  (unnamed virtual subpopulation), depending on whether subpopulation is
	 *  named or if \e subPop is virtual.
	 *  <group>2-subpopname</group>
	 */
	string subPopName(vspID subPop) const;

	/** Return the names of all subpopulations (excluding virtual
	 *  subpopulations). An empty string will be returned for unnamed
	 *  subpopulations.
	 *  <group>2-subpopname</group>
	 */
	vectorstr subPopNames() const;

	/** Assign a name \e name to subpopulation \e subPop. Note that
	 *  subpopulation names do not have to be unique.
	 *  <group>2-subpopname</group>
	 */
	void setSubPopName(const string & name, size_t subPop);

	/** Return the sizes of all subpopulations at the current generation
	 *  (default) or specified ancestral generation \e ancGen. Virtual
	 *  subpopulations are not considered.
	 *  <group>2-subpopsize</group>
	 */
	vectoru subPopSizes(int ancGen = -1) const
	{
		if (ancGen < 0 || ancGen == m_curAncestralGen)
			return m_subPopSize;
		DBG_FAILIF(ancGen > ancestralGens(), IndexError,
			(boost::format("Ancestral generation %1% out of range of 0 ~ %2%") % ancGen %
			 ancestralGens()).str());
		return m_ancestralPops[ancGen - 1].m_subPopSize;
	}


	//@}

	/** @name indexes, conversion between absoluate indexes and relative indexes. return of chromomosome/subpopulation indexes.
	 */
	//@{

	/** Return the total number of individuals in all subpopulations of the
	 *  current generation (default) or the an ancestral generation \e ancGen.
	 *  This function by default returns number of all individuals 
	 *  (<tt>sex=ANY_SEX</tt>), but it will return number of males  (if
	 *  <tt>sex=MALE_ONLY</tt>), number of females (if <tt>sex=MALE_ONLY</tt>),
	 *  and number of male/female pairs (if <tt>sex=PAIR_ONLY</tt>) which is
	 *  essentially less of the number of males and females.
	 *  <group>2-subpopsize</group>
	 */
	size_t popSize(int ancGen = -1, SexChoice sex = ANY_SEX) const;

	/** return the absolute index of an individual \e idx in subpopulation \e subPop.
	 *  <group>2-subpop</group>
	 */
	size_t absIndIndex(size_t idx, size_t subPop) const
	{
		CHECKRANGESUBPOP(subPop);
		CHECKRANGESUBPOPMEMBER(idx, subPop);

		return m_subPopIndex[subPop] + idx;
	}


	/** Return the subpopulation ID and relative index of an individual, given
	 *  its absolute index \c idx.
	 *  <group>2-subpop</group>
	 */
	pairu subPopIndPair(size_t idx) const
	{
		CHECKRANGEIND(idx);

		pairu loc;

		for (size_t i = 1; i <= m_subPopSize.size(); ++i) {
			if (m_subPopIndex[i] > idx) {
				loc.first = i - 1;
				loc.second = idx - m_subPopIndex[i - 1];
				break;
			}
		}
		return loc;
	}


	/** Return the index of the first individual in subpopulation \e subPop.
	 *  <group>2-subpop</group>
	 */
	size_t subPopBegin(size_t subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return m_subPopIndex[subPop];
	}


	/** Return the index of the last individual in subpopulation \e subPop plus
	 *  \c 1, so that <tt>range(subPopBegin(subPop)</tt>,
	 *  <tt>subPopEnd(subPop)</tt> can iterate through the index of all
	 *  individuals in subpopulation \e subPop.
	 *  <group>2-subpop</group>
	 */
	size_t subPopEnd(size_t subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return m_subPopIndex[subPop + 1];
	}


	//@}
	/** @name itertors and accessers, ways to access information, mainly various iterators.
	 */
	//@{


	/** CPPONLY
	 */
	Individual & individual(size_t idx, vspID subPop = vspID())
	{
#ifndef OPTIMIZED
		DBG_FAILIF(subPop.isVirtual(), ValueError,
			"Function Individual currently does not support virtual subpopulation");

		if (!subPop.valid()) {
			CHECKRANGEIND(idx);
		} else {
			CHECKRANGESUBPOPMEMBER(idx, subPop.subPop());
		}
#endif
		return subPop.valid() ? m_inds[subPopBegin(subPop.subPop()) + idx] : m_inds[idx];
	}


	/** CPPONLY
	 */
	const Individual & individual(size_t idx, vspID subPop = vspID()) const
	{
#ifndef OPTIMIZED
		DBG_FAILIF(subPop.isVirtual(), ValueError,
			"Function Individual currently does not support virtual subpopulation");

		if (!subPop.valid()) {
			CHECKRANGEIND(idx);
		} else {
			CHECKRANGESUBPOPMEMBER(idx, subPop.subPop());
		}
#endif
		return subPop.valid() ? m_inds[subPopBegin(subPop.subPop()) + idx] : m_inds[idx];
	}


	/** Return a refernce to individual \e idx in the population
	 * (if <tt>subPop=[]</tt>, default) or a subpopulation (if
	 * <tt>subPop=sp</tt>). Virtual subpopulation is not supported. Note that a
	 * float \e idx is acceptable as long as it rounds closely to an integer.
	 * <group>4-ind</group>
	 */
	Individual & individual(double idx, vspID subPop = vspID())
	{
		size_t intIdx = toID(idx);

		DBG_FAILIF(fabs(idx - double(intIdx)) > 1e-8, ValueError,
			"individual index has to be integer (or a double round to full iteger).");
#ifndef OPTIMIZED
		DBG_FAILIF(subPop.isVirtual(), ValueError,
			"Function Individual currently does not support virtual subpopulation");

		if (!subPop.valid()) {
			CHECKRANGEIND(intIdx);
		} else {
			CHECKRANGESUBPOPMEMBER(intIdx, subPop.subPop());
		}
#endif
		return subPop.valid() ? m_inds[subPopBegin(subPop.subPop()) + intIdx] : m_inds[intIdx];
	}


	/** Return a reference to individual with \e id stored in information
	 *  field \e idField (default to \c ind_id). This function by default
	 *  search the present and all ancestral generations (\c ancGen=ALL_AVAIL),
	 *  but you can limit the search in specific generations if you know which
	 *  generations to search (<tt>ancGens=[0,1]</tt> for present and
	 *  parental generations) or \c UNSPECIFIED to search only the current
	 *  generation. If no individual with \e id is found, an \c IndexError will
	 *  be raised. A float \e id is acceptable as long as it rounds closely to
	 *  an integer. Note that this function uses a dynamic searching algorithm
	 *  which tends to be slow. If you need to look for multiple individuals
	 *  from a static population, you might want to convert a population
	 *  object to a pedigree object and use function <tt>Pedigree.indByID</tt>.
	 *  <group>4-ind</group>
	 */
	Individual & indByID(double id, const uintList & ancGens = uintList(), const string & idField = "ind_id");

	/** CPPONLY: const version of the ind function.
	 */
	const Individual & individual(double idx, vspID subPop = vspID()) const
	{
		size_t intIdx = toID(idx);

		DBG_FAILIF(fabs(idx - double(intIdx)) > 1e-8, ValueError,
			"individual index has to be integer (or a double round to full iteger).");
#ifndef OPTIMIZED
		DBG_FAILIF(subPop.isVirtual(), ValueError,
			"Function Individual currently does not support virtual subpopulation");

		if (!subPop.valid()) {
			CHECKRANGEIND(intIdx);
		} else {
			CHECKRANGESUBPOPMEMBER(intIdx, subPop.subPop());
		}
#endif
		return subPop.valid() ? m_inds[subPopBegin(subPop.subPop()) + intIdx] : m_inds[intIdx];
	}


	/** Return a reference to individual \c idx in ancestral generation \c gen.
	 *  The correct individual will be returned even if the current generation
	 *  is not the present one (see also \c useAncestralGen). If a valid
	 *  \e subPop is specified, \e index is relative to that \e subPop.
	 *  Virtual subpopulation is not supported. Note that a float \e idx is
	 *  acceptable as long as it rounds closely to an integer.
	 *  <group>6-ancestral</group>
	 */
	Individual & ancestor(double idx, ssize_t gen, vspID subPop = vspID());

	/** CPPONLY const version of ancestor().
	 *  <group>6-ancestral</group>
	 */
	const Individual & ancestor(double idx, ssize_t gen, vspID subPop = vspID()) const;

	/** Return an iterator that can be used to iterate through all individuals
	 *  in a population (if <tt>subPop=[]</tt>, default), or a (virtual)
	 *  subpopulation (<tt>subPop=spID</tt> or <tt>(spID, vspID)</tt>). If you
	 *  would like to iterate through multiple subpopulations in multiple
	 *  ancestral generations, please use function \c Population.allIndividuals().
	 *  <group>4-ind</group>
	 */
	pyIndIterator individuals(vspID subPop = vspID());

	/// CPPONLY
	bool indOrdered() const
	{
		return m_indOrdered;
	}


	/// CPPONLY
	void setIndOrdered(bool s) const
	{
		m_indOrdered = s;
	}


	/// CPPONLY Individual iterator: without subPop info
	IndIterator indIterator()
	{
		return IndIterator(m_inds.begin(), m_inds.end(), !hasActivatedVirtualSubPop());
	}


	/** CPPONLY Individual iterator: with subPop info.
	 *  The iterator will skip invisible Individuals
	 */
	IndIterator indIterator(size_t subPop)
	{
		CHECKRANGESUBPOP(subPop);

		return IndIterator(m_inds.begin() + m_subPopIndex[subPop],
			m_inds.begin() + m_subPopIndex[subPop + 1], !hasActivatedVirtualSubPop(subPop));
	}


#ifdef _OPENMP
	/** CPPONLY Individual iterator: with subPop info and Thread ID.
	 *  The iterator will skip invisible Individuals
	 */
	IndIterator indIterator(size_t subPop, size_t threadID)
	{
		CHECKRANGESUBPOP(subPop);
		DBG_FAILIF(threadID >= numThreads(), RuntimeError,
			(boost::format("Thread ID %1% execeed total number of threads %2%") % threadID % numThreads()).str());
		size_t blockSize = m_subPopSize[subPop] / numThreads();
		if (threadID + 1 != numThreads())
			return IndIterator(m_inds.begin() + m_subPopIndex[subPop] + blockSize * threadID,
				m_inds.begin() + m_subPopIndex[subPop] + blockSize * (threadID + 1),
				!hasActivatedVirtualSubPop(subPop));
		else
			return IndIterator(m_inds.begin() + m_subPopIndex[subPop] + blockSize * threadID,
				m_inds.begin() + m_subPopIndex[subPop + 1], !hasActivatedVirtualSubPop(subPop));
	}


#endif

	/** CPPONLY Individual iterator: without subPop info
	 *  The iterator will skip invisible Individuals
	 */
	ConstIndIterator indIterator() const
	{
		return ConstIndIterator(m_inds.begin(), m_inds.end(), !hasActivatedVirtualSubPop());
	}


	/** CPPONLY Individual iterator: with subPop info.
	 *  The iterator will skip invisible Individuals
	 */
	ConstIndIterator indIterator(size_t subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return ConstIndIterator(m_inds.begin() + m_subPopIndex[subPop],
			m_inds.begin() + m_subPopIndex[subPop + 1], !hasActivatedVirtualSubPop(subPop));
	}


	/** CPPONLY Individual iterator: without subPop info
	 */
	RawIndIterator rawIndBegin()
	{
		return m_inds.begin();
	}


	/** CPPONLY Individual iterator: without subPop info
	 */
	RawIndIterator rawIndEnd()
	{
		return m_inds.end();
	}


	/** CPPONLY Individual iterator: with subPop info.
	 * The iterator will skip invisible Individuals
	 */
	RawIndIterator rawIndBegin(size_t subPop)
	{
		CHECKRANGESUBPOP(subPop);

		return m_inds.begin() + m_subPopIndex[subPop];
	}


	/** CPPONLY Individual iterator: with subPop info.
	 */
	RawIndIterator rawIndEnd(size_t subPop)
	{
		CHECKRANGESUBPOP(subPop);

		return m_inds.begin() + m_subPopIndex[subPop + 1];
	}


	/** CPPONLY Individual iterator: without subPop info
	 * The iterator will skip invisible Individuals
	 */
	ConstRawIndIterator rawIndBegin() const
	{
		return m_inds.begin();
	}


	/** CPPONLY Individual iterator: without subPop info
	 */
	ConstRawIndIterator rawIndEnd() const
	{
		return m_inds.end();
	}


	/** CPPONLY Individual iterator: with subPop info.
	 * The iterator will skip invisible Individuals
	 */
	ConstRawIndIterator rawIndBegin(size_t subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return m_inds.begin() + m_subPopIndex[subPop];
	}


	/** CPPONLY Individual iterator: with subPop info.
	 */
	ConstRawIndIterator rawIndEnd(size_t subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return m_inds.begin() + m_subPopIndex[subPop + 1];
	}


	/// CPPONLY allele iterator that access a locus across all copies of chromosomes and Individual
	/**
	   \param locus allele access, given locus, return the first allele. ptr++ go the next one.
	   Default return the beginning of the first subpopulation, also the first of the whole population

	   \param order order = True: indiviudls in order
	       order = false: do not even respect subpops

	   \note The order of alleles DOES NOT HAVE TO match the order of individuals. Only the boundary of
	   subpopulations will be respected.	Therefore, it is possible to access all alleles within an
	   subpopulation	through such iterators.
	 */
	IndAlleleIterator alleleIterator(size_t locus);

	/// CPPONLY allele begin, for given subPop
	IndAlleleIterator alleleIterator(size_t locus, size_t subPop);

	/// CPPONLY
	ConstIndAlleleIterator alleleIterator(size_t locus) const;


	/// CPPONLY allele begin, for given subPop
	ConstIndAlleleIterator alleleIterator(size_t locus, size_t subPop) const;

#ifdef LINEAGE
	/// CPPONLY lineage begin
	IndLineageIterator lineageIterator(size_t locus);

	/// CPPONLY lineage begin, for given subPop
	IndLineageIterator lineageIterator(size_t locus, size_t subPop);

#endif

	///  CPPONLY allele iterator, go through all allels one by one, without subPop info
	/**
	   if order, in order
	   otherwise, do not even respect subpopulation structure
	 */
	GenoIterator genoBegin(bool order)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This function is not valid with an activated virtual subpopulation");

		if (order && !indOrdered())
			syncIndPointers();

		return m_genotype.begin();
	}


	///  CPPONLY allele iterator
	GenoIterator genoEnd(bool order)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This function is not valid with an activated virtual subpopulation");
		if (order && !indOrdered())
			syncIndPointers();

		return m_genotype.end();
	}


#ifdef LINEAGE

	///  CPPONLY allele iterator, go through all lineages one by one, without subPop info
	/**
	   if order, in order
	   otherwise, do not even respect subpopulation structure
	 */
	LineageIterator lineageBegin(bool order)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This function is not valid with an activated virtual subpopulation");

		if (order && !indOrdered())
			syncIndPointers();

		return m_lineage.begin();
	}


	///  CPPONLY lineage iterator
	LineageIterator lineageEnd(bool order)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This function is not valid with an activated virtual subpopulation");
		if (order && !indOrdered())
			syncIndPointers();

		return m_lineage.end();
	}


#endif  // LINEAGE

	///  CPPONLY allele iterator, go through all allels one by one in a subpopulation
	/**
	   if order, keep order
	   if not order, respect subpopulation structure
	 */
	GenoIterator genoBegin(size_t subPop, bool order)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This function is not valid with an activated virtual subpopulation");
		CHECKRANGESUBPOP(subPop);

		syncIndPointers(order);

		return m_genotype.begin() + m_subPopIndex[subPop] * genoSize();
	}


	/// CPPONLY allele iterator in a subpopulation.
	GenoIterator genoEnd(size_t subPop, bool order)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This function is not valid with an activated virtual subpopulation");
		CHECKRANGESUBPOP(subPop);
		syncIndPointers(order);

		return m_genotype.begin() + m_subPopIndex[subPop + 1] * genoSize();
	}


#ifdef LINEAGE

	///  CPPONLY lineage iterator, go through all lineages one by one in a subpopulation
	/**
	   if order, keep order
	   if not order, respect subpopulation structure
	 */
	LineageIterator lineageBegin(size_t subPop, bool order)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This function is not valid with an activated virtual subpopulation");
		CHECKRANGESUBPOP(subPop);

		syncIndPointers(order);

		return m_lineage.begin() + m_subPopIndex[subPop] * genoSize();
	}


	/// CPPONLY lineage iterator in a subpopulation.
	LineageIterator lineageEnd(size_t subPop, bool order)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This function is not valid with an activated virtual subpopulation");
		CHECKRANGESUBPOP(subPop);
		syncIndPointers(order);

		return m_lineage.begin() + m_subPopIndex[subPop + 1] * genoSize();
	}


#endif  // LINEAGE

	/// CPPONLY genoIterator --- beginning of individual ind.
	GenoIterator indGenoBegin(size_t ind) const
	{
		CHECKRANGEIND(ind);
		return m_inds[ind].genoBegin();
	}


	/// CPPONLY genoIterator -- end of individual ind.
	GenoIterator indGenoEnd(size_t ind) const
	{
		CHECKRANGEIND(ind);
		return m_inds[ind].genoEnd();
	}


#ifdef LINEAGE

	/// CPPONLY lineageIterator --- beginning of individual ind.
	LineageIterator indLineageBegin(size_t ind) const
	{
		CHECKRANGEIND(ind);
		return m_inds[ind].lineageBegin();
	}


	/// CPPONLY lineageIterator -- end of individual ind.
	LineageIterator indLineageEnd(size_t ind) const
	{
		CHECKRANGEIND(ind);
		return m_inds[ind].lineageEnd();
	}


#endif

	/// CPPONLY genoIterator --- beginning of individual ind.
	GenoIterator indGenoBegin(size_t ind, size_t subPop) const
	{
		CHECKRANGESUBPOP(subPop);
		CHECKRANGESUBPOPMEMBER(ind, subPop);

		return m_inds[ subPopBegin(subPop) + ind].genoBegin();
	}


	/// CPPONLY genoIterator -- end of individual ind.
	GenoIterator indGenoEnd(size_t ind, size_t subPop) const
	{
		CHECKRANGESUBPOP(subPop);
		CHECKRANGESUBPOPMEMBER(ind, subPop);

		return m_inds[ subPopBegin(subPop) + ind].genoEnd();
	}


#ifdef LINEAGE

	/// CPPONLY lineageIterator --- beginning of individual ind.
	LineageIterator indLineageBegin(size_t ind, size_t subPop) const
	{
		CHECKRANGESUBPOP(subPop);
		CHECKRANGESUBPOPMEMBER(ind, subPop);

		return m_inds[ subPopBegin(subPop) + ind].lineageBegin();
	}


	/// CPPONLY genoIterator -- end of individual ind.
	LineageIterator indLineageEnd(size_t ind, size_t subPop) const
	{
		CHECKRANGESUBPOP(subPop);
		CHECKRANGESUBPOPMEMBER(ind, subPop);

		return m_inds[ subPopBegin(subPop) + ind].lineageEnd();
	}


#endif

	/** Return an editable array of the genotype of all individuals in
	 *  a population (if <tt>subPop=[]</tt>, default), or individuals in a
	 *  subpopulation \e subPop. Virtual subpopulation is unsupported.
	 *  <group>5-genotype</group>
	 */
	PyObject * genotype(vspID subPop = vspID());

	/** Return an iterator that iterate through mutants of all individuals in
	 *  a population (if <tt>subPop=[]</tt>, default), or individuals in a
	 *  subpopulation \e subPop. Virtual subpopulation is unsupported. Each
	 *  mutant is presented as a tuple of (index, value) where index is the
	 *  index of mutant (from 0 to totNumLoci()*ploidy()) so you will have
	 *  to adjust its value to check multiple alleles at a locus. This
	 *  function ignores type of chromosomes so non-zero alleles in unused
	 *  alleles of sex and mitochondrial chromosomes are also iterated.
	 *  <group>5-genotype</group>
	 */
	pyMutantIterator mutants(vspID subPop = vspID());

	/** Return an editable array of the lineage of alleles for all individuals in
	 *  a population (if <tt>subPop=[]</tt>, default), or individuals in a
	 *  subpopulation \e subPop. Virtual subpopulation is unsupported. <bf>
	 *  This function returns \c None for modules without lineage information.</bf>
	 *  <group>5-genotype</group>
	 */
	PyObject * lineage(vspID subPop = vspID());

	/** Fill the genotype of all individuals in a population (if
	 *  <tt>subPop=[]</tt>) or in a (virtual) subpopulation \e subPop (if
	 *  <tt>subPop=sp</tt> or <tt>(sp, vsp)</tt>) using a list of alleles
	 *  \e geno. \e geno will be reused if its length is less than
	 *  <tt>subPopSize(subPop)*totNumLoci()*ploidy()</tt>.
	 *  <group>5-genotype</group>
	 */
	void setGenotype(const uintList & geno, vspID subPop = vspID());


	/** Fill the lineage of all individuals in a population (if
	 *  <tt>subPop=[]</tt>) or in a (virtual) subpopulation \e subPop (if
	 *  <tt>subPop=sp</tt> or <tt>(sp, vsp)</tt>) using a list of IDs
	 *  \e lineage. \e lineage will be reused if its length is less than
	 *  <tt>subPopSize(subPop)*totNumLoci()*ploidy()</tt>. This function
	 *  returns directly for modules without lineage information.
	 *  <group>5-genotype</group>
	 */
	void setLineage(const uintList & geno, vspID subPop = vspID());

	//@}

	/** @name utility functions, set subpopulation, save and load etc.
	 */
	//@{

	/** Sort individuals according to values at specified information
	 *  fields (\e infoFields). Individuals will be sorted at an increasing
	 *  order unless \e reverse is set to \c true.
	 */
	void sortIndividuals(const stringList & infoFields, bool reverse=false);

	/** Rearrange individuals to their new subpopulations according to their
	 *  integer values at information field \e field (value returned by
	 *  <tt>Individual::info(field)</tt>). individuals with negative values
	 *  at this \e field will be removed. Existing subpopulation names are
	 *  kept. New subpopulations will have empty names.
	 *  <group>7-manipulate</group>
	 */
	void setSubPopByIndInfo(const string & field);

	/** Split subpopulation \e subPop into subpopulations of given \e sizes,
	 *  which should add up to the size of subpopulation \e subPop or \e 1,
	 *  in which case \e sizes are treated as proportions. If \e subPop
	 *  is not the last subpopulation, indexes of subpopulations after
	 *  \e subPop are shifted. If \e subPop is named, the same name will be
	 *  given to all new subpopulations unless a new set of \e names are
	 *  specified for these subpopulations. This function returns the IDs of
	 *  split subpopulations.
	 *  <group>7-manipulate</group>
	 */
	vectoru splitSubPop(size_t subPop, const vectorf & sizes, const vectorstr & names = vectorstr());


	/** Remove (virtual) subpopulation(s) \e subPops and all their individuals.
	 *  This function can be used to remove complete subpopulations (with
	 *  shifted subpopulation indexes) or individuals belonging to virtual
	 *  subpopulations of a subpopulation. In the latter case, the
	 *  subpopulations are kept even if all individuals have been removed.
	 *  This function only handles the present generation.
	 *  <group>7-manipulate</group>
	 */
	void removeSubPops(const subPopList & subPops);

	/// CPPONLY
	void removeMarkedIndividuals();

	/** remove individual(s) by absolute indexes (parameter \e index) or
	 *  their IDs (parameter \e IDs), or using a filter function (paramter
	 *  \e filter). If indexes are used, only individuals at the current
	 *  generation will be removed. If IDs are used, all individuals with
	 *  one of the IDs at information field \e idField (default to \c "ind_id")
	 *  will be removed. Although \c "ind_id" usually stores unique IDs of
	 *  individuals, this function is frequently used to remove groups of
	 *  individuals with the same value at an information field. An
	 *  \c IndexError will be raised if an index is out of bound, but no
	 *  error will be given if an invalid ID is specified. In the last
	 *  case, a user-defined function should be provided. This function
	 *  should accept parameter \c "ind" or one or more of the information
	 *  fields. All individuals, including ancestors if there are multiple
	 *  ancestral generations, will be passed to this function. Individuals
	 *  that returns \c True will be removed. This function does not affect
	 *  subpopulation structure in the sense that a subpopulation will be
	 *  kept even if all individuals from it are removed.
	 *  <group>7-manipulate</group>
	 */
	void removeIndividuals(const uintList & indexes = vectoru(),
		const floatList & IDs = vectorf(), const string & idField = "ind_id",
		PyObject * filter = NULL);

	/** Merge subpopulations \e subPops. If \e subPops is \c ALL_AVAIL (default),
	 *  all subpopulations will be merged. \e subPops do not have to be adjacent
	 *  to each other. They will all be merged to the subpopulation with the
	 *  smallest subpopulation ID, unless a subpopulation ID is specified using
	 *  parameter \c toSubPop. Indexes of the rest of the subpopulation may
	 *  be changed. A new name can be assigned to the merged subpopulation
	 *  through parameter \e name (an empty \e name will be ignored). This
	 *  function returns the ID of the merged subpopulation.
	 *  <group>7-manipulate</group>
	 */
	size_t mergeSubPops(const uintList & subPops = uintList(), const string & name = UnnamedSubPop, int toSubPop=-1);

	/** Add all individuals, including ancestors, in \e pop to the current
	 *  population. Two populations should have the same genotypic structures
	 *  and number of ancestral generations. Subpopulations in population
	 *  \e pop are kept.
	 *  <group>7-manipulate</group>
	 */
	void addIndFrom(const Population & pop);

	/** Add chromosomes in population \e pop to the current population.
	 *  population \e pop should have the same number of individuals as the
	 *  current population in the current and all ancestral generations.
	 *  Chromosomes of \e pop, if named, should not conflict with names of
	 *  existing chromosome. This function merges genotypes on the
	 *  new chromosomes from population \c pop individual by individual.
	 *  <group>7-manipulate</group>
	 */
	void addChromFrom(const Population & pop);

	/** Add loci from population \e pop. By default, chromosomes are merged
	 *  by index and names of merged chromosomes of population \e pop will
	 *  be ignored (merge of two chromosomes with different names will yield
	 *  a warning). If \e byName is set to \c True, chromosomes in \e pop
	 *  will be merged to chromosomes with identical names. Added loci
	 *  will be inserted according to their position. Their position and
	 *  names should not overlap with any locus in the current population.
	 *  population \e pop should have the same number of individuals as the
	 *  current population in the current and all ancestral generations.
	   Allele lineages are also copied from \e pop in modules with
	 *  lineage information.
	 *  <group>7-manipulate</group>
	 */
	void addLociFrom(const Population & pop, bool byName = false);

	/** Add chromosome \e chromName with given type \e chromType to a
	 *  population, with loci \e lociNames inserted at position \e lociPos.
	 *  \e lociPos should be ordered. \e lociNames and \e chromName should not
	 *  exist in the current population. Allele names could be specified for
	 *  all loci (a list of names) or differently for each locus (a nested
	 *  list of names), using parameter \e alleleNames. Empty loci names will
	 *  be used if \e lociNames is not specified. The newly added alleles
	 *  will have zero lineage in modules wiht lineage information.
	 *  <group>7-manipulate</group>
	 */
	void addChrom(const floatList & lociPos, const stringList & lociNames = vectorstr(),
		const string & chromName = string(), const stringMatrix & alleleNames = stringMatrix(),
		size_t chromType = AUTOSOME);

	/** Insert loci \e lociNames at positions \e pos on chromosome \e chrom.
	 *  These parameters should be lists of the same length, although
	 *  \e names may be ignored, in which case empty strings will be assumed.
	 *  Single-value input is allowed for parameter \e chrom and \e pos if only
	 *  one locus is added. Alleles at inserted loci are initialized with zero
	 *  alleles. Note that loci have to be added to existing chromosomes. If
	 *  loci on a new chromosome need to be added, function <tt>addChrom</tt>
	 *  should be used. Optionally, allele names could be specified either
	 *  for all loci (a single list) or each loci (a nested list). This
	 *  function returns indexes of the inserted loci. Newly inserted
	 *  alleles will have zero lineage in modules with lineage information.
	 *  <group>7-manipulate</group>
	 */
	vectoru addLoci(const uintList & chrom, const floatList & pos,
		const stringList & lociNames = vectorstr(), const stringMatrix & alleleNames = stringMatrix());

	/** Resize population by giving new subpopulation sizes \e sizes.
	 *  individuals at the end of some subpopulations will be removed if the
	 *  new subpopulation size is smaller than the old one. New individuals
	 *  will be appended to a subpopulation if the new size is larger. Their
	 *  genotypes will be set to zero (default), or be copied from existing
	 *  individuals if \e propagate is set to \c True. More specifically,
	 *  if a subpopulation with \c 3 individuals is expanded to \c 7, the
	 *  added individuals will copy genotypes from individual \c 1, \c 2,
	 *  \c 3, and \c 1 respectively. Note that this function only resizes
	 *  the current generation.
	 *  <group>7-manipulate</group>
	 */
	void resize(const uintList & sizes, bool propagate = false);


	/** Extract a list of (virtual) subpopulations from a population and create
	 *  a new population. If \e rearrange is \c False (default), structure and
	 *  names of extracted subpopulations are kept although extracted
	 *  subpopulations can have fewer individuals if they are created from
	 *  extracted virtual subpopulations. (e.g. it is possible to extract all
	 *  male individuals from a subpopulation using a \c SexSplitter()). If
	 *  \e rearrange is \c True, each (virtual) subpopulation in \e subPops
	 *  becomes a new subpopulation in the extracted population in the order
	 *  at which they are specified. Because each virtual subpopulation becomes
	 *  a subpopulation, this function could be used, for example, to separate
	 *  male and female individuals to two subpopulations (
	 *  <tt>subPops=[(0,0), (0,1)]</tt>). If overlapping (virtual)
	 *  subpopulations are specified, individuals will be copied multiple times.
	 *  This function only extract individuals from the present generation.
	 *  <group>7-manipulate</group>
	 */
	Population & extractSubPops(const subPopList & subPops = subPopList(), bool rearrange = false) const;


	/// CPPONLY
	Population & extractMarkedIndividuals() const;

	/** Extract individuals with given absolute indexes (parameter \e indexes),
	 *  IDs (parameter \e IDs, stored in information field \e idField,
	 *  default to \c ind_id), or a filter function (parameter \e filter). If a
	 *  list of absolute indexes are specified, the present generation will be
	 *  extracted and form a one-generational population. If a list of IDs are
	 *  specified, this function will look through all ancestral generations
	 *  and extract individuals with given ID. Individuals with shared IDs are
	 *  allowed. In the last case, a user-defined Python function should be
	 *  provided. This function should accept parameter \c "ind" or one or more
	 *  of the information fields. All individuals, including ancestors if
	 *  there are multiple ancestral generations, will be passed to this
	 *  function. Individuals that returns \c True will be extracted. Extracted
	 *  individuals will be in their original ancestral generations and
	 *  subpopulations, even if some subpopulations or generations are empty.
	 *  An \c IndexError will be raised if an index is out of bound but no
	 *  error will be given if an invalid ID is encountered.
	 *  <group>7-manipulate</group>
	 */
	Population & extractIndividuals(const uintList & indexes = vectoru(),
		const floatList & IDs = vectorf(), const string & idField = "ind_id",
		PyObject * filter = NULL) const;

	/** Extract subsets of individuals, loci and/or information fields from the
	 *  current population and create a new population. By default, all
	 *  genotypes and information fields for all individuals in all ancestral
	 *  generations are extracted. If a list of (virtual) subpopulations are
	 *  given, only individuals in these subpopulations are extracted.
	 *  Structure and names of extracted subpopulations will be kept although
	 *  extracted subpopulations can have fewer individuals if they are created
	 *  from extracted virtual subpopulations. (e.g. it is possible to extract
	 *  all male individuals from a subpopulation using a \c SexSplitter()).
	 *  If a list of loci is specified, only genotypes at specified loci are
	 *  extracted. If a list of \e infoFields is specified, only these
	 *  information fields are extracted. If \e ancGens is not \c ALL_AVAIL (default,
	 *  meaing all ancestral generations), only specified ancestral generations
	 *  will be extracted.
	 *  CPPONLY
	 *  <group>7-manipulate</group>
	 */
	Population & extract(const lociList & extractedLoci, const stringList & infoFieldList,
		const subPopList & subPops = subPopList(), const uintList & ancGens = uintList()) const;

	/** Remove \e loci (absolute indexes or names) and genotypes at these loci
	 *  from the current population. Alternatively, a parameter \e keep can be
	 *  used to specify loci that will not be removed.
	 *  <group>7-manipulate</group>
	 */
	void removeLoci(const lociList & loci = lociList(NULL), const lociList & keep = lociList(NULL));

	/** Recode alleles at \e loci (can be a list of loci indexes or names, or
	 *  all loci in a population (\c ALL_AVAIL)) to other values according to
	 *  parameter \e alleles. This parameter can a list of new allele numbers
	 *  for alleles \c 0, \c 1, \c 2, ... (allele \c x will be recoded to
	 *  <tt>newAlleles[x]</tt>, \c x outside of the range of \e newAlleles will
	 *  not be recoded, although a warning will be given if \c DBG_WARNING is
	 *  defined) or a Python function, which should accept one or both
	 *  parameters \c allele (existing allele) and \c locus (index of locus).
	 *  The return value will become the new allele. This function is intended
	 *  to recode some alleles without listing all alleles in a list. It will
	 *  be called once for each existing allele so it is not possible to recode
	 *  an allele to different alleles. A new list of allele names could be
	 *  specified for these \e loci. Different sets of names could be specified
	 *  for each locus if a nested list of names are given. This function recode
	 *  alleles for all subpopulations in all ancestral generations.
	 *  <group>7-manipulate</group>
	 */
	void recodeAlleles(const uintListFunc & alleles, const lociList & loci = lociList(),
		const stringMatrix & alleleNames = stringMatrix());

	/** Push population \e pop into the current population. Both populations
	 *  should have the same genotypic structure. The current population is
	 *  discarded if \e ancestralDepth (maximum number of ancestral generations
	 *  to hold) is zero so no ancestral generation can be kept. Otherise, the
	 *  current population will become the parental generation of \e pop.
	 *  If \e ancGen of a population is positive and there are already
	 *  \e ancGen ancestral generations (c.f. <tt>ancestralGens()</tt>), the
	 *  greatest ancestral generation will be discarded. In any case, Population
	 *  \e pop becomes invalid as all its individuals are absorbed by the
	 *  current population.
	 *  <group>6-ancestral</group>
	 */
	void push(Population & pop);

	/** HIDDEN
	 *  Return the current ancestral generation number.
	 *  <group>6-ancestral</group>
	 */
	size_t curAncestralGen() const
	{
		return m_curAncestralGen;
	}


	/** Return the actual number of ancestral generations stored in a
	 *  population, which does not necessarily equal to the number set by
	 *  \c setAncestralDepth().
	 *  <group>6-ancestral</group>
	 */
	int ancestralGens() const
	{
		return static_cast<int>(m_ancestralPops.size());
	}


	/** CPPONLY
	 *  clear all information field.
	 */
	void clearInfo()
	{
		// fill_n is faster than fill
		std::fill_n(m_info.begin(), m_info.size(), 0.0);
		// I am glad that memset is not faster than fill_n.
		//memset(&*m_info.begin(), 0, sizeof(double) * m_info.size());
	}


	/** CPPONLY
	 *  Mark subpopulation
	 */
	void markIndividuals(vspID subPop, bool mark) const;

	/** Set information field \c field (specified by index or name) of
	 *  all individuals (if <tt>subPop=[]</tt>, default), or individuals in
	 *  a (virtual) subpopulation (<tt>subPop=sp</tt> or <tt>(sp, vsp)</tt>)
	 *  to \e values. \e values will be reused if its length is smaller than
	 *  the size of the population or (virtual) subpopulation.
	 *  <group>8-info</group>
	 */
	void setIndInfo(const floatList & values, const uintString & field,
		vspID subPop = vspID());


	/// CPPONLY info iterator
	IndInfoIterator infoBegin(size_t idx)
	{
		CHECKRANGEINFO(idx);

		// if there is virtual subpop, use Individual based iterator
		// or
		// if requires order, but the information is not ordered
		// use Individual based
		if (hasActivatedVirtualSubPop() || !indOrdered())
			return IndInfoIterator(idx, indIterator());
		else
			// if not required order, or if the information is ordered
			return IndInfoIterator(idx, m_info.begin(), infoSize());
	}


	/// CPPONLY
	IndInfoIterator infoEnd(size_t idx)
	{
		CHECKRANGEINFO(idx);
		if (hasActivatedVirtualSubPop() || !indOrdered())
			return IndInfoIterator(idx, IndIterator(m_inds.end(), m_inds.end(), false));
		else
			return IndInfoIterator(idx, m_info.end(), infoSize());
	}


	/// CPPONLY info iterator
	IndInfoIterator infoBegin(size_t idx, vspID vsp)
	{
		size_t subPop = vsp.subPop();

		CHECKRANGEINFO(idx);
		CHECKRANGESUBPOP(subPop);
		CHECKRANGEVIRTUALSUBPOP(vsp.virtualSubPop());

		DBG_FAILIF(!vsp.isVirtual() && hasActivatedVirtualSubPop(subPop), RuntimeError,
			"vsp is not virtual but infoBegin iterates through a virtual subpopulations.");

		DBG_FAILIF(vsp.isVirtual() && !hasActivatedVirtualSubPop(subPop), RuntimeError,
			"Virtual subpopulation is not activated.");
		//
		if (vsp.isVirtual() || !indOrdered())
			return IndInfoIterator(idx, indIterator(subPop));
		else
			return IndInfoIterator(idx, m_info.begin() + m_subPopIndex[subPop] * infoSize(), infoSize());
	}


	/// CPPONLY
	IndInfoIterator infoEnd(size_t idx, vspID vsp)
	{
		size_t subPop = vsp.subPop();

		CHECKRANGESUBPOP(subPop);
		CHECKRANGEVIRTUALSUBPOP(vsp.virtualSubPop());

		// has to adjust order because of parameter subPop
		if (vsp.isVirtual() || hasActivatedVirtualSubPop(subPop) || !indOrdered())
			return IndInfoIterator(idx, IndIterator(rawIndEnd(subPop), rawIndEnd(subPop), true));
		else
			return IndInfoIterator(idx, m_info.begin() + m_subPopIndex[subPop + 1] * infoSize(), infoSize());
	}


	/** Return the values (as a list) of information field \c field (by index
	 *  or name) of all individuals (if <tt>subPop=[]</tt>, default), or
	 *  individuals in a (virtual) subpopulation (if <tt>subPop=sp</tt> or
	 *  <tt>(sp, vsp)</tt>).
	 *  <group>8-info</group>
	 */
	vectorf indInfo(const uintString & field, vspID subPop = vspID());


	/** Add a list of information fields \e fields to a population and
	 *  initialize their values to \e init. If an information field alreay
	 *  exists, it will be re-initialized.
	 * <group>8-info</group>
	 */
	void addInfoFields(const stringList & fields, double init = 0);

	/** Set information fields \e fields to a population and initialize them
	 *  with value \e init. All existing information fields will be removed.
	 *  <group>8-info</group>
	 */
	void setInfoFields(const stringList & fields, double init = 0);

	/** Remove information fields \e fields from a population.
	 *  <group>8-info</group>
	 */
	void removeInfoFields(const stringList & fields);


	/** Update information fields \e fields from \e fromFields of another
	 *  population (or Pedigree) \e pop. Two populations should have the same
	 *  number of individuals. If \e fromFields is not specified, it is assumed
	 *  to be the same as \e fields. If \e ancGens is not \c ALL_AVAIL, only
	 *  the specified ancestral generations are updated.
	 *  <group>8-info</group>
	 */
	void updateInfoFieldsFrom(const stringList & fields, const Population & pop,
		const stringList & fromFields = vectorstr(),
		const uintList & ancGens = uintList());

	/** set the intended ancestral depth of a population to \e depth, which can
	 *  be \c 0 (does not store any ancestral generation), \c -1 (store all
	 *  ancestral generations), and a positive number (store \e depth ancestral
	 *  generations. If there exists more than \e depth ancestral generations
	 *  (if \e depth > 0), extra ancestral generations are removed.
	 *  <group>6-ancestral</group>
	 */
	void setAncestralDepth(int depth);

	/// CPPONLY remove certain ancestral generations
	void keepAncestralGens(const uintList & ancGens);

	/** Making ancestral generation \e idx (\c 0 for current generation, \c 1
	 *  for parental generation, \c 2 for grand-parental generation, etc) the
	 *  current generation. This is an efficient way to access Population
	 *  properties of an ancestral generation. <tt>useAncestralGen(0)</tt>
	 *  should always be called afterward to restore the correct order of
	 *  ancestral generations.
	 *  <group>6-ancestral</group>
	 */
	void useAncestralGen(ssize_t idx);

	//@}

	/// CPPONLY
	/**
	   some iterators requires that genotype information is within
	   each subpopulation. We need to adjust genotypic info to
	   obey this.
	   This function is const because the population is 'not changed'
	   conceptually.
	 */
	void syncIndPointers(bool infoOnly = false) const;

	/** Save population to a file \e filename, which can be loaded by a global
	 *  function <tt>loadPopulation(filename)</tt>.
	 *  <group>8-pop</group>
	 */
	void save(const string & filename) const;

	/** CPPONLY load Population from file \e filename
	 *  <group>8-pop</group>
	 */
	void load(const string & filename);

public:
	/** return variables of a population as a Python dictionary. If a valid
	 *  subpopulation \e subPop is specified, a dictionary
	 *  <tt>vars()["subPop"][subPop]</tt> is returned. A \c ValueError will be
	 *  raised if key \e subPop does not exist in \c vars(), or if key
	 *  \e subPop does not exist in <tt>vars()["subPop"]</tt>.
	 *  <group>9-var</group>
	 */
	PyObject * vars(vspID subPop = vspID());


	/// CPPONLY The same as vars(), but without increasing reference count.
	PyObject * dict(vspID subPop = vspID());

	/// CPPONLY
	SharedVariables & getVars() const
	{
		return m_vars;
	}


	/// CPPONLY
	void setDict(PyObject * dict)
	{
		DBG_ASSERT(dict != NULL, SystemError, "Dictionary is empty");
		m_vars.setDict(dict);
	}


	/// CPPONLY
	string varsAsString(bool use_pickle=false) const
	{
		if (use_pickle)
			return m_vars.to_pickle();
		else
			return m_vars.asString();
	}


	/// CPPONLY
	void varsFromString(const string & vars, bool use_pickle=false)
	{
		if (use_pickle)
			return m_vars.from_pickle(vars);
		else
			return m_vars.fromString(vars);
	}


	/** HIDDEN
	 *  evaluate a Python statment/expression in the population's local namespace
	 *  This function evaluates a Python statment( \c stmts )/expression( \c expr )
	 *  and return its result as a string. Optionally run statement( \c stmts ) first.
	 *  <group>8-var</group>
	 */
	PyObject * evaluate(const string & expr = string(), const string & stmts = string())
	{
		return Expression(expr, stmts, m_vars.dict()).evaluate();
	}


	/** HIDDEN
	 *  execute a statement (can be a multi-line string) in the population's local namespace
	 */
	void execute(const string & stmts = string())
	{
		Expression("", stmts, m_vars.dict()).evaluate();
	}


private:
	friend class boost::serialization::access;

	void save(boost::archive::text_oarchive & ar, const unsigned int /* version */) const;

	void load(boost::archive::text_iarchive & ar, const unsigned int /* version */);

	BOOST_SERIALIZATION_SPLIT_MEMBER();

private:
	/// population size: number of individual
	size_t m_popSize;

	/// size of each subpopulation
	vectoru m_subPopSize;

	/// names of each subpopulation
	vectorstr m_subPopNames;

	/// index to subPop \todo change to vectori
	vectoru m_subPopIndex;

	///
	BaseVspSplitter * m_vspSplitter;

	/// pool of genotypic information
#ifdef MUTANTALLELE
	vectorm m_genotype;
#else
	vectora m_genotype;
#endif

#ifdef LINEAGE
	vectori m_lineage;
#endif

	/// information
	/// only in head node
	vectorf m_info;

	/// individuals.
	/// only in head node?
	vector<Individual> m_inds;

	int m_ancestralGens;

	/// shared variables for this population
	mutable SharedVariables m_vars;

	/// store previous populations
	/// need to store: subPopSize, genotype and m_inds
	struct popData
	{
		vectoru m_subPopSize;
		vectorstr m_subPopNames;
#ifdef MUTANTALLELE
		vectorm m_genotype;
#else
		vectora m_genotype;
#endif

#ifdef LINEAGE
		vectori m_lineage;
#endif

		vectorf m_info;
		vector<Individual> m_inds;
		bool m_indOrdered;

		// swap between a popData and existing data.
		void swap(Population & pop);

	};

	std::deque<popData> m_ancestralPops;

	/// current ancestral depth
	int m_curAncestralGen;

	/// whether or not individual genotype and information are in order
	/// within a population.
	mutable bool m_indOrdered;

	mutable size_t m_gen;
	mutable size_t m_rep;

public:
	/** CPPONLY
	 *  current replicate in a simulator which is not meaningful for a stand-alone population
	 *	<group>evolve</group>
	 */
	size_t rep()
	{
		return m_rep;
	}


	/// CPPONLY  set rep number
	void setRep(size_t rep)
	{
		m_rep = rep;
		m_vars.setVar("rep", rep);
	}


	/** CPPONLY
	 *  Return the current generation number of this population. This is
	 *  faster than getting generation number from population variable.
	 *  <group>evolve</group>
	 */
	size_t gen() const
	{
		return m_gen;
	}


	/// CPPONLY
	void setGen(size_t gen)
	{
		m_gen = gen;
		m_vars.setVar("gen", gen);
	}


};

/** load a population from a file saved by <tt>Population::save()</tt>.
 */
Population & loadPopulation(const string & file);

}


#ifndef SWIG
#  ifndef _NO_SERIALIZATION_
// version 0: base (reset for version 1.0)
// version 1: with lineage information for lineage-aware modules
// version 2: for memory-efficient save/load
// version 3: use pickle to save load population variables
BOOST_CLASS_VERSION(simuPOP::Population, 3)
#  endif
#endif
#endif
