/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu
*                                                                         *
*   $LastChangedDate$
*   $Rev$                                                     *
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*   This program is distributed in the hope that it will be useful,       *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU General Public License     *
*   along with this program; if not, write to the                         *
*   Free Software Foundation, Inc.,                                       *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/
#ifndef _POPULATION_H
#define _POPULATION_H

/**
   \file
   \brief head file of class population
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

#include <fstream>
using std::ifstream;
using std::ofstream;

// used to save history population
// 0 (first parental) 1, ...., n
#include <deque>
using std::deque;

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/version.hpp>

#include "individual.h"
#include "virtualSubPop.h"


namespace simuPOP {

class pedigree;

/**
 *  A simuPOP population consists of individuals of the same genotypic
 *  structure, organized by generations, subpopulations and virtual
 *  subpopulations. It also contains a Python dictionary that is used to
 *  store arbitrary population variables.
 *
 *  In addition to genotypic structured related functions provided by the
 *  \c genoStruTrait class, the population class provides a large number
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
 *    sequentially individual by individual. Please refer to class \c individual
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
class population : public GenoStruTrait
{
public:
#define Haplodiploid 2.5

	/** @name  constructors and destructor */
	//@{

	/** The following parameters are used to create a population object:
	 *
	 *  \param size A list of subpopulation sizes. The length of this list
	 *    determines the number of subpopulations of this population. If
	 *    there is no subpopulation, <em>size</em><tt>=[popSize]</tt> can be
	 *    written as <em>size</em><tt>=popSize</tt>.
	 *  \param ploidy Number of homologous sets of chromosomes. Default to
	 *    \c 2 (diploid). For efficiency considerations, all chromosomes have
	 *    the same number of homologous sets, even if some chromosomes (e.g.
	 *    mitochondrial) or some individuals (e.g. males in a haplodiploid
	 *    population) have different numbers of homologous sets. The first
	 *    case is handled by setting \e chromTypes of each chromosome. Only
	 *    the haplodiploid populations are handled for the second case, for
	 *    which <tt>ploidy=Haplodiploid</tt> should be used.
	 *  \param loci A list of numbers of loci on each chromosome. The length of
	 *    this parameter determines the number of chromosomes. Default to
	 *    <tt>[1]</tt>, meaning one chromosome with a single locus.
	 *  \param chromTypes A list that specifies the type of each chromosome,
	 *    which can be \c Autosome, \c ChromosomeX, \c ChromosomeY, or
	 *    \c Mitochondrial. All chromosomes are assumed to be autosomes if
	 *    this parameter is ignored. Sex chromosome can only be specified in a
	 *    diploid population where the sex of an individual is determined by
	 *    the existence of these chromosomes using the \c XX (\c Female) and
	 *    \c XY (\c Male) convention. Both sex chromosomes have to be available
	 *    and be specified only once. Because chromosomes \c X and \c Y are
	 *    treated as two chromosomes, recombination on the pseudo-autosomal
	 *    regions of the sex chromsomes is not supported. \c Mitochondrial
	 *    chromosomes are inherited maternally. If more than one mitochondrial
	 *    chromosomes are available, they are selected randomly with
	 *    replacement. No recombination is allowed for mitochondrial
	 *    chromosomes.
	 *  \param lociPos Positions of all loci on all chromosome, as a list of
	 *    float numbers. Default to \c 1, \c 2, ... etc on each chromosome.
	 *    Positions on the same chromosome should be ordered. A nested list
	 *    that specifies positions of loci on each chromosome is also
	 *    acceptable.
	 *  \param ancGen Number of the most recent ancestral generations to keep
	 *    during evolution. Default to \c 0, which means only the current
	 *    generation will be kept. If it is set to \c -1, all ancestral
	 *    generations will be kept in this population (and exhaust your computer
	 *    RAM quickly).
	 *  \param chromNames A list of chromosome names. Default to \c chrom1,
	 *    \c chrom2, ... etc.
	 *  \param alleleNames A list of allele names for all markers. For example,
	 *    <em>alleleNames</em><tt>=('A','C','T','G')</tt> names allele \c 0 --
	 *    \c 3 \c 'A', \c 'C', \c 'T', and \c 'G' respectively. Note that
	 *    simuPOP does not yet support locus-specific allele names.
	 *  \param lociNames A list or a matrix (separated by chromosomes) of names
	 *    for each locus. Default to \c "locX-Y" where \c X and \c Y are 1-based
	 *    chromosome and locus indexes, respectively.
	 *  \param infoFields Names of information fields (named float number) that
	 *    will be attached to each individual.
	 */
	population(const vectorlu & size = vectorlu(),
		float ploidy = 2,
		const vectoru & loci = vectoru(),
		const vectoru & chromTypes = vectoru(),
		const vectorf & lociPos = vectorf(),
		int ancGen = 0,
		const vectorstr & chromNames = vectorstr(),
		const vectorstr & alleleNames = vectorstr(),
		const vectorstr & lociNames = vectorstr(),
		const vectorstr & infoFields = vectorstr());

	/// CPPONLY copy constructor
	population(const population & rhs);

	/** Create a cloned copy of a population. Note that Python statement
	 *  <tt>pop1 = pop</tt> only creates a reference to an existing population
	 *  \c pop.
	 *  <group>1-pop</group>
	 */
	population * clone() const;

	/** HIDDEN (do not see a need to expose this function yet.)
	 *  swap the content of two populations
	 *  <group>1-pop</group>
	 */
	void swap(population & rhs)
	{
		GenoStruTrait::swap(rhs);
		std::swap(m_popSize, rhs.m_popSize);
		m_subPopSize.swap(rhs.m_subPopSize);
		m_subPopIndex.swap(rhs.m_subPopIndex);
		m_genotype.swap(rhs.m_genotype);
		m_info.swap(rhs.m_info);
		m_inds.swap(rhs.m_inds);
		std::swap(m_ancestralGens, rhs.m_ancestralGens);
		m_vars.swap(rhs.m_vars);
		m_ancestralPops.swap(rhs.m_ancestralPops);
		std::swap(m_rep, rhs.m_rep);
		std::swap(m_gen, rhs.m_gen);
		std::swap(m_curAncestralGen, rhs.m_curAncestralGen);
		std::swap(m_indOrdered, rhs.m_indOrdered);
		std::swap(m_vspSplitter, rhs.m_vspSplitter);
	}


	/// destroy a population
	~population();


	/// used by Python print function to print out the general information of the population
	string __repr__()
	{
		return "<simuPOP::population of size " + toStr(popSize()) + ">";
	}


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
	void fitSubPopStru(const vectorlu & newSubPopSizes);

	/** if a population has any activated virtual subpopulations
	 *  CPPONLY
	 */
	bool hasActivatedVirtualSubPop() const;

	/** if a subpopulation has any activated virtual subpopulation
	 *  CPPONLY
	 */
	bool hasActivatedVirtualSubPop(SubPopID subPop) const;

	/** CPPONLY because this is simply numVirtualSubPop() != 0.
	 *  Return True if virtual subpopulations are defined for this population.
	 *  <group>3-VSP</group>
	 */
	bool hasVirtualSubPop() const;

	/// CPPONLY
	vspSplitter * virtualSplitter() const { return m_vspSplitter; }

	/** Set a VSP \e splitter to the population, which defines the same VSPs
	 *  for all subpopulations. If different VSPs are needed for different
	 *  subpopulations, a \c combinedSplitter can be used to make these VSPs
	 *  available to all subpopulations.
	 *  <group>3-VSP</group>
	 */
	void setVirtualSplitter(vspSplitter * splitter);

	/** Return the number of virtual subpopulations (VSP) defined by a VSP
	 *  splitter. Return \c 0 if no VSP is defined.
	 *  <group>3-VSP</group>
	 */
	UINT numVirtualSubPop() const;

	/// HIDDEN activate a virtual subpopulation.
	/**
	   \param id subpopulation id
	   \param vid virtual subpopulation id
	 */
	void activateVirtualSubPop(SubPopID subPop, SubPopID virtualSubPop = InvalidSubPopID,
		IterationType type = VisibleInds);

	/** HIDDEN
	 *  deactivate virtual subpopulations in a given
	 *  subpopulation. In another word, all individuals
	 *  will become visible.
	 */
	void deactivateVirtualSubPop(SubPopID subPop);

	// allow compaison of populations in python
	// only equal or unequal, no greater or less than
	/// a python function used to compare the population objects
	int __cmp__(const population & rhs) const;

	/** HIDDEN
	 *  set population/subpopulation structure given subpopulation sizes
	 *  \param newSubPopSizes an array of new subpopulation sizes. The overall
	 *    population size should not changed.
	 *  <group>2-subpop</group>
	 */
	void setSubPopStru(const vectorlu & newSubPopSizes);

	/** Return the number of subpopulations in a population. Return 1 if there
	 *  is no subpopulation structure.
	 *  <group>2-subpop</group>
	 */
	UINT numSubPop() const
	{
		return m_subPopSize.size();
	}


	/** Return the size of a subpopulation (<tt>subPopSize(sp)</tt>) or a
	 *  virtual subpopulation (<tt>subPopSize([sp, vsp])<tt>).
	 *  <group>2-subpop</group>
	 */
	ULONG subPopSize(vspID subPop) const
	{
		CHECKRANGESUBPOP(subPop.subPop());
		CHECKRANGEVIRTUALSUBPOP(subPop.virtualSubPop());
		if (hasActivatedVirtualSubPop() || subPop.isVirtual())
			return m_vspSplitter->size(*this, subPop.subPop(), subPop.virtualSubPop());
		else
			return m_subPopSize[subPop.subPop()];
	}


	/** Return the name of a virtual subpopulation \e subPop (specified by a
	 *  <tt>(sp, vsp)</tt> pair). Because VSP names are the same across all
	 *  subpopulations, a single VSP index is also acceptable.
	 *  <group>3-VSP</group>
	 */
	string virtualSubPopName(vspID subPop) const;

	/** Return the sizes of all subpopulations in a list. Virtual
	 *  subpopulations are not considered.
	 *  <group>2-subpop</group>
	 */
	vectorlu subPopSizes() const
	{
		return m_subPopSize;
	}


	//@}

	/** @name indexes, conversion between absoluate indexes and relative indexes. return of chromomosome/subpopulation indexes.
	 */
	//@{

	/** Return the total number of individuals in all subpopulations.
	 *  <group>2-subpop</group>
	 */
	ULONG popSize() const
	{
		return m_popSize;
	}


	/** return the absolute index of an individual \e idx in subpopulation \e subPop.
	 *  <group>2-subpop</group>
	 */
	ULONG absIndIndex(ULONG idx, UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);
		CHECKRANGESUBPOPMEMBER(idx, subPop);

		return m_subPopIndex[subPop] + idx;
	}


	/** Return the subpopulation ID and relative index of an individual, given
	 *  its absolute index \c idx.
	 *  <group>2-subpop</group>
	 */
	std::pair<UINT, ULONG> subPopIndPair(ULONG idx)
	{
		CHECKRANGEIND(idx);

		pair<UINT, ULONG> loc;

		for (UINT i = 1; i <= m_subPopSize.size(); ++i) {
			if (m_subPopIndex[i] > idx) {
				loc.first = i - 1;
				loc.second = idx - m_subPopIndex[i - 1];
				break;
			}
		}
		return loc;
	}


	/** Return the index of the first individual in subpopulation \e subPop. An
	 *  \c IndexError will be raised if \e subPop is out of range.
	 *  <group>2-subpop</group>
	 */
	ULONG subPopBegin(UINT subPop) const
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
	ULONG subPopEnd(UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return m_subPopIndex[subPop + 1];
	}


	//@}
	/** @name itertors and accessers, ways to access information, mainly various iterators.
	 */
	//@{

	/** Return a refernce to individual \e ind in the population.
	 *  <group>4-ind</group>
	 */
	individual & ind(ULONG idx)
	{
		CHECKRANGEIND(idx);

		return m_inds[idx];
	}


	/** Return a refernce to individual \e ind in subpopulation \e subPop.
	 *  <group>4-ind</group>
	 */
	individual & ind(ULONG idx, UINT subPop)
	{
		CHECKRANGESUBPOPMEMBER(idx, subPop);

		return m_inds[subPopBegin(subPop) + idx];
	}

	/// CPPONLY refernce to individual \c ind in subpopulation \c subPop
	/** CPPONLY
	 *  Return a reference to individual \e ind from subpopulation \subPop.
	 *  <group>4-ind</group>
	 */
	const individual & ind(ULONG idx, UINT subPop = 0) const
	{
#ifndef OPTIMIZED
		if (subPop > 0) {
			CHECKRANGESUBPOPMEMBER(idx, subPop);
		} else {
			CHECKRANGEIND(idx);
		}
#endif

		return m_inds[subPopBegin(subPop) + idx];
	}


	/** Return a reference to individual \c idx in ancestral generation \c gen.
	 *  The correct individual will be returned even if the current generation
	 *  is not the present one (see also \c useAncestralGen).
	 *  <group>6-ancestral</group>
	 */
	individual & ancestor(ULONG idx, UINT gen);

	/// refrence to an individual \c ind in an ancestral generation
	const individual & ancestor(ULONG ind, UINT gen) const;

	/** Return a reference to individual \c idx of subpopulation \e subPop in
	 *   ancestral generation \c gen.
	 *  <group>6-ancestral</group>
	 */
	individual & ancestor(ULONG ind, UINT subPop, UINT gen);

	/// refrence to an individual \c ind in a specified subpopulaton or an ancestral generation
	const individual & ancestor(ULONG ind, UINT subPop, UINT gen) const;

	/** Return a Python iterator that can be used to iterate through all
	 *  individuals in a population.
	 *  <group>4-ind</group>
	 */
	pyIndIterator individuals()
	{
		// if a virtual subpopulation is activated, this will
		// iterate through virtual subpopulation. However,
		// users are not supposed to manually activate subpopulation
		// so this feature is CPPONLY
		return pyIndIterator(m_inds.begin(), m_inds.end(),
			!hasActivatedVirtualSubPop(), true);
	}


	/** Return an iterator that can be used to iterate through all individuals
	 *  in a subpopulation (<tt>subPop=spID</tt>) or a virtual subpopulation
	 *  (<tt>subPop=[spID, vspID]</tt>).
	 *  <group>4-ind</group>
	 */
	pyIndIterator individuals(vspID subPop)
	{
		SubPopID spID = subPop.subPop();
		SubPopID vspID = subPop.virtualSubPop();

#ifndef OPTIMIZED
		CHECKRANGESUBPOP(spID);
		CHECKRANGEVIRTUALSUBPOP(vspID);
		DBG_FAILIF(hasActivatedVirtualSubPop(spID), ValueError,
			"This operation is not allowed for an activated subpopulation");
#endif
		if (subPop.isVirtual()) {
			// this does not need to be deactivated...
			activateVirtualSubPop(spID, vspID, IteratableInds);
			// if there is no virtual subpop
			return pyIndIterator(m_inds.begin() + subPopBegin(spID),
				m_inds.begin() + subPopEnd(spID),
				// allInds will not work at all, because there will be
				// virtual subpopulation
				false,
				// and we count visible, and iteratable individuals.
				false);
		} else
			return pyIndIterator(m_inds.begin() + subPopBegin(spID),
				m_inds.begin() + subPopEnd(spID),
				// if there is no activated virtual subpopualtions
				// iterate through all individuals.
				!hasActivatedVirtualSubPop(spID),
				// otherwise, iterate through all visible individuals.
				true);
	}


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


	/// CPPONLY individual iterator: without subPop info
	IndIterator indBegin(IterationType type = VisibleInds)
	{
		return IndIterator(m_inds.begin(), m_inds.end(),
			(type == VisibleInds && !hasActivatedVirtualSubPop()) ? AllInds : type);
	}


	/// CPPONLY individual iterator: without subPop info
	IndIterator indEnd(IterationType type = VisibleInds)
	{
		return IndIterator(m_inds.end(), m_inds.end(),
			(type == VisibleInds && !hasActivatedVirtualSubPop()) ? AllInds : type);
	}


	/** CPPONLY individual iterator: with subPop info.
	 *  The iterator will skip invisible individuals
	 */
	IndIterator indBegin(UINT subPop, IterationType type = VisibleInds)
	{
		CHECKRANGESUBPOP(subPop);
		DBG_FAILIF(hasActivatedVirtualSubPop(subPop) && type == VisibleInds,
			ValueError, "Can not iterate through an VSP with iteratable iterator");

		return IndIterator(m_inds.begin() + m_subPopIndex[subPop],
			m_inds.begin() + m_subPopIndex[subPop + 1],
			(type == VisibleInds && !hasActivatedVirtualSubPop(subPop)) ? AllInds : type);
	}


	/// CPPONLY individual iterator: with subPop info.
	IndIterator indEnd(UINT subPop, IterationType type = VisibleInds)
	{
		CHECKRANGESUBPOP(subPop);

		return IndIterator(m_inds.begin() + m_subPopIndex[subPop + 1],
			m_inds.begin() + m_subPopIndex[subPop + 1],
			(type == VisibleInds && !hasActivatedVirtualSubPop(subPop)) ? AllInds : type);
	}


	/** CPPONLY individual iterator: without subPop info
	 *  The iterator will skip invisible individuals
	 */
	ConstIndIterator indBegin(IterationType type = VisibleInds) const
	{
		return ConstIndIterator(m_inds.begin(), m_inds.end(),
			(type == VisibleInds && !hasActivatedVirtualSubPop()) ? AllInds : type);
	}


	/** CPPONLY individual iterator: without subPop info
	 *  It is recommended to use it.valid(), instead of it != indEnd()
	 */
	ConstIndIterator indEnd(IterationType type = VisibleInds) const
	{
		return ConstIndIterator(m_inds.end(), m_inds.end(),
			(type == VisibleInds && !hasActivatedVirtualSubPop()) ? AllInds : type);
	}


	/** CPPONLY individual iterator: with subPop info.
	 *  The iterator will skip invisible individuals
	 */
	ConstIndIterator indBegin(UINT subPop, IterationType type) const
	{
		CHECKRANGESUBPOP(subPop);

		return ConstIndIterator(m_inds.begin() + m_subPopIndex[subPop],
			m_inds.begin() + m_subPopIndex[subPop + 1],
			(type == VisibleInds && !hasActivatedVirtualSubPop(subPop)) ? AllInds : type);
	}


	/** CPPONLY individual iterator: with subPop info.
	 * It is recommended to use it.valid(), instead of it != indEnd(sp)
	 */
	ConstIndIterator indEnd(UINT subPop, IterationType type) const
	{
		CHECKRANGESUBPOP(subPop);

		return ConstIndIterator(m_inds.begin() + m_subPopIndex[subPop],
			m_inds.begin() + m_subPopIndex[subPop + 1],
			(type == VisibleInds && !hasActivatedVirtualSubPop(subPop)) ? AllInds : type);
	}


	/** CPPONLY individual iterator: without subPop info
	 */
	RawIndIterator rawIndBegin()
	{
		return m_inds.begin();
	}


	/** CPPONLY individual iterator: without subPop info
	 */
	RawIndIterator rawIndEnd()
	{
		return m_inds.end();
	}


	/** CPPONLY individual iterator: with subPop info.
	 * The iterator will skip invisible individuals
	 */
	RawIndIterator rawIndBegin(UINT subPop)
	{
		CHECKRANGESUBPOP(subPop);

		return m_inds.begin() + m_subPopIndex[subPop];
	}


	/** CPPONLY individual iterator: with subPop info.
	 */
	RawIndIterator rawIndEnd(UINT subPop)
	{
		CHECKRANGESUBPOP(subPop);

		return m_inds.begin() + m_subPopIndex[subPop + 1];
	}


	/** CPPONLY individual iterator: without subPop info
	 * The iterator will skip invisible individuals
	 */
	ConstRawIndIterator rawIndBegin() const
	{
		return m_inds.begin();
	}


	/** CPPONLY individual iterator: without subPop info
	 *  It is recommended to use it.valid(), instead of it != indEnd()
	 */
	ConstRawIndIterator rawIndEnd() const
	{
		return m_inds.end();
	}


	/** CPPONLY individual iterator: with subPop info.
	 * The iterator will skip invisible individuals
	 */
	ConstRawIndIterator rawIndBegin(UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return m_inds.begin() + m_subPopIndex[subPop];
	}


	/** CPPONLY individual iterator: with subPop info.
	 * It is recommended to use it.valid(), instead of it != indEnd(sp)
	 */
	ConstRawIndIterator rawIndEnd(UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return m_inds.begin() + m_subPopIndex[subPop + 1];
	}


	/// CPPONLY allele iterator that access a locus across all copies of chromosomes and individual
	/**
	   \param locus allele access, given locus, return the first allele. ptr++ go the next one.
	   Default return the beginning of the first subpopulation, also the first of the whole population

	   \param order order = True: indiviudls in order
	       order = false: do not even respect subpops

	   \note The order of alleles DOES NOT HAVE TO match the order of individuals. Only the boundary of
	   subpopulations will be respected.	Therefore, it is possible to access all alleles within an
	   subpopulation	through such iterators.
	 */
	IndAlleleIterator alleleBegin(UINT locus)
	{
		CHECKRANGEABSLOCUS(locus);

		// if there is virtual subpop, use individual based iterator
		// or
		// if requires order, but the alleles are not ordered
		// use individual based
		if (hasActivatedVirtualSubPop() || !indOrdered())
			return IndAlleleIterator(locus, indBegin(), ploidy(), totNumLoci());
		else
			return IndAlleleIterator(m_genotype.begin() + locus, totNumLoci());
	}


	/// CPPONLY allele iterator
	IndAlleleIterator alleleEnd(UINT locus)
	{
		CHECKRANGEABSLOCUS(locus);
		if (hasActivatedVirtualSubPop() || !indOrdered())
			return IndAlleleIterator(locus, indEnd(), ploidy(), totNumLoci());
		else
			return IndAlleleIterator(m_genotype.begin() + locus + m_popSize * genoSize(), totNumLoci());
	}


	/// CPPONLY allele begin, for given subPop
	/**
	   order = True: keep order
	   order = false: repect subpop
	 */
	IndAlleleIterator alleleBegin(UINT locus, UINT subPop)
	{
		CHECKRANGEABSLOCUS(locus);
		CHECKRANGESUBPOP(subPop);

		if (hasActivatedVirtualSubPop() || !indOrdered())
			return IndAlleleIterator(locus, indBegin(subPop), ploidy(), totNumLoci());
		else
			return IndAlleleIterator(m_genotype.begin() + m_subPopIndex[subPop] * genoSize() +
				locus, totNumLoci());
	}


	///  CPPONLY allele iterator
	IndAlleleIterator alleleEnd(UINT locus, UINT subPop)
	{
		CHECKRANGEABSLOCUS(locus);
		CHECKRANGESUBPOP(subPop);

		if (hasActivatedVirtualSubPop() || !indOrdered())
			return IndAlleleIterator(locus, indEnd(subPop), ploidy(), totNumLoci());
		else
			return IndAlleleIterator(m_genotype.begin() + m_subPopIndex[subPop + 1] * genoSize() +
				locus, totNumLoci());
	}


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
			sortIndividuals();

		return m_genotype.begin();
	}


	///  CPPONLY allele iterator
	GenoIterator genoEnd(bool order)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This function is not valid with an activated virtual subpopulation");
		if (order && !indOrdered())
			sortIndividuals();

		return m_genotype.end();
	}


	///  CPPONLY allele iterator, go through all allels one by one in a subpopulation
	/**
	   if order, keep order
	   if not order, respect subpopulation structure
	 */
	GenoIterator genoBegin(UINT subPop, bool order)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This function is not valid with an activated virtual subpopulation");
		CHECKRANGESUBPOP(subPop);

		sortIndividuals();

		return m_genotype.begin() + m_subPopIndex[subPop] * genoSize();
	}


	/// CPPONLY allele iterator in a subpopulation.
	GenoIterator genoEnd(UINT subPop, bool order)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This function is not valid with an activated virtual subpopulation");
		CHECKRANGESUBPOP(subPop);
		sortIndividuals(order);

		return m_genotype.begin() + m_subPopIndex[subPop + 1] * genoSize();
	}


	/// CPPONLY genoIterator --- beginning of individual ind.
	GenoIterator indGenoBegin(ULONG ind) const
	{
		CHECKRANGEIND(ind);
		return m_inds[ind].genoBegin();
	}


	/// CPPONLY genoIterator -- end of individual ind.
	GenoIterator indGenoEnd(ULONG ind) const
	{
		CHECKRANGEIND(ind);
		return m_inds[ind].genoEnd();
	}


	/// CPPONLY genoIterator --- beginning of individual ind.
	GenoIterator indGenoBegin(ULONG ind, UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);
		CHECKRANGESUBPOPMEMBER(ind, subPop);

		return m_inds[ subPopBegin(subPop) + ind].genoBegin();
	}


	/// CPPONLY genoIterator -- end of individual ind.
	GenoIterator indGenoEnd(ULONG ind, UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);
		CHECKRANGESUBPOPMEMBER(ind, subPop);

		return m_inds[ subPopBegin(subPop) + ind].genoEnd();
	}


	/// HIDDEN get the whole genotypes
	/**
	   Return an editable array of all genotypes of the population. You need to
	   know how these genotypes are organized to safely read/write genotype
	   directly.
	   \param order if order is \c true, individuals will be ordered such that
	    <tt>pop.individual(x).arrGenotype() == pop.arrGenotype()[x*pop.genoSize():(x+1)*pop.genoSize()]</tt>.
	 */
	PyObject * arrGenotype(bool order);

	/// HIDDEN get the whole genotypes of individuals in a subpopulation
	/**
	   Return an editable array of all genotype in a subpopulation.
	   \param subPop index of subpopulation (start from 0)
	   \param order if order is \c true, individuals will be ordered.
	 */
	PyObject * arrGenotype(UINT subPop, bool order);

	/** Return an editable array of the genotype of all individuals in this
	 *  population.
	 *  <group>5-genotype</group>
	 */
	PyObject * genotype();

	/** Return an editable array of the genotype of all individuals in
	 *  subpopulation \e subPop.
	 *  <group>5-genotype</group>
	 */
	PyObject * genotype(SubPopID subPop);

	/** Fill the genotype of all individuals of a population using a list of
	 *  alleles \e geno. \e geno will be reused if its length is less than
	 *  <tt>popSize()*totNumLoci()*ploidy()</tt>.
	 *  <group>5-genotype</group>
	 */
	void setGenotype(vectora geno);

	/** Fill the genotype of all individuals of in subpopulation \e subPop
	 *  using a list of alleles \e geno. \e geno will be reused if its length
	 *  is less than <tt>subPopSize(subPop)*totNumLoci()*ploidy()</tt>.
	 *  <group>5-genotype</group>
	 */
	void setGenotype(vectora geno, SubPopID subPop);

	//@}

	/** @name utility functions, set subpopulation, save and load etc.
	 */
	//@{

	/** Rearrange individuals to their new subpopulations according to their
	 *  integer values at information field \e field (value returned by
	 *  <tt>individual::indInfo(field)</tt>). Individuals with negative values
	 *  at this \e field will be removed.
	 *  <group>7-manipulate</group>
	 */
	void setSubPopByIndInfo(const string & field);

	/** Split subpopulation \e subPop into subpopulations of given \e sizes,
	 *  which should add up to the size of subpopulation \e subPop.
	 *  Alternatively, \e sizes can be a list of proportions (add up to \c 1)
	 *  from which the sizes of new subpopulations are determined. If \e subPop
	 *  is not the last subpopulation, subpopulation indexes will be changed.
	 *  For example, if you split the second subpopulation into two, a
	 *  population with three subpopulations will have four subpopulations:
	 *  \c 0 (untouched), \c 1.1->1, \c 1.2->2, \c 2->3 (changed).
	 *  <group>7-manipulate</group>
	 */
	void splitSubPop(UINT subPop, vectorf sizes);

	/** remove empty subpopulations by adjusting subpopulation IDs.
	 *  <group>7-manipulate</group>
	 */
	void removeEmptySubPops();

	/** Remove subpopulations \e subPop and all their individuals. Indexes of
	 *  subpopulations after removed subpopulations will be shifted.
	 *  <group>7-manipulate</group>
	 */
	void removeSubPops(const vectoru & subPops);

	/** remove individuals \e inds (absolute indexes) from the current
	 *  population. A subpopulation will be kept even if all individuals from
	 *  it are removed. This function only affects the current generation.
	 *  <group>7-manipulate</group>
	 */
	void removeIndividuals(const vectoru & inds);

	/** Merge subpopulations \e subPops. If \e subPops is empty (default), all
	 *  subpopulations will be merged. \e subPops do not have to be adjacent to
	 *  each other. They will all be merged to the subpopulation with the
	 *  smallest subpopulation ID. Indexes of the rest of the subpopulation may
	 *  be changed.
	 *  <group>7-manipulate</group>
	 */
	void mergeSubPops(const vectoru & subPops = vectoru());

	/** Add all individuals, including ancestors, in \e pop to the current
	 *  population. Two populations should have the same genotypic structures
	 *  and number of ancestral generations. Subpopulations in population
	 *  \e pop are kept.
	 *  <group>7-manipulate</group>
	 */
	void addIndFromPop(const population & pop);

	/** Add chromosomes in population \e pop to the current population.
	 *  Population \e pop should have the same number of individuals as the
	 *  current population in the current and all ancestral generations.
	 *  This function merges genotypes on the
	 *  new chromosomes from population \c pop individual by individual.
	 *  <group>7-manipulate</group>
	 */
	void addChromFromPop(const population & pop);

	/** Add loci from population \e pop, chromosome by chromosome. Added
	 *  loci will be inserted according to their position. Their position
	 *  and names should not overlap with any locus in the current population.
	 *  Population \e pop should have the same number of individuals as the
	 *  current population in the current and all ancestral generations.
	 *  <group>7-manipulate</group>
	 */
	void addLociFromPop(const population & pop);

	/** Add chromosome \e chromName with given type \e chromType to a
	 *  population, with loci \e lociNames inserted at position \e lociPos.
	 *  \e lociPos should be ordered. \e lociNames and \e chromName should not
	 *  exist in the current population. If they are not specified, simuPOP
	 *  will try to assign default names, and raise a \c ValueError if the
	 *  default names have been used.
	 *  <group>7-manipulate</group>
	 */
	void addChrom(const vectorf & lociPos, const vectorstr & lociNames = vectorstr(),
		const string & chromName = "", UINT chromType = Autosome);

	/** Insert loci \e names at positions \e pos on chromosome \e chrom.
	 *  These parameters should be lists of the same length, although
	 *  \e names may be ignored, in which case random names will be given.
	 *  Alleles at inserted loci are initialized with zero alleles. Note that
	 *  loci have to be added to existing chromosomes. If loci on a new
	 *  chromosome need to be added, function <tt>addChrom</tt> should be
	 *  used. This function returns indexes of the inserted loci.
	 *  <group>7-manipulate</group>
	 */
	vectoru addLoci(const vectoru & chrom, const vectorf & pos,
		const vectorstr & names = vectorstr());

	/** Resize population by giving new subpopulation sizes \e newSubPopSizes.
	 *  Individuals at the end of some subpopulations will be removed if the
	 *  new subpopulation size is smaller than the old one. New individuals
	 *  will be appended to a subpopulation if the new size is larger. Their
	 *  genotypes will be set to zero (default), or be copied from existing
	 *  individuals if \e propagate is set to \c True. More specifically,
	 *  if a subpopulation with \c 3 individuals is expanded to \c 7, the
	 *  added individuals will copy genotypes from individual \c 1, \c 2,
	 *  \c 3, and \c 1 respectively. Note that this function only resizes
	 *  the current generation. This function affects only the current
	 *  generation.
	 *  <group>7-manipulate</group>
	 */
	void resize(const vectorlu & newSubPopSizes, bool propagate = false);

	/** Extract subsets of individuals, loci and/or information fields from the
	 *  current population and create a new one. If information field \e field
	 *  is not \c None, individuals with negative values at this information
	 *  field will be removed, and others are put into subpopulations specified
	 *  by this field. If \e loci is not \c None, only genotypes at \e loci
	 *  are extracted. If \e infoFields is not \c None, only these information
	 *  fields will be extracted. If \e ancGen is not \c -1 (default, meaing
	 *  all ancestral generations), only \e ancGen ancestral generations will
	 *  be kept. As an advanced feature, \e field can be information field
	 *  of a pedigree object \e ped. This allows extraction of individuals
	 *  according to pedigrees identified in a pedigree object. This pedigree
	 *  should have the same number of individuals in all generations.
	 *  <group>7-manipulate</group>
	 */
	population & extract(bool removeInd, const string & field,
		bool removeLoci, const vectoru & loci,
		bool removeInfo, const vectorstr & infoFields, int ancGen = -1,
		pedigree * ped = NULL) const;

	/** Remove \e loci (absolute indexes) and genotypes at these loci from the
	 *  current population. Alternatively, a parameter \e keep can be used to
	 *  specify loci that will not be removed.
	 *  <group>7-manipulate</group>
	 */
	void removeLoci(const vectoru & loci = vectoru(), const vectoru & keep = vectoru());

	/** Push population \e pop into the current population. Both populations
	 *  should have the same genotypic structure. The current population is
	 *  discarded if \e ancestralDepth (maximum number of ancestral generations
	 *  to hold) is zero so no ancestral generation can be kept. Otherise, the
	 *  current population will become the parental generation of \e pop,
	 *  advancing the greatness level of all existing ancestral generations by
	 *  one. If \e ancestralDepth is positive and there are already
	 *  \e ancestralDepth ancestral generations (see also:
	 *  <tt>ancestralGens()</tt>), the greatest ancestral generation will be
	 *  discarded. In any case, population \e pop becomes invalid as all its
	 *  individuals are absorbed by the current population.
	 *  <group>6-ancestral</group>
	 */
	void push(population & pop);

	/** Return the actual number of ancestral generations stored in a
	 *  population, which does not necessarily equal to the number set by
	 *  \c setAncestralDepth().
	 *  <group>6-ancestral</group>
	 */
	UINT ancestralGens() const
	{
		return m_ancestralPops.size();
	}

	/** Set information field \c idx (an index) of the current population to
	 *  \e values. \e values will be reused if its length is smaller than
	 *  <tt>popSize()</tt>.
	 *  <group>8-info</group>
	 */
	void setIndInfo(const vectorinfo & values, UINT idx);

	/** Set information field \c name of the current population to \e values.
	 *  \e values will be reused if its length is smaller than
	 *  <tt>popSize()</tt>.
	 *  <group>8-info</group>
	 */
	void setIndInfo(const vectorinfo & values, const string & name)
	{
		setIndInfo(values, infoIdx(name));
	}


	/** Set information field \c idx (an index) of a subpopulation
	 *  (<tt>subPop=sp</tt>) or a virtual subpopulation
	 *  (<tt>subPop=[sp, vsp]</tt>) to \e values. \e values will be reused if
	 *  its length is smaller than <tt>subPopSize(subPop)</tt>.
	 *  <group>8-info</group>
	 */
	void setIndInfo(const vectorinfo & values, UINT idx, vspID subPop);

	/** Set information field \c name of a subpopulation (<tt>subPop=sp</tt>) or
	 *  a virtual subpopulation (<tt>subPop=[sp, vsp]</tt>) to \e values.
	 *  \e values will be reused if its length is smaller than
	 *  <tt>subPopSize(subPop)</tt>.
	 *  <group>8-info</group>
	 */
	void setIndInfo(const vectorinfo & values, const string & name, vspID subPop)
	{
		setIndInfo(values, infoIdx(name), subPop);
	}


	/// CPPONLY info iterator
	IndInfoIterator infoBegin(UINT idx)
	{
		CHECKRANGEINFO(idx);

		// if there is virtual subpop, use individual based iterator
		// or
		// if requires order, but the information is not ordered
		// use individual based
		if (hasActivatedVirtualSubPop() || !indOrdered())
			return IndInfoIterator(idx, indBegin());
		else
			// if not required order, or if the information is ordered
			return IndInfoIterator(idx, m_info.begin() + idx, infoSize());
	}


	/// CPPONLY
	IndInfoIterator infoEnd(UINT idx)
	{
		CHECKRANGEINFO(idx);
		if (hasActivatedVirtualSubPop() || !indOrdered())
			return IndInfoIterator(idx, indEnd());
		else
			return IndInfoIterator(idx, m_info.begin() + idx + m_info.size(), infoSize());
	}


	/// CPPONLY info iterator
	IndInfoIterator infoBegin(UINT idx, vspID vsp)
	{
		SubPopID subPop = vsp.subPop();

		CHECKRANGEINFO(idx);
		CHECKRANGESUBPOP(subPop);
		CHECKRANGEVIRTUALSUBPOP(vsp.virtualSubPop());

		DBG_FAILIF(hasActivatedVirtualSubPop(subPop) && vsp.isVirtual(),
			ValueError, "Can not iterate through an activated subpopulation");

		// has to adjust order because of parameter subPop
		if (vsp.isVirtual()) {
			activateVirtualSubPop(subPop, vsp.virtualSubPop(), IteratableInds);
			return IndInfoIterator(idx, indBegin(subPop, IteratableInds));
		} else if (hasActivatedVirtualSubPop(subPop) || !indOrdered())
			return IndInfoIterator(idx, indBegin(subPop));
		else
			return IndInfoIterator(idx, m_info.begin() + idx + m_subPopIndex[subPop] * infoSize(), infoSize());
	}


	/// CPPONLY
	IndInfoIterator infoEnd(UINT idx, vspID vsp)
	{
		SubPopID subPop = vsp.subPop();

		CHECKRANGESUBPOP(subPop);
		CHECKRANGEVIRTUALSUBPOP(vsp.virtualSubPop());

		// has to adjust order because of parameter subPop
		if (vsp.isVirtual())
			return IndInfoIterator(idx, indEnd(subPop, IteratableInds));
		else if (hasActivatedVirtualSubPop(subPop) || !indOrdered())
			return IndInfoIterator(idx, indEnd(subPop));
		else
			return IndInfoIterator(idx, m_info.begin() + idx + m_subPopIndex[subPop + 1] * infoSize(), infoSize());
	}


	/** Return the information field \c idx (an index) of all individuals as a
	 *  list.
	 *  <group>8-info</group>
	 */
	vectorinfo indInfo(UINT idx)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This operation is not allowed when there is an activated virtual subpopulation");

		return vectorinfo(infoBegin(idx), infoEnd(idx));
	}


	/** Return the information field \c name of all individuals as a list.
	 *  <group>8-info</group>
	 */
	vectorinfo indInfo(const string & name)
	{
		UINT idx = infoIdx(name);

		return vectorinfo(infoBegin(idx), infoEnd(idx));
	}


	/** Return the information field \c idx (an index) of all individuals in
	 *  (virtual) subpopulation \e subPop as a list.
	 *  <group>8-info</group>
	 */
	vectorinfo indInfo(UINT idx, vspID subPop)
	{
		return vectorinfo(infoBegin(idx, subPop), infoEnd(idx, subPop));
	}


	/** Return the information field \c name of all individuals in
	 *  (virtual) subpopulation \e subPop as a list.
	 *  <group>8-info</group>
	 */
	vectorinfo indInfo(const string & name, vspID subPop)
	{
		UINT idx = infoIdx(name);

		return vectorinfo(infoBegin(idx, subPop), infoEnd(idx, subPop));
	}


	/**	Add an information field \e field to a population and initialize its
	 *  values to \e init.
	 *  <group>8-info</group>
	 */
	void addInfoField(const string & field, double init = 0);

	/** Add a list of information fields \e fields to a population and
	 *  initialize their values to \e init. If an information field alreay
	 *  exists, it will be re-initialized.
	 * <group>8-info</group>
	 */
	void addInfoFields(const vectorstr & fields, double init = 0);

	/** Set information fields \e fields to a population and initialize them
	 *  with value \e init. All existing information fields will be removed.
	 * <group>8-info</group>
	 */
	void setInfoFields(const vectorstr & fields, double init = 0);

	/** set the intended ancestral depth of a population to \e depth, which can
	 *  be \c 0 (does not store any ancestral generation), \c -1 (store all
	 *  ancestral generations), and a positive number (store \e depth ancestral
	 *  generations.
	 *  <group>6-ancestral</group>
	 */
	void setAncestralDepth(int depth);

	/** Making ancestral generation \e idx (\c 0 for current generation, \c 1
	 *  for parental generation, \c 2 for grand-parental generation, etc) the
	 *  current generation. This is an efficient way to access population
	 *  properties of an ancestral generation. <tt>useAncestralGen(0)</tt>
	 *  should always be called afterward to restore the correct order of
	 *  ancestral generations.
	 *  <group>6-ancestral</group>
	 */
	void useAncestralGen(UINT idx);

	/// CPPONLY compare two populations
	bool equalTo(const population & rhs)
	{
		return
		    genoStru() == rhs.genoStru() &&
		    m_subPopSize == rhs.m_subPopSize &&
		    m_inds == rhs.m_inds ;
	}


	//@}

	/// CPPONLY
	/**
	   some iterators requires that genotype information is within
	   each subpopulation. We need to adjust genotypic info to
	   obey this.
	   This function is const because the population is 'not changed'
	   conceptually.
	 */
	void sortIndividuals(bool infoOnly = false) const;

	/** Save population to a file \e filename, which can be loaded by a global
	 *  function <tt>LoadPopulation(filename)</tt>.
	 *  <group>1-pop</group>
	 */
	void save(const string & filename) const;

	/// CPPONLY load population from a file
	/**
	   filename load from filename

	   <group>1-pop</group>
	 */
	void load(const string & filename);

public:
	/// CPPONLY selection is on at any subpopulation?
	bool selectionOn() const
	{
		return !m_selectionFlags.empty();
	}


	/// CPPONLY
	bool selectionOn(UINT sp) const
	{
		DBG_ASSERT(m_selectionFlags.empty() || m_selectionFlags.size() == numSubPop(),
			IndexError, "Selection flags are wrong");
		return !m_selectionFlags.empty() && m_selectionFlags[sp];
	}


	/** HIDDEN
	    Turn off selection for all subpopulations
	   This is only used when you would like to apply two selectors. Maybe using two
	   different information fields.
	 */
	void turnOffSelection()
	{
		m_selectionFlags.clear();
	}


	/// CPPONLY
	void turnOnSelection(UINT sp)
	{
		if (m_selectionFlags.empty())
			m_selectionFlags.resize(numSubPop(), false);
		// there is an extreme case
		// selector turn on ...
		// split population...
		DBG_ASSERT(m_selectionFlags.size() == numSubPop(),
			SystemError, "Selection flags are wrong, did you split or merge populations after a selector is applied?");
		DBG_FAILIF(m_selectionFlags[sp], ValueError,
			"\nOnly one selector is allowed because each individual has only one fitness value\n"
			"If you need to select on more than one locus, use a multi-locus selector\n"
			"If you really want to apply another selector on the same population, call \n"
			"population::turnOffSelection() to walk around this restriction.\n");
		m_selectionFlags[sp] = true;
	}


	/// CPPONLY Turn on selection for all subpopulations
	void turnOnSelection()
	{
		if (m_selectionFlags.empty()) {
			m_selectionFlags.resize(numSubPop(), true);
			return;
		}
		// there is an extreme case
		// selector turn on ...
		// split population...
		DBG_ASSERT(m_selectionFlags.size() == numSubPop(),
			SystemError, "Selection flags are wrong, did you split or merge populations after a selector is applied?");
		DBG_FAILIF(true, ValueError,
			"\nOnly one selector is allowed because each individual has only one fitness value\n"
			"If you need to select on more than one locus, use a multi-locus selector\n"
			"If you really want to apply another selector on the same population, call \n"
			"population::turnOffSelection() to walk around this restriction.\n");
	}


public:
	/** CPPONLY
	 *  current replicate in a simulator which is not meaningful for a stand-alone population
	 *	<group>evolve</group>
	 */
	int rep()
	{
		return m_rep;
	}


	/// CPPONLY  set rep number
	void setRep(int rep, bool setVar = true)
	{
		m_rep = rep;
		if (setVar)
			m_vars.setIntVar("rep", rep);
	}


	/** CPPONLY
	    current generation during evolution
	    <group>evolve</group>
	 */
	ULONG gen()
	{
		return m_gen;
	}


	/// CPPONLY
	void setGen(ULONG gen, bool setVar = true)
	{
		m_gen = gen;
		if (setVar)
			m_vars.setIntVar("gen", gen);
	}


	/** return variables of a population as a Python dictionary.
	 *  <group>9-var</group>
	 */
	PyObject * vars();

	/** return a dictionary <tt>vars()["subPop"][subPop]</tt>. \e subPop can be
	 *  a number (<tt>subPop=spID</tt>), or a pair of numbers 
	 *  (<tt>subPop=(spID, vspID)</tt>). A \c ValueError will be raised if key
	 *  \c 'subPop' does not exist in \c vars(), or if key \e subPop does not
	 *  exist in <tt>vars()["subPop"]</tt>.
	 *  <group>9-var</group>
	 */
	PyObject * vars(vspID subPop);


	/// CPPONLY The same as vars(), but without increasing reference count.
	PyObject * dict(int subPop = -1);

	/// CPPONLY
	void setDict(PyObject * dict)
	{
		DBG_ASSERT(dict != NULL, SystemError, "Dictionary is empty");
		m_vars.setDict(dict);
	}


	/// CPPONLY
	bool hasVar(const string & name)
	{
		return m_vars.hasVar(name);
	}


	/// CPPONLY
	void removeVar(const string & name)
	{
		m_vars.removeVar(name);
	}


	/// CPPONLY
	PyObject * setBoolVar(const string & name, const bool val)
	{
		return m_vars.setBoolVar(name, val);
	}


	/// CPPONLY
	PyObject * setIntVar(const string & name, const int val)
	{
		return m_vars.setIntVar(name, val);
	}


	/// CPPONLY
	PyObject * setDoubleVar(const string & name, const double val)
	{
		return m_vars.setDoubleVar(name, val);
	}


	/// CPPONLY
	PyObject * setStringVar(const string & name, const string & val)
	{
		return m_vars.setStringVar(name, val);
	}


	///CPPONLY
	PyObject * setIntVectorVar(const string & name, const vectori & val)
	{
		return m_vars.setIntVectorVar(name, val);
	}


	///CPPONLY
	PyObject * setDoubleVectorVar(const string & name, const vectorf & val)
	{
		return m_vars.setDoubleVectorVar(name, val);
	}


	/// CPPONLY
	PyObject * setStrDictVar(const string & name, const strDict & val)
	{
		return m_vars.setStrDictVar(name, val);
	}


	/// CPPONLY
	PyObject * setIntDictVar(const string & name, const intDict & val)
	{
		return m_vars.setIntDictVar(name, val);
	}


	/// CPPONLY
	PyObject * setVar(const string & name, PyObject * val)
	{
		return m_vars.setVar(name, val);
	}


	/// CPPONLY
	PyObject * getVar(const string & name, bool nameError = true)
	{
		return m_vars.getVar(name, nameError);
	}


	/// CPPONLY
	bool getVarAsBool(const string & name, bool nameError = true)
	{
		return m_vars.getVarAsBool(name, nameError);
	}


	/// CPPONLY
	int getVarAsInt(const string & name, bool nameError = true)
	{
		return m_vars.getVarAsInt(name, nameError);
	}


	/// CPPONLY
	double getVarAsDouble(const string & name, bool nameError = true)
	{
		return m_vars.getVarAsDouble(name, nameError);
	}


	/// CPPONLY
	string getVarAsString(const string & name, bool nameError = true)
	{
		return m_vars.getVarAsString(name, nameError);
	}


	/// CPPONLY
	strDict getVarAsStrDict(const string & name, bool nameError = true)
	{
		return m_vars.getVarAsStrDict(name, nameError);
	}


	/// CPPONLY
	intDict getVarAsIntDict(const string & name, bool nameError = true)
	{
		return m_vars.getVarAsIntDict(name, nameError);
	}


	/// CPPONLY
	string varsAsString() const
	{
		return m_vars.asString();
	}


	/// CPPONLY
	void varsFromString(const string & vars)
	{
		return m_vars.fromString(vars);
	}


	/** HIDDEN
	 *  evaluate a Python statment/expression in the population's local namespace
	 *  This function evaluates a Python statment( \c stmts )/expression( \c expr )
	 *  and return its result as a string. Optionally run statement( \c stmts ) first.
	 *  <group>8-var</group>
	 */
	PyObject * evaluate(const string & expr = "", const string & stmts = "")
	{
		return Expression(expr, stmts, m_vars.dict() ).evaluate();
	}


	/** HIDDEN
	 *  execute a statement (can be a multi-line string) in the population's local namespace
	 */
	void execute(const string & stmts = "")
	{
		Expression("", stmts, m_vars.dict() ).evaluate();
	}


	/** HIDDEN
	 *  This function changes individual genotype pointer and information
	 *  pointer so that the population appear to be after a massive migration.
	 *  It is used to test population functions and see if they behave the
	 *  same under this situation.
	 */
	void scramble();

private:
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const UINT version) const
	{
		// deep adjustment: everyone in order
		const_cast<population *>(this)->sortIndividuals();

		ar & ModuleMaxAllele;

		DBG_DO(DBG_POPULATION, cout << "Handling geno structure" << endl);
		// GenoStructure genoStru = this->genoStru();
		ar & genoStru();

		ar & m_subPopSize;
		DBG_DO(DBG_POPULATION, cout << "Handling genotype" << endl);
#ifdef BINARYALLELE
		size_t size = m_genotype.size();
		ar & size;
		WORDTYPE * ptr = BITPTR(m_genotype.begin());
		size_t blks = size / WORDBIT;
		size_t rest = size - blks * WORDBIT;
		DBG_ASSERT(WORDBIT >= 32, SystemError, "WordBit should be at least 32 bits");

		WORDTYPE tmp, tmp1;
		for (size_t i = 0; i < blks; ++i) {
			tmp = *ptr++;
			for (size_t j = 0; j < WORDBIT / 32; ++j) {
				tmp1 = tmp & 0xFFFFFFFF;
				tmp = tmp >> 32;
				ar & tmp1;
			}
		}
		// last block
		if (rest > 0) {
			tmp = *ptr;
			for (size_t j = 0; j <= (rest - 1) / 32; ++j) {
				tmp1 = tmp & 0xFFFFFFFF;
				tmp = tmp >> 32;
				ar & tmp1;
			}
		}
#else
		ar & m_genotype;
#endif
		DBG_DO(DBG_POPULATION, cout << "Handling information" << endl);
		ar & m_info;
		DBG_DO(DBG_POPULATION, cout << "Handling individuals" << endl);
		ar & m_inds;
		DBG_DO(DBG_POPULATION, cout << "Handling ancestral populations" << endl);
		ar & m_ancestralGens;
		size_t sz = m_ancestralPops.size();
		ar & sz;
		for (size_t i = 0; i < m_ancestralPops.size(); ++i) {
			const_cast<population *>(this)->useAncestralGen(i + 1);
			// need to make sure ancestral pop also in order
			const_cast<population *>(this)->sortIndividuals();
			ar & m_subPopSize;
#ifdef BINARYALLELE
			size_t size = m_genotype.size();
			ar & size;
			WORDTYPE * ptr = BITPTR(m_genotype.begin());
			size_t blks = size / WORDBIT;
			size_t rest = size - blks * WORDBIT;
			DBG_ASSERT(WORDBIT >= 32, SystemError, "WordBit should be at least 32 bits");

			WORDTYPE tmp, tmp1;
			for (size_t i = 0; i < blks; ++i) {
				tmp = *ptr++;
				for (size_t j = 0; j < WORDBIT / 32; ++j) {
					tmp1 = tmp & 0xFFFFFFFF;
					tmp = tmp >> 32;
					ar & tmp1;
				}
			}
			// last block
			if (rest > 0) {
				tmp = *ptr;
				// rest = 1-31: (rest-1)/32=0, j <= rest/32 = 0
				// rest = 32; j <= (rest-1)/32 = 0
				for (size_t j = 0; j <= (rest - 1) / 32; ++j) {
					tmp1 = tmp & 0xFFFFFFFF;
					tmp = tmp >> 32;
					ar & tmp1;
				}
			}
#else
			ar & m_genotype;
#endif
			ar & m_info;
			ar & m_inds;
		}
		const_cast<population *>(this)->useAncestralGen(0);

		// save shared variables as string.
		// note that many format are not supported.
		try {
			DBG_DO(DBG_POPULATION, cout << "Handling shared variables" << endl);
			string vars = varsAsString();
			ar & vars;
		} catch (...) {
			cout << "Warning: shared variable is not saved correctly.\npopulation should still be usable." << endl;
		}
	}


	template<class Archive>
	void load(Archive & ar, const UINT version)
	{
		ULONG ma;
		ar & ma;

		if (ma > ModuleMaxAllele)
			cout << "Warning: The population is saved in library with more allele states. \n"
			     << "Unless all alleles are less than " << ModuleMaxAllele
			     << ", you should use the modules used to save this file. (c.f. simuOpt.setOptions()\n";

		GenoStructure stru;
		DBG_DO(DBG_POPULATION, cout << "Handling geno structure" << endl);
		ar & stru;
		ar & m_subPopSize;
		DBG_DO(DBG_POPULATION, cout << "Handling genotype" << endl);

#ifdef BINARYALLELE
		// binary from binary
		if (ma == 1) {
			size_t size;
			ar & size;
			size_t blks = size / WORDBIT;
			size_t rest = size - blks * WORDBIT;

			m_genotype.resize(size);
			WORDTYPE tmp, tmp1;
			WORDTYPE * ptr = BITPTR(m_genotype.begin());
			for (size_t i = 0; i < blks; ++i) {
				tmp = 0;
				for (size_t j = 0; j < WORDBIT / 32; ++j) {
					ar & tmp1;
					tmp |= tmp1 << (j * 32);
				}
				*ptr++ = tmp;
			}
			// last block
			if (rest > 0) {
				tmp = 0;
				for (size_t j = 0; j <= (rest - 1) / 32; ++j) {
					ar & tmp1;
					tmp |= tmp1 << (j * 32);
				}
				*ptr = tmp;
			}
		}
		// binary from others (long types)
		else {
			DBG_DO(DBG_POPULATION, cout << "Load bin from long. " << endl);
			vector<unsigned char> tmpgeno;
			ar & tmpgeno;
			m_genotype = vectora(tmpgeno.begin(), tmpgeno.end());
		}
#else
		// long from binary
		if (ma == 1) {
			// for version 2 and higher, archive in 32bit blocks.
			size_t size;
			ar & size;
			m_genotype.resize(size);
			size_t blks = size / 32;
			size_t rest = size - blks * 32;
			DBG_DO(DBG_POPULATION, cout << "Load long from bin. " << size << " rest " << rest << endl);
			DBG_ASSERT(WORDBIT >= 32, SystemError, "WordBit should be at least 32 bits");

			GenoIterator ptr = m_genotype.begin();
			WORDTYPE tmp;
			for (size_t i = 0; i < blks; ++i) {
				ar & tmp;
				for (size_t j = 0; j < 32; ++j) {
					*ptr++ = (tmp & 1UL) != 0;
					tmp = tmp >> 1;
				}
			}
			// last block
			if (rest > 0) {
				ar & tmp;
				for (size_t j = 0; j < rest; ++j) {
					*ptr++ = (tmp & 1UL) != 0;
					tmp = tmp >> 1;
				}
			}
		}                                                                               // if ma == 1
		else {                                                                          // for non-binary types, ...
			DBG_DO(DBG_POPULATION, cout << "Load long from long. " << endl);
			// long from long
			ar & m_genotype;
		}
#endif

		DBG_DO(DBG_POPULATION, cout << "Handling info" << endl);
		ar & m_info;

		DBG_DO(DBG_POPULATION, cout << "Handling individuals" << endl);
		ar & m_inds;

		// set genostructure, check duplication
		// we can not use setGenoStruIdx since stru may be new.
		this->setGenoStructure(stru);

		m_popSize = accumulate(m_subPopSize.begin(), m_subPopSize.end(), 0L);

		DBG_FAILIF(m_info.size() != m_popSize * infoSize(), ValueError, "Wgong size of info vector");

		if (m_popSize != m_inds.size() ) {
			cout << "Number of individuals loaded" << m_inds.size() << endl;
			cout << "population size" << m_popSize << endl;
			throw ValueError("Number of individuals does not match population size.\n"
				             "Please use the same (binary, short or long) module to save and load files.");
		}

		DBG_DO(DBG_POPULATION, cout << "Reconstruct individual genotype" << endl);
		m_subPopIndex.resize(m_subPopSize.size() + 1);
		UINT i = 1;
		for (m_subPopIndex[0] = 0; i <= m_subPopSize.size(); ++i)
			m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];

		// assign genotype location and set structure information for individuals
		GenoIterator ptr = m_genotype.begin();
		UINT step = genoSize();
		InfoIterator infoPtr = m_info.begin();
		UINT infoStep = infoSize();
		for (ULONG i = 0; i < m_popSize; ++i, ptr += step, infoPtr += infoStep) {
			m_inds[i].setGenoStruIdx(genoStruIdx());
			m_inds[i].setGenoPtr(ptr);
			m_inds[i].setInfoPtr(infoPtr);
		}
		m_ancestralGens = 0;
		m_ancestralPops.clear();

		// ancestry populations
		DBG_DO(DBG_POPULATION, cout << "Handling ancestral populations" << endl);
		ar & m_ancestralGens;
		size_t na;
		ar & na;
		for (size_t ap = 0; ap < na; ++ap) {
			popData pd;
			ar & pd.m_subPopSize;
#ifdef BINARYALLELE
			// binary from binary
			if (ma == 1) {
				DBG_DO(DBG_POPULATION, cout << "Load bin from bin. " << endl);
				size_t size;
				ar & size;
				size_t blks = size / WORDBIT;
				size_t rest = size - blks * WORDBIT;

				pd.m_genotype.resize(size);
				WORDTYPE * ptr = BITPTR(pd.m_genotype.begin());
				WORDTYPE tmp, tmp1;
				for (size_t i = 0; i < blks; ++i) {
					tmp = 0;
					for (size_t j = 0; j < WORDBIT / 32; ++j) {
						ar & tmp1;
						tmp |= tmp1 << (j * 32);
					}
					*ptr++ = tmp;
				}
				// last block
				if (rest > 0) {
					tmp = 0;
					for (size_t j = 0; j <= (rest - 1) / 32; ++j) {
						ar & tmp1;
						tmp |= tmp1 << (j * 32);
					}
					*ptr = tmp;
				}
			} else {
				DBG_DO(DBG_POPULATION, cout << "Load bin from long. " << endl);
				// binary from long types
				vector<unsigned char> tmpgeno;
				ar & tmpgeno;
				pd.m_genotype = vectora(tmpgeno.begin(), tmpgeno.end());
			}
#else
			if (ma == 1) {
				// long type from binary
				size_t size;
				ar & size;
				pd.m_genotype.resize(size);
				size_t blks = size / 32;
				size_t rest = size - blks * 32;
				DBG_DO(DBG_POPULATION, cout << "Load long from bin. " << size << " rest " << rest << endl);

				ptr = pd.m_genotype.begin();
				WORDTYPE tmp;
				for (size_t i = 0; i < blks; ++i) {
					ar & tmp;
					for (size_t j = 0; j < 32; ++j) {
						*ptr++ = (tmp & 1UL) != 0;
						tmp = tmp >> 1;
					}
				}
				// last block
				if (rest > 0) {
					ar & tmp;
					for (size_t i = 0; i < rest; ++i) {
						*ptr++ = (tmp & 1UL) != 0;
						tmp = tmp >> 1;
					}
				}
			} else {
				DBG_DO(DBG_POPULATION, cout << "Load long from long. " << endl);
				// long type from long type.
				ar & pd.m_genotype;
			}
#endif
			ar & pd.m_info;
			ar & pd.m_inds;
			// set pointer after copy this thing again (push_back)
			m_ancestralPops.push_back(pd);
			// now set pointers
			popData & p = m_ancestralPops.back();
			// set pointers
			vector<individual> & inds = p.m_inds;
			ULONG ps = inds.size();
			ptr = p.m_genotype.begin();
			infoPtr = p.m_info.begin();

			for (ULONG i = 0; i < ps; ++i, ptr += step, infoPtr += infoStep) {
				inds[i].setGenoPtr(ptr);
				inds[i].setInfoPtr(infoPtr);
				// set new genoStructure
				inds[i].setGenoStruIdx(genoStruIdx());
			}
		}

		// load vars from string
		try {
			DBG_DO(DBG_POPULATION, cout << "Handling shared variables" << endl);
			string vars;
			ar & vars;
			varsFromString(vars);
		} catch (...) {
			cout << "Warning: shared variable is not loaded correctly.\npopulation should still be usable." << endl;
		}

		setIndOrdered(true);
	}


	BOOST_SERIALIZATION_SPLIT_MEMBER();

private:
	/// population size: number of individual
	ULONG m_popSize;

	/// size of each subpopulation
	vectorlu m_subPopSize;

	/// index to subPop \todo change to vectorl
	vectorlu m_subPopIndex;

	///
	vspSplitter * m_vspSplitter;

	/// pool of genotypic information
	vectora m_genotype;

	/// information
	/// only in head node
	vectorinfo m_info;

	/// individuals.
	/// only in head node?
	vector<individual> m_inds;

	int m_ancestralGens;

	/// shared variables for this population
	SharedVariables m_vars;

	/// store previous populations
	/// need to store: subPopSize, genotype and m_inds
	struct popData
	{
		vectorlu m_subPopSize;
		vectora m_genotype;
		vectorinfo m_info;
		vector<individual> m_inds;
		bool m_indOrdered;

		// swap between a popData and existing data.
		void swap(population & pop);

	};

	std::deque<popData> m_ancestralPops;

	/// curent replicate number
	int m_rep;

	/// generation
	ULONG m_gen;

	/// current ancestral depth
	UINT m_curAncestralGen;

	/// whether or not individual genotype and information are in order
	/// within a population.
	mutable bool m_indOrdered;

	/// selection flags for each subpopulation.
	/// empty means no selection
	vector<bool> m_selectionFlags;

};

/** load a population from a file.
 */
population & LoadPopulation(const string & file);

}


#ifndef SWIG
#  ifndef _NO_SERIALIZATION_
// version 0: base (reset for version 1.0)
BOOST_CLASS_VERSION(simuPOP::population, 0)
#  endif
#endif
#endif
