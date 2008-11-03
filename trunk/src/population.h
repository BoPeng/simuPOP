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

enum RelativeType {
	REL_None,           // do nothing
	REL_Self,           // individual himself or herself.
	REL_Offspring,      // All offspring with all spouses (if there are more than one spouse)
	REL_Spouse,         // All spouses (with at least one offspring)
	REL_FullSibling,    // Siblings who share two parents
	REL_Sibling,        // Siblings who share at least one parent
};

enum SexChoice {
	AnySex = 0,
	MaleOnly = 1,
	FemaleOnly = 2,
	OppositeSex = 3
};

namespace simuPOP {

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
	 *    regions of the sex chromsomes is not supported. A \c Mitochondrial
	 *    chromosome only exists in females and is inherited maternally.
	 *  \param lociPos Positions of all loci on all chromosome, as a list of
	 *    float numbers. Default to \c 1, \c 2, ... etc on each chromosome.
	 *    Positions on the same chromosome should be ordered. A nested list
	 *    that specifies positions of loci on each chromosome is also
	 *    acceptable.
	 *  \param ancestralGens Number of the most recent ancestral generations
	 *    to keep during evolution. Default to \c 0, which means only the
	 *    current generation will be kept. If it is set to \c -1, all ancestral
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
		int ancestralGens = 0,
		const vectorstr & chromNames = vectorstr(),
		const vectorstr & alleleNames = vectorstr(),
		const vectorstr & lociNames = vectorstr(),
		const vectorstr & infoFields = vectorstr());

	/// CPPONLY copy constructor
	population(const population & rhs);

	/** Copy a population, with the option to keep all (default), no, or a
	 *  given number of ancestral generations (\e keepAncestralPops = \c -1,
	 *  \c 0, or a positive number, respectively). Note that Python statement
	 *  <tt>pop1 = pop</tt> only creates a reference to an existing population
	 *  \c pop.
	 *  <group>1-pop</group>
	 */
	population * clone(int keepAncestralPops = -1) const;

	/** HIDDEN (do not see a need to expose this function yet.)
	 *  swap the content of two populations
	 *  <group>1-pop</group>
	 */
	void swap(population & rhs)
	{
		GenoStruTrait::swap(rhs);
		std::swap(m_popSize, rhs.m_popSize);
		std::swap(m_numSubPop, rhs.m_numSubPop);
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

	/** Set a VSP splitter to the population, which defines the same VSPs for
	 *  all subpopulations. If different VSPs are needed for different
	 *  subpopulations, a \c combinedSplitter can be used to make these VSPs
	 *  available to all subpopulations.
	 *  <group>3-VSP</group>
	 */
	void setVirtualSplitter(vspSplitter * vsp);

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
		vspSplitter::activateType type = vspSplitter::Visible);

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
		return m_numSubPop;
	}


	/** Return the size of a subpopulation (<tt>subPopSize(sp)</tt>) or a
	 *  virtual subpopulation (<tt>subPopSize([sp, vsp])<tt>).
	 *  <group>2-subpop</group>
	 */
	ULONG subPopSize(vspID vsp) const
	{
		CHECKRANGESUBPOP(vsp.subPop());
		CHECKRANGEVIRTUALSUBPOP(vsp.virtualSubPop());
		if (hasActivatedVirtualSubPop() || vsp.isVirtual())
			return m_vspSplitter->size(*this, vsp.subPop(), vsp.virtualSubPop());
		else
			return m_subPopSize[vsp.subPop()];
	}


	/** Return the name of a virtual subpopulation \e vsp (specified by a
	 *  <tt>(sp, vsp)</tt> pair). Because VSP names are the same across all
	 *  subpopulations, a single <tt>vsp</tt> index is also acceptable.
	 *  <group>3-VSP</group>
	 */
	string virtualSubPopName(vspID vsp) const;

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

		for (UINT i = 1; i <= m_numSubPop; ++i) {
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

	/** Return a refernce to individual \e ind in subpopulation \e subPop.
	 *  <group>4-ind</group>
	 */
	individual & ind(ULONG idx, UINT subPop = 0)
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
	 *  is not the present one (see \c useAncestralGen).
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
	 *  in a subpopulation (<tt>vsp=spID</tt>) or a virtual subpopulation
	 *  (<tt>vsp=[spID, vspID]</tt>).
	 *  <group>4-ind</group>
	 */
	pyIndIterator individuals(vspID vsp)
	{
		SubPopID spID = vsp.subPop();
		SubPopID vspID = vsp.virtualSubPop();

#ifndef OPTIMIZED
		CHECKRANGESUBPOP(spID);
		CHECKRANGEVIRTUALSUBPOP(vspID);
		DBG_FAILIF(hasActivatedVirtualSubPop(spID), ValueError,
			"This operation is not allowed for an activated subpopulation");
#endif
		if (vsp.isVirtual()) {
			// this does not need to be deactivated...
			activateVirtualSubPop(spID, vspID, vspSplitter::Iteratable);
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
	bool indOrdered()
	{
		return m_indOrdered;
	}


	/// CPPONLY
	void setIndOrdered(bool s)
	{
		m_indOrdered = s;
	}


	/// CPPONLY individual iterator: without subPop info
	IndIterator indBegin()
	{
		return IndIterator(m_inds.begin(), m_inds.end(),
			!hasActivatedVirtualSubPop());
	}


	/// CPPONLY individual iterator: without subPop info
	IndIterator indEnd()
	{
		return IndIterator(m_inds.end(), m_inds.end(),
			!hasActivatedVirtualSubPop());
	}


	/** CPPONLY individual iterator: with subPop info.
	 *  The iterator will skip invisible individuals
	 */
	IndIterator indBegin(UINT subPop)
	{
		CHECKRANGESUBPOP(subPop);

		return IndIterator(m_inds.begin() + m_subPopIndex[subPop],
			m_inds.begin() + m_subPopIndex[subPop + 1],
			!hasActivatedVirtualSubPop(subPop));
	}


	/// CPPONLY individual iterator: with subPop info.
	IndIterator indEnd(UINT subPop)
	{
		CHECKRANGESUBPOP(subPop);

		return IndIterator(m_inds.begin() + m_subPopIndex[subPop + 1],
			m_inds.begin() + m_subPopIndex[subPop + 1],
			!hasActivatedVirtualSubPop(subPop));
	}


	/** CPPONLY individual iterator: without subPop info
	 *  The iterator will skip invisible individuals
	 */
	ConstIndIterator indBegin() const
	{
		return ConstIndIterator(m_inds.begin(), m_inds.end(),
			!hasActivatedVirtualSubPop());
	}


	/** CPPONLY individual iterator: without subPop info
	 *  It is recommended to use it.valid(), instead of it != indEnd()
	 */
	ConstIndIterator indEnd() const
	{
		return ConstIndIterator(m_inds.end(), m_inds.end(),
			!hasActivatedVirtualSubPop());
	}


	/** CPPONLY individual iterator: with subPop info.
	 *  The iterator will skip invisible individuals
	 */
	ConstIndIterator indBegin(UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return ConstIndIterator(m_inds.begin() + m_subPopIndex[subPop],
			m_inds.begin() + m_subPopIndex[subPop + 1],
			!hasActivatedVirtualSubPop(subPop));
	}


	/** CPPONLY individual iterator: with subPop info.
	 * It is recommended to use it.valid(), instead of it != indEnd(sp)
	 */
	ConstIndIterator indEnd(UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return ConstIndIterator(m_inds.begin() + m_subPopIndex[subPop + 1],
			m_inds.begin() + m_subPopIndex[subPop + 1],
			!hasActivatedVirtualSubPop(subPop));
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

	/// HIDDEN set subpopulation ID with given ID
	/**
	   Set subpopulation ID of each individual with given ID. Individuals
	   can be rearranged afterwards using \c setSubPopByIndID.

	   \param id an array of the same length of population size, resprenting
	   subpopulation ID of each individual. If the length of \id is less
	   than population size, it is repeated to fill the whole population.
	   \param ancestralPops If true (default to False), set subpop id for ancestral
	   generations as well.
	   \sa individual::setSubPopID, individual::subPopID
	 */
	void setIndSubPopID(const vectori & id, bool ancestralPops = false);

	/// HIDDEN set subpopulation ID of each individual with their current subpopulation ID
	/**
	   \param ancestralPops If true (default to False), set subpop id for ancestral
	   generations as well.
	 */
	void setIndSubPopIDWithID(bool ancestralPops = false);

	/// HIDDEN move individuals to subpopulations according to individual subpopulation IDs
	/**
	   Rearrange individuals to their new subpopulations according to their
	   subpopulation ID (or the new given \c ID). Order within each subpopulation is not respected.

	   \param id new subpopulation ID, if given, current individual subpopulation ID
	   will be ignored.
	   \note Individual with negative info will be removed!
	   \sa setIndSubPopID
	 */
	void setSubPopByIndID(vectori id = vectori());

	/// split a subpopulation into subpopulations of given sizes
	/**
	   The sum of given sizes should be equal to the size of the split subpopulation. Subpopulation
	   IDs can be specified. The subpopulation IDs of non-split subpopulations will be kept. For example, if
	   subpopulation 1 of 0 1 2 3 is split into three parts, the new subpop id will be
	   0 (1 4 5) 2 3.
	   \note \c subpop with negative ID will be removed. So, you can shrink one \c subpop by splitting and
	   setting one of the new \c subpop with negative ID.

	 * <group>7-manipulate</group>
	 */
	void splitSubPop(UINT which, vectorlu sizes, vectoru subPopID = vectoru());

	/// split a subpopulation into subpopulations of given proportions
	/**
	   The sum of given proportions should add up to one. Subpopulation IDs can be specified.
	   \c subpop with negative ID will be removed. So, you can shrink one \c subpop by splitting and
	   setting one of the new \c subpop with negative ID.
	 * <group>7-manipulate</group>
	 */
	void splitSubPopByProportion(UINT which, vectorf proportions, vectoru subPopID = vectoru());

	/** remove empty subpopulations by adjusting subpopulation IDs
	 * <group>7-manipulate</group>
	 */
	void removeEmptySubPops();

	/// remove subpopulations and adjust subpopulation IDs so that there will be no \em 'empty' subpopulation left
	/**
	   Remove specified subpopulations (and all individuals within). If \c shiftSubPopID is false, \c subPopID
	   will be kept intactly.
	 * <group>7-manipulate</group>
	 */
	void removeSubPops(const vectoru & subPops = vectoru(), bool shiftSubPopID = true, bool removeEmptySubPops = false);

	/** remove individuals. If a valid \c subPop is given, remove individuals from this subpopulation. Indexes in \c inds will be treated as relative indexes.
	 * <group>7-manipulate</group>
	 */
	void removeIndividuals(const vectoru & inds = vectoru(), int subPop = -1, bool removeEmptySubPops = false);

	/**
	   Merge subpopulations, the first subpopulation ID (the first one in array
	   \c subPops) will be used as the ID of the new subpopulation. That is to
	   say, all merged subpopulations will take the ID of the first one. The
	   subpopulation ID of the empty subpopulations will be kept (so that other
	   subpopulations are unaffected, unless they are removed by <tt>removeEmptySubPops = True</tt>).
	 * <group>7-manipulate</group>
	 */
	void mergeSubPops(vectoru subPops = vectoru(), bool removeEmptySubPops = false);

	/** Add all individuals, including ancestors, in \e pop to the current
	 *  population. Two populations should have the same genotypic structures
	 *  and number of ancestral generations. Subpopulations in population
	 *  \e pop are kept.
	 *  <group>7-manipulate</group>
	 */
	void addIndFrom(const population & pop);

	/**
	 * <group>7-manipulate</group>
	   merge populations by loci
	   Two populations should have the same number of individuals. This also holds for
	   any ancestral generations. By default, chromosomes of \c pop are appended to the current
	   population.	You can change this arrangement in two ways
	   \li specify new chromosome structure using parameter \c newLoci and \c newLociPos. Loci from new and old
	    populations are still in their original order, but chromosome number and positions
	    can be changed in this way.
	   \li specify \c byChromosome=true so that chromosomes will be merged one by one. In this
	    case, loci position of two popualtions are important because loci will be arranged
	    in the order of loci position; and identical loci position of two loci in two
	   populations will lead to error.

	   \param newNumLoci the new number of loci for the combined genotypic structure.
	   \param newLociPos the new loci position if number of loci on each chromosomes are
	   changed with \c newNumLoci. New loci positions should be in order on the new chromosomes.
	   \param byChromosome merge chromosome by chromosome, loci are ordered by loci position
	   Default to \c False.
	   \note \li Information fields are not merged.
	   \li All ancestral generations are merged because all individuals in a
	   population have to have the same genotypic structure.
	 */
	void mergePopulationByLoci(const population & pop, const vectoru & newNumLoci = vectoru(),
		const vectorf & newLociPos = vectorf(), bool byChromosome = false);

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

	/// resize current population
	/**
	   Resize population by giving new subpopulation sizes.
	   \param newSubPopSizes an array of new subpopulation sizes. If there
	   is only one subpopulation, use <tt>[newPopSize]</tt>.
	   \param propagate if \c propagate is \c true, copy individuals to new comers.
	   I.e., 1, 2, 3 ==> 1, 2, 3, 1, 2, 3, 1
	   \note This function only resizes the current generation.

	 * <group>7-manipulate</group>
	 */
	void resize(const vectorlu & newSubPopSizes, bool propagate = false);

	/// reorder subpopulations by \c order or by \c rank
	/**
	   \param order new order of the subpopulations. For examples, 3 2 0 1
	   means \c subpop3, \c subpop2, \c subpop0, \c subpop1 will be the new layout.
	   \param rank you may also specify a new rank for each subpopulation. For example, 3,2,0,1
	   means the original subpopulations will have new IDs 3,2,0,1, respectively. To achive order 3,2,0,1,
	   the rank should be 1 0 2 3.


	 * <group>7-manipulate</group>
	 */
	void reorderSubPops(const vectoru & order = vectoru(), const vectoru & rank = vectoru(),
		bool removeEmptySubPops = false);

	/**
	   Form a new population according to individual subpopulation ID. Individuals with negative subpopulation
	   ID will be removed.
	 * <group>7-manipulate</group>
	 */
	population & newPopByIndID(int keepAncestralPops = -1,
		const vectori & id = vectori(),
		bool removeEmptySubPops = false);

	/** remove some loci from the current population. Only one of the two parameters can be specified.
	 * <group>7-manipulate</group>
	 */
	void removeLoci(const vectoru & remove = vectoru(), const vectoru & keep = vectoru());

	/// obtain a new population with selected loci
	/**
	   Copy current population to a new one with selected loci \c keep or remove specified loci \c remove
	   (no change on the current population), equivalent to \n <tt>y=x.clone</tt> \n <tt>y.removeLoci(remove, keep)</tt>
	 * <group>7-manipulate</group>
	 */
	population & newPopWithPartialLoci(
		const vectoru & remove = vectoru(),
		const vectoru & keep = vectoru());

	/// HIDDEN absorb \c rhs population as the current generation of a population
	/**
	   This function is used by a simulator to push offspring generation \c rhs
	   to the current population, while the current population is pushed
	   back as an ancestral population (if <tt>ancestralDepath() != 0</tt>). Because \c rhs
	   population is swapped in, \c rhs will be empty after this operation.
	 */
	void pushAndDiscard(population & rhs, bool force = false);

	/** Return the actual number of ancestral generations stored in a
	 *  population, which does not necessarily equal to the number set by
	 *  \c setAncestralDepth().
	 *  <group>6-ancestral</group>
	 */
	UINT ancestralGens() const
	{
		return m_ancestralPops.size();
	}


	/** CPPONLY because I do not see a Python level use case of this function.
	 *  Current ancestral population activated by \c useAncestralGen(). There can be
	 *  several ancestral generations in a population. \c 0 (current), \c 1 (parental)
	 *  etc. When \c useAncestralGen(gen) is used, current generation is set to
	 *  one of the parental generations, which is the information returned by this
	 *  function. \c useAncestralGen(0) should always be used to set a population
	 *  to its usual ancestral order after operations to the ancestral generation are done.
	 *
	 * <group>6-ancestral</group>
	 */
	UINT ancestralGen() const
	{
		return m_curAncestralGen;
	}


	/// set individual information for the given information field \c idx
	/**
	 * <group>8-info</group>
	 */
	template<typename T, typename T1>
	void setIndInfo(const T & values, T1 idx)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This operation is not allowed when there is an activated virtual subpopulation");

		CHECKRANGEINFO(idx);
		DBG_ASSERT(values.size() == popSize(), IndexError,
			"Size of values should be the same as population size");
		typename T::const_iterator infoIter = values.begin();
		for (IndInfoIterator ptr = infoBegin(idx); ptr != infoEnd(idx); ++ptr)
			*ptr = static_cast<InfoType>(*infoIter++);
	}


	/// set individual information for the given information field \c name
	/**
	   <tt>x.setIndInfo(values, name)</tt> is
	   equivalent to the \c idx version <tt>x.setIndInfo(values, x.infoIdx(name))</tt>.
	 * <group>8-info</group>
	 */
	template<class T>
	void setIndInfo(const T & values, const string & name)
	{
		// for mpi version , use gloal idx
		int idx = infoIdx(name);

		setIndInfo<T, UINT>(values, idx);
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
	IndInfoIterator infoBegin(UINT idx, UINT subPop)
	{
		CHECKRANGEINFO(idx);
		CHECKRANGESUBPOP(subPop);

		// has to adjust order because of parameter subPop
		if (hasActivatedVirtualSubPop(subPop) || !indOrdered())
			return IndInfoIterator(idx, indBegin(subPop));
		else
			return IndInfoIterator(idx, m_info.begin() + idx + m_subPopIndex[subPop] * infoSize(), infoSize());
	}


	/// CPPONLY
	IndInfoIterator infoEnd(UINT idx, UINT subPop)
	{
		CHECKRANGEINFO(idx);
		CHECKRANGESUBPOP(subPop);

		// has to adjust order because of parameter subPop
		if (hasActivatedVirtualSubPop(subPop) || !indOrdered())
			return IndInfoIterator(idx, indEnd(subPop));
		else
			return IndInfoIterator(idx, m_info.begin() + idx + m_subPopIndex[subPop + 1] * infoSize(), infoSize());
	}


	/// get information field \c idx of all individuals
	/**

	 * <group>8-info</group>
	 */
	vectorinfo indInfo(UINT idx)
	{
		DBG_FAILIF(hasActivatedVirtualSubPop(), ValueError,
			"This operation is not allowed when there is an activated virtual subpopulation");

		return vectorinfo(infoBegin(idx), infoEnd(idx));
	}


	/// get information field \c name of all individuals
	/**
	   \param name name of the information field
	   \return a vector with value of the information field

	 * <group>8-info</group>
	 */
	vectorinfo indInfo(const string & name)
	{
		UINT idx = infoIdx(name);

		return vectorinfo(infoBegin(idx), infoEnd(idx));
	}


	/// get information field \c idx of all individuals in a subpopulation \c subPop
	/**
	 * <group>8-info</group>
	 */
	vectorinfo indInfo(UINT idx, UINT subPop)
	{
		return vectorinfo(infoBegin(idx, subPop),
			infoEnd(idx, subPop));
	}


	/// get information field \c name of all individuals in a subpopulation \c subPop
	/**
	 * <group>8-info</group>
	 */
	vectorinfo indInfo(const string & name, UINT subPop)
	{
		UINT idx = infoIdx(name);

		return vectorinfo(infoBegin(idx, subPop), infoEnd(idx, subPop));
	}


	///	add an information field to a population
	/**
	 * <group>8-info</group>
	 */
	void addInfoField(const string field, double init = 0);

	/// add one or more information fields to a population
	/**
	   fields an array of new information fields. If one or more of the fields
	   alreay exist, they will be re-initialized.
	   init initial value for the new fields.
	 * <group>8-info</group>
	 */
	void addInfoFields(const vectorstr & fields, double init = 0);

	/// set information fields for an existing population. The existing fields will be removed.
	/**
	   fields an array of fields
	   init initial value for the new fields.
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
	 *  should always be called to restore the correct order of ancestral
	 *  generations.
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
	 */
	void sortIndividuals(bool infoOnly = false);

	/** Save population to a file \e filename. The population can be restored
	 *  from this file, using a global function <tt>LoadPopulation(filename)</tt>.
	 *  <group>1-pop</group>
	 */
	void save(const string & filename) const;

	/// CPPONLY load population from a file
	/**
	   filename load from filename

	   <group>1-pop</group>
	 */
	void load(const string & filename);

private:
	population & newPopByIndIDPerGen(const vectori & id = vectori(),
		bool removeEmptySubPops = false);

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


	/** return variables of a population. If \c subPop is given, return a dictionary for specified subpopulation.
	    <group>9-var</group>
	 */
	PyObject * vars(int subPop = -1);

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


	/// CPPONLY rearrange loci on chromosomes, e.g. combine two chromosomes into one
	/**
	   This is used by \c mergeByLoci.
	 */
	void rearrangeLoci(const vectoru & newNumLoci, const vectorf & newLociPos);

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

		if (version <= 1) {
			ar & m_genotype;
		}
		// the new version
		// the binary genotypes are saved in an efficient way
		else {
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
		}                                                                                 // verion >= 2

		if (version > 0) {
			DBG_DO(DBG_POPULATION, cout << "Handling info" << endl);
			ar & m_info;
		}
		DBG_DO(DBG_POPULATION, cout << "Handling individuals" << endl);
		ar & m_inds;

		// set genostructure, check duplication
		// we can not use setGenoStruIdx since stru may be new.
		this->setGenoStructure(stru);

		m_numSubPop = m_subPopSize.size();
		m_popSize = accumulate(m_subPopSize.begin(), m_subPopSize.end(), 0L);

		DBG_FAILIF(m_info.size() != m_popSize * infoSize(), ValueError, "Wgong size of info vector");

		if (m_popSize != m_inds.size() ) {
			cout << "Number of individuals loaded" << m_inds.size() << endl;
			cout << "population size" << m_popSize << endl;
			throw ValueError("Number of individuals does not match population size.\n"
				             "Please use the same (binary, short or long) module to save and load files.");
		}

		DBG_DO(DBG_POPULATION, cout << "Reconstruct individual genotype" << endl);
		m_subPopIndex.resize(m_numSubPop + 1);
		UINT i = 1;
		for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
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
			// version <= 1, direct handling
			if (version <= 1) {
				ar & pd.m_genotype;
			} else {
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
			}
			if (version > 0)
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

	/// number of subpopulations
	UINT m_numSubPop;

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
	bool m_indOrdered;

	/// selection flags for each subpopulation.
	/// empty means no selection
	vector<bool> m_selectionFlags;

};

/** load a population from a file.
 */
population & LoadPopulation(const string & file);

/// get info through ind.info()
vectorf testGetinfoFromInd(population & pop);

/// get info through IndInfoIterator
vectorf testGetinfoFromPop(population & pop, bool order);

}


#ifndef SWIG
#  ifndef _NO_SERIALIZATION_
// version 0: base
// version 1: save info
// version 2: reduce binary file size
BOOST_CLASS_VERSION(simuPOP::population, 2)
#  endif
#endif
#endif
