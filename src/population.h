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
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/version.hpp>
using boost::serialization::make_nvp;

#include "individual.h"
#include "virtualSubPop.h"

#ifdef SIMUMPI

#  define SLAVE_CREATE 1
#  define SLAVE_DESTROY 2
#endif

namespace simuPOP {

class population;

//************Documentation Format*****************
//   /// brief description
//   /**
//   details
//   \param
//   \note
//   \return
//   \sa
//   etc..
//   */
//

/// A collection of individuals with the same genotypic structure.
/**
   A simuPOP population consists of individuals of the same genotypic structure,
   which refers to the number of chromosomes, numbers and positions of loci on each
   chromosome etc. The most important components of a population are:

 \li subpopulations. A population is divided into subpopulations (unstructured
   population has a single	subpopulation, which is the whole population itself).
   Subpopulation structure limits the usually random exchange of genotypes between
   individuals by disallowing mating between individuals from different subpopulations.
   In the presence of subpopualtion structure, exchange of genetic information
   across subpopulations can only be done through migration. Note that simuPOP uses
   one-level population structure, which means there is no sub-subpopulation or
   family in subpopulations.

 \li variables. Every population has its own variable space, or <em>local
   namespace</em> in simuPOP term. This namespace is a Python dictionary that is
   attached to each population and can be exposed to the users through \c vars()
   or \c dvars() function. Many functions and operators work and store their results
   in this namespace. For example, function \c Stat sets variables
   such as <tt>alleleFreq[loc]</tt>, and you can access it via
   <tt>pop.dvars().alleleFreq[loc][allele]</tt>.

 \li ancestral generations. A population can save arbitrary number of ancestral
   generations. During evolution, the latest several (or all) ancestral generations
   are saved. Functions to switch between ancestral generations are
   provided so that one can examine and modify ancestral generations.
 */
class population : public GenoStruTrait
{
public:
	/** @name  constructors and destructor */
	//@{

	///Create a population object with given size and genotypic structure.
	/**
	 \param size population size. Can be ignored if \c subPop is specified.
	   	In that case, \c size is the sum of \c subPop. Default to \c 0.
	 \param ploidy number of sets of homologous copies of chromosomes. Default to \c 2 (diploid).
	 \param loci an array of numbers of loci on each chromosome. The length
	   	of parameter \c loci determines the number of chromosomes. Default
	   	to <tt>[1]</tt>, meaning one chromosome with a single locus. \n
	   	The last chromosome can be sex chromosome. In this case, the maximum
	   	number of loci on X and Y should be provided. I.e., if there are 3
	   loci on Y chromosme and 5 on X chromosome, use \c 5.
	 \param sexChrom Diploid population only. If this parameter is \c True,
	   the last homologous chromosome will be treated as sex chromosome.
	   (XY for male and XX for female.) If X and Y have different numbers of loci,
	   the number of loci of the longer one of the last (sex) chromosome should be
	   specified in \c loci.
	 \param lociPos a 1-d or 2-d array specifying positions of loci on each
	   chromosome. You can use a nested array to specify loci position for
	   each chromosome. For example, you can use <tt>lociPos=[1,2,3]</tt>
	   when <tt>loci=[3]</tt> or <tt>lociPos=[[1,2],[1.5,3,5]]</tt> for
	   <tt>loci=[2,3]</tt>. simuPOP does not assume a unit for these
	   positions, although they are usually intepreted as centiMorgans.
	   The default values are \c 1, \c 2, etc. on each chromosome.
	 \param subPop an array of subpopulation sizes. Default value is
	   <tt>[size]</tt>	which means a single subpopulation of the whole
	   population. If both	\c size and \c subPop are provided,
	 \c subPop should add up to \c size.
	 \param ancestralDepth number of most recent ancestral generations to keep
	   during evolution. Default to \c 0, which means only the current generation
	   will be available. You can set it to a positive number \c m to
	   keep the latest m generations in the population, or \c -1 to keep all ancestral
	   populations. Note that keeping track of all ancestral generations may quickly
	   exhaust your computer RAM. If you really need to do that, using \c savePopulation
	   operator to save each generation to a file is a much better choice.
	 \param chromNames an array of chromosome names.
	 \param alleleNames an array of allele names. For example, for a locus with alleles
	   A, C, T, G, you can specify \c alleleNames as <tt>('A','C','T','G')</tt>.
	 \param lociNames an array or a matrix (separated by chromosomes) of names for
	   each locus. Default to \c "locX-Y" where \c X is the chromosome index and \c Y
	   is the locus number, both starting from 1.
	 \param maxAllele maximum allele number. Default to the maximum allowed allele state
	   of the current library. This will set a cap for all loci. For individual locus, you can
	   specify \c maxAllele in mutation models, which can be smaller than the global \c maxAllele
	   but not larger. Note that this number is the number of allele states minus 1 since
	   allele number starts from 0.
	 \param infoFields names of information fields that will be attached to each
	   individual. For example, if you need to record the parents of each individual
	   using operator
	   <tt>parentTagger()</tt>, you will need two fields \c father_idx and \c mother_idx.
	 \param chromMap For MPI modules, currently unused.
	 \return no return value. Exception will be thrown when wrong parameters are given.
	 \sa simulator, baseOperator, mating schemes
	 \test src_population.log Population initialization and member functions
	 */
	population(ULONG size = 0,
	           UINT ploidy = 2,
	           const vectoru & loci = vectoru(),
	           bool sexChrom = false,
	           const vectorf & lociPos = vectorf(),
	           const vectorlu & subPop = vectorlu(),
	           int ancestralDepth = 0,
	           const vectorstr & chromNames = vectorstr(),
	           const vectorstr & alleleNames = vectorstr(),
	           const vectorstr & lociNames = vectorstr(),
	           UINT maxAllele = ModuleMaxAllele,
	           const vectorstr & infoFields = vectorstr(),
	           const vectori & chromMap = vectori());

	/// CPPONLY copy constructor
	population(const population & rhs);

	/// deep copy of a population. (In python, <tt>pop1 = pop</tt> will only create a reference to \c pop.)
	/**
	   This function by default copies all ancestral generations, but you can copy only one (current,
	   <tt>keepAncestralPops=0</tt>), or specified number of ancestral generations.
	 */
	population * clone(int keepAncestralPops = -1) const;

	/// swap the content of two populations
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
		std::swap(m_ancestralDepth, rhs.m_ancestralDepth);
		m_vars.swap(rhs.m_vars);
		m_ancestralPops.swap(rhs.m_ancestralPops);
		std::swap(m_rep, rhs.m_rep);
		std::swap(m_grp, rhs.m_grp);
		std::swap(m_gen, rhs.m_gen);
		std::swap(m_curAncestralGen, rhs.m_curAncestralGen);
		std::swap(m_shallowCopied, rhs.m_shallowCopied);
		std::swap(m_infoOrdered, rhs.m_infoOrdered);
	}


	/// destroy a population
	~population();


	/// used by Python print function to print out the general information of the population
	string __repr__()
	{
		return "<simuPOP::population of size " + toStr(popSize()) + ">";
	}


	/// if a population has any activated virtual subpopulations
	/// CPPONLY
	bool hasActivatedVirtualSubPop() const;

	/// if a subpopulation has any activated virtual subpopulation
	/// CPPONLY
	bool hasActivatedVirtualSubPop(SubPopID subPop) const;

	/// if a subpopulation has any virtual subpopulation
	bool hasVirtualSubPop(SubPopID subpop) const;

	/// set a virtual splitter to a given subpopulation.
	vspSplitter * setVirtualSplitter(vspSplitter * vsp, SubPopID sp);

	/// number of virtual subpopulation of a given subpopulation.
	UINT numVirtualSubPop(SubPopID subPop) const;

	/// activate a virtual subpopulation. CPPONLY
	/**
	 \param id subpopulation id
	 \param vid virtual subpopulation id
	 */
	void activateVirtualSubPop(SubPopID subPop, SubPopID virtualSubPop = InvalidSubPopID,
		vspSplitter::activateType type = vspSplitter::Visible);

	/// CPPONLY
	/// Copy virtual splitters if two populations have
	/// the same number of subpopulations.
	void copyVirtualSplitters(const population & rhs);

	/// deactivate virtual subpopulations in a given
	/// subpopulation. In another word, all individuals
	/// will become visible.
	/// CPPONLY
	void deactivateVirtualSubPop(SubPopID subPop);

	// allow compaison of populations in python
	// only equal or unequal, no greater or less than
	/// a python function used to compare the population objects
	int __cmp__(const population & rhs) const;

	/// set population/subpopulation structure given subpopulation sizes
	/**
	 \param newSubPopSizes an array of new subpopulation sizes.
	   The population may or may not change according to
	   parameter \c allowPopSizeChange if the sum of \c subPopSize
	   does not match \c popSize.
	 \param allowPopSizeChange if this parameter is \c True, population will be resized.
	 \return none
	 \sa migration, mating
	 */
	void setSubPopStru(const vectorlu & newSubPopSizes, bool allowPopSizeChange = false);

	///  number of subpopulations in a population.
	/**
	 \return number of subpopulations (>=1)
	 */
	UINT numSubPop() const
	{
		return m_numSubPop;
	}


	/// return size of a subpopulation \c subPop.
	/**
	 \param subPop index of subpopulation (start from 0)

	 \return size of subpopulation \c subPop
	 */
	ULONG subPopSize(SubPopID subPop) const
	{
		CHECKRANGESUBPOP(subPop);
		return m_subPopSize[subPop];
	}


	/**
	   return the size of virtual subpopulation subPop.
	   if subPop is activated, and subPop does not specify
	   which virtual subpopulation to count, the currently
	   activated virtual subpop is returned. Therefore,
	   When it is not certain if a subpopulation has activated
	   virtual subpopulation, this function can be used.

	 \param id subpopulation id
	 \param vid virtual subpopulation id. If not given,
	   current subpopulation, or current actived subpopulation size
	   will be returned.
	 */
	ULONG virtualSubPopSize(SubPopID subPop, SubPopID virtualSubPop = InvalidSubPopID) const;

	/// name of the given virtual subpopulation.
	/**
	 \param id subpopulation id
	 \param vid virtual subpopulation id
	 */
	string virtualSubPopName(SubPopID subPop, SubPopID virtualSubPop = InvalidSubPopID) const;

	/// return an array of all subpopulation sizes.
	/**
	 \return an array of size of subpopulations
	 */
	vectorlu subPopSizes() const
	{
		return m_subPopSize;
	}


	//@}

	/** @name indexes, conversion between absoluate indexes and relative indexes. return of chromomosome/subpopulation indexes.
	 */
	//@{

	/// total population size
	/**
	 \return total number of individuals in this population
	 */
	ULONG popSize() const
	{
		return m_popSize;
	}


	///  return the absolute index of an individual in a subpopulation.
	/**
	 \param index index of an individual in a subpopulation \c subPop
	 \param subPop subpopulation index (start from \c 0)
	 \return the absolute index of an individual in a subpopulation
	 \sa subPopIndPair
	 */
	ULONG absIndIndex(ULONG ind, UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);
		CHECKRANGESUBPOPMEMBER(ind, subPop);

		return m_subPopIndex[subPop] + ind;
	}


	/// return the subpopulation ID and relative index of an individual with absolute index \c ind
	/*
	 \param absInd absolute index of an individual
	 \return a pair of values (subPop, index)
	 \sa absIndIndex
	 */
	std::pair<UINT, ULONG> subPopIndPair(ULONG ind)
	{
		CHECKRANGEIND(ind);

		pair<UINT, ULONG> loc;

		for (UINT i = 1; i <= m_numSubPop; ++i) {
			if (m_subPopIndex[i] > ind) {
				loc.first = i - 1;
				loc.second = ind - m_subPopIndex[i - 1];
				break;
			}
		}
		return loc;
	}


	/// index of the first individual of a subpopulation \c subPop
	/**
	 \result beginning index of this subpopulation
	 \sa absIndIndex
	 */
	ULONG subPopBegin(UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return m_subPopIndex[subPop];
	}


	/// return the value of the index of the last individual of a subpopulation \c subPop plus 1
	/**
	 \return ending index of this subpopulation (not in this subpop)
	 \sa absIndIndex
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

	/// refernce to individual \c ind in subpopulation \c subPop
	/**
	   This function is named \c individual in the Python interface.
	 \param ind individual index within \c subPop
	 \param subPop subpopulation index
	 \return reference to an individual
	 */
	individual & ind(ULONG ind, UINT subPop = 0)
	{
#ifndef OPTIMIZED
		if (subPop > 0) {
			CHECKRANGESUBPOPMEMBER(ind, subPop);
		} else {
			CHECKRANGEIND(ind);
		}
#endif

		return m_inds[ subPopBegin(subPop) + ind];
	}


	/// return an iterator that can be used to iterate through all individuals
	/**
	   Typical usage is \n <tt>for ind in pop.individuals():</tt>
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


	/// return an iterator that can be used to iterate through all individuals in subpopulation \c subPop
	pyIndIterator individuals(SubPopID subPop)
	{
#ifndef OPTIMIZED
		CHECKRANGESUBPOP(subPop);
#endif
		// if a virtual subpopulation is activated, this will
		// iterate through virtual subpopulation. However,
		// users are not supposed to manually activate subpopulation
		// so this feature is CPPONLY
		return pyIndIterator(m_inds.begin() + subPopBegin(subPop),
		           m_inds.begin() + subPopEnd(subPop),
				   // if there is no activated virtual subpopualtions
				   // iterate through all individuals.
		           !hasActivatedVirtualSubPop(subPop),
				   // otherwise, iterate through all visible individuals.
				   true);
	}

	pyIndIterator individuals(SubPopID subPop, SubPopID virtualSubPop)
	{
#ifndef OPTIMIZED
		CHECKRANGESUBPOP(subPop);
#endif
		DBG_FAILIF(hasVirtualSubPop(subPop), ValueError,
			"Population does not have any virtual subpopulation");
		
		// this does not need to be deactivated...
		activateVirtualSubPop(subPop, virtualSubPop, vspSplitter::Iteratable);
		
		// if there is no virtual subpop
		return pyIndIterator(m_inds.begin() + subPopBegin(subPop),
		           m_inds.begin() + subPopEnd(subPop),
				   // allInds will not work at all, because there will be
				   // virtual subpopulation
		           false,
				   // and we count visible, and iteratable individuals.
				   false);
	}
	
	/// CPPONLY refernce to individual \c ind in subpopulation \c subPop
	/**
	   Return individual \ind from subpopulation \subPop. This function
	   is named \c individual in the Python interface.
	 */
	const individual & ind(ULONG ind, UINT subPop = 0) const
	{
#ifndef OPTIMIZED
		if (subPop > 0) {
			CHECKRANGESUBPOPMEMBER(ind, subPop);
		} else {
			CHECKRANGEIND(ind);
		}
#endif

		return m_inds[ subPopBegin(subPop) + ind];
	}


	/// CPPONLY
	bool shallowCopied()
	{
		return m_shallowCopied;
	}


	/// CPPONLY
	void setShallowCopied(bool s)
	{
		m_shallowCopied = s;
	}


	/// CPPONLY
	bool infoOrdered() const
	{
		return m_infoOrdered;
	}


	/// CPPONLY
	void setInfoOrdered(bool s)
	{
		m_infoOrdered = s;
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


	/// CPPONLY individual iterator: with subPop info.
	/// The iterator will skip invisible individuals
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


	/// CPPONLY individual iterator: without subPop info
	/// The iterator will skip invisible individuals
	ConstIndIterator indBegin() const
	{
		return ConstIndIterator(m_inds.begin(), m_inds.end(),
		           !hasActivatedVirtualSubPop());
	}


	/// CPPONLY individual iterator: without subPop info
	/// It is recommended to use it.valid(), instead of it != indEnd()
	ConstIndIterator indEnd() const
	{
		return ConstIndIterator(m_inds.end(), m_inds.end(),
		           !hasActivatedVirtualSubPop());
	}


	/// CPPONLY individual iterator: with subPop info.
	/// The iterator will skip invisible individuals
	ConstIndIterator indBegin(UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return ConstIndIterator(m_inds.begin() + m_subPopIndex[subPop],
		           m_inds.begin() + m_subPopIndex[subPop + 1],
		           !hasActivatedVirtualSubPop(subPop));
	}


	/// CPPONLY individual iterator: with subPop info.
	/// It is recommended to use it.valid(), instead of it != indEnd(sp)
	ConstIndIterator indEnd(UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return ConstIndIterator(m_inds.begin() + m_subPopIndex[subPop + 1],
		           m_inds.begin() + m_subPopIndex[subPop + 1],
		           !hasActivatedVirtualSubPop(subPop));
	}


	/// CPPONLY individual iterator: without subPop info
	RawIndIterator rawIndBegin()
	{
		return m_inds.begin();
	}


	/// CPPONLY individual iterator: without subPop info
	RawIndIterator rawIndEnd()
	{
		return m_inds.end();
	}


	/// CPPONLY individual iterator: with subPop info.
	/// The iterator will skip invisible individuals
	RawIndIterator rawIndBegin(UINT subPop)
	{
		CHECKRANGESUBPOP(subPop);

		return m_inds.begin() + m_subPopIndex[subPop];
	}


	/// CPPONLY individual iterator: with subPop info.
	RawIndIterator rawIndEnd(UINT subPop)
	{
		CHECKRANGESUBPOP(subPop);

		return m_inds.begin() + m_subPopIndex[subPop + 1];
	}


	/// CPPONLY individual iterator: without subPop info
	/// The iterator will skip invisible individuals
	ConstRawIndIterator rawIndBegin() const
	{
		return m_inds.begin();
	}


	/// CPPONLY individual iterator: without subPop info
	/// It is recommended to use it.valid(), instead of it != indEnd()
	ConstRawIndIterator rawIndEnd() const
	{
		return m_inds.end();
	}


	/// CPPONLY individual iterator: with subPop info.
	/// The iterator will skip invisible individuals
	ConstRawIndIterator rawIndBegin(UINT subPop) const
	{
		CHECKRANGESUBPOP(subPop);

		return m_inds.begin() + m_subPopIndex[subPop];
	}


	/// CPPONLY individual iterator: with subPop info.
	/// It is recommended to use it.valid(), instead of it != indEnd(sp)
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
	IndAlleleIterator alleleBegin(UINT locus, bool order)
	{
		CHECKRANGEABSLOCUS(locus);

		// if there is virtual subpop, use individual based iterator
		// or
		// if requires order, but the alleles are not ordered
		// use individual based
		if (hasActivatedVirtualSubPop() || (order && shallowCopied()))
			return IndAlleleIterator(locus, indBegin(), ploidy(), totNumLoci());
		else
			return IndAlleleIterator(m_genotype.begin() + locus, totNumLoci());
	}


	/// CPPONLY allele iterator
	IndAlleleIterator alleleEnd(UINT locus, bool order)
	{
		CHECKRANGEABSLOCUS(locus);
		if (hasActivatedVirtualSubPop() || (order && shallowCopied()))
			return IndAlleleIterator(locus, indEnd(), ploidy(), totNumLoci());
		else
			return IndAlleleIterator(m_genotype.begin() + locus + m_popSize * genoSize(), totNumLoci());
	}


	/// CPPONLY allele begin, for given subPop
	/**
	   order = True: keep order
	   order = false: repect subpop
	 */
	IndAlleleIterator alleleBegin(UINT locus, UINT subPop, bool order)
	{
		CHECKRANGEABSLOCUS(locus);
		CHECKRANGESUBPOP(subPop);

		if (hasActivatedVirtualSubPop() || (order && shallowCopied()))
			return IndAlleleIterator(locus, indBegin(subPop), ploidy(), totNumLoci());
		else
			return IndAlleleIterator(m_genotype.begin() + m_subPopIndex[subPop] * genoSize() +
			           locus, totNumLoci());
	}


	///  CPPONLY allele iterator
	IndAlleleIterator alleleEnd(UINT locus, UINT subPop, bool order)
	{
		CHECKRANGEABSLOCUS(locus);
		CHECKRANGESUBPOP(subPop);

		if (hasActivatedVirtualSubPop() || (order && shallowCopied()))
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
		if (order && shallowCopied())
			adjustGenoPosition(true);

		return m_genotype.begin();
	}


	///  CPPONLY allele iterator
	GenoIterator genoEnd(bool order)
	{
		if (order && shallowCopied())
			adjustGenoPosition(true);

		return m_genotype.end();
	}


	///  CPPONLY allele iterator, go through all allels one by one in a subpopulation
	/**
	   if order, keep order
	   if not order, respect subpopulation structure
	 */
	GenoIterator genoBegin(UINT subPop, bool order)
	{
		CHECKRANGESUBPOP(subPop);

		if (shallowCopied())
			adjustGenoPosition(order);

		return m_genotype.begin() + m_subPopIndex[subPop] * genoSize();
	}


	/// CPPONLY allele iterator in a subpopulation.
	GenoIterator genoEnd(UINT subPop, bool order)
	{
		CHECKRANGESUBPOP(subPop);
		if (shallowCopied())
			adjustGenoPosition(order);

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


	/// get the whole genotypes
	/**
	   Return an editable array of all genotypes of the population. You need to
	   know how these genotypes are organized to safely read/write genotype
	   directly.
	 \param order if order is \c true, individuals will be ordered such that
	   	<tt>pop.individual(x).arrGenotype() == pop.arrGenotype()[x*pop.genoSize():(x+1)*pop.genoSize()]</tt>.
	 */
	PyObject * arrGenotype(bool order);

	/// get the whole genotypes of individuals in a subpopulation
	/**
	   Return an editable array of all genotype in a subpopulation.
	 \param subPop index of subpopulation (start from 0)
	 \param order if order is \c true, individuals will be ordered.
	 */
	PyObject * arrGenotype(UINT subPop, bool order);

	//@}

	/** @name utility functions, set subpopulation, save and load etc.
	 */
	//@{

	/// set subpopulation ID with given ID
	/**
	   Set subpopulation ID of each individual with given ID. Individuals
	   can be rearranged afterwards using \c setSubPopByIndID.

	 \param id an array of the same length of population size, resprenting
	   subpopulation ID of each individual.
	 \sa individual::setSubPopID, individual::subPopID
	 */
	void setIndSubPopID(const vectori & id)
	{
		DBG_ASSERT(id.size() == m_popSize, ValueError,
		    "Info should have the same length as pop size");

		for (ULONG it = 0; it < m_popSize; ++it)
			ind(it).setSubPopID(static_cast<SubPopID>(id[it]));
	}


	/// set subpopulation ID of each individual with their current subpopulation ID
	void setIndSubPopIDWithID()
	{
		for (UINT i = 0, iEnd = numSubPop(); i < iEnd;  ++i)
			for (IndIterator it = indBegin(i), itEnd = indEnd(i); it != itEnd;  ++it)
				it->setSubPopID(i);
	}


	/// move individuals to subpopulations according to individual subpopulation IDs
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
	 */
	void splitSubPop(UINT which, vectorlu sizes, vectoru subPopID = vectoru());

	/// split a subpopulation into subpopulations of given proportions
	/**
	   The sum of given proportions should add up to one. Subpopulation IDs can be specified.
	 \note \c subpop with negative ID will be removed. So, you can shrink one \c subpop by splitting and
	   setting one of the new \c subpop with negative ID.
	 */
	void splitSubPopByProportion(UINT which, vectorf proportions, vectoru subPopID = vectoru());

	/// remove empty subpopulations by adjusting subpopulation IDs
	void removeEmptySubPops();

	/// remove subpopulations and adjust subpopulation IDs so that there will be no \em 'empty' subpopulation left
	/**
	   Remove specified subpopulations (and all individuals within). If \c shiftSubPopID is false, \c subPopID
	   will be kept intactly.
	 */
	void removeSubPops(const vectoru & subPops = vectoru(), bool shiftSubPopID = true, bool removeEmptySubPops = false);

	/// remove individuals. If a valid \c subPop is given, remove individuals from this subpopulation. Indexes in \c inds will be treated as relative indexes.
	void removeIndividuals(const vectoru & inds = vectoru(), int subPop = -1, bool removeEmptySubPops = false);

	/// merge given subpopulations
	/**
	   Merge subpopulations, the first subpopulation ID (the first one in array
	 \c subPops) will be used as the ID of the new subpopulation. That is to
	   say, all merged subpopulations will take the ID of the first one. The
	   subpopulation ID of the empty subpopulations will be kept (so that other
	   subpopulations are unaffected, unless they are removed by <tt>removeEmptySubPops = True</tt>).
	 */
	void mergeSubPops(vectoru subPops = vectoru(), bool removeEmptySubPops = false);

	/// merge populations by individuals
	/**
	   Merge individuals from \c pop to the current population.
	   Two populations should have the same genotypic structures. By default, subpopulations
	   of the merged populations are kept. I.e., if you merge two populations with one
	   subpopulation, the resulting population will have two subpopulations. All ancestral
	   generations are also merged by default.
	 \param newSubPopSizes subpopulation sizes can be specified. The overall size should
	   be the combined size of the two populations. Because this parameter will
	   be used for all ancestral generations, it may fail if ancestral generations have
	   different sizes. To avoid this problem, you can run \c mergePopulation without this parameter,
	   and then adjust subpopulation sizes generation by generation.
	 \param keepAncestralPops ancestral populations to merge, default to all (\c -1)
	 \note Population variables are not copied to \c pop.
	 */
	void mergePopulation(const population & pop, const vectorlu & newSubPopSizes = vectorlu(),
	                     int keepAncestralPops = -1);

	/// merge populations by loci
	/**
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

	/// insert loci before given positions
	/** Insert loci at some given locations. Alleles at inserted loci are initialized with zero allele.
	 \param idx an array of locus index. The loci will be inserted \em before each index.
	   If you need to append to the last locus, use \c insertAfterLoci instead. If your index is the first
	   locus of a chromosome, the inserted locus will become the first locus of that chromosome.
	   If you need to insert multiple loci before a locus, repeat that locus number.
	 \param pos an array of locus positions. The positions of the appended loci have to be between adjacent markers.
	 \param names an array of locus names. If this parameter is not given, some unique names such as
	   "insX_Y" will be given.
	 */
	void insertBeforeLoci(const vectoru & idx, const vectorf & pos, const vectorstr & names = vectorstr());

	/// insert an locus before a given position
	/**
	   <tt>insertBeforeLocus(idx, pos, name)</tt> is a shortcut to <tt>insertBeforeLoci([idx], [pos], [name])</tt>
	 */
	void insertBeforeLocus(UINT idx, double pos, const string & name = string());

	/// append loci after given positions
	/**
	   Append loci at some given locations. Alleles at inserted loci are initialized with zero allele.
	 \param idx an array of locus index. The loci will be added \em after each index.
	   If you need to append to the first locus of a chromosome, use \c insertBeforeLoci instead.
	   If your index is the last locus of a chromosome, the appended locus will become the last locus
	   of that chromosome. If you need to append multiple loci after a locus, repeat that locus number.
	 \param pos an array of locus positions. The positions of the appended loci have to be between adjacent markers.
	 \param names an array of locus names. If this parameter is not given, some unique names such as
	   "insX_Y" will be given.
	 */
	void insertAfterLoci(const vectoru & idx, const vectorf & pos, const vectorstr & names = vectorstr());

	/// append an locus after a given position
	/**
	   <tt>insertAfterLocus(idx, pos, name)</tt> is a shortcut to <tt>insertAfterLoci([idx], [pos], [name])</tt>.
	 */
	void insertAfterLocus(UINT idx, double pos, const string & name = string());

	/// resize current population
	/**
	   Resize population by giving new subpopulation sizes.
	 \param newSubPopSizes an array of new subpopulation sizes. If there
	   is only one subpopulation, use <tt>[newPopSize]</tt>.
	 \param propagate if \c propagate is \c true, copy individuals to new comers.
	   I.e., 1, 2, 3 ==> 1, 2, 3, 1, 2, 3, 1
	 \note This function only resizes the current generation.
	 */
	void resize(const vectorlu & newSubPopSizes, bool propagate = false);

	/// reorder subpopulations by \c order or by \c rank
	/**
	 \param order new order of the subpopulations. For examples, 3 2 0 1
	   means \c subpop3, \c subpop2, \c subpop0, \c subpop1 will be the new layout.
	 \param rank you may also specify a new rank for each subpopulation. For example, 3,2,0,1
	   means the original subpopulations will have new IDs 3,2,0,1, respectively. To achive order 3,2,0,1,
	   the rank should be 1 0 2 3.
	 */
	void reorderSubPops(const vectoru & order = vectoru(), const vectoru & rank = vectoru(),
	                    bool removeEmptySubPops = false);

	/**
	   Form a new population according to individual subpopulation ID. Individuals with negative subpopulation
	   ID will be removed.
	 */
	population & newPopByIndID(int keepAncestralPops = -1,
	                           const vectori & id = vectori(),
	                           bool removeEmptySubPops = false);

	/// remove some loci from the current population. Only one of the two parameters can be specified.
	void removeLoci(const vectoru & remove = vectoru(), const vectoru & keep = vectoru());

	/// obtain a new population with selected loci
	/**
	   Copy current population to a new one with selected loci \c keep or remove specified loci \c remove
	   (no change on the current population), equivalent to \n <tt>y=x.clone</tt> \n <tt>y.removeLoci(remove, keep)</tt>
	 */
	population & newPopWithPartialLoci(
	                                   const vectoru & remove = vectoru(),
	                                   const vectoru & keep = vectoru());

	/// absorb \c rhs population as the current generation of a population
	/**
	   This function is used by a simulator to push offspring generation \c rhs
	   to the current population, while the current population is pushed
	   back as an ancestral population (if <tt>ancestralDepath() != 0</tt>). Because \c rhs
	   population is swapped in, \c rhs will be empty after this operation.
	 */
	void pushAndDiscard(population & rhs, bool force = false);

	/// ancestral depth of the current population
	/** \note The return value is the number of ancestral generations
	   	exist in the population, not necessarily equals to the number set by
	 \c setAncestralDepth().
	 */
	UINT ancestralDepth() const
	{
		return m_ancestralPops.size();
	}


	/// currently used ancestral population (\c 0 for the latest generation)
	/**
	   Current ancestral population activated by \c useAncestralPop(). There can be
	   several ancestral generations in a population. \c 0 (current), \c 1 (parental)
	   etc. When \c useAncestralPop(gen) is used, current generation is set to
	   one of the parental generations, which is the information returned by this
	   function. \c useAncestralPop(0) should always be used to set a population
	   to its usual ancestral order after operations to the ancestral generation are done.
	 */
	int ancestralGen() const
	{
		return m_curAncestralGen;
	}


	// int requestInfoField(const string name);

	/// set individual information for the given information field \c idx
	/**
	 \param values an array that has the same length as population size.
	 \param idx index to the information field.
	 \param order if true, info will be in the order of individuals
	 */
	template<typename T, typename T1>
	void setIndInfo(const T & values, T1 idx, bool order = true)
	{
		CHECKRANGEINFO(idx);
		DBG_ASSERT(values.size() == popSize(), IndexError,
		    "Size of values should be the same as population size");
		UINT is = infoSize();
		if (order && !infoOrdered())
			adjustInfoPosition(true);
		typename T::const_iterator infoIter = values.begin();
		for (vectorinfo::iterator ptr = m_info.begin() + idx;
		     ptr != m_info.end() + idx; ptr += is) {
			*ptr = static_cast<InfoType>(*infoIter++);
		}
	}


	/// set individual information for the given information field \c name
	/**
	   <tt>x.setIndInfo(values, name)</tt> is
	   equivalent to the \c idx version <tt>x.setIndInfo(values, x.infoIdx(name))</tt>.
	 */
	template<class T>
	void setIndInfo(const T & values, const string & name, bool order = true)
	{
		// for mpi version , use gloal idx
		int idx = infoIdx(name);

		setIndInfo<T, UINT>(values, idx, order);
	}


	/// CPPONLY info iterator
	IndInfoIterator infoBegin(UINT idx, bool order)
	{
		CHECKRANGEINFO(idx);

		// if there is virtual subpop, use individual based iterator
		// or
		// if requires order, but the information is not ordered
		// use individual based
		if (hasActivatedVirtualSubPop() || (order && !infoOrdered()))
			return IndInfoIterator(idx, indBegin());
		else
			// if not required order, or if the information is ordered
			return IndInfoIterator(idx, m_info.begin() + idx, infoSize());
	}


	/// CPPONLY
	IndInfoIterator infoEnd(UINT idx, bool order)
	{
		CHECKRANGEINFO(idx);
		if (hasActivatedVirtualSubPop() || (order && !infoOrdered()))
			return IndInfoIterator(idx, indEnd());
		else
			return IndInfoIterator(idx, m_info.begin() + idx + m_info.size(), infoSize());
	}


	/// CPPONLY info iterator
	IndInfoIterator infoBegin(UINT idx, UINT subPop, bool order)
	{
		CHECKRANGEINFO(idx);
		CHECKRANGESUBPOP(subPop);
		if (hasActivatedVirtualSubPop(subPop) || (order && !infoOrdered()))
			return IndInfoIterator(idx, indBegin(subPop));
		else
			return IndInfoIterator(idx, m_info.begin() + idx + m_subPopIndex[subPop] * infoSize(), infoSize());
	}


	/// CPPONLY
	IndInfoIterator infoEnd(UINT idx, UINT subPop, bool order)
	{
		CHECKRANGEINFO(idx);
		CHECKRANGESUBPOP(subPop);

		if (hasActivatedVirtualSubPop(subPop) || (order && !infoOrdered()))
			return IndInfoIterator(idx, indEnd(subPop));
		else
			return IndInfoIterator(idx, m_info.begin() + idx + m_subPopIndex[subPop + 1] * infoSize(), infoSize());
	}


	/// get information field \c idx of all individuals
	/**
	 \param idx index of the information field
	 \param order if true, sort returned vector in individual order
	 \return a vector with value of the information field
	 */
	vectorinfo indInfo(UINT idx, bool order)
	{
		return vectorinfo(infoBegin(idx, order), infoEnd(idx, order));
	}


	/// get information field \c name of all individuals
	/**
	 \param name name of the information field
	 \param order if true, sort returned vector in individual order
	 \return a vector with value of the information field
	 */
	vectorinfo indInfo(const string & name, bool order)
	{
		UINT idx = infoIdx(name);

		return vectorinfo(infoBegin(idx, order), infoEnd(idx, order));
	}


	/// get information field \c idx of all individuals in a subpopulation \c subPop
	/**
	 \param idx index of the information field
	 \param subPop subpopulation index
	 \param order if true, sort returned vector in individual order
	 \return a vector with value of the information field
	 */
	vectorinfo indInfo(UINT idx, UINT subPop, bool order)
	{
		return vectorinfo(infoBegin(idx, subPop, order),
		           infoEnd(idx, subPop, order));
	}


	/// get information field \c name of all individuals in a subpopulation \c subPop
	/**
	 \param name name of the information field
	 \param subPop subpopulation index
	 \param order if true, sort returned vector in individual order
	 \return a vector with value of the information field
	 */
	vectorinfo indInfo(const string & name, UINT subPop, bool order)
	{
		UINT idx = infoIdx(name);

		return vectorinfo(infoBegin(idx, subPop, order), infoEnd(idx, subPop, order));
	}


	/// get an editable array (Python list) of all information fields
	/**
	   The length of the array is <tt>infoSize()*popSize()</tt>.
	 \param order whether or not the list has the same order as individuals
	 \return a python list that has all information fields for all individual
	   values of information fields of the same individual are put together.
	 */
	PyObject * arrIndInfo(bool order)
	{
		if (order && !infoOrdered())
			adjustInfoPosition(true);

		return Info_Vec_As_NumArray(m_info.begin(), m_info.end());
	}


	/// get an editable array (Python list) of all information fields in \c subPop
	/**
	 \param order whether or not the list has the same order as individuals
	 \return a python list that has all information fields for all individual
	   values of information fields of the same individual are put together.
	 */
	PyObject * arrIndInfo(UINT subPop, bool order)
	{
		CHECKRANGESUBPOP(subPop);

		if (!infoOrdered())
			adjustInfoPosition(order);

		return Info_Vec_As_NumArray(m_info.begin() + m_subPopIndex[subPop] * infoSize(),
		           m_info.begin() + m_subPopIndex[subPop + 1] * infoSize());
	}


	///	add an information field to a population
	/**
	 \param field new information field. If it already exists, it will
	   be re-initialized.
	 \param init initial value for the new field.
	 */
	void addInfoField(const string field, double init = 0);

	/// add one or more information fields to a population
	/**
	 \param fields an array of new information fields. If one or more of the fields
	   alreay exist, they will be re-initialized.
	 \param init initial value for the new fields.
	 */
	void addInfoFields(const vectorstr & fields, double init = 0);

	/// set information fields for an existing population. The existing fields will be removed.
	/**
	 \param fields an array of fields
	 \param init initial value for the new fields.
	 */
	void setInfoFields(const vectorstr & fields, double init = 0);

	/// set ancestral depth
	/**
	 \param depth \c 0 for none, \c -1 for unlimited, a positive number sets
	   the number of ancestral generations to save.
	 */
	void setAncestralDepth(int depth);

	// idx = 0 (current), 1 (parents), 2 (grandparents...)
	/// use an ancestral generation. \c 0 for the latest generation.
	/**
	 \param idx Index of the ancestral generation. \c 0 for current,
	 \c 1 for parental, etc. idx can not exceed ancestral depth
	   (see setAncestralDepth).
	 */
	void useAncestralPop(UINT idx);

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
	   order=true: make individuals in order
	   order=false: make individuals in each subpopulation
	 */
	void adjustGenoPosition(bool order);

	/// CPPONLY
	void adjustInfoPosition(bool order);

	/// save population to a file
	/**
	 \param filename save to filename
	 \param format format to save. Can be one of the following: 'txt', 'bin', or 'xml',
	   	or 'auto' which is determined by the extension of \c filename.
	 \sa global function loadPopulation
	 */
	void savePopulation(const string & filename, const string & format = "auto", bool compress = true) const;

	/// CPPONLY load population from a file
	/**
	 \param filename load from filename
	 \param format format to load. Can be one of "text", "bin", "xml". It should match
	   the format used to save the population.
	 \sa savePopulation
	 */
	void loadPopulation(const string & filename, const string & format = "auto");

private:
	population & newPopByIndIDPerGen(const vectori & id = vectori(),
	                                 bool removeEmptySubPops = false);

	void mergePopulationPerGen(const population & pop, const vectorlu & newSubPopSizes);

public:
	/// CPPONLY
	/// selection is on at any subpopulation?
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


	/// turn off selection for all subpopulations
	/**
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


	/// CPPONLY
	/// Turn on selection for all subpopulations
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
	/// current replicate in a simulator which is not meaningful for a stand-alone population
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


	/// current group ID in a simulator which is not meaningful for a stand-alone population.
	int grp()
	{
		return m_grp;
	}


	/// CPPONLY
	void setGrp(int grp, bool setVar = true)
	{
		m_grp = grp;
		if (setVar)
			m_vars.setIntVar("grp", grp);
	}


	/// current generation during evolution
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


	/// return variables of a population. If \c subPop is given, return a dictionary for specified subpopulation.
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


	/// evaluate a Python statment/expression in the population's local namespace
	/**
	   This function evaluates a Python statment( \c stmts )/expression( \c expr )
	   and return its result as a string. Optionally run statement( \c stmts ) first.
	 */
	PyObject * evaluate(const string & expr = "", const string & stmts = "")
	{
		return Expression(expr, stmts, m_vars.dict() ).evaluate();
	}


	/// execute a statement (can be a multi-line string) in the population's local namespace
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
		//const_cast<population*>(this)->adjustGenoPosition(true);

		ar & make_nvp("libraryMaxAllele", ModuleMaxAllele);

		DBG_DO(DBG_POPULATION, cout << "Handling geno structure" << endl);
		// GenoStructure genoStru = this->genoStru();
		ar & make_nvp("geno_structure", this->genoStru());
		ar & make_nvp("subPop_sizes", m_subPopSize);
		DBG_DO(DBG_POPULATION, cout << "Handling genotype" << endl);
#ifdef BINARYALLELE
		size_t size = m_genotype.size();
		ar & make_nvp("size", size);
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
				ar & make_nvp("blocks", tmp1);
			}
		}
		// last block
		if (rest > 0) {
			tmp = *ptr;
			for (size_t j = 0; j <= (rest - 1) / 32; ++j) {
				tmp1 = tmp & 0xFFFFFFFF;
				tmp = tmp >> 32;
				ar & make_nvp("blocks", tmp1);
			}
		}
#else
		ar & make_nvp("genotype", m_genotype);
#endif
		DBG_DO(DBG_POPULATION, cout << "Handling information" << endl);
		ar & make_nvp("info", m_info);
		DBG_DO(DBG_POPULATION, cout << "Handling individuals" << endl);
		ar & make_nvp("individuals", m_inds);
		DBG_DO(DBG_POPULATION, cout << "Handling ancestral populations" << endl);
		ar & make_nvp("ancestry", m_ancestralDepth);
		size_t sz = m_ancestralPops.size();
		ar & make_nvp("numOfAncestralPops", sz);
		for (size_t i = 0; i < m_ancestralPops.size(); ++i) {
			const_cast<population *>(this)->useAncestralPop(i + 1);
			// need to make sure ancestral pop also in order
			const_cast<population *>(this)->adjustGenoPosition(true);
			ar & make_nvp("subPop_sizes", m_subPopSize);
#ifdef BINARYALLELE
			size_t size = m_genotype.size();
			ar & make_nvp("size", size);
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
					ar & make_nvp("blocks", tmp1);
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
					ar & make_nvp("blocks", tmp1);
				}
			}
#else
			ar & make_nvp("genotype", m_genotype);
#endif
			ar & make_nvp("info", m_info);
			ar & make_nvp("individuals", m_inds);
		}
		const_cast<population *>(this)->useAncestralPop(0);

		// save shared variables as string.
		// note that many format are not supported.
		try {
			DBG_DO(DBG_POPULATION, cout << "Handling shared variables" << endl);
			string vars = varsAsString();
			ar & make_nvp("vars", vars);
		} catch (...) {
			cout << "Warning: shared variable is not saved correctly.\npopulation should still be usable." << endl;
		}
	}


	template<class Archive>
	void load(Archive & ar, const UINT version)
	{

		ULONG ma;
		ar & make_nvp("libraryMaxAllele", ma);

		if (ma > ModuleMaxAllele)
			cout << "Warning: The population is saved in library with more allele states. \n"
			     << "Unless all alleles are less than " << ModuleMaxAllele
			     << ", you should use the modules used to save this file. (c.f. simuOpt.setOptions()\n";

		GenoStructure stru;
		DBG_DO(DBG_POPULATION, cout << "Handling geno structure" << endl);
		ar & make_nvp("geno_structure", stru);
		ar & make_nvp("subPop_sizes", m_subPopSize);
		DBG_DO(DBG_POPULATION, cout << "Handling genotype" << endl);

		if (version <= 1) {
			ar & make_nvp("genotype", m_genotype);
		}
		// the new version
		// the binary genotypes are saved in an efficient way
		else {
#ifdef BINARYALLELE
			// binary from binary
			if (ma == 1) {
				size_t size;
				ar & make_nvp("size", size);
				size_t blks = size / WORDBIT;
				size_t rest = size - blks * WORDBIT;

				m_genotype.resize(size);
				WORDTYPE tmp, tmp1;
				WORDTYPE * ptr = BITPTR(m_genotype.begin());
				for (size_t i = 0; i < blks; ++i) {
					tmp = 0;
					for (size_t j = 0; j < WORDBIT / 32; ++j) {
						ar & make_nvp("blocks", tmp1);
						tmp |= tmp1 << (j * 32);
					}
					*ptr++ = tmp;
				}
				// last block
				if (rest > 0) {
					tmp = 0;
					for (size_t j = 0; j <= (rest - 1) / 32; ++j) {
						ar & make_nvp("blocks", tmp1);
						tmp |= tmp1 << (j * 32);
					}
					*ptr = tmp;
				}
			}
			// binary from others (long types)
			else {
				DBG_DO(DBG_POPULATION, cout << "Load bin from long. " << endl);
				vector<unsigned char> tmpgeno;
				ar & make_nvp("genotype", tmpgeno);
				m_genotype = vectora(tmpgeno.begin(), tmpgeno.end());
			}
#else
			// long from binary
			if (ma == 1) {
				// for version 2 and higher, archive in 32bit blocks.
				size_t size;
				ar & make_nvp("size", size);
				m_genotype.resize(size);
				size_t blks = size / 32;
				size_t rest = size - blks * 32;
				DBG_DO(DBG_POPULATION, cout << "Load long from bin. " << size << " rest " << rest << endl);
				DBG_ASSERT(WORDBIT >= 32, SystemError, "WordBit should be at least 32 bits");

				GenoIterator ptr = m_genotype.begin();
				WORDTYPE tmp;
				for (size_t i = 0; i < blks; ++i) {
					ar & make_nvp("blocks", tmp);
					for (size_t j = 0; j < 32; ++j) {
						*ptr++ = (tmp & 1UL) != 0;
						tmp = tmp >> 1;
					}
				}
				// last block
				if (rest > 0) {
					ar & make_nvp("blocks", tmp);
					for (size_t j = 0; j < rest; ++j) {
						*ptr++ = (tmp & 1UL) != 0;
						tmp = tmp >> 1;
					}
				}
			}                                                                               // if ma == 1
			else {                                                                          // for non-binary types, ...
				DBG_DO(DBG_POPULATION, cout << "Load long from long. " << endl);
				// long from long
				ar & make_nvp("genotype", m_genotype);
			}
#endif
		}                                                                                 // verion >= 2

		if (version > 0) {
			DBG_DO(DBG_POPULATION, cout << "Handling info" << endl);
			ar & make_nvp("info", m_info);
		}
		DBG_DO(DBG_POPULATION, cout << "Handling individuals" << endl);
		ar & make_nvp("individuals", m_inds);

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
		m_ancestralDepth = 0;
		m_ancestralPops.clear();

		// ancestry populations
		DBG_DO(DBG_POPULATION, cout << "Handling ancestral populations" << endl);
		ar & make_nvp("ancestry", m_ancestralDepth);
		size_t na;
		ar & make_nvp("numOfAncestralPops", na);
		for (size_t ap = 0; ap < na; ++ap) {
			popData pd;
			ar & make_nvp("subPop_sizes", pd.m_subPopSize);
			// version <= 1, direct handling
			if (version <= 1) {
				ar & make_nvp("genotype", pd.m_genotype);
			} else {
#ifdef BINARYALLELE
				// binary from binary
				if (ma == 1) {
					DBG_DO(DBG_POPULATION, cout << "Load bin from bin. " << endl);
					size_t size;
					ar & make_nvp("size", size);
					size_t blks = size / WORDBIT;
					size_t rest = size - blks * WORDBIT;

					pd.m_genotype.resize(size);
					WORDTYPE * ptr = BITPTR(pd.m_genotype.begin());
					WORDTYPE tmp, tmp1;
					for (size_t i = 0; i < blks; ++i) {
						tmp = 0;
						for (size_t j = 0; j < WORDBIT / 32; ++j) {
							ar & make_nvp("blocks", tmp1);
							tmp |= tmp1 << (j * 32);
						}
						*ptr++ = tmp;
					}
					// last block
					if (rest > 0) {
						tmp = 0;
						for (size_t j = 0; j <= (rest - 1) / 32; ++j) {
							ar & make_nvp("blocks", tmp1);
							tmp |= tmp1 << (j * 32);
						}
						*ptr = tmp;
					}
				} else {
					DBG_DO(DBG_POPULATION, cout << "Load bin from long. " << endl);
					// binary from long types
					vector<unsigned char> tmpgeno;
					ar & make_nvp("genotype", tmpgeno);
					pd.m_genotype = vectora(tmpgeno.begin(), tmpgeno.end());
				}
#else
				if (ma == 1) {
					// long type from binary
					size_t size;
					ar & make_nvp("size", size);
					pd.m_genotype.resize(size);
					size_t blks = size / 32;
					size_t rest = size - blks * 32;
					DBG_DO(DBG_POPULATION, cout << "Load long from bin. " << size << " rest " << rest << endl);

					ptr = pd.m_genotype.begin();
					WORDTYPE tmp;
					for (size_t i = 0; i < blks; ++i) {
						ar & make_nvp("blocks", tmp);
						for (size_t j = 0; j < 32; ++j) {
							*ptr++ = (tmp & 1UL) != 0;
							tmp = tmp >> 1;
						}
					}
					// last block
					if (rest > 0) {
						ar & make_nvp("blocks", tmp);
						for (size_t i = 0; i < rest; ++i) {
							*ptr++ = (tmp & 1UL) != 0;
							tmp = tmp >> 1;
						}
					}
				} else {
					DBG_DO(DBG_POPULATION, cout << "Load long from long. " << endl);
					// long type from long type.
					ar & make_nvp("genotype", pd.m_genotype);
				}
#endif
			}
			if (version > 0)
				ar & make_nvp("info", pd.m_info);
			ar & make_nvp("individuals", pd.m_inds);
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
				// fresh copy so clear shallowcopied flag.
				inds[i].setShallowCopied(false);
			}
		}

		// load vars from string
		try {
			DBG_DO(DBG_POPULATION, cout << "Handling shared variables" << endl);
			string vars;
			ar & make_nvp("vars", vars);
			varsFromString(vars);
		} catch (...) {
			cout << "Warning: shared variable is not loaded correctly.\npopulation should still be usable." << endl;
		}

		m_shallowCopied = false;
		m_infoOrdered = true;
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
	vectorvsp m_virtualSubPops;

	/// pool of genotypic information
	vectora m_genotype;

	/// information
	/// only in head node
	vectorinfo m_info;

	/// individuals.
	/// only in head node?
	vector<individual> m_inds;

	int m_ancestralDepth;

	/// shared variables for this population
	SharedVariables m_vars;

	/// store previous populations
	/// need to store: subPopSize, genotype and m_inds
	struct popData
	{
#ifndef OPTIMIZED
		GenoIterator m_startingGenoPtr;
#endif
		vectorlu m_subPopSize;
		vectora m_genotype;
		vectorinfo m_info;
		vector<individual> m_inds;
	};

	std::deque<popData> m_ancestralPops;

	/// curent replicate, group number
	int m_rep, m_grp;

	/// generation
	ULONG m_gen;

	/// current ancestral depth
	int m_curAncestralGen;

	/// whether or not individual genotype is in order
	bool m_shallowCopied;

	/// whether or not information is ordered
	bool m_infoOrdered;

	/// selection flags for each subpopulation.
	/// empty means no selection
	vector<bool> m_selectionFlags;

#ifdef SIMUMPI
	/// unique population id, used by slave nodes to identify a population
	ULONG m_popID;
#endif

};

/// load a population from a file. The file format is by default determined by file extension (<tt>format="auto"</tt>). Otherwise, \c format can be one of \c txt, \c bin, or \c xml.
population & LoadPopulation(const string & file, const string & format = "auto");

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
