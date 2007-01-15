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

namespace simuPOP
{

	class population;

	/// this class implements a Python itertor class that
	/// can be used to iterate through individuals in a
	/// population.
	///
	/// an instance of this class is returned by
	/// population::individuals() and population::individuals(subPop)
	class individualIterator
	{
		public:
			individualIterator(population* pop, ULONG s, ULONG e) :
			m_population(pop),
				m_index(s),
				m_end(e)
			{
			}

			individualIterator __iter__()
			{
				return *this;
			}

			individual& next();

		private:
			// an instance of the population being iterated.
			population* m_population;

			// current (initialized as strt) index
			ULONG m_index;

			// ending index
			ULONG m_end;
	};

	/** NOTE TO THE MPI VERSION OF THE POPULATION
	 *
	 * The MPI version of simuPOP spread the chromosomes
	 * across the nodes of a cluster system.
	 *
	 * The head node (0) keeps all the individual, variables
	 * and information fields, but no genotype.
	 *
	 * Other nodes (1-...) has one chromosome (0-...) of genotype
	 * nothing else.
	 *
	 **/

	/** \brief a collection of individuals with subpopulation structure

	Please refer to user's Guide for details about this object.

	*/
	class population : public GenoStruTrait
	{
		public:

			/// individual itertor, used to iterate all individuals.
			typedef vector<individual>::iterator IndIterator;
			typedef vector<individual>::const_iterator ConstIndIterator;

			/** @name  constructors and destructor */
			//@{

			/// create a population object with given size and genotypic structure
			/**
			\param size population size. Can be ignored if subPop is specified.
			   In that case, size is sum of subPop. Default to 0.
			\param ploidy number of sets of chromosomes. Default to 2 (diploid).
			\param loci an array of numbers of loci on each chromosome. If
			   not specified, assume a single locus on one chromosome. Number
			   of chromosomes is determined by the size of this array.
			\param lociPos an array of loci distance for each locus. You can
			   also use a nested array to specify loci distance for each chromosome.
			   ( [1,2,3,4,5] or [[1,2],[3,4,5]] are both allowed for loci=[2,3])
			   The default values are 1, 2, etc. on each chromosome.
			\param subPop an array of subpopulation sizes. Default value is [size]
			which means a single subpopulation of the whole population. If both
			size and subPop are given, sum of subPop should agree with size.
			\param ancestralDepth number of ancestral populations to keep. Default to 0,
			meaning only current generation will be available. If -1 is given,
			all ancestral populations will be saved. This may exhaust your RAM
			pretty quickly though.
			\param alleleNames an array of allele names. The first element should be
			given to invalid/unknown allele. For example, for a locus with alleles
			A,C,T,G, you can specify alleleName as ('_','A','C','T','G'). Note that
			simuPOP uses 1,2,3,4 internally and these names will only be used for
			display purpose.
			\param lociNames an array or a matrix (separated by chromosome) of names for
			each loci. Default to "locX-X" where
			X-X is chromosome-loci index starting from 1. This info is rarely used.
			\param maxAllele maximum allele number. Default to the max allowed allele states
			of current library (standard or long allele version)
			\param infoFields: name of information fields that will be attached to each
			individual. For example, if you need to record the parents of each individual
			you will need two, if you need to record the age of individual, you need an additional
			one. Other possibilities include offspring ids etc. Note that you have to plan
			this ahead of time since, for example, tagger will need to know what info unit
			to use. Default to none.
			\return no return value. Exception will be thrown is wrong parameters are given.
			\sa simulator, baseOperator, mating schemes
			\test popInit.log \include popInit.log
			*/
			population( ULONG size=0,
				UINT ploidy=2,
				const vectoru& loci=vectoru(),
				bool sexChrom=false,
				const vectorf& lociPos=vectorf(),
				const vectorlu& subPop=vectorlu(),
				int ancestralDepth=0,
				const vectorstr& alleleNames=vectorstr(),
				const vectorstr& lociNames=vectorstr(),
				UINT maxAllele = MaxAllele,
				const vectorstr& infoFields = vectorstr(),
				const vectori& chromMap = vectori());

			/// CPPONLY copy constructor
			population(const population& rhs);

			///
			population * clone(int keepAncestralPops=-1) const;

			/// SWAP population
			/// swap the content of two populations
			void swap(population& rhs)
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
				std::swap(m_curAncestralPop, rhs.m_curAncestralPop);
				std::swap(m_shallowCopied, rhs.m_shallowCopied);
				std::swap(m_infoOrdered, rhs.m_infoOrdered);
			}

			/// destroy a population
			~population()
			{
				DBG_DO(DBG_POPULATION,
					cout << "Destructor of population is called" << endl);
			}

			// allow str(population) to get something better looking
			string __repr__()
			{
				return "<simuPOP::population of size " + toStr(popSize()) + ">";
			}

			// allow compaison of populations in python
			// only equal or unequal, no greater or less than
			int __cmp__(const population& rhs) const;

			/// \brief set population/subpopulation given subpopulation sizes
			///
			/// \param subPopSize an array of subpopulation sizes
			///    the population may or may not change according to
			///    parameter allowPopSizeChange if sum of subPopSize
			///    does not match popSize.
			/// \param allowPopSizeChange if true, popSize can change to sum of
			///    subPopSize.
			///
			/// \return none
			/// \sa migration, mating
			void setSubPopStru(const vectorlu& newSubPopSizes, bool allowPopSizeChange=false);

			///  number of sub populations.
			/**
			 \return number of subpopulations (>=1)
			 */
			UINT numSubPop() const
			{
				return m_numSubPop;
			}

			/// get size of subpopulation subPop
			/** \param subPop index of subpopulation (start from 0)
				\return size of subpopulation subPop
			*/
			ULONG subPopSize(UINT subPop) const
			{
				CHECKRANGESUBPOP(subPop);

				return m_subPopSize[subPop];
			}

			/// get size of all subpopulations
			/** \return an array of size of subpopulations
			 */
			vectorlu subPopSizes()
			{
				return m_subPopSize;
			}
			//@}

			/** @name indices
			  conversion between absoluate indices and relative indices.
			  return of chromomosome/subpopulation indices.
			*/
			//@{

			/// get population size
			/**
			 \return total number of individuals in this population
			 */
			ULONG popSize() const
			{
				return m_popSize;
			}

			///  absolute index of individual at a subpopulation
			/**
			  \param index index of individual at subpopulation subPop
			  \param subPop subpopulation index
			  \return absolute index of individual at subpopulation subPop
			  \sa subPopIndPair
			*/
			ULONG absIndIndex(ULONG ind, UINT subPop) const
			{
				CHECKRANGESUBPOP(subPop);
				CHECKRANGESUBPOPMEMBER(ind, subPop);

				return( m_subPopIndex[subPop] + ind);
			}

			/// subPop and relative index of an individual
			/*
			  \param absInd absolute index of an individual
			  \return a pair of values (subPop, index)
			  \sa absIndIndex
			*/
			std::pair<UINT, ULONG> subPopIndPair(ULONG ind)
			{
				CHECKRANGEIND(ind);

				pair<UINT, ULONG> loc;

				for(UINT i=1; i<=m_numSubPop; ++i)
				{
					if( m_subPopIndex[i] > ind )
					{
						loc.first = i-1;
						loc.second = ind - m_subPopIndex[i-1];
						break;
					}
				}
				return loc;
			}

			/// beginning index for subpopulation subPop
			/**
			\param subPop subpopulation index
			\result beginning index of this subpopulation
			\sa absIndIndex

			absIndIndex(index, subPop) = subPopBegin(subPop) + index
			*/
			ULONG subPopBegin(UINT subPop) const
			{
				CHECKRANGESUBPOP(subPop);

				return m_subPopIndex[subPop];
			}

			/// ending index for subpopulation subPop
			/**
			\param subPop subpopulation index
			\result ending index of this subpopulation (not in this subpop)
			\sa absIndIndex
			 \note as with all ...End functions, the returning index is out of the range
			 so the actual range is [xxxBegin, xxxEnd). This agrees with all STL
			 conventions.
			*/
			ULONG subPopEnd(UINT subPop) const
			{
				CHECKRANGESUBPOP(subPop);

				return m_subPopIndex[subPop+1];
			}

			//@}
			/** @name itertors and accessers
			  ways to access information, mainly various iterators.
			*/
			//@{

			/// refernce to individual ind in subpopulation subPop
			/** \param ind individual index within subPop
				\param subPop subpopulation index
				\return reference to an individual
			*/
			individual& ind(ULONG ind, UINT subPop=0)
			{
#ifndef OPTIMIZED
				if( subPop > 0 )
				{
					CHECKRANGESUBPOPMEMBER(ind, subPop);
				}
				else
				{
					CHECKRANGEIND(ind);
				}
#endif

				return m_inds[ subPopBegin(subPop) + ind];
			}

			individualIterator individuals()
			{
				return individualIterator(this, 0, popSize());
			}

			individualIterator individuals(UINT subPop)
			{
#ifndef OPTIMIZED
				CHECKRANGESUBPOP(subPop);
#endif

				return individualIterator(this, subPopBegin(subPop), subPopEnd(subPop));
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
			bool infoOrdered()
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
				return m_inds.begin();
			}

			/// CPPONLY individual iterator: without subPop info
			IndIterator indEnd()
			{
				return m_inds.end();
			}

			/// CPPONLY individual iterator: with subPop info.
			IndIterator indBegin(UINT subPop)
			{
				CHECKRANGESUBPOP(subPop);

				return m_inds.begin() + absIndIndex(0, subPop);
			}

			/// CPPONLY individual iterator: with subPop info.
			IndIterator indEnd(UINT subPop)
			{
				CHECKRANGESUBPOP(subPop);

				return m_inds.begin() + m_subPopIndex[subPop+1];
			}

			/// CPPONLY individual iterator: without subPop info
			ConstIndIterator indBegin() const
			{
				return m_inds.begin();
			}

			/// CPPONLY individual iterator: without subPop info
			ConstIndIterator indEnd() const
			{
				return m_inds.end();
			}

			/// CPPONLY individual iterator: with subPop info.
			ConstIndIterator indBegin(UINT subPop) const
			{
				CHECKRANGESUBPOP(subPop);

				return m_inds.begin() + absIndIndex(0, subPop);
			}

			/// CPPONLY individual iterator: with subPop info.
			ConstIndIterator indEnd(UINT subPop) const
			{
				CHECKRANGESUBPOP(subPop);

				return m_inds.begin() + m_subPopIndex[subPop+1];
			}

			/// allele iterator that access a locus across all copies of chromosomes and individual
			/** CPPONLY
			\param locus
			  allele access, given locus, return
			 the first allele. ptr++ go the next one.
			 default return the beginning of the first subpopulation,
			 also the first of the whole population

			\note The order of alleles DOES NOT HAVE TO match the order of
			individuals. Only the boundary of subpopulations will be respected.
			Therefore, it is possible to access all alleles within an subpopulation
			through such iterators.

			order = True: indiviudls in order
			order = false: do not even respect subpops
			*/
			GappedAlleleIterator alleleBegin(UINT locus, bool order)
			{
				CHECKRANGEABSLOCUS(locus);
				if(order && shallowCopied())
					adjustGenoPosition(true);
#ifdef SIMUMPI
				UINT rank = rankOfLocus(locus);
				if (mpiRank() == rank)
					return GappedAlleleIterator(m_genotype.begin() + locus - beginLocus(),
						localNumLoci());
				else
					// this iterator is invalid
					return GappedAlleleIterator(m_genotype.begin(), 0);
#else
				return GappedAlleleIterator( m_genotype.begin()+locus, totNumLoci());
#endif
			}

			/// CPPONLY allele iterator
			GappedAlleleIterator alleleEnd(UINT locus, bool order)
			{
				CHECKRANGEABSLOCUS( locus);
				if(order && shallowCopied())
					adjustGenoPosition(true);
#ifdef SIMUMPI
				UINT rank = rankOfLocus(locus);
				if (mpiRank() == rank)
					return GappedAlleleIterator(m_genotype.begin() + locus -
						beginLocus() + m_popSize*localGenoSize(), localNumLoci());
				else
					return GappedAlleleIterator(m_genotype.begin(), 0);
#else
				return GappedAlleleIterator(m_genotype.begin() + locus + m_popSize*genoSize(), totNumLoci());
#endif
			}

			///  CPPONLY allele begin, for given subPopu
			/// order = True: keep order
			/// order = false; repect subpop
			GappedAlleleIterator alleleBegin(UINT locus, UINT subPop, bool order)
			{
				CHECKRANGEABSLOCUS(locus);
				CHECKRANGESUBPOP(subPop);

				if(shallowCopied())
					adjustGenoPosition(order);
#ifdef SIMUMPI
				UINT rank = rankOfLocus(locus);
				if (mpiRank() == rank)
					return GappedAlleleIterator(m_genotype.begin() + m_subPopIndex[subPop]*localGenoSize() +
						locus - beginLocus(), localNumLoci());
				else
					return GappedAlleleIterator( m_genotype.begin(), 0);

#else
				return GappedAlleleIterator( m_genotype.begin() + m_subPopIndex[subPop]*genoSize() +
					locus, totNumLoci());
#endif
			}

			///  CPPONLY allele iterator
			GappedAlleleIterator alleleEnd( UINT locus, UINT subPop, bool order)
			{
				CHECKRANGEABSLOCUS(locus);
				CHECKRANGESUBPOP(subPop);

				if(shallowCopied())
					adjustGenoPosition(order);

#ifdef SIMUMPI
				UINT rank = rankOfLocus(locus);
				if (mpiRank() == rank)
					return GappedAlleleIterator(m_genotype.begin() + m_subPopIndex[subPop+1]*localGenoSize() +
						locus - beginLocus(), localNumLoci());
				else
					return GappedAlleleIterator(m_genotype.begin(), 0);
#else
				return GappedAlleleIterator( m_genotype.begin() + m_subPopIndex[subPop+1]*genoSize() +
					locus, totNumLoci());
#endif
			}

			///  CPPONLY allele iterator, go through all allels one by one, without subPop info
			/// if order, in order
			/// otherwise, do not even respect subpopulation structure
			GenoIterator genoBegin(bool order)
			{
				if(order && shallowCopied())
					adjustGenoPosition(true);

				return m_genotype.begin();
			}

			///  CPPONLY allele iterator
			GenoIterator genoEnd(bool order)
			{
				if(order && shallowCopied())
					adjustGenoPosition(true);

				return m_genotype.end();
			}

			///  CPPONLY allele iterator, go through all allels one by one in a subpopulation
			/// if order, keep order
			/// if not order, respect subpopulation structure
			GenoIterator genoBegin(UINT subPop, bool order)
			{
				CHECKRANGESUBPOP(subPop);

				if(shallowCopied())
					adjustGenoPosition(order);

#ifdef SIMUMPI
				return m_genotype.begin() + m_subPopIndex[subPop]*localGenoSize();
#else
				return m_genotype.begin() + m_subPopIndex[subPop]*genoSize();
#endif
			}

			///  CPPONLY allele iterator in a subpopulation.
			GenoIterator genoEnd(UINT subPop, bool order)
			{
				CHECKRANGESUBPOP(subPop);
				if(shallowCopied())
					adjustGenoPosition(order);

#ifdef SIMUMPI
				return m_genotype.begin() + m_subPopIndex[subPop+1]*localGenoSize();
#else
				return m_genotype.begin() + m_subPopIndex[subPop+1]*genoSize();
#endif
			}

			///  CPPONLY genoIterator --- beginning of individual ind.
			GenoIterator indGenoBegin(ULONG ind) const
			{
				CHECKRANGEIND(ind);
				return m_inds[ind].genoBegin();
			}

			///  CPPONLY genoIterator -- end of individual ind.
			GenoIterator indGenoEnd(ULONG ind) const
			{
				CHECKRANGEIND(ind);
				return m_inds[ind].genoEnd();
			}

			///  CPPONLY genoIterator --- beginning of individual ind.
			GenoIterator indGenoBegin(ULONG ind, UINT subPop) const
			{
				CHECKRANGESUBPOP(subPop);
				CHECKRANGESUBPOPMEMBER(ind, subPop);

				return m_inds[ subPopBegin(subPop) + ind].genoBegin();
			}

			///  CPPONLY genoIterator -- end of individual ind.
			GenoIterator indGenoEnd(ULONG ind, UINT subPop) const
			{
				CHECKRANGESUBPOP(subPop);
				CHECKRANGESUBPOPMEMBER(ind, subPop);

				return m_inds[ subPopBegin(subPop) + ind].genoEnd();
			}

			/// get the whole genotype.
			/// individuals will be in order before exposing
			/// their genotypes.
			///
			/// if order, respect order,
			/// if false, do not repect population structure
			PyObject* arrGenotype(bool order);

			/// get the whole genotype.
			/// individuals will be in order before exposing
			/// their genotypes.
			///
			/// if order: keep order
			/// otherwise: respect subpop structure
			PyObject* arrGenotype(UINT subPop, bool order);

			//@}

			/** @name utility functions
			 set subpopulation, save and load etc.
			*/
			//@{

			/** brief return individual affected status in pop namespace
			 */
			PyObject* exposeAffectedness(string name="affected")
			{
				// regardness of info type, treat it as int.
				vectori val(m_popSize);
				for(ULONG it=0; it < m_popSize; ++it)
					val[it] = static_cast<int>(ind(it).affected());
				// this has to be changed according to info type.
				PyObject* var = setIntVectorVar(name, val);
				Py_INCREF(var);
				return var;
			}

			// set info field of all individuals
			/**
			This function set info field of all individuals. Info can be used to
			sort individuals and set subpopulations. Therefore, the following code
			will do a migration:

			setIndSubPopID([ an_array_of_info_value ])

			setSubPopByIndID()

			\param info an array of info values, should have length of pop size
			\sa individual::setSubPopID, individual::subPopID
			*/
			void setIndSubPopID( const vectori& id)
			{
				DBG_ASSERT( id.size() == m_popSize, ValueError,
					"Info should have the same length as pop size");

				for(ULONG it=0; it < m_popSize; ++it)
					ind(it).setSubPopID( static_cast<SubPopID>(id[it]) );
			}

			/// set individual info with their subpopulation id.
			/**
			set individual info by subpop id.
			*/
			void setIndSubPopIDWithID()
			{
				for( UINT i=0, iEnd = numSubPop(); i < iEnd;  ++i)
					for(IndIterator it = indBegin(i), itEnd = indEnd(i); it != itEnd;  ++it)
						it -> setSubPopID(i);
			}

			/// adjust subpopulation according to individual info values
			/** assume individual has subpop index in their info value and put
			them into corresponding subpopulations.
			\param info optional info that will be used to set sub pop
			\note individual with negative info will be removed!
			\sa setIndSubPopID,
			*/
			void setSubPopByIndID(vectori id=vectori());

			/// split population
			/** split subpopulation 'which' into subpopulations with size specified in subPops,
			optionally given subPopID.
			The subPOP ID of Other subpop will be kept. For example, if
			subpop 1 of 0 1 2 3 is split into three parts, the new subpop id will be
			0 (1 4 5) 2 3.
			\note subpop with negative id will be removed. So, you can shrink
			one subpop by split and set one of the new subpop with negative id.

			*/
			void splitSubPop(UINT which, vectorlu sizes, vectoru subPopID=vectoru());

			/// split population
			/** split subpopulation 'which' into subpopulations with specified 'proportions',
			optionally given subPopID.
			\note subpop with negative id will be removed. So, you can shrink
			one subpop by split and set one of the new subpop with negative id.
			*/
			void splitSubPopByProportion(UINT which, vectorf proportions, vectoru subPopID=vectoru());

			/** remove empty subpops, this will adjust subPOP ID of other subpops */
			void removeEmptySubPops();

			/**  remove subpop, adjust subpop numbers so that there will be no 'empty'
			subpops left */
			void removeSubPops(const vectoru& subPops=vectoru(), bool shiftSubPopID=true, bool removeEmptySubPops=false);

			/**  remove subpop, adjust subpop numbers so that there will be no 'empty'
			subpops left */
			void removeIndividuals(const vectoru& inds=vectoru(), int subPop=-1, bool removeEmptySubPops=false);

			/// merge population
			/** merge subpopulations, subpop id will be the ID of the first in array subPops
			  all subpopulation will take the id of the first one.
			*/
			void mergeSubPops(vectoru subPops=vectoru(), bool removeEmptySubPops=false);

			/// \brief reorder subpopulations
			/**
			\param order new order of the subpopulations. For examples, 3 2 0 1
			means subpop3, subpop2, subpop0, subpop1 will be the new layout.
			\param rank you can also specify new rank for each subpop. For example, 3,2,0,1
			means the original subpopulations will have new ID 3,2,0,1. To achive order 3,2,0,1.
			the rank should be 1 0 2 3.
			*/
			void reorderSubPops(const vectoru& order=vectoru(), const vectoru& rank=vectoru(),
				bool removeEmptySubPops=false);

			/** form a new population according to info, info can be given directly
				keepAncestralPops=-1: keep all
				0: only current
				1: keep one ...
			*/
			population& newPopByIndID(int keepAncestralPops=-1,
				const vectori& id=vectori(),
				bool removeEmptySubPops=false);

			void removeLoci( const vectoru& remove=vectoru(), const vectoru& keep=vectoru());

			/** get a new population with selected loci */
			population& newPopWithPartialLoci(
				const vectoru& remove=vectoru(),
				const vectoru& keep=vectoru());

			// swap in rhs. (usually a scratch population
			// if fixRhs is false, leave rhs in broken status
			// if force is true, push current population
			// regardless of m_ancestray settings.
			void pushAndDiscard(population& rhs, bool force=false);

			UINT ancestralDepth()
			{
				return m_ancestralPops.size();
			}

			/// return the current ancestral gen index.
			int ancestralGen()
			{
				return m_curAncestralPop;
			}

			// int requestInfoField(const string name);

			/// set info for all individuals
			template<typename T, typename T1>
				void setIndInfo(const T& values, T1 idx)
			{
#ifdef SIMUMPI
				if(mpiRank()==0)
				{
#endif
					CHECKRANGEINFO(idx);
					DBG_ASSERT(values.size() == popSize(), IndexError,
						"Size of values should be the same as population size");
					UINT is = infoSize();
					typename T::const_iterator infoIter = values.begin();
					for(vectorinfo::iterator ptr=m_info.begin() + idx;
						ptr != m_info.end() + idx; ptr += is)
					*ptr = static_cast<InfoType>(*infoIter++);
#ifdef SIMUMPI
				}
#endif
			}

			template<class T>
				void setIndInfo(const T& values, const string& name)
			{
#ifdef SIMUMPI
				if(mpiRank()==0)
				{
#endif
					int idx = infoIdx(name);
					DBG_ASSERT(idx>=0, IndexError,
						"Info name " + name + " is not a valid values field name");
					setIndInfo<T, UINT>(values, idx);
#ifdef SIMUMPI
				}
#endif
			}

#ifdef SIMUMPI
			GappedInfoIterator infoBegin(UINT idx, bool order);
			GappedInfoIterator infoEnd(UINT idx, bool order);
			GappedInfoIterator infoBegin(UINT index, UINT subPop, bool order);
			GappedInfoIterator infoEnd(UINT index, UINT subPop, bool order);
			vectorinfo indInfo(UINT idx, bool order);
			vectorinfo indInfo(const string& name, bool order);
			vectorinfo indInfo(UINT idx, UINT subPop, bool order);
			vectorinfo indInfo(const string& name, UINT subPop, bool order);
			PyObject* arrIndInfo(bool order);
			PyObject* arrIndInfo(UINT subPop, bool order);
#else
			/// CPPONLY
			/// info iterator
			/// if order=true, keep order,
			/// if flase, do not respect pop structure
			GappedInfoIterator infoBegin(UINT idx, bool order)
			{
				CHECKRANGEINFO(idx);
				if(order && !infoOrdered())
					adjustInfoPosition(true);
				return GappedInfoIterator(m_info.begin()+idx, infoSize());
			}

			/// CPPONLY
			GappedInfoIterator infoEnd(UINT idx, bool order)
			{
				CHECKRANGEINFO(idx);
				if(order && !infoOrdered())
					adjustInfoPosition(true);
				return GappedInfoIterator(m_info.begin()+idx+m_info.size(), infoSize());
			}

			/// info iterator
			/// oder = true: keep order
			/// otherwise, respect subpop structure
			GappedInfoIterator infoBegin(UINT index, UINT subPop, bool order)
			{
				CHECKRANGEINFO(index);
				CHECKRANGESUBPOP(subPop);

				if(!infoOrdered())
					adjustInfoPosition(order);

				return GappedInfoIterator(m_info.begin()+index+m_subPopIndex[subPop]*infoSize(), infoSize());
			}

			///
			GappedInfoIterator infoEnd(UINT index, UINT subPop, bool order)
			{
				CHECKRANGEINFO(index);
				CHECKRANGESUBPOP(subPop);

				if(!infoOrdered())
					adjustInfoPosition(order);

				return GappedInfoIterator(m_info.begin()+index+m_subPopIndex[subPop+1]*infoSize(), infoSize());
			}

			vectorinfo indInfo(UINT idx, bool order)
			{
				return vectorinfo(infoBegin(idx, order), infoEnd(idx, order));
			}

			vectorinfo indInfo(const string& name, bool order)
			{
				UINT idx = infoIdx(name);
				return vectorinfo(infoBegin(idx, order), infoEnd(idx, order));
			}

			vectorinfo indInfo(UINT idx, UINT subPop, bool order)
			{
				return vectorinfo(infoBegin(idx, subPop, order), infoEnd(idx, subPop, order));
			}

			vectorinfo indInfo(const string& name, UINT subPop, bool order)
			{
				UINT idx = infoIdx(name);
				return vectorinfo(infoBegin(idx, subPop, order), infoEnd(idx, subPop, order));
			}

			/// if order: keep order
			/// otherwise: do not respect subpop info
			PyObject* arrIndInfo(bool order)
			{
				if(order && !infoOrdered())
					adjustInfoPosition(true);

				return Info_Vec_As_NumArray(m_info.begin(), m_info.end());
			}

			/// if order: keep order
			/// otherwise: respect subpop info
			PyObject* arrIndInfo(UINT subPop, bool order)
			{
				CHECKRANGESUBPOP(subPop);

				if(!infoOrdered())
					adjustInfoPosition(order);

				return Info_Vec_As_NumArray(m_info.begin() + m_subPopIndex[subPop]*infoSize(),
					m_info.begin() + m_subPopIndex[subPop+1]*infoSize());
			}
#endif

			/// add information field to a population
			/// if the field already exists in the geno structure
			/// match local info fields with the structure.
			///
			/// the difficult part is for ancestral generations.
			int addInfoField(const string field, double init=0);
			void addInfoFields(const vectorstr& fields, double init=0);

			/// set information fields, remove the old one
			void setInfoFields(const vectorstr& fields, double init=0);

			/// set ancestral depth, can be -1
			void setAncestralDepth(int depth);

			int ancestralPop()
			{
				return m_curAncestralPop;
			}

			// idx = 0 (current, do nothing) -1, -2, ..
			//
			void useAncestralPop(int idx);

			/// compare two populations
			bool equalTo(const population& rhs)
			{
				return(
					genoStru() == rhs.genoStru() &&
					m_subPopSize == rhs.m_subPopSize &&
					m_inds == rhs.m_inds );
			}

			//@}

			/// CPPONLY
			/// some iterators requires that genotype information is within
			/// each subpopulation. We need to adjust genotypic info to
			/// obey this.
			/// order=true: make individuals in order
			/// order=false: make individuals in each subpopulation
			void adjustGenoPosition(bool order);
			void adjustInfoPosition(bool order);

			/// save population to a file
			/**
			\param filename save to filename
			\param format format to save. Can be one of 'text', 'bin', 'xml'

			The default format is 'text' but the output is not suppored to be read.
			'bin' has smaller size and should be used for large populations.
			'xml' format is most readable and should be used when you would
			like to convert simuPOP populations to other formats.
			\sa global function loadPopulation
			*/
			void savePopulation(const string& filename, const string& format="auto", bool compress=true) const;

			/// CPPONLY load population from a file
			/**
			\param filename load from filename
			\param format format to load. Can be one of "text", "bin", "xml". It should match
			  the format used to save the population.

			\sa savePopulation
			*/
			void loadPopulation(const string& filename, const string& format="auto");

		public:

			int rep()
			{
				return m_rep;
			}

			/// CPPONLY  set rep number
			void setRep(int rep, bool setVar=true)
			{
				m_rep = rep;
				if(setVar)
					m_vars.setIntVar("rep", rep);
			}

			int grp()
			{
				return m_grp;
			}

			/// CPPONLY
			void setGrp(int grp, bool setVar=true)
			{
				m_grp = grp;
				if(setVar)
					m_vars.setIntVar("grp", grp);
			}

			ULONG gen()
			{
				return m_gen;
			}

			/// CPPONLY
			void setGen(ULONG gen, bool setVar=true)
			{
				m_gen = gen;
				if(setVar)
					m_vars.setIntVar("gen", gen);
			}

			/// return variables of this population
			/// if subPop is given, return dictionary for
			/// specified subpopulation.
			PyObject* vars(int subPop=-1);

			/// CPPONLY
			/// The same as vars(), but without increasing
			/// reference count.
			PyObject* dict(int subPop=-1);

			/// CPPONLY
			void setDict(PyObject* dict)
			{
				DBG_ASSERT(dict != NULL, SystemError, "Dictionary is empty");
				m_vars.setDict(dict);
			}

			///
			bool hasVar(const string& name)
			{
				return m_vars.hasVar(name);
			}

			/// CPPNLY
			void removeVar(const string& name)
			{
				m_vars.removeVar(name);
			}

			/// CPPONLY
			PyObject* setBoolVar(const string& name, const bool val)
			{
				return  m_vars.setBoolVar(name, val);
			}

			/// CPPONLY
			PyObject* setIntVar(const string& name, const int val)
			{
				return m_vars.setIntVar(name, val);
			}

			/// CPPONLY
			PyObject* setDoubleVar(const string& name, const double val)
			{
				return m_vars.setDoubleVar(name, val);
			}

			/// CPPONLY
			PyObject* setStringVar(const string& name, const string& val)
			{
				return m_vars.setStringVar(name, val);
			}

			///CPPONLY
			PyObject* setIntVectorVar(const string& name, const vectori& val)
			{
				return m_vars.setIntVectorVar(name, val);
			}

			///CPPONLY
			PyObject* setDoubleVectorVar(const string& name, const vectorf& val)
			{
				return m_vars.setDoubleVectorVar(name, val);
			}

			/// CPPONLY
			PyObject* setStrDictVar(const string& name, const strDict& val)
			{
				return m_vars.setStrDictVar(name, val);
			}

			/// CPPONLY
			PyObject* setIntDictVar(const string& name, const intDict& val)
			{
				return m_vars.setIntDictVar(name, val);
			}

			/// CPPONLY
			PyObject* setVar(const string& name, PyObject* val)
			{
				return m_vars.setVar(name, val);
			}

			/// CPPONLY
			PyObject* getVar(const string&name, bool nameError=true)
			{
				return m_vars.getVar(name, nameError);
			}

			/// CPPONLY
			bool getVarAsBool(const string& name, bool nameError=true)
			{
				return m_vars.getVarAsBool(name, nameError);
			}

			/// CPPONLY
			int getVarAsInt(const string& name, bool nameError=true)
			{
				return m_vars.getVarAsInt(name, nameError);
			}

			/// CPPONLY
			double getVarAsDouble(const string& name, bool nameError=true)
			{
				return m_vars.getVarAsDouble(name, nameError);
			}

			/// CPPONLY
			string getVarAsString(const string& name, bool nameError=true)
			{
				return m_vars.getVarAsString(name, nameError);
			}

			/// CPPONLY
			strDict getVarAsStrDict(const string& name, bool nameError=true)
			{
				return m_vars.getVarAsStrDict(name, nameError);
			}

			/// CPPONLY
			intDict getVarAsIntDict(const string& name, bool nameError=true)
			{
				return m_vars.getVarAsIntDict(name, nameError);
			}

			/// CPPONLY
			string varsAsString() const
			{
				return m_vars.asString();
			}

			/// CPPONLY
			void varsFromString(const string& vars)
			{
				return m_vars.fromString(vars);
			}

			/// evaluate python statment/expressions
			/**
			  this function evaluate python expressions
			  and return as string representing the result
				*/
			PyObject* evaluate(const string& expr="", const string& stmts="")
			{
				return Expression(expr, stmts, m_vars.dict() ).evaluate();
			}

			///
			void execute(const string& stmts="")
			{
				Expression("", stmts, m_vars.dict() ).evaluate();
			}

		private:

			friend class boost::serialization::access;

			template<class Archive>
				void save(Archive &ar, const UINT version) const
			{
#ifdef SIMUMPI
				PENDING_MPI;
#else

				// deep adjustment: everyone in order
				//const_cast<population*>(this)->adjustGenoPosition(true);

				ar & make_nvp("libraryMaxAllele", MaxAllele);

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
				for(size_t i=0; i < blks; ++i)
				{
					tmp = *ptr++;
					for(size_t j = 0; j < WORDBIT/32; ++j)
					{
						tmp1 = tmp & 0xFFFFFFFF;
						tmp = tmp >> 32;
						ar & make_nvp("blocks", tmp1);
					}
				}
				// last block
				if (rest > 0)
				{
					tmp = *ptr;
					for(size_t j = 0; j <= rest/32; ++j)
					{
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
				for(size_t i=0; i< m_ancestralPops.size(); ++i)
				{
					const_cast<population*>(this)->useAncestralPop(i+1);
					// need to make sure ancestral pop also in order
					const_cast<population*>(this)->adjustGenoPosition(true);
					ar & make_nvp("subPop_sizes", m_subPopSize);
#ifdef BINARYALLELE
					size_t size = m_genotype.size();
					ar & make_nvp("size", size);
					WORDTYPE * ptr = BITPTR(m_genotype.begin());
					size_t blks = size / WORDBIT;
					size_t rest = size - blks * WORDBIT;
					DBG_ASSERT(WORDBIT >= 32, SystemError, "WordBit should be at least 32 bits");

					WORDTYPE tmp, tmp1;
					for(size_t i=0; i < blks; ++i)
					{
						tmp = *ptr++;
						for(size_t j = 0; j < WORDBIT/32; ++j)
						{
							tmp1 = tmp & 0xFFFFFFFF;
							tmp = tmp >> 32;
							ar & make_nvp("blocks", tmp1);
						}
					}
					// last block
					if (rest > 0)
					{
						tmp = *ptr;
						for(size_t j = 0; j <= rest/32; ++j)
						{
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
				const_cast<population*>(this)->useAncestralPop(0);

				// save shared variables as string.
				// note that many format are not supported.
				try
				{
					DBG_DO(DBG_POPULATION, cout << "Handling shared variables" << endl);
					string vars = varsAsString();
					ar & make_nvp("vars", vars);
				}
				catch(...)
				{
					cout << "Warning: shared variable is not saved correctly.\npopulation should still be usable." << endl;
				}
#endif
			}

			template<class Archive>
				void load(Archive &ar, const UINT version)
			{
#ifdef SIMUMPI
				PENDING_MPI;
#else

				ULONG ma;
				ar & make_nvp("libraryMaxAllele", ma);

				if( ma > MaxAllele)
					cout << "Warning: The population is saved in library with more allele states. \n"
						<< "Unless all alleles are less than "  << MaxAllele
						<< ", you should use the modules used to save this file. (c.f. simuOpt.setOptions()\n";

				GenoStructure stru;
				DBG_DO(DBG_POPULATION, cout << "Handling geno structure" << endl);
				ar & make_nvp("geno_structure", stru);
				ar & make_nvp("subPop_sizes", m_subPopSize);
				DBG_DO(DBG_POPULATION, cout << "Handling genotype" << endl);

				if (version <= 1)
				{
					ar & make_nvp("genotype", m_genotype);
				}
				// the new version
				// the binary genotypes are saved in an efficient way
				else
				{
#ifdef BINARYALLELE
					// binary from binary
					if (ma == 1)
					{
						DBG_DO(DBG_POPULATION, cout << "Load bin from bin. " << endl);
						size_t size;
						ar & make_nvp("size", size);
						size_t blks = size / WORDBIT;
						size_t rest = size - blks * WORDBIT;

						m_genotype.resize(size);
						WORDTYPE tmp, tmp1;
						WORDTYPE * ptr = BITPTR(m_genotype.begin());
						for(size_t i=0; i < blks; ++i)
						{
							tmp = 0;
							for(size_t j = 0; j < WORDBIT/32; ++j)
							{
								ar & make_nvp("blocks", tmp1);
								tmp |= tmp1 << (j * 32);
							}
							*ptr++ = tmp;
						}
						// last block
						if (rest > 0)
						{
							tmp = 0;
							for(size_t j = 0; j <= rest/32; ++j)
							{
								ar & make_nvp("blocks", tmp1);
								tmp |= tmp1 << (j * 32);
							}
							*ptr = tmp;
						}
					}
					// binary from others (long types)
					else
					{
						DBG_DO(DBG_POPULATION, cout << "Load bin from long. " << endl);
						vector<unsigned char> tmpgeno;
						ar & make_nvp("genotype", tmpgeno);
						m_genotype = vectora(tmpgeno.begin(), tmpgeno.end());
					}
#else
					// long from binary
					if (ma == 1)
					{
						DBG_DO(DBG_POPULATION, cout << "Load long from bin. " << endl);
						// for version 2 and higher, archive in 32bit blocks.
						size_t size;
						ar & make_nvp("size", size);
						m_genotype.resize(size);
						size_t blks = size / 32;
						size_t rest = size - blks * 32;
						DBG_ASSERT(WORDBIT >= 32, SystemError, "WordBit should be at least 32 bits");

						GenoIterator ptr = m_genotype.begin();
						WORDTYPE tmp;
						for(size_t i=0; i < blks; ++i)
						{
							ar & make_nvp("blocks", tmp);
							for(size_t j = 0; j < 32; ++j)
							{
								*ptr++ = (tmp & 1UL) != 0;
								tmp = tmp >> 1;
							}
						}
						// last block
						if (rest > 0)
						{
							ar & make_nvp("blocks", tmp);
							for(size_t j = 0; j < rest; ++j)
							{
								*ptr++ = (tmp & 1UL) != 0;
								tmp = tmp >> 1;
							}
						}
					}							  // if ma == 1
					else						  // for non-binary types, ...
					{
						DBG_DO(DBG_POPULATION, cout << "Load long from long. " << endl);
						// long from long
						ar & make_nvp("genotype", m_genotype);
					}
#endif
				}								  // verion >= 2

				if ( version > 0)
				{
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

				DBG_FAILIF(m_info.size() != m_popSize*infoSize(), ValueError, "Wgong size of info vector");

				if( m_popSize != m_inds.size() )
				{
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
				for(ULONG i=0; i< m_popSize; ++i, ptr += step, infoPtr += infoStep)
				{
					m_inds[i].setGenoStruIdx(genoStruIdx());
					m_inds[i].setGenoPtr( ptr );
					m_inds[i].setInfoPtr( infoPtr );
				}
				m_ancestralDepth = 0;
				m_ancestralPops.clear();

				// ancestry populations
				DBG_DO(DBG_POPULATION, cout << "Handling ancestral populations" << endl);
				ar & make_nvp("ancestry", m_ancestralDepth);
				size_t na;
				ar & make_nvp("numOfAncestralPops", na);
				for(size_t ap=0; ap< na; ++ap)
				{
					popData pd;
					ar & make_nvp("subPop_sizes", pd.m_subPopSize);
					// version <= 1, direct handling
					if (version <= 1)
					{
						ar & make_nvp("genotype", pd.m_genotype);
					}
					else
					{
#ifdef BINARYALLELE
						// binary from binary
						if(ma == 1)
						{
							DBG_DO(DBG_POPULATION, cout << "Load bin from bin. " << endl);
							size_t size;
							ar & make_nvp("size", size);
							size_t blks = size / WORDBIT;
							size_t rest = size - blks * WORDBIT;

							pd.m_genotype.resize(size);
							WORDTYPE * ptr = BITPTR(pd.m_genotype.begin());
							WORDTYPE tmp, tmp1;
							for(size_t i=0; i < blks; ++i)
							{
								tmp = 0;
								for(size_t j = 0; j < WORDBIT/32; ++j)
								{
									ar & make_nvp("blocks", tmp1);
									tmp |= tmp1 << (j * 32);
								}
								*ptr++ = tmp;
							}
							// last block
							if (rest > 0)
							{
								tmp = 0;
								for(size_t j = 0; j <= rest/32; ++j)
								{
									ar & make_nvp("blocks", tmp1);
									tmp |= tmp1 << (j * 32);
								}
								*ptr = tmp;
							}
						}
						else
						{
							DBG_DO(DBG_POPULATION, cout << "Load bin from long. " << endl);
							// binary from long types
							vector<unsigned char> tmpgeno;
							ar & make_nvp("genotype", tmpgeno);
							pd.m_genotype = vectora(tmpgeno.begin(), tmpgeno.end());
						}
#else
						if(ma == 1)
						{
							DBG_DO(DBG_POPULATION, cout << "Load long from bin. " << endl);
							// long type from binary
							size_t size;
							ar & make_nvp("size", size);
							size_t blks = size / 32;
							size_t rest = size - blks * 32;
							pd.m_genotype.resize(size);

							ptr = pd.m_genotype.begin();
							WORDTYPE tmp;
							for(size_t i=0; i < blks; ++i)
							{
								ar & make_nvp("blocks", tmp);
								for(size_t j = 0; j < 32; ++j)
								{
									*ptr++ = (tmp & 1UL) != 0;
									tmp = tmp >> 1;
								}
							}
							// last block
							if (rest > 0)
							{
								ar & make_nvp("blocks", tmp);
								for(size_t i=0; i < rest; ++i)
								{
									*ptr++ = (tmp & 1UL) != 0;
									tmp = tmp >> 1;
								}
							}
						}
						else
						{
							DBG_DO(DBG_POPULATION, cout << "Load long from long. " << endl);
							// long type from long type.
							ar & make_nvp("genotype", pd.m_genotype);
						}
#endif
					}
					if ( version > 0)
						ar & make_nvp("info", pd.m_info);
					ar & make_nvp("individuals", pd.m_inds);
					// set pointer after copy this thing again (push_back)
					m_ancestralPops.push_back(pd);
					// now set pointers
					popData& p = m_ancestralPops.back();
					// set pointers
					vector<individual>& inds = p.m_inds;
					ULONG ps = inds.size();
					ptr = p.m_genotype.begin();
					infoPtr = p.m_info.begin();

					for(ULONG i=0; i< ps; ++i, ptr+=step, infoPtr += infoStep)
					{
						inds[i].setGenoPtr( ptr );
						inds[i].setInfoPtr( infoPtr );
						// set new genoStructure
						inds[i].setGenoStruIdx(genoStruIdx());
						// fresh copy so clear shallowcopied flag.
						inds[i].setShallowCopied(false);
					}
				}

				// load vars from string
				try
				{
					DBG_DO(DBG_POPULATION, cout << "Handling shared variables" << endl);
					string vars;
					ar & make_nvp("vars", vars);
					varsFromString(vars);
				}
				catch(...)
				{
					cout << "Warning: shared variable is not loaded correctly.\npopulation should still be usable." << endl;
				}

				m_shallowCopied = false;
				m_infoOrdered = true;
#endif
			}

			BOOST_SERIALIZATION_SPLIT_MEMBER();

		private:
			/// population size: number of individual
			/// MPI: in all nodes
			ULONG m_popSize;

			/// number of subpopulations
			/// MPI: in all nodes
			UINT m_numSubPop;

			/// size of each subpopulation
			/// MPI: in all nodes
			vectorlu m_subPopSize;

			/// index to subPop \todo change to vectorl
			/// MPI: in all nodes
			vectorlu m_subPopIndex;

			/// pool of genotypic information
			/// MPI: empty in head node, one chromosome for others.
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
			int m_curAncestralPop;

			/// whether or not individual genotype is in order
			bool m_shallowCopied;

			/// whether or not information is ordered
			bool m_infoOrdered;
	};

	population& LoadPopulation(const string& file, const string& format="auto");

	/// get info through ind.info()
	vectorf testGetinfoFromInd(population& pop);

	/// get info through GappedInfoIterator
	vectorf testGetinfoFromPop(population& pop, bool order);

}



#ifndef SWIG
#ifndef _NO_SERIALIZATION_
// version 0: base
// version 1: save info
// version 2: reduce binary file size
BOOST_CLASS_VERSION(simuPOP::population, 2)
#endif
#endif
#endif
