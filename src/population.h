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
#ifndef __NO_XML_SUPPORT__
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#endif
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

	/** \brief a collection of individuals with subpopulation structure

	Please refer to user's Guide for details about this object.

	*/
	class population : public GenoStruTrait
	{
		public:

			/// individual itertor, used to iterate all individuals.
			typedef vector<individual>::iterator IndIterator;

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
			\param infoName: name of information fields that will be attached to each
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
				const vectorstr& infoName = vectorstr());

			/// CPPONLY copy constructor
			population(const population& rhs);

			///
			population& clone(bool keepAncestralPops=true) const
			{
				population& p = *new population(*this);
				if( keepAncestralPops == false)
					p.m_ancestralPops.clear();
				return p;
			}

			/// SWAP population
			/// swap the content of two populations
			void swap(population& rhs)
			{
				GenoStruTrait::swap(rhs);
				std::swap(m_popSize, rhs.m_popSize);
				std::swap(m_numSubPop, rhs.m_numSubPop);
				m_subPopSize.swap(rhs.m_subPopSize);
				std::swap(m_popGenoSize, rhs.m_popGenoSize);
				m_subPopIndex.swap(rhs.m_subPopIndex);
				m_genotype.swap(rhs.m_genotype);
				m_inds.swap(rhs.m_inds);
				std::swap(m_ancestralDepth, rhs.m_ancestralDepth);
				m_vars.swap(rhs.m_vars);
				m_ancestralPops.swap(rhs.m_ancestralPops);
				std::swap(m_rep, rhs.m_rep);
				std::swap(m_grp, rhs.m_grp);
				std::swap(m_gen, rhs.m_gen);
				std::swap(m_curAncestralPop, rhs.m_curAncestralPop);
				// FIXME: remove m_fitness, use generalized info field.
				m_fitness.swap(rhs.m_fitness);
				std::swap(m_shallowCopied, rhs.m_shallowCopied);
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
			*/
			GappedAlleleIterator alleleBegin(UINT locus)
			{
				CHECKRANGEABSLOCUS(locus);

				return GappedAlleleIterator( m_genotype.begin()+locus, totNumLoci());
			}

			/// CPPONLY allele iterator
			GappedAlleleIterator alleleEnd(UINT locus)
			{
				CHECKRANGEABSLOCUS( locus);

				return GappedAlleleIterator( m_genotype.begin() + locus + m_popGenoSize , totNumLoci());
			}

			///  CPPONLY allele begin, for given subPopu
			GappedAlleleIterator alleleBegin(UINT locus, UINT subPop)
			{
				CHECKRANGEABSLOCUS(locus);
				CHECKRANGESUBPOP(subPop);

				if(shallowCopied())
					adjustGenoPosition();

				return GappedAlleleIterator( m_genotype.begin() + m_subPopIndex[subPop]*genoSize() +
					locus, totNumLoci());
			}

			///  CPPONLY allele iterator
			GappedAlleleIterator alleleEnd( UINT locus, UINT subPop)
			{
				CHECKRANGEABSLOCUS(locus);
				CHECKRANGESUBPOP(subPop);

				if(shallowCopied())
					adjustGenoPosition();

				return GappedAlleleIterator( m_genotype.begin() + m_subPopIndex[subPop+1]*genoSize() +
					locus, totNumLoci());
			}

			///  CPPONLY allele iterator, go through all allels one by one, without subPop info
			GenoIterator genoBegin()
			{
				return m_genotype.begin();
			}

			///  CPPONLY allele iterator
			GenoIterator genoEnd()
			{
				return m_genotype.end();
			}

			///  CPPONLY allele iterator, go through all allels one by one in a subpopulation
			GenoIterator genoBegin(UINT subPop)
			{
				CHECKRANGESUBPOP(subPop);

				if(shallowCopied())
					adjustGenoPosition();

				return m_genotype.begin() + m_subPopIndex[subPop]*genoSize();
			}

			///  CPPONLY allele iterator in a subpopulation.
			GenoIterator genoEnd(UINT subPop)
			{
				CHECKRANGESUBPOP(subPop);

				if(shallowCopied())
					adjustGenoPosition();
				return m_genotype.begin() + m_subPopIndex[subPop+1]*genoSize();
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
			PyObject* arrGenotype()
			{
				if(shallowCopied())
					// adjust position. deep=true
					adjustGenoPosition(true);

				// false: directly expose values. Do not copy data over.
				return Allele_Vec_As_NumArray(m_genotype.begin(), m_genotype.end());
			}

			/// get the whole genotype.
			/// individuals will be in order before exposing
			/// their genotypes.
			PyObject* arrGenotype(UINT subPop)
			{
				CHECKRANGESUBPOP(subPop);
				if(shallowCopied())
					// adjust position. deep=true
					adjustGenoPosition(true);

				return Allele_Vec_As_NumArray( genoBegin(subPop), genoEnd(subPop));
			}

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
			void setIndSubPopID( const vectori& info)
			{
				DBG_ASSERT( info.size() == m_popSize, ValueError,
					"Info should have the same length as pop size");

				for(ULONG it=0; it < m_popSize; ++it)
					ind(it).setSubPopID( static_cast<SubPop_ID>(info[it]) );
			}

			/// set individual info with their subpopulation id.
			/**
			set individual info by subpop id.
			*/
			void setIndSubPopIDWithID()
			{
				for( UINT i=0, iEnd = numSubPop(); i < iEnd;  ++i)
					for( IndIterator it = indBegin(i), itEnd = indEnd(i); it != itEnd;  ++it)
						it -> setSubPopID(i);
			}

			/// adjust subpopulation according to individual info values
			/** assume individual has subpop index in their info value and put
			them into corresponding subpopulations.
			\param info optional info that will be used to set sub pop
			\note individual with negative info will be removed!
			\sa setIndSubPopID,
			*/
			void setSubPopByIndID(vectori info=vectori());

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

			/** form a new population according to info, info can be given directly */
			population& newPopByIndID(bool keepAncestralPops=true,
				vectori info=vectori(), bool removeEmptySubPops=false);

			void removeLoci( const vectoru& remove=vectoru(), const vectoru& keep=vectoru());

			/** get a new population with selected loci */
			population& newPopWithPartialLoci(
				const vectoru& remove=vectoru(),
				const vectoru& keep=vectoru());

			// return reference to fitness vector
			// CPPONLY
			vectorf& fitness()
			{
				return m_fitness;
			}

			// this funciton allows fitness to be
			// accessed from python.
			PyObject* arrFitness()
			{
				if(m_fitness.size() != m_popSize)
					m_fitness.resize(m_popSize);
				return Double_Vec_As_NumArray(m_fitness.begin(), m_fitness.end());
			}

			// swap in rhs. (usually a scratch population
			// if fixRhs is false, leave rhs in broken status
			// if force is true, push current population
			// regardless of m_ancestray settings.
			void pushAndDiscard(population& rhs, bool force=false);

			UINT ancestralDepth()
			{
				return m_ancestralPops.size();
			}

			int requestInfoField(const string name);

			/// set info for all individuals
			template<class T>
			void setIndInfo(const T& info, UINT index)
			{
				CHECKRANGEINFO(index);
				DBG_ASSERT(info.size() == popSize(), IndexError,
					"Size of info should be the same as population size");
				UINT is = infoSize();
				typename T::const_iterator infoIter = info.begin();
				for(vector<InfoType>::iterator ptr=m_info.begin() + index;
					ptr != m_info.end() + index; ptr += is)
					*ptr = static_cast<InfoType>(*infoIter++);
			}

			/// info iterator
			GappedInfoIterator infoBegin(UINT index)
			{
				CHECKRANGEINFO(index);
				return GappedInfoIterator(m_info.begin()+index, infoSize());
			}
			
			GappedInfoIterator infoEnd(UINT index)
			{
				CHECKRANGEINFO(index);
				return GappedInfoIterator(m_info.begin()+index+m_info.size(), infoSize());
			}
			
			/// info iterator
			GappedInfoIterator infoBegin(UINT index, UINT subPop)
			{
				CHECKRANGEINFO(index);
				CHECKRANGESUBPOP(subPop);

				if(shallowCopied())
					adjustGenoPosition();

				return GappedInfoIterator(m_info.begin()+index+m_subPopIndex[subPop]*infoSize(), infoSize());
			}
		
			PyObject* arrIndInfo()
			{
				return Info_Vec_As_NumArray(m_info.begin(), m_info.end());
			}	

			PyObject* arrIndInfo(UINT subPop)
			{
				CHECKRANGESUBPOP(subPop);

				if(shallowCopied())
					adjustGenoPosition();

				return Info_Vec_As_NumArray(m_info.begin() + m_subPopIndex[subPop]*infoSize(), 
					m_info.begin() + m_subPopIndex[subPop+1]*infoSize());
			}
			
			GappedInfoIterator infoEnd(UINT index, UINT subPop)
			{
				CHECKRANGEINFO(index);
				CHECKRANGESUBPOP(subPop);

				if(shallowCopied())
					adjustGenoPosition();

				return GappedInfoIterator(m_info.begin()+index+m_subPopIndex[subPop+1]*infoSize(), infoSize());
			}
			
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
			void adjustGenoPosition(bool deep=false);

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

#ifndef OPTIMIZED
			void checkRefCount()
			{
				m_vars.checkRefCount();
			}
#endif

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
				// deep adjustment: everyone in order
				//const_cast<population*>(this)->adjustGenoPosition(true);

				ar & make_nvp("libraryMaxAllele", MaxAllele);

				DBG_DO(DBG_POPULATION, cout << "Handling geno structure" << endl);
				// GenoStructure genoStru = this->genoStru();
				ar & make_nvp("geno_structure", this->genoStru());
				ar & make_nvp("subPop_sizes", m_subPopSize);
				DBG_DO(DBG_POPULATION, cout << "Handling genotype" << endl);
				ar & make_nvp("genotype", m_genotype);
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
					ar & make_nvp("genotype", m_genotype);
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
			}

			template<class Archive>
				void load(Archive &ar, const UINT version)
			{
				ULONG ma;
				ar & make_nvp("libraryMaxAllele", ma);

				if( ma > MaxAllele)
					cout << "Warning: The population is saved in library with more allele states. \n"
						<< "Unless all alleles are less than "  << MaxAllele
						<< ", you should use the long allele library. (c.f. simuOpt.setOptions()\n";

				GenoStructure  stru;
				DBG_DO(DBG_POPULATION, cout << "Handling geno structure" << endl);
				ar & make_nvp("geno_structure", stru);
				ar & make_nvp("subPop_sizes", m_subPopSize);
				DBG_DO(DBG_POPULATION, cout << "Handling genotype" << endl);
				ar & make_nvp("genotype", m_genotype);
				DBG_DO(DBG_POPULATION, cout << "Handling individuals" << endl);
				ar & make_nvp("individuals", m_inds);

				// set genostructure, check duplication
				// we can not use setGenoStruIdx since stru may be new.
				this->setGenoStructure(stru);

				m_numSubPop = m_subPopSize.size();
				m_popSize = accumulate(m_subPopSize.begin(), m_subPopSize.end(), 0L);

				if( m_popSize != m_inds.size() )
				{
					cout << "Number of individuals loaded" << m_inds.size() << endl;
					cout << "population size" << m_popSize << endl;
					throw ValueError("Number of individuals does not match population size.\n"
						"Please use the same (binary, short or long) module to save and load files.");
				}

				m_popGenoSize = totNumLoci() * GenoStruTrait::ploidy()
					* m_popSize;

				DBG_DO(DBG_POPULATION, cout << "Reconstruct individual genotype" << endl);
				m_subPopIndex.resize(m_numSubPop + 1);
				UINT i = 1;
				for (m_subPopIndex[0] = 0; i <= m_numSubPop; ++i)
					m_subPopIndex[i] = m_subPopIndex[i - 1] + m_subPopSize[i - 1];

				// assign genotype location and set structure information for individuals
				GenoIterator ptr = m_genotype.begin();
				UINT step = genoSize();
				for(ULONG i=0; i< m_popSize; ++i, ptr+=step)
				{
					m_inds[i].setGenoStruIdx(genoStruIdx());
					m_inds[i].setGenoPtr( ptr );
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
					ar & make_nvp("genotype", pd.m_genotype);
					ar & make_nvp("individuals", pd.m_inds);
					// set pointer after copy this thing again (push_back)
					m_ancestralPops.push_back(pd);
					// now set pointers
					popData& p = m_ancestralPops.back();
					// set pointers
					vector<individual>& inds = p.m_inds;
					ULONG ps = inds.size();
					ptr = p.m_genotype.begin();

					for(ULONG i=0; i< ps; ++i, ptr+=step)
					{
						inds[i].setGenoPtr( ptr );
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

				// m_fitness was not saved
				m_fitness = vectorf();
				m_shallowCopied = false;
			}

			BOOST_SERIALIZATION_SPLIT_MEMBER();

		private:
			/// population size: number of individual
			ULONG m_popSize;

			/// number of subpopulations
			UINT m_numSubPop;

			/// size of each subpopulation
			vectorlu m_subPopSize;

			///  size of genotypic information of the whole poopulation
			LONG m_popGenoSize;

			/// index to subPop \todo change to vectorl
			vectorlu m_subPopIndex;

			/// pool of genotypic information
			vectora m_genotype;

			/// information
			vector<InfoType> m_info;

			/// individuals.
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
				vector<InfoType> m_info;
				vector<individual> m_inds;
			};

			std::deque<popData> m_ancestralPops;

			/// curent replicate, group number
			int m_rep, m_grp;

			/// generation
			ULONG m_gen;

			/// current ancestral depth
			int m_curAncestralPop;

			/// fitness
			vectorf m_fitness;

			/// whether or not individual genotype is in order
			bool m_shallowCopied;
	};

	population& LoadPopulation(const string& file, const string& format="auto");

}


#ifndef SWIG
#ifndef _NO_SERIALIZATION_
// version 0: base
BOOST_CLASS_VERSION(simuPOP::population, 0)
#endif
#endif
#endif
