/************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu                                                        *
 *                                                                         *
 *   $LastChangedDate: 2007-01-08 22:49:35 -0600 (Mon, 08 Jan 2007) $
 *   $Rev: 726 $
 *
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

#ifndef _GENOSTRU_H
#define _GENOSTRU_H

/**
\file
\brief class genoStru and genoTrait
*/

#include "utility.h"
#include "simuPOP_cfg.h"

//
// the following is required by a vc7.1 bug.
#if  defined(_WIN32) || defined(__WIN32__)
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <fstream>
using std::ofstream;
using std::ifstream;
#endif											  // win32

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>
using boost::serialization::make_nvp;

#include <iterator>
using std::ostream;
using std::ostream_iterator;

#include <algorithm>
using std::copy;

#include <iostream>
using std::cout;
using std::endl;
using std::hex;
using std::dec;

#include <numeric>
using std::pair;

namespace simuPOP
{

	/** \brief CPPONLY genetic structure. Shared by individuals of one population

	populations create a copy of GenoStrcture and assign its pointer to each individual.
	This strcuture will be destroyed when population is destroyed.

	population with the same geneotype structure as an old one will use that,
	instead of creating a new one. This is ensured by GenoStructureTrait.

	Different populations will have different individuals but comparison, copy etc
	are forbidden even if they do have the same genotypic structure.
	 */
	class GenoStructure
	{

		public:

			/// CPPONLY serialization library requires a default constructor
			GenoStructure():m_ploidy(2), m_totNumLoci(0), m_genoSize(0), m_numChrom(0),
				m_numLoci(0), m_sexChrom(false), m_lociPos(0), m_chromIndex(0),
				m_alleleNames(), m_lociNames(), m_maxAllele(), m_infoFields(0),
				m_chromMap()
	#ifdef SIMUMPI
				, m_beginChrom(0), m_endChrom(0), m_beginLocus(0), m_endLocus(0),
				m_localNumLoci(0), m_localGenoSize(0), m_localChromIndex(0)
	#endif
				{}

			/** \brief constructor. The ONLY way to construct this strucuture. There is not set... functions
			CPPONLY
			\param ploidy number of sets of chromosomes
			\param loci number of loci on each chromosome.
			\param lociPos loci distance on each chromosome. the default values
			   are 1,2,etc.
			\param alleleNames allele names
			\param lociNames name of loci
			\param maxAllele maximum possible allele number for all alleles.
			\param length of info field
			*/
			GenoStructure(UINT ploidy, const vectoru& loci, bool sexChrom,
				const vectorf& lociPos, const vectorstr& alleleNames,
				const vectorstr& lociNames, UINT maxAllele, const vectorstr& infoFields,
				const vectori& chromMap);

			bool operator== (const GenoStructure& rhs);

			bool operator!= (const GenoStructure& rhs)
			{
				return !( *this == rhs);
			}

			/// destructor, do nothing.
			~GenoStructure()
			{
			}

#if  defined(_WIN32) || defined(__WIN32__)

			// due to an weird compiling error fo vc7.1,
			// if I do not specify these two functions, the ar & process
			// will fail to compile.
			// This will only be defined for win32 system
			/// CPPONLY
			void saveStru(string filename)
			{
				ofstream ofs(filename.c_str());
				boost::archive::binary_oarchive oa(ofs);
				oa << boost::serialization::make_nvp("geno_structure",*this);
			}

			/// CPPONLY
			void loadStru(string filename)
			{
				ifstream ifs(filename.c_str());

				boost::archive::binary_iarchive ia(ifs);
				ia >> boost::serialization::make_nvp("geno_structure",*this);
			}
#endif									  // win32

		private:

			friend class boost::serialization::access;

			template<class Archive>
				void save(Archive &ar, const UINT version) const
			{
#ifdef SIMUMPI
				if(mpiRank() == 0)
				{
#endif
					ar & make_nvp("ploidy", m_ploidy);
					ar & make_nvp("num_of_chrom", m_numChrom);
					ar & make_nvp("num_of_loci_on_each_chrom", m_numLoci);
					ar & make_nvp("sex_chromosome", m_sexChrom);
					ar & make_nvp("loci_distance_on_chrom", m_lociPos);
					ar & make_nvp("allele_name", m_alleleNames);
					ar & make_nvp("loci_name", m_lociNames);
					ar & make_nvp("max_allele", m_maxAllele);
					ar & make_nvp("info_name", m_infoFields);
					/// do not save load chromosome map
#ifdef SIMUMPI
				}
#endif
			}

			template<class Archive>
				void load(Archive &ar, const UINT version)
			{
#ifdef SIMUMPI
				if(mpiRank() == 0)
				{
#endif

					ar & make_nvp("ploidy", m_ploidy);
					ar & make_nvp("num_of_chrom", m_numChrom);
					ar & make_nvp("num_of_loci_on_each_chrom", m_numLoci);
					// after simuPOP 0.6.8, we have m_sexChrom
					// before that, there is no sex chromosome
					if(version > 0)
						ar & make_nvp("sex_chromosome", m_sexChrom);
					else
						m_sexChrom = false;
					ar & make_nvp("loci_distance_on_chrom", m_lociPos);
					ar & make_nvp("allele_name", m_alleleNames);
					ar & make_nvp("loci_name", m_lociNames);
					ar & make_nvp("max_allele", m_maxAllele);
					if(version > 1)
						ar & make_nvp("info_name", m_infoFields);

					// build chromosome index
					m_chromIndex.resize(m_numLoci.size()+1);
					ULONG i;
					for(m_chromIndex[0] = 0, i = 1; i <= m_numChrom; ++i)
						m_chromIndex[i] = m_chromIndex[i - 1] + m_numLoci[i - 1];

					m_totNumLoci = m_chromIndex[m_numChrom];
					m_genoSize = m_totNumLoci*m_ploidy;
					/// do not save load chromosome map
#ifdef SIMUMPI
				}
#endif
			}

			BOOST_SERIALIZATION_SPLIT_MEMBER();

			/// ploidy
			UINT m_ploidy;

			/// total number of loci
			UINT m_totNumLoci;

			/// total number of loci times ploidy
			UINT m_genoSize;

			/// number of chrom
			UINT m_numChrom;

			/// number of loci
			vectoru m_numLoci;

			/// whether or not the last chromosome is sex chromosome
			bool m_sexChrom;

			/// position of loci on chromosome, recommended with unit cM
			vectorf m_lociPos;

			/// loci index
			vectoru m_chromIndex;

			/// allele names
			vectorstr m_alleleNames;

			/// loci names
			vectorstr m_lociNames;

			/// max allele
			UINT m_maxAllele;

			/// name of the information field
			vectorstr m_infoFields;

			/// chromosome map for mpi modules
			/// This field is not saved/restored
			vectori m_chromMap;

#ifdef SIMUMPI
			/// begin chromosome for this node
			/// This field is not saved/restored
			UINT m_beginChrom;

			/// end chromosome for this node
			/// This field is not saved/restored
			UINT m_endChrom;

			/// begin locus for this node
			UINT m_beginLocus;

			/// end locus for this node
			UINT m_endLocus;

			/// total number of loci at this node
			UINT m_localNumLoci;

			/// local genotype size
			UINT m_localGenoSize;

			/// loci index
			vectoru m_localChromIndex;
#endif

			friend class GenoStruTrait;
	};
}



#ifndef SWIG
// set version for GenoStructure class
// version 0: base
// version 1: add sexChrom indicator
// version 2: add infoSize
BOOST_CLASS_VERSION(simuPOP::GenoStructure, 2)
#endif

namespace simuPOP
{

	/** \brief genoStruTrait

	A trait class that maintain a static array of geno structure,
	and provide interfaces around a GenoStructure Index.
	*/
	class GenoStruTrait
	{
		private:

#define TraitIndexType unsigned char
#define TraitMaxIndex 0xFF

		public:
			/// constructor, but m_genoStruIdx will be set later.
			GenoStruTrait():m_genoStruIdx(TraitMaxIndex)
			{
			}

			/// set genotypic structure
			/// CPPONLY
			void setGenoStructure(UINT ploidy, const vectoru& loci, bool sexChrom,
				const vectorf& lociPos, const vectorstr& alleleNames,
				const vectorstr& lociNames, UINT maxAllele, const vectorstr& infoFields,
				const vectori& chromMap);

			/// set an existing geno structure, simply use it
			/// This is NOT efficient! (but has to be used when, for example,
			/// loading a structure from file
			void setGenoStructure(GenoStructure& rhs);

			/// CPPONLY set index directly
			void setGenoStruIdx(size_t idx)
			{
				DBG_FAILIF( idx >= s_genoStruRepository.size(), IndexError,
					"Index " + toStr(idx) + " to geno structure repository should be less than " +
					toStr( s_genoStruRepository.size() ) );
				m_genoStruIdx = static_cast<TraitIndexType>(idx);
			}

			/// return the GenoStructure
			/// CPPONLY
			GenoStructure& genoStru() const
			{
				return s_genoStruRepository[m_genoStruIdx];
			}

			/// return the GenoStructure index
			/// CPPONLY
			size_t genoStruIdx() const
			{
				return static_cast<size_t>(m_genoStruIdx);
			}

			/// return ploidy
			UINT ploidy() const
			{

				DBG_FAILIF( m_genoStruIdx == TraitMaxIndex, SystemError,
					"Ploidy: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

				return s_genoStruRepository[m_genoStruIdx].m_ploidy;
			}

			/// return ploidy
			string ploidyName() const;

			/// number of loci on chromosome \c chrom
			UINT numLoci(UINT chrom) const
			{
				DBG_FAILIF( m_genoStruIdx == TraitMaxIndex, SystemError,
					"numLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

				CHECKRANGECHROM(chrom);
				return s_genoStruRepository[m_genoStruIdx].m_numLoci[chrom];
			}

			/// whether or not the last chromosome is sex chromosome
			bool sexChrom() const
			{
				DBG_FAILIF( m_genoStruIdx == TraitMaxIndex, SystemError,
					"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

				return s_genoStruRepository[m_genoStruIdx].m_sexChrom;
			}

			/// return totNumLoci (STATIC)
			UINT totNumLoci() const
			{

				DBG_FAILIF( m_genoStruIdx == TraitMaxIndex, SystemError,
					"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

				return s_genoStruRepository[m_genoStruIdx].m_totNumLoci;
			}

			/// return totNumLoci * ploidy
			UINT genoSize() const
			{
				DBG_FAILIF( m_genoStruIdx == TraitMaxIndex, SystemError,
					"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

				return s_genoStruRepository[m_genoStruIdx].m_genoSize;
			}

			/// locus distance.
			double locusPos(UINT locus) const
			{
				DBG_FAILIF( m_genoStruIdx == TraitMaxIndex, SystemError,
					"locusPos: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

				CHECKRANGEABSLOCUS(locus);
				return s_genoStruRepository[m_genoStruIdx].m_lociPos[locus];
			}

			/// expose loci distance
			PyObject* arrLociPos()
			{
				return Double_Vec_As_NumArray( s_genoStruRepository[m_genoStruIdx].m_lociPos.begin(),
					s_genoStruRepository[m_genoStruIdx].m_lociPos.end() );
			}

			/// expose loci distance of a chromosome
			PyObject* arrLociPos(UINT chrom)
			{
				CHECKRANGECHROM(chrom);

				return Double_Vec_As_NumArray(
					s_genoStruRepository[m_genoStruIdx].m_lociPos.begin() + chromBegin(chrom),
					s_genoStruRepository[m_genoStruIdx].m_lociPos.begin() + chromEnd(chrom) );
			}

			/// number of chromosome
			UINT numChrom() const
			{
				DBG_FAILIF( m_genoStruIdx == TraitMaxIndex, SystemError,
					"numChrom: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

				return s_genoStruRepository[m_genoStruIdx].m_numChrom;
			}

			/// chromosome index
			const vectoru& chromIndex() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_chromIndex;
			}

			/// chromosome index of chromosome \c chrom
			UINT chromBegin(UINT chrom) const
			{
				DBG_FAILIF( m_genoStruIdx == TraitMaxIndex, SystemError,
					"chromBegin: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

				CHECKRANGECHROM(chrom);

				return s_genoStruRepository[m_genoStruIdx].m_chromIndex[chrom];
			}

			/// chromosome index of chromosome \c chrom
			UINT chromEnd(UINT chrom) const
			{
				DBG_FAILIF( m_genoStruIdx == TraitMaxIndex, SystemError,
					"chromEnd: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

				CHECKRANGECHROM(chrom);

				return s_genoStruRepository[m_genoStruIdx].m_chromIndex[chrom+1];
			}

			/// convert from relative locus (on chromsome) to absolute locus (no chromosome structure)
			UINT absLocusIndex(UINT chrom, UINT locus)
			{
				CHECKRANGECHROM(chrom);
				CHECKRANGELOCUS(chrom, locus);

				return( s_genoStruRepository[m_genoStruIdx].m_chromIndex[chrom] + locus );
			}

			/// return chrom, locus pair from an absolute locus position.
			std::pair<UINT, UINT> chromLocusPair(UINT locus) const;

			/// return allele name
			string alleleName(const Allele allele) const;

			/// allele names
			vectorstr alleleNames() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_alleleNames;
			}

			/// return locus name
			string locusName(const UINT loc) const
			{
				DBG_FAILIF( loc >= s_genoStruRepository[m_genoStruIdx].m_totNumLoci, IndexError,
					"Locus index " + toStr(loc) + " out of range of 0 ~ " +
					toStr(s_genoStruRepository[m_genoStruIdx].m_totNumLoci));

				return s_genoStruRepository[m_genoStruIdx].m_lociNames[loc];
			}

			UINT maxAllele() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_maxAllele;
			}

			void setMaxAllele(UINT maxAllele)
			{
#ifdef BINARYALLELE
				DBG_ASSERT(maxAllele == 1,  ValueError,
					"max allele must be 1 for binary modules");
#else
				s_genoStruRepository[m_genoStruIdx].m_maxAllele = maxAllele;
#endif
			}

			/// get info length
#ifdef SIMUMPI
			UINT localInfoSize() const;
			UINT infoSize() const;
			vectorstr infoFields() const;
			string infoField(UINT idx) const;
			UINT infoIdx(const string& name) const;
			void struAddInfoField(const string& field);
			void struSetInfoFields(const vectorstr& fields);
#else
			UINT infoSize() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_infoFields.size();
			}

			vectorstr infoFields() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_infoFields;
			}

			string infoField(UINT idx) const
			{
				CHECKRANGEINFO(idx);
				return s_genoStruRepository[m_genoStruIdx].m_infoFields[idx];
			}

			/// return the index of field name, return -1 if not found.
			UINT infoIdx(const string& name) const;

			/// add a new information field
			/// NOTE: should only be called by population::requestInfoField
			/// Right now, do not allow dynamic addition of these fields.
			/// CPPONLY
			void struAddInfoField(const string& field)
			{
				vectorstr& fields = s_genoStruRepository[m_genoStruIdx].m_infoFields;
				fields.push_back(field);
				return fields.size()-1;
			}

			/// should should only be called from population
			/// CPPONLY
			void struSetInfoFields(const vectorstr& fields)
			{
				s_genoStruRepository[m_genoStruIdx].m_infoFields = fields;
			}
#endif

			void swap(GenoStruTrait& rhs)
			{
				std::swap(m_genoStruIdx, rhs.m_genoStruIdx);
			}

			///
			vectori chromMap() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_chromMap;
			}

#ifdef SIMUMPI
			/// return node rank by chromosome number, according to map on setChromMap
			UINT rankOfChrom(UINT chrom) const;

			/// return node rank by locus id, according to map on setChromMap
			UINT rankOfLocus(UINT locus) const
			{
				return rankOfChrom(chromLocusPair(locus).first);
			}

			/// begin chromosome for a given rank
			UINT beginChromOfRank(UINT rank) const;

			/// end chromosome for a given rank (actually begin chromosome for the next rank)
			UINT endChromOfRank(UINT rank) const;

			/// begin locus for a given rank
			UINT beginLocusOfRank(UINT rank) const
			{
				return chromBegin(beginChromOfRank(rank));
			}

			/// end locus for a given rank
			UINT endLocusOfRank(UINT rank) const
			{
				return rank==0?0:chromEnd(endChromOfRank(rank)-1);
			}

			/// begin chromosome for current node
			UINT beginChrom() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_beginChrom;
			}

			/// end chromosome for current node
			UINT endChrom() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_endChrom;
			}

			UINT localChromBegin(UINT chrom) const
			{
				return s_genoStruRepository[m_genoStruIdx].m_localChromIndex[chrom-beginChrom()];
			}

			UINT localChromEnd(UINT chrom) const
			{
				return s_genoStruRepository[m_genoStruIdx].m_localChromIndex[chrom+1-beginChrom()];
			}

			/// begin locus of current rank
			UINT beginLocus() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_beginLocus;
			}

			/// end locus of current rank
			UINT endLocus() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_endLocus;
			}

			UINT localNumLoci() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_localNumLoci;
			}

			UINT localGenoSize() const
			{
				return s_genoStruRepository[m_genoStruIdx].m_localGenoSize;
			}

			bool hasChrom(UINT chrom) const
			{
				return mpiRank() == rankOfChrom(chrom);
			}

			bool hasLocus(UINT locus) const
			{
				return mpiRank() == rankOfLocus(locus);
			}
#endif

		private:

			friend class boost::serialization::access;

			template<class Archive>
				void serialize(Archive & ar, const UINT version)
			{
				// do not archive index.
			}

		private:
			/// m_genoStru is originally a pointer,
			/// I am using a short index now to save a few RAM (4 vs 1)
			/// This may become significant since this info is avaiable for
			/// all individuals.
			TraitIndexType m_genoStruIdx;

			/// glocal genotypic strcuture repository
			/// only unique structure will be saved
			/// store pointers instead of object to avoid relocation of
			/// objects themselves by vector
			static vector<GenoStructure> s_genoStruRepository;
	};
}



#ifndef SWIG
BOOST_CLASS_TRACKING(simuPOP::GenoStruTrait, track_never)
#endif
#endif
