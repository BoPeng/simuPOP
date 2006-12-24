/*************************enumerate(pop.individuals()):
# find spose.**************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu                                                        *
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$
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

#ifndef _INDIVIDUAL_H
#define _INDIVIDUAL_H

/**
\file
\brief class individual, individualWithAge etc.
*/

#include "utility.h"
#include "simupop_cfg.h"

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
				m_alleleNames(), m_lociNames(), m_maxAllele(), m_infoFields(0)
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
				const vectorstr& lociNames, UINT maxAllele, const vectorstr& infoFields);

			/// copy constructor
			/// CPPONLY
			GenoStructure(const GenoStructure& rhs) :
			m_ploidy(rhs.m_ploidy),
				m_totNumLoci(rhs.m_totNumLoci),
				m_genoSize(rhs.m_genoSize),
				m_numChrom(rhs.m_numChrom),
				m_numLoci(rhs.m_numLoci),
				m_sexChrom(rhs.m_sexChrom),
				m_lociPos(rhs.m_lociPos),
				m_chromIndex(rhs.m_chromIndex),
				m_alleleNames(rhs.m_alleleNames),
				m_lociNames(rhs.m_lociNames),
				m_maxAllele(rhs.m_maxAllele),
				m_infoFields(rhs.m_infoFields)
			{
			}

			bool operator== (const GenoStructure& rhs)
			{
				// compare pointer directly will be fastest
				if(this == &rhs || (
					( m_ploidy == rhs.m_ploidy) &&
					( m_numLoci == rhs.m_numLoci) &&
					( m_sexChrom == rhs.m_sexChrom) &&
					( m_lociPos == rhs.m_lociPos) &&
					( m_alleleNames == rhs.m_alleleNames) &&
					( m_lociNames == rhs.m_lociNames) &&
					( m_maxAllele == rhs.m_maxAllele) &&
					( m_infoFields == rhs.m_infoFields) ))
					return true;
				else
					return false;
			}

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
				ar & make_nvp("ploidy", m_ploidy);
				ar & make_nvp("num_of_chrom", m_numChrom);
				ar & make_nvp("num_of_loci_on_each_chrom", m_numLoci);
				ar & make_nvp("sex_chromosome", m_sexChrom);
				ar & make_nvp("loci_distance_on_chrom", m_lociPos);
				ar & make_nvp("allele_name", m_alleleNames);
				ar & make_nvp("loci_name", m_lociNames);
				ar & make_nvp("max_allele", m_maxAllele);
				ar & make_nvp("info_name", m_infoFields);
			}

			template<class Archive>
				void load(Archive &ar, const UINT version)
			{
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
				const vectorstr& lociNames, UINT maxAllele, const vectorstr& infoFields);

			/// set an existing geno structure, simply use it
			/// This is NOT efficient! (but has to be used when, for example,
			/// loading a structure from file
			void setGenoStructure(GenoStructure& rhs)
			{
				for(TraitIndexType it = 0; it < s_genoStruRepository.size();
					++it)
				{
												  // object comparison
					if( s_genoStruRepository[it] == rhs )
					{
						m_genoStruIdx = it;
						return;
					}
				}

				// if not found, make a copy and store it.
				s_genoStruRepository.push_back( rhs );
				m_genoStruIdx = s_genoStruRepository.size() - 1;
			}

			/// CPPONLY set index directly
			void setGenoStruIdx(size_t idx)
			{
				DBG_FAILIF( idx >= s_genoStruRepository.size(), IndexError,
					"Index " + toStr(idx) + " to geno structure repository should be less than " +
					toStr( s_genoStruRepository.size() ) );
				m_genoStruIdx = static_cast<TraitIndexType>(idx);
			}

			/// no destructure since a pointer will be shared by all indiviudals and a population
			/// only population will call destroyGenoStructure in its destructor.
			/// CPPONLY
			// void destroyGenoStructure()
			//{
			//  delete m_genoStruIdx;
			// }

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
			string ploidyName() const
			{
				DBG_FAILIF( m_genoStruIdx == TraitMaxIndex, SystemError,
					"PloidyName: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

				if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 1)
					return "haploid";
				else if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 2)
					return "diploid";
				else if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 3)
					return "triploid";
				else if (s_genoStruRepository[m_genoStruIdx].m_ploidy == 4)
					return "tetraploid";
				else
					return toStr(s_genoStruRepository[m_genoStruIdx].m_ploidy) + "-polid";
			}

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
			std::pair<UINT, UINT> chromLocusPair(UINT locus) const
			{
				CHECKRANGEABSLOCUS(locus);

				pair<UINT, UINT> loc;

				for(UINT i=1, iEnd =numChrom(); i <= iEnd;  ++i)
				{
					if( s_genoStruRepository[m_genoStruIdx].m_chromIndex[i] > locus)
					{
						loc.first = i-1;
						loc.second = locus - s_genoStruRepository[m_genoStruIdx].m_chromIndex[i-1];
						break;
					}
				}
				return loc;
			}

			/// return allele name
			string alleleName(const Allele allele) const
			{
#ifndef BINARYALLELE                
				DBG_FAILIF(allele > s_genoStruRepository[m_genoStruIdx].m_maxAllele,
					IndexError, "Allele out of range of 0 ~ " +
					toStr(s_genoStruRepository[m_genoStruIdx].m_maxAllele));
#endif
				if( allele < s_genoStruRepository[m_genoStruIdx].m_alleleNames.size() )
				{
					DBG_FAILIF( allele >= s_genoStruRepository[m_genoStruIdx].m_alleleNames.size() ,
						IndexError, "No name for allele " + toStr(static_cast<UINT>(allele)));

					return s_genoStruRepository[m_genoStruIdx].m_alleleNames[allele];
				}
				else
					return toStr(static_cast<int>(allele));
			}

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
			UINT infoIdx(const string& name) const
			{
				vectorstr& names = s_genoStruRepository[m_genoStruIdx].m_infoFields;

				for(UINT i=0; i< names.size(); ++i)
				{
					if(names[i] == name)
						return i;
				}
				throw IndexError("Info field '" + name + "' is not found. "
					"Plese use infoFields=['" + name + "'] option of population() during construction\n"
					"or use addInfoField('" + name + "') to add to an existing population.");
				// this should never be reached.
				return 0;
			}

			/// add a new information field
			/// NOTE: should only be called by population::requestInfoField
			/// return the index of the newly added field
			/// Right now, do not allow dynamic addition of these fields.
			/// CPPONLY
			int struAddInfoField(const string& field)
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

			void swap(GenoStruTrait& rhs)
			{
				std::swap(m_genoStruIdx, rhs.m_genoStruIdx);
			}

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

namespace simuPOP
{
	/** \brief Basic individual class

	class individual with
	- genotypic information
	- shared genotypic structure info (through a GenoStructure pointer)
	- flags about sex, affected status
	- an internal info field
	.

	other individuals will be derived from this class, adding age info etc.

	\b Note that
	- individual DOES NOT manage memory. It will use a pointer passed from
	class population. This causes A LOT of trouble and I have not
	evaluated how much benefic I get.
	- operator = uses shallow copy. This is required by
	sort algorithm since otherwise individuals are
	non-copiable. However, in population memory management,
	it is showtimes required that genotypic information within
	one subPop should go together. This is done by
	a shollow_copied flag for each individual and for all
	individuals. population might have to re-arrange
	individuals to solve this problem.
	- output of individual can be adjusted by setOutputDelimeter.
	.

	Usage info: (for population classes developers)
	- for individuals are created, you are responsible to set its genotypic
	pointer and genotypic information. This is done by
	\code
	setGenoStructure(GenoStructure gs)
	\endcode
	- \c setSubPopID() and \c subPopID() can be used for any \e temporary purpose.

	*/
	class individual : public GenoStruTrait
	{

		protected:
			/// 0: male, 1: female regardless of outside coding
			static const size_t m_flagFemale          = 1;

			/// if this individual is affect
			static const size_t m_flagAffected        = 2;

			/// if this individual is the result of a shoallow copy
			static const size_t m_flagShallowCopied   = 4;

		public:

			///  @name constructor, destructor etc
			//@{
			/// default constructor,
			individual():m_flags(0),m_subPopID(0)
			{
			}

			/// CPPONLY
			/// copy constructor will be a shallow copied one
			individual(const individual& ind) :
			GenoStruTrait(ind), m_flags(ind.m_flags),
				m_subPopID(ind.m_subPopID),
				m_genoPtr(ind.m_genoPtr),
				m_infoPtr(ind.m_infoPtr)
			{
				setShallowCopied(true);
			}

			/// destructor. Do nothing.
			~individual()
			{
			}

			/// CPPONLY
			/// set genotype pointer (use if Allele*pos can not
			/// be determined during construction.
			void setGenoPtr(GenoIterator pos)
			{
				m_genoPtr = pos;
			}

			/// CPPONLY
			/// set pointer to individual info
			void setInfoPtr(InfoIterator pos)
			{
				m_infoPtr = pos;
			}

			/// shallow copy of an object.
			individual& operator= (const individual& rhs)
			{
				setShallowCopied(true);

				/// when the source is being moved. Its relative position may change
				/// so it also becomes shallowCopied.
				const_cast<individual&>(rhs).setShallowCopied(true);

				m_flags = rhs.m_flags;
				setSubPopID(rhs.subPopID());
				setGenoPtr(rhs.genoPtr());
				setInfoPtr(rhs.infoPtr());
				// also copy genoStru pointer...
				this->setGenoStruIdx(rhs.genoStruIdx());
				return *this;
			}

			/// Deep copy! Important!
			individual& copyFrom( const individual& rhs)
			{
				m_flags = rhs.m_flags;
				setSubPopID(rhs.subPopID());
				copy(rhs.genoBegin(), rhs.genoEnd(), genoBegin());
				copy(rhs.infoBegin(), rhs.infoEnd(), infoBegin());
				// also copy genoStru pointer...
				this->setGenoStruIdx(rhs.genoStruIdx());
				setShallowCopied(false);
				return *this;
			}

			//@}
			/// @name readonly structural info
			//@{

			/// pointer to alleles
			/// CPPONLY
			GenoIterator genoPtr() const
			{
				return m_genoPtr;
			}

			/// CPPONLY
			InfoIterator infoPtr() const
			{
				return m_infoPtr;
			}

			//@}
			/// @name allele, info get/set functions
			//@{

			/// return genotype as python Numeric.array object
			/// This is the whole genotype (all)
			PyObject* arrGenotype()
			{
				// this &* is to avoid any possible type mismatch thing.
				return Allele_Vec_As_NumArray( m_genoPtr, m_genoPtr + genoSize() );
			}

			/// return genotype as python Numeric.array object
			/// This is the p'th copy of chromosomes
			PyObject* arrGenotype(UINT p)
			{
				CHECKRANGEPLOIDY(p);

				return Allele_Vec_As_NumArray( m_genoPtr + p*totNumLoci(),
					m_genoPtr + (p+1)*totNumLoci() );
			}

			/// return genotype as python Numeric.array object
			/// This is the ch chromosome of the pth copy of chromosome
			PyObject* arrGenotype(UINT p, UINT ch)
			{
				CHECKRANGEPLOIDY(p);

				return Allele_Vec_As_NumArray( m_genoPtr + p*totNumLoci() + chromBegin(ch),
					m_genoPtr + p*totNumLoci() +chromEnd(ch));
			}

			PyObject* arrInfo()
			{
				return Info_Vec_As_NumArray(m_infoPtr, m_infoPtr + infoSize() );
			}

			/// get allele from an index
			/** \param index index from the beginning of genotypic info
			 */
			Allele allele(UINT index) const
			{
				CHECKRANGEGENOSIZE(index);

				return *(m_genoPtr+index);
			}

			/// get allele from an index, on the pth set of chromosome
			/** \param index index from the begining of the p'th set of chromosomes.
				\param p on p'th set of chromosomes, default to 0
			*/
			Allele allele(UINT index, UINT p) const
			{
				CHECKRANGEABSLOCUS(index);
				CHECKRANGEPLOIDY(p);

				return *(m_genoPtr+index + p* totNumLoci() );
			}

			Allele allele(UINT index, UINT p, UINT ch) const
			{
				CHECKRANGELOCUS(ch, index);
				CHECKRANGEPLOIDY(p);
				CHECKRANGECHROM(ch);

				return *(m_genoPtr + index + p* totNumLoci() + chromBegin(ch));
			}

			string alleleChar(UINT index) const
			{
				CHECKRANGEGENOSIZE(index);

				return this->alleleName(*(m_genoPtr + index));
			}

			/// get allele from an index, on the pth set of chromosome
			/** \param index index from the begining of the p'th set of chromosomes.
				\param p on p'th set of chromosomes, p=0 by default
			*/
			string alleleChar(UINT index, UINT p) const
			{
				CHECKRANGEABSLOCUS(index);
				CHECKRANGEPLOIDY(p);

				return this->alleleName(*(m_genoPtr + index + p* totNumLoci() ));
			}

			/// get allele from an index, on the pth set of chromosome
			/** \param index index from the begining of the p'th set of chromosomes.
				\param p on p'th set of chromosomes, p=0 by default
			*/
			string alleleChar(UINT index, UINT p, UINT ch) const
			{
				CHECKRANGELOCUS(ch, index);
				CHECKRANGEPLOIDY(p);
				CHECKRANGECHROM(ch);

				return this->alleleName(*(m_genoPtr + index + p* totNumLoci()
					+ chromBegin(ch) ) );
			}

			/// set allele from an index.
			/** \param index index from the begining of genotype
			 */
			void setAllele(Allele allele, UINT index)
			{
				CHECKRANGEGENOSIZE(index);

				*(m_genoPtr+index) = allele;
			}

			/// set allele from an index.
			/** \param allele allele to set
			\param index index from the begining of genotype
			\param p on p'th set of chromosome, p=0 by default
			 */
			void setAllele(Allele allele, UINT index, UINT p)
			{
				CHECKRANGEABSLOCUS(index);
				CHECKRANGEPLOIDY(p);

				*(m_genoPtr + index+p*totNumLoci()) = allele;
			}

			void setAllele(Allele allele, UINT index, UINT p, UINT ch)
			{
				CHECKRANGELOCUS(ch, index);
				CHECKRANGEPLOIDY(p);
				CHECKRANGECHROM(ch);

				*(m_genoPtr + index + p*totNumLoci() + chromBegin(ch) ) = allele;
			}

			/// sex?
			Sex sex() const
			{
				if( ISSETFLAG(m_flags, m_flagFemale) ) return Female;
				else return Male;
			}

			/// return M or F for sex, for display purpose
			char sexChar() const
			{
				if( ISSETFLAG(m_flags, m_flagFemale) ) return 'F';
				else return 'M';
			}

			/// set sex
			void setSex(Sex sex)
			{
				CHECKRANGESEX(sex);

				if( sex == Male ) RESETFLAG( m_flags, m_flagFemale);
				else SETFLAG(m_flags, m_flagFemale);
			}

			/// affected?
			bool affected() const
			{
				return( ISSETFLAG(m_flags, m_flagAffected));
			}

			/// unaffected?
			bool unaffected() const
			{
				return( ! ISSETFLAG(m_flags, m_flagAffected));
			}

			/// return A or U for affected/Unaffected, for display purpose
			char affectedChar() const
			{
				if( ISSETFLAG( m_flags, m_flagAffected))
					return 'A';
				else
					return 'U';
			}

			/// set affected status
			void setAffected(bool affected)
			{
				if(affected)
					SETFLAG(m_flags, m_flagAffected);
				else
					RESETFLAG(m_flags, m_flagAffected);
			}

			/// get subpop id
			SubPop_ID subPopID() const
			{
				return m_subPopID;
			}

			/// set subpop if
			void setSubPopID(SubPop_ID id)
			{
				m_subPopID = id;
			}

			/// get info
			InfoType info(UINT idx) const
			{
				CHECKRANGEINFO(idx);
				return m_infoPtr[idx];
			}

			/// set info
			void setInfo(InfoType value, UINT idx)
			{
				CHECKRANGEINFO(idx);
				m_infoPtr[idx] = value;
			}

			/// get info
			InfoType info(const string& name) const
			{
				int idx = infoIdx(name);
				DBG_ASSERT(idx>=0, IndexError,
					"Info name " + name + " is not a valid info field name");
				return m_infoPtr[idx];
			}

			/// set info
			void setInfo(InfoType value, const string& name)
			{
				int idx = infoIdx(name);
				DBG_ASSERT(idx>=0, IndexError,
					"Info name " + name + " is not a valid info field name");
				m_infoPtr[idx] = value;
			}

			/// start of alleles
			/// CPPONLY
			GenoIterator genoBegin() const
			{
				return m_genoPtr;
			}

			/// end of allele
			/// CPPONLY
			GenoIterator genoEnd() const
			{
				return m_genoPtr + genoSize();
			}

			/// start of allele of the pth set of chromosome
			/// CPPONLY
			GenoIterator genoBegin(UINT p) const
			{
				CHECKRANGEPLOIDY(p);

				return m_genoPtr + p*totNumLoci();
			}

			/// end of allele of the pth set of chromosome
			/// CPPONLY
			GenoIterator genoEnd(UINT p) const
			{
				CHECKRANGEPLOIDY(p);

				return m_genoPtr + (p+1)*totNumLoci();
			}

			/// start of allele of the pth set of chromosome, chrom ch
			/// CPPONLY
			GenoIterator genoBegin(UINT p, UINT chrom) const
			{
				CHECKRANGEPLOIDY(p);
				CHECKRANGECHROM(chrom);

				return m_genoPtr + p*totNumLoci() + chromBegin(chrom);
			}

			/// end of allele of the pth set of chromosome
			/// CPPONLY
			GenoIterator genoEnd(UINT p, UINT chrom) const
			{
				CHECKRANGEPLOIDY(p);
				CHECKRANGECHROM(chrom);
				return m_genoPtr + p*totNumLoci() + chromEnd(chrom);
			}

			/// start of info
			/// CPPONLY
			InfoIterator infoBegin() const
			{
				return m_infoPtr;
			}

			/// end of info
			/// CPPONLY
			InfoIterator infoEnd() const
			{
				return m_infoPtr + infoSize();
			}

			//@}
			/// @name copy, comparison, swap operations to objects.
			//@{
			/// compare if two individuals are the same used in case of serialization etc
			/** Note that we do not compare info because
			   m_subPopID is considered temporary.
			*/
			bool operator== (const individual& rhs) const
			{
				if( genoStruIdx() != rhs.genoStruIdx() )
					return false;

				if(ISSETFLAG(m_flags, m_flagFemale) != ISSETFLAG(rhs.m_flags, m_flagFemale)
					|| ISSETFLAG(m_flags, m_flagAffected) != ISSETFLAG(rhs.m_flags, m_flagAffected) )
					return false;

				for( UINT i=0, iEnd = infoSize(); i < iEnd;  ++i)
					if( info(i) != rhs.info(i) )
						return false;

				for( UINT i=0, iEnd = genoSize(); i < iEnd;  ++i)
					if( allele(i) != rhs.allele(i) )
						return false;

				return true;
			}

			///
			bool operator!= (const individual& rhs) const
			{
				return ! (*this == rhs);
			}

			// allow compaison of individuals in python
			// only equal or unequal, no greater or less than
			int __cmp__(const individual& rhs) const
			{
				if( genoStruIdx() != rhs.genoStruIdx() )
					return 1;

				if( m_flags != rhs.m_flags )
					return 1;

				for( UINT i=0, iEnd = infoSize(); i < iEnd;  ++i)
					if( info(i) != rhs.info(i) )
						return 1;

				for( UINT i=0, iEnd = genoSize(); i < iEnd;  ++i)
					if( allele(i) != rhs.allele(i) )
						return 1;

				return 0;
			}

			/// there is usally no >, < comparison for individuals
			/// if order is required, it is a comparison of info.
			/// this behavior is used in migration.
			bool operator< (const individual& rhs) const
			{
				return subPopID() < rhs.subPopID();
			}

			// allow str(population) to get something better looking
			string __repr__()
			{
				return "<simuPOP::individual>";
			}

			/// swap individuals
			/**
			The default behavior is swapping all info, but not the
			position of genotypic info. If swapContent is false,
			pointer to genotypic info is swapped instead. This
			will lead to better performance for swapping but
			may affected performance of allele counting.

			\param ind individual to be swapped in
			\param swapContent swapContent or only the pointers.

			The guideline is that if we swap individuals across
			subpopulation, we should swap content. Otherwise,
			swap pointers. (There is no order right now within
			subpopulation so the later case is rare, at best.
			*/
			void swap(individual& ind, bool swapContent=true)
			{
				if( genoStruIdx() != ind.genoStruIdx() )
					throw SystemError("Can only swap individuals with different geno structure.");

				std::swap(m_subPopID, ind.m_subPopID);
				std::swap(m_infoPtr, ind.m_infoPtr);

				if(swapContent)
				{
					Allele tmp;
					for(UINT i=0, iEnd = genoSize(); i < iEnd;  i++)
					{
						tmp = m_genoPtr[i];
						m_genoPtr[i] = ind.m_genoPtr[i];
						ind.m_genoPtr[i] = tmp;
					}
				}
				else
				{
					setShallowCopied(true);
					ind.setShallowCopied(true);
					std::swap(m_genoPtr, ind.m_genoPtr);
				}
			}

			//@}
			/// @name misc (only relevant to developers.
			//@{

			/// is this individual a result of shallow copy?
			/// CPPONLY
			bool shallowCopied() const
			{
				return ISSETFLAG(m_flags, m_flagShallowCopied);
			}

			/// set shallowCopied flag.
			/// CPPONLY
			void setShallowCopied(bool shallowCopied)
			{
				if( shallowCopied )
					SETFLAG(m_flags, m_flagShallowCopied);
				else
					RESETFLAG(m_flags, m_flagShallowCopied);
			}

			/// CPPONLY
			void display( ostream& out, int width=1, const vectori& chrom=vectori(), const vectori& loci=vectori() )
			{
				out << sexChar() << affectedChar() << " ";
				DBG_DO(DBG_POPULATION, out <<  subPopID() << " ");
				for(UINT p=0, pEnd = ploidy(); p < pEnd;  ++p)
				{
					//      copy( genoBegin()+i, genoBegin()+i+totNumLoci(),
					//        std::ostream_iterator<string>(out, outputSeparator()) );
					if(chrom.empty() && loci.empty())
					{
						for(UINT ch=0, chEnd=numChrom(); ch<chEnd; ++ch)
						{
							for(UINT j = 0, jEnd = numLoci(ch); j < jEnd;  ++j)
								out << setw(width) << alleleChar(j, p, ch);
							out << " ";
						}
					}
					else if(! chrom.empty() && loci.empty())
					{
						for(vectori::const_iterator ch=chrom.begin(); ch != chrom.end(); ++ch)
						{
							for(UINT j = 0, jEnd = numLoci(*ch); j < jEnd;  ++j)
								out << setw(width) << alleleChar(j, p, *ch);
							out << " ";
						}
					}
					else if( chrom.empty() && ! loci.empty())
					{
						for(vectori::const_iterator loc=loci.begin(); loc != loci.end(); ++loc)
							out << setw(width) << alleleChar(*loc, p);
						out << " ";
					}
					else						  // both specified
						throw ValueError("Please specify only one of chrom and loci.");

					if( p != pEnd-1)
						out << "| ";
				}
			}

			//@}
		private:

			friend class boost::serialization::access;

			template<class Archive>
				void save(Archive &ar, const UINT version) const
			{
				// ar & boost::serialization::make_nvp("base ptr",
				//  boost::serialization::base_object<GenoStruTrait>(*this));
				bool b;
				b= ISSETFLAG(m_flags, m_flagFemale);
				ar & boost::serialization::make_nvp("sex",b);

				b= ISSETFLAG(m_flags, m_flagAffected);
				ar & boost::serialization::make_nvp("affected",b);
			}

			template<class Archive>
				void load(Archive &ar, const UINT version)
			{
				bool b;
				m_flags = 0;
				ar & boost::serialization::make_nvp("sex",b);
				if(b) SETFLAG(m_flags, m_flagFemale);
				ar & boost::serialization::make_nvp("affected",b);
				if(b) SETFLAG(m_flags, m_flagAffected);

				RESETFLAG(m_flags, m_flagShallowCopied);

				if (version < 1)
				{
					std::pair<int, int> tag;
					ar & make_nvp("tag", tag);
					ar & make_nvp("info", m_subPopID);
				}
			}

			BOOST_SERIALIZATION_SPLIT_MEMBER();

		protected:

			/// internal flag. Can be used to perform many things.
			/// bitset<3> was previously used but that will take 4 bytes.
			unsigned char m_flags;

			/// temporary information
			SubPop_ID m_subPopID;

			/// pointer to genotype.
			GenoIterator m_genoPtr;

			/// pointer to info
			InfoIterator m_infoPtr;
	};

}



#ifndef SWIG
// set version for GenoStructure class
// version 0: base
// version 1: add sexChrom indicator
// version 2: add infoSize
BOOST_CLASS_VERSION(simuPOP::individual, 1)
#endif
#endif
