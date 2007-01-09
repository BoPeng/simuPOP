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
#include "simuPOP_cfg.h"
#include "genoStru.h"

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
			SubPopID subPopID() const
			{
				return m_subPopID;
			}

			/// set subpop if
			void setSubPopID(SubPopID id)
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
			SubPopID m_subPopID;

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
