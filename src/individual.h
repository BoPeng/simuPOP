/***************************************************************************
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
#if  defined (_WIN32) || defined (__WIN32__)
#  include <boost/archive/binary_iarchive.hpp>
#  include <boost/archive/binary_oarchive.hpp>
#  include <fstream>
using std::ofstream;
using std::ifstream;
#endif                                                                                    // win32

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

namespace simuPOP {
/// individuals with genotype, affection status, sex etc.
/**
   Individuals are the building blocks of populations, each having
   the following individual information:
 \li shared genotypic structure information
 \li genotype
 \li sex, affection status, subpopulation ID
 \li optional information fields

   Individual genotypes are arranged by locus, chromosome, ploidy, in that order,
   and can be accessed from a single index. For example, for a diploid individual with
   two loci on the first chromosome, one locus on the second, its genotype is arranged
   as <tt> 1-1-1 1-1-2 1-2-1 2-1-1 2-1-2 2-2-1 </tt> where \c x-y-z represents ploidy \c x
   chromosome \c y and locus \c z. An allele \c 2-1-2 can be accessed by
 \c allele(4) (by absolute index), \c allele(1, 1) (by index and ploidy) or \c allele(1, 1, 0)
   (by index, ploidy and chromosome). Individuals are created by populations automatically.
   Do not call the constructor function directly.
 */
/*
   Usage information: (for population class developers)
 \li for individuals created, you are responsible for setting their genotypic
   pointer and genotypic information by using
   <tt>setGenoStructure(GenoStructure gs)</tt>.
 \li \c setSubPopID() and \c subPopID() can be used for any \em temporary purpose.

 \note
 \li \c individual does \em not manage memory. Instead, it use a pointer passed
   from \c population class. This may cause a \em lot of troubles.
 \li \c operator= uses shallow copy. This is required by \em sort algorithm since
   otherwise individuals will not be able to be copied. However, in population
   memory management, it is sometimes required that genotypic information within
   one subpopulation should go together. This is done by using a \c shollow_copied
   flag for each individual and for all individuals. Population might have to rearrange
   individuals to solve this problem.
 \li Output of \c individual can be adjusted by \c setOutputDelimeter.
 */
class individual : public GenoStruTrait
{

protected:
	/// 0: male, 1: female regardless of outside coding
	static const size_t m_flagFemale = 1;

	/// if this individual is affect
	static const size_t m_flagAffected = 2;

	/// if this individual is the result of a shoallow copy
	static const size_t m_flagShallowCopied = 4;

	/// if this individual is visible. This is used
	/// to implement virtual subpopulations
	static const size_t m_flagVisible = 8;

public:
	///  @name constructor, destructor etc
	//@{
	///
	/**
	 \test src_individual.log Individual member functions
	 */
	individual() : m_flags(m_flagVisible), m_subPopID(0)
	{
	}


	/// CPPONLY copy constructor will be a shallow copied one
	individual(const individual & ind) :
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


	/// CPPONLY set genotype pointer (use if Allele*pos can not be determined during construction)
	void setGenoPtr(GenoIterator pos)
	{
		m_genoPtr = pos;
	}


	/// CPPONLY set pointer to individual info
	void setInfoPtr(InfoIterator pos)
	{
		m_infoPtr = pos;
	}


	/// shallow copy of an individual class
	individual & operator=(const individual & rhs);

	/// CPPONLY deep copy of an individual class
	individual & copyFrom(const individual & rhs);

	//@}
	/// @name readonly structural info
	//@{

	/// CPPONLY pointer to alleles
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

	/// return an editable array (a carray of length <tt>totNumLoci()*ploidy()</tt>) of genotypes of an individual
	/**
	   This function returns the whole genotype. Although this function is
	   not as easy to use as other functions that access alleles,
	   it is the fastest one since you can read/write genotype directly.
	 */
	PyObject * arrGenotype();

	/// return a carray with the genotypes of the \c p-th copy of the chromosomes
	PyObject * arrGenotype(UINT p);

	/// return a carray with the genotypes of the \c ch-th chromosome in the \c p-th chromosome set
	PyObject * arrGenotype(UINT p, UINT ch);

	/// return a carray of all information fields (of size \c infoSize()) of this individual
	PyObject * arrInfo();

	/// return the allele at locus \c index
	/**
	 \param index absolute index from the beginning of the genotype, ranging from \c 0
	   	to <tt> totNumLoci()*ploidy() </tt>
	 */
#ifdef SIMUMPI
	Allele allele(UINT index) const;

#else
	Allele allele(UINT index) const
	{
		CHECKRANGEGENOSIZE(index);
		return *(m_genoPtr + index);
	}


#endif

	/// return the allele at locus \c index of the \c p-th copy of the chromosomes
	/**
	 \param index index from the begining of the \c p-th set of the chromosomes, ranging from
	 \c 0 to <tt> totNumLoci() </tt>
	 \param p index of the ploidy
	 */
#ifdef SIMUMPI
	Allele allele(UINT index, UINT p) const;

#else
	Allele allele(UINT index, UINT p) const
	{
		CHECKRANGEABSLOCUS(index);
		CHECKRANGEPLOIDY(p);
		return *(m_genoPtr + index + p * totNumLoci() );
	}


#endif

	/// return the allele at locus \c index of the \c ch-th chromosome in the \c p-th chromosome set
	/**
	 \param index index from the begining of chromosome \c ch of ploidy \c p,
	   	ranging from \c 0 to <tt> numLoci(ch) </tt>
	 \param p index of the polidy
	 \param ch index of the chromosome in the \c p-th chromosome set
	 */
#ifdef SIMUMPI
	Allele allele(UINT index, UINT p, UINT ch) const;

#else
	Allele allele(UINT index, UINT p, UINT ch) const
	{
		CHECKRANGELOCUS(ch, index);
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(ch);
		return *(m_genoPtr + index + p * totNumLoci() + chromBegin(ch));
	}


#endif

	/// return the name of \c allele(index)
	string alleleChar(UINT index) const
	{
		CHECKRANGEGENOSIZE(index);

		return this->alleleName(allele(index));
	}


	/// return the name of <tt>allele(index, p)</tt>
	string alleleChar(UINT index, UINT p) const
	{
		CHECKRANGEABSLOCUS(index);
		CHECKRANGEPLOIDY(p);

		return this->alleleName(allele(index, p));
	}


	/// return the name of <tt>allele(idx, p, ch)</tt>
	string alleleChar(UINT index, UINT p, UINT ch) const
	{
		CHECKRANGELOCUS(ch, index);
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(ch);

		return this->alleleName(allele(index, p, ch));
	}


	/// set the allele at locus \c index
	/**
	 \param allele allele to be set
	 \param index index from the begining of genotype, ranging from \c 0
	   	to <tt> totNumLoci()*ploidy() </tt>
	 */
#ifdef SIMUMPI
	void setAllele(Allele allele, UINT index);

#else
	void setAllele(Allele allele, UINT index)
	{
		CHECKRANGEGENOSIZE(index);
		*(m_genoPtr + index) = allele;
	}


#endif

	/// set the allele at locus \c index of the \c p-th copy of the chromosomes
	/**
	 \param allele allele to be set
	 \param index index from the begining of the poloidy \c p, ranging from \c 0 to <tt> totNumLoci(p) </tt>
	 \param p index of the poloidy
	 */
#ifdef SIMUMPI
	void setAllele(Allele allele, UINT index, UINT p);

#else
	void setAllele(Allele allele, UINT index, UINT p)
	{
		CHECKRANGEABSLOCUS(index);
		CHECKRANGEPLOIDY(p);
		*(m_genoPtr + index + p * totNumLoci()) = allele;
	}


#endif

	/// set the allele at locus \c index of the \c ch-th chromosome in the \c p-th chromosome set
	/**
	 \param allele allele to be set
	 \param index index from the begining of the chromosome, ranging from \c 0 to \c numLoci(ch)
	 \param p index of the ploidy
	 \param ch index of the chromosome in ploidy \c p
	 */
#ifdef SIMUMPI
	void setAllele(Allele allele, UINT index, UINT p, UINT ch);

#else
	void setAllele(Allele allele, UINT index, UINT p, UINT ch)
	{
		CHECKRANGELOCUS(ch, index);
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(ch);
		*(m_genoPtr + index + p * totNumLoci() + chromBegin(ch) ) = allele;
	}


#endif

	/// return the sex of an individual, \c 1 for males and \c 2 for females.
	Sex sex() const
	{
		if (ISSETFLAG(m_flags, m_flagFemale) )
			return Female;
		else
			return Male;
	}


	/// return the sex of an individual, \c M or \c F
	char sexChar() const
	{
		return sex() == Female ? 'F' : 'M';
	}


	/// set sex. \c sex can be \c Male of \c Female.
	void setSex(Sex sex)
	{
		CHECKRANGESEX(sex);

		if (sex == Male)
			RESETFLAG(m_flags, m_flagFemale);
		else
			SETFLAG(m_flags, m_flagFemale);
	}


	/// whether or not an individual is affected
	bool affected() const
	{
		return ISSETFLAG(m_flags, m_flagAffected);
	}


	/// equals to <tt>not affected()</tt>
	bool unaffected() const
	{
		return !affected();
	}


	/// return \c A (affected) or \c U (unaffected) for affection status
	char affectedChar() const
	{
		return affected() ? 'A' : 'U';
	}


	/// set affection status
	void setAffected(bool affected)
	{
		if (affected)
			SETFLAG(m_flags, m_flagAffected);
		else
			RESETFLAG(m_flags, m_flagAffected);
	}


	bool visible() const
	{
		return ISSETFLAG(m_flags, m_flagVisible);
	}


	void setVisible(bool visible)
	{
		if (visible)
			SETFLAG(m_flags, m_flagVisible);
		else
			RESETFLAG(m_flags, m_flagVisible);
	}


	/// return the ID of the subpopulation to which this individual blongs
	/**
	 \note \c subPopID is not set by default. It only corresponds to the subpopulation
	   	in which this individual resides after \c pop::setIndSubPopID is called.
	 */
	SubPopID subPopID() const
	{
		return m_subPopID;
	}


	/// set new subpopulation ID, \c pop.rearrangeByIndID will move this individual to that population
	void setSubPopID(SubPopID id)
	{
		m_subPopID = id;
	}


	/// get information field \c idx
	/**
	 \param idx index of the information field
	 */
	InfoType info(UINT idx) const
	{
		CHECKRANGEINFO(idx);

		return m_infoPtr[idx];
	}


	/// get information field \c name
	/**
	   Equivalent to <tt>info(infoIdx(name))</tt>.
	 \param name name of the information field
	 */
	InfoType info(const string & name) const
	{
		int idx = infoIdx(name);

		DBG_ASSERT(idx >= 0, IndexError,
		    "Info name " + name + " is not a valid info field name");
		return m_infoPtr[idx];
	}


	/// set information field by \c idx
	void setInfo(InfoType value, UINT idx)
	{
		CHECKRANGEINFO(idx);
		m_infoPtr[idx] = value;
	}


	/// set information field by \c name
	void setInfo(InfoType value, const string & name)
	{
		int idx = infoIdx(name);

		DBG_ASSERT(idx >= 0, IndexError,
		    "Info name " + name + " is not a valid info field name");
		m_infoPtr[idx] = value;
	}


	/// CPPONLY start of alleles
	GenoIterator genoBegin() const
	{
		return m_genoPtr;
	}


	/// CPPONLY end of allele
	GenoIterator genoEnd() const
	{
#ifdef SIMUMPI
		return m_genoPtr + localGenoSize();
#else
		return m_genoPtr + genoSize();
#endif
	}


	/// CPPONLY start of allele of the pth set of chromosome
	GenoIterator genoBegin(UINT p) const
	{
		CHECKRANGEPLOIDY(p);
#ifdef SIMUMPI
		return m_genoPtr + p * localNumLoci();

#else
		return m_genoPtr + p * totNumLoci();

#endif
	}


	/// CPPONLY end of allele of the pth set of chromosome
	GenoIterator genoEnd(UINT p) const
	{
		CHECKRANGEPLOIDY(p);
#ifdef SIMUMPI
		return m_genoPtr + (p + 1) * localNumLoci();
#else
		return m_genoPtr + (p + 1) * totNumLoci();
#endif
	}


	/// CPPONLY start of allele of the pth set of chromosome, chrom ch
	GenoIterator genoBegin(UINT p, UINT chrom) const
	{
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(chrom);
		return m_genoPtr + p * totNumLoci() + chromBegin(chrom);

	}


	/// CPPONLY end of allele of the pth set of chromosome
	GenoIterator genoEnd(UINT p, UINT chrom) const
	{
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(chrom);
		return m_genoPtr + p * totNumLoci() + chromEnd(chrom);

	}


	/// CPPONLY start of info
	InfoIterator infoBegin() const
	{
		return m_infoPtr;
	}


	/// CPPONLY end of info
	InfoIterator infoEnd() const
	{
		return m_infoPtr + infoSize();
	}


	//@}
	/// @name copy, comparison, swap operations to objects.
	//@{
	/// compare if two individuals are the same used in case of serialization etc.
	/**
	 \note We do not compare info because \c m_subPopID is considered temporary.
	 */
	bool operator==(const individual & rhs) const;

	/// compare if two individuals are not the same used in case of serialization etc.
	bool operator!=(const individual & rhs) const
	{
		return !(*this == rhs);
	}


	// allow compaison of individuals in python
	// only equal or unequal, no greater or less than
	/// a python function used to compare the individual objects
	int __cmp__(const individual & rhs) const;

	/**
	   There is usally no '>', '<' comparisons for individuals.
	   If order is required, it is a comparison of \c info.
	   This behavior is used in migration.
	 */
	bool operator<(const individual & rhs) const
	{
		return subPopID() < rhs.subPopID();
	}


	// allow str(population) to get something better looking
	/// used by Python print function to print out the general information of the individual
	string __repr__();

	/// CPPONLY swap individuals
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
	   subpopulation so the later case is rare, at best.)
	 */
	void swap(individual & ind, bool swapContent = true);

	//@}
	/// @name misc (only relevant to developers)
	//@{

	/// CPPONLY is this individual a result of shallow copy?
	bool shallowCopied() const
	{
		return ISSETFLAG(m_flags, m_flagShallowCopied);
	}


	/// CPPONLY set shallowCopied flag
	void setShallowCopied(bool shallowCopied)
	{
		if (shallowCopied)
			SETFLAG(m_flags, m_flagShallowCopied);
		else
			RESETFLAG(m_flags, m_flagShallowCopied);
	}


	/// CPPONLY
	void display(ostream & out, int width = 1, const vectori & chrom = vectori(), const vectori & loci = vectori() );

	//@}

private:
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const UINT version) const
	{
		// ar & boost::serialization::make_nvp("base ptr",
		//  boost::serialization::base_object<GenoStruTrait>(*this));
		bool b;

		b = ISSETFLAG(m_flags, m_flagFemale);
		ar & boost::serialization::make_nvp("sex", b);

		b = ISSETFLAG(m_flags, m_flagAffected);
		ar & boost::serialization::make_nvp("affected", b);
	}


	template<class Archive>
	void load(Archive & ar, const UINT version)
	{
		bool b;

		m_flags = 0;
		ar & boost::serialization::make_nvp("sex", b);
		if (b) SETFLAG(m_flags, m_flagFemale);
		ar & boost::serialization::make_nvp("affected", b);
		if (b) SETFLAG(m_flags, m_flagAffected);

		RESETFLAG(m_flags, m_flagShallowCopied);

		if (version < 1) {
			std::pair<int, int> tag;
			ar & make_nvp("tag", tag);
			ar & make_nvp("info", m_subPopID);
		}
	}


	BOOST_SERIALIZATION_SPLIT_MEMBER();

protected:
	// internal flag. Can be used to perform many things.
	// bitset<3> was previously used but that will take 4 bytes.
	unsigned char m_flags;

	/// temporary information
	SubPopID m_subPopID;

#ifdef SIMUMPI
	/// pointer to genotype.
	union {
		// used by slave node
		GenoIterator m_genoPtr;
		// used by master node
		ULONG m_indIndex;
	};

	// population ID.
	ULONG m_popID;
#else
	/// pointer to genotype.
	GenoIterator m_genoPtr;
#endif

	/// pointer to info
	InfoIterator m_infoPtr;
};


/**
   	this class implements a Python itertor class that can be used to iterate
   	through individuals in a (sub)population. If allInds are true,
   	visiblility of individuals will not be checked. Note that
   	individualIterator *will* iterate through only visible individuals, and
   	allInds is only provided when we know in advance that all individuals are
   	visible. This is a way to obtain better performance in simple cases.

   	An instance of this class is returned by
   	population::individuals() and population::individuals(subPop)
 */
class pyIndIterator
{
public:
	pyIndIterator(vector<individual>::iterator const begin,
	              vector<individual>::iterator const end, bool allInds) :
		m_index(begin),
		m_end(end),
		m_allInds(allInds)
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


	individual & next();

private:
	// current (initial individual)
	vector<individual>::iterator m_index;

	// ending index
	vector<individual>::iterator m_end;

	//
	bool m_allInds;
};

/**
   	this class implements a C++ iterator class that iterate through
   	individuals in a (sub)population. If allInds are true, the
   	visiblility of individuals will not be checked. Note that
   	individualIterator *will* iterate through only visible individuals, and
   	allInds is only provided when we know in advance that all individuals are
   	visible. This is a way to obtain better performance in simple cases.
 */
template <typename T>
class IndividualIterator
{
public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef typename T::value_type value_type;
	typedef long int difference_type;
	typedef typename T::reference reference;
	typedef typename T::pointer pointer;

	IndividualIterator() : m_it(), m_end(), m_allInds(true)
	{
	}


	IndividualIterator(T it, T end, bool allInds = true)
		: m_it(it), m_end(end), m_allInds(allInds)
	{
		// m_it need to point to the first valid 
		// individual. otherwise *it will fail.
		while (!m_it->visible() && m_it < m_end)
			++m_it;
	}


	bool valid()
	{
		return m_it < m_end;
	}


	T rawIter()
	{
		return m_it;
	}


	reference operator*() const
	{
		return *m_it;
	}


	pointer operator->() const
	{
		return & * m_it;
	}


	// return, then advance.
	IndividualIterator operator++(int)
	{
		DBG_ASSERT(m_it < m_end, ValueError,
		    "Can not advance invalid iterator");

		if (m_allInds)
			return IndividualIterator(m_it++, m_end, m_allInds);

		// save current state
		IndividualIterator tmp(*this);
		// move forward
		while (m_it < m_end)
			if ((++m_it)->visible())
				break;
		// return the original one
		return tmp;
	}


	IndividualIterator operator++()
	{
		DBG_ASSERT(m_it < m_end, ValueError,
		    "Can not advance invalid iterator");
		if (m_allInds) {
			++m_it;
			return *this;
		}

		while (m_it < m_end)
			if ((++m_it)->visible())
				return *this;
		DBG_ASSERT(m_it == m_end, IndexError,
		    "Something wrong with operator++ here");
		return *this;
	}


	// NOTE: This is very slow, and it is here only
	// because sort() requires this. This means that
	// sorting in virtual populations will be slow
	//
	IndividualIterator operator+(difference_type diff)
	{
		if (!m_allInds)
			return IndividualIterator(m_it + diff, m_end, m_allInds);
		IndividualIterator tmp(*this);
		DBG_ASSERT(tmp.m_it < tmp.m_end, ValueError,
		    "Can not advance invalid iterator");
		difference_type i = 0;
		while (i < diff && tmp.m_it < tmp.m_end)
			if ((++tmp.m_it)->visible())
				++i;
		DBG_FAILIF(i != diff, ValueError,
		    "Can not add to IndIterator");
		return tmp;
	}


	IndividualIterator operator-(difference_type diff)
	{
		if (m_allInds)
			return IndividualIterator(m_it - diff, m_end, m_allInds);
		else {
			IndividualIterator tmp(*this);
			// can not check. Possible problem
			for (difference_type i = 0; i < diff; ++i)
				while (!(--tmp.m_it)->visible()) ;
			return tmp;
		}
	}


	difference_type operator-(IndividualIterator rhs)
	{
		if (m_allInds)
			return m_it - rhs.m_it;
		else {
			difference_type i = 0;
			for (T it = rhs.m_it; it != m_it; ++it)
				if (it->visible())
					++i;
			return i;
		}
	}


	IndividualIterator operator--(int)
	{
		if (m_allInds)
			return IndividualIterator(m_it--, m_end, m_allInds);
		IndividualIterator tmp(*this);
		while (!(--m_it)->visible()) ;
		return tmp;
	}


	IndividualIterator operator--()
	{
		if (m_allInds) {
			--m_it;
			return *this;
		}
		while (!(--m_it)->visible()) ;
		return *this;
	}


	bool operator!=(const IndividualIterator & rhs)
	{
		return m_it != rhs.m_it;
	}


	bool operator==(const IndividualIterator & rhs)
	{
		return m_it == rhs.m_it;
	}


	bool operator<(const IndividualIterator & rhs)
	{
		return m_it < rhs.m_it;
	}


	bool operator>(const IndividualIterator & rhs)
	{
		return m_it > rhs.m_it;
	}


private:
	/// The current individual iterator
	T m_it;

	/// the ending iterator. Used to make sure that m_it will not
	/// go past m_end beceause m_it->visible() will be invalid otherwise.
	T m_end;

	// a shortcut. If m_allInds is set, using a simpler algorithm.
	bool m_allInds;
};

//
typedef vector<individual>::iterator RawIndIterator;
typedef vector<individual>::const_iterator ConstRawIndIterator;

typedef IndividualIterator<RawIndIterator> IndIterator;
typedef IndividualIterator<ConstRawIndIterator> ConstIndIterator;

}


#ifndef SWIG
// set version for individual class
// version 0: base
// version 1: add sexChrom indicator
// version 2: add infoSize
BOOST_CLASS_VERSION(simuPOP::individual, 1)
#endif

#endif
