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

#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>

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
/**
 *  A \c population consists of individuals with the same genotypic structure.
 *  An \c individual object cannot be created independently, but refences to
 *  inidividuals can be retrieved using member functions of a \c population
 *  object. In addition to structural information shared by all individuals in
 *  a population (provided by class \c genoStruTrait), a \c individual class
 *  provides member functions to get and set \e genotype, \e sex, <em>affection
 *  status</em> and <em>information fields</em> of an individual.
 * 
 *  \par
 *  Genotypes of an individual are stored sequentially and can be accessed
 *  locus by locus, or in batch. The alleles are arranged by position,
 *  chromosome and ploidy. That is to say, the first allele on the first
 *  chromosome of the first homologous set is followed by alleles at other loci
 *  on the same chromsome, then markers on the second and later chromosomes,
 *  followed by alleles on the second homologous set of the chromosomes for
 *  a diploid individual. A consequence of this memory layout is that alleles
 *  at the same locus of a non-haploid individual are separated by
 *  <tt>individual::totNumLoci()</tt> loci.
 */
class individual : public GenoStruTrait
{

protected:
	/// 0: male, 1: female regardless of outside coding
	static const size_t m_flagFemale = 1;

	/// if this individual is affect
	static const size_t m_flagAffected = 2;

	/// if this individual is visible. This is used
	/// to implement virtual subpopulations
	static const size_t m_flagVisible = 4;

	/// if this individual is iteratable. This will
	/// not affect how activated virtual subpoulations
	/// behave, but will affect how pyIndIterator
	/// iterate through the population.
	///
	/// In short, this is supposed to be a temporary, light
	/// weight flag that help iterators go through virtual
	/// subpopulation.
	static const size_t m_flagIteratable = 8;

public:
	///  @name constructor, destructor etc
	//@{
	///
	/**
	 * An \c individual object cannot be created directly. It has to be accessed
	 * from a \c population object using functions such as
	 * <tt>population::individual(</tt><em>idx</em><tt>)</tt>.
	 */
	individual() : m_flags(m_flagVisible), m_subPopID(0)
	{
	}


	/// CPPONLY
	individual(const individual & ind) :
		GenoStruTrait(ind), m_flags(ind.m_flags),
		m_subPopID(ind.m_subPopID),
		m_genoPtr(ind.m_genoPtr),
		m_infoPtr(ind.m_infoPtr)
	{
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

	/// HIDDEN return an editable array (a carray of length <tt>totNumLoci()*ploidy()</tt>) of genotypes of an individual
	/**
	   This function returns the whole genotype. Although this function is
	   not as easy to use as other functions that access alleles,
	   it is the fastest one since you can read/write genotype directly.
	 */
	PyObject * arrGenotype();

	/// HIDDEN return a carray with the genotypes of the \c p-th copy of the chromosomes
	PyObject * arrGenotype(UINT p);

	/// HIDDEN return a carray with the genotypes of the \c ch-th chromosome in the \c p-th chromosome set
	PyObject * arrGenotype(UINT p, UINT ch);

	/// HIDDEN return a carray of all information fields (of size \c infoSize()) of this individual
	PyObject * arrInfo();

	/** return the current allele at a locus, using its absolute index \e idx.
	 *
	 * <group>allele</group>
	 */
	UINT allele(UINT idx) const
	{
		CHECKRANGEGENOSIZE(idx);
		return static_cast<UINT>(*(m_genoPtr + idx));
	}


	/** return the current allele at locus \e idx on the <em>p</em>-th set of
	 *  homologous chromosomes.
	 * <group>allele</group>
	 */
	UINT allele(UINT idx, UINT p) const
	{
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(p);
		return static_cast<UINT>(*(m_genoPtr + idx + p * totNumLoci() ));
	}


	/** return the current allele at locus \e idx on chromosome \e chrom of
	 *  the <em>p</em>-th set of homologous chromosomes.
	 * <group>allele</group>
	 */
	UINT allele(UINT idx, UINT p, UINT chrom) const
	{
		CHECKRANGELOCUS(chrom, idx);
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(chrom);
		return static_cast<UINT>(*(m_genoPtr + idx + p * totNumLoci() + chromBegin(chrom)));
	}


	/** return the name of \c allele(idx)
	 *  HIDDEN
	 * <group>allele</group>
	 */
	string alleleChar(UINT idx) const
	{
		CHECKRANGEGENOSIZE(idx);

		return this->alleleName(allele(idx));
	}


	/** HIDDEN
	 * return the name of <tt>allele(idx, p)</tt>
	 * <group>allele</group>
	 */
	string alleleChar(UINT idx, UINT p) const
	{
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(p);

		return this->alleleName(allele(idx, p));
	}


	/** HIDDEN
	 * return the name of <tt>allele(idx, p, ch)</tt>
	 * <group>allele</group>
	 */
	string alleleChar(UINT idx, UINT p, UINT ch) const
	{
		CHECKRANGELOCUS(ch, idx);
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(ch);

		return this->alleleName(allele(idx, p, ch));
	}


	/** set allele \e allele to a locus, using its absolute index \e idx.
	 * <group>allele</group>
	 */
	void setAllele(Allele allele, UINT idx)
	{
		CHECKRANGEGENOSIZE(idx);
		*(m_genoPtr + idx) = allele;
	}


	/** set allele \e allele to locus \e idx on the <em>p</em>-th homologous
	 *  set of chromosomes.
	 * <group>allele</group>
	 */
	void setAllele(Allele allele, UINT idx, UINT p)
	{
		CHECKRANGEABSLOCUS(idx);
		CHECKRANGEPLOIDY(p);
		*(m_genoPtr + idx + p * totNumLoci()) = allele;
	}


	/** set allele \e allele to locus \e idx on chromosome \e chrom of the
	 *  <em>p</em>-th homologous set of chromosomes.
	 *  <group>allele</group>
	 */
	void setAllele(Allele allele, UINT idx, UINT p, UINT chrom)
	{
		CHECKRANGELOCUS(chrom, idx);
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(chrom);
		*(m_genoPtr + idx + p * totNumLoci() + chromBegin(chrom) ) = allele;
	}


	/** return an editable array (a \c carray of length <tt>totNumLoci()*ploidy()</tt>)
	 *  that represents all alleles of an individual.
	 * <group>genotype</group>
	 */
	PyObject * genotype()
	{
		// The following implementation has comparable performance as
		// the direct memory access implementation, but it has two 
		// problems:
		// 1. genotype has to be copied out, which requires additional
		//    allocation of memory.
		// 2. the return value is a tuple which lacks some functions of
		//    a list, mostly notably the count() function.
		/*
		UINT num = ploidy() * totNumLoci();
		vectora geno(num);	
		for(UINT i =0; i < num; i++)
			geno[i] = *(m_genoPtr + i);
		return geno;
		*/
		return Allele_Vec_As_NumArray(m_genoPtr, m_genoPtr + genoSize());
	}


    /** return an editable array (a \c carray of length <tt>totNumLoci()</tt>)
	 *  that represents all alleles on the <em>p</em>-th homologous set of
	 *  chromosomes.
	 * <group>genotype</group>
	 */
	PyObject * genotype(UINT p)
	{
		CHECKRANGEPLOIDY(p);
	    return Allele_Vec_As_NumArray(m_genoPtr + p * totNumLoci(),
			m_genoPtr + (p + 1) * totNumLoci() );
	}


	/** return an editable array (a \c carrary of legnth 
	 *  <tt>numLoci(</tt><em>chrom</em><tt>)</tt>) that represents all alleles
	 *  on chromosome \e chrom of the <em>p</em>-th homologous set of
	 *  chromosomes.
	 *  <group>genotype</group>
	 */
	PyObject * genotype(UINT p, UINT chrom)
	{
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(chrom);
		return Allele_Vec_As_NumArray(m_genoPtr + p * totNumLoci() + chromBegin(chrom),
			m_genoPtr + p * totNumLoci() + chromEnd(chrom));
	}


	/** Fill the genotype of an individual using a list of alleles \e geno.
	 *  \c geno will be reused if its length is less than 
	 *  <tt>totNumLoci()*ploidy()</tt>.
	 *  <group>genotype</group>
	 */
	void setGenotype(vectora geno) 
	{
		UINT sz = geno.size();
		for(UINT i =0; i < totNumLoci()*ploidy(); i++)
			*(m_genoPtr + i) = geno[i % sz];
	}


	/** Fill the genotype of the <em>p</em>-th homologous set of chromosomes
	 *  using a list of alleles \e geno. \c geno will be reused if its length
	 *  is less than <tt>totNumLoci()</tt>.
	 *  <group>genotype</group>
	 */
	void setGenotype(vectora geno, UINT p) 
	{
		CHECKRANGEPLOIDY(p);
		GenoIterator ptr = m_genoPtr + p*totNumLoci();
		UINT sz = geno.size();
		for(UINT i =0; i < totNumLoci(); i++)
			*(ptr + i) = geno[i % sz];
	}


	/** Fill the genotype of chromosome \e chrom on the <em>p</em>-th
	 *  homologous set of chromosomes using a list of alleles \e geno.
	 *  \c geno will be reused if its length is less than 
	 *  <tt>mumLoci(</tt><em>chrom</em><tt>)</tt>.
	 * <group>genotype</group>
	 */
	void setGenotype(vectora geno, UINT p, UINT chrom) 
	{
		CHECKRANGEPLOIDY(p);
	    CHECKRANGECHROM(chrom);
		GenoIterator ptr = m_genoPtr + p*totNumLoci()+chromBegin(chrom);
		UINT sz = geno.size();
		for(UINT i =0; i < numLoci(ch); i++)
			*(ptr+i) = geno[i % sz];
	}


	/** return the sex of an individual, \c 1 for male and \c 2 for female.
	 * <group>sex</group>
	 */
	Sex sex() const
	{
		if (ISSETFLAG(m_flags, m_flagFemale) )
			return Female;
		else
			return Male;
	}


	/** return the sex of an individual, \c M for male or \c F for \c female.
	 * <group>sex</group>
	 */
	char sexChar() const
	{
		return sex() == Female ? 'F' : 'M';
	}


	/** set individual sex to \c Male or \c Female.
	 * <group>sex</group>
	 */
	void setSex(Sex sex)
	{
		CHECKRANGESEX(sex);

		if (sex == Male)
			RESETFLAG(m_flags, m_flagFemale);
		else
			SETFLAG(m_flags, m_flagFemale);
	}


	/** Return \c True if this individual is affected.
	 * <group>affection</group>
	 */
	bool affected() const
	{
		return ISSETFLAG(m_flags, m_flagAffected);
	}


	/** HIDDEN
	 */
	bool unaffected() const
	{
		return !affected();
	}


	/** Return \c A if this individual is affected, or \c U otherwise.
	 * <group>affection</group>
	 */
	char affectedChar() const
	{
		return affected() ? 'A' : 'U';
	}


	/** set affection status to \e affected (\c True or \c False).
	 * <group>affection</group>
	 */
	void setAffected(bool affected)
	{
		if (affected)
			SETFLAG(m_flags, m_flagAffected);
		else
			RESETFLAG(m_flags, m_flagAffected);
	}


	/// CPPONLY
	bool iteratable() const
	{
		return ISSETFLAG(m_flags, m_flagIteratable);
	}


	/// CPPONLY
	void setIteratable(bool iteratable)
	{
		if (iteratable)
			SETFLAG(m_flags, m_flagIteratable);
		else
			RESETFLAG(m_flags, m_flagIteratable);
	}


	/// CPPONLY
	bool visible() const
	{
		return ISSETFLAG(m_flags, m_flagVisible);
	}


	/// CPPONLY
	void setVisible(bool visible)
	{
		if (visible)
			SETFLAG(m_flags, m_flagVisible);
		else
			RESETFLAG(m_flags, m_flagVisible);
	}


	/// return the ID of the subpopulation to which this individual blongs
	/** HIDDEN
	 \note \c subPopID is not set by default. It only corresponds to the subpopulation
	   	in which this individual resides after \c pop::setIndSubPopID is called.
	 */
	SubPopID subPopID() const
	{
		return m_subPopID;
	}


	/** HIDDEN
	 *set new subpopulation ID, \c pop.rearrangeByIndID will move this individual to that population
	 */
	void setSubPopID(SubPopID id)
	{
		m_subPopID = id;
	}


	/** Return the value of an information field \e idx (an index).
	 * <group>info</group>
	 */
	InfoType info(UINT idx) const
	{
		CHECKRANGEINFO(idx);

		return m_infoPtr[idx];
	}


	/** Return the value of an information field \e idx (an index) as an integer number.
	 * <group>info</group>
	 */
	int intInfo(UINT idx) const
	{
		CHECKRANGEINFO(idx);
		return static_cast<int>(m_infoPtr[idx]);
	}


	/** Return the value of an information field \e name.
	 *  <group>info</group>
	 */
	InfoType info(const string & name) const
	{
		int idx = infoIdx(name);

		DBG_ASSERT(idx >= 0, IndexError,
			"Info name " + name + " is not a valid info field name");
		return m_infoPtr[idx];
	}


	/** Return the value of an information field \e name as an integer number.
	 * <group>info</group>
	 */
	int intInfo(const string & name) const
	{
		int idx = infoIdx(name);

		DBG_ASSERT(idx >= 0, IndexError,
			"Info name " + name + " is not a valid info field name");
		return static_cast<int>(m_infoPtr[idx]);
	}


	/** set the value of an information field \e idx (an index) to \e value.
	 *  <group>info</group>
	 */
	void setInfo(InfoType value, UINT idx)
	{
		CHECKRANGEINFO(idx);
		m_infoPtr[idx] = value;
	}


	/** set the value of an information field \e name to \e value.
	 *  <group>info</group>
	 */
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
		return m_genoPtr + genoSize();
	}


	/// CPPONLY start of allele of the pth set of chromosome
	GenoIterator genoBegin(UINT p) const
	{
		CHECKRANGEPLOIDY(p);
		return m_genoPtr + p * totNumLoci();

	}


	/// CPPONLY end of allele of the pth set of chromosome
	GenoIterator genoEnd(UINT p) const
	{
		CHECKRANGEPLOIDY(p);
		return m_genoPtr + (p + 1) * totNumLoci();
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

	   The guideline is that if we swap individuals across
	   subpopulation, we should swap content. Otherwise,
	   swap pointers. (There is no order right now within
	   subpopulation so the later case is rare, at best.)
	 */
	void swap(individual & ind, bool swapContent = true);

	//@}
	/// @name misc (only relevant to developers)
	//@{

	/// CPPONLY
	void display(ostream & out, int width = 1, const vectori & chrom = vectori(), const vectori & loci = vectori() );

	//@}

private:
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const UINT version) const
	{
		//  boost::serialization::base_object<GenoStruTrait>(*this));
		bool b;

		b = ISSETFLAG(m_flags, m_flagFemale);
		ar & b;

		b = ISSETFLAG(m_flags, m_flagAffected);
		ar & b;
	}


	template<class Archive>
	void load(Archive & ar, const UINT version)
	{
		bool b;

		m_flags = 0;
		ar & b;
		if (b) SETFLAG(m_flags, m_flagFemale);
		ar & b;
		if (b) SETFLAG(m_flags, m_flagAffected);
		SETFLAG(m_flags, m_flagVisible);
		SETFLAG(m_flags, m_flagIteratable);

		if (version < 1) {
			std::pair<int, int> tag;
			ar & tag;
			ar & m_subPopID;
		}
	}


	BOOST_SERIALIZATION_SPLIT_MEMBER();

protected:
	// internal flag. Can be used to perform many things.
	// bitset<3> was previously used but that will take 4 bytes.
	unsigned char m_flags;

	/// temporary information
	SubPopID m_subPopID;

	/// pointer to genotype.
	GenoIterator m_genoPtr;

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
	              vector<individual>::iterator const end,
	              bool allInds, bool allVisibles) :
		m_index(begin),
		m_end(end),
		m_allInds(allInds),
		m_allVisibles(allVisibles)
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

	// ending idx
	vector<individual>::iterator m_end;

	//
	bool m_allInds;
	//
	bool m_allVisibles;
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
		while (m_it < m_end && !m_it->visible())
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
		if (m_allInds)
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


	IndividualIterator operator+=(difference_type diff)
	{
		if (m_allInds) {
			m_it += diff;
			return *this;
		}
		DBG_ASSERT(m_it < m_end, ValueError,
			"Can not advance invalid iterator");
		difference_type i = 0;
		while (i < diff && m_it < m_end)
			if ((++m_it)->visible())
				++i;
		DBG_FAILIF(i != diff, ValueError, "Can not add to IndIterator");
		return *this;
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

/**
   	this class implements a C++ iterator class that iterate through
   	infomation fields in a (sub)population using
   	1. an IndIterator that	will skip invisible individuals, or
   	2. a gapped iterator that will run faster.
   	Note that 1, 2 should yield identical result, and 2 should be used
   	when there is no virtual subpopulation.q
 */
template <typename T>
class InformationIterator
{
public:
	typedef std::forward_iterator_tag iterator_category;
	typedef typename T::value_type value_type;
	typedef long int difference_type;
	typedef typename T::reference reference;
	typedef typename T::pointer pointer;


	InformationIterator()
	{
	}


	InformationIterator(UINT info, InfoIterator ptr, UINT size)
		: m_info(info), m_useGappedIterator(true),
		m_it(), m_ptr(& * ptr), m_step(size)
	{
	}


	InformationIterator(UINT info, IndividualIterator<T> it)
		: m_info(info), m_useGappedIterator(false),
		m_it(it), m_ptr(), m_step()
	{
	}


	// this is the most important part!
	InfoType & operator *() const
	{
		if (m_useGappedIterator)
			return *m_ptr;
		else
			return *(m_it->infoPtr() + m_info);
	}


	// return, then advance.
	InformationIterator operator++(int)
	{
		// save current state
		InformationIterator tmp(*this);

		if (m_useGappedIterator)
			m_ptr += m_step;
		else
			++m_it;
		return tmp;
	}


	InformationIterator operator++()
	{
		if (m_useGappedIterator)
			m_ptr += m_step;
		else
			++m_it;
		return *this;
	}


	bool operator!=(const InformationIterator & rhs)
	{
		if (m_useGappedIterator) {
			DBG_FAILIF(m_step != rhs.m_step, ValueError,
				"Iterator comparison failed");
			return m_ptr != rhs.m_ptr || m_info != rhs.m_info;
		} else
			return m_it != rhs.m_it || m_info != rhs.m_info;
	}


private:
	// idx of the information field
	UINT m_info;
	///
	bool m_useGappedIterator;
	// individual iterator
	IndividualIterator<T> m_it;
	//
	InfoType * m_ptr;
	//
	UINT m_step;
};


typedef InformationIterator<RawIndIterator> IndInfoIterator;
typedef InformationIterator<ConstRawIndIterator> ConstIndInfoIterator;

/**
   	this class implements a C++ iterator class that iterate through
   	infomation fields in a (sub)population using
   	1. an IndIterator that	will skip invisible individuals, or
   	2. a gapped iterator that will run faster.
   	Note that 1, 2 should yield identical result, and 2 should be used
   	when there is no virtual subpopulation.q
 */
template <typename T>
class CombinedAlleleIterator
{
public:
	typedef std::forward_iterator_tag iterator_category;
	typedef Allele value_type;
	typedef long int difference_type;
	typedef AlleleRef reference;
	typedef GenoIterator pointer;


	CombinedAlleleIterator()
	{
	}


	CombinedAlleleIterator(GenoIterator ptr, UINT size)
		: m_useGappedIterator(true), m_it(), m_ptr(ptr), m_size(size)
	{
	}


	CombinedAlleleIterator(UINT idx, IndividualIterator<T> it,
	                       UINT ploidy, UINT size)
		: m_index(idx), m_useGappedIterator(false),
		m_it(it), m_ptr(), m_p(0), m_ploidy(ploidy), m_size(size)
	{
	}


	// this is the most important part!
	AlleleRef operator *() const
	{
		if (m_useGappedIterator)
			return *m_ptr;
		else
			return *(m_it->genoBegin() + m_index + m_p * m_size);
	}


	GenoIterator ptr()
	{
		if (m_useGappedIterator)
			return m_ptr;
		else
			return m_it->genoBegin() + m_index + m_p * m_size;
	}


	// return, then advance.
	CombinedAlleleIterator operator++(int)
	{
		// save current state
		CombinedAlleleIterator tmp(*this);

		if (m_useGappedIterator)
			m_ptr += m_size;
		else {
			++m_p;
			if (m_p == m_ploidy) {
				m_p = 0;
				++m_it;
			}
		}
		return tmp;
	}


	CombinedAlleleIterator operator++()
	{
		if (m_useGappedIterator)
			m_ptr += m_size;
		else {
			++m_p;
			if (m_p == m_ploidy) {
				m_p = 0;
				++m_it;
			}
		}
		return *this;
	}


	CombinedAlleleIterator & operator+=(difference_type diff)
	{
		if (m_useGappedIterator)
			m_ptr += diff * m_size;
		else {
			m_p += diff;
			m_it += m_p / m_ploidy;
			m_p %= m_ploidy;
		}
		return *this;
	}


	CombinedAlleleIterator operator+(difference_type diff)
	{
		CombinedAlleleIterator tmp(*this);

		if (m_useGappedIterator)
			tmp.m_ptr += diff * m_size;
		else {
			tmp.m_p += diff;
			tmp.m_it += tmp.m_p / m_ploidy;
			tmp.m_p %= m_ploidy;
		}
		return tmp;
	}


	bool operator!=(const CombinedAlleleIterator & rhs)
	{
		if (m_useGappedIterator)
			return m_ptr != rhs.m_ptr;
		else {
			DBG_FAILIF(m_ploidy != rhs.m_ploidy || m_size != rhs.m_size, ValueError,
				"Iterator comparison fails");
			return m_it != rhs.m_it || m_index != rhs.m_index || m_p != rhs.m_p;
		}
	}


private:
	// index of the information field
	UINT m_index;
	///
	bool m_useGappedIterator;
	// individual iterator
	IndividualIterator<T> m_it;
	//
	GenoIterator m_ptr;
	// current ploidy, used in individualiterator
	UINT m_p;
	// overall ploidy
	UINT m_ploidy;
	// genosize
	UINT m_size;
};


typedef CombinedAlleleIterator<RawIndIterator> IndAlleleIterator;
typedef CombinedAlleleIterator<ConstRawIndIterator> ConstIndAlleleIterator;


}

#ifndef SWIG
// set version for individual class
// version 0: base
// version 1: add sexChrom indicator
BOOST_CLASS_VERSION(simuPOP::individual, 1)
#endif

#endif
