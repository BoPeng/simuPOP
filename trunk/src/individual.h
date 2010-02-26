/**
 *  $File: individual.h $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2009 Bo Peng (bpeng@mdanderson.org)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _INDIVIDUAL_H
#define _INDIVIDUAL_H

/**
   \file
   \brief class Individual, IndividualWithAge etc.
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
using std::endl;
using std::hex;
using std::dec;

#include <numeric>
using std::pair;

namespace simuPOP {
/**
 *  A \c Population consists of individuals with the same genotypic structure.
 *  An \c Individual object cannot be created independently, but refences to
 *  inidividuals can be retrieved using member functions of a \c Population
 *  object. In addition to structural information shared by all individuals in
 *  a population (provided by class \c GenoStruTrait), the \c Individual class
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
 *  <tt>Individual::totNumLoci()</tt> loci. It is worth noting that access to
 *  invalid chromosomes, such as the Y chromosomes of female individuals, is
 *  not restricted.
 */
class Individual : public GenoStruTrait
{

protected:
	/// 0: male, 1: female regardless of outside coding
	static const size_t m_flagFemale = 1;

	/// if this individual is affect
	static const size_t m_flagAffected = 2;

	/// if this individual is visible. This is used
	/// to implement virtual subpopulations
	static const size_t m_flagVisible = 4;

	/// a temporary mark to mark individuals for deletion
	/// or extraction.
	static const size_t m_flagMarked = 8;

public:
	///  @name constructor, destructor etc
	//@{
	///
	/**
	 * An \c Individual object cannot be created directly. It has to be accessed
	 * from a \c Population object using functions such as
	 * <tt>Population::Individual(idx)</tt>.
	 */
	Individual() : m_flags(m_flagVisible)
	{
	}


	/// CPPONLY
	Individual(const Individual & ind) :
		GenoStruTrait(ind), m_flags(ind.m_flags),
		m_genoPtr(ind.m_genoPtr),
		m_infoPtr(ind.m_infoPtr)
	{
	}


	/// destructor. Do nothing.
	~Individual()
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
	Individual & operator=(const Individual & rhs);

	/// CPPONLY deep copy of an individual class
	Individual & copyFrom(const Individual & rhs);

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

	/** return the current allele at a locus, using its absolute index \e idx.
	 *  If a ploidy \e ploidy and/or a chromosome indexes is given, \e idx is
	 *  relative to the beginning of specified homologous copy of chromosomes
	 *  (if \e chrom=-1) or the beginning of the specified homologous copy of
	 *  specified chromosome (if \e chrom >= 0).
	 *
	 * <group>1-allele</group>
	 */
	UINT allele(UINT idx, int ploidy = -1, int chrom = -1) const;


	/** return the name of <tt>allele(idx, ploidy, chrom)</tt>. If idx is
	 *  invalid (e.g. second homologus copy of chromosome Y), '_' is returned.
	 *  <group>1-allele</group>
	 */
	string alleleChar(UINT idx, int ploidy = -1, int chrom = -1) const;

	/** set allele \e allele to a locus, using its absolute index \e idx.
	 *  If a ploidy \e ploidy and/or a chromosome indexes are given, \e idx is
	 *  relative to the beginning of specified homologous copy of chromosomes
	 *  (if \e chrom=-1) or the beginning of the specified homologous copy of
	 *  specified chromosome (if \e chrom >= 0).
	 *  <group>1-allele</group>
	 */
	void setAllele(Allele allele, UINT idx, int ploidy = -1, int chrom = -1);

	/** return an editable array (a \c carray object) that represents all
	 *  alleles of an individual. If \e ploidy or \e chroms is given, only
	 *  alleles on the specified chromosomes and homologous copy of chromosomes
	 *  will be returned. If multiple chromosomes are specified, there should
	 *  not be gaps between chromosomes.
	 *  <group>2-genotype</group>
	 */
	PyObject * genotype(const uintList & ploidy = uintList(), const uintList & chroms = uintList());


	/** CPPONLY
	 *  Return a Python object with alleles at specified loci. This function
	 *  is usually used to collect alleles to send to a user-provided function.
	 */
	PyObject * genoAtLoci(const vectoru & loci);

	/** Fill the genotype of an individual using a list of alleles \e geno.
	 *  If parameters \e ploidy and/or \e chroms are specified, alleles will
	 *  be copied to only all or specified chromosomes on selected homologous
	 *  copies of chromosomes. \c geno will be reused if its length is less than
	 *  number of alleles to be filled.
	 *  <group>2-genotype</group>
	 */
	void setGenotype(const vectoru & geno, const uintList & ploidy = uintList(), const uintList & chroms = uintList());

	/** return the sex of an individual, \c 1 for male and \c 2 for female.
	 * <group>3-sex</group>
	 */
	Sex sex() const
	{
		if (ISSETFLAG(m_flags, m_flagFemale))
			return FEMALE;
		else
			return MALE;
	}


	/** set individual sex to \c MALE or \c FEMALE.
	 * <group>3-sex</group>
	 */
	void setSex(Sex sex)
	{
		CHECKRANGESEX(sex);

		if (sex == MALE)
			RESETFLAG(m_flags, m_flagFemale);
		else
			SETFLAG(m_flags, m_flagFemale);
	}


	/** Return \c True if this individual is affected.
	 * <group>4-affection</group>
	 */
	bool affected() const
	{
		return ISSETFLAG(m_flags, m_flagAffected);
	}


	/** set affection status to \e affected (\c True or \c False).
	 * <group>4-affection</group>
	 */
	void setAffected(bool affected)
	{
		if (affected)
			SETFLAG(m_flags, m_flagAffected);
		else
			RESETFLAG(m_flags, m_flagAffected);
	}


	/// CPPONLY
	bool visible() const
	{
		return ISSETFLAG(m_flags, m_flagVisible);
	}


	/// CPPONLY
	void setVisible(bool visible) const
	{
		if (visible)
			SETFLAG(m_flags, m_flagVisible);
		else
			RESETFLAG(m_flags, m_flagVisible);
	}


	/** CPPONLY
	 *  check if an individual is marked. This is a temporary flag that is
	 *  usually used to mark individuals for removal or extraction. It might be
	 *  be reset with any Population operation.
	 *  <group>5-mark</group>
	 */
	bool marked() const
	{
		return ISSETFLAG(m_flags, m_flagMarked);
	}


	/** CPPONLY
	 *  mark (default) or unmark (if \e mark=false) an individual.
	 *  <group>5-mark</group>
	 */
	void setMarked(bool mark = true) const
	{
		if (mark)
			SETFLAG(m_flags, m_flagMarked);
		else
			RESETFLAG(m_flags, m_flagMarked);
	}


	/** Return the value of an information field \e filed (by index or name).
	 *  <tt>ind.info(name)</tt> is equivalent to <tt>ind.name</tt> although the
	 *  function form allows the use of indexes of information fieldes.
	 * <group>5-info</group>
	 */
	double info(const uintString & field) const
	{
		UINT idx = field.empty() ? field.value() : infoIdx(field.name());

		CHECKRANGEINFO(idx);
		return m_infoPtr[idx];
	}


	/** Return the value of an information field \e field (by index or name)
	 *  as an integer number.
	 *  CPPONLY
	 */
	int intInfo(const uintString & field) const
	{
		UINT idx = field.empty() ? field.value() : infoIdx(field.name());

		CHECKRANGEINFO(idx);
		return static_cast<int>(m_infoPtr[idx]);
	}


	/** read info as attribute
	 */
	double __getattr__(const string & field) const
	{
		return m_infoPtr[infoIdx(field)];
	}


	/** write info as attribute
	 */
	void __setattr__(const string & field, double value) const
	{
		m_infoPtr[infoIdx(field)] = value;
	}


	/** set the value of an information field \e field (by index or name) to
	 *  \e value. <tt>ind.setInfo(value, field)</tt> is equivalent to
	 *  <tt>ind.field = value</tt> although the function form allows the use
	 *  of indexes of information fieldes.
	 *  <group>5-info</group>
	 */
	void setInfo(double value, const uintString & field)
	{
		UINT idx = field.empty() ? field.value() : infoIdx(field.name());

		CHECKRANGEINFO(idx);
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
	bool operator==(const Individual & rhs) const;

	/// compare if two individuals are not the same used in case of serialization etc.
	bool operator!=(const Individual & rhs) const
	{
		return !(*this == rhs);
	}


	// allow compaison of individuals in python
	// only equal or unequal, no greater or less than
	/// a python function used to compare the individual objects
	int __cmp__(const Individual & rhs) const;


	/// CPPONLY swap Individuals
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
	void swap(Individual & ind, bool swapContent = true);

	//@}
	/// @name misc (only relevant to developers)
	//@{

	/// CPPONLY
	void display(ostream & out, int width = 1, const vectoru & loci = vectoru());

	//@}

private:
	bool validIndex(UINT idx) const;

	bool validIndex(UINT idx, UINT p) const;

	bool validIndex(UINT idx, UINT p, UINT ch) const;

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
	}


	BOOST_SERIALIZATION_SPLIT_MEMBER();

protected:
	// internal flag. Can be used to perform many things.
	// bitset<3> was previously used but that will take 4 bytes.
	mutable unsigned char m_flags;

	/// pointer to genotype.
	GenoIterator m_genoPtr;

	/// pointer to info
	InfoIterator m_infoPtr;
};


/** CPPONLY
 *  A class used to compare two individuals by an information field.
 */
class indCompare
{
public:
	// accept the index to an information field
	indCompare(UINT idx) : m_field(idx) {}

	bool operator()(const Individual & lhs, const Individual & rhs)
	{
		return lhs.info(m_field) < rhs.info(m_field);
	}


private:
	UINT m_field;
};


/**
    this class implements a C++ iterator class that iterate through
    individuals in a (sub)population. If allInds are true, the
    visiblility of individuals will not be checked. Note that
    IndividualIterator *will* iterate through only visible individuals, and
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


	IndividualIterator(T it, T end, bool allInds)
		: m_it(it), m_end(end), m_allInds(allInds)
	{
		// m_it need to point to the first valid
		// Individual. otherwise *it will fail.
		if (!allInds)
			while (m_it < m_end && !m_it->visible())
				++m_it;
	}


	bool valid() const
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
		return &*m_it;
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
		IndividualIterator tmp(*this);
		// can not check. Possible problem
		for (difference_type i = 0; i < diff; ++i)
			while (!(--tmp.m_it)->visible()) ;
		return tmp;
	}


	difference_type operator-(IndividualIterator rhs)
	{
		if (m_allInds)
			return m_it - rhs.m_it;
		difference_type i = 0;
		for (T it = rhs.m_it; it != m_it; ++it)
			if (it->visible())
				++i;
		return i;
	}


	IndividualIterator operator--(int)
	{
		if (m_allInds)
			return IndividualIterator(m_it--, m_end, m_allInds);
		IndividualIterator tmp(*this);
		while (!(--tmp.m_it)->visible()) ;
		return tmp;
	}


	IndividualIterator operator--()
	{
		if (m_allInds)
			--m_it;
		else
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

	//
	bool m_allInds;
};

//
typedef vector<Individual>::iterator RawIndIterator;
typedef vector<Individual>::const_iterator ConstRawIndIterator;

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
		m_it(), m_ptr(&*ptr), m_step(size)
	{
	}


	InformationIterator(UINT info, IndividualIterator<T> it)
		: m_info(info), m_useGappedIterator(false),
		m_it(it), m_ptr(), m_step()
	{
	}


	bool valid()
	{
		return m_useGappedIterator || m_it.valid();
	}


	// this is the most important part!
	double & operator *() const
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
	// Individual iterator
	IndividualIterator<T> m_it;
	//
	double * m_ptr;
	//
	UINT m_step;
};


typedef InformationIterator<RawIndIterator> IndInfoIterator;
typedef InformationIterator<ConstRawIndIterator> ConstIndInfoIterator;

/** This class implements a C++ iterator class that iterate through
    all alleles in a (virtual) (sub)population using
    1. an IndIterator that will skip invisible individuals and invalid
        alleles, or
    2. a gapped iterator that will run faster, in the case that
      a): no virtual subpopulation
      b): not sex chromosomes
      c): not haplodiploid

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


	CombinedAlleleIterator(GenoIterator ptr, GenoIterator ptrEnd, UINT size)
		: m_useGappedIterator(true), m_valid(true),
		m_ptr(ptr), m_ptrEnd(ptrEnd), m_size(size),
		// ignored
		m_it(), m_index(0), m_ploidy(0), m_chromType(0),
		m_haplodiploid(false), m_p(0)
	{
		m_valid = m_ptr != m_ptrEnd;
	}


	CombinedAlleleIterator(UINT idx, IndividualIterator<T> it)
		: m_useGappedIterator(false), m_valid(true),
		m_ptr(), m_ptrEnd(), m_size(0), // belong to a previous one
		m_it(it), m_index(idx), m_ploidy(0), m_chromType(0),
		m_haplodiploid(false), m_p(0)
	{
		if (!it.valid()) {
			m_valid = false;
			return;
		}

		m_size = it->totNumLoci();
		m_ploidy = it->ploidy();
		m_haplodiploid = it->isHaplodiploid();
		m_chromType = it->chromType(it->chromLocusPair(idx).first);
		// we do not know anything about customized chromosome
		// so we just assume it is autosome.
		if (m_chromType == CUSTOMIZED)
			m_chromType = AUTOSOME;
		//
		if (m_chromType == CHROMOSOME_Y) {
			if (m_it->sex() == FEMALE) {
				while (m_it.valid())
					if ((++m_it)->sex() == MALE)
						break;
				m_valid = m_it.valid();
			}
			m_p = 1;
		}
	}


	bool valid()
	{
		return m_valid;
	}


	// this is the most important part!
	AlleleRef operator *() const
	{
		if (m_useGappedIterator)
			return *m_ptr;
		else {
			DBG_ASSERT(m_it.valid(), SystemError, "Cannot refer to an invalid individual iterator");
			return *(m_it->genoBegin() + m_index + m_p * m_size);
		}
	}


	GenoIterator ptr()
	{
		if (m_useGappedIterator)
			return m_ptr;
		else
			return m_it->genoBegin() + m_index + m_p * m_size;
	}


	void advance(IndividualIterator<T> & it, UINT & p, bool & valid)
	{
		DBG_ASSERT(valid, RuntimeError, "Can not advance invalid allele iterator");
		if (m_chromType == AUTOSOME) {
			++p;
			if (p == m_ploidy) {
				p = 0;
				++it;
				valid = it.valid();
			}
		} else if (m_chromType == CHROMOSOME_X) {
			if (it->sex() == FEMALE) {
				// X0 -> X1
				if (p == 0)
					++p;
				// X1 -> X0 of next ind (male or female)
				else {
					p = 0;
					++it;
					valid = it.valid();
				}
			} else {
				// male, no X1.
				DBG_ASSERT(p == 0, SystemError,
					"Male Individual only has the first homologous copy of chromosome X");
				// next Individual, ploidy 0, sex does not matter.
				++it;
				valid = it.valid();
			}
		} else if (m_chromType == CHROMOSOME_Y) {
			DBG_ASSERT(it->sex() == MALE, SystemError,
				"There is no chromosome Y for FEMALE Individuals");
			while (it.valid())
				if ((++it)->sex() == MALE)
					break;
			p = 1;
			valid = it.valid();
		}
	}


	// return, then advance.
	CombinedAlleleIterator operator++(int)
	{
		// save current state
		CombinedAlleleIterator tmp(*this);

		if (!valid())
			return tmp;

		if (m_useGappedIterator) {
			m_ptr += m_size;
			m_valid = m_ptr != m_ptrEnd;
		} else {
			DBG_ASSERT(m_it.valid(), SystemError, "Cannot refer to an invalid individual iterator");
			advance(m_it, m_p, m_valid);
		}
		return tmp;
	}


	CombinedAlleleIterator operator++()
	{
		if (!valid())
			return *this;
		if (m_useGappedIterator) {
			m_ptr += m_size;
			m_valid = m_ptr != m_ptrEnd;
		} else {
			DBG_ASSERT(m_it.valid(), SystemError, "Cannot refer to an invalid individual iterator");
			advance(m_it, m_p, m_valid);
		}
		return *this;
	}


	CombinedAlleleIterator & operator+=(difference_type diff)
	{
		if (!valid())
			return *this;
		if (m_useGappedIterator) {
			m_ptr += diff * m_size;
			if (m_ptr > m_ptrEnd)
				m_ptr = m_ptrEnd;
			m_valid = m_ptr != m_ptrEnd;
		} else {
			DBG_ASSERT(m_it.valid(), SystemError, "Cannot refer to an invalid individual iterator");
			for (int i = 0; i < diff && m_valid; ++i)
				advance(m_it, m_p, m_valid);
		}
		return *this;
	}


	CombinedAlleleIterator operator+(difference_type diff)
	{
		CombinedAlleleIterator tmp(*this);

		if (!valid())
			return tmp;

		if (m_useGappedIterator) {
			tmp.m_ptr += diff * m_size;
			if (tmp.m_ptr > tmp.m_ptrEnd)
				tmp.m_ptr = tmp.m_ptrEnd;
			tmp.m_valid = tmp.m_ptr != tmp.m_ptrEnd;
		} else {
			DBG_ASSERT(m_it.valid(), SystemError, "Cannot refer to an invalid individual iterator");
			for (int i = 0; i < diff && tmp.m_valid; ++i)
				advance(tmp.m_it, tmp.m_p, tmp.m_valid);
		}
		return tmp;
	}


	bool operator!=(const CombinedAlleleIterator & rhs)
	{
		if (m_useGappedIterator)
			return m_ptr != rhs.m_ptr;
		else {
			//DBG_FAILIF(m_it.valid() && rhs.m_it.valid() &&
			//	(m_ploidy != rhs.m_ploidy || m_size != rhs.m_size
			//	|| m_chromType != rhs.m_chromType), ValueError,
			//	"Iterator comparison fails");
			return !(m_index == rhs.m_index && m_it == rhs.m_it &&
			         (m_p == rhs.m_p || !m_it.valid() || !rhs.m_it.valid()));
		}
	}


private:
	///
	bool m_useGappedIterator;
	///
	bool m_valid;
	//
	GenoIterator m_ptr;
	//
	GenoIterator m_ptrEnd;
	// genosize
	UINT m_size;
	//
	// The second iteration method
	// Individual iterator
	IndividualIterator<T> m_it;
	// index of the locus
	UINT m_index;
	// overall ploidy
	UINT m_ploidy;
	// chromosome type
	int m_chromType;
	//
	bool m_haplodiploid;
	// current ploidy, used in individualiterator
	UINT m_p;
};


typedef CombinedAlleleIterator<RawIndIterator> IndAlleleIterator;
typedef CombinedAlleleIterator<ConstRawIndIterator> ConstIndAlleleIterator;


}

#ifndef SWIG
// set version for individual class
// version 0: base (reset for version 1.0)
BOOST_CLASS_VERSION(simuPOP::Individual, 0)
#endif

#endif
