/**
 *  $File: individual.h $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
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


/** this class implements a Python itertor class that can be used to iterate
 *  through individuals in a (sub)population. If allInds are true,
 *  visiblility of individuals will not be checked. Otherwise, a functor
 *  will be used to check if indiviudals belong to a specified virtual
 *  subpopulation.
 *
 *  An instance of this class is returned by
 *  population::Individuals() and Population::Individuals(subPop)
 */
class pyMutantIterator
{
public:
	pyMutantIterator(GenoIterator base, size_t begin,
		size_t end, size_t step) :
#ifdef MUTANTALLELE
		m_ptr((base + begin).get_val_iterator()),
		m_end((base + end).get_val_iterator()),
#else
		m_base(base),
		m_ptr(begin),
		m_end(end),
#endif
		m_step(step)
	{
	}


	~pyMutantIterator()
	{
	}


	pyMutantIterator __iter__()
	{
		return *this;
	}


	// python 2.x uses next()
	pairu next()
	{
#ifdef MUTANTALLELE
		do {
			if (m_ptr == m_end)
				throw StopIteration("");
			// the mutants can be mutated back to zero, in this case the mutant is
			// not removed (because removal is an expensive operation)
			else if (m_ptr->second != 0) {
				vectorm::const_val_iterator tmp = m_ptr++;
				return pairu(tmp->first % m_step, tmp->second);
			} else
				++m_ptr;
		} while (true);
#else
		do {
			if (m_ptr == m_end)
				throw StopIteration("");
			else if (*(m_base + m_ptr) != 0) {
				size_t tmp = m_ptr++;
				return pairu(tmp % m_step, *(m_base + tmp));
			} else
				++m_ptr;
		} while (true);
#endif
	}


	// python 3.x uses __next__ instead of next.
	pairu __next__()
	{
		return next();
	}


private:
#ifdef MUTANTALLELE
	vectorm::const_val_iterator m_ptr;
	vectorm::const_val_iterator m_end;
#else
	GenoIterator m_base;
	size_t m_ptr;
	size_t m_end;
#endif
	size_t m_step;
};


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
	static const unsigned char m_flagFemale = 1;

	/// if this individual is affect
	static const unsigned char m_flagAffected = 2;

	/// if this individual is visible. This is used
	/// to implement virtual subpopulations
	static const unsigned char m_flagVisible = 4;

	/// a temporary mark to mark individuals for deletion
	/// or extraction.
	static const unsigned char m_flagMarked = 8;

	/// a temporary mark to mark individuals as the first
	/// offspring in the family
	static const unsigned char m_flagFirstOffspring = 16;

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
#ifdef LINEAGE
		m_lineagePtr(ind.m_lineagePtr),
#endif
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


#ifdef LINEAGE
	/// CPPONLY set lineage pointer
	void setLineagePtr(LineageIterator pos)
	{
		m_lineagePtr = pos;
	}


#endif

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


#ifdef LINEAGE
	/// CPPONLY pointer to lineage
	LineageIterator lineagePtr() const
	{
		return m_lineagePtr;
	}


#endif

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
	ULONG allele(size_t idx, ssize_t ploidy = -1, ssize_t chrom = -1) const;

	/** return the name of <tt>allele(idx, ploidy, chrom)</tt>. If idx is
	 *  invalid (e.g. second homologus copy of chromosome Y), '_' is returned.
	 *  <group>1-allele</group>
	 */
	string alleleChar(size_t idx, ssize_t ploidy = -1, ssize_t chrom = -1) const;

	/** set allele \e allele to a locus, using its absolute index \e idx.
	 *  If a ploidy \e ploidy and/or a chromosome indexes are given, \e idx is
	 *  relative to the beginning of specified homologous copy of chromosomes
	 *  (if \e chrom=-1) or the beginning of the specified homologous copy of
	 *  specified chromosome (if \e chrom >= 0).
	 *  <group>1-allele</group>
	 */
	void setAllele(ULONG allele, size_t idx, int ploidy = -1, int chrom = -1);

	/** return the lineage of the allele at a locus, using its absolute index
	 *  \e idx. If a ploidy \e ploidy and/or a chromosome indexes is given,
	 *  \e idx is relative to the beginning of specified homologous copy of
	 *  chromosomes (if \e chrom=-1) or the beginning of the specified
	 *  homologous copy of specified chromosome (if \e chrom >= 0). <bf>
	 *  This function returns 0 for modules without lineage information.</bf>
	 *
	 * <group>1-allele</group>
	 */
	long alleleLineage(size_t idx, ssize_t ploidy = -1, ssize_t chrom = -1) const;


	/** set lineage \e lineage to an allele, using its absolute index \e idx.
	 *  If a ploidy \e ploidy and/or a chromosome indexes are given, \e idx is
	 *  relative to the beginning of specified homologous copy of chromosomes
	 *  (if \e chrom=-1) or the beginning of the specified homologous copy of
	 *  specified chromosome (if \e chrom >= 0). This function does nothing
	 *  for modules without lineage information.
	 *  <group>1-allele</group>
	 */
	void setAlleleLineage(long lineage, size_t idx, int ploidy = -1, int chrom = -1);

	/** return an editable array (a \c carray object) that represents all
	 *  alleles of an individual. If \e ploidy or \e chroms is given, only
	 *  alleles on the specified chromosomes and homologous copy of chromosomes
	 *  will be returned. If multiple chromosomes are specified, there should
	 *  not be gaps between chromosomes. This function ignores type of
	 *  chromosomes so it will return unused alleles for sex and mitochondrial
	 *  chromosomes.
	 *  <group>2-genotype</group>
	 */
	PyObject * genotype(const uintList & ploidy = uintList(), const uintList & chroms = uintList());

	/** return an itertor that iterate through all mutants (non-zero alleles) of
	 *  an individual. Each mutant is presented as a tuple of (index, value)
	 *  where index is the index of mutant ranging from zero to totNumLoci() *
	 *  ploidy() - 1, so you will have to adjust indexes to check multiple alleles
	 *  at a locus. If \e ploidy or \e chroms is given, only alleles on the specified
	 *  chromosomes and homologous copy of chromosomes will be iterated. If multiple
	 *  chromosomes are specified, there should not be gaps between chromosomes. This
	 *  function ignores type of chromosomes so it will return unused alleles for sex
	 *  and mitochondrial chromosomes.
	 *  <group>2-genotype</group>
	 */
	pyMutantIterator mutants(const uintList & ploidy = uintList(), const uintList & chroms = uintList());


	/** return an editable array (a \c carray_lineage object) that represents the
	 *  lineages of all alleles of an individual. If \e ploidy or \e chroms is
	 *  given, only lineages on the specified chromosomes and homologous copy of
	 *  chromosomes will be returned. If multiple chromosomes are specified,
	 *  there should not be gaps between chromosomes. This function ignores
	 *  type of chromosomes so it will return lineage of unused alleles for sex
	 *  and mitochondrial chromosomes. A \c None object will be returned for
	 *  modules without lineage information.
	 *  <group>2-genotype</group>
	 */
	PyObject * lineage(const uintList & ploidy = uintList(), const uintList & chroms = uintList());

	/** CPPONLY
	 *  Return a Python dictionary with alleles at specified loci. This function
	 *  is usually used to collect alleles to send to a user-provided function.
	 */
	PyObject * mutAtLoci(const lociList & loci);

	/** CPPONLY
	 *  Return a Python tuple with alleles at specified loci. This function
	 *  is usually used to collect alleles to send to a user-provided function.
	 */
	PyObject * genoAtLoci(const lociList & loci);

	/** Fill the genotype of an individual using a list of alleles \e geno.
	 *  If parameters \e ploidy and/or \e chroms are specified, alleles will
	 *  be copied to only all or specified chromosomes on selected homologous
	 *  copies of chromosomes. \c geno will be reused if its length is less than
	 *  number of alleles to be filled. This function ignores type of
	 *  chromosomes so it will set genotype for unused alleles for sex and
	 *  mitochondrial chromosomes.
	 *  <group>2-genotype</group>
	 */
	void setGenotype(const uintList & geno, const uintList & ploidy = uintList(), const uintList & chroms = uintList());

	/** Fill the lineage of an individual using a list of IDs \e lineage.
	 *  If parameters \e ploidy and/or \e chroms are specified, lineages will
	 *  be copied to only all or specified chromosomes on selected homologous
	 *  copies of chromosomes. \c lineage will be reused if its length is less than
	 *  number of allelic lineage to be filled. This function ignores type of
	 *  chromosomes so it will set lineage to unused alleles for sex and
	 *  mitochondrial chromosomes. It does nothing for modules without lineage
	 *  information.
	 *  <group>2-genotype</group>
	 */
	void setLineage(const uintList & lineage, const uintList & ploidy = uintList(), const uintList & chroms = uintList());


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


	/// CPPONLY
	bool firstOffspring() const
	{
		return ISSETFLAG(m_flags, m_flagFirstOffspring);
	}


	/// CPPONLY
	void setFirstOffspring(bool first) const
	{
		if (first)
			SETFLAG(m_flags, m_flagFirstOffspring);
		else
			RESETFLAG(m_flags, m_flagFirstOffspring);
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
		size_t idx = field.empty() ? field.value() : infoIdx(field.name());

		CHECKRANGEINFO(idx);
		return m_infoPtr[idx];
	}


	/** Return the value of an information field \e field (by index or name)
	 *  as an integer number.
	 *  CPPONLY
	 */
	int intInfo(const uintString & field) const
	{
		size_t idx = field.empty() ? field.value() : infoIdx(field.name());

		CHECKRANGEINFO(idx);
		return static_cast<int>(m_infoPtr[idx]);
	}


	/** set the value of an information field \e field (by index or name) to
	 *  \e value. <tt>ind.setInfo(value, field)</tt> is equivalent to
	 *  <tt>ind.field = value</tt> although the function form allows the use
	 *  of indexes of information fieldes.
	 *  <group>5-info</group>
	 */
	void setInfo(double value, const uintString & field)
	{
		size_t idx = field.empty() ? field.value() : infoIdx(field.name());

		CHECKRANGEINFO(idx);
		m_infoPtr[idx] = value;
	}


	/// CPPONLY start of alleles
	GenoIterator genoBegin() const
	{
#ifdef MUTANTALLELE
		// Call + operator in order to find the correct m_com_index;
		// more information see operator+ function in vectorm.h
		return m_genoPtr + 0;
#else
		return m_genoPtr;
#endif
	}


	/// CPPONLY end of allele
	GenoIterator genoEnd() const
	{
		return m_genoPtr + genoSize();
	}


#ifdef LINEAGE
	/// CPPONLY start of lineage
	LineageIterator lineageBegin() const
	{
		return m_lineagePtr;
	}


	/// CPPONLY end of lineage
	LineageIterator lineageEnd() const
	{
		return m_lineagePtr + genoSize();
	}


#endif

	/// CPPONLY start of allele of the pth set of chromosome
	GenoIterator genoBegin(size_t p) const
	{
		CHECKRANGEPLOIDY(p);
		return m_genoPtr + p * totNumLoci();

	}


	/// CPPONLY end of allele of the pth set of chromosome
	GenoIterator genoEnd(size_t p) const
	{
		CHECKRANGEPLOIDY(p);
		return m_genoPtr + (p + 1) * totNumLoci();
	}


#ifdef LINEAGE
	/// CPPONLY start of lineage of the pth set of chromosome
	LineageIterator lineageBegin(size_t p) const
	{
		CHECKRANGEPLOIDY(p);
		return m_lineagePtr + p * totNumLoci();

	}


	/// CPPONLY end of lineage of the pth set of chromosome
	LineageIterator lineageEnd(size_t p) const
	{
		CHECKRANGEPLOIDY(p);
		return m_lineagePtr + (p + 1) * totNumLoci();
	}


#endif

	/// CPPONLY start of allele of the pth set of chromosome, chrom ch
	GenoIterator genoBegin(size_t p, size_t chrom) const
	{
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(chrom);
		return m_genoPtr + p * totNumLoci() + chromBegin(chrom);

	}


	/// CPPONLY end of allele of the pth set of chromosome
	GenoIterator genoEnd(size_t p, size_t chrom) const
	{
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(chrom);
		return m_genoPtr + p * totNumLoci() + chromEnd(chrom);

	}


#ifdef LINEAGE
	/// CPPONLY start of lineage of the pth set of chromosome, chrom ch
	LineageIterator lineageBegin(size_t p, size_t chrom) const
	{
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(chrom);
		return m_lineagePtr + p * totNumLoci() + chromBegin(chrom);

	}


	/// CPPONLY end of lineage of the pth set of chromosome
	LineageIterator lineageEnd(size_t p, size_t chrom) const
	{
		CHECKRANGEPLOIDY(p);
		CHECKRANGECHROM(chrom);
		return m_lineagePtr + p * totNumLoci() + chromEnd(chrom);

	}


#endif

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
	void display(ostream & out, int width = 1, const vectoru & loci = vectoru(),
		const vectoru & infoIdx = vectoru());

	//@}

private:
	bool validIndex(size_t idx) const;

	bool validIndex(size_t idx, size_t p) const;

	bool validIndex(size_t idx, size_t p, size_t ch) const;

	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const size_t /* version */) const
	{
		//  boost::serialization::base_object<GenoStruTrait>(*this));
		bool b;

		b = ISSETFLAG(m_flags, m_flagFemale);
		ar & b;

		b = ISSETFLAG(m_flags, m_flagAffected);
		ar & b;
	}


	template<class Archive>
	void load(Archive & ar, const size_t /* version */)
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

#ifdef LINEAGE
	/// pointer to lineage.
	LineageIterator m_lineagePtr;
#endif

	/// pointer to info
	InfoIterator m_infoPtr;
};


/** CPPONLY
 *  A class used to compare two individuals by one or more information fields.
 */
class indCompare
{
public:
	// accept the index to an information field
	indCompare(size_t idx, bool reverse=false) : m_fields(1, idx), m_reverse(reverse) {}

	indCompare(const vectoru & idx, bool reverse=false) : m_fields(idx), m_reverse(reverse) {}

	bool operator()(const Individual & lhs, const Individual & rhs) const
	{
		for (size_t i = 0; i < m_fields.size(); ++i) {
			double v1 = lhs.info(m_fields[i]);
			double v2 = rhs.info(m_fields[i]);
			if (v1 == v2)
				// if equal, look for the next field
				continue;
			else
				// if < or >, we have a comparison
				return m_reverse ? v1 > v2 : v1 < v2;
		}
		// they are actually equal, but we are inplementing operator <
		return false;
	}


private:
	vectoru m_fields;

	bool m_reverse;
};


/**
    this class implements a C++ iterator class that iterate through
    individuals in a (sub)population. If allInds are true, the
    visiblility of individuals will not be checked. Note that
    IndividualIterator *will* iterate through only visible individuals, and
    allInds is only provided when we know in advance that all individuals are
    visible. This is a way to obtain better performance in simple cases.
 */
template <typename T, typename PTR, typename REF>
class IndividualIterator
{
public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef typename T::value_type value_type;
	typedef long int difference_type;
	typedef REF reference;
	typedef PTR pointer;

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


	REF operator*() const
	{
		DBG_ASSERT(m_it < m_end, ValueError,
			"Can not refer to an invalid iterator");

		return *m_it;
	}


	PTR operator->() const
	{
		DBG_ASSERT(m_it < m_end, ValueError,
			"Can not refer to an invalid iterator");

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
		while (++m_it < m_end)
			if (m_it->visible())
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

		while (++m_it < m_end)
			if (m_it->visible())
				return *this;
		// in this case, m_it must equals m_end
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
			return IndividualIterator(m_end - m_it >= diff ? m_it + diff : m_end, m_end, m_allInds);
		IndividualIterator tmp(*this);
		DBG_ASSERT(tmp.m_it < tmp.m_end, ValueError,
			"Can not advance invalid iterator");
		difference_type i = 0;
		while (i < diff && ++tmp.m_it < tmp.m_end)
			if (tmp.m_it->visible())
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
		if (diff == 0)
			return *this;
		difference_type i = 0;
		while (++m_it < m_end) {
			if (m_it->visible()) {
				++i;
				if (i == diff)
					break;
			}
		}
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
typedef vector<Individual>::pointer RawIndPointer;
typedef vector<Individual>::const_pointer ConstRawIndPointer;
typedef vector<Individual>::reference RawIndReference;
typedef vector<Individual>::const_reference ConstRawIndReference;

typedef IndividualIterator<RawIndIterator, RawIndPointer, RawIndReference> IndIterator;
typedef IndividualIterator<ConstRawIndIterator, ConstRawIndPointer, ConstRawIndReference> ConstIndIterator;

/**
 *  this class implements a C++ iterator class that iterate through
 *  infomation fields in a (sub)population using
 *  1. an IndIterator that	will skip invisible individuals, or
 *  2. a gapped iterator that will run faster.
 *  Note that 1, 2 should yield identical result, and 2 should be used
 *  when there is no virtual subpopulation.q
 *
 *  CPPONLY
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


	InformationIterator(size_t info, InfoIterator ptr, size_t size)
		: m_info(info), m_useGappedIterator(true),
		m_it(), m_ptr(ptr), m_step(size)
	{
	}


	InformationIterator(size_t info, IndividualIterator<T, pointer, reference> it)
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
			return *(m_ptr + m_info);
		else
			return *(m_it->infoPtr() + m_info);
	}


	//
	pointer operator->() const
	{
		if (m_useGappedIterator)
			return m_ptr + m_info;
		else
			return m_it->infoPtr() + m_info;
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
	size_t m_info;
	///
	bool m_useGappedIterator;
	// Individual iterator
	IndividualIterator<T, pointer, reference> m_it;
	//
	InfoIterator m_ptr;
	//
	size_t m_step;
};


typedef InformationIterator<RawIndIterator> IndInfoIterator;
typedef InformationIterator<ConstRawIndIterator> ConstIndInfoIterator;

/** This class implements a C++ iterator class that iterate through
 *  all alleles in a (virtual) (sub)population using
 *  1. an IndIterator that will skip invisible individuals and invalid
 *      alleles, or
 *  2. a gapped iterator that will run faster, in the case that
 *    a): no virtual subpopulation
 *    b): not sex chromosomes
 *    c): not haplodiploid
 *
 *   CPPONLY
 */
template <typename T, typename ITER, typename REF>
class CombinedAlleleIterator
{
public:
	typedef std::forward_iterator_tag iterator_category;
	typedef Allele value_type;
	typedef long int difference_type;
	typedef REF reference;
	//typedef GenoIterator pointer;

	CombinedAlleleIterator()
	{
	}


	CombinedAlleleIterator(size_t shift, IndividualIterator<T, typename T::pointer, typename T::reference> it, ITER ptr, ITER ptrEnd, size_t size)
		: m_useGappedIterator(true), m_valid(true), m_shift(shift),
		m_ptr(ptr), m_ptrBegin(ptr), m_ptrEnd(ptrEnd), m_size(size),
		// ignored
		m_it(it), m_index(0), m_ploidy(0), m_chromType(0),
		m_haplodiploid(false), m_p(0)
	{
		m_valid = (m_ptr != m_ptrEnd) && it.valid();
		if (m_valid)
			m_ploidy = it->ploidy();
	}


	CombinedAlleleIterator(size_t idx, IndividualIterator<T, typename T::pointer, typename T::reference> it)
		: m_useGappedIterator(false), m_valid(true), m_shift(),
		m_ptr(), m_ptrBegin(), m_ptrEnd(), m_size(0), // belong to a previous one
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


	size_t currentPloidy()
	{
		if (m_useGappedIterator) {
			// NOTE: this iterator is used only when indOrdered() is set to true for
			// the whole population (see Population.lineageIterator()). It is therefore
			// possible to get the index of individual from index of m_ptr.
			//
			// There is a conversion from size_t to long (difference_type), a possible
			// loss of data
			difference_type offset = static_cast<difference_type>((m_ptr - m_ptrBegin) / m_size);
			return static_cast<size_t>(offset % m_ploidy);
		} else {
			return m_p;
		}
	}


	IndividualIterator<T, typename T::pointer, typename T::reference> individual()
	{
		if (m_useGappedIterator) {
			// NOTE: this iterator is used only when indOrdered() is set to true for
			// the whole population (see Population.lineageIterator()). It is therefore
			// possible to get the index of individual from index of m_ptr.
			//
			// There is a conversion from size_t to long (difference_type), a possible
			// loss of data
			difference_type offset = static_cast<difference_type>((m_ptr - m_ptrBegin) / (m_size * m_ploidy));
			return(m_it + offset);
		} else {
			return(m_it);
		}
	}


	void assignIfDiffer(Allele a)
	{
		if (m_useGappedIterator)
			return REF_ASSIGN_ALLELE(m_ptr + m_shift, a);
		else {
			DBG_ASSERT(m_it.valid(), SystemError, "Cannot refer to an invalid individual iterator");
			return REF_ASSIGN_ALLELE(m_it->genoBegin() + m_index + m_p * m_size, a);
		}
	}


#ifdef MUTANTALLELE
	// do not dereference which will create something ...
	Allele value()
	{
		if (m_useGappedIterator)
			return (m_ptr + m_shift).value();
		else {
			DBG_ASSERT(m_it.valid(), SystemError, "Cannot refer to an invalid individual iterator");
			return (m_it->genoBegin() + m_index + m_p * m_size).value();
		}
	}


#else
	Allele value()
	{
		if (m_useGappedIterator)
			return *(m_ptr + m_shift);
		else {
			DBG_ASSERT(m_it.valid(), SystemError, "Cannot refer to an invalid individual iterator");
			return *(m_it->genoBegin() + m_index + m_p * m_size);
		}
	}


	// this is the most important part!
	REF operator *()
	{
		if (m_useGappedIterator)
			return *(m_ptr + m_shift);
		else {
			DBG_ASSERT(m_it.valid(), SystemError, "Cannot refer to an invalid individual iterator");
			return *(m_it->genoBegin() + m_index + m_p * m_size);
		}
	}


#endif

	ITER ptr()
	{
		if (m_useGappedIterator)
			return m_ptr + m_shift;
		else
			return m_it->genoBegin() + m_index + m_p * m_size;
	}


	void advance(IndividualIterator<T, typename T::pointer, typename T::reference> & it, size_t & p, bool & valid)
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
			while ((++it).valid())
				if (it->sex() == MALE)
					break;
			p = 1;
			valid = it.valid();
		} else if (m_chromType == MITOCHONDRIAL) {
			// only the first homologous copy is valid
			DBG_ASSERT(p == 0, SystemError, "Only the first homologous copy of mitochondrial DNA can be iterated.");
			++it;
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
			if (m_ptrEnd > m_ptr + diff * m_size)
				m_ptr += diff * m_size;
			else {
				m_ptr = m_ptrEnd;
				m_valid = false;
			}
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
			return m_ptr != rhs.m_ptr || m_shift != rhs.m_shift;
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
	size_t m_shift;
	//
	ITER m_ptr;
	//
	ITER m_ptrBegin;
	//
	ITER m_ptrEnd;
	// genosize
	size_t m_size;
	//
	// The second iteration method
	// Individual iterator
	IndividualIterator<T, typename T::pointer, typename T::reference> m_it;
	// index of the locus
	size_t m_index;
	// overall ploidy
	size_t m_ploidy;
	// chromosome type
	size_t m_chromType;
	//
	bool m_haplodiploid;
	// current ploidy, used in individualiterator
	size_t m_p;
};


typedef CombinedAlleleIterator<RawIndIterator, GenoIterator, AlleleRef> IndAlleleIterator;
typedef CombinedAlleleIterator<ConstRawIndIterator, ConstGenoIterator, ConstAlleleRef> ConstIndAlleleIterator;

#ifdef LINEAGE
/* This is the lineageIterator version of the CombinedAlleleIterator. It cannot
 * be derived from the previous template because of the use of lineage specific
 * functions such as LineageBegin().
 */
template <typename T>
class CombinedLineageIterator
{
public:
	typedef std::forward_iterator_tag iterator_category;
	typedef long value_type;
	typedef long int difference_type;
	typedef long & reference;
	//typedef GenoIterator pointer;

	CombinedLineageIterator()
	{
	}


	CombinedLineageIterator(size_t shift, IndividualIterator<T, typename T::pointer, typename T::reference> it,
		LineageIterator ptr, LineageIterator ptrEnd, size_t size)
		: m_useGappedIterator(true), m_valid(true), m_shift(shift),
		m_ptr(ptr), m_ptrBegin(ptr), m_ptrEnd(ptrEnd), m_size(size), m_it(it),
		// ignored
		m_index(0), m_ploidy(0), m_chromType(0),
		m_haplodiploid(false), m_p(0)
	{
		m_valid = m_ptr != m_ptrEnd;
		m_ploidy = it->ploidy();
	}


	CombinedLineageIterator(size_t idx, IndividualIterator<T, typename T::pointer, typename T::reference> it)
		: m_useGappedIterator(false), m_valid(true), m_shift(),
		m_ptr(), m_ptrBegin(), m_ptrEnd(), m_size(0), // belong to a previous one
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


	size_t currentPloidy()
	{
		if (m_useGappedIterator) {
			// NOTE: this iterator is used only when indOrdered() is set to true for
			// the whole population (see Population.lineageIterator()). It is therefore
			// possible to get the index of individual from index of m_ptr.
			//
			// There is a conversion from size_t to long (difference_type), a possible
			// loss of data
			difference_type offset = static_cast<difference_type>((m_ptr - m_ptrBegin) / m_size);
			return static_cast<size_t>(offset % m_ploidy);
		} else {
			return m_p;
		}
	}


	IndividualIterator<T, typename T::pointer, typename T::reference> individual()
	{
		if (m_useGappedIterator) {
			// NOTE: this iterator is used only when indOrdered() is set to true for
			// the whole population (see Population.lineageIterator()). It is therefore
			// possible to get the index of individual from index of m_ptr.
			//
			// There is a conversion from size_t to long (difference_type), a possible
			// loss of data
			difference_type offset = static_cast<difference_type>((m_ptr - m_ptrBegin) / (m_size * m_ploidy));
			return(m_it + offset);
		} else {
			return(m_it);
		}
	}


	// this is the most important part!
	long & operator *() const
	{
		if (m_useGappedIterator)
			return *(m_ptr + m_shift);
		else {
			DBG_ASSERT(m_it.valid(), SystemError, "Cannot refer to an invalid individual iterator");
			return *(m_it->lineageBegin() + m_index + m_p * m_size);
		}
	}


	LineageIterator ptr()
	{
		if (m_useGappedIterator)
			return m_ptr + m_shift;
		else
			return m_it->lineageBegin() + m_index + m_p * m_size;
	}


	void advance(IndividualIterator<T, typename T::pointer, typename T::reference> & it, size_t & p, bool & valid)
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
			while ((++it).valid())
				if (it->sex() == MALE)
					break;
			p = 1;
			valid = it.valid();
		} else if (m_chromType == MITOCHONDRIAL) {
			// only the first homologous copy is valid
			DBG_ASSERT(p == 0, SystemError, "Only the first homologous copy of mitochondrial DNA can be iterated.");
			++it;
			valid = it.valid();
		}
	}


	// return, then advance.
	CombinedLineageIterator operator++(int)
	{
		// save current state
		CombinedLineageIterator tmp(*this);

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


	CombinedLineageIterator operator++()
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


	CombinedLineageIterator & operator+=(difference_type diff)
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


	CombinedLineageIterator operator+(difference_type diff)
	{
		CombinedLineageIterator tmp(*this);

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


	bool operator!=(const CombinedLineageIterator & rhs)
	{
		if (m_useGappedIterator)
			return m_ptr != rhs.m_ptr || m_shift != rhs.m_shift;
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
	size_t m_shift;
	//
	LineageIterator m_ptr;
	//
	LineageIterator m_ptrBegin;
	//
	LineageIterator m_ptrEnd;
	// genosize
	size_t m_size;
	//
	// The second iteration method
	// Individual iterator
	IndividualIterator<T, typename T::pointer, typename T::reference> m_it;
	// index of the locus
	size_t m_index;
	// overall ploidy
	size_t m_ploidy;
	// chromosome type
	size_t m_chromType;
	//
	bool m_haplodiploid;
	// current ploidy, used in individualiterator
	size_t m_p;
};

typedef CombinedLineageIterator<RawIndIterator> IndLineageIterator;
#endif


}

#ifndef SWIG
// set version for individual class
// version 0: base (reset for version 1.0)
BOOST_CLASS_VERSION(simuPOP::Individual, 0)
#endif

#endif
