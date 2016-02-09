/**
 *  $File: genoStru.h $
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

#ifndef _GENOSTRU_H
#define _GENOSTRU_H

/**
   \file
   \brief class genoStru and genoTrait
 */

#include "utility.h"
#include "simuPOP_cfg.h"

#include "boost_pch.hpp"
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
/// CPPONLY
/** genetic structure. Shared by individuals of one population

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
	GenoStructure() : m_ploidy(2), m_totNumLoci(0),
		m_numLoci(0), m_chromTypes(), m_chromX(-1), m_chromY(-1), m_mitochondrial(-1), m_customized(),
		m_haplodiploid(false), m_lociPos(0), m_chromIndex(0),
		m_chromNames(), m_alleleNames(), m_lociNames(), m_lociNameMap(), m_infoFields(0), m_lociPosMap(),
		m_refCount(0)
	{
	}


	/** CPPONLY \brief constructor. The ONLY way to construct this strucuture. There is not set... functions

	   \param ploidy number of sets of chromosomes
	   \param loci number of loci on each chromosome.
	   \param lociPos loci distance on each chromosome. the default values
	   are 1,2,etc.
	   \param chromNames chromosome names
	   \param alleleNames allele names
	   \param lociNames name of loci
	   \param length of info field
	 */
	GenoStructure(UINT ploidy, const vectoru & loci, const vectoru & chromTypes, bool haplodiploid,
		const vectorf & lociPos, const vectorstr & chromNames, const matrixstr & alleleNames,
		const vectorstr & lociNames, const vectorstr & infoFields);

	bool operator==(const GenoStructure & rhs);

	bool operator!=(const GenoStructure & rhs)
	{
		return !(*this == rhs);
	}


	/// destructor, do nothing.
	~GenoStructure()
	{
	}


	/// CPPONLY
	// Population::mergePopulationByLoci needs to hanle this
	double locusPos(UINT locus) const
	{
		return m_lociPos[locus];
	}


	/// CPPONLY
	// Population::mergePopulationByLoci needs to hanle this
	size_t chromIndex(UINT ch) const
	{
		return m_chromIndex[ch];
	}


	/// CPPONLY
	void buildLociPosMap() const;

	/// CPPONLY
	void setChromTypes(const vectoru & chromTypes);

private:
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const UINT /* version */) const
	{
		ar & m_ploidy;
		ar & m_numLoci;
		ar & m_chromTypes;
		ar & m_haplodiploid;
		ar & m_lociPos;
		ar & m_chromNames;
		ar & m_alleleNames;
		ar & m_lociNames;
		ar & m_infoFields;
		/// do not save load chromosome map
	}


	template<class Archive>
	void load(Archive & ar, const UINT /* version */)
	{

		ar & m_ploidy;
		ar & m_numLoci;

		ar & m_chromTypes;

		// set m_chromX etc.
		setChromTypes(m_chromTypes);
		// haplodiploid flag is introduced in 0.8.5
		ar & m_haplodiploid;
		//
		ar & m_lociPos;
		ar & m_chromNames;
		ar & m_alleleNames;
		ar & m_lociNames;
		ar & m_infoFields;

		m_lociNameMap.clear();
		if (!m_lociNames.empty()) {
			for (size_t i = 0; i < m_lociNames.size(); ++i) {
				if (!m_lociNames[i].empty())
					m_lociNameMap[m_lociNames[i]] = i;
			}
		}

		// build chromosome index
		m_chromIndex.resize(m_numLoci.size() + 1);
		m_chromIndex[0] = 0;
		for (size_t i = 1; i <= m_numLoci.size(); ++i)
			m_chromIndex[i] = m_chromIndex[i - 1] + m_numLoci[i - 1];

		m_totNumLoci = m_chromIndex[m_numLoci.size()];
		/// do not save load chromosome map
	}


	BOOST_SERIALIZATION_SPLIT_MEMBER();

	/// ploidy
	UINT m_ploidy;

	/// total number of loci
	size_t m_totNumLoci;

	/// number of loci
	vectoru m_numLoci;

	/// Type of each chromosome.
	vectoru m_chromTypes;

	/// index of chromosome X, -1 if not exist
	int m_chromX;

	/// index of chromosome Y, -1 if not exist
	int m_chromY;

	/// indexes of mitochondrial chromosomes
	int m_mitochondrial;

	/// indexes of customized chromosomes
	vectoru m_customized;

	/// whether or not this population is haplodiploid.
	bool m_haplodiploid;

	/// position of loci on chromosome, recommended with unit cM
	vectorf m_lociPos;

	/// loci index
	vectoru m_chromIndex;

	/// chromosome names
	vectorstr m_chromNames;

	/// allele names
	matrixstr m_alleleNames;

	/// loci names
	vectorstr m_lociNames;

	/// map of locinames
	map<string, size_t> m_lociNameMap;

	/// name of the information field
	vectorstr m_infoFields;

	mutable map<genomic_pos, size_t> m_lociPosMap;

	mutable UINT m_refCount;

	friend class GenoStruTrait;
};
}


#ifndef SWIG
// set version for GenoStructure class
// version 0: base (reset for 1.0)
BOOST_CLASS_VERSION(simuPOP::GenoStructure, 0)
#endif

namespace simuPOP {
/**
 *  All individuals in a population share the same genotypic properties such as
 *  number of chromosomes, number and position of loci, names of markers,
 *  chromosomes, and information fields. These properties are stored in this
 *  \c GenoStruTrait class and are accessible from both \c Individual and
 *  \c Population classes. Currently, a genotypic structure consists of
 *
 *  \li Ploidy, namely the number of homologous sets of chromosomes, of a
 *      population. Haplodiploid population is also supported.
 *  \li Number of chromosomes and number of loci on each chromosome.
 *  \li Positions of loci, which determine the relative distance between loci
 *      on the same chromosome. No unit is assumed so these positions can be
 *      ordinal (\c 1, \c 2, \c 3, ..., the default), in physical distance
 *      (\c bp, \c kb or \c mb), or in map distance (e.g. \c centiMorgan)
 *      depending on applications.
 *  \li Names of alleles, which can either be shared by all loci or be
 *      specified for each locus.
 *  \li Names of loci and chromosomes.
 *  \li Names of information fields attached to each individual.
 *
 *  In addition to basic property access functions, this class provides
 *  some utility functions such as \c locusByName, which looks up a locus by
 *  its name.
 */
class GenoStruTrait
{
public:
	/**
	 * A \c GenoStruTrait object is created with the construction of a
	 * \c Population object and cannot be initialized directly.
	 */
	GenoStruTrait() : m_genoStruIdx(MaxTraitIndex)
	{
	}


	/// CPPONLY set genotypic structure
	void setGenoStructure(UINT ploidy, const vectoru & loci, const vectoru & chromTypes, bool haplodiploid,
		const vectorf & lociPos, const vectorstr & chromNames, const matrixstr & alleleNames,
		const vectorstr & lociNames, const vectorstr & infoFields);

	/// CPPONLY set an existing geno structure
	/**
	   \note This is \em NOT efficient! However, this function has to be used when, for example,
	   loading a structure from a file.
	 */
	void setGenoStructure(const GenoStructure & rhs);

	/// CPPONLY set index directly
	void setGenoStruIdx(size_t idx)
	{
		DBG_FAILIF(idx >= s_genoStruRepository.size(), IndexError,
			(boost::format("Index %1%  to geno structure repository should be less than %2%") % idx %
			 s_genoStruRepository.size()).str());
		m_genoStruIdx = static_cast<TraitIndexType>(idx);
	}


	/// CPPONLY set geno stru
	void swapGenoStru(GenoStruTrait & rhs)
	{
		std::swap(m_genoStruIdx, rhs.m_genoStruIdx);
	}


	/** Return the distance between loci \e locus1 and \e locus2 on the
	 *  same chromosome. A negative value will be returned if \e locus1 is
	 *  after \e locus2.
	 *  <group>3-locus</group>
	 */
	double lociDist(size_t locus1, size_t locus2) const;

	/** HIDDEN
	 * return the number of loci left on that chromosome, including locus \c loc
	 *  <group>3-locus</group>
	 */
	size_t lociLeft(size_t locus) const;

	/** HIDDEN
	 *  Distance between locus \c locus and the last locus that is on the same
	 *  chromsome as \c locus.
	 *  <group>3-locus</group>
	 */
	double distLeft(size_t locus) const;

	/** HIDDEN
	 *  starting from \c locus, how many markers are covered by distance \c dist (>=0)
	 *  the result will be at least 1, even if dist = 0.
	 *  <group>3-locus</group>
	 */
	size_t lociCovered(size_t locus, double dist) const;

	/** CPPONLY
	 *  Add chromosomes from another genotypic structure and
	 *  create a new structure.
	 */
	const GenoStructure gsAddChromFromStru(size_t idx) const;

	/** CPPONLY
	 *  Add loci (merge loci on the same chromsoomes) from another genotypic
	 *  structure and create a new structure. index1 and index2 are used
	 *  to return the indexes of old loci in the new structure.
	 */
	const GenoStructure gsAddLociFromStru(size_t idx, vectoru & index1, vectoru & index2) const;

	/** CPPONLY
	 *  Add loci (merge loci on the same chromsoomes with identical name) from another genotypic
	 *  structure and create a new structure. index1 and index2 are used
	 *  to return the indexes of old loci in the new structure.
	 */
	const GenoStructure gsAddLociByNameFromStru(size_t idx, vectoru & index1, vectoru & index2) const;


	/** CPPONLY
	 *  Remove a list of loci from the current genotypic structure
	 *  and create a new structure.
	 */
	const GenoStructure gsRemoveLoci(const vectoru & kept);

	/** CPPONLY
	 *  add a new chromosome to genotype structure and create a new structure.
	 */
	const GenoStructure gsAddChrom(const vectorf & lociPos,
		const vectorstr & lociNames, const string & chromName,
		const matrixstr & alleleNames, size_t chromType) const;

	/** CPPONLY
	 *  Create a geno structure using new allele names.
	 */
	const GenoStructure gsSetAlleleNames(const lociList & loci, const matrixstr & alleleNames);

	/** CPPONLY
	 *  add some loci to genotype structure, newIndex
	 *  is used to return the indexes of these loci in the new
	 *  structure
	 */
	const GenoStructure gsAddLoci(const vectoru & chrom, const vectorf & pos,
		const vectorstr & lociNames, const matrixstr & alleleNames, vectoru & newIndex) const;

	/// CPPONLY return the GenoStructure
	GenoStructure & genoStru() const
	{
		return s_genoStruRepository[m_genoStruIdx];
	}


	/// CPPONLY return the GenoStructure index
	TraitIndexType genoStruIdx() const
	{
		return m_genoStruIdx;
	}


	/** return the number of homologous sets of chromosomes, specified by the
	 *  \e ploidy parameter of the \c Population function. Return 2 for a
	 *  haplodiploid population because two sets of chromosomes are stored
	 *  for both males and females in such a population.
	 *  <group>1-ploidy</group>
	 */
	UINT ploidy() const
	{

		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"Ploidy: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_ploidy;
	}


	/** return the ploidy name of this population, can be one of \c haploid,
	 *  \c diploid, \c haplodiploid, \c triploid, \c tetraploid or \c #-ploid
	 *  where \c # is the ploidy number.
	 *  <group>1-ploidy</group>
	 */
	string ploidyName() const;

	/** return the number of loci on chromosome \e chrom.
	 *  <group>3-locus</group>
	 */
	size_t numLoci(size_t chrom) const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"numLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGECHROM(chrom);
		return s_genoStruRepository[m_genoStruIdx].m_numLoci[chrom];
	}


	/** return a list of the number of loci on all chromosomes.
	 *  <group>3-locus</group>
	 */
	vectoru numLoci() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_numLoci;
	}


	/** CPPONLY
	 *  Return the index of chromosome X, \c -1 if there is no X chromosome.
	 *  <group>2-chromosome</group>
	 */
	int chromX() const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_chromX;
	}


	/** CPPONLY
	 *  Return the index of chromosome Y, \c -1 if there is no Y chromosome.
	 *  <group>2-chromosome</group>
	 */
	int chromY() const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_chromY;
	}


	/** CPPONLY
	 *  Return the indexes of the customized chromosomes.
	 *  <group>2-chromosome</group>
	 */
	vectoru customizedChroms() const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_customized;
	}


	/** CPPONLY
	 *  Return the indexes of mitochondrial chromosomes.
	 *  <group>2-chromosome</group>
	 */
	int mitochondrial() const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_mitochondrial;
	}


	/// HIDDEN
	bool sexChrom() const
	{
		return chromX() != -1;
	}


	/** HIDDEN becuase it can be replaced by ploidyName() == 'haplodiploid'
	 *  Return \c True if this population is haplodiploid.
	 *  <group>1-ploidy</group>
	 */
	bool isHaplodiploid() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_haplodiploid;
	}


	/** return the total number of loci on all chromosomes.
	 *  <group>3-locus</group>
	 */
	size_t totNumLoci() const
	{

		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_totNumLoci;
	}


	/** HIDDEN
	 *  return the total number of loci on all homologous chromosomes, which
	 *  is <tt>totNumLoci()*ploidy()</tt>.
	 */
	size_t genoSize() const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_ploidy * s_genoStruRepository[m_genoStruIdx].m_totNumLoci;
	}


	/** return the position of locus \e locus specified by the \e lociPos
	 *  parameter of the \c Population function.
	 *  <group>3-locus</group>
	 */
	double locusPos(size_t locus) const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"locusPos: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGEABSLOCUS(locus);
		return s_genoStruRepository[m_genoStruIdx].m_lociPos[locus];
	}


	/** return the positions of all loci, specified by the \e lociPos prameter
	 *  of the \c Population function. The default positions are 1, 2, 3, 4, ...
	 *  on each chromosome.
	 *  <group>3-locus</group>
	 */
	vectorf lociPos() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_lociPos;
	}


	/** return the number of chromosomes.
	 *  <group>2-chromosome</group>
	 */
	size_t numChrom() const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"numChrom: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_numLoci.size();
	}


	/// CPPONLY return an array of chromosome indexes
	const vectoru & chromIndex() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_chromIndex;
	}


	/** return the index of the first locus on chromosome \e chrom.
	 *  <group>2-chromosome</group>
	 */
	size_t chromBegin(size_t chrom) const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"chromBegin: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGECHROM(chrom);

		return s_genoStruRepository[m_genoStruIdx].m_chromIndex[chrom];
	}


	/** return the index of the last locus on chromosome \e chrom plus 1.
	 *  <group>2-chromosome</group>
	 */
	size_t chromEnd(size_t chrom) const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"chromEnd: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGECHROM(chrom);

		return s_genoStruRepository[m_genoStruIdx].m_chromIndex[chrom + 1];
	}


	/** return the absolute index of locus \e locus on chromosome \e chrom.
	 *  c.f. \c chromLocusPair.
	 *  <group>3-locus</group>
	 */
	size_t absLocusIndex(UINT chrom, UINT locus) const
	{
		CHECKRANGECHROM(chrom);
		CHECKRANGELOCUS(chrom, locus);

		return s_genoStruRepository[m_genoStruIdx].m_chromIndex[chrom] + locus ;
	}


	/**
	 * return the chromosome and relative index of a locus using its absolute
	 * index \e locus. c.f. \c absLocusIndex.
	 *  <group>3-locus</group>
	 */
	pairu chromLocusPair(size_t locus) const;


	/**
	 * return the name of a chromosome \e chrom.
	 *  <group>2-chromosome</group>
	 */
	string chromName(const size_t chrom) const
	{
		DBG_FAILIF(chrom >= s_genoStruRepository[m_genoStruIdx].m_numLoci.size(), IndexError,
			(boost::format("Chromosome index %1% out of range of 0 ~ %2%") % chrom %
			 s_genoStruRepository[m_genoStruIdx].m_numLoci.size()).str());

		return s_genoStruRepository[m_genoStruIdx].m_chromNames[chrom];
	}


	/** return a list of the names of all chromosomes.
	 *  <group>2-chromosome</group>
	 */
	vectorstr chromNames() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_chromNames;
	}


	/** return the type of a chromosome \e chrom (\c CUSTOMIZED, \c AUTOSOME,
	 *  \c CHROMOSOME_X, \c CHROMOSOME_Y or \c MITOCHONDRIAL.
	 *  <group>2-chromosome</group>
	 */
	size_t chromType(const size_t chrom) const
	{
		DBG_FAILIF(chrom >= s_genoStruRepository[m_genoStruIdx].m_numLoci.size(), IndexError,
			(boost::format("Chromosome index %1% out of range of 0 ~ %2%") % chrom %
			 s_genoStruRepository[m_genoStruIdx].m_numLoci.size()).str());

		return s_genoStruRepository[m_genoStruIdx].m_chromTypes[chrom];
	}


	/** return the type of all chromosomes (\c CUSTOMIZED, \c AUTOSOME,
	 *  \c CHROMOSOME_X, \c CHROMOSOME_Y, or \c MITOCHONDRIAL).
	 *  <group>2-chromosome</group>
	 */
	vectoru chromTypes() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_chromTypes;
	}


	/** return the index of a chromosome by its \e name.
	 *  <group>2-chromosome</group>
	 */
	size_t chromByName(const string name) const
	{
		const vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_chromNames;
		vectorstr::const_iterator it = std::find(names.begin(), names.end(), name);

		if (it == names.end())
			throw ValueError("Failed to find chrom with name " + name);
		return it - names.begin();
	}


	/** return the name of allele \e allele at \e lcous specified by the
	 *  \e alleleNames parameter of the \c Population function. \e locus could
	 *  be ignored if alleles at all loci share the same names. If the name of
	 *  an allele is unspecified, its value (\c '0', \c '1', \c '2', etc) is
	 *  returned.
	 *  <group>4-allele</group>
	 */
	string alleleName(const ULONG allele, const size_t locus = 0) const;

	/** CPPONLY
	 *  Return all allele names
	 */
	matrixstr allAlleleNames() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_alleleNames;
	}


	/** return a list of allele names at \locus given by the \e alleleNames
	 *  parameter of the \c Population function. \e locus could be ignored if
	 *  alleles at all loci share the same names. This list does not have to
	 *  cover all possible allele states of a population so
	 *  <tt>alleleNames()[</tt><em>allele</em><tt>]</tt> might fail
	 *  (use <tt>alleleNames(</tt><em>allele</em><tt>)</tt> instead).
	 *  <group>4-allele</group>
	 */
	vectorstr alleleNames(const size_t locus = 0) const;

	/** return the name of locus \e locus specified by the \e lociNames parameter of
	 *  the \c Population function. An empty string will be returned if no name
	 *  has been given to locus \e locus.
	 *  <group>3-locus</group>
	 */
	string locusName(const size_t locus) const
	{
		DBG_FAILIF(locus >= s_genoStruRepository[m_genoStruIdx].m_totNumLoci, IndexError,
			(boost::format("Locus index %1% out of range of 0 ~ %2%") % locus %
			 s_genoStruRepository[m_genoStruIdx].m_totNumLoci).str());

		const vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_lociNames;
		return names.empty() ? string() : names[locus];
	}


	/** return the names of all loci specified by the \e lociNames parameter of
	 *  the \c Population function. An empty list will be returned if
	 *  \e lociNames was not specified.
	 *  <group>3-locus</group>
	 */
	vectorstr lociNames() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_lociNames;
	}


	/** return the index of a locus with name \e name. Raise a \c ValueError
	 *  if no locus is found. Note that empty strings are used for loci without
	 *  name but you cannot lookup such loci using this function.
	 *  <group>3-locus</group>
	 */
	size_t locusByName(const string name) const
	{
		const map<string, size_t> & names = s_genoStruRepository[m_genoStruIdx].m_lociNameMap;

		map<string, size_t>::const_iterator it = names.find(name);

		if (it == names.end())
			throw ValueError("Failed to find locus with name " + name);
		return it->second;
	}


	/** return the indexes of loci with names \e names. Raise a \c ValueError
	 *  if any of the loci cannot be found.
	 *  <group>3-locus</group>
	 */
	vectoru lociByNames(const vectorstr & names) const;

	/** return the indexes of loci with positions \e positions (list of (chr, pos)
	 *  pairs). Raise a \c ValueError if any of the loci cannot be found.
	 */
	vectoru indexesOfLoci(const lociList & loci = lociList()) const;

	/** CPPONLY return the indexes of loci with positions \e positions (list of (chr, pos)
	 *  pairs). Raise a \c ValueError if any of the loci cannot be found.
	 */
	vectoru lociByPos(const vectorpos & positions) const;

	/** HIDDEN
	    Return \c True if \c name is one of the information fields of this population.
	 */
	bool hasInfoField(const string & name) const
	{
		vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_infoFields;

		return std::find(names.begin(), names.end(), name) != names.end();
	}


	/// HIDDEN obtain the number of information fields
	size_t infoSize() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_infoFields.size();
	}


	/** return a list of the names of all information fields of the population.
	 *  <group>5-info</group>
	 */
	vectorstr infoFields() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_infoFields;
	}


	/** return the name of information field \e idx.
	 *  <group>5-info</group>
	 */
	string infoField(size_t idx) const
	{
		CHECKRANGEINFO(idx);
		return s_genoStruRepository[m_genoStruIdx].m_infoFields[idx];
	}


	/** return the index of information field \e name. Raise an \c IndexError
	 * if \e name is not one of the information fields.
	 *  <group>5-info</group>
	 */
	size_t infoIdx(const string & name) const;

	/// CPPONLY add a new information field
	/**
	   \note Should only be called by Population::requestInfoField.
	   Right now, do not allow dynamic addition of these fields.
	 */
	const GenoStructure gsAddInfoFields(const vectorstr & fields);

	/// CPPONLY should only be called from population
	const GenoStructure gsSetInfoFields(const vectorstr & fields);

	/// CPPONLY swap a geno structure with the current one
	void swap(GenoStruTrait & rhs)
	{
		std::swap(m_genoStruIdx, rhs.m_genoStruIdx);
	}


	/// CPPONLY Increase the reference count of this structure
	void incGenoStruRef() const
	{
		++s_genoStruRepository[m_genoStruIdx].m_refCount;
		DBG_DO(DBG_POPULATION, cerr << "Inc ref of " << int(m_genoStruIdx) << " to "
			                        << s_genoStruRepository[m_genoStruIdx].m_refCount << endl);
	}


	/// CPPONLY Decrease the reference count of this structure
	void decGenoStruRef() const
	{
		DBG_FAILIF(s_genoStruRepository[m_genoStruIdx].m_refCount == 0, SystemError,
			(boost::format("Unknow error for reference counting of genotypic structure %1%.") % int(m_genoStruIdx)).str());
		--s_genoStruRepository[m_genoStruIdx].m_refCount;
		DBG_DO(DBG_POPULATION, cerr << "Dec ref of " << int(m_genoStruIdx) << " to "
			                        << s_genoStruRepository[m_genoStruIdx].m_refCount << endl);
	}


private:
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive &, const UINT)
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
