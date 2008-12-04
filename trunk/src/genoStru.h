/************************************************************************
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
#if  defined (_WIN32) || defined (__WIN32__)
#  include <boost/archive/text_iarchive.hpp>
#  include <boost/archive/text_oarchive.hpp>
#  include <fstream>
using std::ofstream;
using std::ifstream;
#endif                                                                                    // win32

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
		m_numLoci(0), m_chromTypes(), m_chromX(-1), m_chromY(-1), m_mitochondrial(),
		m_haplodiploid(false), m_lociPos(0), m_chromIndex(0),
		m_chromNames(), m_alleleNames(), m_lociNames(), m_infoFields(0)
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
		const vectorf & lociPos, const vectorstr & chromNames, const vectorstr & alleleNames,
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
	// population::mergePopulationByLoci needs to hanle this
	double locusPos(UINT locus) const
	{
		return m_lociPos[locus];
	}


	/// CPPONLY
	// population::mergePopulationByLoci needs to hanle this
	const size_t chromIndex(UINT ch) const
	{
		return m_chromIndex[ch];
	}


#if  defined (_WIN32) || defined (__WIN32__)

	// due to an weird compiling error fo vc7.1,
	// if I do not specify these two functions, the ar & process
	// will fail to compile.
	// This will only be defined for win32 system
	/// CPPONLY
	void saveStru(string filename)
	{
		ofstream ofs(filename.c_str());
		boost::archive::text_oarchive oa(ofs);

		oa << *this;
	}


	/// CPPONLY
	void loadStru(string filename)
	{
		ifstream ifs(filename.c_str());

		boost::archive::text_iarchive ia(ifs);

		ia >> *this;
	}


#endif                                                                    // win32

	/// CPPONLY
	void setChromTypes(const vectoru & chromTypes);

private:
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const UINT version) const
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
	void load(Archive & ar, const UINT version)
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

		// build chromosome index
		m_chromIndex.resize(m_numLoci.size() + 1);
		ULONG i;
		for (m_chromIndex[0] = 0, i = 1; i <= m_numLoci.size(); ++i)
			m_chromIndex[i] = m_chromIndex[i - 1] + m_numLoci[i - 1];

		m_totNumLoci = m_chromIndex[m_numLoci.size()];
		/// do not save load chromosome map
	}


	BOOST_SERIALIZATION_SPLIT_MEMBER();

	/// ploidy
	UINT m_ploidy;

	/// total number of loci
	UINT m_totNumLoci;

	/// number of loci
	vectoru m_numLoci;

	/// Type of each chromosome.
	vectoru m_chromTypes;

	/// index of chromosome X, -1 if not exist
	int m_chromX;

	/// index of chromosome Y, -1 if not exist
	int m_chromY;

	/// indexes of mitochondrial chromosome
	vectoru m_mitochondrial;

	/// whether or not this population is haplodiploid.
	bool m_haplodiploid;

	/// position of loci on chromosome, recommended with unit cM
	vectorf m_lociPos;

	/// loci index
	vectoru m_chromIndex;

	/// chromosome names
	vectorstr m_chromNames;

	/// allele names
	vectorstr m_alleleNames;

	/// loci names
	vectorstr m_lociNames;

	/// name of the information field
	vectorstr m_infoFields;

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
 *  \c GenoStruTrait class and are accessible from \c individual,
 *  \c population, and \c simulator classes. Currently, a genotypic structure
 *  consists of
 *
 *  \li Ploidy, namely the number of homologous sets of chromosomes, of a
 *      population. Haplodiploid population is also supported.
 *  \li Number of chromosomes and number of loci on each chromosome.
 *  \li Positions of loci, which determine the relative distance between loci
 *      on the same chromosome. No unit is assumed so these positions can be
 *      ordinal (\c 1, \c 2, \c 3, ..., the default), in physical distance
 *      (\c bp, \c kb or \c mb), or in map distance (e.g. centiMorgan)
 *      depending on applications.
 *  \li Names of alleles. Although alleles at different loci usually have
 *      different names, simuPOP uses the same names for alleles across loci
 *      for simplicity.
 *  \li Names of loci and chromosomes.
 *  \li Names of information fields attached to each individual.
 *
 *  In addition to basic property access functions, this class also provides
 *  some utility functions such as \c locusByName, which looks up a locus by
 *  its name.
 */
class GenoStruTrait
{
public:
	/**
	 * A \c GenoStruTrait object is created with the creation of a \c population
	 * so it cannot be initialized directly.
	 */
	GenoStruTrait() : m_genoStruIdx(MaxTraitIndex)
	{
	}


	/// CPPONLY set genotypic structure
	void setGenoStructure(UINT ploidy, const vectoru & loci, const vectoru & chromTypes, bool haplodiploid,
		const vectorf & lociPos, const vectorstr & chromNames, const vectorstr & alleleNames,
		const vectorstr & lociNames, const vectorstr & infoFields);

	/// CPPONLY set an existing geno structure
	/**
	   \note This is \em NOT efficient! However, this function has to be used when, for example,
	   loading a structure from a file.
	 */
	void setGenoStructure(GenoStructure & rhs);

	/// CPPONLY set index directly
	void setGenoStruIdx(size_t idx)
	{
		DBG_FAILIF(idx >= s_genoStruRepository.size(), IndexError,
			"Index " + toStr(idx) + " to geno structure repository should be less than " +
			toStr(s_genoStruRepository.size() ) );
		m_genoStruIdx = static_cast<TraitIndexType>(idx);
	}


	/** Return the distance between loci \e loc1 and \e loc2 on the same
	 *  chromosome. A negative value will be returned if \e loc1 is after
	 *  \e loc2.
	 *  <group>3-locus</group>
	 */
	double lociDist(UINT loc1, UINT loc2) const;

	/** HIDDEN
	 * return the number of loci left on that chromosome, including locus \c loc
	 *  <group>3-locus</group>
	 */
	UINT lociLeft(UINT loc) const;

	/** HIDDEN
	 *  Distance between locus \c loc and the last locus that is on the same
	 *  chromsome as \c loc.
	 *  <group>3-locus</group>
	 */
	double distLeft(UINT loc) const;

	/** HIDDEN
	 *  starting from \c loc, how many markers are covered by distance \c dist (>=0)
	 *  the result will be at least 1, even if dist = 0.
	 *  <group>3-locus</group>
	 */
	UINT lociCovered(UINT loc, double dist) const;

	/** CPPONLY
	 *  Add chromosomes from another genotypic structure and
	 *  create a new structure.
	 */
	GenoStructure & gsAddChromFromStru(size_t idx) const;

	/** CPPONLY
	 *  Add loci (merge loci on the same chromsoomes) from another genotypic
	 *  structure and create a new structure.
	 */
	GenoStructure & gsAddLociFromStru(size_t idx) const;

	/** CPPONLY
	 *  Remove a list of loci from the current genotypic structure
	 *  and create a new structure. Use \e kept to return indexes
	 *  of the remaining loci. However, if a non-empty kept is passed
	 *  it will be used directly (loci is ignored).
	 */
	GenoStructure & gsRemoveLoci(const vectoru & loci, vectoru & kept);

	/** CPPONLY
	 *  add a new chromosome to genotype structure and create a new structure.
	 */
	GenoStructure & gsAddChrom(const vectorf & lociPos,
		const vectorstr & lociNames, const string & chromName, UINT chromType) const;

	/** CPPONLY
	 *  add some loci to genotype structure, newIndex
	 *  is used to return the indexes of these loci in the new
	 *  structure
	 */
	GenoStructure & gsAddLoci(const vectoru & chrom, const vectorf & pos,
		const vectorstr & names, vectoru & newIndex) const;

	/// CPPONLY return the GenoStructure
	GenoStructure & genoStru() const
	{
		return s_genoStruRepository[m_genoStruIdx];
	}


	/// CPPONLY return the GenoStructure index
	size_t genoStruIdx() const
	{
		return static_cast<size_t>(m_genoStruIdx);
	}


	/** return the number of homologous sets of chromosomes, specified by the
	 *  \e ploidy parameter of the \c population function. Return 2 for a
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

	/** return the number of loci on chromosome \e chrom, equivalent to
	 *  <tt>numLoci()[</tt><em>chrom</em></tt>]</tt>.
	 *  <group>3-locus</group>
	 */
	UINT numLoci(UINT chrom) const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"numLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGECHROM(chrom);
		return s_genoStruRepository[m_genoStruIdx].m_numLoci[chrom];
	}


	/** return the number of loci on all chromosomes.
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
	 *  Return the indexes of the mitochondrial chromosomes.
	 *  <group>2-chromosome</group>
	 */
	vectoru mitochondrial() const
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
	UINT totNumLoci() const
	{

		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_totNumLoci;
	}


	/** HIDDEN
	 *  return the total number of loci on all homologous chromosomes, which
	 *  is <tt>totNumLoci()*ploidy()</tt>.
	 */
	UINT genoSize() const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_ploidy * s_genoStruRepository[m_genoStruIdx].m_totNumLoci;
	}


	/** return the position of locus \e loc specified by the \e lociPos
	 *  parameter of the \c population function. An \c IndexError will be
	 *  raised if the absolute index \e loc is greater than or equal to
	 *  the total number of loci.
	 *  <group>3-locus</group>
	 */
	double locusPos(UINT loc) const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"locusPos: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGEABSLOCUS(loc);
		return s_genoStruRepository[m_genoStruIdx].m_lociPos[loc];
	}


	/** return the positions of all loci, specified by the \e lociPos prameter
	 *  of the \c population function. The default positions are 1, 2, 3, 4, ...
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
	UINT numChrom() const
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
	UINT chromBegin(UINT chrom) const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"chromBegin: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGECHROM(chrom);

		return s_genoStruRepository[m_genoStruIdx].m_chromIndex[chrom];
	}


	/** return the index of the last locus on chromosome \e chrom plus 1.
	 *  <group>2-chromosome</group>
	 */
	UINT chromEnd(UINT chrom) const
	{
		DBG_FAILIF(m_genoStruIdx == MaxTraitIndex, SystemError,
			"chromEnd: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGECHROM(chrom);

		return s_genoStruRepository[m_genoStruIdx].m_chromIndex[chrom + 1];
	}


	/** return the absolute index of locus \e locus on chromosome \e chrom.
	 *  An \c IndexError will be raised if \e chrom or \e locus is out of
	 *  range. c.f. \c chromLocusPair.
	 *  <group>3-locus</group>
	 */
	UINT absLocusIndex(UINT chrom, UINT locus)
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
	std::pair<UINT, UINT> chromLocusPair(UINT locus) const;


	/**
	 * return the name of a chromosome \e chrom. Default to \c chrom# where #
	 * is the 1-based index of the chromosome.
	 *  <group>2-chromosome</group>
	 */
	string chromName(const UINT chrom) const
	{
		DBG_FAILIF(chrom >= s_genoStruRepository[m_genoStruIdx].m_numLoci.size(), IndexError,
			"Chromosome index " + toStr(chrom) + " out of range of 0 ~ " +
			toStr(s_genoStruRepository[m_genoStruIdx].m_numLoci.size()));

		return s_genoStruRepository[m_genoStruIdx].m_chromNames[chrom];
	}


	/** return a list of the names of all chromosomes.
	 *  <group>2-chromosome</group>
	 */
	vectorstr chromNames() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_chromNames;
	}


	/** return the type of a chromosome \e chrom (\c 1 for \c Autosome, \c 2 for
	 *  \c ChromosomeX, \c 3 for \c ChromosomeY, and \c 4 for \c Mitochondrial).
	 *  <group>2-chromosome</group>
	 */
	int chromType(const UINT chrom) const
	{
		DBG_FAILIF(chrom >= s_genoStruRepository[m_genoStruIdx].m_numLoci.size(), IndexError,
			"Chromosome index " + toStr(chrom) + " out of range of 0 ~ " +
			toStr(s_genoStruRepository[m_genoStruIdx].m_numLoci.size()));

		return s_genoStruRepository[m_genoStruIdx].m_chromTypes[chrom];
	}


	/** return the type of all chromosomes (\c 1 for \c Autosome, \c 2 for
	 *  \c ChromosomeX, \c 3 for \c ChromosomeY, and \c 4 for \c Mitochondrial).
	 *  <group>2-chromosome</group>
	 */
	vectoru chromTypes() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_chromTypes;
	}


	/** return the index of a chromosome by its \e name.
	 *  <group>2-chromosome</group>
	 */
	UINT chromByName(const string name) const
	{
		const vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_chromNames;
		vectorstr::const_iterator it = std::find(names.begin(), names.end(), name);

		if (it == names.end())
			throw ValueError("Failed to find chrom with name " + name);
		return it - names.begin();
	}


	/** return the name of allele \e allele specified by the \e alleleNames parameter of
	 *  the \c population function. If the name of an allele is not specified, its
	 *  index (\c '0', \c '1', \c '2', etc) is returned. An \c IndexError will be
	 *  raised if \e allele is larger than the maximum allowed allele state of this
	 *  module (<tt>MaxAllele()</tt>).
	 *  <group>4-allele</group>
	 */
	string alleleName(const UINT allele) const;

	/** return a list of allele names given by the \e alleleNames parameter of the
	 *  \c population function. This list does not have to cover all possible allele
	 *  states of a population so <tt>alleleNames()[</tt><em>allele</em><tt>]</tt>
	 *  might fail (use <tt>alleleNames(</tt><em>allele</em><tt>)</tt> instead).
	 *  <group>4-allele</group>
	 */
	vectorstr alleleNames() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_alleleNames;
	}


	/** return the name of locus \e loc specified by the \e lociNames parameter of
	 *  the \c population function. Default to \c locX-Y where \c X and \c Y
	 *  are 1-based chromosome and locus indexes (\c loc1-1, \c loc1-2, ... etc)
	 *  <group>3-locus</group>
	 */
	string locusName(const UINT loc) const
	{
		DBG_FAILIF(loc >= s_genoStruRepository[m_genoStruIdx].m_totNumLoci, IndexError,
			"Locus index " + toStr(loc) + " out of range of 0 ~ " +
			toStr(s_genoStruRepository[m_genoStruIdx].m_totNumLoci));

		return s_genoStruRepository[m_genoStruIdx].m_lociNames[loc];
	}


	/** return the names of all loci specified by the \e lociNames parameter of
	 *  the \c population function.
	 *  <group>3-locus</group>
	 */
	vectorstr lociNames() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_lociNames;
	}


	/** return the index of a locus with name \e name. Raise a \c ValueError
	 *  if no locus is found.
	 *  <group>3-locus</group>
	 */
	UINT locusByName(const string name) const
	{
		const vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_lociNames;
		vectorstr::const_iterator it = std::find(names.begin(), names.end(), name);

		if (it == names.end())
			throw ValueError("Failed to find locus with name " + name);
		return it - names.begin();
	}


	/** return the indexes of loci with names \e names. Raise a \c ValueError
	 *  if any of the loci cannot be found.
	 *  <group>3-locus</group>
	 */
	vectoru lociByNames(const vectorstr & names) const
	{
		vectoru indexes;

		for (vectorstr::const_iterator name = names.begin(); name != names.end(); ++name)
			indexes.push_back(locusByName(*name));
		return indexes;
	}


	/** HIDDEN
	    Return \c True if \c name is one of the information fields of this population.
	 */
	bool hasInfoField(const string & name) const
	{
		vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_infoFields;

		return std::find(names.begin(), names.end(), name) != names.end();
	}


	/// HIDDEN obtain the number of information fields
	UINT infoSize() const
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
	string infoField(UINT idx) const
	{
		CHECKRANGEINFO(idx);
		return s_genoStruRepository[m_genoStruIdx].m_infoFields[idx];
	}


	/** return the index of information field \e name. Raise an \c IndexError
	 * if \e name is not one of the information fields.
	 *  <group>5-info</group>
	 */
	UINT infoIdx(const string & name) const;

	/// CPPONLY add a new information field
	/**
	   \note Should only be called by population::requestInfoField.
	   Right now, do not allow dynamic addition of these fields.
	 */
	GenoStructure & struAddInfoFields(const vectorstr & fields);

	/// CPPONLY should should only be called from population
	GenoStructure & struSetInfoFields(const vectorstr & fields);

	/// CPPONLY swap a geno structure with the current one
	void swap(GenoStruTrait & rhs)
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
#endif
