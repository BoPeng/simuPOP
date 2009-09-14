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
/// CPPONLY
/** \brief genetic structure. Shared by individuals of one population

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
	GenoStructure() : m_ploidy(2), m_totNumLoci(0), m_genoSize(0), m_numChrom(0),
		m_numLoci(0), m_sexChrom(false), m_haplodiploid(false), m_lociPos(0), m_chromIndex(0),
		m_chromNames(), m_alleleNames(), m_lociNames(), m_maxAllele(), m_infoFields(0)
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
	 \param maxAllele maximum possible allele number for all alleles.
	 \param length of info field
	 */
	GenoStructure(UINT ploidy, const vectoru & loci, bool sexChrom, bool haplodiploid,
	              const vectorf & lociPos, const vectorstr & chromNames, const vectorstr & alleleNames,
	              const vectorstr & lociNames, UINT maxAllele, const vectorstr & infoFields);

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
		boost::archive::binary_oarchive oa(ofs);

		oa << boost::serialization::make_nvp("geno_structure", *this);
	}


	/// CPPONLY
	void loadStru(string filename)
	{
		ifstream ifs(filename.c_str());

		boost::archive::binary_iarchive ia(ifs);

		ia >> boost::serialization::make_nvp("geno_structure", *this);
	}


#endif                                                                    // win32

private:
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const UINT version) const
	{
		ar & make_nvp("ploidy", m_ploidy);
		ar & make_nvp("num_of_chrom", m_numChrom);
		ar & make_nvp("num_of_loci_on_each_chrom", m_numLoci);
		ar & make_nvp("sex_chromosome", m_sexChrom);
		ar & make_nvp("haplodiploid", m_haplodiploid);
		ar & make_nvp("loci_distance_on_chrom", m_lociPos);
		ar & make_nvp("chrom_name", m_chromNames);
		ar & make_nvp("allele_name", m_alleleNames);
		ar & make_nvp("loci_name", m_lociNames);
		ar & make_nvp("max_allele", m_maxAllele);
		ar & make_nvp("info_name", m_infoFields);
		/// do not save load chromosome map
	}


	template<class Archive>
	void load(Archive & ar, const UINT version)
	{

		ar & make_nvp("ploidy", m_ploidy);
		ar & make_nvp("num_of_chrom", m_numChrom);
		ar & make_nvp("num_of_loci_on_each_chrom", m_numLoci);

		// after simuPOP 0.6.8, we have m_sexChrom
		// before that, there is no sex chromosome
		if (version >= 1)
			ar & make_nvp("sex_chromosome", m_sexChrom);
		else
			m_sexChrom = false;
		// haplodiploid flag is introduced in 0.8.5
		if (version >= 4)
			ar & make_nvp("haplodiploid", m_haplodiploid);
		else
			m_haplodiploid = false;
		//
		ar & make_nvp("loci_distance_on_chrom", m_lociPos);
		if (version >= 3)
			ar & make_nvp("chrom_name", m_chromNames);
		else {
			m_chromNames.resize(m_numChrom);
			for (size_t i = 0; i < m_numChrom; ++i)
				m_chromNames[i] = "chrom" + toStr(i + 1);
		}
		ar & make_nvp("allele_name", m_alleleNames);
		ar & make_nvp("loci_name", m_lociNames);
		ar & make_nvp("max_allele", m_maxAllele);
		if (version >= 2)
			ar & make_nvp("info_name", m_infoFields);

		// build chromosome index
		m_chromIndex.resize(m_numLoci.size() + 1);
		ULONG i;
		for (m_chromIndex[0] = 0, i = 1; i <= m_numChrom; ++i)
			m_chromIndex[i] = m_chromIndex[i - 1] + m_numLoci[i - 1];

		m_totNumLoci = m_chromIndex[m_numChrom];
		m_genoSize = m_totNumLoci * m_ploidy;
		/// do not save load chromosome map
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
// version 2: add info name
// version 3: add chromName
// version 4: add haplodiploid
BOOST_CLASS_VERSION(simuPOP::GenoStructure, 4)
#endif

namespace simuPOP {
/// genotypic structure related functions, can be accessed from individuals, populations and simulator levels.
/**
   Genotypic structure refers to the number of chromosomes, the
   number and position of loci on each chromosome, and allele and locus names etc. All individuals
   in a population share the same genotypic structure. Because class \c GenoStruTrait
   is inherited by class \c population, class \c individual, and class \c simulator,
   functions provided in this class can be accessed at the individual, population and
   simulator levels. This object can not be created directly. It is created by a population.
 */
class GenoStruTrait
{
private:
#define TraitIndexType unsigned char
#define TraitMaxIndex 0xFF

public:
	///
	/**
	 \test src_genoStruTrait.log Genotypic structure
	 */
	GenoStruTrait() : m_genoStruIdx(TraitMaxIndex)
	{
	}


	/// CPPONLY set genotypic structure
	void setGenoStructure(UINT ploidy, const vectoru & loci, bool sexChrom, bool haplodiploid,
		const vectorf & lociPos, const vectorstr & chromNames, const vectorstr & alleleNames,
		const vectorstr & lociNames, UINT maxAllele, const vectorstr & infoFields);

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


	/// distance between loci \c loc1 and \c loc2. These two loci should be
	/// on the same chromosome. The distance will be negative if \c loc1 is after
	/// \c loc2.
	double lociDist(UINT loc1, UINT loc2) const;

	/// return the number of loci left on that chromosome, including locus \c loc
	UINT lociLeft(UINT loc) const;

	/// distance left to the right of the loc, till the end of chromosome
	double distLeft(UINT loc) const;

	/// starting from \c loc, how many markers are covered by distance \c dist (>=0)
	/// the result will be at least 1, even if dist = 0.
	UINT lociCovered(UINT loc, double dist) const;


	/// CPPONLY merge two genotype structure
	GenoStructure & mergeGenoStru(size_t idx, bool byChromosome) const;

	/// CPPONLY
	GenoStructure & removeLociFromGenoStru(const vectoru & remove = vectoru(), const vectoru & keep = vectoru());

	/// CPPONLY add some loci to genotype structure
	GenoStructure & insertBeforeLociToGenoStru(const vectoru & idx, const vectorf & pos, const vectorstr & names) const;

	/// CPPONLY append some loci to genotype structure
	GenoStructure & insertAfterLociToGenoStru(const vectoru & idx, const vectorf & pos, const vectorstr & names) const;

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


	/// return ploidy, the number of homologous sets of chromosomes
	UINT ploidy() const
	{

		DBG_FAILIF(m_genoStruIdx == TraitMaxIndex, SystemError,
			"Ploidy: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_ploidy;
	}


	/// return ploidy name, \c haploid, \c diploid, or \c triploid etc.
	string ploidyName() const;

	/// return the number of loci on chromosome \c chrom, equivalent to <tt> numLoci()[chrom] </tt>
	UINT numLoci(UINT chrom) const
	{
		DBG_FAILIF(m_genoStruIdx == TraitMaxIndex, SystemError,
			"numLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGECHROM(chrom);
		return s_genoStruRepository[m_genoStruIdx].m_numLoci[chrom];
	}


	/// return the number of loci on all chromosomes
	vectoru numLoci() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_numLoci;
	}


	/// determine whether or not the last chromosome is sex chromosome
	bool sexChrom() const
	{
		DBG_FAILIF(m_genoStruIdx == TraitMaxIndex, SystemError,
			"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_sexChrom;
	}


	bool haplodiploid() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_haplodiploid;
	}


	/// return the total number of loci on all chromosomes
	UINT totNumLoci() const
	{

		DBG_FAILIF(m_genoStruIdx == TraitMaxIndex, SystemError,
			"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_totNumLoci;
	}


	/// return the total number of loci times ploidy
	UINT genoSize() const
	{
		DBG_FAILIF(m_genoStruIdx == TraitMaxIndex, SystemError,
			"totNumLoci: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_genoSize;
	}


	/// return the position of a locus
	double locusPos(UINT locus) const
	{
		DBG_FAILIF(m_genoStruIdx == TraitMaxIndex, SystemError,
			"locusPos: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGEABSLOCUS(locus);
		return s_genoStruRepository[m_genoStruIdx].m_lociPos[locus];
	}


	/// return loci positions
	vectorf lociPos() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_lociPos;
	}


	/// return a \c carray of loci positions of all loci
	/**
	 \note Modifying loci position directly using this function is strongly discouraged.
	 */
	PyObject * arrLociPos()
	{
		return Double_Vec_As_NumArray(s_genoStruRepository[m_genoStruIdx].m_lociPos.begin(),
			s_genoStruRepository[m_genoStruIdx].m_lociPos.end() );
	}


	/// return a \c carray of loci positions on a given chromosome
	/**
	 \note Modifying loci position directly using this function is strongly discouraged.
	 */
	PyObject * arrLociPos(UINT chrom)
	{
		CHECKRANGECHROM(chrom);

		return Double_Vec_As_NumArray(
			s_genoStruRepository[m_genoStruIdx].m_lociPos.begin() + chromBegin(chrom),
			s_genoStruRepository[m_genoStruIdx].m_lociPos.begin() + chromEnd(chrom) );
	}


	/// return the number of chromosomes
	UINT numChrom() const
	{
		DBG_FAILIF(m_genoStruIdx == TraitMaxIndex, SystemError,
			"numChrom: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		return s_genoStruRepository[m_genoStruIdx].m_numChrom;
	}


	/// CPPONLY return an array of chromosome indexes
	const vectoru & chromIndex() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_chromIndex;
	}


	/// return the index of the first locus on a chromosome
	UINT chromBegin(UINT chrom) const
	{
		DBG_FAILIF(m_genoStruIdx == TraitMaxIndex, SystemError,
			"chromBegin: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGECHROM(chrom);

		return s_genoStruRepository[m_genoStruIdx].m_chromIndex[chrom];
	}


	/// return the index of the last locus on a chromosome plus 1
	UINT chromEnd(UINT chrom) const
	{
		DBG_FAILIF(m_genoStruIdx == TraitMaxIndex, SystemError,
			"chromEnd: You have not set genoStructure. Please use setGenoStrucutre to set such info.");

		CHECKRANGECHROM(chrom);

		return s_genoStruRepository[m_genoStruIdx].m_chromIndex[chrom + 1];
	}


	/// return the absolute index of a locus on a chromosome. c.f. \c chromLocusPair
	UINT absLocusIndex(UINT chrom, UINT locus)
	{
		CHECKRANGECHROM(chrom);
		CHECKRANGELOCUS(chrom, locus);

		return s_genoStruRepository[m_genoStruIdx].m_chromIndex[chrom] + locus ;
	}


	/// return a <tt>(chrom, locus)</tt> pair of an absolute locus index, c.f. \c absLocusIndex
	std::pair<UINT, UINT> chromLocusPair(UINT locus) const;

	/// return the name of an chrom
	string chromName(const UINT chrom) const
	{
		DBG_FAILIF(chrom >= s_genoStruRepository[m_genoStruIdx].m_numChrom, IndexError,
			"Chromosome index " + toStr(chrom) + " out of range of 0 ~ " +
			toStr(s_genoStruRepository[m_genoStruIdx].m_numChrom));

		return s_genoStruRepository[m_genoStruIdx].m_chromNames[chrom];
	}


	/// return an array of chrom names
	vectorstr chromNames() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_chromNames;
	}


	/// return the index of a chromosome by its name
	UINT chromByName(const string name) const
	{
		const vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_chromNames;
		vectorstr::const_iterator it = std::find(names.begin(), names.end(), name);

		if (it == names.end())
			throw ValueError("Failed to find chrom with name " + name);
		return it - names.begin();
	}


	/// return the name of an allele (if previously specified). Default to allele index.
	string alleleName(const UINT allele) const;

	/// return an array of allele names
	vectorstr alleleNames() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_alleleNames;
	}


	/// return the name of a locus
	string locusName(const UINT loc) const
	{
		DBG_FAILIF(loc >= s_genoStruRepository[m_genoStruIdx].m_totNumLoci, IndexError,
			"Locus index " + toStr(loc) + " out of range of 0 ~ " +
			toStr(s_genoStruRepository[m_genoStruIdx].m_totNumLoci));

		return s_genoStruRepository[m_genoStruIdx].m_lociNames[loc];
	}


	/// return names of all loci
	vectorstr lociNames() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_lociNames;
	}


	/// return the index of a locus by its locus name
	UINT locusByName(const string name) const
	{
		const vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_lociNames;
		vectorstr::const_iterator it = std::find(names.begin(), names.end(), name);

		if (it == names.end())
			throw ValueError("Failed to find locus with name " + name);
		return it - names.begin();
	}


	/// return an array of locus indexes by locus names
	vectoru lociByNames(const vectorstr & names) const
	{
		vectoru indexes;

		for (vectorstr::const_iterator name = names.begin(); name != names.end(); ++name)
			indexes.push_back(locusByName(*name));
		return indexes;
	}


	/// return the maximum allele value for all loci. Default to maximum allowed allele state.
	/**
	   Maximum allele value has to be \c 1 for binary modules. \c maxAllele is
	   the maximum possible allele value, which allows <tt>maxAllele+1</tt> alleles
	   <tt>0, 1, ..., maxAllele</tt>.
	 \sa setMaxAllele
	 */
	UINT maxAllele() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_maxAllele;
	}


	/// CPPONLY set the maximum allele value for all loci
	void setMaxAllele(UINT maxAllele)
	{
#ifdef BINARYALLELE
		DBG_ASSERT(maxAllele == 1,  ValueError,
			"max allele must be 1 for binary modules");
#else
		s_genoStruRepository[m_genoStruIdx].m_maxAllele = maxAllele;
#endif
	}


	/// determine if an information field exists
	bool hasInfoField(const string & name) const
	{
		vectorstr & names = s_genoStruRepository[m_genoStruIdx].m_infoFields;

		return std::find(names.begin(), names.end(), name) != names.end();
	}


	/// obtain the number of information fields
	UINT infoSize() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_infoFields.size();
	}


	/// return an array of all information fields
	vectorstr infoFields() const
	{
		return s_genoStruRepository[m_genoStruIdx].m_infoFields;
	}


	/// obtain the name of information field \c idx
	string infoField(UINT idx) const
	{
		CHECKRANGEINFO(idx);
		return s_genoStruRepository[m_genoStruIdx].m_infoFields[idx];
	}


	/// return the index of the field \c name, return \c -1 if not found
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