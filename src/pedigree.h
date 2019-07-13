/**
 *  $File: pedigree.h $
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

#ifndef _PEDIGREE_H
#define _PEDIGREE_H
/**
   \file
   \brief head file of class Pedigree
 */
#include "population.h"

#include "boost_pch.hpp"


namespace simuPOP {

/** The pedigree class is derived from the population class. Unlike a
 *  population class that emphasizes on individual properties, the pedigree
 *  class emphasizes on relationship between individuals. An unique ID for
 *  all individuals is needed to create a pedigree object from a population
 *  object. Compared to the \c Population class, a \c Pedigree object is
 *  optimized for access individuals by their IDs, regardless of population
 *  structure and ancestral generations. Note that the implementation of
 *  some algorithms rely on the fact that parental IDs are smaller than their
 *  offspring because individual IDs are assigned sequentially during
 *  evolution. Pedigrees with manually assigned IDs should try to obey such
 *  a rule.
 */
class Pedigree : public Population
{
public:
	/** Create a pedigree object from a population, using a subset of loci
	 *  (parameter \e loci, can be a list of loci indexes, names, or
	 *  \c ALL_AVAIL, default to no locus), information fields
	 *  (parameter \e infoFields, default to no information field besides
	 *  \e idField, \e fatherField and \e motherField), and ancestral
	 *  generations (parameter \e ancGens, default to all ancestral generations).
	 *  By default, information field \c father_id (parameter \e fatherField)
	 *  and \c mother_id (parameter \e motherField) are used to locate parents
	 *  identified by \c ind_id (parameter \e idField), which should store an
	 *  unique ID for all individuals. Multiple individuls with the same ID are
	 *  allowed and will be considered as the same individual, but a warning
	 *  will be given if they actually differ in genotype or information fields.
	 *  Operators \c IdTagger and  \c PedigreeTagger are usually used to assign
	 *  such IDs, although function \c sampling.indexToID could be used to
	 *  assign unique IDs and construct parental IDs from index based
	 *  relationship recorded by operator \c ParentsTagger. A pedigree object
	 *  could be constructed with one or no parent but certain functions such
	 *  as relative tracking will not be available for such pedigrees. In case
	 *  that your are no longer using your population object, you could steal
	 *  the content from the population by setting \e stealPop to \c True.
	 */
	Pedigree(const Population & pop, const lociList & loci = vectoru(),
		const stringList & infoFields = vectorstr(),
		const uintList & ancGens = uintList(),
		const string & idField = "ind_id", const string & fatherField = "father_id",
		const string & motherField = "mother_id", bool stealPop = false);

	/// CPPONLY copy constructor
	Pedigree(const Pedigree & rhs);

	/// CPPONLY
	size_t idIdx() const
	{
		return m_idIdx;
	}


	/// CPPONLY Return the ID of the father of individual id.
	/// return 0 if id is zero or invalid, or father_idx is -1.
	size_t fatherOf(size_t id) const
	{
		if (id == 0 || m_fatherIdx == -1)
			return 0;
		IdMap::iterator it = m_idMap.find(id);
		if (it == m_idMap.end())
			return 0;
		return toID(it->second->info(m_fatherIdx));
	}


	/// CPPONLY Return the ID of the mother of individual id.
	/// return 0 if id is zero or invalid, or mother_idx is -1.
	size_t motherOf(size_t id) const
	{
		if (id == 0 || m_motherIdx == -1)
			return 0;
		IdMap::iterator it = m_idMap.find(id);
		if (it == m_idMap.end())
			return 0;
		return toID(it->second->info(m_motherIdx));
	}


	/** Create a cloned copy of a Pedigree.
	 *  <group>1-ped</group>
	 */
	Pedigree * clone() const;

	/** Save a pedigree to file \e filename. This function goes through all
	 *  individuals of a pedigree and outputs in each line the ID of individual,
	 *  IDs of his or her parents, sex (\c 'M' or \c 'F'), affection status
	 *  (\c 'A' or \c 'U'), values of specified information fields
	 *  \e infoFields and genotypes at specified loci (parameter \c loci, which
	 *  can be a list of loci indexes, names, or \c ALL_AVAIL). Allele numbers,
	 *  instead of their names are outputed. Two columns are used for each
	 *  locus if the population is diploid. This file can be loaded using
	 *  function \c loadPedigree although additional information such as names
	 *  of information fields need to be specified. This format differs from a
	 *  \c .ped file used in some genetic analysis software in that there is
	 *  no family ID and IDs of all individuals have to be unique. Note that
	 *  parental IDs will be set to zero if the parent is not in the pedigree
	 *  object. Therefore, the parents of individuals in the top-most ancestral
	 *  generation will always be zero.
	 *  <group>1-ped</group>
	 */
	void save(const string & filename, const stringList & infoFields = vectorstr(),
		const lociList & loci = vectoru()) const;

	/** Return a reference to individual with \e id. An \c IndexError will be
	 *  raised if no individual with \e id is found. An float \e id is
	 *  acceptable as long as it rounds closely to an integer.
	 *  <group>4-ind</group>
	 */
	Individual & indByID(double id) const;

	/** CPPONLY */
	Individual & indByID(size_t id) const
	{
		IdMap::iterator it = m_idMap.find(id);

		// if still cannot be found, raise an IndexError.
		if (it == m_idMap.end())
			throw IndexError((boost::format("No individual with ID %1% could be found.") % id).str());
		return *it->second;
	}


	/** Return the number of parents each individual has. This function returns
	 *  the number of information fields used to store parental indexes, even
	 *  if one of the fields are unused.
	 *  HIDDEN
	 *  <group>2-info</group>
	 */
	size_t numParents() const;

	/** This function locates relatives (of type \e relType) of each individual
	 *  and store their IDs in information fields \e relFields. The length of
	 *  \e relFields determines how many relatives an individual can have.
	 *
	 *  Parameter \e relType specifies what type of relative to locate, which
	 *  can be
	 *  \li \c SPOUSE locate spouses with whom an individual has at least one
	 *       common offspring.
	 *  \li \c OUTBRED_SPOUSE locate non-slibling spouses, namely spouses with
	 *       no shared parent.
	 *  \li \c OFFSPRING all offspring of each individual.
	 *  \li \c COMMON_OFFSPRING common offspring between each individual and its
	 *       spouse (located by \c SPOUSE or \c OUTBRED_SPOUSE). \e relFields
	 *       should consist of an information field for spouse and \c m-1 fields
	 *       for offspring where \c m is the number of fields.
	 *  \li \c FULLSIBLING siblings with common father and mother,
	 *  \li \c SIBLING siblings with at least one common parent.
	 *
	 *  Optionally, you can specify the sex and affection status of relatives
	 *  you would like to locate, using parameters \e sex and
	 *  \e affectionStatus. \e sex can be \c ANY_SEX (default), \c MALE_ONLY,
	 *  \c FEMALE_ONLY, \c SAME_SEX or \c OPPOSITE_SEX, and \e affectionStatus can
	 *  be \c AFFECTED, \c UNAFFECTED or \c ANY_AFFECTION_STATUS (default). Only
	 *  relatives with specified properties will be located.
	 *
	 *  This function will by default go through all ancestral generations and
	 *  locate relatives for all individuals. This can be changed by setting
	 *  parameter \e ancGens to certain ancestral generations you would like
	 *  to process.
	 *  <group>4-locate</group>
	 */
	void locateRelatives(RelativeType relType, const vectorstr & resultFields = vectorstr(),
		SexChoice sex = ANY_SEX, AffectionStatus affectionStatus = ANY_AFFECTION_STATUS,
		const uintList & ancGens = uintList());

	/** Trace a relative path in a population and record the result in the
	 *  given information fields \e resultFields. This function is used to
	 *  locate more distant relatives based on the relatives located by
	 *  function \c locateRelatives. For example, after siblings and offspring
	 *  of all individuals are located, you can locate mother's sibling's
	 *  offspring using a <em>relative path</em>, and save their indexes
	 *  in each individuals information fields \e resultFields.
	 *
	 *  A <em>relative path</em> consits of a \e fieldPath that specifies
	 *  which information fields to look for at each step, a \e sex
	 *  specifies sex choices at each generation, and a \e affectionStatus
	 *  that specifies affection status at each generation. \e fieldPath
	 *  should be a list of information fields, \e sex and
	 *  \e affectionStatus are optional. If specified, they should be a list of
	 *  \c ANY_SEX, \c MALE_ONLY, \c FEMALE_ONLY, \c SAME_SEX and \c OppsiteSex
	 *  for parameter \e sex, and a list of \c UNAFFECTED, \c AFFECTED
	 *  and \c ANY_AFFECTION_STATUS for parameter \e affectionStatus.
	 *
	 *  For example, if <tt>fieldPath = [['father_id', 'mother_id'],
	 *  ['sib1', 'sib2'], ['off1', 'off2']]</tt>, and <tt>sex = [ANY_SEX,
	 *  MALE_ONLY, FEMALE_ONLY]</tt>, this function will locate \c father_id
	 *  and \c mother_id for each individual, find all individuals referred
	 *  by \c father_id and \c mother_id, find informaton fields \c sib1 and
	 *  \c sib2 from these parents and locate male individuals referred by
	 *  these two information fields. Finally, the information fields \c off1
	 *  and \c off2 from these siblings are located and are used to locate
	 *  their female offspring. The results are father or mother's brother's
	 *  daughters. Their indexes will be saved in each individuals information
	 *  fields \e resultFields. If a list of ancestral generations is given in
	 *  parameter \e ancGens is given, only individuals in these ancestral
	 *  generations will be processed.
	 *  <group>4-locate</group>
	 */
	bool traceRelatives(const stringMatrix & fieldPath,
		const uintList & sex = vectoru(),
		const uintList & affectionStatus = vectoru(),
		const stringList & resultFields = vectorstr(),
		const uintList & ancGens = uintList());

	/** Return a list of IDs of individuals who have non-negative values at
	 *  information fields \e infoFields. Additional requirements could be
	 *  specified by parameters \e sex and \e affectionStatus.
	 *  \e sex can be \c ANY_SEX (default), \c MALE_ONLY, \c FEMALE_ONLY,
	 *  \c SAME_SEX or \c OPPOSITE_SEX, and \e affectionStatus can be
	 *  \c AFFECTED, \c UNAFFECTED or \c ANY_AFFECTION_STATUS (default). This
	 *  function by default check all individuals in all ancestral generations,
	 *  but you could limit the search using parameter \e subPops (a list of
	 *  (virtual) subpopulations) and ancestral generations \e ancGens.
	 *  Relatives fall out of specified subpopulations and ancestral generaions
	 *  will be considered invalid.
	 *  <group>4-locate</group>
	 */
	vectoru individualsWithRelatives(const stringList & infoFields, const uintList & sex = vectoru(),
		const uintList & affectionStatus = vectoru(), const subPopList & subPops = subPopList(),
		const uintList & ancGens = uintList());

	/** This function goes through all individuals in a pedigree and group
	 *  related individuals into families. If an information field \e pedField
	 *  is given, indexes of families will be assigned to this field of each
	 *  family member. The return value is a list of family sizes corresponding
	 *  to families 0, 1, 2, ... etc. If a list of (virtual) subpopulations
	 *  (parameter \e subPops) or ancestral generations are specified
	 *  (parameter \e ancGens), the search will be limited to individuals in
	 *  these subpopulations and generations.
	 *  <group>4-locate</group>
	 */
	vectoru identifyFamilies(const string & pedField = string(),
		const subPopList & subPops = subPopList(),
		const uintList & ancGens = uintList());

	/** If a list of individuals (\e IDs) is given, this function traces
	 *  backward in time and find all ancestors of these individuals. If \e IDs
	 *  is \c ALL_AVAIL, ancestors of all individuals in the present generation
	 *  will be located. If a list of (virtual) subpopulations (\e subPops) or
	 *  ancestral geneartions (\e ancGens) is given, the search will be limited
	 *  to individuals in these subpopulations and generations. This could be
	 *  used to, for example, find all fathers of \e IDs. This function returns
	 *  a list of IDs, which includes valid specified IDs. Invalid IDs will be
	 *  silently ignored. Note that parameters \e subPops and \e ancGens will
	 *  limit starting IDs if \c IDs is set to \c ALL_AVAIL, but specified
	 *  IDs will not be trimmed according to these parameters.
	 *  <group>4-locate</group>
	 */
	vectoru identifyAncestors(const uintList & IDs = uintList(),
		const subPopList & subPops = subPopList(),
		const uintList & ancGens = uintList());

	/** This function traces forward in time and find all offspring of
	 *  individuals specified in parameter \e IDs. If a list of (virtual)
	 *  subpopulations (\e subPops) or ancestral geneartions (\e ancGens) is
	 *  given, the search will be limited to individuals in these
	 *  subpopulations and generations. This could be used to, for example,
	 *  find all male offspring of \e IDs. This function returns a list of IDs,
	 *  which includes valid starting \e IDs. Invalid IDs are silently ignored.
	 *  Note that parameters \e subPops and \e ancGens will limit search result
	 *  but will not be used to trim specified \e IDs.
	 *  <group>4-locate</group>
	 */
	vectoru identifyOffspring(const uintList & IDs = vectoru(),
		const subPopList & subPops = subPopList(),
		const uintList & ancGens = uintList());

	/** HIDDEN This function has the potential to change individuals in a
	 *  population so the ID map needs to be rebuilt.
	 */
	void removeIndividuals(const uintList & indexes = vectoru(),
		const floatList & IDs = floatList(), const string & idField = "ind_id",
		PyObject * filter = NULL);

	/** HIDDEN This function has the potential to change individuals in a
	 *  population so the ID map needs to be rebuilt.
	 */
	void removeSubPops(const subPopList & subPops);

	/** HIDDEN This function has the potential to change individuals in a
	 *  population so the ID map needs to be rebuilt.
	 */
	void push(Population & pop);

	/** HIDDEN This function has the potential to change individuals in a
	 *  population so the ID map needs to be rebuilt.
	 */
	void addChrom(const vectorf & lociPos, const vectorstr & lociNames = vectorstr(),
		const string & chromName = string(), const stringMatrix & alleleNames = stringMatrix(),
		size_t chromType = AUTOSOME);

	/** HIDDEN This function has the potential to change individuals in a
	 *  population so the ID map needs to be rebuilt.
	 */
	void addChromFrom(const Population & pop);

	/** HIDDEN This function has the potential to change individuals in a
	 *  population so the ID map needs to be rebuilt.
	 */
	void addIndFrom(const Population & pop);

	/** HIDDEN This function has the potential to change individuals in a
	 *  population so the ID map needs to be rebuilt.
	 */
	size_t mergeSubPops(const uintList & subPops = uintList(), const string & name = UnnamedSubPop);

	/** HIDDEN This function has the potential to change individuals in a
	 *  population so the ID map needs to be rebuilt.
	 */
	void resize(const uintList & sizes, bool propagate = false);

	/** HIDDEN This function has the potential to change individuals in a
	 *  population so the ID map needs to be rebuilt.
	 */
	void setSubPopByIndInfo(const string & field);

private:
	void buildIDMap();

	bool acceptableSex(Sex mySex, Sex relSex, SexChoice choice);

	bool acceptableAffectionStatus(bool affected, AffectionStatus choice);

	// a list of functions that will be used in locateRelatives.
	// they are called only once. The reason this is separated is because
	// they are too long when putting in one function.

	void locateSelf(SexChoice relSex, AffectionStatus relAffection,
		const vectorstr & relFields, const vectoru & ancGens);

	void locateSpouse(SexChoice relSex, AffectionStatus relAffection,
		const vectorstr & relFields, const vectoru & ancGens, bool excludeOutbred);

	void locateSibling(SexChoice relSex, AffectionStatus relAffection,
		const vectorstr & relFields, const vectoru & ancGens);

	void locateFullSibling(SexChoice relSex, AffectionStatus relAffection,
		const vectorstr & relFields, const vectoru & ancGens);

	void locateOffspring(SexChoice relSex, AffectionStatus relAffection,
		const vectorstr & relFields, const vectoru & ancGens);

	void locateCommonOffspring(SexChoice relSex, AffectionStatus relAffection,
		const vectorstr & relFields, const vectoru & ancGens);

private:
	const string m_idField;
	const string m_fatherField;
	const string m_motherField;

	int m_idIdx;
	int m_fatherIdx;
	int m_motherIdx;

#if TR1_SUPPORT == 0
	typedef std::map<size_t, Individual *> IdMap;
#elif TR1_SUPPORT == 1
	typedef std::unordered_map<size_t, Individual *> IdMap;
#else
	typedef std::tr1::unordered_map<size_t, Individual *> IdMap;
#endif
	mutable IdMap m_idMap;
};


/** Load a pedigree from a file saved by operator \c PedigreeTagger or function
 *  \c Pedigree.save. This file contains the ID of each offspring and their
 *  parent(s) and optionally sex ('M' or 'F'), affection status ('A' or 'U'),
 *  values of information fields and genotype at some loci. IDs of each
 *  individual and their parents are loaded to information fields \e idField,
 *  \e fatherField and \e motherField. Only numeric IDs are allowed, and
 *  individual IDs must be unique across all generations.
 *
 *  Because this file does not contain generation information, generations to
 *  which offspring belong are determined by the parent-offspring relationships.
 *  Individuals without parents are assumed to be in the top-most ancestral
 *  generation. This is the case for individuals in the top-most ancestral
 *  generation if the file is saved by function ``Pedigree.save()``, and for
 *  individuals who only appear as another individual's parent, if the file is
 *  saved by operator ``PedigreeTagger``. The order at which offsprng is
 *  specified is not important because this function essentially creates a
 *  top-most ancestral generation using IDs without parents, and creates the
 *  next generation using offspring of these parents, and so on until all
 *  generations are recreated. That is to say, if you have a mixture of
 *  pedigrees with different generations, they will be lined up from the top
 *  most ancestral generation.
 *
 *  If individual sex is not specified, sex of of parents are determined by
 *  their parental roles (father or mother) but the sex of individuals in
 *  the last generation can not be determined so they will all be males. If
 *  additional information fields are given, their names have to be specified
 *  using parameter \e infoFields. The rest of the columns are assued to be
 *  alleles, arranged \e ploidy consecutive columns for each locus. If
 *  paraemter \e loci is not specified, the number of loci is calculated by
 *  number of columns divided by \e ploidy (default to 2). All loci are assumed
 *  to be on one chromosome unless parameter \e loci is used to specified number
 *  of loci on each chromosome. Additional parameters such as \e ploidy,
 *  \e chromTypes, \e lociPos, \e chromNames, \e alleleNames, \e lociNames
 *  could be used to specified the genotype structured of the loaded pedigree.
 *  Please refer to class \c Population for details about these parameters.
 */
Pedigree loadPedigree(const string & file,
	const string & idField = "ind_id",
	const string & fatherField = "father_id",
	const string & motherField = "mother_id",
	float ploidy = 2,
	const uintList & loci = vectoru(),
	const uintList & chromTypes = vectoru(),
	const floatList & lociPos = floatList(),
	const stringList & chromNames = vectorstr(),
	const stringMatrix & alleleNames = stringMatrix(),
	const stringList & lociNames = vectorstr(),
	const stringList & subPopNames = vectorstr(),
	const stringList & infoFields = vectorstr());

}
#endif
