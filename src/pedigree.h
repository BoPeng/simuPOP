/**
 *  $File: pedigree.h $
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

#ifndef _PEDIGREE_H
#define _PEDIGREE_H
/**
   \file
   \brief head file of class Pedigree
 */
#include "population.h"


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
	 *  (parameter \e loci, default to no locus), information fields
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
	Pedigree(const Population & pop, const uintList & loci = vectoru(),
		const stringList & infoFields = vectorstr(),
		const uintList & ancGens = uintList(),
		const string & idField = "ind_id", const string & fatherField = "father_id",
		const string & motherField = "mother_id", bool stealPop = false);

	/// CPPONLY copy constructor
	Pedigree(const Pedigree & rhs);

	/** Create a cloned copy of a Pedigree.
	 *  <group>1-ped</group>
	 */
	Pedigree * clone() const;

	/** Return a reference to individual with \e id. An \c IndexError will be
	 *  raised if no individual with \e id is found. An float \e id is
	 *  acceptable as long as it rounds closely to an integer.
	 *  <group>4-ind</group>
	 */
	Individual & indByID(double id);


	/** Return the number of parents each individual has. This function returns
	 *  the number of information fields used to store parental indexes, even
	 *  if one of the fields are unused.
	 *  HIDDEN
	 *  <group>2-info</group>
	 */
	UINT numParents();

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

	/** HIDDEN This function has the potential to change individuals in a
	 *  population so the ID map needs to be rebuilt.
	 */
	void removeIndividuals(const uintList & indexes = vectoru(),
		const floatList & IDs = vectorf(), const string & idField = "ind_id",
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
		UINT chromType = AUTOSOME);

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
	UINT mergeSubPops(const uintList & subPops = uintList(), const string & name = UnnamedSubPop);

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

	mutable std::map<ULONG, Individual *> m_idMap;
};

// /// A pedigree manipulation class
// /**
//    A pedigree.has all the pedigree information that is needed to look at parent
//    offspring relationship in a multi-generation Population.
//
//    Conceptually, there are n generations with the latest generation being
//    generation 0. The number of generations (c.f. gen()) is the number of
//    parental generations plus 1. Therefore, each individual can be identified
//    by (gen, idx).
//
//    Each Individual can have a few properties
//    1. mother (c.f. mother())
//    2. father (c.f. father(), optional because a pedigree can have only one sex)
//    3. subpopulation (if subpopulation structure is given)
//    4. sex (c.f. info('sex'))
//    5. affection (c.f. info('affection')
//    6. arbitrary information fields (c.f. info())
//
//  */
// class pedigree
// {
// private:
//  /// The reason why I do not want to use a
//  /// vector<vector< struct<dad, mom, vector<info>>>
//  /// structure is because mom and info are optional.
//  /// Because pedigree can be large, it makes sense to treat
//  /// them separately.
//
//  typedef vector<vectoru> Pedigree;
//  // pedigree subpopulation size information ...
//  typedef pedigree PedSize;
//  typedef vector<vector<vectorf> > PedInfo;
//
// public:
//  /// create a Pedigree. If a filename \c pedfile is given, the pedgree
//  /// will be loaded from this file.
//  Pedigree(int numParents = 2, const string & pedfile = string());
//
//  /// population size at generation \c gen
//  ULONG popSize(ULONG gen)
//  {
//      DBG_FAILIF(gen >= m_paternal.size(), IndexError,
//          "Generation number " + toStr(gen) + " out of bound (<"
//          + toStr(m_paternal.size()) + ")");
//      return m_paternal[gen].size();
//  }
//
//
//  /// return the number of parents for each individual
//  UINT numParents()
//  {
//      return m_numParents;
//  }
//
//
//  /// Return the number of generations of this pedigree
//  ULONG gen()
//  {
//      return m_paternal.size();
//  }
//
////  /// Add a generation to the existing Pedigree, with given subpopulation sizes
//  /// \c subPopSize . All parental indexes and information will be set to zero
//  /// for the new generation.
//  void addGen(const vectoru & subPopSize);
//
//  /// Return the subpopulation sizes of generation \c gen
//  vectoru subPopSizes(ULONG gen);
//
//  /// Return the subpopulation size of subpopulation \c subPop of generation \c gen
//  ULONG subPopSize(ULONG gen, SubPopID subPop);
//
//  /// Make a copy of this pedigree
//  pedigree * clone()
//  {
//      return new Pedigree();
//  }
//
//
//  /// load pedigree from a file, the file is usually saved by \c parentTagger or
//  /// \c ParentsTagger. The format is described in the simuPOP reference manual
//  void load(const string & filename);
//
//  /// load information \c name from a information pedigree file and add to this pedigree
//  /// Information \c name should not have existed in the Pedigree. The information
//  /// pedigree file \c filename is usually produced by taggers such as \c sexTagger
//  /// \c affectionTagger, \c PyTagger and \c infoTagger.
//  void loadInfo(const string & filename, const string & name);
//
//  /// load information \c names from a information pedigree file and add to this pedigree
//  /// Information names in \c names should not have existed in the Pedigree.
//  /// pedigree file \c filename is usually produced by taggers such as \c sexTagger
//  /// \c affectionTagger, \c PyTagger and \c infoTagger.
//  void loadInfo(const string & filename, const vectorstr & names);
//
//  /// add an information field to the Pedigree, with given initial value
//  void addInfo(const string & name, double init = 0);
//
//  /// Write the pedigree to a file
//  void save(const string & filename);
//
//  /// Save auxiliary information \c name to an information pedigree file
//  void saveInfo(const string & filename, const string & name);
//
//  /// save auxiliary information \c names to an information pedigree file
//  void saveInfo(const string & filename, const vectorstr & names);
//
//  /// Mark individuals \c inds as 'unrelated' in the last generation
//  /// \c removeUnrelated function will remove these individuals from the pdeigree
//  void selectIndividuals(const vectoru & inds);
//
//  /// mark individuals that are unrelated to the visible individuals at the
//  /// last generation (not marked by selectIndividuals) from the Pedigree.
//  void markUnrelated();
//
//  /// remove individuals that are unrelated to the last generation from the pedigree
//  /// indexes will be adjusted.
//  /// WARNING: if adjust_index=false, an invalid pedigree will be generated
//  void removeUnrelated(bool adjust_index = true);
//
// private:
//  int m_numParents;
//
//  pedigree m_paternal;
//  ///
//  pedigree m_maternal;
//  ///
//  PedSize m_pedSize;
//  ///
//  PedInfo m_info;
//  ///
//  vectorstr m_infoNames;
// };
//


}
#endif
