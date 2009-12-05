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
   \brief head file of class pedigree
 */
#include "population.h"


namespace simuPOP {

/** The pedigree class is derived from the population class. Unlike a
 *  population class that emphasizes on individual properties, the pedigree
 *  class emphasizes on relationship between individuals. An unique ID for
 *  all individuals is needed to create a pedigree object from a population
 *  object.
 *
 *  A pedigree object can be created from a population, or loaded from
 *  a disk file, which is usually saved by an operator during a previous
 *  evolutionary process. Depending on how a pedigree is saved, sex and
 *  affection status information may be missing.
 */
class pedigree : public population
{
public:
	/** Create a pedigree object from a population, using a subset of loci
	 *  (parameter \e loci, default to no loci), information fields
	 *  (parameter \e infoFields, default to no information field except for
	 *  \e idField, \e fatherField and \e motherField), and ancestral
	 *  generations (parameter \e ancGen, default to all ancestral generations).
	 *  By default, information field \c father_id (parameter \e fatherField)
	 *  and \c mother_id (parameter \e motherField) are used to locate parents
	 *  with \c ind_id (parameter \e idField) as an ID field storing an unique
	 *  ID for every individual. Multiple individuls with the same ID are
	 *  allowed and will be considered as the same individual, but a warning
	 *  will be given if they actually differ in genotype or information fields.
	 *  Operators \c idTagger and  \c pedigreeTagger are usually used to assign
	 *  such IDs, although function \c sampling.indexToID could be used to
	 *  assign unique IDs and construct parental IDs from index based
	 *  relationship recorded by operator \c parentsTagger. This pedigree
	 *  object works with one or no parents but certain functions such as
	 *  relative tracking will not be available for such pedigrees.
	 */
	pedigree(const population & pop, const uintList & loci = vectoru(),
		const stringList & infoFields = vectorstr(), int ancGen = -1,
		const string & idField = "", const string & fatherField = "",
		const string & motherField = "");

	/// CPPONLY copy constructor
	pedigree(const pedigree & rhs);

	/** Create a cloned copy of a pedigree.
	 *  <group>1-ped</group>
	 */
	pedigree * clone() const;

	/** Return the ID of the father of individual \e ID.
	 *  <group>3-rel</group>
	 */
	ULONG father(ULONG ID);

	/** Return the ID of the mather of individual \e ID.
	 *  <group>3-rel</group>
	 */
	ULONG mother(ULONG ID);


	/** Return a reference to individual with \e id. An \c IndexError will be
	 *  raised if no individual with \e id is found. Note that a float \e id
	 *  is acceptable as long as it rounds closely to an integer.
	 *  <group>4-ind</group>
	 */
	individual & indByID(double id);


	/** Return the number of parents each individual has. This function returns
	 *  the number of information fields used to store parental indexes, even
	 *  if one of the fields are unused.
	 *  <group>2-info</group>
	 */
	UINT numParents();

	/** This function locates relatives (of type \e relType) of each individual
	 *  and store their ID (if \c idField is specified) or indexes (if
	 *  \c idField is not specified) in information fields \e relFields. The
	 *  length of \e relFields determines how many relatives an individual can
	 *  have.
	 *
	 *  Parameter \e relType specifies what type of relative to locate. It can
	 *  be
	 *  \li \c Spouse locate spouses with whom an individual has at least one
	 *       common offspring.
	 *  \li \c OutbredSpouse locate non-slibling spouses, namely spouses with
	 *       no shared parent.
	 *  \li \c Offspring all offspring of each individual.
	 *  \li \c CommonOffspring common offspring between each individual and its
	 *       spouse (located by Spouse or OutbredSpouse). \e relFields should
	 *       consist of an information field for spouse and \n m-1 fields for
	 *       offspring where \n m is the number of fields.
	 *  \li \c FullSibling siblings with common father and mother,
	 *  \li \c Sibling siblings with at least one common parent.
	 *
	 *  Optionally, you can specify the sex of relatives you would like to
	 *  locate, in the form of <tt>relType=(type, sexChoice)</tt>. \e sexChoice
	 *  can be \c AnySex (default), \c MaleOnly, \c FemaleOnly, \c SameSex or
	 *  \c OppositeSex.
	 *
	 *  This function will by default go through all ancestral generations and
	 *  locate relatives for all individuals. This can be changed by setting
	 *  parameter \e ancGen to the greatest ancestral generation you would like
	 *  to process.
	 *  <group>4-locate</group>
	 */
	void locateRelatives(uintList relType = vectoru(), const vectorstr & relFields = vectorstr(),
		int ancGen = -1);

	/** Trace a relative path in a population and record the result in the
	 *  given information fields \e resultFields. This function is used to
	 *  locate more distant relatives based on the relatives located by
	 *  function \c locateRelatives. For example, after siblings and offspring
	 *  of all individuals are located, you can locate mother's sibling's
	 *  offspring using a <em>relative path</em>, and save their indexes
	 *  in each individuals information fields \e resultFields.
	 *
	 *  A <em>relative path</em> consits of a \e pathFields that specifies
	 *  which information fields to look for at each step, and a \e pathSex
	 *  specifies sex choices at each generation, which should be a list of
	 *  \c AnySex, \c MaleOnly, \c FemaleOnly, \c SameSex and \c OppsiteSex.
	 *  The default value for this paramter is \c AnySex at all steps.
	 *
	 *  For example, if <tt>pathFields = [['father_ID', 'mother_ID'],
	 *  ['sib1', 'sib2'], ['off1', 'off2']]</tt>, and <tt>pathSex = [AnySex,
	 *  MaleOnly, FemaleOnly]</tt>, this function will locate \c father_ID
	 *  and \c mother_ID for each individual, find all individuals referred
	 *  by \c father_ID and \c mother_ID, find informaton fields \c sib1 and
	 *  \c sib2 from these parents and locate male individuals referred by
	 *  these two information fields. Finally, the information fields \c off1
	 *  and \c off2 from these siblings are located and are used to locate
	 *  their female offspring at the present geneartion. The results are
	 *  father or mother's brother's daughters. Their indexes will be saved in
	 *  each individuals information fields \e resultFields. If a non-negative
	 *  \e ancGen is given, only individuals in these ancestral generations
	 *  will be processed.
	 *  <group>4-locate</group>
	 */
	bool traceRelatives(const stringMatrix & pathFields,
		const vectori & pathSex = vectori(),
		const vectorstr & resultFields = vectorstr(), int ancGen = -1);

private:
	// a list of functions that will be used in locateRelatives.
	// they are called only once. The reason this is separated is because
	// they are too long when putting in one function.

	void locateSelf(SexChoice relSex, const vectorstr & relFields, int ancGen);

	void locateSpouse(SexChoice relSex, const vectorstr & relFields, int ancGen, bool excludeOutbred);

	void locateSibling(SexChoice relSex, const vectorstr & relFields, int ancGen);

	void locateFullSibling(SexChoice relSex, const vectorstr & relFields, int ancGen);

	void locateOffspring(SexChoice relSex, const vectorstr & relFields, int ancGen);

	void locateCommonOffspring(SexChoice relSex, const vectorstr & relFields, int ancGen);

private:
	string m_idField;
	string m_fatherField;
	string m_motherField;

	int m_idIdx;
	int m_fatherIdx;
	int m_motherIdx;

	std::map<ULONG, individual *> m_idMap;
};

// /// A pedigree manipulation class
// /**
//    A pedigree has all the pedigree information that is needed to look at parent
//    offspring relationship in a multi-generation population.
//
//    Conceptually, there are n generations with the latest generation being
//    generation 0. The number of generations (c.f. gen()) is the number of
//    parental generations plus 1. Therefore, each individual can be identified
//    by (gen, idx).
//
//    Each individual can have a few properties
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
//  typedef Pedigree PedSize;
//  typedef vector<vector<vectorf> > PedInfo;
//
// public:
//  /// create a pedigree. If a filename \c pedfile is given, the pedgree
//  /// will be loaded from this file.
//  pedigree(int numParents = 2, const string & pedfile = string());
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
////  /// Add a generation to the existing pedigree, with given subpopulation sizes
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
//      return new pedigree();
//  }
//
//
//  /// load pedigree from a file, the file is usually saved by \c parentTagger or
//  /// \c parentsTagger. The format is described in the simuPOP reference manual
//  void load(const string & filename);
//
//  /// load information \c name from a information pedigree file and add to this pedigree
//  /// Information \c name should not have existed in the pedigree. The information
//  /// pedigree file \c filename is usually produced by taggers such as \c sexTagger
//  /// \c affectionTagger, \c pyTagger and \c infoTagger.
//  void loadInfo(const string & filename, const string & name);
//
//  /// load information \c names from a information pedigree file and add to this pedigree
//  /// Information names in \c names should not have existed in the pedigree.
//  /// pedigree file \c filename is usually produced by taggers such as \c sexTagger
//  /// \c affectionTagger, \c pyTagger and \c infoTagger.
//  void loadInfo(const string & filename, const vectorstr & names);
//
//  /// add an information field to the pedigree, with given initial value
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
//  /// last generation (not marked by selectIndividuals) from the pedigree.
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
//  Pedigree m_paternal;
//  ///
//  Pedigree m_maternal;
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
