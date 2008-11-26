/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu
*                                                                         *
*   $LastChangedDate$
*   $Rev$                                                       *
*                                                                         *
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

#ifndef _PEDIGREE_H
#define _PEDIGREE_H
/**
 \file
 \brief head file of class pedigree
 */
#include "population.h"

enum RelativeType {
	REL_None,           // do nothing
	REL_Self,           // individual himself or herself.
	REL_Offspring,      // All offspring with all spouses (if there are more than one spouse)
	REL_Spouse,         // All spouses (with at least one offspring)
	REL_FullSibling,    // Siblings who share two parents
	REL_Sibling,        // Siblings who share at least one parent
};

enum SexChoice {
	AnySex = 0,
	MaleOnly = 1,
	FemaleOnly = 2,
	OppositeSex = 3
};

namespace simuPOP {

class pedigree : public population
{
public:
	/** Create a pedigree object from a population, using a subset of loci,
	 *  information fields and ancestral generations.
	 */
	pedigree(const population & pop, const vectoru & loci = vectoru(),
		const vectorstr & infoFields = vectorstr(), int ancGen = -1);

	/** This function locates relatives (of type \e relType, and sex \e relSex) 
	 *  of each individual and store their indexes in specified information
	 *  fields \e relFields. The indexes of parents in the parental generation
	 *  should be available in information fields \e parentFields (default to
	 *  <tt>['father_idx', 'mother_idx']</tt> which are the information fields
	 *  used by operator \c parentsTagger. This function currently only work
	 *  for diploid populations. \n
	 *
	 *  \param relType Relative type, which can be
	 *  \li REL_Self set indexes of individual themselves.
	 *  \li REL_Spouse locate spouses of individuals in the current generation.
	 *    A spouse is defined as two individuals
	   			having an offspring with shared \c parentFields. If more than one \c infoFields is given,
	   			multiple spouses can be identified.
	 \li REL_Offspring index of offspring in the offspring generation. If only one
	   			parent is given, only paternal or maternal relationship is considered. For example,
	   			<tt>parentFields=['father_idx']</tt> will locate offspring for all fathers.
	 \li REL_FullSibling all siblings with the same parents
	 \li REL_Sibling all sibs with at least one shared parent
	 * \param relFields information fields to hold relatives. The number of these fields
	 *		limits the number of relatives to locate.
	 * \param gen Find relatives for individuals for how many generations. Default to -1,
	 *      meaning for all generations. If a non-negative number is given, up till generation
	 *      gen will be processed.
	 * \param relSex Whether or not only locate relative or certain sex. It can be
	 *		AnySex (do not care, default), MaleOnly, FemaleOnly, or OppositeSex (only locate
	 *      relatives of opposite sex.
	 * \return if relatives are successfully located. Possible problems are
	 *      father and mother indexes are not available, or insufficient parental generations.
	 *
	 * <group>6-ancestral</group>
	 */
	void locateRelatives(RelativeType relType, const vectorstr & relFields,
		int gen = -1, SexChoice relSex = AnySex,
		const vectorstr & parentFields = vectorstr());

	/// Trace a relative path in a population and record the result in the given information fields.
	/**
	 \param pathGen A list of generations that form a relative path. This array is one element longer
	   	than \c pathFields, with gen_i, gen_i+1 indicating the current and destinating generation
	   	of information fields path_i.
	 \param pathFields A list of list of information fields forming a path to trace a certain
	   	type of relative.
	 \param resultFields Where to store located relatives. Note that the result will be saved
	   	in the starting generation specified in \c pathGen[0], which is usually 0.
	 \param pathSex (Optional) A list of sex choices, AnySex, Male, Female or OppositeSex,
	   	that is used to choose individuals at each step. Default to AnySex.

	   For example,
	   <tt>
	   	setInfoWithRelatives(pathGen = [0, 1, 1, 0],
	   		pathFields = [['father_idx', 'mother_idx'], ['sib1', 'sib2'],
	   			['off1', 'off2']],
	   		pathSex = [AnySex, MaleOnly, FemaleOnly],
	   		resultFields = ['cousin1', 'cousin2'])
	   </tt>
	   This function will
	   1. locate father_idx and mother_idx for each individual at generation 0 (\c pathGen[0])
	   2. find AnySex individuals referred by father_idx and mother_idx at generation 1 (\c pathGen[1])
	   3. find informaton fields \c sib1 and \c sib2 from these parents
	   4. locate MaleOnly individuals referred by \c sib1 and \c sib2 from generation 1 (\c pathGen[2])
	   5. find information fields \c off1 and \c off2 from these individuals, and
	   6. locate FemaleOnly indiviudals referred by \c off1 and \of2 from geneartion 0 (\c pathGen[3])
	   7. Save index of these individuals to information fields \c cousin1 and \c cousin2 at
	   	genearation \c pathGen[0].

	   In short, this function locates father or mother's brother's daughters.
	 * <group>6-ancestral</group>
	 */
	bool setIndexesOfRelatives(const vectoru & pathGen,
		const stringMatrix & pathFields,
		const vectori & pathSex = vectori(),
		const vectorstr & resultFields = vectorstr());
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
// 	/// The reason why I do not want to use a
// 	/// vector<vector< struct<dad, mom, vector<info>>>
// 	/// structure is because mom and info are optional.
// 	/// Because pedigree can be large, it makes sense to treat
// 	/// them separately.
// 
// 	typedef vector<vectorlu> Pedigree;
// 	// pedigree subpopulation size information ...
// 	typedef Pedigree PedSize;
// 	typedef vector<vector<vectorf> > PedInfo;
// 
// public:
// 	/// create a pedigree. If a filename \c pedfile is given, the pedgree
// 	/// will be loaded from this file.
// 	pedigree(int numParents = 2, const string & pedfile = string());
// 
// 	/// population size at generation \c gen
// 	ULONG popSize(ULONG gen)
// 	{
// 		DBG_FAILIF(gen >= m_paternal.size(), IndexError,
// 			"Generation number " + toStr(gen) + " out of bound (<"
// 			+ toStr(m_paternal.size()) + ")");
// 		return m_paternal[gen].size();
// 	}
// 
// 
// 	/// return the number of parents for each individual
// 	UINT numParents()
// 	{
// 		return m_numParents;
// 	}
// 
// 
// 	/// Return the number of generations of this pedigree
// 	ULONG gen()
// 	{
// 		return m_paternal.size();
// 	}
// 
// 
// 	/// Return the index of the father of individual \c idx at generation \c gen
// 	/// The returned index is the absolute index of father in the parental generation
// 	ULONG father(ULONG gen, ULONG idx);
// 
// 	/// Return the index of the mother of individual \c idx at generation \c gen
// 	/// The returned index is the absolute index of mother in the parental generation
// 	ULONG mother(ULONG gen, ULONG idx);
// 
// 	/// Return the index of the father of individual \c idx of subpopulation
// 	/// \c subPop at generation \c gen
// 	/// The returned index is the absolute index of father in the parental generation
// 	ULONG father(ULONG gen, SubPopID subPop, ULONG idx);
// 
// 	/// Return the index of the mother of individual \c idx of subpopulation
// 	/// \c subPop at generation \c gen
// 	/// The returned index is the absolute index of mother in the parental generation
// 	ULONG mother(ULONG gen, SubPopID subPop, ULONG idx);
// 
// 	/// Set the index of the father of individual  \c idx at generation \c gen
// 	void setFather(ULONG parent, ULONG gen, ULONG idx);
// 
// 	/// Set the index of the mother of individual  \c idx at generation \c gen
// 	void setMother(ULONG parent, ULONG gen, ULONG idx);
// 
// 	/// Set the index of the father of individual \c idx of subpopulation \c subPop at generation \c gen
// 	void setFather(ULONG parent, ULONG gen, SubPopID subPop, ULONG idx);
// 
// 	/// Set the index of the mother of individual \c idx of subpopulation \c subPop at generation \c gen
// 	void setMother(ULONG parent, ULONG gen, SubPopID subPop, ULONG idx);
// 
// 	/// Return information \c name of individual \c idx at generation \c gen
// 	double info(ULONG gen, ULONG idx, const string & name);
// 
// 	/// Return information \c name of all individuals at generation \c gen
// 	vectorf info(ULONG gen, const string & name);
// 
// 	/// Return information \c name of individual \c idx of subpopulation \c subPop at generation \c gen
// 	double info(ULONG gen, SubPopID subPop, ULONG idx, const string & name);
// 
// 	/// Set information \c name  of individual \c idx at generation \c gen
// 	void setInfo(double info, ULONG gen, ULONG idx, const string & name);
// 
// 	/// Set information \c name of individual \c idx of subpopulation \c subPop at generation \c gen
// 	void setInfo(double info, ULONG gen, SubPopID subPop, ULONG idx, const string & name);
// 
// 	/// Add a generation to the existing pedigree, with given subpopulation sizes
// 	/// \c subPopSize . All parental indexes and information will be set to zero
// 	/// for the new generation.
// 	void addGen(const vectorlu & subPopSize);
// 
// 	/// Return the subpopulation sizes of generation \c gen
// 	vectorlu subPopSizes(ULONG gen);
// 
// 	/// Return the subpopulation size of subpopulation \c subPop of generation \c gen
// 	ULONG subPopSize(ULONG gen, SubPopID subPop);
// 
// 	/// Make a copy of this pedigree
// 	pedigree * clone()
// 	{
// 		return new pedigree();
// 	}
// 
// 
// 	/// load pedigree from a file, the file is usually saved by \c parentTagger or
// 	/// \c parentsTagger. The format is described in the simuPOP reference manual
// 	void load(const string & filename);
// 
// 	/// load information \c name from a information pedigree file and add to this pedigree
// 	/// Information \c name should not have existed in the pedigree. The information
// 	/// pedigree file \c filename is usually produced by taggers such as \c sexTagger
// 	/// \c affectionTagger, \c pyTagger and \c infoTagger.
// 	void loadInfo(const string & filename, const string & name);
// 
// 	/// load information \c names from a information pedigree file and add to this pedigree
// 	/// Information names in \c names should not have existed in the pedigree.
// 	/// pedigree file \c filename is usually produced by taggers such as \c sexTagger
// 	/// \c affectionTagger, \c pyTagger and \c infoTagger.
// 	void loadInfo(const string & filename, const vectorstr & names);
// 
// 	/// add an information field to the pedigree, with given initial value
// 	void addInfo(const string & name, double init = 0);
// 
// 	/// Write the pedigree to a file
// 	void save(const string & filename);
// 
// 	/// Save auxiliary information \c name to an information pedigree file
// 	void saveInfo(const string & filename, const string & name);
// 
// 	/// save auxiliary information \c names to an information pedigree file
// 	void saveInfo(const string & filename, const vectorstr & names);
// 
// 	/// Mark individuals \c inds as 'unrelated' in the last generation
// 	/// \c removeUnrelated function will remove these individuals from the pdeigree
// 	void selectIndividuals(const vectorlu & inds);
// 
// 	/// mark individuals that are unrelated to the visible individuals at the
// 	/// last generation (not marked by selectIndividuals) from the pedigree.
// 	void markUnrelated();
// 
// 	/// remove individuals that are unrelated to the last generation from the pedigree
// 	/// indexes will be adjusted.
// 	/// WARNING: if adjust_index=false, an invalid pedigree will be generated
// 	void removeUnrelated(bool adjust_index = true);
// 
// private:
// 	int m_numParents;
// 
// 	Pedigree m_paternal;
// 	///
// 	Pedigree m_maternal;
// 	///
// 	PedSize m_pedSize;
// 	///
// 	PedInfo m_info;
// 	///
// 	vectorstr m_infoNames;
// };
// 


}
#endif
