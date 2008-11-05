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

#ifndef _MATING_H
#define _MATING_H
/**
 \file
 \brief head file of class mating and its subclasses
 */
#include "utility.h"
// for trajectory simulation functions
#include "simuPOP_cfg.h"
#include "misc.h"
#include "pedigree.h"
#include "population.h"
#include "operator.h"

#include <string>
#include <algorithm>
using std::max_element;
using std::string;

#include <stack>
using std::stack;

typedef std::vector<PyObject * >           vectorobj;

namespace simuPOP {

/**
   Offspring generators generate offspring from given parents. Generators differ
   from each other by how and how many offspring is generated at each mating event.

   Parameters \c mode, \c numOffspring, \c maxNumOffspring and \c numOffspringFunc
   are used to specify how many offspring will be produced at each mating event.
 \c mode can be one of
 \li \c MATE_NumOffspring: a fixed number of offspring will be produced
   	at all mating events .
 \li \c MATE_PyNumOffspring: A python function, specified by parameter
 \c numOffspringFunc, is called at each mating event to determine the number
   	of offspring to produce.
 \li \c MATE_GeometricDistribution: a Geometric distribution with parameter \c numOffspring
   	   is used to determine the number of offspring of each family.
 \li \c MATE_PoissonDistribution: a Poisson distribution with parameter \c numOffspring
   	   is used to determine the number of offspring of each family.
 \li \c MATE_BinomialDistribution: a Binomial distribution with parameter \c numOffspring
   	   is used to determine the number of offspring of each family.
 \li \c MATE_UniformDistribution: a Uniform <tt> [a, b] </tt> distribution with parameter
 \c numOffspring (a) and \c maxNumOffspring (b) is used to determine the number of offspring of each family.

   This is the base class of all offspring generators, and should not
   be used directly.
   <applicability>all ploidy</applicability>
 */
class offspringGenerator
{
public:
	// numOffspring: constant, numOffspringFunc: call each time before mating
#define MATE_NumOffspring           1
	// call numOffspringFunc each time during mating.
#define MATE_NumOffspringEachFamily 2 // This name is obsolete
#define MATE_PyNumOffspring         2
	// numOffspring and numOffsrpingsFunc call each time before mating is
	// the p for a geometric distribution
#define MATE_GeometricDistribution   3
#define MATE_PoissonDistribution     4
#define MATE_BinomialDistribution    5
	// uniform between numOffspring and maxNumOffspring
#define MATE_UniformDistribution     6


#define MATE_RandomSex               1
#define MATE_ProbOfMale              2
#define MATE_NumOfMale               3
#define MATE_NumOfFemale             4

public:
	/**
	 \param numOffspring Depending on \mode, this paramter can be the number of offspring to
	   	produce, or a paremter of a random distribution.
	 \param numOffspringFunc a Python function that returns the number of offspring at each
	   	mating event. The setting of this parameter implies \mode=MATE_PyNumOffspring.
	 \param maxNumOffspring used when \c numOffspring is generated from a binomial or random
	   	distribution.
	 \param mode can be one of <tt>MATE_NumOffspring, MATE_PyNumOffspring,
	   	MATE_GeometricDistribution, MATE_PoissonDistribution, MATE_BinomialDistribution,</tt>
	   	or <tt>MATE_UniformDistribution</tt>.
	 \param sexParam parameter that controls the sex distribution among offspring. Its exact
	   	meaning is determined by parameter sexMode. Default to 0.5.
	 \param sexMode can be one of
	 \li MATE_RandomSex  Set sex to Male or Female with probability 0.5. Parameter
	 \c sexParam is ignored in this case. This is the default mode.
	 \li MATE_ProbOfMale Set an offspring to Male with probability \c sexParam (default to 0.5)
	 \li MATE_NumOfMale Set \c sexParam offspring to Male
	 \li MATE_NumOfFemale Set \c sexParam offspring to Female.
	 \note: Parameter \c sexMode and \c sexParam are ignored if sex chromosome is defined.
	   	Offspring sex is defined by genotype in this case.
	 */
	offspringGenerator(double numOffspring, PyObject * numOffspringFunc,
	                   UINT maxNumOffspring, UINT mode,
	                   double sexParam, UINT sexMode);


	virtual ~offspringGenerator()
	{
		if (m_numOffspringFunc != NULL)
			Py_DECREF(m_numOffspringFunc);
	}


	/// CPPONLY
	offspringGenerator(const offspringGenerator & rhs);

	virtual offspringGenerator * clone() const = 0;


	/// generate \c numOff offspring
	/// CPPONLY
	virtual UINT generateOffspring(population & pop, individual * dad, individual * mom,
		RawIndIterator & offBegin,
		RawIndIterator & offEnd,
		vector<baseOperator *> & ops) = 0;

	/// CPPONLY
	bool fixedFamilySize() const
	{
		return m_mode == MATE_NumOffspring;
	}


	/// create an offspring generator, save information from \c pop and \c ops to
	/// speed up the calls to \c generateOffspring
	/// CPPONLY
	virtual void initialize(const population & pop, vector<baseOperator *> const & ops);

	virtual void finalize(const population & pop)
	{
		m_initialized = false;
	}


	/// CPPONLY the number of offspring of a genaration \c gen
	/**
	   This is called whenever a family size is needed.
	   Its actual meaning depending on \c mode.
	 */
	ULONG numOffspring(int gen);

	/// return sex according to m_sexParam and m_sexMode
	/// \param count the index of offspring
	Sex getSex(int count);

	/// CPPONLY
	bool initialized() const
	{
		return m_initialized;
	}


	/// CPPONLY
	void setNumParents(int numParents)
	{
		m_numParents = numParents;
	}


	/// CPPONLY
	int numParents() const
	{
		return m_numParents;
	}


public:
	/// number of offspring each mate
	double m_numOffspring;

	/// number of offspring func
	PyObject * m_numOffspringFunc;

	///
	UINT m_maxNumOffspring;

	/// whether or not call m_numOffspringFunc each time
	UINT m_mode;

	///
	double m_sexParam;

	///
	UINT m_sexMode;

protected:
	/// check if any of the during mating operators will set genotype
	/// Note that such an operator will also set individual sex
	/// if sexChromosome is present.
	bool checkFormOffspringGenotype(vector<baseOperator *> const & ops);

	/// see if who will generate offspring genotype
	bool m_formOffGenotype;

#ifndef OPTIMIZED
	size_t m_genoStruIdx;
#endif

private:
	// number of parents needed
	int m_numParents;

	bool m_initialized;
};

/** clone offspring generator copies parental geneotype to a number
   of offspring. Only one parent is accepted. The number of offspring
   produced is controled by parameters \c numOffspring, \c numOffspringFunc,
 \c maxNumOffspring and \c mode. Parameters \c sexParam and \c sexMode is
   ignored.
   <applicability>all ploidy</applicability>
 */
class cloneOffspringGenerator : public offspringGenerator
{
public:
	/**
	 \param sexParam ignored because sex is copied from the parent.
	 \param sexMode ignored because sex is copied from the parent.
	 */
	cloneOffspringGenerator(double numOffspring = 1,
	                        PyObject * numOffspringFunc = NULL,
	                        UINT maxNumOffspring = 1,
	                        UINT mode = MATE_NumOffspring,
	                        double sexParam = 0.5,
	                        UINT sexMode = MATE_RandomSex
	                        ) :
		offspringGenerator(numOffspring, numOffspringFunc, maxNumOffspring,
		                   mode, sexParam, sexMode)
	{
		setNumParents(1);
	}


	offspringGenerator * clone() const
	{
		return new cloneOffspringGenerator(*this);
	}


	/// CPPONLY
	UINT generateOffspring(population & pop, individual * dad, individual * mom,
		RawIndIterator & offBegin,
		RawIndIterator & offEnd,
		vector<baseOperator *> & ops);

};

/** Mendelian offspring generator accepts two parents and pass their
   genotype to a number of offspring following Mendelian's law. Basically,
   one of the paternal chromosomes is chosen randomly to form the paternal
   copy of the offspring, and one of the maternal chromosome is chosen
   randomly to form the maternal copy of the offspring. The number of offspring
   produced is controled by parameters \c numOffspring, \c numOffspringFunc,
 \c maxNumOffspring and \c mode. Recombination will not happen unless
   a during-mating operator recombinator is used.

   <applicability>diploid only</applicability>
 */
class mendelianOffspringGenerator : public offspringGenerator
{
public:
	mendelianOffspringGenerator(double numOffspring = 1,
	                            PyObject * numOffspringFunc = NULL,
	                            UINT maxNumOffspring = 1,
	                            UINT mode = MATE_NumOffspring,
	                            double sexParam = 0.5,
	                            UINT sexMode = MATE_RandomSex
	                            ) :
		offspringGenerator(numOffspring, numOffspringFunc, maxNumOffspring,
		                   mode, sexParam, sexMode),
		m_bt(rng())
	{
		setNumParents(2);
	}


	offspringGenerator * clone() const
	{
		return new mendelianOffspringGenerator(*this);
	}


	/// CPPONLY
	virtual void initialize(const population & pop, vector<baseOperator *> const & ops);

	// the default method to produce offspring
	/// CPPONLY
	/// \param count index of offspring, used to set offspring sex
	/// does not set sex if count == -1.
	void formOffspringGenotype(individual * parent,
		RawIndIterator & it, int ploidy, int count);

	/// CPPONLY
	UINT generateOffspring(population & pop, individual * dad, individual * mom,
		RawIndIterator & offBegin,
		RawIndIterator & offEnd,
		vector<baseOperator *> & ops);

protected:
	// use bernullitrisls with p=0.5 for free recombination
	BernulliTrials m_bt;

	/// sex chromosome handling
	bool m_hasSexChrom;

	// cache chromBegin, chromEnd for better performance.
	vectoru m_chIdx;

};


/** selfing offspring generator works similarly as a mendelian offspring
   generator but a single parent produces both the paternal and maternal
   copy of the offspring chromosomes. This offspring generator accepts a
   dipload parent. A random copy of the parental chromosomes is chosen
   randomly to form the parental copy of the offspring chromosome, and
   is chosen randomly again to form the maternal copy of the offspring
   chromosome.
   <applicability>diploid only</applicability>
 */
class selfingOffspringGenerator : public mendelianOffspringGenerator
{
public:
	selfingOffspringGenerator(double numOffspring = 1,
	                          PyObject * numOffspringFunc = NULL,
	                          UINT maxNumOffspring = 1,
	                          UINT mode = MATE_NumOffspring,
	                          double sexParam = 0.5,
	                          UINT sexMode = MATE_RandomSex
	                          )
		: mendelianOffspringGenerator(numOffspring, numOffspringFunc, maxNumOffspring,
		                              mode, sexParam, sexMode)
	{
		setNumParents(1);
	}


	offspringGenerator * clone() const
	{
		return new selfingOffspringGenerator(*this);
	}


	/// CPPONLY
	UINT generateOffspring(population & pop, individual * parent, individual *,
		RawIndIterator & offBegin,
		RawIndIterator & offEnd,
		vector<baseOperator *> & ops);

};


/** haplodiploid offspring generator mimics sex-determination in honey bees.
   Given a female (queen) parent and a male parent, the female is considered
   as diploid with two set of chromosomes, and the male is condiered as haploid.
   Actually, the first set of male chromosomes are used. During mating,
   female produce eggs, subject to potential recombination and gene conversion,
   while male sperm is identical to the parental chromosome.

   Female offspring has two sets of chromosomes, one from mother and one from
   father. Male offspring has one set of chromosomes from his mother.
   <applicability>haplodiploid only</applicability>
 */
class haplodiploidOffspringGenerator : public mendelianOffspringGenerator
{
public:
	haplodiploidOffspringGenerator(double numOffspring = 1,
	                               PyObject * numOffspringFunc = NULL,
	                               UINT maxNumOffspring = 1,
	                               UINT mode = MATE_NumOffspring,
	                               double sexParam = 0.5,
	                               UINT sexMode = MATE_RandomSex
	                               )
		: mendelianOffspringGenerator(numOffspring, numOffspringFunc, maxNumOffspring,
		                              mode, sexParam, sexMode)
	{
		setNumParents(2);
	}


	void copyParentalGenotype(individual * parent,
		RawIndIterator & it, int ploidy, int count);

	offspringGenerator * clone() const
	{
		return new haplodiploidOffspringGenerator(*this);
	}


	/// CPPONLY
	UINT generateOffspring(population & pop, individual * dad, individual * mom,
		RawIndIterator & offBegin,
		RawIndIterator & offEnd,
		vector<baseOperator *> & ops);

};


/** Parent choosers repeatedly choose parent(s) from a parental
   population, and pass them to offspring generators. A parent
   chooser can select one or two parents, which should match what is
   required by the offspring generator used.

   This is the base class of all parent choosers, and should not
   be used directly.
 */
class parentChooser
{
public:
	typedef std::pair<individual *, individual *> individualPair;

public:
	// CPPONLY
	// numParents can be 0 (undetermined, can be 1 or 2)
	// 1 (one parent), or 2 (two parents)
	parentChooser(int numParents) : m_initialized(false)
	{
		m_numParents = numParents;
	}


	virtual parentChooser * clone() const
	{
		return new parentChooser(*this);
	}


	/// CPPONLY
	virtual void initialize(population & pop, SubPopID subPop) { }

	virtual void finalize(population & pop, SubPopID subPop)
	{
		m_initialized = false;
	}


	/// CPPONLY
	bool initialized() const
	{
		return m_initialized;
	}


	/// CPPONLY
	int numParents()
	{
		return m_numParents;
	}


	/// CPPONLY
	virtual individual * chooseParent(RawIndIterator basePtr)
	{
		return NULL;
	}


	/// CPPONLY
	virtual individualPair chooseParents(RawIndIterator basePtr)
	{
		return individualPair(NULL, NULL);
	}


	virtual ~parentChooser() { }

protected:
	int m_numParents;
	bool m_initialized;
};


/** This parent chooser chooses a parent linearly, regardless of sex
   or fitness values (selection is not considered).
   <applicability>all ploidy</applicability>
 */
class sequentialParentChooser : public parentChooser
{
public:
	sequentialParentChooser() : parentChooser(1)
	{
	}


	parentChooser * clone() const
	{
		return new sequentialParentChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	/// CPPONLY
	individual * chooseParent(RawIndIterator basePtr);

private:
	bool m_selection;
	/// starting individual
	IndIterator m_begin;
	/// ending individual
	IndIterator m_end;
	/// current individual
	IndIterator m_ind;
};


/** This parents chooser chooses two parents sequentially. The
   parents are chosen from their respective sex groups. Selection
   is not considered.
   <applicability>all ploidy</applicability>
 */
class sequentialParentsChooser : public parentChooser
{
public:
	sequentialParentsChooser() :
		parentChooser(2), m_maleIndex(0), m_femaleIndex(0),
		m_numMale(0), m_numFemale(0),
		m_curMale(0), m_curFemale(0)
	{
	}


	parentChooser * clone() const
	{
		return new sequentialParentsChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	/// CPPONLY
	individualPair chooseParents(RawIndIterator basePtr);

	/// CPPONLY
	ULONG numMale() { return m_numMale; }

	/// CPPONLY
	ULONG numFemale() { return m_numFemale; }

private:
	/// internal index to female/males.
	vector<RawIndIterator> m_maleIndex;
	vector<RawIndIterator> m_femaleIndex;

	ULONG m_numMale;
	ULONG m_numFemale;

	ULONG m_curMale;
	ULONG m_curFemale;
};


/** This parents chooser chooses one or two parents from
   	a given pedigree. It works even when only one parent
   	is needed.
   <applicability>all ploidy</applicability>
 */
// class pedigreeParentsChooser : public parentChooser
// {
// public:
// 	pedigreeParentsChooser(const pedigree & ped) :
// 		parentChooser(0), m_pedigree(ped), m_index(0)
// 	{
// 	}
// 
// 
// 	parentChooser * clone() const
// 	{
// 		return new pedigreeParentsChooser(*this);
// 	}
// 
// 
// 	vectorlu subPopSizes(ULONG gen)
// 	{
// 		return m_pedigree.subPopSizes(gen);
// 	}
// 
// 
// 	/// CPPONLY
// 	void initialize(population & pop, SubPopID sp);
// 
// 	/// CPPONLY
// 	individualPair chooseParents(RawIndIterator basePtr);
// 
// private:
// 	pedigree m_pedigree;
// 	UINT m_gen;
// 	SubPopID m_subPop;
// 	RawIndIterator m_begin;
// 	ULONG m_index;
// };

/** This parent chooser chooses a parent randomly from the
   parental generation. If selection is turned on, parents are
   chosen with probabilities that are proportional to their
   fitness values. Sex is not considered. Parameter \c replacement
   determines if a parent can be chosen multiple times. In case
   that \c replacement=false, paremeter \c replenish=true allows
   restart of the process if all parents are exhausted.
   Note that selection is not allowed when \c replacement=false
   because this poses a particular order on individuals in the
   offspring generation.
   <applicability>all ploidy</applicability>
 */
class randomParentChooser : public parentChooser
{
public:
	/**
	 \param replacement if replacement is false, a parent can not
	   		be chosen more than once.
	 \param replenish if all parent has been chosen, choose from
	   		the whole parental population again.
	 */
	randomParentChooser(bool replacement = true, bool replenish = false) :
		parentChooser(1),
		m_replacement(replacement), m_replenish(replenish),
		m_index(0), m_chosen(0), m_sampler(rng())
	{
	}


	parentChooser * clone() const
	{
		return new randomParentChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	/// CPPONLY
	individual * chooseParent(RawIndIterator basePtr);

protected:
	bool m_replacement;
	bool m_replenish;

	bool m_selection;
	///
	vector<RawIndIterator> m_index;
	vector<RawIndIterator> m_chosen;
	/// accumulative fitness
	Weightedsampler m_sampler;
	/// individuals to choose
	size_t m_size;
};


/** This parent chooser chooses two parents randomly, a male
   and a female, from their respective sex groups randomly.
   If selection is turned on, parents are chosen from their sex
   groups with probabilities that are proportional to their
   fitness values. If parameter \c replacement is false,
   a chosen pair of parents can no longer be selected. This feature
   can be used to simulate monopoly. If \c replenish is true,
   a sex group can be replenished when it is exhausted. Note that
   selection is not allowed in the case of monopoly because this
   poses a particular order on individuals in the offspring generation.
   This parents chooser also allows polygamous mating by reusing
   a parent multiple times when returning parents, and allows
   specification of a few alpha individuals who will be the only
   mating individuals in their sex group.
   <applicability>all ploidy</applicability>
 */
class randomParentsChooser : public parentChooser
{
public:
	/**
	 \param replacement choose with (\c True, default) or without (\c False)
	   		replacement. When choosing without replacement, parents
	   		will be paired and can only mate once.
	 \param replenish if set to true, one or both sex groups will
	   		be replenished if they are exhausted.
	 \param polySex Male (polygyny) or Female (polyandry) parent that
	   		will have \c polyNum sex partners.
	 \param polyNum Number of sex partners.
	 \param alphaSex the sex of the alpha individual, i.e. alpha male
	   	or alpha female who be the only mating individuals in their
	   	sex group.
	 \param alphaNum Number of alpha individuals. If \c infoField is
	   	not given, \c alphaNum random individuals with \c alphaSex
	   	will be chosen. If selection is enabled, individuals with higher
	   	fitness values have higher probability to be selected. There is
	   	by default no alpha individual (\c alphaNum = 0).
	 \param alphaField if an information field is given, individuals
	   	with non-zero values at this information field are alpha individuals.
	   	Note that these individuals must have \c alphaSex.

	   Note: If selection is enabled, it works regularly on on-alpha sex, but
	   	works twice on alpha sex. That is to say, \c alphaNum alpha indiviudals
	   	are chosen selectively, and selected again during mating.
	 */
	randomParentsChooser(bool replacement = true, bool replenish = false,
	                     Sex polySex = Male, UINT polyNum = 1,
	                     Sex alphaSex = Male, UINT alphaNum = 0, string alphaField = string()) :
		parentChooser(2),
		m_replacement(replacement), m_replenish(replenish),
		m_polySex(polySex), m_polyNum(polyNum), m_polyCount(0),
		m_alphaSex(alphaSex), m_alphaNum(alphaNum), m_alphaField(alphaField),
		m_lastParent(NULL), m_maleIndex(0), m_femaleIndex(0),
		m_chosenMale(0), m_chosenFemale(0),
		m_maleFitness(0), m_femaleFitness(0),
		m_malesampler(rng()), m_femalesampler(rng())
	{
		DBG_FAILIF(polyNum < 1, ValueError,
			"Number of sex partners has to be at least one");
	}


	parentChooser * clone() const
	{
		return new randomParentsChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	/// CPPONLY
	individualPair chooseParents(RawIndIterator basePtr);

	/// CPPONLY
	ULONG numMale() { return m_numMale; }
	/// CPPONLY
	ULONG numFemale() { return m_numFemale; }

private:
	bool m_replacement;
	bool m_replenish;
	Sex m_polySex;
	UINT m_polyNum;

	UINT m_polyCount;

	Sex m_alphaSex;
	UINT m_alphaNum;
	string m_alphaField;


	individual * m_lastParent;

	bool m_selection;

	ULONG m_numMale;
	ULONG m_numFemale;

	/// internal index to female/males.
	vector<RawIndIterator> m_maleIndex;
	vector<RawIndIterator> m_femaleIndex;
	vector<RawIndIterator> m_chosenMale;
	vector<RawIndIterator> m_chosenFemale;

	vectorf m_maleFitness;
	vectorf m_femaleFitness;

	// weighted sampler
	Weightedsampler m_malesampler;
	Weightedsampler m_femalesampler;

};


/** This parents chooser choose an individual randomly, but choose
   	his/her spouse from a given set of information fields, which stores
   	indexes of individuals in the same generation. A field will be ignored
   	if its value is negative, or if sex is compatible.

   	Depending on what indexes are stored in these information fields,
   	this parent chooser can be used to implement consanguineous mating
   	where close relatives are located for each individual, or certain
   	non-random mating schemes where each individual can only mate with a small
   	number of pre-determinable individuals.

   	This parent chooser (currently) uses \c randomParentChooser to choose
   	one parent and randomly choose another one from the information fields.
   	Because of potentially non-even distribution of valid information
   	fields, the overall process may not be as random as expected, especially
   	when selection is applied.

   	Note: if there is no valid individual, this parents chooser works like
   	a double parentChooser.
   <applicability>all ploidy</applicability>
 */
class infoParentsChooser : public randomParentChooser
{
public:
	/**
	 \param infoFields information fields that store index of matable
	   		individuals.
	 \param replacement if replacement is false, a parent can not
	   		be chosen more than once.
	 \param replenish if all parent has been chosen, choose from
	   		the whole parental population again.
	 */
	infoParentsChooser(const vectorstr & infoFields = vectorstr(),
	                   bool replacement = true, bool replenish = false) :
		randomParentChooser(replacement, replenish),
		m_infoFields(infoFields), m_degenerate(false)
	{
		m_numParents = 2;
		DBG_FAILIF(m_infoFields.empty(), ValueError,
			"At least one information field should be provided for this infoParentsChooser");
	}


	parentChooser * clone() const
	{
		return new infoParentsChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	/// CPPONLY
	individualPair chooseParents(RawIndIterator basePtr);

private:
	vectorstr m_infoFields;
	vectori m_infoIdx;
	// if there is no valid individual, this mating schemes
	// works like a double parentChooser.
	bool m_degenerate;
};


/** This parents chooser accept a Python generator function that yields
   repeatedly an index (relative to each subpopulation) of a parent, or
   indexes of two parents as a Python list of tuple. The generator function
   is responsible for handling sex or selection if needed.
   <applicability>all ploidy</applicability>
 */
class pyParentsChooser : public parentChooser
{
public:
	/**
	 \param parentsGenerator A Python generator function
	 */
	pyParentsChooser(PyObject * parentsGenerator);

	/// CPPONLY
	pyParentsChooser(const pyParentsChooser & rhs)
		: parentChooser(rhs),
		m_func(rhs.m_func),
		m_generator(rhs.m_generator),
		m_parIterator(rhs.m_parIterator)
	{
		Py_INCREF(m_func);
		if (m_generator != NULL)
			Py_INCREF(m_generator);
		if (m_parIterator != NULL)
			Py_INCREF(m_parIterator);
	}


	parentChooser * clone() const
	{
		return new pyParentsChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	void finalize(population & pop, SubPopID sp)
	{
		m_initialized = false;
	}


	/// destructor
	~pyParentsChooser()
	{
		if (m_func != NULL)
			Py_XDECREF(m_func);
		if (m_generator != NULL)
			Py_XDECREF(m_generator);
		if (m_parIterator != NULL)
			Py_XDECREF(m_parIterator);
	}


	/// CPPONLY
	individualPair chooseParents(RawIndIterator basePtr);

private:
#ifndef OPTIMIZED
	ULONG m_size;
#endif
	IndIterator m_begin;

	PyObject * m_func;
	PyObject * m_generator;
	PyObject * m_parIterator;
};


/// the base class of all mating schemes - a required parameter of \c simulator
/**
   Mating schemes specify how to generate offspring from the current population.
   It must be provided when a simulator is created. Mating can perform the following tasks:
 \li change population/subpopulation sizes;
 \li randomly select parent(s) to generate offspring to populate the offspring generation;
 \li apply \em during-mating operators;
 \li apply selection if applicable.
 */
class mating
{

public:
	/// CPPONLY check if the mating type is compatible with the population structure
	/**
	   possible things to check:
	 \li need certain types of individuals (age, sex etc)
	 \li need resizeable population...
	 */
	virtual bool isCompatible(const population & pop) const
	{
		return true;
	}


	/// create a mating scheme (do not use this base mating scheme, use one of its derived classes instead)
	/**
	   By default, a mating scheme keeps a constant population size, generates
	   one offspring per mating event. These can be changed using certain
	   parameters. \c newSubPopSize, \c newSubPopSizeExpr and \c newSubPopSizeFunc
	   can be used to specify subpopulation sizes of the offspring generation.

	 \param newSubPopSize an array of subpopulations sizes, should have the same
	   	number of subpopulations as the current population
	 \param newSubPopSizeExpr an expression that will be evaluated as an array of new subpopulation sizes
	 \param newSubPopSizeFunc a function that takes parameters \c gen (generation number) and \c oldsize
	   (an array of current population size) and return an array of subpopulation sizes of the next generation.
	   This is usually easier to use than its expression version of this parameter.
	 \param subPop if this parameter is given, the mating scheme
	   	will be applied only to the given (virtual) subpopulation.
	   	This is only used in heteroMating where mating schemes
	   	are passed to.
	 \param weight When subPop is virtual, this is used to detemine
	   	the number of offspring for this mating scheme. Weight can be
	 \li 0 (default) the weight will be proportional to the
	   		current (virtual) subpopulation size. If other virutal
	   		subpopulation has non-zero weight, this virtual subpopulation
	   		will produce no offspring (weight 0).
	 \li any negative number -n: the size will be n*m where m is the size
	   	of the (virtual) subpopulation of the parental generation.
	 \li any positive number n: the size will be determined by
	   		weights from all (virtual) subpopulations.

	 \test src_mating.log Demographic models and control of number of offspring per mating event
	 */
	mating(vectorlu newSubPopSize = vectorlu(),
	       string newSubPopSizeExpr = "",
	       PyObject * newSubPopSizeFunc = NULL,
	       vspID subPop = vspID(),
	       double weight = 0);

	/// CPPONLY
	mating(const mating & rhs);

	/// destructor
	virtual ~mating()
	{
		if (m_subPopSizeFunc != NULL)
			Py_DECREF(m_subPopSizeFunc);
	}


	/// CPPONLY
	SubPopID subPop() const
	{
		return m_subPop.subPop();
	}


	/// CPPONLY
	SubPopID virtualSubPop() const
	{
		return m_subPop.virtualSubPop();
	}


	/// CPPONLY
	double weight() const
	{
		return m_weight;
	}


	/// deep copy of a mating scheme
	virtual mating * clone() const
	{
		return NULL;
	}


	/// used by Python print function to print out the general information of the mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::generic mating scheme>";
	}


	/// CPPONLY
	virtual void preparePopulation(population & pop) { }

	/// CPPONLY
	/// a common submit procedure is defined.
	virtual void submitScratch(population & pop, population & scratch);

	/// CPPONLY
	virtual bool mateSubPop(population & pop, SubPopID subPop,
	                        RawIndIterator offBegin, RawIndIterator offEnd,
	                        vector<baseOperator * > & ops)
	{
		return true;
	}


	/// CPPONLY this is not supposed to be called for a base mating class
	/**
	 \param pop population
	 \param scratch scratch population
	 \param ops during-mating operators
	 \return return false when mating fail
	 */
	virtual bool mate(population & pop, population & scratch, vector<baseOperator * > & ops, bool submit);

public:
	/// CPPONLY dealing with the \c pop/subPop size change, copy of structure etc.
	void prepareScratchPop(population & pop, population & scratch);

protected:
	/// new subpopulation size. mostly used to 'keep' subPopsize
	/// after migration.
	vectorlu m_subPopSize;

	/// expression to evaluate subPopSize
	/// population size can change as a result of this
	/// e.g. "%popSize*1.3" whereas %popSize is predefined
	Expression m_subPopSizeExpr;

	/// the function version of the parameter.
	/// the python function should take one parameter (gen) and
	/// return a vector of subpop size.
	PyObject * m_subPopSizeFunc;

	///
	vspID m_subPop;

	///
	double m_weight;

#ifndef OPTIMIZED

public:
	/// record family sizes
	vectori m_famSize;
#endif
};

/// a mating scheme that does nothing
/**
   In this scheme, there is
 \li no mating. Parent generation will be considered as offspring generation.
 \li no subpopulation change. \em During-mating operators will be applied, but
   	the return values are not checked. I.e., subpopulation size parameters will be ignored
   	although some during-mating operators might be applied.

   Note that because the offspring population is the same as parental population,
   this mating scheme can not be used with other mating schemes in a heterogeneous
   mating scheme. cloneMating is recommended for that purpose.
   <applicability>all ploidy</applicability>
 */
class noMating : public mating
{
public:
	/// creat a scheme with no mating
	/**
	 \note All parameters are ignored!
	 */
	noMating(double numOffspring = 1.0,
	         PyObject * numOffspringFunc = NULL,
	         UINT maxNumOffspring = 0,
	         UINT mode = MATE_NumOffspring,
	         double sexParam = 0.5,
	         UINT sexMode = MATE_RandomSex,
	         vectorlu newSubPopSize = vectorlu(),
	         string newSubPopSizeExpr = "",
	         PyObject * newSubPopSizeFunc = NULL,
	         vspID subPop = vspID(),
	         double weight = 0)
	{
		DBG_FAILIF(subPop.isVirtual(), ValueError,
			"noMating can not be used in virtual subpopulations");
	}


	/// destructor
	~noMating()
	{
	}


	/// deep copy of a scheme with no mating
	/**
	 \sa mating::clone() const
	 */
	virtual mating * clone() const
	{
		return new noMating(*this);
	}


	/// used by Python print function to print out the general information of the scheme with no mating
	virtual string __repr__()
	{
		return "<simuPOP::no mating>";
	}


	/// CPPONLY
	void submitScratch(population & pop, population & scratch)
	{
	}


	/// CPPONLY perform mating, but no mating in this scheme
	/**
	   All individuals will be passed to during-mating operators but
	   no one will die (ignore the during-mating failing signal).
	 */
	virtual bool mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit);

};

/// a clone mating that copy everyone from parental to offspring generation.
/**
   Note that
 \li selection is not considered (fitness is ignored)
 \li sequentialParentMating is used. If offspring (virtual) subpopulation size
   is smaller than parental subpopulation size, not all parents will be cloned.
   If offspring (virtual) subpopulation size is larger, some parents will be
   cloned more than once.
 \li numOffspring interface is respected.
 \li during mating operators are applied.
   <applicability>all ploidy</applicability>
 */
class cloneMating : public mating
{
public:
	/// create a binomial selection mating scheme
	/**
	   Please refer to class \c mating for parameter descriptions.
	 */
	cloneMating(double numOffspring = 1.,
	            PyObject * numOffspringFunc = NULL,
	            UINT maxNumOffspring = 0,
	            UINT mode = MATE_NumOffspring,
	            double sexParam = 0.5,
	            UINT sexMode = MATE_RandomSex,
	            vectorlu newSubPopSize = vectorlu(),
	            string newSubPopSizeExpr = "",
	            PyObject * newSubPopSizeFunc = NULL,
	            vspID subPop = vspID(),
	            double weight = 0)
		: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc, subPop, weight),
		m_offGenerator(numOffspring, numOffspringFunc, maxNumOffspring,
		               mode, sexParam, sexMode)
	{
	}


	/// destructor
	~cloneMating()
	{
	}


	/// deep copy of a binomial selection mating scheme
	/**
	 \sa mating::clone() const
	 */
	virtual mating * clone() const
	{
		return new cloneMating(*this);
	}


	/// used by Python print function to print out the general information of the binomial selection mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::binomial random selection>";
	}


	/// CPPONLY perform sexless random mating
	/**
	 \param pop population
	 \param scratch scratch population
	 \param ops during-mating operators
	 \return return false when mating fails
	 */
	virtual bool mateSubPop(population & pop, SubPopID subPop,
		RawIndIterator offBegin, RawIndIterator offEnd,
		vector<baseOperator * > & ops);

protected:
	cloneOffspringGenerator m_offGenerator;
};


/// a mating scheme that uses binomial selection, regardless of sex
/**
   No sex information is involved (binomial random selection). Offspring is chosen from parental generation
   by random or according to the fitness values.
   In this mating scheme,
 \li \c numOffspring protocol is honored;
 \li population size changes are allowed;
 \li selection is possible;
 \li haploid population is allowed.
   <applicability>all ploidy</applicability>
 */
class binomialSelection : public mating
{
public:
	/// create a binomial selection mating scheme
	/**
	   Please refer to class \c mating for parameter descriptions.
	 */
	binomialSelection(double numOffspring = 1.,
	                  PyObject * numOffspringFunc = NULL,
	                  UINT maxNumOffspring = 0,
	                  UINT mode = MATE_NumOffspring,
	                  double sexParam = 0.5,
	                  UINT sexMode = MATE_RandomSex,
	                  vectorlu newSubPopSize = vectorlu(),
	                  string newSubPopSizeExpr = "",
	                  PyObject * newSubPopSizeFunc = NULL,
					  vspID subPop = vspID(),
	                  double weight = 0)
		: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc, subPop, weight),
		m_offGenerator(numOffspring, numOffspringFunc, maxNumOffspring,
		               mode, sexParam, sexMode)
	{
	}


	/// destructor
	~binomialSelection()
	{
	}


	/// deep copy of a binomial selection mating scheme
	/**
	 \sa mating::clone() const
	 */
	virtual mating * clone() const
	{
		return new binomialSelection(*this);
	}


	/// used by Python print function to print out the general information of the binomial selection mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::binomial random selection>";
	}


	/// CPPONLY perform sexless random mating
	/**
	 \param pop population
	 \param scratch scratch population
	 \param ops during-mating operators
	 \return return false when mating fails
	 */
	virtual bool mateSubPop(population & pop, SubPopID subPop,
		RawIndIterator offBegin, RawIndIterator offEnd,
		vector<baseOperator * > & ops);

protected:
	cloneOffspringGenerator m_offGenerator;
};


/**
   This base class defines a general random mating scheme that
   makes full use of a general random parents chooser, and a Mendelian
   offspring generator. A general random parents chooser allows
   selection without replacement, polygemous parents selection
   (a parent with more than one partners), and the definition of
   several alpha individuals.

   Direct use of this mating scheme is not recommended.
 \c randomMating, \c monogemousMating, \c polygemousMating,
 \c alphaMating are all special cases of this mating scheme.
   They should be used whenever possible.
   <applicability>diploid only</applicability>
 */
class baseRandomMating : public mating
{
public:
	/**
	 \param replacement If set to \c True, a parent can be chosen to mate again.
	   		Default to False.
	 \param replenish In case that \c replacement=True, whether or not replenish
	   		 a sex group when it is exhausted.
	 \param polySex sex of polygamous mating. Male for polygyny, Female for polyandry.
	 \param polyNum Number of sex partners.
	 \param alphaSex the sex of the alpha individual, i.e. alpha male
	   		   or alpha female who be the only mating individuals in their
	   		   sex group.
	 \param alphaNum Number of alpha individuals. If \c infoField is
	   		   not given, \c alphaNum random individuals with \c alphaSex
	   		   will be chosen. If selection is enabled, individuals with higher+               fitness values have higher probability to be selected. There is
	   		   by default no alpha individual (\c alphaNum = 0).
	 \param alphaField if an information field is given, individuals
	   		   with non-zero values at this information field are alpha individuals.
	   		   Note that these individuals must have \c alphaSex.

	 */
	baseRandomMating(bool replacement = true,
	                 bool replenish = false,
	                 Sex polySex = Male,
	                 UINT polyNum = 1,
	                 Sex alphaSex = Male,
	                 UINT alphaNum = 0,
	                 string alphaField = string(),
	                 double numOffspring = 1.,
	                 PyObject * numOffspringFunc = NULL,
	                 UINT maxNumOffspring = 0,
	                 UINT mode = MATE_NumOffspring,
	                 double sexParam = 0.5,
	                 UINT sexMode = MATE_RandomSex,
	                 vectorlu newSubPopSize = vectorlu(),
	                 string newSubPopSizeExpr = "",
	                 PyObject * newSubPopSizeFunc = NULL,
	                 bool contWhenUniSex = true,
					 vspID subPop = vspID(),
	                 double weight = 0)
		: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc, subPop, weight),
		m_offspringGenerator(numOffspring, numOffspringFunc,
		                     maxNumOffspring, mode, sexParam, sexMode),
		m_replacement(replacement), m_replenish(replenish),
		m_polySex(polySex), m_polyNum(polyNum),
		m_alphaSex(alphaSex), m_alphaNum(alphaNum), m_alphaField(alphaField),
		m_contWhenUniSex(contWhenUniSex)
	{
	}


	/// destructor
	~baseRandomMating()
	{
	}


	/// deep copy of a random mating scheme
	virtual mating * clone() const
	{
		return new baseRandomMating(*this);
	}


	/// CPPONLY
	virtual bool isCompatible(const population & pop) const
	{
		// test if individual has sex
		// if not, will yield compile time error.
		pop.indBegin()->sex();

#ifndef OPTIMIZED
		if (pop.ploidy() != 2)
			cout << "Warning: This mating type only works with diploid population." << endl;
#endif

		return true;
	}


	/// used by Python print function to print out the general information of the random mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::sexual random mating>";
	}


	/// CPPONLY perform random mating
	virtual bool mateSubPop(population & pop, SubPopID subPop,
		RawIndIterator offBegin, RawIndIterator offEnd,
		vector<baseOperator * > & ops);

protected:
	mendelianOffspringGenerator m_offspringGenerator;

	/// random parents chooser parameters
	bool m_replacement;
	bool m_replenish;
	Sex m_polySex;
	UINT m_polyNum;
	Sex m_alphaSex;
	UINT m_alphaNum;
	string m_alphaField;

	/// if no other sex exist in a subpopulation,
	/// same sex mating will occur if m_contWhenUniSex is set.
	/// otherwise, an exception will be thrown.
	bool m_contWhenUniSex;
};


/// a mating scheme of basic sexually random mating
/**
   In this scheme, sex information is considered for each individual,
   and ploidy is always 2. Within each subpopulation, males and females
   are randomly chosen. Then randomly get one copy of chromosomes from
   father and mother. If only one sex exists in a subpopulation, a
   parameter (\c contWhenUniSex) can be set to determine the behavior.
   Default to continuing without warning.
   <applicability>diploid only</applicability>
 */
class randomMating : public baseRandomMating
{
public:
	/**
	 \param contWhenUniSex continue when there is only one sex in the population. Default to \c True.
	 \n

	   Please refer to class \c mating for descriptions of other parameters.
	 */
	randomMating(double numOffspring = 1.,
	             PyObject * numOffspringFunc = NULL,
	             UINT maxNumOffspring = 0,
	             UINT mode = MATE_NumOffspring,
	             double sexParam = 0.5,
	             UINT sexMode = MATE_RandomSex,
	             vectorlu newSubPopSize = vectorlu(),
	             PyObject * newSubPopSizeFunc = NULL,
	             string newSubPopSizeExpr = "",
	             bool contWhenUniSex = true,
				 vspID subPop = vspID(),
	             double weight = 0)
		:  baseRandomMating(true, false, Male, 1, Male, 0, string(),
		                    numOffspring, numOffspringFunc, maxNumOffspring, mode,
		                    sexParam, sexMode, newSubPopSize, newSubPopSizeExpr,
		                    newSubPopSizeFunc, contWhenUniSex, subPop, weight)
	{
	}


	/// deep copy of a random mating scheme
	virtual mating * clone() const
	{
		return new randomMating(*this);
	}


	/// used by Python print function to print out the general information of the random mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::sexual random mating>";
	}


};


/// a mating scheme of monogamy
/**
   This mating scheme is identical to random mating except that parents
   are chosen without replacement. Under this mating scheme, offspring share
   the same mother must share the same father. In case that all parental
   pairs are exhausted, parameter \c replenish=True allows for the replenishment
   of one or both sex groups.
   <applicability>diploid only</applicability>
 */
class monogamousMating : public baseRandomMating
{
public:
	/**
	 \c replenish This parameter allows replenishment of one or both
	   	parental sex groups in case that they are are exhausted. Default to False.
	   Please refer to class \c mating for descriptions of other parameters.
	 */
	monogamousMating(bool replenish = false,
	                 double numOffspring = 1.,
	                 PyObject * numOffspringFunc = NULL,
	                 UINT maxNumOffspring = 0,
	                 UINT mode = MATE_NumOffspring,
	                 double sexParam = 0.5,
	                 UINT sexMode = MATE_RandomSex,
	                 vectorlu newSubPopSize = vectorlu(),
	                 PyObject * newSubPopSizeFunc = NULL,
	                 string newSubPopSizeExpr = "",
	                 bool contWhenUniSex = true,
					 vspID subPop = vspID(),
	                 double weight = 0)
		:  baseRandomMating(false, replenish, Male, 1, Male, 0, string(),
		                    numOffspring, numOffspringFunc, maxNumOffspring, mode,
		                    sexParam, sexMode, newSubPopSize, newSubPopSizeExpr,
		                    newSubPopSizeFunc, contWhenUniSex, subPop, weight)
	{
	}


	/// deep copy of a random mating scheme
	virtual mating * clone() const
	{
		return new monogamousMating(*this);
	}


	/// used by Python print function to print out the general information of the random mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::monogamous random mating>";
	}


};


/// a mating scheme of polygymy or polyandry
/**
   This mating scheme is composed of a random parents chooser that allows for
   polygamous mating, and a mendelian offspring generator. In this mating scheme,
   a male (or female) parent will have more than one sex partner (\c numPartner).
   Parents returned from this parents chooser will yield the same male (or female)
   parents, each with varying partners.
   <applicability>diploid only</applicability>
 */
class polygamousMating : public baseRandomMating
{
public:
	/**
	 \param polySex sex of polygamous mating. Male for polygyny, Female for polyandry.
	 \param polyNum Number of sex partners.
	 \param replacement If set to \c True, a parent can be chosen to mate again.
	   		Default to False.
	 \param replenish In case that \c replacement=True, whether or not replenish
	   		 a sex group when it is exhausted.
	   Please refer to class \c mating for descriptions of other parameters.
	 */
	polygamousMating(Sex polySex = Male,
	                 UINT polyNum = 1,
	                 bool replacement = false,
	                 bool replenish = false,
	                 double numOffspring = 1.,
	                 PyObject * numOffspringFunc = NULL,
	                 UINT maxNumOffspring = 0,
	                 UINT mode = MATE_NumOffspring,
	                 double sexParam = 0.5,
	                 UINT sexMode = MATE_RandomSex,
	                 vectorlu newSubPopSize = vectorlu(),
	                 PyObject * newSubPopSizeFunc = NULL,
	                 string newSubPopSizeExpr = "",
	                 bool contWhenUniSex = true,
					 vspID subPop = vspID(),
	                 double weight = 0)
		:  baseRandomMating(replacement, replenish, polySex, polyNum, Male, 0, string(),
		                    numOffspring, numOffspringFunc, maxNumOffspring, mode,
		                    sexParam, sexMode, newSubPopSize, newSubPopSizeExpr,
		                    newSubPopSizeFunc, contWhenUniSex, subPop, weight)
	{
	}


	/// deep copy of a random mating scheme
	virtual mating * clone() const
	{
		return new polygamousMating(*this);
	}


	/// used by Python print function to print out the general information of the random mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::sexual random mating>";
	}


};


/// Only a number of alpha individuals can mate with individuals of opposite sex.
/**
   This mating scheme is composed of an random parents chooser with alpha individuals,
   and a Mendelian offspring generator. That is to say, a certain number of alpha
   individual (male or female) are determined by \c alphaNum or an information field. Then,
   only these alpha individuals are able to mate with random individuals of
   opposite sex.
   <applicability>diploid only</applicability>
 */
class alphaMating : public baseRandomMating
{
public:
	/**
	 \param alphaSex the sex of the alpha individual, i.e. alpha male
	   		   or alpha female who be the only mating individuals in their
	   		   sex group.
	 \param alphaNum Number of alpha individuals. If \c infoField is
	   		   not given, \c alphaNum random individuals with \c alphaSex
	   		   will be chosen. If selection is enabled, individuals with higher+               fitness values have higher probability to be selected. There is
	   		   by default no alpha individual (\c alphaNum = 0).
	 \param alphaField if an information field is given, individuals
	   		   with non-zero values at this information field are alpha individuals.
	   		   Note that these individuals must have \c alphaSex.

	   Please refer to class \c mating for descriptions of other parameters.
	   Note: If selection is enabled, it works regularly on on-alpha sex, but
	   		   works twice on alpha sex. That is to say, \c alphaNum alpha indiviudals
	   		   are chosen selectively, and selected again during mating.

	 */
	alphaMating(Sex alphaSex = Male,
	            UINT alphaNum = 0,
	            string alphaField = string(),
	            double numOffspring = 1.,
	            PyObject * numOffspringFunc = NULL,
	            UINT maxNumOffspring = 0,
	            UINT mode = MATE_NumOffspring,
	            double sexParam = 0.5,
	            UINT sexMode = MATE_RandomSex,
	            vectorlu newSubPopSize = vectorlu(),
	            PyObject * newSubPopSizeFunc = NULL,
	            string newSubPopSizeExpr = "",
				vspID subPop = vspID(),
	            double weight = 0)
		:  baseRandomMating(true, false, Male, 1, alphaSex, alphaNum, alphaField,
		                    numOffspring, numOffspringFunc, maxNumOffspring, mode,
		                    sexParam, sexMode, newSubPopSize, newSubPopSizeExpr,
		                    newSubPopSizeFunc, false, subPop, weight)
	{
	}


	/// deep copy of a random mating scheme
	virtual mating * clone() const
	{
		return new alphaMating(*this);
	}


	/// used by Python print function to print out the general information of the random mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::alpha random mating>";
	}


};


/// haplodiploid mating scheme of many hymemopterans
/**
   This mating scheme is composed of an alphaParentChooser and a
   haplodiploidOffspringGenerator. The alphaParentChooser chooses a single
   Female randomly or from a given information field. This female will
   mate with random males from the colony. The offspring will have one of the
   two copies of chromosomes from the female parent, and the first copy
   of chromosomes from the male parent. Note that if a recombinator
   is used, it should disable recombination of male parent.
   <applicability>haplodiploid only</applicability>
 */
class haplodiploidMating : public mating
{
public:
	/**
	 \param alphaSex sex of the alpha individual. Default to Female.
	 \param alphaNum Number of alpha individual. Default to one.
	 \param alphaField information field that identifies the queen of the colony.
	   	By default, a random female will be chosen.

	   Please refer to class \c mating for descriptions of other parameters.
	 */
	haplodiploidMating(Sex alphaSex = Female,
	                   UINT alphaNum = 1,
	                   string alphaField = string(),
	                   double numOffspring = 1.,
	                   PyObject * numOffspringFunc = NULL,
	                   UINT maxNumOffspring = 0,
	                   UINT mode = MATE_NumOffspring,
	                   double sexParam = 0.5,
	                   UINT sexMode = MATE_RandomSex,
	                   vectorlu newSubPopSize = vectorlu(),
	                   PyObject * newSubPopSizeFunc = NULL,
	                   string newSubPopSizeExpr = "",
					   vspID subPop = vspID(),
	                   double weight = 0)
		: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc, subPop, weight),
		m_offspringGenerator(numOffspring, numOffspringFunc,
		                     maxNumOffspring, mode, sexParam, sexMode),
		m_alphaSex(alphaSex), m_alphaNum(alphaNum), m_alphaField(alphaField)
	{
	}


	/// destructor
	~haplodiploidMating()
	{
	}


	/// deep copy of a random mating scheme
	virtual mating * clone() const
	{
		return new haplodiploidMating(*this);
	}


	/// CPPONLY
	virtual bool isCompatible(const population & pop) const
	{
#ifndef OPTIMIZED
		if (pop.ploidy() != 2)
			cout << "Warning: This mating type only works with diploid population." << endl;
#endif
		return true;
	}


	/// used by Python print function to print out the general information of the random mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::sexual random mating>";
	}


	/// CPPONLY perform random mating
	virtual bool mateSubPop(population & pop, SubPopID subPop,
		RawIndIterator offBegin, RawIndIterator offEnd,
		vector<baseOperator * > & ops);

protected:
	haplodiploidOffspringGenerator m_offspringGenerator;

	Sex m_alphaSex;
	UINT m_alphaNum;
	string m_alphaField;
};


/// a mating scheme that follows a given pedigree
// /**
//    In this scheme, a pedigree is given and the mating scheme will
//    choose parents and produce offspring strictly following the pedigree.
//    Parameters setting number of offspring per mating event, and
//    size of the offspring generations are ignored.
// 
//    To implement this mating scheme in pyMating,
//    1.) a newSubPopSizeFunc should be given to return the exact subpopulation
//    	size, returned from pedigree.subPopSizes(gen).
//    2.) use pedigreeChooser to choose parents
//    3.) use a suitable offspring generator to generate offspring.
// 
//    This pedigreeMating helps you do 1 and 2, and use a mendelianOffspringGenerator
//    as the default offspring generator. You can use another offspring generator
//    by setting the generator parameter. Note that the offspring generator can
//    generate one and only one offspring each time.
//    <applicability>all ploid</applicability>
//  */
// class pedigreeMating : public mating
// {
// public:
// 	/**
// 	   Please refer to class \c mating for descriptions of other parameters.
// 	 */
// 	pedigreeMating(pedigree & ped,
//                    offspringGenerator & generator,
// 	               vectorlu newSubPopSize = vectorlu(),
// 	               PyObject * newSubPopSizeFunc = NULL,
// 	               string newSubPopSizeExpr = "",
//					vspID subPop = vspID(),
// 	               double weight = 0);
// 
// 	/// destructor
// 	~pedigreeMating()
// 	{
// 	}
// 
// 
// 	/// deep copy of a random mating scheme
// 	virtual mating * clone() const
// 	{
// 		return new pedigreeMating(*this);
// 	}
// 
// 
// 	/// CPPONLY
// 	virtual bool isCompatible(const population & pop) const
// 	{
// 		return true;
// 	}
// 
// 
// 	/// used by Python print function to print out the general information of the random mating scheme
// 	virtual string __repr__()
// 	{
// 		return "<simuPOP::pedigree-following mating>";
// 	}
// 
// 
// 	/// CPPONLY perform random mating
// 	bool mateSubPop(population & pop, SubPopID subPop,
// 		RawIndIterator offBegin, RawIndIterator offEnd,
// 		vector<baseOperator * > & ops);
// 
// 	bool mate(population & pop, population & scratch, vector<baseOperator * > & ops, bool submit);
// 
// protected:
// 	pedigree * m_pedigree;
// 	offspringGenerator * m_offspringGenerator;
// 	pedigreeParentsChooser m_pedParentsChooser;
// };

/// a mating scheme of selfing
/**
   In this mating scheme, a parent is choosen randomly, acts
   both as father and mother in the usual random mating. The parent
   is chosen randomly, regardless of sex. If selection is turned on,
   the probability that an individual is chosen is proportional to
   his/her fitness.
   <applicability>diploid only</applicability>
 */
class selfMating : public mating
{
public:
	/// create a self mating scheme
	/**
	 \param contWhenUniSex continue when there is only one sex in the population. Default to \c True.
	 \n

	   Please refer to class \c mating for descriptions of other parameters.
	 */
	selfMating(double numOffspring = 1.,
	           PyObject * numOffspringFunc = NULL,
	           UINT maxNumOffspring = 0,
	           UINT mode = MATE_NumOffspring,
	           double sexParam = 0.5,
	           UINT sexMode = MATE_RandomSex,
	           vectorlu newSubPopSize = vectorlu(),
	           PyObject * newSubPopSizeFunc = NULL,
	           string newSubPopSizeExpr = "",
	           bool contWhenUniSex = true,
	           vspID subPop = vspID(),
	           double weight = 0)
		: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc, subPop, weight),
		m_offspringGenerator(numOffspring, numOffspringFunc,
		                     maxNumOffspring, mode, sexParam, sexMode)
	{
	}


	/// destructor
	~selfMating()
	{
	}


	/// deep copy of a self mating scheme
	virtual mating * clone() const
	{
		return new selfMating(*this);
	}


	/// CPPONLY
	virtual bool isCompatible(const population & pop) const
	{
#ifndef OPTIMIZED
		if (pop.ploidy() != 2)
			cout << "Warning: This mating type only works with diploid population." << endl;
#endif
		return true;
	}


	/// used by Python print function to print out the general information of the self mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::self mating>";
	}


	/// CPPONLY perform self mating
	virtual bool mateSubPop(population & pop, SubPopID subPop,
		RawIndIterator offBegin, RawIndIterator offEnd,
		vector<baseOperator *> & ops);

protected:
	selfingOffspringGenerator m_offspringGenerator;
};

/// a mating scheme of consanguineous mating
/**
   In this mating scheme, a parent is choosen randomly and mate with a
   relative that has been located and written to a number of information
   fields.
   <applicability>diploid only</applicability>
 */
class consanguineousMating : public mating
{
public:
	/// create a consanguineous mating scheme
	/**
	   This mating scheme randomly choose a parent and then choose his/her spouse from indexes
	   stored in \c infoFields.

	 \param relativeFields The information fields that stores indexes to other individuals
	   	in a population. If more than one valid (positive value) indexes exist, a random
	   	index will be chosen. (c.f. \c infoParentsChooser ) If there is no individual
	   	having any valid index, the second parent will be chosen randomly from the
	   	whole population.

	 \param func A python function that can be used to prepare the indexes of these
	   	information fields. For example, functions population::locateRelatives and/or
	   	population::setIndexesOfRelatives can be used to locate certain types of relatives
	   	of each individual.

	 \param param An optional parameter that can be passed to \c func.

	   Please refer to \c infoParentsChooser and \c mendelianOffspringGenerator for
	   other parameters.
	 */
	consanguineousMating(
	                     const vectorstr & relativeFields = vectorstr(),
	                     PyObject * func = NULL,
	                     PyObject * param = NULL,
	                     bool replacement = false,
	                     bool replenish = true,
	                     double numOffspring = 1.,
	                     PyObject * numOffspringFunc = NULL,
	                     UINT maxNumOffspring = 0,
	                     UINT mode = MATE_NumOffspring,
	                     double sexParam = 0.5,
	                     UINT sexMode = MATE_RandomSex,
	                     vectorlu newSubPopSize = vectorlu(),
	                     PyObject * newSubPopSizeFunc = NULL,
	                     string newSubPopSizeExpr = "",
	                     bool contWhenUniSex = true,
						 vspID subPop = vspID(),
	                     double weight = 0);

	/// destructor
	~consanguineousMating()
	{
		if (m_func)
			Py_DECREF(m_func);
		if (m_param)
			Py_DECREF(m_param);
	}


	/// CPPONLY
	consanguineousMating(const consanguineousMating & rhs) :
		mating(rhs),
		m_offspringGenerator(rhs.m_offspringGenerator),
		m_relativeFields(rhs.m_relativeFields),
		m_func(rhs.m_func),
		m_param(rhs.m_param),
		m_replacement(rhs.m_replacement),
		m_replenish(rhs.m_replenish)
	{
		if (m_func)
			Py_INCREF(m_func);
		if (m_param)
			Py_INCREF(m_param);
	}


	/// deep copy of a consanguineous mating scheme
	virtual mating * clone() const
	{
		return new consanguineousMating(*this);
	}


	/// CPPONLY
	virtual bool isCompatible(const population & pop) const
	{
#ifndef OPTIMIZED
		if (pop.ploidy() != 2)
			cout << "Warning: This mating type only works with diploid population." << endl;
#endif
		return true;
	}


	/// used by Python print function to print out the general information of the consanguineous mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::consanguineous mating>";
	}


	/// CPPONLY
	virtual void preparePopulation(population & pop);

	/// CPPONLY perform consanguineous mating
	virtual bool mateSubPop(population & pop, SubPopID subPop,
		RawIndIterator offBegin, RawIndIterator offEnd,
		vector<baseOperator *> & ops);

private:
	mendelianOffspringGenerator m_offspringGenerator;
	vectorstr m_relativeFields;
	PyObject * m_func;
	PyObject * m_param;
	bool m_replacement;
	bool m_replenish;
};


/// CPPONLY
void countAlleles(population & pop, int subpop, const vectori & loci, const vectori & alleles,
	vectorlu & numAllele);

/// a controlled mating scheme
/**
   This is an experimental mating scheme that uses a frequency range to control
   the allele frequency of the offspring generation at given loci. When allele frequencies
   at the offspring generation does not fall into the given range, the offspring generation
   is regenerated. Any mating scheme can be used with this mating scheme by passing through
   parameter \c matingScheme.
   <applicability>diploid only</applicability>
 */
class controlledMating : public mating
{
public:
	/// control allele frequencies at a locus
	/**
	 \param matingScheme a mating scheme
	 \param loci loci at which allele frequency is controlled. Note that controlling
	   	the allele frequencies at several loci may take a long time.
	 \param alleles alleles to control at each locus. Should have the same length as \c loci.
	 \param freqFunc frequency boundaries. If the length of the returned value equals the size
	   	of \c loci, the range for loci will be <tt>[value0, value0+range]</tt>,
	   	<tt>[value1, value1+range]</tt> etc. If the length of the returned value is 2 times
	   	the size of \c loci, it will be interpreted as <tt>[low1, high1, low2, high2, ...]</tt>.
	 */
	controlledMating(
	                 mating & matingScheme,
	                 vectori loci,
	                 vectori alleles,
	                 PyObject * freqFunc,
	                 double range = 0.01
	                 )
		: m_loci(loci), m_alleles(alleles),
		m_freqFunc(freqFunc),
		m_range(range)
	{
		DBG_FAILIF(m_loci.empty(), ValueError, "Have to specify a locus (or a loci) to control");

		DBG_FAILIF(m_alleles.empty(), ValueError, "Have to specify allele at each locus");

		DBG_FAILIF(m_loci.size() != m_alleles.size(), ValueError, "Should specify allele for each locus");

		if (m_freqFunc == NULL || !PyCallable_Check(m_freqFunc))
			throw ValueError("Please specify a valid frequency function");
		else
			Py_INCREF(m_freqFunc);

		m_matingScheme = matingScheme.clone();
	}


	/// CPPONLY
	controlledMating(const controlledMating & rhs)
		: m_loci(rhs.m_loci),
		m_alleles(rhs.m_alleles),
		m_freqFunc(rhs.m_freqFunc),
		m_range(rhs.m_range)
	{
		Py_INCREF(m_freqFunc);
		m_matingScheme = rhs.m_matingScheme->clone();
	}


	/// destructor
	~controlledMating()
	{
		if (m_freqFunc != NULL)
			Py_DECREF(m_freqFunc);
		delete m_matingScheme;
	}


	/// CPPONLY
	void submitScratch(population & pop, population & scratch)
	{
	}


	/// deep copy of a controlled mating scheme
	virtual mating * clone() const
	{
		return new controlledMating(*this);
	}


	/// CPPONLY
	virtual bool isCompatible(const population & pop) const
	{
		return m_matingScheme->isCompatible(pop);
	}


	/// used by Python print function to print out the general information of the controlled mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::controlled mating>";
	}


	/// CPPONLY perform controlled mating
	virtual bool mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit);

private:
	/// mating scheme
	mating * m_matingScheme;

	/// loci at which mating is controlled.
	vectori m_loci;

	/// allele to be controlled at each locus
	vectori m_alleles;

	/// function that return an array of frquency range
	PyObject * m_freqFunc;

	/// range, used when m_freqFunc returns a vector of the same length as m_loci
	double m_range;
};

/// CPPONLY
void getExpectedAlleles(population & pop, vectorf & expFreq, const vectori & loci, const vectori & alleles,
	vectoru & expAlleles);

/// a controlled random mating scheme
/** This is the controlled random mating scheme described in
   <em> Peng 2007 (PLoS Genetics) </em>. Basically, a \c freqFunc
   is passed to this mating scheme and set the allele frequencies of given
   alleles at given loci at the offspring generation.
 \n
   The offspring generation is conceptually populated in two steps.
   At the first step, only families with disease alleles are accepted
   until the expected number of disease alleles are met. At the second
   step, only families with wide type alleles are accepted to populate
   the rest of the offspring generation.
   <applicability>diploid only</applicability>
 */
class controlledRandomMating : public randomMating
{
public:
	/// create a controlled random mating scheme
	/**
	 \param loci loci at which allele frequencies are monitored (controlled)
	 \param alleles alleles at given loci. It should have the same length as \c loci
	 \param freqFunc a Python function that accepts a generation number and returns
	   	expected allele frequencies at given loci
	 \param acceptScheme internal use only

	   Please refer to class \c mating for descriptions of other parameters.
	 */
	controlledRandomMating(
	                       vectori loci,
	                       vectori alleles,
	                       PyObject * freqFunc,
	                       int acceptScheme = 0,
	                       double numOffspring = 1.,
	                       double sexParam = 0.5,
	                       UINT sexMode = MATE_RandomSex,
	                       PyObject * numOffspringFunc = NULL,
	                       UINT maxNumOffspring = 0,
	                       UINT mode = MATE_NumOffspring,
	                       vectorlu newSubPopSize = vectorlu(),
	                       PyObject * newSubPopSizeFunc = NULL,
	                       string newSubPopSizeExpr = "",
	                       bool contWhenUniSex = true,
						   vspID subPop = vspID(),
	                       double weight = 0)
		: randomMating(numOffspring,
		               numOffspringFunc, maxNumOffspring, mode,
		               sexParam, sexMode,
		               newSubPopSize,
		               newSubPopSizeFunc,
		               newSubPopSizeExpr,
		               contWhenUniSex,
		               subPop, weight),
		m_loci(loci),
		m_alleles(alleles),
		m_freqFunc(freqFunc),
		m_stack()
	{
		if (m_freqFunc == NULL || !PyCallable_Check(m_freqFunc))
			throw ValueError("Please specify a valid frequency function");
		else
			Py_INCREF(m_freqFunc);
	}


	/// CPPONLY
	controlledRandomMating(const controlledRandomMating & rhs)
		: randomMating(rhs),
		m_loci(rhs.m_loci),
		m_alleles(rhs.m_alleles),
		m_freqFunc(rhs.m_freqFunc),
		m_stack()
	{
		Py_INCREF(m_freqFunc);
	}


	/// destructor
	~controlledRandomMating()
	{
		if (m_freqFunc != NULL)
			Py_DECREF(m_freqFunc);
	}


	/// deep copy of a controlled random mating scheme
	virtual mating * clone() const
	{
		return new controlledRandomMating(*this);
	}


	/// CPPONLY
	virtual bool isCompatible(const population & pop) const
	{
#ifndef OPTIMIZED
		if (pop.ploidy() != 2)
			cout << "Warning: This mating type only works with diploid population." << endl;
#endif

		return true;
	}


	/// used by Python print function to print out the general information of the controlled random mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::controlled sexual random mating>";
	}


	/// CPPONLY
	void submitScratch(population & pop, population & scratch)
	{
		pop.turnOffSelection();
		// use scratch population,
		pop.pushAndDiscard(scratch);
		DBG_DO(DBG_MATING, pop.setIntVectorVar("famSizes", m_famSize));
	}


	/// CPPONLY perform controlled random mating
	virtual bool mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit);

private:
	/// locus at which mating is controlled.
	vectori m_loci;

	/// allele to be controlled at each locus
	vectori m_alleles;

	/// function that return an array of frquency range
	PyObject * m_freqFunc;

	///
	stack<RawIndIterator> m_stack;
};

/// a Python mating scheme
/**
   	This hybrid mating scheme does not have to involve a python function.
   	It requires a parent chooser, and an offspring generator. The parent
   	chooser chooses parent(s) and pass them to the offspring generator to
   	produce offspring.
   <applicability>all ploidy</applicability>
 */
class pyMating : public mating
{
public:
	/// create a Python mating scheme
	/**
	 \param chooser a parent chooser that chooses parent(s) from the parental
	   	generation.
	 \param generator an offspring generator that produce offspring of given
	   	parents.

	 */
	pyMating(parentChooser & chooser,
	         offspringGenerator & generator,
	         vectorlu newSubPopSize = vectorlu(),
	         string newSubPopSizeExpr = "",
	         PyObject * newSubPopSizeFunc = NULL,
	         vspID subPop = vspID(),
	         double weight = 0);

	/// destructor
	~pyMating()
	{
		delete m_parentChooser;
		delete m_offspringGenerator;
	}


	/// CPPONLY
	pyMating(const pyMating & rhs) :
		mating(rhs)
	{
		m_offspringGenerator = rhs.m_offspringGenerator->clone();
		m_parentChooser = rhs.m_parentChooser->clone();
	}


	/// deep copy of a Python mating scheme
	virtual mating * clone() const
	{
		return new pyMating(*this);
	}


	/// used by Python print function to print out the general information of the Python mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::pyMating>";
	}


	/// CPPONLY perform Python mating
	/**
	   All individuals will be passed to during mating operators but
	   no one will die (ignore during mating failing signal).
	 */
	virtual bool mateSubPop(population & pop, SubPopID subPop,
		RawIndIterator offBegin, RawIndIterator offEnd,
		vector<baseOperator * > & ops);

private:
	parentChooser * m_parentChooser;
	offspringGenerator * m_offspringGenerator;

};

typedef std::vector<mating *> vectormating;

/** a heterogeneous mating scheme that applies a list of mating
   schemes to different (virtual) subpopulations.
   <applicability>diploid only</applicability>
 */
class heteroMating : public mating
{
public:
	/// create a heterogeneous Python mating scheme
	/**
	 \param matingSchemes A list of mating schemes. If parameter \c subPop of an
	   	mating scheme is specified, it will be applied to specific subpopulation.
	   	If \c virtualSubPop if specified, it will be applied to specifc virtual
	   	subpopulations.

	   Parameter subpop, virtualSubPOp and weight of this mating scheme is ignored.
	 */
	heteroMating(const vectormating & matingSchemes,
	             vectorlu newSubPopSize = vectorlu(),
	             string newSubPopSizeExpr = "",
	             PyObject * newSubPopSizeFunc = NULL,
	             bool shuffleOffspring = true,
	             vspID subPop = vspID(),
	             double weight = 0);

	/// destructor
	~heteroMating();

	/// CPPONLY
	heteroMating(const heteroMating & rhs);

	/// deep copy of a Python mating scheme
	virtual mating * clone() const
	{
		return new heteroMating(*this);
	}


	/// used by Python print function to print out the general information of the Python mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::heteroMating>";
	}


	/*
	   mateSubPop is not redefined. They are called by each
	   mating schemes.
	 */
	bool mate(population & pop, population & scratch, vector<baseOperator * > & ops, bool submit);

private:
	vectormating m_matingSchemes;
	///
	bool m_shuffleOffspring;
};


}
#endif
