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
#include "recombinator.h"

#include <string>
#include <algorithm>
using std::max_element;
using std::string;

#include <stack>
using std::stack;

namespace simuPOP {

/**
   Offspring generators generate offspring from given parents. Generators differ
   from each other by how and how many offspring is generated at each mating event.

   Parameters \c mode, \c numOffspring, \c numOffspringParam and \c numOffspringFunc
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
   \c numOffspring (a) and \c numOffspringParam (b) is used to determine the number of offspring of each family.

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
	// uniform between numOffspring and numOffspringParam
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
	   \param numOffspringParam used when \c numOffspring is generated from a binomial or random
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
	   If there are sex chromosomes, sex is determined by sex chromosomes when \c sexMode
	    id \c MATE_RandomSex. Otherwise, some offspring will be rejected so that offspring
	    sex match what is specified in other modes.
	    \param transmitter is an during mating operator, that will be used if
	        no during mating operator is used to produce offspring.
	 */
	offspringGenerator(const vectorop & ops, UINT numParents = 0,
		double numOffspring = 1., PyObject * numOffspringFunc = NULL,
		UINT numOffspringParam = 1, UINT mode = MATE_NumOffspring,
		double sexParam = 0.5, UINT sexMode = MATE_RandomSex);

	virtual ~offspringGenerator()
	{
		if (m_numOffspringFunc != NULL)
			Py_DECREF(m_numOffspringFunc);
		for (size_t i = 0; i < m_transmitters.size(); ++i)
			delete m_transmitters[i];
	}


	/// CPPONLY
	offspringGenerator(const offspringGenerator & rhs);

	virtual offspringGenerator * clone() const
	{
		return new offspringGenerator(*this);
	}


	/** create an offspring generator, save information from \c pop and \c ops to
	 *  speed up the calls to \c generateOffspring
	 *  CPPONLY
	 */
	virtual void initialize(const population & pop, SubPopID subPop, vector<baseOperator *> const & ops);

	/// CPPONLY
	virtual UINT generateOffspring(population & pop, individual * dad, individual * mom,
		RawIndIterator & offBegin,
		RawIndIterator & offEnd,
		vector<baseOperator *> & ops);

	virtual void finalize(const population & pop)
	{
		m_initialized = false;
	}


	bool initialized()
	{
		return m_initialized;
	}


	/** CPPONLY
	 *  Return the number of offspring of a genaration \e gen
	 *  This is called whenever a family size is needed.
	 *  Its actual meaning depending on \c mode.
	 */
	ULONG numOffspring(int gen);

	/** CPPONLY
	 *  return sex according to m_sexParam, m_sexMode and
	 *  \e count, which is the index of offspring
	 */
	Sex getSex(int count);

	/// CPPONLY
	int numParents() const
	{
		return m_numParents;
	}


protected:
	/// number of offspring each mate
	double m_numOffspring;

	/// number of offspring func
	PyObject * m_numOffspringFunc;

	///
	UINT m_numOffspringParam;

	/// whether or not call m_numOffspringFunc each time
	UINT m_mode;

	///
	double m_sexParam;

	///
	UINT m_sexMode;

	// number of parents needed
	int m_numParents;

	/// default transmitter
	vectorop m_transmitters;

protected:
	/// see if who will generate offspring genotype
	bool m_formOffGenotype;

	bool m_initialized;
};

/** The offspring generation is conceptually populated in two steps.
   At the first step, only families with disease alleles are accepted
   until the expected number of disease alleles are met. At the second
   step, only families with wide type alleles are accepted to populate
   the rest of the offspring generation.
 */
class controlledOffspringGenerator : public offspringGenerator
{
public:
	/**
	   \param loci loci at which allele frequencies are monitored (controlled)
	   \param alleles alleles at given loci. It should have the same length as \c loci
	   \param freqFunc a Python function that accepts a generation number and returns
	    expected allele frequencies at given loci
	   \param acceptScheme internal use only
	 */
	controlledOffspringGenerator(
		vectori loci,
		vectori alleles,
		PyObject * freqFunc,
		//
		const vectorop & ops = vectorop(),
		UINT numParents = 0,
		//
		double numOffspring = 1.,
		PyObject * numOffspringFunc = NULL,
		UINT numOffspringParam = NULL,
		UINT mode = MATE_NumOffspring,
		//
		double sexParam = 0.5,
		UINT sexMode = MATE_RandomSex);


	/// CPPONLY
	controlledOffspringGenerator(const controlledOffspringGenerator & rhs);

	/// destructor
	~controlledOffspringGenerator()
	{
		if (m_freqFunc != NULL)
			Py_DECREF(m_freqFunc);
	}


	void initialize(const population & pop, SubPopID subPop, vector<baseOperator *> const & ops);

	/// CPPONLY
	virtual UINT generateOffspring(population & pop, individual * dad, individual * mom,
		RawIndIterator & offBegin,
		RawIndIterator & offEnd,
		vector<baseOperator *> & ops);

	//
	//
	/// deep copy of a controlled random mating scheme
	virtual offspringGenerator * clone() const
	{
		return new controlledOffspringGenerator(*this);
	}


private:
	void getExpectedAlleles(const population & pop, vectorf & expFreq);

private:
	/// locus at which mating is controlled.
	vectori m_loci;
	//
	/// allele to be controlled at each locus
	vectori m_alleles;

	/// function that return an array of frquency range
	PyObject * m_freqFunc;
	//
	// expected alleles
	vectoru m_expAlleles;
	vectoru m_flip;         // in a subpop
	vectoru m_totAllele;    // in a subpop
	vectoru m_curAllele;    // in a subpop
	//
	int m_AAattempt;
	bool m_freqRequMet;
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


/** This parent chooser chooses a parent randomly from the
   parental generation. If selection is turned on, parents are
   chosen with probabilities that are proportional to their
   fitness values. Sex is not considered. Parameter \c replacement
   determines if a parent can be chosen multiple times.
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
	 */
	randomParentChooser(bool replacement = true) :
		parentChooser(1), m_replacement(replacement),
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
   fitness values. If replacement = False, each parent can
   only be used once.
   <applicability>all ploidy</applicability>
 */
class randomParentsChooser : public parentChooser
{
public:
	/**
	    Note: If selection is enabled, it works regularly on on-alpha sex, but
	    works twice on alpha sex. That is to say, \c alphaNum alpha indiviudals
	    are chosen selectively, and selected again during mating.
	 */
	randomParentsChooser(bool replacement = true) :
		parentChooser(2), m_replacement(replacement),
		m_maleIndex(0), m_femaleIndex(0),
		m_maleFitness(0), m_femaleFitness(0),
		m_malesampler(rng()), m_femalesampler(rng())
	{
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

	individual * m_lastParent;

	bool m_selection;

	ULONG m_numMale;
	ULONG m_numFemale;

	/// internal index to female/males.
	vector<RawIndIterator> m_maleIndex;
	vector<RawIndIterator> m_femaleIndex;

	vectorf m_maleFitness;
	vectorf m_femaleFitness;

	// weighted sampler
	Weightedsampler m_malesampler;
	Weightedsampler m_femalesampler;

};


/** This parent chooser chooses two parents randomly, a male
   and a female, from their respective sex groups randomly.
   If selection is turned on, parents are chosen from their sex
   groups with probabilities that are proportional to their
   fitness values. Note that
   selection is not allowed in the case of monopoly because this
   poses a particular order on individuals in the offspring generation.
   This parents chooser also allows polygamous mating by reusing
   a parent multiple times when returning parents, and allows
   specification of a few alpha individuals who will be the only
   mating individuals in their sex group.
   <applicability>all ploidy</applicability>
 */
class polyParentsChooser : public parentChooser
{
public:
	/**
	   \param polySex Male (polygyny) or Female (polyandry) parent that
	        will have \c polyNum sex partners.
	   \param polyNum Number of sex partners.

	   Note: If selection is enabled, it works regularly on on-alpha sex, but
	    works twice on alpha sex. That is to say, \c alphaNum alpha indiviudals
	    are chosen selectively, and selected again during mating.
	 */
	polyParentsChooser(Sex polySex = Male, UINT polyNum = 1) :
		parentChooser(2),
		m_polySex(polySex), m_polyNum(polyNum), m_polyCount(0),
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
		return new polyParentsChooser(*this);
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
	Sex m_polySex;
	UINT m_polyNum;

	UINT m_polyCount;

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


/** This parent chooser chooses two parents randomly, a male
   and a female, from their respective sex groups randomly.
   If selection is turned on, parents are chosen from their sex
   groups with probabilities that are proportional to their
   fitness values.  This parents chooser also allows polygamous mating by reusing
   a parent multiple times when returning parents, and allows
   specification of a few alpha individuals who will be the only
   mating individuals in their sex group.
   <applicability>all ploidy</applicability>
 */
class alphaParentsChooser : public parentChooser
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
	alphaParentsChooser(Sex alphaSex = Male, UINT alphaNum = 0, string alphaField = string()) :
		parentChooser(2),
		m_alphaSex(alphaSex), m_alphaNum(alphaNum), m_alphaField(alphaField),
		m_maleIndex(0), m_femaleIndex(0),
		m_maleFitness(0), m_femaleFitness(0),
		m_malesampler(rng()), m_femalesampler(rng())
	{
	}


	parentChooser * clone() const
	{
		return new alphaParentsChooser(*this);
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
	Sex m_alphaSex;
	UINT m_alphaNum;
	string m_alphaField;

	bool m_selection;

	ULONG m_numMale;
	ULONG m_numFemale;

	/// internal index to female/males.
	vector<RawIndIterator> m_maleIndex;
	vector<RawIndIterator> m_femaleIndex;

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
	 */
	infoParentsChooser(const vectorstr & infoFields = vectorstr(),
		PyObject * func = NULL,
		PyObject * param = NULL,
		bool replacement = true) :
		randomParentChooser(replacement),
		m_infoFields(infoFields), m_func(NULL), m_param(NULL),
		m_infoIdx(0), m_degenerate(false)
	{
		m_numParents = 2;
		DBG_FAILIF(m_infoFields.empty(), ValueError,
			"At least one information field should be provided for this infoParentsChooser");
		if (func && func != Py_None) {
			if (!PyCallable_Check(func))
				throw ValueError("Passed variable is not a callable Python function.");
			m_func = func;
			Py_XINCREF(func);
		}
		if (param && param != Py_None) {
			m_param = param;
			Py_XINCREF(param);
		}
	}


	~infoParentsChooser()
	{
		if (m_func)
			Py_DECREF(m_func);
		if (m_param)
			Py_DECREF(m_param);
	}


	infoParentsChooser(const infoParentsChooser & rhs) :
		randomParentChooser(rhs),
		m_infoFields(rhs.m_infoFields),
		m_func(rhs.m_func),
		m_param(rhs.m_param)
	{
		if (m_func)
			Py_INCREF(m_func);
		if (m_param)
			Py_INCREF(m_param);
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

	PyObject * m_func;
	PyObject * m_param;

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
	   parameters. \c subPopSize, and \c subPopSizeFunc
	   can be used to specify subpopulation sizes of the offspring generation.

	   \param subPopSize an array of subpopulations sizes, should have the same
	    number of subpopulations as the current population
	   \param subPopSizeFunc a function that takes parameters \c gen (generation number) and \c oldsize
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

	 */
	mating(vectorlu subPopSize = vectorlu(),
		PyObject * subPopSizeFunc = NULL,
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


	/** CPPONLY
	 * a common submit procedure is defined.
	 */
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

	/// CPPONLY dealing with the \c pop/subPop size change, copy of structure etc.
	bool prepareScratchPop(population & pop, population & scratch);

protected:
	/// new subpopulation size. mostly used to 'keep' subPopsize
	/// after migration.
	vectorlu m_subPopSize;

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


class pedigreeMating : public mating
{
public:
	pedigreeMating(const pedigree & ped, const offspringGenerator & generator,
		const string & fatherField = "father_idx", const string & motherField = "mother_idx");

	/// destructor
	~pedigreeMating();

	/// CPPONLY
	pedigreeMating(const pedigreeMating & rhs);


	/// deep copy of a Python mating scheme
	virtual mating * clone() const
	{
		return new pedigreeMating(*this);
	}


	/// used by Python print function to print out the general information of the Python mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::pedigreeMating>";
	}


	bool prepareScratchPop(population & pop, population & scratch);

	virtual bool mate(population & pop, population & scratch, vector<baseOperator * > & ops, bool submit);

private:
	pedigree m_ped;

	offspringGenerator * m_generator;

	string m_fatherField;
	string m_motherField;
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
		vectorlu subPopSize = vectorlu(),
		PyObject * subPopSizeFunc = NULL,
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
		vectorlu subPopSize = vectorlu(),
		PyObject * subPopSizeFunc = NULL,
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
