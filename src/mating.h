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
#include "population.h"
#include "operator.h"

#include <string>
#include <algorithm>
using std::max_element;
using std::string;

#include <stack>
using std::stack;

namespace simuPOP {

/**
   A mating scheme is a complex object that does many things.
   1. create offspring generation
   2. how to choose parents
   3. how to produce offspring per mating event

 */
/// CPPONLY
/** How to generate offspring.
 */
class offspringGenerator
{
public:
	/// create an offspring generator, save information from \c pop and \c ops to speed up the calls to \c generateOffspring
	offspringGenerator(const population & pop, vector<baseOperator *> & ops);

	/// generate \c numOff offspring
	virtual void generateOffspring(population & pop, individual * dad, individual * mom, UINT numOff,
	                               population::IndIterator & offBegin) = 0;

	virtual ~offspringGenerator()
	{
	}

protected:
	bool checkFormOffspringGenotype();

	/// see if who will generate offspring genotype
	bool m_formOffGenotype;

	/// cache during-mating operators
	/// we do not cache pop since it may be changed during mating.
	vector<baseOperator *> & m_ops;

#ifndef OPTIMIZED
	size_t m_genoStruIdx;
#endif
};


class cloneOffspringGenerator : public offspringGenerator
{
public:
	cloneOffspringGenerator(const population & pop, vector<baseOperator *> & ops)
		: offspringGenerator(pop, ops)
	{
	}


	// the default method to produce offspring
	void generateOffspring(population & pop, individual * dad, individual * mom, UINT numOff,
	                       population::IndIterator & offBegin);

};


class mendelianOffspringGenerator : public offspringGenerator
{
public:
	mendelianOffspringGenerator(const population & pop,
	                            vector<baseOperator *> & ops);

	void generateOffspring(population & pop, individual * dad, individual * mom, UINT numOff,
	                       population::IndIterator & offBegin);

	// the default method to produce offspring
	void formOffspring(individual * dad, individual * mom,
	                   population::IndIterator & it);

protected:
	// use bernullitrisls with p=0.5 for free recombination
	BernulliTrials m_bt;

	/// sex chromosome handling
	bool m_hasSexChrom;

	// cache chromBegin, chromEnd for better performance.
	vectoru m_chIdx;

};


class selfingOffspringGenerator : public mendelianOffspringGenerator
{
public:
	selfingOffspringGenerator(const population & pop, vector<baseOperator *> & ops)
		: mendelianOffspringGenerator(pop, ops)
	{
	}


	// the default method to produce offspring
	void generateOffspring(population & pop, individual * parent, individual *, UINT numOff,
	                       population::IndIterator & offBegin);

	void formOffspring(individual * parent, population::IndIterator & it);

};


// a base class for parent choosers
class parentChooser
{
public:
	typedef std::pair<individual *, individual *> individualPair;

public:
	parentChooser(int numParents)
	{
		m_numParents = numParents;
	}


	int numParents()
	{
		return m_numParents;
	}


	virtual individual * chooseParent()
	{
		return NULL;
	}


	virtual individualPair chooseParents()
	{
		return individualPair(NULL, NULL);
	}


	virtual ~parentChooser() { }

private:
	int m_numParents;

};


/// choose a parent randomly
class randomParentChooser : public parentChooser
{
public:
	randomParentChooser(population & pop, size_t sp);

	individual * chooseParent();

private:
	bool m_selection;
	/// accumulative fitness
	Weightedsampler m_sampler;
	/// starting individual
	population::IndIterator m_begin;
	/// individuals to choose
	size_t m_size;
};


/// choose two parents randomly, considering selection
class randomParentsChooser : public parentChooser
{
public:
	randomParentsChooser(population & pop, size_t sp);

	individualPair chooseParents();

	ULONG numMale() { return m_numMale; }
	ULONG numFemale() { return m_numFemale; }

private:
	bool m_selection;

	ULONG m_numMale;
	ULONG m_numFemale;

	/// internal index to female/males.
	vectorlu m_maleIndex;
	vectorlu m_femaleIndex;

	vectorf m_maleFitness;
	vectorf m_femaleFitness;

	// weighted sampler
	Weightedsampler m_malesampler;
	Weightedsampler m_femalesampler;

	/// starting individual of the whole population
	/// because m_maleIndex is relative to the whole population
	population::IndIterator m_begin;
};


/// choose a parent using a Python generator,
/// does not consider selection
class pyParentChooser : public parentChooser
{
public:
	pyParentChooser(population & pop, size_t sp, PyObject * parentGenerator);

	/// destructor
	~pyParentChooser()
	{
	}


	individual * chooseParent();

private:
	PyObject * m_parentGenerator;

#ifndef OPTIMIZED
	ULONG m_size;
#endif
	population::IndIterator m_begin;
};


/// choose two parents using a Python generator
/// does not consider selection
class pyParentsChooser : public parentChooser
{
public:
	pyParentsChooser(population & pop, size_t sp, PyObject * parentsGenerator);

	/// destructor
	~pyParentsChooser()
	{
	}


	individualPair chooseParents();

private:
	PyObject * m_parentsGenerator;

#ifndef OPTIMIZED
	ULONG m_size;
#endif
	population::IndIterator m_begin;
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
	// numOffspring: constant, numOffspringFunc: call each time before mating
#define MATE_NumOffspring           1
	// call numOffspringFunc each time during mating.
#define MATE_NumOffspringEachFamily 2
	// numOffspring and numOffsrpingsFunc call each time before mating is
	// the p for a geometric distribution
#define MATE_GeometricDistribution   3
#define MATE_PoissonDistribution     4
#define MATE_BinomialDistribution    5
	// uniform between numOffspring and maxNumOffspring
#define MATE_UniformDistribution     6

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
	 \c mode, \c numOffspring, \c maxNumOffspring can be used to specify how
	   many offspring will be produced at each mating event. This \c mode parameter
	   can be one of
	 \li \c MATE_NumOffspring: a fixed number of offspring at all mating events at this generation.
	   	If \c numOffspring is given, all generations use this fixed number. If \c numOffspringFunc
	   	is given, the number of offspring at each generation is determined by the value
	   returned from this function.
	 \li \c MATE_NumOffspringEachFamily: each family can have its own number of offspring.
	   Usually, \c numOffspringFunc is used to determine the number of offspring of each
	   family. If \c numOffspring is used, \c MATE_NumOffspringEachFamily is equivalent to
	 \c MATE_NumOffspring.
	 \li \c MATE_GeometricDistribution: a Geometric distribution with parameter \c numOffspring
	   is used to determine the number of offspring of each family.
	 \li \c MATE_PoissonDistribution: a Poisson distribution with parameter \c numOffspring
	   is used to determine the number of offspring of each family.
	 \li \c MATE_BinomialDistribution: a Binomial distribution with parameter \c numOffspring
	   is used to determine the number of offspring of each family.
	 \li \c MATE_UniformDistribution: a Uniform <tt> [a, b] </tt> distribution with parameter
	 \c numOffspring (a) and \c maxNumOffspring (b) is used to determine the number of offspring of each family.

	 \param numOffspring the number of offspring or \em p for a random distribution.
	   Default to \c 1. This parameter determines the number of offspring that a
	   mating event will produce. Therefore, it determines the family size.
	 \param numOffspringFunc a Python function that returns the number of offspring or \em p
	 \param maxNumOffspring used when \c numOffspring is generated from a binomial distribution
	 \param mode can be one of <tt>MATE_NumOffspring, MATE_NumOffspringEachFamily,
	   MATE_GeometricDistribution, MATE_PoissonDistribution, MATE_BinomialDistribution,</tt>
	   or <tt>MATE_UniformDistribution</tt>.
	 \param newSubPopSize an array of subpopulations sizes, should have the same
	   	number of subpopulations as the current population
	 \param newSubPopSizeExpr an expression that will be evaluated as an array of new subpopulation sizes
	 \param newSubPopSizeFunc a function that takes parameters \c gen (generation number) and \c oldsize
	   (an array of current population size) and return an array of subpopulation sizes of the next generation.
	   This is usually easier to use than its expression version of this parameter.

	 \test src_mating.log Demographic models and control of number of offspring per mating event
	 */
	mating(double numOffspring = 1.0,
	       PyObject * numOffspringFunc = NULL,
	       UINT maxNumOffspring = 0,
	       UINT mode = MATE_NumOffspring,
	       vectorlu newSubPopSize = vectorlu(),
	       string newSubPopSizeExpr = "",
	       PyObject * newSubPopSizeFunc = NULL);

	/// CPPONLY
	mating(const mating & rhs)
		: m_numOffspring(rhs.m_numOffspring),
		m_numOffspringFunc(rhs.m_numOffspringFunc),
		m_maxNumOffspring(rhs.m_maxNumOffspring),
		m_mode(rhs.m_mode), m_firstOffspring(true),
		m_subPopSize(rhs.m_subPopSize),
		m_subPopSizeExpr(rhs.m_subPopSizeExpr),
		m_subPopSizeFunc(rhs.m_subPopSizeFunc)
	{
		if (m_subPopSizeFunc != NULL)
			Py_INCREF(m_subPopSizeFunc);
		if (m_numOffspringFunc != NULL)
			Py_INCREF(m_numOffspringFunc);
	}


	/// destructor
	virtual ~mating()
	{
		if (m_subPopSizeFunc != NULL)
			Py_DECREF(m_subPopSizeFunc);
		if (m_numOffspringFunc != NULL)
			Py_DECREF(m_numOffspringFunc);
	}


	/// deep copy of a mating scheme
	virtual mating * clone() const
	{
		return new mating(*this);
	}


	/// used by Python print function to print out the general information of the mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::generic mating scheme>";
	}


	/// CPPONLY
	virtual void submitScratch(population & pop, population & scratch)
	{
	}


	/// CPPONLY this is not supposed to be called for a base mating class
	/**
	 \param pop population
	 \param scratch scratch population
	 \param ops during-mating operators
	 \return return false when mating fail
	 */
	virtual bool mate(population & pop, population & scratch, vector<baseOperator * > & ops, bool submit)
	{
		throw SystemError("You are not supposed to call base mating scheme.");
	}


	/// CPPONLY
	bool fixedFamilySize();

	/// CPPONLY the number of offspring of a genaration \c gen
	/**
	   This is called whenever a family size is needed.
	   Its actual meaning depending on \c mode.
	 */
	ULONG numOffspring(int gen);

	/// CPPONLY reset the number of offspring
	/**
	   This is used to tell a mating scheme that the current number
	   of offspring need to be recaculated. e.g., when a new generation
	   starts and mode=MATE_numOffspring and numOffspringFunc is given.
	 */
	void resetNumOffspring()
	{
		m_firstOffspring =  true;
	}


public:
	/// CPPONLY dealing with the \c pop/subPop size change, copy of structure etc.
	void prepareScratchPop(population & pop, population & scratch);

protected:
	/// number of offspring each mate
	double m_numOffspring;

	/// number of offspring func
	PyObject * m_numOffspringFunc;

	///
	UINT m_maxNumOffspring;

	/// whether or not call m_numOffspringFunc each time
	UINT m_mode;

	///
	bool m_firstOffspring;

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
};

/// a mating scheme that does nothing
/**
   In this scheme, there is
 \li no mating. Parent generation will be considered as offspring generation.
 \li no subpopulation change. \em During-mating operators will be applied, but
   	the return values are not checked. I.e., subpopulation size parameters will be ignored
   	although some during-mating operators might be applied.
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
	         vectorlu newSubPopSize = vectorlu(),
	         string newSubPopSizeExpr = "",
	         PyObject * newSubPopSizeFunc = NULL)
	{

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
	virtual void submitScratch(population & pop, population & scratch)
	{
	}


	/// CPPONLY perform mating, but no mating in this scheme
	/**
	   All individuals will be passed to during-mating operators but
	   no one will die (ignore the during-mating failing signal).
	 */
	virtual bool mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit);

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
	                  vectorlu newSubPopSize = vectorlu(),
	                  string newSubPopSizeExpr = "",
	                  PyObject * newSubPopSizeFunc = NULL)
		: mating(numOffspring,
		         numOffspringFunc, maxNumOffspring, mode,
		         newSubPopSize, newSubPopSizeExpr,
		         newSubPopSizeFunc)
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


	///CPPONLY
	virtual void submitScratch(population & pop, population & scratch)
	{
		pop.turnOffSelection();
		// use scratch population,
		pop.pushAndDiscard(scratch);
		DBG_DO(DBG_MATING, pop.setIntVectorVar("famSizes", m_famSize));
	}


	/// CPPONLY perform sexless random mating
	/**
	 \param pop population
	 \param scratch scratch population
	 \param ops during-mating operators
	 \return return false when mating fails
	 */
	virtual bool mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit);

protected:
#ifndef OPTIMIZED
	///
	vectori m_famSize;
#endif
};

/// a mating scheme of basic sexually random mating
/**
   In this scheme, sex information is considered for each individual,
   and ploidy is always 2. Within each subpopulation, males and females
   are randomly chosen. Then randomly get one copy of chromosomes from
   father and mother. If only one sex exists in a subpopulation, a
   parameter (\c contWhenUniSex) can be set to determine the behavior.
   Default to continuing without warning.
 */
class randomMating : public mating
{
public:
	/// create a random mating scheme
	/**
	 \param contWhenUniSex continue when there is only one sex in the population. Default to \c True.
	 \n

	   Please refer to class \c mating for descriptions of other parameters.
	 */
	randomMating(double numOffspring = 1.,
	             PyObject * numOffspringFunc = NULL,
	             UINT maxNumOffspring = 0,
	             UINT mode = MATE_NumOffspring,
	             vectorlu newSubPopSize = vectorlu(),
	             PyObject * newSubPopSizeFunc = NULL,
	             string newSubPopSizeExpr = "",
	             bool contWhenUniSex = true)
		: mating(numOffspring,
		         numOffspringFunc, maxNumOffspring, mode,
		         newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc),
		m_contWhenUniSex(contWhenUniSex)
	{
	}


	/// destructor
	~randomMating()
	{
	}


	/// deep copy of a random mating scheme
	virtual mating * clone() const
	{
		return new randomMating(*this);
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


	/// CPPONLY
	virtual void submitScratch(population & pop, population & scratch)
	{
		pop.turnOffSelection();
		// use scratch population,
		pop.pushAndDiscard(scratch);
		DBG_DO(DBG_MATING, pop.setIntVectorVar("famSizes", m_famSize));
	}


	/// CPPONLY perform random mating
	virtual bool mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit);

protected:
	/// if no other sex exist in a subpopulation,
	/// same sex mating will occur if m_contWhenUniSex is set.
	/// otherwise, an exception will be thrown.
	bool m_contWhenUniSex;
#ifndef OPTIMIZED
	///
	vectori m_famSize;
#endif
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
	                       PyObject * numOffspringFunc = NULL,
	                       UINT maxNumOffspring = 0,
	                       UINT mode = MATE_NumOffspring,
	                       vectorlu newSubPopSize = vectorlu(),
	                       PyObject * newSubPopSizeFunc = NULL,
	                       string newSubPopSizeExpr = "",
	                       bool contWhenUniSex = true)
		: randomMating(numOffspring,
		               numOffspringFunc, maxNumOffspring, mode,
		               newSubPopSize,
		               newSubPopSizeFunc,
		               newSubPopSizeExpr,
		               contWhenUniSex),
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
		// test if individual has sex
		// if not, will yield compile time error.
		pop.indBegin()->sex();

#ifndef OPTIMIZED
		if (pop.ploidy() != 2)
			cout << "Warning: This mating type only works with diploid population." << endl;
#endif

		return true;
	}


	/// used by Python print function to print out the general information of the controlled random mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::sexual random mating>";
	}


	/// CPPONLY
	virtual void submitScratch(population & pop, population & scratch)
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
	stack<population::IndIterator> m_stack;
};

/// a Python mating scheme
/**
   Hybird mating scheme. This mating scheme takes a Python generator
   that generate parents that will be mated by the mating scheme.
   The mating scheme will generate offspring population (controlled by
   newSubPopSize etc), call this function repeatedly to get parents,
   perform the mating, produce a number of offspring  (controlled by
   numOffspring etc), and apply given during mating operators.

   The parentsGenerator is not a usually Python function, rather a
   Python generator (use of yield keyword). Please refer to simuPOP
   user's guide for an example of how to use this mating scheme.

 */
class pyMating : public mating
{
public:
#define MATE_RandomParentChooser  1
#define MATE_RandomParentsChooser 2
#define MATE_PyParentChooser      3
#define MATE_PyParentsChooser     4

#define MATE_CloneOffspringGenerator     1
#define MATE_MendelianOffspringGenerator 2
#define MATE_SelfingOffspringGenerator   3

	/// create a Python mating scheme
	/**
	 \param parentsGenerator a Python generator that accepts the parental
	   	population, and yield parents.
	 \n

	   Please refer to class \c mating for descriptions of other parameters.
	 */
	pyMating(
	         vector<int> parentChooser = vector<int>(1, 2),
	         vector<PyObject *> pyChoosers = std::vector<PyObject *>(),
	         vector<int> offspringGenerator = vector<int>(1, 1),
	         double numOffspring = 1.,
	         PyObject * numOffspringFunc = NULL,
	         UINT maxNumOffspring = 0,
	         UINT mode = MATE_NumOffspring,
	         vectorlu newSubPopSize = vectorlu(),
	         string newSubPopSizeExpr = "",
	         PyObject * newSubPopSizeFunc = NULL
	         )
		: mating(numOffspring,
		         numOffspringFunc, maxNumOffspring, mode,
		         newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc),
		m_parentChooser(parentChooser),
		m_pyChoosers(),
		m_offspringGenerator(offspringGenerator)
	{
		vector<PyObject *>::iterator it = pyChoosers.begin();
		vector<PyObject *>::iterator it_end = pyChoosers.end();
		for (; it != it_end; ++it) {
#if PY_VERSION_HEX < 0x02040000
			throw SystemError("Your Python version does not have good support for generator"
			    " so operator pyMating can not be used.");
#else
			if (!PyGen_Check(*it))
				throw ValueError("Passed variable is not a Python generator.");

			Py_XINCREF(*it);
			m_pyChoosers.push_back(*it);
			Py_XINCREF(*it);
#endif
		}
	}


	/// destructor
	~pyMating()
	{
		vector<PyObject *>::iterator it = m_pyChoosers.begin();
		vector<PyObject *>::iterator it_end = m_pyChoosers.end();
		for (; it != it_end; ++it)
			Py_XDECREF(*it);
	}


	/// CPPONLY
	pyMating(const pyMating & rhs)
	{
		vector<PyObject *>::iterator it = m_pyChoosers.begin();
		vector<PyObject *>::iterator it_end = m_pyChoosers.end();
		for (; it != it_end; ++it)
			Py_INCREF(*it);
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
	virtual bool mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit);

private:
	vector<int> m_parentChooser;

	vector<PyObject *> m_pyChoosers;

	vector<int> m_offspringGenerator;

#ifndef OPTIMIZED
	///
	vectori m_famSize;
#endif
};

}
#endif
