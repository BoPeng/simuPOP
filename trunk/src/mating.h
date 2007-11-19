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

   This is the base class of all parent choosers, and should not
   be used directly.
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
	 */
	offspringGenerator(double numOffspring, PyObject * numOffspringFunc,
	                   UINT maxNumOffspring, UINT mode);

	virtual ~offspringGenerator()
	{
		if (m_numOffspringFunc != NULL)
			Py_DECREF(m_numOffspringFunc);
	}


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

	/// CPPONLY the number of offspring of a genaration \c gen
	/**
	   This is called whenever a family size is needed.
	   Its actual meaning depending on \c mode.
	 */
	ULONG numOffspring(int gen);


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

protected:
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
 \c maxNumOffspring and \c mode.
 */
class cloneOffspringGenerator : public offspringGenerator
{
public:
	cloneOffspringGenerator(double numOffspring = 1,
	                        PyObject * numOffspringFunc = NULL,
	                        UINT maxNumOffspring = 1,
	                        UINT mode = MATE_NumOffspring
	                        ) :
		offspringGenerator(numOffspring, numOffspringFunc, maxNumOffspring, mode)
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

   This offspring generator only works for diploid populations.
 */
class mendelianOffspringGenerator : public offspringGenerator
{
public:
	mendelianOffspringGenerator(double numOffspring = 1,
	                            PyObject * numOffspringFunc = NULL,
	                            UINT maxNumOffspring = 1,
	                            UINT mode = MATE_NumOffspring
	                            ) :
		offspringGenerator(numOffspring, numOffspringFunc, maxNumOffspring, mode),
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
	void formOffspringGenotype(individual * parent,
	                           RawIndIterator & it, int ploidy, bool setSex);

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
 */
class selfingOffspringGenerator : public mendelianOffspringGenerator
{
public:
	selfingOffspringGenerator(double numOffspring = 1,
	                          PyObject * numOffspringFunc = NULL,
	                          UINT maxNumOffspring = 1,
	                          UINT mode = MATE_NumOffspring
	                          )
		: mendelianOffspringGenerator(numOffspring, numOffspringFunc, maxNumOffspring, mode)
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
	virtual individual * chooseParent()
	{
		return NULL;
	}


	/// CPPONLY
	virtual individualPair chooseParents()
	{
		return individualPair(NULL, NULL);
	}


	virtual ~parentChooser() { }

private:
	int m_numParents;

protected:
	bool m_initialized;
};


/** This parent chooser chooses a parent linearly, regardless of sex
   or fitness values (selection is not considered).
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
	individual * chooseParent();

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
	individualPair chooseParents();

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
   fitness values. Sex is not considered.
 */
class randomParentChooser : public parentChooser
{
public:
	randomParentChooser() : parentChooser(1), m_sampler(rng())
	{
	}


	parentChooser * clone() const
	{
		return new randomParentChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	/// CPPONLY
	individual * chooseParent();

private:
	bool m_selection;
	/// accumulative fitness
	Weightedsampler m_sampler;
	/// starting individual
	IndIterator m_begin;
	/// individuals to choose
	size_t m_size;
};


/** This parent chooser chooses two parents randomly, a male
   and a female, from their respective sex groups randomly.
   If selection is turned on, parents are chosen from their sex
   groups with probabilities that are proportional to their
   fitness values.
 */
class randomParentsChooser : public parentChooser
{
public:
	randomParentsChooser() :
		parentChooser(2), m_maleIndex(0), m_femaleIndex(0),
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
	individualPair chooseParents();

	/// CPPONLY
	ULONG numMale() { return m_numMale; }
	/// CPPONLY
	ULONG numFemale() { return m_numFemale; }

private:
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


/** This parents chooser accept a Python generator function that yields
   repeatedly an index (relative to each subpopulation) of a parent, or
   indexes of two parents as a Python list of tuple. The generator function
   is responsible for handling sex or selection if needed.
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
	individualPair chooseParents();

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
	       SubPopID subPop = InvalidSubPopID,
	       SubPopID virtualSubPop = InvalidSubPopID,
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
		return m_subPop;
	}


	/// CPPONLY
	SubPopID virtualSubPop() const
	{
		return m_virtualSubPop;
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
	SubPopID m_subPop;

	///
	SubPopID m_virtualSubPop;

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
	         PyObject * newSubPopSizeFunc = NULL,
	         SubPopID subPop = InvalidSubPopID,
	         SubPopID virtualSubPop = InvalidSubPopID,
	         double weight = 0)
	{
		DBG_FAILIF(virtualSubPop != InvalidSubPopID, ValueError,
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
	            vectorlu newSubPopSize = vectorlu(),
	            string newSubPopSizeExpr = "",
	            PyObject * newSubPopSizeFunc = NULL,
	            SubPopID subPop = InvalidSubPopID,
	            SubPopID virtualSubPop = InvalidSubPopID,
	            double weight = 0)
		: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc, subPop, virtualSubPop, weight),
		m_offGenerator(numOffspring, numOffspringFunc, maxNumOffspring, mode)
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
	                  PyObject * newSubPopSizeFunc = NULL,
	                  SubPopID subPop = InvalidSubPopID,
	                  SubPopID virtualSubPop = InvalidSubPopID,
	                  double weight = 0)
		: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc, subPop, virtualSubPop, weight),
		m_offGenerator(numOffspring, numOffspringFunc, maxNumOffspring, mode)
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
	             bool contWhenUniSex = true,
	             SubPopID subPop = InvalidSubPopID,
	             SubPopID virtualSubPop = InvalidSubPopID,
	             double weight = 0)
		: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc, subPop, virtualSubPop, weight),
		m_offspringGenerator(numOffspring,
		                     numOffspringFunc, maxNumOffspring, mode),
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


	/// CPPONLY perform random mating
	virtual bool mateSubPop(population & pop, SubPopID subPop,
	                        RawIndIterator offBegin, RawIndIterator offEnd,
	                        vector<baseOperator * > & ops);

protected:
	mendelianOffspringGenerator m_offspringGenerator;

	/// if no other sex exist in a subpopulation,
	/// same sex mating will occur if m_contWhenUniSex is set.
	/// otherwise, an exception will be thrown.
	bool m_contWhenUniSex;

};


/// a mating scheme of selfing
/**
   In this mating scheme, a parent is choosen randomly, acts
   both as father and mother in the usual random mating. The parent
   is chosen randomly, regardless of sex. If selection is turned on,
   the probability that an individual is chosen is proportional to
   his/her fitness.
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
	           vectorlu newSubPopSize = vectorlu(),
	           PyObject * newSubPopSizeFunc = NULL,
	           string newSubPopSizeExpr = "",
	           bool contWhenUniSex = true,
	           SubPopID subPop = InvalidSubPopID,
	           SubPopID virtualSubPop = InvalidSubPopID,
	           double weight = 0)
		: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc, subPop, virtualSubPop, weight),
		m_offspringGenerator(numOffspring,
		                     numOffspringFunc, maxNumOffspring, mode)
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
	                       bool contWhenUniSex = true,
	                       SubPopID subPop = InvalidSubPopID,
	                       SubPopID virtualSubPop = InvalidSubPopID,
	                       double weight = 0)
		: randomMating(numOffspring,
		               numOffspringFunc, maxNumOffspring, mode,
		               newSubPopSize,
		               newSubPopSizeFunc,
		               newSubPopSizeExpr,
		               contWhenUniSex,
		               subPop, virtualSubPop, weight),
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
	         SubPopID subPop = InvalidSubPopID,
	         SubPopID virtualSubPop = InvalidSubPopID,
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
 */
class heteroMating : public mating
{
public:
	/// create a heterogeneous Python mating scheme
	/**
	 \param matingSchemes A list of mating schemes. It parameter \c subPop of an
	   	mating scheme is specified, it will be applied to specific subpopulation.
	   	If \c virtualSubPop if specified, it will be applied to specifc virtual
	   	subpopulations. The \c weight parameter is used to control how many
	   	offspring to produce in case that more than one mating schemes are applied
	   	to the same subpopulation.
	 */
	heteroMating(const vectormating & matingSchemes,
	             vectorlu newSubPopSize = vectorlu(),
	             string newSubPopSizeExpr = "",
	             PyObject * newSubPopSizeFunc = NULL,
	             bool shuffleOffspring = true,
	             SubPopID subPop = InvalidSubPopID,
	             SubPopID virtualSubPop = InvalidSubPopID,
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
