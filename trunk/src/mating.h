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

/** An <em>offspring generator</em> generates offspring from parents chosen by
 *  a parent chooser. It is responsible for creating a certain number of
 *  offspring, determinning their sex, and transmitting genotypes from parents
 *  to offspring.\n
 */
class offspringGenerator
{
public:
	/** Create a basic offspring generator. This offspring generator uses
	 *  \e ops genotype transmitters to transmit genotypes from parents to
	 *  offspring. It expects \e numParents from an upstream parents chooser
	 *  and raises an \c RuntimeError if incorrect number of parents are
	 *  passed. If both one and two parents can be handled, \c 0 should be
	 *  specified for this parameter.\n
	 *
	 *  A number of <em>genotype transmitters</em> can be used to transmit
	 *  genotype from parents to offspring. Additional during-mating operators
	 *  can be passed from the \c evolve() function of a \e simulator, but the
	 *  \e ops operators will be applied before them. An exception is that if
	 *  one of the passed operators is set to form offspring genotype (a flag
	 *  \c setOffGenotype), operators in \e ops with the same flag will not be
	 *  applied. For example, a \c recombinator will override a
	 *  \c mendelianGenoTransmitter used in \c randomMating if it is used in
	 *  the \c ops parameter of the \c evolve function.\n
	 *
	 *  A number of derived offspring generators are available with a default
	 *  transmitter. For example, a \c mendelianOffspringGenerator uses
	 *  a \c mendelianGenoTransmitter to transmit genotypes.\n
	 *
	 *  Parameter \e numOffspring is used to control the number of offspring
	 *  per mating event, or in another word the number of offspring in each
	 *  family. It can be a number, a function, or a mode parameter followed by
	 *  some optional arguments. If a number is given, given number of
	 *  offspring will be generated at each mating event. If a Python function
	 *  is given, it will be called each time when a mating event happens.
	 *  Current generation number will be passed to this function, and its
	 *  return value will be considered the number of offspring. In the last
	 *  case, a tuple (or a list) in one of the following forms:
	 *  <tt>(GeometricDistribution, p)</tt>, <tt>(PoissonDistribution, p)</tt>,
	 *  <tt>(BinomialDistribution, p, N)</tt>, or
	 *  <tt>(UniformDistribution, a, b)</tt> can be given. The number of
	 *  offspring will be determined randomly following these statistical
	 *  distributions. Please refer to the simuPOP user's guide for a detailed
	 *  description of these distribution and their parameters.\n
	 *
	 *  Parameter \e sexMode is used to control the sex of each offspring. Its
	 *  default value is usually \e RandomSex which assign \c Male or \c Female
	 *  to each individual randomly, with equal probabilities. If \c NoSex is
	 *  given, all individuals will be \c Male. \e sexMode can also be one of
	 *  <tt>(ProbOfMale, p)</tt>, <tt>(NumOfMale, n)</tt>, and
	 *  <tt>(NumOfFemale, n)</tt>. The first case specifies the probability
	 *  of male for each offspring. The next two cases specifies the number of
	 *  male or female individuals in each family, respectively. If \c n is
	 *  greater than or equal to the number of offspring in this family, all
	 *  offspring in this family will be \c Male or \c Female.
	 */
	offspringGenerator(const vectorop & ops, const floatListFunc & numOffspring = 1,
		const floatList & sexMode = RandomSex);

	virtual ~offspringGenerator()
	{
		for (size_t i = 0; i < m_transmitters.size(); ++i)
			delete m_transmitters[i];
	}


	/// CPPONLY
	offspringGenerator(const offspringGenerator & rhs);

	/// Make a deep copy of this offspring generator.
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

	/// CPPONLY
	virtual void finalize(const population & pop)
	{
		m_initialized = false;
	}


	/// CPPONLY
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

protected:
	/// number of offspring
	floatListFunc m_numOffspring;

	/// paramter to determine offspring sex
	floatList m_sexMode;

	/// default transmitter
	vectorop m_transmitters;

protected:
	/// see if who will generate offspring genotype
	bool m_formOffGenotype;

	bool m_initialized;
};


/** This offspring generator populates an offspring population and controls
 *  allele frequencies at specified loci. At each generation, expected allele
 *  frequencies at these loci are passed from a user defined allele frequency
 *  \e trajectory function. The offspring population is populated in two steps.
 *  At the first step, only families with disease alleles are accepted until
 *  until the expected number of disease alleles are met. At the second step,
 *  only families with wide type alleles are accepted to populate the rest of
 *  the offspring generation. This method is described in detail in
 *  "Peng et al, (2007) <em>Forward-time simulations of populations with
 *  complex human diseases</em>, PLoS Genetics".
 */
class controlledOffspringGenerator : public offspringGenerator
{
public:
	/** Create an offspring generator that selects offspring so that allele
	 *  frequency at specified loci in the offspring generation reaches specified
	 *  allele frequency. At the beginning of each generation, expected allele
	 *  frequency of \e alleles at \e loci is returned from a user-defined
	 *  trajectory function \e freqFunc. If there are multiple subpopulations,
	 *  \e freqFunc can return a list of allele frequencies for each subpopulation,
	 *  or a list of allele frequencies in the whole population. In the latter
	 *  case, overall expected number of alleles are scattered to each
	 *  subpopulation in proportion to existing number of alleles in each
	 *  subpopulation, using a multi-nomial distribution.\n
	 *
	 *  After the expected alleles are calculated, this offspring generator
	 *  accept and reject families according to their genotype at \e loci until
	 *  allele frequecies reach their expected values. The rest of the
	 *  offspring generation is then filled with families without only wild
	 *  type alleles at these \e loci.\n
	 *
	 *  This offspring generator is derived from class \e offspringGenerator.
	 *  Please refer to class \e offspringGenerator for a detailed description
	 *  of parameters \e ops, \e numOffspring and \e sexMode.
	 */
	controlledOffspringGenerator(const vectori & loci, const vectori & alleles,
		PyObject * freqFunc, const vectorop & ops = vectorop(), 
		const floatListFunc & numOffspring = 1, const floatList & sexMode = RandomSex);


	/// CPPONLY
	controlledOffspringGenerator(const controlledOffspringGenerator & rhs);

	/// CPPONLY
	void initialize(const population & pop, SubPopID subPop, vector<baseOperator *> const & ops);

	/// CPPONLY
	virtual UINT generateOffspring(population & pop, individual * dad, individual * mom,
		RawIndIterator & offBegin,
		RawIndIterator & offEnd,
		vector<baseOperator *> & ops);

	/// Deep copy of a controlled random mating scheme
	virtual offspringGenerator * clone() const
	{
		return new controlledOffspringGenerator(*this);
	}


private:
	void getExpectedAlleles(const population & pop, vectorf & expFreq);

	/// locus at which mating is controlled.
	vectori m_loci;
	//
	/// allele to be controlled at each locus
	vectori m_alleles;

	/// function that return an array of frquency range
	pyFunc m_freqFunc;
	//
	// expected alleles
	vectoru m_expAlleles;
	vectoru m_flip;         // in a subpop
	vectoru m_totAllele;    // in a subpop
	vectoru m_curAllele;    // in a subpop
	//
	int m_AAattempt;
	int m_aaAttempt;
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
	parentChooser() : m_initialized(false)
	{
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
	virtual individualPair chooseParents(RawIndIterator basePtr)
	{
		return individualPair(NULL, NULL);
	}


	virtual ~parentChooser() { }

protected:
	bool m_initialized;
};


/** This parent chooser chooses a parent linearly, regardless of sex
   or fitness values (selection is not considered).
   <applicability>all ploidy</applicability>
 */
class sequentialParentChooser : public parentChooser
{
public:
	sequentialParentChooser() : parentChooser()
	{
	}


	parentChooser * clone() const
	{
		return new sequentialParentChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	/// CPPONLY
	individualPair chooseParents(RawIndIterator basePtr);

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
		parentChooser(), m_maleIndex(0), m_femaleIndex(0),
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
		parentChooser(), m_replacement(replacement),
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
	individualPair chooseParents(RawIndIterator basePtr);

protected:
	bool m_replacement;

	bool m_selection;
	///
	vector<RawIndIterator> m_index;
	vector<RawIndIterator> m_chosen;
	/// accumulative fitness
	weightedSampler m_sampler;
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
		parentChooser(), m_replacement(replacement),
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
	weightedSampler m_malesampler;
	weightedSampler m_femalesampler;

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
		parentChooser(),
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
	weightedSampler m_malesampler;
	weightedSampler m_femalesampler;
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
		parentChooser(),
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
	weightedSampler m_malesampler;
	weightedSampler m_femalesampler;
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
		m_infoFields(infoFields), m_func(func), m_param(param),
		m_infoIdx(0), m_degenerate(false)
	{
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

	pyFunc m_func;
	pyObject m_param;

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

	pyFunc m_func;
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
	   parameters. \c subPopSize,
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
	mating(const uintListFunc & subPopSize = uintListFunc(),
		vspID subPop = vspID(),
		double weight = 0);

	/// CPPONLY
	mating(const mating & rhs);

	/// destructor
	virtual ~mating()
	{
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
	uintListFunc m_subPopSize;

	///
	vspID m_subPop;

	///
	double m_weight;

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
		const uintListFunc & subPopSize = uintListFunc(),
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
		const uintListFunc & subPopSize = uintListFunc(),
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
