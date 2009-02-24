/**
 *  $File: mating.h $
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
#include "transmitter.h"

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
 *  to offspring.
 */
class offspringGenerator
{
public:
	/** Create a basic offspring generator. This offspring generator uses
	 *  \e ops genotype transmitters to transmit genotypes from parents to
	 *  offspring. It expects \e numParents from an upstream parents chooser
	 *  and raises an \c RuntimeError if incorrect number of parents are
	 *  passed. If both one and two parents can be handled, \c 0 should be
	 *  specified for this parameter.
	 *
	 *  A number of <em>genotype transmitters</em> can be used to transmit
	 *  genotype from parents to offspring. Additional during-mating operators
	 *  can be passed from the \c evolve() function of a \e simulator, but the
	 *  \e ops operators will be applied before them. An exception is that if
	 *  one of the passed operators is set to form offspring genotype (a flag
	 *  \c setOffGenotype), operators in \e ops with the same flag will not be
	 *  applied. For example, a \c recombinator will override a
	 *  \c mendelianGenoTransmitter used in \c randomMating if it is used in
	 *  the \c ops parameter of the \c evolve function. This general offspring
	 *  generator does not use any genotype transmitter. A number of derived
	 *  offspring generators are available with a default transmitter. For
	 *  example, a \c mendelianOffspringGenerator uses a
	 *  \c mendelianGenoTransmitter to transmit genotypes.
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
	 *  description of these distribution and their parameters.
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

	virtual ~offspringGenerator();

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
	vectorf m_sexMode;

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
	 *  subpopulation, using a multi-nomial distribution.
	 *
	 *  After the expected alleles are calculated, this offspring generator
	 *  accept and reject families according to their genotype at \e loci until
	 *  allele frequecies reach their expected values. The rest of the
	 *  offspring generation is then filled with families without only wild
	 *  type alleles at these \e loci.
	 *
	 *  This offspring generator is derived from class \e offspringGenerator.
	 *  Please refer to class \e offspringGenerator for a detailed description
	 *  of parameters \e ops, \e numOffspring and \e sexMode.
	 */
	controlledOffspringGenerator(const uintList & loci, const uintList & alleles,
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
	vectorlu m_loci;
	//
	/// allele to be controlled at each locus
	vectorlu m_alleles;

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


/** A parent chooser repeatedly chooses parent(s) from a parental population
 *  and pass them to an offspring generator. A parent chooser can select one or
 *  two parents, which should be matched by the offspring generator. This class
 *  is the base class of all parent choosers, and should not be used directly.
 */
class parentChooser
{
public:
	typedef std::pair<individual *, individual *> individualPair;

public:
	// CPPONLY
	parentChooser(const string & selectionField = string()) : m_initialized(false),
		m_selectionField(selectionField)
	{
	}


	/// Deep copy of a parent chooser
	virtual parentChooser * clone() const
	{
		return new parentChooser(*this);
	}


	/// CPPONLY
	virtual void initialize(population & pop, SubPopID subPop) { }

	/// CPPONLY
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


	/// destructor
	virtual ~parentChooser() { }

protected:
	bool m_initialized;

	string m_selectionField;
};


/** This parent chooser chooses a parent from a parental (virtual) subpopulation
 *  sequentially. Sex and selection is not considered. If the last parent is
 *  reached, this parent chooser will restart from the beginning of the
 *  (virtual) subpopulation.
 */
class sequentialParentChooser : public parentChooser
{
public:
	/** Create a parent chooser that chooses a parent from a parental (virtual)
	 *  subpopulation sequentially.
	 */
	sequentialParentChooser() : parentChooser()
	{
	}


	/// Deep copy of a sequential parent chooser.
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


/** This parent chooser chooses two parents (a father and a mother)
 *  sequentially from their respective sex groups. Selection is not considered.
 *  If all fathers (or mothers) are exhausted, this parent chooser will choose
 *  fathers (or mothers) from the beginning of the (virtual) subpopulation
 *  again.
 */
class sequentialParentsChooser : public parentChooser
{
public:
	/** Create a parent chooser that chooses two parents sequentially from a
	 *  parental (virtual) subpopulation.
	 */
	sequentialParentsChooser() :
		parentChooser(), m_maleIndex(0), m_femaleIndex(0),
		m_numMale(0), m_numFemale(0),
		m_curMale(0), m_curFemale(0)
	{
	}


	/// Deep copy of a sequential parents chooser.
	parentChooser * clone() const
	{
		return new sequentialParentsChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	/// CPPONLY
	individualPair chooseParents(RawIndIterator basePtr);

private:
	/// internal index to female/males.
	vector<RawIndIterator> m_maleIndex;
	vector<RawIndIterator> m_femaleIndex;

	ULONG m_numMale;
	ULONG m_numFemale;

	ULONG m_curMale;
	ULONG m_curFemale;
};


/** This parent chooser chooses a parent randomly from a (virtual) parental
 *  subpopulation. Parents are chosen with or without replacement.
 *  If parents are chosen with replacement, a parent can be selected multiple
 *  times. If natural selection is enabled, the probability that an individual
 *  is chosen is proportional to his/her fitness value stored in an information
 *  field \e selectionField (default to \c "fitness"). If parents are chosen
 *  without replacement, a parent can be chosen only once. An \c RuntimeError
 *  will be raised if all parents are exhausted. Selection is disabled in the
 *  without-replacement case.
 */
class randomParentChooser : public parentChooser
{
public:
	/** Create a random parent chooser that choose parents with or without
	 *  replacement (parameter \e replacement, default to \c True). If selection
	 *  is enabled and information field \e selectionField exists in the passed
	 *  population, the probability that a parent is chosen is proportional to
	 *  his/her fitness value stored in \e selectionField.
	 */
	randomParentChooser(bool replacement = true,
		const string & selectionField = "fitness") :
		parentChooser(selectionField), m_replacement(replacement),
		m_index(0), m_chosen(0), m_sampler(rng())
	{
	}


	/// Deep copy of a random parent chooser.
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


/** This parent chooser chooses two parents, a male and a female, randomly from
 *  a (virtual) parental subpopulation. Parents are chosen with or without
 *  replacement from their respective sex group. If parents are chosen with
 *  replacement, a parent can be selected multiple times. If natural selection
 *  is enabled, the probability that an individual is chosen is proportional to
 *  his/her fitness value among all individuals with the same sex. Selection
 *  will be disabled if specified information field \e selectionField (default
 *  to \c "fitness") does not exist.If parents are chosen without replacement,
 *  a parent can be chosen only once. An \c RuntimeError will be raised if all
 *  males or females are exhausted. Selection is disabled in the
 *  without-replacement case.
 */
class randomParentsChooser : public parentChooser
{
public:
	/** Create a random parents chooser that choose two parents with or without
	 *  replacement (parameter \e replacement, default to \c True). If selection
	 *  is enabled and information field \e selectionField exists in the passed
	 *  population, the probability that a parent is chosen is proportional to
	 *  his/her fitness value stored in \e selectionField.
	 */
	randomParentsChooser(bool replacement = true, const string & selectionField = "fitness") :
		parentChooser(selectionField), m_replacement(replacement),
		m_maleIndex(0), m_femaleIndex(0), m_maleFitness(0), m_femaleFitness(0),
		m_malesampler(rng()), m_femalesampler(rng())
	{
	}


	/// Deep copy of a random parents chooser.
	parentChooser * clone() const
	{
		return new randomParentsChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	/// CPPONLY
	individualPair chooseParents(RawIndIterator basePtr);

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


/** This parent chooser is similar to random parents chooser but instead of
 *  selecting a new pair of parents each time, one of the parents in this
 *  parent chooser will mate with several spouses before he/she is replaced.
 *  This mimicks multi-spouse mating schemes such as polygyny or polyandry
 *  in some populations. Natural selection is supported for both sexes.
 */
class polyParentsChooser : public parentChooser
{
public:
	/** Create a multi-spouse parents chooser where each father (if \e polySex
	 *  is Male) or mother (if \e polySex is Female) has \e polyNum spouses.
	 *  The parents are chosen with replacement. If natural selection
	 *  is enabled, the probability that an individual is chosen is
	 *  proportional to his/her fitness value among all individuals with the
	 *  same sex. Selection will be disabled if specified information field
	 *  \e selectionField (default to \c "fitness") does not exist.
	 */
	polyParentsChooser(Sex polySex = Male, UINT polyNum = 1,
		const string & selectionField = "fitness") :
		parentChooser(selectionField),
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


/** This parent chooser mimicks some animal populations where only certain
 *  individuals (usually males) can mate. Alpha individuals can be chosen
 *  either randomly (with natural selection) or according to an information
 *  field. After the alpha individuals are selected, the parent chooser works
 *  identical to a random mating scheme, except that one of the parents are
 *  chosen from these alpha individuals.
 */
class alphaParentsChooser : public parentChooser
{
public:
	/** Create a parent chooser that chooses father (if \e alphaSex is \c Male)
	 *  or mother (if \e alphaSex is \c Female) from a selected group of alpha
	 *  individuals. If \e alphaNum is given, alpha individuals are chosen
	 *  randomly or according to individual fitness if natural selection is
	 *  enabled. If \e alphaField is given, individuals with non-zero values
	 *  at this information field are considered as alpha individuals. After
	 *  alpha individuals are selected, \e alphaSex parent will be chosen
	 *  from the alpha individuals randomly or according to individual fitness.
	 *  The other parents are chosen randomly.
	 */
	alphaParentsChooser(Sex alphaSex = Male, UINT alphaNum = 0, string alphaField = string(),
		const string & selectionField = "fitness") :
		parentChooser(selectionField),
		m_alphaSex(alphaSex), m_alphaNum(alphaNum), m_alphaField(alphaField),
		m_maleIndex(0), m_femaleIndex(0),
		m_maleFitness(0), m_femaleFitness(0),
		m_malesampler(rng()), m_femalesampler(rng())
	{
	}


	/// Deep copy of an alpha parents chooser.
	parentChooser * clone() const
	{
		return new alphaParentsChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	/// CPPONLY
	individualPair chooseParents(RawIndIterator basePtr);

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


/** This parent chooser chooses an individual randomly, and then his/her spouse
 *  his/her spouse from a given set of information fields, which stores indexes
 *  of individuals in the same generation. An information field will be ignored
 *  if its value is negative, or if sex is incompatible.
 *
 *  Depending on what indexes are stored in these information fields, this
 *  parent chooser can be used to implement different types of mating schemes
 *  where selection of spouse is limited. For example, a consanguineous mating
 *  scheme can be implemeneted using this mating scheme if certain type of
 *  relatives are located for each individual, and are used for mating.
 *
 *  This parent chooser uses \c randomParentChooser to choose one parent and
 *  randomly choose another one from the information fields. Natural selection
 *  is supported during the selection of the first parent. Because of
 *  potentially uneven distribution of valid information fields, the overall
 *  process may not be as random as expected.
 */
class infoParentsChooser : public randomParentChooser
{
public:
	/** Create a information parent chooser a parent randomly (with replacement,
	 *  and with selection if natural selection is enabled), and then his/her
	 *  spouse from indexes stored in \e infoFields. If a Python function
	 *  \e func is specified, it will be called before parents are chosen. This
	 *  function accepts the parental population and an optional parameter
	 *  \e param and is usually used to locate qualified spouse for each parent.
	 *  The return value of this function is ignored.
	 */
	infoParentsChooser(const vectorstr & infoFields = vectorstr(),
		PyObject * func = NULL, PyObject * param = NULL,
		const string & selectionField = "fitness") :
		randomParentChooser(true, selectionField),
		m_infoFields(infoFields), m_func(func), m_param(param),
		m_infoIdx(0), m_degenerate(false)
	{
		DBG_FAILIF(m_infoFields.empty(), ValueError,
			"At least one information field should be provided for this infoParentsChooser");
	}


	/// Deep copy of a infomation parent chooser.
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


/** This parents chooser accept a Python generator function that repeatedly
 *  yields an index (relative to each subpopulation) of a parent, or indexes
 *  of two parents as a Python list of tuple. The parent chooser calls the
 *  generator function with parental population and a subpopulation index
 *  for each subpopulation and retrieves indexes of parents repeatedly using
 *  the iterator interface of the generator function.
 *
 *  This parent chooser does not support virtual subpopulation directly. A
 *  \c ValueError will be raised if this parent chooser is applied to a
 *  virtual subpopulation. However, because virtual subpopulations are defined
 *  in the passed parental population, it is easy to return parents from a
 *  particular virtual subpopulation using virtual subpopulation related
 *  functions.
 */
class pyParentsChooser : public parentChooser
{
public:
	/** Create a Python parent chooser using a Python generator function
	 *  \e parentsGenerator. This function should accept a population object
	 *  (the parental population) and a subpopulation number and return
	 *  the index of a parent or a pair of parents repeatedly using the
	 *  iterator interface of the generator function.
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


	/// Deep copy of a python parent chooser.
	parentChooser * clone() const
	{
		return new pyParentsChooser(*this);
	}


	/// CPPONLY
	void initialize(population & pop, SubPopID sp);

	/// CPPONLY
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


/** This class is the base class of all mating schemes. It should not be used
 * directly.
 * HIDDEN
 */
class mating
{

public:
	/** CPPONLY
	 *  check if the mating type is compatible with the population structure.
	 */
	virtual bool isCompatible(const population & pop) const
	{
		return true;
	}


	/** Create a mating scheme. \e subPopSize can be used to determine
	 * subpopulatio sizes of an offspring generation.
	 */
	mating(const uintListFunc & subPopSize = uintListFunc());

	/// destructor
	virtual ~mating()
	{
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
	 */
	virtual void submitScratch(population & pop, population & scratch);

	/** CPPONLY
	 *  mate a subpopulation, called by mate().
	 */
	virtual bool mateSubPop(population & pop, SubPopID subPop,
	                        RawIndIterator offBegin, RawIndIterator offEnd,
	                        vector<baseOperator * > & ops)
	{
		return true;
	}


	/** CPPONLY
	 *  Generate an offspring population \e scratch from parental population
	 *  \e pop. During mating operators \e ops will be applied after each
	 *  offspring is generated.
	 */
	virtual bool mate(population & pop, population & scratch, vector<baseOperator * > & ops);

	/** CPPONLY
	 *  Prepare a scratch population \e scratch.
	 */
	bool prepareScratchPop(population & pop, population & scratch);

protected:
	/** Specify subpopulation size of the offspring generation. Can be a
	 *  list of subpopulation sizes or a function.
	 */
	uintListFunc m_subPopSize;
};

/** A pedigree mating scheme that evolves a population following a
 *  pedigree object.
 */
class pedigreeMating : public mating
{
public:
	/** Creates a mating scheme that evolve a population following a pedigree
	 *  object \e ped. Considering this pedigree as a population with \c N
	 *  ancestral generations, the starting population is the greatest ancestral
	 *  generation of \e ped. The mating scheme creates an offspring generation
	 *  that match the size of generation \c N-1 and chooses parents according
	 *  to the parents of individuals at this generation. Depending on the \e gen
	 *  parameter of the simulator, the process continues generation by
	 *  generation for \c N generations if \c gen >= N), or \c gen generations
	 *  if \c gen < \c N. During the evolution, an offspring generator
	 *  \e generator is used to produce one offspring at a time, regardless of
	 *  the \e numOffspring setting of this offspring generator. If individuals
	 *  in pedigree \e ped has only one parent, the offspring generator should
	 *  be compatible.
	 *
	 *  By default, the pedigree mating scheme does not set offspring sex and
	 *  affection status using sex and affection status of corresponding
	 *  individuals in the pedigree. However, if such information is valid
	 *  in the pedigree object \e ped, you can set parameters \e setSex and/or
	 *  \e setAffection to \c True to set sex and/of affection status to
	 *  offspring during the evolutionary process. Similarly, you can specify
	 *  some information fields in \e copyFields to copy some information
	 *  fields from pedigree to the evolving population. Note that these
	 *  information will be copied also to the starting population (from the
	 *  greatest ancestral generation in \e ped).
	 */
	pedigreeMating(const pedigree & ped, const offspringGenerator & generator,
		bool setSex = false, bool setAffection = false,
		const vectorstr & copyFields = vectorstr());

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


	/// CPPONLY
	bool prepareScratchPop(population & pop, population & scratch);

	/// CPPONLY
	virtual bool mate(population & pop, population & scratch, vector<baseOperator * > & ops);

private:
	pedigree m_ped;

	offspringGenerator * m_generator;

	int m_parentalPopSize;

	bool m_setSex;
	bool m_setAffection;

	vectorstr m_copyFields;
};


/** A homogeneous mating scheme that uses a parent chooser to choose parents
 *  from a prental generation, and an offspring generator to generate offspring
 *  from chosen parents. It can be either used directly, or within a
 *  heterogeneous mating scheme. In the latter case, it can be applied to a
 *  (virtual) subpopulation.
 */
class homoMating : public mating
{
public:
	/** Create a homogeneous mating scheme using a parent chooser \e chooser
	 *  and an offspring generator \e generator.
	 *
	 *  If this mating scheme is used directly in a simulator, it will be
	 *  responsible for creating an offspring population according to parameter
	 *  \e subPopSize. This parameter can be a list of subpopulation sizes
	 *  (or a number if there is only one subpopulation) or a Python function.
	 *  The function should take two parameters, a generation number and a
	 *  list of subpopulation sizes before mating, and return a list of
	 *  subpopulation sizes for the offspring generation. A single number can
	 *  be returned if there is only one subpopulation. If latter form is used,
	 *  the specified function will be called at each generation to determine
	 *  the size of the offspring generation. Parameters \e subPop and
	 *  \e weight are ignored in this case.
	 *
	 *  If this mating shcme is used within a heterogeneous mating scheme.
	 *  Parameters \e subPop and \e weight are used to determine which (virtual)
	 *  subpopulation this mating scheme will be applied to, and how many
	 *  offspring this mating scheme will produce. Please refer to mating scheme
	 *  \c heteroMating for the use of these two parameters.
	 */
	homoMating(parentChooser & chooser,
		offspringGenerator & generator,
		const uintListFunc & subPopSize = uintListFunc(),
		vspID subPop = vspID(),
		double weight = 0);

	/// destructor
	~homoMating()
	{
		delete m_parentChooser;
		delete m_offspringGenerator;
	}


	/// CPPONLY
	homoMating(const homoMating & rhs) :
		mating(rhs), m_subPop(rhs.m_subPop), m_weight(rhs.m_weight)
	{
		m_offspringGenerator = rhs.m_offspringGenerator->clone();
		m_parentChooser = rhs.m_parentChooser->clone();
	}


	/// Deep copy of a homogeneous mating scheme
	virtual mating * clone() const
	{
		return new homoMating(*this);
	}


	/// used by Python print function to print out the general information of the Python mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::homoMating>";
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

	///
	vspID m_subPop;

	///
	double m_weight;

};

typedef std::vector<homoMating *> vectormating;

/** A heterogeneous mating scheme that applies a list of mating
 *  schemes to different (virtual) subpopulations.
 */
class heteroMating : public mating
{
public:
	/** Create a heterogeneous mating scheme that will apply a list of
	 *  homogeneous mating schemes \e matingSchemes to different (virtual)
	 *  subpopulations. The size of the offspring generation is determined
	 *  by parameter \e subPopSize, which can be a list of subpopulation sizes
	 *  or a Python function that returns a list of subpopulation sizes at
	 *  each generation. Please refer to \c homoMating for a detailed
	 *  explanation of this parameter.
	 *
	 *  Each mating scheme defined in \e matingSchemes can be applied to
	 *  one or more (virtual) subpopulation. If parameter \e subPop is not
	 *  specified, a mating scheme will be applied to all subpopulations.
	 *  If a (virtual) subpopulation is specified, a mating scheme will be
	 *  applied to a specific (virtual) subpopulation. A special case is
	 *  when \e subPop is given as <tt>(-1, vsp)</tt>. In this case, the
	 *  mating scheme will be applied to virtual subpopulation \c vsp in
	 *  all subpopulations.
	 *
	 *  If multiple mating schemes are applied to the same subpopulation, a
	 *  weight (parameter \e weight) can be given to each mating scheme to
	 *  determine how many offspring it will produce. The default \weight for
	 *  all mating schemes are \c 0. In this case, the number of offspring each
	 *  mating scheme produces is proportional to the size of its parental
	 *  (virtual) subpopulation. If all weights are negative, the numbers of
	 *  offspring are determined by the multiplication of the absolute values
	 *  of the weights and their respective parental (virtual) subpopulation
	 *  sizes. If all weights are positive, the number of offspring produced by
	 *  each mating scheme is proportional to these weights. Mating schemes
	 *  with zero weight in this case will produce no offspring. If both
	 *  negative and positive weights are present, negative weights are
	 *  processed before positive ones.
	 *
	 *  If multiple mating schemes are applied to the same subpopulation,
	 *  offspring produced by these mating schemes are shuffled randomly. If this
	 *  is not desired, you can turn off offspring shuffling by setting parameter
	 *  \e shuffleOffspring to \c False.
	 */
	heteroMating(const vectormating & matingSchemes,
		const uintListFunc & subPopSize = uintListFunc(),
		bool shuffleOffspring = true);

	/// destructor
	~heteroMating();

	/// CPPONLY
	heteroMating(const heteroMating & rhs);

	/// deep copy of a heterogeneous mating scheme
	virtual mating * clone() const
	{
		return new heteroMating(*this);
	}


	/// used by Python print function to print out the general information of the Python mating scheme
	virtual string __repr__()
	{
		return "<simuPOP::heteroMating>";
	}


	/** CPPONLY Call each homogeneous mating scheme to populate offspring
	 *  generation.
	 */
	bool mate(population & pop, population & scratch, vector<baseOperator * > & ops);

private:
	vectormating m_matingSchemes;
	///
	bool m_shuffleOffspring;
};


}
#endif
