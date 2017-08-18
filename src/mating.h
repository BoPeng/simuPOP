/**
 *  $File: mating.h $
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

#ifndef _MATING_H
#define _MATING_H
/**
   \file
   \brief head file of class mating and its subclasses
 */
#include "utility.h"
// for trajectory simulation functions
#include "simuPOP_cfg.h"
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

/// CPPONLY
class SexModel
{
public:
	SexModel() {}
	virtual ~SexModel() {}
	virtual Sex getSex(UINT count) = 0;

	virtual void reset() {}
	virtual SexModel * clone() = 0;

	virtual bool parallelizable() const
	{
		return false;
	}


};

/// CPPONLY
class NoSexModel : public SexModel
{
public:
	NoSexModel() {}

	SexModel * clone()
	{
		return new NoSexModel(*this);
	}


	Sex getSex(UINT)
	{
		return MALE;
	}


	bool parallelizable() const
	{
		return true;
	}


};

/// CPPONLY
class RandomSexModel : public SexModel
{
public:
	RandomSexModel() {}

	SexModel * clone()
	{
		return new RandomSexModel(*this);
	}


	Sex getSex(UINT)
	{
		return getRNG().randBit() ? MALE : FEMALE;
	}


	bool parallelizable() const
	{
		return true;
	}


};


/// CPPONLY
class ProbOfMalesSexModel : public SexModel
{
public:
	ProbOfMalesSexModel(double probOfMales) : m_probOfMales(probOfMales)
	{
	}


	SexModel * clone()
	{
		return new ProbOfMalesSexModel(*this);
	}


	Sex getSex(UINT)
	{
		return getRNG().randUniform() < m_probOfMales ? MALE : FEMALE;
	}


	bool parallelizable() const
	{
		return true;
	}


private:
	double m_probOfMales;
};


/// CPPONLY
class NumOfMalesSexModel : public SexModel
{
public:
	NumOfMalesSexModel(UINT numOfMales) : m_numOfMales(numOfMales)
	{
	}


	SexModel * clone()
	{
		return new NumOfMalesSexModel(*this);
	}


	Sex getSex(UINT count)
	{
		return count < m_numOfMales ? MALE : FEMALE;
	}


	bool parallelizable() const
	{
		return true;
	}


private:
	UINT m_numOfMales;
};

/// CPPONLY
class NumOfFemalesSexModel : public SexModel
{
public:
	NumOfFemalesSexModel(UINT numOfFemales) : m_numOfFemales(numOfFemales)
	{
	}


	SexModel * clone()
	{
		return new NumOfFemalesSexModel(*this);
	}


	Sex getSex(UINT count)
	{
		return count < m_numOfFemales ? FEMALE : MALE;
	}


	bool parallelizable() const
	{
		return true;
	}


private:
	UINT m_numOfFemales;
};


/// CPPONLY
class SeqSexModel : public SexModel
{
public:
	SeqSexModel(const vectorf & sex);

	SexModel * clone()
	{
		return new SeqSexModel(*this);
	}


	Sex getSex(UINT count)
	{
		return m_sex[count % m_sex.size()];
	}


	bool parallelizable() const
	{
		return true;
	}


private:
	vector<Sex> m_sex;
};


/// CPPONLY
class GlobalSeqSexModel : public SexModel
{
public:
	GlobalSeqSexModel(const vectorf & sex);

	SexModel * clone()
	{
		return new GlobalSeqSexModel(*this);
	}


	Sex getSex(UINT)
	{
		return m_sex[m_index++ % m_sex.size()];
	}


	void reset()
	{
		m_index = 0;
	}


	bool parallelizable() const
	{
		return false;
	}


private:
	vector<Sex> m_sex;
	int m_index;
};


/// CPPONLY
class FuncSexModel : public SexModel
{
public:
	FuncSexModel(const pyFunc & func) : m_func(func), m_generator(NULL)
	{
	}


	SexModel * clone()
	{
		return new FuncSexModel(*this);
	}


	Sex getSex(UINT count);


	void reset()
	{
		m_generator.set(NULL);
	}


	bool parallelizable() const
	{
		return false;
	}


private:
	pyFunc m_func;
	pyGenerator m_generator;
};

/// CPPONLY
class NumOffModel
{
public:
	NumOffModel() {}
	virtual ~NumOffModel() {}
	virtual UINT getNumOff(ssize_t gen) = 0;

	virtual void reset() {}
	virtual NumOffModel * clone() = 0;

	virtual bool parallelizable() const
	{
		return false;
	}


};

/// CPPONLY
class ConstNumOffModel : public NumOffModel
{
public:
	ConstNumOffModel(UINT numOff) : m_numOff(numOff)
	{
	}


	NumOffModel * clone()
	{
		return new ConstNumOffModel(*this);
	}


	UINT getNumOff(ssize_t /* gen */)
	{
		return m_numOff;
	}


	bool parallelizable() const
	{
		return true;
	}


private:
	UINT m_numOff;
};


/// CPPONLY
class GeometricNumOffModel : public NumOffModel
{
public:
	GeometricNumOffModel(double p) : m_p(p)
	{
	}


	NumOffModel * clone()
	{
		return new GeometricNumOffModel(*this);
	}


	UINT getNumOff(ssize_t /* gen */)
	{
		// possible 64 bit -> 32 bit conversion
		return static_cast<UINT>(getRNG().randGeometric(m_p));
	}


private:
	double m_p;
};


/// CPPONLY
class PoissonNumOffModel : public NumOffModel
{
public:
	PoissonNumOffModel(double mu) : m_mu(mu)
	{
	}


	NumOffModel * clone()
	{
		return new PoissonNumOffModel(*this);
	}


	UINT getNumOff(ssize_t /* gen */)
	{
		// possible 64 bit -> 32 bit conversion
		return static_cast<UINT>(getRNG().randTruncatedPoisson(m_mu));
	}


private:
	double m_mu;
};


/// CPPONLY
class BinomialNumOffModel : public NumOffModel
{
public:
	BinomialNumOffModel(ULONG N, double mu) : m_N(N), m_mu(mu)
	{
	}


	NumOffModel * clone()
	{
		return new BinomialNumOffModel(*this);
	}


	UINT getNumOff(ssize_t /* gen */)
	{
		return static_cast<UINT>(getRNG().randTruncatedBinomial(static_cast<ULONG>(m_N), m_mu));
	}


private:
	size_t m_N;
	double m_mu;
};


/// CPPONLY
class UniformNumOffModel : public NumOffModel
{
public:
	UniformNumOffModel(UINT low, UINT high) : m_low(low), m_high(high)
	{
	}


	NumOffModel * clone()
	{
		return new UniformNumOffModel(*this);
	}


	UINT getNumOff(ssize_t /* gen */)
	{
		// max: 5
		// num: 2
		// randint(4)  ==> 0, 1, 2, 3
		// + 2 ==> 2, 3, 4, 5
		// returning 0 is meaningless so we shift low from 0 to 1. This is actually a
		// truncated uniform distribution with U(X) = P(X) / (1-P(0)) = 1/n / (1-1/n) = 1/(n-1)
		if (m_low == 0)
			m_low = 1;
		return getRNG().randInt(m_high - m_low + 1) + m_low;
	}


private:
	UINT m_low;
	UINT m_high;
};


/// CPPONLY
class FuncNumOffModel : public NumOffModel
{
public:
	FuncNumOffModel(const pyFunc & func) : m_func(func), m_generator(NULL)
	{
	}


	NumOffModel * clone()
	{
		return new FuncNumOffModel(*this);
	}


	UINT getNumOff(ssize_t gen);

	void reset()
	{
		m_generator.set(NULL);
	}


private:
	pyFunc m_func;
	pyGenerator m_generator;
};


/** An <em>offspring generator</em> generates offspring from parents chosen by
 *  a parent chooser. It is responsible for creating a certain number of
 *  offspring, determinning their sex, and transmitting genotypes from parents
 *  to offspring.
 */
class OffspringGenerator
{
public:
	/** Create a basic offspring generator. This offspring generator uses
	 *  \e ops genotype transmitters to transmit genotypes from parents to
	 *  offspring.
	 *
	 *  A number of <em>during-mating operators</em> (parameter \e ops) can be
	 *  used to, among other possible duties such as setting information fields
	 *  of offspring, transmit genotype from parents to offspring. This general
	 *  offspring generator does not have any default during-mating operator
	 *  but all stock mating schemes use an offspring generator with a default
	 *  operator. For example, a \c mendelianOffspringGenerator is used by
	 *  \c RandomMating to trasmit genotypes. Note that applicability
	 *  parameters \c begin, \c step, \c end, \c at and \c reps could be used
	 *  in these operators but negative population and generation indexes are
	 *  unsupported.
	 *
	 *  Parameter \e numOffspring is used to control the number of offspring
	 *  per mating event, or in another word the number of offspring in each
	 *  family. It can be a number, a Python function or generator, or a mode
	 *  parameter followed by some optional arguments. If a number is given,
	 *  given number of offspring will be generated at each mating event. If a
	 *  Python function is given, it will be called each time when a mating
	 *  event happens. When a generator function is specified, it will be
	 *  called for each subpopulation to provide number of offspring for all
	 *  mating events during the populating of this subpopulation. Current
	 *  generation number will be passed to this function or generator function
	 *  if parameter "gen" is used in this function. In the last case, a tuple
	 *  (or a list) in one of the following forms can be given:
	 *  \li <tt>(GEOMETRIC_DISTRIBUTION, p)</tt>
	 *  \li <tt>(POISSON_DISTRIBUTION, p)</tt>, p > 0
	 *  \li <tt>(BINOMIAL_DISTRIBUTION, p, N)</tt>, 0 < p <=1, N > 0
	 *  \li <tt>(UNIFORM_DISTRIBUTION, a, b)</tt>, 0 <= a <= b.
	 *
	 *  In this case, the number of offspring will be determined randomly
	 *  following the specified statistical distributions. Because families
	 *  with zero offspring are silently ignored, the distribution of the
	 *  observed number of offspring per mating event (excluding zero)
	 *  follows zero-truncated versions of these distributions.
	 *
	 *  Parameter \e numOffspring specifies the number of offspring per mating
	 *  event but the actual surviving offspring can be less than specified.
	 *  More spefically, if any during-mating operator returns \c False, an
	 *  offspring will be discarded so the actually number of offspring of a
	 *  mating event will be reduced. This is essentially how during-mating
	 *  selector works.
	 *
	 *  Parameter \e sexMode is used to control the sex of each offspring. Its
	 *  default value is usually \e RANDOM_SEX which assign \c MALE or \c FEMALE
	 *  to each individual randomly, with equal probabilities. If \c NO_SEX is
	 *  given, offspring sex will not be changed. \e sexMode can also be one of
	 *  \li <tt>(PROB_OF_MALES, p)</tt> where \c p is the probability of male
	 *       for each offspring,
	 *  \li <tt>(NUM_OF_MALES, n)</tt> where \c n is the number of males in a
	 *       mating event. If \c n is greater than or equal to the number of
	 *       offspring in this family, all offspring in this family will be
	 *       \c MALE.
	 *  \li <tt>(NUM_OF_FEMALES, n)</tt> where \c n is the number of females in
	 *       a mating event,
	 *  \li <tt>(SEQUENCE_OF_SEX, s1, s2 ...)</tt> where \c s1, \c s2 etc are
	 *       MALE or FEMALE. The sequence will be used for each mating event.
	 *       It will be reused if the number of offspring in a mating event
	 *       is greater than the length of sequence.
	 *  \li <tt>(GLOBAL_SEQUENCE_OF_SEX, s1, s2, ...)</tt> where \c s1, \c s2
	 *       etc are MALE or FEMALE. The sequence will be used across mating
	 *       events. It will be reused if the number of offspring in a
	 *       subpopulation is greater than the length of sequence.
	 *
	 *  Finally, parameter \e sexMode accepts a function or a generator function.
	 *  A function will be called whenever an offspring is produced. A generator
	 *  will be created at each subpopulation and will be used to produce sex
	 *  for all offspring in this subpopulation. No parameter is accepted.
	 */
	OffspringGenerator(const opList & ops, const floatListFunc & numOffspring = 1,
		const floatListFunc & sexMode = RANDOM_SEX);

	/// CPPONLY
	OffspringGenerator(const OffspringGenerator & rhs)
		: m_transmitters(rhs.m_transmitters), m_initialized(rhs.m_initialized)
	{
		m_numOffModel = rhs.m_numOffModel->clone();
		m_sexModel = rhs.m_sexModel->clone();
	}


	virtual ~OffspringGenerator()
	{
		delete m_numOffModel;
		delete m_sexModel;
	}


	/// HIDDEN Make a deep copy of this offspring generator.
	virtual OffspringGenerator * clone() const
	{
		return new OffspringGenerator(*this);
	}


	/** create an offspring generator, save information from \c pop and \c ops to
	 *  speed up the calls to \c generateOffspring
	 *  CPPONLY
	 */
	virtual void initialize(const Population & pop, size_t subPop);

	/// CPPONLY
	virtual UINT generateOffspring(Population & pop, Population & offPop, Individual * dad, Individual * mom,
		RawIndIterator & offBegin, RawIndIterator & offEnd);

	/// CPPONLY
	virtual void finalize(const Population & /* pop */)
	{
		m_numOffModel->reset();
		m_sexModel->reset();
		m_initialized = false;
	}


	/// HIDDEN describe an offspring generator
	virtual string describe(bool format = true) const;

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
	UINT numOffspring(ssize_t gen);


	/** CPPONLY
	 *  return sex according to m_sexParam, m_sexMode and
	 *  \e count, which is the index of offspring
	 */
	Sex getSex(UINT count);


	/// CPPONLY
	bool parallelizable() const;

protected:
	/// number of offspring
	NumOffModel * m_numOffModel;

	SexModel * m_sexModel;

	/// default transmitter
	opList m_transmitters;

protected:
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
 *  "Peng et al, (2007) PLoS Genetics".
 */
class ControlledOffspringGenerator : public OffspringGenerator
{
public:
	/** Create an offspring generator that selects offspring so that allele
	 *  frequency at specified loci in the offspring generation reaches
	 *  specified allele frequency. At the beginning of each generation,
	 *  expected allele frequency of \e alleles at \e loci is returned from a
	 *  user-defined trajectory function \e freqFunc. Parameter \e loci can
	 *  be a list of loci indexes, names, or ALL_AVAIL. If there is no
	 *  subpopulation, this function should return a list of frequencies for
	 *  each locus. If there are multiple subpopulations, \e freqFunc can
	 *  return a list of allele frequencies for all subpopulations or
	 *  combined frequencies that ignore population structure. In the former
	 *  case, allele frequencies should be arranged by loc0_sp0, loc1_sp0, ...
	 *  loc0_sp1, loc1_sp1, ..., and so on. In the latter case, overall expected
	 *  number of alleles are scattered to each subpopulation in proportion to
	 *  existing number of alleles in each subpopulation, using a multinomial
	 *  distribution.
	 *
	 *  After the expected alleles are calculated, this offspring generator
	 *  accept and reject families according to their genotype at \e loci until
	 *  allele frequecies reach their expected values. The rest of the
	 *  offspring generation is then filled with families without only wild
	 *  type alleles at these \e loci.
	 *
	 *  This offspring generator is derived from class \e OffspringGenerator.
	 *  Please refer to class \e OffspringGenerator for a detailed description
	 *  of parameters \e ops, \e numOffspring and \e sexMode.
	 */
	ControlledOffspringGenerator(const lociList & loci, const uintList & alleles,
		PyObject * freqFunc, const opList & ops = vectorop(),
		const floatListFunc & numOffspring = 1, const floatListFunc & sexMode = RANDOM_SEX);


	/// CPPONLY
	ControlledOffspringGenerator(const ControlledOffspringGenerator & rhs);

	/// CPPONLY
	void initialize(const Population & pop, size_t subPop);

	/// CPPONLY
	virtual UINT generateOffspring(Population & pop, Population & offPop, Individual * dad, Individual * mom,
		RawIndIterator & offBegin,
		RawIndIterator & offEnd);

	/// HIDDEN Deep copy of a controlled random mating scheme
	virtual OffspringGenerator * clone() const
	{
		return new ControlledOffspringGenerator(*this);
	}


	/// HIDDEN describe a controlled offspring generator
	virtual string describe(bool format = true) const;

private:
	void getExpectedAlleles(const Population & pop, vectorf & expFreq);

	/// locus at which mating is controlled.
	lociList m_loci;
	//
	/// allele to be controlled at each locus
	vectoru m_alleles;

	/// function that return an array of frquency range
	pyFunc m_freqFunc;
	//
	// expected alleles
	vectoru m_expAlleles;
	vectoru m_flip;                 // in a subpop
	vectoru m_totAllele;            // in a subpop
	vectoru m_curAllele;            // in a subpop
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
class ParentChooser
{
public:
	typedef std::pair<Individual *, Individual *> IndividualPair;

public:
	// CPPONLY
	ParentChooser(const string & selectionField = string())
		: m_initialized(false), m_selectionField(selectionField)
	{
	}


	/// HIDDEN Deep copy of a parent chooser
	virtual ParentChooser * clone() const
	{
		return new ParentChooser(*this);
	}


	/// Initialize a parent chooser for subpopulation \e subPop of \e population pop
	virtual void initialize(Population & pop, size_t subPop)
	{
	}


	/// CPPONLY
	virtual void finalize()
	{
		m_initialized = false;
	}


	/// HIDDEN describe a general parent chooser
	virtual string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.ParentChooser> (base class)";
	}


	/// CPPONLY
	virtual bool parallelizable() const
	{
		return false;
	}


	/// CPPONLY
	bool initialized() const
	{
		return m_initialized;
	}


	/// Return chosen parents from a population if the parent chooser object is created with a population
	virtual IndividualPair chooseParents()
	{
		// do not use NULL because VC2010 uses nullptr and who knows what gcc uses
		return IndividualPair((Individual *)(0), (Individual *)(0));
	}


	/// destructor
	virtual ~ParentChooser()
	{
		finalize();
	}


protected:
	bool m_initialized;

	string m_selectionField;
};


/** This parent chooser chooses a parent from a parental (virtual) subpopulation
 *  sequentially. Natural selection is not considered. If the last parent is
 *  reached, this parent chooser will restart from the beginning of the
 *  (virtual) subpopulation.
 */
class SequentialParentChooser : public ParentChooser
{
public:
	/** Create a parent chooser that chooses a parent from a parental (virtual)
	 *  subpopulation sequentially. Parameter \e choice can be \c ANY_SEX
	 *  (default), \c MALE_ONLY and \c FEMALE_ONLY. In the latter two cases,
	 *  only male or female individuals are selected. A \c RuntimeError will be
	 *  raised if there is no male or female individual from the population.
	 */
	SequentialParentChooser(SexChoice sexChoice = ANY_SEX)
		: ParentChooser(), m_choice(sexChoice)
	{
		DBG_FAILIF(m_choice != ANY_SEX && m_choice != MALE_ONLY && m_choice != FEMALE_ONLY,
			ValueError, "choice of sex can only be ANY_SEX, MALE_ONLY or FEMALE_ONLY.");
	}


	/// HIDDEN Deep copy of a sequential parent chooser.
	ParentChooser * clone() const
	{
		return new SequentialParentChooser(*this);
	}


	/// HIDDEN describe a sequential parent chooser
	virtual string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.SequentialParentChooser> chooses a parent sequentially";
	}


	/// Initialize a parent chooser for subpopulation \e subPop of \e population pop
	void initialize(Population & pop, size_t subPop);

	/// Return chosen parents from a population if the parent chooser object is created with a population
	IndividualPair chooseParents();

private:
	///
	SexChoice m_choice;
	/// starting individual
	IndIterator m_begin;
	/// current individual
	IndIterator m_ind;
	/// for sexual selection of parents
	vector<RawIndIterator> m_index;
	size_t m_curInd;
};


/** This parent chooser chooses a parent randomly from a (virtual) parental
 *  subpopulation. Parents are chosen with or without replacement.
 *  If parents are chosen with replacement, a parent can be selected multiple
 *  times. If individual fitness values are assigned to individuals (
 *  stored in an information  field \e selectionField (default to \c "fitness"),
 *  individuals will be chosen at a probability proportional to his or her
 *  fitness value. If parents are chosen without replacement, a parent can be
 *  chosen only once. An \c RuntimeError will be raised if all parents are
 *  exhausted. Natural selection is disabled in the without-replacement case.
 */
class RandomParentChooser : public ParentChooser
{
public:
	/** Create a random parent chooser that choose parents with or without
	 *  replacement (parameter \e replacement, default to \c True). If selection
	 *  is enabled and information field \e selectionField exists in the passed
	 *  population, the probability that a parent is chosen is proportional to
	 *  his/her fitness value stored in \e selectionField. This parent chooser
	 *  by default chooses parent from all individuals (\c ANY_SEX), but it
	 *  can be made to select only male (\c MALE_ONLY) or female
	 *  (\c FEMALE_ONLY) individuals by setting parameter \e sexChoice.
	 */
	RandomParentChooser(bool replacement = true,
		const string & selectionField = "fitness", SexChoice sexChoice = ANY_SEX)
		: ParentChooser(selectionField), m_basePtr(), m_replacement(replacement), m_choice(sexChoice),
		m_index(0), m_chosen(0), m_sampler(), m_size(0), m_shift(0)
	{
	}


	/// HIDDEN Deep copy of a random parent chooser.
	ParentChooser * clone() const
	{
		return new RandomParentChooser(*this);
	}


	/// HIDDEN describe a random parent chooser
	virtual string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.RandomParentChooser> chooses one parent randomly";
	}


	/// CPPONLY
	virtual bool parallelizable() const
	{
		// if we choose parents with replacment, we can do that
		// in parallel. Otherwise, because a parent can only be
		// chosen once, we can not allow independent selection of parents.
		return m_replacement;
	}


	/// Initialize a parent chooser for subpopulation \e subPop of \e population pop
	void initialize(Population & pop, size_t subPop);

	/// Return chosen parents from a population if the parent chooser object is created with a population
	IndividualPair chooseParents();

protected:
	RawIndIterator m_basePtr;

	bool m_replacement;

	bool m_selection;
	///
	SexChoice m_choice;
	///
	vector<RawIndIterator> m_index;
	vector<RawIndIterator> m_chosen;
	/// accumulative fitness
	WeightedSampler m_sampler;
	/// individuals to choose
	size_t m_size;
	/// index to the subpopulation
	size_t m_shift;
};


/** This parent chooser chooses two parents, a male and a female, randomly from
 *  a (virtual) parental subpopulation. Parents are chosen with or without
 *  replacement from their respective sex group. If parents are chosen with
 *  replacement, a parent can be selected multiple times. If individual fitness
 *  values are assigned (stored in information field \e selectionField, default
 *  to \c "fitness", the probability that an individual is chosen is
 *  proportional to his/her fitness value among all individuals with the same
 *  sex. If parents are chosen without replacement, a parent can be chosen only
 *  once. An \c RuntimeError will be raised if all males or females are
 *  exhausted. Natural selection is disabled in the without-replacement case.
 */
class RandomParentsChooser : public ParentChooser
{
public:
	/** Create a random parents chooser that choose two parents with or without
	 *  replacement (parameter \e replacement, default to \c True). If selection
	 *  is enabled and information field \e selectionField exists in the passed
	 *  population, the probability that a parent is chosen is proportional to
	 *  his/her fitness value stored in \e selectionField.
	 */
	RandomParentsChooser(bool replacement = true, const string & selectionField = "fitness")
		: ParentChooser(selectionField), m_replacement(replacement),
		m_index(0), m_fitness(0), m_malesampler(), m_femalesampler()
	{
	}


	/// HIDDEN Deep copy of a random parents chooser.
	ParentChooser * clone() const
	{
		return new RandomParentsChooser(*this);
	}


	/// HIDDEN describe a random parents chooser
	virtual string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.RandomParentsChooser> chooses two parents randomly";
	}


	/// CPPONLY
	virtual bool parallelizable() const
	{
		// if we choose parents with replacment, we can do that
		// in parallel. Otherwise, because a parent can only be
		// chosen once, we can not allow independent selection of parents.
		return m_replacement;
	}


	/// Initialize a parent chooser for subpopulation \e subPop of \e population pop
	void initialize(Population & pop, size_t subPop);

	/// Return chosen parents from a population if the parent chooser object is created with a population
	IndividualPair chooseParents();

private:
	bool m_replacement;

	Individual * m_lastParent;

	bool m_selection;

	size_t m_numMale;
	size_t m_numFemale;

	/// internal index to female/males.
	vector<RawIndIterator> m_index;

	size_t m_femaleIndex;
	size_t m_maleIndex;

	vectorf m_fitness;

	// weighted sampler
	WeightedSampler m_malesampler;
	WeightedSampler m_femalesampler;
};


/** This parent chooser is similar to random parents chooser but instead of
 *  selecting a new pair of parents each time, one of the parents in this
 *  parent chooser will mate with several spouses before he/she is replaced.
 *  This mimicks multi-spouse mating schemes such as polygyny or polyandry
 *  in some populations. Natural selection is supported for both sexes.
 */
class PolyParentsChooser : public ParentChooser
{
public:
	/** Create a multi-spouse parents chooser where each father (if \e polySex
	 *  is MALE) or mother (if \e polySex is FEMALE) has \e polyNum spouses.
	 *  The parents are chosen with replacement. If individual fitness values
	 *  are assigned (stored to information field \c selectionField, default
	 *  to \c "fitness"), the probability that an individual is chosen is
	 *  proportional to his/her fitness value among all individuals with the
	 *  same sex.
	 */
	PolyParentsChooser(Sex polySex = MALE, UINT polyNum = 1,
		const string & selectionField = "fitness") :
		ParentChooser(selectionField),
		m_polySex(polySex), m_polyNum(polyNum),
		m_maleIndex(0), m_femaleIndex(0),
		m_chosenMale(0), m_chosenFemale(0),
		m_maleFitness(0), m_femaleFitness(0),
		m_malesampler(), m_femalesampler()
	{
#ifdef _OPENMP
		m_polyCount.resize(numThreads());
		m_lastParent.resize(numThreads());
		for (size_t i = 0; i < numThreads(); i++) {
			m_polyCount[i] = 0;
			m_lastParent[i] = NULL;
		}
#else
		m_polyCount = 0;
		m_lastParent = NULL;
#endif

		DBG_FAILIF(polyNum < 1, ValueError,
			"Number of sex partners has to be at least one");
	}


	ParentChooser * clone() const
	{
		return new PolyParentsChooser(*this);
	}


	/// HIDDEN describe a polygenic parents chooser
	virtual string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.PolyParentsChooser> chooses parents with several spouses";
	}


	/// CPPONLY
	virtual bool parallelizable() const
	{
		return true;
	}


	/// Initialize a parent chooser for subpopulation \e subPop of \e population pop
	void initialize(Population & pop, size_t subPop);

	/// Return chosen parents from a population if the parent chooser object is created with a population
	IndividualPair chooseParents();

private:
	Sex m_polySex;
	UINT m_polyNum;
#ifdef _OPENMP
	vector<UINT> m_polyCount;

	vector<Individual *> m_lastParent;
#else
	UINT m_polyCount;

	Individual * m_lastParent;
#endif

	bool m_selection;

	size_t m_numMale;
	size_t m_numFemale;

	/// internal index to female/males.
	vector<RawIndIterator> m_maleIndex;
	vector<RawIndIterator> m_femaleIndex;
	vector<RawIndIterator> m_chosenMale;
	vector<RawIndIterator> m_chosenFemale;

	vectorf m_maleFitness;
	vectorf m_femaleFitness;

	// weighted sampler
	WeightedSampler m_malesampler;
	WeightedSampler m_femalesampler;
};


/*  This parent chooser chooses an individual randomly, and then his/her spouse
 *  his/her spouse from a given set of information fields, which stores indexes
 *  or IDs (if an idField is given) of individuals in the parental population.
 *  If the parental population has multiple ancestral generations and IDs are
 *  used, this parents chooser allows cross-generational mating. An information
 *  field will be ignored if its value is negative, or if sex is incompatible.
 *
 *  Depending on what indexes are stored in these information fields, this
 *  parent chooser can be used to implement different types of mating schemes
 *  where selection of spouse is limited. For example, a consanguineous mating
 *  scheme can be implemeneted using this mating scheme if certain type of
 *  relatives are located for each individual, and are used for mating.
 *
 *  This parent chooser uses \c RandomParentChooser to choose one parent and
 *  randomly choose another one from the information fields. Natural selection
 *  is supported during the selection of the first parent. Because of
 *  potentially uneven distribution of valid information fields, the overall
 *  process may not be as random as expected.
 */
/*
   class infoParentsChooser : public RandomParentChooser
   {
   public:
 */
/*  Create a information parent chooser a parent randomly (with replacement,
 *  and with selection if natural selection is enabled), and then his/her
 *  spouse from indexes or IDs (if \e idFiels is not empty) stored in
 *  \e infoFields. If a Python function \e func is specified, it will be
 *  called before parents are chosen. This function accepts the parental
 *  population and an optional parameter \e param and is usually used to
 *  locate qualified spouse for each parent. The return value of this
 *  function is ignored.
 */
/*
   infoParentsChooser(const stringList & infoFields = vectorstr(),
   PyObject * func = NULL, PyObject * param = NULL,
   const string & idField = '',
   const string & selectionField = "fitness") :
   RandomParentChooser(true, selectionField),
   m_infoFields(infoFields.elems()), m_idField(idField), m_func(func), m_param(param),
   m_infoIdx(0), m_degenerate(false)
   {
   DBG_FAILIF(m_infoFields.empty(), ValueError,
       "At least one information field should be provided for this infoParentsChooser");
   }


   /// HIDDEN Deep copy of a infomation parent chooser.
   ParentChooser * clone() const
   {
   return new infoParentsChooser(*this);
   }


   /// HIDDEN describe a parents chooser using information fields
   virtual string describe(bool format = true) const
   {
   return "<simuPOP.infoParentsChooser> chooses parents from specified information fields";
   }


   /// CPPONLY
   void initialize(Population & pop, size_t sp);

    /// Return chosen parents
   IndividualPair chooseParents();

   private:
   vectorstr m_infoFields;
   string m_idField;

   pyFunc m_func;
   pyObject m_param;

   vectori m_infoIdx;
   // if there is no valid individual, this mating schemes
   // works like a double ParentChooser.
   bool m_degenerate;
   }; */


/** This parent chooser accepts two parent choosers. It takes one parent from
 *  each parent chooser and return them as father and mother. Because two
 *  parent choosers do not have to choose parents from the same virtual
 *  subpopulation, this parent chooser allows you to choose parents from
 *  different subpopulations.
 */
class CombinedParentsChooser : public ParentChooser
{
public:
	/** Create a Python parent chooser using two parent choosers
	 *  \e fatherChooser and \e motherChooser. It takes one parent from each
	 *  parent chooser and return them as father and mother. If two valid
	 *  parents are returned, the first valid parent (father) will be used
	 *  for \e fatherChooser, the second valid parent (mother) will be used
	 *  for \e motherChooser. Although these two parent choosers are supposed
	 *  to return a father and a mother respectively, the sex of returned
	 *  parents are not checked so it is possible to return parents with the
	 *  same sex using this parents chooser. This choose by default allows
	 *  the selection of the same parents as father and mother
	 *  (self-fertilization), unless a parameter \e allowSelfing is used to
	 *  disable it.
	 */
	CombinedParentsChooser(const ParentChooser & fatherChooser,
		const ParentChooser & motherChooser, bool allowSelfing=true);

	/// CPPONLY
	CombinedParentsChooser(const CombinedParentsChooser & rhs)
		: ParentChooser("fitness"), m_fatherChooser(rhs.m_fatherChooser->clone()),
		m_motherChooser(rhs.m_motherChooser->clone()),
		m_allowSelfing(rhs.m_allowSelfing)
	{
		m_initialized = false;
	}


	~CombinedParentsChooser()
	{
		delete m_fatherChooser;
		delete m_motherChooser;
	}


	/// HIDDEN Deep copy of a python parent chooser.
	ParentChooser * clone() const
	{
		return new CombinedParentsChooser(*this);
	}


	/// HIDDEN describe a hybrid parent chooser
	virtual string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.CombinedParentsChooser> chooses parents using two parent choosers";
	}


	/// CPPONLY
	virtual bool parallelizable() const
	{
		return m_fatherChooser->parallelizable() && m_motherChooser->parallelizable();
	}


	/// Initialize a parent chooser for subpopulation \e subPop of \e population pop
	void initialize(Population & pop, size_t subPop);

	/// CPPONLY
	void finalize();

	/// Return chosen parents from a population if the parent chooser object is created with a population
	IndividualPair chooseParents();

private:
	ParentChooser * m_fatherChooser;
	ParentChooser * m_motherChooser;
	bool m_allowSelfing;
};

/** This parent chooser accepts a Python generator function that repeatedly
 *  yields one or two parents, which can be references to individual objects
 *  or indexes relative to each subpopulation. The parent chooser calls the
 *  generator function with parental population and a subpopulation index
 *  for each subpopulation and retrieves parents repeatedly using the iterator
 *  interface of the generator function.
 *
 *  This parent chooser does not support virtual subpopulation directly.
 *  However, because virtual subpopulations are defined in the passed parental
 *  population, it is easy to return parents from a particular virtual
 *  subpopulation using virtual subpopulation related functions.
 */
class PyParentsChooser : public ParentChooser
{
public:
	/** Create a Python parent chooser using a Python generator function
	 *  \e parentsGenerator. This function should accept one or both of
	 *  parameters \e pop (the parental population) and \e subPop (index
	 *  of subpopulation) and return the reference or index (relative to
	 *  subpopulation) of a parent or a pair of parents repeatedly using
	 *  the iterator interface of the generator function.
	 */
	PyParentsChooser(PyObject * generator);

	/// CPPONLY
	PyParentsChooser(const PyParentsChooser & rhs)
		: ParentChooser(rhs), m_func(rhs.m_func),
		m_popObj(NULL), m_generator(NULL)
	{
		m_initialized = false;
	}


	/// HIDDEN Deep copy of a python parent chooser.
	ParentChooser * clone() const
	{
		return new PyParentsChooser(*this);
	}


	/// HIDDEN describe a hybrid parent chooser
	virtual string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.PyParentsChooser> chooses parents according to a user-provided Python function";
	}


	/** CPPONLY
	 *  PyParentsChooser can't be parallelizable because it call external functions
	 *  from Python which don't support multi-thread.
	 */
	virtual bool parallelizable() const
	{
		return false;
	}


	/// Initialize a parent chooser for subpopulation \e subPop of \e population pop
	void initialize(Population & pop, size_t subPop);

	/// CPPONLY
	void finalize();

	/// destructor
	~PyParentsChooser()
	{
		// this is a bad idea because PyParentsChooser will not be properly released when
		// some error happens, and this exception actually blocks some previous exception.
		// DBG_FAILIF(m_popObj != NULL, SystemError, "Python parents chooser is not properly destroyed.");
	}


	/// Return chosen parents from a population if the parent chooser object is created with a population
	IndividualPair chooseParents();

private:
#ifndef OPTIMIZED
	size_t m_size;
#endif
	IndIterator m_begin;

	pyFunc m_func;
	PyObject * m_popObj;
	pyGenerator m_generator;
};


/** This mating scheme is the base class of all mating schemes. It evolves
 *  a population generation by generation but does not actually transmit
 *  genotype.
 */
class MatingScheme
{
public:
	/** Create a base mating scheme that evolves a population without
	 *  transmitting genotypes. At each generation, this mating scheme
	 *  creates an offspring generation according to parameter \e subPopSize,
	 *  which can be a list of subpopulation sizes (or a number if there is
	 *  only one subpopulation) or a Python function which will be called at
	 *  each generation, just before mating, to determine the subpopulation
	 *  sizes of the offspring generation. The function should be defined
	 *  with one or both parameters of \c gen and \c pop where \c gen is the
	 *  current generation number and \c pop is the parental population just
	 *  before mating. The return value of this function should be a list of
	 *  subpopulation sizes for the offspring generation. A single number can
	 *  be returned if there is only one subpopulation. The passed parental
	 *  population is usually used to determine offspring population size from
	 *  parental population size but you can also modify this population to
	 *  prepare for mating. A common practice is to split and merge parental
	 *  populations in this function so that you demographic related
	 *  information and actions could be implemented in the same function.
	 */
	MatingScheme(const uintListFunc & subPopSize = uintListFunc());

	/// destructor
	virtual ~MatingScheme()
	{
	}


	/// HIDDEN Deep copy of a mating scheme
	virtual MatingScheme * clone() const
	{
		return new MatingScheme(*this);
	}


	/// HIDDEN describe a general mating scheme.
	virtual string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.mating> A mating scheme";
	}


	/** CPPONLY
	 */
	virtual void submitScratch(Population & pop, Population & scratch);

	/** CPPONLY
	 *  mate a subpopulation, called by mate().
	 */
	virtual bool mateSubPop(Population & /* pop */, Population & /* offPop */,
	                        size_t /* subPop */,
	                        RawIndIterator /* offBegin */, RawIndIterator /* offEnd */)
	{
		return true;
	}


	/** CPPONLY
	 *  Generate an offspring population \e scratch from parental population
	 *  \e pop.
	 */
	virtual bool mate(Population & pop, Population & scratch);

	/** CPPONLY
	 *  Prepare a scratch population \e scratch.
	 */
	bool prepareScratchPop(Population & pop, Population & scratch);

	/** CPPONLY
	 *  Use to generate a warning when subPopSize is specified in a homogeneous
	 *  mating scheme called in a heterogeneous mating scheme.
	 */
	bool subPopSizeSpecified()
	{
		return !m_subPopSize.empty() || m_subPopSize.func().isValid();
	}


protected:
	/** Specify subpopulation size of the offspring generation. Can be a
	 *  list of subpopulation sizes or a function.
	 */
	uintListFunc m_subPopSize;
};


/** A homogeneous mating scheme that uses a parent chooser to choose parents
 *  from a prental generation, and an offspring generator to generate offspring
 *  from chosen parents. It can be either used directly, or within a
 *  heterogeneous mating scheme. In the latter case, it can be applied to a
 *  (virtual) subpopulation.
 */
class HomoMating : public MatingScheme
{
public:
	/** Create a homogeneous mating scheme using a parent chooser \e chooser
	 *  and an offspring generator \e generator.
	 *
	 *  If this mating scheme is used directly in a simulator, it will be
	 *  responsible for creating an offspring population according to parameter
	 *  \e subPopSize. This parameter can be a list of subpopulation sizes (or
	 *  a number if there is only one subpopulation) or a Python function which
	 *  will be called at each generation to determine the subpopulation sizes
	 *  of the offspring generation. Please refer to class \c MatingScheme
	 *  for details about this parameter.
	 *
	 *  If this mating shcme is used within a heterogeneous mating scheme.
	 *  Parameters \e subPops and \e weight are used to determine which (virtual)
	 *  subpopulations this mating scheme will be applied to, and how many
	 *  offspring this mating scheme will produce. Please refer to mating scheme
	 *  \c HeteroMating for the use of these two parameters.
	 */
	HomoMating(ParentChooser & chooser,
		OffspringGenerator & generator,
		const uintListFunc & subPopSize = uintListFunc(),
		subPopList subPops = subPopList(),
		double weight = 0);

	/// destructor
	~HomoMating()
	{
		delete m_ParentChooser;
		delete m_OffspringGenerator;
	}


	/// CPPONLY
	HomoMating(const HomoMating & rhs) :
		MatingScheme(rhs), m_subPops(rhs.m_subPops), m_weight(rhs.m_weight)
	{
		m_OffspringGenerator = rhs.m_OffspringGenerator->clone();
		m_ParentChooser = rhs.m_ParentChooser->clone();
	}


	/// HIDDEN Deep copy of a homogeneous mating scheme
	virtual MatingScheme * clone() const
	{
		return new HomoMating(*this);
	}


	/// HIDDEN describe a homogeneous mating scheme.
	virtual string describe(bool format = true) const;


	/// CPPONLY
	subPopList subPops() const
	{
		return m_subPops;
	}


	/// CPPONLY
	double weight() const
	{
		return m_weight;
	}


	/// CPPONLY
	virtual bool mateSubPop(Population & pop, Population & offPop, size_t subPop,
		RawIndIterator offBegin, RawIndIterator offEnd);

private:
	ParentChooser * m_ParentChooser;
	OffspringGenerator * m_OffspringGenerator;

	///
	subPopList m_subPops;

	///
	double m_weight;

};


/** This mating scheme evolves a population following an existing pedigree
 *  structure. If the \c Pedigree object has \c N ancestral generations and a
 *  present generation, it can be used to evolve a population for \c N
 *  generations, starting from the topmost ancestral generation. At the \e k-th
 *  generation, this mating scheme produces an offspring generation according
 *  to subpopulation structure of the <tt>N-k-1</tt> ancestral generation in
 *  the pedigree object (e.g. producing the offspring population of generation
 *  0 according to the <tt>N-1</tt> ancestral generation of the pedigree object
 *  ). For each offspring, this mating scheme copies individual ID and sex from
 *  the corresponing individual in the pedigree object. It then locates the
 *  parents of each offspring using their IDs in the pedigree object. A list of
 *  during mating operators are then used to transmit parental genotype to
 *  the offspring. The population being evolved must have an information field
 *  \c 'ind_id'.
 */
class PedigreeMating : public MatingScheme
{
public:
	/** Creates a pedigree mating scheme that evolves a population according to
	 *  \c Pedigree object \e ped. The evolved population should contain
	 *  individuals with ID (at information field \e idField, default to
	 *  \c 'ind_id') that match those individual in the topmost ancestral
	 *  generation who have offspring. After parents of each individuals are
	 *  determined from their IDs, a list of during-mating operators
	 *  \e ops are applied to transmit genotypes. The return value of these
	 *  operators are not checked.
	 */
	PedigreeMating(const Pedigree & ped, const opList & ops,
		const string & idField = "ind_id") :
		m_ped(ped), m_transmitters(ops), m_idField(idField), m_gen(ped.ancestralGens() - 1)
	{
	}


	~PedigreeMating()
	{
	}


	/// HIDDEN Deep copy of a Python mating scheme
	virtual MatingScheme * clone() const
	{
		return new PedigreeMating(*this);
	}


	/// HIDDEN describe a pedigree mating scheme.
	virtual string describe(bool format = true) const;

	/** CPPONLY
	 *  Generate an offspring population \e scratch from parental population
	 *  \e pop.
	 */
	virtual bool mate(Population & pop, Population & scratch);

	bool parallelizable() const;

private:
	const Pedigree & m_ped;

	opList m_transmitters;

	const string m_idField;

	mutable ssize_t m_gen;
};


typedef std::vector<HomoMating *> vectormating;

/** A heterogeneous mating scheme that applies a list of homogeneous mating
 *  schemes to different (virtual) subpopulations.
 */
class HeteroMating : public MatingScheme
{
public:
	/** Create a heterogeneous mating scheme that will apply a list of
	 *  homogeneous mating schemes \e matingSchemes to different (virtual)
	 *  subpopulations. The size of the offspring generation is determined
	 *  by parameter \e subPopSize, which can be a list of subpopulation sizes
	 *  or a Python function that returns a list of subpopulation sizes at
	 *  each generation. Please refer to class \c MatingScheme for a detailed
	 *  explanation of this parameter.
	 *
	 *  Each mating scheme defined in \e matingSchemes can be applied to
	 *  one or more (virtual) subpopulation. If parameter \e subPops is not
	 *  specified, a mating scheme will be applied to all subpopulations.
	 *  If a list of (virtual) subpopulation is specified, the mating scheme
	 *  will be applied to specific (virtual) subpopulations.
	 *
	 *  If multiple mating schemes are applied to the same subpopulation, a
	 *  weight (parameter \e weight) can be given to each mating scheme to
	 *  determine how many offspring it will produce. The default \weight for
	 *  all mating schemes are \c 0. In this case, the number of offspring each
	 *  mating scheme produces is proportional to the number of individuals in its
	 *  parental (virtual) subpopulation (default to all parents, but can be
	 *  father for \c weightBy=MALE_ONLY, mother for \c weightBy=FEMALE_ONLY,
	 *  or father mother pairs (less of number of father and mothers) for
	 *  \c weightBy=PAIR_ONLY). If all weights are negative, the numbers of
	 *  offspring are determined by the multiplication of the absolute values
	 *  of the weights and their respective parental (virtual) subpopulation
	 *  sizes. If all weights are positive, the number of offspring produced by
	 *  each mating scheme is proportional to these weights, except for mating
	 *  schemes with zero parental population size (or no father, no mother, 
	 *  or no pairs, depending on value of parameter \c weightBy). Mating schemes
	 *  with zero weight in this case will produce no offspring. If both
	 *  negative and positive weights are present, negative weights are
	 *  processed before positive ones.
	 *
	 *  A sexual mating scheme might fail if a parental (virtual) subpopulation
	 *  has no father or mother. In this case, you can set \c weightBy to
	 *  \c PAIR_ONLY so a (virtual) subpopulation will appear to have zero size,
	 *  and will thus contribute no offspring to the offspring population. Note
	 *  that the perceived parental (virtual) subpopulation size in this mode (and
	 *  in modes of \c MALE_ONLY, \c FEMALE_ONLY) during the calculation of the
	 *  size of the offspring subpopulation will be roughly half of the actual
	 *  population size so you might have to use \c weight=-2 if you would like
	 *  to have an offspring subpopulation that is roughly the same size of the
	 *  parental (virtual) subpopulation.
	 *
	 *  If multiple mating schemes are applied to the same subpopulation,
	 *  offspring produced by these mating schemes are shuffled randomly. If this
	 *  is not desired, you can turn off offspring shuffling by setting parameter
	 *  \e shuffleOffspring to \c False.
	 */
	HeteroMating(const vectormating & matingSchemes,
		const uintListFunc & subPopSize = uintListFunc(),
		bool shuffleOffspring = true, SexChoice weightBy = ANY_SEX);

	/// destructor
	~HeteroMating();

	/// CPPONLY
	HeteroMating(const HeteroMating & rhs);

	/// HIDDEN Deep copy of a heterogeneous mating scheme
	virtual MatingScheme * clone() const
	{
		return new HeteroMating(*this);
	}


	/// HIDDEN describe a heterogeneous mating scheme.
	virtual string describe(bool format = true) const;

	/** CPPONLY Call each homogeneous mating scheme to populate offspring
	 *  generation.
	 */
	bool mate(Population & pop, Population & scratch);

private:
	vectormating m_matingSchemes;
	///
	bool m_shuffleOffspring;
	///
	SexChoice m_weightBy;
};


/** A conditional mating scheme that applies different mating schemes
 *  according to a condition (similar to operator IfElse). The condition can
 *  be a fixed condition, an expression or a user-defined function, to
 *  determine which mating scheme to be used.
 */
class ConditionalMating : public MatingScheme
{
public:
	/** Create a conditional mating scheme that applies mating scheme \e
	 *  ifMatingScheme if the condition \e cond is \c True, or \e
	 *  elseMatingScheme if \e cond is \c False. If a Python expression
	 *  (a string) is given to parameter \e cond, the expression will be
	 *  evalulated in parental population's local namespace. When a Python
	 *  function is specified, it accepts parameter \c pop for the parental
	 *  population. The return value of this function should be \c True or
	 *  \c False. Otherwise, parameter \e cond will be treated as a fixed
	 *  condition (converted to \c True or \c False) upon which \e ifMatingScheme
	 *  or \e elseMatingScheme will alway be applied.
	 */
	ConditionalMating(PyObject * cond, const MatingScheme & ifMatingScheme,
		const MatingScheme & elseMatingScheme);

	/// destructor
	~ConditionalMating();

	/// CPPONLY
	ConditionalMating(const ConditionalMating & rhs);

	/// HIDDEN Deep copy of a heterogeneous mating scheme
	virtual MatingScheme * clone() const
	{
		return new ConditionalMating(*this);
	}


	/// HIDDEN describe a heterogeneous mating scheme.
	virtual string describe(bool format = true) const;

	/** CPPONLY Call each homogeneous mating scheme to populate offspring
	 *  generation.
	 */
	bool mate(Population & pop, Population & scratch);

private:
	Expression m_cond;
	pyFunc m_func;
	int m_fixedCond;
	MatingScheme * m_ifMS;
	MatingScheme * m_elseMS;
};


}
#endif
