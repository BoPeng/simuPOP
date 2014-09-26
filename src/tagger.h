/**
 *  $File: tagger.h $
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

#ifndef _TAGGER_H
#define _TAGGER_H
/**
   \file
   \brief head file of class tagger: public BaseOperator
 */
#include "operator.h"

namespace simuPOP {


/** An IdTagger gives a unique ID for each individual it is applies to. These
 *  ID can be used to uniquely identify an individual in a multi-generational
 *  population and be used to reliably reconstruct a Pedigree.
 *
 *  To ensure uniqueness across populations, a single source of ID is used for
 *  this operator. individual IDs are assigned consecutively starting from 1.
 *  Value 1 instead of 0 is used because most software applications use 0 as
 *  missing values for parentship. If you would like to reset the sequence or
 *  start from a different number, you can call the \c reset(startID) function
 *  of any \c IdTagger.
 *
 *  An \c IdTagger is usually used during-mating to assign ID to each offspring.
 *  However, if it is applied directly to a population, it will assign unique
 *  IDs to all individuals in this population. This property is usually used
 *  in the \c preOps parameter of function \c Simulator.evolve to assign
 *  initial ID to a population.
 */
class IdTagger : public BaseOperator
{
public:
	/** Create an \c IdTagger that assign an unique ID for each individual it
	 *  is applied to. The IDs are created sequentially and are stored in an
	 *  information field specified in parameter \e infoFields (default to
	 *  \c ind_id). This operator is considered a during-mating operator but it
	 *  can be used to set ID for all individuals of a population when it is
	 *  directly applied to the population.
	 */
	IdTagger(int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringFunc & output = "",
		const stringList & infoFields = vectorstr(1, "ind_id")) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields)
	{
		DBG_FAILIF(infoFields.elems().size() != 1, ValueError,
			"One and only one information field is needed for IdTagger.");
	};

	virtual ~IdTagger()
	{
	}


	/// HIDDEN
	string describe(bool format = true) const;


	/** Reset the global individual ID number so that IdTaggers will start
	 *  from id (default to 1) again.
	 */
	void reset(ULONG startID = 1);

	/** HIDDEN Set an unique ID to all individuals with zero ID. */
	virtual bool apply(Population & pop) const;

	/** CPPONLY
	 *  apply the \c IdTagger
	 */
	bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// HIDDEN Deep copy of an \c IdTagger
	virtual BaseOperator * clone() const
	{
		return new IdTagger(*this);
	}


	/// CPPONLY
	bool parallelizable() const
	{
		return true;
	}


};


/** An inheritance tagger passes values of parental information field(s) to the
 *  corresponding fields of offspring. If there are two parental values from
 *  parents of a sexual mating event, a parameter \e mode is used to specify
 *  how to assign offspring information fields.
 */
class InheritTagger : public BaseOperator
{
public:
	/** Creates an inheritance tagger that passes values of parental
	 *  information fields (parameter \e infoFields) to the corresponding
	 *  fields of offspring. If there is only one parent, values at the
	 *  specified information fields are copied directly. If there are two
	 *  parents, parameter \e mode specifies how to pass them to an offspring.
	 *  More specifically,
	 *  \li \c mode=MATERNAL Passing the value from mother.
	 *  \li \c mode=PATERNAL Passing the value from father.
	 *  \li \c mode=MEAN Passing the average of two values.
	 *  \li \c mode=MAXIMUM Passing the maximum value of two values.
	 *  \li \c mode=MINIMUM Passing the minimum value of two values.
	 *  \li \c mode=SUMMATION Passing the summation of two values.
	 *  \li \c mode=MULTIPLICATION Passing the multiplication of two values.
	 *
	 *  An \c RuntimeError will be raised if any of the parents does not exist.
	 *  This operator does not support parameter \e subPops and does not output
	 *  any information.
	 */
	InheritTagger(InheritanceType mode = PATERNAL,
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringFunc & output = "",
		const stringList & infoFields = vectorstr()) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields), m_mode(mode)
	{
	};

	virtual ~InheritTagger()
	{
	}


	/// HIDDEN
	string describe(bool format = true) const;

	/** CPPONLY
	 *  apply the \c InheritTagger
	 */
	bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// HIDDEN Deep copy of a \c InheritTagger
	virtual BaseOperator * clone() const
	{
		return new InheritTagger(*this);
	}


	/// CPPONLY
	bool parallelizable() const
	{
		return true;
	}


private:
	const int m_mode;
};

/** A summary tagger summarize values of one or more parental information field
 *  to another information field of an offspring. If mating is sexual, two sets
 *  of parental values will be involved.
 */
class SummaryTagger : public BaseOperator
{
public:
	/** Creates a summary tagger that summarize values of one or more parental
	 *  information field (\e infoFields[:-1]) to an offspring information
	 *  field (\e infoFields[-1]). A parameter \e mode specifies how to pass
	 *  summarize parental values. More specifically,
	 *  \li \c mode=MEAN Passing the average of values.
	 *  \li \c mode=MAXIMUM Passing the maximum value of values.
	 *  \li \c mode=Minumum Passing the minimum value of values.
	 *  \li \c mode=SUMMATION Passing the sum of values.
	 *  \li \c mode=MULTIPLICATION Passing the multiplication of values.
	 *
	 *  This operator does not support parameter \e subPops and does not output
	 *  any information.
	 */
	SummaryTagger(InheritanceType mode = MEAN,
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringFunc & output = "",
		const stringList & infoFields = vectorstr()) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields), m_mode(mode)
	{
		DBG_FAILIF(infoFields.elems().size() < 2, ValueError,
			"Please specify at least one parental field and one offspring field.");
	};

	virtual ~SummaryTagger()
	{
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.SummaryTagger>" ;
	}


	/** CPPONLY
	 *  apply the \c SummaryTagger
	 */
	bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// HIDDEN Deep copy of a \c SummaryTagger
	virtual BaseOperator * clone() const
	{
		return new SummaryTagger(*this);
	}


private:
	const int m_mode;
};


/** This tagging operator records the indexes of parents (relative to the
 *  parental generation) of each offspring in specified information fields (
 *  default to \c father_idx and \c mother_idx). Only one information field
 *  should be specified if an asexsual mating scheme is used so there is one
 *  parent for each offspring. Information recorded by this operator is
 *  intended to be used to look up parents of each individual in
 *  multi-generational Population.
 */
class ParentsTagger : public BaseOperator
{
public:
	/** Create a parents tagger that records the indexes of parents of each
	 *  offspring when it is applied to an offspring during-mating. If two
	 *  information fields are specified (parameter \e infoFields, with default
	 *  value <tt>['father_idx', 'mother_idx']</tt>), they are used to record
	 *  the indexes of each individual's father and mother. Value \c -1 will be
	 *  assigned if any of the parent is missing. If only one information field
	 *  is given, it will be used to record the index of the first valid parent
	 *  (father if both parents are valid). This operator ignores parameters
	 *  \e output and \e subPops.
	 */
	ParentsTagger(int begin = 0, int end = -1,
		int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringFunc & output = "",
		const stringList & infoFields = stringList("father_idx", "mother_idx")) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields)
	{
	}


	virtual ~ParentsTagger()
	{
	}


	/// HIDDEN Deep copy of a \c ParentsTagger
	virtual BaseOperator * clone() const
	{
		return new ParentsTagger(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;


	/** CPPONLY
	 * apply the \c ParentsTagger
	 */
	bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// CPPONLY
	bool parallelizable() const
	{
		return true;
	}


};


/** This tagging operator records the indexes of offspring within a family
 *  (sharing the same parent or parents) in specified information field
 *  (default to \c offspring_idx). This tagger can be used to control the
 *  number of offspring during mating.
 */
class OffspringTagger : public BaseOperator
{
public:
	/** Create an offspring tagger that records the indexes of offspring within
	 *  a family. The index is determined by successful production of offspring
	 *  during a mating events so the it does not increase the index if a
	 *  previous offspring is discarded, and it resets index even if adjacent
	 *  families share the same parents. This operator ignores parameters
	 *  \e stage, \e output, and \e subPops.
	 */
	OffspringTagger(int begin = 0, int end = -1,
		int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringFunc & output = "",
		const stringList & infoFields = stringList("offspring_idx")) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields)
	{
		DBG_ASSERT(infoSize() == 1, ValueError,
			"A single field to record offspring index is required for operator OffspringTagger.");
	}


	virtual ~OffspringTagger()
	{
	}


	/// HIDDEN Deep copy of a \c OffspringTagger
	virtual BaseOperator * clone() const
	{
		return new OffspringTagger(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;


	/** CPPONLY
	 * apply the \c OffspringTagger
	 */
	bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// CPPONLY
	bool parallelizable() const
	{
		return true;
	}


};


/** This tagging operator records the ID of parents of each offspring in
 *  specified information fields (default to \c father_id and \c mother_id).
 *  Only one information field should be specified if an asexsual mating
 *  scheme is used so there is one parent for each offspring. Information
 *  recorded by this operator is intended to be used to record full pedigree
 *  information of an evolutionary process.
 */
class PedigreeTagger : public BaseOperator
{
public:
	/** Create a pedigree tagger that records the ID of parents of each
	 *  offspring when it is applied to an offspring during-mating. If two
	 *  information fields are specified (parameter \e infoFields, with default
	 *  value <tt>['father_id', 'mother_id']</tt>), they are used to record
	 *  the ID of each individual's father and mother stored in the \e idField
	 *  (default to \c ind_id) field of the parents. Value \c -1 will be
	 *  assigned if any of the parent is missing. If only one information field
	 *  is given, it will be used to record the ID of the first valid parent
	 *  (father if both pedigree are valid).
	 *
	 *  This operator by default does not send any output. If a valid output
	 *  stream is given (should be in the form of \c '>>filename' so that
	 *  output will be concatenated), this operator will output the ID of
	 *  offspring, IDs of his or her parent(s), sex and affection status of
	 *  offspring, and values at specified information fields (\e outputFields)
	 *  and loci (\e outputLoci) in the format of  <tt>off_id father_id
	 *  mother_id M/F A/U fields genotype</tt>. \c father_id or \c mother_id
	 *  will be ignored if only one parent is involved. This file format
	 *  can be loaded using function \c loadPedigree.
	 *
	 *  Because only offspring will be outputed, individuals in the top-most
	 *  ancestral generation will not be outputed. This is usually not a
	 *  problem because individuals who have offspring in the next generation
	 *  will be constructed by function \c loadPedigree, although their
	 *  information fields and genotype will be missing. If you would like
	 *  to create a file with complete pedigree information, you can apply
	 *  this operator before evolution in the \e initOps parameter of functions
	 *  \c Population.evolve or \c Simulator.evolve. This will output all
	 *  individuals in the initial population (the top-most ancestral
	 *  population after evolution) in the same format. Note that sex,
	 *  affection status and genotype can be changed by other operators so this
	 *  operator should usually be applied after all other operators are
	 *  applied.
	 */
	PedigreeTagger(const string & idField = "ind_id", const stringFunc & output = "",
		const stringList & outputFields = vectorstr(), const uintList & outputLoci = vectoru(),
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("father_id", "mother_id")) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_idField(idField), m_outputFields(outputFields), m_outputLoci(outputLoci)
	{
	}


	virtual ~PedigreeTagger()
	{
	}


	/// HIDDEN Deep copy of a \c PedigreeTagger
	virtual BaseOperator * clone() const
	{
		return new PedigreeTagger(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;

	/** HIDDEN If an output stream is specified, output current population. */
	virtual bool apply(Population & pop) const;

	/** CPPONLY
	 * apply the \c PedigreeTagger
	 */
	bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;


	/// CPPONLY
	bool parallelizable() const
	{
		return noOutput();
	}


private:
	void outputIndividual(ostream & out, const Individual * ind,
		const vectorf & IDs) const;

private:
	const string m_idField;
	stringList m_outputFields;
	uintList m_outputLoci;
};


/** A Python tagger takes some information fields from both parents, pass them
 *  to a user provided Python function and set the offspring individual fields
 *  with the return values.
 */
class PyTagger : public BaseOperator
{
public:
	/** Create a hybrid tagger that provides an user provided function \e func
	 *  with values of specified information fields (determined by parameter
	 *  names of this function) of parents and assign corresponding information
	 *  fields of offspring with its return value. If more than one parent are
	 *  available, maternal values are passed after paternal values. For
	 *  example, if a function <tt>func(A, B)</tt> is passed, this operator
	 *  will send two tuples with parental values of information fields \c 'A'
	 *  and \c 'B' to this function and assign its return values to fields
	 *  \c 'A' and \c 'B' of each offspring. The return value of this function
	 *  should be a list, although a single value will be accepted if only one
	 *  information field is specified. This operator ignores parameters
	 *  \e stage, \e output and \e subPops.
	 */
	PyTagger(PyObject * func = NULL, int begin = 0, int end = -1,
		int step = 1, const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringFunc & output = "", const stringList & infoFields = vectorstr()) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_func(func)
	{
		DBG_ASSERT(infoSize() == 0, ValueError,
			"Parameter infoFields of this operator is not used.");

		DBG_ASSERT(m_func.isValid(), ValueError,
			"Passed variable is not a callable python function.");
	};


	/// HIDDEN Deep copy of a \c PyTagger
	virtual BaseOperator * clone() const
	{
		return new PyTagger(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.PyTagger>" ;
	}


	/** CPPONLY
	 *  apply the \c PyTagger
	 */
	virtual bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/** CPPONLY
	 *  PyTagger can't be parallelizable because it call external functions
	 *  from Python which don't support multi-thread
	 */
	bool parallelizable() const
	{
		return false;
	}


private:
	const pyFunc m_func;
};


}
#endif
