/**
 *  $File: tagger.h $
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

#ifndef _TAGGER_H
#define _TAGGER_H
/**
   \file
   \brief head file of class tagger: public baseOperator
 */
#include "operator.h"

namespace simuPOP {


/** An idTagger gives a unique ID for each individual it is applies to. These
 *  ID can be used to uniquely identify an individual in a multi-generational
 *  population and be used to reliably reconstruct a pedigree.
 *
 *  To ensure uniqueness across populations, a single source of ID is used for
 *  this operator. Individual IDs are assigned consecutively starting from 0.
 *  If you would like to reset the sequence or start from a different number,
 *  you can call the \c reset(startID) function of any \c idTagger.
 *
 *  An \c idTagger is usually used during-mating to assign ID to each offspring.
 *  However, if it is applied directly to a population, it will assign unique
 *  IDs to all individuals in this population. This property is usually used
 *  in the \c preOps parameter of function \c simulator.evolve to assign
 *  initial ID to a population.
 *  <funcForm>TagID</funcForm>
 */
class idTagger : public baseOperator
{
public:
	/** Create an \c idTagger that assign an unique ID for each individual it
	 *  is applied to. The IDs are created sequentially and are stored in an
	 *  information field specified in parameter \e infoFields (default to
	 *  \c ind_id). This operator is considered a during-mating operator but it
	 *  can be used to set ID for all individuals of a population when it is
	 *  directly applied to the population. Because the information field is
	 *  supposed to record a unique ID for the whole population, and because
	 *  the IDs are increasingly assigned, this operator will raise a
	 *  \c RuntimeError if parental IDs are the same, or are larger than the ID
	 *  to be assigned to an offspring.
	 */
	idTagger(int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringFunc & output = "",
		const stringList & infoFields = vectorstr(1, "ind_id")) :
		baseOperator(output, begin, end, step, at, reps, subPops, infoFields)
	{
		DBG_FAILIF(infoFields.elems().size() != 1, ValueError,
			"One and only one information field is needed for idTagger.");
	};

	virtual ~idTagger()
	{
	}


	/// HIDDEN
	string describe();


	/** Reset the global individual ID number so that idTaggers will start
	 *  from id (default to 0) again.
	 */
	void reset(ULONG startID = 0);

	/** Set an unique ID to all individuals with zero ID. */
	virtual bool apply(population & pop);

	/** CPPONLY
	 *  apply the \c idTagger
	 */
	bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// deep copy of an \c idTagger
	virtual baseOperator * clone() const
	{
		return new idTagger(*this);
	}


};


/** An inheritance tagger passes values of parental information field(s) to the
 *  corresponding fields of offspring. If there are two parental values from
 *  parents of a sexual mating event, a parameter \e mode is used to specify
 *  how to assign offspring information fields.
 */
class inheritTagger : public baseOperator
{
public:
	/** Creates an inheritance tagger that passes values of parental
	 *  information fields (parameter \e infoFields) to the corresponding
	 *  fields of offspring. If there is only one parent, values at the
	 *  specified information fields are copied directly. If there are two
	 *  parents, parameter \e mode specifies how to pass them to an offspring.
	 *  More specifically,
	 *  \li \c mode=Maternal Passing the value from mother.
	 *  \li \c mode=Paternal Passing the value from father.
	 *  \li \c mode=Mean Passing the average of two values.
	 *  \li \c mode=Maximum Passing the maximum value of two values.
	 *  \li \c mode=Minumum Passing the minimum value of two values.
	 *  \li \c mode=Summation Passing the summation of two values.
	 *  \li \c mode=Multiplication Passing the multiplication of two values.
	 *
	 *  An \c RuntimeError will be raised if any of the parents does not exist.
	 *  This operator does not support parameter \e subPops and does not output
	 *  any information.
	 */
	inheritTagger(InheritanceType mode = Paternal,
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringFunc & output = "",
		const stringList & infoFields = vectorstr()) :
		baseOperator(output, begin, end, step, at, reps, subPops, infoFields), m_mode(mode)
	{
	};

	virtual ~inheritTagger()
	{
	}


	/// HIDDEN
	string describe()
	{
		return "<simuPOP.inheritTagger>" ;
	}


	/** CPPONLY
	 *  apply the \c inheritTagger
	 */
	bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// deep copy of a \c inheritTagger
	virtual baseOperator * clone() const
	{
		return new inheritTagger(*this);
	}


private:
	int m_mode;
};

/** A summary tagger summarize values of one or more parental information field
 *  to another information field of an offspring. If mating is sexual, two sets
 *  of parental values will be involved.
 */
class summaryTagger : public baseOperator
{
public:
	/** Creates a summary tagger that summarize values of one or more parental
	 *  information field (\e infoFields[:-1]) to an offspring information
	 *  field (\e infoFields[-1]). A parameter \e mode specifies how to pass
	 *  summarize parental values. More specifically,
	 *  \li \c mode=Mean Passing the average of values.
	 *  \li \c mode=Maximum Passing the maximum value of values.
	 *  \li \c mode=Minumum Passing the minimum value of values.
	 *  \li \c mode=Summation Passing the sum of values.
	 *  \li \c mode=Multiplication Passing the multiplication of values.
	 *
	 *  This operator does not support parameter \e subPops and does not output
	 *  any information.
	 */
	summaryTagger(InheritanceType mode = Mean,
		int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringFunc & output = "",
		const stringList & infoFields = vectorstr()) :
		baseOperator(output, begin, end, step, at, reps, subPops, infoFields), m_mode(mode)
	{
		DBG_FAILIF(infoFields.elems().size() < 2, ValueError,
			"Please specify at least one parental field and one offspring field.");
	};

	virtual ~summaryTagger()
	{
	}


	/// HIDDEN
	string describe()
	{
		return "<simuPOP.summaryTagger>" ;
	}


	/** CPPONLY
	 *  apply the \c summaryTagger
	 */
	bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// deep copy of a \c summaryTagger
	virtual baseOperator * clone() const
	{
		return new summaryTagger(*this);
	}


private:
	int m_mode;
};


/** This tagging operator records the indexes of parents (relative to the
 *  parental generation) of each offspring in specified information fields (
 *  default to \c father_idx and \c mother_idx). Only one information field
 *  should be specified if an asexsual mating scheme is used so there is one
 *  parent for each offspring. Information recorded by this operator is
 *  intended to be used to look up parents of each individual in
 *  multi-generational population.
 */
class parentsTagger : public baseOperator
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
	 *  \e stage, \e output, and \e subPops.
	 */
	parentsTagger(int begin = 0, int end = -1,
		int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringFunc & output = "",
		const stringList & infoFields = stringList("father_idx", "mother_idx")) :
		baseOperator(output, begin, end, step, at, reps, subPops, infoFields)
	{
	};

	virtual ~parentsTagger()
	{
	}


	/// deep copy of a \c parentsTagger
	virtual baseOperator * clone() const
	{
		return new parentsTagger(*this);
	}


	/// HIDDEN
	string describe();


	/** CPPONLY
	 * apply the \c parentsTagger
	 */
	bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

};


/** This tagging operator records the ID of parents of each offspring in
 *  specified information fields (default to \c father_id and \c mother_id).
 *  Only one information field should be specified if an asexsual mating
 *  scheme is used so there is one parent for each offspring. Information
 *  recorded by this operator is intended to be used to record full pedigree
 *  information of an evolutionary process.
 */
class pedigreeTagger : public baseOperator
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
	 *  This operator by default does not send any output but will output the
	 *  ID of offspring, father, and mother if a valid output stream is
	 *  specified. The output will be in the format of
	 *  <tt>off_id father_id mother_id</tt>. \c father_id or \c mother_id will
	 *  be ignored if only one parent is involved. This operator ignores
	 *  parameter \e stage, and \e subPops.
	 */
	pedigreeTagger(const string & idField = "ind_id", const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList("father_id", "mother_id")) :
		baseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_idField(idField)
	{
	};

	virtual ~pedigreeTagger()
	{
	}


	/// deep copy of a \c pedigreeTagger
	virtual baseOperator * clone() const
	{
		return new pedigreeTagger(*this);
	}


	/// HIDDEN
	string describe();

	/** CPPONLY
	 * apply the \c pedigreeTagger
	 */
	bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

private:
	string m_idField;
};


/** A Python tagger takes some information fields from both parents, pass them
 *  to a user provided Python function and set the offspring individual fields
 *  with the return values.
 */
class pyTagger : public baseOperator
{
public:
	/** Create a hybrid tagger that passes parental information fields
	 *  (parameter \e infoFields) to an user provided function \e func and use
	 *  its return values to assign corresponding information fields of
	 *  offspring. If more than one parent are available, maternal values are
	 *  passed after paternal values. For example, if
	 *  <tt>infoFields=['A', 'B']</tt>, the user-defined function should expect
	 *  an array of size 4, with paternal values at fields \c 'A', \c 'B',
	 *  followed by maternal values at these fields. The return value of this
	 *  function should be a list, although a single value will be accepted
	 *  if only one information field is specified. This operator ignores
	 *  parameters \e stage, \e output and \e subPops.
	 */
	pyTagger(PyObject * func = NULL, int begin = 0, int end = -1,
		int step = 1, const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringFunc & output = "", const stringList & infoFields = vectorstr()) :
		baseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_func(func)
	{
		DBG_FAILIF(infoSize() == 0, ValueError,
			"infoFields can not be empty.");

		DBG_ASSERT(m_func.isValid(), ValueError,
			"Passed variable is not a callable python function.");
	};


	/// deep copy of a \c pyTagger
	virtual baseOperator * clone() const
	{
		return new pyTagger(*this);
	}


	/// HIDDEN
	string describe()
	{
		return "<simuPOP.pyTagger>" ;
	}


	/** CPPONLY
	 *  apply the \c pyTagger
	 */
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

private:
	pyFunc m_func;
};
}
#endif
