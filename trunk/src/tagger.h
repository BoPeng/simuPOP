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
 */
class idTagger : public baseOperator
{
public:
	/** Create an \c idTagger that assign a unique ID for each individual it is
	 *  applied to. The IDs are created sequentially and are stored in an
	 *  information field specified in parameter \e infoFields (default to
	 *  \c ind_id). A \e startID can be specified.
	 */
	idTagger(ULONG startID = 1, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringFunc & output = "",
		const stringList & infoFields = vectorstr(1, "ind_id")) :
		baseOperator(output, DuringMating, begin, end, step, at, reps, subPops, infoFields),
		m_id(startID)
	{
		DBG_FAILIF(infoFields.elems().size() != 1, ValueError,
			"One and only one information field is needed for idTagger.");
	};

	virtual ~idTagger()
	{
	}


	/// used by Python print function to print out the general information of the \c inheritTagger
	virtual string __repr__()
	{
		return "<simuPOP::id tagger>" ;
	}


	/** CPPONLY
	 *  apply the \c idTagger
	 */
	bool applyDuringMating(population & pop, RawIndIterator offspring,
	                       individual * dad = NULL, individual * mom = NULL)
	{
		offspring->setInfo(m_id++, infoField(0));
		return true;
	}


	/// deep copy of an \c idTagger
	virtual baseOperator * clone() const
	{
		return new idTagger(*this);
	}


private:
	int m_id;
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
	 *  \li \c mode=Average Passing the average of two values.
	 *  \li \c mode=Maximum Passing the maximum value of two values.
	 *  \li \c mode=Minumum Passing the minimum value of two values.
	 *  \li \c mode=Summation Passing the summation of two values.
	 *
	 *  An \c RuntimeError will be raised if any of the parents does not exist.
	 *  This operator does not support parameter \e subPops and does not output
	 *  any information.
	 */
	inheritTagger(InheritanceType mode = Maternal, int begin = 0, int end = -1, int step = 1,
		const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringFunc & output = "",
		const stringList & infoFields = vectorstr()) :
		baseOperator(output, DuringMating, begin, end, step, at, reps, subPops, infoFields), m_mode(mode)
	{
	};

	virtual ~inheritTagger()
	{
	}


	/// used by Python print function to print out the general information of the \c inheritTagger
	virtual string __repr__()
	{
		return "<simuPOP::inherit tagger>" ;
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

/// tagging according to parental indexes
/**
   This during-mating operator
   set \c tag() each  individual with indexes of his/her parent in the parental population.
   Because only one parent is recorded, this is recommended to be used for mating schemes
   that requires only one parent (such as selfMating).

   This tagger record indexes to information field parent_idx, and/or a given file. The usage
   is similar to parentsTagger.
 */
class parentTagger : public baseOperator
{
public:
	/// create a \c parentTagger
	// string can be any string (m_Delimiter will be ignored for this class.)
	//  %r will be replicate number %g will be generation number.
	parentTagger(int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringFunc & output = "",
		const stringList & infoFields = stringList("parent_idx")) :
		baseOperator(output, DuringMating, begin, end, step, at, reps, subPops, infoFields),
		m_subPopSize(1, 0)
	{
	};

	virtual ~parentTagger()
	{
	}


	/// deep copy of a \c parentTagger
	virtual baseOperator * clone() const
	{
		return new parentTagger(*this);
	}


	/// used by Python print function to print out the general information of the \c parentTagger
	virtual string __repr__()
	{
		return "<simuPOP::parenttagger>" ;
	}


	/** CPPONLY
	 * apply the \c parentTagger
	 */
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/** at the end of a generation, write population structure information to a file
	 * with a newline.
	 */
	bool apply(population & pop);

private:
	/// number of offspring from each subpopulation, counted
	/// from the origin of parent
	vectoru m_subPopSize;
};


/// tagging according to parents' indexes
/**
   This during-mating operator
   set \c tag(), currently a pair of numbers, of each
   individual with indexes of his/her parents in the parental population. This information
   will be used by pedigree-related operators like \c affectedSibpairSample to track
   the pedigree information. Because parental population will be discarded or stored after
   mating, these index will not be affected by post-mating operators.

   This tagger record parental index to one or both
   \li one or two information fields. Default to father_idx and mother_idx.
   If only one parent is passed in a mating scheme (such as selfing), only the first
   information field is used. If two parents are passed, the first information
   field records paternal index, and the second records maternal index.
   \li a file. Indexes will be written to this file. This tagger will also
   act as a post-mating operator to add a new-line to this file.
 */
class parentsTagger : public baseOperator
{
public:
	/// create a \c parentsTagger
	// string can be any string (m_Delimiter will be ignored for this class.)
	//  %r will be replicate number %g will be generation number.
	parentsTagger(int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringFunc & output = "",
		const stringList & infoFields = stringList(ParentsFields[0], ParentsFields[1])) :
		baseOperator(output, DuringMating, begin, end, step, at, reps, subPops, infoFields),
		m_subPopSize(1, 0)
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


	/// used by Python print function to print out the general information of the \c parentsTagger
	virtual string __repr__()
	{
		return "<simuPOP::parentstagger>" ;
	}


	/** CPPONLY
	 * apply the \c parentsTagger
	 */
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/** at the end of a generation, write population structure information to a file
	 * with a newline.
	 */
	bool apply(population & pop);

private:
	/// number of offspring from each subpopulation, counted
	/// from the origin of parents
	vectoru m_subPopSize;
};


/** Pedigree tagger is used to save a complete pedigree to a pedigree file
 *  during an evolution process.
 *  Because
 *  is destroyedof record individuals involved in an evolutioary process.
   This is a simple post-mating tagger that write given
   information fields to a file (or standard output).
 */
class pedigreeTagger : public baseOperator
{
public:
	pedigreeTagger(int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		int stage = PostMating, const stringFunc & output = ">",
		const stringList & pedigreeFields = stringList());

	bool apply(population & pop);

};


/// Python tagger
/**
   This tagger takes some information fields from both parents, pass to a Python
   function and set the individual field with the return value. This operator can
   be used to trace the inheritance of trait values.
 */
class pyTagger : public baseOperator
{
public:
	/// creates a \c pyTagger that works on specified information fields
	/**
	   \param infoFields information fields. The user should gurantee the existence
	    of these fields.
	   \param func a Pyton function that returns a list to assign the information fields.
	    e.g., if <tt>fields=['A', 'B']</tt>, the function will pass values of fields
	   \c 'A' and \c 'B' of father, followed by mother if there is one, to this
	    function. The return value is assigned to fields \c 'A' and \c 'B' of the
	    offspring. The return value has to be a list even if only one field is given.
	 */
	pyTagger(PyObject * func = NULL, int begin = 0, int end = -1,
		int step = 1, const intList & at = vectori(), const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringFunc & output = "",
		const stringList & infoFields = vectorstr()) :
		baseOperator(output, DuringMating, begin, end, step, at, reps, subPops, infoFields),
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


	/// used by Python print function to print out the general information of the \c pyTagger
	virtual string __repr__()
	{
		return "<simuPOP::pyTagger>" ;
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
