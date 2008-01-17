/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu
*                                                                         *
*   $LastChangedDate$
*   $Rev$                                                      *
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

#ifndef _TAGGER_H
#define _TAGGER_H
/**
 \file
 \brief head file of class tagger: public baseOperator
 */
#include "operator.h"

const string TAG_InheritFields[2] = { "paternal_tag", "maternal_tag" };
const string TAG_ParentsFields[2] = { "father_idx", "mother_idx" };

namespace simuPOP {
/// base class of tagging individuals
/**
   This is a during-mating operator that tags individuals with various information.
   Potential usages are:
 \li recording the parental information to track pedigree;
 \li tagging an individual/allele and monitoring its spread in the population etc.
 */
class tagger : public baseOperator
{

public:
	/// create a \c tagger, default to be always active but no output
	tagger(string output = "", string outputExpr = "",
	       int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
	       int rep = REP_ALL, int grp = GRP_ALL,
	       const vectorstr & infoFields = vectorstr());

	/// destructor
	virtual ~tagger()
	{
	};

	/// deep copy of a \ tagger
	virtual baseOperator * clone() const
	{
		return new tagger(*this);
	}


	/// CPPONLY
	/// add a newline
	bool apply(population & pop);

};

/// inherite tag from parents
/**
   This during-mating operator will copy the tag (information field) from his/her parents.
   Depending on \c mode parameter, this tagger will obtain tag from his/her father
   (two tag fields), mother (two tag fields) or both (first tag field from both
   father and mother). \n

   An example may be tagging one or a few parents and examining, at the last generation,
   how many offspring they have.
 */
class inheritTagger : public tagger
{
public:
#define TAG_Paternal   0
#define TAG_Maternal   1
#define TAG_Both       2

public:
	/// create an \c inheritTagger that inherits a tag from one or both parents
	/**
	 \param mode can be one of \c TAG_Paternal, \c TAG_Maternal, and \c TAG_Both
	 */
	inheritTagger(int mode = TAG_Paternal, int begin = 0, int end = -1, int step = 1,
	              vectorl at = vectorl(), int rep = REP_ALL, int grp = GRP_ALL,
	              string output = "", string outputExpr = "",
	              const vectorstr & infoFields = vectorstr (TAG_InheritFields, TAG_InheritFields + 2)) :
		tagger(output, outputExpr, begin, end, step, at, rep, grp, infoFields), m_mode(mode)
	{
		DBG_ASSERT(infoSize() == 2, ValueError,
			"Inherit tagger needs to know the information fields of both parents");
	};

	virtual ~inheritTagger()
	{
	}


	/// used by Python print function to print out the general information of the \c inheritTagger
	virtual string __repr__()
	{
		return "<simuPOP::inherittagger>" ;
	}


	/// CPPONLY
	/// apply the \c inheritTagger
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// deep copy of a \c inheritTagger
	virtual baseOperator * clone() const
	{
		return new inheritTagger(*this);
	}


private:
	/// mode can be
	/// TAG_Paternal: get dad's info
	/// TAG_Maternal: get mon's info
	/// TAG_BOTH:     get parents' first field
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
class parentTagger : public tagger
{
public:
	/// create a \c parentTagger
	// string can be any string (m_Delimiter will be ignored for this class.)
	//  %r will be replicate number %g will be generation number.
	parentTagger(int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(), int rep = REP_ALL, int grp = GRP_ALL,
	             string output = "", string outputExpr = "",
	             const vectorstr & infoFields = vectorstr(1, "parent_idx")) :
		tagger(output, outputExpr, begin, end, step, at, rep, grp, infoFields),
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


	/// CPPONLY
	/// apply the \c parentTagger
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// at the end of a generation, write population structure information to a file
	/// with a newline.
	bool apply(population & pop);

private:
	/// number of offspring from each subpopulation, counted
	/// from the origin of parent
	vectorlu m_subPopSize;
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
class parentsTagger : public tagger
{
public:
	/// create a \c parentsTagger
	// string can be any string (m_Delimiter will be ignored for this class.)
	//  %r will be replicate number %g will be generation number.
	parentsTagger(int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(), int rep = REP_ALL, int grp = GRP_ALL,
	              string output = "", string outputExpr = "",
	              const vectorstr & infoFields = vectorstr (TAG_ParentsFields, TAG_ParentsFields + 2)) :
		tagger(output, outputExpr, begin, end, step, at, rep, grp, infoFields),
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


	/// CPPONLY
	/// apply the \c parentsTagger
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// at the end of a generation, write population structure information to a file
	/// with a newline.
	bool apply(population & pop);

private:
	/// number of offspring from each subpopulation, counted
	/// from the origin of parents
	vectorlu m_subPopSize;
};


/// Tagging information fields
/** This is a simple post-mating tagger that write given
   information fields to a file (or standard output).
 */
class infoTagger : public tagger
{
public:
	infoTagger(int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
	           int rep = REP_ALL, int grp = GRP_ALL,
	           int stage = PostMating, string output = "", string outputExpr = "",
	           const vectorstr & infoFields = vectorstr()) :
		tagger(output, outputExpr, begin, end, step, at, rep, grp, infoFields)
	{
		setApplicableStage(stage);
	}


	bool apply(population & pop);

};

/// Tagging sex status
/** This is a simple post-mating tagger that write sex
   	status to a file. By default, 1 for Male, 2 for Female.
 */
class sexTagger : public tagger
{
public:
	/**
	 \param code  code for Male and Female, default to 1 and 2, respectively.
	   This is used by Linkage format.
	 */
	sexTagger(const vectori & code = vectori(), int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(), int rep = REP_ALL, int grp = GRP_ALL,
	          int stage = PostMating, string output = "", string outputExpr = "",
	          const vectorstr & infoFields = vectorstr());

	bool apply(population & pop);

private:
	vectori m_code;
};

/// Tagging affection status
/** This is a simple post-mating tagger that write affection status
   	to a file. By default, 1 for unaffected, 2 for affected.
 */
class affectionTagger : public tagger
{
public:
	/**
	 \param code  code for Male and Female, default to 1 and 2, respectively.
	   This is used by Linkage format.
	 */
	affectionTagger(const vectori & code = vectori(), int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(), int rep = REP_ALL, int grp = GRP_ALL,
	                int stage = PostMating, string output = "", string outputExpr = "",
	                const vectorstr & infoFields = vectorstr());

	bool apply(population & pop);

private:
	vectori m_code;
};


/// Python tagger
/**
   This tagger takes some information fields from both parents, pass to a Python
   function and set the individual field with the return value. This operator can
   be used to trace the inheritance of trait values.
 */
class pyTagger : public tagger
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
	         int step = 1, vectorl at = vectorl(), int rep = REP_ALL, int grp = GRP_ALL,
	         string output = "", string outputExpr = "",
	         const vectorstr & infoFields = vectorstr()) :
		tagger(output, outputExpr, begin, end, step, at, rep, grp, infoFields)
	{
		DBG_FAILIF(infoSize() == 0, ValueError,
			"infoFields can not be empty.");

		DBG_ASSERT(PyCallable_Check(func), ValueError,
			"Passed variable is not a callable python function.");

		Py_XINCREF(func);
		m_func = func;
	};

	virtual ~pyTagger()
	{
		if (m_func != NULL)
			Py_DECREF(m_func);
	}


	/// CPPONLY
	pyTagger(const pyTagger & rhs) :
		tagger(rhs), m_func(rhs.m_func)
	{
		if (m_func != NULL)
			Py_INCREF(m_func);
	}


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


	/// CPPONLY
	/// apply the \c pyTagger
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

private:
	PyObject * m_func;
};
}
#endif
