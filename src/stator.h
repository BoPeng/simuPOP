/**************************************************************************
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

#ifndef _STATOR_H
#define _STATOR_H
/**
   \file
   \brief head file of class baseOperator:public baseOperator
 */
#include "utility.h"
#include "population.h"
#include "operator.h"

#include <string>
#include <numeric>
using std::string;
using std::count_if;
using std::bind2nd;
using std::equal;
using std::greater;
using std::min;
using std::max;

#include <iomanip>
using std::setprecision;

namespace simuPOP {


/// evaluate an expression
/**
   Python expressions/statements will be executed when \c pyEval is applied
   to a population by using parameters <tt>expr/stmts</tt>. Statements can
   also been executed when \c pyEval is created and destroyed or before \c expr
   is executed. The corresponding parameters are \c preStmts, \c postStmts
   and \c stmts. For example, operator \c varPlotter
   uses this feature to initialize R plots and save plots to a file when finished.
   <funcForm>PyEval</funcForm>
 */
class pyEval : public baseOperator
{
public:
	/// evaluate expressions/statments in the local namespace of a replicate
	/**
	   \param expr the expression to be evaluated. The result will be sent to \c output.
	   \param stmts the statement that will be executed before the expression
	   \param preStmts the statement that will be executed when the operator is constructed
	   \param postStmts the statement that will be executed when the operator is destroyed
	   \param exposePop if \c True, expose the current population as a variable named \c pop
	   \param name used to let pure Python operator to identify themselves
	   \param output default to \c >. I.e., output to standard output.
	 */
	pyEval(const string & expr = "", const string & stmts = "", const string & preStmts = "",
		const string & postStmts = "", bool exposePop = false, const string & name = "",
		string output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr())
		: baseOperator(output, stage, begin, end, step, at, rep, subPops, infoFields),
		m_expr(expr, stmts), m_postExpr("", postStmts), m_exposePop(exposePop), m_name(name)
	{
		if (preStmts != "")
			Expression("", preStmts).evaluate();
	}


	~pyEval()
	{
		m_postExpr.evaluate();
	}


	/// deep copy of a \c pyEval operator
	virtual baseOperator * clone() const
	{
		return new pyEval(*this);
	}


	// check all alleles in vector allele if they are fixed.
	/// apply the \c pyEval operator
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c pyEval operator
	virtual string __repr__()
	{
		return "<simuPOP::pyEval " + m_name + ">";
	}


	/// return the name of an expression
	/** The name of a \c pyEval operator is given by an optional parameter \c name.
	   It can be used to identify this \c pyEval operator in debug output, or in the dryrun
	   mode of simulator::evolve.
	 */
	string & name()
	{
		return m_name;
	}


private:
	/// expression to evaluate
	Expression m_expr, m_postExpr;

	/// if expose pop
	bool m_exposePop;

	string m_name;
};

/// execute a Python statement
/**
   This operator takes a list of statements and executes them. No value will be returned or outputted.
   <funcForm>PyExec</funcForm>
 */
class pyExec : public pyEval
{
public:
	/// evaluate statments in the local replicate namespace, no return value
	/**
	   Please refer to class \c pyEval for parameter descriptions.
	 */
	pyExec(const string & stmts = "", const string & preStmts = "", const string & postStmts = "",
		bool exposePop = false, const string & name = "",
		string output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr())
		: pyEval("", stmts, preStmts, postStmts, exposePop, name, "",
		         stage, begin, end, step, at, rep, subPops, infoFields)
	{
	}


	~pyExec()
	{
	}


	/// deep copy of a \c pyExec operator
	virtual baseOperator * clone() const
	{
		return new pyExec(*this);
	}


	/// used by Python print function to print out the general information of the \c pyExec operator
	virtual string __repr__()
	{
		return "<simuPOP::pyExec " + this->name() + ">";
	}


};


/** Unlike operator pyEval and pyExec that work at the population level, in
   its local namespace, infoEval works at the individual level, working
   with individual information fields.
   is statement can
   change the value of existing information fields. Optionally, variables in
   population's local namespace can be used in the statement, but this should
   be used with caution.
   <funcForm>infoEval</funcForm>
 */
class infoEval : public baseOperator
{
public:
	/// evaluate Python statements with variables being an individual's information fields
	/**
	   The expression and statements will be executed for each individual, in a
	   Python namespace (dictionary) where individual information fields are made
	   available as variables. Population dictionary can be made avaialbe with option
	   usePopVars. Changes to these variables will change the corresponding information
	   fields of individuals.

	   Please note that,
	   1. If population variables are used, and there are name conflicts between
	   information fields and variables, population variables will be overridden
	   by information fields, without any warning.
	   2. Information fields are float numbers. An exceptions will raise if an information
	   field can not be converted to a float number.
	   3. This operator can be used in all stages. When it is used during-mating,
	   it will act on each offspring.

	   \param expr the expression to be evaluated. The result will be sent to \c output.
	   \param stmts the statement that will be executed before the expression
	   \param subPop a shortcut to <tt>subPops=[subPop]</tt>
	   \param subPops subpopulations this operator will apply to. Default to all.
	   \param usePopVars if \c True, import variables from expose the current population as a variable named \c pop
	   \param exposePop if \c True, expose the current population as a variable named \c pop
	   \param name used to let pure Python operator to identify themselves
	   \param output default to \c >. I.e., output to standard output. Note that
	    because the expression will be executed for each individual, the output
	    can be large.
	 */
	infoEval(const string & expr = "", const string & stmts = "",
		bool usePopVars = false,  bool exposePop = false, const string & name = "",
		string output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr())
		: baseOperator(output, stage, begin, end, step, at, rep, subPops, infoFields),
		m_expr(expr, stmts), m_usePopVars(usePopVars), m_exposePop(exposePop),
		m_name(name), m_dict(NULL)
	{
	}


	~infoEval()
	{
		if (!m_usePopVars && m_dict != NULL)
			Py_DECREF(m_dict);
	}


	/// deep copy of a \c infoEval operator
	virtual baseOperator * clone() const
	{
		return new infoEval(*this);
	}


	// check all alleles in vector allele if they are fixed.
	/// apply the \c infoEval operator
	virtual bool apply(population & pop);

	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// used by Python print function to print out the general information of the \c infoEval operator
	virtual string __repr__()
	{
		return "<simuPOP::infoEval " + m_name + ">";
	}


	/// return the name of an expression
	/** The name of a \c infoEval operator is given by an optional parameter \c name.
	   It can be used to identify this \c infoEval operator in debug output, or in the dryrun
	   mode of simulator::evolve.
	 */
	string & name()
	{
		return m_name;
	}


private:
	void prepareDict(population & pop);

	string evalInfo(individual *);

private:
	/// expression to evaluate
	Expression m_expr;

	bool m_usePopVars;

	bool m_exposePop;

	string m_name;

	PyObject * m_dict;
};


/// execute a Python statement for each individual, using information fields
/**
   This operator takes a list of statements and executes them. No value will be returned or outputted.
   <funcForm>infoExec</funcForm>
 */
class infoExec : public infoEval
{
public:
	/// evaluate statments in the a namespace consists of individual information
	/// fields, optionally with variable in population's local namespace
	/**
	   Please refer to class \c infoEval for parameter descriptions.
	 */
	infoExec(const string & stmts = "",  bool usePopVars = false,
		bool exposePop = false, const string & name = "",
		string output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(), const vectorstr & infoFields = vectorstr())
		: infoEval("", stmts, usePopVars, exposePop, name, output,
		           stage, begin, end, step, at, rep, subPops, infoFields)
	{
	}


	~infoExec()
	{
	}


	/// deep copy of a \c infoExec operator
	virtual baseOperator * clone() const
	{
		return new infoExec(*this);
	}


	/// used by Python print function to print out the general information of the \c infoExec operator
	virtual string __repr__()
	{
		return "<simuPOP::infoExec " + this->name() + ">";
	}


};


// The following classes apply various statistics
// and stat class will provide an opearator interface
// to all of them.

//
// each class defines how to apply the statistics and
// provide interface to allow others to retrieve the value.
// NOTE: the values are population dependent so these
// values are only meaningful within the same stat.apply
// call;
//
// The design separate the calculation of statistics
// as much as possible and allow easier debug and writting
// new statistics.

/// CPPONLY return {'a-b-b'} for a b c
string haploKey(const vectori & seq);

/// CPPONLY post population sizes etc.
class statPopSize
{
private:
#define  numSubPop_String   "numSubPop"
#define  popSize_String     "popSize"
#define  virtualPopSize_String "virtualPopSize"
#define  subPopSize_String  "subPopSize"

public:
	statPopSize(bool popSize = false)
		: m_isActive(popSize)
	{
	}


	void activate()
	{
		m_isActive = true;
	}


	bool apply(population & pop);

private:
	bool m_isActive;
};

/// CPPONLY
class statNumOfMale
{
private:
#define  numOfMale_String    "numOfMale"
#define  propOfMale_String   "propOfMale"
#define  numOfFemale_String  "numOfFemale"
#define  propOfFemale_String   "propOfFemale"

public:
	statNumOfMale(bool numOfMale = false, const strDict & param = strDict())
		: m_numOfMale(numOfMale ? 1 : 0), m_numOfFemale(0),
		m_evalInSubPop(true),
		m_output_numOfMale(true),
		m_output_propOfMale(true),
		m_output_numOfFemale(true),
		m_output_propOfFemale(true)
	{
		if (!param.empty()) {
			strDict::const_iterator it;
			strDict::const_iterator itEnd = param.end();
			if ((it = param.find("subPop")) != itEnd)
				m_evalInSubPop = it->second != 0.;
			if (param.find(numOfMale_String) != itEnd ||
			    param.find(propOfMale_String) != itEnd ||
			    param.find(numOfFemale_String) != itEnd ||
			    param.find(propOfFemale_String) != itEnd) {
				m_output_numOfMale = false;
				m_output_propOfMale = false;
				m_output_numOfFemale = false;
				m_output_propOfFemale = false;
				if ((it = param.find(numOfMale_String)) != itEnd)
					m_output_numOfMale = it->second != 0.;
				if ((it = param.find(propOfMale_String)) != itEnd)
					m_output_propOfMale = it->second != 0.;
				if ((it = param.find(numOfFemale_String)) != itEnd)
					m_output_numOfFemale = it->second != 0.;
				if ((it = param.find(propOfFemale_String)) != itEnd)
					m_output_propOfFemale = it->second != 0.;
			}
		}
	}


	void activate(bool yes = true)
	{
		m_numOfMale.resize(yes ? 1 : 0);
		m_numOfFemale.resize(yes ? 1 : 0);
	}


	ULONG numOfMale()
	{
		DBG_ASSERT(m_numOfMale.size() > 1, ValueError,
			"num of male has not been counted.");

		return m_numOfMale[ m_numOfMale.size() - 1 ];
	}


	ULONG numOfFemale()
	{
		DBG_ASSERT(m_numOfFemale.size() > 1, ValueError,
			"num of female has not been counted.");

		return m_numOfFemale[ m_numOfFemale.size() - 1 ];
	}


	ULONG numOfMale(UINT subPop)
	{
		DBG_ASSERT(m_numOfMale.size() > 1, ValueError,
			"num of male has not been counted.");

		DBG_ASSERT(subPop >= m_numOfMale.size() - 1, ValueError,
			"subPop index out of range.");

		return m_numOfMale[ subPop ];
	}


	ULONG numOfFemale(UINT subPop)
	{
		DBG_ASSERT(m_numOfFemale.size() > 1, ValueError,
			"num of male has not been counted.");

		DBG_ASSERT(subPop >= m_numOfFemale.size() - 1, ValueError,
			"subPop index out of range.");

		return m_numOfFemale[ subPop ];
	}


	bool apply(population & pop);

private:
	/// whether or not to apply number of male/female
	vectorlu m_numOfMale, m_numOfFemale;
	bool m_evalInSubPop;
	bool m_output_numOfMale;
	bool m_output_propOfMale;
	bool m_output_numOfFemale;
	bool m_output_propOfFemale;
};

/// CPPONLY
class statNumOfAffected
{
private:
#define  numOfAffected_String    "numOfAffected"
#define  propOfAffected_String   "propOfAffected"
#define  numOfUnaffected_String  "numOfUnaffected"
#define  propOfUnaffected_String  "propOfUnaffected"

public:
	statNumOfAffected(bool numOfAffected = false, const strDict & param = strDict())
		: m_numOfAffected(numOfAffected ? 1 : 0), m_numOfUnaffected(0),
		m_evalInSubPop(true),
		m_output_numOfAffected(true),
		m_output_propOfAffected(true),
		m_output_numOfUnaffected(true),
		m_output_propOfUnaffected(true)
	{
		if (!param.empty()) {
			strDict::const_iterator it;
			strDict::const_iterator itEnd = param.end();
			if ((it = param.find("subPop")) != itEnd)
				m_evalInSubPop = it->second != 0.;
			if (param.find(numOfAffected_String) != itEnd ||
			    param.find(propOfAffected_String) != itEnd ||
			    param.find(numOfUnaffected_String) != itEnd ||
			    param.find(propOfUnaffected_String) != itEnd) {
				m_output_numOfAffected = false;
				m_output_propOfAffected = false;
				m_output_numOfUnaffected = false;
				m_output_propOfUnaffected = false;
				if ((it = param.find(numOfAffected_String)) != itEnd)
					m_output_numOfAffected = it->second != 0.;
				if ((it = param.find(propOfAffected_String)) != itEnd)
					m_output_propOfAffected = it->second != 0.;
				if ((it = param.find(numOfUnaffected_String)) != itEnd)
					m_output_numOfUnaffected = it->second != 0.;
				if ((it = param.find(propOfUnaffected_String)) != itEnd)
					m_output_propOfUnaffected = it->second != 0.;
			}
		}
	}


	~statNumOfAffected()
	{
	}


	void activate(bool yes = true)
	{
		m_numOfAffected.resize(yes ? 1 : 0);
		m_numOfUnaffected.resize(yes ? 1 : 0);
	}


	ULONG numOfAffected()
	{
		DBG_ASSERT(m_numOfAffected.size() > 1, ValueError,
			"num of affected has not been counted.");

		return m_numOfAffected[ m_numOfAffected.size() - 1 ];
	}


	ULONG numOfUnaffected()
	{
		DBG_ASSERT(m_numOfUnaffected.size() > 1, ValueError,
			"num of unaffected has not been counted.");

		return m_numOfUnaffected[ m_numOfUnaffected.size() - 1 ];
	}


	ULONG numOfAffected(UINT subPop)
	{
		DBG_ASSERT(m_numOfAffected.size() > 1, ValueError,
			"num of affected has not been counted.");

		DBG_ASSERT(subPop >= m_numOfAffected.size() - 1, ValueError,
			"subPop index out of range.");

		return m_numOfAffected[ subPop ];
	}


	ULONG numOfUnaffected(UINT subPop)
	{
		DBG_ASSERT(m_numOfUnaffected.size() > 1, ValueError,
			"num of unaffected has not been counted.");

		DBG_ASSERT(subPop >= m_numOfUnaffected.size() - 1, ValueError,
			"subPop index out of range.");

		return m_numOfUnaffected[ subPop ];
	}


	bool apply(population & pop);

private:
	/// record the result.
	vectorlu m_numOfAffected, m_numOfUnaffected;
	bool m_evalInSubPop;
	bool m_output_numOfAffected;
	bool m_output_propOfAffected;
	bool m_output_numOfUnaffected;
	bool m_output_propOfUnaffected;
};

/// CPPONLY
class statAlleleFreq
{
private:
#define  NumOfAlleles_String  "numOfAlleles"
#define  AlleleNum_String     "alleleNum"
#define  AlleleFreq_String    "alleleFreq"

public:
	statAlleleFreq(const vectori & atLoci = vectori(), const strDict & param = strDict())
		: m_atLoci(atLoci), m_ifPost(atLoci.size()), m_numOfAlleles(0),
		m_alleleNum(0), m_alleleFreq(0),
		m_evalInSubPop(true),
		m_output_alleleNum(true),
		m_output_alleleFreq(true),
		m_output_numOfAlleles(false)
	{
		for (size_t i = 0; i < atLoci.size(); ++i)
			m_ifPost[i] = 1; // true, post result
		if (!param.empty()) {
			strDict::const_iterator it;
			strDict::const_iterator itEnd = param.end();
			if ((it = param.find("subPop")) != itEnd)
				m_evalInSubPop = it->second != 0.;
			if (param.find(AlleleNum_String) != itEnd ||
			    param.find(AlleleFreq_String) != itEnd ||
			    param.find(NumOfAlleles_String) != itEnd) {
				m_output_alleleNum = false;
				m_output_alleleFreq = false;
				m_output_numOfAlleles = false;
				if ((it = param.find(AlleleNum_String)) != itEnd)
					m_output_alleleNum = it->second != 0.;
				if ((it = param.find(AlleleFreq_String)) != itEnd)
					m_output_alleleFreq = it->second != 0.;
				if ((it = param.find(NumOfAlleles_String)) != itEnd)
					m_output_numOfAlleles = it->second != 0.;
			}
		}
	}


	/// destructor, nested vectors have to be cleared manually
	~statAlleleFreq();

	void addLocus(int locus, bool post, bool subPop, bool numOfAlleles);

	intMatrix & alleleNumAll()
	{
		return m_alleleNum.back();
	}


	vectori & alleleNumVec(int loc)
	{
		return m_alleleNum.back()[loc];
	}


	int alleleNum(UINT allele, int loc)
	{
		// test for size?
		vectori & an = m_alleleNum.back()[loc];

		return allele < an.size() ? an[allele] : 0;
	}


	matrix & alleleFreqAll()
	{
		return m_alleleFreq.back();
	}


	vectorf & alleleFreqVec(int loc)
	{
		return m_alleleFreq.back()[loc];
	}


	double alleleFreq(UINT allele, int loc)
	{
		vectorf & af = m_alleleFreq.back()[loc];

		return allele < af.size() ? af[allele] : 0;
	}


	intMatrix & alleleNumAll(UINT subPop)
	{
		return m_alleleNum[subPop];
	}


	vectori & alleleNumVec(int loc, UINT subPop)
	{
		return m_alleleNum[subPop][loc];
	}


	int alleleNum(UINT allele, int loc, UINT subPop)
	{
		vectori & an = m_alleleNum[subPop][loc];

		return allele < an.size() ? an[allele] : 0;
	}


	matrix & alleleFreqAll(UINT subPop)
	{
		return m_alleleFreq[subPop];
	}


	vectorf & alleleFreqVec(int loc, UINT subPop)
	{
		return m_alleleFreq[subPop][loc];
	}


	double alleleFreq(UINT allele, int loc, UINT subPop)
	{
		vectorf & af = m_alleleFreq[subPop][loc];

		return allele < af.size() ? af[allele] : 0.;
	}


	vectori & numOfAlleles()
	{
		UINT idx = m_numOfAlleles.size() == 2 ? 0 : (m_numOfAlleles.size() - 1);

		return m_numOfAlleles[idx];
	}


	vectori & numOfAlleles(UINT subPop)
	{
		DBG_ASSERT(subPop < m_numOfAlleles.size() - 1, IndexError,
			"Subpop index " + toStr(subPop) + " out of range of 0 ~ "
			+ toStr(m_numOfAlleles.size() - 2));
		return m_numOfAlleles[subPop];
	}


	vectori alleles(int loc)
	{
		vectori al;

		for (size_t j = 0; j < m_alleleNum.back()[loc].size(); ++j)
			if (m_alleleNum.back()[loc][j] > 0)
				al.push_back(j);

		DBG_ASSERT(al.size() == static_cast<UINT>(numOfAlleles()[loc]),
			SystemError, "Number of alleles at locus " + toStr(loc)
			+ " does not match.Observed "
			+ toStr(al.size()) + " previous count: "
			+ toStr(numOfAlleles()[loc]));

		return al;
	}


	vectori alleles(int loc, UINT subPop)
	{
		vectori al;

		// whether or not count 0 allele?
		// consider NA as an allele is reasonable.
		for (size_t j = 0; j < m_alleleNum[subPop][loc].size(); ++j)
			if (m_alleleNum[subPop][loc][j] != 0)
				al.push_back(j);

#ifndef BINARYALLELE
		DBG_WARNING(m_alleleNum[subPop][loc][0] != 0,
			"Having zero (NA) allele, counted as one allele.");
#endif

		DBG_ASSERT(al.size() == static_cast<UINT>(numOfAlleles(subPop)[loc]),
			SystemError, "Number of alleles at locus " + toStr(loc)
			+ " at subpop " + toStr(subPop) + " does not match. Observed "
			+ toStr(al.size()) + " previous count: "
			+ toStr(numOfAlleles(subPop)[loc]));

		return al;
	}


	bool apply(population & pop);

private:
	/// which alleles?
	vectori m_atLoci;

	/// whether or not post result
	vectori m_ifPost;

	/// count number of alleles
	intMatrix m_numOfAlleles;

	/// allele counter! use this matrix to avoid frequent
	/// allocation of memory
	vector<intMatrix> m_alleleNum;

	/// allele Freq
	vector<matrix> m_alleleFreq;

	///
	bool m_evalInSubPop;
	bool m_output_alleleNum;
	bool m_output_alleleFreq;
	bool m_output_numOfAlleles;
};

// use alleleFreq to get number of alleles
// so statNumOfAlleles is just a proxy class
/// CPPONLY
class statNumOfAlleles
{
public:
	statNumOfAlleles(statAlleleFreq & calc, const vectori & atLoci = vectori(),
		const strDict & param = strDict())
		: m_calc(calc), m_evalInSubPop(true)
	{
		if (!param.empty()) {
			strDict::const_iterator it;
			strDict::const_iterator itEnd = param.end();
			if ((it = param.find("subPop")) != itEnd)
				m_evalInSubPop = it->second != 0.;
		}
		for (vectori::const_iterator it = atLoci.begin(); it != atLoci.end(); ++it)
			m_calc.addLocus(*it, true, m_evalInSubPop, true);
	}


	~statNumOfAlleles()
	{
	}


	// do nothing. m_calc.spply will be called by stat.
	bool apply(population & pop)
	{
		return true;
	}


private:
	/// a reference to an existing allelefreq calculator
	statAlleleFreq & m_calc;

	bool m_evalInSubPop;
};

/// CPPONLY
class statHeteroFreq
{
private:
#define HeteroNum_String        "heteroNum"
#define HeteroFreq_String       "heteroFreq"
#define AllHeteroNum_String     "HeteroNum"
#define AllHeteroFreq_String    "HeteroFreq"
#define HomoNum_String          "homoNum"
#define HomoFreq_String         "homoFreq"

	int locusIdx(int loc)
	{
		UINT idx = 0;

		while (m_atLoci[idx] != loc && idx < m_atLoci.size() )
			idx++;
		DBG_ASSERT(m_atLoci[idx] == loc, ValueError,
			"Can not find allele freq for locus " + toStr(loc));
		return idx;
	}


	// the result layout is
	//  loc 1, sup1
	//  loc 2, sup1
	//  ...
	//  loc 1 , sup2
	//  loc 2, sup2
	//  ...
	//  loc 1, summary
	//  loc 2, summary
	//  ...

	int resIdx(int idx)                                                       // for summary
	{
		UINT numSP = m_heteroNum.size() / m_atLoci.size() - 1;

		if (numSP == 1)
			return idx;
		else
			return idx + numSP * m_atLoci.size();
	}


	int resIdx(int idx, UINT subPop)
	{
		DBG_ASSERT(subPop < m_heteroNum.size() / m_atLoci.size() - 1,
			IndexError,
			"Subpop index " + toStr(subPop) + " out of range of 0 ~ "
			+ toStr(m_heteroNum.size() / m_atLoci.size() - 2));

		return idx + subPop * m_atLoci.size();
	}


public:
	statHeteroFreq(const vectori & heteroFreq = vectori(),
		const vectori & homoFreq = vectori())
		: m_atLoci(heteroFreq), m_ifPost(0),
		m_postHetero(!heteroFreq.empty()), m_postHomo(!homoFreq.empty()),
		m_heteroNum(0), m_heteroFreq(0), m_homoNum(0), m_homoFreq(0)
	{
		// add homofreq to m_atLoci
		for (size_t i = 0; i < homoFreq.size(); ++i)
			if (find(m_atLoci.begin(), m_atLoci.end(), homoFreq[i]) == m_atLoci.end())
				m_atLoci.push_back(homoFreq[i]);

		// post result
		m_ifPost.resize(m_atLoci.size(), 1);
	}


	void addLocus(int locus, bool post = false)
	{
		if (find(m_atLoci.begin(), m_atLoci.end(), locus) == m_atLoci.end() ) {
			m_atLoci.push_back(locus);
			// default not post result
			m_ifPost.push_back(static_cast<int>(post));
		}
		m_postHetero = true;
	}


	int heteroNum(UINT allele, int loc)
	{
		UINT idx = locusIdx(loc);

		vectori & hn = m_heteroNum[resIdx(idx)];

		return allele < hn.size() ? hn[allele] : 0;
	}


	double heteroFreq(UINT allele, int loc)
	{
		UINT idx = locusIdx(loc);
		vectorf & hf = m_heteroFreq[resIdx(idx)];

		return allele < hf.size() ? hf[allele] : 0.;
	}


	int heteroNum(UINT allele, int loc, UINT subPop)
	{
		UINT idx = locusIdx(loc);
		vectori & hn = m_heteroNum[resIdx(idx, subPop)];

		return allele < hn.size() ? hn[allele] : 0;
	}


	double heteroFreq(UINT allele, int loc, UINT subPop)
	{
		UINT idx = locusIdx(loc);
		vectorf & hf = m_heteroFreq[resIdx(idx, subPop)];

		return allele < hf.size() ? hf[allele] : 0;
	}


	bool apply(population & pop);

private:
	/// heteroFreq
	vectori m_atLoci;

	///
	vectori m_ifPost;

	bool m_postHetero, m_postHomo;

	/// hetero counter
	intMatrix m_heteroNum;

	/// hetero Freq
	matrix m_heteroFreq;

	/// expected heterozygosity
	matrix m_expHetero;

	/// homozygosity number and freq
	intMatrix m_homoNum;

	matrix m_homoFreq;
};

/// CPPONLY
class statExpHetero
{
private:
#define ExpHetero_String "expHetero"

public:
	statExpHetero(statAlleleFreq & alleleFreq, const vectori & expHetero = vectori(),
		const strDict & param = strDict())
		: m_alleleFreq(alleleFreq), m_atLoci(expHetero), m_expHetero(0),
		m_midValues(false), m_evalInSubPop(true)
	{
		if (!param.empty()) {
			strDict::const_iterator it;
			strDict::const_iterator itEnd = param.end();
			if ((it = param.find("subPop")) != itEnd)
				m_evalInSubPop = it->second != 0.;
			if ((it = param.find("midValues")) != itEnd)
				m_midValues = it->second != 0.;
		}
		// add expected hetero to m_alleleFreq
		for (size_t i = 0; i < expHetero.size(); ++i)
			m_alleleFreq.addLocus(expHetero[i], m_midValues, m_evalInSubPop, false);
	}


	bool apply(population & pop);

private:
	/// need this to apply alleleFreq
	statAlleleFreq & m_alleleFreq;

	/// heteroFreq
	vectori m_atLoci;

	/// expected heterozygosity
	matrix m_expHetero;

	/// whether or not keep intermediate values
	bool m_midValues;

	/// whether or not calculate statistics for subpopulations
	bool m_evalInSubPop;

};

// currently there is no need to expose the result.
// may add that later.
/// CPPONLY
class statGenoFreq
{
private:
#define  GenotypeNum_String   "genoNum"
#define  GenotypeFreq_String  "genoFreq"

public:
	statGenoFreq(const vectori & genoFreq = vectori(),
		const strDict & param = strDict());

	bool apply(population & pop);

private:
	/// which genotypes
	vectori m_atLoci;

	/// phase
	bool m_phase;
};

/// CPPONLY
class statHaploFreq
{
private:
#define  HaplotypeNum_String    "haploNum"
#define  HaplotypeFreq_String   "haploFreq"
	int haploIndex(const vectori & haplo)
	{
		// first locate haplo
		UINT idx = 0;

		while (m_haplotypes[idx] != haplo && idx < m_haplotypes.size())
			idx++;

		DBG_ASSERT(idx != m_haplotypes.size(), ValueError,
			"Can not find haplotype." + toStr(haplo[0]) + ", " + toStr(haplo[1]));

		return idx;
	}


public:
	statHaploFreq(const intMatrix & haploFreq = intMatrix())
		: m_haplotypes(haploFreq), m_ifPost(haploFreq.size())
	{
		for (size_t i = 0; i < haploFreq.size(); ++i)
			m_ifPost[i] = 1;
	}


	~statHaploFreq()
	{
	}


	void addHaplotype(const vectori & haplo, bool post = false);

	int numOfHaplotypes(const vectori & haplo)
	{
		// first locate haplo
		UINT idx = haploIndex(haplo);

		return m_haploNum[idx + m_haploNum.size() -
		                  m_haplotypes.size()].size();
	}


	int numOfHaplotypes(const vectori & haplo, UINT subPop)
	{
		// first locate haplo
		UINT idx = haploIndex(haplo);

		return m_haploNum[idx + subPop * m_haplotypes.size()].size();
	}


	map<vectori, UINT> & haploNum(const vectori & haplo)
	{
		UINT idx = haploIndex(haplo);

		return m_haploNum[idx + m_haploNum.size() - m_haplotypes.size() ];
	}


	map<vectori, double> & haploFreq(const vectori & haplo)
	{
		UINT idx = haploIndex(haplo);

		return m_haploFreq[idx + m_haploNum.size() - m_haplotypes.size() ];
	}


	map<vectori, UINT> & haploNum(const vectori & haplo, UINT subPop)
	{
		UINT idx = haploIndex(haplo);

		return m_haploNum[idx + subPop * m_haplotypes.size() ];
	}


	map<vectori, double> & haploFreq(const vectori & haplo, UINT subPop)
	{
		UINT idx = haploIndex(haplo);

		return m_haploFreq[idx + subPop * m_haplotypes.size() ];
	}


	bool apply(population & pop);

private:
	/// which haplotypes
	intMatrix m_haplotypes;

	vectori m_ifPost;

	/// keep results
	vector< map< vectori, UINT> > m_haploNum;

	/// keep result
	vector< map< vectori, double> > m_haploFreq;
};

/// CPPONLY
class statLD
{
private:
	// these are names of calcualted statistics, will be accessed like
	// pop.dvars().r2 or pop.vars()['r2']
#define   LD_String           "ld"
#define   LDPRIME_String      "ld_prime"
#define   R2_String           "r2"
#define   DELTA2_String       "delta2"
	// these are LD averaged across all alleles
	// for diallelic loci, they are the same as single-allele values
#define   AvgLD_String        "LD"
#define   AvgLDPRIME_String   "LD_prime"
#define   AvgR2_String        "R2"
#define   AvgDELTA2_String    "Delta2"

public:
	// alleleFreq and haploFreq is required to calculate LD
	// needed allele and halplotype are added to alleleFreq and haploFreq
	// objects during the initialization of statLD, as well as stat.
	// In stat::apply(), alleleFreq.apply() and haploFreq.apply()
	// is called before statLD.apply() and ensures that allele frequencies
	// are calculated when statLD needs them.
	statLD(statAlleleFreq & alleleFreq, statHaploFreq & haploFreq,
		const intMatrix & LD = intMatrix(), const strDict & LD_param = strDict());

	// calculate, right now,  do not tempt to save values
	bool apply(population & pop);

private:
	// calculate single allele LD values
	void calculateLD(const vectori & hapLoci, const vectori & hapAlleles, UINT sp, bool subPop,
		double & P_A, double & P_B, double & D, double & D_prime, double & r2, double & delta2);

	// output statistics
	void outputLD(population & pop, const vectori & hapLoci, const string & allele_string, UINT sp, bool subPop,
		bool valid_delta2, double D, double D_prime, double r2, double delta2);

private:
	/// need to get allele freq
	statAlleleFreq & m_alleleFreq;

	/// need to get haplofreq
	statHaploFreq & m_haploFreq;

	/// LD
	intMatrix m_LD;

	/// whether or not keep intermediate values
	bool m_midValues;

	/// whether or not calculate statistics for subpopulations
	bool m_evalInSubPop;

	/// whether or not calculate the following statistics
	bool m_output_ld;
	bool m_output_ld_prime;
	bool m_output_r2;
	bool m_output_delta2;
	bool m_output_LD;
	bool m_output_LD_prime;
	bool m_output_R2;
	bool m_output_Delta2;
};

/// CPPONLY
class statAssociation
{
private:
	// these are names of calcualted statistics, will be accessed like
	// pop.dvars().Chisq or pop.vars()['Chisq']
#define   ChiSq_String      "ChiSq"
#define   ChiSq_P_String    "ChiSq_P"
#define   UCU_String        "UC_U"
#define   CramerV_String    "CramerV"

public:
	// alleleFreq and haploFreq is required to calculate Chisq
	// needed allele and halplotype are added to alleleFreq and haploFreq
	// objects during the initialization of statAssociation, as well as stat.
	// In stat::apply(), alleleFreq.apply() and haploFreq.apply()
	// is called before statAssociation.apply() and ensures that allele frequencies
	// are calculated when statAssociation needs them.
	statAssociation(statAlleleFreq & alleleFreq, statHaploFreq & haploFreq,
		const intMatrix & Association = intMatrix(), const strDict & param = strDict());

	// calculate, right now,  do not tempt to save values
	bool apply(population & pop);

private:
	/// need to get allele freq
	statAlleleFreq & m_alleleFreq;

	/// need to get haplofreq
	statHaploFreq & m_haploFreq;

	/// Association
	intMatrix m_association;

	///
	bool m_midValues;
	bool m_evalInSubPop;
	bool m_output_ChiSq;
	bool m_output_UCU;
	bool m_output_CramerV;

};

/// CPPONLY currently there is no need to retrieve calculated value
class statFst
{

private:
#define  Fst_String  "Fst"
#define  Fis_String  "Fis"
#define  Fit_String  "Fit"
#define  AvgFst_String  "AvgFst"
#define  AvgFis_String  "AvgFis"
#define  AvgFit_String  "AvgFit"

public:
	statFst(statAlleleFreq & alleleFreq, statHeteroFreq & heteroFreq,
		const vectori & Fst = vectori(), const strDict & param = strDict());

	double Fst()
	{
		return m_avgFst;
	}


	double Fis()
	{
		return m_avgFis;
	}


	double Fit()
	{
		return m_avgFit;
	}


	double Fst(UINT loc)
	{
		return m_Fst[loc];
	}


	double Fis(UINT loc)
	{
		return m_Fis[loc];
	}


	double Fit(UINT loc)
	{
		return m_Fit[loc];
	}


	bool apply(population & pop);

private:
	statAlleleFreq & m_alleleFreq;
	statHeteroFreq & m_heteroFreq;

	/// Fst
	vectori m_atLoci;

	/// result
	vectorf m_Fst, m_Fit, m_Fis;

	double m_avgFst, m_avgFit, m_avgFis;

	bool m_midValues;
	bool m_output_Fst;
	bool m_output_Fis;
	bool m_output_Fit;
	bool m_output_AvgFst;
	bool m_output_AvgFis;
	bool m_output_AvgFit;
};

#define REL_Queller             1
#define REL_Lynch               2
#define REL_IR                  3
#define REL_D2                  4
#define REL_Rel                 5

// the relatedness measure between two individuals/families
// using Queller and Goodnight or Lynch's method.
// or internal relatedness values
/// CPPONLY
class statRelatedness
{

public:
#define Rel_Queller_String "relQueller"
#define Rel_Lynch_String   "relLynch"
#define Rel_IR_String      "relIR"
#define Rel_D2_String      "relD2"
#define Rel_Rel_String     "relRel"

public:
	typedef std::pair<double, double> fraction;

public:
	/// \brief calculate relatedness measures between elements in groups
	/**
	   \param groups can be [ [1,2,3],[4,5,6],[7,8,9]] as three groups of
	   individuals; or [ 1 3 4] as three subpopulations. To specify between
	   individual relatedness, use [[1],[2],[3]] (the first form). If this
	   parameter is ignored, this operator calculate relatedness between
	   all subpopulations.

	   \param method can be REL_Queller, REL_Lynch, REL_IR, REL_D2
	   or REL_Rel. Please refer to the manual for details.
	 */
	statRelatedness(statAlleleFreq & alleleFreq, const intMatrix & groups = intMatrix(),
		bool useSubPop = false, const vectori & loci = vectori(), vectori method = vectori(),
		int minScored = 10, const strDict & param = strDict());

	// relatedness between individuals
	fraction relQueller(individual ind1,
		individual ind2);

	fraction relLynch(individual ind1,
		individual ind2);

	// IR measure for individual ind at specified locus
	fraction relIR(individual ind1, int locus);

	// D2 measure for individual ind at specified locus
	fraction relD2(individual ind1, int locus);

	// REL measure for individual ind at specified locus
	fraction relRel(individual ind1,
		individual ind2,  int locus);

	// between group i and j if method=REL_Queller and REL_Lynch
	/// for group i and locus j otherwise
	double groupRelatedness(population & pop, int i, int j, int method);

	bool apply(population & pop);

private:
	/// need to get allele freq
	statAlleleFreq & m_alleleFreq;

	/// method to use
	intMatrix m_groups;

	///
	bool m_useSubPop;

	/// loci used
	vectori m_atLoci;

	/// method
	vectori m_method;

	/// control number of scored loci
	int m_minScored;

	// save result
	matrix m_relQueller, m_relLynch, m_relIR, m_relD2, m_relRel;

	bool m_midValues;
};

/// calculate statistics
/**
   Operator \c stat calculates various basic statistics for the population
   and sets variables in the local namespace. Other operators or functions can
   refer to the results from the namespace after \c stat is applied. \c Stat
   is the function form of the operator. \n

   Note that these statistics are dependent to each other. For example,
   heterotype and allele frequencies of related loci will be automatically
   calculated if linkage diseqilibrium is requested.

   <funcForm>Stat</funcForm>
 */
class stat : public baseOperator
{
public:
	/// create an \c stat operator
	/**

	   \param popSize whether or not calculate population and virtual subpopulation
	    sizes. This parameter will set the following variables:
	   \li \c numSubPop the number of subpopulations.
	   \li \c subPopSize an array of subpopulation sizes.
	   \li \c virtualSubPopSize (optional) an array of virtual subpopulation sizes.
	    If a subpopulation does not have any virtual subpopulation, the
	    subpopulation size is returned.
	   \li \c popSize, <tt>subPop[sp]['popSize']</tt> the population/subpopulation size.

	   \param numOfMale whether or not count the numbers or proportions of males and females.
	    This parameter can set the following variables by user's specification:
	   \li \c numOfMale, <tt>subPop[sp]['numOfMale']</tt> the number of males in the population/subpopulation.
	   \li \c numOfFemale, <tt>subPop[sp]['numOfFemale']</tt> the number of females in the population/subpopulation.
	   \li \c propOfMale, <tt>subPop[sp]['propOfMale']</tt> the proportion of
	   males in the population/subpopulation.
	   \li \c propOfFemale, <tt>subPop[sp]['propOfFemale']</tt> the proportion of
	   females in the population/subpopulation.

	   \param numOfMale_param a dictionary of parameters of \c numOfMale statistics.
	   Can be one or more items choosen from the following options: \c numOfMale,
	   \c propOfMale, \c numOfFemale, and \c propOfFemale.

	   \param numOfAffected whether or not count the numbers or proportions of affected and unaffected individuals.
	   This parameter can set the following variables by user's specification:
	   \li \c numOfAffected, <tt>subPop[sp]['numOfAffected']</tt> the number of
	   affected individuals in the population/subpopulation.
	   \li \c numOfUnaffected, <tt>subPop[sp]['numOfUnAffected']</tt> the number of
	   unaffected individuals in the population/subpopulation.
	   \li \c propOfAffected, <tt>subPop[sp]['propOfAffected']</tt> the proportion of
	   affected individuals in the population/subpopulation.
	   \li \c propOfUnaffected, <tt>subPop[sp]['propOfUnAffected']</tt> the proportion of
	   unaffected individuals in the population/subpopulation.

	   \param numOfAffected_param a dictionary of parameters of \c numOfAffected statistics.
	   Can be one or more items choosen from the following options: \c numOfAffected,
	   \c propOfAffected, \c numOfUnaffected, \c propOfUnaffected.

	   \param numOfAlleles an array of loci at which the numbers of distinct alleles
	   will be counted (<tt>numOfAlleles=[loc1, loc2, ...]</tt> where \c loc1 etc.
	   are absolute locus indexes). This is done through the calculation of allele
	   frequencies. Therefore, allele frequencies will also be calculated if this
	   statistics is requested. This parameter will set the following variables
	   (\c carray objects of the numbers of alleles for <em>all loci</em>). Unrequested loci will
	   have \c 0 distinct alleles.
	   \li \c numOfAlleles, <tt>subPop[sp]['numOfAlleles']</tt> the number of distinct
	   alleles at each locus. (Calculated only at requested loci.)

	   \param numOfAlleles_param a dictionary of parameters of \c numOfAlleles statistics.
	   Can be one or more items choosen from the following options: \c numOfAffected,
	   \c propOfAffected, \c numOfUnaffected, \c propOfUnaffected.

	   \param alleleFreq an array of loci at which all allele frequencies will be
	   calculated (<tt>alleleFreq=[loc1, loc2, ...]</tt> where \c loc1 etc. are
	   loci where allele frequencies will be calculated). This parameter will set
	   the following variables (\c carray objects); for example, <tt>alleleNum[1][2]</tt>
	   will be the number of allele \c 2 at locus \c 1:
	   \li <tt>alleleNum[a]</tt>, <tt>subPop[sp]['alleleNum'][a]</tt>
	   \li <tt>alleleFreq[a]</tt>, <tt>subPop[sp]['alleleFreq'][a]</tt>.

	   \param alleleFreq_param a dictionary of parameters of \c alleleFreq statistics.
	   Can be one or more items choosen from the following options: \c numOfAlleles,
	   \c alleleNum, and \c alleleFreq.

	   \param genoFreq an array of loci at which all genotype frequencies will be
	   calculated (<tt>genoFreq=[loc1, loc2, ...]</tt>. You may use parameter
	   \c genoFreq_param to control if <tt>a/b</tt> and <tt>b/a</tt> are the same
	   genotype. This parameter will set the following
	   dictionary variables. Note that unlike list used for \c alleleFreq etc.,
	   the indexes \c a, \c b of <tt>genoFreq[loc][a][b]</tt> are dictionary keys,
	   so you will get a \em KeyError when you used a wrong key. You can get around
	   this problem by using expressions like <tt>genoNum[loc].setDefault(a,{})</tt>.
	   \li <tt>genoNum[loc][allele1][allele2]</tt> and <tt>subPop[sp]['genoNum'][loc][allele1][allele2]</tt>,
	   the number of genotype \c allele1-allele2 at locus \c loc.
	   \li <tt>genoFreq[loc][allele1][allele2]</tt> and <tt>subPop[sp]['genoFreq'][loc][allele1][allele2]</tt>,
	   the frequency of genotype \c allele1-allele2 at locus \c loc.
	   \li genoFreq_param a dictionary of parameters of \c phase = 0 or 1.

	   \param heteroFreq an array of loci at which observed heterozygosities will be calculated
	   (<tt>heteroFreq=[loc1, loc2, ...]</tt>). For each locus, the number and frequency of
	   allele specific and overall heterozygotes will be calcuated and stored in four population
	   variables. For example, <tt>heteroNum[loc][1]</tt> stores number of heterozygotes
	   at locus \c loc, with respect to allele \c 1, which is the number of all genotype
	   \c 1x or \c x1 where \x does not equal to \c 1. All other genotypes such as \c 02 are
	   considered as homozygotes when <tt>heteroFreq[loc][1]</tt> is calculated.
	   The overall number of heterozygotes (<tt>HeteroNum[loc]</tt>) is the number of
	   genotype \c xy if \c x does not equal to \c y.
	   \li <tt>HeteroNum[loc]</tt>, <tt>subPop[sp]['HeteroNum'][loc]</tt>, the overall heterozygote count.
	   \li <tt>HeteroFreq[loc]</tt>, <tt>subPop[sp]['HeteroFreq'][loc]</tt>, the overall heterozygote frequency.
	   \li <tt>heteroNum[loc][allele]</tt>, <tt>subPop[sp]['heteroNum'][loc][allele]</tt>, allele-specific heterozygote counts.
	   \li <tt>heteroFreq[loc][allele]</tt>, <tt>subPop[sp]['heteroFreq'][loc][allele]</tt>, allele-specific heterozygote frequency.

	   \param homoFreq an array of loci to calculate observed homozygosities
	   and expected homozygosities (<tt>homoFreq=[loc1, loc2, ...]</tt>).
	   This parameter will calculate the numbers and frequencies of homozygotes
	   \b xx and set the following variables:
	   \li <tt>homoNum[loc]</tt>, <tt>subPop[sp]['homoNum'][loc]</tt>.
	   \li <tt>homoFreq[loc]</tt>, <tt>subPop[sp]['homoFreq'][loc]</tt>.

	   \param expHetero an array of loci at which the expected heterozygosities will
	   be calculated (<tt>expHetero=[loc1, loc2, ...]</tt>). The expected heterozygosity
	   is calculated by \f[ h_{exp}=1-p_{i}^{2}, \f] where \f$ p_i \f$ is the allele
	   frequency of allele \f$ i \f$. The following variables will be set:
	   \li <tt>expHetero[loc]</tt>, <tt>subPop[sp]['expHetero'][loc]</tt>.

	   \param expHetero_param a dictionary of parameters of \c expHetero statistics.
	   Can be one or more items choosen from the following options: \c subpop and
	   \c midValues.

	   \param haploFreq a matrix of haplotypes (allele sequences on different loci) to
	   count. For example, <tt>haploFreq = [ [ 0,1,2 ], [1,2] ]</tt>  will count
	   all haplotypes on loci 0, 1 and 2; and all haplotypes on loci 1, 2.
	   If only one haplotype is specified, the outer <tt>[]</tt> can be omitted. I.e.,
	   <tt>haploFreq=[0,1]</tt> is acceptable. The following dictionary variables
	   will be set with keys <tt>0-1-2</tt> etc. For example, <tt>haploNum['1-2']['5-6']</tt>
	   is the number of allele pair 5, 6 (on loci 1 and 2 respectively) in the population.
	   \li <tt>haploNum[haplo]</tt> and <tt>subPop[sp]['haploNum'][haplo]</tt>, the number
	   of allele sequencies on loci \c haplo.
	   \li <tt>haploFreq[haplo]</tt>, <tt>subPop[sp]['haploFreq'][haplo]</tt>, the frequency
	   of allele sequencies on loci \c haplo.

	   \param LD calculate linkage disequilibria \f$ LD \f$, \f$ LD' \f$
	   and \f$ r^{2} \f$, given
	   <tt>LD=[ [loc1, loc2], [ loc1, loc2, allele1, allele2], ... ]</tt>.
	   For each item <tt>[loc1, loc2, allele1, allele2]</tt>, \f$ D \f$, \f$ D' \f$
	   and \f$ r^{2} \f$ will be calculated
	   based on \c allele1 at \c loc1 and \c allele2 at \c loc2. If only two loci are given,
	   the LD values are averaged over all allele pairs. For example, for allele \f$ A \f$ at
	   locus \c 1 and allele \f$ B \f$ at locus \c 2,
	   \f[ D = P_{AB}-P_{A}P_{B} \f]
	   \f[ D' = D/D_{max} \f]
	   \f[ D_{max} =
	   \min\left(P_{A}\left(1-P_{B}\right),\left(1-P_{A}\right)P_{B}\right)  \textrm{if }D>0 \ \
	   \min\left(P_{A}P_{B},\left(1-P_{A}\right)\left(1-P_{B}\right)\right)  \textrm{if }D<0 \f]
	   \f[ r^{2} = \frac{D^{2}}{P_{A}\left(1-P_{A}\right)P_{B}\left(1-P_{B}\right)} \f]
	   If only one item is specified, the
	   outer <tt>[]</tt> can be ignored. I.e., <tt>LD=[loc1, loc2]</tt> is acceptable.
	   This parameter will set the following variables. Please note that the difference between
	   the data structures used for \c ld and \c LD.
	   \li <tt>ld['loc1-loc2']['allele1-allele2']</tt>, <tt>subPop[sp]['ld']['loc1-loc2']['allele1-allele2']</tt>.
	   \li <tt>ld_prime['loc1-loc2']['allele1-allele2']</tt>, <tt>subPop[sp]['ld_prime']['loc1-loc2']['allele1-allele2']</tt>.
	   \li <tt>r2['loc1-loc2']['allele1-allele2']</tt>, <tt>subPop[sp]['r2']['loc1-loc2']['allele1-allele2']</tt>.
	   \li <tt>LD[loc1][loc2]</tt>, <tt>subPop[sp]['LD'][loc1][loc2]</tt>.
	   \li <tt>LD_prime[loc1][loc2]</tt>, <tt>subPop[sp]['LD_prime'][loc1][loc2]</tt>.
	   \li <tt>R2[loc1][loc2]</tt>, <tt>subPop[sp]['R2'][loc1][loc2]</tt>.

	   \param LD_param a dictionary of parameters of \c LD statistics. Can have key \c stat which is
	   a list of statistics to calculate. Default to all. If any statistics is specified,
	   only those specified will be calculated. For example, you may use <tt>LD_param={LD_prime}</tt>
	   to calculate D' only, where <tt>LD_prime</tt> is a shortcut for <tt>'stat':['LD_prime']</tt>.
	   Other parameters that you may use are:
	   \li \c subPop whether or not calculate statistics for subpopulations.
	   \li \c midValues whether or not keep intermediate results.

	   \param association association measures

	   \param association_param a dictionary of parameters of \c association statistics.
	   Can be one or more items choosen from the following options: \c ChiSq,
	   \c ChiSq_P, \c UC_U, and \c CramerV.

	   \param Fst calculate \f$ F_{st} \f$, \f$ F_{is} \f$, \f$ F_{it} \f$.
	   For example, <tt>Fst = [0,1,2]</tt> will calculate \f$ F_{st} \f$, \f$ F_{is} \f$,
	   \f$ F_{it} \f$ based on alleles at loci \c 0, \c 1, \c 2. The locus-specific values will be used
	   to calculate \c AvgFst, which is an average value over all alleles (Weir
	   & Cockerham, 1984). Terms and values that match Weir & Cockerham are:
	   \li \f$ F \f$ (\f$ F_{IT} \f$) the correlation of genes within individuals (inbreeding);
	   \li \f$ \theta \f$ (\f$ F_{ST} \f$) the correlation of genes of difference individuals
	   in the same population (will evaluate for each subpopulation and the whole population)
	   \li \f$ f \f$ (\f$ F_{IS} \f$) the correlation of genes within individuals within
	   populations.

	   This parameter will set the following variables:
	   \li <tt>Fst[loc]</tt>, <tt>Fis[loc]</tt>, <tt>Fit[loc]</tt>
	   \li <tt>AvgFst</tt>, <tt>AvgFis</tt>, <tt>AvgFit</tt>.

	   \param Fst_param a dictionary of parameters of \c Fst statistics.
	   Can be one or more items choosen from the following options: \c Fst,
	   \c Fis, \c Fit, \c AvgFst, \c AvgFis, and \c AvgFit.

	   \param relMethod method used to calculate relatedness. Can be either \c  REL_Queller or \c REL_Lynch.
	   The relatedness values between two individuals, or two groups of individuals are calculated according
	   to Queller & Goodnight (1989) (<tt>method=REL_Queller</tt>) and Lynch et al. (1999) (<tt>method=REL_Lynch</tt>).
	   The results are pairwise relatedness values, in the form of a matrix. Original group or subpopulation
	   numbers are discarded. There is no subpopulation level relatedness value.

	   \param relGroups calculate pairwise relatedness between groups. Can be in the form of either
	   <tt>[[1,2,3],[5,6,7],[8,9]]</tt> or <tt>[2,3,4]</tt>. The first one specifies groups of
	   individuals, while the second specifies subpopulations. By default, relatedness between
	   subpopulations is calculated.

	   \param relLoci loci on which relatedness values are calculated

	   \param rel_param a dictionary of parameters of relatedness statistics.
	   Can be one or more items choosen from the following options: \c Fst,
	   \c Fis, \c Fit, \c AvgFst, \c AvgFis, and \c AvgFit.

	   \param hasPhase if a/b and b/a are the same genotype. Default to \c False.

	   \param midValues whether or not post intermediate results. Default to \c False.
	   For example, \c Fst will need to calculate allele frequencise. If \c midValues
	   is set to \c True, allele frequencies will be posted as well. This will be
	   helpful in debugging and sometimes in deriving statistics.
	 **/
	stat(bool popSize = false,
		//
		bool numOfMale = false,
		strDict numOfMale_param = strDict(),
		//
		bool numOfAffected = false,
		strDict numOfAffected_param = strDict(),
		//
		vectori numOfAlleles = vectori(),
		strDict numOfAlleles_param = strDict(),
		//
		vectori alleleFreq = vectori(),
		strDict alleleFreq_param = strDict(),
		//
		vectori heteroFreq = vectori(),
		vectori expHetero = vectori(),
		strDict expHetero_param = strDict(),
		//
		vectori homoFreq = vectori(),
		vectori genoFreq = vectori(),
		strDict genoFreq_param = strDict(),
		intMatrix haploFreq = intMatrix(),
		//
		intMatrix LD = intMatrix(),
		strDict LD_param = strDict(),
		//
		intMatrix association = intMatrix(),
		strDict association_param = strDict(),
		//
		vectori Fst = vectori(),
		strDict Fst_param = strDict(),
		//
		intMatrix relGroups = intMatrix(),
		vectori relLoci = vectori(),
		strDict rel_param = strDict(),
		//
		bool relBySubPop = false,                                           // internal use
		vectori relMethod = vectori(),
		int relMinScored = 10,                                              // minimal number of loci required.
		bool hasPhase = false,
		bool midValues = false,                                             // this parameter will be removed after all _param parameter is given.
		// regular parameters
		string output = "",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPops = subPopList(),
		const vectorstr & infoFields = vectorstr());

	~stat()
	{
	}


	/// deep copy of a \c stat operator
	virtual baseOperator * clone() const
	{
		return new stat(*this);
	}


	// count various statistics.
	// use m_alleles etc to save (potentially) time to
	// resize all these variables.
	/// apply the \c stat operator
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c stat operator
	virtual string __repr__()
	{
		return "<simuPOP::statistics>";
	}


private:
	statPopSize m_popSize;
	statNumOfMale m_numOfMale;
	statNumOfAffected m_numOfAffected;
	statAlleleFreq m_alleleFreq;
	statNumOfAlleles m_numOfAlleles;
	statHeteroFreq m_heteroFreq;
	statExpHetero m_expHetero;
	statGenoFreq m_genoFreq;
	statHaploFreq m_haploFreq;
	statLD m_LD;
	statAssociation m_association;
	statFst m_Fst;
	statRelatedness m_relatedness;
};
}
#endif
