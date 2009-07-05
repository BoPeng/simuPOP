/**
 *  $File: stator.h $
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


/** A \c pyEval operator evaluates a Python expression in a population's local
 *  namespace when it is applied to this population. The result is written to
 *  an output specified by parameter \e output.
 *  <funcForm>PyEval</funcForm>
 */
class pyEval : public baseOperator
{
public:
	/** Crete a \c pyEval operator that evaluates a Python expression \e expr
	 *  in a population's local namespace when it is applied to this population.
	 *  If Python statements \e stmts is given (a single or multi-line string),
	 *  the statement will be executed before \e expr. If \e exposePop is set
	 *  to an non-empty string, the current population will be exposed in its
	 *  own local namespace as a variable with this name. This allows the
	 *  execution of expressions such as <tt>'pop.individual(0).allele(0)'</tt>.
	 *  The result of \e expr will be sent to an output stream specified by
	 *  parameter \c output. The exposed population variable will be removed
	 *  after \e expr is evaluated. Please refer to class \c baseOperator for
	 *  other parameters.
	 *
	 *  \note Although the statements and expressions are evaluated in a
	 *  population's local namespace, they have access to a **global**
	 *  namespace which is the module global namespace. It is therefore
	 *  possible to refer to any module variable in these expressions. Such
	 *  mixed use of local and global variables is, however, strongly
	 *  discouraged.
	 */
	pyEval(const string & expr = string(), const string & stmts = string(),
		const string & exposePop = string(), const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops, const stringList & infoFields = stringList())
		: baseOperator(output, stage, begin, end, step, at, reps, subPops, infoFields),
		m_expr(expr, stmts), m_exposePop(exposePop)
	{
	}


	~pyEval()
	{
	}


	/// deep copy of a \c pyEval operator
	virtual baseOperator * clone() const
	{
		return new pyEval(*this);
	}


	/** Evaluate the expression and optional statements in the local namespace
	 *  of population \e pop and return its result as a string.
	 */
	string evaluate(population & pop);

	/// Apply the \c pyEval operator to population \e pop.
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c pyEval operator
	virtual string __repr__()
	{
		return "<simuPOP::pyEval>";
	}


private:
	/// expression to evaluate
	Expression m_expr;

	/// if expose pop
	string m_exposePop;
};


/** This operator executes given Python statements in a population's local
 *  namespace when it is applied to this population.
 *  <funcForm>PyExec</funcForm>
 */
class pyExec : public pyEval
{
public:
	/** Create a \c pyExec operator that executes statements \e stmts in a
	 *  population's local namespace when it is applied to this population.
	 *  If \e exposePop is given, current population will be exposed in
	 *  its local namespace as a variable named by \e exposePop. Although
	 *  multiple statements can be executed, it is recommended that you use
	 *  this operator to execute short statements and use \c pyOperator for
	 *  more complex once. Note that exposed population variable will be
	 *  removed after the statements are executed.
	 */
	pyExec(const string & stmts = string(), const string & exposePop = string(),
		const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops, const stringList & infoFields = stringList())
		: pyEval("", stmts, exposePop, "", stage, begin, end, step, at, reps, subPops, infoFields)
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
		return "<simuPOP::pyExec>";
	}


};


/** Unlike operator \c pyEval and \c pyExec that work at the population level,
 *  in a population's local namespace, operator \c infoEval works at the
 *  individual level, working with individual information fields. When this
 *  operator is applied to a population, information fields of eligible
 *  individuals are put into either a temporary dictionary or in the local
 *  namespace of the population. A Python expression is then evaluated for
 *  each individual. The result is written to an output.
 *
 *  \note Unlike operator ``infoExec``, individual information fields are not
 *  updated after this operator is applied to a population.
 *
 *  \note This operator tends to generate a large amount of output so use it
 *  is with caution.
 *
 *  <funcForm>InfoEval</funcForm>
 */
class infoEval : public baseOperator
{
public:
	/** Create an operator that evaluate a Python expression \e expr using
	 *  individual information fields as variables. For each eligible
	 *  individual (individuals in (virtual) subpopulations specified by
	 *  parameter \e subPops, default to all individuals), its information
	 *  fields are copied either to a temporary namespace (default) or the
	 *  population's local namespace (if \e usePopVars is \c True). If
	 *  \e exposeInd is not empty, the individual itself will be exposed in
	 *  this namespace as a variable with name specified by \e exposeInd. In
	 *  the <tt>usePopVars=True</tt> case, any population variable whose name
	 *  matches an information field or \e exposeInd will be silently
	 *  overridden.
	 *
	 *  A Python expression (\e expr) is evaluated for each individual. The
	 *  results are converted to strings and are written to an output specified
	 *  by parameter \e output. Optionally, a statement (or several statements
	 *  separated by newline) can be executed before \e expr is evaluated.
	 *
	 *  This operator is by default applied post-mating. If it stage is set to
	 *  \c DuringMating, it will be applied to all offspring, regardless of
	 *  \c subPops settings.
	 *
	 *  \note Although \e expr is evaluated in individual or population level
	 *  local namespaces, it can also access a global namespace which is
	 *  the module namespace of your script. However, using module level
	 *  variables and functions in this operator is discouraged.
	 */
	infoEval(const string & expr = string(), const string & stmts = string(), bool usePopVars = false,
		const string & exposeInd = string(),
		const stringFunc & output = ">", int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops, const stringList & infoFields = stringList())
		: baseOperator(output, stage, begin, end, step, at, reps, subPops, infoFields),
		m_expr(expr, stmts), m_usePopVars(usePopVars), m_exposeInd(exposeInd), m_dict(NULL)
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
	bool apply(population & pop);

	bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// used by Python print function to print out the general information of the \c infoEval operator
	virtual string __repr__()
	{
		return "<simuPOP::infoEval>";
	}


protected:
	string evalInfo(individual *, bool update);

	/// expression to evaluate
	Expression m_expr;

	bool m_usePopVars;

	string m_exposeInd;

	PyObject * m_dict;
};

/** Operator \c infoExec is similar to \c infoEval in that it works at the
 *  individual level, using individual information fields as variables. The
 *  difference is that instead of evaluating an expression and outputing its
 *  result, this operator execute one or more statements and <bf>update
 *  individual information fields</bf> from the namespace after the
 *  specified statements are execuated.
 *
 *  <funcForm>InfoExec</funcForm>
 */
class infoExec : public infoEval
{
public:
	/** Create an operator that executes Python statements \e stmts using
	 *  individual information fields as variables. For each eligible
	 *  individual (individuals in (virtual) subpopulations specified by
	 *  parameter \e subPops, default to all individuals), its information
	 *  fields are copied either to a temporary namespace (default) or the
	 *  population's local namespace (if \e usePopVars is \c True). If
	 *  \e exposeInd is not empty, the individual itself will be exposed in
	 *  this namespace as a variable with name specified by \e exposeInd. In
	 *  the <tt>usePopVars=True</tt> case, any population variable whose name
	 *  matches an information field or \e exposeInd will be silently
	 *  overridden.
	 *
	 *  One or more python statements (\e stmts) are executed for each
	 *  individual. Information fields of these individuals are then updated
	 *  from the corresponding variables. For example, <tt>a=1</tt> will set
	 *  information field \e a of all individuals to \c 1, <tt>a=b</tt> will
	 *  set information field \e a of all individuals to information field
	 *  \c b or a population variable \c b if \c b is not an information field
	 *  but a population variable (needs <tt>usePopVars=True</tt>), and
	 *  <tt>a=ind.sex()</tt> will set information field \e a of all individuals
	 *  to its sex (needs <tt>exposeInd='ind'</tt>.
	 *
	 *  This operator is by default applied post-mating. If it stage is set to
	 *  \c DuringMating, it will be applied to all offspring, regardless of
	 *  \c subPops settings.
	 *
	 *  \note Although \e stmts are executed in individual or population level
	 *  local namespaces, they also have access to a global namespace which is
	 *  the module namespace of your script. However, using module level
	 *  variables and functions in \e stmts is discouraged.
	 */
	infoExec(const string & stmts = string(), bool usePopVars = false,  const string & exposeInd = string(),
		const stringFunc & output = "", int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops, const stringList & infoFields = stringList())
		: infoEval(string(), stmts, usePopVars, exposeInd, output, stage, begin, end, step, at, reps, subPops, infoFields),
		m_simpleStmt(stmts)
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


	// check all alleles in vector allele if they are fixed.
	/// apply the \c infoExec operator
	bool apply(population & pop);

	bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// used by Python print function to print out the general information of the \c infoExec operator
	virtual string __repr__()
	{
		return "<simuPOP::infoExec>";
	}


private:
	simpleStmt m_simpleStmt;
};


// The following classes apply various statistics
// and stat class will provide an opearator interface
// to all of them.


/// CPPONLY return {'a-b-b'} for a b c
string haploKey(const vectori & seq);

/// CPPONLY post population sizes etc.
class statPopSize
{
private:
#define  popSize_String     "popSize"
#define  subPopSize_String  "subPopSize"

public:
	statPopSize(bool popSize, const subPopList & subPops, const stringList & vars)
		: m_isActive(popSize), m_subPops(subPops), m_vars(vars)
	{
	}

	bool apply(population & pop);

private:
	bool m_isActive;
	subPopList m_subPops;
	stringList m_vars;
};

/// CPPONLY
class statNumOfMale
{
private:
#define  numOfMale_String      "numOfMale"
#define  propOfMale_String     "propOfMale"
#define  numOfFemale_String    "numOfFemale"
#define  propOfFemale_String   "propOfFemale"

public:
	statNumOfMale(bool numOfMale, const subPopList & subPops, const stringList & vars)
		: m_isActive(numOfMale), m_subPops(subPops), m_vars(vars)
	{
	}

	bool apply(population & pop);

private:
	/// whether or not to apply number of male/female
	bool m_isActive;
	subPopList m_subPops;
	stringList m_vars;
};

/// CPPONLY
class statNumOfAffected
{
private:
#define  numOfAffected_String     "numOfAffected"
#define  propOfAffected_String    "propOfAffected"
#define  numOfUnaffected_String   "numOfUnaffected"
#define  propOfUnaffected_String  "propOfUnaffected"

public:
	statNumOfAffected(bool numOfAffected, const subPopList & subPops, const stringList & vars)
		: m_isActive(numOfAffected), m_subPops(subPops), m_vars(vars)
	{
	}

	bool apply(population & pop);

private:
	/// whether or not to apply number of affected
	bool m_isActive;
	subPopList m_subPops;
	stringList m_vars;
};


/// CPPONLY
class statAlleleFreq
{
private:
#define  AlleleNum_String     "alleleNum"
#define  AlleleFreq_String    "alleleFreq"

public:
	statAlleleFreq(const vectorlu & loci, const subPopList & subPops, const stringList & vars)
		: m_loci(loci), m_subPops(subPops), m_vars(vars)
	{
	}


	/// destructor, nested vectors have to be cleared manually
	~statAlleleFreq()
	{
	}


	void addLocus(UINT locus, const subPopList & subPops = AllSubPops,
		const stringList & vars = stringList());

	vectori numOfAlleles(population & pop);

	vectorf alleleFreqVec(population & pop, int loc);

	double alleleFreq(population & pop, UINT allele, int loc);

	vectorf alleleFreqVec(population & pop, int loc, UINT subPop);

	double alleleFreq(population & pop, UINT allele, int loc, UINT subPop);

	vectori alleles(population & pop, int loc);

	bool apply(population & pop);

private:
	/// which alleles?
	vectorlu m_loci;

	subPopList m_subPops;
	stringList m_vars;
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

	int locusIdx(UINT loc)
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
	statHeteroFreq(const vectorlu & heteroFreq = vectorlu(),
		const vectorlu & homoFreq = vectorlu())
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


	void addLocus(UINT locus, bool post = false)
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
	vectorlu m_atLoci;

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
	statExpHetero(statAlleleFreq & alleleFreq, const vectorlu & expHetero = vectorlu(),
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
			m_alleleFreq.addLocus(expHetero[i]);
	}


	/// CPPONLY: used for copy constructor
	statExpHetero(statAlleleFreq & alleleFreq, const statExpHetero & rhs)
		: m_alleleFreq(alleleFreq), m_atLoci(rhs.m_atLoci),
		m_expHetero(rhs.m_expHetero), m_midValues(rhs.m_midValues),
		m_evalInSubPop(rhs.m_evalInSubPop)
	{

	}


	bool apply(population & pop);

private:
	/// need this to apply alleleFreq
	statAlleleFreq & m_alleleFreq;

	/// heteroFreq
	vectorlu m_atLoci;

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
	statGenoFreq(const vectorlu & genoFreq = vectorlu(),
		const strDict & param = strDict());

	// Return AA, Aa and aa, wild type is A.
	vectorlu countGenotype(population & pop, UINT loc, UINT wildtype);

	vectorlu countGenotype(population & pop, UINT loc, SubPopID subPop, UINT wildtype);

	bool apply(population & pop);

private:
	/// which genotypes
	vectorlu m_atLoci;

	/// phase
	bool m_phase;
};

/// CPPONLY
class statHaploFreq
{
private:
#define  HaplotypeNum_String    "haploNum"
#define  HaplotypeFreq_String   "haploFreq"
	int haploIndex(const vectori & haplo);

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

	// association tests
#define   ChiSq_String      "LD_ChiSq"
#define   ChiSq_P_String    "LD_ChiSq_P"
#define   UCU_String        "UC_U"
#define   CramerV_String    "CramerV"

public:
	// alleleFreq and haploFreq is required to calculate LD
	// needed allele and halplotype are added to alleleFreq and haploFreq
	// objects during the initialization of statLD, as well as stat.
	// In stat::apply(), alleleFreq.apply() and haploFreq.apply()
	// is called before statLD.apply() and ensures that allele frequencies
	// are calculated when statLD needs them.
	statLD(statAlleleFreq & alleleFreq, statHaploFreq & haploFreq,
		const intMatrix & LD = intMatrix(), const strDict & LD_param = strDict());

	/// CPPONLY: for copy constructor
	statLD(statAlleleFreq & alleleFreq, statHaploFreq & haploFreq,
		const statLD & rhs) :
		m_alleleFreq(alleleFreq), m_haploFreq(haploFreq),
		m_LD(rhs.m_LD), m_midValues(rhs.m_midValues),
		m_evalInSubPop(rhs.m_evalInSubPop),
		m_output_ld(rhs.m_output_ld),
		m_output_ld_prime(rhs.m_output_ld_prime),
		m_output_r2(rhs.m_output_r2),
		m_output_delta2(rhs.m_output_delta2),
		m_output_LD(rhs.m_output_LD),
		m_output_LD_prime(rhs.m_output_LD_prime),
		m_output_R2(rhs.m_output_R2),
		m_output_Delta2(rhs.m_output_Delta2),
		m_output_ChiSq(rhs.m_output_ChiSq),
		m_output_UCU(rhs.m_output_UCU),
		m_output_CramerV(rhs.m_output_CramerV)
	{
	}


	// calculate, right now,  do not tempt to save values
	bool apply(population & pop);

private:
	// calculate single allele LD values
	void calculateLD(population & pop, const vectori & hapLoci, const vectori & hapAlleles, UINT sp, bool subPop,
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
	bool m_output_ChiSq;
	bool m_output_UCU;
	bool m_output_CramerV;
};

/// CPPONLY
class statAssociation
{
private:
#define Asso_ChiSq_String   "ChiSq"
#define Asso_ChiSq_P_String "ChiSq_P"

public:
	statAssociation(const vectorlu & loci = vectorlu(),
		const subPopList & subPops = AllSubPops);

	// calculate, right now,  do not tempt to save values
	bool apply(population & pop);

private:
	void calcChiSq(ULONG aff_0, ULONG aff_1, ULONG unaff_0, ULONG unaff_1,
		double & chisq, double & pvalue);

	void countAlleles(IndIterator & it, UINT loc,
		ULONG & aff_0, ULONG & aff_1, ULONG & unaff_0, ULONG & unaff_1);

private:
	/// Association
	vectorlu m_loci;

	subPopList m_subPops;
};

/// CPPONLY
class statNeutrality
{
private:
#define Neutra_Pi_String   "Pi"

public:
	statNeutrality(const vectorlu & loci = vectorlu(),
		const subPopList & subPops = AllSubPops);

	// calculate, right now,  do not tempt to save values
	bool apply(population & pop);

private:
	double calcPi(IndIterator & it);

private:
	/// Neutrality
	vectorlu m_loci;

	subPopList m_subPops;
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
		const vectorlu & Fst = vectorlu(), const strDict & param = strDict());

	/// CPPONLY: for copy constructor
	statFst(statAlleleFreq & alleleFreq, statHeteroFreq & heteroFreq,
		const statFst & rhs) :
		m_alleleFreq(alleleFreq), m_heteroFreq(heteroFreq),
		m_atLoci(rhs.m_atLoci), m_Fst(rhs.m_Fst),
		m_Fit(rhs.m_Fit), m_Fis(rhs.m_Fis),
		m_avgFst(rhs.m_avgFst), m_avgFit(rhs.m_avgFit), m_avgFis(rhs.m_avgFis),
		m_midValues(rhs.m_midValues),
		m_output_Fst(rhs.m_output_Fst), m_output_Fis(rhs.m_output_Fis), m_output_Fit(rhs.m_output_Fit),
		m_output_AvgFst(rhs.m_output_AvgFst), m_output_AvgFis(rhs.m_output_AvgFis),
		m_output_AvgFit(rhs.m_output_AvgFit)
	{
	}


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
	vectorlu m_atLoci;

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


/// CPPONLY
class statHWE
{
private:
#define  HWE_String  "HWE"

public:
	statHWE(statGenoFreq & genoFreq, const vectorlu & loci = vectorlu());

	/// CPPONLY: for copy constructor
	statHWE(statGenoFreq & genoFreq, const statHWE & rhs) :
		m_genoFreq(genoFreq), m_loci(rhs.m_loci)
	{
	}


	bool apply(population & pop);

private:
	double calcHWE(const vectorlu & cnt);

	statGenoFreq & m_genoFreq;

	vectorlu m_loci;
};


/** Operator \c stat calculates various statistics of the population being
 *  applied and sets variables in its local namespace. Other operators or
 *  functions can retrieve results from or evalulate expressions in this local
 *  namespace after \c stat is applied.
 *  <funcForm>Stat</funcForm>
 */
class stat : public baseOperator
{
public:
	/** Create a \c stat operator that calculates specified statistics of a
	 *  population when it is applied to this population. This operator is
	 *  by default applied after mating (parameter \e stage) and can be applied
	 *  to specified replicates (parameter \e rep) at specified generations
	 *  (parameter \e begin, \e end, \e step, and \e at). This operator does
	 *  not produce any output (ignore parameter \e output) after statistics
	 *  are calculated. Instead, it stores results in the local namespace of
	 *  the population being applied. Other operators can retrieve these
	 *  variables or evalulate expression directly in this local namespace.
	 *  Please refer to operator \c baseOperator for a detailed explanation of
	 *  these common operator parameters. 
	 *
	 *  \c stat supports parameter \e subPops. It usually calculate the same
	 *  set of statistics for all subpopulations (<tt>subPops=[]</tt>). If
	 *  a list of (virtual) subpopulations are specified, statistics for only
	 *  specified subpopulations will be calculated. If an invalid list of
	 *  subpopulations is given (<tt>subPops=-1</tt>), statistics will not
	 *  be calculated for any subpopulation. However, different statistics
	 *  treat this parameter differently and it is very important to check its
	 *  reference before you use \e subPops for any statistics.
	 *
	 *  Calculated statistics are saved as variables in a population's local
	 *  namespace. These variables can be numbers, lists or dictionaries and
	 *  can be retrieved using functions <tt>population.vars()</tt> or
	 *  <tt>population.dvars()</tt>. If the same variable is calculated for
	 *  one or more (virtual) subpopulation, the variables are stored in
	 *  <tt>vars()['subPop'][sp]['var']</tt> where sp is a subpopulation ID
	 *  (\c sp) or a tuple of virtual subpopulation ID (<tt>(sp, vsp)</tt>).
	 *  <tt>population.vars(sp)</tt> and <tt>population.dvars(sp)</tt> provide
	 *  shortcuts to these variables.
	 *
	 *  Operator \e stat outputs a number of most useful variables for each
	 *  type of statistic. For example, <tt>alleleFreq</tt> calculates both
	 *  allele counts and allele frequencies and it by default sets
	 *  variable \c alleleFreq for each (virtual) subpopulation
	 *  (<tt>dvars(sp).alleleFreq<tt>)and for all subpopulations
	 *  (<tt>dvars().alleleFreq</tt>). If this does not fit your need, you can
	 *  use parameter \e vars to output additional parameters, or limit the
	 *  output of existing parameters. More specifically, for this particular
	 *  statistic, the available variables are \c 'alleleFreq', \c 'alleleNum',
	 *  \c 'alleleFreq_sp' (\c 'alleleFreq' in each subpopulation), and
	 *  \c 'alleleNum_sp' (\c 'alleleNum' in each subpopulation). You can set
	 *  <tt>vars=['alleleNum']</tt> to only output overall allele count.

	 Statistic calculator can output a large number of variables caParameter \e vars 
	 *  Operator \c stat supports the following statistics:
	 *
	 *  <b>Population size</b>: If \e popSize=True, population size of all or
	 *  specified subpopulations (parameter \e subPops) will be set to the
	 *  following variables:
	 *  \li \c popSize: Number of individuals in all or specified
	 *       subpopulations. Because \e subPops does not have to cover all
	 *       individuals, it may not be the actual population size.
	 *  \li \c popSize_sp: Size of (virtual) subpopulation \c sp.
	 *  \li \c subPopSize: A list of subpopulation sizes.
	 *       <tt>sum(subPopSize)</tt> is the total population size.
	 *
	 *

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
	   \li <tt>LD_ChiSq[loc1][loc2]</tt>, <tt>subPop[s]['LD_ChiSq'][loc1][loc2]</tt>
	   \li <tt>LD_ChiSq_P[loc1][loc2]</tt>, <tt>subPop[s]['LD_ChiSq_P'][loc1][loc2]</tt>
	   \li <tt>LD_UC_U[loc1][loc2]</tt>, <tt>subPop[s]['LD_UC_U'][loc1][loc2]</tt>

	   \param LD_param a dictionary of parameters of \c LD statistics. Can have key \c stat which is
	   a list of statistics to calculate. Default to LD, D' and R2. If any statistics is specified,
	   only those specified will be calculated. For example, you may use <tt>LD_param={LD_prime}</tt>
	   to calculate D' only, where <tt>LD_prime</tt> is a shortcut for <tt>'stat':['LD_prime']</tt>.
	   Other parameters that you may use are:
	   \li \c subPop whether or not calculate statistics for subpopulations.
	   \li \c midValues whether or not keep intermediate results.

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

	   \param hasPhase if a/b and b/a are the same genotype. Default to \c False.

	   \param midValues whether or not post intermediate results. Default to \c False.
	   For example, \c Fst will need to calculate allele frequencise. If \c midValues
	   is set to \c True, allele frequencies will be posted as well. This will be
	   helpful in debugging and sometimes in deriving statistics.
	 **/
	stat(bool popSize = false,
		//
		bool numOfMale = false,
		//
		bool numOfAffected = false,
		//
		const uintList & alleleFreq = uintList(),
		//
		const uintList & heteroFreq = uintList(),
		const uintList & expHetero = uintList(),
		const strDict & expHetero_param = strDict(),
		const uintList & homoFreq = uintList(),
		//
		const uintList & genoFreq = uintList(),
		const strDict & genoFreq_param = strDict(),
		const intMatrix & haploFreq = intMatrix(),
		//
		const intMatrix & LD = intMatrix(),
		const strDict & LD_param = strDict(),
		//
		const uintList & association = uintList(),
		//const strDict & association_param = strDict(),
		//
		const uintList & neutrality = uintList(),
		//
		const uintList & Fst = uintList(),
		const strDict & Fst_param = strDict(),
		//
		const uintList & HWE = uintList(),
		//
		const stringList & vars = stringList(),
		// regular parameters
		const stringFunc & output = "",
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & reps = AllReps, const subPopList & subPops = AllSubPops,
		const stringList & infoFields = stringList());

	~stat()
	{
	}


	/// CPPONLY
	stat(const stat & rhs);

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
	statHeteroFreq m_heteroFreq;
	statExpHetero m_expHetero;
	statGenoFreq m_genoFreq;
	statHaploFreq m_haploFreq;
	statLD m_LD;
	statAssociation m_association;
	statNeutrality m_neutrality;
	statFst m_Fst;
	statHWE m_HWE;
};

}
#endif
