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
#define  popSize_String        "popSize"
#define  subPopSize_String     "subPopSize"
#define  popSize_sp_String     "popSize_sp"

public:
	statPopSize(bool popSize, const subPopList & subPops, const stringList & vars)
		: m_isActive(popSize), m_subPops(subPops), m_vars()
	{
		const char * allowedVars[] = { popSize_String,    popSize_sp_String,
			                           subPopSize_String, "" };
		const char * defaultVars[] = { popSize_String, subPopSize_String, "" };

		m_vars.obtainFrom(vars, allowedVars, defaultVars);
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
#define  numOfMale_String         "numOfMale"
#define  propOfMale_String        "propOfMale"
#define  numOfFemale_String       "numOfFemale"
#define  propOfFemale_String      "propOfFemale"
#define  numOfMale_sp_String      "numOfMale_sp"
#define  propOfMale_sp_String     "propOfMale_sp"
#define  numOfFemale_sp_String    "numOfFemale_sp"
#define  propOfFemale_sp_String   "propOfFemale_sp"

public:
	statNumOfMale(bool numOfMale, const subPopList & subPops, const stringList & vars)
		: m_isActive(numOfMale), m_subPops(subPops), m_vars()
	{
		const char * allowedVars[] = {
			numOfMale_String,      propOfMale_String,
			numOfFemale_String,    propOfFemale_String,
			numOfMale_sp_String,   propOfMale_sp_String,
			numOfFemale_sp_String, propOfFemale_sp_String,""
		};
		const char * defaultVars[] = { numOfMale_String, numOfFemale_String, "" };

		m_vars.obtainFrom(vars, allowedVars, defaultVars);
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
#define  numOfAffected_String        "numOfAffected"
#define  propOfAffected_String       "propOfAffected"
#define  numOfUnaffected_String      "numOfUnaffected"
#define  propOfUnaffected_String     "propOfUnaffected"
#define  numOfAffected_sp_String     "numOfAffected_sp"
#define  propOfAffected_sp_String    "propOfAffected_sp"
#define  numOfUnaffected_sp_String   "numOfUnaffected_sp"
#define  propOfUnaffected_sp_String  "propOfUnaffected_sp"

public:
	statNumOfAffected(bool numOfAffected, const subPopList & subPops, const stringList & vars)
		: m_isActive(numOfAffected), m_subPops(subPops), m_vars()
	{
		const char * allowedVars[] = {
			numOfAffected_String,      propOfAffected_String,
			numOfUnaffected_String,    propOfUnaffected_String,
			numOfAffected_sp_String,   propOfAffected_sp_String,
			numOfUnaffected_sp_String, propOfUnaffected_sp_String,""
		};
		const char * defaultVars[] = { numOfAffected_String, numOfUnaffected_String, "" };

		m_vars.obtainFrom(vars, allowedVars, defaultVars);
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
#define  AlleleNum_String        "alleleNum"
#define  AlleleFreq_String       "alleleFreq"
#define  AlleleNum_sp_String     "alleleNum_sp"
#define  AlleleFreq_sp_String    "alleleFreq_sp"

public:
	statAlleleFreq(const vectorlu & loci, const subPopList & subPops, const stringList & vars)
		: m_loci(loci), m_subPops(subPops), m_vars()
	{
		const char * allowedVars[] = {
			AlleleNum_String,    AlleleFreq_String,
			AlleleNum_sp_String, AlleleFreq_sp_String,""
		};
		const char * defaultVars[] = { AlleleFreq_String, AlleleNum_String, "" };

		m_vars.obtainFrom(vars, allowedVars, defaultVars);
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
#define HomoNum_String          "homoNum"
#define HomoFreq_String         "homoFreq"
#define HeteroNum_sp_String     "heteroNum_sp"
#define HeteroFreq_sp_String    "heteroFreq_sp"
#define HomoNum_sp_String       "homoNum_sp"
#define HomoFreq_sp_String      "homoFreq_sp"

public:
	statHeteroFreq(const vectorlu & heteroFreq, const vectorlu & homoFreq,
		const subPopList & subPops, const stringList & vars);

	void addLocus(UINT locus, const subPopList & subPops = AllSubPops,
		const stringList & vars = stringList());

	double heteroFreq(population & pop, UINT allele, int loc, vspID subPop);

	bool apply(population & pop);

private:
	/// heteroFreq
	vectorlu m_loci;

	subPopList m_subPops;
	stringList m_vars;
};


/// CPPONLY
class statGenoFreq
{
private:
#define  GenotypeNum_String      "genoNum"
#define  GenotypeFreq_String     "genoFreq"
#define  GenotypeNum_sp_String   "genoNum_sp"
#define  GenotypeFreq_sp_String  "genoFreq_sp"

public:
	statGenoFreq(const vectorlu & genoFreq,
		const subPopList & subPops, const stringList & vars);

	// Return AA, Aa and aa, wild type is A.
	vectorlu countGenotype(population & pop, UINT loc, UINT wildtype);

	vectorlu countGenotype(population & pop, UINT loc, SubPopID subPop, UINT wildtype);

	bool apply(population & pop);

private:
	/// which genotypes
	vectorlu m_loci;

	subPopList m_subPops;
	stringList m_vars;
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
		m_loci(rhs.m_loci), m_Fst(rhs.m_Fst),
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
	vectorlu m_loci;

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
	 *  set of statistics for all subpopulations (<tt>subPops=AllSubPops</tt>).
	 *  If a list of (virtual) subpopulations are specified, statistics for
	 *  only specified subpopulations will be calculated. However, different
	 *  statistics treat this parameter differently and it is very important
	 *  to check its reference before you use \e subPops for any statistics.
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
	 *  variable \c alleleFreq (<tt>dvars().alleleFreq</tt>) for all or
	 *  specified subpopulations. If this does not fit your need, you can
	 *  use parameter \e vars to output additional parameters, or limit the
	 *  output of existing parameters. More specifically, for this particular
	 *  statistic, the available variables are \c 'alleleFreq', \c 'alleleNum',
	 *  \c 'alleleFreq_sp' (\c 'alleleFreq' in each subpopulation), and
	 *  \c 'alleleNum_sp' (\c 'alleleNum' in each subpopulation). You can set
	 *  <tt>vars=['alleleNum_sp']</tt> to output only subpopulation specific
	 *  allele count.
	 *
	 *  Operator \c stat supports the following statistics:
	 *
	 *  <b>popSize</b>: If \e popSize=True, number of individuals in all or
	 *  specified subpopulations (parameter \e subPops) will be set to the
	 *  following variables:
	 *  \li \c popSize (default): Number of individuals in all or specified
	 *       subpopulations. Because \e subPops does not have to cover all
	 *       individuals, it may not be the actual population size.
	 *  \li \c popSize_sp: Size of (virtual) subpopulation \c sp.
	 *  \li \c subPopSize (default): A list of subpopulation sizes.
	 *       <tt>sum(subPopSize)</tt> is the total population size.
	 *
	 *  <b>numOfMale</b>: If \e numOfMale=True, number of male individuals in
	 *  all or specified (virtual) subpopulations will be set to the following
	 *  variables:
	 *  \li \c numOfMale (default): Total number of male individuals in all
	 *       or specified (virtual) subpopulations.
	 *  \li \c numOfMale (default): Total number of female individuals in all
	 *       or specified (virtual) subpopulations.
	 *  \li \c propOfMale: Proportion of male individuals.
	 *  \li \c propOfFemale: Proportion of female individuals.
	 *  \li \c numOfMale_sp: Number of male individuals in each (virtual)
	 *       subpopulation.
	 *  \li \c numOfFemale_sp: Number of female individuals in each (virtual)
	 *       subpopulation.
	 *  \li \c propOfMale_sp: Proportion of male individuals in each (virtual)
	 *       subpopulation.
	 *  \li \c propOfFemale_sp: Proportion of female individuals in each
	 *       (virtual) subpopulation.
	 *
	 *  <b>numOfAffected</b>: If \e numOfAffected=True, number of affected
	 *  individuals in all or specified (virtual) subpopulations will be set to the
	 *  following variables:
	 *  \li \c numOfAffected (default): Total number of affected individuals in
	 *       all or specified (virtual) subpopulations.
	 *  \li \c numOfAffected (default): Total number of unaffected individuals
	 *       in all or specified (virtual) subpopulations.
	 *  \li \c propOfAffected: Proportion of affected individuals.
	 *  \li \c propOfUnaffected: Proportion of unaffected individuals.
	 *  \li \c numOfAffected_sp: Number of affected individuals in each (virtual)
	 *       subpopulation.
	 *  \li \c numOfUnaffected_sp: Number of unaffected individuals in each
	 *       (virtual) subpopulation.
	 *  \li \c propOfAffected_sp: Proportion of affected individuals in each
	 *       (virtual) subpopulation.
	 *  \li \c propOfUnaffected_sp: Proportion of unaffected individuals in
	 *       each (virtual) subpopulation.
	 *
	 *  <b>alleleFreq</b>: This parameter accepts a list of loci (by indexes),
	 *  at which allele frequencies will be calculated. This statistic outputs
	 *  the following variables:
	 *  \li \c alleleFreq (default): <tt>alleleFreq[loc][a]</tt> is the
	 *       frequency of allele \c a at locus \loc for all or specified
	 *       (virtual) subpopulations. Variable \c alleleFreq is a dictionary
	 *       (with loci indexes as keys) of lists. The length of the frequency
	 *       list is determined by the maximum allele, but is guaranteed to
	 *       have at least 2 items (frequency for alleles 0 and 1).
	 *  \li \c alleleNum (default): <tt>alleleNum[loc][a]</tt> is the number of
	 *       allele \c a at locus \loc for all or specified (virtual)
	 *       subpopulations.
	 *  \li \c alleleFreq_sp: Allele frequency in each (virtual) subpopulation.
	 *  \li \c alleleNum_sp: Allele count in each (virtual) subpopulation.
	 *
	 *  <b>heteroFreq</b> and <b>homoFreq</b>: These parameters accept a list
	 *  of loci (by indexes), at which the number and frequency of homozygotes
	 *  and/or heterozygotes will be calculated. These statistics are only
	 *  available for diploid populations. The following variables will be
	 *  outputted:
	 *  \li \c heteroFreq (default for parameter \e heteroFreq): A dictionary
	 *       of proportion of heterozygotes in all or specified (virtual)
	 *       subpopulations, with loci indexes as dictionary keys.
	 *  \li \c homoFreq (default for parameter \e homoFreq): A dictionary of
	 *       proportion of homozygotes in all or specified (virtual)
	 *       subpopulations.
	 *  \li \c heteroNum: A dictionary of number of heterozygotes in all or
	 *       specified (virtual) subpopulations.
	 *  \li \c homoNum: A dictionary of number of homozygotes in all or
	 *       specified (virtual) subpopulations.
	 *  \li \c heteroFreq_sp: A dictionary of proportion of heterozygotes in
	 *       each (virtual) subpopulation.
	 *  \li \c homoFreq_sp: A dictionary of proportion of homozygotes in each
	 *       (virtual) subpopulation.
	 *  \li \c heteroNum_sp: A dictionary of number of heterozygotes in each
	 *       (virtual) subpopulation.
	 *  \li \c homoNum_sp: A dictionary of number of homozygotes in each
	 *       (virtual) subpopulation.
	 *
	 *  <b>genoFreq</b>: This parameter accept a list of loci (by index) at
	 *  which number and frequency of all genotypes are outputed as a
	 *  dictionary (indexed by loci indexes) of dictionaries (indexed by tuples
	 *  of possible indexes). This statistic is available for all population
	 *  types with genotype defined as ordered alleles at a locus. The length of
	 *  genotype is determined by population ploidy. Because genotypes are
	 *  ordered, <tt>(1, 0)</tt> and <tt>(0, 1)</tt> (two possible genotypes in
	 *  a diploid population) are considered as different genotypes. This
	 *  statistic outputs the following variables:
	 *  \li \c genoFreq (default): A dictionary (by loci indexes) of
	 *       dictionaries (by genotype) of genotype frequencies. For example,
	 *       <tt>genoFreq[1][(1, 0)]</tt> is the frequency of genotype (1, 0)
	 *       at locus 1.
	 *  \li \c genoNum (default): A dictionary of dictionaries of genotype
	 *       counts of all or specified (virtual) subpopulations.
	 *  \li \c genoFreq_sp: genotype frequency in each specified (virtual)
	 *      subpopulation.
	 *  \li \c genoFreq_sp: genotype count in each specified (virtual)
	 *      subpopulation.
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
		const uintList & homoFreq = uintList(),
		//
		const uintList & genoFreq = uintList(),
		//
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
