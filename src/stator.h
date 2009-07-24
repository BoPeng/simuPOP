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


/// CPPONLY post population sizes etc.
class statPopSize
{
private:
#define  popSize_String        "popSize"
#define  subPopSize_String     "subPopSize"
#define  popSize_sp_String     "popSize_sp"

public:
	statPopSize(bool popSize, const subPopList & subPops, const stringList & vars);

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
	statNumOfMale(bool numOfMale, const subPopList & subPops, const stringList & vars);

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
	statNumOfAffected(bool numOfAffected, const subPopList & subPops, const stringList & vars);

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
	statAlleleFreq(const vectorlu & loci, const subPopList & subPops, const stringList & vars);

	/// destructor, nested vectors have to be cleared manually
	~statAlleleFreq()
	{
	}


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
#define  HaplotypeNum_String       "haploNum"
#define  HaplotypeFreq_String      "haploFreq"
#define  HaplotypeNum_sp_String    "haploNum_sp"
#define  HaplotypeFreq_sp_String   "haploFreq_sp"

public:
	statHaploFreq(const intMatrix & haploFreq, const subPopList & subPops, const stringList & vars);

	~statHaploFreq()
	{
	}


	bool apply(population & pop);

private:
	// key string (in the format of a tuple)
	string dictKey(const vectori & loci);

private:
	/// haplotype at which loci
	intMatrix m_loci;

	subPopList m_subPops;
	stringList m_vars;
};


/// CPPONLY
class statInfo
{
private:
#define  SumOfInfo_String         "sumOfInfo"
#define  MeanOfInfo_String        "meanOfInfo"
#define  VarOfInfo_String         "varOfInfo"
#define  MaxOfInfo_String         "maxOfInfo"
#define  MinOfInfo_String         "minOfInfo"
#define  SumOfInfo_sp_String      "sumOfInfo_sp"
#define  MeanOfInfo_sp_String     "meanOfInfo_sp"
#define  VarOfInfo_sp_String      "varOfInfo_sp"
#define  MaxOfInfo_sp_String      "maxOfInfo_sp"
#define  MinOfInfo_sp_String      "minOfInfo_sp"

public:
	statInfo(const vectorstr & sumOfInfo, const vectorstr & meanOfInfo,
		const vectorstr & varOfInfo, const vectorstr & maxOfInfo,
		const vectorstr & minOfInfo,
		const subPopList & subPops, const stringList & vars);

	~statInfo()
	{
	}


	bool apply(population & pop);

private:
	vectorstr m_sumOfInfo;
	vectorstr m_meanOfInfo;
	vectorstr m_varOfInfo;
	vectorstr m_maxOfInfo;
	vectorstr m_minOfInfo;

	subPopList m_subPops;
	stringList m_vars;
};


/// CPPONLY
class statLD
{
private:
#define   LD_String            "LD"
#define   LD_prime_String      "LD_prime"
#define   R2_String            "R2"
#define   ChiSq_String         "LD_ChiSq"
#define   ChiSq_p_String       "LD_ChiSq_p"
#define   CramerV_String       "CramerV"

#define   LD_sp_String         "LD_sp"
#define   LD_prime_sp_String   "LD_prime_sp"
#define   R2_sp_String         "R2_sp"
#define   ChiSq_sp_String      "LD_ChiSq_sp"
#define   ChiSq_p_sp_String    "LD_ChiSq_p_sp"
#define   CramerV_sp_String    "CramerV_sp"

public:
	// In the previous versions (< 0.9.6), statLD relies statAlleleFreq
	// and statHaploFreq to obtain allele and haplotype frequencies. This
	// complicates the structure of these two statistics (because they need
	// to decide how to save and provide statistics for statLD), and make
	// statLD less efficient. The only advanced was that if both allele
	// or haplotype frequencies and LD are requested, the frequencies will be
	// calculated only once. However, this appear to be a rare case that does
	// not worth special optimization. The newer version calculates allele and
	// haplotype frequencies locally and in a more readable way.
	statLD(const intMatrix & LD, const subPopList & subPops,
		const stringList & vars);

	// calculate, right now,  do not tempt to save values
	bool apply(population & pop);

private:
	typedef map<UINT, UINT> ALLELECNT;
	typedef vector<ALLELECNT> ALLELECNTLIST;
	typedef map<std::pair<UINT, UINT>, UINT> HAPLOCNT;
	typedef vector<HAPLOCNT> HAPLOCNTLIST;

	// calculate single allele LD values
	void calculateLD(const vectoru & lociMap,
		const ALLELECNTLIST & alleleCnt, const HAPLOCNTLIST & haploCnt,
		vectorf & LD, vectorf & D_prime, vectorf & R2, vectorf & ChiSq, vectorf & ChiSq_p,
		vectorf & CramerV);

	void outputVar(population & pop, const string & name, const vectorf & value);

private:
	/// LD
	intMatrix m_LD;

	subPopList m_subPops;
	stringList m_vars;
};

/// CPPONLY
class statAssociation
{
private:
#define Allele_ChiSq_String      "Allele_ChiSq"
#define Allele_ChiSq_p_String    "Allele_ChiSq_p"
#define Geno_ChiSq_String        "Geno_ChiSq"
#define Geno_ChiSq_p_String      "Geno_ChiSq_p"
#define Armitage_p_String        "Armitage_p"

#define Allele_ChiSq_sp_String   "Allele_ChiSq_sp"
#define Allele_ChiSq_p_sp_String "Allele_ChiSq_p_sp"
#define Geno_ChiSq_sp_String     "Geno_ChiSq_sp"
#define Geno_ChiSq_p_sp_String   "Geno_ChiSq_p_sp"
#define Armitage_p_sp_String     "Armitage_p_sp"

private:
	typedef map<UINT, ULONG>  ALLELECNT;
	typedef vector<ALLELECNT> ALLELECNTLIST;
	typedef map<std::pair<UINT, UINT>, ULONG>  GENOCNT;
	typedef vector<GENOCNT> GENOCNTLIST;

public:
	statAssociation(const vectorlu & loci,
		const subPopList & subPops, const stringList & vars);

	// calculate, right now,  do not tempt to save values
	bool apply(population & pop);

private:
	void alleleChiSqTest(const ALLELECNT & caseCnt,
		const ALLELECNT & controlCnt, double & chisq,
		double & chisq_p);

	void genoChiSqTest(const GENOCNT & caseCnt,
		const GENOCNT & controlCnt, double & chisq,
		double & chisq_p);

	double armitageTest(const GENOCNT & caseCnt,
		const GENOCNT & controlCnt);

private:
	/// Association
	vectorlu m_loci;

	subPopList m_subPops;
	stringList m_vars;
};

/// CPPONLY
class statNeutrality
{
private:
#define Neutra_Pi_String      "Pi"
#define Neutra_Pi_sp_String   "Pi_sp"

public:
	statNeutrality(const vectorlu & loci, const subPopList & subPops,
		const stringList & vars);

	// calculate, right now,  do not tempt to save values
	bool apply(population & pop);

private:
	typedef vector<vectora> HAPLOLIST;
	double calcPi(HAPLOLIST::const_iterator begin, HAPLOLIST::const_iterator end);

private:
	/// Neutrality
	vectorlu m_loci;

	subPopList m_subPops;
	stringList m_vars;
};

/// CPPONLY currently there is no need to retrieve calculated value
class statFst
{

private:
#define  Fst_String     "Fst"
#define  Fis_String     "Fis"
#define  Fit_String     "Fit"
#define  AvgFst_String  "AvgFst"
#define  AvgFis_String  "AvgFis"
#define  AvgFit_String  "AvgFit"

public:
	statFst(const vectorlu & Fst, const subPopList & subPops, const stringList & vars);

	bool apply(population & pop);

private:
	/// Fst
	vectorlu m_loci;

	subPopList m_subPops;
	stringList m_vars;
};


/// CPPONLY
class statHWE
{
private:
#define  HWE_String     "HWE"
#define  HWE_sp_String  "HWE_sp"

private:
	typedef map<std::pair<UINT, UINT>, ULONG>  GENOCNT;
	typedef vector<GENOCNT> GENOCNTLIST;

public:
	statHWE(const vectorlu & loci, const subPopList & subPops,
		const stringList & vars);

	bool apply(population & pop);

private:
	vectorlu mapToCount(const GENOCNT & cnt);

private:
	vectorlu m_loci;
	subPopList m_subPops;
	stringList m_vars;
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
	 *  \li \c subPopSize (default): A list of (virtual) subpopulation sizes.
	 *      This variable is easier to use than accessing popSize from each
	 *      (virtual) subpopulation.
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
	 *
	 *  <b>haploFreq</b>: This parameter accepts one or more lists of loci (by
	 *  index) at which number and frequency of haplotypes are outputted as
	 *  dictionaries. <tt>[(1,2)]</tt> can be abbreviated to <tt>(1,2)</tt>.
	 *  For example, using parameter <tt>haploFreq=(1,2,4)</tt>, all haplotypes
	 *  at these three loci are counted. Results are saved in variables such as
	 *  <tt>haploFreq[(1,2,4)][(1,1,0)]</tt> (frequency of haplotype (1,1,0)).
	 *  This statistic works for all population types. Number of haplotypes for
	 *  each individual equals to his/her ploidy number. Haplodiploid
	 *  populations are supported in the sense that the second homologous copy
	 *  of the haplotype is not counted for male individuals. This statistic
	 *  outputs the following variables:
	 *  \li \c haploFreq (default): A dictionary (with tuples of loci indexes
	 *       as keys) of dictionaries of haplotype frequencies. For example,
	 *       <tt>haploFreq[(0, 1)][(1,1)]</tt> records the frequency of
	 *       haplotype <tt>(1,1)</tt> at loci <tt>(0, 1)</tt> in all or
	 *       specified (virtual) subpopulations.
	 *  \li \c haploNum (default): A dictionary of dictionaries of haplotype
	 *       counts in all or specified (virtual) subpopulations.
	 *  \li \c haploFreq_sp: Halptype frequencies in each (virtual)
	 *       subpopulation.
	 *  \li \c haploNum_sp: Halptype count in each (virtual) subpopulation.
	 *
	 *  <b>sumOfinfo</b>, <b>meanOfInfo</b>, <b>varOfInfo</b>, <b>maxOfInfo</b>
	 *  and <b>minOfInfo</b>: Each of these five parameters accepts a list of
	 *  information fields. For each information field, the sum, mean, variance,
	 *  maximum or minimal (depending on the specified parameter(s)) of this
	 *  information field at iddividuals in all or specified (virtual)
	 *  subpopulations will be calculated. The results will be put into the
	 *  following population variables:
	 *  \li \c sumOfInfo (default for \e sumOfInfo): A dictionary of the sum of
	 *       specified information fields of individuals in all or specified
	 *       (virtual) subpopulations. This dictionary is indexed by names of
	 *       information fields.
	 *  \li \c meanOfInfo (default for \e meanOfInfo): A dictionary of the mean
	 *       of information fields of all individuals.
	 *  \li \c varOfInfo (default for \e varOfInfo): A dictionary of the sample
	 *       variance of information fields of all individuals.
	 *  \li \c maxOfInfo (default for \e maxOfInfo): A dictionary of the
	 *       maximum value of information fields of all individuals.
	 *  \li \c minOfInfo (default for \e minOfInfo): A dictionary of the
	 *       minimal value of information fields of all individuals.
	 *  \li \c sumOfInfo_sp: A dictionary of the sum of information fields of
	 *       individuals in each subpopulation.
	 *  \li \c meanOfInfo_sp: A dictionary of the mean of information fields of
	 *       individuals in each subpopulation.
	 *  \li \c varOfInfo_sp: A dictionary of the sample variance of information
	 *       fields of individuals in each subpopulation.
	 *  \li \c maxOfInfo_sp: A dictionary of the maximum value of information
	 *       fields of individuals in each subpopulation.
	 *  \li \c minOfInfo_sp: A dictionary of the minimal value of information
	 *       fields of individuals in each subpopulation.
	 *
	 *  <b>LD</b>: Parameter \c LD accepts one or a list of loci pairs (e.g.
	 *  <tt>LD=[[0,1], [2,3]]</tt>) with optional primary alleles at both loci
	 *  (e.g. <tt>LD=[0,1,0,0]</tt>). For each pair of loci, this operator
	 *  calculates linkage disequilibrium and optional association statistics
	 *  between two loci. When primary alleles are specified, signed linkage
	 *  disequilibrium values are calculated with non-primary alleles are
	 *  combined. Otherwise, absolute values of diallelic measures are combined
	 *  to yield positive measure of LD. Association measures are calculated
	 *  from a \c m by \c n contigency of haplotype counts (<tt>m=n=2</tt>
	 *  if primary alleles are specified). Please refer to the simuPOP user's
	 *  guide for detailed information. This statistic sets the following
	 *  variables:
	 *  \li \c LD (default) Basic LD measure for haplotypes in all or specified
	 *       (virtual) subpopulations. Signed if primary alleles are specified.
	 *  \li \c LD_prime (default) Lewontin's D' measure for haplotypes in all
	 *       or specified (virtual) subpopulations. Signed if primary alleles
	 *       are specified.
	 *  \li \c R2 (default) Correlation LD measure for haplotypes in all or
	 *       specified (virtual) subpopulations.
	 *  \li \c LD_ChiSq ChiSq statistics for a contigency table with
	 *       frequencies of haplotypes in all or specified (virtual)
	 *       subpopulations.
	 *  \li \c LD_ChiSq_p Single side p-value for the ChiSq statistic. Degrees
	 *       of freedom is determined by number of alleles at both loci and the
	 *       specification of primary alleles.
	 *  \li \c CramerV Normalized ChiSq statistics.
	 *  \li \c LD_sp Basic LD measure for haplotypes in each (virtual)
	 *       subpopulation.
	 *  \li \c LD_prime_sp Lewontin's D' measure for haplotypes in each
	 *       (virtual) subpopulation.
	 *  \li \c R2_sp R2 measure for haplotypes in each (virtual) subpopulation.
	 *  \li \c LD_ChiSq_sp ChiSq statistics for each (virtual) subpopulation.
	 *  \li \c LD_ChiSq_p_sp p value for the ChiSq statistics for each
	 *       (virtual) subpopulation.
	 *  \li \c CramerV_sp Cramer V statistics for each (virtual) subpopulation.
	 *
	 *  <b>association</b>: Parameter \c association accepts a list of loci.
	 *  At each locus, one or more statistical tests will be performed to test
	 *  association between this locus and individual affection status.
	 *  Currently, simuPOP provides the following tests:
	 *  \li An allele-based Chi-square test using alleles counts. This test
	 *       can be applied to loci with more than two alleles, and to haploid
	 *       populations.
	 *  \li A genotype-based Chi-square test using genotype counts. This test
	 *       can be applied to loci with more than two alleles (more than 3
	 *       genotypes) in diploid populations. \c aA and \c Aa are considered
	 *       to be the same genotype.
	 *  \li A genotype-based Cochran-Armitage trend test. This test can only
	 *       be applied to diallelic loci in diploid populations. A codominant
	 *       model is assumed.
	 *
	 *  This statistic sets the following variables:
	 *  \li \c Allele_ChiSq A dictionary of allele-based Chi-Square statistics
	 *       for each locus, using cases and controls in all or specified
	 *       (virtual) subpopulations.
	 *  \li \c Allele_ChiSq_p (default) A dictionary of \e p-values of the
	 *       corresponding Chi-square statistics.
	 *  \li \c Geno_ChiSq A dictionary of genotype-based Chi-Square statistics
	 *       for each locus, using cases and controls in all or specified
	 *       (virtual) subpopulations.
	 *  \li \c Geno_ChiSq_p A dictionary of \e p-values of the corresponding
	 *       genotype-based Chi-square test.
	 *  \li \c Armitage_p A dictionary of \e p-values of the Cochran-Armitage
	 *       tests, using cases and controls in all or specified (virtual)
	 *       subpopulations.
	 *  \li \c Allele_ChiSq_sp A dictionary of allele-based Chi-Square
	 *       statistics for each locus, using cases and controls from each
	 *       subpopulation.
	 *  \li \c Allele_ChiSq_p_sp A dictionary of p-values of allele-based
	 *       Chi-square tests, using cases and controls from each
	 *       (virtual) subpopulation.
	 *  \li \c Geno_ChiSq_sp A dictionary of genotype-based Chi-Square tests for
	 *       each locus, using cases and controls from each subpopulation.
	 *  \li \c Geno_ChiSq_p_sp A dictionary of p-values of genotype-based
	 *       Chi-Square tests, using cases and controls from each
	 *       subpopulation.
	 *  \li \c Armitage_p_sp A dictionary of \e p-values of the Cochran-
	 *       Armitage tests, using cases and controls from each subpopulation.
	 *
	 *  <b>neutrality</b>: This parameter performs neutrality tests (detection
	 *  of natural selection) on specified loci. It currently only outputs
	 *  \e Pi, which is the average number of pairwise difference between loci.
	 *  This statistic outputs the following variables:
	 *  \li \c Pi Mean pairwise difference between all sequences from all or
	 *       specified (virtual) subpopulations.
	 *  \li \c Pi_sp Mean paiewise difference between all sequences in each
	 *       (virtual) subpopulation.
	 *
	 *  <b>Fst</b>: Parameter \c Fst accepts a list of loci at which level of
	 *  population structure is measured by statistic \e Fst using an algorithm
	 *  developed in Cockerham & Weir, 1984. \e Fst is by default calculated
	 *  for all subpopulation but a subset of subpopulations could be specified
	 *  using parameter \e subPops. Virtual subpopulations are supported which
	 *  makes it possible to estimate population structure from groups of
	 *  individuals from the same subpopulation. This statistic calulates
	 *  \e Fis, \e Fit and \e Fst for each locus and for all loci and set
	 *  variables:
	 *  \li \c AvgFst (default) The \e Fst statistic estimated for all
	 *       specified loci.
	 *  \li \c AvgFis The \e Fis statistic estimated for all specified loci.
	 *  \li \c AvgFit The \e Fit statistic estimated for all specified loci.
	 *  \li \c Fst A dictionary of locus level \e Fst values.
	 *  \li \c Fis A dictionary of locus level \e Fis values.
	 *  \li \c Fit A dictionary of locus level \e Fit values.
	 *
	 *  <b>HWE</b>: Parameter \c HWE accepts a list of loci at which exact
	 *  two-side tests for Hardy-Weinberg equilibrium will be performed. This
	 *  statistic is only available for diallelic loci in diploid populations.
	 *  It outputs the following variables:
	 *  \li \c HWE (default) A dictionary of p-values of HWE tests using
	 *       genotypes in all or specified (virtual) subpopulations.
	 *  \li \c HWE_sp A dictionary of p-values of HWS tests using genotypes
	 *       in each (virtual) subpopulation.
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
		const stringList & sumOfInfo = stringList(),
		const stringList & meanOfInfo = stringList(),
		const stringList & varOfInfo = stringList(),
		const stringList & maxOfInfo = stringList(),
		const stringList & minOfInfo = stringList(),
		//
		const intMatrix & LD = intMatrix(),
		//
		const uintList & association = uintList(),
		//
		const uintList & neutrality = uintList(),
		//
		const uintList & Fst = uintList(),
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
	statInfo m_info;
	statLD m_LD;
	statAssociation m_association;
	statNeutrality m_neutrality;
	statFst m_Fst;
	statHWE m_HWE;
};

}
#endif
