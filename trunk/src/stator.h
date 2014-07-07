/**
 *  $File: stator.h $
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

#ifndef _STATOR_H
#define _STATOR_H
/**
   \file
   \brief head file of class BaseOperator:public BaseOperator
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


/** A \c PyEval operator evaluates a Python expression in a population's local
 *  namespace when it is applied to this population. The result is written to
 *  an output specified by parameter \e output.
 */
class PyEval : public BaseOperator
{
public:
	/** Create a \c PyEval operator that evaluates a Python expression \e expr
	 *  in a population's local namespaces when it is applied to this
	 *  population. This namespace can either be the population's local
	 *  namespace (<tt>pop.vars()</tt>), or namespaces <tt>subPop[sp]</tt> for
	 *  (virtual) subpop (<tt>pop.vars(subpop)</tt>) in specified \e subPops.
	 *  If Python statements \e stmts is given (a single or multi-line string),
	 *  the statement will be executed before \e expr. If \e exposePop is set
	 *  to an non-empty string, the current population will be exposed in its
	 *  own local namespace as a variable with this name. This allows the
	 *  execution of expressions such as <tt>'pop.individual(0).allele(0)'</tt>.
	 *  The result of \e expr will be sent to an output stream specified by
	 *  parameter \c output. The exposed population variable will be removed
	 *  after \e expr is evaluated. Please refer to class \c BaseOperator for
	 *  other parameters.
	 *
	 *  \note Although the statements and expressions are evaluated in a
	 *  population's local namespace, they have access to a global
	 *  namespace which is the module global namespace. It is therefore
	 *  possible to refer to any module variable in these expressions. Such
	 *  mixed use of local and global variables is, however, strongly
	 *  discouraged.
	 */
	PyEval(const string & expr = string(), const string & stmts = string(),
		const string & exposePop = string(), const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = Py_False,
		const stringList & infoFields = vectorstr())
		: BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_expr(expr, stmts), m_exposePop(exposePop)
	{
	}


	~PyEval()
	{
	}


	/// HIDDEN Deep copy of a \c PyEval operator
	virtual BaseOperator * clone() const
	{
		return new PyEval(*this);
	}


	/** HIDDEN Evaluate the expression and optional statements in the local namespace
	 *  of population \e pop and return its result as a string.
	 */
	string evaluate(Population & pop) const;

	/// HIDDEN Apply the \c PyEval operator to population \e pop.
	virtual bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const;

protected:
	/// expression to evaluate
	const Expression m_expr;

	/// if expose pop
	const string m_exposePop;
};


/** This operator executes given Python statements in a population's local
 *  namespace when it is applied to this population.
 */
class PyExec : public PyEval
{
public:
	/** Create a \c PyExec operator that executes statements \e stmts in a
	 *  population's local namespace when it is applied to this population.
	 *  This namespace can either be the population's local namespace
	 *  (<tt>pop.vars()</tt>), or namespaces <tt>subPop[sp]</tt> for each
	 *  (virtual) subpop (<tt>pop.vars(subpop)</tt>) in specified \e subPops.
	 *  If \e exposePop is given, current population will be exposed in
	 *  its local namespace as a variable named by \e exposePop. Although
	 *  multiple statements can be executed, it is recommended that you use
	 *  this operator to execute short statements and use \c PyOperator for
	 *  more complex once. Note that exposed population variables will be
	 *  removed after the statements are executed.
	 */
	PyExec(const string & stmts = string(), const string & exposePop = string(),
		const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = Py_False,
		const stringList & infoFields = vectorstr())
		: PyEval("", stmts, exposePop, "", begin, end, step, at, reps, subPops, infoFields)
	{
		(void)output;  // avoid warning about unused parameter
	}


	~PyExec()
	{
	}


	/// HIDDEN Deep copy of a \c PyExec operator
	virtual BaseOperator * clone() const
	{
		return new PyExec(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;

};


/** Unlike operator \c PyEval and \c PyExec that work at the population level,
 *  in a population's local namespace, operator \c InfoEval works at the
 *  individual level, working with individual information fields. When this
 *  operator is applied to a population, information fields of eligible
 *  individuals are put into the local namespace of the population. A
 *  Python expression is then evaluated for each individual. The result is
 *  written to an output.
 */
class InfoEval : public BaseOperator
{
public:
	/** Create an operator that evaluate a Python expression \e expr using
	 *  individual information fields and population variables as variables.
	 *  If \e exposeInd is not empty, the individual itself will be exposed in
	 *  the population's local namespace as a variable with name specified by
	 *  \e exposeInd.
	 *
	 *  A Python expression (\e expr) is evaluated for each individual. The
	 *  results are converted to strings and are written to an output specified
	 *  by parameter \e output. Optionally, a statement (or several statements
	 *  separated by newline) can be executed before \e expr is evaluated. The
	 *  evaluation of this statement may change the value of information
	 *  fields.
	 *
	 *  Parameter \e usePopVars is obsolete because population variables are
	 *  always usable in such expressions.
	 */
	InfoEval(const string & expr = string(), const string & stmts = string(), bool usePopVars = false,
		const string & exposeInd = string(),
		const stringFunc & output = ">", int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr())
		: BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_expr(expr, stmts), m_exposeInd(exposeInd), m_lastValues()
	{
		(void)usePopVars;  // this parameter is obsolete, use (void) to avoid a warning message
		DBG_WARNIF(debug(DBG_COMPATIBILITY) && usePopVars, "WARNING: parameter usePopVars is obsolete.");
	}


	~InfoEval()
	{
	}


	/// HIDDEN Deep copy of a \c InfoEval operator
	virtual BaseOperator * clone() const
	{
		return new InfoEval(*this);
	}


	// check all alleles in vector allele if they are fixed.
	/// HIDDEN apply the \c InfoEval operator
	bool apply(Population & pop) const;

	bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// HIDDEN
	string describe(bool format = true) const;

protected:
	string evalInfo(Individual * ind, PyObject * dict) const;

	void clearVars(Population & pop) const;

	/// expression to evaluate
	const Expression m_expr;

	const string m_exposeInd;
	/// cache last values to speed up evaluation, more specifically,
	/// if the next individual holds the same value at an information field
	/// existing variable will not be set again.
	mutable vectorf m_lastValues;
};

/** Operator \c InfoExec is similar to \c InfoEval in that it works at the
 *  individual level, using individual information fields as variables.
 *  This is usually used to change the value of information fields. For
 *  example, <tt>"b=a*2"</tt> will set the value of information field \c b
 *  to <tt>a*a</tt> for all individuals.
 */
class InfoExec : public InfoEval
{
public:
	/** Create an operator that executes Python statements \e stmts using
	 *  individual information fields and population variables as variables.
	 *  If \e exposeInd is not empty, the individual itself will be exposed in
	 *  the population's local namespace as a variable with name specified by
	 *  \e exposeInd.
	 *
	 *  One or more python statements (\e stmts) are executed for each
	 *  individual. Information fields of these individuals are then updated
	 *  from the corresponding variables. For example, <tt>a=1</tt> will set
	 *  information field \e a of all individuals to \c 1, <tt>a=b</tt> will
	 *  set information field \e a of all individuals to information field
	 *  \c b or a population variable \c b if \c b is not an information field
	 *  but a population variable, and <tt>a=ind.sex()</tt> will set
	 *  information field \e a of all individuals to its sex (needs
	 *  <tt>exposeInd='ind'</tt>.
	 *
	 *  Parameter \e usePopVars is obsolete because population variables will
	 *  always be usable.
	 */
	InfoExec(const string & stmts = string(), bool usePopVars = false,  const string & exposeInd = string(),
		const stringFunc & output = "", int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr())
		: InfoEval(string(), stmts, usePopVars, exposeInd, output, begin, end, step, at, reps, subPops, infoFields),
		m_simpleStmt(stmts, exposeInd)
	{
	}


	~InfoExec()
	{
	}


	/// HIDDEN Deep copy of a \c InfoExec operator
	virtual BaseOperator * clone() const
	{
		return new InfoExec(*this);
	}


	// check all alleles in vector allele if they are fixed.
	/// HIDDEN apply the \c InfoExec operator
	bool apply(Population & pop) const;

	bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// HIDDEN
	string describe(bool format = true) const;

private:
	const simpleStmt m_simpleStmt;
};


/// CPPONLY post population sizes etc.
class statPopSize
{
private:
#define  popSize_String        "popSize"
#define  subPopSize_String     "subPopSize"
#define  popSize_sp_String     "popSize_sp"

public:
	statPopSize(bool popSize, const subPopList & subPops,
		const stringList & vars, const string & suffix);

	bool apply(Population & pop) const;

	string describe(bool format = true) const;

private:
	bool m_isActive;

	subPopList m_subPops;

	stringList m_vars;
	string m_suffix;
};

/// CPPONLY
class statNumOfMales
{
private:
#define  numOfMales_String         "numOfMales"
#define  propOfMales_String        "propOfMales"
#define  numOfFemales_String       "numOfFemales"
#define  propOfFemales_String      "propOfFemales"
#define  numOfMales_sp_String      "numOfMales_sp"
#define  propOfMales_sp_String     "propOfMales_sp"
#define  numOfFemales_sp_String    "numOfFemales_sp"
#define  propOfFemales_sp_String   "propOfFemales_sp"

public:
	statNumOfMales(bool numOfMales, const subPopList & subPops,
		const stringList & vars, const string & suffix);

	bool apply(Population & pop) const;

	string describe(bool format = true) const;

private:
	/// whether or not to apply number of male/female
	bool m_isActive;
	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
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
	statNumOfAffected(bool numOfAffected, const subPopList & subPops,
		const stringList & vars, const string & suffix);

	string describe(bool format = true) const;

	bool apply(Population & pop) const;

private:
	/// whether or not to apply number of affected
	bool m_isActive;
	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
};


/// CPPONLY
class statNumOfSegSites
{
private:
#define  numOfSegSites_String        "numOfSegSites"
#define  numOfFixedSites_String      "numOfFixedSites"
#define  numOfSegSites_sp_String     "numOfSegSites_sp"
#define  numOfFixedSites_sp_String   "numOfFixedSites_sp"
#define  segSites_String             "segSites"
#define  fixedSites_String           "fixedSites"
#define  segSites_sp_String          "segSites_sp"
#define  fixedSites_sp_String        "fixedSites_sp"

public:
	statNumOfSegSites(const lociList & loci, const subPopList & subPops,
		const stringList & vars, const string & suffix);

	string describe(bool format = true) const;

	bool apply(Population & pop) const;

private:
	/// whether or not to apply number of affected
	lociList m_loci;
	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
};


/// CPPONLY
class statNumOfMutants

{
private:
#define  numOfMutants_String        "numOfMutants"
#define  numOfMutants_sp_String     "numOfMutants_sp"

public:
	statNumOfMutants(const lociList & loci, const subPopList & subPops,
		const stringList & vars, const string & suffix);

	string describe(bool format = true) const;

	bool apply(Population & pop) const;

private:
	/// whether or not to apply number of affected
	lociList m_loci;
	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
};


/// CPPONLY
class statAlleleFreq
{
private:
#define  AlleleNum_String        "alleleNum"
#define  AlleleFreq_String       "alleleFreq"
#define  AlleleNum_sp_String     "alleleNum_sp"
#define  AlleleFreq_sp_String    "alleleFreq_sp"

private:
	typedef uintDict ALLELECNT;
	typedef vector<ALLELECNT> ALLELECNTLIST;

public:
	statAlleleFreq(const lociList & loci, const subPopList & subPops,
		const stringList & vars, const string & suffix);

	string describe(bool format = true) const;

	/// destructor, nested vectors have to be cleared manually
	~statAlleleFreq()
	{
	}


	bool apply(Population & pop) const;

private:
	/// which alleles?
	lociList m_loci;

	subPopList m_subPops;

	stringList m_vars;
	string m_suffix;
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
	statHeteroFreq(const lociList & heteroFreq, const lociList & homoFreq,
		const subPopList & subPops, const stringList & vars, const string & suffix);

	string describe(bool format = true) const;

	bool apply(Population & pop) const;

private:
	/// heteroFreq
	lociList m_loci;

	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
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
	statGenoFreq(const lociList & genoFreq,  const subPopList & subPops,
		const stringList & vars, const string & suffix);

	string describe(bool format = true) const;

	bool apply(Population & pop) const;

private:
	/// which genotypes
	lociList m_loci;

	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
};


/// CPPONLY
class statHaploFreq
{
private:
#define HaplotypeNum_String         "haploNum"
#define HaplotypeFreq_String        "haploFreq"
#define HaplotypeNum_sp_String      "haploNum_sp"
#define HaplotypeFreq_sp_String     "haploFreq_sp"

public:
	statHaploFreq(const intMatrix & haploFreq, const subPopList & subPops,
		const stringList & vars, const string & suffix);

	string describe(bool format = true) const;

	~statHaploFreq()
	{
	}


	bool apply(Population & pop) const;

private:
	// key string (in the format of a tuple)
	string dictKey(const vectori & loci) const;

private:
	/// haplotype at which loci
	matrixi m_loci;

	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
};


/// CPPONLY
class statHaploHomoFreq
{
private:
#define HaploHeteroNum_String        "haploHeteroNum"
#define HaploHeteroFreq_String       "haploHeteroFreq"
#define HaploHomoNum_String          "haploHomoNum"
#define HaploHomoFreq_String         "haploHomoFreq"
#define HaploHeteroNum_sp_String     "haploHeteroNum_sp"
#define HaploHeteroFreq_sp_String    "haploHeteroFreq_sp"
#define HaploHomoNum_sp_String       "haploHomoNum_sp"
#define HaploHomoFreq_sp_String      "haploHomoFreq_sp"

public:
	statHaploHomoFreq(const intMatrix & haploHeteroFreq, const intMatrix & haploHomoFreq,
		const subPopList & subPops, const stringList & vars, const string & suffix);

	~statHaploHomoFreq()
	{
	}


	string describe(bool format = true) const;

	bool apply(Population & pop) const;

private:
	/// heteroFreq
	matrixi m_loci;

	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
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
		const vectorstr & minOfInfo, const subPopList & subPops,
		const stringList & vars, const string & suffix);


	string describe(bool format = true) const;

	~statInfo()
	{
	}


	bool apply(Population & pop) const;

private:
	vectorstr m_sumOfInfo;
	vectorstr m_meanOfInfo;
	vectorstr m_varOfInfo;
	vectorstr m_maxOfInfo;
	vectorstr m_minOfInfo;

	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
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
		const stringList & vars, const string & suffix);

	string describe(bool format = true) const;

	// calculate, right now,  do not tempt to save values
	bool apply(Population & pop) const;

private:
	typedef map<size_t, size_t> ALLELECNT;
	typedef vector<ALLELECNT> ALLELECNTLIST;
	typedef map<pairu, size_t> HAPLOCNT;
	typedef vector<HAPLOCNT> HAPLOCNTLIST;

	// calculate single allele LD values
	void calculateLD(const vectoru & lociMap,
		const ALLELECNTLIST & alleleCnt, const HAPLOCNTLIST & haploCnt,
		vectorf & LD, vectorf & D_prime, vectorf & R2, vectorf & ChiSq, vectorf & ChiSq_p,
		vectorf & CramerV) const;

	void outputVar(Population & pop, const string & name, const vectorf & value) const;

private:
	/// LD
	matrixi m_LD;

	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
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
	typedef map<Allele, size_t>  ALLELECNT;
	typedef vector<ALLELECNT> ALLELECNTLIST;
	typedef map<std::pair<Allele, Allele>, size_t>  GENOCNT;
	typedef vector<GENOCNT> GENOCNTLIST;

public:
	statAssociation(const lociList & loci, const subPopList & subPops,
		const stringList & vars, const string & suffix);


	string describe(bool format = true) const;

	// calculate, right now,  do not tempt to save values
	bool apply(Population & pop) const;

private:
	void alleleChiSqTest(const ALLELECNT & caseCnt,
		const ALLELECNT & controlCnt, double & chisq,
		double & chisq_p) const;

	void genoChiSqTest(const GENOCNT & caseCnt,
		const GENOCNT & controlCnt, double & chisq,
		double & chisq_p) const;

	double armitageTest(const GENOCNT & caseCnt,
		const GENOCNT & controlCnt) const;

private:
	/// Association
	lociList m_loci;

	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
};

/// CPPONLY
class statNeutrality
{
private:
#define Neutra_Pi_String      "Pi"
#define Neutra_Pi_sp_String   "Pi_sp"

public:
	statNeutrality(const lociList & loci, const subPopList & subPops,
		const stringList & vars, const string & suffix);

	string describe(bool format = true) const;

	// calculate, right now,  do not tempt to save values
	bool apply(Population & pop) const;

private:
	typedef vector<vectora> HAPLOLIST;
	double calcPi(HAPLOLIST::const_iterator begin, HAPLOLIST::const_iterator end) const;

private:
	/// Neutrality
	lociList m_loci;

	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
};

/// CPPONLY currently there is no need to retrieve calculated value
class statStructure
{

private:
#define  fst_String     "f_st"
#define  fis_String     "f_is"
#define  fit_String     "f_it"
#define  Fst_String     "F_st"
#define  Fis_String     "F_is"
#define  Fit_String     "F_it"

#define  Gst_String     "G_st"
#define  gst_String     "g_st"

public:
	statStructure(const lociList & Fst, const subPopList & subPops,
		const stringList & vars, const string & suffix);

	string describe(bool format = true) const;

	bool apply(Population & pop) const;

private:
	typedef map<size_t, float> FREQ;
	typedef map<size_t, FREQ> LOCIFREQ;
	typedef vector<LOCIFREQ> LOCIFREQLIST;
	typedef map<Allele, bool> ALLELES;
	typedef vector<ALLELES> ALLELELIST;

	void calcGst_Nei73(const vectoru & loci, const vectoru & n_i, LOCIFREQLIST & alleleFreq,
		const ALLELELIST & alleles, double & Gst, uintDict & gst) const;

	void calcFst_WC84(const vectoru & loci, const vectoru & n_i, LOCIFREQLIST & alleleFreq, LOCIFREQLIST & heteroFreq,
		const ALLELELIST & alleles, double & Fst, double & Fis, double & Fit,
		uintDict & fst, uintDict & fis, uintDict & fit) const;

private:
	/// Fst
	lociList m_loci;

	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
};


/// CPPONLY
class statHWE
{
private:
#define  HWE_String     "HWE"
#define  HWE_sp_String  "HWE_sp"

private:
	typedef map<pairu, size_t>  GENOCNT;
	typedef vector<GENOCNT> GENOCNTLIST;

public:
	statHWE(const lociList & loci, const subPopList & subPops,
		const stringList & vars, const string & suffix);


	string describe(bool format = true) const;

	bool apply(Population & pop) const;

private:
	vectoru mapToCount(const GENOCNT & cnt) const;

private:
	lociList m_loci;
	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
};


/// CPPONLY
class statInbreeding
{
private:
#define  IBD_freq_String     "IBD_freq"
#define  IBD_freq_sp_String  "IBD_freq_sp"
#define  IBS_freq_String     "IBS_freq"
#define  IBS_freq_sp_String  "IBS_freq_sp"

public:
	statInbreeding(const lociList & loci, const subPopList & subPops,
		const stringList & vars, const string & suffix);

	string describe(bool format = true) const;

	bool apply(Population & pop) const;

private:
	lociList m_loci;
	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
};


/// CPPONLY
class statEffectiveSize
{
private:
#define  Ne_demo_base_String         "Ne_demo_base"
#define  Ne_demo_base_sp_String      "Ne_demo_base_sp"
#define  Ne_demo_String              "Ne_demo"
#define  Ne_demo_sp_String           "Ne_demo_sp"

#define  Ne_temporal_base_String     "Ne_temporal_base"
#define  Ne_temporal_base_sp_String  "Ne_temporal_base_sp"

	// deprecated, use Ne_waples89_P1 etc instead
#define  Ne_waples89_String       "Ne_waples89"
#define  Ne_waples89_sp_String    "Ne_waples89_sp"
#define  Ne_tempoFS_String        "Ne_tempoFS"
#define  Ne_tempoFS_sp_String     "Ne_tempoFS_sp"

#define  Ne_waples89_P1_String       "Ne_waples89_P1"
#define  Ne_waples89_P1_sp_String    "Ne_waples89_P1_sp"
#define  Ne_tempoFS_P1_String        "Ne_tempoFS_P1"
#define  Ne_tempoFS_P1_sp_String     "Ne_tempoFS_P1_sp"
#define  Ne_waples89_P2_String       "Ne_waples89_P2"
#define  Ne_waples89_P2_sp_String    "Ne_waples89_P2_sp"
#define  Ne_tempoFS_P2_String        "Ne_tempoFS_P2"
#define  Ne_tempoFS_P2_sp_String     "Ne_tempoFS_P2_sp"

#define  Ne_LD_String             "Ne_LD"
#define  Ne_LD_sp_String          "Ne_LD_sp"
#define  Ne_LD_mono_String        "Ne_LD_mono"
#define  Ne_LD_mono_sp_String     "Ne_LD_mono_sp"

private:
	typedef uintDict ALLELECNT;
	typedef vector<ALLELECNT> ALLELECNTLIST;


	// calculate moment based estimate of Ne based on Waples 89
	void Waples89(size_t N, size_t S0, size_t St, size_t t,
		const ALLELECNTLIST & P0, const ALLELECNTLIST & Pt,
		vectorf & res1, vectorf & res2) const;

	// calculate moment based estimate of Ne based on Jorde & Ryman's (2007)
	void TempoFS(size_t N, size_t S0, size_t St, size_t t,
		const ALLELECNTLIST & P0, const ALLELECNTLIST & Pt,
		vectorf & res1, vectorf & res2) const;

	// calculate LD based on genotype counts
	// genotype at two loci i,j, k, l
	// no phase is assumed so i < j and k <l is assumed
	typedef std::vector<size_t> GENOTYPE;
	typedef std::map<GENOTYPE, size_t> GENOTYPECNT;
	typedef uintDict HOMOCNT;

	// the first list for R2, the second list for weight, for
	// 0.01, 0.02 and 0.05 cutoff values
	typedef std::pair<vectorf, vectoru> R2WEIGHT;
	typedef std::vector<R2WEIGHT> LDLIST;

	R2WEIGHT Burrows(size_t N, const ALLELECNT & a1, const ALLELECNT & a2,
		const HOMOCNT & h1, const HOMOCNT & h2, const GENOTYPECNT & g) const;

	void LDNe(const LDLIST & ld, int cutoff, size_t S, vectorf & res, vectorf & res_mono) const;

public:
	statEffectiveSize(const lociList & loci, const subPopList & subPops,
		const stringList & vars, const string & suffix);


	string describe(bool format = true) const;

	bool apply(Population & pop) const;

	bool demographicEffectiveSize(Population & pop) const;

	bool temporalEffectiveSize(Population & pop) const;

	bool LDEffectiveSize(Population & pop) const;

private:
	lociList m_loci;
	subPopList m_subPops;
	stringList m_vars;
	string m_suffix;
};


/** Operator \c Stat calculates various statistics of the population being
 *  applied and sets variables in its local namespace. Other operators or
 *  functions can retrieve results from or evalulate expressions in this local
 *  namespace after \c Stat is applied.
 */
class Stat : public BaseOperator
{
public:
	/** Create a \c Stat operator that calculates specified statistics of a
	 *  population when it is applied to this population. This operator can
	 *  be applied to specified replicates (parameter \e rep) at specified
	 *  generations (parameter \e begin, \e end, \e step, and \e at). This
	 *  operator does not produce any output (ignore parameter \e output)
	 *  after statistics are calculated. Instead, it stores results in the
	 *  local namespace of the population being applied. Other operators can
	 *  retrieve these variables or evalulate expression directly in this
	 *  local namespace. Please refer to operator \c BaseOperator for a
	 *  detailed explanation of these common operator parameters.
	 *
	 *  \c Stat supports parameter \e subPops. It usually calculate the same
	 *  set of statistics for all subpopulations (<tt>subPops=subPopList()</tt>).
	 *  If a list of (virtual) subpopulations are specified, statistics for
	 *  only specified subpopulations will be calculated. However, different
	 *  statistics treat this parameter differently and it is very important
	 *  to check its reference before you use \e subPops for any statistics.
	 *
	 *  Calculated statistics are saved as variables in a population's local
	 *  namespace. These variables can be numbers, lists or dictionaries and
	 *  can be retrieved using functions <tt>Population.vars()</tt> or
	 *  <tt>Population.dvars()</tt>. A special default dictionary (\c defdict)
	 *  is used for dictionaries whose keys are determined dynamically.
	 *  Accessing elements of such a dictionary with an invalid key will yield
	 *  value 0 instead of a \c KeyError. If the same variables are calculated
	 *  for one or more (virtual) subpopulation, the variables are stored in
	 *  <tt>vars()['subPop'][sp]['var']</tt> where sp is a subpopulation ID
	 *  (\c sp) or a tuple of virtual subpopulation ID (<tt>(sp, vsp)</tt>).
	 *  <tt>Population.vars(sp)</tt> and <tt>Population.dvars(sp)</tt> provide
	 *  shortcuts to these variables.
	 *
	 *  Operator \e Stat outputs a number of most useful variables for each
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
	 *  allele count. An optional suffix (parameter \e suffix) can be used
	 *  to append a suffix to default parameter names. This parameter can be
	 *  used, for example, to calculate and store the same statistics for
	 *  different subpopulations (e.g. pairwise \c Fst).
	 *
	 *  Operator \c Stat supports the following statistics:
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
	 *  <b>numOfMales</b>: If \e numOfMales=True, number of male individuals in
	 *  all or specified (virtual) subpopulations will be set to the following
	 *  variables:
	 *  \li \c numOfMales (default): Total number of male individuals in all
	 *       or specified (virtual) subpopulations.
	 *  \li \c numOfFemales (default): Total number of female individuals in all
	 *       or specified (virtual) subpopulations.
	 *  \li \c propOfMales: Proportion of male individuals.
	 *  \li \c propOfFemales: Proportion of female individuals.
	 *  \li \c numOfMales_sp: Number of male individuals in each (virtual)
	 *       subpopulation.
	 *  \li \c numOfFemales_sp: Number of female individuals in each (virtual)
	 *       subpopulation.
	 *  \li \c propOfMales_sp: Proportion of male individuals in each (virtual)
	 *       subpopulation.
	 *  \li \c propOfFemales_sp: Proportion of female individuals in each
	 *       (virtual) subpopulation.
	 *
	 *  <b>numOfAffected</b>: If \e numOfAffected=True, number of affected
	 *  individuals in all or specified (virtual) subpopulations will be set to the
	 *  following variables:
	 *  \li \c numOfAffected (default): Total number of affected individuals in
	 *       all or specified (virtual) subpopulations.
	 *  \li \c numOfUnaffected (default): Total number of unaffected individuals
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
	 *  <b>numOfSegSites</b>: Parameter \e numOfSegSites accepts a list of loci
	 *  (loci indexes, names, or \c ALL_AVAIL) and count the number of loci with
	 *  at least two different alleles (segregating sites) or loci with only one
	 *  non-zero allele (no zero allele, not segragating) for individuals in all
	 *  or specified (virtual) subpopulations. This parameter sets variables
	 *  \li \c numOfSegSites (default): Number of segregating sites in all or
	 *      specified (virtual) subpopulations.
	 *  \li \c numOfSegSites_sp: Number of segregating sites in each (virtual)
	 *      subpopulation.
	 *  \li \c numOfFixedSites: Number of sites with one non-zero allele in all
	 *      or specified (virtual) subpopulations.
	 *  \li \c numOfFixedSites_sp: Number of sites with one non-zero allele in
	 *      in each (virtual) subpopulations.
	 *  \li \c segSites: A list of segregating sites in all or specified
	 *      (virtual) subpopulations.
	 *  \li \c segSites_sp: A list of segregating sites in each (virtual)
	 *      subpopulation.
	 *  \li \c fixedSites: A list of sites with one non-zero allele in all
	 *      or specified (virtual) subpopulations.
	 *  \li \c fixedSites_sp: A list of sites with one non-zero allele in
	 *      in each (virtual) subpopulations.
	 *
	 *  <b>numOfMutants</b>: Parameter \e numOfMutants accepts a list of loci
	 *  (loci indexes, names, or \c ALL_AVAIL) and count the number of mutants
	 *  (non-zero alleles) for individuals in all or specified (virtual)
	 *  subpopulations. It sets variables
	 *  \li \c numOfMutants (default): Number of mutants in all or specified
	 *      (virtual) subpopulations.
	 *  \li \c numOfMutants_sp: Number of mutants in each (virtual)
	 *      subpopulations.
	 *
	 *  <b>alleleFreq</b>: This parameter accepts a list of loci (loci indexes,
	 *  names, or \c ALL_AVAIL), at which allele frequencies will be calculated.
	 *  This statistic outputs the following variables, all of which are
	 *  dictionary (with loci indexes as keys) of default dictionaries (with
	 *  alleles as keys). For example, <tt>alleleFreq[loc][a]</tt> returns 0
	 *  if allele \c a does not exist.
	 *  \li \c alleleFreq (default): <tt>alleleFreq[loc][a]</tt> is the
	 *       frequency of allele \c a at locus \loc for all or specified
	 *       (virtual) subpopulations.
	 *  \li \c alleleNum (default): <tt>alleleNum[loc][a]</tt> is the number of
	 *       allele \c a at locus \loc for all or specified (virtual)
	 *       subpopulations.
	 *  \li \c alleleFreq_sp: Allele frequency in each (virtual) subpopulation.
	 *  \li \c alleleNum_sp: Allele count in each (virtual) subpopulation.
	 *
	 *  <b>heteroFreq</b> and <b>homoFreq</b>: These parameters accept a list
	 *  of loci (by indexes or names), at which the number and frequency of
	 *  homozygotes and/or heterozygotes will be calculated. These statistics
	 *  are only available for diploid populations. The following variables
	 *  will be outputted:
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
	 *  <b>genoFreq</b>: This parameter accept a list of loci (by indexes or
	 *  names) at which number and frequency of all genotypes are outputed as a
	 *  dictionary (indexed by loci indexes) of default dictionaries (indexed
	 *  by tuples of possible indexes). This statistic is available for all
	 *  population types with genotype defined as ordered alleles at a locus.
	 *  The length of genotype equals the number of homologous copies of
	 *  chromosomes (ploidy) of a population. Genotypes for males or females
	 *  on sex chromosomes or in haplodiploid populations will have different
	 *  length. Because genotypes are ordered, <tt>(1, 0)</tt> and
	 *  <tt>(0, 1)</tt> (two possible genotypes in a diploid population) are
	 *  considered as different genotypes. This statistic outputs the following
	 *  variables:
	 *  \li \c genoFreq (default): A dictionary (by loci indexes) of
	 *       default dictionaries (by genotype) of genotype frequencies. For
	 *       example, <tt>genoFreq[1][(1, 0)]</tt> is the frequency of genotype
	 *       (1, 0) at locus 1.
	 *  \li \c genoNum (default): A dictionary of default dictionaries of
	 *       genotype counts of all or specified (virtual) subpopulations.
	 *  \li \c genoFreq_sp: genotype frequency in each specified (virtual)
	 *      subpopulation.
	 *  \li \c genoFreq_sp: genotype count in each specified (virtual)
	 *      subpopulation.
	 *
	 *  <b>haploFreq</b>: This parameter accepts one or more lists of loci (by
	 *  index) at which number and frequency of haplotypes are outputted as
	 *  default dictionaries. <tt>[(1,2)]</tt> can be abbreviated to
	 *  <tt>(1,2)</tt>. For example, using parameter <tt>haploFreq=(1,2,4)</tt>,
	 *  all haplotypes at loci \c 1, \c 2 and \c 4 are counted. This statistic
	 *  saves results to dictionary (with loci index as keys) of default
	 *  dictionaries (with haplotypes as keys) such as
	 *  <tt>haploFreq[(1,2,4)][(1,1,0)]</tt> (frequency of haplotype
	 *  <tt>(1,1,0)</tt> at loci <tt>(1,2,3)</tt>). This statistic works for
	 *  all population types. Number of haplotypes for each individual equals
	 *  to his/her ploidy number. Haplodiploid populations are supported in
	 *  the sense that the second homologous copy of the haplotype is not
	 *  counted for male individuals. This statistic outputs the following
	 *  variables:
	 *  \li \c haploFreq (default): A dictionary (with tuples of loci indexes
	 *       as keys) of default dictionaries of haplotype frequencies. For
	 *       example, <tt>haploFreq[(0, 1)][(1,1)]</tt> records the frequency
	 *       of haplotype <tt>(1,1)</tt> at loci <tt>(0, 1)</tt> in all or
	 *       specified (virtual) subpopulations.
	 *  \li \c haploNum (default): A dictionary of default dictionaries of
	 *       haplotype counts in all or specified (virtual) subpopulations.
	 *  \li \c haploFreq_sp: Halptype frequencies in each (virtual)
	 *       subpopulation.
	 *  \li \c haploNum_sp: Halptype count in each (virtual) subpopulation.
	 *
	 *  <b>haploHeteroFreq</b> and <b>haploHomoFreq</b>: These parameters accept
	 *  a list of haplotypes (list of loci), at which the number and frequency of
	 *  haplotype homozygotes and/or heterozygotes will be calculated. Note that
	 *  these statistics are \b observed count of haplotype heterozygote. The
	 *  following variables will be outputted:
	 *  \li \c haploHeteroFreq (default for parameter \e haploHeteroFreq): A
	 *       dictionary of proportion of haplotype heterozygotes in all or
	 *       specified (virtual) subpopulations, with haplotype indexes as
	 *       dictionary keys.
	 *  \li \c haploHomoFreq (default for parameter \e haploHomoFreq): A
	 *       dictionary of proportion of homozygotes in all or specified
	 *       (virtual) subpopulations.
	 *  \li \c haploHeteroNum: A dictionary of number of heterozygotes in all
	 *       or specified (virtual) subpopulations.
	 *  \li \c haploHomoNum: A dictionary of number of homozygotes in all or
	 *       specified (virtual) subpopulations.
	 *  \li \c haploHeteroFreq_sp: A dictionary of proportion of heterozygotes
	 *       in each (virtual) subpopulation.
	 *  \li \c haploHomoFreq_sp: A dictionary of proportion of homozygotes in
	 *       each (virtual) subpopulation.
	 *  \li \c haploHeteroNum_sp: A dictionary of number of heterozygotes in
	 *       each (virtual) subpopulation.
	 *  \li \c haploHomoNum_sp: A dictionary of number of homozygotes in each
	 *       (virtual) subpopulation.
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
	 *  <b>association</b>: Parameter \c association accepts a list of loci,
	 *  which can be a list of indexes, names, or \c ALL_AVAIL. At each locus,
	 *  one or more statistical tests will be performed to test association
	 *  between this locus and individual affection status. Currently,
	 *  simuPOP provides the following tests:
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
	 *  of natural selection) on specified loci, which can be a list of loci
	 *  indexes, names or \c ALL_AVAIL. It currently only outputs \e Pi, which
	 *  is the average number of pairwise difference between loci. This
	 *  statistic outputs the following variables:
	 *  \li \c Pi Mean pairwise difference between all sequences from all or
	 *       specified (virtual) subpopulations.
	 *  \li \c Pi_sp Mean paiewise difference between all sequences in each
	 *       (virtual) subpopulation.
	 *
	 *  <b>structure</b>: Parameter \c structure accepts a list of loci at
	 *  which statistics that measure population structure are calculated.
	 *  \e structure accepts a list of loci indexes, names or \c ALL_AVAIL.
	 *  This parameter currently supports the following statistics:
	 *  \li Weir and Cockerham's Fst (1984). This is the most widely used
	 *       estimator of Wright's fixation index and can be used to measure
	 *       Population differentiation. However, this method is designed to
	 *       estimate Fst from samples of larger populations and might not be
	 *       appropriate for the calculation of Fst of large populations.
	 *  \li Nei's Gst (1973). The Gst estimator is another estimator for
	 *       Wright's fixation index but it is extended for multi-allele (more
	 *       than two alleles) and multi-loci cases. This statistics should
	 *       be used if you would like to obtain a \e true Fst value of a large
	 *       Population.
	 *  Nei's Gst uses only allele frequency information so it is available
	 *  for all population type (haploid, diploid etc). Weir and Cockerham's
	 *  Fst uses heterozygosity frequency so it is best for autosome of
	 *  diploid populations. For non-diploid population, sex, and mitochondrial
	 *  DNAs, simuPOP uses expected heterozygosity (1 - sum p_i^2) when
	 *  heterozygosity is needed. These statistics output the following
	 *  variables:
	 *
	 *  \li \c F_st (default) The WC84 \e Fst statistic estimated for all *       specified loci.
	 *  \li \c F_is The WC84 \e Fis statistic estimated for all specified loci.
	 *  \li \c F_it The WC84 \e Fit statistic estimated for all specified loci.
	 *  \li \c f_st A dictionary of locus level WC84 \e Fst values.
	 *  \li \c f_is A dictionary of locus level WC84 \e Fis values.
	 *  \li \c f_it A dictionary of locus level WC84 \e Fit values.
	 *  \li \c G_st Nei's Gst statistic estimated for all specified loci.
	 *  \li \c g_st A dictionary of Nei's Gst statistic estimated for each
	 *	     locus.
	 *
	 *  <b>HWE</b>: Parameter \c HWE accepts a list of loci at which exact
	 *  two-side tests for Hardy-Weinberg equilibrium will be performed. This
	 *  statistic is only available for diallelic loci in diploid populations.
	 *  \e HWE can be a list of loci indexes, names or \c ALL_AVAIL. This
	 *  statistic outputs the following variables:
	 *  \li \c HWE (default) A dictionary of p-values of HWE tests using
	 *       genotypes in all or specified (virtual) subpopulations.
	 *  \li \c HWE_sp A dictionary of p-values of HWS tests using genotypes
	 *       in each (virtual) subpopulation.
	 *
	 *  <b>inbreeding</b>: Inbreeding measured by Identitcal by Decent (and by
	 *  State). This statistics go through all loci of individuals in a diploid
	 *  population and calculate the number and proportions of alleles that are
	 *  identitcal by decent and by state. Because ancestral information is only
	 *  available in lineage module, variables IBD_freq are always set to zero
	 *  in other modules. Loci on sex and mitochondrial chromosomes, and
	 *  non-diploid populations are currently not supported. This statistic
	 *  outputs the following variables:
	 *  \li \c IBD_freq (default) The frequency of IBD pairs among all allele pairs.
	 *       To use this statistic, the population must be initialized by
	 *       operator InitLineage() to assign each ancestral allele an unique
	 *       identify.
	 *  \li \c IBS_freq (default) The proportion of IBS pairs among all allele pairs.
	 *  \li \c IBD_freq_sp frequency of IBD in each (virtual) subpopulations.
	 *  \li \c IBS_freq_sp frequency of IBS in each (virtual) subpopulations.
	 *
	 *  <b>effectiveSize</b>: Parameter \c effectiveSize accepts a list of loci
	 *  at which the effective population size for the whole or specified
	 *  (virtual) subpopulations is calculated. \e effectiveSize can be a list
	 *  of loci indexes, names or \c ALL_AVAIL. Parameter \e subPops is usually
	 *  used to define samples from which effective sizes are estimated. This
	 *  statistic allows the calculation of true effective size based on
	 *  number of gametes each parents transmit to the offspring population
	 *  (per-locus before and after mating), and estimated effective size
	 *  based on sample genotypes. Due to the temporal natural of some methods,
	 *  more than one Stat operators might be needed to calculate effective
	 *  size. The \e vars parameter specified which method to use and which
	 *  variable to set. Acceptable values include:
	 *  \li \c Ne_demo_base When this variable is set before mating, it stores
	 *       IDs of breeding parents and, more importantly, assign an unique
	 *       lineage value to alleles at specified loci of each individual.
	 *       <b>This feature is only available for lineage modules and will
	 *       change lineage values at specified loci of all individuals</b>.
	 *  \li \c Ne_demo_base_sp Pre-mating information for each (virtual)
	 *       subpopulation, used by variable \c Ne_demo_sp.
	 *  \li \c Ne_demo A dictionary of locus-specific demographic effective
	 *       population size, calculated using number of gemetes each parent
	 *       transmits to the offspring population. The method is vased on
	 *       Crow & Denniston 1988 (Ne = KN-1/k-1+Vk/k) and need variable
	 *       \c Ne_demo_base set before mating. <b>Effective size estimated
	 *       from this formula is model dependent and might not be applicable
	 *       to your mating schemes.</b>
	 *  \li \c Ne_demo_sp Calculate subpopulation-specific effective size.
	 *  \li \c Ne_temporal_base When this variable is set in parameter \e vars,
	 *       the Stat operator saves baseline allele frequencies and other
	 *       information in this variable, which are used by temporary methods
	 *       to estimate effective population size according to changes in
	 *       allele frequency between the baseline and present generations.
	 *       This variable could be set repeatedly to change baselines.
	 *  \li \c Ne_temporal_base_sp Set baseline information for each (virtual)
	 *       subpopulation specified.
	 *  \li \c Ne_tempoFS_P1 Effective population size, 2.5% and 97.5%
	 *       confidence interval for sampling plan 1 as a list of size 3,
	 *       estimated using a temporal method as described in Jorde & Ryman
	 *       (2007), and as implemented by software tempoFS
	 *       (http://www.zoologi.su.se/~ryman/). This variable is set to census
	 *       population size if no baseline has been set, and to the temporal
	 *       effective size between the present and the baseline generation
	 *       otherwise. This method uses population size or sum of
	 *       subpopulation sizes of specified (virtual) subpopulations as
	 *       census population size for the calculation based on plan 1.
	 *  \li \c Ne_tempoFS_P2 Effective population size, 2.5% and 97.5%
	 *       confidence interval for sampling plan 2 as a list of size 6,
	 *       estimated using a temporal method as described in Jorde & Ryman
	 *       (2007). This variable is set to census population size no baseline
	 *       has been set, and to the temporal effective size between the present
	 *       and the baseline generation otherwise. This method assumes that
	 *       the sample is drawn from an infinitely-sized population.
	 *  \li \c Ne_tempoFS deprecated, use \c Ne_tempoFS_P2 instead.
	 *  \li \c Ne_tempoFS_P1_sp Estimate effective size of each (virtual)
	 *       subpopulation using method Jorde & Ryman 2007, assuming sampling
	 *       plan 1. The census population sizes for sampling plan 1 are the
	 *       sizes for each subpopulation that contain the specified (virtual)
	 *       subpopulations.
	 *  \li \c Ne_tempoFS_P2_sp Estimate effective size of each (virtual)
	 *       subpopulation using method Jorde & Ryman 2007, assuming sampling
	 *       plan 2.
	 *  \li \c Ne_tempoFS_sp deprecated, use \c Ne_tempoFS_P2_sp instead.
	 *  \li \c Ne_waples89_P1 Effective population size, 2.5% and 97.5%
	 *       confidence interval for sampling plan 1 as a list of size 6,
	 *       estimated using a temporal method as described in Waples 1989,
	 *       Genetics. Because this is a temporal method, Ne_waples89 estimates
	 *       effective size between the present and the baseline generation set
	 *       by variable \c Ne_temporal_base. Census population size will be
	 *       resutned if no baseline has been set. This method uses population
	 *       size or sum of subpopulation sizes of specified (virtual)
	 *       subpopulations as census population size for the calculation
	 *       based on plan 1.
	 *  \li \c Ne_waples89_P2 Effective population size, 2.5% and 97.5%
	 *       confidence interval for sampling plan 2 as a list of size 6,
	 *       estimated using a temporal method as described in Waples 1989,
	 *       Genetics. Because this is a temporal method, Ne_waples89 estimates
	 *       effective size between the present and the baseline generation
	 *       set by variable \c Ne_temporal_base. Census population size will
	 *       be returned if no baseline has been set.
	 *  \li \c Ne_waples89_P1_sp Estimate effective size for each (virtual)
	 *       subpopulation using method Waples 89, assuming sampling plan 1.
	 *       The census population sizes are the sizes for each subpopulation
	 *       that contain the specified (virtual) subpopulation.
	 *  \li \c Ne_waples89_P2_sp Estimate effective size for each (virtual)
	 *       subpopulation using method Waples 89, assuming sampling plan 2.
	 *  \li \c Ne_waples89_sp deprecated, use \c Ne_waples89_P2_sp instead.
	 *  \li \c Ne_LD Lists of length three for effective population size, 2.5%
	 *       and 97.% confidence interval for cutoff allele frequency 0., 0.01,
	 *       0.02 and 0.05 (as dictionary keys), using a parametric method,
	 *       estimated from linkage disequilibrim information of one sample,
	 *       using LD method developed by Waples & Do 2006 (LDNe). This method
	 *       assumes unlinked loci and uses LD measured from genotypes at loci.
	 *       Because this is a sample based method, it should better be applied
	 *       to a random sample of the population. 95% CI is calculated using a
	 *       Jackknife estimated effective number of independent alleles. Please
	 *       refer to relevant papers and the LDNe user's guide for details.
	 *  \li \c Ne_LD_sp Estimate LD-based effective population size for each
	 *       specified (virtual) subpopulation.
	 *  \li \c Ne_LD_mono A version of Ne_LD that assumes monogamy (see Waples 2006
	 *        for details.
	 *  \li \c Ne_LD_mono_sp Ne_LD_mono calculated for each (virtual) subpopulation.
	 **/
	Stat(bool popSize = false,
		//
		bool numOfMales = false,
		//
		bool numOfAffected = false,
		//
		const lociList & numOfSegSites = vectoru(),
		//
		const lociList & numOfMutants = vectoru(),
		//
		const lociList & alleleFreq = vectoru(),
		//
		const lociList & heteroFreq = vectoru(),
		const lociList & homoFreq = vectoru(),
		//
		const lociList & genoFreq = vectoru(),
		//
		const intMatrix & haploFreq = intMatrix(),
		const intMatrix & haploHeteroFreq = intMatrix(),
		const intMatrix & haploHomoFreq = intMatrix(),
		//
		const stringList & sumOfInfo = vectorstr(),
		const stringList & meanOfInfo = vectorstr(),
		const stringList & varOfInfo = vectorstr(),
		const stringList & maxOfInfo = vectorstr(),
		const stringList & minOfInfo = vectorstr(),
		//
		const intMatrix & LD = intMatrix(),
		//
		const lociList & association = vectoru(),
		//
		const lociList & neutrality = vectoru(),
		//
		const lociList & structure = vectoru(),
		//
		const lociList & HWE = vectoru(),
		//
		const lociList & inbreeding = vectoru(),
		//
		const lociList & effectiveSize = vectoru(),
		//
		const stringList & vars = stringList(),
		const string & suffix = string(),
		// regular parameters
		const stringFunc & output = "",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr());

	~Stat()
	{
	}


	/// HIDDEN
	string describe(bool format = true) const;


	/// HIDDEN Deep copy of a \c Stat operator
	virtual BaseOperator * clone() const
	{
		return new Stat(*this);
	}


	// count various statistics.
	// use m_alleles etc to save (potentially) time to
	// resize all these variables.
	/// HIDDEN apply the \c Stat operator
	virtual bool apply(Population & pop) const;

private:
	const statPopSize m_popSize;
	const statNumOfMales m_numOfMales;
	const statNumOfAffected m_numOfAffected;
	const statNumOfSegSites m_numOfSegSites;
	const statNumOfMutants m_numOfMutants;
	const statAlleleFreq m_alleleFreq;
	const statHeteroFreq m_heteroFreq;
	const statGenoFreq m_genoFreq;
	const statHaploFreq m_haploFreq;
	const statHaploHomoFreq m_haploHomoFreq;
	const statInfo m_info;
	const statLD m_LD;
	const statAssociation m_association;
	const statNeutrality m_neutrality;
	const statStructure m_structure;
	const statHWE m_HWE;
	const statInbreeding m_Inbreeding;
	const statEffectiveSize m_effectiveSize;
};

}
#endif
