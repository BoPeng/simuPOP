/**
 *  $File: operator.h $
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

#ifndef _OPERATOR_H
#define _OPERATOR_H
/**
   \file
   \brief head file of class Operator
 */
#include "utility.h"
#include "simuPOP_cfg.h"

#include <iostream>
#include <iomanip>
using std::endl;

#include <fstream>
using std::ofstream;

#include <string>
using std::string;

// for operator TicToc
#include <time.h>

#include "individual.h"
#include "population.h"

namespace simuPOP {

/** Operators are objects that act on populations. They can be applied to
 *  populations directly using their function forms, but they are usually
 *  managed and applied by a simulator. In the latter case, operators are
 *  passed to the \c evolve function of a simulator, and are applied repeatedly
 *  during the evolution of the simulator.
 *
 *  The \e BaseOperator class is the base class for all operators. It defines
 *  a common user interface that specifies at which generations, at which stage
 *  of a life cycle, to which populations and subpopulations an operator is
 *  applied. These are achieved by a common set of parameters such as \c begin,
 *  \c end, \c step, \c at, \c stage for all operators. Note that a specific
 *  operator does not have to honor all these parameters. For example, a
 *  Recombinator can only be applied during mating so it ignores the \c stage
 *  parameter.
 *
 *  An operator can be applied to all or part of the generations during the
 *  evolution of a simulator. At the beginning of an evolution, a simulator
 *  is usually at the beginning of generation \c 0. If it evolves \c 10
 *  generations, it evolves generations \c 0, \c 1, ,,,., and \c 9 (\c 10
 *  generations) and stops at the begging of generation \c 10. A negative
 *  generation number \c a has generation number <tt>10 + a</tt>, with -1
 *  referring to the last evolved generation \c 9. Note that the starting
 *  generation number of a simulator can be changed by its \c setGen()
 *  member function.
 *
 *  Output from an operator is usually directed to the standard output
 *  (\c sys.stdout). This can be configured using a output specification
 *  string, which can be <tt>''</tt> for no output, <tt>'>'</tt> standard
 *  terminal output (default), a filename prefixed by one or more
 *  <tt>'>'</tt> characters or a Python expression indicated by a leading
 *  exclamation mark (<tt>'!expr'</tt>). In the case of \c '>filename' (or
 *  equivalently \c 'filename'), the output from an operator is written to this
 *  file. However, if two operators write to the same file \c filename, or if
 *  an operator writes to this file more than once, only the last write
 *  operation will succeed. In the case of <tt>'>>filename'</tt>, file
 *  \c filename will be opened at the beginning of the evolution and closed at
 *  the end. Outputs from multiple operators are appended. <tt>>>>filename</tt>
 *  works similar to <tt>>>filename</tt> but \c filename, if it already exists
 *  at the beginning of an evolutionary process, will not be cleared. If the
 *  output specification is prefixed by an exclamation mark, the string after
 *  the mark is considered as a Python expression. When an operator is applied
 *  to a population, this expression will be evaluated within the population's
 *  local namespace to obtain a population specific output specification.
 *  As an advanced feature, a Python function can be assigned to this
 *  parameter. Output strings will be sent to this function for processing.
 *  Lastly, if the output stream only accept a binary output (e.g.
 *  a gzip stream), <tt>WithMode(output, 'b')</tt> should be used to
 *  let simuPOP convert string to bytes before writing to the output.
 */
class BaseOperator
{
public:
	/** @name constructor and destructor */
	//@{

	/** The following parameters can be specified by all operators. However,
	 *  an operator can ignore some parameters and the exact meaning of a
	 *  parameter can vary.
	 *
	 *  \param output A string that specifies how output from an operator is
	 *    written, which can be \c '' (no output), \c '>' (standard output),
	 *    \c 'filename' prefixed by one or more '>', or an Python expression
	 *    prefixed by an exclamation mark (\c '!expr'). If a \c file object, or
	 *    any Python object with a \c write function is provided, the output
	 *    will be write to this file. Alternatively, a Python function or a
	 *    file object (any Python object with a \c write function) can be
	 *    given which will be called with a string of output content. A
	 *    global function \c WithMode can be used to let simuPOP output
	 *    bytes instead of string.
	 *  \param begin The starting generation at which an operator will be
	 *    applied. Default to \c 0. A negative number is interpreted as a
	 *    generation counted from the end of an evolution (-1 being the last
	 *    evolved generation).
	 *  \param end The last generation at which an operator will be applied.
	 *    Default to \c -1, namely the last generation.
	 *  \param step The number of generations between applicable generations.
	 *    Default to \c 1.
	 *  \param at A list of applicable generations. Parameters \c begin,
	 *    \c end, and \c step will be ignored if this parameter is specified.
	 *    A single generation number is also acceptable.
	 *  \param reps A list of applicable replicates. A common default value
	 *    \c ALL_AVAIL is interpreted as all replicates in a simulator. Negative
	 *    indexes such as \c -1 (last replicate) is acceptable. <tt>rep=idx</tt>
	 *    can be used as a shortcut for <tt>rep=[idx]</tt>.
	 *  \param subPops A list of applicable (virtual) subpopulations, such as
	 *    <tt>subPops=[sp1, sp2, (sp2, vsp1)]</tt>. <tt>subPops=[sp1]</tt>
	 *    can be simplied as <tt>subPops=sp1</tt>. Negative indexes are not
	 *    supported. A common default value (\c ALL_AVAIL) of this parameter
	 *    reprents all subpopulations of the population being aplied. Suport
	 *    for this parameter vary from operator to operator and some operators
	 *    do not support virtual subpopulations at all. Please refer to the
	 *    reference manual of individual operators for their support for this
	 *    parameter.
	 *  \param infoFields A list of information fields that will be used by an
	 *    operator. You usually do not need to specify this parameter because
	 *    operators that use information fields usually have default values for
	 *    this parameter.
	 */
	BaseOperator(const stringFunc & output, int begin, int end, int step, const intList & at,
	             const intList & reps, const subPopList & subPops, const stringList & infoFields) :
		m_beginGen(begin), m_endGen(end), m_stepGen(step), m_atGen(at.elems()),
		m_flags(0), m_reps(reps), m_subPops(subPops),
		m_ostream(output.value(), output.func(), output.mode()), m_infoFields(infoFields)
	{
		PARAM_FAILIF(step <= 0, ValueError, "step need to be at least one");

		setFlags();
	}


	/// destroy an operator
	virtual ~BaseOperator()
	{
	}


	/** Return a cloned copy of an operator. This function is available to all
	 *  operators.
	 */
	virtual BaseOperator * clone() const
	{
		return new BaseOperator(*this);
	}


	//@}

	/** @name  applicable generations (also judge from rep). use of parameter start, end, every, at, rep
	 */
	//@{

	/** CPPONLY determine if this operator is active
	   Determine if this operator is active under the conditions such as the current
	   replicate, current generation, ending generation etc.
	   \note This function will be called by Simulators before applying.
	 */
	bool isActive(size_t rep, ssize_t gen, ssize_t end, const vector<bool> & activeRep, bool repOnly = false) const;

	/** CPPONLY Another version of isActive when negative gen is not considered.
	 */
	bool isActive(ssize_t rep, ssize_t gen) const;

	//@}

	/** @name applicable stages	pre, during, post-mating methods */
	//@{


	/// get the length of information fields for this operator
	/// CPPONLY
	size_t infoSize(const GenoStruTrait * trait = NULL) const
	{
		return m_infoFields.elems(trait).size();
	}


	/// get the information field specified by user (or by default)
	/// CPPONLY
	string infoField(size_t idx, const GenoStruTrait * trait = NULL) const
	{
		// suppress a warning about unused variable
		(void)trait;
		DBG_ASSERT(idx < m_infoFields.elems(trait).size(), IndexError,
			(boost::format("Given info index %1% is out of range of 0 ~ %2%") % idx % m_infoFields.elems().size()).str());
		return m_infoFields.elems()[idx];
	}


	/// CPPONLY
	const stringList & infoFields() const
	{
		return m_infoFields;
	}


	/** Apply an operator to population \e pop directly, without checking its
	 *  applicability.
	 */
	virtual bool apply(Population & pop) const;


	/// CPPONLY apply during mating, given \c pop, \c offspring, \c dad and \c mom
	virtual bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;


	//@}
	/** @name dealing with output separator, persistant files, $gen etc substitution.
	 */
	//@{


	/// CPPONLY get output stream. This function is not exposed to user.
	ostream & getOstream(PyObject * dict = NULL, bool readable = false) const
	{
		return m_ostream.getOstream(dict, readable);
	}


	/// CPPONLY close output stream and delete output stream pointer, if it is a output stream
	void closeOstream() const
	{
		m_ostream.closeOstream();
	}


	/// CPPONLY say something about the applicability of this operator.
	string applicability(bool subPops = true, bool gen = true) const;


	/// HIDDEN
	virtual string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.operator> a based operator that should not be used directly." ;
	}


	//@}

	/// CPPONLY determine if \c output=">". Used internally.
	bool noOutput() const
	{
		return m_ostream.noOutput();
	}


	/// CPPONLY
	subPopList applicableSubPops(const Population & pop) const
	{
		return m_subPops.expandFrom(pop);
	}


	/// CPPONLY
	bool applicableToAllOffspring() const
	{
		return ISSETFLAG(m_flags, m_flagAllSubPops);
	}


	/// CPPONLY
	bool applicableToOffspring(const Population & pop, RawIndIterator offspring) const;


	/// CPPONLY
	virtual bool parallelizable() const
	{
		return false;
	}


	/// CPPONLY
	virtual void initialize(const Individual & ind) const
	{
		(void)ind;  // avoid warning about unused parameter
	}


	/// CPPONLY
	virtual void initializeIfNeeded(const Individual & ind) const
	{
		(void)ind;  // avoid warning about unused parameter
	}


protected:
	/// analyze active generations: set m_flagAtAllGen etc
	void setFlags();

private:
	/// internal m_flags of the operator. They are set during initialization for
	/// performance considerations.
	static const unsigned char m_flagAtAllGen = 1;
	static const unsigned char m_flagOnlyAtBegin = 2;
	static const unsigned char m_flagOnlyAtEnd = 4;
	static const unsigned char m_flagAllSubPops = 8;

private:
	/// starting generation, default to 0
	const int m_beginGen;

	/// ending generation, default to -1 (last genrration)
	const int m_endGen;

	/// interval between active states, default to 1 (active at every generation)
	const int m_stepGen;

	/// a list of generations that this oeprator will be active.
	/// typical usage is m_atGen=-1 to apply at the last generation.
	const vectori m_atGen;

	/// various m_flags of Operator for faster processing.
	mutable unsigned char m_flags;

	/// apply to all (-1) or one of the replicates.
	const intList m_reps;

	/// apply to some of the (virtual) subpops.
	const subPopList m_subPops;

	/// the output stream
	mutable StreamProvider m_ostream;

	/// information fields that will be used by this operator
	const stringList m_infoFields;
};

typedef std::vector< BaseOperator * > vectorop;

class opList
{
public:
	typedef vectorop::const_iterator const_iterator;

public:
	opList(PyObject * obj = NULL);

	/// CPPONLY
	opList(const vectorop & ops);

	/// CPPONLY
	opList(const opList & rhs);

	const BaseOperator * operator[](size_t idx) const
	{
		DBG_ASSERT(idx < m_elems.size(), IndexError,
			"Operator index out of range");
		return m_elems[idx];
	}


	~opList();


	/// CPPONLY
	const_iterator begin() const
	{
		return m_elems.begin();
	}


	/// CPPONLY
	const_iterator end() const
	{
		return m_elems.end();
	}


	/// CPPONLY
	size_t size() const
	{
		return m_elems.size();
	}


	/// CPPONLY
	bool empty() const
	{
		return m_elems.empty();
	}


	/// CPPONLY
	const vectorop & elems() const
	{
		return m_elems;
	}


protected:
	vectorop m_elems;
};

/** This operator pauses the evolution of a simulator at given generations or
 *  at a key stroke. When a simulator is stopped, you can go to a Python
 *  shell to examine the status of an evolutionary process, resume or stop the
 *  evolution.
 */
class Pause : public BaseOperator
{

public:
	/** Create an operator that pause the evolution of a population when it is
	 *  applied to this population. If \e stopOnKeyStroke is \c False (default),
	 *  it will always pause a population when it is applied, if this parameter
	 *  is set to \c True, the operator will pause a population if *any* key
	 *  has been pressed. If a specific character is set, the operator will stop
	 *  when this key has been pressed. This allows, for example, the use of
	 *  several pause operators to pause different populations.
	 *
	 *  After a population has been paused, a message will be displayed (unless
	 *  \e prompt is set to \c False) and tells you how to proceed. You can
	 *  press \c 's' to stop the evolution of this population, \c 'S' to
	 *  stop the evolution of all populations, or \c 'p' to enter a Python
	 *  shell. The current population will be available in this Python shell
	 *  as \c "pop_X_Y" when \c X is generation number and \c Y is replicate
	 *  number. The evolution will continue after you exit this interactive
	 *  Python shell.
	 *
	 *  \note Ctrl-C will be intercepted even if a specific character is
	 *  specified in parameter \e stopOnKeyStroke.
	 */
	Pause(char stopOnKeyStroke = false, bool prompt = true,
	      const stringFunc & output = ">",
	      int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
	      const intList & reps = intList(), const subPopList & subPops = subPopList(),
	      const stringList & infoFields = vectorstr()) :
		BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
		m_prompt(prompt), m_stopOnKeyStroke(stopOnKeyStroke)
	{
		(void)output;  // avoid warning about unused parameter
	}


	/// destructor
	virtual ~Pause()
	{
	}


	/// HIDDEN
	virtual BaseOperator * clone() const
	{
		return new Pause(*this);
	}


	/// HIDDEN apply the \c Pause operator to one population
	bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const;

private:
	const bool m_prompt;

	const char m_stopOnKeyStroke;

	static vectori s_cachedKeys;
};

/** This operator does nothing when it is applied to a population. It is
 *  usually used as a placeholder when an operator is needed syntactically.
 */
class NoneOp : public BaseOperator
{

public:
	/** Create a \c NoneOp.
	 */
	NoneOp(const stringFunc & output = ">",
	       int begin = 0, int end = 0, int step = 1, const intList & at = vectori(),
	       const intList & reps = intList(), const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr()) :
		BaseOperator("", begin, end, step, at, reps, subPops, infoFields)
	{
		(void)output;  // avoid warning about unused parameter
	}


	/// destructor
	virtual ~NoneOp()
	{
	}


	/// HIDDEN
	virtual BaseOperator * clone() const
	{
		return new NoneOp(*this);
	}


	/// CPPONLY apply during mating, given \c pop, \c offspring, \c dad and \c mom
	virtual bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
	                               Individual * dad = NULL, Individual * mom = NULL) const
	{
		(void)pop;          // avoid warning about unused parameter
		(void)offPop;       // avoid warning about unused parameter
		(void)offspring;    // avoid warning about unused parameter
		(void)dad;          // avoid warning about unused parameter
		(void)mom;          // avoid warning about unused parameter
		return true;
	}


	/// HIDDEN apply the \c NoneOp operator to one population
	virtual bool apply(Population & /* pop */) const
	{
		return true;
	}


	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.None> an operator that does nothing" ;
	}


};

/** This operator uses a condition, which can be a fixed condition, an
 *  expression or a user-defined function, to determine which operators to be
 *  applied when this operator is applied. A list of if-operators will be
 *  applied when the condition is \c True. Otherwise a list of else-operators
 *  will be applied.
 */
class IfElse : public BaseOperator
{

public:
	/** Create a conditional operator that will apply operators \e ifOps if
	 *  condition \e cond is met and \e elseOps otherwise. If a Python
	 *  expression (a string) is given to parameter \e cond, the expression
	 *  will be evalulated in each population's local namespace when this
	 *  operator is applied. When a Python function is specified, it accepts
	 *  parameter \c pop when it is applied to a population, and one or more
	 *  parameters \c pop, \c off, \c dad or \c mom when it is applied during
	 *  mating. The return value of this function should be \c True or \c False.
	 *  Otherwise, parameter \e cond will be treated as a fixed condition
	 *  (converted to \c True or \c False) upon which one set of operators is
	 *  always applied. The applicability of \e ifOps and \e elseOps are
	 *  controlled by parameters \e begin, \e end, \e step,  \e at and \e rep
	 *  of both the \c IfElse operator and individual operators but \e ifOps
	 *  and \e elseOps opeartors does not support negative indexes
	 *  for replicate and generation numbers.
	 */
	IfElse(PyObject * cond, const opList & ifOps = opList(), const opList & elseOps = opList(),
		const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr());

	/// destructor
	virtual ~IfElse()
	{
	}


	/// HIDDEN
	virtual BaseOperator * clone() const
	{
		return new IfElse(*this);
	}


	/// CPPONLY apply during mating, given \c pop, \c offspring, \c dad and \c mom
	virtual bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// HIDDEN apply the \c IfElse operator to population \e pop.
	virtual bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const;

private:
	/// These will be kept constant (they are set in constructor only)
	Expression m_cond;
	pyFunc m_func;
	int m_fixedCond;

	const opList m_ifOps;
	const opList m_elseOps;
};


/** This operator evaluates an expression in a population's local namespace
 *  and terminate the evolution of this population, or the whole simulator,
 *  if the return value of this expression is \c True. Termination caused by
 *  an operator will stop the execution of all operators after it. The
 *  generation at which the population is terminated will be counted in the
 *  <em>evolved generations</em> (return value from <tt>Simulator::evolve</tt>)
 *  if termination happens after mating.
 */
class TerminateIf : public BaseOperator
{

public:
	/** Create a terminator with an expression \e condition, which will be evalulated
	 *  in a population's local namespace when the operator is applied to this
	 *  population. If the return value of \e condition is \c True, the evolution
	 *  of the population will be terminated. If \e stopAll is set to \c True, the
	 *  evolution of all replicates of the simulator will be terminated. If this
	 *  operator is allowed to write to an \e output (default to ""), the generation
	 *  number, proceeded with an optional \e message.
	 */
	TerminateIf(string condition = string(), bool stopAll = false, string message = string(),
	            const stringFunc & output = "", int begin = 0, int end = -1,
	            int step = 1, const intList & at = vectori(), const intList & reps = intList(),
	            const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr()) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_expr(condition), m_stopAll(stopAll), m_message(message)
	{
	}


	/// HIDDEN Deep copy of a \c TerminateIf terminator
	virtual BaseOperator * clone() const
	{
		return new TerminateIf(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;


	/// HIDDEN check all alleles in vector allele if they are fixed.
	bool apply(Population & pop) const;

	virtual ~TerminateIf()
	{
	}


private:
	/// alleles to check. If empty, check all alleles.
	const Expression m_expr;

	///
	const bool m_stopAll;

	/// message to print when terminated
	const string m_message;
};


/** This operator replaces the current evolving population by a population
 *  loaded from a specified filename if certain condition is met. It is mostly
 *  used to return to a previously saved state if the simulation process
 *  fails to met a condition (e.g. a disease allele is lost).
 */
class RevertIf : public BaseOperator
{

public:
	/** Replaces the current evolving population by a population loaded from
	 *  \e fromPop, which should be a file saved by function \c Population.save()
	 *  or operator \c SavePopulation. If a Python expression (a string) is
	 *  given to parameter \e cond, the expression will be evalulated in each
	 *  population's local namespace when this operator is applied. When a Python
	 *  function with optional parameter \c pop is specified, it should accept
	 *  the current population (to parameter \c pop) and converts and return
	 *  \c True or \c False. Otherwise, parameter \e cond will be treated as a
	 *  fixed condition (converted to \c True or \c False) upon which the population
	 *  is reverted. After the population is reverted, an optional list of
	 *  operators \e ops could be applied to the population.
	 */
	RevertIf(PyObject * cond, const string & fromPop = string(), const opList & ops = opList(),
		const stringFunc & output = "", int begin = 0, int end = -1,
		int step = 1, const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr());

	/// HIDDEN Deep copy of a \c RevertIf terminator
	virtual BaseOperator * clone() const
	{
		return new RevertIf(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;


	/// HIDDEN check all alleles in vector allele if they are fixed.
	bool apply(Population & pop) const;

	virtual ~RevertIf()
	{
	}


private:
	///
	Expression m_cond;
	pyFunc m_func;
	int m_fixedCond;

	/// message to print when terminated
	const string m_fromPop;

	const opList m_ops;
};


/** This operator discards individuals according to either an expression that
 *  evaluates according to individual information field, or a Python function
 *  that accepts individual and its information fields.
 */
class DiscardIf : public BaseOperator
{

public:
	/** Create an operator that discard individuals according to an expression
	 *  or the return value of a Python function (parameter \e cond). This
	 *  operator can be applied to a population before or after mating, or
	 *  to offspring during mating. If an expression is passed to \e cond, it
	 *  will be evalulated with each individual's information fields (see
	 *  operator \c InfoEval for details). If \e exposeInd is non-empty,
	 *  individuals will be available for evaluation in the expression as an
	 *  variable with name spacied by \e exposeInd. If the expression is
	 *  evaluated to be \c True, individuals (if applied before or after
	 *  mating) or offspring (if applied during mating) will be removed or
	 *  discard. Otherwise the return value should be either \c False (not
	 *  discard), or a float number between \c 0 and \c 1 as the probability
	 *  that the individual is removed. If a function is passed to \e cond, it
	 *  should accept paramters \e ind and \e pop or names of information
	 *  fields when it is applied to a population (pre or post mating), or
	 *  parameters \e off, \e dad, \e mom, \e pop (parental population), or
	 *  names of information fields if the operator is applied during mating.
	 *  Individuals will be discarded if this function returns \c True or
	 *  at a probability if a float number between 0 and 1 is returned. A
	 *  constant expression (e.g. \c True, \c False, \c 0.4) is also
	 *  acceptable, with the last example (\c cond=0.1) that removes
	 *  10% of individuals at randomly. This operator supports parameter
	 *  \e subPops and will remove only individuals belonging to specified
	 *  (virtual) subpopulations.
	 */
	DiscardIf(PyObject * cond, const string & exposeInd = string(),
		const stringFunc & output = "", int begin = 0, int end = -1,
		int step = 1, const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr());


	/// HIDDEN Deep copy of a \c DiscardIf terminator
	virtual BaseOperator * clone() const
	{
		return new DiscardIf(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;

	/// HIDDEN
	virtual bool apply(Population & pop) const;

	/// CPPONLY
	bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad, Individual * mom) const;

	virtual ~DiscardIf()
	{
		Py_XDECREF(m_dict);
	}


private:
	Expression m_cond;
	pyFunc m_func;
	double m_fixedCond;
	string m_fixedField;

	const string m_exposeInd;
	mutable PyObject * m_dict;
	mutable vectorf m_lastValues;
};


/** This operator, when called, output the difference between current and the
 *  last called clock time. This can be used to estimate execution time of
 *  each generation. Similar information can also be obtained from
 *  <tt>turnOnDebug("DBG_PROFILE")</tt>, but this operator has the advantage
 *  of measuring the duration between several generations by setting \c step
 *  parameter. As an advanced feature that mainly used for performance testing,
 *  this operator accepts a parameter \e stopAfter (seconds), and will stop the
 *  evolution of a population if the overall time exceeds \e stopAfter. Note
 *  that elapsed time is only checked when this operator is applied to a
 *  population so it might not be able to stop the evolution process right
 *  after \e stopAfter seconds. This operator can also be applied during
 *  mating. Note that to avoid excessive time checking, this operator does
 *  not always check system time accurately.
 */
class TicToc : public BaseOperator
{
public:
	/** Create a \c TicToc operator that outputs the elapsed since the last
	 *  time it was applied, and the overall time since the first time this
	 *  operator is applied.
	 */
	TicToc(const stringFunc & output = ">", const double stopAfter = 0,
	       int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
	       const intList & reps = intList(), const subPopList & subPops = subPopList(),
	       const stringList & infoFields = vectorstr()) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_stopAfter(stopAfter), m_startTime(0), m_lastTime(0),
		m_totalCount(0), m_countPerSec(0), m_counter(0)
	{
	}


	/// destructor
	virtual ~TicToc()
	{
	}


	/// HIDDEN
	virtual BaseOperator * clone() const
	{
		return new TicToc(*this);
	}


	/// HIDDEN
	virtual bool apply(Population & pop) const;

	/// CPPONLY
	virtual bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// HIDDEN
	string describe(bool format = true) const;


	/// CPPONLY
	bool parallelizable() const
	{
		return true;
	}


private:
	const double m_stopAfter;
	mutable time_t m_startTime;
	mutable time_t m_lastTime;
	mutable ULONG m_totalCount;
	mutable ULONG m_countPerSec;
	mutable ULONG m_counter;
};


/** An operator that calls a user-defined function when it is applied to a
 *  population (pre- or post-mating) or offsprings (during-mating). The
 *  function can have have parameters \c pop when the operator is applied
 *  pre- or post-mating, <tt>pop, off, dad, mom</tt> when the operator is
 *  applied during-mating. An optional parameter can be passed if parameter
 *  \e param is given. In the during-mating case, parameters \c pop,
 *  \c dad and \c mom can be ignored if \e offspringOnly is set to \c True.
 */
class PyOperator : public BaseOperator
{
public:
	/** Create a pure-Python operator that calls a user-defined function when
	 *  it is applied. If this operator is applied before or after mating,
	 *  your function should have form <tt>func(pop)</tt> or
	 *  <tt>func(pop, param)</tt> where \c pop is the population to which
	 *  the operator is applied, \c param is the value specified in parameter
	 *  \e param. \c param will be ignored if your function only accepts
	 *  one parameter. Althernatively, the function should have form
	 *  <tt>func(ind)</tt> with optional parameters \c pop and \c param.
	 *  In this case, the function will be called for all individuals, or
	 *  individuals in subpopulations \e subPops. Individuals for which the
	 *  function returns \c False will be removed from the population. This
	 *  operator can therefore perform similar functions as operator
	 *  \c DiscardIf.
	 *
	 *  If this operator is applied during mating, your function should accept
	 *  parameters \c pop, \c off (or \c ind), \c dad, \c mom and \c param where
	 *  \c pop is the parental population, and \c off or \c ind, \c dad, and
	 *  \c mom are offspring and their parents for each mating event, and
	 *  \c param is an optional parameter. If \e subPops are provided, only
	 *  offspring in specified (virtual) subpopulations are acceptable.
	 *
	 *  This operator does not support parameters \e output, and
	 *  \e infoFields. If certain output is needed, it should be handled in the
	 *  user defined function \e func. Because the status of files used by
	 *  other operators through parameter \e output is undetermined during
	 *  evolution, they should not be open or closed in this Python operator.
	 */
	PyOperator(PyObject * func, PyObject * param = NULL,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr());

	/// HIDDEN
	virtual BaseOperator * clone() const
	{
		return new PyOperator(*this);
	}


	/** HIDDEN Apply the \c PyOperator operator to population \e pop. Calling this
	 *  function is equivalent to call \e func with parameter \e pop and
	 *  optional parameter \e param.
	 */
	virtual bool apply(Population & pop) const;

	/// CPPONLY
	virtual bool applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
		Individual * dad = NULL, Individual * mom = NULL) const;

	/// HIDDEN
	string describe(bool format = true) const;

private:
	/// the function
	const pyFunc m_func;

	/// parammeters
	const pyObject m_param;
};


/** HIDDEN
 *  This function is used to test during mating operators. It simply apply
 *  operator \e op to \e dad, \e mom to offspring \e off through \e off1
 *  of population \e pop. If index of dad or mom is negative, NULL will be passed.
 */
void applyDuringMatingOperator(const BaseOperator & op,
	Population * pop, Population * offPop, ssize_t dad, ssize_t mom, const pairu & off);

}
#endif
