/**
 *  $File: operator.h $
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

// for operator ticToc
#include <time.h>

#include "individual.h"
#include "population.h"

namespace simuPOP {

/** A class to specify (virtual) subpopulation list. Using a dedicated class
 *  allows users to specify a single subpopulation, or a list of (virutal)
 *  subpoulations easily.
 */
class subPopList
{
public:
	// for some unknown reason, std:: is required for this type to be recognized
	// by swig.
	typedef std::vector<vspID> vectorvsp;
	typedef vectorvsp::iterator iterator;
	typedef vectorvsp::const_iterator const_iterator;

public:
	///
	subPopList(PyObject * obj = NULL);

	/// CPPONLY
	subPopList(const vectorvsp & subPops);

	/// CPPONLY
	bool allAvail()
	{
		return m_allAvail;
	}


	/// CPPONLY
	bool empty() const
	{
		return m_subPops.empty();
	}


	/// CPPONLY
	size_t size() const
	{
		return m_subPops.size();
	}


	int __len__() const
	{
		return m_subPops.size();
	}


	/// CPPONLY
	vspID operator[](unsigned int idx) const
	{
		DBG_FAILIF(idx >= m_subPops.size(), IndexError,
			"Index out of range.");
		return m_subPops[idx];
	}


	/// CPPONLY
	void push_back(const vspID subPop)
	{
		m_subPops.push_back(subPop);
	}


	/// CPPONLY
	bool contains(const vspID subPop) const
	{
		return find(m_subPops.begin(), m_subPops.end(), subPop) != m_subPops.end();
	}


	/// CPPONLY
	const_iterator begin() const
	{
		return m_subPops.begin();
	}


	/// CPPONLY
	const_iterator end() const
	{
		return m_subPops.end();
	}


	/// CPPONLY
	iterator begin()
	{
		return m_subPops.begin();
	}


	/// CPPONLY
	iterator end()
	{
		return m_subPops.end();
	}


	///  CPPONLY If a subPopList is invalid (none), it will not be expanded.
	void useSubPopsFrom(const population & pop)
	{
		DBG_FAILIF(m_allAvail && !m_subPops.empty(), SystemError,
			"Only when no subpopulation is specified can this function be called."
			"This is likely caused by the use of persistent subPops for different populations.");
		if (m_allAvail)
			for (size_t sp = 0; sp < pop.numSubPop(); ++sp)
				m_subPops.push_back(vspID(sp));
	}


private:
	vectorvsp m_subPops;
	bool m_allAvail;
};


/** Operators are objects that act on populations. They can be applied to
 *  populations directly using their function forms, but they are usually
 *  managed and applied by a simulator. In the latter case, operators are
 *  passed to the \c evolve function of a simulator, and are applied repeatedly
 *  during the evolution of the simulator.
 *
 *  The \e baseOperator class is the base class for all operators. It defines
 *  a common user interface that specifies at which generations, at which stage
 *  of a life cycle, to which populations and subpopulations an operator is
 *  applied. These are achieved by a common set of parameters such as \c begin,
 *  \c end, \c step, \c at, \c stage for all operators. Note that a specific
 *  operator does not have to honor all these parameters. For example, a
 *  recombinator can only be applied during mating so it ignores the \c stage
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
 */
class baseOperator
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
	 *    prefixed by an exclamation mark (\c '!expr'). Alternatively, a
	 *    Python function can be given to handle outputs.
	 *  \param stage Stage(s) of a life cycle at which an operator will be
	 *    applied. It can be \c PreMating, \c DuringMating, \c PostMating or
	 *    any of their combined stages \c PrePostMating, \c PreDuringMating,
	 *    \c DuringPostMating and \c PreDuringPostMating. Note that all
	 *    operators have their default stage parameter and some of them ignore
	 *    this parameter because they can only be applied at certain stage(s)
	 *    of a life cycle.
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
	 *    \c AllAvail is interpreted as all replicates in a simulator. Negative
	 *    indexes such as \c -1 (last replicate) is acceptable. <tt>rep=idx</tt>
	 *    can be used as a shortcut for <tt>rep=[idx]</tt>.
	 *  \param subPops A list of applicable (virtual) subpopulations, such as
	 *    <tt>subPops=[sp1, sp2, (sp2, vsp1)]</tt>. <tt>subPops=[sp1]</tt>
	 *    can be simplied as <tt>subPops=sp1</tt>. Negative indexes are not
	 *    supported. A common default value (\c AllAvail) of this parameter
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
	baseOperator(const stringFunc & output, int stage, int begin, int end, int step, const intList & at,
		const intList & rep, const subPopList & subPops, const stringList & infoFields) :
		m_beginGen(begin), m_endGen(end), m_stepGen(step), m_atGen(at.elems()),
		m_flags(0), m_rep(rep), m_subPops(subPops),
		m_ostream(output.value(), output.func()), m_infoFields(infoFields.elems()),
		m_lastPop(MaxTraitIndex)
	{
		DBG_FAILIF(step <= 0, ValueError, "step need to be at least one");

		setApplicableStage(stage);
		setFlags();
	}


	/// destroy an operator
	virtual ~baseOperator()
	{
	}


	/** Return a cloned copy of an operator. This function is available to all
	 *  operators.
	 */
	virtual baseOperator * clone() const
	{
		return new baseOperator(*this);
	}


	//@}

	/** @name  applicable generations (also judge from rep). use of parameter start, end, every, at, rep
	 */
	//@{

	/// CPPONLY determine if this operator is active
	/**
	   Determine if this operator is active under the conditions such as the current
	   replicate, current generation, ending generation etc.
	   \note This function will be called by simulators before applying.
	 */
	bool isActive(UINT rep, long gen, long end, const vector<bool> & activeRep, bool repOnly = false);


	//@}

	/** @name applicable stages	pre, during, post-mating methods */
	//@{

	/// set applicable stage. Another way to set \c stage parameter.
	/// CPPONLY
	void setApplicableStage(int stage)
	{
		RESETFLAG(m_flags, PreDuringPostMating);
		SETFLAG(m_flags, stage);
	}


	/// set if this operator can be applied \em pre-mating
	/// CPPONLY
	bool canApplyPreMating() const
	{
		return ISSETFLAG(m_flags, m_flagPreMating);
	}


	/// set if this operator can be applied \em during-mating
	/// CPPONLY
	bool canApplyDuringMating() const
	{
		return ISSETFLAG(m_flags, m_flagDuringMating);
	}


	/// set if this operator can be applied \em post-mating
	/// CPPONLY
	bool canApplyPostMating() const
	{
		return ISSETFLAG(m_flags, m_flagPostMating);
	}


	/// CPPONLY
	virtual bool isCompatible(const population & pop)
	{
		return true;
	}


	/// determine if the operator can be applied only for haploid population
	/// CPPONLY
	bool haploidOnly()
	{
		return ISSETFLAG(m_flags, m_flagHaploid);
	}


	/// determine if the operator can be applied only for diploid population
	/// CPPONLY
	bool diploidOnly()
	{
		return ISSETFLAG(m_flags, m_flagDiploid);
	}


	/// CPPONLY set that the operator can be applied only for haploid populations
	void setHaploidOnly()
	{
		SETFLAG(m_flags, m_flagHaploid);
	}


	/// CPPONLY set that the operator can be applied only for diploid populations
	void setDiploidOnly()
	{
		SETFLAG(m_flags, m_flagDiploid);
	}


	/// get the length of information fields for this operator
	/// CPPONLY
	UINT infoSize()
	{
		return m_infoFields.size();
	}


	/// get the information field specified by user (or by default)
	/// CPPONLY
	string infoField(UINT idx)
	{
		DBG_ASSERT(idx < m_infoFields.size(), IndexError, "Given info index " + toStr(idx) +
			" is out of range of 0 ~ " + toStr(m_infoFields.size()));
		return m_infoFields[idx];
	}


	/** Apply an operator to population \e pop directly, without checking its
	 *  applicability.
	 */
	virtual bool apply(population & pop);


	/// CPPONLY apply during mating, given \c pop, \c offspring, \c dad and \c mom
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);


	//@}
	/** @name dealing with output separator, persistant files, $gen etc substitution.
	 */
	//@{


	/// CPPONLY get output stream. This function is not exposed to user.
	ostream & getOstream(PyObject * dict = NULL, bool readable = false)
	{
		return m_ostream.getOstream(dict, readable);
	}


	/// CPPONLY close output stream and delete output stream pointer, if it is a output stream
	void closeOstream()
	{
		m_ostream.closeOstream();
	}


	/// CPPONLY say something about active states
	virtual string atRepr()
	{
		if (ISSETFLAG(m_flags, m_flagAtAllGen))
			return " at all generations";

		if (ISSETFLAG(m_flags, m_flagOnlyAtBegin) )
			return " at generation 0";

		if (ISSETFLAG(m_flags, m_flagOnlyAtEnd) )
			return " at ending generation";

		if (!m_atGen.empty() ) {
			string atStr = " at generation(s)";
			for (size_t i = 0; i < m_atGen.size(); ++i)
				atStr += " " + toStr(m_atGen[i]);
			return atStr;
		}

		string repr = " ";
		if (m_beginGen != 0)
			repr += "begin at " + toStr(m_beginGen) + " ";
		if (m_endGen != -1)
			repr += "end at " + toStr(m_endGen) + " ";
		if (m_stepGen != 1)
			repr += "at interval " + toStr(m_stepGen);
		return repr;
	}


	/// used by Python print function to print out the general information of the operator
	virtual string __repr__()
	{
		return "<simuPOP::operator> " ;
	}


	//@}

	/// CPPONLY determine if \c output=">". Used internally.
	bool noOutput()
	{
		return m_ostream.noOutput();
	}


	/// HIDDEN Initialize an operator against a population.
	virtual void initialize(const population & pop) {}

	/// CPPONLY
	subPopList applicableSubPops() const { return m_subPops; }

protected:
	/// analyze active generations: set m_flagAtAllGen etc
	void setFlags();

private:
	/// internal m_flags of the operator. They are set during initialization for
	/// performance considerations.
	static const size_t m_flagPreMating = PreMating;
	static const size_t m_flagDuringMating = DuringMating;
	static const size_t m_flagPostMating = PostMating;
	static const size_t m_flagAtAllGen = 8;
	static const size_t m_flagOnlyAtBegin = 16;
	static const size_t m_flagOnlyAtEnd = 32;
	// limited to haploid?
	static const size_t m_flagHaploid = 64;
	// limited to diploid?
	static const size_t m_flagDiploid = 128;

private:
	/// starting generation, default to 0
	int m_beginGen;

	/// ending generation, default to -1 (last genrration)
	int m_endGen;

	/// interval between active states, default to 1 (active at every generation)
	int m_stepGen;

	/// a list of generations that this oeprator will be active.
	/// typical usage is m_atGen=-1 to apply at the last generation.
	vectori m_atGen;

	/// various m_flags of Operator for faster processing.
	unsigned char m_flags;

	/// apply to all (-1) or one of the replicates.
	intList m_rep;

	/// apply to some of the (virtual) subpops.
	subPopList m_subPops;

	/// the output stream
	StreamProvider m_ostream;

	/// information fields that will be used by this operator
	vectorstr m_infoFields;

	/// last population to which this operator is applied.
	/// If the population is changed, the operator needs to be
	/// re-initialized.
	size_t m_lastPop;

};

typedef std::vector< baseOperator * > vectorop;

class opList
{
public:
	typedef vectorop::iterator iterator;
	typedef vectorop::const_iterator const_iterator;

public:
	opList(const vectorop & ops = vectorop());

	opList(const baseOperator & op);

	/// CPPONLY
	opList(const opList & rhs);

	baseOperator * operator[](size_t idx) const
	{
		DBG_ASSERT(idx < m_elems.size(), IndexError,
			"Operator index out of range");
		return m_elems[idx];
	}


	~opList();

	/// CPPONLY
	iterator begin()
	{
		return m_elems.begin();
	}


	/// CPPONLY
	iterator end()
	{
		return m_elems.end();
	}


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
class pause : public baseOperator
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
	pause(char stopOnKeyStroke = false, bool prompt = true,
		const stringFunc & output = ">", int stage = PostMating,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		baseOperator("", stage, begin, end, step, at, reps, subPops, infoFields),
		m_prompt(prompt), m_stopOnKeyStroke(stopOnKeyStroke)
	{
	}


	/// destructor
	virtual ~pause()
	{
	}


	/// HIDDEN
	virtual baseOperator * clone() const
	{
		return new pause(*this);
	}


	/// apply the \c pause operator to one population
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c pause operator
	virtual string __repr__()
	{
		return "<simuPOP::pause simulation>" ;
	}


private:
	bool m_prompt;

	char m_stopOnKeyStroke;

	static vectori s_cachedKeys;
};

/** This operator does nothing when it is applied to a population. It is
 *  usually used as a placeholder when an operator is needed syntactically.
 */
class noneOp : public baseOperator
{

public:
	/** Create a \c noneOp.
	 */
	noneOp(const stringFunc & output = ">",
		int stage = PostMating, int begin = 0, int end = 0, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr()) :
		baseOperator("", stage, begin, end, step, at, reps, subPops, infoFields)
	{
	}


	/// destructor
	virtual ~noneOp()
	{
	}


	/// HIDDEN
	virtual baseOperator * clone() const
	{
		return new noneOp(*this);
	}


	/// CPPONLY apply during mating, given \c pop, \c offspring, \c dad and \c mom
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
	                               individual * dad = NULL, individual * mom = NULL)
	{
		return true;
	}


	/// apply the \c noneOp operator to one population
	virtual bool apply(population & pop)
	{
		return true;
	}


	/// used by Python print function to print out the general information of the \c noneOp operator
	virtual string __repr__()
	{
		return "<simuPOP::None>" ;
	}


};

/** This operator accepts an expression that will be evaluated when this
 *  operator is applied. A list of if-operators will be applied when the
 *  expression returns \c True. Otherwise a list of else-operators will be
 *  applied.
 */
class ifElse : public baseOperator
{

public:
	/** Create a conditional operator that will apply operators \e ifOps if
	 *  condition \e cond is met and \e elseOps otherwise. The replicate and
	 *  generation applicability parameters (\e begin, \e end, \e step, \e at
	 *  and \e rep) of the \e ifOps and \e elseOps are ignored because their
	 *  applicability is determined by the \c ifElse operator.
	 */
	ifElse(const string & cond, const opList & ifOps = opList(), const opList & elseOps = opList(),
		const stringFunc & output = ">", int stage = PostMating,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		baseOperator("", stage, begin, end, step, at, reps, subPops, infoFields),
		m_cond(cond, ""), m_ifOps(ifOps), m_elseOps(elseOps)
	{
	};

	/// destructor
	virtual ~ifElse()
	{
	}


	/// HIDDEN
	virtual baseOperator * clone() const
	{
		return new ifElse(*this);
	}


	/// CPPONLY apply during mating, given \c pop, \c offspring, \c dad and \c mom
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// apply the \c ifElse operator to population \e pop.
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c ifElse operator
	virtual string __repr__()
	{
		return "<simuPOP::if else operator >";
	}


private:
	Expression m_cond;

	opList m_ifOps;
	opList m_elseOps;
};

/** This operator, when called, output the difference between current and the
 *  last called clock time. This can be used to estimate execution time of
 *  each generation. Similar information can also be obtained from
 *  <tt>turnOnDebug("DBG_PROFILE")</tt>, but this operator has the advantage
 *  of measuring the duration between several generations by setting \c step
 *  parameter.
 */
class ticToc : public baseOperator
{
public:
	/** Create a \c ticToc operator that outputs the elapsed since the last
	 *  time it was applied, and the overall time since it was created.
	 */
	ticToc(const stringFunc & output = ">",
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr()) :
		baseOperator(">", stage, begin, end, step, at, reps, subPops, infoFields)
	{
		time(&m_startTime);
		m_lastTime = m_startTime;
	};

	/// destructor
	virtual ~ticToc()
	{
	};

	/// HIDDEN
	virtual baseOperator * clone() const
	{
		return new ticToc(*this);
	}


	/// HIDDEN
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c ticToc operator
	virtual string __repr__()
	{
		return "<simuPOP::tic toc performance monitor>" ;
	}


private:
	time_t m_startTime, m_lastTime;
};

/** This operator sets the number of ancestral generations to keep during the
 *  evolution of a population. This is usually used to start storing ancestral
 *  generations at the end of an evolutionary process. A typical usage is
 *  <tt>setAncestralDepth(1, at=-1)</tt> which will cause the parental
 *  generation of the present population to be stored at the last generation
 *  of an evolutionary process.
 */
class setAncestralDepth : public baseOperator
{

public:
	/** Create a \c setAncestralDepth operator that sets the ancestral depth of
	 *  an population. It basically calls the
	 *  <tt>population.setAncestralDepth</tt> member function of a population.
	 */
	setAncestralDepth(int depth, const stringFunc & output = ">",
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		baseOperator(">", stage, begin, end, step, at, reps, subPops, infoFields),
		m_depth(depth)
	{
	};

	/// destructor
	virtual ~setAncestralDepth()
	{
	};

	/// HIDDEN
	virtual baseOperator * clone() const
	{
		return new setAncestralDepth(*this);
	}


	/// apply the \c setAncestralDepth operator to population \e pop.
	virtual bool apply(population & pop)
	{
		pop.setAncestralDepth(m_depth);
		return true;
	}


	/// used by Python print function to print out the general information of the \c setAncestralDepth operator
	virtual string __repr__()
	{
		return "<simuPOP::setAncestralDepth>";
	}


private:
	int m_depth;
};

/** Turn on debug. There are several ways to turn on debug information for
 *  non-optimized modules, namely
 *  \li set environment variable \c SIMUDEBUG.
 *  \li use <tt>simuOpt.setOptions(debug)</tt> function.
 *  \li use function \c TurnOnDebug
 *  \li use the \c turnOnDebug operator
 *
 *  The advantage of using an operator is that you can turn on debug at
 *  given generations.
 *  <funcForm>TurnOnDebug</funcForm>
 */
class turnOnDebug : public baseOperator
{
public:
	/** create a \c turnOnDebug operator that turns on debug information \e code
	 *  when it is applied to a population.
	 */
	turnOnDebug(const string & code,
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr()) :
		baseOperator(">", stage, begin, end, step, at, reps, subPops, infoFields),
		m_code(code)
	{
	};

	/// destructor
	virtual ~turnOnDebug()
	{
	};

	/// HIDDEN
	virtual baseOperator * clone() const
	{
		return new turnOnDebug(*this);
	}


	/// HIDDEN
	virtual bool apply(population & pop)
	{
		TurnOnDebug(m_code);
		return true;
	}


	/// used by Python print function to print out the general information of the \c turnOnDebug operator
	virtual string __repr__()
	{
		return "<simuPOP::turnOnDebug>";
	}


private:
	string m_code;
};


/** Turn off certain debug information. Please refer to operator \c turnOnDebug
 *  for detailed usages.
 *  <funcForm>TurnOffDebug</funcForm>
 */
class turnOffDebug : public baseOperator
{

public:
	/** create a \c turnOffDebug operator that turns off debug information
	 *  \e code when it is applied to a population.
	 */
	turnOffDebug(const string & code = "DBG_ALL",
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr()) :
		baseOperator(">", stage, begin, end, step, at, reps, subPops, infoFields),
		m_code(code)
	{
	};

	/// destructor
	virtual ~turnOffDebug()
	{
	};

	/// HIDDEN
	virtual baseOperator * clone() const
	{
		return new turnOffDebug(*this);
	}


	/// HIDDEN
	virtual bool apply(population & pop)
	{
		TurnOffDebug(m_code);
		return true;
	}


	/// used by Python print function to print out the general information of the \c turnOffDebug operator
	virtual string __repr__()
	{
		return "<simuPOP::turnOffDebug>";
	}


private:
	string m_code;
};


/** An operator that calls a user-defined function when it is applied to a
 *  population (pre- or post-mating) or offsprings (during-mating). The
 *  function can have have parameters \c pop when the operator is applied
 *  pre- or post-mating, <tt>pop, off, dad, mom</tt> when the operator is
 *  applied during-mating. An optional parameter can be passed if parameter
 *  \e param is given. In the during-mating case, parameters \c pop,
 *  \c dad and \c mom can be ignored if \e offspringOnly is set to \c True.
 */
class pyOperator : public baseOperator
{
public:
	/** Create a pure-Python operator that calls a user-defined function when
	 *  it is applied. Depending on parameters \e stage, \e param, and
	 *  \e offspringOnly, the function should have one of the following forms:
	 *  \li <tt>func(pop)</tt> if <tt>stage=PreMating</tt> or \c PostMating,
	 *    and without \c param.
	 *  \li <tt>func(pop, param)</tt> if <tt>stage=PreMating</tt> or
	 *    \c PostMating, and with \c param.
	 *  \li <tt>func(pop, off, dad, mom)</tt> if <tt>stage=DuringMating</tt> and
	 *    without \c param.
	 *  \li <tt>func(pop, off, dad, mom, param)</tt> if <tt>stage=DuringMating</tt>,
	 *    and with \c param.
	 *  \li <tt>func(off)</tt> if <tt>stage=DuringMating</tt>,
	 *    <tt>offspringOnly=True</tt> and without \c param.
	 *  \li <tt>func(off, param)</tt> if <tt>stage=DuringMating</tt>,
	 *    <tt>offspringOnly=True</tt> and with \c param.
	 *
	 *  where \c pop is the population to which the operator is applied, \c off
	 *  is the offspring of \c dad and \c mom, and \c param is the parameter
	 *  \e param specified when the operator is created. When this operator is
	 *  applied during mating, it can be used in the \c ops parameter of a
	 *  mating scheme, or used in the \c ops parameter of \c simulator.evolve
	 *  and be applied after an offspring has been created. Please refer to the
	 *  simuPOP user's guide for a detailed explanation.
	 *
	 *  This operator does not support parameters \e output, \e subPops and
	 *  \e infoFields. If certain output is needed, it should be handled in the
	 *  user defined function \e func. Because the status of files used by
	 *  other operators through parameter \e output is undetermined during
	 *  evolution, they should not be open or closed in this Python operator.
	 */
	pyOperator(PyObject * func, PyObject * param = NULL,
		int stage = PostMating, bool offspringOnly = false,
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr());

	/// HIDDEN
	virtual baseOperator * clone() const
	{
		return new pyOperator(*this);
	}


	/** Apply the \c pyOperator operator to population \e pop. Calling this
	 *  function is equivalent to call \e func with parameter \e pop and
	 *  optional parameter \e param.
	 */
	virtual bool apply(population & pop);

	/// CPPONLY
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// used by Python print function to print out the general information of the \c pyOperator operator
	virtual string __repr__()
	{
		return "<simuPOP::pyOperator>";
	}


private:
	/// the function
	pyFunc m_func;

	/// parammeters
	pyObject m_param;

	// whether or not pass pop, dad, mon in a duringMating py function.
	bool m_passOffspringOnly;
};


/** HIDDEN
 *  This function is used to test during mating operators. It simply apply
 *  operator \e op to \e dad, \e mom and \e off of population \e pop.
 *  If index of dad or mom is negative, NULL will be passed.
 */
void ApplyDuringMatingOperator(const baseOperator & op,
	population * pop, int dad, int mom, ULONG off);

}
#endif
