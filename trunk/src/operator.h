/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu
*                                                                         *
*   $LastChangedDate$
*   $Rev$                                                   *
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
using std::cout;
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
	typedef vectorvsp::const_iterator iterator;

public:
	subPopList(const vectorvsp & subPops = vectorvsp()) : m_subPops(subPops)
	{
	}


	bool empty() const
	{
		return m_subPops.empty();
	}


	size_t size() const
	{
		return m_subPops.size();
	}


	vspID operator[](unsigned int idx)
	{
		DBG_FAILIF(idx >= m_subPops.size(), IndexError,
			"Index out of range.");
		return m_subPops[idx];
	}


	void push_back(const vspID subPop)
	{
		m_subPops.push_back(subPop);
	}


	vectorvsp::const_iterator begin() const
	{
		return m_subPops.begin();
	}


	vectorvsp::const_iterator end() const
	{
		return m_subPops.end();
	}


private:
	vectorvsp m_subPops;
};


/** Operators are objects that act on populations. They can be applied to
 *  populations directly using their function forms, but they are usually
 *  managed and applied by a simulator. In the latter case, operators are
 *  passed to the \c evolve function of a simulator, and are applied repeatedly
 *  during the evolution of the simulator.
 *
 *  The \e baseOperator class is the base class for all operators. It defines
 *  a common user interface that specifies at which generations, at which stage
 *  of a life cycle, to which populations and subpopulation an operator will be
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
 *  terminal output (default), or a filename prefixed by one or more
 *  <tt>'>'</tt> characters. In the case of \c '>filename' (or equivalently
 *  \c 'filename'), the output from an operator is written to this file.
 *  However, if two operators write to the same file \c filename, or if an
 *  operator write to this file more than once, only the last write operation
 *  will succeed. In the case of <tt>'>>filename'</tt>, file \c filename
 *  will be opened at the beginning of the evolution and closed at the end.
 *  Outputs from multiple operators are appended. <tt>>>>filename</tt> works
 *  similar to <tt>>>filename</tt> but \c filename, if it already exists at the
 *  beginning of an evolutionary process, will not be cleared.
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
	 *    or \c 'filename' prefixed by one or more '>'.
	 *  \param stage Stage(s) of a life cycle at which an operator will be
	 *    applied. It can be \c PreMating, \c DuringMating, \c PostMating and
	 *    any of their combined stages \c PrePostMating, \c PreDuringMating
	 *    \c DuringPostMating and \c PreDuringPostMating. Note that all
	 *    operators have their default stage parameter and some of them ignores
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
	 *  \param rep A list of applicable replicates. An empty list (default) is
	 *    interpreted as all replicates in a simulator. Negative indexes such
	 *    as \c -1 (last replicate) is acceptable. <tt>rep=idx</tt> can be used
	 *    as a shortcut for <tt>rep=[idx]</tt>.
	 *  \param subPop A list of applicable (virtual) subpopulations, such as
	 *    <tt>subPop=[sp1, sp2, (sp2, vsp1)]</tt>. An empty list (default) is
	 *    interpreted as all subpopulations. <tt>subPop=[sp1]</tt> can be
	 *    simplied as <tt>subPop=sp1</tt>. Negative indexes are not supported.
	 *    Suport for this parameter vary from operator to operator. Some
	 *    operators do not support virtual subpopulations and some operators
	 *    do not support this parameter at all. Please refer to the reference
	 *    manual of individual operators for their support for this parameter.
	 *  \param infoFields A list of information fields that will be used by an
	 *    operator. You usually do not need to specify this parameter because
	 *    operators that use information fields usually have default values for
	 *    this parameter.
	 */
	baseOperator(string output, int stage, int begin, int end, int step, const intList & at,
		const repList & rep, const subPopList & subPop, const vectorstr & infoFields) :
		m_beginGen(begin), m_endGen(end), m_stepGen(step), m_atGen(at),
		m_flags(0), m_rep(rep), m_subPop(subPop),
		m_ostream(output), m_infoFields(infoFields),
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
	bool canApplyPreMating()
	{
		return ISSETFLAG(m_flags, m_flagPreMating);
	}


	/// set if this operator can be applied \em during-mating
	/// CPPONLY
	bool canApplyDuringMating()
	{
		return ISSETFLAG(m_flags, m_flagDuringMating);
	}


	/// set if this operator can be applied \em post-mating
	/// CPPONLY
	bool canApplyPostMating()
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


	/// CPPONLY if the operator will form genotype of offspring
	/**
	   If none of the during mating operator can form offspring, default will be used.
	 */
	bool formOffGenotype()
	{
		return ISSETFLAG(m_flags, m_flagFormOffGenotype);
	}


	/// CPPONLY
	void setFormOffGenotype(bool flag = true)
	{
		if (flag)
			SETFLAG(m_flags, m_flagFormOffGenotype);
		else
			RESETFLAG(m_flags, m_flagFormOffGenotype);
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


	virtual void initialize(const population & pop) {}

	/// CPPONLY
	subPopList applicableSubPops() const { return m_subPop; }

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
	static const size_t m_flagFormOffGenotype = 64;
	// limited to haploid?
	static const size_t m_flagHaploid = 128;
	// limited to diploid?
	static const size_t m_flagDiploid = 256;

private:
	/// starting generation, default to 0
	int m_beginGen;

	/// ending generation, default to -1 (last genrration)
	int m_endGen;

	/// interval between active states, default to 1 (active at every generation)
	int m_stepGen;

	/// a list of generations that this oeprator will be active.
	/// typical usage is m_atGen=-1 to apply at the last generation.
	intList m_atGen;

	/// various m_flags of Operator for faster processing.
	unsigned char m_flags;

	/// apply to all (-1) or one of the replicates.
	repList m_rep;

	/// apply to some of the (virtual) subpops.
	subPopList m_subPop;

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

/// pause a simulator
/**
   This operator pauses the evolution of a simulator at given generations or
   at a key stroke, using <tt>stopOnKeyStroke=True</tt> option. Users can
   use \c 'q' to stop an evolution. When a simulator is stopped, press any
   other key to resume	the simulation or escape to a Python shell to examine
   the status of the simulation by pressing \c 's'.

   There are two ways to use this operator, the first one is to pause the
   simulation at specified generations, using the usual operator parameters
   such as \c at. Another way is to pause a simulation with any key stroke,
   using the \c stopOnKeyStroke parameter. This feature is useful for a
   presentation or an interactive simulation. When \c 's' is pressed, this operator
   expose the current population to the main Python dictionary as variable \c pop
   and enter an interactive Python session. The way current population is exposed
   can be controlled by parameter \c exposePop and \c popName. This feature is
   useful when you want to examine the properties of a population during evolution.
 */
class pause : public baseOperator
{

public:
	/// stop a simulation. Press \c 'q' to exit or any other key to continue.
	/**
	   \param prompt if \c True (default), print prompt message.
	   \param stopOnKeyStroke if \c True, stop only when a key was pressed.
	   \param exposePop whether or not expose \c pop to user namespace, only
	    useful when user choose \c 's' at pause. Default to \c True.
	   \param popName by which name the population is exposed. Default to \c pop.
	 */
	pause(bool prompt = true, bool stopOnKeyStroke = false,
		bool exposePop = true, string popName = "pop",
		string output = ">", 
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = -1, const subPopList & subPop = subPopList(), const vectorstr & infoFields = vectorstr()) :
		baseOperator("", stage, begin, end, step, at, rep, subPop, infoFields),
		m_prompt(prompt), m_stopOnKeyStroke(stopOnKeyStroke),
		m_exposePop(exposePop), m_popName(popName)
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

	bool m_stopOnKeyStroke;

	/// whether or not expose population to user namespace
	bool m_exposePop;

	/// name of the population object
	string m_popName;
};

/// none operator
/**
   This operator does nothing.
 */
class noneOp : public baseOperator
{

public:
	/// create a none operator
	/**
	 */
	noneOp(string output = ">", 
		int stage = PostMating, int begin = 0, int end = 0, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(), const vectorstr & infoFields = vectorstr()) :
		baseOperator("", stage, begin, end, step, at, rep, subPop, infoFields)
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

/// conditional operator
/**
   This operator accepts
   \li an expression that will be evaluated when this operator is applied.
   \li an operator that will be applied if the expression is <tt>True</tt> (default to null).
   \li an operator that will be applied if the expression is <tt>False</tt> (default to null).

   When this operator is applied to a population, it will evaluate the expression and
   depending on its value, apply the supplied operator. Note that the \c begin, \c end,
   \c step, and \c at parameters of \c ifOp and \c elseOp will be ignored.
   For example, you can mimic the \c at parameter of an operator by
   <tt>ifElse('rep in [2,5,9]' operator)</tt>. The real use of this machanism is
   to monitor the population statistics and act accordingly.
 */
class ifElse : public baseOperator
{

public:
	/// create a conditional operator
	/**
	   \param cond expression that will be treated as a boolean variable
	   \param ifOp an operator that will be applied when \c cond is \c True
	   \param elseOp an operator that will be applied when \c cond is \c False

	 */
	ifElse(const string & cond, baseOperator * ifOp = NULL, baseOperator * elseOp = NULL,
		string output = ">", 
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(), const vectorstr & infoFields = vectorstr()) :
		baseOperator("", stage, begin, end, step, at, rep, subPop, infoFields),
		m_cond(cond, ""), m_ifOp(NULL), m_elseOp(NULL)
	{
		if (ifOp != NULL)
			m_ifOp = ifOp->clone();
		if (elseOp != NULL)
			m_elseOp = elseOp->clone();
	};

	/// destructor
	virtual ~ifElse()
	{
		if (m_ifOp != NULL)
			delete m_ifOp;
		if (m_elseOp != NULL)
			delete m_elseOp;
	};

	/// CPPONLY copy constructor
	ifElse(const ifElse & rhs)
		: baseOperator(rhs), m_cond(rhs.m_cond), m_ifOp(NULL), m_elseOp(NULL)
	{
		if (rhs.m_ifOp != NULL)
			m_ifOp = rhs.m_ifOp->clone();
		if (rhs.m_elseOp != NULL)
			m_elseOp = rhs.m_elseOp->clone();
	}


	/// HIDDEN
	virtual baseOperator * clone() const
	{
		return new ifElse(*this);
	}


	/// CPPONLY apply during mating, given \c pop, \c offspring, \c dad and \c mom
	virtual bool applyDuringMating(population & pop, RawIndIterator offspring,
		individual * dad = NULL, individual * mom = NULL);

	/// apply the \c ifElse operator to one population
	virtual bool apply(population & pop);

	/// used by Python print function to print out the general information of the \c ifElse operator
	virtual string __repr__()
	{
		return "<simuPOP::if else operator >" + m_ifOp->__repr__() ;
	}


private:
	Expression m_cond;

	baseOperator * m_ifOp;
	baseOperator * m_elseOp;
};

/// timer operator
/**
   This operator, when called, output the difference between current and the
   last called clock time. This can be used to estimate execution time of
   each generation. Similar information can also be obtained from
   <tt>turnOnDebug(DBG_PROFILE)</tt>, but this operator has the advantage
   of measuring the duration between several generations by setting \c step
   parameter.
   <funcForm>TicToc</funcForm>
 */
class ticToc : public baseOperator
{
public:
	/// create a timer
	ticToc(string output = ">", 
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(), const vectorstr & infoFields = vectorstr()) :
		baseOperator(">", stage, begin, end, step, at, rep, subPop, infoFields)
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

/// set ancestral depth
/**
   This operator set the number of ancestral generations to keep in a population.
   It is usually called like <tt>setAncestral(at=[-2])</tt> to start recording
   ancestral generations to a population at the end of the evolution. This is
   useful when constructing pedigree trees from a population.
 */
class setAncestralDepth : public baseOperator
{

public:
	/// create a \c setAncestralDepth operator
	setAncestralDepth(int depth, string output = ">", 
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(), const vectorstr & infoFields = vectorstr()) :
		baseOperator(">", stage, begin, end, step, at, rep, subPop, infoFields),
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


	/// apply the \c setAncestralDepth operator to one population
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

/// set debug on
/**
   Turn on debug. There are several ways to turn on debug information for
   non-optimized modules, namely
   \li set environment variable \c SIMUDEBUG.
   \li use <tt>simuOpt.setOptions(debug)</tt> function.
   \li use \c TurnOnDebug or \c TurnOnDebugByName function.
   \li use this \c turnOnDebug operator

   The advantage of using this operator is that you can turn on debug at
   given generations.
   <funcForm>TurnOnDebug</funcForm>
 */
class turnOnDebug : public baseOperator
{
public:
	/// create a \c turnOnDebug operator
	turnOnDebug(DBG_CODE code,
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(), const vectorstr & infoFields = vectorstr()) :
		baseOperator(">", stage, begin, end, step, at, rep, subPop, infoFields),
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
		return "<simuPOP::turnOnDebug " + dbgString(m_code) + ">";
	}


private:
	DBG_CODE m_code;
};

/// set debug off
/**
   Turn off debug.
   <funcForm>TurnOffDebug</funcForm>
 */
class turnOffDebug : public baseOperator
{

public:
	/// create a \c turnOffDebug operator
	turnOffDebug(DBG_CODE code,
		int stage = PreMating, int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(), const vectorstr & infoFields = vectorstr()) :
		baseOperator(">", stage, begin, end, step, at, rep, subPop, infoFields),
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


	/// apply the \c turnOffDebug operator to one population
	virtual bool apply(population & pop)
	{
		TurnOffDebug(m_code);
		return true;
	}


	/// used by Python print function to print out the general information of the \c turnOffDebug operator
	virtual string __repr__()
	{
		return "<simuPOP::turnOffDebug " + dbgString(m_code) + ">";
	}


private:
	DBG_CODE m_code;
};

/// A python operator that directly operate a population
/**
   This operator accepts a function that can take the form of
   \li <tt>func(pop)</tt> when <tt>stage=PreMating</tt> or \c PostMating, without setting \c param;
   \li <tt>func(pop, param)</tt> when <tt>stage=PreMating</tt> or \c PostMating, with \c param;
   \li <tt>func(pop, off, dad, mom)</tt> when <tt>stage=DuringMating</tt> and <tt>passOffspringOnly=False</tt>, without setting \c param;
   \li <tt>func(off)</tt> when <tt>stage=DuringMating</tt> and <tt>passOffspringOnly=True</tt>, and without setting \c param;
   \li <tt>func(pop, off, dad, mom, param)</tt> when <tt>stage=DuringMating</tt> and <tt>passOffspringOnly=False</tt>, with \c param;
   \li <tt>func(off, param)</tt> when <tt>stage=DuringMating</tt> and <tt>passOffspringOnly=True</tt>, with \c param.

   For \c Pre- and \c PostMating usages, a population and an optional parameter is passed to the given
   function. For \c DuringMating usages, population, offspring, its parents and an optional parameter
   are passed to the given function. Arbitrary operations can be applied to the population and
   offspring (if <tt>stage=DuringMating</tt>).

 */
class pyOperator : public baseOperator
{
public:
	/// Python operator, using a function that accepts a population object
	/**
	   \param func a Python function. Its form is determined by other parameters.
	   \param param any Python object that will be passed to \c func after \c pop parameter.
	    Multiple parameters can be passed as a tuple.
	   \param formOffGenotype This option tells the mating scheme this operator will set
	    the genotype of offspring (valid only for <tt>stage=DuringMating</tt>). By default
	    (<tt>formOffGenotype=False</tt>), a mating scheme will set the genotype of offspring before it is
	    passed to the given Python function. Otherwise, a 'blank' offspring will be passed.
	   \param passOffspringOnly if \c True, \c pyOperator will expect a function of form <tt>func(off [,param])</tt>,
	    instead of <tt>func(pop, off, dad, mom [, param])</tt> which is used when \c passOffspringOnly
	    is \c False. Because many during-mating \c pyOperator only need access to offspring,
	   this will improve efficiency. Default to \c False.

	   \note
	   \li Output to \c output is not supported. That is to say,
	   you have to open/close/append to files explicitly in the Python function.
	   Because files specified by \c output are controlled (opened/closed) by
	   simulators, they should not be manipulated in a \c pyOperator operator.
	   \li This operator can be applied \c Pre-, \c During- or <tt>Post- Mating</tt> and is applied \c PostMating
	   by default. For example, if you would like to examine the fitness values set by
	   a selector, a \c PreMating Python operator should be used.
	 */
	pyOperator(PyObject * func, PyObject * param = NULL,
		int stage = PostMating, bool formOffGenotype = false, bool offspringOnly = false,
		int begin = 0, int end = -1, int step = 1, const intList & at = intList(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(), const vectorstr & infoFields = vectorstr());

	/// HIDDEN
	virtual baseOperator * clone() const
	{
		return new pyOperator(*this);
	}


	/// apply the \c pyOperator operator to one population
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
