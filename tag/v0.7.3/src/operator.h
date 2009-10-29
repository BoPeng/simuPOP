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
#include "simupop_cfg.h"

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

namespace simuPOP
{

	/** \brief base class of all classes that manipulate populations.

	Operators are object that act on populations. They can be applied
	to populations directly using apply() member function, but most of
	the time they are managed and applied by a simulator.

	Operators can be applied at different stage(s) of a life cycle. More
	specifically, at pre-, duing- or post mating stage(s). Note that it
	is possible for an operator to apply multiple times in a life cycle.
	For example, an save-to-file operator might be applied before and
	after mating to trace parental information.

	Operators do not have to be applied at all generations. You can specify
	starting genertion, ending generation, gaps between applicable generations,
	or even specific generations to apply. For example, you might want to
	start applying migrations after certain heat-up generation; or you want
	to calculate every 10 generations.

	Operators can have outputs. Output can be standard output (terminal)
	or a file, which can be constant, or change with generation or replicate.
	Different operators can append to the same file to form table-like
	outputs.

	filename can have the following format:

	1. 'filename' this file will be closed after each use. I.e., if several
	operators output to the same file, only the last one will succeed.

	2. '>filename' the same as 'filaname'

	3. '>>filename' The file will be created at the beginning of evolution
	(simulator::evolve) and close at the end. Several operators can
	output to this file to form a table.

	4. '>>>filename' The same as '>>filename' except that the file will not
	be cleared at the beginning of evolution if it is not empty.

	5. '>' out put to standard output.

	6. '' supress output.

	Most operators are applied to every replicate of a simulator during
	evolution. However, you can apply opertors to one or a group of
	replicates only. For example, you can initialize different replicates
	with different initial values and then start evolution. c.f.
	simulator::setGroup .

	Please refer to help(baseOperator) and help(baseOperator.__init__) for
	detailed information about member functions and parameters.

	@author Bo Peng
	*/

	class Operator
	{
		public:

			/** @name constructor and destructor */
			//@{

			/// create an operator (this function is not supposed to be called directly)
			/**
			\param output a string of output filename. Different operators will have
			   different default output (most commonly '>' or '')
			\param outputExpr an expression that determines output filename dynamically.
			\param begin start generation. default to 1. negative number is interpreted
			as endGeneration + begin
			\param end stop applying after this generation. negative number is allowed
			\param step number of generations between active generations. default to 1
			\param at an array of active generations. If given, stage, begin, end,
			step will be ignored.
			\param rep applicable replicate. It can be replicate number 0 ~ (number of replicate -1),
			REP_ALL (all replicates) or REP_LAST (only to the last replicate).
			Usually default to REP_ALL.
			\param grp applicable group, default to GRP_ALL.  A group number for each
			replicate is set by simulator.__init__ or simulator::setGroup(). grp, if not
			GRP_ALL, will be compared to the group number of this replicate before applying.

			DEVONLY{ DO NOT SET DEFAULT PARAMETE. This will force all
			derived classes to pay attention to parameter number. }
			*/
			Operator(string output, string outputExpr, int stage,
				int begin, int end, int step, vectorl at,
				int rep, int grp, const vectorstr& infoFields):
			m_beginGen(begin), m_endGen(end), m_stepGen(step), m_atGen(at),
				m_flags(0), m_rep(rep), m_grp(grp),
				m_ostream(output, outputExpr), m_infoFields(infoFields)
			{
				DBG_FAILIF(step<=0, ValueError, "step need to be at least one");

				setApplicableStage(stage);
				setFlags();
			}

			/// destroy an operator
			virtual ~Operator()
			{
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new Operator(*this);
			}

			//@}

			/** @name  applicable generations (also judge from rep and group.
			  use of parameter start, end, every, at, group, rep
			*/
			//@{

			/// judge if this operator is active
			/** Judge if this operator is active under the conditions like current
			 replicate, group number of current replicate, current generation.
			ending generation etc.

			\note will be called by simulator before applying.
			*/
			bool isActive(UINT rep, UINT numRep, long gen, long end, int grp, bool repOnly =false);

			/// return applicable group
			int applicableGroup()
			{
				return m_grp;
			}

			/// set applicable group.
			/** GRP_ALL is the default value (applicable to
			all groups. )
				  Otherwise, the operator is applicable to ONE group of replicates.
				  groups can be set in simulator::setGroup()
			*/
			void setApplicableGroup(int grp=GRP_ALL)
			{
				m_grp = grp;
			}

			/// return applicable replicate
			int applicableReplicate()
			{
				return m_rep;
			}

			/// set applicable replicate.
			void setApplicableReplicate(int rep)
			{
				m_rep = rep;
			}

			/// set applicable generation parrameters: stage, begin, end, step and at
			void setActiveGenerations(int begin=0, int end=-1, int step=1, vectorl at=vectorl())
			{
				m_beginGen = begin;
				m_endGen = end;
				m_stepGen = step;
				m_atGen = at;

				DBG_FAILIF(step<=0, ValueError, "step need to be at least one");

				/// set certain m_flags to speed up using this machanism.
				setFlags();
			}

			//@}

			/** @name applicable stages
					pre, during, post-mating methods */
			//@{

			/// set m_stage settings. This is usually not usable since
			/// the m_stage setting are set by default for each Operator.
			void setApplicableStage(int stage)
			{
				RESETFLAG(m_flags, PreDuringPostMating);
				SETFLAG(m_flags, stage);
			}

			/// Can this operator be applied pre-mating?
			bool canApplyPreMating()
			{
				return ISSETFLAG(m_flags, m_flagPreMating);
			}

			/// Can this operator be applied uring mating?
			bool canApplyDuringMating()
			{
				return ISSETFLAG(m_flags, m_flagDuringMating);
			}

			/// Can this operator be applied post-mating?
			bool canApplyPostMating()
			{
				return ISSETFLAG(m_flags, m_flagPostMating);
			}

			/// can be applied pre or post mating.
			bool canApplyPreOrPostMating()
			{
				return ISSETFLAG(m_flags, m_flagPreMating) || ISSETFLAG(m_flags, m_flagPostMating);
			}

			/// get the number of information fields for this operator
			UINT infoSize()
			{
				return m_infoFields.size();
			}
			
			/// get the information field specified by user (or by default)
			string infoField(UINT idx)
			{
				DBG_ASSERT(idx < m_infoFields.size(), IndexError, "Given info index " + toStr(idx) + 
					" is out of range of 0 ~ " + toStr(m_infoFields.size()));
				return m_infoFields[idx];
			}
			
			/// CPPONLY
			/// if the operator will form genotype of offspring.
			/// if none of the during mating operator can form offspring, default will be used
			bool formOffGenotype()
			{
				return ISSETFLAG(m_flags, m_flagFormOffGenotype);
			}

			/// CPPONLY
			void setFormOffGenotype(bool flag=true)
			{
				if(flag)
					SETFLAG(m_flags, m_flagFormOffGenotype);
				else
					RESETFLAG(m_flags, m_flagFormOffGenotype);
			}

			///  providing interface to apply operator before during or after
			///  mating.
			virtual bool applyWithScratch(population& pop, population& scratch, int stage)
			{
				return apply(pop);
			}

			/// apply to one population, does not check if the oeprator is activated.
			virtual bool apply(population& pop)
			{
				throw SystemError("This function Operator::apply() is not supposed to be called.");
			}

			/// give pop, offspring, pop and mom.
			virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
				individual* dad=NULL, individual* mom=NULL)
			{
				throw SystemError("ApplyDuringMating: This function is not supposed to be called.");
			}

			//@}
			/** @name dealing with output
			  separator, persistant files, $gen etc substitution.
			*/
			//@{

			/// ostream, if not set during construction.
			void setOutput(string output="", string outputExpr="")
			{
				m_ostream.setOutput(output, outputExpr);
			}

			/// get output stream. This function is not exposed to user.
			ostream& getOstream( PyObject* dict=NULL, bool readable=false)
			{
				return m_ostream.getOstream( dict, readable);
			}

			/// close ostream and delete ostream pointer.. if it is a ofstream.
			void closeOstream()
			{
				m_ostream.closeOstream();
			}

			/// CPPONLY, say something about active states
			virtual string atRepr()
			{
				if( ISSETFLAG(m_flags, m_flagAtAllGen))
					return(" at all generations");

				if( ISSETFLAG( m_flags, m_flagOnlyAtBegin) )
					return(" at generation 0");

				if( ISSETFLAG(m_flags, m_flagOnlyAtEnd) )
					return(" at ending generation");

				if( !m_atGen.empty() )
				{
					string atStr = " at generation(s)";
					for(size_t i=0; i<m_atGen.size(); ++i)
						atStr += " " + toStr(m_atGen[i]);
					return atStr;
				}

				string repr = " ";
				if( m_beginGen != 0 )
					repr += "begin at " + toStr(m_beginGen) + " ";
				if( m_endGen != -1 )
					repr += "end at " + toStr(m_endGen) + " ";
				if( m_stepGen != 1)
					repr += "at interval " + toStr(m_stepGen);
				return repr;
			}

			virtual string __repr__()
			{
				return "<simuPOP::operator> " ;
			}

			//@}

			/// if output=">". used internally
			bool noOutput(){ return m_ostream.noOutput(); }

		protected:

			/// analyze active generations: set m_flagAtAllGen etc
			void setFlags();

		private:

			/// internal m_flags of the operator. They are set during initialization for
			/// performance considerations.
			static const size_t m_flagPreMating      = PreMating;
			static const size_t m_flagDuringMating   = DuringMating;
			static const size_t m_flagPostMating     = PostMating;
			static const size_t m_flagAtAllGen       = 8;
			static const size_t m_flagOnlyAtBegin    = 16;
			static const size_t m_flagOnlyAtEnd      = 32;
			static const size_t m_flagFormOffGenotype= 64;

		private:

			/// starting generation, default to 0
			int m_beginGen;

			/// ending generation, default to -1 (last genrration)
			int m_endGen;

			/// interval between active states, default to 1 (active at every generation)
			int m_stepGen;

			/// a list of generations that this oeprator will be active.
			/// typical usage is m_atGen=[-1] to apply at the last generation.
			vectorl m_atGen;

			/// various m_flags of Operator for faster processing.
			unsigned char m_flags;

			/// apply to all (-1) or one of the replicates.
			int m_rep;

			/// apply to all (-1) or one of the groups of simulatior replicates.
			int m_grp;

			/// the output stream
			StreamProvider m_ostream;

			/// information fields that will be used by this operator
			vectorstr m_infoFields;
	};

	/* 
	pause simulation, press any key to stop
	*/

	class pause: public Operator
	{

		public:
			/** \brief
			stop simulation. press q to exit and any other key to continue

			\param prompt if true (default), print prompt message.
			\param stopOnKeyStroke if true, goon if no key was pressed
			\param exposePop whether or not expose pop to user namespace, only
			  useful when user choose 's' at pause. Default to true.
			\param popName by which name the population is exposed? default to 'pop'

			*/
			pause(bool prompt=true, bool stopOnKeyStroke=false,
				bool exposePop=true, string popName="pop",
				string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_LAST, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			Operator("", "", stage, begin, end, step, at, rep, grp, infoFields),
				m_prompt(prompt), m_stopOnKeyStroke(stopOnKeyStroke),
				m_exposePop(exposePop), m_popName(popName)
			{
			}

			/// destructor
			virtual ~pause(){}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new pause(*this);
			}

			/// simply output some info
			virtual bool apply(population& pop)
			{
				char a;

				if(m_stopOnKeyStroke)
				{
					// check if key is already pressed
					if( !simuPOP_kbhit() )
						return true;
					else
						simuPOP_getch();
				}
				// clear input and wait for user input
				// std::cin.clear();
				// std::cin.ignore(std::numeric_limits<int>::max());

				if(m_prompt)
				{
					cout << "Simulation paused. " << endl
						<< " Press " << endl
						<< "   q to stop evolution, " << endl
						<< "   s to start an interative shell, (current population is ";
					if(m_exposePop)
						cout << "exported as " << m_popName << ')' << endl;
					else
						cout << "not exported)" << endl;
					cout << "   or any other key to continue...." << endl;
				}
				a = simuPOP_getch();			  // std::cin.get(a);

				if( a == 'q' || a=='Q' )
					throw SystemError("Terminated by user");
				else if( a == 's' || a == 'S' )
				{
					// export current population
					PyObject* popObj;
					if(m_exposePop)
					{
						popObj = pyPopObj(static_cast<void*>(&pop));
						if( popObj == NULL)
							throw SystemError("Could not expose population pointer. Compiled with the wrong version of SWIG? ");

						// get global dictionary
						mainVars().setVar(m_popName, popObj);
					}
					PyRun_InteractiveLoop(stdin, "<stdin>");
					// if expose pop, release it.
					if(m_exposePop)
						mainVars().removeVar(m_popName);
				}

				// clear input and wait for user input
				// std::cin.clear();
				// std::cin.ignore(std::numeric_limits<int>::max());

				return true;
			}

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

	/* 
	  None operator, do nothing.
	  */

	class noneOp: public Operator
	{

		public:
			/** \brief do nothing
			 */
			noneOp( string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=0, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			Operator("", "", stage, begin, end, step, at, rep, grp, infoFields)
			{
			}

			/// destructor
			virtual ~noneOp()
			{
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new noneOp(*this);
			}

			/// simply output some info
			///  providing interface to apply operator before during or after
			///  mating.
			virtual bool applyWithScratch(population& pop, population& scratch, int stage)
			{
				return true;
			}

			/// give pop, offspring, pop and mom.
			virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
				individual* dad=NULL, individual* mom=NULL)
			{
				return true;
			}

			virtual bool apply(population& pop)
			{
				return true;
			}

			virtual string __repr__()
			{
				return "<simuPOP::None>" ;
			}
	};

	/* 
	 evaluate an expression, call if_op or else_op
	 */

	class ifElse: public Operator
	{

		public:
			/** \brief conditional operator
			\param cond expression, will be treated as bool variable.
			\param ifOp if operator, be called when expr is true
			\param elseOp else operator, be called when expr is false
			*/
			ifElse(const string& cond, Operator* ifOp=NULL, Operator* elseOp = NULL,
				string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			Operator("", "", stage, begin, end, step, at, rep, grp, infoFields),
				m_cond(cond,""), m_ifOp(NULL), m_elseOp(NULL)
			{
				if( ifOp != NULL)
					m_ifOp = ifOp->clone();
				if( elseOp != NULL)
					m_elseOp = elseOp->clone();
			};

			/// destructor
			virtual ~ifElse()
			{
				if( m_ifOp != NULL)
					delete m_ifOp;
				if( m_elseOp != NULL)
					delete m_elseOp;
			};

			/// copy constructor
			/// CPPONLY
			ifElse(const ifElse& rhs)
				: Operator(rhs), m_cond(rhs.m_cond), m_ifOp(NULL), m_elseOp(NULL)
			{
				if( rhs.m_ifOp != NULL)
					m_ifOp = rhs.m_ifOp->clone();
				if( rhs.m_elseOp != NULL)
					m_elseOp = rhs.m_elseOp->clone();
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new ifElse(*this);
			}

			/// simply output some info
			///  providing interface to apply operator before during or after
			///  mating.
			virtual bool applyWithScratch(population& pop, population& scratch, int stage)
			{
				m_cond.setLocalDict(pop.dict());
				bool res = m_cond.valueAsBool();

				if(res && m_ifOp != NULL)
					return m_ifOp->applyWithScratch(pop, scratch, stage);
				else if( ! res && m_elseOp != NULL)
					return m_elseOp->applyWithScratch(pop, scratch, stage);
				return true;
			}

			/// give pop, offspring, pop and mom.
			virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
				individual* dad=NULL, individual* mom=NULL)
			{
				m_cond.setLocalDict(pop.dict());
				bool res = m_cond.valueAsBool();

				if(res && m_ifOp != NULL)
					return m_ifOp->applyDuringMating(pop, offspring, dad, mom);
				else if(!res && m_elseOp != NULL)
					return m_elseOp->applyDuringMating(pop, offspring, dad, mom);
				return true;
			}

			virtual bool apply(population& pop)
			{
				m_cond.setLocalDict(pop.dict());
				bool res = m_cond.valueAsBool();

				if(res && m_ifOp != NULL)
					return m_ifOp->apply(pop);
				else if( !res && m_elseOp != NULL)
					return m_elseOp->apply(pop);
				return true;
			}

			virtual string __repr__()
			{
				return "<simuPOP::if else operator >" + m_ifOp->__repr__() ;
			}

		private:
			Expression m_cond;

			Operator* m_ifOp;
			Operator *m_elseOp;
	};

	/* 
	  evaluate an expression, call if_op or else_op
	  */

	class ticToc: public Operator
	{

		public:
			/** \brief timer
			if called, output time passed since last calling time.
			*/
			ticToc( string output=">", string outputExpr="",
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			Operator(">", "", stage, begin, end, step, at, rep, grp, infoFields)
			{
				time(&m_startTime);
				m_lastTime = m_startTime;
			};

			/// destructor
			virtual ~ticToc()
			{
			};

			/// this function is very important
			virtual Operator* clone() const
			{
				return new ticToc(*this);
			}

			virtual bool apply(population& pop)
			{
				time_t tmpTime;

				// this may not be correct but wrap is a possible problem.
				if( ! this->noOutput() )
				{
					ostream & out = this->getOstream(pop.dict());
					// since last time
					double timeDiff = difftime( time(&tmpTime), m_lastTime);
					out << "Elapsed Time: " << int(timeDiff*100)/100. ;
					// since beginning
					timeDiff = difftime( time(&tmpTime), m_startTime);
					int h = int(timeDiff/3600);
					int m = int((timeDiff - h*3600)/60);
					int s = int(timeDiff - h*3600 - m*60);
					out << "s  Overall Time: " << std::setw(2) << std::setfill('0') << h
						<< ":" << std::setw(2) << std::setfill('0') << m << ":" << std::setw(2)
						<< std::setfill('0') << s << endl;
					this->closeOstream();
				}

				time(&m_lastTime);
				return true;
			}

			virtual string __repr__()
			{
				return "<simuPOP::tic toc performance monitor>" ;
			}

		private:

			time_t m_startTime, m_lastTime;
	};

	/* 
	   set ancestral depth
	   */

	class setAncestralDepth: public Operator
	{

		public:
			/** \brief timer
			if called, output time passed since last calling time.
			*/
			setAncestralDepth( int depth, string output=">", string outputExpr="",
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			Operator(">", "", stage, begin, end, step, at, rep, grp, infoFields),
				m_depth(depth)
			{
			};

			/// destructor
			virtual ~setAncestralDepth()
			{
			};

			/// this function is very important
			virtual Operator* clone() const
			{
				return new setAncestralDepth(*this);
			}

			virtual bool apply(population& pop)
			{
				pop.setAncestralDepth(m_depth);
				return true;
			}

			virtual string __repr__()
			{
				return "<simuPOP::setAncestralDepth>";
			}

		private:
			int m_depth;
	};

	/* 
	   set debug
	   */

	class turnOnDebugOp: public Operator
	{

		public:
			/** \brief turn on debug

			 */
			turnOnDebugOp(DBG_CODE code,
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			Operator(">", "", stage, begin, end, step, at, rep, grp, infoFields),
				m_code(code)
			{
			};

			/// destructor
			virtual ~turnOnDebugOp()
			{
			};

			/// this function is very important
			virtual Operator* clone() const
			{
				return new turnOnDebugOp(*this);
			}

			virtual bool apply(population& pop)
			{
				TurnOnDebug(m_code);
				return true;
			}

			virtual string __repr__()
			{
				return "<simuPOP::turnOnDebug " + dbgString(m_code) + ">";
			}

		private:
			DBG_CODE m_code;
	};

	/* 
	   set debug
	   */

	class turnOffDebugOp: public Operator
	{

		public:
			/** \brief turn on debug

			 */
			turnOffDebugOp(DBG_CODE code,
				int stage=PreMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			Operator(">", "", stage, begin, end, step, at, rep, grp, infoFields),
				m_code(code)
			{
			};

			/// destructor
			virtual ~turnOffDebugOp()
			{
			};

			/// this function is very important
			virtual Operator* clone() const
			{
				return new turnOffDebugOp(*this);
			}

			virtual bool apply(population& pop)
			{
				TurnOffDebug(m_code);
				return true;
			}

			virtual string __repr__()
			{
				return "<simuPOP::turnOffDebug " + dbgString(m_code) + ">";
			}

		private:
			DBG_CODE m_code;
	};

	/* 
	   pyOperator, the ultimate python operator
	   */

	class pyOperator: public Operator
	{

		public:
			/** \brief python operator, using a function that accept a population object

			\param func a python function that accept a population and perform whatever
			operation it wants to.
			\param para any python object that will be passed to func after pop parameter.
			  Multiple parameter can be passed as a tuple.
			\param formOffGenotype if stage=DuringMating, set this parameter to false
			  will disallow random mating to set genotype.
			\param passOffspringOnly Default to false. If true, pyOperator will expect
			  a function of form func(off, param), instead of func(pop, off, dad, mon, param)
			  when passOffspringOnly is false. Since many duringMating pyOperator only need
			access to offspring, this will imporve efficiency.

			Note: (FIXME) output to output or outputExpr is not yet supported. Ideally,
			this func will take two parameters with pop and then a filehandle to output,
			however, differentiating output, append etc is too troublesome right now.
			*/
			pyOperator(PyObject* func, PyObject* param=NULL,
				int stage=PostMating, bool formOffGenotype=false, bool passOffspringOnly=false,
				int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			Operator(">", "", stage, begin, end, step, at, rep, grp, infoFields),
				m_func(func), m_param(param), m_passOffspringOnly(passOffspringOnly)
			{
				if( !PyCallable_Check(func))
					throw ValueError("Passed variable is not a callable python function.");

				// inc reference count
				Py_XINCREF(func);

				if(param != NULL)
					Py_XINCREF(param);

				this->setFormOffGenotype(formOffGenotype);
			};

			/// destructor
			virtual ~pyOperator()
			{
				Py_DECREF(m_func);

				if( m_param != NULL)
					Py_DECREF(m_param);
			};

			/// CPPONLY need a copy operator because of m_func
			pyOperator(const pyOperator& rhs):
			Operator(rhs),
				m_func(rhs.m_func),
				m_param(rhs.m_param)
			{
				Py_INCREF( m_func);

				if( m_param != NULL)
					Py_INCREF( m_param);
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new pyOperator(*this);
			}

			virtual bool apply(population& pop)
			{
				// call the python function, pass the whole population in it.
				// get pop object
				PyObject* popObj = pyPopObj(static_cast<void*>(&pop));
				// if pop is valid?
				if(popObj == NULL)
					throw SystemError("Could not pass population to the provided function. \n"
						"Compiled with the wrong version of SWIG?");

				// parammeter list, ref count increased
				bool resBool;
				// parenthesis is needed since PyCallFuncX are macros.
				if( m_param == NULL)
				{
					PyCallFunc(m_func, "(O)", popObj, resBool, PyObj_As_Bool);
				}
				else
				{
					PyCallFunc2(m_func, "(OO)", popObj, m_param, resBool, PyObj_As_Bool);
				}

				Py_DECREF(popObj);
				return resBool;
			}

			virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
				individual* dad=NULL, individual* mom=NULL)
			{
				// get offspring object
				PyObject* offObj = pyIndObj(static_cast<void*>(&(*offspring)));
				DBG_FAILIF(offObj == NULL, SystemError,
					"Could not pass offspring to the provided function. \n"
					"Compiled with the wrong version of SWIG?");

				PyObject* arglist, *result;
				if( m_passOffspringOnly )
				{
					// parammeter list, ref count increased
					if( m_param == NULL)
						arglist = Py_BuildValue("(O)", offObj);
					else
						arglist = Py_BuildValue("(OO)", offObj, m_param);

					// we do not need to catch exceptions here,
					// our wrapper will do that
					result = PyEval_CallObject(m_func, arglist);
					Py_DECREF(offObj);
				}
				else
				{
					// call the python function, pass all the parameters to it.
					// get pop object
					PyObject* popObj = pyPopObj(static_cast<void*>(&pop));

					// get dad object
					PyObject* dadObj, *momObj;
					if(dad == NULL)
					{
						Py_INCREF(Py_None);
						dadObj = Py_None;
					}
					else
						dadObj = pyIndObj(static_cast<void*>(dad));

					if(mom == NULL)
					{
						Py_INCREF(Py_None);
						momObj = Py_None;
					}
					else
						momObj = pyIndObj(static_cast<void*>(mom));

					// if pop is valid?
					DBG_FAILIF(popObj == NULL || dadObj == NULL || momObj == NULL, SystemError,
						"Could not pass population or parental individuals to the provided function. \n"
						"Compiled with the wrong version of SWIG?");

					// parammeter list, ref count increased
					if( m_param == NULL)
						arglist = Py_BuildValue("(OOOO)", popObj, offObj, dadObj, momObj);
					else
						arglist = Py_BuildValue("(OOOOO)", popObj, offObj, dadObj, momObj, m_param);

					// we do not need to catch exceptions here,
					// our wrapper will do that
					result = PyEval_CallObject(m_func, arglist);

					Py_DECREF(offObj);
					Py_DECREF(popObj);
					Py_DECREF(dadObj);
					Py_DECREF(momObj);
				}
				// release arglist
				Py_DECREF(arglist);

				if( result == NULL)
				{
					PyErr_Print();
					throw ValueError("Invalid return from provided function. (Be sure to return True or False)");
				}
				// result should be a boolean value
				bool resBool;
				// defined in utility.h
				PyObj_As_Bool(result, resBool);
				Py_DECREF(result);
				return resBool;
			}

			virtual string __repr__()
			{
				return "<simuPOP::pyOperator>";
			}

		private:

			/// the function
			PyObject * m_func;

			/// parammeters
			PyObject * m_param;

			// whether or not pass pop, dad, mon in a duringMating py function.
			bool m_passOffspringOnly;
	};


	/* 
	   infoOperator that manipulate individuals
    */
	class pyIndOperator: public Operator
	{

		public:
			/** \brief individual operator that apply a function to each individual

			\param func a python function that accept a individual and perform whatever
			operation it wants to.
			\param para any python object that will be passed to func after pop parameter.
			  Multiple parameter can be passed as a tuple.
			\param formOffGenotype if stage=DuringMating, set this parameter to false
			  will disallow random mating to set genotype.

			Note: (FIXME) output to output or outputExpr is not yet supported. Ideally,
			this func will take two parameters with pop and then a filehandle to output,
			however, differentiating output, append etc is too troublesome right now.
			*/
			pyIndOperator(PyObject* func, PyObject* param=NULL,
				int stage=PostMating, bool formOffGenotype=false, 
				int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			Operator(">", "", stage, begin, end, step, at, rep, grp, infoFields),
				m_func(func), m_param(param)
			{
				if( !PyCallable_Check(func))
					throw ValueError("Passed variable is not a callable python function.");

				// inc reference count
				Py_XINCREF(func);

				if(param != NULL)
					Py_XINCREF(param);

				this->setFormOffGenotype(formOffGenotype);
			};

			/// destructor
			virtual ~pyIndOperator()
			{
				Py_DECREF(m_func);

				if( m_param != NULL)
					Py_DECREF(m_param);
			};

			/// CPPONLY need a copy operator because of m_func
			pyIndOperator(const pyIndOperator& rhs):
			Operator(rhs),
				m_func(rhs.m_func),
				m_param(rhs.m_param)
			{
				Py_INCREF( m_func);

				if( m_param != NULL)
					Py_INCREF( m_param);
			}

			/// this function is very important
			virtual Operator* clone() const
			{
				return new pyIndOperator(*this);
			}

			virtual bool apply(population& pop)
			{
				// call the python function, pass the each individual to it.
				// get pop object
				for(size_t i=0; i<pop.popSize(); ++i) 
				{
					PyObject* indObj = pyIndObj(static_cast<void*>(&pop.ind(i)));
					// if pop is valid?
					if(indObj == NULL)
						throw SystemError("Could not pass population to the provided function. \n"
							"Compiled with the wrong version of SWIG?");

					// parammeter list, ref count increased
					bool resBool;
					// parenthesis is needed since PyCallFuncX are macros.
					if( m_param == NULL)
					{
						PyCallFunc(m_func, "(O)", indObj, resBool, PyObj_As_Bool);
					}
					else
					{
						PyCallFunc2(m_func, "(OO)", indObj, m_param, resBool, PyObj_As_Bool);
					}

					Py_DECREF(indObj);

					if(!resBool)
						return false;
				}
				return true;
			}

			virtual string __repr__()
			{
				return "<simuPOP::pyIndOperator>";
			}

		private:

			/// the function
			PyObject * m_func;

			/// parammeters
			PyObject * m_param;
	};

}
#endif