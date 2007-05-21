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

namespace simuPOP
{

    /// base class of all classes that manipulate populations
	/** 
	Operators are objects that act on populations. They can be
    applied	to populations directly using their function forms,
    but they are usually managed and applied by a simulator. \n

	Operators can be applied at different stages of the life cycle of
    a generation. More specifically, they can be applied at \em pre-,
    \em during-, \em post-mating, or a combination of these stages. Applicable stages
    are usually set by default but you can change it by setting
    <tt>stage=(PreMating|PostMating|DuringMating|PrePostMating)</tt> parameter.
    Some operators ignore \c stage parameter because they only
    work at one stage. \n

	Operators do not have to be applied at all generations. You can specify
	starting/ending generation, gaps between applicable generations,
	or even specific generations. For example, you might want to
	start applying migrations after certain burn-in generations, or
    calculate certain statistics only sparsely. \n

	Operators can have outputs. Output can be standard (terminal)
	or a file, which can vary with replicates and/or generations.
    Outputs from different operators can be accumulated to the same file to form
    table-like outputs. \n

    Operators are applied to every replicate of a simulator by default.
    However, you can apply operators to one or a group of
	replicates using parameter \c rep or \c grp. \n
    
	Filenames can have the following format:

	\li \c 'filename' this file will be overwritten each time. If two operators
      output to the same file, only the last one will succeed;

	\li \c '>filename' the same as \c 'filename';

	\li <tt>'>>filename'</tt> the file will be created at the beginning of evolution
	  (\c simulator::evolve) and closed at the end. Output from several operators
      is allowed;

	\li <tt>'>>>filename'</tt> the same as <tt>'>>filename'</tt> except that the file will not
	  be cleared at the beginning of evolution if it is not empty;

	\li \c '>' standard output (terminal);

	\li \c '' supress output.
	*/
	class Operator
	{
		public:

			/** @name constructor and destructor */
			//@{

			/// common interface for all operators (this base operator does nothing by itself.)
			/**

			\param begin the starting generation. Default to \c 0. Negative numbers are allowed.
			\param end stop applying after this generation. Negative numbers are allowed.
			\param step the number of generations between active generations. Default to \c 1.
			\param at an array of active generations. If given, \c stage, \c begin, \c end,
			  and \c step will be ignored.
			\param rep applicable replicates. It can be a valid replicate number, \c REP_ALL
              (all replicates, default), or \c REP_LAST (only the last replicate). \c REP_LAST
              is useful in adding newlines to a table output.
			\param grp applicable group. Default to \c GRP_ALL.  A group number for each
			  replicate is set by <tt>simulator.__init__</tt> or <tt>simulator::setGroup()</tt>.
			\param output a string of the output filename. Different operators will have
			  different default \c output (most commonly \c '>' or \c ''). 
			\param outputExpr an expression that determines the output filename dynamically. This
              expression will be evaluated against a population's local namespace each time when
              an output filename is required. For example, <tt> "'>>out%s_%s.xml' % (gen, rep)" </tt>
              will output to <tt> >>>out1_1.xml </tt> for replicate \c 1 at generation \c 1.
            
            \note Negative generation numbers are allowed for \c begin, \c end and \c at. They are
              intepretted as <tt>endGen + gen + 1</tt>. For example, <tt>begin = -2</tt> in
              <tt>simu.evolve(..., end=20)</tt> starts at generation \c 19.
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

			/// deep copy of an operator
			virtual Operator* clone() const
			{
				return new Operator(*this);
			}

			//@}

			/** @name  applicable generations (also judge from rep and group). use of parameter start, end, every, at, group, rep
			*/
			//@{

			/// CPPONLY determine if this operator is active
			/**
            Determine if this operator is active under the conditions such as the current
            replicate, group number of the current replicate, current generation, ending generation etc.
			\note This function will be called by simulators before applying.
			*/
			bool isActive(UINT rep, UINT numRep, long gen, long end, int grp, bool repOnly =false);

			/// return applicable group
			int applicableGroup()
			{
				return m_grp;
			}

			/// set applicable group
			/**
			Default to \c GRP_ALL (applicable to all groups).
			Otherwise, the operator is applicable to only \em one group of replicates.
			Groups can be set in \c simulator::setGroup().
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

			/// set applicable replicate
			void setApplicableReplicate(int rep)
			{
				m_rep = rep;
			}

			/// set applicable generation parameters: \c begin, \c end, \c step and \c at
			void setActiveGenerations(int begin=0, int end=-1, int step=1, vectorl at=vectorl())
			{
				m_beginGen = begin;
				m_endGen = end;
				m_stepGen = step;
				m_atGen = at;

				DBG_FAILIF(step<=0, ValueError, "step need to be at least one");

				// set certain m_flags to speed up using this machanism
				setFlags();
			}

			//@}

			/** @name applicable stages	pre, during, post-mating methods */
			//@{

			/// set applicable stage. Another way to set \c stage parameter.
			void setApplicableStage(int stage)
			{
				RESETFLAG(m_flags, PreDuringPostMating);
				SETFLAG(m_flags, stage);
			}

			/// set if this operator can be applied \em pre-mating
			bool canApplyPreMating()
			{
				return ISSETFLAG(m_flags, m_flagPreMating);
			}

			/// set if this operator can be applied \em during-mating
			bool canApplyDuringMating()
			{
				return ISSETFLAG(m_flags, m_flagDuringMating);
			}

			/// set if this operator can be applied \em post-mating
			bool canApplyPostMating()
			{
				return ISSETFLAG(m_flags, m_flagPostMating);
			}

			/// set of this operator can be applied \em pre- or \em post-mating
			bool canApplyPreOrPostMating()
			{
				return ISSETFLAG(m_flags, m_flagPreMating) || ISSETFLAG(m_flags, m_flagPostMating);
			}

            /// CPPONLY
			virtual bool isCompatible(const population & pop)
			{
#ifdef SIMUMPI
				DBG_ASSERT(MPIReady(), ValueError,
					"Operator " + __repr__() + " is not MPI ready. ");
				return true;
#else
				return true;
#endif
			}

			/// determine if the operator can be applied only for haploid population
			bool haploidOnly()
			{
				return ISSETFLAG(m_flags, m_flagHaploid);
			}

			/// determine if the operator can be applied only for diploid population
			bool diploidOnly()
			{
				return ISSETFLAG(m_flags, m_flagDiploid);
			}

			/// determine if this operator can be used in a MPI module
			bool MPIReady()
			{
				return ISSETFLAG(m_flags, m_flagMPI);
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

			/// CPPONLY set the operator is ready for MPI version
			void setMPIReady()
			{
				SETFLAG(m_flags, m_flagMPI);
			}

			/// get the length of information fields for this operator
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


			/// CPPONLY if the operator will form genotype of offspring
			/**
			If none of the during mating operator can form offspring, default will be used.
			*/
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

			/// CPPONLY provide interface to apply operator before during or after mating
			virtual bool applyWithScratch(population& pop, population& scratch, int stage)
			{
				return apply(pop);
			}

			/// apply to one population. It does not check if the operator is activated.
			virtual bool apply(population& pop)
			{
				throw SystemError("This function Operator::apply() is not supposed to be called.");
			}

			/// CPPONLY apply during mating, given \c pop, \c offspring, \c dad and \c mom
			virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
				individual* dad=NULL, individual* mom=NULL)
			{
				throw SystemError("ApplyDuringMating: This function is not supposed to be called.");
			}

			//@}
			/** @name dealing with output separator, persistant files, $gen etc substitution.
			*/
			//@{

			/// set ouput stream, if was not set during construction
			void setOutput(string output="", string outputExpr="")
			{
				m_ostream.setOutput(output, outputExpr);
			}

			/// CPPONLY get output stream. This function is not exposed to user.
			ostream& getOstream( PyObject* dict=NULL, bool readable=false)
			{
				return m_ostream.getOstream( dict, readable);
			}

			/// CPPONLY close output stream and delete output stream pointer, if it is a output stream
			void closeOstream()
			{
				m_ostream.closeOstream();
			}

			/// CPPONLY say something about active states
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

            /// used by Python print function to print out the general information of the operator
			virtual string __repr__()
			{
				return "<simuPOP::operator> " ;
			}

			//@}

			/// CPPONLY determine if \c output=">". Used internally.
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
			// limited to haploid?
			static const size_t m_flagHaploid        = 128;
			// limited to diploid?
			static const size_t m_flagDiploid        = 256;
			// can be used for MPI version?
			static const size_t m_flagMPI            = 512;

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

	/**
	Pause the evolution of a simulator at given generations or at
	key stroke, using \c stopOnKeyStroke=True option. When a simulator
	is stopped, user can resume simulation by pressing '/??' or escape
	to a python shell to examine the status of the simulation.
	*/
	class pause: public Operator
	{

		public:
			/// stop simulation. press q to exit and any other key to continue
            /**
			\param prompt if true (default), print prompt message. findpause
			\param stopOnKeyStroke if true, go on if no key was pressed
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
			virtual bool apply(population& pop);

            /// used by Python print function to print out the general information of the pause
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

			///  providing interface to apply operator before during or after mating
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
			virtual bool applyWithScratch(population& pop, population& scratch, int stage);

			/// give pop, offspring, pop and mom.
			virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
				individual* dad=NULL, individual* mom=NULL);

			virtual bool apply(population& pop);

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

			virtual bool apply(population& pop);

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

			virtual bool apply(population& pop);

			virtual bool applyDuringMating(population& pop, population::IndIterator offspring,
				individual* dad=NULL, individual* mom=NULL);

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

			virtual bool apply(population& pop);

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
