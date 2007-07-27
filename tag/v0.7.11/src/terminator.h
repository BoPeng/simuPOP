/***************************************************************************
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

#ifndef _TERMINATOR_H
#define _TERMINATOR_H
/**
\file
\brief head file of class terminator: public Operator
*/
#include "population.h"
#include "operator.h"

#include <iostream>
#include <iomanip>

namespace simuPOP
{

    /// terminate the evolution
    /**
    These operators are used to see if an evolution is running as expected, and
    terminate the evolution if a certain condition fails.
    */
	class terminator: public Operator
	{

		public:
			/// create a terminator, default to be always active
			terminator(string message = "", string output=">", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(), int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			Operator(output, outputExpr, stage, begin, end, step, at, rep, grp, infoFields),
				m_message(message)
			{
			};

			/// destructor
			virtual ~terminator(){};

            /// deep copy of a terminator
			virtual Operator* clone() const
			{
				return new terminator(*this);
			}

            /// return the message to print when terminated???
			string message()
			{
				return m_message;
			}

		private:
			/// message to print when terminated
			string m_message;
	};

	/// terminate according to a condition
	/**
	This operator terminates the evolution under certain conditions. For example,
    <tt>terminateIf(condition='alleleFreq[0][1]<0.05', begin=100)</tt>
    terminates the evolution if the allele frequency of allele \c 1 at locus \c 0
    is less than 0.05. Of course, to make this opertor work, you will need to use
    a \c stat operator before it so that variable \c alleleFreq exists in the local namespace. \n

    When the condition is true, a shared variable <tt>var="terminate"</tt> will be
	set to the current generation.
	*/
	class terminateIf: public terminator
	{

		public:
            /// create a \c terminateIf terminator
			terminateIf(string condition="", string message="", string var="terminate",
				string output="", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			terminator(message, output, outputExpr, stage, begin, end, step, at,
				rep, grp), m_expr(condition ), m_var(var)
			{
			}

            /// deep copy of a \c terminateIf terminator
			virtual Operator* clone() const
			{
				return new terminateIf(*this);
			}

            /// used by Python print function to print out the general information of the \c terminateIf terminator
			virtual string __repr__()
			{
				return "<simuPOP::terminateIf>";
			}

			// check all alleles in vector allele if they are fixed.
			/// apply the \c terminateIf terminator
			virtual bool apply(population& pop)
			{
				// experssion return true
				m_expr.setLocalDict(pop.dict());

				if( m_expr.valueAsBool() == true )
				{
					// store the generations this replicate stops
					int gen = pop.gen();		  // mainVars().getVarAsInt("gen");

					pop.setIntVar(m_var, gen);
					if( !this->noOutput() )
					{
						ostream & out = this->getOstream(pop.dict());
						out << gen << endl;
						this->closeOstream();
					}
					if(this->message() != "")
						cout << this->message() << endl;
					return false;				  // return false, this replicate will be stopped
				}
				else
					return true;
			}

			virtual ~terminateIf(){};

		private:
			/// alleles to check. If empty, check all alleles.
			Expression m_expr;

			/// variable to set when terminated
			string m_var;
	};

	/// terminate according to a condition failure
	/** 
	The same as \c terminateIf but continue if the condition is \c True.
    */
	class continueIf: public terminator
	{

		public:
            /// create a \c continueIf terminator
			continueIf(string condition="", string message="", string var="terminate",
				string output="", string outputExpr="",
				int stage=PostMating, int begin=0, int end=-1, int step=1, vectorl at=vectorl(),
				int rep=REP_ALL, int grp=GRP_ALL, const vectorstr& infoFields=vectorstr()):
			terminator(message, output, outputExpr, stage, begin, end, step, at,
				rep, grp), m_expr(condition ), m_var(var)
			{
			}

            /// deep copy of a \c continueIf terminator
			virtual Operator* clone() const
			{
				return new continueIf(*this);
			}

            /// used by Python print function to print out the general information of the \c continueIf terminator
			virtual string __repr__()
			{
				return "<simuPOP::terminateIf>";
			}

			// check all alleles in vector allele if they are fixed.
			/// apply the \c continueIf terminator???
			virtual bool apply(population& pop)
			{
				// experssion return true
				m_expr.setLocalDict(pop.dict());

				if( m_expr.valueAsBool() == false )
				{
					// store the generations this replicate stops
					int gen = pop.gen();		  // mainVars().getVarAsInt("gen");

					pop.setIntVar(m_var, gen);
					if( !this->noOutput() )
					{
						ostream & out = this->getOstream(pop.dict());
						out << gen << endl;
						this->closeOstream();
					}
					if(this->message() != "")
						cout << this->message() << endl;
					return false;				  // return false, this replicate will be stopped
				}
				else
					return true;
			}

			virtual ~continueIf(){};

		private:
			/// alleles to check. If empty, check all alleles.
			Expression m_expr;

			/// variable to set when terminated
			string m_var;
	};

}
#endif
