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
 \brief head file of class terminator: public baseOperator
 */
#include "population.h"
#include "operator.h"

#include <iostream>
#include <iomanip>

namespace simuPOP {

/// Base class of all terminators
/**
   Teminators are used to see if an evolution is running as expected, and
   terminate the evolution if a certain condition fails.
 */
class terminator : public baseOperator
{

public:
	/// create a terminator
	/**
	 \param message a message that will be displayed when the evolution is terminated.
     \param stopAll stop all replicates if this replicate is stopped.
	 */
	terminator(string message = "", bool stopAll=false,
                string output = ">", string outputExpr = "",
	           int stage = PostMating, int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(), int rep = REP_ALL, const vectorstr & infoFields = vectorstr()) :
		baseOperator(output, outputExpr, stage, begin, end, step, at, rep, infoFields),
		m_message(message), m_stopAll(stopAll)
	{
	};

	/// destructor
	virtual ~terminator()
	{
	};

	/// deep copy of a terminator
	virtual baseOperator * clone() const
	{
		return new terminator(*this);
	}


protected:
	/// message to print when terminated
	string m_message;

    /// 
    bool m_stopAll;
};

/// terminate according to a condition
/**
   This operator terminates the evolution under certain conditions. For example,
   <tt>terminateIf(condition='alleleFreq[0][1]<0.05', begin=100)</tt>
   terminates the evolution if the allele frequency of allele \c 1 at locus \c 0
   is less than 0.05. Of course, to make this opertor work, you will need to use
   a \c stat operator before it so that variable \c alleleFreq exists in the local namespace. \n

   When the value of condition is \c True, a shared variable <tt>var="terminate"</tt> will be
   set to the current generation.
 */
class terminateIf : public terminator
{

public:
	/// create a \c terminateIf terminator
	terminateIf(string condition = "", bool stopAll=false, string message = "", string var = "terminate",
	            string output = "", string outputExpr = "",
	            int stage = PostMating, int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
	            int rep = REP_ALL, const vectorstr & infoFields = vectorstr()) :
		terminator(message, stopAll, output, outputExpr, stage, begin, end, step, at,
		           rep), m_expr(condition), m_var(var)
	{
	}


	/// deep copy of a \c terminateIf terminator
	virtual baseOperator * clone() const
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
	virtual bool apply(population & pop);

	virtual ~terminateIf()
	{
	};

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
class continueIf : public terminator
{

public:
	/// create a \c continueIf terminator
	continueIf(string condition = "", bool stopAll=false, string message = "", string var = "terminate",
	           string output = "", string outputExpr = "",
	           int stage = PostMating, int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
	           int rep = REP_ALL,const vectorstr & infoFields = vectorstr()) :
		terminator(message, stopAll, output, outputExpr, stage, begin, end, step, at,
		           rep), m_expr(condition), m_var(var)
	{
	}


	/// deep copy of a \c continueIf terminator
	virtual baseOperator * clone() const
	{
		return new continueIf(*this);
	}


	/// used by Python print function to print out the general information of the \c continueIf terminator
	virtual string __repr__()
	{
		return "<simuPOP::terminateIf>";
	}


	virtual bool apply(population & pop);

	virtual ~continueIf()
	{
	};

private:
	/// alleles to check. If empty, check all alleles.
	Expression m_expr;

	/// variable to set when terminated
	string m_var;
};

}
#endif
