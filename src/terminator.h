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


namespace simuPOP {

/** This operator evaluates an expression in a population's local namespace
 *  and terminate the evolution of this population, or the whole simulator,
 *  if the return value of this expression is \c True. Termination caused by
 *  an operator will stop the execution of all operators after it. Because
 *  a life-cycle is considered to be complete if mating is complete, the
 *  <em>evolved generations</em> (return value from <tt>simulator::evolve</tt>)
 *  of a terminated replicate is determined by when the last evolution cycle is
 *  terminated.
 */
class terminateIf : public baseOperator
{

public:
	/** Create a terminator with an expression \e condition, which will be evalulated
	 *  in a population's local namespace when the operator is applied to this
	 *  population. If the return value of \e condition is \c True, the evolution
	 *  of the population will be terminated. If \e stopAll is set to \c True, the
	 *  evolution of all replicates of the simulator will be terminated. If this
	 *  operator is allowed to write to an \e output (default to ""), the generation
	 *  number, preceeded with an optional \e message will be
	 *  written to it.
	 */
	terminateIf(string condition = "", bool stopAll = false, string message = "",
		string output = "", 
		int stage = PostMating, int begin = 0, int end = -1, int step = 1, vectorl at = vectorl(),
		const repList & rep = repList(), const subPopList & subPop = subPopList(), const vectorstr & infoFields = vectorstr()) :
		baseOperator(output, stage, begin, end, step, at, rep, subPop, infoFields),
		m_expr(condition), m_stopAll(stopAll), m_message(message)
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
	virtual bool apply(population & pop);

	virtual ~terminateIf()
	{
	}


private:
	/// alleles to check. If empty, check all alleles.
	Expression m_expr;

	///
	bool m_stopAll;

	/// message to print when terminated
	string m_message;
};

}

#endif
