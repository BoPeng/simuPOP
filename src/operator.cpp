/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu                                                        *
*                                                                         *
*   $LastChangedDate$
*   $Rev$
*
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

#include "operator.h"

namespace simuPOP {

bool baseOperator::isActive(UINT rep, UINT numRep, long gen, long end, bool repOnly)
{
	// rep does not match
	if (!m_rep.match(rep, numRep))
		return false;

	// only check for rep value.
	if (repOnly)
		return true;

	if (gen < 0)
		return true;

	// if all active? (begin=0, end=-1)
	if (ISSETFLAG(m_flags, m_flagAtAllGen))
		return true;

	DBG_FAILIF(end > 0 && gen > end, IndexError,
		"Current generation can not be bigger than ending generation.");

	if (ISSETFLAG(m_flags, m_flagOnlyAtBegin) ) {
		if (gen == 0) return true;
		else return false;
	}

	if (ISSETFLAG(m_flags, m_flagOnlyAtEnd) && end > 0) {
		if (gen == end) return true;
		else return false;
	}

	// at Gen has higher priority.
	if (!m_atGen.empty() ) {
		// chech atGen.
		for (size_t i = 0, iEnd = m_atGen.size(); i < iEnd;  ++i) {
			int atGen = m_atGen[i];

			if (atGen > 0) {
				if (gen == atGen) return true;
				else continue;
			}                                                                               // atGen <=0
			else {
				if (end < 0)                                                                // can not determine.
					continue;
				// now end >= 0 atGen <=0
				if (end + atGen + 1 == gen)
					return true;
				else
					continue;
			}
		}
		// Do not check other parameters
		return false;
	}

	// finally, check start, end, every
	if (end < 0) {                                                            // if we do not know ending generation.
		// can not determine starting gen.
		if (m_beginGen < 0 || m_beginGen > gen)
			return false;

		// can only check positive begin + every
		if ( ( ((gen - m_beginGen) % m_stepGen) == 0) && (m_endGen < 0 || m_endGen >= gen) )
			return true;
		else
			return false;
	}                                                                                         // know ending generation
	else {
		int realStartGen = m_beginGen >= 0 ? m_beginGen : m_beginGen + end + 1;
		int realEndGen = m_endGen >= 0 ? m_endGen : m_endGen + end + 1;

		if (realStartGen > realEndGen)
			return false;

		return gen >= realStartGen && gen <= realEndGen && (gen - realStartGen) % m_stepGen == 0  ;
	}
	return false;
}


void baseOperator::setFlags()
{
	RESETFLAG(m_flags, m_flagAtAllGen);
	RESETFLAG(m_flags, m_flagOnlyAtBegin);
	RESETFLAG(m_flags, m_flagOnlyAtEnd);

	// atGen has higher priority: if it is not empty, use it.
	if (m_atGen.empty() ) {
		if (m_beginGen == 0 && m_endGen == 0)
			SETFLAG(m_flags, m_flagOnlyAtBegin);

		if (m_endGen == -1 && m_beginGen == -1)
			SETFLAG(m_flags, m_flagOnlyAtEnd);

		if (m_beginGen == 0 && m_endGen == -1 && m_stepGen == 1)
			SETFLAG(m_flags, m_flagAtAllGen);

	} else if (m_atGen.size() == 1) {
		if (m_stepGen < 1)
			throw IndexError("active generation interval should be at least 1.");

		if (m_atGen[0] == 0) SETFLAG(m_flags, m_flagOnlyAtBegin);
		if (m_atGen[0] == -1) SETFLAG(m_flags, m_flagOnlyAtEnd);
	}
}


bool baseOperator::apply(population & pop)
{
	if (m_lastPop != pop.genoStruIdx()) {
		initialize(pop);
		m_lastPop = pop.genoStruIdx();
	}
	return true;
}


bool baseOperator::applyDuringMating(population & pop, RawIndIterator offspring,
                                     individual * dad, individual * mom)
{
	if (m_lastPop != pop.genoStruIdx()) {
		initialize(pop);
		m_lastPop = pop.genoStruIdx();
	}
	return true;
}


bool pause::apply(population & pop)
{
	// call initialize if needed.
	baseOperator::apply(pop);

	char a;

	if (m_stopOnKeyStroke) {
		// check if key is already pressed
		if (!simuPOP_kbhit() )
			return true;
		else
			simuPOP_getch();
	}
	// clear input and wait for user input
	// std::cin.clear();
	// std::cin.ignore(std::numeric_limits<int>::max());

	if (m_prompt) {
		cout << "Simulation paused. " << endl
		     << " Press " << endl
		     << "   q to stop evolution, " << endl
		     << "   s to start an interative shell, (current population is ";
		if (m_exposePop)
			cout << "exported as " << m_popName << ')' << endl;
		else
			cout << "not exported)" << endl;
		cout << "   or any other key to continue...." << endl;
	}
	a = simuPOP_getch();                                              // std::cin.get(a);

	if (a == 'q' || a == 'Q')
		throw SystemError("Terminated by user");
	else if (a == 's' || a == 'S') {
		// export current population
		PyObject * popObj;
		if (m_exposePop) {
			popObj = pyPopObj(static_cast<void *>(&pop));
			if (popObj == NULL)
				throw SystemError("Could not expose population pointer. Compiled with the wrong version of SWIG? ");

			// get global dictionary
			mainVars().setVar(m_popName, popObj);
		}
		PyRun_InteractiveLoop(stdin, "<stdin>");
		// if expose pop, release it.
		if (m_exposePop)
			mainVars().removeVar(m_popName);
	}

	// clear input and wait for user input
	// std::cin.clear();
	// std::cin.ignore(std::numeric_limits<int>::max());

	return true;
}


bool ifElse::applyDuringMating(population & pop, RawIndIterator offspring,
                               individual * dad, individual * mom)
{
	m_cond.setLocalDict(pop.dict());
	bool res = m_cond.valueAsBool();

	if (res && m_ifOp != NULL)
		return m_ifOp->applyDuringMating(pop, offspring, dad, mom);
	else if (!res && m_elseOp != NULL)
		return m_elseOp->applyDuringMating(pop, offspring, dad, mom);
	return true;
}


bool ifElse::apply(population & pop)
{
	m_cond.setLocalDict(pop.dict());
	bool res = m_cond.valueAsBool();

	if (res && m_ifOp != NULL)
		return m_ifOp->apply(pop);
	else if (!res && m_elseOp != NULL)
		return m_elseOp->apply(pop);
	return true;
}


bool ticToc::apply(population & pop)
{
	time_t tmpTime;

	// this may not be correct but wrap is a possible problem.
	if (!this->noOutput() ) {
		ostream & out = this->getOstream(pop.dict());
		// since last time
		double timeDiff = difftime(time(&tmpTime), m_lastTime);
		out << "Elapsed Time: " << int(timeDiff * 100) / 100. ;
		// since beginning
		timeDiff = difftime(time(&tmpTime), m_startTime);
		int h = int(timeDiff / 3600);
		int m = int((timeDiff - h * 3600) / 60);
		int s = int(timeDiff - h * 3600 - m * 60);
		out << "s  Overall Time: " << std::setw(2) << std::setfill('0') << h
		    << ":" << std::setw(2) << std::setfill('0') << m << ":" << std::setw(2)
		    << std::setfill('0') << s << endl;
		this->closeOstream();
	}

	time(&m_lastTime);
	return true;
}


pyOperator::pyOperator(PyObject * func, PyObject * param,
	int stage, bool formOffGenotype, bool passOffspringOnly,
	int begin, int end, int step, vectorl at,
	const repList & rep, const subPopList & subPop, const vectorstr & infoFields) :
	baseOperator(">", "", stage, begin, end, step, at, rep, subPop, infoFields),
	m_func(func), m_param(param), m_passOffspringOnly(passOffspringOnly)
{
	if (!m_func.isValid())
		throw ValueError("Passed variable is not a callable Python function.");

	this->setFormOffGenotype(formOffGenotype);
}

bool pyOperator::apply(population & pop)
{
	// call the python function, pass the whole population in it.
	// get pop object
	PyObject * popObj = pyPopObj(static_cast<void *>(&pop));

	// if pop is valid?
	if (popObj == NULL)
		throw SystemError("Could not pass population to the provided function. \n"
			              "Compiled with the wrong version of SWIG?");

	// parammeter list, ref count increased
	bool resBool;
	// parenthesis is needed since PyCallFuncX are macros.
	if (m_param.isValid())
		resBool = m_func.call("(OO)", PyObj_As_Bool, popObj, m_param.object());
	else
		resBool = m_func.call("(O)", PyObj_As_Bool, popObj);

	Py_DECREF(popObj);
	return resBool;
}


bool pyOperator::applyDuringMating(population & pop, RawIndIterator offspring,
                                   individual * dad, individual * mom)
{
	// get offspring object
	PyObject * offObj = pyIndObj(static_cast<void *>(&(*offspring)));

	DBG_FAILIF(offObj == NULL, SystemError,
		"Could not pass offspring to the provided function. \n"
		"Compiled with the wrong version of SWIG?");

	bool res;
	if (m_passOffspringOnly) {
		// parammeter list, ref count increased
		if (m_param.isValid())
			res = m_func.call("(OO)", PyObj_As_Bool, offObj, m_param.object());
		else
			res = m_func.call("(O)", PyObj_As_Bool, offObj);

	} else {
		// call the python function, pass all the parameters to it.
		// get pop object
		PyObject * popObj = pyPopObj(static_cast<void *>(&pop));

		// get dad object
		PyObject * dadObj, * momObj;
		if (dad == NULL) {
			Py_INCREF(Py_None);
			dadObj = Py_None;
		} else
			dadObj = pyIndObj(static_cast<void *>(dad));

		if (mom == NULL) {
			Py_INCREF(Py_None);
			momObj = Py_None;
		} else
			momObj = pyIndObj(static_cast<void *>(mom));

		// if pop is valid?
		DBG_FAILIF(popObj == NULL || dadObj == NULL || momObj == NULL, SystemError,
			"Could not pass population or parental individuals to the provided function. \n"
			"Compiled with the wrong version of SWIG?");

		// parammeter list, ref count increased
		if (m_param.isValid())
			res  = m_func.call("(OOOOO)", PyObj_As_Bool, popObj, offObj, dadObj, momObj, m_param.object());
		else
			res = m_func.call("(OOOO)", PyObj_As_Bool, popObj, offObj, dadObj, momObj);

		Py_DECREF(popObj);
		Py_DECREF(dadObj);
		Py_DECREF(momObj);
	}
	Py_DECREF(offObj);

	return res;
}


void ApplyDuringMatingOperator(const baseOperator & op,
                               population * pop, int dad, int mom, ULONG off)
{
	baseOperator * opPtr = op.clone();

	opPtr->applyDuringMating(*pop, pop->rawIndBegin() + off,
		dad < 0 ? NULL : &pop->ind(dad),
		mom < 0 ? NULL : &pop->ind(mom));
}


}
