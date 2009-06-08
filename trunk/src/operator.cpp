/**
 *  $File: operator.cpp $
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

#include "operator.h"

namespace simuPOP {

bool baseOperator::isActive(UINT rep, long gen, long end,
                            const vector<bool> & activeRep, bool repOnly)
{
	// rep does not match
	if (!m_rep.match(rep, activeRep))
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


opList::opList(const vectorop & values) : m_elems(0)
{
	vectorop::const_iterator it = values.begin();
	vectorop::const_iterator itEnd = values.end();

	for (; it != itEnd; ++it)
		m_elems.push_back((*it)->clone());
}


opList::opList(const baseOperator & op) : m_elems(0)
{
	m_elems.push_back(op.clone());
}


opList::opList(const opList & ops) : m_elems(0)
{
	vectorop::const_iterator it = ops.m_elems.begin();
	vectorop::const_iterator itEnd = ops.m_elems.end();

	for (; it != itEnd; ++it)
		m_elems.push_back((*it)->clone());
}


opList::~opList()
{
	vectorop::const_iterator it = m_elems.begin();
	vectorop::const_iterator itEnd = m_elems.end();

	for (; it != itEnd; ++it)
		delete *it;
}


vectori pause::s_cachedKeys = vectori();

bool pause::apply(population & pop)
{
	// call initialize if needed.
	baseOperator::apply(pop);

	char a;

	if (m_stopOnKeyStroke != static_cast<char>(false)) {
		// stop on any key
		if (m_stopOnKeyStroke == static_cast<char>(true)) {
			// check if key is already pressed
			if (!simuPOP_kbhit() )
				return true;
		} else {
			// get all keys
			while (simuPOP_kbhit())
				s_cachedKeys.push_back(simuPOP_getch());
			// if the required key is pressed?
			vectori::iterator it = find(s_cachedKeys.begin(), s_cachedKeys.end(),
				static_cast<int>(m_stopOnKeyStroke));
			if (it == s_cachedKeys.end())
				// also look for C-c
				it = find(s_cachedKeys.begin(), s_cachedKeys.end(), 3);
			if (it == s_cachedKeys.end())
				return true;
			s_cachedKeys.erase(it);
		}
	}
	// clear input and wait for user input
	std::cin.clear();

	string popName = "pop_" + toStr(pop.gen()) + "_" + toStr(pop.rep());
	if (m_prompt) {
		cout << "Simulation paused for population " << pop.rep() << "\n"
		     << "Press\n"
		     << "   's' to (s)top the evolution of this population,\n"
		     << "   'q' to quit (stop the evolution of all populations),\n"
		     << "   'p' to start an interative (P)ython shell, (current population will be exported as "
		     << popName << ")\n"
		     << "   'r' or any other key to (r)esume evolution...." << endl;
	}
	a = simuPOP_getch();

	if (a == 's' || a == 'S') {
		cout << "Evolution of population " << pop.rep() << " is stopped." << endl;
		return false;
	} else if (a == 'q' || a == 'Q') {
		cout << "Evolution of all populations are terminated." << endl;
		throw StopEvolution(string());
	}
	if (a == 'p' || a == 'P') {
		// export current population
		PyObject * popObj;
		popObj = pyPopObj(static_cast<void *>(&pop));
		if (popObj == NULL)
			throw SystemError("Could not expose population pointer. Compiled with the wrong version of SWIG? ");

		// get global dictionary
		mainVars().setVar(popName, popObj);
		PyRun_InteractiveLoop(stdin, "<stdin>");
		// if expose pop, release it.
		if (mainVars().hasVar(popName))
			mainVars().removeVar(popName);
	}
	cout << "Resume evolution of population " << pop.rep() << "..." << endl;

	return true;
}


bool ifElse::applyDuringMating(population & pop, RawIndIterator offspring,
                               individual * dad, individual * mom)
{
	m_cond.setLocalDict(pop.dict());
	bool res = m_cond.valueAsBool();

	if (res && !m_ifOps.empty()) {
		const vectorop & ops = m_ifOps.elems();
		vectorop::const_iterator it = ops.begin();
		vectorop::const_iterator itEnd = ops.end();
		for (; it != itEnd; ++it) {
			bool ret = (*it)->applyDuringMating(pop, offspring, dad, mom);
			if (!ret)
				return false;
		}
		return true;
	} else if (!res && !m_elseOps.empty()) {
		const vectorop & ops = m_elseOps.elems();
		vectorop::const_iterator it = ops.begin();
		vectorop::const_iterator itEnd = ops.end();
		for (; it != itEnd; ++it) {
			bool ret = (*it)->applyDuringMating(pop, offspring, dad, mom);
			if (!ret)
				return false;
		}
		return true;
	}
	return true;
}


bool ifElse::apply(population & pop)
{
	m_cond.setLocalDict(pop.dict());
	bool res = m_cond.valueAsBool();

	if (res && !m_ifOps.empty()) {
		const vectorop & ops = m_ifOps.elems();
		vectorop::const_iterator it = ops.begin();
		vectorop::const_iterator itEnd = ops.end();
		for (; it != itEnd; ++it) {
			bool ret = (*it)->apply(pop);
			if (!ret)
				return false;
		}
		return true;
	} else if (!res && !m_elseOps.empty()) {
		const vectorop & ops = m_elseOps.elems();
		vectorop::const_iterator it = ops.begin();
		vectorop::const_iterator itEnd = ops.end();
		for (; it != itEnd; ++it) {
			bool ret = (*it)->apply(pop);
			if (!ret)
				return false;
		}
		return true;
	}
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
	int stage, bool isTransmitter, bool passOffspringOnly,
	int begin, int end, int step, const intList & at,
	const repList & rep, const subPopList & subPops,
	const stringList & infoFields) :
	baseOperator(">", stage, begin, end, step, at, rep, subPops, infoFields),
	m_func(func), m_param(param), m_passOffspringOnly(passOffspringOnly)
{
	if (!m_func.isValid())
		throw ValueError("Passed variable is not a callable Python function.");

	this->setTransmitter(isTransmitter);
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
		resBool = m_func(PyObj_As_Bool, "(OO)", popObj, m_param.object());
	else
		resBool = m_func(PyObj_As_Bool, "(O)", popObj);

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
			res = m_func(PyObj_As_Bool, "(OO)", offObj, m_param.object());
		else
			res = m_func(PyObj_As_Bool, "(O)", offObj);

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
			res = m_func(PyObj_As_Bool, "(OOOOO)", popObj, offObj, dadObj, momObj, m_param.object());
		else
			res = m_func(PyObj_As_Bool, "(OOOO)", popObj, offObj, dadObj, momObj);

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
