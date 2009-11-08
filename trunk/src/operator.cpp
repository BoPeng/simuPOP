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

subPopList::subPopList(PyObject * obj) : m_subPops(), m_allAvail(false)
{
	if (obj == NULL || obj == Py_None)
		// accept NULL
		m_allAvail = true;
	else if (PyBool_Check(obj)) {
		// accept True/False
		m_allAvail = obj == Py_True;
		return;
	} else if (PyNumber_Check(obj)) {
		// accept a number
		m_allAvail = false;
		m_subPops.push_back(vspID(PyInt_AS_LONG(obj)));
	} else if (PySequence_Check(obj)) {
		m_subPops.resize(PySequence_Size(obj));
		// assign values
		for (size_t i = 0, iEnd = m_subPops.size(); i < iEnd; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			if (PyNumber_Check(item))           // subpopulation
				m_subPops[i] = vspID(PyInt_AS_LONG(item));
			else if (PySequence_Check(item)) {  // virtual subpopulation
				DBG_ASSERT(PySequence_Size(item) == 2, ValueError,
					"Invalid virtual subpopulation ID");
				PyObject * sp = PySequence_GetItem(item, 0);
				PyObject * vsp = PySequence_GetItem(item, 1);
				DBG_ASSERT(PyNumber_Check(sp), ValueError, "Invalid input for a list of (virtual) subpopulations.");
				DBG_ASSERT(PyNumber_Check(vsp), ValueError, "Invalid input for a list of (virtual) subpopulations.");
				m_subPops[i] = vspID(PyInt_AS_LONG(sp), PyInt_AS_LONG(vsp));
				Py_DECREF(sp);
				Py_DECREF(vsp);
			} else {
				DBG_FAILIF(true, ValueError, "Invalid input for a list of (virtual) subpopulations.");
			}
			Py_DECREF(item);
		}
	} else {
		DBG_FAILIF(true, ValueError, "Invalid input for a list of (virtual) subpopulations.");
	}
}


subPopList::subPopList(const vectorvsp & subPops)
	: m_subPops(subPops), m_allAvail(false)
{
	for (size_t i = 0; i < m_subPops.size(); ++i) {
		DBG_ASSERT(m_subPops[i].valid(), ValueError,
			"Invalid subpopulation ID");
	}
}


bool baseOperator::isActive(UINT rep, long gen)
{
	if (!m_reps.match(rep))
		return false;

	DBG_FAILIF(gen < 0, ValueError, "Negative generation number.");

	// if all active? (begin=0, end=-1)
	if (ISSETFLAG(m_flags, m_flagAtAllGen))
		return true;

	if (ISSETFLAG(m_flags, m_flagOnlyAtBegin) )
		return gen == 0;

	// at gen has higher priority.
	if (!m_atGen.empty() ) {
		// chech atGen.
		for (size_t i = 0, iEnd = m_atGen.size(); i < iEnd;  ++i) {
			int atGen = m_atGen[i];
			DBG_FAILIF(atGen < 0, ValueError,
				"During-mating operator used in a mating scheme does not support negative generation index.");
			if (gen == atGen)
				return true;
			else
				continue;
		}
		// Do not check other parameters
		return false;
	}

	// finally, check start, end, every
	DBG_FAILIF(m_beginGen < 0, ValueError,
		"During-mating operator used in a mating scheme does not support negative generation index.");

	if (m_endGen >= 0 && m_beginGen > m_endGen)
		return false;

	return gen >= m_beginGen && (m_endGen < 0 || gen <= m_endGen) && (gen - m_beginGen) % m_stepGen == 0;
}


bool baseOperator::isActive(UINT rep, long gen, long end,
                            const vector<bool> & activeRep, bool repOnly)
{
	// rep does not match
	if (!m_reps.match(rep, &activeRep))
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

	if (ISSETFLAG(m_flags, m_flagOnlyAtBegin) )
		return gen == 0;

	if (ISSETFLAG(m_flags, m_flagOnlyAtEnd) && end > 0)
		return gen == end;

	// at Gen has higher priority.
	if (!m_atGen.empty() ) {
		// chech atGen.
		for (size_t i = 0, iEnd = m_atGen.size(); i < iEnd;  ++i) {
			int atGen = m_atGen[i];

			if (atGen >= 0) {
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


void baseOperator::initializeIfNeeded(const population & pop)
{
	if (m_lastPop != pop.genoStruIdx()) {
		initialize(pop);
		m_lastPop = pop.genoStruIdx();
	}
}


bool baseOperator::apply(population & pop)
{
	DBG_FAILIF(true, RuntimeError,
		"This operator can only be applied during mating.");
	return true;
}


bool baseOperator::applyDuringMating(population & pop, RawIndIterator offspring,
                                     individual * dad, individual * mom)
{
	DBG_FAILIF(true, RuntimeError,
		"This operator cannot be applied during mating.");
	return true;
}


string baseOperator::applicability(bool subPops, bool gen)
{
	string desc;

	if (gen) {
		if (ISSETFLAG(m_flags, m_flagAtAllGen))
			desc += "";
		else if (ISSETFLAG(m_flags, m_flagOnlyAtBegin) )
			desc += "at generation 0";
		else if (ISSETFLAG(m_flags, m_flagOnlyAtEnd) )
			desc += "at ending generation";
		else if (!m_atGen.empty() ) {
			if (m_atGen.size() == 1)
				desc += "at generation";
			else
				desc += "at generations";
			for (size_t i = 0; i < m_atGen.size(); ++i) {
				if (i == 0)
					desc += ", ";
				desc += " " + toStr(m_atGen[i]);
			}
		} else {
			if (m_beginGen != 0)
				desc += "begin at " + toStr(m_beginGen) + " ";
			if (m_endGen != -1)
				desc += "end at " + toStr(m_endGen) + " ";
			if (m_stepGen != 1)
				desc += "at interval " + toStr(m_stepGen);
		}
	}
	if (subPops) {
		if (m_subPops.allAvail())
			desc += "";
		else {
			if (desc.size() != 1)
				desc += ", ";
			desc += "to subpopulations ";
			for (size_t i = 0; i < m_subPops.size(); ++i) {
				vspID sp = m_subPops[i];
				if (i != 0)
					desc += ", ";
				if (sp.isVirtual())
					desc += "(" + toStr(sp.subPop()) + ", " + toStr(sp.virtualSubPop()) + ")";
				else
					desc += toStr(sp.subPop());
			}
		}
	}
	if (desc.empty())
		return desc;
	return "(" + desc + ")";
}


opList::opList(const vectorop & ops) : m_elems(0)
{
	vectorop::const_iterator it = ops.begin();
	vectorop::const_iterator itEnd = ops.end();

	for (; it != itEnd; ++it)
		m_elems.push_back((*it)->clone());
}


opList::opList(const baseOperator & op) : m_elems(0)
{
	m_elems.push_back(op.clone());
}


opList::opList(const opList & rhs) : m_elems(0)
{
	vectorop::const_iterator it = rhs.m_elems.begin();
	vectorop::const_iterator itEnd = rhs.m_elems.end();

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
	initializeIfNeeded(pop);

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
		cerr << "Simulation paused for population " << pop.rep() << "\n"
		     << "Press\n"
		     << "   's' to (s)top the evolution of this population,\n"
		     << "   'q' to quit (stop the evolution of all populations),\n"
		     << "   'p' to start an interative (P)ython shell, (current population will be exported as "
		     << popName << ")\n"
		     << "   'r' or any other key to (r)esume evolution...." << endl;
	}
	a = simuPOP_getch();

	if (a == 's' || a == 'S') {
		cerr << "Evolution of population " << pop.rep() << " is stopped." << endl;
		return false;
	} else if (a == 'q' || a == 'Q') {
		cerr << "Evolution of all populations are terminated." << endl;
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
	cerr << "Resume evolution of population " << pop.rep() << "..." << endl;

	return true;
}


ifElse::ifElse(PyObject * cond, const opList & ifOps, const opList & elseOps,
	const stringFunc & output, int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	const stringList & infoFields) :
	baseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_cond(), m_fixedCond(-1), m_ifOps(ifOps), m_elseOps(elseOps)
{
	if (PyString_Check(cond))
		m_cond.setExpr(PyString_AsString(cond));
	else {
		bool c;
		PyObj_As_Bool(cond, c);
		m_fixedCond = c ? 1 : 0;
	}
}


bool ifElse::applyDuringMating(population & pop, RawIndIterator offspring,
                               individual * dad, individual * mom)
{
	m_cond.setLocalDict(pop.dict());
	bool res = m_fixedCond == -1 ? m_cond.valueAsBool() : m_fixedCond == 1;

	if (res && !m_ifOps.empty()) {
		opList::const_iterator it = m_ifOps.begin();
		opList::const_iterator itEnd = m_ifOps.end();
		for (; it != itEnd; ++it) {
			if (!(*it)->isActive(pop.rep(), pop.gen()))
				continue;
			bool ret = (*it)->applyDuringMating(pop, offspring, dad, mom);
			if (!ret)
				return false;
		}
		return true;
	} else if (!res && !m_elseOps.empty()) {
		opList::const_iterator it = m_elseOps.begin();
		opList::const_iterator itEnd = m_elseOps.end();
		for (; it != itEnd; ++it) {
			if (!(*it)->isActive(pop.rep(), pop.gen()))
				continue;
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
	bool res = m_fixedCond == -1 ? m_cond.valueAsBool() : m_fixedCond == 1;

	if (res && !m_ifOps.empty()) {
		const vectorop & ops = m_ifOps.elems();
		vectorop::const_iterator it = ops.begin();
		vectorop::const_iterator itEnd = ops.end();
		for (; it != itEnd; ++it) {
			if (!(*it)->isActive(pop.rep(), pop.gen()))
				continue;
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
			if (!(*it)->isActive(pop.rep(), pop.gen()))
				continue;
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
	bool passOffspringOnly,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	const stringList & infoFields) :
	baseOperator(">", begin, end, step, at, reps, subPops, infoFields),
	m_func(func), m_param(param), m_passOffspringOnly(passOffspringOnly)
{
	if (!m_func.isValid())
		throw ValueError("Passed variable is not a callable Python function.");
}


string pyOperator::describe(bool format)
{
	PyObject * name = PyObject_GetAttrString(m_func.func(), "__name__");
	DBG_FAILIF(name == NULL, RuntimeError, "Passwd object does not have attribute __name__.");
	return "<simuPOP.pyOperator> calling a Python function " + string(PyString_AsString(name));
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
