/**
 *  $File: operator.cpp $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of simuPOP, a forward-time population genetics
 *  simulation environment. Please visit http://simupop.sourceforge.net
 *  for details.
 *
 *  Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
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
#include <sstream>
#include <locale>
using std::ostringstream;

#if PY_VERSION_HEX >= 0x03000000
#  define PyString_Check PyUnicode_Check
#endif

#include "utility.h"

namespace simuPOP {

bool BaseOperator::isActive(ssize_t rep, ssize_t gen) const
{
	if (!m_reps.match(rep))
		return false;

	DBG_FAILIF(gen < 0, ValueError, "Negative generation number.");

	// if all active? (begin=0, end=-1)
	if (ISSETFLAG(m_flags, m_flagAtAllGen))
		return true;

	if (ISSETFLAG(m_flags, m_flagOnlyAtBegin))
		return gen == 0;

	// at gen has higher priority.
	if (!m_atGen.empty()) {
		// chech atGen.
		for (size_t i = 0, iEnd = m_atGen.size(); i < iEnd; ++i) {
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


bool BaseOperator::isActive(size_t rep, ssize_t gen, ssize_t end,
                            const vector<bool> & activeRep, bool repOnly) const
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

	if (ISSETFLAG(m_flags, m_flagOnlyAtBegin))
		return gen == 0;

	if (ISSETFLAG(m_flags, m_flagOnlyAtEnd) && end > 0)
		return gen == end;

	// at Gen has higher priority.
	if (!m_atGen.empty()) {
		// chech atGen.
		for (size_t i = 0, iEnd = m_atGen.size(); i < iEnd; ++i) {
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
	if (end < 0) {                                                          // if we do not know ending generation.
		// can not determine starting gen.
		if (m_beginGen < 0 || m_beginGen > gen)
			return false;

		// can only check positive begin + every
		if ( ( ((gen - m_beginGen) % m_stepGen) == 0) && (m_endGen < 0 || m_endGen >= gen))
			return true;
		else
			return false;
	}                                                                                        // know ending generation
	else {
		ssize_t realStartGen = m_beginGen >= 0 ? m_beginGen : m_beginGen + end + 1;
		ssize_t realEndGen = m_endGen >= 0 ? m_endGen : m_endGen + end + 1;

		if (realStartGen > realEndGen)
			return false;

		return gen >= realStartGen && gen <= realEndGen && (gen - realStartGen) % m_stepGen == 0  ;
	}
	return false;
}


void BaseOperator::setFlags()
{
	RESETFLAG(m_flags, m_flagAtAllGen);
	RESETFLAG(m_flags, m_flagOnlyAtBegin);
	RESETFLAG(m_flags, m_flagOnlyAtEnd);
	RESETFLAG(m_flags, m_flagAllSubPops);

	// atGen has higher priority: if it is not empty, use it.
	if (m_atGen.empty()) {
		if (m_beginGen == 0 && m_endGen == 0)
			SETFLAG(m_flags, m_flagOnlyAtBegin);

		if (m_endGen == -1 && m_beginGen == -1)
			SETFLAG(m_flags, m_flagOnlyAtEnd);

		if (m_beginGen == 0 && m_endGen == -1 && m_stepGen == 1)
			SETFLAG(m_flags, m_flagAtAllGen);

	} else if (m_atGen.size() == 1) {
		if (m_stepGen < 1)
			throw IndexError("active generation interval should be at least 1.");

		if (m_atGen[0] == 0)
			SETFLAG(m_flags, m_flagOnlyAtBegin);
		if (m_atGen[0] == -1)
			SETFLAG(m_flags, m_flagOnlyAtEnd);
	}
	if (m_subPops.allAvail())
		SETFLAG(m_flags, m_flagAllSubPops);
}


bool BaseOperator::apply(Population & /* pop */) const
{
	DBG_FAILIF(true, RuntimeError,
		"This operator can only be applied during mating.");
	return true;
}


bool BaseOperator::applyDuringMating(Population & /* pop */, Population & /* offPop */,
                                     RawIndIterator /* offspring */,
                                     Individual * /* dad */, Individual * /* mom */) const
{
	DBG_FAILIF(true, RuntimeError,
		"This operator cannot be applied during mating.");
	return true;
}


string BaseOperator::applicability(bool subPops, bool gen) const
{
	ostringstream desc;

	if (gen) {
		if (ISSETFLAG(m_flags, m_flagAtAllGen))
			desc << "";
		else if (ISSETFLAG(m_flags, m_flagOnlyAtBegin))
			desc << "at generation 0";
		else if (ISSETFLAG(m_flags, m_flagOnlyAtEnd))
			desc << "at ending generation";
		else if (!m_atGen.empty()) {
			if (m_atGen.size() == 1)
				desc << "at generation";
			else
				desc << "at generations";
			for (size_t i = 0; i < m_atGen.size(); ++i) {
				if (i == 0)
					desc << ", ";
				desc << " " << m_atGen[i];
			}
		} else {
			if (m_beginGen != 0)
				desc << "begin at " << m_beginGen << " ";
			if (m_endGen != -1)
				desc << "end at " << m_endGen << " ";
			if (m_stepGen != 1)
				desc << "at interval " << m_stepGen;
		}
	}
	if (subPops) {
		if (m_subPops.allAvail())
			desc << "";
		else {
			if (desc.str().size() != 1)
				desc << ", ";
			desc << "to subpopulations ";
			for (size_t i = 0; i < m_subPops.size(); ++i) {
				vspID sp = m_subPops[i];
				if (i != 0)
					desc << ", ";
				if (sp.isVirtual())
					desc << "(" << sp.subPop() << ", " << sp.virtualSubPop() << ")";
				else
					desc << sp.subPop();
			}
		}
	}
	if (desc.str().empty())
		return desc.str();
	return "(" + desc.str() + ")";
}


bool BaseOperator::applicableToOffspring(const Population & pop, RawIndIterator offspring) const
{
	pairu pp = pop.subPopIndPair(offspring - pop.rawIndBegin());

	subPopList subPops = m_subPops.expandFrom(pop);

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();

	for (; sp != spEnd; ++sp)
		if (static_cast<ULONG>(sp->subPop()) == pp.first &&
		    (!sp->isVirtual() || pop.virtualSplitter()->contains(pop, pp.second, *sp))) {
			return true;
		}
	return false;
}


opList::opList(PyObject * obj) : m_elems(0)
{
	if (!obj)
		return;
	void * op = pyOpPointer(obj);
	// if this is a valid operator
	if (op)
		m_elems.push_back(reinterpret_cast<BaseOperator *>(op)->clone());
	else if (PySequence_Check(obj)) {
		m_elems.resize(PySequence_Size(obj));
		// assign values
		for (size_t i = 0, iEnd = m_elems.size(); i < iEnd; ++i) {
			PyObject * item = PySequence_GetItem(obj, i);
			void * op = pyOpPointer(item);
			DBG_ASSERT(op, ValueError, "Invalid operator passed.");
			m_elems[i] = reinterpret_cast<BaseOperator *>(op)->clone();
			Py_DECREF(item);
		}
	} else {
		// anything like set, dict_values, iterator can be retrieved this way.
		PyObject * iter = PyObject_GetIter(obj);
		DBG_ASSERT(iter, ValueError, "Invalid operator passed.");
		PyObject * item = NULL;
		while (item == PyIter_Next(iter)) {
			void * op = pyOpPointer(item);
			DBG_ASSERT(op, ValueError, "Invalid operator passed.");
			m_elems.push_back(reinterpret_cast<BaseOperator *>(op)->clone());
			Py_DECREF(item);
		}
		Py_DECREF(iter);
	}
}


opList::opList(const vectorop & ops) : m_elems(0)
{
	vectorop::const_iterator it = ops.begin();
	vectorop::const_iterator itEnd = ops.end();

	for (; it != itEnd; ++it)
		m_elems.push_back((*it)->clone());
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


vectori Pause::s_cachedKeys = vectori();

string Pause::describe(bool /* format */) const
{
	string desc = "<simuPOP.Pause> Pause an evolutionary process";

	if (m_stopOnKeyStroke)
		desc += m_stopOnKeyStroke == static_cast<char>(true) ?
		        string(" with any key stroke") : " with key " + string(1, m_stopOnKeyStroke);
	return desc;
}


bool Pause::apply(Population & pop) const
{
	char a;

	if (m_stopOnKeyStroke != static_cast<char>(false)) {
		// stop on any key
		if (m_stopOnKeyStroke == static_cast<char>(true)) {
			// check if key is already pressed
			if (!simuPOP_kbhit())
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

	string popName = (boost::format("pop_%1%_%2%") % pop.gen() % pop.rep()).str();
	if (m_prompt) {
		cerr	<< "Simulation Paused for population " << pop.rep() << "\n"
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


IfElse::IfElse(PyObject * cond, const opList & ifOps, const opList & elseOps,
               const stringFunc & output, int begin, int end, int step, const intList & at,
               const intList & reps, const subPopList & subPops,
               const stringList & infoFields) :
	BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_cond(PyString_Check(cond) ? PyObj_AsString(cond) : string()),
	m_func(PyCallable_Check(cond) ? cond : NULL),
	m_fixedCond(-1), m_ifOps(ifOps), m_elseOps(elseOps)
{
	(void)output;  // avoid warning about unused parameter
	if (!PyString_Check(cond) && !PyCallable_Check(cond)) {
		bool c;
		PyObj_As_Bool(cond, c);
		const_cast<IfElse *>(this)->m_fixedCond = c ? 1 : 0;
	}
}


string IfElse::describe(bool format) const
{
	string desc = "<simuPOP.IfElse>";
	string ifDesc;
	string elseDesc;

	if (!m_ifOps.empty()) {
		ifDesc = "<ul>\n";
		opList::const_iterator it = m_ifOps.begin();
		opList::const_iterator itEnd = m_ifOps.end();
		for (; it != itEnd; ++it)
			ifDesc += "<li>" + (*it)->describe(false) + " " + (*it)->applicability() + "\n";
		ifDesc += "</ul>";
	}
	if (!m_elseOps.empty()) {
		elseDesc = "<ul>\n";
		opList::const_iterator it = m_elseOps.begin();
		opList::const_iterator itEnd = m_elseOps.end();
		for (; it != itEnd; ++it)
			elseDesc += "<li>" + (*it)->describe(false) + " " + (*it)->applicability() + "\n";
		elseDesc += "</ul>";
	}
	if (m_fixedCond != -1)
		desc += " always apply opertors\n" + (m_fixedCond == 1 ? ifDesc : elseDesc);
	else if (m_func.isValid())
		desc += " apply operators\n" + ifDesc + "\n<indent>if function " + m_func.name() + " returns True";
	else {
		desc += " apply operators \n" + ifDesc + "\n<indent>if " + m_cond.expr();
		if (!m_elseOps.empty())
			desc += ", and otherwise apply operators \n" + elseDesc;
	}
	return format ? formatDescription(desc) : desc;
}


bool IfElse::applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
                               Individual * dad, Individual * mom) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;

	bool res = true;

	if (m_fixedCond != -1)
		res = m_fixedCond == 1;
	else if (m_func.isValid()) {
		PyObject * args = PyTuple_New(m_func.numArgs());

		DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

		for (size_t i = 0; i < m_func.numArgs(); ++i) {
			const string & arg = m_func.arg(i);
			if (arg == "pop")
				PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
			else if (arg == "off")
				PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(&*offspring)));
			else if (arg == "dad")
				PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(dad)));
			else if (arg == "mom")
				PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(mom)));
			else {
				DBG_FAILIF(true, ValueError, "Only parameters 'pop', 'off', 'dad', and 'mom' are acceptable in "
					                         "function" + m_func.name());
			}
		}
		res = m_func(PyObj_As_Bool, args);
		Py_XDECREF(args);
	} else {
		m_cond.setLocalDict(pop.dict());
		res = m_cond.valueAsBool();
	}

	if (res && !m_ifOps.empty()) {
		opList::const_iterator it = m_ifOps.begin();
		opList::const_iterator itEnd = m_ifOps.end();
		for (; it != itEnd; ++it) {
			if (!(*it)->isActive(pop.rep(), pop.gen()))
				continue;
			bool ret = (*it)->applyDuringMating(pop, offPop, offspring, dad, mom);
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
			bool ret = (*it)->applyDuringMating(pop, offPop, offspring, dad, mom);
			if (!ret)
				return false;
		}
		return true;
	}
	return true;
}


bool IfElse::apply(Population & pop) const
{
	bool res = true;

	if (m_fixedCond != -1)
		res = m_fixedCond == 1;
	else if (m_func.isValid()) {
		PyObject * args = PyTuple_New(m_func.numArgs());

		DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

		for (size_t i = 0; i < m_func.numArgs(); ++i) {
			const string & arg = m_func.arg(i);
			if (arg == "pop")
				PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
			else {
				DBG_FAILIF(true, ValueError, "Only parameters 'pop', 'off', 'dad', and 'mom' are acceptable in "
					                         "function" + m_func.name());
			}
		}
		res = m_func(PyObj_As_Bool, args);
		Py_XDECREF(args);
	} else {
		m_cond.setLocalDict(pop.dict());
		res = m_cond.valueAsBool();
	}

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


string TerminateIf::describe(bool /* format */) const
{
	return string("<simuPOP.TerminateIf> terminate the evolution of ") +
	       (m_stopAll ? "all populations" : "the current population") +
	       " if expression \"" + m_expr.expr() + "\" is evalated to be True";
}


bool TerminateIf::apply(Population & pop) const
{
	// experssion return true
	m_expr.setLocalDict(pop.dict());

	if (m_expr.valueAsBool()) {
		if (!noOutput()) {
			ostream & out = getOstream(pop.dict());
			out << m_message << pop.gen() << endl;
			closeOstream();
		}
		if (m_stopAll)
			throw StopEvolution(m_message);
		return false;                                            // return false, this replicate will be stopped
	} else
		return true;
}


RevertIf::RevertIf(PyObject * cond, const string & fromPop, const opList & ops,
                   const stringFunc & output, int begin, int end, int step, const intList & at,
                   const intList & reps, const subPopList & subPops,
                   const stringList & infoFields) :
	BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
	m_cond(PyString_Check(cond) ? PyObj_AsString(cond) : string()),
	m_func(PyCallable_Check(cond) ? cond : NULL),
	m_fixedCond(-1), m_fromPop(fromPop), m_ops(ops)
{
	(void)output;  // avoid warning about unused parameter
	if (!PyString_Check(cond) && !PyCallable_Check(cond)) {
		bool c;
		PyObj_As_Bool(cond, c);
		const_cast<RevertIf *>(this)->m_fixedCond = c ? 1 : 0;
	}
}


string RevertIf::describe(bool /* format */) const
{
	return string("<simuPOP.RevertIf> revert the evolution of population to ") +
	       m_fromPop + " if a condition is met";
}


bool RevertIf::apply(Population & pop) const
{
	bool res = true;

	if (m_fixedCond != -1)
		res = m_fixedCond == 1;
	else if (m_func.isValid()) {
		PyObject * args = PyTuple_New(m_func.numArgs());

		DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

		for (size_t i = 0; i < m_func.numArgs(); ++i) {
			const string & arg = m_func.arg(i);
			if (arg == "pop")
				PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
			else {
				DBG_FAILIF(true, ValueError, "Only parameters 'pop', 'off', 'dad', and 'mom' are acceptable in "
					                         "function" + m_func.name());
			}
		}
		res = m_func(PyObj_As_Bool, args);
		Py_XDECREF(args);
	} else {
		m_cond.setLocalDict(pop.dict());
		res = m_cond.valueAsBool();
	}

	if (res) {
		if (!noOutput()) {
			ostream & out = getOstream(pop.dict());
			out << "Revert to population " << m_fromPop << " at " << pop.gen() << endl;
			closeOstream();
		}
		size_t rep = pop.rep();
		pop.load(m_fromPop);
		if (!pop.getVars().hasVar("gen"))
			pop.setGen(0);
		if (!pop.getVars().hasVar("rep"))
			pop.setRep(rep);
		//
		const vectorop & ops = m_ops.elems();
		vectorop::const_iterator it = ops.begin();
		vectorop::const_iterator itEnd = ops.end();
		for (; it != itEnd; ++it) {
			if (!(*it)->isActive(pop.rep(), pop.gen()))
				continue;
			bool ret = (*it)->apply(pop);
			if (!ret)
				return false;
		}
		throw RevertEvolution("");
	} else
		return true;
}


DiscardIf::DiscardIf(PyObject * cond, const string & exposeInd,
                     const stringFunc & output, int begin, int end, int step, const intList & at,
                     const intList & reps, const subPopList & subPops,
                     const stringList & infoFields) :
	BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_cond(PyString_Check(cond) ? PyObj_AsString(cond) : string()),
	m_func(PyCallable_Check(cond) ? cond : NULL),
	m_fixedCond(-1), m_fixedField(), m_exposeInd(exposeInd), m_dict(NULL),
	m_lastValues()
{
	(void)output;  // avoid warning about unused parameter
	if (!PyString_Check(cond) && !PyCallable_Check(cond)) {
		double res;
		if (cond == Py_True)
			res = 1.;
		else if (cond == Py_False)
			res = 0.;
		else {
			PyObj_As_Double(cond, res);
			DBG_FAILIF(res < 0 || res > 1, RuntimeError,
				"A fixed condition should be either True or False, or a float number between 0 and 1");
		}
		const_cast<DiscardIf *>(this)->m_fixedCond = res;
	} else if (PyString_Check(cond)) {
		string res;
		PyObj_As_String(cond, res);
		//
		std::locale loc;
		bool simple_name = true;
		for (size_t i = 0; i < res.length(); ++i)
			if (!std::isalnum(res[i])) {
				simple_name = false;
				break;
			}
		if (simple_name)
			const_cast<DiscardIf *>(this)->m_fixedField = res;
	}
}


string DiscardIf::describe(bool format) const
{
	string desc;

	if (m_fixedCond != -1) {
		if (m_fixedCond == 1)
			desc += "<simuPOP.DiscardIf> discard all individuals";
		else
			desc += "<simuPOP.DiscardIf> keep all individuals";
	} else if (m_func.isValid())
		desc += "<simuPOP.DiscardIf> discard individuals if function " + m_func.name() + " returns True";
	else
		desc += "<simuPOP.DiscardIf> discard individuals if expression " + m_cond.expr();
	desc += applicability();
	return format ? formatDescription(desc) : desc;
}


bool DiscardIf::apply(Population & pop) const
{
	// mark as not remove for everyone
	pop.markIndividuals(vspID(), false);

	subPopList subPops = applicableSubPops(pop);
	for (UINT idx = 0; idx < subPops.size(); ++idx) {
		if (subPops[idx].isVirtual())
			pop.activateVirtualSubPop(subPops[idx]);

		IndIterator it = pop.indIterator(subPops[idx].subPop());
		for (; it.valid(); ++it) {
			double res = 0;
			if (m_fixedCond != -1)
				res = m_fixedCond;
			else if (!m_fixedField.empty())
				res = it->info(m_fixedField);
			else if (m_func.isValid()) {
				PyObject * args = PyTuple_New(m_func.numArgs());

				DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

				for (size_t i = 0; i < m_func.numArgs(); ++i) {
					const string & arg = m_func.arg(i);
					if (arg == "pop")
						PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
					else if (arg == "ind")
						PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(&*it)));
					else {
						DBG_FAILIF(!it->hasInfoField(arg), ValueError,
							"Only parameters 'ind', 'pop', and names of information fields are "
							"acceptable in function " + m_func.name());
						PyTuple_SET_ITEM(args, i, PyFloat_FromDouble(it->info(arg)));
					}
				}
				PyObject * resObj = m_func(args);
				if (resObj == Py_True)
					res = 1.;
				else if (resObj == Py_False)
					res = 0.;
				else {
					PyObj_As_Double(resObj, res);
					DBG_FAILIF(res < 0 || res > 1, RuntimeError,
						"A callback function for operator DiscardIf has to return either True or False, or a float number between 0 and 1");
				}
				Py_XDECREF(args);
			} else {
				if (!m_dict)
					m_dict = PyDict_New();

				vectorstr infos = it->infoFields();

				for (UINT idx = 0; idx < infos.size(); ++idx) {
					string name = infos[idx];
					double val = it->info(idx);
					// if the value is unchanged, do not set new value
					if (m_lastValues.size() <= idx || m_lastValues[idx] != val) {
						if (m_lastValues.size() <= idx)
							m_lastValues.push_back(val);
						else
							m_lastValues[idx] = val;
						PyObject * obj = PyFloat_FromDouble(val);
						int err = PyDict_SetItemString(m_dict, name.c_str(), obj);
						Py_DECREF(obj);
						if (err != 0) {
#ifndef OPTIMIZED
							if (debug(DBG_GENERAL)) {
								PyErr_Print();
								PyErr_Clear();
							}
#endif
							throw RuntimeError("Setting information fields as variables failed");
						}
					}
				}

				if (!m_exposeInd.empty()) {
					PyObject * indObj = pyIndObj(static_cast<void *>(&*it));
					if (indObj == NULL)
						throw SystemError("Could not expose individual pointer. Compiled with the wrong version of SWIG? ");

					// set dictionary variable pop to this object
					PyDict_SetItemString(m_dict, m_exposeInd.c_str(), indObj);
					Py_DECREF(indObj);
				}

				m_cond.setLocalDict(m_dict);
				// evaluate
				PyObject * resObj = m_cond.evaluate();
				if (resObj == NULL)
					throw RuntimeError("Evaluation of expression failed");

				if (resObj == Py_True)
					res = 1.;
				else if (resObj == Py_False)
					res = 0.;
				else {
					PyObj_As_Double(resObj, res);
					DBG_FAILIF(res < 0 || res > 1, RuntimeError,
						"A callback function for operator DiscardIf has to return either True or False, or a float number between 0 and 1");
				}
				Py_XDECREF(resObj);
			}
			it->setMarked(res == 1. || getRNG().randUniform() < res);
		}
		if (subPops[idx].isVirtual())
			pop.deactivateVirtualSubPop(subPops[idx].subPop());
	}
	pop.removeMarkedIndividuals();
	return true;
}


bool DiscardIf::applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
                                  Individual * dad, Individual * mom) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;

	double res = 0;

	if (m_fixedCond != -1)
		res = m_fixedCond;
	else if (!m_fixedField.empty())
		res = offspring->info(m_fixedField);
	else if (m_func.isValid()) {
		PyObject * args = PyTuple_New(m_func.numArgs());

		DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

		for (size_t i = 0; i < m_func.numArgs(); ++i) {
			const string & arg = m_func.arg(i);
			if (arg == "pop")
				PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
			else if (arg == "off")
				PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(&*offspring)));
			else if (arg == "dad")
				PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(dad)));
			else if (arg == "mom")
				PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(mom)));
			else {
				DBG_FAILIF(!offspring->hasInfoField(arg), ValueError,
					"Only parameters 'off', 'dad', 'mom', 'pop', and names of information fields are "
					"acceptable in function " + m_func.name());
				PyTuple_SET_ITEM(args, i, PyFloat_FromDouble(offspring->info(arg)));
			}
		}
		PyObject * resObj = m_func(args);
		if (resObj == Py_True)
			res = 1.;
		else if (resObj == Py_False)
			res = 0.;
		else {
			PyObj_As_Double(resObj, res);
			DBG_FAILIF(res < 0 || res > 1, RuntimeError,
				"A callback function for operator DiscardIf has to return either True or False, or a float number between 0 and 1");
		}
		Py_XDECREF(args);
	} else {
		if (!m_dict)
			m_dict = PyDict_New();

		vectorstr infos = offspring->infoFields();

		for (UINT idx = 0; idx < infos.size(); ++idx) {
			string name = infos[idx];
			double val = offspring->info(idx);
			// if the value is unchanged, do not set new value
			if (m_lastValues.size() <= idx || m_lastValues[idx] != val) {
				if (m_lastValues.size() <= idx)
					m_lastValues.push_back(val);
				else
					m_lastValues[idx] = val;
				PyObject * obj = PyFloat_FromDouble(val);
				int err = PyDict_SetItemString(m_dict, name.c_str(), obj);
				Py_DECREF(obj);
				if (err != 0) {
#ifndef OPTIMIZED
					if (debug(DBG_GENERAL)) {
						PyErr_Print();
						PyErr_Clear();
					}
#endif
					throw RuntimeError("Setting information fields as variables failed");
				}
			}
		}

		if (!m_exposeInd.empty()) {
			PyObject * indObj = pyIndObj(static_cast<void *>(&*offspring));
			if (indObj == NULL)
				throw SystemError("Could not expose individual pointer. Compiled with the wrong version of SWIG? ");

			// set dictionary variable pop to this object
			PyDict_SetItemString(m_dict, m_exposeInd.c_str(), indObj);
			Py_DECREF(indObj);
		}

		m_cond.setLocalDict(m_dict);
		// evaluate
		PyObject * resObj = m_cond.evaluate();
		if (resObj == NULL)
			throw RuntimeError("Evaluation of expression failed");

		if (resObj == Py_True)
			res = 1.;
		else if (resObj == Py_False)
			res = 0.;
		else {
			PyObj_As_Double(resObj, res);
			DBG_FAILIF(res < 0 || res > 1, RuntimeError,
				"A callback function for operator DiscardIf has to return either True or False, or a float number between 0 and 1");
		}
	}
	// discard if result is True
	return !res;
}


string TicToc::describe(bool /* format */) const
{
	return "<simuPOP.TicToc> output performance monitor>" ;
}


bool TicToc::apply(Population & pop) const
{
#ifdef _OPENMP
	if (omp_get_thread_num() != 0)
		return true;
#endif

	if (m_startTime == 0)
		m_startTime = time(NULL);

	m_totalCount++;
	time_t lastTime = m_lastTime;
	if (m_countPerSec == 0) {
		m_lastTime = time(NULL);
	} else {
		if (m_counter == m_countPerSec) {
			m_lastTime = time(NULL);
			DBG_DO(DBG_OPERATOR, cerr << "Skipped clock check at " <<
				    difftime(m_lastTime, m_startTime) << endl);
			m_counter = 0;
		} else {
			m_lastTime += static_cast<time_t>(1. / m_countPerSec);
			++m_counter;
		}
	}

	double overallTime = difftime(m_lastTime, m_startTime);
	if (overallTime > 5 && m_counter == 0) { // over 5s
		m_countPerSec = static_cast<ULONG>(m_totalCount / overallTime);
		DBG_DO(DBG_OPERATOR, cerr	<< m_totalCount << " hits after " << overallTime << " seconds. Using "
			                        << m_countPerSec << " counts per second." << endl);
	}

	if (!noOutput()) {
		ostream & out = getOstream(pop.dict());
		if (lastTime == 0)
			out << "Start stopwatch." << endl;
		else {
			// since last time
			out << "Elapsed time: " << std::fixed << std::setprecision(2)
			    << difftime(m_lastTime, lastTime)
			    << "s\t Overall time: " << overallTime << "s"
			    << std::resetiosflags(std::ios::fixed) << std::setprecision(-1) << endl;
		}
		closeOstream();
	}
	if (m_stopAfter != 0 && overallTime > m_stopAfter)
		return false;
	return true;
}


bool TicToc::applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
                               Individual * /* dad */, Individual * /* mom */) const
{
#ifdef _OPENMP
	if (omp_get_thread_num() != 0)
		return true;
#endif
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;

	if (m_startTime == 0)
		m_startTime = time(NULL);

	m_totalCount++;
	//
	time_t lastTime = m_lastTime;
	if (m_countPerSec == 0) {
		m_lastTime = time(NULL);
	} else {
		if (m_counter == m_countPerSec) {
			m_lastTime = time(NULL);
			DBG_DO(DBG_OPERATOR, cerr << "Skipped clock check at " <<
				    difftime(m_lastTime, m_startTime) << endl);
			m_counter = 0;
		} else {
			m_lastTime += static_cast<time_t>(1. / m_countPerSec);
			++m_counter;
		}
	}

	double overallTime = difftime(m_lastTime, m_startTime);
	if (overallTime > 5 && m_counter == 0) { // over 5s
		m_countPerSec = static_cast<ULONG>(m_totalCount / overallTime);
		DBG_DO(DBG_OPERATOR, cerr	<< m_totalCount << " hits after " << overallTime << " seconds. Using "
			                        << m_countPerSec << " counts per second." << endl);
	}

	if (!noOutput()) {
		ostream & out = getOstream(pop.dict());
		if (lastTime == 0)
			out << "Start stopwatch." << endl;
		else {
			// since last time
			out << "Elapsed time: " << std::fixed << std::setprecision(2)
			    << difftime(m_lastTime, lastTime)
			    << "s\t Overall time: " << overallTime << "s"
			    << std::resetiosflags(std::ios::fixed) << std::setprecision(-1) << endl;
		}
		closeOstream();
	}
	if (m_stopAfter != 0 && overallTime > m_stopAfter)
		throw StopEvolution("TimeOut by operator TicToc");
	return true;
}


PyOperator::PyOperator(PyObject * func, PyObject * param,
                       int begin, int end, int step, const intList & at,
                       const intList & reps, const subPopList & subPops,
                       const stringList & infoFields) :
	BaseOperator(">", begin, end, step, at, reps, subPops, infoFields),
	m_func(func), m_param(param, true)
{
	if (!m_func.isValid())
		throw ValueError("Passed variable is not a callable Python function.");
}


string PyOperator::describe(bool /* format */) const
{
	return "<simuPOP.PyOperator> calling a Python function " + m_func.name();
}


bool PyOperator::apply(Population & pop) const
{
	PyObject * args = PyTuple_New(m_func.numArgs());

	DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

	bool popMode = true;
	for (size_t i = 0; i < m_func.numArgs(); ++i) {
		if (m_func.arg(i) == "ind") {
			popMode = false;
			break;
		}
	}
	// when the operator is applied to the whole population.
	if (popMode) {
		for (size_t i = 0; i < m_func.numArgs(); ++i) {
			const string & arg = m_func.arg(i);
			if (arg == "pop")
				PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
			else if (arg == "param") {
				Py_INCREF(m_param.object());
				PyTuple_SET_ITEM(args, i, m_param.object());
			} else {
				DBG_FAILIF(true, ValueError,
					"Only parameters 'pop' and 'param' are acceptable in function " + m_func.name());
			}
		}

		PyObject * res = m_func(args);
		Py_XDECREF(args);
		DBG_FAILIF(res != Py_True && res != Py_False, RuntimeError,
			"A callback function for operator PyOperator has to return either True or False");
		return res == Py_True;
	}
	//
	// apply to qualified individual
	pop.markIndividuals(vspID(), false);

	subPopList subPops = applicableSubPops(pop);
	for (UINT idx = 0; idx < subPops.size(); ++idx) {
		if (subPops[idx].isVirtual())
			pop.activateVirtualSubPop(subPops[idx]);

		IndIterator it = pop.indIterator(subPops[idx].subPop());
		for (; it.valid(); ++it) {
			PyObject * args = PyTuple_New(m_func.numArgs());

			DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

			for (size_t i = 0; i < m_func.numArgs(); ++i) {
				const string & arg = m_func.arg(i);
				if (arg == "ind")
					PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(&*it)));
				else if (arg == "param") {
					Py_INCREF(m_param.object());
					PyTuple_SET_ITEM(args, i, m_param.object());
				} else if (arg == "pop")
					PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
				else {
					DBG_FAILIF(true, ValueError,
						"Only parameters 'ind', 'pop', and 'param' are "
						"acceptable in function " + m_func.name());
				}
			}
			PyObject * resObj = m_func(args);
			DBG_FAILIF(resObj != Py_True && resObj != Py_False, RuntimeError,
				"A callback function for operator DiscardIf has to return either True or False");
			Py_XDECREF(args);
			it->setMarked(resObj == Py_False);
		}
		if (subPops[idx].isVirtual())
			pop.deactivateVirtualSubPop(subPops[idx].subPop());
	}
	pop.removeMarkedIndividuals();
	return true;
}


bool PyOperator::applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
                                   Individual * dad, Individual * mom) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;

	PyObject * args = PyTuple_New(m_func.numArgs());

	DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

	for (size_t i = 0; i < m_func.numArgs(); ++i) {
		const string & arg = m_func.arg(i);
		if (arg == "pop")
			PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
		else if (arg == "off" || arg == "ind")
			PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(&*offspring)));
		else if (arg == "dad")
			PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(dad)));
		else if (arg == "mom")
			PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(mom)));
		else if (arg == "param") {
			Py_INCREF(m_param.object());
			PyTuple_SET_ITEM(args, i, m_param.object());
		} else {
			DBG_FAILIF(true, ValueError, "Only parameters 'pop', 'off', 'ind', 'dad', 'mom', and 'param' "
				                         "are acceptable in function" + m_func.name());
		}
	}

	PyObject * res = m_func(args);
	Py_XDECREF(args);
	DBG_FAILIF(res != Py_True && res != Py_False, RuntimeError,
		"A callback function for operator PyOperator has to return either True or False");
	return res == Py_True;
}


void applyDuringMatingOperator(const BaseOperator & op,
                               Population * pop, Population * offPop, ssize_t dad, ssize_t mom,
                               const pairu & off)
{
	BaseOperator * opPtr = op.clone();

	opPtr->initializeIfNeeded(*pop->rawIndBegin());
	// NOTE: putting d and m in the applyDuringMating function will cause
	// compiling error gor gcc/llvm under mac.
	Individual * d = dad < 0 ? NULL : &pop->individual(static_cast<size_t>(dad));
	Individual * m = mom < 0 ? NULL : &pop->individual(static_cast<size_t>(mom));
#ifndef MUTANTALLELE
#  pragma omp parallel for
#endif
	// i needs to be int since some openMP implementation does not handle unsigned index
	for (int i = static_cast<int>(off.first); i < static_cast<int>(off.second); ++i)
		opPtr->applyDuringMating(*pop, *offPop, pop->rawIndBegin() + i, d, m);
}


}
