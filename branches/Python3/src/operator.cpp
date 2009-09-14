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
	if ( (m_rep >= 0 && static_cast<UINT>(m_rep) != rep) ||
	    (m_rep == REP_LAST && rep != numRep - 1) )
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


bool pause::apply(population & pop)
{
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


bool ifElse::applyWithScratch(population & pop, population & scratch, int stage)
{
	m_cond.setLocalDict(pop.dict());
	bool res = m_cond.valueAsBool();

	if (res && m_ifOp != NULL)
		return m_ifOp->applyWithScratch(pop, scratch, stage);
	else if (!res && m_elseOp != NULL)
		return m_elseOp->applyWithScratch(pop, scratch, stage);
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
		out << "Elapsed Time: " << int (timeDiff * 100) / 100. ;
		// since beginning
		timeDiff = difftime(time(&tmpTime), m_startTime);
		int h = int (timeDiff / 3600);
		int m = int ((timeDiff - h * 3600) / 60);
		int s = int (timeDiff - h * 3600 - m * 60);
		out << "s  Overall Time: " << std::setw(2) << std::setfill('0') << h
		    << ":" << std::setw(2) << std::setfill('0') << m << ":" << std::setw(2)
		    << std::setfill('0') << s << endl;
		this->closeOstream();
	}

	time(&m_lastTime);
	return true;
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
	if (m_param == NULL) {
		PyCallFunc(m_func, "(O)", popObj, resBool, PyObj_As_Bool);
	} else {
		PyCallFunc2(m_func, "(OO)", popObj, m_param, resBool, PyObj_As_Bool);
	}

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

	PyObject * arglist, * result;
	if (m_passOffspringOnly) {
		// parammeter list, ref count increased
		if (m_param == NULL)
			arglist = Py_BuildValue("(O)", offObj);
		else
			arglist = Py_BuildValue("(OO)", offObj, m_param);

		// we do not need to catch exceptions here,
		// our wrapper will do that
		result = PyEval_CallObject(m_func, arglist);
		Py_DECREF(offObj);
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
		if (m_param == NULL)
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

	if (result == NULL) {
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


bool pyIndOperator::apply(population & pop)
{
	// if loci is given
	vectora alleles;
	PyObject * numArray = NULL;
	UINT pld = pop.ploidy();

	//
	if (!m_loci.empty()) {
		alleles.resize(m_loci.size() * pld);
		numArray = Allele_Vec_As_NumArray(alleles.begin(), alleles.end());
	}
	vectoru infoIdx(infoSize());
	if (!infoIdx.empty()) {
		for (size_t i = 0; i < infoIdx.size(); ++i)
			infoIdx[i] = pop.infoIdx(infoField(i));
	}
	// call the python function, pass the each individual to it.
	// get pop object
	for (IndIterator it = pop.indBegin(); it.valid(); ++it) {
		PyObject * indObj = pyIndObj(static_cast<void *>(& * it));
		// if pop is valid?
		if (indObj == NULL)
			throw SystemError("Could not pass population to the provided function. \n"
				              "Compiled with the wrong version of SWIG?");

		// loci
		if (!m_loci.empty()) {
			for (size_t i = 0, iEnd = m_loci.size(), j = 0; i < iEnd; ++i)
				for (UINT p = 0; p < pld; ++p)
					alleles[j++] = it->allele(m_loci[i], p);
		}
		// parammeter list, ref count increased
		bool resBool;
		vectorf resArray;

		// if infoFields is not given.
		if (infoSize() == 0) {
			// parenthesis is needed since PyCallFuncX are macros.
			if (m_param == NULL) {
				if (m_loci.empty()) {
					PyCallFunc(m_func, "(O)", indObj, resBool, PyObj_As_Bool);
				} else {
					PyCallFunc2(m_func, "(OO)", indObj, numArray, resBool, PyObj_As_Bool);
				}
			} else {
				if (m_loci.empty()) {
					PyCallFunc2(m_func, "(OO)", indObj, m_param, resBool, PyObj_As_Bool);
				} else {
					PyCallFunc3(m_func, "(OOO)", indObj, numArray, m_param, resBool, PyObj_As_Bool);
				}
			}
			if (!resBool)
				return false;
		} else {
			// parenthesis is needed since PyCallFuncX are macros.
			if (m_param == NULL) {
				if (m_loci.empty()) {
					PyCallFunc(m_func, "(O)", indObj, resArray, PyObj_As_Array);
				} else {
					PyCallFunc2(m_func, "(OO)", indObj, numArray, resArray, PyObj_As_Array);
				}
			} else {
				if (m_loci.empty()) {
					PyCallFunc2(m_func, "(OO)", indObj, m_param, resArray, PyObj_As_Array);
				} else {
					PyCallFunc3(m_func, "(OOO)", indObj, numArray, m_param, resArray, PyObj_As_Array);
				}
			}
			DBG_FAILIF(resArray.size() != infoIdx.size(), ValueError,
				"Returned array should have the same size as given information fields");
			for (size_t i = 0; i < infoSize(); ++i)
				it->setInfo(resArray[i], infoIdx[i]);
		}
		Py_DECREF(indObj);
	}
	return true;
}


}