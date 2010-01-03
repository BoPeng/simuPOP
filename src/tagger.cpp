/**
 *  $File: tagger.cpp $
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

#include "tagger.h"

namespace simuPOP {

// ID starts from 0 to avoid trouble with other programs that
// treats 0 as missing.
ULONG g_indID = 1;


void IdTagger::reset(ULONG startID)
{
	DBG_FAILIF(startID == 0, ValueError, "Individual ID must start from 1.");
	g_indID = startID;
}


string IdTagger::describe(bool format)
{
	return "<simuPOP.IdTagger> assign an unique ID to individuals" ;
}


bool IdTagger::apply(Population & pop)
{
	DBG_DO(DBG_TAGGER, cerr << "Applying IdTagger with current ID " << g_indID << endl);

	UINT idx = pop.infoIdx(infoField(0));

	int curGen = pop.curAncestralGen();
	for (int depth = pop.ancestralGens(); depth >= 0; --depth) {
		pop.useAncestralGen(depth);
		for (ULONG i = 0, iEnd = pop.popSize(); i < iEnd; ++i)
			pop.ind(i).setInfo(g_indID++, idx);
	}
	pop.useAncestralGen(curGen);
	return true;
}


bool IdTagger::applyDuringMating(Population & pop, RawIndIterator offspring,
                                 individual * dad, individual * mom)
{
	UINT idx = pop.infoIdx(infoField(0));

	DBG_FAILIF(dad != NULL && dad->info(idx) >= g_indID, RuntimeError,
		"Paternal ID is larger than or equal to offspring ID (wrong startID?).");
	DBG_FAILIF(mom != NULL && mom->info(idx) >= g_indID, RuntimeError,
		"Matental ID is larger than or equal to offspring ID (wrong startID?).");
	DBG_FAILIF(mom != NULL && dad != NULL && mom->info(idx) == dad->info(idx), RuntimeError,
		"Parental IDs are not unique (forgot InitInfo?)");

	offspring->setInfo(g_indID++, idx);
	return true;
}


bool InheritTagger::applyDuringMating(Population & pop, RawIndIterator offspring,
                                      individual * dad, individual * mom)
{
	UINT sz = infoSize();

	if (sz == 0)
		return true;

	for (size_t i = 0; i < sz; ++i) {
		UINT idx = pop.infoIdx(infoField(i));

		if (m_mode == PATERNAL) {
			DBG_FAILIF(dad == NULL, RuntimeError,
				"Invalid father for paternal inheritance");
			offspring->setInfo(dad->info(idx), idx);
		} else if (m_mode == MATERNAL) {
			DBG_FAILIF(mom == NULL, RuntimeError,
				"Invalid mother for maternal inheritance");
			offspring->setInfo(mom->info(idx), idx);
		} else if (m_mode == MEAN) {
			DBG_FAILIF(mom == NULL || dad == NULL, RuntimeError,
				"Invalid father or mother for average inheritance");
			offspring->setInfo((mom->info(idx) + dad->info(idx)) / 2, idx);
		}  else if (m_mode == MAXIMUM) {
			DBG_FAILIF(mom == NULL || dad == NULL, RuntimeError,
				"Invalid father or mother for maximum inheritance");
			offspring->setInfo(std::max(mom->info(idx), dad->info(idx)), idx);
		}  else if (m_mode == MINIMUM) {
			DBG_FAILIF(mom == NULL || dad == NULL, RuntimeError,
				"Invalid father or mother for minimum inheritance");
			offspring->setInfo(std::min(mom->info(idx), dad->info(idx)), idx);
		}  else if (m_mode == SUMMATION) {
			DBG_FAILIF(mom == NULL || dad == NULL, RuntimeError,
				"Invalid father or mother for summation inheritance");
			offspring->setInfo(mom->info(idx) + dad->info(idx), idx);
		}  else if (m_mode == MULTIPLICATION) {
			DBG_FAILIF(mom == NULL || dad == NULL, RuntimeError,
				"Invalid father or mother for summation inheritance");
			offspring->setInfo(mom->info(idx) * dad->info(idx), idx);
		} else {
			DBG_FAILIF(true, ValueError, "Invalid inheritance mode");
		}
	}
	return true;
}


bool SummaryTagger::applyDuringMating(Population & pop, RawIndIterator offspring,
                                      individual * dad, individual * mom)
{
	DBG_FAILIF(mom == NULL && dad == NULL, RuntimeError,
		"Invalid father and mother for SummaryTagger.");

	UINT sz = infoSize();

	if (m_mode == MEAN) {
		double all = 0;
		UINT cnt = 0;
		for (size_t i = 0; i < sz - 1; ++i) {
			if (dad != 0) {
				all += dad->info(infoField(i));
				cnt += 1;
			}
			if (mom != 0) {
				all += mom->info(infoField(i));
				cnt += 1;
			}
		}
		offspring->setInfo(all / cnt, infoField(sz - 1));
	}  else if (m_mode == MAXIMUM) {
		double dadMax = 0;
		double momMax = 0;
		if (dad != NULL) {
			double dadMax = dad->info(infoField(0));
			for (size_t i = 1; i < sz - 1; ++i)
				dadMax = std::max(dadMax, dad->info(infoField(i)));
		}
		if (mom != NULL) {
			double momMax = mom->info(infoField(0));
			for (size_t i = 1; i < sz - 1; ++i)
				momMax = std::max(momMax, mom->info(infoField(i)));
		}
		double allMax = dad == NULL ? momMax : (mom == NULL ? dadMax : std::max(dadMax, momMax));
		offspring->setInfo(allMax, infoField(sz - 1));
	}  else if (m_mode == MINIMUM) {
		double dadMin = 0;
		double momMin = 0;
		if (dad != NULL) {
			double dadMin = dad->info(infoField(0));
			for (size_t i = 1; i < sz - 1; ++i)
				dadMin = std::min(dadMin, dad->info(infoField(i)));
		}
		if (mom != NULL) {
			double momMin = mom->info(infoField(0));
			for (size_t i = 1; i < sz - 1; ++i)
				momMin = std::min(momMin, mom->info(infoField(i)));
		}
		double allMin = dad == NULL ? momMin : (mom == NULL ? dadMin : std::min(dadMin, momMin));
		offspring->setInfo(allMin, infoField(sz - 1));
	}  else if (m_mode == SUMMATION) {
		double all = 0;
		for (size_t i = 0; i < sz - 1; ++i) {
			if (dad != 0)
				all += dad->info(infoField(i));
			if (mom != 0)
				all += mom->info(infoField(i));
		}
		offspring->setInfo(all, infoField(sz - 1));
	}  else if (m_mode == MULTIPLICATION) {
		double all = 1.;
		for (size_t i = 0; i < sz - 1; ++i) {
			if (dad != 0)
				all *= dad->info(infoField(i));
			if (mom != 0)
				all *= mom->info(infoField(i));
		}
		offspring->setInfo(all, infoField(sz - 1));
	} else {
		DBG_FAILIF(true, ValueError, "Invalid inheritance mode");
	}

	return true;
}


string ParentsTagger::describe(bool format)
{
	if (infoSize() == 1)
		return "<simuPOP.ParentsTagger> record index of parent in the parental generation "
		       " to information field " + infoField(0) + " of each offspring.";
	else
		return "<simuPOP.ParentsTagger> record indexes of parents in the parental generation "
		       " to information fields " + infoField(0) + " and " + infoField(1) + " of each offspring.";
	// avoid warning
	return "";
}


bool ParentsTagger::applyDuringMating(Population & pop, RawIndIterator offspring,
                                      individual * dad, individual * mom)
{
	DBG_FAILIF(mom == NULL && dad == NULL, ValueError,
		"Both parents are invalid");

	// record to one or two information fields
	size_t is = infoSize();
	if (is == 1) {
		if (dad != NULL)
			offspring->setInfo(dad - &*pop.indIterator(), infoField(0));
		else if (mom != NULL)
			offspring->setInfo(mom - &*pop.indIterator(), infoField(0));
	} else if (is == 2) {
		offspring->setInfo(dad == NULL ? -1 : dad - &*pop.indIterator(), infoField(0));
		offspring->setInfo(mom == NULL ? -1 : mom - &*pop.indIterator(), infoField(1));
	}
	return true;
}


string PedigreeTagger::describe(bool format)
{
	return "<simuPOP.PedigreeTagger> record parental IDs (" + infoField(0) + " and "
	       + infoField(1) + ") to field " + m_idField + " of each offspring.";
}


bool PedigreeTagger::applyDuringMating(Population & pop, RawIndIterator offspring,
                                       individual * dad, individual * mom)
{
	DBG_FAILIF(mom == NULL && dad == NULL, ValueError,
		"Both parents are invalid");

	// record to one or two information fields
	size_t is = infoSize();
	if (is == 1) {
		if (dad != NULL)
			offspring->setInfo(dad->info(m_idField), infoField(0));
		else if (mom != NULL)
			offspring->setInfo(mom->info(m_idField), infoField(0));
	} else if (is == 2) {
		offspring->setInfo(dad == NULL ? -1 : dad->info(m_idField), infoField(0));
		offspring->setInfo(mom == NULL ? -1 : mom->info(m_idField), infoField(1));
	}

	if (noOutput())
		return true;

	ostream & out = getOstream(pop.dict());
	out << offspring->info(m_idField);
	if (dad != NULL)
		out << ' ' << dad->info(m_idField);
	if (mom != NULL)
		out << ' ' << mom->info(m_idField);
	out << '\n';
	return true;
}


bool PyTagger::applyDuringMating(Population & pop, RawIndIterator offspring,
                                 individual * dad, individual * mom)
{
	PyObject * args = PyTuple_New(m_func.numArgs());

	DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

	for (int i = 0; i < m_func.numArgs(); ++i) {
		const string & arg = m_func.arg(i);

		PyObject * item = PyTuple_New((dad != NULL) + (mom != NULL));
		int idx = 0;
		if (dad != NULL) {
			PyTuple_SET_ITEM(item, idx, PyFloat_FromDouble(dad->info(arg)));
			++idx;
		}
		if (mom != NULL)
			PyTuple_SET_ITEM(item, idx, PyFloat_FromDouble(mom->info(arg)));

		PyTuple_SET_ITEM(args, i, item);
	}

	//
	vectorf res = m_func(PyObj_As_Array, args);

	DBG_FAILIF(res.size() != static_cast<size_t>(m_func.numArgs()), ValueError,
		"Please return a value for each information field");

	// assign return values to offspring
	for (size_t i = 0; i < res.size(); ++i)
		offspring->setInfo(res[i], m_func.arg(i));

	Py_DECREF(args);
	return true;
}


}


