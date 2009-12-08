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


void idTagger::reset(ULONG startID)
{
	DBG_FAILIF(startID == 0, ValueError, "Individual ID must start from 1.");
	g_indID = startID;
}


string idTagger::describe(bool format)
{
	return "<simuPOP.idTagger> assign an unique ID to individuals" ;
}


bool idTagger::apply(population & pop)
{
	DBG_DO(DBG_TAGGER, cerr << "Applying idTagger with current ID " << g_indID << endl);

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


bool idTagger::applyDuringMating(population & pop, RawIndIterator offspring,
                                 individual * dad, individual * mom)
{
	UINT idx = pop.infoIdx(infoField(0));

	DBG_FAILIF(dad != NULL && dad->info(idx) >= g_indID, RuntimeError,
		"Paternal ID is larger than or equal to offspring ID (wrong startID?).");
	DBG_FAILIF(mom != NULL && mom->info(idx) >= g_indID, RuntimeError,
		"Matental ID is larger than or equal to offspring ID (wrong startID?).");
	DBG_FAILIF(mom != NULL && dad != NULL && mom->info(idx) == dad->info(idx), RuntimeError,
		"Parental IDs are not unique (forgot initInfo?)");

	offspring->setInfo(g_indID++, idx);
	return true;
}


bool inheritTagger::applyDuringMating(population & pop, RawIndIterator offspring,
                                      individual * dad, individual * mom)
{
	UINT sz = infoSize();

	if (sz == 0)
		return true;

	for (size_t i = 0; i < sz; ++i) {
		UINT idx = pop.infoIdx(infoField(i));

		if (m_mode == Paternal) {
			DBG_FAILIF(dad == NULL, RuntimeError,
				"Invalid father for paternal inheritance");
			offspring->setInfo(dad->info(idx), idx);
		} else if (m_mode == Maternal) {
			DBG_FAILIF(mom == NULL, RuntimeError,
				"Invalid mother for maternal inheritance");
			offspring->setInfo(mom->info(idx), idx);
		} else if (m_mode == Mean) {
			DBG_FAILIF(mom == NULL || dad == NULL, RuntimeError,
				"Invalid father or mother for average inheritance");
			offspring->setInfo((mom->info(idx) + dad->info(idx)) / 2, idx);
		}  else if (m_mode == Maximum) {
			DBG_FAILIF(mom == NULL || dad == NULL, RuntimeError,
				"Invalid father or mother for maximum inheritance");
			offspring->setInfo(std::max(mom->info(idx), dad->info(idx)), idx);
		}  else if (m_mode == Minimum) {
			DBG_FAILIF(mom == NULL || dad == NULL, RuntimeError,
				"Invalid father or mother for minimum inheritance");
			offspring->setInfo(std::min(mom->info(idx), dad->info(idx)), idx);
		}  else if (m_mode == Summation) {
			DBG_FAILIF(mom == NULL || dad == NULL, RuntimeError,
				"Invalid father or mother for summation inheritance");
			offspring->setInfo(mom->info(idx) + dad->info(idx), idx);
		}  else if (m_mode == Multiplication) {
			DBG_FAILIF(mom == NULL || dad == NULL, RuntimeError,
				"Invalid father or mother for summation inheritance");
			offspring->setInfo(mom->info(idx) * dad->info(idx), idx);
		} else {
			DBG_FAILIF(true, ValueError, "Invalid inheritance mode");
		}
	}
	return true;
}


bool summaryTagger::applyDuringMating(population & pop, RawIndIterator offspring,
                                      individual * dad, individual * mom)
{
	DBG_FAILIF(mom == NULL && dad == NULL, RuntimeError,
		"Invalid father and mother for summaryTagger.");

	UINT sz = infoSize();

	if (m_mode == Mean) {
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
	}  else if (m_mode == Maximum) {
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
	}  else if (m_mode == Minimum) {
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
	}  else if (m_mode == Summation) {
		double all = 0;
		for (size_t i = 0; i < sz - 1; ++i) {
			if (dad != 0)
				all += dad->info(infoField(i));
			if (mom != 0)
				all += mom->info(infoField(i));
		}
		offspring->setInfo(all, infoField(sz - 1));
	}  else if (m_mode == Multiplication) {
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


string parentsTagger::describe(bool format)
{
	if (infoSize() == 1)
		return "<simuPOP.parentsTagger> record index of parent in the parental generation "
		       " to information field " + infoField(0) + " of each offspring.";
	else
		return "<simuPOP.parentsTagger> record indexes of parents in the parental generation "
		       " to information fields " + infoField(0) + " and " + infoField(1) + " of each offspring.";
	// avoid warning
	return "";
}


bool parentsTagger::applyDuringMating(population & pop, RawIndIterator offspring,
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


string pedigreeTagger::describe(bool format)
{
	return "<simuPOP.pedigreeTagger> record parental IDs (" + infoField(0) + " and "
	       + infoField(1) + ") to field " + m_idField + " of each offspring.";
}


bool pedigreeTagger::applyDuringMating(population & pop, RawIndIterator offspring,
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


bool pyTagger::applyDuringMating(population & pop, RawIndIterator offspring,
                                 individual * dad, individual * mom)
{
	UINT numFields = infoSize();

	vectoru idx(numFields);

	for (size_t i = 0; i < numFields; ++i)
		idx[i] = pop.infoIdx(infoField(i));

	vectorf values;
	if (dad != NULL) {
		for (size_t i = 0; i < numFields; ++i)
			values.push_back(dad->info(idx[i]));
	}
	if (mom != NULL) {
		for (size_t i = 0; i < numFields; ++i)
			values.push_back(mom->info(idx[i]));
	}
	//
	vectorf res = m_func(PyObj_As_Array, "(O)", Double_Vec_As_NumArray(values.begin(), values.end()));

	DBG_FAILIF(res.size() != numFields, ValueError, "Please return a value for each information field");

	// assign return values to offspring
	for (size_t i = 0; i < numFields; ++i)
		offspring->setInfo(res[i], idx[i]);

	return true;
}


}


