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
		} else if (m_mode == Average) {
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
		} else {
			DBG_FAILIF(true, ValueError, "Invalid inheritance mode");
		}
	}
	return true;
}


bool parentTagger::applyDuringMating(population & pop, RawIndIterator offspring,
                                     individual * dad, individual * mom)
{
	DBG_FAILIF(mom == NULL && dad == NULL, ValueError,
		"Both parents are invalid");

	// record to one or two information fields
	size_t is = infoSize();
	if (is >= 1) {
		if (dad != NULL)
			offspring->setInfo(dad - & * pop.indIterator(), infoField(0));
		else if (mom != NULL)
			offspring->setInfo(mom - & * pop.indIterator(), infoField(0));
	}
	// output to a file?
	if (noOutput())
		return true;

	// output one number
	ostream & out = getOstream(pop.dict());
	ULONG dadIdx = dad == NULL ? 0 : dad - & * pop.indIterator();
	ULONG momIdx = mom == NULL ? 0 : mom - & * pop.indIterator();
	UINT spID = pop.subPopIndPair(std::max(dadIdx, momIdx)).first;
	// record subpopulation count
	if (m_subPopSize.size() < spID + 1)
		m_subPopSize.resize(spID + 1, 0);
	m_subPopSize[spID]++;
	if (dad == NULL)
		out << momIdx << '\t';
	else
		out << dadIdx << '\t';
	closeOstream();
	return true;
}


bool parentTagger::apply(population & pop)
{
	ostream & out = this->getOstream(pop.dict());

	out << "#\t";
	for (size_t i = 0; i < m_subPopSize.size(); ++i)
		out << m_subPopSize[i] << '\t';
	out << '\n';
	closeOstream();
	m_subPopSize.clear();
	return true;
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
			offspring->setInfo(dad - & * pop.indIterator(), infoField(0));
		else if (mom != NULL)
			offspring->setInfo(mom - & * pop.indIterator(), infoField(0));
	} else if (is == 2) {
		offspring->setInfo(dad == NULL ? -1 : dad - & * pop.indIterator(), infoField(0));
		offspring->setInfo(mom == NULL ? -1 : mom - & * pop.indIterator(), infoField(1));
	}
	// output to a file?
	if (noOutput())
		return true;

	// always output two numbers. Because it is possible to have heteroMating
	// with selfing + random mating.
	ostream & out = getOstream(pop.dict());
	ULONG dadIdx = dad == NULL ? -1 : dad - & * pop.indIterator();
	ULONG momIdx = mom == NULL ? -1 : mom - & * pop.indIterator();
	UINT spID = pop.subPopIndPair(std::max(dadIdx, momIdx)).first;
	// record subpopulation count
	if (m_subPopSize.size() < spID + 1)
		m_subPopSize.resize(spID + 1, 0);
	m_subPopSize[spID]++;
	out << dadIdx << '\t';
	out << momIdx << '\t';
	closeOstream();
	return true;
}


bool parentsTagger::apply(population & pop)
{
	ostream & out = this->getOstream(pop.dict());

	out << "#\t";
	for (size_t i = 0; i < m_subPopSize.size(); ++i)
		out << m_subPopSize[i] << '\t';
	out << '\n';
	closeOstream();
	m_subPopSize.clear();
	return true;
}


pedigreeTagger::pedigreeTagger(int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	int stage, const stringFunc & output,
	const stringList & pedigreeFields) :
	baseOperator(output, DuringMating, begin, end, step, at, reps, subPops, pedigreeFields)
{
	setApplicableStage(stage);
}


bool pedigreeTagger::apply(population & pop)
{
	if (infoSize() == 0)
		return true;

	ostream & out = this->getOstream(pop.dict());

	vectoru idx;
	UINT is = infoSize();
	for (size_t i = 0; i < is; ++i)
		idx.push_back(pop.infoIdx(infoField(0)));

	IndIterator it = pop.indIterator();
	for (; it.valid(); ++it)
		for (size_t i = 0; i < is; ++i)
			out << it->info(idx[i]) << '\t';
	out << '\n';
	closeOstream();
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

	// output to a file?
	if (noOutput())
		return true;
	ostream & out = getOstream(pop.dict());
	for (size_t i = 0; i < numFields; ++i)
		out << res[i] << '\t';
	closeOstream();

	return true;
}


}
