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

#include "tagger.h"

namespace simuPOP {
// create a \c tagger, default to be always active but no output
tagger::tagger(string output, string outputExpr, int stage,
	int begin, int end, int step, vectorl at, const repList & rep, const subPopList & subPop,
	const vectorstr & infoFields) :
	// stage is automatically determined.
	baseOperator(output, outputExpr, stage, begin, end, step, at, rep, subPop, infoFields)
{
	if (!noOutput())
		setApplicableStage(DuringPostMating);
};


bool tagger::apply(population & pop)
{
	ostream & out = this->getOstream(pop.dict());

	out << '\n';
	closeOstream();
	return true;
}


bool inheritTagger::applyDuringMating(population & pop, RawIndIterator offspring,
                                      individual * dad, individual * mom)
{
	UINT id1 = 0, id2 = 0;

	if (m_mode == TAG_Paternal || m_mode == TAG_Maternal)
		id1 = pop.infoIdx(infoField(0));
	else {
		id1 = pop.infoIdx(infoField(0));
		id2 = pop.infoIdx(infoField(1));
	}

	if (m_mode == TAG_Paternal)
		offspring->setInfo(dad == NULL ? 0 : dad->info(id1), id1);
	else if (m_mode == TAG_Maternal)
		offspring->setInfo(mom == NULL ? 0 : mom->info(id1), id1);
	else {
		offspring->setInfo(dad == NULL ? 0 : dad->info(id1), id1);
		offspring->setInfo(mom == NULL ? 0 : mom->info(id2), id2);
	}
	// output to a file?
	if (noOutput())
		return true;
	ostream & out = getOstream(pop.dict());
	if (m_mode == TAG_Paternal || m_mode == TAG_Maternal)
		out << offspring->info(id1) << '\t';
	else
		out << offspring->info(id1) << '\t' << offspring->info(id2) << '\t';

	closeOstream();
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
			offspring->setInfo(dad - & * pop.indBegin(), infoField(0));
		else if (mom != NULL)
			offspring->setInfo(mom - & * pop.indBegin(), infoField(0));
	}
	// output to a file?
	if (noOutput())
		return true;

	// output one number
	ostream & out = getOstream(pop.dict());
	ULONG dadIdx = dad == NULL ? 0 : dad - & * pop.indBegin();
	ULONG momIdx = mom == NULL ? 0 : mom - & * pop.indBegin();
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
			offspring->setInfo(dad - & * pop.indBegin(), infoField(0));
		else if (mom != NULL)
			offspring->setInfo(mom - & * pop.indBegin(), infoField(0));
	} else if (is == 2) {
		offspring->setInfo(dad == NULL ? 0 : dad - & * pop.indBegin(), infoField(0));
		offspring->setInfo(mom == NULL ? 0 : mom - & * pop.indBegin(), infoField(1));
	}
	// output to a file?
	if (noOutput())
		return true;

	// always output two numbers. Because it is possible to have heteroMating
	// with selfing + random mating.
	ostream & out = getOstream(pop.dict());
	ULONG dadIdx = dad == NULL ? 0 : dad - & * pop.indBegin();
	ULONG momIdx = mom == NULL ? 0 : mom - & * pop.indBegin();
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


pedigreeTagger::pedigreeTagger(int begin, int end, int step, vectorl at,
	const repList & rep, const subPopList & subPop,
	int stage, string output, string outputExpr,
	const vectorstr & pedigreeFields) :
	tagger(output, outputExpr, DuringMating, begin, end, step, at, rep, subPop, pedigreeFields)
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

	IndIterator it = pop.indBegin();
	IndIterator it_end = pop.indEnd();
	for (; it != it_end; ++it)
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
	vectorf res = m_func.call("(O)", Double_Vec_As_NumArray(values.begin(), values.end()),
		PyObj_As_Array);

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
