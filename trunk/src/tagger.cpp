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
/// create a \c tagger, default to be always active but no output
tagger::tagger(string output, string outputExpr,
               int begin, int end, int step, vectorl at, int rep, int grp,
               const vectorstr & infoFields) :
	// stage is automatically determined.
	baseOperator(output, outputExpr, DuringMating, begin, end, step, at, rep, grp, infoFields)
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

	if (m_mode == TAG_Paternal)
		id1 = pop.infoIdx(infoField(0));
	else if (m_mode == TAG_Maternal)
		id1 = pop.infoIdx(infoField(1));
	else {
		id1 = pop.infoIdx(infoField(0));
		id2 = pop.infoIdx(infoField(1));
	}

	if (m_mode == TAG_Paternal) {
		if (dad == NULL)
			offspring->setInfo(0, id1);
		else
			offspring->setInfo(dad->info(id1), id1);
	} else if (m_mode == TAG_Maternal) {
		if (mom == NULL)
			offspring->setInfo(0, id1);
		else
			offspring->setInfo(mom->info(id1), id1);
	} else {
		if (dad == NULL)
			offspring->setInfo(0, id1);
		else
			offspring->setInfo(dad->info(id1), id1);
		if (mom == NULL)
			offspring->setInfo(0, id2);
		else
			offspring->setInfo(mom->info(id2), id2);
	}
	// output to a file?
	if (noOutput())
		return true;
	ostream & out = getOstream(pop.dict());
	if (m_mode == TAG_Paternal) {
		if (dad == NULL)
			out << 0 << '\t';
		else
			out << dad->info(id1) << '\t';
	} else if (m_mode == TAG_Maternal) {
		if (mom == NULL)
			out << 0 << '\t';
		else
			out << mom->info(id1) << '\t';
	} else {
		if (dad == NULL)
			out << 0 << '\t';
		else
			out << dad->info(id1) << '\t';
		if (mom == NULL)
			out << 0 << '\t';
		else
			out << mom->info(id2) << '\t';
	}

	closeOstream();
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
		UINT id1 = pop.infoIdx(infoField(0));

		if (dad != NULL)
			offspring->setInfo(dad - & * pop.indBegin(), id1);
		else if (mom != NULL)
			offspring->setInfo(mom - & * pop.indBegin(), id1);
		else
			offspring->setInfo(0, id1);
	} else if (is == 2) {
		UINT id1 = pop.infoIdx(infoField(0));
		UINT id2 = pop.infoIdx(infoField(1));

		if (dad == NULL)
			offspring->setInfo(0, id1);
		else
			offspring->setInfo(dad - & * pop.indBegin(), id1);

		if (mom == NULL)
			offspring->setInfo(0, id2);
		else
			offspring->setInfo(mom - & * pop.indBegin(), id2);
	}
	// output to a file?
	if (noOutput())
		return true;

	// always output two numbers. Because it is possible to have heteroMating
	// with selfing + random mating.
	ostream & out = getOstream(pop.dict());
	if (dad == NULL)
		out << 0 << '\t';
	else
		out << dad - & * pop.indBegin() << '\t';

	if (mom == NULL)
		out << 0 << '\t';
	else
		out << mom - & * pop.indBegin() << '\t';
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
	vectorf res;
	PyCallFunc(m_func, "(O)", Double_Vec_As_NumArray(values.begin(), values.end()),
	    res, PyObj_As_Array);

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
