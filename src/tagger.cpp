/**
 *  $File: tagger.cpp $
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

#include "tagger.h"
#include "pedigree.h"

#include <sstream>
using std::ostringstream;

namespace simuPOP {

// ID starts from 0 to avoid trouble with other programs that
// treats 0 as missing.
ATOMICLONG g_indID = 1;


void IdTagger::reset(ULONG startID)
{
	DBG_FAILIF(startID == 0, ValueError, "individual ID must start from 1.");
	g_indID = startID;
}


string IdTagger::describe(bool /* format */) const
{
	return "<simuPOP.IdTagger> assign an unique ID to individuals" ;
}


bool IdTagger::apply(Population & pop) const
{
	DBG_DO(DBG_TAGGER, cerr << "Applying IdTagger with current ID " << g_indID << endl);

	size_t idx = pop.infoIdx(infoField(0));

	size_t curGen = pop.curAncestralGen();
	for (int depth = pop.ancestralGens(); depth >= 0; --depth) {
		pop.useAncestralGen(depth);
		for (size_t i = 0, iEnd = pop.popSize(); i < iEnd; ++i)
			pop.individual(i).setInfo(static_cast<double>(g_indID++), idx);
	}
	pop.useAncestralGen(curGen);
	return true;
}


bool IdTagger::applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
                                 Individual * dad, Individual * mom) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;
	size_t idx = pop.infoIdx(infoField(0));

	(void)dad;  // avoid a warning message in optimized modules
	(void)mom;  // avoid a warning message in optimized modules
	DBG_FAILIF(dad != NULL && dad->info(idx) >= g_indID, RuntimeError,
		"Paternal ID is larger than or equal to offspring ID (wrong startID?).");
	DBG_FAILIF(mom != NULL && mom->info(idx) >= g_indID, RuntimeError,
		"Matental ID is larger than or equal to offspring ID (wrong startID?).");
#ifdef _OPENMP
	ATOMICLONG id = fetchAndIncrement(&g_indID);
	offspring->setInfo(static_cast<long>(id), idx);
#else
	offspring->setInfo(static_cast<long>(g_indID++), idx);
#endif
	return true;
}


bool InheritTagger::applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
                                      Individual * dad, Individual * mom) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;
	size_t sz = infoSize();

	if (sz == 0)
		return true;

	for (size_t i = 0; i < sz; ++i) {
		size_t idx = pop.infoIdx(infoField(i));

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


string InheritTagger::describe(bool /* format */) const
{
	ostringstream desc;

	desc << "<simuPOP.InheritTagger> ";
	if (m_mode == PATERNAL)
		desc << "pass information fields " << infoFields().elems() << " from father to offspring";
	else if (m_mode == MATERNAL)
		desc << "pass information fields " << infoFields().elems() << " from mother to offspring";
	else if (m_mode == MEAN)
		desc << "pass the mean value of information fields " << infoFields().elems() << " of parents to offspring";
	else if (m_mode == MAXIMUM)
		desc << "pass the maximum value of information fields " << infoFields().elems() << " of parents to offspring";
	else if (m_mode == MINIMUM)
		desc << "pass the minimum value of information fields " << infoFields().elems() << " of parents to offspring";
	else if (m_mode == SUMMATION)
		desc << "pass the sum of information fields " << infoFields().elems() << " of parents to offspring";
	else if (m_mode == MULTIPLICATION)
		desc << "pass the product of information fields " << infoFields().elems() << " of parents to offspring";
	return desc.str();
}


bool SummaryTagger::applyDuringMating(Population & /* pop */, Population & offPop, RawIndIterator offspring,
                                      Individual * dad, Individual * mom) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;
	DBG_FAILIF(mom == NULL && dad == NULL, RuntimeError,
		"Invalid father and mother for SummaryTagger.");

	size_t sz = infoSize();

	if (m_mode == MEAN) {
		double all = 0;
		size_t cnt = 0;
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


string ParentsTagger::describe(bool /* format */) const
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


bool ParentsTagger::applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
                                      Individual * dad, Individual * mom) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;
	DBG_FAILIF(mom == NULL && dad == NULL, ValueError,
		"Both parents are invalid");

	// record to one or two information fields
	size_t is = infoSize();
	if (is == 1) {
		if (dad != NULL)
			offspring->setInfo(static_cast<double>(dad - &*pop.indIterator()), infoField(0));
		else if (mom != NULL)
			offspring->setInfo(static_cast<double>(mom - &*pop.indIterator()), infoField(0));
	} else if (is == 2) {
		offspring->setInfo(static_cast<double>(dad == NULL ? -1 : dad - &*pop.indIterator()), infoField(0));
		offspring->setInfo(static_cast<double>(mom == NULL ? -1 : mom - &*pop.indIterator()), infoField(1));
	}
	return true;
}


string OffspringTagger::describe(bool /* format */) const
{
	return "<simuPOP.OffspringTagger> records indexes of offspring within family in the offspring population "
	       " to information fields " + infoField(0) + " of each offspring.";
}


bool OffspringTagger::applyDuringMating(Population &, Population & offPop, RawIndIterator offspring,
                                        Individual *, Individual *) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;

	// record to one or two information fields
	if (offspring->firstOffspring())
		offspring->setInfo(0, infoField(0));
	else
		offspring->setInfo((*(offspring - 1)).intInfo(infoField(0)) + 1, infoField(0));
	return true;
}


string PedigreeTagger::describe(bool /* format */) const
{
	return "<simuPOP.PedigreeTagger> record parental IDs (" + infoField(0) + " and "
	       + infoField(1) + ") to field " + m_idField + " of each offspring.";
}


void PedigreeTagger::outputIndividual(ostream & out, const Individual * ind,
                                      const vectorf & IDs) const
{
	// out << .... is very slow compared to the sprintf implementation.
	//
	// three numbers (maximum 20 charameters) + M F, the buffer should be long enough
	char buffer[96];
	char sexChar = ind->sex() == MALE ? 'M' : 'F';
	char affChar = ind->affected() ? 'A' : 'U';

	if (IDs.empty())
		sprintf(buffer, SIZE_T_FORMAT " %c %c", toID(ind->info(m_idField)), sexChar, affChar);
	else if (IDs.size() == 1)
		sprintf(buffer, SIZE_T_FORMAT " " SIZE_T_FORMAT " %c %c", (unsigned long)toID(ind->info(m_idField)), toID(IDs[0]), sexChar, affChar);
	else
		sprintf(buffer, SIZE_T_FORMAT " " SIZE_T_FORMAT " " SIZE_T_FORMAT " %c %c", toID(ind->info(m_idField)), toID(IDs[0]), toID(IDs[1]), sexChar, affChar);
	out << buffer;
	// it is difficult to create buffers for the following, but we do not really care
	// because writing information fields and genotype is rare.
	if (m_outputFields.allAvail())
		for (size_t i = 0; i < ind->infoSize(); ++i)
			out << ' ' << ind->info(i);
	else if (!m_outputFields.elems().empty()) {
		const vectorstr & fields = m_outputFields.elems();
		for (size_t i = 0; i < fields.size(); ++i)
			out << ' ' << ind->info(fields[i]);
	}
	if (m_outputLoci.allAvail()) {
		size_t pldy = ind->ploidy();
		for (size_t i = 0; i < ind->totNumLoci(); ++i)
			for (size_t p = 0; p < pldy; ++p)
				out << ' ' << ind->allele(i, p);
	} else if (!m_outputLoci.elems().empty()) {
		size_t pldy = ind->ploidy();
		const vectoru & loci = m_outputLoci.elems();
		for (size_t i = 0; i < loci.size(); ++i)
			for (size_t p = 0; p < pldy; ++p)
				out << ' ' << ind->allele(loci[i], p);
	}
	out << '\n';
}


bool PedigreeTagger::apply(Population & pop) const
{
	if (noOutput())
		return true;

	//an ID map
	std::map<size_t, int> idMap;

	ostream & out = getOstream(pop.dict());
	size_t is = infoSize();
	vectorf IDs(is);
	vectoru idx(is);
	for (size_t i = 0; i < infoSize(); ++i)
		idx[i] = pop.infoIdx(infoField(i));

	size_t idIdx = pop.infoIdx(m_idField);
	size_t curGen = pop.curAncestralGen();
	for (int depth = pop.ancestralGens(); depth >= 0; --depth) {
		pop.useAncestralGen(depth);
		ConstRawIndIterator it = pop.rawIndBegin();
		ConstRawIndIterator it_end = pop.rawIndEnd();
		for (; it != it_end; ++it) {
			size_t myID = toID(it->info(idIdx));
			idMap[myID] = 1;
			for (size_t i = 0; i < is; ++i) {
				IDs[i] = it->info(idx[i]);
				if (idMap.find(toID(IDs[i])) == idMap.end())
					IDs[i] = 0;
			}
			outputIndividual(out, &*it, IDs);
		}
	}
	pop.useAncestralGen(curGen);
	return true;
}


bool PedigreeTagger::applyDuringMating(Population & pop, Population & offPop, RawIndIterator offspring,
                                       Individual * dad, Individual * mom) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;
	DBG_FAILIF(mom == NULL && dad == NULL, ValueError,
		"Both parents are invalid");

	size_t idIdx = pop.infoIdx(m_idField);
	// record to one or two information fields
	size_t is = infoSize();
	vectorf IDs(is);
	if (is == 1) {
		if (dad != NULL)
			IDs[0] = dad->info(idIdx);
		else if (mom != NULL)
			IDs[0] = mom->info(idIdx);
		offspring->setInfo(IDs[0], pop.infoIdx(infoField(0)));
	} else if (is == 2) {
		IDs[0] = dad == NULL ? 0 : dad->info(idIdx);
		IDs[1] = mom == NULL ? 0 : mom->info(idIdx);
		offspring->setInfo(IDs[0], pop.infoIdx(infoField(0)));
		offspring->setInfo(IDs[1], pop.infoIdx(infoField(1)));
	}

	if (noOutput())
		return true;

	ostream & out = getOstream(pop.dict());
	outputIndividual(out, &*offspring, IDs);
	return true;
}


bool PyTagger::applyDuringMating(Population & /* pop */, Population & offPop, RawIndIterator offspring,
                                 Individual * dad, Individual * mom) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;
	PyObject * args = PyTuple_New(m_func.numArgs());

	DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

	for (size_t i = 0; i < m_func.numArgs(); ++i) {
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


