/**
 *  $File: stator.cpp $
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

#include "stator.h"

#include <sstream>
using std::ostringstream;

namespace simuPOP {

string pyEval::describe(bool format)
{
	string desc = m_expr.expr();

	if (desc.size() > 50)
		desc = desc.substr(0, 50) + "... ";
	return "<simuPOP.pyEval> evalulate expression \"" + desc + "\"";
}


string pyEval::evaluate(population & pop)
{
	if (!m_exposePop.empty()) {
		PyObject * popObj = pyPopObj(static_cast<void *>(&pop));
		if (popObj == NULL)
			throw SystemError("Could not expose population pointer. Compiled with the wrong version of SWIG? ");

		// set dictionary variable pop to this object
		pop.getVars().setVar(m_exposePop, popObj);
	}

	m_expr.setLocalDict(pop.dict());
	string res = m_expr.valueAsString();
	if (!m_exposePop.empty())
		pop.getVars().removeVar(m_exposePop);
	return res;
}


string pyExec::describe(bool format)
{
	string desc = m_expr.stmts();

	if (desc.size() > 50)
		desc = desc.substr(0, 50) + "... ";
	return "<simuPOP.pyExec> execute statements " + desc;
}


string infoEval::describe(bool format)
{
	return "<simuPOP.infoEval> evaluate expression " + m_expr.expr() + " using information fields as variables.";
}


bool pyEval::apply(population & pop)
{
	string res = evaluate(pop);

	if (!this->noOutput() ) {
		ostream & out = this->getOstream(pop.dict());
		out << res;
		this->closeOstream();
	}
	return true;
}


string infoEval::evalInfo(individual * ind, bool update)
{
	vectorstr infos = ind->infoFields();

	for (UINT idx = 0; idx < infos.size(); ++idx) {
		string name = infos[idx];
		double val = ind->info(idx);
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

	if (!m_exposeInd.empty()) {
		PyObject * indObj = pyIndObj(static_cast<void *>(ind));
		if (indObj == NULL)
			throw SystemError("Could not expose individual pointer. Compiled with the wrong version of SWIG? ");

		// set dictionary variable pop to this object
		PyDict_SetItemString(m_dict, m_exposeInd.c_str(), indObj);
		Py_DECREF(indObj);
	}

	m_expr.setLocalDict(m_dict);
	// evaluate
	string res = m_expr.valueAsString();
	// update
	if (update) {
		for (UINT idx = 0; idx < infos.size(); ++idx) {
			double info = 0;
			string name = infos[idx];
			try {
				PyObject * var = PyDict_GetItemString(m_dict, name.c_str());
				PyObj_As_Double(var, info);
				ind->setInfo(info, idx);
			} catch (...) {
				DBG_WARNING(true, "Failed to update information field " + name +
					" from a dictionary of information fields.");
			}
		}
	}

	//
	for (UINT idx = 0; idx < infos.size(); ++idx) {
		string name = infos[idx];
		int err = PyDict_DelItemString(m_dict, name.c_str());
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
	if (!m_exposeInd.empty())
		PyDict_DelItemString(m_dict, m_exposeInd.c_str());
	return res;
}


bool infoEval::apply(population & pop)
{
	m_dict = m_usePopVars ? pop.dict() : PyDict_New();

	subPopList subPops = applicableSubPops();
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for ( ; sp != spEnd; ++sp) {
		pop.activateVirtualSubPop(*sp);
		IndIterator ind = const_cast<population &>(pop).indIterator(sp->subPop());
		for (; ind.valid(); ++ind) {
			string res = evalInfo(&*ind, false) ;
			if (!this->noOutput() ) {
				ostream & out = this->getOstream(pop.dict());
				out << res;
				this->closeOstream();
			}
		}
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


bool infoEval::applyDuringMating(population & pop, RawIndIterator offspring,
                                 individual * dad, individual * mom)
{
	m_dict = m_usePopVars ? pop.dict() : PyDict_New();

	string res = evalInfo(&*offspring, false);

	if (!this->noOutput() ) {
		ostream & out = this->getOstream(pop.dict());
		out << res;
		this->closeOstream();
	}
	return true;
}


string infoExec::describe(bool format)
{
	return "<simuPOP.infoExec> execute statement " + m_expr.stmts() + " using information fields as variables.";
}


bool infoExec::apply(population & pop)
{
	m_dict = m_usePopVars ? pop.dict() : PyDict_New();

	subPopList subPops = applicableSubPops();
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);

	simpleStmt::OperationType oType = m_simpleStmt.operation();
	string oVar = m_simpleStmt.var();
	double oValue = m_simpleStmt.value();
	UINT oVarIdx = 0;
	if (oType != simpleStmt::NoOperation)
		oVarIdx = pop.infoIdx(oVar);

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for ( ; sp != spEnd; ++sp) {
		pop.activateVirtualSubPop(*sp);
		IndIterator ind = const_cast<population &>(pop).indIterator(sp->subPop());
		for (; ind.valid(); ++ind) {
			switch (m_simpleStmt.operation()) {
			case simpleStmt::NoOperation:
				evalInfo(&*ind, true);
				break;
			case simpleStmt::Assignment:
				ind->setInfo(oValue, oVarIdx);
				break;
			case simpleStmt::Increment:
				ind->setInfo(ind->info(oVarIdx) + oValue, oVarIdx);
				break;
			case simpleStmt::Decrement:
				ind->setInfo(ind->info(oVarIdx) - oValue, oVarIdx);
				break;
			case simpleStmt::MultipliedBy:
				ind->setInfo(ind->info(oVarIdx) * oValue, oVarIdx);
				break;
			case simpleStmt::SetSex:
				ind->setInfo(ind->sex(), oVarIdx);
				break;
			case simpleStmt::SetAffected:
				ind->setInfo(ind->affected(), oVarIdx);
				break;
			case simpleStmt::SetUnaffected:
				ind->setInfo(!ind->affected(), oVarIdx);
				break;
			default:
				throw RuntimeError("Incorrect operation type");
			}
		}
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


bool infoExec::applyDuringMating(population & pop, RawIndIterator offspring,
                                 individual * dad, individual * mom)
{
	m_dict = m_usePopVars ? pop.dict() : PyDict_New();

	evalInfo(&*offspring, true);
	return true;
}


string subPopVar_String(vspID sp, const string & var)
{
	if (sp.isVirtual())
		return "subPop{(" + toStr(sp.subPop()) + "," + toStr(sp.virtualSubPop()) + ")}{\'" + var + "\'}";
	else
		return "subPop{" + toStr(sp.subPop()) + "}{\'" + var + "\'}";
}


string haploKey(const vectori & seq)
{
	string key = "{'" + toStr(seq[0]);

	for (size_t i = 1; i < seq.size(); ++i)
		key += toStr("-") + toStr(seq[i]);

	return key + "'}";
}


stat::stat(
	bool popSize,
	//
	bool numOfMales,
	//
	bool numOfAffected,
	//
	const uintList & alleleFreq,
	//
	const uintList & heteroFreq,
	const uintList & homoFreq,
	//
	const uintList & genoFreq,
	//
	const intMatrix & haploFreq,
	//
	const stringList & sumOfInfo,
	const stringList & meanOfInfo,
	const stringList & varOfInfo,
	const stringList & maxOfInfo,
	const stringList & minOfInfo,
	//
	const intMatrix & LD,
	//
	const uintList & association,
	//
	const uintList & neutrality,
	//
	const uintList & structure,
	//
	const uintList & HWE,
	//
	const stringList & vars,
	const string & suffix,
	// regular parameters
	const stringFunc & output,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops, const stringList & infoFields)
	: baseOperator("", begin, end, step, at, reps, subPops, infoFields),
	// the order of initialization is meaningful since they may depend on each other
	m_popSize(popSize, subPops, vars, suffix),
	m_numOfMales(numOfMales, subPops, vars, suffix),
	m_numOfAffected(numOfAffected, subPops, vars, suffix),
	m_alleleFreq(alleleFreq.elems(), subPops, vars, suffix),
	m_heteroFreq(heteroFreq.elems(), homoFreq.elems(), subPops, vars, suffix),
	m_genoFreq(genoFreq.elems(), subPops, vars, suffix),
	m_haploFreq(haploFreq, subPops, vars, suffix),
	m_info(sumOfInfo.elems(), meanOfInfo.elems(), varOfInfo.elems(), maxOfInfo.elems(), minOfInfo.elems(), subPops, vars, suffix),
	m_LD(LD, subPops, vars, suffix),
	m_association(association.elems(), subPops, vars, suffix),
	m_neutrality(neutrality.elems(), subPops, vars, suffix),
	m_structure(structure.elems(), subPops, vars, suffix),
	m_HWE(HWE.elems(), subPops, vars, suffix)
{
}


string stat::describe(bool format)
{
	string desc = "<simuPOP.stat> Calculate statistics\n<ul>\n";
	vectorstr descs;

	descs.push_back(m_popSize.describe(false));
	descs.push_back(m_numOfMales.describe(false));
	descs.push_back(m_numOfAffected.describe(false));
	descs.push_back(m_alleleFreq.describe(false));
	descs.push_back(m_heteroFreq.describe(false));
	descs.push_back(m_genoFreq.describe(false));
	descs.push_back(m_haploFreq.describe(false));
	descs.push_back(m_info.describe(false));
	descs.push_back(m_LD.describe(false));
	descs.push_back(m_association.describe(false));
	descs.push_back(m_neutrality.describe(false));
	descs.push_back(m_structure.describe(false));
	descs.push_back(m_HWE.describe(false));
	for (size_t i = 0; i < descs.size(); ++i) {
		if (!descs[i].empty())
			desc += "<li>" + descs[i] + "\n";
	}
	desc += "</ul>";
	return desc;
}


bool stat::apply(population & pop)
{
	return m_popSize.apply(pop) &&
	       m_numOfMales.apply(pop) &&
	       m_numOfAffected.apply(pop) &&
	       m_alleleFreq.apply(pop) &&
	       m_heteroFreq.apply(pop) &&
	       m_genoFreq.apply(pop) &&
	       m_haploFreq.apply(pop) &&
	       m_info.apply(pop) &&
	       m_LD.apply(pop) &&
	       m_association.apply(pop) &&
	       m_neutrality.apply(pop) &&
	       m_structure.apply(pop) &&
	       m_HWE.apply(pop);
}


statPopSize::statPopSize(bool popSize, const subPopList & subPops,
	const stringList & vars, const string & suffix)
	: m_isActive(popSize), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = { popSize_String,	  popSize_sp_String,
		                           subPopSize_String, "" };
	const char * defaultVars[] = { popSize_String, subPopSize_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);
}


string statPopSize::describe(bool format)
{
	return m_isActive ? "calculate population size" : "";
}


bool statPopSize::apply(population & pop)
{
	if (!m_isActive)
		return true;

	// popSize = ...
	ULONG popSize = 0;
	vectori spSize;
	// for each (virtual) subpopulation
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	for (; it != itEnd; ++it) {
		ULONG spPopSize = pop.subPopSize(*it);
		popSize += spPopSize;
		spSize.push_back(spPopSize);
		if (m_vars.contains(popSize_sp_String))
			pop.getVars().setIntVar(subPopVar_String(*it, popSize_String) + m_suffix, spPopSize);
	}
	// NOTE: popSize does not have to be the total population size
	if (m_vars.contains(popSize_String))
		pop.getVars().setIntVar(popSize_String + m_suffix, popSize);
	// subPopSize = ...
	if (m_vars.contains(subPopSize_String))
		pop.getVars().setIntVectorVar(subPopSize_String + m_suffix, spSize);
	return true;
}


statNumOfMales::statNumOfMales(bool numOfMales, const subPopList & subPops, const stringList & vars, const string & suffix)
	: m_isActive(numOfMales), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = {
		numOfMales_String,		propOfMales_String,
		numOfFemales_String,	propOfFemales_String,
		numOfMales_sp_String,	propOfMales_sp_String,
		numOfFemales_sp_String, propOfFemales_sp_String,   ""
	};
	const char * defaultVars[] = { numOfMales_String, numOfFemales_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);
}


string statNumOfMales::describe(bool format)
{
	return m_isActive ? "count number of male individuals" : "";
}


bool statNumOfMales::apply(population & pop)
{
	if (!m_isActive)
		return true;

	ULONG allMaleCnt = 0;
	ULONG allFemaleCnt = 0;
	ULONG allTotalCnt = 0;
	// for each subpopulation.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for (; sp != spEnd; ++sp) {
		ULONG maleCnt = 0;
		ULONG femaleCnt = 0;
		ULONG totalCnt = 0;
		pop.activateVirtualSubPop(*sp);

		IndIterator it = pop.indIterator(sp->subPop());
		for (; it.valid(); ++it)
			if (it->sex() == MALE)
				maleCnt++;
			else
				femaleCnt++;

		pop.deactivateVirtualSubPop(sp->subPop());

		totalCnt = maleCnt + femaleCnt;

		if (m_vars.contains(numOfMales_sp_String))
			pop.getVars().setIntVar(subPopVar_String(*sp, numOfMales_String) + m_suffix, maleCnt);
		if (m_vars.contains(propOfMales_sp_String))
			pop.getVars().setDoubleVar(subPopVar_String(*sp, propOfMales_String) + m_suffix,
				totalCnt == 0 ? 0 : static_cast<double>(maleCnt) / totalCnt);
		if (m_vars.contains(numOfFemales_sp_String))
			pop.getVars().setIntVar(subPopVar_String(*sp, numOfFemales_String) + m_suffix, femaleCnt);
		if (m_vars.contains(propOfFemales_sp_String))
			pop.getVars().setDoubleVar(subPopVar_String(*sp, propOfFemales_String) + m_suffix,
				totalCnt == 0 ? 0 : static_cast<double>(femaleCnt) / totalCnt);

		allMaleCnt += maleCnt;
		allFemaleCnt += femaleCnt;
		allTotalCnt += totalCnt;
	}

	// output whole population
	if (m_vars.contains(numOfMales_String))
		pop.getVars().setIntVar(numOfMales_String + m_suffix, allMaleCnt);
	if (m_vars.contains(propOfMales_String))
		pop.getVars().setDoubleVar(propOfMales_String + m_suffix, allTotalCnt == 0 ? 0. : static_cast<double>(allMaleCnt) / allTotalCnt);
	if (m_vars.contains(numOfFemales_String))
		pop.getVars().setIntVar(numOfFemales_String + m_suffix, allFemaleCnt);
	if (m_vars.contains(propOfFemales_String))
		pop.getVars().setDoubleVar(propOfFemales_String + m_suffix, allTotalCnt == 0 ? 0 : static_cast<double>(allFemaleCnt) / allTotalCnt);
	return true;
}


statNumOfAffected::statNumOfAffected(bool numOfAffected, const subPopList & subPops,
	const stringList & vars, const string & suffix)
	: m_isActive(numOfAffected), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = {
		numOfAffected_String,	   propOfAffected_String,
		numOfUnaffected_String,	   propOfUnaffected_String,
		numOfAffected_sp_String,   propOfAffected_sp_String,
		numOfUnaffected_sp_String, propOfUnaffected_sp_String,""
	};
	const char * defaultVars[] = { numOfAffected_String, numOfUnaffected_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);
}


string statNumOfAffected::describe(bool format)
{
	return m_isActive ? "count number of affected individuals" : "";
}


bool statNumOfAffected::apply(population & pop)
{
	if (!m_isActive)
		return true;

	ULONG allAffectedCnt = 0;
	ULONG allUnaffectedCnt = 0;
	ULONG allTotalCnt = 0;
	// for each subpopulation.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for (; sp != spEnd; ++sp) {
		ULONG affectedCnt = 0;
		ULONG unaffectedCnt = 0;
		ULONG totalCnt = 0;
		pop.activateVirtualSubPop(*sp);

		IndIterator it = pop.indIterator(sp->subPop());
		for (; it.valid(); ++it)
			if (it->affected())
				affectedCnt++;
			else
				unaffectedCnt++;

		pop.deactivateVirtualSubPop(sp->subPop());

		totalCnt = affectedCnt + unaffectedCnt;

		if (m_vars.contains(numOfAffected_sp_String))
			pop.getVars().setIntVar(subPopVar_String(*sp, numOfAffected_String) + m_suffix, affectedCnt);
		if (m_vars.contains(propOfAffected_sp_String))
			pop.getVars().setDoubleVar(subPopVar_String(*sp, propOfAffected_String) + m_suffix,
				totalCnt == 0 ? 0 : static_cast<double>(affectedCnt) / totalCnt);
		if (m_vars.contains(numOfUnaffected_sp_String))
			pop.getVars().setIntVar(subPopVar_String(*sp, numOfUnaffected_String) + m_suffix, unaffectedCnt);
		if (m_vars.contains(propOfUnaffected_sp_String))
			pop.getVars().setDoubleVar(subPopVar_String(*sp, propOfUnaffected_String) + m_suffix,
				totalCnt == 0 ? 0 : static_cast<double>(unaffectedCnt) / totalCnt);

		allAffectedCnt += affectedCnt;
		allUnaffectedCnt += unaffectedCnt;
		allTotalCnt += totalCnt;
	}

	// output whole population
	if (m_vars.contains(numOfAffected_String))
		pop.getVars().setIntVar(numOfAffected_String + m_suffix, allAffectedCnt);
	if (m_vars.contains(propOfAffected_String))
		pop.getVars().setDoubleVar(propOfAffected_String + m_suffix,
			allTotalCnt == 0 ? 0. : static_cast<double>(allAffectedCnt) / allTotalCnt);
	if (m_vars.contains(numOfUnaffected_String))
		pop.getVars().setIntVar(numOfUnaffected_String + m_suffix, allUnaffectedCnt);
	if (m_vars.contains(propOfUnaffected_String))
		pop.getVars().setDoubleVar(propOfUnaffected_String + m_suffix,
			allTotalCnt == 0 ? 0 : static_cast<double>(allUnaffectedCnt) / allTotalCnt);

	return true;
}


statAlleleFreq::statAlleleFreq(const vectoru & loci, const subPopList & subPops,
	const stringList & vars, const string & suffix)
	: m_loci(loci), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = {
		AlleleNum_String,	 AlleleFreq_String,
		AlleleNum_sp_String, AlleleFreq_sp_String,""
	};
	const char * defaultVars[] = { AlleleFreq_String, AlleleNum_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);
}


string statAlleleFreq::describe(bool format)
{
	ostringstream desc;

	if (!m_loci.empty())
		desc << "calculate allele frequency at " << (m_loci.size() == 1 ? "locus " : "loci ") << m_loci;
	return desc.str();
}


bool statAlleleFreq::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	DBG_DO(DBG_STATOR, cerr << "Calculated allele frequency for loci " << m_loci << endl);

	// count for all specified subpopulations
	ALLELECNTLIST alleleCnt(m_loci.size());
	vectoru allAllelesCnt(m_loci.size(), 0);
	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	for (; it != itEnd; ++it) {
		if (m_vars.contains(AlleleNum_sp_String))
			pop.getVars().removeVar(subPopVar_String(*it, AlleleNum_String) + m_suffix);
		if (m_vars.contains(AlleleFreq_sp_String))
			pop.getVars().removeVar(subPopVar_String(*it, AlleleFreq_String) + m_suffix);

		pop.activateVirtualSubPop(*it);

		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];

#ifdef LONGALLELE
			intDict alleles;
#else
			vectoru alleles(2, 0);
#endif
			size_t allAlleles = 0;

			// go through all alleles
			IndAlleleIterator a = pop.alleleIterator(loc, it->subPop());
			// use allAllelel here because some marker does not have full number
			// of alleles (e.g. markers on chromosome X and Y).
			for (; a.valid(); ++a) {
#ifndef BINARYALLELE
#  ifndef LONGALLELE
				if (*a >= alleles.size())
					alleles.resize(*a + 1, 0);
#  endif
#endif
				alleles[*a]++;
				allAlleles++;
			}
			// total allele count
#ifdef LONGALLELE
			intDict::iterator cnt = alleles.begin();
			intDict::iterator cntEnd = alleles.end();
			for ( ; cnt != cntEnd; ++cnt)
				alleleCnt[idx][cnt->first] += cnt->second;
#else
			for (size_t i = 0; i < alleles.size(); ++i)
				if (alleles[i] != 0)
					alleleCnt[idx][i] += alleles[i];
#endif
			allAllelesCnt[idx] += allAlleles;
			// output variable.
#ifdef LONGALLELE
			if (m_vars.contains(AlleleNum_sp_String))
				pop.getVars().setIntDefDictVar(subPopVar_String(*it, AlleleNum_String) + m_suffix + "{" + toStr(loc) + "}", alleles);
			if (m_vars.contains(AlleleFreq_sp_String)) {
				intDict::iterator cnt = alleles.begin();
				intDict::iterator cntEnd = alleles.end();
				for ( ; cnt != cntEnd; ++cnt)
					cnt->second /= static_cast<double>(allAlleles);
				pop.getVars().setIntDefDictVar(subPopVar_String(*it, AlleleFreq_String) + m_suffix + "{" + toStr(loc) + "}", alleles);
			}
#else
			if (m_vars.contains(AlleleNum_sp_String)) {
				intDict d;
				for (size_t i = 0; i < alleles.size(); ++i)
					if (alleles[i] != 0)
						d[i] = alleles[i];
				pop.getVars().setIntDefDictVar(subPopVar_String(*it, AlleleNum_String) + m_suffix + "{" + toStr(loc) + "}", d);
			}
			if (m_vars.contains(AlleleFreq_sp_String)) {
				intDict d;
				for (size_t i = 0; i < alleles.size(); ++i)
					if (alleles[i] != 0)
						d[i] = alleles[i] / static_cast<double>(allAlleles);
				pop.getVars().setIntDefDictVar(subPopVar_String(*it, AlleleFreq_String) + m_suffix + "{" + toStr(loc) + "}", d);
			}
#endif
		}
		pop.deactivateVirtualSubPop(it->subPop());
	}

	if (m_vars.contains(AlleleNum_String)) {
		pop.getVars().removeVar(AlleleNum_String + m_suffix);
		for (size_t idx = 0; idx < m_loci.size(); ++idx)
			pop.getVars().setIntDefDictVar(AlleleNum_String + m_suffix + "{" + toStr(m_loci[idx]) + "}",
				alleleCnt[idx]);
	}
	if (m_vars.contains(AlleleFreq_String)) {
		pop.getVars().removeVar(AlleleFreq_String + m_suffix);
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			if (allAllelesCnt[idx] != 0) {
				intDict::iterator cnt = alleleCnt[idx].begin();
				intDict::iterator cntEnd = alleleCnt[idx].end();
				for ( ; cnt != cntEnd; ++cnt)
					cnt->second /= static_cast<double>(allAllelesCnt[idx]);
			}
			pop.getVars().setIntDefDictVar(AlleleFreq_String + m_suffix + "{" + toStr(m_loci[idx]) + "}",
				alleleCnt[idx]);
		}
	}
	return true;
}


statHeteroFreq::statHeteroFreq(const vectoru & heteroFreq, const vectoru & homoFreq,
	const subPopList & subPops, const stringList & vars, const string & suffix)
	: m_loci(heteroFreq), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	// add homofreq to m_loci
	for (size_t i = 0; i < homoFreq.size(); ++i)
		if (find(m_loci.begin(), m_loci.end(), homoFreq[i]) == m_loci.end())
			m_loci.push_back(homoFreq[i]);
	//
	const char * allowedVars[] = {
		HeteroNum_String,	 HeteroFreq_String,
		HeteroNum_sp_String, HeteroFreq_sp_String,
		HomoNum_String,		 HomoFreq_String,
		HomoNum_sp_String,	 HomoFreq_sp_String,
		""
	};

	const char * defaultVars[] = { "" };
	m_vars.obtainFrom(vars, allowedVars, defaultVars);
	// add default variables
	if (m_vars.empty()) {
		if (!heteroFreq.empty())
			m_vars.push_back(HeteroFreq_String);
		if (!homoFreq.empty())
			m_vars.push_back(HomoFreq_String);
	}
}


string statHeteroFreq::describe(bool format)
{
	string desc;

	if (!m_loci.empty())
		desc += "calculate heterozygotes frequency";
	return desc;
}


bool statHeteroFreq::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	DBG_FAILIF(pop.ploidy() != 2, ValueError,
		"Heterozygote frequency can only be calculated for diploid populations.");

	DBG_DO(DBG_STATOR, cerr << "Calculated heterozygote frequency for loci " << m_loci << endl);

	// count for all specified subpopulations
	intDict allHeteroCnt;
	intDict allHomoCnt;

	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	for (; it != itEnd; ++it) {
		pop.activateVirtualSubPop(*it);

		intDict heteroCnt;
		intDict homoCnt;

		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];

#ifndef OPTIMIZED
			int chromType = pop.chromType(pop.chromLocusPair(loc).first);
			DBG_FAILIF(chromType == CHROMOSOME_X || chromType == CHROMOSOME_Y,
				ValueError, "Heterozygosity count for sex chromosomes is not supported.");
#endif
			size_t hetero = 0;
			size_t homo = 0;

			// go through all alleles
			IndAlleleIterator a = pop.alleleIterator(loc, it->subPop());
			for (; a.valid(); a += 2) {
				if (*a != *(a + 1))
					hetero += 1;
				else
					homo += 1;
			}
			heteroCnt[loc] = hetero;
			homoCnt[loc] = homo;
			//
			allHeteroCnt[loc] += heteroCnt[loc];
			allHomoCnt[loc] += homoCnt[loc];
		}
		pop.deactivateVirtualSubPop(it->subPop());
		// output subpopulation variable?
		if (m_vars.contains(HeteroNum_sp_String))
			pop.getVars().setIntDictVar(subPopVar_String(*it, HeteroNum_String) + m_suffix, heteroCnt);
		if (m_vars.contains(HomoNum_sp_String))
			pop.getVars().setIntDictVar(subPopVar_String(*it, HomoNum_String) + m_suffix, homoCnt);
		if (m_vars.contains(HeteroFreq_sp_String)) {
			intDict freq;
			for (size_t idx = 0; idx < m_loci.size(); ++idx) {
				UINT loc = m_loci[idx];
				double all = heteroCnt[loc] + homoCnt[loc];
				freq[loc] = all == 0. ? 0 : heteroCnt[loc] / all;
			}
			pop.getVars().setIntDictVar(subPopVar_String(*it, HeteroFreq_String) + m_suffix, freq);
		}
		if (m_vars.contains(HomoFreq_sp_String)) {
			intDict freq;
			for (size_t idx = 0; idx < m_loci.size(); ++idx) {
				UINT loc = m_loci[idx];
				double all = heteroCnt[loc] + homoCnt[loc];
				freq[loc] = all == 0. ? 0 : homoCnt[loc] / all;
			}
			pop.getVars().setIntDictVar(subPopVar_String(*it, HomoFreq_String) + m_suffix, freq);
		}
	}
	// for whole population.
	if (m_vars.contains(HeteroNum_String))
		pop.getVars().setIntDictVar(HeteroNum_String + m_suffix, allHeteroCnt);
	if (m_vars.contains(HomoNum_String))
		pop.getVars().setIntDictVar(HomoNum_String + m_suffix, allHomoCnt);
	if (m_vars.contains(HeteroFreq_String)) {
		intDict freq;
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];
			double all = allHeteroCnt[loc] + allHomoCnt[loc];
			freq[loc] = all == 0. ? 0 : allHeteroCnt[loc] / all;
		}
		pop.getVars().setIntDictVar(HeteroFreq_String + m_suffix, freq);
	}
	if (m_vars.contains(HomoFreq_String)) {
		intDict freq;
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];
			double all = allHeteroCnt[loc] + allHomoCnt[loc];
			freq[loc] = all == 0. ? 0 : allHomoCnt[loc] / all;
		}
		pop.getVars().setIntDictVar(HomoFreq_String + m_suffix, freq);
	}

	return true;
}


statGenoFreq::statGenoFreq(const vectoru & genoFreq, const subPopList & subPops,
	const stringList & vars, const string & suffix)
	: m_loci(genoFreq), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = {
		GenotypeNum_String,	   GenotypeFreq_String,
		GenotypeNum_sp_String, GenotypeFreq_sp_String,""
	};
	const char * defaultVars[] = { GenotypeFreq_String, GenotypeNum_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);
}


string statGenoFreq::describe(bool format)
{
	string desc;

	if (!m_loci.empty())
		desc += "calculate genotype frequency";
	return desc;
}


bool statGenoFreq::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	vectoru chromTypes;
	for (size_t i = 0; i < m_loci.size(); ++i)
		chromTypes.push_back(pop.chromType(pop.chromLocusPair(m_loci[i]).first));

	DBG_DO(DBG_STATOR, cerr << "Calculated genotype frequency for loci " << m_loci << endl);

	// count for all specified subpopulations
	vector<tupleDict> genotypeCnt(m_loci.size());
	vectoru allGenotypeCnt(m_loci.size(), 0);
	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	UINT ply = pop.ploidy();
	for (; it != itEnd; ++it) {
		if (m_vars.contains(GenotypeNum_sp_String))
			pop.getVars().removeVar(subPopVar_String(*it, GenotypeNum_String) + m_suffix);
		if (m_vars.contains(GenotypeFreq_sp_String))
			pop.getVars().removeVar(subPopVar_String(*it, GenotypeFreq_String) + m_suffix);

		pop.activateVirtualSubPop(*it);

		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];

			tupleDict genotypes;
			size_t allGenotypes = 0;

			// go through all alleles
			IndIterator ind = pop.indIterator(it->subPop());
			// the simple case, the speed is potentially faster
			if (!pop.isHaplodiploid() && (chromTypes[idx] == AUTOSOME || chromTypes[idx] == CUSTOMIZED)) {
				for (; ind.valid(); ++ind) {
					vectori genotype(ply);
					for (size_t p = 0; p < ply; ++p)
						genotype[p] = ind->allele(loc, p);
					genotypes[genotype]++;
					allGenotypes++;
				}
			} else {
				for (; ind.valid(); ++ind) {
					vectori genotype;
					for (size_t p = 0; p < ply; ++p) {
						if (p == 1 && ind->sex() == MALE && pop.isHaplodiploid())
							continue;
						if (chromTypes[idx] == CHROMOSOME_Y && ind->sex() == FEMALE)
							continue;
						if (((chromTypes[idx] == CHROMOSOME_X && p == 1) ||
						     (chromTypes[idx] == CHROMOSOME_Y && p == 0)) && ind->sex() == MALE)
							continue;
						genotype.push_back(ind->allele(loc, p));
					}
					genotypes[genotype]++;
					allGenotypes++;
				}
			}
			// total allele count
			tupleDict::iterator dct = genotypes.begin();
			tupleDict::iterator dctEnd = genotypes.end();
			for (; dct != dctEnd; ++dct)
				genotypeCnt[idx][dct->first] += dct->second;
			allGenotypeCnt[idx] += allGenotypes;
			// output variable.
			if (m_vars.contains(GenotypeNum_sp_String))
				pop.getVars().setTupleDefDictVar(subPopVar_String(*it, GenotypeNum_String)
					+ m_suffix + "{" + toStr(loc) + "}", genotypes);
			// note that genotyeps is changed in place.
			if (m_vars.contains(GenotypeFreq_sp_String)) {
				if (allGenotypes != 0) {
					tupleDict::iterator dct = genotypes.begin();
					tupleDict::iterator dctEnd = genotypes.end();
					for (; dct != dctEnd; ++dct)
						dct->second /= allGenotypes;
				}
				pop.getVars().setTupleDefDictVar(subPopVar_String(*it, GenotypeFreq_String)
					+ m_suffix + "{" + toStr(loc) + "}", genotypes);
			}
		}
		pop.deactivateVirtualSubPop(it->subPop());
	}

	if (m_vars.contains(GenotypeNum_String)) {
		pop.getVars().removeVar(GenotypeNum_String + m_suffix);
		for (size_t idx = 0; idx < m_loci.size(); ++idx)
			pop.getVars().setTupleDefDictVar(GenotypeNum_String + m_suffix + "{" + toStr(m_loci[idx]) + "}",
				genotypeCnt[idx]);
	}
	// note that genotyeCnt[idx] is changed in place.
	if (m_vars.contains(GenotypeFreq_String)) {
		pop.getVars().removeVar(GenotypeFreq_String + m_suffix);
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];
			if (allGenotypeCnt[idx] != 0) {
				tupleDict::iterator dct = genotypeCnt[idx].begin();
				tupleDict::iterator dctEnd = genotypeCnt[idx].end();
				for (; dct != dctEnd; ++dct)
					dct->second /= allGenotypeCnt[idx];
			}
			pop.getVars().setTupleDefDictVar(GenotypeFreq_String + m_suffix + "{" + toStr(loc) + "}",
				genotypeCnt[idx]);
		}
	}

	return true;
}


statHaploFreq::statHaploFreq(const intMatrix & haploFreq, const subPopList & subPops,
	const stringList & vars, const string & suffix)
	: m_loci(haploFreq), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = {
		HaplotypeNum_String,	HaplotypeFreq_String,
		HaplotypeNum_sp_String, HaplotypeFreq_sp_String,""
	};
	const char * defaultVars[] = { HaplotypeFreq_String, HaplotypeNum_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);
}


string statHaploFreq::describe(bool format)
{
	string desc;

	if (!m_loci.empty())
		desc += "calculate haplotype frequency";
	return desc;
}


string statHaploFreq::dictKey(const vectori & loci)
{
	string key = "(";

	for (size_t i = 0; i < loci.size(); ++i) {
		if (i != 0)
			key += ",";
		key += toStr(loci[i]);
	}
	key += ")";
	return key;
}


bool statHaploFreq::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	DBG_DO(DBG_STATOR, cerr << "Calculated haplotype frequency for loci " << m_loci << endl);

	// count for all specified subpopulations
	vector<tupleDict> haplotypeCnt(m_loci.size());
	vectoru allHaplotypeCnt(m_loci.size());
	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	UINT ply = pop.ploidy();
	for (; it != itEnd; ++it) {
		if (m_vars.contains(HaplotypeNum_sp_String))
			pop.getVars().removeVar(subPopVar_String(*it, HaplotypeNum_String) + m_suffix);
		if (m_vars.contains(HaplotypeFreq_sp_String))
			pop.getVars().removeVar(subPopVar_String(*it, HaplotypeFreq_String) + m_suffix);

		pop.activateVirtualSubPop(*it);

		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			const vectori & loci = m_loci[idx];
			UINT nLoci = loci.size();
			if (nLoci == 0)
				continue;

			int chromType = pop.chromType(pop.chromLocusPair(loci[0]).first);
#ifndef OPTIMIZED
			for (size_t i = 1; i < nLoci; ++i) {
				DBG_FAILIF(pop.chromType(pop.chromLocusPair(loci[i]).first) != chromType, ValueError,
					"Haplotype must be on the chromosomes of the same type");
			}
#endif
			string key = dictKey(loci);

			tupleDict haplotypes;
			size_t allHaplotypes = 0;

			// go through all individual
			IndIterator ind = pop.indIterator(it->subPop());
			for (; ind.valid(); ++ind) {
				vectori haplotype(loci.size());
				for (size_t p = 0; p < ply; ++p) {
					if (p == 1 && ind->sex() == MALE && pop.isHaplodiploid())
						continue;
					if (chromType == CHROMOSOME_Y && ind->sex() == FEMALE)
						continue;
					if (((chromType == CHROMOSOME_X && p == 1) ||
					     (chromType == CHROMOSOME_Y && p == 0)) && ind->sex() == MALE)
						continue;
					for (size_t idx = 0; idx < nLoci; ++idx)
						haplotype[idx] = ind->allele(loci[idx], p);
					haplotypes[haplotype]++;
					allHaplotypes++;
				}
			}
			// total haplotype count
			tupleDict::iterator dct = haplotypes.begin();
			tupleDict::iterator dctEnd = haplotypes.end();
			for (; dct != dctEnd; ++dct)
				haplotypeCnt[idx][dct->first] += dct->second;
			allHaplotypeCnt[idx] += allHaplotypes;
			// output variable.
			if (m_vars.contains(HaplotypeNum_sp_String))
				pop.getVars().setTupleDefDictVar(subPopVar_String(*it, HaplotypeNum_String) + m_suffix + "{"
					+ key + "}", haplotypes);
			// note that genotyeps is changed in place.
			if (m_vars.contains(HaplotypeFreq_sp_String)) {
				if (allHaplotypes != 0) {
					tupleDict::iterator dct = haplotypes.begin();
					tupleDict::iterator dctEnd = haplotypes.end();
					for (; dct != dctEnd; ++dct)
						dct->second /= allHaplotypes;
				}
				pop.getVars().setTupleDefDictVar(subPopVar_String(*it, HaplotypeFreq_String) + m_suffix + "{"
					+ key + "}", haplotypes);
			}
		}
		pop.deactivateVirtualSubPop(it->subPop());
	}

	if (m_vars.contains(HaplotypeNum_String)) {
		pop.getVars().removeVar(HaplotypeNum_String + m_suffix);
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			string key = dictKey(m_loci[idx]);
			pop.getVars().setTupleDefDictVar(string(HaplotypeNum_String) + m_suffix + "{" + key + "}",
				haplotypeCnt[idx]);
		}
	}
	// note that genotyeCnt[idx] is changed in place.
	if (m_vars.contains(HaplotypeFreq_String)) {
		pop.getVars().removeVar(HaplotypeFreq_String + m_suffix);
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			string key = dictKey(m_loci[idx]);
			if (allHaplotypeCnt[idx] != 0) {
				tupleDict::iterator dct = haplotypeCnt[idx].begin();
				tupleDict::iterator dctEnd = haplotypeCnt[idx].end();
				for (; dct != dctEnd; ++dct)
					dct->second /= allHaplotypeCnt[idx];
			}
			pop.getVars().setTupleDefDictVar(string(HaplotypeFreq_String) + m_suffix + "{" + key + "}",
				haplotypeCnt[idx]);
		}
	}

	return true;
}


statInfo::statInfo(const vectorstr & sumOfInfo, const vectorstr & meanOfInfo,
	const vectorstr & varOfInfo, const vectorstr & maxOfInfo,
	const vectorstr & minOfInfo,
	const subPopList & subPops, const stringList & vars, const string & suffix)
	: m_sumOfInfo(sumOfInfo), m_meanOfInfo(meanOfInfo), m_varOfInfo(varOfInfo),
	m_maxOfInfo(maxOfInfo), m_minOfInfo(minOfInfo),
	m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = {
		SumOfInfo_String,	 MeanOfInfo_String,	   VarOfInfo_String,
		MaxOfInfo_String,	 MinOfInfo_String,
		SumOfInfo_sp_String, MeanOfInfo_sp_String, VarOfInfo_sp_String,
		MaxOfInfo_sp_String, MinOfInfo_sp_String,
		""
	};
	const char * defaultVars[] = { "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);
	if (m_vars.empty()) {
		if (!m_sumOfInfo.empty())
			m_vars.push_back(SumOfInfo_String);
		if (!m_meanOfInfo.empty())
			m_vars.push_back(MeanOfInfo_String);
		if (!m_varOfInfo.empty())
			m_vars.push_back(VarOfInfo_String);
		if (!m_maxOfInfo.empty())
			m_vars.push_back(MaxOfInfo_String);
		if (!m_minOfInfo.empty())
			m_vars.push_back(MinOfInfo_String);
	}
}


string statInfo::describe(bool format)
{
	if (m_sumOfInfo.empty() && m_meanOfInfo.empty() && m_varOfInfo.empty()
	    && m_maxOfInfo.empty() && m_minOfInfo.empty())
		return "";

	string desc = "calculate\n<ul>/n";
	if (!m_sumOfInfo.empty()) {
		desc += "<li>Sum of information fields ";
		for (size_t i = 0; i < m_sumOfInfo.size(); ++i) {
			if (i != 0)
				desc += ", ";
			desc += m_sumOfInfo[i];
		}
		desc += "\n";
	}
	if (!m_meanOfInfo.empty()) {
		desc += "<li>Mean of information fields ";
		for (size_t i = 0; i < m_meanOfInfo.size(); ++i) {
			if (i != 0)
				desc += ", ";
			desc += m_meanOfInfo[i];
		}
		desc += "\n";
	}
	if (!m_varOfInfo.empty()) {
		desc += "<li>Variance of information fields ";
		for (size_t i = 0; i < m_varOfInfo.size(); ++i) {
			if (i != 0)
				desc += ", ";
			desc += m_varOfInfo[i];
		}
		desc += "\n";
	}
	if (!m_maxOfInfo.empty()) {
		desc += "<li>Maximum of information fields ";
		for (size_t i = 0; i < m_maxOfInfo.size(); ++i) {
			if (i != 0)
				desc += ", ";
			desc += m_maxOfInfo[i];
		}
		desc += "\n";
	}
	if (!m_minOfInfo.empty()) {
		desc += "<li>Minimal of information fields ";
		for (size_t i = 0; i < m_minOfInfo.size(); ++i) {
			if (i != 0)
				desc += ", ";
			desc += m_minOfInfo[i];
		}
		desc += "\n";
	}
	desc += "</ul>";
	return desc;
}


bool statInfo::apply(population & pop)
{
	if (m_sumOfInfo.empty() && m_meanOfInfo.empty() && m_varOfInfo.empty()
	    && m_maxOfInfo.empty() && m_minOfInfo.empty())
		return true;

	// field indexes
	UINT numSumFld = m_sumOfInfo.size();
	UINT numMeanFld = m_meanOfInfo.size();
	UINT numVarFld = m_varOfInfo.size();
	UINT numMaxFld = m_maxOfInfo.size();
	UINT numMinFld = m_minOfInfo.size();
	//
	vectoru sumOfInfo(m_sumOfInfo.size());
	vectoru meanOfInfo(m_meanOfInfo.size());
	vectoru varOfInfo(m_varOfInfo.size());
	vectoru maxOfInfo(m_maxOfInfo.size());
	vectoru minOfInfo(m_minOfInfo.size());
	for (size_t i = 0; i < numSumFld; ++i)
		sumOfInfo[i] = pop.infoIdx(m_sumOfInfo[i]);
	for (size_t i = 0; i < numMeanFld; ++i)
		meanOfInfo[i] = pop.infoIdx(m_meanOfInfo[i]);
	for (size_t i = 0; i < numVarFld; ++i)
		varOfInfo[i] = pop.infoIdx(m_varOfInfo[i]);
	for (size_t i = 0; i < numMaxFld; ++i)
		maxOfInfo[i] = pop.infoIdx(m_maxOfInfo[i]);
	for (size_t i = 0; i < numMinFld; ++i)
		minOfInfo[i] = pop.infoIdx(m_minOfInfo[i]);

	vectorf allSumVal(numSumFld);
	vectorf allMeanSumVal(numMeanFld);
	vectoru allMeanNumVal(numMeanFld);
	vectorf allVarSumVal(numVarFld);
	vectorf allVarSum2Val(numVarFld);
	vectoru allVarNumVal(numVarFld);
	vectorf allMaxVal(0);
	vectorf allMinVal(0);
	// for each subpopulation.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for (; sp != spEnd; ++sp) {
		vectorf sumVal(numSumFld, 0.);
		vectorf meanSumVal(numMeanFld, 0.);
		vectoru meanNumVal(numMeanFld, 0);
		vectorf varSumVal(numVarFld, 0.);
		vectorf varSum2Val(numVarFld, 0.);
		vectoru varNumVal(numVarFld, 0);
		vectorf maxVal(0);
		vectorf minVal(0);

		pop.activateVirtualSubPop(*sp);

		IndIterator it = pop.indIterator(sp->subPop());
		for (; it.valid(); ++it) {
			for (size_t i = 0; i < numSumFld; ++i)
				sumVal[i] += it->info(sumOfInfo[i]);
			for (size_t i = 0; i < numMeanFld; ++i) {
				meanSumVal[i] += it->info(meanOfInfo[i]);
				meanNumVal[i]++;
			}
			for (size_t i = 0; i < numVarFld; ++i) {
				double val = it->info(varOfInfo[i]);
				varSumVal[i] += val;
				varSum2Val[i] += val * val;
				varNumVal[i]++;
			}
			if (maxVal.empty()) {
				for (size_t i = 0; i < numMaxFld; ++i)
					maxVal.push_back(it->info(maxOfInfo[i]));
			} else {
				for (size_t i = 0; i < numMaxFld; ++i) {
					if (maxVal[i] < it->info(maxOfInfo[i]))
						maxVal[i] = it->info(maxOfInfo[i]);
				}
			}
			if (minVal.empty()) {
				for (size_t i = 0; i < numMinFld; ++i)
					minVal.push_back(it->info(minOfInfo[i]));
			} else {
				for (size_t i = 0; i < numMinFld; ++i) {
					if (minVal[i] > it->info(minOfInfo[i]))
						minVal[i] = it->info(minOfInfo[i]);
				}
			}
		}

		pop.deactivateVirtualSubPop(sp->subPop());

		for (size_t i = 0; i < numSumFld; ++i)
			allSumVal[i] += sumVal[i];
		for (size_t i = 0; i < numMeanFld; ++i) {
			allMeanSumVal[i] += meanSumVal[i];
			allMeanNumVal[i] += meanNumVal[i];
		}
		for (size_t i = 0; i < numVarFld; ++i) {
			allVarSumVal[i] += varSumVal[i];
			allVarSum2Val[i] += varSum2Val[i];
			allVarNumVal[i] += varNumVal[i];
		}
		if (allMaxVal.empty()) {
			for (size_t i = 0; i < numMaxFld; ++i)
				allMaxVal.push_back(maxVal[i]);
		} else {
			for (size_t i = 0; i < numMaxFld; ++i)
				if (allMaxVal[i] < maxVal[i])
					allMaxVal[i] = maxVal[i];
		}
		if (allMinVal.empty()) {
			for (size_t i = 0; i < numMinFld; ++i)
				allMinVal.push_back(minVal[i]);
		} else {
			for (size_t i = 0; i < numMinFld; ++i)
				if (allMinVal[i] > minVal[i])
					allMinVal[i] = minVal[i];
		}
		// output variable
		if (m_vars.contains(SumOfInfo_sp_String)) {
			strDict dct;
			for (size_t i = 0; i < numSumFld; ++i)
				dct[m_sumOfInfo[i]] = sumVal[i];
			pop.getVars().setStrDictVar(subPopVar_String(*sp, SumOfInfo_String) + m_suffix, dct);
		}
		if (m_vars.contains(MeanOfInfo_sp_String)) {
			strDict dct;
			for (size_t i = 0; i < numMeanFld; ++i)
				dct[m_meanOfInfo[i]] = meanNumVal[i] == 0 ? 0 : meanSumVal[i] / meanNumVal[i];
			pop.getVars().setStrDictVar(subPopVar_String(*sp, MeanOfInfo_String) + m_suffix, dct);
		}
		if (m_vars.contains(VarOfInfo_sp_String)) {
			strDict dct;
			for (size_t i = 0; i < numVarFld; ++i)
				dct[m_varOfInfo[i]] = varNumVal[i] <= 1 ? 0 :
				                      (varSum2Val[i] - varSumVal[i] * varSumVal[i] / varNumVal[i]) / (varNumVal[i] - 1);
			pop.getVars().setStrDictVar(subPopVar_String(*sp, VarOfInfo_String) + m_suffix, dct);
		}
		if (m_vars.contains(MaxOfInfo_sp_String)) {
			strDict dct;
			for (size_t i = 0; i < numMaxFld; ++i)
				dct[m_maxOfInfo[i]] = maxVal[i];
			pop.getVars().setStrDictVar(subPopVar_String(*sp, MaxOfInfo_String) + m_suffix, dct);
		}
		if (m_vars.contains(MinOfInfo_sp_String)) {
			strDict dct;
			for (size_t i = 0; i < numMinFld; ++i)
				dct[m_minOfInfo[i]] = minVal[i];
			pop.getVars().setStrDictVar(subPopVar_String(*sp, MinOfInfo_String) + m_suffix, dct);
		}
	}
	if (m_vars.contains(SumOfInfo_String)) {
		strDict dct;
		for (size_t i = 0; i < m_sumOfInfo.size(); ++i)
			dct[m_sumOfInfo[i]] = allSumVal[i];
		pop.getVars().setStrDictVar(SumOfInfo_String + m_suffix, dct);
	}
	if (m_vars.contains(MeanOfInfo_String)) {
		strDict dct;
		for (size_t i = 0; i < numMeanFld; ++i)
			dct[m_meanOfInfo[i]] = allMeanNumVal[i] == 0 ? 0 : allMeanSumVal[i] / allMeanNumVal[i];
		pop.getVars().setStrDictVar(MeanOfInfo_String + m_suffix, dct);
	}
	if (m_vars.contains(VarOfInfo_String)) {
		strDict dct;
		for (size_t i = 0; i < numVarFld; ++i)
			dct[m_varOfInfo[i]] = allVarNumVal[i] <= 1 ? 0 :
			                      (allVarSum2Val[i] - allVarSumVal[i] * allVarSumVal[i] / allVarNumVal[i]) / (allVarNumVal[i] - 1);
		pop.getVars().setStrDictVar(VarOfInfo_String + m_suffix, dct);
	}
	if (m_vars.contains(MaxOfInfo_String)) {
		strDict dct;
		for (size_t i = 0; i < numMaxFld; ++i)
			dct[m_maxOfInfo[i]] = allMaxVal[i];
		pop.getVars().setStrDictVar(MaxOfInfo_String + m_suffix, dct);
	}
	if (m_vars.contains(MinOfInfo_String)) {
		strDict dct;
		for (size_t i = 0; i < numMinFld; ++i)
			dct[m_minOfInfo[i]] = allMinVal[i];
		pop.getVars().setStrDictVar(MinOfInfo_String + m_suffix, dct);
	}
	return true;
}


statLD::statLD(const intMatrix & LD,  const subPopList & subPops,
	const stringList & vars, const string & suffix)
	: m_LD(LD), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = {
		LD_String,		 LD_prime_String,		R2_String,
		ChiSq_String,	 ChiSq_p_String,		CramerV_String,
		LD_sp_String,	 LD_prime_sp_String,	R2_sp_String,
		ChiSq_sp_String, ChiSq_p_sp_String,		CramerV_sp_String,
		""
	};
	const char * defaultVars[] = { LD_String, LD_prime_String, R2_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);

	for (size_t i = 0; i < m_LD.size(); ++i) {
		DBG_FAILIF(m_LD[i].size() != 2 && m_LD[i].size() != 4, ValueError,
			"Parameter LD should be a list of loci pairs with optional primary alleles.");
	}
}


string statLD::describe(bool format)
{
	string desc;

	if (!m_LD.empty())
		desc += "calculate Linkage disequilibrium";
	return desc;
}


void statLD::calculateLD(const vectoru & lociMap, const ALLELECNTLIST & alleleCnt,
                         const HAPLOCNTLIST & haploCnt, vectorf & LD, vectorf & D_prime,
                         vectorf & R2, vectorf & ChiSq, vectorf & ChiSq_p, vectorf & CramerV)
{
	for (size_t idx = 0; idx < m_LD.size(); ++idx) {
		UINT loc1 = m_LD[idx][0];
		UINT loc2 = m_LD[idx][1];
		const ALLELECNT & alleleCnt1 = alleleCnt[lociMap[loc1]];
		const ALLELECNT & alleleCnt2 = alleleCnt[lociMap[loc2]];
		vectora alleles1;
		vectorf freq1;
		vectora alleles2;
		vectorf freq2;
		ALLELECNT::const_iterator cnt = alleleCnt1.begin();
		ALLELECNT::const_iterator cntEnd = alleleCnt1.end();
		double allAlleles = 0;
		for (; cnt != cntEnd; ++cnt) {
			alleles1.push_back(ToAllele(cnt->first));
			freq1.push_back(cnt->second);
			allAlleles += cnt->second;
		}
		for (size_t i = 0; i < freq1.size(); ++i)
			freq1[i] /= allAlleles;
		cnt = alleleCnt2.begin();
		cntEnd = alleleCnt2.end();
		allAlleles = 0;
		for (; cnt != cntEnd; ++cnt) {
			alleles2.push_back(ToAllele(cnt->first));
			freq2.push_back(cnt->second);
			allAlleles += cnt->second;
		}
		for (size_t i = 0; i < freq2.size(); ++i)
			freq2[i] /= allAlleles;
		DBG_DO(DBG_STATOR, cerr << "Loc " << loc1 << " freq " << freq1 << endl
			                    << "Loc " << loc2 << " freq " << freq2 << endl);
		UINT nAllele1 = alleles1.size();
		UINT nAllele2 = alleles2.size();
		// get haplotype count
		const HAPLOCNT & haplos = haploCnt[idx];
		// get total haplotype count (used to calculate haplotype frequency)
		double allHaplo = 0;
		HAPLOCNT::const_iterator hCnt = haplos.begin();
		HAPLOCNT::const_iterator hCntEnd = haplos.end();
		for (; hCnt != hCntEnd; ++hCnt) {
			allHaplo += hCnt->second;
			DBG_DO(DBG_STATOR, cerr << "Haplotype " << hCnt->first.first << ", " << hCnt->first.second
				                    << " Cnt " << hCnt->second << endl);
		}
		if (allHaplo == 0)
			return;
		//
		// calculate LD
		if (!LD.empty()) {
			// when primary alleles are specified
			if (m_LD[idx].size() == 4) {
				UINT A = m_LD[idx][2];
				UINT B = m_LD[idx][3];
				double P_AB = 0;
				// haplotype might not exist.
				HAPLOCNT::const_iterator h = haplos.find(HAPLOCNT::key_type(A, B));
				if (h != haplos.end())
					P_AB = h->second / allHaplo;
				// get allele freq from the m_alleleFreq object
				double P_A = 0;
				ALLELECNT::const_iterator cnt = alleleCnt1.find(A);
				if (cnt != alleleCnt1.end())
					P_A = cnt->second / allAlleles;
				double P_B = 0;
				cnt = alleleCnt2.find(B);
				if (cnt != alleleCnt2.end())
					P_B = cnt->second / allAlleles;

				// calculate LD
				double D = P_AB - P_A * P_B;
				// calculate LD'
				double D_max = D > 0 ? std::min(P_A * (1 - P_B), (1 - P_A) * P_B) : std::min(P_A * P_B, (1 - P_A) * (1 - P_B));
				// fcmp_eq is the float comparison operator, which treat (-1e-10, 1e-10) or so as 0 (platform dependent)
				double Dp = fcmp_eq(D_max, 0.) ? 0. : D / D_max;
				double r2 = (fcmp_eq(P_A, 0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1)) ? 0. : D * D / P_A / (1 - P_A) / P_B / (1 - P_B);
				// if TurnOnDebug(DBG_STATOR) is called in python, the following will be printed.
				DBG_DO(DBG_STATOR, cerr << "P_AB: " << P_AB
					                    << " P_A: " << P_A << " P_B: " << P_B << " D_max: " << D_max <<
					" LD: " << D << " LD': " << Dp << " r2: " << r2 << endl);

				// LD and D' can be negative in this case
				LD[idx] = D;
				DBG_ASSERT(fcmp_ge(LD[idx], -0.25) && fcmp_le(LD[idx], 0.25), SystemError,
					"Calculated LD value " + toStr(LD[idx]) + " out of range of [-1/4, 1/4]");

				D_prime[idx] = Dp;
				DBG_ASSERT(fcmp_ge(D_prime[idx], -1) && fcmp_le(D_prime[idx], 1), SystemError,
					"Calculated D' value " + toStr(D_prime[idx]) + " out of range of [-1, 1]");

				R2[idx] = r2;
				DBG_ASSERT(fcmp_ge(R2[idx], 0) && fcmp_le(R2[idx], 1), SystemError,
					"Calculated R2 value " + toStr(R2[idx]) + " out of range of [0, 1]");
			} else {
				for (size_t i = 0; i < nAllele1; ++i) {
					for (size_t j = 0; j < nAllele2; ++j) {
						UINT A = alleles1[i];
						UINT B = alleles2[j];
						// no particular haplotype
						double P_AB = 0;
						HAPLOCNT::const_iterator h = haplos.find(HAPLOCNT::key_type(A, B));
						if (h != haplos.end())
							P_AB = h->second / allHaplo;
						// get allele freq from the m_alleleFreq object
						double P_A = freq1[i];
						double P_B = freq2[j];

						// calculate LD
						double D = P_AB - P_A * P_B;
						// calculate LD'
						double D_max = D > 0 ? std::min(P_A * (1 - P_B), (1 - P_A) * P_B) : std::min(P_A * P_B, (1 - P_A) * (1 - P_B));
						// fcmp_eq is the float comparison operator, which treat (-1e-10, 1e-10) or so as 0 (platform dependent)
						double Dp = fcmp_eq(D_max, 0.) ? 0. : D / D_max;
						double r2 = (fcmp_eq(P_A, 0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1)) ? 0. : D * D / P_A / (1 - P_A) / P_B / (1 - P_B);
						// if TurnOnDebug(DBG_STATOR) is called in python, the following will be printed.
						DBG_DO(DBG_STATOR, cerr << "P_AB: " << P_AB
							                    << " P_A: " << P_A << " P_B: " << P_B << " D_max: " << D_max <<
							" LD: " << D << " LD': " << Dp << " r2: " << r2 << endl);

						if (nAllele1 <= 2 && nAllele2 <= 2) {
							// for the monomorphic or diallelic case, there is no need to do an average.
							LD[idx] = fabs(D);
							DBG_ASSERT(fcmp_ge(LD[idx], 0) && fcmp_le(LD[idx], 0.25), SystemError,
								"Calculated LD value " + toStr(LD[idx]) + " out of range of [0, 1/4]");

							D_prime[idx] = fabs(Dp);
							DBG_ASSERT(fcmp_ge(D_prime[idx], 0) && fcmp_le(D_prime[idx], 1), SystemError,
								"Calculated D' value " + toStr(D_prime[idx]) + " out of range of [0, 1]");

							R2[idx] = r2;
							DBG_ASSERT(fcmp_ge(R2[idx], 0) && fcmp_le(R2[idx], 1), SystemError,
								"Calculated R2 value " + toStr(R2[idx]) + " out of range of [0, 1]");
							break;
						} else {
							// for the monomorphic or diallelic case, there is no need to do an average.
							LD[idx] += P_A * P_B * fabs(D);
							DBG_ASSERT(fcmp_ge(LD[idx], 0) && fcmp_le(LD[idx], 0.25), SystemError,
								"Calculated LD value " + toStr(LD[idx]) + " out of range of [0, 1/4]");

							D_prime[idx] += P_A * P_B * fabs(Dp);
							DBG_ASSERT(fcmp_ge(D_prime[idx], 0) && fcmp_le(D_prime[idx], 1), SystemError,
								"Calculated D' value " + toStr(D_prime[idx]) + " out of range of [0, 1]");

							R2[idx] += P_A * P_B * r2;
							DBG_ASSERT(fcmp_ge(R2[idx], 0) && fcmp_le(R2[idx], 1), SystemError,
								"Calculated R2 value " + toStr(R2[idx]) + " out of range of [0, 1]");
						}
					}
					if (nAllele1 <= 2 && nAllele2 <= 2)
						break;
				}
			}       // with/without primary alleles
		}           // LD measures
		if (!ChiSq.empty()) {
			// if ChiSq is empty, do not calculate any association stuff.
			if (m_LD[idx].size() == 4) {
				UINT A = m_LD[idx][2];
				UINT B = m_LD[idx][3];

				// calculate association
				vector<vectoru> cont_table(2);
				for (size_t i = 0; i < 2; ++i)
					cont_table[i].resize(2);
				// get P_ij
				HAPLOCNT::const_iterator hCnt = haplos.begin();
				HAPLOCNT::const_iterator hCntEnd = haplos.end();
				for (; hCnt != hCntEnd; ++hCnt) {
					UINT i = hCnt->first.first != A;
					UINT j = hCnt->first.second != B ;
					cont_table[i][j] += hCnt->second;
				}
				// calculate statistics
				chisqTest(cont_table, ChiSq[idx], ChiSq_p[idx]);
				CramerV[idx] = sqrt(ChiSq[idx] / allHaplo);
			} else {
				// calculate association
				vector<vectoru> cont_table(nAllele1);
				for (size_t i = 0; i < nAllele1; ++i)
					cont_table[i].resize(nAllele2);
				// get P_ij
				for (size_t i = 0; i < nAllele1; ++i)
					for (size_t j = 0; j < nAllele2; ++j) {
						HAPLOCNT::const_iterator it = haplos.find(HAPLOCNT::key_type(alleles1[i], alleles2[j]));
						if (it != haplos.end())
							cont_table[i][j] = it->second;
					}
				// calculate statistics
				chisqTest(cont_table, ChiSq[idx], ChiSq_p[idx]);
				CramerV[idx] = sqrt(ChiSq[idx] / (allHaplo * std::min(nAllele1 - 1, nAllele2 - 1)));
			}   // with/without primary alleles
		}       // association
	}           // for all m_LD
}


void statLD::outputVar(population & pop, const string & name, const vectorf & value)
{
	if (value.empty())
		return;

	DBG_FAILIF(value.size() != m_LD.size(), RuntimeError,
		"Return result has incorrect value");

	map<UINT, intDict> res;
	for (size_t i = 0; i < m_LD.size(); ++i)
		res[m_LD[i][0]][m_LD[i][1]] = value[i];

	pop.getVars().removeVar(name);
	map<UINT, intDict>::const_iterator it = res.begin();
	map<UINT, intDict>::const_iterator itEnd = res.end();
	for (; it != itEnd; ++it)
		pop.getVars().setIntDictVar(name + "{" + toStr(it->first) + "}", it->second);
}


bool statLD::apply(population & pop)
{
	if (m_LD.empty())
		return true;

	UINT nLD = m_LD.size();
	// determine involved LD.
	vectoru loci;
	vectoru lociMap(pop.totNumLoci());
	vectoru chromTypes;
	for (size_t i = 0; i < nLD; ++i) {
		for (size_t j = 0; j < 2; ++j) {
			DBG_FAILIF(m_LD[i][j] < 0 || static_cast<UINT>(m_LD[i][j]) >= pop.totNumLoci(),
				IndexError, "Locus index " + toStr(m_LD[i][j]) + " out of range.");
			if (find(loci.begin(), loci.end(), static_cast<UINT>(m_LD[i][j])) == loci.end()) {
				loci.push_back(m_LD[i][j]);
				lociMap[m_LD[i][j]] = loci.size() - 1;
				chromTypes.push_back(pop.chromType(pop.chromLocusPair(m_LD[i][j]).first));
			}
		}
		DBG_FAILIF(pop.chromType(pop.chromLocusPair(m_LD[i][0]).first) !=
			pop.chromType(pop.chromLocusPair(m_LD[i][1]).first),
			ValueError, "Two loci must be on chromosome(s) of the same type");
	}
	UINT nLoci = loci.size();

	ALLELECNTLIST allAlleleCnt(loci.size());
	HAPLOCNTLIST allHaploCnt(m_LD.size());

	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	UINT ply = pop.ploidy();
	for (; it != itEnd; ++it) {
		const char * spVars[] = {
			LD_sp_String,	 LD_prime_sp_String,	R2_sp_String,
			ChiSq_sp_String, ChiSq_p_sp_String,		CramerV_sp_String,
			""
		};
		for (size_t i = 0; spVars[i][0]; ++i) {
			if (m_vars.contains(spVars[i]))
				pop.getVars().removeVar(subPopVar_String(*it, spVars[i]));
		}

		pop.activateVirtualSubPop(*it);

		ALLELECNTLIST alleleCnt(loci.size());
		HAPLOCNTLIST haploCnt(m_LD.size());

		// count allele and genotype
		IndIterator ind = pop.indIterator(it->subPop());
		for (; ind.valid(); ++ind) {
			for (size_t p = 0; p < ply; ++p) {
				if (ply == 2 && p == 1 && ind->sex() == MALE && pop.isHaplodiploid())
					continue;
				GenoIterator geno = ind->genoBegin(p);
				// allele frequency
				for (size_t idx = 0; idx < nLoci; ++idx) {
					if (ply == 2 && chromTypes[idx] == CHROMOSOME_Y && ind->sex() == FEMALE)
						continue;
					if (ply == 2 && ((chromTypes[idx] == CHROMOSOME_X && p == 1) ||
					                 (chromTypes[idx] == CHROMOSOME_Y && p == 0)) && ind->sex() == MALE)
						continue;
					alleleCnt[idx][*(geno + loci[idx])]++;
				}
				// haplotype frequency
				for (size_t idx = 0; idx < nLD; ++idx) {
					UINT chromType = chromTypes[lociMap[m_LD[idx][0]]];
					if (chromType == CHROMOSOME_Y && ind->sex() == FEMALE)
						continue;
					if (((chromType == CHROMOSOME_X && p == 1) ||
					     (chromType == CHROMOSOME_Y && p == 0)) && ind->sex() == MALE)
						continue;
					haploCnt[idx][HAPLOCNT::key_type(*(geno + m_LD[idx][0]), *(geno + m_LD[idx][1]))]++;
				}
			}
		}
		pop.deactivateVirtualSubPop(it->subPop());
		// add to all count
		for (size_t idx = 0; idx < nLoci; ++idx) {
			ALLELECNT::iterator cnt = alleleCnt[idx].begin();
			ALLELECNT::iterator cntEnd = alleleCnt[idx].end();
			for (; cnt != cntEnd; ++cnt)
				allAlleleCnt[idx][cnt->first] += cnt->second;
		}
		//
		for (size_t idx = 0; idx < m_LD.size(); ++idx) {
			HAPLOCNT::iterator cnt = haploCnt[idx].begin();
			HAPLOCNT::iterator cntEnd = haploCnt[idx].end();
			for (; cnt != cntEnd; ++cnt)
				allHaploCnt[idx][cnt->first] += cnt->second;
		}
		// calculate statistics
		UINT ldSize = 0;
		if (m_vars.contains(LD_sp_String) || m_vars.contains(LD_prime_sp_String) || m_vars.contains(R2_sp_String))
			ldSize = m_LD.size();
		vectorf LD(ldSize);
		vectorf D_prime(ldSize);
		vectorf R2(ldSize);
		UINT assoSize = 0;
		if (m_vars.contains(ChiSq_sp_String) || m_vars.contains(ChiSq_p_sp_String) || m_vars.contains(CramerV_sp_String))
			assoSize = m_LD.size();
		vectorf ChiSq(assoSize);
		vectorf ChiSq_p(assoSize);
		vectorf CramerV(assoSize);
		calculateLD(lociMap, alleleCnt, haploCnt, LD, D_prime, R2,
			ChiSq, ChiSq_p, CramerV);

		// output statistics for subpopulation
		if (m_vars.contains(LD_sp_String))
			outputVar(pop, subPopVar_String(*it, LD_String), LD);
		if (m_vars.contains(LD_prime_sp_String))
			outputVar(pop, subPopVar_String(*it, LD_prime_String), D_prime);
		if (m_vars.contains(R2_sp_String))
			outputVar(pop, subPopVar_String(*it, R2_String), R2);
		if (m_vars.contains(ChiSq_sp_String))
			outputVar(pop, subPopVar_String(*it, ChiSq_String), ChiSq);
		if (m_vars.contains(ChiSq_p_sp_String))
			outputVar(pop, subPopVar_String(*it, ChiSq_p_String), ChiSq_p);
		if (m_vars.contains(CramerV_sp_String))
			outputVar(pop, subPopVar_String(*it, CramerV_String), CramerV);
	}
	// output statistics for all (virtual) subpopulations
	// calculate statistics
	UINT ldSize = 0;
	if (m_vars.contains(LD_String) || m_vars.contains(LD_prime_String) || m_vars.contains(R2_String))
		ldSize = m_LD.size();
	vectorf LD(ldSize);
	vectorf D_prime(ldSize);
	vectorf R2(ldSize);
	UINT assoSize = 0;
	if (m_vars.contains(ChiSq_String) || m_vars.contains(ChiSq_p_String) || m_vars.contains(CramerV_String))
		assoSize = m_LD.size();
	vectorf ChiSq(assoSize);
	vectorf ChiSq_p(assoSize);
	vectorf CramerV(assoSize);
	calculateLD(lociMap, allAlleleCnt, allHaploCnt, LD, D_prime, R2,
		ChiSq, ChiSq_p, CramerV);

	// output statistics for subpopulation
	if (m_vars.contains(LD_String))
		outputVar(pop, LD_String, LD);
	if (m_vars.contains(LD_prime_String))
		outputVar(pop, LD_prime_String, D_prime);
	if (m_vars.contains(R2_String))
		outputVar(pop, R2_String, R2);
	if (m_vars.contains(ChiSq_String))
		outputVar(pop, ChiSq_String, ChiSq);
	if (m_vars.contains(ChiSq_p_String))
		outputVar(pop, ChiSq_p_String, ChiSq_p);
	if (m_vars.contains(CramerV_String))
		outputVar(pop, CramerV_String, CramerV);
	return true;
}


statAssociation::statAssociation(const vectoru & loci,
	const subPopList & subPops, const stringList & vars, const string & suffix)
	: m_loci(loci), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = {
		Allele_ChiSq_String,	Allele_ChiSq_p_String,
		Geno_ChiSq_String,		Geno_ChiSq_p_String,
		Armitage_p_String,
		Allele_ChiSq_sp_String, Allele_ChiSq_p_sp_String,
		Geno_ChiSq_sp_String,	Geno_ChiSq_p_sp_String,
		Armitage_p_sp_String,
		""
	};
	const char * defaultVars[] = { Allele_ChiSq_p_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);
}


string statAssociation::describe(bool format)
{
	string desc;

	if (!m_loci.empty())
		desc += "perform association tests";
	return desc;
}


void statAssociation::alleleChiSqTest(const ALLELECNT & caseCnt,
                                      const ALLELECNT & controlCnt, double & chisq,
                                      double & chisq_p)
{
	// figure out alleles
	map<Allele, int> alleles;
	ALLELECNT::const_iterator cnt = caseCnt.begin();
	ALLELECNT::const_iterator cntEnd = caseCnt.end();
	for (; cnt != cntEnd; ++cnt)
		if (alleles.find(cnt->first) == alleles.end()) {
			UINT idx = alleles.size();
			alleles[cnt->first] = idx;
		}
	cnt = controlCnt.begin();
	cntEnd = controlCnt.end();
	for (; cnt != cntEnd; ++cnt)
		if (alleles.find(cnt->first) == alleles.end()) {
			UINT idx = alleles.size();
			alleles[cnt->first] = idx;
		}

	vector<vectoru> table(2);
	for (size_t i = 0; i < 2; ++i)
		table[i].resize(alleles.size(), 0);
	// fill the table
	cnt = caseCnt.begin();
	cntEnd = caseCnt.end();
	for (; cnt != cntEnd; ++cnt)
		table[0][alleles[cnt->first]] = cnt->second;
	cnt = controlCnt.begin();
	cntEnd = controlCnt.end();
	for (; cnt != cntEnd; ++cnt)
		table[1][alleles[cnt->first]] = cnt->second;
	chisqTest(table, chisq, chisq_p);
}


void statAssociation::genoChiSqTest(const GENOCNT & caseCnt,
                                    const GENOCNT & controlCnt, double & chisq,
                                    double & chisq_p)
{
	// figure out alleles
	map<std::pair<UINT, UINT>, int> genotypes;
	GENOCNT::const_iterator cnt = caseCnt.begin();
	GENOCNT::const_iterator cntEnd = caseCnt.end();
	for (; cnt != cntEnd; ++cnt)
		if (genotypes.find(cnt->first) == genotypes.end()) {
			UINT idx = genotypes.size();
			genotypes[cnt->first] = idx;
		}
	cnt = controlCnt.begin();
	cntEnd = controlCnt.end();
	for (; cnt != cntEnd; ++cnt)
		if (genotypes.find(cnt->first) == genotypes.end()) {
			UINT idx = genotypes.size();
			genotypes[cnt->first] = idx;
		}

	vector<vectoru> table(2);
	for (size_t i = 0; i < 2; ++i)
		table[i].resize(genotypes.size(), 0);
	// fill the table
	cnt = caseCnt.begin();
	cntEnd = caseCnt.end();
	for (; cnt != cntEnd; ++cnt)
		table[0][genotypes[cnt->first]] = cnt->second;
	cnt = controlCnt.begin();
	cntEnd = controlCnt.end();
	for (; cnt != cntEnd; ++cnt)
		table[1][genotypes[cnt->first]] = cnt->second;
	chisqTest(table, chisq, chisq_p);
}


double statAssociation::armitageTest(const GENOCNT & caseCnt,
                                     const GENOCNT & ctrlCnt)
{
	// figure out alleles and their frequencies
	map<Allele, int> alleles;
	GENOCNT::const_iterator cnt = caseCnt.begin();
	GENOCNT::const_iterator cntEnd = caseCnt.end();
	for (; cnt != cntEnd; ++cnt) {
		alleles[cnt->first.first] += cnt->second;
		alleles[cnt->first.second] += cnt->second;
	}
	cnt = ctrlCnt.begin();
	cntEnd = ctrlCnt.end();
	for (; cnt != cntEnd; ++cnt) {
		alleles[cnt->first.first] += cnt->second;
		alleles[cnt->first.second] += cnt->second;
	}
	// figure out major and minor allele
	DBG_FAILIF(alleles.size() > 2, ValueError,
		"Armitage trend test can only be applied to diallelic markers.");
	if (alleles.size() != 2)
		return 1.;
	map<Allele, int>::const_iterator first = alleles.begin();
	map<Allele, int>::const_iterator second = first;
	++second;
	Allele major = 0;
	Allele minor = 0;
	if (first->second > second->second) {
		major = first->first;
		minor = second->first;
	} else {
		major = second->first;
		minor = first->first;
	}
	vector<vectoru> table(2);
	for (size_t i = 0; i < 2; ++i)
		table[i].resize(3, 0);
	// fill the table
	GENOCNT::const_iterator it = caseCnt.find(GENOCNT::key_type(major, major));
	if (it != caseCnt.end())
		table[1][0] = it->second;
	it = caseCnt.find(major > minor ? GENOCNT::key_type(minor, major) : GENOCNT::key_type(major, minor));
	if (it != caseCnt.end())
		table[1][1] = it->second;
	it = caseCnt.find(GENOCNT::key_type(minor, minor));
	if (it != caseCnt.end())
		table[1][2] = it->second;
	it = ctrlCnt.find(GENOCNT::key_type(major, major));
	if (it != ctrlCnt.end())
		table[0][0] = it->second;
	it = ctrlCnt.find(major > minor ? GENOCNT::key_type(minor, major) : GENOCNT::key_type(major, minor));
	if (it != ctrlCnt.end())
		table[0][1] = it->second;
	it = ctrlCnt.find(GENOCNT::key_type(minor, minor));
	if (it != ctrlCnt.end())
		table[0][2] = it->second;
	vectorf weight(3);
	for (size_t i = 0; i < 3; ++i)
		weight[i] = i;
	return armitageTrendTest(table, weight);
}


bool statAssociation::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	vectoru chromTypes;
	for (size_t i = 0; i < m_loci.size(); ++i)
		chromTypes.push_back(pop.chromType(pop.chromLocusPair(m_loci[i]).first));

	UINT ply = pop.ploidy();
	bool hasAlleleTest = false;
	bool hasGenoTest = false;
	for (size_t i = 0; i < m_vars.elems().size(); ++i) {
		if (m_vars.elems()[i].compare(0, 6, "Allele") == 0)
			hasAlleleTest = true;
		if (m_vars.elems()[i].compare(0, 4, "Geno") == 0 || m_vars.elems()[i].compare(0, 8, "Armitage") == 0) {
			DBG_FAILIF(ply != 2, ValueError, "Genotype test can only be performed for diploid populations.");
			hasGenoTest = true;
		}
	}
	// count for all specified subpopulations
	UINT nLoci = m_loci.size();
	ALLELECNTLIST allCaseAlleleCnt(nLoci);
	ALLELECNTLIST allCtrlAlleleCnt(nLoci);
	GENOCNTLIST allCaseGenoCnt(nLoci);
	GENOCNTLIST allCtrlGenoCnt(nLoci);
	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	for (; it != itEnd; ++it) {
		ALLELECNTLIST caseAlleleCnt(nLoci);
		ALLELECNTLIST ctrlAlleleCnt(nLoci);
		GENOCNTLIST caseGenoCnt(nLoci);
		GENOCNTLIST ctrlGenoCnt(nLoci);

		pop.activateVirtualSubPop(*it);

		IndIterator ind = pop.indIterator(it->subPop());
		for (; ind.valid(); ++ind) {
			if (hasAlleleTest) {
				for (UINT p = 0; p < ply; ++p) {
					if (ply == 2 && p == 1 && ind->sex() == MALE && pop.isHaplodiploid())
						continue;
					GenoIterator geno = ind->genoBegin(p);
					// allele count
					for (size_t idx = 0; idx < nLoci; ++idx) {
						if (ply == 2 && chromTypes[idx] == CHROMOSOME_Y && ind->sex() == FEMALE)
							continue;
						if (ply == 2 && ((chromTypes[idx] == CHROMOSOME_X && p == 1) ||
						                 (chromTypes[idx] == CHROMOSOME_Y && p == 0)) && ind->sex() == MALE)
							continue;
						if (ind->affected())
							caseAlleleCnt[idx][*(geno + m_loci[idx])]++;
						else
							ctrlAlleleCnt[idx][*(geno + m_loci[idx])]++;
					}
				}
			}
			// genotype
			if (hasGenoTest) {
				GenoIterator geno1 = ind->genoBegin(0);
				GenoIterator geno2 = ind->genoBegin(1);
				for (size_t idx = 0; idx < nLoci; ++idx) {
					if (chromTypes[idx] == CHROMOSOME_X || chromTypes[idx] == CHROMOSOME_Y)
						continue;
					Allele a1 = *(geno1 + m_loci[idx]);
					Allele a2 = *(geno2 + m_loci[idx]);
					if (a1 > a2)
						std::swap(a1, a2);
					if (ind->affected())
						caseGenoCnt[idx][GENOCNT::key_type(a1, a2)]++;
					else
						ctrlGenoCnt[idx][GENOCNT::key_type(a1, a2)]++;
				}
			}
		}
		pop.deactivateVirtualSubPop(it->subPop());
		//
		// output variable.
		if (m_vars.contains(Allele_ChiSq_sp_String) || m_vars.contains(Allele_ChiSq_p_sp_String)) {
			intDict chisq;
			intDict chisq_p;
			for (size_t i = 0; i < m_loci.size(); ++i)
				alleleChiSqTest(caseAlleleCnt[i], ctrlAlleleCnt[i], chisq[m_loci[i]], chisq_p[m_loci[i]]);
			if (m_vars.contains(Allele_ChiSq_sp_String))
				pop.getVars().setIntDictVar(subPopVar_String(*it, Allele_ChiSq_String) + m_suffix, chisq);
			if (m_vars.contains(Allele_ChiSq_p_sp_String))
				pop.getVars().setIntDictVar(subPopVar_String(*it, Allele_ChiSq_p_String) + m_suffix, chisq_p);
		}
		if (m_vars.contains(Geno_ChiSq_sp_String) || m_vars.contains(Geno_ChiSq_p_sp_String)) {
			intDict chisq;
			intDict chisq_p;
			for (size_t i = 0; i < m_loci.size(); ++i)
				genoChiSqTest(caseGenoCnt[i], ctrlGenoCnt[i], chisq[m_loci[i]], chisq_p[m_loci[i]]);
			if (m_vars.contains(Geno_ChiSq_sp_String))
				pop.getVars().setIntDictVar(subPopVar_String(*it, Geno_ChiSq_String) + m_suffix, chisq);
			if (m_vars.contains(Geno_ChiSq_p_sp_String))
				pop.getVars().setIntDictVar(subPopVar_String(*it, Geno_ChiSq_p_String) + m_suffix, chisq_p);
		}
		if (m_vars.contains(Armitage_p_sp_String)) {
			intDict pvalues;
			for (size_t i = 0; i < m_loci.size(); ++i)
				pvalues[m_loci[i]] = armitageTest(caseGenoCnt[i], ctrlGenoCnt[i]);
			pop.getVars().setIntDictVar(subPopVar_String(*it, Armitage_p_String) + m_suffix, pvalues);
		}
		// total allele count
		if (hasAlleleTest) {
			for (size_t idx = 0; idx < nLoci; ++idx) {
				ALLELECNT::const_iterator cnt = caseAlleleCnt[idx].begin();
				ALLELECNT::const_iterator cntEnd = caseAlleleCnt[idx].end();
				for (; cnt != cntEnd; ++cnt)
					allCaseAlleleCnt[idx][cnt->first] += cnt->second;
				cnt = ctrlAlleleCnt[idx].begin();
				cntEnd = ctrlAlleleCnt[idx].end();
				for (; cnt != cntEnd; ++cnt)
					allCtrlAlleleCnt[idx][cnt->first] += cnt->second;
			}
		}
		if (hasGenoTest) {
			for (size_t idx = 0; idx < nLoci; ++idx) {
				GENOCNT::const_iterator cnt = caseGenoCnt[idx].begin();
				GENOCNT::const_iterator cntEnd = caseGenoCnt[idx].end();
				for (; cnt != cntEnd; ++cnt)
					allCaseGenoCnt[idx][cnt->first] += cnt->second;
				cnt = ctrlGenoCnt[idx].begin();
				cntEnd = ctrlGenoCnt[idx].end();
				for (; cnt != cntEnd; ++cnt)
					allCtrlGenoCnt[idx][cnt->first] += cnt->second;
			}
		}
	}
	//
	if (m_vars.contains(Allele_ChiSq_String) || m_vars.contains(Allele_ChiSq_p_String)) {
		intDict chisq;
		intDict chisq_p;
		for (size_t i = 0; i < m_loci.size(); ++i)
			alleleChiSqTest(allCaseAlleleCnt[i], allCtrlAlleleCnt[i],
				chisq[m_loci[i]], chisq_p[m_loci[i]]);
		// output variable.
		if (m_vars.contains(Allele_ChiSq_String))
			pop.getVars().setIntDictVar(Allele_ChiSq_String + m_suffix, chisq);
		if (m_vars.contains(Allele_ChiSq_p_String))
			pop.getVars().setIntDictVar(Allele_ChiSq_p_String + m_suffix, chisq_p);
	}
	if (m_vars.contains(Geno_ChiSq_String) || m_vars.contains(Geno_ChiSq_p_String)) {
		intDict chisq;
		intDict chisq_p;
		for (size_t i = 0; i < m_loci.size(); ++i)
			genoChiSqTest(allCaseGenoCnt[i], allCtrlGenoCnt[i],
				chisq[m_loci[i]], chisq_p[m_loci[i]]);
		if (m_vars.contains(Geno_ChiSq_String))
			pop.getVars().setIntDictVar(Geno_ChiSq_String + m_suffix, chisq);
		if (m_vars.contains(Geno_ChiSq_p_String))
			pop.getVars().setIntDictVar(Geno_ChiSq_p_String + m_suffix, chisq_p);
	}
	if (m_vars.contains(Armitage_p_String)) {
		intDict pvalues;
		for (size_t i = 0; i < m_loci.size(); ++i)
			pvalues[m_loci[i]] = armitageTest(allCaseGenoCnt[i], allCtrlGenoCnt[i]);
		pop.getVars().setIntDictVar(Armitage_p_String + m_suffix, pvalues);
	}
	return true;
}


statNeutrality::statNeutrality(const vectoru & loci, const subPopList & subPops,
	const stringList & vars, const string & suffix) :
	m_loci(loci), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = {
		Neutra_Pi_String, Neutra_Pi_sp_String, ""
	};
	const char * defaultVars[] = { Neutra_Pi_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);

	sort(m_loci.begin(), m_loci.end());
	for (UINT i = 1; i < m_loci.size(); ++i) {
		DBG_FAILIF(m_loci[i - 1] == m_loci[i], ValueError,
			"Duplicated loci occurred in the loci list.");
	}
}


string statNeutrality::describe(bool format)
{
	string desc;

	if (!m_loci.empty())
		desc += "perform neutrality tests";
	return desc;
}


double statNeutrality::calcPi(HAPLOLIST::const_iterator begin, HAPLOLIST::const_iterator end)
{
	double diffCnt = 0;
	int numComparison = 0;

	HAPLOLIST::const_iterator it = begin;

	for (; it != end; ++it) {
		HAPLOLIST::const_iterator it1 = it;
		for (++it1; it1 != end; ++it1) {
			const vectora & seq1 = *it;
			const vectora & seq2 = *it1;
			size_t sz = seq1.size();
			for (size_t i = 0; i < sz; ++i)
				diffCnt += seq1[i] != seq2[i];
			++numComparison;
		}
	}
	// return 0 if there is only one sequence
	return numComparison == 0 ? 0 : diffCnt / numComparison;
}


bool statNeutrality::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	UINT nLoci = m_loci.size();
	int chromType = pop.chromType(pop.chromLocusPair(m_loci[0]).first);
#ifndef OPTIMIZED
	for (size_t i = 1; i < nLoci; ++i) {
		DBG_ASSERT(chromType == pop.chromType(pop.chromLocusPair(m_loci[i]).first),
			ValueError, "All loci must be from chromosomes of the same type.");
	}
#endif
	// count for all specified subpopulations
	HAPLOLIST allHaplotypes;
	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	UINT ply = pop.ploidy();
	for (; it != itEnd; ++it) {
		pop.activateVirtualSubPop(*it);

		UINT spBegin = allHaplotypes.size();
		// go through all individual
		IndIterator ind = pop.indIterator(it->subPop());
		for (; ind.valid(); ++ind) {
			vectora haplotype(nLoci);
			for (size_t p = 0; p < ply; ++p) {
				if (p == 1 && ind->sex() == MALE && pop.isHaplodiploid())
					continue;
				if (chromType == CHROMOSOME_Y && ind->sex() == FEMALE)
					continue;
				if (((chromType == CHROMOSOME_X && p == 1) ||
				     (chromType == CHROMOSOME_Y && p == 0)) && ind->sex() == MALE)
					continue;
				for (size_t idx = 0; idx < nLoci; ++idx)
					haplotype[idx] = ToAllele(ind->allele(m_loci[idx], p));
				allHaplotypes.push_back(haplotype);
			}
		}
		// output variable.
		if (m_vars.contains(Neutra_Pi_sp_String))
			pop.getVars().setDoubleVar(subPopVar_String(*it, Neutra_Pi_String) + m_suffix,
				calcPi(allHaplotypes.begin() + spBegin, allHaplotypes.end()));
		pop.deactivateVirtualSubPop(it->subPop());
	}

	if (m_vars.contains(Neutra_Pi_String))
		pop.getVars().setDoubleVar(Neutra_Pi_String + m_suffix, calcPi(allHaplotypes.begin(), allHaplotypes.end()));
	return true;
}


statStructure::statStructure(const vectoru & Fst, const subPopList & subPops, const stringList & vars, const string & suffix)
	: m_loci(Fst), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = {
		fst_String, fis_String, fit_String,
		Fst_String, Fis_String, Fit_String,
		Gst_String, gst_String, ""
	};
	const char * defaultVars[] = { Fst_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);
}


string statStructure::describe(bool format)
{
	string desc;

	if (!m_loci.empty())
		desc += "evaluate population structure";
	return desc;
}


void statStructure::calcGst_Nei73(const vectoru & n_i, LOCIFREQLIST & alleleFreq,
                                  const ALLELELIST & alleles, double & Gst, intDict & gst)
{
	double H_t_all = 0.;
	double D_st_all = 0.;

	UINT numSP = n_i.size();
	double n = accumulate(n_i.begin(), n_i.end(), 0);

	// for each locus
	for (size_t st = 0; st < m_loci.size(); ++st) {
		UINT loc = m_loci[st];
		// D_st = Sum_i,j D_ij / s^2
		double D_st = 0;
		// i, j for subpopulation
		for (size_t i = 0; i < numSP; ++i) {
			for (size_t j = 0; j < numSP; ++j) {
				// D_ij = Sum_k (x_ik - x_jk)^2 /2 (i,j are subpops, k is allele)
				double D_ij = 0;
				// for each allele
				ALLELES::const_iterator aIt = alleles[st].begin();
				ALLELES::const_iterator aEnd = alleles[st].end();
				for (; aIt != aEnd; ++aIt) {
					Allele allele = aIt->first;
					D_ij += pow(alleleFreq[i][loc][allele] - alleleFreq[j][loc][allele], 2);
				}
				D_ij /= 2.;
				D_st += D_ij;
			}
		}
		D_st /= (numSP * numSP);
		//
		// J_t = Sum_k (x_dot_k^2)
		//
		// NOTE: w_i is chosen as n_i/n instead of 1/numSP as
		// used in the paper.
		double J_t = 0;
		ALLELES::const_iterator aIt = alleles[st].begin();
		ALLELES::const_iterator aEnd = alleles[st].end();
		for (; aIt != aEnd; ++aIt) {
			Allele allele = aIt->first;
			double x_dotk = 0;
			for (size_t i = 0; i < numSP; ++i)
				x_dotk += n_i[i] * alleleFreq[i][loc][allele];
			x_dotk /= n;
			J_t += pow(x_dotk, 2);
		}
		gst[loc] = fcmp_eq(J_t, 1.0) ? 0 : D_st / (1.0 - J_t);
		H_t_all += 1 - J_t;
		D_st_all += D_st;
	}
	Gst = fcmp_eq(H_t_all, 0.) ? 0 : D_st_all / H_t_all;
}


void statStructure::calcFst_WC84(const vectoru & n_i, LOCIFREQLIST & alleleFreq, LOCIFREQLIST & heteroFreq,
                                 const ALLELELIST & alleles, double & Fst, double & Fis, double & Fit,
                                 intDict & fst, intDict & fis, intDict & fit)
{
	double aa = 0.;
	double bb = 0.;
	double cc = 0.;

	// vector to store p[i]
	vectorf p_i = vectorf(n_i.size());
	double n = accumulate(n_i.begin(), n_i.end(), 0);

	// calculate Fst for each locus
	for (size_t st = 0; st < m_loci.size(); ++st) {
		int loc = m_loci[st];

		// n_bar
		double r = n_i.size();
		double n_bar = n / r;

		// n_c
		double n_c = n;
		for (int i = 0; i < r; ++i)
			n_c -= n_i[i] * n_i[i] / n;
		n_c /= (r - 1);

		double a = 0.0, b = 0.0, c = 0.0;

		ALLELES::const_iterator aIt = alleles[st].begin();
		ALLELES::const_iterator aEnd = alleles[st].end();
		for (; aIt != aEnd; ++aIt) {
			Allele allele = aIt->first;
			// p_i
			for (int sp = 0; sp < r; ++sp)
				p_i[sp] = alleleFreq[sp][loc][allele];

			// p_bar (there are 2n alleles, but this does not affect the result)
			double p_bar = 0;
			for (int sp = 0; sp < r; ++sp)
				p_bar += n_i[sp] * p_i[sp];
			p_bar /= n;

			// s^2
			double s_2 = 0;
			for (int sp = 0; sp < r; ++sp)
				s_2 += n_i[sp] * (p_i[sp] - p_bar) * (p_i[sp] - p_bar);
			s_2 /= (r - 1) * n_bar;

			// h_bar
			double h_bar = 0;
			for (int sp = 0; sp < r; ++sp)
				h_bar += heteroFreq[sp][loc][allele] * n_i[sp];
			h_bar /= n;

			// a, b, c
			a += n_bar / n_c * (s_2 - (p_bar * (1 - p_bar) - (r - 1.) / r * s_2 - h_bar / 4.) / (n_bar - 1.) );
			b += n_bar / (n_bar - 1) * (p_bar * (1 - p_bar) - (r - 1) / r * s_2 - (2 * n_bar - 1) / (4. * n_bar) * h_bar);
			c += h_bar / 2.;

			DBG_DO(DBG_STATOR, cerr << "allele " << allele << "\tn_c: " << n_c
				                    << "\tp_i: " << p_i << "\tp_bar: " << p_bar << "\ts^2: " << s_2 << "\th_bar:"
				                    << h_bar << "\ta: " << a << "\tb: " << b << "\tc: " << c << endl);
		}                                                                                 // each allele

		DBG_DO(DBG_STATOR, cerr << "Fst= " << a / (a + b + c) << endl);

		fst[loc] = fcmp_eq(a + b + c, 0.) ? 0. : (a / (a + b + c));
		fit[loc] = fcmp_eq(a + b + c, 0.) ? 1. : (1 - c / (a + b + c));
		fis[loc] = fcmp_eq(b + c, 0.) ? 1. : (1 - c / (b + c));

		aa += a;
		bb += b;
		cc += c;
	}
	Fst = fcmp_eq(aa + bb + cc, 0.) ? 0 : (aa / (aa + bb + cc));
	Fit = fcmp_eq(aa + bb + cc, 0.) ? 1. : (1 - cc / (aa + bb + cc));
	Fis = fcmp_eq(aa + bb + cc, 0) ? 1. : (1 - cc / (bb + cc));
}


bool statStructure::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	DBG_ASSERT(pop.ploidy() == 2, ValueError,
		"Fst statistics is available only for diploid populations.");

#ifndef OPTIMIZED
	for (size_t idx = 0; idx < m_loci.size(); ++idx) {
		int chromType = pop.chromType(pop.chromLocusPair(m_loci[idx]).first);
		DBG_FAILIF(chromType == CHROMOSOME_X || chromType == CHROMOSOME_Y, ValueError,
			"Fst can not be esimated from markers on sex chromosomes");
	}
#endif

	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	// count for all specified subpopulations
	vectoru n_i(0);
	ALLELELIST allAlleles(m_loci.size());
	LOCIFREQLIST alleleFreq(subPops.size());
	LOCIFREQLIST heteroFreq(subPops.size());
	for (size_t spIdx = 0; it != itEnd; ++it, ++spIdx) {
		if (m_vars.contains(AlleleNum_sp_String))
			pop.getVars().removeVar(subPopVar_String(*it, AlleleNum_String) + m_suffix);
		if (m_vars.contains(AlleleFreq_sp_String))
			pop.getVars().removeVar(subPopVar_String(*it, AlleleFreq_String) + m_suffix);

		pop.activateVirtualSubPop(*it);

		UINT spSize = 0;
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];
			FREQ & af = alleleFreq[spIdx][loc];
			FREQ & hf = heteroFreq[spIdx][loc];
			ALLELES & alleles = allAlleles[idx];
			UINT cnt = 0;

			// go through all alleles
			IndAlleleIterator a = pop.alleleIterator(loc, it->subPop());
			for (; a.valid(); ++cnt) {
				Allele a1 = *a++;
				Allele a2 = *a++;
				++af[a1];
				++af[a2];
				hf[a1] += a1 != a2;
				hf[a2] += a1 != a2;
				alleles[a1] = true;
				alleles[a2] = true;
			}
			// heterozygote frequency
			map<UINT, float>::iterator it = af.begin();
			map<UINT, float>::iterator itEnd = af.end();
			for (; it != itEnd; ++it)
				it->second /= 2 * cnt;
			// heterozygote frequency
			it = hf.begin();
			itEnd = hf.end();
			for (; it != itEnd; ++it)
				it->second /= cnt;
			//
			DBG_FAILIF(spSize != 0 && spSize != cnt, SystemError,
				"Subpopulation size counts are inconsistent. (locus not on autosome?)");
			spSize = cnt;
		}
		// (virtual) subpopulation size
		n_i.push_back(spSize);
		pop.deactivateVirtualSubPop(it->subPop());
	}

	// Nei's Gst
	double Gst = 0;
	intDict gst;
	calcGst_Nei73(n_i, alleleFreq, allAlleles, Gst, gst);
	if (m_vars.contains(Gst_String))
		pop.getVars().setDoubleVar(Gst_String + m_suffix, Gst);
	if (m_vars.contains(gst_String))
		pop.getVars().setIntDictVar(gst_String + m_suffix, gst);

	// Weir and Cockerham 1984 Fst
	double Fst = 0;
	double Fis = 0;
	double Fit = 0;
	intDict fst;
	intDict fis;
	intDict fit;
	calcFst_WC84(n_i, alleleFreq, heteroFreq, allAlleles, Fst, Fis, Fit, fst, fis, fit);
	// post results
	if (m_vars.contains(Fst_String))
		pop.getVars().setDoubleVar(Fst_String + m_suffix, Fst);
	if (m_vars.contains(Fis_String))
		pop.getVars().setDoubleVar(Fis_String + m_suffix, Fis);
	if (m_vars.contains(Fit_String))
		pop.getVars().setDoubleVar(Fit_String + m_suffix, Fit);
	if (m_vars.contains(fst_String))
		pop.getVars().setIntDictVar(fst_String + m_suffix, fst);
	if (m_vars.contains(fis_String))
		pop.getVars().setIntDictVar(fis_String + m_suffix, fis);
	if (m_vars.contains(fit_String))
		pop.getVars().setIntDictVar(fit_String + m_suffix, fit);
	return true;
}


statHWE::statHWE(const vectoru & loci,  const subPopList & subPops,
	const stringList & vars, const string & suffix)
	: m_loci(loci), m_subPops(subPops), m_vars(), m_suffix(suffix)
{
	const char * allowedVars[] = {
		HWE_String, HWE_sp_String, ""
	};
	const char * defaultVars[] = { HWE_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);
}


string statHWE::describe(bool format)
{
	string desc;

	if (!m_loci.empty())
		desc += "perform Hardy Weinberg tests";
	return desc;
}


vectoru statHWE::mapToCount(const GENOCNT & cnt)
{
	vectoru res(3, 0);

	DBG_FAILIF(cnt.size() > 3, ValueError,
		"More than three genotypes are detected for HWE test");

	GENOCNT::const_iterator it = cnt.begin();
	GENOCNT::const_iterator itEnd = cnt.end();
	for (; it != itEnd; ++it) {
		if (it->first.first != it->first.second)
			res[1] = it->second;
		else if (res[0] == 0)
			res[0] = it->second;
		else
			res[2] = it->second;
	}
	return res;
}


bool statHWE::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	DBG_FAILIF(pop.ploidy() != 2, ValueError,
		"HWE test is only available for diploid populations.");

#ifndef OPTIMIZED
	for (size_t i = 0; i < m_loci.size(); ++i) {
		int chromType = pop.chromType(pop.chromLocusPair(m_loci[i]).first);
		DBG_FAILIF(chromType == CHROMOSOME_X || chromType == CHROMOSOME_Y, ValueError,
			"Can not run HWE test on loci on sex chromosomes.");
	}
#endif
	// count for all specified subpopulations
	UINT nLoci = m_loci.size();
	GENOCNTLIST allGenoCnt(nLoci);
	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	if (subPops.allAvail())
		subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	for (; it != itEnd; ++it) {
		GENOCNTLIST genoCnt(nLoci);

		pop.activateVirtualSubPop(*it);

		IndIterator ind = pop.indIterator(it->subPop());
		for (; ind.valid(); ++ind) {
			GenoIterator geno1 = ind->genoBegin(0);
			GenoIterator geno2 = ind->genoBegin(1);
			for (size_t idx = 0; idx < nLoci; ++idx) {
				Allele a1 = *(geno1 + m_loci[idx]);
				Allele a2 = *(geno2 + m_loci[idx]);
				if (a1 > a2)
					std::swap(a1, a2);
				genoCnt[idx][GENOCNT::key_type(a1, a2)]++;
			}
		}
		pop.deactivateVirtualSubPop(it->subPop());
		//
		// output variable.
		if (m_vars.contains(HWE_sp_String)) {
			intDict hwe;
			for (size_t i = 0; i < m_loci.size(); ++i) {
				hwe[m_loci[i]] = hweTest(mapToCount(genoCnt[i]));
			}
			pop.getVars().setIntDictVar(subPopVar_String(*it, HWE_String) + m_suffix, hwe);
		}
		for (size_t idx = 0; idx < nLoci; ++idx) {
			GENOCNT::const_iterator cnt = genoCnt[idx].begin();
			GENOCNT::const_iterator cntEnd = genoCnt[idx].end();
			for (; cnt != cntEnd; ++cnt)
				allGenoCnt[idx][cnt->first] += cnt->second;
		}
	}
	//
	if (m_vars.contains(HWE_String)) {
		intDict hwe;
		for (size_t i = 0; i < m_loci.size(); ++i)
			hwe[m_loci[i]] = hweTest(mapToCount(allGenoCnt[i]));
		pop.getVars().setIntDictVar(HWE_String + m_suffix, hwe);
	}
	return true;
}


}


