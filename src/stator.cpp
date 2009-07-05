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

namespace simuPOP {

string pyEval::evaluate(population & pop)
{
	if (!m_exposePop.empty()) {
		PyObject * popObj = pyPopObj(static_cast<void *>(&pop));
		if (popObj == NULL)
			throw SystemError("Could not expose population pointer. Compiled with the wrong version of SWIG? ");

		// set dictionary variable pop to this object
		pop.setVar(m_exposePop, popObj);
	}

	m_expr.setLocalDict(pop.dict());
	string res = m_expr.valueAsString();
	if (!m_exposePop.empty())
		pop.removeVar(m_exposePop);
	return res;
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
	if (subPops.empty())
		subPops.useSubPopsFrom(pop);

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for ( ; sp != spEnd; ++sp) {
		if (sp->isVirtual())
			pop.activateVirtualSubPop(*sp, IteratableInds);
		IndIterator ind = const_cast<population &>(pop).indIterator(sp->subPop(), sp->isVirtual() ? IteratableInds : AllInds);
		for (; ind.valid(); ++ind) {
			string res = evalInfo(& * ind, false) ;
			if (!this->noOutput() ) {
				ostream & out = this->getOstream(pop.dict());
				out << res;
				this->closeOstream();
			}
		}
	}
	return true;
}


bool infoEval::applyDuringMating(population & pop, RawIndIterator offspring,
                                 individual * dad, individual * mom)
{
	m_dict = m_usePopVars ? pop.dict() : PyDict_New();

	string res = evalInfo(& * offspring, false);

	if (!this->noOutput() ) {
		ostream & out = this->getOstream(pop.dict());
		out << res;
		this->closeOstream();
	}
	return true;
}


bool infoExec::apply(population & pop)
{
	m_dict = m_usePopVars ? pop.dict() : PyDict_New();

	subPopList subPops = applicableSubPops();
	if (subPops.empty())
		subPops.useSubPopsFrom(pop);

	OperationType oType = m_simpleStmt.operation();
	string oVar = m_simpleStmt.var();
	double oValue = m_simpleStmt.value();
	UINT oVarIdx = 0;
	if (oType != NoOperation)
		oVarIdx = pop.infoIdx(oVar);

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for ( ; sp != spEnd; ++sp) {
		if (sp->isVirtual())
			pop.activateVirtualSubPop(*sp, IteratableInds);
		IndIterator ind = const_cast<population &>(pop).indIterator(sp->subPop(), sp->isVirtual() ? IteratableInds : AllInds);
		for (; ind.valid(); ++ind) {
			switch (m_simpleStmt.operation()) {
			case NoOperation:
				evalInfo(& * ind, true);
				break;
			case Assignment:
				ind->setInfo(oValue, oVarIdx);
				break;
			case Increment:
				ind->setInfo(ind->info(oVarIdx) + oValue, oVarIdx);
				break;
			case Decrement:
				ind->setInfo(ind->info(oVarIdx) - oValue, oVarIdx);
				break;
			case MultipliedBy:
				ind->setInfo(ind->info(oVarIdx) * oValue, oVarIdx);
				break;
			default:
				throw RuntimeError("Incorrect operation type");
			}
		}
	}
	return true;
}


bool infoExec::applyDuringMating(population & pop, RawIndIterator offspring,
                                 individual * dad, individual * mom)
{
	m_dict = m_usePopVars ? pop.dict() : PyDict_New();

	evalInfo(& * offspring, true);
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
	bool numOfMale,
	//
	bool numOfAffected,
	//
	const uintList & alleleFreq,
	//
	const uintList & heteroFreq,
	const uintList & expHetero,
	const strDict & expHetero_param,
	const uintList & homoFreq,
	//
	const uintList & genoFreq,
	const strDict & genoFreq_param,
	const intMatrix & haploFreq,
	//
	const intMatrix & LD,
	const strDict & LD_param,
	//
	const uintList & association,
	//
	const uintList & neutrality,
	//
	const uintList & Fst,
	const strDict & Fst_param,
	//
	const uintList & HWE,
	//
	const stringList & vars,
	// regular parameters
	const stringFunc & output,
	int stage, int begin, int end, int step, const intList & at,
	const repList & reps, const subPopList & subPops, const stringList & infoFields)
	: baseOperator("", stage, begin, end, step, at, reps, subPops, infoFields),
	// the order of initialization is meaningful since they may depend on each other
	m_popSize(popSize, subPops, vars),
	m_numOfMale(numOfMale, subPops, vars),
	m_numOfAffected(numOfAffected, subPops, vars),
	m_alleleFreq(alleleFreq.elems(), subPops, vars),
	m_heteroFreq(heteroFreq.elems(), homoFreq.elems()),
	m_expHetero(m_alleleFreq, expHetero.elems(), expHetero_param),
	m_genoFreq(genoFreq.elems(), genoFreq_param),
	m_haploFreq(haploFreq),
	m_LD(m_alleleFreq, m_haploFreq, LD, LD_param),
	m_association(association.elems(), subPops),
	m_neutrality(neutrality.elems(), subPops),
	m_Fst(m_alleleFreq, m_heteroFreq, Fst.elems(), Fst_param),
	m_HWE(m_genoFreq, HWE.elems())
{
}


stat::stat(const stat & rhs) :
	baseOperator(rhs),
	// the order of initialization is meaningful since they may depend on each other
	m_popSize(rhs.m_popSize),
	m_numOfMale(rhs.m_numOfMale),
	m_numOfAffected(rhs.m_numOfAffected),
	m_alleleFreq(rhs.m_alleleFreq),
	m_heteroFreq(rhs.m_heteroFreq),
	m_expHetero(m_alleleFreq, rhs.m_expHetero),     // and this one
	m_genoFreq(rhs.m_genoFreq),
	m_haploFreq(rhs.m_haploFreq),
	m_LD(m_alleleFreq, m_haploFreq, rhs.m_LD), // and this one
	m_association(rhs.m_association),
	m_neutrality(rhs.m_neutrality),
	m_Fst(m_alleleFreq, m_heteroFreq, rhs.m_Fst),
	m_HWE(m_genoFreq, rhs.m_HWE)
{
}


bool stat::apply(population & pop)
{
	return m_popSize.apply(pop) &&
	       m_numOfMale.apply(pop) &&
	       m_numOfAffected.apply(pop) &&
	       m_alleleFreq.apply(pop) &&
	       m_heteroFreq.apply(pop) &&
	       m_expHetero.apply(pop) &&
	       m_genoFreq.apply(pop) &&
	       m_haploFreq.apply(pop) &&
	       m_LD.apply(pop) &&
	       m_association.apply(pop) &&
	       m_neutrality.apply(pop) &&
	       m_Fst.apply(pop) &&
	       m_HWE.apply(pop);
}


bool statPopSize::apply(population & pop)
{
	if (!m_isActive)
		return true;

	// popSize = ...
	ULONG popSize = 0;
	// for each (virtual) subpopulation
	subPopList subPops = m_subPops;
	subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	for (; it != itEnd; ++it) {
		ULONG spPopSize = pop.subPopSize(*it);
		popSize += spPopSize;
		if (m_vars.contains(popSize_sp_String))
			pop.setIntVar(subPopVar_String(*it, popSize_String), spPopSize);
	}
	// NOTE: popSize does not have to be the total population size
	if (m_vars.contains(popSize_String))
		pop.setIntVar(popSize_String, popSize);
	// subPopSize = ...
	// type mismatch, can not use subPopSizes() directly.
	if (m_vars.contains(subPopSize_String)) {
		vectori spSize(pop.numSubPop());
		for (size_t sp = 0; sp < spSize.size(); ++sp)
			spSize[sp] = pop.subPopSize(sp);
		pop.setIntVectorVar(subPopSize_String, spSize);
	}
	return true;
}


bool statNumOfMale::apply(population & pop)
{
	if (!m_isActive)
		return true;

	ULONG allMaleCnt = 0;
	ULONG allFemaleCnt = 0;
	ULONG allTotalCnt = 0;
	// for each subpopulation.
	subPopList subPops = m_subPops;
	subPops.useSubPopsFrom(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for (; sp != spEnd; ++sp) {
		ULONG maleCnt = 0;
		ULONG femaleCnt = 0;
		ULONG totalCnt = 0;
		if (sp->isVirtual())
			pop.activateVirtualSubPop(*sp);

		IndIterator it = pop.indIterator(sp->subPop());
		for (; it.valid(); ++it)
			if (it->sex() == Male)
				maleCnt++;
			else
				femaleCnt++;
		
		if (sp->isVirtual())
			pop.deactivateVirtualSubPop(sp->subPop());

		totalCnt = maleCnt + femaleCnt;
			
		if (m_vars.contains(string(numOfMale_String) + "_sp"))
			pop.setIntVar(subPopVar_String(*sp, numOfMale_String), maleCnt);
		if (m_vars.contains(string(propOfMale_String) + "_sp"))
			pop.setDoubleVar(subPopVar_String(*sp, propOfMale_String),
				totalCnt == 0 ? 0 : static_cast<double>(maleCnt) / totalCnt);
		if (m_vars.contains(string(numOfFemale_String) + "_sp"))
			pop.setIntVar(subPopVar_String(*sp, numOfFemale_String), femaleCnt);
		if (m_vars.contains(string(propOfFemale_String) + "_sp"))
			pop.setDoubleVar(subPopVar_String(*sp, propOfFemale_String),
				totalCnt == 0 ? 0 : static_cast<double>(femaleCnt) / totalCnt);

		allMaleCnt += maleCnt;
		allFemaleCnt += femaleCnt;
		allTotalCnt += totalCnt;
	}

	// output whole population
	if (m_vars.contains(numOfMale_String))
		pop.setIntVar(numOfMale_String, allMaleCnt);
	if (m_vars.contains(propOfMale_String))
		pop.setDoubleVar(propOfMale_String, allTotalCnt == 0 ? 0. : static_cast<double>(allMaleCnt) / allTotalCnt);
	if (m_vars.contains(numOfFemale_String))
		pop.setIntVar(numOfFemale_String, allFemaleCnt);
	if (m_vars.contains(propOfFemale_String))
		pop.setDoubleVar(propOfFemale_String, allTotalCnt == 0 ? 0 : static_cast<double>(allFemaleCnt) / allTotalCnt);
	return true;
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
	subPops.useSubPopsFrom(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for (; sp != spEnd; ++sp) {
		ULONG affectedCnt = 0;
		ULONG unaffectedCnt = 0;
		ULONG totalCnt = 0;
		if (sp->isVirtual())
			pop.activateVirtualSubPop(*sp);

		IndIterator it = pop.indIterator(sp->subPop());
		for (; it.valid(); ++it)
			if (it->affected())
				affectedCnt++;
			else
				unaffectedCnt++;
		
		if (sp->isVirtual())
			pop.deactivateVirtualSubPop(sp->subPop());

		totalCnt = affectedCnt + unaffectedCnt;
			
		if (m_vars.contains(string(numOfAffected_String) + "_sp"))
			pop.setIntVar(subPopVar_String(*sp, numOfAffected_String), affectedCnt);
		if (m_vars.contains(string(propOfAffected_String) + "_sp"))
			pop.setDoubleVar(subPopVar_String(*sp, propOfAffected_String),
				totalCnt == 0 ? 0 : static_cast<double>(affectedCnt) / totalCnt);
		if (m_vars.contains(string(numOfUnaffected_String) + "_sp"))
			pop.setIntVar(subPopVar_String(*sp, numOfUnaffected_String), unaffectedCnt);
		if (m_vars.contains(string(propOfUnaffected_String) + "_sp"))
			pop.setDoubleVar(subPopVar_String(*sp, propOfUnaffected_String),
				totalCnt == 0 ? 0 : static_cast<double>(unaffectedCnt) / totalCnt);

		allAffectedCnt += affectedCnt;
		allUnaffectedCnt += unaffectedCnt;
		allTotalCnt += totalCnt;
	}

	// output whole population
	if (m_vars.contains(numOfAffected_String))
		pop.setIntVar(numOfAffected_String, allAffectedCnt);
	if (m_vars.contains(propOfAffected_String))
		pop.setDoubleVar(propOfAffected_String, allTotalCnt == 0 ? 0. : static_cast<double>(allAffectedCnt) / allTotalCnt);
	if (m_vars.contains(numOfUnaffected_String))
		pop.setIntVar(numOfUnaffected_String, allUnaffectedCnt);
	if (m_vars.contains(propOfUnaffected_String))
		pop.setDoubleVar(propOfUnaffected_String, allTotalCnt == 0 ? 0 : static_cast<double>(allUnaffectedCnt) / allTotalCnt);

	return true;
}


void statAlleleFreq::addLocus(UINT locus, const subPopList & subPops,
                              const stringList & vars)
{
	if (find(m_loci.begin(), m_loci.end(), locus) == m_loci.end())
		m_loci.push_back(locus);
	if (!subPops.empty()) {
		subPopList::const_iterator it = subPops.begin();
		subPopList::const_iterator itEnd = subPops.end();
		for (; it != itEnd; ++it)
			if (!m_subPops.contains(*it))
				m_subPops.push_back(*it);
	}
	if (!vars.empty() && m_vars.empty()) {
		vectorstr::const_iterator it = vars.elems().begin();
		vectorstr::const_iterator itEnd = vars.elems().end();
		for (; it != itEnd; ++it)
			if (!m_vars.contains(*it))
				m_vars.push_back(*it);
	}
}


bool statAlleleFreq::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	DBG_DO(DBG_STATOR, cout << "Calculated allele frequency for loci " << m_loci << endl);

	// count for all specified subpopulations
	vector<vectori> alleleCnt(m_loci.size(), vectori(2, 0));
	vectorlu allAllelesCnt(m_loci.size(), 0);
	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	for (; it != itEnd; ++it) {
		if (m_vars.contains(string(AlleleNum_String) + "_sp"))
			pop.removeVar(subPopVar_String(*it, AlleleNum_String));
		if (m_vars.contains(string(AlleleFreq_String) + "_sp"))
			pop.removeVar(subPopVar_String(*it, AlleleFreq_String));

		if (it->isVirtual())
			pop.activateVirtualSubPop(*it);

		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];

			vectori alleles(2, 0);
			size_t allAlleles = 0;

			// go through all alleles
			IndAlleleIterator a = pop.alleleIterator(loc, it->subPop());
			// use allAllelel here because some marker does not have full number
			// of alleles (e.g. markers on chromosome X and Y).
			for (; a.valid(); ++a) {
				if (AlleleUnsigned(*a) >= alleles.size())
					alleles.resize(*a + 1, 0);
				alleles[*a]++;
				allAlleles++;
			}
			// output variable.
			if (m_vars.contains(string(AlleleNum_String) + "_sp"))
				pop.setIntVectorVar(subPopVar_String(*it, AlleleNum_String) + "{" + toStr(loc) + "}", alleles);
			if (m_vars.contains(string(AlleleFreq_String) + "_sp")) {
				vectorf freq(alleles.size(), 0.);
				if (allAlleles != 0) {
					for (size_t i = 0; i < alleles.size(); ++i)
						freq[i] = static_cast<double>(alleles[i]) / allAlleles;
				}
				pop.setDoubleVectorVar(subPopVar_String(*it, AlleleFreq_String) + "{" + toStr(loc) + "}", freq);
			}
			// total allele count
			if (alleleCnt[idx].size() < alleles.size())
				alleleCnt[idx].resize(alleles.size(), 0);
			for (size_t i = 0; i < alleles.size(); ++i)
				alleleCnt[idx][i] += alleles[i];
			allAllelesCnt[idx] += allAlleles;
		}
		if (it->isVirtual())
			pop.deactivateVirtualSubPop(it->subPop());
	}

	if (m_vars.contains(AlleleNum_String)) {
		pop.removeVar(AlleleNum_String);
		for (size_t idx = 0; idx < m_loci.size(); ++idx)
			pop.setIntVectorVar(string(AlleleNum_String) + "{" + toStr(m_loci[idx]) + "}",
				alleleCnt[idx]);
	}
	if (m_vars.contains(AlleleFreq_String)) {
		pop.removeVar(AlleleFreq_String);
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			vectorf freq(alleleCnt[idx].size(), 0.);
			if (allAllelesCnt[idx] != 0) {
				for (size_t i = 0; i < alleleCnt[idx].size(); ++i)
					freq[i] = static_cast<double>(alleleCnt[idx][i]) / allAllelesCnt[idx];
			}
			pop.setDoubleVectorVar(string(AlleleFreq_String) + "{" + toStr(m_loci[idx]) + "}",
				freq);
		}
	}
	return true;
}


vectori statAlleleFreq::numOfAlleles(population & pop)
{
	UINT maxLocus = 0;

	for (size_t loc = 0; loc < m_loci.size(); ++loc)
		if (maxLocus < m_loci[loc])
			maxLocus = m_loci[loc];
	vectori res(maxLocus + 1, 0);
	for (size_t loc = 0; loc < m_loci.size(); ++loc) {
		string varname = string(AlleleNum_String) + "{" + toStr(m_loci[loc]) + "}";
		PyObject * d = pop.getVar(varname);
		vectori num;
		PyObj_As_IntArray(d, num);

		int cnt = 0;
		for (size_t j = 0; j < num.size(); ++j)
			if (num[j] != 0)
				cnt += 1;
		res[m_loci[loc]] = cnt;
	}
	return res;
}


vectorf statAlleleFreq::alleleFreqVec(population & pop, int loc)
{
	string varname = string(AlleleFreq_String) + "{" + toStr(loc) + "}";
	PyObject * d = pop.getVar(varname);
	vectorf res;

	PyObj_As_Array(d, res);
	return res;
}


double statAlleleFreq::alleleFreq(population & pop, UINT allele, int loc)
{
	string varname = string(AlleleFreq_String) + "{" + toStr(loc) + "}[" + toStr(allele) + "]";
	PyObject * d = pop.getVar(varname);
	double af;

	PyObj_As_Double(d, af);
	return af;
}


vectorf statAlleleFreq::alleleFreqVec(population & pop, int loc, UINT subPop)
{
	string varname = subPopVar_String(subPop, AlleleFreq_String) + "{" + toStr(loc) + "}";
	PyObject * d = pop.getVar(varname);
	vectorf res;

	PyObj_As_Array(d, res);
	return res;
}


double statAlleleFreq::alleleFreq(population & pop, UINT allele, int loc, UINT subPop)
{
	string varname = subPopVar_String(subPop, AlleleFreq_String) + "{" + toStr(loc) + "}[" + toStr(allele) + "]";
	PyObject * d = pop.getVar(varname);
	double af;

	PyObj_As_Double(d, af);
	return af;
}


vectori statAlleleFreq::alleles(population & pop, int loc)
{
	string varname = string(AlleleNum_String) + "{" + toStr(loc) + "}";
	PyObject * d = pop.getVar(varname);
	vectori res;

	PyObj_As_IntArray(d, res);

	vectori al;

	for (size_t j = 0; j < res.size(); ++j) {
		if (res[j] > 0)
			al.push_back(j);
	}
	return al;
}


bool statHeteroFreq::apply(population & pop)
{
	if (m_atLoci.empty())
		return true;

	UINT numSP = pop.numSubPop();
	UINT numLoci = m_atLoci.size();

	pop.removeVar(HeteroNum_String);
	pop.removeVar(HeteroFreq_String);
	pop.removeVar(AllHeteroNum_String);
	pop.removeVar(AllHeteroFreq_String);
	pop.removeVar(HomoNum_String);
	pop.removeVar(HomoFreq_String);
	for (UINT sp = 0; sp < numSP; ++sp) {
		pop.removeVar(subPopVar_String(sp, HeteroNum_String));
		pop.removeVar(subPopVar_String(sp, HeteroFreq_String));
		pop.removeVar(subPopVar_String(sp, AllHeteroNum_String));
		pop.removeVar(subPopVar_String(sp, AllHeteroFreq_String));
		pop.removeVar(subPopVar_String(sp, HomoNum_String));
		pop.removeVar(subPopVar_String(sp, HomoFreq_String));
	}

	// may be resizing for different replicate of populations.
	// if not initialized or m_atLoci/numSP changes
	if (m_heteroNum.size() != (numSP + 1) * numLoci) {
		m_heteroNum.resize((numSP + 1) * numLoci);
		m_heteroFreq.resize((numSP + 1) * numLoci);
		m_homoNum.resize(numSP + 1);
		m_homoFreq.resize(numSP + 1);
	}

	string varname;
	ULONG popSize = pop.popSize();

	for (size_t i = 0; i < numLoci; ++i) {
		UINT loc = m_atLoci[i];
		DBG_DO(DBG_STATOR, cout << "Counting heterozygotes at locus " << loc << endl);

		vectori & sum = m_heteroNum[resIdx(i)];
		fill(sum.begin(), sum.end(), 0);
		int sumAll = 0;

		// for each subpopulation
		for (UINT sp = 0; sp < numSP;  ++sp) {
			vectori & num = m_heteroNum[resIdx(i, sp)];
			fill(num.begin(), num.end(), 0);
			int numAll = 0;

			// go through all alleles
			//?>> \todo here we assume diploid population
			// FIXME: alleleIterator will jump over chromosome x and Y etc.
			IndAlleleIterator a = pop.alleleIterator(loc, sp);
			for (; a.valid(); a += 2) {
				if (AlleleUnsigned(*a) >= num.size() )
					num.resize(*a + 1);

				if (AlleleUnsigned(*(a + 1)) >= num.size() )
					num.resize(*(a + 1) + 1);

				if (*a != *(a + 1) ) {
					num[*a]++;
					num[*(a + 1)]++;
					numAll++;
				}
			}

			sumAll += numAll;
			// add this number to overall num
			// calculate frequency
			if (numSP > 1) {
				if (sum.size() < num.size())
					sum.resize(num.size());
				for (size_t e = 0, eEnd = num.size(); e < eEnd; ++e)
					sum[e] += num[e];
			}

			vectorf & freq = m_heteroFreq[resIdx(i, sp)];
			freq.resize(num.size());
			for (size_t e = 0, eEnd = num.size(); e < eEnd; ++e)
				freq[e] = static_cast<double>(num[e]) / pop.subPopSize(sp);

			// set variable.
			if (m_ifPost[i] && m_postHetero) {
				varname = subPopVar_String(sp, HeteroNum_String) + "[" + toStr(loc) + "]";
				PyObject * d = pop.setIntVectorVar(varname, num);
				if (numSP == 1) {
					Py_INCREF(d);
					varname = toStr(HeteroNum_String) + "[" + toStr(loc) + "]";
					pop.setVar(varname, d);
				}

				varname = subPopVar_String(sp, HeteroFreq_String) + "[" + toStr(loc) + "]";
				d = pop.setDoubleVectorVar(varname, freq);

				if (numSP == 1) {
					Py_INCREF(d);
					varname = toStr(HeteroFreq_String) + "[" + toStr(loc) + "]";
					pop.setVar(varname, d);
				}

				// overall hetero
				varname = subPopVar_String(sp, AllHeteroNum_String) + "[" + toStr(loc) + "]";
				d = pop.setIntVar(varname, numAll);
				if (numSP == 1) {
					Py_INCREF(d);
					varname = toStr(AllHeteroNum_String) + "[" + toStr(loc) + "]";
					pop.setVar(varname, d);
				}

				varname = subPopVar_String(sp, AllHeteroFreq_String) + "[" + toStr(loc) + "]";
				d = pop.setDoubleVar(varname, static_cast<double>(numAll) / pop.subPopSize(sp));

				if (numSP == 1) {
					Py_INCREF(d);
					varname = toStr(AllHeteroFreq_String) + "[" + toStr(loc) + "]";
					pop.setVar(varname, d);
				}
			}
		}

		if (numSP > 1 && m_postHetero) {
			vectorf & freq = m_heteroFreq[resIdx(i)];
			freq.resize(sum.size());
			for (size_t e = 0, eEnd = sum.size(); e < eEnd; ++e)
				freq[e] = static_cast<double>(sum[e]) / popSize;

			if (m_ifPost[i]) {                                                        //
				varname = string(HeteroNum_String) + "[" + toStr(loc) + "]";
				pop.setIntVectorVar(varname, sum);

				varname = string(HeteroFreq_String) + "[" + toStr(loc) + "]";
				pop.setDoubleVectorVar(varname, freq);

				varname = string(AllHeteroNum_String) + "[" + toStr(loc) + "]";
				pop.setIntVar(varname, sumAll);

				varname = string(AllHeteroFreq_String) + "[" + toStr(loc) + "]";
				pop.setDoubleVar(varname, static_cast<double>(sumAll) / pop.popSize());
			}
		}                                                                                           // whole population
	}                                                                                               // for all loci

	if (m_postHomo) {
		for (size_t i = 0; i < numLoci; ++i) {
			UINT loc = m_atLoci[i];

			if (loc + 1 >= m_homoFreq[0].size()) {
				for (UINT sp = 0; sp < numSP + 1;  ++sp) {
					m_homoFreq[sp].resize(loc + 1, 0.0);
					m_homoNum[sp].resize(loc + 1, 0);
				}
			}

			// calculate homoNum
			for (UINT sp = 0; sp < numSP;  ++sp) {
				m_homoNum[sp][loc] = pop.subPopSize(sp) - m_heteroNum[resIdx(i, sp)][0];
				m_homoFreq[sp][loc] = (double)(m_homoNum[sp][loc]) / pop.subPopSize(sp);
			}
			m_homoNum[numSP][loc] = pop.popSize() - m_heteroNum[resIdx(i)][0];
			m_homoFreq[numSP][loc] = (double)(m_homoNum[numSP][loc]) / pop.popSize();
		}                                                                                 // all loci
		// post result
		for (UINT sp = 0; sp < numSP; ++sp) {
			pop.setIntVectorVar(subPopVar_String(sp, HomoNum_String),
				m_homoNum[sp]);
			pop.setDoubleVectorVar(subPopVar_String(sp, HomoFreq_String),
				m_homoFreq[sp]);
		}
		pop.setIntVectorVar(HomoNum_String, m_homoNum[numSP]);
		pop.setDoubleVectorVar(HomoFreq_String, m_homoFreq[numSP]);
	}
	return true;
}


bool statExpHetero::apply(population & pop)
{
	if (m_atLoci.empty())
		return true;

	UINT numSP = pop.numSubPop();
	UINT numLoci = m_atLoci.size();

	pop.removeVar(ExpHetero_String);
	for (UINT sp = 0; sp < numSP; ++sp)
		pop.removeVar(subPopVar_String(sp, ExpHetero_String));

	if (m_expHetero.size() != numSP + 1)
		m_expHetero.resize(numSP + 1);

	for (size_t i = 0; i < numLoci; ++i) {
		UINT loc = m_atLoci[i];

		for (UINT sp = 0; sp < numSP + 1;  ++sp) {
			if (m_expHetero[sp].size() < loc + 1)
				m_expHetero[sp].resize(loc + 1, 0.0);
		}

		// for each subpopulation
		for (UINT sp = 0; sp < numSP;  ++sp) {
			// calculate expected heterozygosity
			// get allele frequency
			vectorf af = m_alleleFreq.alleleFreqVec(pop, loc, sp);
			double expHeter = 1;
			// 1-sum pi^2
			for (int al = 0, alEnd = af.size() ; al < alEnd; al++)
				expHeter -= af[al] * af[al];

			m_expHetero[sp][loc] = expHeter;
		}

		vectorf af = m_alleleFreq.alleleFreqVec(pop, loc);
		double expHeter = 1;
		// 1-sum pi^2
		for (int al = 0, alEnd = af.size(); al < alEnd; al++)
			expHeter -= af[al] * af[al];

		m_expHetero[numSP][loc] = expHeter;
	}

	// post result
	for (UINT sp = 0; sp < numSP; ++sp)
		pop.setDoubleVectorVar(subPopVar_String(sp, ExpHetero_String),
			m_expHetero[sp]);
	pop.setDoubleVectorVar(ExpHetero_String,
		m_expHetero[numSP]);

	return true;
}


statGenoFreq::statGenoFreq(const vectorlu & genoFreq,
	const strDict & param)
	: m_atLoci(genoFreq), m_phase(false)
{
	if (!param.empty()) {
		strDict::const_iterator it;
		strDict::const_iterator itEnd = param.end();
		if ((it = param.find("phase")) != itEnd)
			m_phase = it->second != 0.;
	}
}


bool statGenoFreq::apply(population & pop)
{
	if (m_atLoci.empty())
		return true;

	DBG_ASSERT(pop.ploidy() == 2, ValueError,
		"Genotype can only be calcualted for diploid populations.");

	UINT numSP = pop.numSubPop();
	ULONG popSize = pop.popSize();

	for (size_t i = 0, iEnd = m_atLoci.size(); i < iEnd;  ++i) {
		if (static_cast<UINT>(m_atLoci[i]) >= pop.totNumLoci() )
			throw IndexError("Absolute locus index "
				+ toStr(m_atLoci[i]) + " is out of range of 0 ~ "
				+ toStr(pop.totNumLoci() - 1));
	}

	// first remove genoNum that may be set by previous count.
	// for example if genoNum[a][10] was set but this time there
	// is no allele 10, we will not reset genoNum[a][10] ...
	pop.removeVar(GenotypeNum_String);
	pop.removeVar(GenotypeFreq_String);
	for (UINT sp = 0; sp < numSP;  ++sp) {
		// remove genoNum, genoFreq that may be set by previous run.
		pop.removeVar(subPopVar_String(sp, GenotypeNum_String));
		pop.removeVar(subPopVar_String(sp, GenotypeFreq_String));
	}

	string varname;

	// deal with genotype
	for (size_t i = 0, iEnd = m_atLoci.size(); i < iEnd;  ++i) {
		// for each locus, we need to use a vector of dictionaries.
		vector<intDict> sum;

		UINT loc = m_atLoci[i];

#ifndef BINARYALLELE
		Allele a, b;
#else
		unsigned short a, b;
#endif

		// for each subpopulation
		for (UINT sp = 0; sp < numSP;  ++sp) {
			DBG_DO(DBG_STATOR, cout << "Counting genotypes at locus " <<
				loc << " subPop " << sp << endl);

			vector<intDict> num;

			// go through a single allele for all individual, all diploid
			IndIterator it = pop.indIterator(sp);
			for (; it.valid(); ++it) {
				a = it->allele(loc, 0);
				b = it->allele(loc, 1);
				if (!m_phase && a > b)
					std::swap(a, b);

				if (a >= num.size() )
					num.resize(a + 1);

				num[a][b]++;

				if (a >= sum.size() )
					sum.resize(a + 1);

				sum[a][b]++;
			}

			// register values for this subpopulation
			for (a = 0; a < num.size(); ++a) {
				// need to replace previous values
				// if( num[a].empty() )
				//   continue;

				// empty dictionary should be allowed
				varname = subPopVar_String(sp, GenotypeNum_String) +
				          + "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
				pop.setIntDictVar(varname, num[a]);

				// apply frequency
				for (intDict::iterator it = num[a].begin(), itEnd = num[a].end(); it != itEnd; ++it)
					it->second = it->second / pop.subPopSize(sp);

				varname = subPopVar_String(sp, GenotypeFreq_String) +
				          + "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
				pop.setIntDictVar(varname, num[a]);
			}
		}

		for (a = 0; a < sum.size(); ++a) {
			// if( sum[a].empty() )
			//  continue;

			// empty dictionary should be allowed
			varname = toStr(GenotypeNum_String) + "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
			pop.setIntDictVar(varname, sum[a]);

			// apply frequency
			for (intDict::iterator it = sum[a].begin(), itEnd = sum[a].end(); it != itEnd; ++it)
				it->second = it->second / popSize;

			varname = toStr(GenotypeFreq_String) + "[" + toStr(loc) + "][" + toStr(int(a)) + "]";
			pop.setIntDictVar(varname, sum[a]);
		}
	}
	return true;
}


vectorlu statGenoFreq::countGenotype(population & pop, UINT loc, UINT wildtype)
{
	DBG_ASSERT(pop.ploidy() == 2, ValueError,
		"Genotype can only be calcualted for diploid populations.");

	DBG_FAILIF(loc >= pop.totNumLoci(), IndexError,
		"Locus index out of range");

	vectorlu res(3, 0);

	IndIterator it = pop.indIterator();
	for (; it.valid(); ++it) {
		Allele a = ToAllele(it->allele(loc, 0));
		Allele b = ToAllele(it->allele(loc, 1));
		if (a == wildtype) {
			if (b == wildtype)
				res[0]++;
			else
				res[1]++;
		} else {
			if (b == wildtype)
				res[1]++;
			else
				res[2]++;
		}
	}
	return res;
}


vectorlu statGenoFreq::countGenotype(population & pop, UINT loc, SubPopID subPop, UINT wildtype)
{
	DBG_ASSERT(pop.ploidy() == 2, ValueError,
		"Genotype can only be calcualted for diploid populations.");

	DBG_FAILIF(loc >= pop.totNumLoci(), IndexError,
		"Locus index out of range");

	vectorlu res(3, 0);

	IndIterator it = pop.indIterator(subPop);
	for (; it.valid(); ++it) {
		Allele a = ToAllele(it->allele(loc, 0));
		Allele b = ToAllele(it->allele(loc, 1));
		if (a == wildtype) {
			if (b == wildtype)
				res[0]++;
			else
				res[1]++;
		} else {
			if (b == wildtype)
				res[1]++;
			else
				res[2]++;
		}
	}
	return res;
}


int statHaploFreq::haploIndex(const vectori & haplo)
{
	// first locate haplo
	UINT idx = 0;

	while (m_haplotypes[idx] != haplo && idx < m_haplotypes.size())
		idx++;

	DBG_ASSERT(idx != m_haplotypes.size(), ValueError,
		"Can not find haplotype " + toStr(haplo[0]) + ", " + toStr(haplo[1]));
	return idx;
}


void statHaploFreq::addHaplotype(const vectori & haplo, bool post)
{
	intMatrix::iterator it;

	if ((it = find(m_haplotypes.begin(), m_haplotypes.end(), haplo)) == m_haplotypes.end()) {
		m_haplotypes.push_back(haplo);
		m_ifPost.push_back(static_cast<int>(post));
	} else
		m_ifPost[it - m_haplotypes.begin()] |= static_cast<int>(post);
}


bool statHaploFreq::apply(population & pop)
{
	if (m_haplotypes.empty())
		return true;

	pop.removeVar(HaplotypeNum_String);
	pop.removeVar(HaplotypeFreq_String);
	for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
		pop.removeVar(subPopVar_String(sp, HaplotypeNum_String));
		pop.removeVar(subPopVar_String(sp, HaplotypeFreq_String));
	}

	UINT nHap = m_haplotypes.size();

	// first time?
	if (m_haploNum.size() != (pop.numSubPop() + 1) * nHap) {
		for (size_t i = 0; i < nHap;  ++i) {
			vectori & haplotype = m_haplotypes[i];
			size_t sz = haplotype.size();

			if (sz == 0)
				throw ValueError("has to specify some haplotype.");
			else if (sz == 1)
				throw ValueError("Haplotype must contain alleles from more than one locus.");
		}
		m_haploNum.resize( (pop.numSubPop() + 1) * nHap);
		m_haploFreq.resize( (pop.numSubPop() + 1) * nHap);
	}

	UINT numSP = pop.numSubPop();

	DBG_DO(DBG_STATOR, cout << "Counting haplotypes" << endl);

	// clear all statistics
	for (size_t h = 0; h < nHap * (numSP + 1); ++h) {
		m_haploNum[h].clear();
		m_haploFreq[h].clear();
	}

	// for each subpopulation
	for (UINT sp = 0; sp < numSP;  ++sp) {
		for (size_t h = 0; h < nHap; ++h) {
			vectori & haplotype = m_haplotypes[h];
			map< vectori, UINT> & count = m_haploNum[h + sp * nHap];
			map< vectori, UINT> & sum = m_haploNum[h + numSP * nHap];

			size_t sz = haplotype.size();

			vectori sampleHap(sz);

			IndAlleleIterator it = pop.alleleIterator(0, sp);
			for (; it.valid(); ++it) {
				for (size_t hap = 0; hap < sz; ++hap)
					sampleHap[hap] = *(it.ptr() + haplotype[hap]);

				// add sampleHap count
				count[sampleHap]++;
				sum[sampleHap]++;
			}
		}
	}
	// finish count

	// calculate haploFreq,
	// for each subpopulation
	for (UINT sp = 0; sp < numSP;  ++sp) {
		// record both num and freq
		string varNumName = subPopVar_String(sp, HaplotypeNum_String);
		string varFreqName = subPopVar_String(sp, HaplotypeFreq_String);
		for (size_t h = 0; h < nHap; ++h) {
			vectori & haplotype = m_haplotypes[h];
			map< vectori, UINT> & count = m_haploNum[h + sp * nHap];
			map< vectori, double> & freq = m_haploFreq[h + sp * nHap];

			for (map<vectori, UINT>::iterator it = count.begin(); it != count.end(); ++it) {
				freq[ it->first] = double(it->second) / (pop.subPopSize(sp) * pop.ploidy());
				if (m_ifPost[h]) {
					pop.setIntVar(varNumName + haploKey(haplotype) + haploKey(it->first), it->second);
					pop.setDoubleVar(varFreqName + haploKey(haplotype) + haploKey(it->first),
						double(it->second) / (pop.subPopSize(sp) * pop.ploidy()));
				}
			}
		}
	}
	// whole population
	for (size_t h = 0; h < nHap; ++h) {
		vectori & haplotype = m_haplotypes[h];
		map< vectori, UINT> & count = m_haploNum[h + numSP * nHap];
		map< vectori, double> & freq = m_haploFreq[h + numSP * nHap];

		for (map<vectori, UINT>::iterator it = count.begin(); it != count.end(); ++it) {
			freq[ it->first] = double(it->second) / (pop.popSize() * pop.ploidy());
			if (m_ifPost[h]) {
				pop.setIntVar(HaplotypeNum_String + haploKey(haplotype) + haploKey(it->first), it->second);
				pop.setDoubleVar(HaplotypeFreq_String + haploKey(haplotype) + haploKey(it->first),
					double(it->second) / (pop.popSize() * pop.ploidy()));
			}
		}
	}
	return true;
}


statLD::statLD(statAlleleFreq & alleleFreq, statHaploFreq & haploFreq,
	const intMatrix & LD, const strDict & param)
	: m_alleleFreq(alleleFreq), m_haploFreq(haploFreq),
	m_LD(LD),
	// default values,
	m_midValues(false),
	m_evalInSubPop(true),
	m_output_ld(true),
	m_output_ld_prime(true),
	m_output_r2(true),
	m_output_delta2(true),
	m_output_LD(true),
	m_output_LD_prime(true),
	m_output_R2(true),
	m_output_Delta2(true),
	m_output_ChiSq(false),
	m_output_UCU(false),
	m_output_CramerV(false)
{
	// parameters
	if (!param.empty()) {
		strDict::const_iterator it;
		strDict::const_iterator itEnd = param.end();
		if ((it = param.find("subPop")) != itEnd)
			m_evalInSubPop = it->second != 0.;
		if ((it = param.find("midValues")) != itEnd)
			m_midValues = it->second != 0.;
		// if any statistics is specified, other unspecified ones are not calculated
		if (param.find(LD_String) != itEnd ||
		    param.find(LDPRIME_String) != itEnd ||
		    param.find(R2_String) != itEnd ||
		    param.find(DELTA2_String) != itEnd ||
		    param.find(AvgLD_String) != itEnd ||
		    param.find(AvgLDPRIME_String) != itEnd ||
		    param.find(AvgR2_String) != itEnd ||
		    param.find(AvgDELTA2_String) != itEnd ||
		    param.find(ChiSq_String) != itEnd ||
		    param.find(UCU_String) != itEnd ||
		    param.find(CramerV_String) != itEnd) {
			m_output_ld = false;
			m_output_ld_prime = false;
			m_output_r2 = false;
			m_output_delta2 = false;
			m_output_LD = false;
			m_output_LD_prime = false;
			m_output_R2 = false;
			m_output_Delta2 = false;
			m_output_ChiSq = false;
			m_output_UCU = false;
			m_output_CramerV = false;
			// if has key, and is True or 1
			if ((it = param.find(LD_String)) != itEnd)
				m_output_ld = it->second != 0.;
			if ((it = param.find(LDPRIME_String)) != itEnd)
				m_output_ld_prime = it->second != 0.;
			if ((it = param.find(R2_String)) != itEnd)
				m_output_r2 = it->second != 0.;
			if ((it = param.find(DELTA2_String)) != itEnd)
				m_output_delta2 = it->second != 0.;
			if ((it = param.find(AvgLD_String)) != itEnd)
				m_output_LD = it->second != 0.;
			if ((it = param.find(AvgLDPRIME_String)) != itEnd)
				m_output_LD_prime = it->second != 0.;
			if ((it = param.find(AvgR2_String)) != itEnd)
				m_output_R2 = it->second != 0.;
			if ((it = param.find(AvgDELTA2_String)) != itEnd)
				m_output_Delta2 = it->second != 0.;
			if ((it = param.find(ChiSq_String)) != itEnd)
				m_output_ChiSq = it->second != 0.;
			if ((it = param.find(UCU_String)) != itEnd)
				m_output_UCU = it->second != 0.;
			if ((it = param.find(CramerV_String)) != itEnd)
				m_output_CramerV = it->second != 0.;
		}
	}
	//
	for (size_t i = 0, iEnd = m_LD.size(); i < iEnd; ++i) {
		// these asserts will only be checked in non-optimized modules
		DBG_FAILIF(m_LD[i].size() != 2 && m_LD[i].size() != 4,
			ValueError, "Expecting [locus locus [allele allele ]] items");

		// midValues is used to tell alleleFreq that the calculated allele
		// frequency values should not be posted to pop.dvars()
		//
		// That is to say,
		//     stat(LD=[0,1])
		// will not generate
		//     pop.dvars().alleleFreq
		// unless stat() is called as
		//     stat(LD=[0,1], midValues=True)
		//
		m_alleleFreq.addLocus(m_LD[i][0], AllSubPops, stringList("alleleNum"));
		m_alleleFreq.addLocus(m_LD[i][1], AllSubPops, stringList("alleleNum"));
		// also need haplotype.
		if (m_LD[i][0] != m_LD[i][1]) {
			vectori hap(2);
			hap[0] = m_LD[i][0];
			hap[1] = m_LD[i][1];
			m_haploFreq.addHaplotype(hap, m_midValues);
		}
	}
}


// this function calculate single-allele LD measures
// D, D_p and r2 are used to return calculated values.
// LD for subpopulation sp is calculated if subPop is true
void statLD::calculateLD(population & pop, const vectori & hapLoci, const vectori & hapAlleles, UINT sp, bool subPop,
                         double & P_A, double & P_B, double & D, double & D_prime, double & r2, double & delta2)
{
	if (subPop) {
		// get haplotype freq from the m_haploFreq object
		double P_AB;
		if (hapLoci[0] == hapLoci[1])
			P_AB = 0;
		else
			P_AB = m_haploFreq.haploFreq(hapLoci, sp)[hapAlleles];
		// get allele freq from the m_alleleFreq object
		P_A = m_alleleFreq.alleleFreq(pop, hapAlleles[0], hapLoci[0], sp);
		P_B = m_alleleFreq.alleleFreq(pop, hapAlleles[1], hapLoci[1], sp);

		// calculate LD
		D = P_AB - P_A * P_B;
		// calculate LD'
		double D_max = D > 0 ? std::min(P_A * (1 - P_B), (1 - P_A) * P_B) : std::min(P_A * P_B, (1 - P_A) * (1 - P_B));
		// fcmp_eq is the float comparison operator, which treat (-1e-10, 1e-10) or so as 0 (platform dependent)
		D_prime = fcmp_eq(D_max, 0.) ? 0. : D / D_max;
		r2 = (fcmp_eq(P_A, 0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1)) ? 0. : D * D / P_A / (1 - P_A) / P_B / (1 - P_B);
		// calculate delta2
		delta2 = (fcmp_eq(P_A, 0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1)) ? 0. : pow((P_AB * ((1 - P_A) - (P_B - P_AB)) - (P_A - P_AB) * (P_B - P_AB)), 2) / (P_A * (1 - P_A) * P_B * (1 - P_B));
		// if environmental variable SIMUDEBUG is set to DBG_STATOR, or
		// if TurnOnDebug(DBG_STATOR) is called in python, the following will be printed.
		DBG_DO(DBG_STATOR, cout << "LD: subpop " << sp << " : P_AB: " << P_AB
			                    << " P_A: " << P_A << " P_B: " << P_B << " D_max: " << D_max <<
			" LD: " << D << " LD': " << D_prime << " r2: " << r2 << " delta2: " << delta2 << endl);
	} else {
		// whole population
		// get haplotype freq
		double P_AB;
		if (hapLoci[0] == hapLoci[1])
			P_AB = 0;
		else
			P_AB = m_haploFreq.haploFreq(hapLoci)[hapAlleles];
		P_A = m_alleleFreq.alleleFreq(pop, hapAlleles[0], hapLoci[0]);
		P_B = m_alleleFreq.alleleFreq(pop, hapAlleles[1], hapLoci[1]);

		// calculate LD
		D = P_AB - P_A * P_B;
		// calculate LD'
		double D_max = D > 0 ? std::min(P_A * (1 - P_B), (1 - P_A) * P_B) : std::min(P_A * P_B, (1 - P_A) * (1 - P_B));
		D_prime = fcmp_eq(D_max, 0) ? 0 : D / D_max;
		r2 = (fcmp_eq(P_A, 0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1)) ? 0 : D * D / P_A / (1 - P_A) / P_B / (1 - P_B);
		delta2 = (fcmp_eq(P_A, 0) || fcmp_eq(P_B, 0) || fcmp_eq(P_A, 1) || fcmp_eq(P_B, 1)) ? 0. : pow((P_AB * ((1 - P_A) - (P_B - P_AB)) - (P_A - P_AB) * (P_B - P_AB)), 2) / (P_A * (1 - P_A) * P_B * (1 - P_B));

		DBG_DO(DBG_STATOR, cout << "LD: P_AB: " << P_AB
			                    << " P_A: " << P_A << " P_B: " << P_B << " D_max: " << D_max <<
			" LD: " << D << " LD': " << D_prime << " r2: " << r2 << " delta2: " << delta2 << endl);
	}
}


// try to shorten statLD::apply
void statLD::outputLD(population & pop, const vectori & hapLoci, const string & allele_string, UINT sp, bool subPop,
                      bool valid_delta2, double D, double D_prime, double r2, double delta2)
{
	string ld_name, ldp_name, r2_name, d2_name, key_name;
	bool ld_cond, ldp_cond, r2_cond, d2_cond;

	// single allele cases
	if (allele_string != "") {
		ld_name = LD_String;
		ldp_name = LDPRIME_String;
		r2_name = R2_String;
		d2_name = DELTA2_String;
		key_name = haploKey(hapLoci) + allele_string;
		ld_cond = m_output_ld;
		ldp_cond = m_output_ld_prime;
		r2_cond = m_output_r2;
		d2_cond = m_output_delta2;
	} else {
		ld_name = AvgLD_String;
		ldp_name = AvgLDPRIME_String;
		r2_name = AvgR2_String;
		d2_name = AvgDELTA2_String;
		key_name = "{" + toStr(hapLoci[0]) + "}{" + toStr(hapLoci[1]) + '}';
		ld_cond = m_output_LD;
		ldp_cond = m_output_LD_prime;
		r2_cond = m_output_R2;
		d2_cond = m_output_Delta2;
	}
	if (subPop) {
		ld_name = subPopVar_String(sp, ld_name);
		ldp_name = subPopVar_String(sp, ldp_name);
		r2_name = subPopVar_String(sp, r2_name);
		d2_name = subPopVar_String(sp, d2_name);
	}
	DBG_DO(DBG_STATOR, cout << "Output statistics " << ldp_name + key_name << endl);
	if (ld_cond)
		pop.setDoubleVar(ld_name + key_name, D);
	if (ldp_cond)
		pop.setDoubleVar(ldp_name + key_name, D_prime);
	if (r2_cond)
		pop.setDoubleVar(r2_name + key_name, r2);
	if (valid_delta2 && d2_cond)
		pop.setDoubleVar(d2_name + key_name, delta2);
}


// this function is called by stat::apply(pop). It is called
// after m_alleleFreq.apply(pop) and m_haploFreq.apply(pop) so
// allele and haplotype frequencies should be available.
bool statLD::apply(population & pop)
{
	if (m_LD.empty())
		return true;

	UINT numSP = pop.numSubPop();
	UINT nLD = m_LD.size();
	// used for delta2 which can only be computed for 2 alleles
	vectori numofalleles = m_alleleFreq.numOfAlleles(pop);
	bool valid_delta2 = false;

	// remove previous values.
	pop.removeVar(LD_String);
	pop.removeVar(LDPRIME_String);
	pop.removeVar(R2_String);
	pop.removeVar(DELTA2_String);
	pop.removeVar(AvgLD_String);
	pop.removeVar(AvgLDPRIME_String);
	pop.removeVar(AvgR2_String);
	pop.removeVar(AvgDELTA2_String);
	pop.removeVar(ChiSq_String);
	pop.removeVar(ChiSq_P_String);
	pop.removeVar(UCU_String);
	pop.removeVar(CramerV_String);
	// also vars at each subpopulations
	for (UINT sp = 0; sp < numSP;  ++sp) {
		// subPopVar_String is nothing but subPop[sp]['string']
		pop.removeVar(subPopVar_String(sp, LD_String));
		pop.removeVar(subPopVar_String(sp, LDPRIME_String));
		pop.removeVar(subPopVar_String(sp, R2_String));
		pop.removeVar(subPopVar_String(sp, DELTA2_String));
		pop.removeVar(subPopVar_String(sp, AvgLD_String));
		pop.removeVar(subPopVar_String(sp, AvgLDPRIME_String));
		pop.removeVar(subPopVar_String(sp, AvgR2_String));
		pop.removeVar(subPopVar_String(sp, AvgDELTA2_String));
		pop.removeVar(subPopVar_String(sp, ChiSq_String));
		pop.removeVar(subPopVar_String(sp, ChiSq_P_String));
		pop.removeVar(subPopVar_String(sp, UCU_String));
		pop.removeVar(subPopVar_String(sp, CramerV_String));
	}
	for (size_t i = 0; i < nLD; ++i) {
		// specifying alleles
		if (m_LD[i].size() == 4) {
			vectori hapLoci(2);
			vectori hapAlleles(2);

			hapLoci[0] = m_LD[i][0];
			hapLoci[1] = m_LD[i][1];
			if (numofalleles[hapLoci[0]] <= 2 && numofalleles[hapLoci[1]] <= 2)
				valid_delta2 = true;
			hapAlleles[0] = m_LD[i][2];
			hapAlleles[1] = m_LD[i][3];

			// whole population, P_A, P_B is ignored
			double D = 0;
			double D_prime = 0;
			double r2 = 0;
			double delta2 = 0;
			double P_A = 0;
			double P_B = 0;
			calculateLD(pop, hapLoci, hapAlleles, 0, false, P_A, P_B, D, D_prime, r2, delta2);
			outputLD(pop, hapLoci, haploKey(hapAlleles), 0, false, valid_delta2, D, D_prime, r2, delta2);

			if (m_evalInSubPop) {
				if (numSP == 1)     // use the whole population result
					outputLD(pop, hapLoci, haploKey(hapAlleles), 0, true, valid_delta2, D, D_prime, r2, delta2);
				else {
					for (UINT sp = 0; sp < numSP;  ++sp) {
						// get LD values, P_A, P_B is ignored.
						double D = 0;
						double D_prime = 0;
						double r2 = 0;
						double delta2 = 0;
						double P_A = 0;
						double P_B = 0;
						calculateLD(pop, hapLoci, hapAlleles, sp, true, P_A, P_B, D, D_prime, r2, delta2);
						outputLD(pop, hapLoci, haploKey(hapAlleles), sp, true, valid_delta2, D, D_prime, r2, delta2);
					}
				}
			}                                                                         // eval in subpop
		} else {
			// No alleles specified, average over all available alleles
			vectori hapLoci(2);
			vectori hapAlleles(2);

			hapLoci[0] = m_LD[i][0];
			hapLoci[1] = m_LD[i][1];
			if (numofalleles[hapLoci[0]] <= 2 && numofalleles[hapLoci[1]] <= 2)
				valid_delta2 = true;

			// find out all alleles
			vectori A_alleles = m_alleleFreq.alleles(pop, hapLoci[0]);
			vectori B_alleles = m_alleleFreq.alleles(pop, hapLoci[1]);

			// whole population
			double D = 0.0, D_prime = 0.0, r2 = 0.0, delta2 = 0.0;
			for (vectori::iterator A_ale = A_alleles.begin();
			     A_ale != A_alleles.end(); ++A_ale) {
				for (vectori::iterator B_ale = B_alleles.begin();
				     B_ale != B_alleles.end(); ++B_ale) {
					hapAlleles[0] = *A_ale;
					hapAlleles[1] = *B_ale;
					double D_ = 0;
					double D_prime_ = 0;
					double r2_ = 0;
					double delta2_ = 0;
					double P_A = 0;
					double P_B = 0;
					calculateLD(pop, hapLoci, hapAlleles, 0, false, P_A, P_B, D_, D_prime_, r2_, delta2_);
					if (m_midValues)
						outputLD(pop, hapLoci, haploKey(hapAlleles), 0, false, valid_delta2, D_, D_prime_, r2_, delta2_);

					D += P_A * P_B * fabs(D_);
					D_prime += P_A * P_B * fabs(D_prime_);
					r2 += P_A * P_B * r2_;
					if (valid_delta2)
						delta2 += P_A * P_B * delta2_;
					DBG_DO(DBG_STATOR, cout << "Sum D " << D << " D' " << D_prime << " r2 " << r2 << endl);
				}                                                                 // all haplotypes
			}
			outputLD(pop, hapLoci, "", 0, false, valid_delta2, D, D_prime, r2, delta2);

			if (m_evalInSubPop) {
				if (numSP == 1)     // use the whole population result
					outputLD(pop, hapLoci, "", 0, true, valid_delta2, D, D_prime, r2, delta2);
				else {
					for (UINT sp = 0; sp < numSP;  ++sp) {
						double D = 0.0, D_prime = 0.0, r2 = 0.0, delta2 = 0.0;
						// iterate through all alleles at locus A and B
						for (vectori::iterator A_ale = A_alleles.begin();
						     A_ale != A_alleles.end(); ++A_ale) {
							for (vectori::iterator B_ale = B_alleles.begin();
							     B_ale != B_alleles.end(); ++B_ale) {
								// this is now single allele ...
								hapAlleles[0] = *A_ale;
								hapAlleles[1] = *B_ale;
								double D_ = 0;
								double D_prime_ = 0;
								double r2_ = 0;
								double delta2_ = 0;
								double P_A = 0;
								double P_B = 0;
								calculateLD(pop, hapLoci, hapAlleles, sp, true, P_A, P_B, D_, D_prime_, r2_, delta2_);

								// store allele-specific LD values as well.
								if (m_midValues)
									outputLD(pop, hapLoci, haploKey(hapAlleles), sp, true, valid_delta2, D_, D_prime_, r2_, delta2_);

								D += P_A * P_B * fabs(D_);
								D_prime += P_A * P_B * fabs(D_prime_);
								r2 += P_A * P_B * r2_;
								if (valid_delta2)
									delta2 += P_A * P_B * delta2_;
								DBG_DO(DBG_STATOR, cout << "Sum D " << D << " D' " << D_prime << " r2 " << r2 << "delta2" << delta2 << endl);
							}
						}
						outputLD(pop, hapLoci, "", sp, true, valid_delta2, D, D_prime, r2, delta2);
					}                                                                               // each subpop
				}                                                                                   // numSP
			}                                                                                       // eval in subpop
		}                                                                                           // size = 2, 4
	}                                                                                               // for all LD
	/* FIXME: the following is moved from the old association class... they should
	 * be merged to the above calculation
	 */
	if (!m_output_ChiSq && !m_output_UCU && !m_output_CramerV)
		return true;

	for (size_t i = 0; i < m_LD.size(); ++i) {
		//
		vectori hapLoci = m_LD[i];
		hapLoci.resize(2);

		// find out all alleles
		vectori A_alleles = m_alleleFreq.alleles(pop, hapLoci[0]);
		vectori B_alleles = m_alleleFreq.alleles(pop, hapLoci[1]);
		string hapLociStr = '[' + toStr(hapLoci[0]) + "][" +
		                    toStr(hapLoci[1]) + ']';

		UINT as = A_alleles.size();
		UINT bs = B_alleles.size();

		// the whole population
		vector<vectorf> cont_table(as + 1);
		for (size_t i = 0; i <= as; ++i)
			cont_table[i].resize(bs + 1);
		double n = static_cast<double>(pop.popSize());
		double ChiSq = 0.0, ChiSq_P = 0.0, UC_U = 0.0, CramerV = 0.0;
		double HA = 0.0, HB = 0.0, HAB = 0.0;
		// initialize last line/column
		for (size_t i = 0; i < as; ++i)
			cont_table[i][bs] = 0;
		for (size_t j = 0; j <= bs; ++j)
			cont_table[as][j] = 0;
		// get P_ij
		for (size_t i = 0; i < as; ++i) {
			for (size_t j = 0; j < bs; ++j) {
				vectori hapAlleles(2);
				hapAlleles[0] = A_alleles[i];
				hapAlleles[1] = B_alleles[j];
				cont_table[i][j] = m_haploFreq.haploFreq(hapLoci)[hapAlleles];
				cont_table[i][bs] += cont_table[i][j];
				cont_table[as][j] += cont_table[i][j];
				cont_table[as][bs] += cont_table[i][j];
			}
		}
		DBG_ASSERT(fcmp_eq(cont_table[as][bs], 1.), ValueError,
			"Sum of haplotype frequencies is not 1. Association will not be computed.");
		//DBG_DO(DBG_STATOR, for(size_t i=0; i <= as; ++i) cout << cont_table[i] << endl);
		// calculate statistics
		//ChiSq
		if (m_output_ChiSq || m_output_CramerV) {
			for (size_t i = 0; i < as; ++i)
				for (size_t j = 0; j < bs; ++j)
					ChiSq += pow((n * cont_table[i][j] - n * cont_table[i][bs] * cont_table[as][j]), 2)
					         / (n * cont_table[i][bs] * cont_table[as][j]);
			ChiSq_P = GetRNG().pvalChiSq(ChiSq, (as - 1) * (bs - 1));
		}
		if (m_output_ChiSq) {
			pop.setDoubleVar(ChiSq_String + hapLociStr, ChiSq);
			pop.setDoubleVar(ChiSq_P_String + hapLociStr, ChiSq_P);
		}
		//UC_U
		if (m_output_UCU) {
			for (size_t i = 0; i < as; ++i)
				HA += -cont_table[i][bs] * log(cont_table[i][bs]);
			for (size_t j = 0; j <= bs; ++j)
				HB += -cont_table[as][j] * log(cont_table[as][j]);
			for (size_t i = 0; i < as; ++i)
				for (size_t j = 0; j < bs; ++j)
					HAB += -cont_table[i][j] * log(cont_table[i][j]);
			UC_U = 2 * ((HA + HB - HAB) / (HA + HB));
			pop.setDoubleVar(UCU_String + hapLociStr, UC_U);
		}
		//CramerV
		if (m_output_CramerV) {
			CramerV = sqrt(ChiSq / (n * std::min(as - 1, bs - 1)));
			pop.setDoubleVar(CramerV_String + hapLociStr, CramerV);
		}

		if (!m_evalInSubPop)
			continue;
		// subpopulations ...
		if (numSP == 1) {
			pop.setDoubleVar(subPopVar_String(0, ChiSq_String) + hapLociStr, ChiSq);
			pop.setDoubleVar(subPopVar_String(0, ChiSq_P_String) + hapLociStr, ChiSq_P);
			pop.setDoubleVar(subPopVar_String(0, UCU_String) + hapLociStr, UC_U);
			pop.setDoubleVar(subPopVar_String(0, CramerV_String) + hapLociStr, CramerV);
		} else {
			for (UINT sp = 0; sp < numSP;  ++sp) {
				UINT as = A_alleles.size();
				UINT bs = B_alleles.size();
				vector<vectorf> cont_table(as + 1);
				for (size_t i = 0; i <= as; ++i)
					cont_table[i].resize(bs + 1);
				double n = static_cast<double>(pop.subPopSize(sp));
				double ChiSq = 0.0, ChiSq_P = 0.0, UC_U = 0.0, CramerV = 0.0;
				double HA = 0.0, HB = 0.0, HAB = 0.0;
				// initialize last line/column
				for (size_t i = 0; i < as; ++i)
					cont_table[i][bs] = 0;
				for (size_t j = 0; j <= bs; ++j)
					cont_table[as][j] = 0;
				// get P_ij
				for (size_t i = 0; i < as; ++i)
					for (size_t j = 0; j < bs; ++j) {
						vectori hapAlleles(2);
						hapAlleles[0] = A_alleles[i];
						hapAlleles[1] = B_alleles[j];
						cont_table[i][j] = m_haploFreq.haploFreq(hapLoci, sp)[hapAlleles];
						cont_table[i][bs] += cont_table[i][j];
						cont_table[as][j] += cont_table[i][j];
						cont_table[as][bs] += cont_table[i][j];
					}
				DBG_ASSERT(fcmp_eq(cont_table[as][bs], 1.), ValueError,
					"Sum of haplotype frequencies is not 1. Association will not be computed.");
				//DBG_DO(DBG_STATOR, for(size_t i=0; i <= as; ++i) cout << cont_table[i] << endl);
				// calculate statistics
				//ChiSq
				if (m_output_ChiSq || m_output_CramerV) {
					for (size_t i = 0; i < as; ++i)
						for (size_t j = 0; j < bs; ++j)
							ChiSq += pow((n * cont_table[i][j] - n * cont_table[i][bs] * cont_table[as][j]), 2)
							         / (n * cont_table[i][bs] * cont_table[as][j]);
					ChiSq_P = GetRNG().pvalChiSq(ChiSq, (as - 1) * (bs - 1));
				}
				if (m_output_ChiSq) {
					pop.setDoubleVar(subPopVar_String(sp, ChiSq_String) + hapLociStr, ChiSq);
					pop.setDoubleVar(subPopVar_String(sp, ChiSq_P_String) + hapLociStr, ChiSq_P);
				}
				//UC_U
				if (m_output_UCU) {
					for (size_t i = 0; i < as; ++i)
						HA += -cont_table[i][bs] * log(cont_table[i][bs]);
					for (size_t j = 0; j <= bs; ++j)
						HB += -cont_table[as][j] * log(cont_table[as][j]);
					for (size_t i = 0; i < as; ++i)
						for (size_t j = 0; j < bs; ++j)
							HAB += -cont_table[i][j] * log(cont_table[i][j]);
					UC_U = 2 * ((HA + HB - HAB) / (HA + HB));
					pop.setDoubleVar(subPopVar_String(sp, UCU_String) + hapLociStr, UC_U);
				}
				//CramerV
				if (m_output_CramerV) {
					CramerV = sqrt(ChiSq / (n * std::min(as - 1, bs - 1)));
					pop.setDoubleVar(subPopVar_String(sp, CramerV_String) + hapLociStr, CramerV);
				}
			}                                                                                       // each sp
		}                                                                                           // numSP > 1
	}
	return true;
}


statAssociation::statAssociation(const vectorlu & loci, const subPopList & subPops) :
	m_loci(loci), m_subPops(subPops)
{
}


void statAssociation::calcChiSq(ULONG aff_0, ULONG aff_1, ULONG unaff_0, ULONG unaff_1,
                                double & chisq, double & pvalue)
{
	double cases = aff_0 + aff_1;
	double controls = unaff_0 + unaff_1;
	double total = cases + controls;
	double allele0 = aff_0 + unaff_0;
	double allele1 = aff_1 + unaff_1;

	if (cases == 0 || controls == 0 || allele0 == 0 || allele1 == 0) {
		chisq = 0;
		pvalue = 1;
		return;
	}
	chisq = 0;
	double nij = cases * allele0 / total;
	chisq += pow(aff_0 - nij, 2) / nij;
	nij = cases * allele1 / total;
	chisq += pow(aff_1 - nij, 2) / nij;
	nij = controls * allele0 / total;
	chisq += pow(unaff_0 - nij, 2) / nij;
	nij = controls * allele1 / total;
	chisq += pow(unaff_1 - nij, 2) / nij;
	DBG_DO(DBG_STATOR, cout << " counts: "
		                    << aff_0 << " " << aff_1 << " " << unaff_0 << " " << unaff_1 << " ChiSq: " << chisq << endl);

	pvalue = GetRNG().pvalChiSq(chisq, 1);
}


void statAssociation::countAlleles(IndIterator & it, UINT loc,
                                   ULONG & aff_0, ULONG & aff_1, ULONG & unaff_0, ULONG & unaff_1)
{
	UINT ploidy = it->ploidy();

	aff_0 = 0;
	aff_1 = 0;
	unaff_0 = 0;
	unaff_1 = 0;
	for (; it.valid(); ++it) {
		for (UINT p = 0; p < ploidy; ++p) {
			if (it->affected()) {
				if (ToAllele(it->allele(loc, p)) == 0)
					++aff_0;
				else
					++aff_1;
			} else {
				if (ToAllele(it->allele(loc, p)) == 0)
					++unaff_0;
				else
					++unaff_1;
			}
		}
	}
}


bool statAssociation::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	pop.removeVar(Asso_ChiSq_String);
	pop.removeVar(Asso_ChiSq_P_String);
	pop.removeVar(Fit_String);

	vectorf chisq(pop.totNumLoci());
	vectorf pvalue(pop.totNumLoci());
	ULONG aff_0, aff_1, unaff_0, unaff_1;

	for (UINT i = 0; i < m_loci.size(); ++i) {
		UINT loc = m_loci[i];
		IndIterator it = pop.indIterator();
		countAlleles(it, loc, aff_0, aff_1, unaff_0, unaff_1);
		calcChiSq(aff_0, aff_1, unaff_0, unaff_1, chisq[loc], pvalue[loc]);
		DBG_DO(DBG_STATOR, cout << "Locus: " << loc << " ChiSq: " << chisq[loc]
			                    << " p-value: " << pvalue[loc] << endl);
	}
	pop.setDoubleVectorVar(Asso_ChiSq_String, chisq);
	pop.setDoubleVectorVar(Asso_ChiSq_P_String, pvalue);
	return true;
}


statNeutrality::statNeutrality(const vectorlu & loci, const subPopList & subPops) :
	m_loci(loci), m_subPops(subPops)
{
	sort(m_loci.begin(), m_loci.end());
	for (UINT i = 1; i < m_loci.size(); ++i) {
		DBG_FAILIF(m_loci[i - 1] == m_loci[i], ValueError, "Duplicated loci occurred in the loci list.");
	}
}


double statNeutrality::calcPi(IndIterator & it)
{
	UINT ploidy = it->ploidy();
	double diffCnt = 0;
	int numComparison = 0;
	GenoIterator gt, gt1;

	for (; it.valid(); ++it)
		for (IndIterator it1 = it; it1.valid(); ++it1)
			for (UINT p = 0; p < ploidy ; p++)
				for (UINT p1 = 0; p1 < ploidy; p1++) {
					// count numbers of the same person
					if (it == it1) {
						if (p < p1) {
							gt = it->genoBegin(p);
							gt1 = it->genoBegin(p1);
							for (vectorlu::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc)
								diffCnt += *(gt + *loc) != *(gt1 + *loc);
							++numComparison;
						}
					}
					// count numbers of different persons
					else {
						gt = it->genoBegin(p);
						gt1 = it1->genoBegin(p1);
						for (vectorlu::iterator loc = m_loci.begin(); loc != m_loci.end(); ++loc)
							diffCnt += *(gt + *loc) != *(gt1 + *loc);
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

	DBG_FAILIF(m_loci.size() > pop.totNumLoci(), IndexError,
		"Locus index out of range");

	pop.removeVar(Neutra_Pi_String);
	IndIterator it = pop.indIterator();
	pop.setDoubleVar(Neutra_Pi_String, calcPi(it));
	return true;
}


statFst::statFst(statAlleleFreq & alleleFreq, statHeteroFreq & heteroFreq,
	const vectorlu & Fst, const strDict & param)
	: m_alleleFreq(alleleFreq), m_heteroFreq(heteroFreq), m_atLoci(Fst),
	m_midValues(false),
	m_output_Fst(true),
	m_output_Fis(true),
	m_output_Fit(true),
	m_output_AvgFst(true),
	m_output_AvgFis(true),
	m_output_AvgFit(true)
{
	strDict::const_iterator it;
	strDict::const_iterator itEnd = param.end();

	if ((it = param.find("midValues")) != itEnd)
		m_midValues = it->second != 0.;
	// if any statistics is specified, other unspecified ones are not calculated
	if (param.find(Fst_String) != itEnd ||
	    param.find(Fis_String) != itEnd ||
	    param.find(Fit_String) != itEnd ||
	    param.find(AvgFst_String) != itEnd ||
	    param.find(AvgFis_String) != itEnd ||
	    param.find(AvgFit_String) != itEnd) {
		m_output_Fst = false;
		m_output_Fis = false;
		m_output_Fit = false;
		m_output_AvgFst = false;
		m_output_AvgFis = false;
		m_output_AvgFit = false;
		// if has key, and is True or 1
		if ((it = param.find(Fst_String)) != itEnd)
			m_output_Fst = it->second != 0.;
		if ((it = param.find(Fis_String)) != itEnd)
			m_output_Fis = it->second != 0.;
		if ((it = param.find(Fit_String)) != itEnd)
			m_output_Fit = it->second != 0.;
		if ((it = param.find(AvgFst_String)) != itEnd)
			m_output_AvgFst = it->second != 0.;
		if ((it = param.find(AvgFis_String)) != itEnd)
			m_output_AvgFis = it->second != 0.;
		if ((it = param.find(AvgFit_String)) != itEnd)
			m_output_AvgFit = it->second != 0.;
	}

	for (size_t i = 0; i < m_atLoci.size(); ++i) {
		// need to get allele frequency at this locus
		m_alleleFreq.addLocus(m_atLoci[i]);

		// need to get heterozygous proportion  at this locus
		m_heteroFreq.addLocus(m_atLoci[i]);
	}
}


bool statFst::apply(population & pop)
{
	if (m_atLoci.empty())
		return true;

	pop.removeVar(Fst_String);
	pop.removeVar(Fis_String);
	pop.removeVar(Fit_String);

	m_Fst.clear();
	m_Fit.clear();
	m_Fis.clear();

	// dicitonary to save values.
	UINT numSP = pop.numSubPop();
	ULONG popSize = pop.popSize();

	// do not save these values now
	double aa = 0., bb = 0., cc = 0.;

	// vector to store p[i]
	vectorf p_i = vectorf(numSP);

	// calculate Fst for each locus
	for (size_t st = 0; st < m_atLoci.size(); ++st) {
		int loc = m_atLoci[st];

		DBG_ASSERT(static_cast<size_t>(loc) < pop.totNumLoci(), IndexError,
			"Index out of range of 0 ~ " + toStr(pop.totNumLoci() - 1));

		// get all available alleles
		vectori alleles = m_alleleFreq.alleles(pop, loc);

		DBG_DO(DBG_STATOR, cout << "Using alleles " << alleles << endl);

		// n_bar
		double r = numSP;
		double n = popSize;
		double n_bar = n / r;
		vectorlu n_i = pop.subPopSizes();

		// n_c
		double n_c = n;
		for (int i = 0; i < r; ++i)
			n_c -= n_i[i] * n_i[i] / n;
		n_c /= (r - 1);

		double a = 0.0, b = 0.0, c = 0.0;

		for (vectori::iterator ale = alleles.begin(); ale != alleles.end(); ++ale) {
			// p_i
			for (int sp = 0; sp < r; ++sp)
				p_i[sp] = m_alleleFreq.alleleFreq(pop, *ale, loc, sp);

			// p_bar
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
				h_bar += m_heteroFreq.heteroFreq(*ale, loc, sp) * n_i[sp];
			h_bar /= n;

			// a, b, c
			a += n_bar / n_c * (s_2 - (p_bar * (1 - p_bar) - (r - 1.) / r * s_2 - h_bar / 4.) / (n_bar - 1.) );
			b += n_bar / (n_bar - 1) * (p_bar * (1 - p_bar) - (r - 1) / r * s_2 - (2 * n_bar - 1) / (4. * n_bar) * h_bar);
			c += h_bar / 2.;

			DBG_DO(DBG_STATOR, cout << "allele " << *ale << "\tn_c: " << n_c
				                    << "\tp_i: " << p_i << "\tp_bar: " << p_bar << "\ts^2: " << s_2 << "\th_bar:"
				                    << h_bar << "\ta: " << a << "\tb: " << b << "\tc: " << c << endl);
		}                                                                                 // each allele

		DBG_DO(DBG_STATOR, cout << "Fst= " << a / (a + b + c) << endl);

		if (static_cast<size_t>(loc) >= m_Fst.size()) {
			m_Fst.resize(loc + 1, 0.);
			m_Fit.resize(loc + 1, 0.);
			m_Fis.resize(loc + 1, 0.);
		}
		m_Fst[loc] = fcmp_eq(a + b + c, 0.) ? 0. : (a / (a + b + c));
		m_Fit[loc] = fcmp_eq(a + b + c, 0.) ? 1. : (1 - c / (a + b + c));
		m_Fis[loc] = fcmp_eq(b + c, 0.) ? 1. : (1 - c / (b + c));

		aa += a;
		bb += b;
		cc += c;
	}
	m_avgFst = fcmp_eq(aa + bb + cc, 0.) ? 0 : (aa / (aa + bb + cc));
	m_avgFit = fcmp_eq(aa + bb + cc, 0.) ? 1. : (1 - cc / (aa + bb + cc));
	m_avgFis = fcmp_eq(aa + bb + cc, 0) ? 1. : (1 - cc / (bb + cc));

	// post results
	if (m_output_Fst)
		pop.setDoubleVectorVar(Fst_String, m_Fst);
	if (m_output_Fit)
		pop.setDoubleVectorVar(Fit_String, m_Fit);
	if (m_output_Fis)
		pop.setDoubleVectorVar(Fis_String, m_Fis);
	if (m_output_AvgFst)
		pop.setDoubleVar(AvgFst_String, m_avgFst);
	if (m_output_AvgFit)
		pop.setDoubleVar(AvgFit_String, m_avgFit);
	if (m_output_AvgFis)
		pop.setDoubleVar(AvgFis_String, m_avgFis);
	return true;
}


statHWE::statHWE(statGenoFreq & genoFreq, const vectorlu & loci)
	: m_genoFreq(genoFreq), m_loci(loci)
{
}


double statHWE::calcHWE(const vectorlu & cnt)
{
	// Calculates exact two-sided hardy-weinberg p-value. Parameters
	// are number of genotypes, number of rare alleles observed and
	// number of heterozygotes observed.
	//
	// (c) 2003 Jan Wigginton, Goncalo Abecasis
	int obsAA = cnt[2];                                             // in this algorithm, AA is rare.
	int obsAB = cnt[1];
	int obsBB = cnt[0];

	int diplotypes = obsAA + obsAB + obsBB;
	int rare = (obsAA * 2) + obsAB;
	int hets = obsAB;


	//make sure "rare" allele is really the rare allele
	if (rare > diplotypes)
		rare = 2 * diplotypes - rare;

	//make sure numbers aren't screwy
	if (hets > rare)
		throw ValueError("HW test: " + toStr(hets) + "heterozygotes but only " + toStr(rare) + "rare alleles.");

	vectorf tailProbs(rare + 1, 0.0);

	//start at midpoint
	//all the casting is to make sure we don't overflow ints if there are 10's of 1000's of inds
	int mid = (int)((double)rare * (double)(2 * diplotypes - rare) / (double)(2 * diplotypes));

	//check to ensure that midpoint and rare alleles have same parity
	if (((rare & 1) ^ (mid & 1)) != 0) {
		mid++;
	}
	int het = mid;
	int hom_r = (rare - mid) / 2;
	int hom_c = diplotypes - het - hom_r;

	//Calculate probability for each possible observed heterozygote
	//count up to a scaling constant, to avoid underflow and overflow
	tailProbs[mid] = 1.0;
	double sum = tailProbs[mid];
	for (het = mid; het > 1; het -= 2) {
		tailProbs[het - 2] = (tailProbs[het] * het * (het - 1.0)) / (4.0 * (hom_r + 1.0) * (hom_c + 1.0));
		sum += tailProbs[het - 2];
		//2 fewer hets for next iteration -> add one rare and one common homozygote
		hom_r++;
		hom_c++;
	}

	het = mid;
	hom_r = (rare - mid) / 2;
	hom_c = diplotypes - het - hom_r;
	for (het = mid; het <= rare - 2; het += 2) {
		tailProbs[het + 2] = (tailProbs[het] * 4.0 * hom_r * hom_c) / ((het + 2.0) * (het + 1.0));
		sum += tailProbs[het + 2];
		//2 more hets for next iteration -> subtract one rare and one common homozygote
		hom_r--;
		hom_c--;
	}

	for (size_t z = 0; z < tailProbs.size(); z++)
		tailProbs[z] /= sum;

	double top = tailProbs[hets];
	for (int i = hets + 1; i <= rare; i++)
		top += tailProbs[i];

	double otherSide = tailProbs[hets];
	for (int i = hets - 1; i >= 0; i--)
		otherSide += tailProbs[i];

	if (top > 0.5 && otherSide > 0.5) {
		return 1.0;
	} else {
		if (top < otherSide)
			return top * 2;
		else
			return otherSide * 2;
	}

	return 1.0;
}


bool statHWE::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	UINT numSP = pop.numSubPop();
	UINT numLoci = m_loci.size();

	pop.removeVar(HWE_String);
	for (UINT sp = 0; sp < numSP; ++sp)
		pop.removeVar(subPopVar_String(sp, HWE_String));

	// calculate HWE for the whole population
	vectorf HWE(numLoci);
	for (size_t i = 0; i < numLoci; ++i) {
		UINT loc = m_loci[i];

		DBG_FAILIF(loc >= pop.totNumLoci(), IndexError,
			"Locus index out of range.");

		if (HWE.size() <= loc)
			HWE.resize(loc + 1, 0);

		HWE[loc] = calcHWE(m_genoFreq.countGenotype(pop, loc, 0));
	}
	pop.setDoubleVectorVar(HWE_String, HWE);

	if (numSP == 1) {
		string varname = subPopVar_String(0, HWE_String);
		pop.setDoubleVectorVar(varname, HWE);
		return true;
	}

	for (UINT sp = 0; sp < numSP; ++sp) {
		vectorf HWE(numLoci);
		for (size_t i = 0; i < numLoci; ++i) {
			UINT loc = m_loci[i];
			if (HWE.size() <= loc)
				HWE.resize(loc + 1);
			HWE[loc] = calcHWE(m_genoFreq.countGenotype(pop, loc, sp, 0));
		}
		string varname = subPopVar_String(sp, HWE_String);
		pop.setDoubleVectorVar(varname, HWE);
	}

	return true;
}


}

