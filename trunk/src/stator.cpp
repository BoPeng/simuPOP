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
	m_heteroFreq(heteroFreq.elems(), homoFreq.elems(), subPops, vars),
	m_genoFreq(genoFreq.elems(), subPops, vars),
	m_haploFreq(haploFreq, subPops, vars),
	m_info(sumOfInfo.elems(), meanOfInfo.elems(), varOfInfo.elems(), maxOfInfo.elems(), minOfInfo.elems(), subPops, vars),
	m_LD(LD, subPops, vars),
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
	m_genoFreq(rhs.m_genoFreq),
	m_haploFreq(rhs.m_haploFreq),
	m_info(rhs.m_info),
	m_LD(rhs.m_LD),
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
	       m_genoFreq.apply(pop) &&
	       m_haploFreq.apply(pop) &&
	       m_info.apply(pop) &&
	       m_LD.apply(pop) &&
	       m_association.apply(pop) &&
	       m_neutrality.apply(pop) &&
	       m_Fst.apply(pop) &&
	       m_HWE.apply(pop);
}


statPopSize::statPopSize(bool popSize, const subPopList & subPops, const stringList & vars)
		: m_isActive(popSize), m_subPops(subPops), m_vars()
	{
		const char * allowedVars[] = { popSize_String,    popSize_sp_String,
			                           subPopSize_String, "" };
		const char * defaultVars[] = { popSize_String, subPopSize_String, "" };

		m_vars.obtainFrom(vars, allowedVars, defaultVars);
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
	subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	for (; it != itEnd; ++it) {
		ULONG spPopSize = pop.subPopSize(*it);
		popSize += spPopSize;
		spSize.push_back(spPopSize);
		if (m_vars.contains(popSize_sp_String))
			pop.setIntVar(subPopVar_String(*it, popSize_String), spPopSize);
	}
	// NOTE: popSize does not have to be the total population size
	if (m_vars.contains(popSize_String))
		pop.setIntVar(popSize_String, popSize);
	// subPopSize = ...
	if (m_vars.contains(subPopSize_String))
		pop.setIntVectorVar(subPopSize_String, spSize);
	return true;
}


statNumOfMale::statNumOfMale(bool numOfMale, const subPopList & subPops, const stringList & vars)
		: m_isActive(numOfMale), m_subPops(subPops), m_vars()
	{
		const char * allowedVars[] = {
			numOfMale_String,      propOfMale_String,
			numOfFemale_String,    propOfFemale_String,
			numOfMale_sp_String,   propOfMale_sp_String,
			numOfFemale_sp_String, propOfFemale_sp_String,""
		};
		const char * defaultVars[] = { numOfMale_String, numOfFemale_String, "" };

		m_vars.obtainFrom(vars, allowedVars, defaultVars);
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

		if (m_vars.contains(numOfMale_sp_String))
			pop.setIntVar(subPopVar_String(*sp, numOfMale_String), maleCnt);
		if (m_vars.contains(propOfMale_sp_String))
			pop.setDoubleVar(subPopVar_String(*sp, propOfMale_String),
				totalCnt == 0 ? 0 : static_cast<double>(maleCnt) / totalCnt);
		if (m_vars.contains(numOfFemale_sp_String))
			pop.setIntVar(subPopVar_String(*sp, numOfFemale_String), femaleCnt);
		if (m_vars.contains(propOfFemale_sp_String))
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


statNumOfAffected::statNumOfAffected(bool numOfAffected, const subPopList & subPops, const stringList & vars)
		: m_isActive(numOfAffected), m_subPops(subPops), m_vars()
	{
		const char * allowedVars[] = {
			numOfAffected_String,      propOfAffected_String,
			numOfUnaffected_String,    propOfUnaffected_String,
			numOfAffected_sp_String,   propOfAffected_sp_String,
			numOfUnaffected_sp_String, propOfUnaffected_sp_String,""
		};
		const char * defaultVars[] = { numOfAffected_String, numOfUnaffected_String, "" };

		m_vars.obtainFrom(vars, allowedVars, defaultVars);
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

		if (m_vars.contains(numOfAffected_sp_String))
			pop.setIntVar(subPopVar_String(*sp, numOfAffected_String), affectedCnt);
		if (m_vars.contains(propOfAffected_sp_String))
			pop.setDoubleVar(subPopVar_String(*sp, propOfAffected_String),
				totalCnt == 0 ? 0 : static_cast<double>(affectedCnt) / totalCnt);
		if (m_vars.contains(numOfUnaffected_sp_String))
			pop.setIntVar(subPopVar_String(*sp, numOfUnaffected_String), unaffectedCnt);
		if (m_vars.contains(propOfUnaffected_sp_String))
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
		pop.setDoubleVar(propOfAffected_String,
			allTotalCnt == 0 ? 0. : static_cast<double>(allAffectedCnt) / allTotalCnt);
	if (m_vars.contains(numOfUnaffected_String))
		pop.setIntVar(numOfUnaffected_String, allUnaffectedCnt);
	if (m_vars.contains(propOfUnaffected_String))
		pop.setDoubleVar(propOfUnaffected_String,
			allTotalCnt == 0 ? 0 : static_cast<double>(allUnaffectedCnt) / allTotalCnt);

	return true;
}


statAlleleFreq::statAlleleFreq(const vectorlu & loci, const subPopList & subPops, const stringList & vars)
		: m_loci(loci), m_subPops(subPops), m_vars()
	{
		const char * allowedVars[] = {
			AlleleNum_String,    AlleleFreq_String,
			AlleleNum_sp_String, AlleleFreq_sp_String,""
		};
		const char * defaultVars[] = { AlleleFreq_String, AlleleNum_String, "" };

		m_vars.obtainFrom(vars, allowedVars, defaultVars);
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
		if (m_vars.contains(AlleleNum_sp_String))
			pop.removeVar(subPopVar_String(*it, AlleleNum_String));
		if (m_vars.contains(AlleleFreq_sp_String))
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
			if (m_vars.contains(AlleleNum_sp_String))
				pop.setIntVectorVar(subPopVar_String(*it, AlleleNum_String) + "{" + toStr(loc) + "}", alleles);
			if (m_vars.contains(AlleleFreq_sp_String)) {
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


statHeteroFreq::statHeteroFreq(const vectorlu & heteroFreq, const vectorlu & homoFreq,
	const subPopList & subPops, const stringList & vars)
	: m_loci(heteroFreq), m_subPops(subPops), m_vars()
{
	// add homofreq to m_loci
	for (size_t i = 0; i < homoFreq.size(); ++i)
		if (find(m_loci.begin(), m_loci.end(), homoFreq[i]) == m_loci.end())
			m_loci.push_back(homoFreq[i]);
	//
	const char * allowedVars[] = {
		HeteroNum_String,    HeteroFreq_String,
		HeteroNum_sp_String, HeteroFreq_sp_String,
		HomoNum_String,      HomoFreq_String,
		HomoNum_sp_String,   HomoFreq_sp_String,
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



bool statHeteroFreq::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	DBG_FAILIF(pop.ploidy() != 2, ValueError,
		"Heterozygote frequency can only be calculated for diploid populations.");

	DBG_DO(DBG_STATOR, cout << "Calculated heterozygote frequency for loci " << m_loci << endl);

	// count for all specified subpopulations
	intDict allHeteroCnt;
	intDict allHomoCnt;

	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	for (; it != itEnd; ++it) {
		if (it->isVirtual())
			pop.activateVirtualSubPop(*it);

		intDict heteroCnt;
		intDict homoCnt;

		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];

#ifndef OPTIMIZED
			int chromType = pop.chromType(pop.chromLocusPair(loc).first);
			DBG_FAILIF(chromType == ChromosomeX || chromType == ChromosomeY,
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
		if (it->isVirtual())
			pop.deactivateVirtualSubPop(it->subPop());
		// output subpopulation variable?
		if (m_vars.contains(HeteroNum_sp_String))
			pop.setIntDictVar(subPopVar_String(*it, HeteroNum_String), heteroCnt);
		if (m_vars.contains(HomoNum_sp_String))
			pop.setIntDictVar(subPopVar_String(*it, HomoNum_String), homoCnt);
		if (m_vars.contains(HeteroFreq_sp_String)) {
			intDict freq;
			for (size_t idx = 0; idx < m_loci.size(); ++idx) {
				UINT loc = m_loci[idx];
				double all = heteroCnt[loc] + homoCnt[loc];
				freq[loc] = all == 0. ? 0 : heteroCnt[loc] / all;
			}
			pop.setIntDictVar(subPopVar_String(*it, HeteroFreq_String), freq);
		}
		if (m_vars.contains(HomoFreq_sp_String)) {
			intDict freq;
			for (size_t idx = 0; idx < m_loci.size(); ++idx) {
				UINT loc = m_loci[idx];
				double all = heteroCnt[loc] + homoCnt[loc];
				freq[loc] = all == 0. ? 0 : homoCnt[loc] / all;
			}
			pop.setIntDictVar(subPopVar_String(*it, HomoFreq_String), freq);
		}
	}
	// for whole population.
	if (m_vars.contains(HeteroNum_String))
		pop.setIntDictVar(HeteroNum_String, allHeteroCnt);
	if (m_vars.contains(HomoNum_String))
		pop.setIntDictVar(HomoNum_String, allHomoCnt);
	if (m_vars.contains(HeteroFreq_String)) {
		intDict freq;
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];
			double all = allHeteroCnt[loc] + allHomoCnt[loc];
			freq[loc] = all == 0. ? 0 : allHeteroCnt[loc] / all;
		}
		pop.setIntDictVar(HeteroFreq_String, freq);
	}
	if (m_vars.contains(HomoFreq_String)) {
		intDict freq;
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];
			double all = allHeteroCnt[loc] + allHomoCnt[loc];
			freq[loc] = all == 0. ? 0 : allHomoCnt[loc] / all;
		}
		pop.setIntDictVar(HomoFreq_String, freq);
	}

	return true;
}


statGenoFreq::statGenoFreq(const vectorlu & genoFreq, const subPopList & subPops, const stringList & vars)
	: m_loci(genoFreq), m_subPops(subPops), m_vars()
{
	const char * allowedVars[] = {
		GenotypeNum_String,    GenotypeFreq_String,
		GenotypeNum_sp_String, GenotypeFreq_sp_String,""
	};
	const char * defaultVars[] = { GenotypeFreq_String, GenotypeNum_String, "" };

	m_vars.obtainFrom(vars, allowedVars, defaultVars);
}


bool statGenoFreq::apply(population & pop)
{
	if (m_loci.empty())
		return true;

	DBG_DO(DBG_STATOR, cout << "Calculated genotype frequency for loci " << m_loci << endl);

	// count for all specified subpopulations
	vector<tupleDict> genotypeCnt(m_loci.size());
	vectorlu allGenotypeCnt(m_loci.size(), 0);
	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	UINT ply = pop.ploidy();
	for (; it != itEnd; ++it) {
		if (m_vars.contains(GenotypeNum_sp_String))
			pop.removeVar(subPopVar_String(*it, GenotypeNum_String));
		if (m_vars.contains(GenotypeFreq_sp_String))
			pop.removeVar(subPopVar_String(*it, GenotypeFreq_String));

		if (it->isVirtual())
			pop.activateVirtualSubPop(*it);

		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];

			tupleDict genotypes;
			size_t allGenotypes = 0;

			// go through all alleles
			IndIterator ind = pop.indIterator(it->subPop());
			for (; ind.valid(); ++ind) {
				vectori genotype(ply);
				for (size_t p = 0; p < ply; ++p)
					genotype[p] = ind->allele(loc, p);
				genotypes[genotype]++;
				allGenotypes++;
			}
			// total allele count
			tupleDict::iterator dct = genotypes.begin();
			tupleDict::iterator dctEnd = genotypes.end();
			for (; dct != dctEnd; ++dct)
				genotypeCnt[idx][dct->first] += dct->second;
			allGenotypeCnt[idx] += allGenotypes;
			// output variable.
			if (m_vars.contains(GenotypeNum_sp_String))
				pop.setTupleDictVar(subPopVar_String(*it, GenotypeNum_String) + "{" + toStr(loc) + "}", genotypes);
			// note that genotyeps is changed in place.
			if (m_vars.contains(GenotypeFreq_sp_String)) {
				if (allGenotypes != 0) {
					tupleDict::iterator dct = genotypes.begin();
					tupleDict::iterator dctEnd = genotypes.end();
					for (; dct != dctEnd; ++dct)
						dct->second /= allGenotypes;
				}
				pop.setTupleDictVar(subPopVar_String(*it, GenotypeFreq_String) + "{" + toStr(loc) + "}", genotypes);
			}
		}
		if (it->isVirtual())
			pop.deactivateVirtualSubPop(it->subPop());
	}

	if (m_vars.contains(GenotypeNum_String)) {
		pop.removeVar(GenotypeNum_String);
		for (size_t idx = 0; idx < m_loci.size(); ++idx)
			pop.setTupleDictVar(string(GenotypeNum_String) + "{" + toStr(m_loci[idx]) + "}",
				genotypeCnt[idx]);
	}
	// note that genotyeCnt[idx] is changed in place.
	if (m_vars.contains(GenotypeFreq_String)) {
		pop.removeVar(GenotypeFreq_String);
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			UINT loc = m_loci[idx];
			if (allGenotypeCnt[idx] != 0) {
				tupleDict::iterator dct = genotypeCnt[idx].begin();
				tupleDict::iterator dctEnd = genotypeCnt[idx].end();
				for (; dct != dctEnd; ++dct)
					dct->second /= allGenotypeCnt[idx];
			}
			pop.setTupleDictVar(string(GenotypeFreq_String) + "{" + toStr(loc) + "}",
				genotypeCnt[idx]);
		}
	}

	return true;
}


statHaploFreq::statHaploFreq(const intMatrix & haploFreq, const subPopList & subPops, const stringList & vars)
		: m_loci(haploFreq), m_subPops(subPops), m_vars()
	{
		const char * allowedVars[] = {
			HaplotypeNum_String,	HaplotypeFreq_String,
			HaplotypeNum_sp_String, HaplotypeFreq_sp_String,""
		};
		const char * defaultVars[] = { HaplotypeFreq_String, HaplotypeNum_String, "" };

		m_vars.obtainFrom(vars, allowedVars, defaultVars);
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

	DBG_DO(DBG_STATOR, cout << "Calculated haplotype frequency for loci " << m_loci << endl);

	// count for all specified subpopulations
	vector<tupleDict> haplotypeCnt(m_loci.size());
	vectorlu allHaplotypeCnt(m_loci.size());
	// selected (virtual) subpopulatons.
	subPopList subPops = m_subPops;
	subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	UINT ply = pop.ploidy();
	for (; it != itEnd; ++it) {
		if (m_vars.contains(HaplotypeNum_sp_String))
			pop.removeVar(subPopVar_String(*it, HaplotypeNum_String));
		if (m_vars.contains(HaplotypeFreq_sp_String))
			pop.removeVar(subPopVar_String(*it, HaplotypeFreq_String));

		if (it->isVirtual())
			pop.activateVirtualSubPop(*it);

		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			const vectori & loci = m_loci[idx];
			UINT nLoci = loci.size();
			if (nLoci == 0)
				continue;

			UINT chromType = pop.chromType(pop.chromLocusPair(loci[0]).first);
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
					if (p == 1 && ind->sex() == Male && pop.isHaplodiploid())
						continue;
					if (chromType == ChromosomeY && ind->sex() == Female)
						continue;
					if (((chromType == ChromosomeX && p == 1) || 
						(chromType == ChromosomeY && p == 0)) && ind->sex() == Male)
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
				pop.setTupleDictVar(subPopVar_String(*it, HaplotypeNum_String) + "{"
					+ key + "}", haplotypes);
			// note that genotyeps is changed in place.
			if (m_vars.contains(HaplotypeFreq_sp_String)) {
				if (allHaplotypes != 0) {
					tupleDict::iterator dct = haplotypes.begin();
					tupleDict::iterator dctEnd = haplotypes.end();
					for (; dct != dctEnd; ++dct)
						dct->second /= allHaplotypes;
				}
				pop.setTupleDictVar(subPopVar_String(*it, HaplotypeFreq_String) + "{"
					+ key + "}", haplotypes);
			}
		}
		if (it->isVirtual())
			pop.deactivateVirtualSubPop(it->subPop());
	}

	if (m_vars.contains(HaplotypeNum_String)) {
		pop.removeVar(HaplotypeNum_String);
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			string key = dictKey(m_loci[idx]);
			pop.setTupleDictVar(string(HaplotypeNum_String) + "{" + key + "}",
				haplotypeCnt[idx]);
		}
	}
	// note that genotyeCnt[idx] is changed in place.
	if (m_vars.contains(HaplotypeFreq_String)) {
		pop.removeVar(HaplotypeFreq_String);
		for (size_t idx = 0; idx < m_loci.size(); ++idx) {
			string key = dictKey(m_loci[idx]);
			if (allHaplotypeCnt[idx] != 0) {
				tupleDict::iterator dct = haplotypeCnt[idx].begin();
				tupleDict::iterator dctEnd = haplotypeCnt[idx].end();
				for (; dct != dctEnd; ++dct)
					dct->second /= allHaplotypeCnt[idx];
			}
			pop.setTupleDictVar(string(HaplotypeFreq_String) + "{" + key + "}",
				haplotypeCnt[idx]);
		}
	}

	return true;
}


statInfo::statInfo(const vectorstr & sumOfInfo, const vectorstr & meanOfInfo,
		const vectorstr & varOfInfo, const vectorstr & maxOfInfo,
		const vectorstr & minOfInfo,
		const subPopList & subPops, const stringList & vars)
		: m_sumOfInfo(sumOfInfo), m_meanOfInfo(meanOfInfo), m_varOfInfo(varOfInfo),
		m_maxOfInfo(maxOfInfo), m_minOfInfo(minOfInfo),
		m_subPops(subPops), m_vars()
	{
		const char * allowedVars[] = {
			SumOfInfo_String,    MeanOfInfo_String,    VarOfInfo_String,	MaxOfInfo_String,    MinOfInfo_String,
			SumOfInfo_sp_String, MeanOfInfo_sp_String, VarOfInfo_sp_String, MaxOfInfo_sp_String, MinOfInfo_sp_String,
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


bool statInfo::apply(population & pop)
{
	if (m_sumOfInfo.empty() && m_meanOfInfo.empty() && m_maxOfInfo.empty() && m_minOfInfo.empty())
		return true;

	// field indexes
	UINT numSumFld = m_sumOfInfo.size();
	UINT numMeanFld = m_meanOfInfo.size();
	UINT numVarFld = m_meanOfInfo.size();
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
	vectorlu allMeanNumVal(numMeanFld);
	vectorf allVarSumVal(numVarFld);
	vectorf allVarSum2Val(numVarFld);
	vectorlu allVarNumVal(numVarFld);
	vectorf allMaxVal(0);
	vectorf allMinVal(0);
	// for each subpopulation.
	subPopList subPops = m_subPops;
	subPops.useSubPopsFrom(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for (; sp != spEnd; ++sp) {
		vectorf sumVal(numSumFld, 0.);
		vectorf meanSumVal(numMeanFld, 0.);
		vectorlu meanNumVal(numMeanFld, 0);
		vectorf varSumVal(numVarFld, 0.);
		vectorf varSum2Val(numVarFld, 0.);
		vectorlu varNumVal(numVarFld, 0);
		vectorf maxVal(0);
		vectorf minVal(0);

		if (sp->isVirtual())
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

		if (sp->isVirtual())
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
			pop.setStrDictVar(subPopVar_String(*sp, SumOfInfo_String), dct);
		}
		if (m_vars.contains(MeanOfInfo_sp_String)) {
			strDict dct;
			for (size_t i = 0; i < numMeanFld; ++i)
				dct[m_meanOfInfo[i]] = meanNumVal[i] == 0 ? 0 : meanSumVal[i] / meanNumVal[i];
			pop.setStrDictVar(subPopVar_String(*sp, MeanOfInfo_String), dct);
		}
		if (m_vars.contains(VarOfInfo_sp_String)) {
			strDict dct;
			for (size_t i = 0; i < numVarFld; ++i)
				dct[m_varOfInfo[i]] = varNumVal[i] <= 1 ? 0 :
				                      (varSum2Val[i] - varSumVal[i] * varSumVal[i] / varNumVal[i]) / (varNumVal[i] - 1);
			pop.setStrDictVar(subPopVar_String(*sp, VarOfInfo_String), dct);
		}
		if (m_vars.contains(MaxOfInfo_sp_String)) {
			strDict dct;
			for (size_t i = 0; i < numMaxFld; ++i)
				dct[m_maxOfInfo[i]] = maxVal[i];
			pop.setStrDictVar(subPopVar_String(*sp, MaxOfInfo_String), dct);
		}
		if (m_vars.contains(MinOfInfo_sp_String)) {
			strDict dct;
			for (size_t i = 0; i < numMinFld; ++i)
				dct[m_minOfInfo[i]] = minVal[i];
			pop.setStrDictVar(subPopVar_String(*sp, MinOfInfo_String), dct);
		}
	}
	if (m_vars.contains(SumOfInfo_String)) {
		strDict dct;
		for (size_t i = 0; i < m_sumOfInfo.size(); ++i)
			dct[m_sumOfInfo[i]] = allSumVal[i];
		pop.setStrDictVar(SumOfInfo_String, dct);
	}
	if (m_vars.contains(MeanOfInfo_String)) {
		strDict dct;
		for (size_t i = 0; i < numMeanFld; ++i)
			dct[m_meanOfInfo[i]] = allMeanNumVal[i] == 0 ? 0 : allMeanSumVal[i] / allMeanNumVal[i];
		pop.setStrDictVar(MeanOfInfo_String, dct);
	}
	if (m_vars.contains(VarOfInfo_String)) {
		strDict dct;
		for (size_t i = 0; i < numVarFld; ++i)
			dct[m_varOfInfo[i]] = allVarNumVal[i] <= 1 ? 0 :
			                      (allVarSum2Val[i] - allVarSumVal[i] * allVarSumVal[i] / allVarNumVal[i]) / (allVarNumVal[i] - 1);
		pop.setStrDictVar(VarOfInfo_String, dct);
	}
	if (m_vars.contains(MaxOfInfo_String)) {
		strDict dct;
		for (size_t i = 0; i < numMaxFld; ++i)
			dct[m_maxOfInfo[i]] = allMaxVal[i];
		pop.setStrDictVar(MaxOfInfo_String, dct);
	}
	if (m_vars.contains(MinOfInfo_String)) {
		strDict dct;
		for (size_t i = 0; i < numMinFld; ++i)
			dct[m_minOfInfo[i]] = allMinVal[i];
		pop.setStrDictVar(MinOfInfo_String, dct);
	}
	return true;
}

	
statLD::statLD(const intMatrix & LD,  const subPopList & subPops,
		const stringList & vars)
	: m_LD(LD), m_subPops(subPops), m_vars()
{
	const char * allowedVars[] = {
		LD_String, LD_prime_String, R2_String,
		ChiSq_String, ChiSq_p_String, CramerV_String,
		LD_sp_String, LD_prime_sp_String, R2_sp_String,
		ChiSq_sp_String, ChiSq_p_sp_String, CramerV_sp_String,
		""
	};
	const char * defaultVars[] = { LD_String, LD_prime_String, R2_String, "" };
	m_vars.obtainFrom(vars, allowedVars, defaultVars);

	for (size_t i = 0; i < m_LD.size(); ++i) {
		DBG_FAILIF(m_LD[i].size() != 2, ValueError,
			"Parameter LD should be a list of loci pairs.");
	}
}


void statLD::calculateLD(const vectoru & lociMap,
		const ALLELECNTLIST & alleleCnt, const HAPLOCNTLIST & haploCnt,
		vectorf & LD, vectorf & D_prime, vectorf & R2, vectorf & ChiSq, vectorf & ChiSq_p,
		vectorf & CramerV)
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
		for (; cnt != cntEnd; ++cnt) {
			alleles1.push_back(cnt->first);
			freq1.push_back(cnt->second);
		}
		cnt = alleleCnt2.begin();
		cntEnd = alleleCnt2.end();
		for (; cnt != cntEnd; ++cnt) {
			alleles2.push_back(cnt->first);
			freq2.push_back(cnt->second);
		}
		UINT nAllele1 = alleles1.size();
		UINT nAllele2 = alleles2.size();
		// get haplotype count
		const HAPLOCNT & haplos = haploCnt[idx];
		// get total haplotype count (used to calculate haplotype frequency)
		double allHaplo = 0;
		HAPLOCNT::const_iterator hCnt = haplos.begin();
		HAPLOCNT::const_iterator hCntEnd = haplos.end();
		for (; hCnt != hCntEnd; ++hCnt)
			allHaplo += hCnt->second;
		//
		// calculate LD
		if (!LD.empty()) {
			for (size_t i = 0; i < nAllele1; ++i) {
				for (size_t j = 0; j < nAllele2; ++i) {
					UINT A = alleles1[i];
					UINT B = alleles2[j];
					double P_AB = haplos.find(HAPLOCNT::key_type(A, B))->second / allHaplo;
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
					DBG_DO(DBG_STATOR, cout << "P_AB: " << P_AB
								<< " P_A: " << P_A << " P_B: " << P_B << " D_max: " << D_max <<
								" LD: " << D << " LD': " << D_prime << " r2: " << r2 << endl);

					if (nAllele1 <= 2 && nAllele2 <= 2) {
						// for the monomorphic or diallelic case, there is no need to do an average.
						LD[idx] = fabs(D);
						D_prime[idx] = Dp;
						R2[idx] = r2;
						break;
					} else {
						// for the monomorphic or diallelic case, there is no need to do an average.
						LD[idx] += P_A * P_B * fabs(D);
						D_prime[idx] += P_A * P_B * Dp;
						R2[idx] += P_A * P_B * r2;
					}
				}
				if (nAllele1 <= 2 && nAllele2 <= 2)
					break;
			}
		}
		// if ChiSq is empty, do not calculate any association stuff.
		if (!ChiSq.empty()) {
			// calculate association
			vector<vectorf> cont_table(nAllele1 + 1);
			for (size_t i = 0; i <= nAllele1; ++i)
				cont_table[i].resize(nAllele2 + 1);
			// initialize last line/column
			for (size_t i = 0; i < nAllele1; ++i)
				cont_table[i][nAllele2] = 0;
			for (size_t j = 0; j <= nAllele2; ++j)
				cont_table[nAllele1][j] = 0;
			// get P_ij
			for (size_t i = 0; i < nAllele1; ++i) {
				for (size_t j = 0; j < nAllele2; ++j) {
					cont_table[i][j] = haplos.find(HAPLOCNT::key_type(alleles1[i], alleles2[j]))->second / allHaplo;
					cont_table[i][nAllele2] += cont_table[i][j];
					cont_table[nAllele1][j] += cont_table[i][j];
					cont_table[nAllele2][nAllele1] += cont_table[i][j];
				}
			}
			DBG_ASSERT(fcmp_eq(cont_table[nAllele1][nAllele2], 1.), ValueError,
				"Sum of haplotype frequencies is not 1. Association will not be computed.");
			// calculate statistics
			//ChiSq
			for (size_t i = 0; i < nAllele1; ++i)
				for (size_t j = 0; j < nAllele2; ++j)
					ChiSq[idx] += pow((allHaplo * cont_table[i][j] - allHaplo * cont_table[i][nAllele2] * cont_table[nAllele1][j]), 2)
					         / (allHaplo * cont_table[i][nAllele2] * cont_table[nAllele1][j]);
			ChiSq_p[idx] = GetRNG().pvalChiSq(ChiSq[idx], (nAllele1 - 1) * (nAllele2 - 1));
			CramerV[idx] = sqrt(ChiSq[idx] / (allHaplo * std::min(nAllele1 - 1, nAllele2 - 1)));
		}
	}
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
	
	pop.removeVar(name);
	map<UINT, intDict>::const_iterator it = res.begin();
	map<UINT, intDict>::const_iterator itEnd = res.end();
	for (; it != itEnd; ++it)
		pop.setIntDictVar(name + "{" + toStr(it->first) + "}", it->second);
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
			if (find(loci.begin(), loci.end(), m_LD[i][j]) == loci.end()) {
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
	subPops.useSubPopsFrom(pop);
	subPopList::const_iterator it = subPops.begin();
	subPopList::const_iterator itEnd = subPops.end();
	UINT ply = pop.ploidy();
	for (; it != itEnd; ++it) {
		const char * spVars[] = {
			LD_sp_String, LD_prime_sp_String, R2_sp_String,
			ChiSq_sp_String, ChiSq_p_sp_String, CramerV_sp_String,
			""
		};
		for (size_t i = 0; spVars[i][0]; ++i) {
			if (m_vars.contains(spVars[i]))
				pop.removeVar(subPopVar_String(*it, spVars[i]));
		}

		if (it->isVirtual())
			pop.activateVirtualSubPop(*it);

		ALLELECNTLIST alleleCnt(loci.size());
		HAPLOCNTLIST haploCnt(m_LD.size());
	
		// count allele and genotype
		IndIterator ind = pop.indIterator(it->subPop());
		for (; ind.valid(); ++ind) {
			for (size_t p = 0; p < ply; ++p) {
				if (p == 1 && ind->sex() == Male && pop.isHaplodiploid())
					continue;
				GenoIterator geno = ind->genoBegin(p);
				// allele frequency
				for (size_t idx = 0; idx < nLoci; ++idx) {
					if (chromTypes[idx] == ChromosomeY && ind->sex() == Female)
						continue;
					if (((chromTypes[idx] == ChromosomeX && p == 1) || 
						(chromTypes[idx] == ChromosomeY && p == 0)) && ind->sex() == Male)
						continue;
					alleleCnt[idx][*(geno + loci[idx])] ++;
				}
				// haplotype frequency
				for (size_t idx = 0; idx < nLD; ++idx) {
					UINT chromType = chromTypes[lociMap[m_LD[idx][0]]];
					if (chromType == ChromosomeY && ind->sex() == Female)
						continue;
					if (((chromType == ChromosomeX && p == 1) || 
						(chromType == ChromosomeY && p == 0)) && ind->sex() == Male)
						continue;
					haploCnt[idx][HAPLOCNT::key_type(*(geno + m_LD[idx][0]), *(geno + m_LD[idx][1]))] ++;
				}
			}
		}
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
	: m_alleleFreq(alleleFreq), m_heteroFreq(heteroFreq), m_loci(Fst),
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

	for (size_t i = 0; i < m_loci.size(); ++i) {
		// need to get allele frequency at this locus
		//m_alleleFreq.addLocus(m_loci[i]);

		// need to get heterozygous proportion  at this locus
		//m_heteroFreq.addLocus(m_loci[i]);
	}
}


bool statFst::apply(population & pop)
{
	if (m_loci.empty())
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
	for (size_t st = 0; st < m_loci.size(); ++st) {
		int loc = m_loci[st];

		DBG_ASSERT(static_cast<size_t>(loc) < pop.totNumLoci(), IndexError,
			"Index out of range of 0 ~ " + toStr(pop.totNumLoci() - 1));

		// get all available alleles
		vectori alleles; // = m_alleleFreq.alleles(pop, loc);

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
				p_i[sp] = 0; //m_alleleFreq.alleleFreq(pop, *ale, loc, sp);

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
				h_bar += 0; //m_heteroFreq.heteroFreq(pop, *ale, loc, sp) * n_i[sp];
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

		HWE[loc] = calcHWE(vectorlu(0)); //m_genoFreq.countGenotype(pop, loc, 0));
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
			HWE[loc] = calcHWE(vectorlu(0)); //m_genoFreq.countGenotype(pop, loc, sp, 0));
		}
		string varname = subPopVar_String(sp, HWE_String);
		pop.setDoubleVectorVar(varname, HWE);
	}

	return true;
}


}

