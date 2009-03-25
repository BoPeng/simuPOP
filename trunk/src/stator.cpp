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
		for (UINT i = 0; i < pop.numSubPop(); ++i)
			subPops.push_back(i);

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for ( ; sp != spEnd; ++sp) {
		if (sp->isVirtual())
			pop.activateVirtualSubPop(*sp, IteratableInds);
		IndIterator ind = const_cast<population &>(pop).indBegin(sp->subPop(), sp->isVirtual() ? IteratableInds : AllInds);
		IndIterator indEnd = const_cast<population &>(pop).indEnd(sp->subPop(), sp->isVirtual() ? IteratableInds : AllInds);
		for (; ind != indEnd; ++ind) {
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
		for (UINT i = 0; i < pop.numSubPop(); ++i)
			subPops.push_back(i);

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
		IndIterator ind = const_cast<population &>(pop).indBegin(sp->subPop(), sp->isVirtual() ? IteratableInds : AllInds);
		IndIterator indEnd = const_cast<population &>(pop).indEnd(sp->subPop(), sp->isVirtual() ? IteratableInds : AllInds);
		for (; ind != indEnd; ++ind) {
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
	strDict numOfMale_param,
	//
	bool numOfAffected,
	strDict numOfAffected_param,
	//
	const uintList & numOfAlleles,
	strDict numOfAlleles_param,
	//
	vectori alleleFreq,
	strDict alleleFreq_param,
	//
	vectori heteroFreq,
	vectori expHetero,
	strDict expHetero_param,
	//
	vectori homoFreq,
	vectori genoFreq,
	strDict genoFreq_param,
	intMatrix haploFreq,
	//
	intMatrix LD,
	strDict LD_param,
	//
	//intMatrix association,
	//strDict association_param,
	//
	vectori Fst,
	strDict Fst_param,
	//
	const uintList & HWE,
	//
	bool hasPhase,
	bool midValues,                                                 // this parameter will be removed after all _param parameter is given.
	// regular parameters
	const stringFunc & output,
	int stage, int begin, int end, int step, const intList & at,
	const repList & rep, const subPopList & subPops, const vectorstr & infoFields)
	: baseOperator("", stage, begin, end, step, at, rep, subPops, infoFields),
	// the order of initialization is meaningful since they may depend on each other
	m_popSize(popSize),
	m_numOfMale(numOfMale, numOfMale_param),
	m_numOfAffected(numOfAffected, numOfAffected_param),
	m_alleleFreq(alleleFreq, alleleFreq_param),
	m_numOfAlleles(m_alleleFreq, numOfAlleles.elems(), numOfAlleles_param),
	m_heteroFreq(heteroFreq, homoFreq),
	m_expHetero(m_alleleFreq, expHetero, expHetero_param),
	m_genoFreq(genoFreq, genoFreq_param),
	m_haploFreq(haploFreq),
	m_LD(m_alleleFreq, m_haploFreq, LD, LD_param),
	//m_association(m_alleleFreq, m_haploFreq, association, association_param),
	m_Fst(m_alleleFreq, m_heteroFreq, Fst, Fst_param),
	m_HWE(m_genoFreq, HWE.elems())
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
	       //m_association.apply(pop) &&
	       m_Fst.apply(pop) &&
	       m_HWE.apply(pop);
}


bool statPopSize::apply(population & pop)
{
	if (!m_isActive)
		return true;

	UINT numSP = pop.numSubPop();
	ULONG popSize = pop.popSize();

	pop.setIntVar(numSubPop_String, numSP);
	pop.setIntVar(popSize_String, popSize);

	// type mismatch, can not use subPopSizes() directly.
	vectori spSize(numSP);
	vectori vspSize;
	for (size_t sp = 0; sp < numSP; ++sp) {
		spSize[sp] = pop.subPopSize(sp);
		size_t numVSP = pop.numVirtualSubPop();
		if (numVSP == 0)
			pop.setIntVar(virtualPopSize_String + string("[") + toStr(sp) + "]", spSize[sp]);
		else {
			vspSize.clear();
			for (size_t vsp = 0; vsp < numVSP; ++vsp)
				vspSize.push_back(pop.subPopSize(vspID(sp, vsp)));
			pop.setIntVectorVar(virtualPopSize_String + string("[") + toStr(sp) + "]", vspSize);
		}
	}

	pop.setIntVectorVar(subPopSize_String, spSize);

	for (size_t sp = 0; sp < numSP; ++sp)
		pop.setIntVar(subPopVar_String(sp, popSize_String), spSize[sp]);
	return true;
}


bool statNumOfMale::apply(population & pop)
{
	if (m_numOfMale.empty())
		return true;

	UINT numSP = pop.numSubPop();
	m_numOfMale.resize(numSP + 1);
	m_numOfFemale.resize(numSP + 1);

	ULONG numOfMale = 0;
	for (size_t sp = 0; sp < numSP; ++sp) {
		ULONG n = 0;
		for (IndIterator it = pop.indBegin(sp); it.valid(); ++it) {
			if (it->sex() == Male)
				n++;
		}
		numOfMale += n;
		m_numOfMale[sp] = n;

		if (m_evalInSubPop) {
			if (m_output_numOfMale)
				pop.setIntVar(subPopVar_String(sp, numOfMale_String), n);
			if (m_output_propOfMale)
				pop.setDoubleVar(subPopVar_String(sp, propOfMale_String),
					(double)(n) / pop.subPopSize(sp));
		}

		n = pop.subPopSize(sp) - n;
		m_numOfFemale[sp] = n;

		if (m_evalInSubPop) {
			if (m_output_numOfFemale)
				pop.setIntVar(subPopVar_String(sp, numOfFemale_String), n);
			if (m_output_propOfFemale)
				pop.setDoubleVar(subPopVar_String(sp, propOfFemale_String),
					(double)(n) / pop.subPopSize(sp));
		}
	}
	if (m_output_numOfMale)
		pop.setIntVar(numOfMale_String, numOfMale);
	if (m_output_numOfFemale)
		pop.setIntVar(numOfFemale_String, pop.popSize() - numOfMale);
	if (m_output_propOfMale)
		pop.setDoubleVar(propOfMale_String, (double)(numOfMale) / pop.popSize());
	if (m_output_propOfFemale)
		pop.setDoubleVar(propOfFemale_String, (double)(pop.popSize() - numOfMale) / pop.popSize());
	m_numOfMale[numSP] = numOfMale;
	m_numOfFemale[numSP] = pop.popSize() - numOfMale;
	return true;
}


bool statNumOfAffected::apply(population & pop)
{
	if (m_numOfAffected.empty() )
		return true;

	ULONG numOfAffected = 0;
	UINT numSP = pop.numSubPop();
	m_numOfAffected.resize(numSP + 1);
	m_numOfUnaffected.resize(numSP + 1);
	if (m_evalInSubPop) {
		for (size_t sp = 0; sp < numSP; ++sp) {
			ULONG n = count_if(pop.indBegin(sp), pop.indEnd(sp),
				isAffected<individual>());
			numOfAffected += n;
			m_numOfAffected[sp] = n;
			if (m_output_numOfAffected)
				pop.setIntVar(subPopVar_String(sp, numOfAffected_String), n);
			if (m_output_propOfAffected)
				pop.setDoubleVar(subPopVar_String(sp, propOfAffected_String),
					(double)(n) / pop.subPopSize(sp));
			if (m_output_numOfUnaffected || m_output_propOfUnaffected) {
				n = pop.subPopSize(sp) - n;
				m_numOfUnaffected[sp] = n;
			}
			if (m_output_numOfUnaffected)
				pop.setIntVar(subPopVar_String(sp, numOfUnaffected_String), n);
			if (m_output_propOfUnaffected)
				pop.setDoubleVar(subPopVar_String(sp, propOfUnaffected_String),
					(double)(n) / pop.subPopSize(sp));
		}
	}
	if (m_output_numOfAffected)
		pop.setIntVar(numOfAffected_String, numOfAffected);
	if (m_output_numOfUnaffected)
		pop.setIntVar(numOfUnaffected_String, pop.popSize() - numOfAffected);
	if (m_output_propOfAffected)
		pop.setDoubleVar(propOfAffected_String, (double)(numOfAffected) / pop.popSize());
	if (m_output_propOfUnaffected)
		pop.setDoubleVar(propOfUnaffected_String, (double)(pop.popSize() - numOfAffected) / pop.popSize());
	m_numOfAffected[numSP] = numOfAffected;
	m_numOfUnaffected[numSP] = pop.popSize() - numOfAffected;
	return true;
}


statAlleleFreq::~statAlleleFreq()
{
}


void statAlleleFreq::addLocus(int locus, bool post, bool subPop, bool numOfAlleles)
{
	vectori::const_iterator it;

	// a new one
	if ( (it = find(m_atLoci.begin(), m_atLoci.end(), locus)) == m_atLoci.end() ) {
		m_atLoci.push_back(locus);
		m_ifPost.push_back(static_cast<int>(post));
	}
	// existing one
	else
		m_ifPost[ it - m_atLoci.begin() ] |= static_cast<int>(post);

	m_evalInSubPop |= subPop;
	m_output_numOfAlleles |= numOfAlleles;
}


bool statAlleleFreq::apply(population & pop)
{
	if (m_atLoci.empty())
		return true;

	UINT numSP = pop.numSubPop();
	UINT numLoci = m_atLoci.size();

	pop.removeVar(NumOfAlleles_String);
	pop.removeVar(AlleleNum_String);
	pop.removeVar(AlleleFreq_String);
	for (UINT sp = 0; sp < numSP; ++sp) {
		pop.removeVar(subPopVar_String(sp, NumOfAlleles_String));
		pop.removeVar(subPopVar_String(sp, AlleleNum_String));
		pop.removeVar(subPopVar_String(sp, AlleleFreq_String));
	}

	UINT len = numSP == 1 ? 1 : (numSP + 1);
	// if not initialized or m_atLoci/numSP changes
	if (m_alleleNum.size() != len) {
		for (size_t i = 0; i < numLoci;  ++i)
			DBG_FAILIF(static_cast<UINT>(m_atLoci[i]) >= pop.totNumLoci(),
				IndexError, "locus index (" + toStr(m_atLoci[i])
				+ ") out of range of 0 - " + toStr(pop.totNumLoci() - 1));

		m_alleleNum.resize(len);
		m_alleleFreq.resize(len);

		m_numOfAlleles.resize(len);
		for (size_t i = 0; i < len; ++i) {
			m_alleleNum[i].resize(pop.totNumLoci());
			m_alleleFreq[i].resize(pop.totNumLoci());
			m_numOfAlleles[i].resize(pop.totNumLoci());
			fill(m_numOfAlleles[i].begin(), m_numOfAlleles[i].end(), 0);
		}
	}

	string varname;

	for (size_t i = 0; i < numLoci; ++i) {
		UINT loc = m_atLoci[i];

		vectori & sum = m_alleleNum.back()[loc];
		fill(sum.begin(), sum.end(), 0);
		ULONG sumAll = 0;

		// for each subpopulation
		for (UINT sp = 0; sp < numSP;  ++sp) {
			vectori & num = m_alleleNum[sp][loc];
			// clear all current values
			fill(num.begin(), num.end(), 0);
			// for convenience, gurantees the existence
			// of num for 0 and 1...
			if (num.size() < 2)
				num.resize(2, 0);

			// go through all alleles
			IndAlleleIterator a = pop.alleleBegin(loc, sp);
			IndAlleleIterator aEnd = pop.alleleEnd(loc, sp);
			// use allAllelel here because some marker does not have full number
			// of alleles (e.g. markers on chromosome X and Y).
			ULONG allAllele = 0;
			for (; a != aEnd; ++a) {
				if (AlleleUnsigned(*a) >= num.size())
					num.resize(*a + 1, 0);
				num[*a]++;
				allAllele++;
			}

			// add this number to overall num
			// calculate frequency
			// if there is only one sp, no need to do so.
			if (numSP > 1) {
				if (sum.size() < num.size())
					sum.resize(num.size(), 0);
				for (size_t e = 0, eEnd = num.size(); e < eEnd; ++e)
					sum[e] += num[e];
				sumAll += allAllele;
			}

			vectorf & freq = m_alleleFreq[sp][loc];
			freq.resize(num.size(), 0.);
			for (size_t e = 0, eEnd = num.size(); e < eEnd; ++e)
				freq[e] = allAllele == 0 ? 0 : static_cast<double>(num[e]) / allAllele;

			// post result at this locus
			if (m_ifPost[i]) {
				if (m_output_alleleNum) {
					varname = subPopVar_String(sp, AlleleNum_String) + "[" + toStr(loc) + "]";
					PyObject * d = pop.setIntVectorVar(varname, num);

					// do not need a separate result
					if (numSP == 1 && sp == 0) {
						varname = toStr(AlleleNum_String) + "[" + toStr(loc) + "]";
						Py_INCREF(d);
						pop.setVar(varname, d);
					}
				}

				if (m_output_alleleFreq) {
					varname = subPopVar_String(sp, AlleleFreq_String) + "[" + toStr(loc) + "]";
					PyObject * d = pop.setDoubleVectorVar(varname, freq);

					// do not need a separate result
					if (numSP == 1 && sp == 0) {
						varname = toStr(AlleleFreq_String) + "[" + toStr(loc) + "]";
						Py_INCREF(d);
						pop.setVar(varname, d);
					}
				}
			}                                                                         // post

			// set numOfAlleles if necessary
			if (m_output_numOfAlleles)
				m_numOfAlleles[sp][loc] = count_if(num.begin(), num.end(),
					bind2nd(std::greater<int>(), 0));
		}                                                                                       // subpop

		if (numSP > 1) {                                                                        // calculate sum and post overall result
			// summary?
			vectorf & freq = m_alleleFreq.back()[loc];
			freq.resize(sum.size());
			for (size_t e = 0, eEnd = sum.size(); e < eEnd; ++e)
				freq[e] = sumAll == 0 ? 0 : static_cast<double>(sum[e]) / sumAll;

			if (m_ifPost[i]) {
				if (m_output_alleleNum) {
					varname = string(AlleleNum_String) + "[" + toStr(loc) + "]";
					pop.setIntVectorVar(varname, sum);
				}
				if (m_output_alleleFreq) {
					varname = string(AlleleFreq_String) + "[" + toStr(loc) + "]";
					pop.setDoubleVectorVar(varname, freq);
				}
			}

			// set numOfAlleles if necessary
			if (m_output_numOfAlleles)
				m_numOfAlleles.back()[loc] = count_if(sum.begin(), sum.end(),
					bind2nd(std::greater<int>(), 0));
		}
	}                                                                                         // all loci

	// post number of alleles
	if (m_output_numOfAlleles && accumulate(m_ifPost.begin(), m_ifPost.end(), 0) > 0) {
		// post number of alleles
		for (UINT sp = 0; sp < numSP; ++sp) {
			PyObject * d = pop.setIntVectorVar(subPopVar_String(sp, NumOfAlleles_String),
				m_numOfAlleles[sp]);
			if (numSP == 1) {
				Py_INCREF(d);
				pop.setVar(NumOfAlleles_String, d);
			}
		}
		if (numSP > 1) {
			pop.setIntVectorVar(NumOfAlleles_String, m_numOfAlleles.back());
		}
	}
	return true;
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
			IndAlleleIterator a = pop.alleleBegin(loc, sp);
			IndAlleleIterator aEnd = pop.alleleEnd(loc, sp);
			for (; a != aEnd; a += 2) {
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
			vectorf & af = m_alleleFreq.alleleFreqVec(loc, sp);
			double expHeter = 1;
			// 1-sum pi^2
			for (int al = 0, alEnd = af.size() ; al < alEnd; al++)
				expHeter -= af[al] * af[al];

			m_expHetero[sp][loc] = expHeter;
		}

		vectorf & af = m_alleleFreq.alleleFreqVec(loc);
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


statGenoFreq::statGenoFreq(const vectori & genoFreq,
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

		int loc = m_atLoci[i];

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
			IndIterator it = pop.indBegin(sp);
			IndIterator itEnd = pop.indEnd(sp);
			for (; it != itEnd; ++it) {
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

	IndIterator it = pop.indBegin();
	IndIterator itEnd = pop.indEnd();
	for (; it != itEnd; ++it) {
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

	IndIterator it = pop.indBegin(subPop);
	IndIterator itEnd = pop.indEnd(subPop);
	for (; it != itEnd; ++it) {
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

			IndAlleleIterator it = pop.alleleBegin(0, sp);
			IndAlleleIterator itEnd = pop.alleleEnd(0, sp);
			for (; it != itEnd; ++it) {
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
		m_alleleFreq.addLocus(m_LD[i][0], m_midValues, true, true);
		m_alleleFreq.addLocus(m_LD[i][1], m_midValues, true, true);
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
void statLD::calculateLD(const vectori & hapLoci, const vectori & hapAlleles, UINT sp, bool subPop,
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
		P_A = m_alleleFreq.alleleFreq(hapAlleles[0], hapLoci[0], sp);
		P_B = m_alleleFreq.alleleFreq(hapAlleles[1], hapLoci[1], sp);

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
		P_A = m_alleleFreq.alleleFreq(hapAlleles[0], hapLoci[0]);
		P_B = m_alleleFreq.alleleFreq(hapAlleles[1], hapLoci[1]);

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
	vectori numofalleles = m_alleleFreq.numOfAlleles();
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
			calculateLD(hapLoci, hapAlleles, 0, false, P_A, P_B, D, D_prime, r2, delta2);
			outputLD(pop, hapLoci, haploKey(hapAlleles), 0, false, valid_delta2, D, D_prime, r2, delta2);

			if (m_evalInSubPop) {
				if (numSP == 1)  // use the whole population result
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
						calculateLD(hapLoci, hapAlleles, sp, true, P_A, P_B, D, D_prime, r2, delta2);
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
			vectori A_alleles = m_alleleFreq.alleles(hapLoci[0]);
			vectori B_alleles = m_alleleFreq.alleles(hapLoci[1]);

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
					calculateLD(hapLoci, hapAlleles, 0, false, P_A, P_B, D_, D_prime_, r2_, delta2_);
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
				if (numSP == 1)  // use the whole population result
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
								calculateLD(hapLoci, hapAlleles, sp, true, P_A, P_B, D_, D_prime_, r2_, delta2_);

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
	/* FIXME: the following is moved from the old statAssociation... they should
	 * be merged to the above calculation
	*/
	if (!m_output_ChiSq && !m_output_UCU && !m_output_CramerV)
		return true;

	for (size_t i = 0; i < m_LD.size(); ++i) {
		//
		vectori hapLoci = m_LD[i];
		hapLoci.resize(2);

		// find out all alleles
		vectori A_alleles = m_alleleFreq.alleles(hapLoci[0]);
		vectori B_alleles = m_alleleFreq.alleles(hapLoci[1]);
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
			ChiSq_P = rng().pvalChiSq(ChiSq, (as - 1) * (bs - 1));
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
					ChiSq_P = rng().pvalChiSq(ChiSq, (as - 1) * (bs - 1));
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


statFst::statFst(statAlleleFreq & alleleFreq, statHeteroFreq & heteroFreq,
	const vectori & Fst, const strDict & param)
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
		m_alleleFreq.addLocus(m_atLoci[i], m_midValues, true, true);

		// need to get heterozygous proportion  at this locus
		m_heteroFreq.addLocus(m_atLoci[i], m_midValues);
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
		vectori alleles = m_alleleFreq.alleles(loc);

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
				p_i[sp] = m_alleleFreq.alleleFreq(*ale, loc, sp);

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
	int obsAA = cnt[2]; // in this algorithm, AA is rare.
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

