/**
 *  $File: BaseMutator.cpp $
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

#include "mutator.h"

#if PY_VERSION_HEX >= 0x03000000
#  define PyInt_FromLong(x) PyLong_FromLong(x)
#endif

namespace simuPOP {

double BaseMutator::mutRate(size_t loc) const
{
	DBG_FAILIF(m_rates.empty(), ValueError, "Please specify mutation rate.");
	if (m_loci.allAvail()) {
		DBG_ASSERT(m_rates.size() == 1 || m_rates.size() > loc,
			IndexError, "Locus index out of range when retrieving mutation rate.");
		return m_rates.size() == 1 ? m_rates[0] : m_rates[loc];
	}
	// without a population reference, we assume that m_loci has already been filled
	// with correct loci.
	const vectoru & loci = m_loci.elems(NULL);
	vectoru::const_iterator it = find(loci.begin(), loci.end(), loc);

	DBG_ASSERT(it != loci.end(), RuntimeError,
		"Failed to find mutation rate for locus " + toStr(loc));
	DBG_ASSERT(m_rates.size() == loci.size(), SystemError,
		"Incorrect rate and loci size");
	return m_rates[it - loci.begin()];
}


void BaseMutator::fillContext(const Population & pop, IndAlleleIterator ptr, size_t locus) const
{
	// chromosome?
	size_t chrom = pop.chromLocusPair(locus).first;
	size_t beg = pop.chromBegin(chrom);
	size_t end = pop.chromEnd(chrom);
	size_t cnt = m_context.size() / 2;

	for (size_t i = 0; i < cnt; ++i) {
		if (locus >= beg + i)
			m_context[i] = *(ptr.ptr() - (cnt - i));
		else
			m_context[i] = InvalidValue;
	}
	for (size_t i = 0; i < cnt; ++i) {
		if (locus + i < end)
			m_context[cnt + i] = *(ptr.ptr() + i + 1);
		else
			m_context[cnt + i] = InvalidValue;
	}
	if (!m_mapIn.empty() || m_mapIn.func().isValid()) {
		for (size_t i = 0; i < m_context.size(); ++i) {
			if (m_context[i] == InvalidValue)
				continue;
			vectoru const & mapInList = m_mapIn.elems();
			if (mapInList.size() > 0) {
				if (static_cast<UINT>(m_context[i]) < mapInList.size())
					m_context[i] = mapInList[m_context[i]];
			} else {
				m_context[i] = m_mapIn.func() (PyObj_As_Int, "(i)",
					m_context[i]);
			}
		}
	}
}


bool BaseMutator::apply(Population & pop) const
{
	DBG_DO(DBG_MUTATOR, cerr << "Mutate replicate " << pop.rep() << endl);

#ifdef LINEAGE
	DBG_WARNIF(infoSize() > 0 && !pop.hasInfoField(infoField(0)),
		"Specified information field " + infoField(0) + " does not exist.");
	bool assignLineage = infoSize() > 0 && pop.hasInfoField(infoField(0));
	size_t lineageIdx = assignLineage ? pop.infoIdx(infoField(0)) : 0;
	DBG_DO(DBG_MUTATOR, cerr << (assignLineage ? "Assign lineage using field " + infoField(0) :
		                         "Not assigning lineage (number of info fields: " + toStr(infoSize()) + ")") << endl);
#endif

	// mapIn and mapOut
	bool mapIn = !m_mapIn.empty() || m_mapIn.func().isValid();
	vectoru const & mapInList = m_mapIn.elems();
	pyFunc mapInFunc = m_mapIn.func();
	size_t numMapInAllele = mapInList.size();
	bool mapOut = !m_mapOut.empty() || m_mapOut.func().isValid();
	vectoru const & mapOutList = m_mapOut.elems();
	size_t numMapOutAllele = mapOutList.size();
	pyFunc mapOutFunc = m_mapOut.func();
	// mutate each mutable locus

	subPopList subPops = applicableSubPops(pop);

	Bernullitrials bt(getRNG());

	DBG_FAILIF(m_rates.empty(), ValueError, "Please specify mutation rate or rates.");
	// all use the same rate
	vectorf rates = m_rates;
	bool rare = true;
	if (rates.size() == 1) {
		if (rates[0] > 1e02)
			rare = false;
		rates.resize(m_loci.allAvail() ? pop.totNumLoci() : m_loci.size());
		fill(rates.begin() + 1, rates.end(), rates[0]);
	} else {
		for (size_t i = 0; i < rates.size(); ++i)
			if (rates[i] > 1e-2) {
				rare = false;
				break;
			}
	}
	// multiple (virtual) subpopulations
	for (size_t idx = 0; idx < subPops.size(); ++idx) {
		size_t sp = subPops[idx].subPop();

		// fromSubPops out of range....
		DBG_FAILIF(sp >= pop.numSubPop(), IndexError,
			"Subpopulation index " + toStr(sp) + " out of range");

		size_t popSize = pop.subPopSize(subPops[idx]);
		DBG_DO(DBG_MUTATOR, cerr << "SP " << subPops[idx] << " size: " << popSize << endl);
		if (popSize == 0)
			continue;

		if (subPops[idx].isVirtual())
			pop.activateVirtualSubPop(subPops[idx]);

		size_t max_pos = pop.ploidy() * popSize;
		if (!rare) {
			bt.setParameter(rates, max_pos);
			bt.doTrial();
		}
		const vectoru & loci = m_loci.elems(&pop);
		size_t iEnd = m_loci.allAvail() ? pop.totNumLoci() : loci.size();
		for (size_t i = 0; i < iEnd; ++i) {
			size_t locus = loci[i];
			DBG_DO(DBG_MUTATOR, cerr << "Mutate at locus " << locus << endl);
			size_t pos = 0;
			if (rare) {
				size_t step = getRNG().randGeometric(rates[i]);
				pos = step == 0 ? Bernullitrials::npos : (step - 1);
			} else
				pos = bt.trialFirstSucc(i);
			size_t lastPos = 0;
			IndAlleleIterator ptr = pop.alleleIterator(locus, sp);
			LINEAGE_EXPR(IndLineageIterator lineagePtr = pop.lineageIterator(locus, sp));
			if (pos != Bernullitrials::npos) {
				do {
#ifdef LINEAGE
					long lineage = 0;
					if (assignLineage) {
						lineagePtr += static_cast<IndLineageIterator::difference_type>(pos - lastPos);
						int sign = m_lineageMode == FROM_INFO ? 1 : (lineagePtr.currentPloidy() % 2 == 0 ? 1 : -1);
						lineage = toLineage(lineagePtr.individual()->info(lineageIdx) * sign);
					}
#endif
					ptr += static_cast<IndAlleleIterator::difference_type>(pos - lastPos);
					lastPos = pos;
					if (!ptr.valid())
						break;
#ifdef MUTANTALLELE
					Allele oldAllele = ptr.value();
#else
					Allele oldAllele = *ptr;
#endif
					(void)oldAllele;  // suppress a warning for unused variable
					DBG_DO(DBG_MUTATOR, cerr << "Allele " << int(oldAllele) << " at locus " << locus);
					if (mapIn) {
						if (numMapInAllele > 0) {
							if (static_cast<size_t>(oldAllele) < numMapInAllele)
								oldAllele = ToAllele(mapInList[oldAllele]);
						} else {
							oldAllele = ToAllele(mapInFunc(PyObj_As_Int, "(i)",
									static_cast<int>(oldAllele)));
						}
					}
					if (!m_context.empty())
						fillContext(pop, ptr, locus);
					// The virtual mutate functions in derived operators will be called.
					Allele newAllele = mutate(oldAllele, locus);
					if (mapOut) {
						if (numMapOutAllele > 0) {
							if (static_cast<size_t>(newAllele) < numMapOutAllele)
								newAllele = ToAllele(mapOutList[newAllele]);
						} else {
							newAllele = ToAllele(mapOutFunc(PyObj_As_Int, "(i)",
									static_cast<int>(newAllele)));
						}
					}
					if (oldAllele != newAllele)
						RefAssign(ptr, newAllele);

					DBG_DO(DBG_MUTATOR, cerr << " is mutated from ");
					DBG_DO(DBG_MUTATOR, cerr << int(oldAllele) << " to " << int(newAllele) << endl);
#ifdef LINEAGE
					if (assignLineage && oldAllele != newAllele) {
						DBG_DO(DBG_MUTATOR, cerr << "Lineage updated from " << *lineagePtr);
						DBG_DO(DBG_MUTATOR, cerr << " to " << lineage << endl);
						*lineagePtr = lineage;
					}
#endif
					if (rare) {
						size_t step = getRNG().randGeometric(rates[i]);
						pos = (step == 0 || step + pos >= max_pos) ? Bernullitrials::npos : (pos + step);
					} else
						pos = bt.trialNextSucc(i, pos);
				} while (pos != Bernullitrials::npos);
			}                                                                                           // succ.any
		}

		if (subPops[idx].isVirtual())
			pop.deactivateVirtualSubPop(sp);
	}   // each subpopulation
	return true;
}


MatrixMutator::MatrixMutator(const floatMatrix & rate,
	const lociList & loci, const uintListFunc & mapIn, const uintListFunc & mapOut,
	const stringFunc & output,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	const stringList & infoFields,
	int lineageMode)
	: BaseMutator(vectorf(1, 0), loci, mapIn, mapOut, 0, output, begin, end, step,
	              at, reps, subPops, infoFields, lineageMode)
{
	matrixf rateMatrix = rate.elems();
	// step 0, determine mu
	double mu = 0;

	for (size_t i = 0; i < rateMatrix.size(); ++i) {
		DBG_ASSERT(rateMatrix[i].size() == rateMatrix.size(), ValueError,
			"A n by n matrix is required.");
		double sum = 0;
		for (size_t j = 0; j < rateMatrix[i].size(); ++j) {
			// ignore p_ii
			if (i == j)
				continue;
			DBG_FAILIF(rateMatrix[i][j] < 0 || rateMatrix[i][j] > 1, ValueError,
				"Elements in a mutation matrix must be between 0 and 1. " + toStr(rateMatrix[i][j]) + " observed.");
			sum += rateMatrix[i][j];
		}
		DBG_FAILIF(sum > 1, ValueError, "Sum of P_ij should not exceed 1");
		if (mu < sum)
			mu = sum;
	}
	DBG_DO(DBG_MUTATOR, cerr << "Mu " << mu << endl);
	setRate(vectorf(1, mu), loci);
	if (mu == 0.)
		return;
	// re-calculate probability
	m_sampler.clear();
	for (size_t i = 0; i < rateMatrix.size(); ++i) {
		double sum = 0;
		for (size_t j = 0; j < rateMatrix[i].size(); ++j) {
			if (i == j)
				continue;
			sum += rateMatrix[i][j];
			rateMatrix[i][j] /= mu;
		}
		rateMatrix[i][i] = 1 - sum / mu;
		DBG_DO(DBG_MUTATOR, cerr << "Setting weight for allele " << i << " to " << rateMatrix[i] << endl);
		m_sampler.push_back(WeightedSampler(rateMatrix[i]));
	}
}


Allele MatrixMutator::mutate(Allele allele, size_t) const
{
	if (static_cast<size_t>(allele) >= m_sampler.size()) {
		DBG_WARNIF(true, "Allele " + toStr(static_cast<size_t>(allele))
			+ " will not be mutated because mutate rates are only defined for alleles 0 ... "
			+ toStr(m_sampler.size() - 1));
		return allele;
	}
	return ToAllele(m_sampler[allele].draw());
}


// mutate to a state other than current state with equal probability
Allele KAlleleMutator::mutate(Allele allele, size_t) const
{
	if (static_cast<size_t>(allele) >= m_k) {
		DBG_WARNIF(true, "Allele " + toStr(static_cast<size_t>(allele))
			+ " will not be mutated because mutate rates are only defined for alleles 0 ... "
			+ toStr(m_k - 1));
		return allele;
	}
#ifdef BINARYALLELE
	return !allele;
#else
	Allele new_allele = static_cast<Allele>(getRNG().randInt(m_k - 1));
	return (new_allele >= allele) ? new_allele + 1 : new_allele;
#endif
}


StepwiseMutator::StepwiseMutator(const floatList & rates, const lociList & loci,
	double incProb, UINT maxAllele, const floatListFunc & mutStep,
	const uintListFunc & mapIn, const uintListFunc & mapOut, const stringFunc & output,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	const stringList & infoFields, int lineageMode)
	: BaseMutator(rates, loci, mapIn, mapOut, 0, output, begin, end,
	              step, at, reps, subPops, infoFields, lineageMode),
	m_incProb(incProb), m_maxAllele(maxAllele), m_mutStep(mutStep)
{
#ifdef BINARYALLELE
	DBG_WARNIF(true, "Symetric stepwise mutation does not work well on two state alleles.");
#endif
	DBG_ASSERT(fcmp_ge(m_incProb, 0.) && fcmp_le(m_incProb, 1.),
		ValueError, "Inc probability should be between [0,1], given " + toStr(m_incProb));

	if (m_maxAllele == 0)
		m_maxAllele = ModuleMaxAllele;
	if (static_cast<ULONG>(m_maxAllele) > ModuleMaxAllele)
		throw ValueError("maxAllele exceeds maximum allowed allele in this module.");

	DBG_FAILIF(!m_mutStep.func().isValid() && m_mutStep.empty(), ValueError,
		"Parameter mutStep must be a number, a list or a valid function.");

	DBG_FAILIF(m_mutStep.size() > 1 &&
		(fcmp_lt(m_mutStep[1], 0) || fcmp_gt(m_mutStep[1], 1)), ValueError,
		"Probability for the geometric distribution has to be between 0 and 1");
}


Allele StepwiseMutator::mutate(Allele allele, size_t) const
{
#ifdef BINARYALLELE
	if (getRNG().randUniform() < m_incProb)
		return 1;
	// decrease
	else
		return 0;
#else
	UINT step = 1;

	if (m_mutStep.size() == 1)
		step = static_cast<UINT>(m_mutStep[0]);
	else if (m_mutStep.size() == 2) {
		DBG_ASSERT(static_cast<int>(m_mutStep[0]) == GEOMETRIC_DISTRIBUTION, ValueError,
			"Incorrect mode for generating mutation step.");
		step = getRNG().randGeometric(m_mutStep[1]);
	} else {
		DBG_ASSERT(m_mutStep.func().isValid(), ValueError,
			"Invalid Python function for StepwiseMutator");
		step = m_mutStep.func() (PyObj_As_Int, "(i)", static_cast<int>(allele));
	}
	// increase
	if (getRNG().randUniform() < m_incProb) {
		if (static_cast<UINT>(allele + step) < m_maxAllele)
			return ToAllele(allele + step);
		else
			return ToAllele(m_maxAllele);
	}
	// decrease
	else {
		if (allele > step)
			return ToAllele(allele - step);
		else
			return 0;
	}
#endif
	return allele;
}


Allele PyMutator::mutate(Allele allele, size_t) const
{
	int resInt = 0;

	PyObject * args = PyTuple_New(m_func.numArgs());

	for (size_t i = 0; i < m_func.numArgs(); ++i) {
		const string & arg = m_func.arg(i);
		if (arg == "allele")
			PyTuple_SET_ITEM(args, i, PyInt_FromLong(static_cast<int>(allele)));
		else if (arg == "context") {
			const vectoru & cnt = context();
			PyObject * c = PyTuple_New(cnt.size());
			for (size_t j = 0; j < cnt.size(); ++j)
				PyTuple_SET_ITEM(c, j, PyInt_FromLong(static_cast<int>(cnt[j])));
			PyTuple_SET_ITEM(args, i, c);
		} else {
			DBG_FAILIF(true, ValueError,
				"Only parameters 'allele' and 'context' are acceptable in a user-provided mutation function.");
		}
	}
	resInt = m_func(PyObj_As_Int, args);
	Py_DECREF(args);

#ifdef BINARYALLELE
	DBG_ASSERT(resInt == 0 || resInt == 1, ValueError,
		"Can only mutate to 0 or 1 in binary mode.");
	return resInt != 0;
#else
	DBG_ASSERT(static_cast<unsigned>(resInt) <= ModuleMaxAllele, ValueError,
		"Mutated to an allele greater than maximum allowed allele value");
	return static_cast<Allele>(resInt);
#endif
}


Allele MixedMutator::mutate(Allele allele, size_t locus) const
{
	size_t idx = m_sampler.draw();
	const BaseMutator * mut = dynamic_cast<const BaseMutator *>(m_mutators[idx]);
	double mu = mut->mutRate(locus);

	if (mu == 1.0 || getRNG().randUniform() < mu)
		return mut->mutate(allele, locus);
	return allele;
}


Allele ContextMutator::mutate(Allele allele, size_t locus) const
{
	const vectoru & alleles = context();

	for (size_t i = 0; i < m_contexts.size(); ++i) {
		bool match = true;
		for (size_t j = 0; j < alleles.size(); ++j) {
			if (ToAllele(m_contexts[i][j]) != alleles[j]) {
				match = false;
				break;
			}
		}
		if (match) {
			DBG_DO(DBG_MUTATOR, cerr << "Context " << alleles << " mutator " << i << endl);
			const BaseMutator * mut = dynamic_cast<const BaseMutator *>(m_mutators[i]);
			if (getRNG().randUniform() < mut->mutRate(locus))
				return mut->mutate(allele, locus);
			return allele;
		}
	}
	if (m_contexts.size() + 1 == m_mutators.size()) {
		DBG_DO(DBG_MUTATOR, cerr << "No context found. Use last mutator." << endl);
		const BaseMutator * mut = dynamic_cast<const BaseMutator *>(m_mutators[m_contexts.size()]);
		if (getRNG().randUniform() < mut->mutRate(locus))
			return mut->mutate(allele, locus);
	} else {
		cerr << "Failed to find context " << alleles << endl;
		throw RuntimeError("No match context is found and there is no default mutator");
	}
	return allele;
}


bool PointMutator::apply(Population & pop) const
{
	subPopList subPops = applicableSubPops(pop);

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();

#ifdef LINEAGE
	bool assignLineage = infoSize() > 0 && pop.hasInfoField(infoField(0));
	size_t lineageIdx = assignLineage ? pop.infoIdx(infoField(0)) : 0;
#endif

	for (; sp != spEnd; ++sp) {
		if (sp->isVirtual())
			pop.activateVirtualSubPop(*sp);

		vectoru::const_iterator it = m_inds.begin();
		vectoru::const_iterator itEnd = m_inds.end();
		for (; it != itEnd; ++it) {
			IndIterator ind = pop.indIterator(sp->subPop()) + static_cast<IndIterator::difference_type>(*it);
			if (!ind.valid())
				continue;

			for (size_t p = 0; p < m_ploidy.size(); ++p) {

				int currentPloidy = static_cast<int>(m_ploidy[p]);

#ifdef LINEAGE
				long lineage = 0;
				if (assignLineage) {
					int sign = m_lineageMode == FROM_INFO ? 1 : (currentPloidy % 2 == 0 ? 1 : -1);
					lineage = toLineage(ind->info(lineageIdx) * sign);
				}
#endif
				if (m_loci.allAvail()) {
					for (size_t i = 0; i < pop.totNumLoci(); ++i) {
#ifdef LINEAGE
						if (assignLineage && m_allele != ind->allele(i, currentPloidy))
							ind->setAlleleLineage(lineage, i, currentPloidy);
#endif
						ind->setAllele(m_allele, i, currentPloidy);
						DBG_DO(DBG_MUTATOR,
							cerr	<< "Mutate locus " << i << " at ploidy " << currentPloidy
							        << " to allele " << int(m_allele) << " at generation "
							        << pop.gen() << endl);
					}
				} else {
					const vectoru & loci = m_loci.elems(&pop);
					for (size_t i = 0; i < loci.size(); ++i) {
#ifdef LINEAGE
						if (assignLineage && m_allele != ind->allele(loci[i], currentPloidy))
							ind->setAlleleLineage(lineage, loci[i], currentPloidy);
#endif
						ind->setAllele(m_allele, loci[i], currentPloidy);
						DBG_DO(DBG_MUTATOR,
							cerr	<< "Mutate locus " << loci[i] << " at ploidy " << currentPloidy
							        << " to allele " << int(m_allele) << " at generation "
							        << pop.gen() << endl);
					}
				}
			}   // ploidy
		}       // Individual
		if (sp->isVirtual())
			pop.deactivateVirtualSubPop(sp->subPop());
	}       // subpopulation

	return true;
}


bool RevertFixedSites::apply(Population & pop) const
{
	if (pop.popSize() == 0 || pop.totNumLoci() == 0)
		return true;

	const vectoru & loci = m_loci.elems(&pop);
	if (loci.size() == 0)
		return true;

	const subPopList subPops = applicableSubPops(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();

	for (; sp != spEnd; ++sp) {
		if (sp->isVirtual())
			pop.activateVirtualSubPop(*sp);

#ifdef MUTANTALLELE
		size_t nLoci = pop.totNumLoci();
		// find the first guy ...
		IndIterator ind = pop.indIterator(sp->subPop());
		// index on the first ploidy, need to get loci index
		GenoIterator s = ind->genoBegin();
		GenoIterator e = s + pop.totNumLoci();
		vectorm::val_iterator s_it = s.get_val_iterator();
		vectorm::val_iterator e_it = e.get_val_iterator();
		vectoru idx;
		for (; s_it != e_it; ++s_it)
			idx.push_back(s_it->first % nLoci);
		// now, find the intersection of loci
		vectoru myloci;
		std::set_intersection(idx.begin(), idx.end(),
			loci.begin(), loci.end(), std::back_inserter(myloci));
		for (size_t idx = 0; idx < myloci.size(); ++idx) {
			size_t loc = myloci[idx];
#else
		for (size_t idx = 0; idx < loci.size(); ++idx) {
			size_t loc = loci[idx];
#endif
			bool fixed = true;
			ConstIndAlleleIterator a = const_cast<const Population &>(pop).alleleIterator(loc, sp->subPop());
			for (; a.valid(); ++a) {
				if (!DerefAllele(a)) {
					fixed = false;
					break;
				}
			}
			// revert fixed allele
			if (fixed) {
				IndAlleleIterator a = pop.alleleIterator(loc, sp->subPop());
				for (; a.valid(); ++a)
					RefAssign(a, 0);
			}
		}

		if (sp->isVirtual())
			pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


size_t FiniteSitesMutator::locateVacantLocus(Population & /* pop */, size_t beg, size_t end, std::set<size_t> & mutants) const
{
	// FIXME: IGNORE this for now
	// this get a random locations
	size_t loc = getRNG().randInt(static_cast<ULONG>(end - beg)) + beg;

	// see if it exists in existing mutants.
	std::set<size_t>::iterator it = std::find(mutants.begin(), mutants.end(), loc);
	if (it == mutants.end())
		return loc;
	//
	// FIXME:
	//
	//   assume mutants have   100 105 106 107
	//   and loc = 106,
	//  goes back and 104, and return
	//
	// look forward and backward
	size_t loc1 = loc + 1;
	std::set<size_t>::iterator it1(it);
	++it1;
	for (; it1 != mutants.end() && loc1 != end; ++it1, ++loc1) {
		if (*it1 != loc1)
			return loc1;
	}
	size_t loc2 = loc - 1;
	std::set<size_t>::reverse_iterator it2(it);
	--it2;
	for (; it2 != mutants.rend() && loc2 != beg; --it2, --loc2) {
		if (*it2 != loc2)
			return loc2;
	}
	// still cannot find
	return 0;
}


bool FiniteSitesMutator::apply(Population & pop) const
{
	// FIXME:
	//
	const matrixi & ranges = m_ranges.elems();
	vectoru width(ranges.size());

	width[0] = ranges[0][1] - ranges[0][0];
	for (size_t i = 1; i < width.size(); ++i)
		width[i] = ranges[i][1] - ranges[i][0] + width[i - 1];
	//
	// FIXME: width[i] is the accumulative number of loci on each chromosome
	//
	// for three ranges(chromsomes) with loci: 500, 1000, 2000
	// width[i] = 0, 500, 1500, 35000
	//
	//  ploidy with = totalNumLoci()
	//size_t ploidyWidth = width.back();
	//
	//size_t indWidth = pop.ploidy() * ploidyWidth;

	size_t indWidth = pop.genoSize();
	size_t ploidyWidth = pop.totNumLoci();

	ostream * out = NULL;
	if (!noOutput())
		out = &getOstream(pop.dict());

	// build a set of existing mutants
	std::set<size_t> mutants;
	//bool saturated = false;

	subPopList subPops = applicableSubPops(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator spEnd = subPops.end();
	for (; sp != spEnd; ++sp) {
		DBG_FAILIF(sp->isVirtual(), ValueError, "This operator does not support virtual subpopulation.");
		for (size_t indIndex = 0; indIndex < pop.subPopSize(sp->subPop()); ++indIndex) {
			size_t loc = 0;
			while (true) {
				// using a geometric distribution to determine mutants
				loc += getRNG().randGeometric(m_rate);
				if (loc > indWidth)
					break;
				Individual & ind = pop.individual(indIndex);
				size_t p = (loc - 1) / ploidyWidth;
				// chromosome and position on chromosome?
				size_t mutLoc = (loc - 1) - p * ploidyWidth;
				size_t ch = 0;
				for (size_t reg = 0; reg < width.size(); ++reg) {
					if (mutLoc < width[reg]) {
						ch = reg;
						break;
					}
				}
				mutLoc += ranges[ch][0];
				if (ch > 0)
					mutLoc -= width[ch - 1];

				// FIXME: as the first step, ignore model 2, so you do not have to figure out
				// the function locateVacantLocus...
				/*
				            // mutLoci is the location of the mutant, for individual ..., chromsome ...
				            if (m_model == 2) {
				                // under an infinite-site model
				                if (saturated) {
				                    if (out)
				                        (*out)	<< pop.gen() << '\t' << mutLoc << '\t' << indIndex
				                                << "\t3\n";
				                    continue;
				                }
				                bool ok = false;
				                //
				                // if the first time
				                if (mutants.empty()) {
				                    // first try our luck...
				                    // FIXME:
				                    // Here we are checing all genotypes (mutant), if mutLoc exists.
				                    // for the new module, yo uneed to check pop.mutBegin()??? (iterate through all
				                    // loci with non-zero allele)
				                    ok = find(pop.genoBegin(false), pop.genoEnd(false), ToAllele(mutLoc)) == pop.genoEnd(false);
				                    if (!ok) {
				                        std::set<size_t> existing(pop.genoBegin(false), pop.genoEnd(false));
				                        mutants.swap(existing);
				                        mutants.erase(0);
				                        saturated = mutants.size() == ploidyWidth;
				                        if (saturated)
				                            cerr << "Failed to introduce new mutants at generation " << pop.gen() << " because all loci have existing mutants." << endl;
				                    }
				                }
				                if (!ok && mutants.find(mutLoc) != mutants.end()) {

				                    size_t newLoc = locateVacantLocus(pop, ranges[ch][0], ranges[ch][1], mutants);
				                    // nothing is found
				                    if (out)
				                        (*out)	<< pop.gen() << '\t' << mutLoc << '\t' << indIndex
				                                << (newLoc == 0 ? "\t3\n" : "\t2\n");
				                    if (newLoc != 0)
				                        mutLoc = newLoc;
				                    else {
				                        cerr << "Failed to introduce a new mutant at generation " << pop.gen() << " because all loci have existing mutants." << endl;
				                        // ignore this mutation, and subsequent mutations...
				                        saturated = true;
				                        continue;
				                    }
				                    // if there is no existing mutant, new mutant is allowed
				                }
				                mutants.insert(mutLoc);
				            }
				 */
				ind.setAllele(1, mutLoc, int(p), int(ch));
			}   // while
		}       // each individual
	}           // each subpopulation
	if (out)
		closeOstream();
	return true;
}


}
