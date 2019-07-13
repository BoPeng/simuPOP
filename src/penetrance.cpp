/**
 *  $File: penetrance.cpp $
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

#include "penetrance.h"

#include <sstream>
using std::ostringstream;

#if PY_VERSION_HEX >= 0x03000000
#  define PyInt_FromLong(x) PyLong_FromLong(x)
#endif


namespace simuPOP {

// set pentrance to all individuals and record penetrance if requested.
bool BasePenetrance::apply(Population & pop) const
{
	bool savePene = infoSize() > 0;
	size_t infoIdx = 0;

	if (savePene)
		infoIdx = pop.infoIdx(infoField(0));

	vectoru gens = m_ancGens.elems();
	if (m_ancGens.allAvail())
		for (int gen = 0; gen <= pop.ancestralGens(); ++gen)
			gens.push_back(gen);
	else if (m_ancGens.unspecified())
		gens.push_back(pop.curAncestralGen());

	size_t oldGen = pop.curAncestralGen();
	for (unsigned genIdx = 0; genIdx < gens.size(); ++genIdx) {
		pop.useAncestralGen(gens[genIdx]);

		subPopList subPops = applicableSubPops(pop);

		subPopList::const_iterator sp = subPops.begin();
		subPopList::const_iterator spEnd = subPops.end();
		for (; sp != spEnd; ++sp) {
			if (sp->isVirtual())
				pop.activateVirtualSubPop(*sp);

			if (numThreads() > 1 && parallelizable()) {
#pragma omp parallel
				{
#ifdef _OPENMP
					IndIterator ind = pop.indIterator(sp->subPop(), omp_get_thread_num());
					for (; ind.valid(); ++ind) {
						double p = penet(&pop, ind.rawIter());

						if (savePene)
							ind->setInfo(p, infoIdx);

						if (getRNG().randUniform() < p)
							ind->setAffected(true);
						else
							ind->setAffected(false);
					}
#endif
				}
			} else {
				IndIterator ind = pop.indIterator(sp->subPop());
				for (; ind.valid(); ++ind) {
					double p = penet(&pop, ind.rawIter());

					if (savePene)
						ind->setInfo(p, infoIdx);

					if (getRNG().randUniform() < p)
						ind->setAffected(true);
					else
						ind->setAffected(false);
				}
			}
			if (sp->isVirtual())
				pop.deactivateVirtualSubPop(sp->subPop());
		}
	}
	pop.useAncestralGen(oldGen);

	return true;
}


bool BasePenetrance::applyToIndividual(Individual * ind, Population * pop)
{
	double p = penet(pop, pop->rawIndBegin() + (ind - &*pop->rawIndBegin()));

	if (infoSize() > 0)
		ind->setInfo(p, infoField(0));
	bool affected = getRNG().randUniform() < p;
	ind->setAffected(affected);
	return affected;
}


bool BasePenetrance::applyDuringMating(Population & /* pop */, Population & offPop, RawIndIterator offspring,
                                       Individual * /* dad */, Individual * /* mom */) const
{
	// if offspring does not belong to subPops, do nothing, but does not fail.
	if (!applicableToAllOffspring() && !applicableToOffspring(offPop, offspring))
		return true;
	double p = penet(&offPop, offspring);

	if (infoSize() > 0)
		offspring->setInfo(p, infoField(0));
	offspring->setAffected(getRNG().randUniform() < p);
	return true;
}


// this function is the same as MapPenetrance.
double MapPenetrance::penet(Population * /* pop */, RawIndIterator ind) const
{
	vectoru chromTypes;

	const vectoru & loci = m_loci.elems(&*ind);

	for (size_t i = 0; i < loci.size(); ++i)
		chromTypes.push_back(ind->chromType(ind->chromLocusPair(loci[i]).first));

	size_t ply = ind->ploidy();
	if (ind->isHaplodiploid() && ind->sex() == MALE)
		ply = 1;

	vectori alleles;
	alleles.reserve(ply * loci.size());

	for (size_t idx = 0; idx < loci.size(); ++idx) {
		for (size_t p = 0; p < ply; ++p) {
			if (chromTypes[idx] == CHROMOSOME_Y && ind->sex() == FEMALE)
				continue;
			if (((chromTypes[idx] == CHROMOSOME_X && p == 1) ||
			     (chromTypes[idx] == CHROMOSOME_Y && p == 0)) && ind->sex() == MALE)
				continue;
			if (chromTypes[idx] == MITOCHONDRIAL && p > 0)
				continue;
			alleles.push_back(ind->allele(loci[idx], p));
		}
	}

	tupleDict::const_iterator pos = m_dict.find(alleles);

	if (pos != m_dict.end())
		return pos->second;

	if (ply > 1) {
		// try to look up the key without phase
		tupleDict::const_iterator it = m_dict.begin();
		tupleDict::const_iterator itEnd = m_dict.end();
		for (; it != itEnd; ++it) {
			bool ok = true;
			const tupleDict::key_type & key = it->first;
			size_t begin_idx = 0;
			size_t end_idx = 0;
			for (size_t i = 0; i < loci.size(); ++i) {
				if (chromTypes[i] == CHROMOSOME_Y) {
					if (ind->sex() == FEMALE)
						continue;
					else
						++end_idx;
				} else if (chromTypes[i] == CHROMOSOME_X && ind->sex() == MALE)
					++end_idx;
				else if (chromTypes[i] == MITOCHONDRIAL)
					++end_idx;
				else
					end_idx += ply;
				if (key.size() != end_idx - begin_idx) {
					ok = false;
					break;
				}
				if (ply == 2) {
					if ((alleles[begin_idx] != key[0] || alleles[end_idx - 1] != key[1]) &&
					    (alleles[begin_idx] != key[1] || alleles[end_idx - 1] != key[0])) {
						ok = false;
						break;
					}
				} else {
					std::sort(alleles.begin() + begin_idx, alleles.begin() + end_idx);
					tupleDict::key_type sorted_key = it->first;
					std::sort(sorted_key.begin(), sorted_key.end());
					for (size_t j = 0; j < sorted_key.size(); ++j) {
						if (alleles[ply * i + j] != sorted_key[j]) {
							ok = false;
							break;
						}
					}
				}
				begin_idx = end_idx;
			}
			if (ok)
				return it->second;
		}
	}
	// no match
	ostringstream allele_string;
	allele_string << "(";
	for (size_t i = 0; i < alleles.size(); ++i) {
		if (i != 0)
			allele_string << ", ";
		allele_string << alleles[i];
	}
	allele_string << ")";
	throw ValueError("No penetrance value for genotype " + allele_string.str());
	// this line should not be reached.
	return 0;
}


string MaPenetrance::describe(bool /* format */) const
{
	return "<simuPOP.MaPenetrance> multiple-alleles penetrance" ;
}


// this function is the same as MaPenetrance.
double MaPenetrance::penet(Population * /* pop */, RawIndIterator ind) const
{
	UINT index = 0;
	bool singleST = m_wildtype.size() == 1;
	const vectoru & loci = m_loci.elems(&*ind);

	DBG_FAILIF((ind->ploidy() == 2 && m_penetrance.size() != static_cast<UINT>(pow(3., static_cast<double>(loci.size())))) ||
		(ind->ploidy() == 1 && m_penetrance.size() != static_cast<UINT>(pow(2., static_cast<double>(loci.size())))),
		ValueError, "Please specify penetrance for each combination of genotype.");

	for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc) {
		if (ind->ploidy() == 1) {
			Allele a = TO_ALLELE(ind->allele(*loc));
			if (singleST)
				index = index * 2 + (ALLELE_AS_UNSINGED(a) != m_wildtype[0]);
			else
				index = index * 2 + (find(m_wildtype.begin(), m_wildtype.end(), ALLELE_AS_UNSINGED(a)) == m_wildtype.end());
		} else if (ind->ploidy() == 2) {
			// get genotype of ind
			Allele a = TO_ALLELE(ind->allele(*loc, 0));
			Allele b = TO_ALLELE(ind->allele(*loc, 1));

			int numWildtype = 0;

			// count number of wildtype
			// this improve the performance a little bit
			if (singleST) {
				numWildtype = (ALLELE_AS_UNSINGED(a) == m_wildtype[0]) + (ALLELE_AS_UNSINGED(b) == m_wildtype[0]);
			} else {
				if (find(m_wildtype.begin(), m_wildtype.end(), ALLELE_AS_UNSINGED(a)) != m_wildtype.end())
					numWildtype++;

				if (find(m_wildtype.begin(), m_wildtype.end(), ALLELE_AS_UNSINGED(b)) != m_wildtype.end())
					numWildtype++;
			}

			index = index * 3 + 2 - numWildtype;
		} else {
			DBG_FAILIF(true, ValueError, "The MaPenetrance only supports haploid and diploid populations.");
		}
	}
	return m_penetrance[index];
}


/** CPPONLY
 *  An accumulator class to generate overall fitness values from existing
 *  ones.
 */
class PenetranceAccumulator
{
public:
	PenetranceAccumulator(int mode) : m_mode(mode)
	{
		m_value = m_mode == ADDITIVE ? 0 : 1;
	}


	void push(double val)
	{
		if (m_mode == MULTIPLICATIVE)
			m_value *= val;
		else if (m_mode == ADDITIVE)
			m_value += val;
		else if (m_mode == HETEROGENEITY)
			m_value *= 1 - val;
	}


	double value()
	{
		if (m_mode == MULTIPLICATIVE)
			return m_value;
		else if (m_mode == ADDITIVE)
			return m_value > 1 ? 1 : m_value;
		else if (m_mode == HETEROGENEITY)
			return m_value > 1 ? 0 : 1 - m_value;
		else
			throw ValueError("Unrecognized accumulation mode");
		return 0;
	}


private:
	int m_mode;

	double m_value;
};


double MlPenetrance::penet(Population * pop, RawIndIterator ind) const
{
	PenetranceAccumulator p(m_mode);

	vectorop::const_iterator s = m_peneOps.begin();
	vectorop::const_iterator sEnd = m_peneOps.end();

	for (; s != sEnd; ++s) {
		if (pop && (
		            (!(*s)->isActive(pop->rep(), pop->gen())) ||
		            (!(*s)->applicableToAllOffspring() && !(*s)->applicableToOffspring(*pop, ind))))
			continue;
		p.push(dynamic_cast<const BasePenetrance *>(*s)->penet(pop, ind));
	}
	return p.value();
}


// the same as PyPenetrance
double PyPenetrance::penet(Population * pop, RawIndIterator ind) const
{
	PyObject * args = PyTuple_New(m_func.numArgs());

	DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

	for (size_t i = 0; i < m_func.numArgs(); ++i) {
		const string & arg = m_func.arg(i);
		if (arg == "ind")
			PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(&*ind)));
		else if (arg == "geno")
			PyTuple_SET_ITEM(args, i, ind->genoAtLoci(m_loci));
		else if (arg == "mut")
			PyTuple_SET_ITEM(args, i, ind->mutAtLoci(m_loci));
		else if (arg == "gen") {
			DBG_FAILIF(pop == NULL, ValueError, "No valid population reference is passed.");
			PyTuple_SET_ITEM(args, i, PyInt_FromLong(static_cast<long>(pop->gen())));
		} else if (arg == "pop") {
			DBG_FAILIF(pop == NULL, ValueError, "No valid population reference is passed.");
			PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(pop)));
		} else {
			DBG_FAILIF(!ind->hasInfoField(arg), ValueError,
				"Only parameters 'ind', 'geno', 'gen', 'pop' and names of information fields are "
				"acceptable in function " + m_func.name());
			PyTuple_SET_ITEM(args, i, PyFloat_FromDouble(ind->info(arg)));
		}
	}

	double penetrance = m_func(PyObj_As_Double, args);
	Py_XDECREF(args);
	return penetrance;
}


PyMlPenetrance::PyMlPenetrance(PyObject * func, int mode, const lociList & loci,
	const uintList & ancGens,
	const stringFunc & /* output */, int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops, const stringList & infoFields) :
	BasePenetrance(ancGens, begin, end, step, at, reps, subPops, infoFields),
	m_func(func), m_mode(mode), m_loci(loci), m_searchMode(0),
	m_penetFactory()
{
	DBG_ASSERT(m_func.isValid(), ValueError, "Passed variable is not a callable python function.");

	if (m_func.numArgs() == 0)
		m_searchMode = 10;
	else if (m_func.numArgs() == 1) {
		if (m_func.arg(0) == "loc")
			m_searchMode = 11;
		else if (m_func.arg(0) == "alleles")
			m_searchMode = 12;
		else
			throw ValueError("Invalid parameter for passed function " + m_func.name());
	} else if (m_func.numArgs() == 2) {
		if (m_func.arg(0) == "alleles" && m_func.arg(1) == "loc")
			m_searchMode = 13;
		else if (m_func.arg(0) == "loc" && m_func.arg(1) == "alleles")
			m_searchMode = 14;
		else
			throw ValueError("Invalid parameters for passed function " + m_func.name() +
				" (allele1 and allele2 have to be specified together, in that order)");
	} else {
		DBG_FAILIF(true, ValueError,
			"Python function for m_func only accepts paramters loc and alleles");
	}
}


double PyMlPenetrance::penet(Population * /* pop */, RawIndIterator ind) const
{
	PenetranceAccumulator pnt(m_mode);
	const vectoru & loci = m_loci.elems(&*ind);

	size_t ply = ind->ploidy();

	if (ind->isHaplodiploid() && ind->sex() == MALE)
		ply = 1;

#ifdef MUTANTALLELE
	if (ply == 1) {
		vectora geno(1, TO_ALLELE(0));
		size_t numLoci = ind->totNumLoci();
		GenoIterator it = ind->genoBegin();
		GenoIterator it_end = ind->genoEnd();
		vectorm::val_iterator iit = it.get_val_iterator();
		vectorm::val_iterator iit_end = it_end.get_val_iterator();
		for (; iit != iit_end; ) {
			if (iit->second != 0 &&
			    (m_loci.allAvail() || find(loci.begin(), loci.end(), (iit->first) % numLoci) != loci.end())) {
				geno[0] = iit->second;
				pnt.push(getPenetranceValue(LocGenotype((iit->first) % numLoci, geno)));
				++iit;
			}
		}
	} else if (ply == 2) {
		size_t numLoci = ind->totNumLoci();
		GenoIterator it0 = ind->genoBegin(0);
		GenoIterator it0_end = ind->genoEnd(0);
		GenoIterator it1 = ind->genoBegin(1);
		GenoIterator it1_end = ind->genoEnd(1);
		vectorm::val_iterator iit0 = it0.get_val_iterator();
		vectorm::val_iterator iit0_end = it0_end.get_val_iterator();
		vectorm::val_iterator iit1 = it1.get_val_iterator();
		vectorm::val_iterator iit1_end = it1_end.get_val_iterator();
		vectora geno(2, TO_ALLELE(0));
		for (; iit0 != iit0_end || iit1 != iit1_end; ) {
			if (iit1 == iit1_end || iit0->first + numLoci < iit1->first) {
				if (iit0->second != 0 &&
				    (m_loci.allAvail() || find(loci.begin(), loci.end(), (iit0->first) % numLoci) != loci.end())) {
					geno[0] = iit0->second;
					geno[1] = 0;
					pnt.push(getPenetranceValue(LocGenotype((iit0->first) % numLoci, geno)));
				}
				++iit0;
			} else if (iit0 == iit0_end || iit0->first + numLoci > iit1->first) {
				if (iit1->second != 0 &&
				    (m_loci.allAvail() || find(loci.begin(), loci.end(), (iit1->first) % numLoci) != loci.end())) {
					geno[0] = 0;
					geno[1] = iit1->second;
					pnt.push(getPenetranceValue(LocGenotype((iit1->first) % numLoci, geno)));
				}
				++iit1;
			} else {
				if ((iit0->second != 0 || iit1->second != 0) &&
				    (m_loci.allAvail() || find(loci.begin(), loci.end(), (iit1->first) % numLoci) != loci.end())) {
					geno[0] = iit0->second;
					geno[1] = iit1->second;
					pnt.push(getPenetranceValue(LocGenotype((iit1->first) % numLoci, geno)));
				}
				++iit0;
				++iit1;
			}
		}
	} else {
		ValueError("Operator PyMlSelector currently only supports haploid and diploid populations");
	}
#else
	if (ply == 1) {
		vectora geno(1, TO_ALLELE(0));
		if (m_loci.allAvail()) {
			GenoIterator it = ind->genoBegin();
			GenoIterator it_end = ind->genoEnd();
			for (size_t index = 0; it != it_end; ++it, ++index)
				if (*it != 0) {
					geno[0] = *it;
					pnt.push(getPenetranceValue(LocGenotype(index, geno)));
				}
		} else {
			GenoIterator it = ind->genoBegin();
			vectoru::const_iterator loc = loci.begin();
			vectoru::const_iterator locEnd = loci.end();
			for (; loc != locEnd; ++loc)
				if (*(it + *loc) != 0) {
					geno[0] = *(it + *loc);
					pnt.push(getPenetranceValue(LocGenotype((*loc), geno)));
				}
		}
	} else if (ply == 2) {
		vectora geno(2, TO_ALLELE(0));
		if (m_loci.allAvail()) {
			GenoIterator it0 = ind->genoBegin(0);
			GenoIterator it0_end = ind->genoEnd(0);
			GenoIterator it1 = ind->genoBegin(1);
			for (size_t index = 0; it0 != it0_end; ++it0, ++it1, ++index)
				if (*it0 != 0 || *it1 != 0) {
					geno[0] = *it0;
					geno[1] = *it1;
					pnt.push(getPenetranceValue(LocGenotype(index, geno)));
				}
		} else {
			GenoIterator it0 = ind->genoBegin(0);
			GenoIterator it1 = ind->genoBegin(1);
			vectoru::const_iterator loc = loci.begin();
			vectoru::const_iterator locEnd = loci.end();
			for (; loc != locEnd; ++loc)
				if (*(it0 + *loc) != 0 || *(it1 + *loc) != 0) {
					geno[0] = *(it0 + *loc);
					geno[1] = *(it1 + *loc);
					pnt.push(getPenetranceValue(LocGenotype((*loc), geno)));
				}
		}
	} else {
		ValueError("Operator PyMlSelector currently only supports haploid and diploid populations");
	}
#endif
	return pnt.value();
}


double PyMlPenetrance::getPenetranceValue(const LocGenotype & geno) const
{
	LocGenotype tmp(geno);

	if (geno.second.size() == 2 && TO_ALLELE(geno.second[0]) > TO_ALLELE(geno.second[1])) {
		tmp.second[0] = TO_ALLELE(geno.second[1]);
		tmp.second[1] = TO_ALLELE(geno.second[0]);
	}
	GenoPenetranceMap::iterator sit = m_penetFactory.find(tmp);
	if (sit != m_penetFactory.end())
		return sit->second;

	size_t nGeno = geno.second.size();
	double penet = 0;
	PyObject * res = NULL;

	if (m_searchMode == 10)
		res = m_func("()");
	else if (m_searchMode == 11)
		res = m_func("(i)", geno.first);
	else if (m_searchMode == 12) {
		if (nGeno == 1)
			res = m_func("((i))", TO_ALLELE(geno.second[0]));
		else
			res = m_func("((ii))", TO_ALLELE(geno.second[0]), TO_ALLELE(geno.second[1]));
	} else if (m_searchMode == 13) {
		if (nGeno == 1)
			res = m_func("((i)i)", TO_ALLELE(geno.second[0]), geno.first);
		else
			res = m_func("((ii)i)", TO_ALLELE(geno.second[0]), TO_ALLELE(geno.second[1]), geno.first);
	} else if (m_searchMode == 14) {
		if (nGeno == 1)
			res = m_func("(i(i))", geno.first, TO_ALLELE(geno.second[0]));
		else
			res = m_func("(i(ii))", geno.first, TO_ALLELE(geno.second[0]), TO_ALLELE(geno.second[1]));
	}


	penet = PyFloat_AsDouble(res);
	Py_DECREF(res);
	m_penetFactory[tmp] = penet;
	return penet;
}


}


