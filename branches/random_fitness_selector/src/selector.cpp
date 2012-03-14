/**
 *  $File: selector.cpp $
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

#include "selector.h"

#if PY_VERSION_HEX >= 0x03000000
#  define PyInt_FromLong(x) PyLong_FromLong(x)
#endif

namespace simuPOP {
bool BaseSelector::apply(Population & pop) const
{
	size_t fit_id = pop.infoIdx(this->infoField(0));

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
				for (; ind.valid(); ++ind)
					ind->setInfo(indFitness(pop, &*ind), fit_id);
#endif
			}

		} else {
			IndIterator ind = pop.indIterator(sp->subPop());
			for (; ind.valid(); ++ind)
				ind->setInfo(indFitness(pop, &*ind), fit_id);
		}
		if (sp->isVirtual())
			pop.deactivateVirtualSubPop(sp->subPop());
	}

	return true;
}


double MapSelector::indFitness(Population & /* pop */, Individual * ind) const
{
	vectoru chromTypes;
	const vectoru & loci = m_loci.elems(ind);

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
	string allele_string = "(";
	for (size_t i = 0; i < alleles.size(); ++i) {
		if (i != 0)
			allele_string += ", ";
		allele_string += toStr(alleles[i]);
	}
	allele_string += ")";
	throw ValueError("No fitness value for genotype " + allele_string);
	// this line should not be reached.
	return 0;
}


// currently assuming diploid
double MaSelector::indFitness(Population & /* pop */, Individual * ind) const
{
	UINT index = 0;
	bool singleST = m_wildtype.size() == 1;
	const vectoru & loci = m_loci.elems(ind);

	DBG_FAILIF((ind->ploidy() == 2 && m_fitness.size() != static_cast<UINT>(pow(3., static_cast<double>(loci.size())))) ||
		(ind->ploidy() == 1 && m_fitness.size() != static_cast<UINT>(pow(2., static_cast<double>(loci.size())))),
		ValueError, "Please specify fitness for each combination of genotype.");

	for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc) {
		if (ind->ploidy() == 1) {
			Allele a = ToAllele(ind->allele(*loc));
			if (singleST)
				index = index * 2 + (AlleleUnsigned(a) != m_wildtype[0]);
			else
				index = index * 2 + (find(m_wildtype.begin(), m_wildtype.end(), AlleleUnsigned(a)) == m_wildtype.end());
		} else if (ind->ploidy() == 2) {
			// get genotype of ind
			Allele a = ToAllele(ind->allele(*loc, 0));
			Allele b = ToAllele(ind->allele(*loc, 1));

			int numWildtype = 0;

			// count number of wildtype
			// this improve the performance a little bit
			if (singleST) {
				numWildtype = (AlleleUnsigned(a) == m_wildtype[0]) + (AlleleUnsigned(b) == m_wildtype[0]);
			} else {
				if (find(m_wildtype.begin(), m_wildtype.end(), AlleleUnsigned(a)) != m_wildtype.end())
					numWildtype++;

				if (find(m_wildtype.begin(), m_wildtype.end(), AlleleUnsigned(b)) != m_wildtype.end())
					numWildtype++;
			}

			index = index * 3 + 2 - numWildtype;
		} else {
			DBG_FAILIF(true, ValueError, "The MaSelector only supports haploid and diploid populations.");
		}
	}
	return m_fitness[index];
}


double MlSelector::indFitness(Population & pop, Individual * ind) const
{
	if (m_mode == MULTIPLICATIVE) {
		double fit = 1;
		for (opList::const_iterator s = m_selectors.begin(), sEnd = m_selectors.end();
		     s != sEnd; ++s)
			fit *= dynamic_cast<const BaseSelector * >(*s)->indFitness(pop, ind);
		return fit;
	} else if (m_mode == ADDITIVE) {
		double fit = 1;
		for (opList::const_iterator s = m_selectors.begin(), sEnd = m_selectors.end();
		     s != sEnd; ++s)
			fit -= 1 - dynamic_cast<const BaseSelector * >(*s)->indFitness(pop, ind);
		return fit < 0 ? 0. : fit;
	} else if (m_mode == HETEROGENEITY) {
		double fit = 1;
		for (opList::const_iterator s = m_selectors.begin(), sEnd = m_selectors.end();
		     s != sEnd; ++s)
			fit *= 1 - dynamic_cast<const BaseSelector * >(*s)->indFitness(pop, ind);
		return fit < 1 ? 1 - fit : 0;
	} else if (m_mode == EXPONENTIAL) {
		double fit = 0;
		for (opList::const_iterator s = m_selectors.begin(), sEnd = m_selectors.end();
		     s != sEnd; ++s)
			fit += 1 - dynamic_cast<const BaseSelector * >(*s)->indFitness(pop, ind);
		return exp(-fit);
	}
	// this is the case for none.
	return 1.0;
}


double PySelector::indFitness(Population & pop, Individual * ind) const
{
	PyObject * args = PyTuple_New(m_func.numArgs());

	DBG_ASSERT(args, RuntimeError, "Failed to create a parameter tuple");

	for (size_t i = 0; i < m_func.numArgs(); ++i) {
		const string & arg = m_func.arg(i);
		if (arg == "ind")
			PyTuple_SET_ITEM(args, i, pyIndObj(static_cast<void *>(ind)));
		else if (arg == "geno")
			PyTuple_SET_ITEM(args, i, ind->genoAtLoci(m_loci));
		else if (arg == "gen")
			PyTuple_SET_ITEM(args, i, PyInt_FromLong(static_cast<long>(pop.gen())));
		else if (arg == "pop")
			PyTuple_SET_ITEM(args, i, pyPopObj(static_cast<void *>(&pop)));
		else {
			DBG_FAILIF(!ind->hasInfoField(arg), ValueError,
				"Only parameters 'ind', 'geno', 'gen', 'pop' and names of information fields are "
				"acceptable in function " + m_func.name());
			PyTuple_SET_ITEM(args, i, PyFloat_FromDouble(ind->info(arg)));
		}
	}

	double fitness = m_func(PyObj_As_Double, args);
	Py_XDECREF(args);
	return fitness;
}


double RandomFitnessSelector::indFitness(Population & /* pop */, Individual * ind) const
{
	if (m_mode == MULTIPLICATIVE) {
		return randomSelMulFitnessExt(ind->genoBegin(), ind->genoEnd());
	} else if (m_mode == ADDITIVE) {
		if (m_additive)
			return randomSelAddFitness(ind->genoBegin(), ind->genoEnd());
		else
			return randomSelAddFitnessExt(ind->genoBegin(), ind->genoEnd());
	} else if (m_mode == EXPONENTIAL) {
		if (m_additive)
			return randomSelExpFitness(ind->genoBegin(), ind->genoEnd());
		else
			return randomSelExpFitnessExt(ind->genoBegin(), ind->genoEnd());
	}
	return 0;
}


bool RandomFitnessSelector::apply(Population & pop) const
{
	DBG_ASSERT(pop.ploidy() != 2, ValueError, "Operator RandomFitnessSelector is applicable "
		"only to diploid populations.");
	m_numLoci = pop.totNumLoci();
	m_newMutants.clear();
	if (!BaseSelector::apply(pop))
		return false;
	// output NEW mutant...
	if (!m_newMutants.empty() && !noOutput()) {
		ostream & out = getOstream(pop.dict());
		vectoru::const_iterator it = m_newMutants.begin();
		vectoru::const_iterator it_end = m_newMutants.end();
		for (; it != it_end; ++it) {
			SelCoef s = m_selFactory[*it];
			out << it->first << '\t' << it->second << '\t' << s.first << '\t' << s.second << '\n';
		}
		closeOstream();
	}
	return true;
}


RandomFitnessSelector::SelCoef RandomFitnessSelector::getFitnessValue(const Mutant & mutant) const
{
	size_t sz = m_selDist.size();
	double s = 0;
	double h = 0.5;

	if (sz == 0) {
		// call a function
		const pyFunc & func = m_selDist.func();
		PyObject * res;
		if (func.numArgs() == 0)
			res = func("()");
		else {
			DBG_FAILIF(func.arg(0) != "mut", ValueError,
				"Only parameter loc is accepted for this user-defined function.");
			res = func("(iii)", mutant[0], mutant[1], mutant[2]);
		}
		if (PyNumber_Check(res)) {
			s = PyFloat_AsDouble(res);
		} else if (PySequence_Check(res)) {
			size_t sz = PySequence_Size(res);
			DBG_FAILIF(sz == 0, RuntimeError, "Function return an empty list.");
			PyObject * item = PySequence_GetItem(res, 0);
			s = PyFloat_AsDouble(item);
			Py_DECREF(item);
			if (sz > 1) {
				item = PySequence_GetItem(res, 1);
				h = PyFloat_AsDouble(item);
				Py_DECREF(item);
			}
		}
		Py_DECREF(res);
	} else {
		int mode = static_cast<int>(m_selDist[0]);
		if (mode == CONSTANT) {
			// constant
			s = m_selDist[1];
			if (m_selDist.size() > 2)
				h = m_selDist[2];
		} else {
			// a gamma distribution
			s = getRNG().randGamma(m_selDist[1], m_selDist[2]);
			if (m_selDist.size() > 3)
				h = m_selDist[3];
		}
	}
	m_selFactory[mutant] = SelCoef(s, h);
	m_newMutants.push_back(mutant);
	if (m_additive && h != 0.5)
		m_additive = false;
	return SelCoef(s, h);
}


double RandomFitnessSelector::randomSelAddFitness(GenoIterator it, GenoIterator it_end) const
{
#ifdef MUTANTALLELE
	double s = 0;
	compressed_vector<size_t>::index_array_type::iterator index_it = it.getIndexIterator();
	compressed_vector<size_t>::index_value_type::iterator value_it = it.getValueIterator();
	compressed_vector<size_t>::index_array_type::iterator index_it_end = it_end.getIndexIterator();
	for (; index_it != index_it_end; ++index_it, ++value_it) {
		if (*value_it == 0)
			continue;
		Mutant mut(*index_it, *value_it);
		SelMap::iterator sit = m_selFactory.find(mut);
		if (sit == m_selFactory.end())
			s += getFitnessValue(mut).first / 2.;
		else
			s += sit->second.first / 2;
	}
	return 1 - s > 0 ? 1 - s : 0;
#else
	size_t idx = 0;
	for (; it != it_end; ++it, ++idx) {
		if (*it != 0) {
			Mutant mut(idx, *it);
			SelMap::iterator sit = m_selFactory.find(mut);
			if (sit == m_selFactory.end())
				s += getFitnessValue(mut).first / 2.;
			else
				s += sit->second.first / 2;
		}
	}
	return 1 - s > 0 ? 1 - s : 0;
#endif
}


double RandomFitnessSelector::randomSelExpFitness(GenoIterator it, GenoIterator it_end) const
{
#ifdef MUTANTALLELE
	double s = 0;

	compressed_vector<size_t>::index_array_type::iterator index_it = it.getIndexIterator();
	compressed_vector<size_t>::index_value_type::iterator value_it = it.getValueIterator();
	compressed_vector<size_t>::index_array_type::iterator index_it_end = it_end.getIndexIterator();
	for (; index_it != index_it_end; ++index_it, ++value_it) {
		if (*index_it == 0)
			continue;
		Mutant mut(*index_it, *value_it);
		SelMap::iterator sit = m_selFactory.find(mut);
		if (sit == m_selFactory.end())
			s += getFitnessValue(mut).first / 2.;
		else
			s += sit->second.first / 2;
	}
	return exp(-s);
#else
	size_t idx = 0;
	for (; it != it_end; ++it, ++idx) {
		if (*it != 0) {
			Mutant mut(idx, *it);
			SelMap::iterator sit = m_selFactory.find(mut);
			if (sit == m_selFactory.end())
				s += getFitnessValue(mut).first / 2.;
			else
				s += sit->second.first / 2;
		}
	}
	return exp(-s);
#endif
}


double RandomFitnessSelector::randomSelMulFitnessExt(GenoIterator it, GenoIterator it_end) const
{
#ifdef MUTANTALLELE
	MutCounter cnt;
	compressed_vector<size_t>::index_array_type::iterator index_it = it.getIndexIterator();
	compressed_vector<size_t>::index_value_type::iterator value_it = it.getValueIterator();
	compressed_vector<size_t>::index_array_type::iterator index_it_end = it_end.getIndexIterator();
	for (; index_it != index_it_end; ++index_it, ++value_it) {
		if (*value_it == 0)
			continue;
		Mutant mut(*index_it, *value_it);
		MutCounter::iterator mit = cnt.find(mut);
		if (mit == cnt.end())
			cnt[mut] = 1;
		else
			++mit->second;
	}

	double s = 1;
	MutCounter::iterator mit = cnt.begin();
	MutCounter::iterator mit_end = cnt.end();
	for (; mit != mit_end; ++mit) {
		SelMap::iterator sit = m_selFactory.find(mit->first);
		if (sit == m_selFactory.end()) {
			SelCoef sf = getFitnessValue(mit->first);
			if (mit->second == 1)
				s *= 1 - sf.first * sf.second;
			else
				s *= 1 - sf.first;
		} else {
			if (mit->second == 1)
				s *= 1 - sit->second.first * sit->second.second;
			else
				s *= 1 - sit->second.first;
		}
	}
	return s;
#else
	size_t idx = 0;
	for (; it != it_end; ++it, ++idx) {
		if (*it != 0) {
			Mutant mut(idx, *it);
			SelMap::iterator sit = m_selFactory.find(mut);
			if (sit == m_selFactory.end())
				s += getFitnessValue(mut).first / 2.;
			else
				s += sit->second.first / 2;
		}
	}
	return exp(-s);
#endif
}


double RandomFitnessSelector::randomSelAddFitnessExt(GenoIterator it, GenoIterator it_end) const
{
#ifdef MUTANTALLELE
	MutCounter cnt;

	compressed_vector<size_t>::index_array_type::iterator index_it = it.getIndexIterator();
	compressed_vector<size_t>::index_array_type::iterator index_it_end = it_end.getIndexIterator();
	for (; index_it != index_it_end; ++index_it) {
		if (*index_it == 0)
			continue;
		MutCounter::iterator mit = cnt.find(*index_it);
		if (mit == cnt.end())
			cnt[*index_it] = 1;
		else
			++mit->second;
	}

	double s = 0;
	MutCounter::iterator mit = cnt.begin();
	MutCounter::iterator mit_end = cnt.end();
	for (; mit != mit_end; ++mit) {
		SelMap::iterator sit = m_selFactory.find(mit->first);
		if (sit == m_selFactory.end()) {
			SelCoef sf = getFitnessValue(mit->first);
			if (mit->second == 1)
				s += sf.first * sf.second;
			else
				s += sf.first;
		} else {
			if (mit->second == 1)
				s += sit->second.first * sit->second.second;
			else
				s += sit->second.first;
		}
	}
	return 1 - s > 0 ? 1 - s : 0;
#else
	(void)it;
	(void)it_end;
	return 0;
#endif
}


double RandomFitnessSelector::randomSelExpFitnessExt(GenoIterator it, GenoIterator it_end) const
{
#ifdef MUTANTALLELE
	MutCounter cnt;

	compressed_vector<size_t>::index_array_type::iterator index_it = it.getIndexIterator();
	compressed_vector<size_t>::index_array_type::iterator index_it_end = it_end.getIndexIterator();
	for (; index_it != index_it_end; ++index_it) {
		if (*index_it == 0)
			continue;
		MutCounter::iterator mit = cnt.find(*index_it);
		if (mit == cnt.end())
			cnt[*index_it] = 1;
		else
			++mit->second;
	}

	double s = 0;
	MutCounter::iterator mit = cnt.begin();
	MutCounter::iterator mit_end = cnt.end();
	for (; mit != mit_end; ++mit) {
		SelMap::iterator sit = m_selFactory.find(mit->first);
		if (sit == m_selFactory.end()) {
			SelCoef sf = getFitnessValue(mit->first);
			if (mit->second == 1)
				s += sf.first * sf.second;
			else
				s += sf.first;
		} else {
			if (mit->second == 1)
				s += sit->second.first * sit->second.second;
			else
				s += sit->second.first;
		}
	}
	return exp(-s);
#else
	(void)it;
	(void)it_end;
	return 0;
#endif
}


}


