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


/** CPPONLY
 *  An accumulator class to generate overall fitness values from existing
 *  ones.
 */
class FitnessAccumulator
{
public:
	FitnessAccumulator(int mode) : m_mode(mode)
	{
		m_value = m_mode == EXPONENTIAL ? 0 : 1;
	}


	void push(double val)
	{
		if (m_mode == MULTIPLICATIVE)
			m_value *= val;
		else if (m_mode == ADDITIVE)
			m_value -= 1 - val;
		else if (m_mode == HETEROGENEITY)
			m_value *= 1 - val;
		else if (m_mode == EXPONENTIAL)
			m_value += 1 - val;
	}


	double value()
	{
		if (m_mode == MULTIPLICATIVE)
			return m_value;
		else if (m_mode == ADDITIVE)
			return m_value < 0 ? 0. : m_value;
		else if (m_mode == HETEROGENEITY)
			return m_value < 1 ? 1 - m_value : 0;
		else if (m_mode == EXPONENTIAL)
			return exp(-m_value);
		else
			throw ValueError("Unrecognized accumulation mode");
		return 0;
	}


private:
	int m_mode;

	double m_value;
};


double MlSelector::indFitness(Population & pop, Individual * ind) const
{
	FitnessAccumulator fit(m_mode);

	opList::const_iterator s = m_selectors.begin();
	opList::const_iterator sEnd = m_selectors.end();

	for (; s != sEnd; ++s)
		fit.push(dynamic_cast<const BaseSelector * >(*s)->indFitness(pop, ind));
	return fit.value();
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


RandomFitnessSelector::RandomFitnessSelector(const floatListFunc & selDist, int mode,
	const lociList & loci, const stringFunc & output, int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops, const stringList & infoFields) :
	m_selDist(selDist), m_mode(mode), m_loci(loci), m_searchMode(0),
	m_newGenotypes(), m_fitnessFactory()
{
	size_t sz = selDist.size();

	if (sz == 0) {
		// a function ?
		DBG_ASSERT(selDist.func().isValid(), ValueError,
			"A Python function or a list is expected for parameter selDist.");
		const pyFunc & func = m_selDist.func();
		if (func.numArgs() == 0)
			m_searchMode = 10;
		else if (func.numArgs() == 1) {
			if (func.arg(0) == "loc")
				m_searchMode = 11;
			else
				throw ValueError("Invalid parameter for passed function " + func.name() +
					" (allele1 and allele2 have to be specified together)");
		} else if (func.numArgs() == 2) {
			if (func.arg(0) == "allele1" and func.arg(1) == "allele2")
				m_searchMode = 12;
			else
				throw ValueError("Invalid parameters for passed function " + func.name() +
					" (allele1 and allele2 have to be specified together, in that order)");
		} else if (func.numArgs() == 3) {
			if (func.arg(0) == "loc" and func.arg(1) == "allele1" and func.arg(2) == "allele2")
				m_searchMode = 13;
			else if (func.arg(0) == "allele1" and func.arg(1) == "allele2" and func.arg(2) == "loc")
				m_searchMode = 14;
			else
				throw ValueError("Invalid parameters for passed function " + func.name());
		} else {
			DBG_FAILIF(true, ValueError,
				"Python function for selDist only accepts paramters loc, allele, allele1, or allele2");
		}
	} else {
		int mode = static_cast<int>(m_selDist[0]);
		if (mode == CONSTANT) {
			if (m_selDist.size() == 2)
				m_searchMode = 1;
			else if (m_selDist.size() == 3)
				m_searchMode = 2;
			else {
				DBG_FAILIF(true, ValueError,
					"CONSTANT mode can only be followed by one or two parameters");
			}
		} else if (mode == GAMMA_DISTRIBUTION) {
			if (m_selDist.size() == 3)
				m_searchMode = 3;
			else if (m_selDist.size() == 4)
				m_searchMode = 4;
			else {
				DBG_FAILIF(true, ValueError,
					"GAMMA_DISTRIBUTION mode can only be followed by two or three parameters");
			}
		} else {
			DBG_FAILIF(true, ValueError,
				"selDist only accept CONSTANT or GAMMA_DISTRIBUTION modes");
		}
	}
}


double RandomFitnessSelector::indFitness(Population & /* pop */, Individual * ind) const
{
	FitnessAccumulator fit(m_mode);
	const vectoru & loci = m_loci.elems(ind);

#ifdef MUTANTALLELE
	size_t numLoci = ind->totNumLoci();
	GenoIterator it0 = ind->genoBegin(0);
	GenoIterator it0_end = ind->genoEnd(0);
	GenoIterator it1 = ind->genoBegin(1);
	GenoIterator it1_end = ind->genoEnd(1);
	compressed_vector<size_t>::index_array_type::iterator index_it0 = it0.getIndexIterator();
	compressed_vector<size_t>::index_array_type::iterator index_it0_end = it0_end.getIndexIterator();
	compressed_vector<Allele>::value_array_type::iterator value_it0 = it0.getValueIterator();
	compressed_vector<size_t>::index_array_type::iterator index_it1 = it1.getIndexIterator();
	compressed_vector<size_t>::index_array_type::iterator index_it1_end = it1_end.getIndexIterator();
	compressed_vector<Allele>::value_array_type::iterator value_it1 = it1.getValueIterator();
	for (; index_it0 != index_it0_end || index_it1 != index_it1_end; ) {
		if (index_it1 == index_it1_end || *index_it0 + numLoci < *index_it1) {
			if (*value_it0 != 0 &&
			    (m_loci.allAvail() || find(loci.begin(), loci.end(), (*index_it0) % numLoci) != loci.end()))
				fit.push(getGenotypeFitnessValue(LocGenotype((*index_it0) % numLoci, std::pair<Allele, Allele>(*value_it0, ToAllele(0)))));
			++index_it0;
			++value_it0;
		} else if (index_it0 == index_it0_end || *index_it0 + numLoci > *index_it1) {
			if (*value_it1 != 0 &&
			    (m_loci.allAvail() || find(loci.begin(), loci.end(), (*index_it1) % numLoci) != loci.end()))
				fit.push(getGenotypeFitnessValue(LocGenotype((*index_it1) % numLoci, std::pair<Allele, Allele>(ToAllele(0), *value_it1))));
			++index_it1;
			++value_it1;
		} else {
			if ((*value_it0 != 0 || *value_it1 != 0) &&
			    (m_loci.allAvail() || find(loci.begin(), loci.end(), (*index_it1) % numLoci) != loci.end()))
				fit.push(getGenotypeFitnessValue(LocGenotype((*index_it1) % numLoci, std::pair<Allele, Allele>(*value_it0, *value_it1))));
			++index_it0;
			++value_it0;
			++index_it1;
			++value_it1;
		}
	}
#else
	if (m_loci.allAvail()) {
		GenoIterator it0 = ind->genoBegin(0);
		GenoIterator it0_end = ind->genoEnd(0);
		GenoIterator it1 = ind->genoBegin(1);
		for (size_t index = 0; it0 != it0_end; ++it0, ++it1, ++index)
			if (*it0 != 0 || *it1 != 0)
				fit.push(getGenotypeFitnessValue(LocGenotype(index, std::pair<Allele, Allele>(*it0, *it1))));
	} else {
		GenoIterator it0 = ind->genoBegin(0);
		GenoIterator it1 = ind->genoBegin(1);
		vectoru::const_iterator loc = loci.begin();
		vectoru::const_iterator locEnd = loci.end();
		for (; loc != locEnd; ++loc)
			if (*(it0 + *loc) != 0 || *(it1 + *loc) != 0)
				fit.push(getGenotypeFitnessValue(LocGenotype((*loc), std::pair<Allele, Allele>(*(it0 + *loc), *(it1 + *loc)))));
	}
#endif
	return fit.value();
}


bool RandomFitnessSelector::apply(Population & pop) const
{
	DBG_ASSERT(pop.ploidy() == 2, ValueError, "Operator RandomFitnessSelector only works with diploid populations");

	m_newGenotypes.clear();
	if (!BaseSelector::apply(pop))
		return false;
	// output NEW mutant...
	if (!m_newGenotypes.empty() && !noOutput()) {
		ostream & out = getOstream(pop.dict());
		vector<LocGenotype>::const_iterator it = m_newGenotypes.begin();
		vector<LocGenotype>::const_iterator it_end = m_newGenotypes.end();
		for (; it != it_end; ++it) {
			double s = m_fitnessFactory[*it];
			out << it->first << '\t' << it->second.first << '\t' << it->second.second << s << '\n';
		}
		closeOstream();
	}
	return true;
}


double RandomFitnessSelector::getGenotypeFitnessValue(const LocGenotype & geno) const
{
	GenoSelMap::iterator sit = m_fitnessFactory.find(geno);
	if (sit != m_fitnessFactory.end())
		return sit->second;
	
	double fitness = 0;
	if (m_searchMode == 1 or m_searchMode == 2) {
		if (geno.second.first != 0 && geno.second.second != 0)
			fitness = 1 - m_selDist[1];
		else
			fitness = 1 - m_selDist[1] * (m_searchMode == 1 ? 0.5 : m_selDist[2]);
	} else if (m_searchMode == 3 or m_searchMode == 4) {
		double s = getRNG().randGamma(m_selDist[1], m_selDist[2]);
		if (geno.second.first != 0 && geno.second.second != 0)
			fitness = 1 - s;
		else
			fitness = 1 - s * (m_searchMode == 3 ? 0.5 : m_selDist[3]);
	} else {
		const pyFunc & func = m_selDist.func();
		PyObject * res = NULL;
		if (m_searchMode == 10)
			res = func("()");
		else if (m_searchMode == 11)
			res = func("(i)", geno.first);
		else if (m_searchMode == 12)
			res = func("(ii)", geno.second.first, geno.second.second);
		else if (m_searchMode == 13)
			res = func("(iii)", geno.first, geno.second.first, geno.second.second);
		else if (m_searchMode == 14)
			res = func("(iii)", geno.second.first, geno.second.second, geno.first);
		fitness = 1 - PyFloat_AsDouble(res);
		Py_DECREF(res);
	}
	m_fitnessFactory[geno] = fitness;
	if (!noOutput())
		m_newGenotypes.push_back(geno);
	return fitness;
}


}


