/**
 *  $File: initializer.cpp $
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

#include "initializer.h"

namespace simuPOP {

string InitSex::describe(bool /* format */) const
{
	string desc = "<simuPOP.InitSex> initialize sex ";

	if (!m_sex.empty())
		desc += "using a list";
	else if (m_maleProp < 0) {
		desc += "randomly";
		if (m_maleFreq != 0.5)
			desc += (boost::format(" with probability %1% being a male") % m_maleFreq).str();
	} else
		desc += (boost::format("randomly with %1% percent of males") % (m_maleProp * 100)).str();
	return desc;
}


bool InitSex::apply(Population & pop) const
{
	const subPopList subPops = applicableSubPops(pop);

	size_t idx = 0;
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator sp_end = subPops.end();

	for (; sp != sp_end; ++sp) {
		WeightedSampler ws;
		if (m_maleProp >= 0) {
			vectorf prop(2, m_maleProp);
			prop[1] = 1 - prop[0];
			ws.set(prop.begin(), prop.end(), pop.subPopSize(*sp));
		} else if (m_sex.empty()) {
			vectorf prop(2, m_maleFreq);
			prop[1] = 1 - prop[0];
			ws.set(prop.begin(), prop.end());
		}
		pop.activateVirtualSubPop(*sp);
		IndIterator ind = pop.indIterator(sp->subPop());
		size_t sexSz = m_sex.size();
		if (!m_sex.empty())
			for (; ind.valid(); ++ind, ++idx)
				ind->setSex(m_sex[idx % sexSz] == 1 ? MALE : FEMALE);
		else {
			if (numThreads() > 1) {
#ifdef _OPENMP
#  pragma omp parallel private(ind)
				{
					ind = pop.indIterator(sp->subPop(), omp_get_thread_num());
					for (; ind.valid(); ++ind)
						ind->setSex(ws.draw() == 0 ? MALE : FEMALE);
				}
#endif
			} else
				for (; ind.valid(); ++ind)
					ind->setSex(ws.draw() == 0 ? MALE : FEMALE);
		}
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


string InitInfo::describe(bool /* format */) const
{
	string desc = "<simuPOP.InitInfo> initialize information field";

	if (infoSize() > 1)
		desc += "s";
	for (size_t i = 0; i < infoSize(); ++i) {
		desc += i == 0 ? " " : ", ";
		desc += infoField(i);
	}
	if (m_values.empty())
		desc += " using a Python function " + m_values.func().name();
	else
		desc += " using a list of values";
	return desc;
}


bool InitInfo::apply(Population & pop) const
{
	vectoru infoIdx(infoSize());

	if (infoIdx.empty())
		return true;

	for (size_t i = 0; i < infoIdx.size(); ++i)
		infoIdx[i] = pop.infoIdx(infoField(i));

	const subPopList subPops = applicableSubPops(pop);
	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator sp_end = subPops.end();
	size_t idx = 0;
	const vectorf & values = m_values.elems();

	for (; sp != sp_end; ++sp) {
		pop.activateVirtualSubPop(*sp);
		size_t numValues = m_values.size();
		if (numThreads() > 1 && !values.empty()) {
#ifdef _OPENMP
#  pragma omp parallel firstprivate (idx)
			{
				size_t id = omp_get_thread_num();
				IndIterator ind = pop.indIterator(sp->subPop(), id);
				idx = idx + id * (pop.subPopSize(sp->subPop()) / numThreads());
				for (; ind.valid(); ++ind, ++idx)
					for (size_t i = 0; i < infoIdx.size(); ++i) {
						ind->setInfo(values[idx % numValues], infoIdx[i]);
					}
			}
			idx = idx + pop.subPopSize(sp->subPop());
#endif
		} else {
			IndIterator ind = pop.indIterator(sp->subPop());
			for (; ind.valid(); ++ind, ++idx) {
				for (size_t i = 0; i < infoIdx.size(); ++i) {
					if (values.empty())
						ind->setInfo(m_values.func() (PyObj_As_Double, "()"), infoIdx[i]);
					else
						ind->setInfo(values[idx % numValues], infoIdx[i]);
				}
			}
		}
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


InitGenotype::InitGenotype(const floatListFunc & freq,
	const uintList & genotype, const floatListFunc & prop,
	const intMatrix & haplotypes,
	const intMatrix & genotypes,
	const lociList & loci,
	const uintList & ploidy,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	const stringList & infoFields)
	: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_freq(freq), m_genotype(genotype.elems()), m_prop(prop),
	m_haplotypes(haplotypes.elems()), m_genotypes(genotypes.elems()),
	m_loci(loci), m_ploidy(ploidy)
{
	for (size_t i = 0; i < m_freq.elems().size(); ++i)
		PARAM_FAILIF(fcmp_lt(m_freq.elems()[i], 0) || fcmp_gt(m_freq.elems()[i], 1), ValueError,
			"Allele frequency should be between 0 and 1");
	for (size_t i = 0; i < m_prop.elems().size(); ++i)
		PARAM_FAILIF(fcmp_lt(m_prop.elems()[i], 0) || fcmp_gt(m_prop.elems()[i], 1), ValueError,
			"Allele proportion should be between 0 and 1");
	PARAM_FAILIF(!m_haplotypes.empty() && !m_freq.elems().empty() && m_haplotypes.size() != m_freq.size(),
		ValueError, "Haplotype frequency, if specified, should be specified for each haplotype.")
	PARAM_FAILIF(!m_haplotypes.empty() && !m_prop.elems().empty() && m_haplotypes.size() != m_prop.size(),
		ValueError, "Haplotype proportion, if specified, should be specified for each haplotype.")
	PARAM_FAILIF(!m_genotypes.empty() && !m_freq.elems().empty() && m_genotypes.size() != m_freq.size(),
		ValueError, "Genotype frequency, if specified, should be specified for each genotype.")
	PARAM_FAILIF(!m_genotypes.empty() && !m_prop.elems().empty() && m_genotypes.size() != m_prop.size(),
		ValueError, "Genotype proportion, if specified, should be specified for each genotype.")
	PARAM_FAILIF(!m_haplotypes.empty() && !m_genotypes.empty(),
		ValueError, "Please specify only one of parameters genotypes and haplotypes")
	PARAM_FAILIF(!m_freq.elems().empty() && fcmp_ne(accumulate(m_freq.elems().begin(), m_freq.elems().end(), 0.), 1.0), ValueError,
		"Allele frequencies should add up to one.");
	PARAM_FAILIF(!m_prop.elems().empty() && fcmp_ne(accumulate(m_prop.elems().begin(), m_prop.elems().end(), 0.), 1.0), ValueError,
		"Allele proportions should add up to one.");
	PARAM_ASSERT(!m_haplotypes.empty() || !m_genotypes.empty() || (m_freq.valid() + !m_genotype.empty() + m_prop.valid() == 1), ValueError,
		"Please specify one and only one of parameters freq, prop (or haplotypes) and genotype.");
}


string InitGenotype::describe(bool /* format */) const
{
	string desc = "<simuPOP.InitGenotype> initialize individual genotype ";

	if (m_freq.valid())
		desc += string("acccording to ") + (m_haplotypes.empty() ? " allele " : " haplotype ") + "frequencies.";
	else if (m_prop.valid())
		desc += string("according to proportion of ") + (m_haplotypes.empty() ? (m_genotypes.empty() ? " alleles." : "genotypes.") : " haplotypes.");
	else
		desc += "using specified genotype.";
	return desc;
}


#if PY_VERSION_HEX >= 0x03000000
#  define PyInt_FromLong(x) PyLong_FromLong(x)
#endif

vectorf InitGenotype::getFreqOrProp(size_t locus, const vspID & vsp) const
{
	pyFunc func = m_freq.func().isValid() ? m_freq.func() : m_prop.func();

	PyObject * args = PyTuple_New(func.numArgs());
	for (size_t i = 0; i < func.numArgs(); ++i) {
		const string & arg = func.arg(i);
		if (arg == "loc")
			PyTuple_SET_ITEM(args, i, PyInt_FromLong(locus));
		else if (arg == "subPop")
			PyTuple_SET_ITEM(args, i, PyInt_FromLong(vsp.subPop()));
		else if (arg == "vsp") {
			PyObject * v = PyTuple_New(2);
			PyTuple_SET_ITEM(v, 0, PyInt_FromLong(vsp.subPop()));
			if (vsp.isVirtual())
				PyTuple_SET_ITEM(v, 1, PyInt_FromLong(vsp.virtualSubPop()));
			else {
				Py_INCREF(Py_None);
				PyTuple_SET_ITEM(v, 1, Py_None);
			}
			// it steals reference count of v
			PyTuple_SET_ITEM(args, i, v);
		}
	}
	vectorf freq = func(PyObj_As_Array, args);
	if (freq.size() == 1 && freq.back() != 1.0)
		freq.push_back(1.0 - freq[0]);
	for (size_t i = 0; i < freq.size(); ++i) {
		if (fcmp_lt(freq[i], 0) || fcmp_gt(freq[i], 1))
			throw RuntimeError("Returned allele frequency should be between 0 and 1.");
	}
	if (fcmp_ne(std::accumulate(freq.begin(), freq.end(), 0.0), 1.0))
		throw RuntimeError("Returned allele frequency should add up to 1.");
	Py_XDECREF(args);
	return freq;
}

bool InitGenotype::apply(Population & pop) const
{
	const subPopList subPops = applicableSubPops(pop);

	const vectoru & loci = m_loci.elems(&pop);

	for (size_t i = 0; i < m_haplotypes.size(); ++i) {
		DBG_WARNIF(m_haplotypes[i].size() != loci.size(),
			(boost::format("Haplotype [%1%] specified in operator InitGenotype has %2% alleles but %3%"
				           " alleles are needed. This haplotype will be truncated or repeated.") % m_haplotypes[i] % m_haplotypes[i].size() % loci.size()).str());
	}
	for (size_t i = 0; i < m_genotypes.size(); ++i) {
		DBG_FAILIF(m_genotypes[i].size() != pop.ploidy(), ValueError,
			(boost::format("Genotype [%1%] specified in operator InitGenotype should have size %2% (ploidy)")
				% m_genotypes[i] % pop.ploidy()).str());
	}

	vectoru ploidy = m_ploidy.elems();
	if (m_ploidy.allAvail())
		for (size_t i = 0 ; i < pop.ploidy(); ++i)
			ploidy.push_back(i);

	pop.syncIndPointers();

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator sp_end = subPops.end();
	size_t sz = m_genotype.size();

	// check function
	int by_locus = false;
	if (m_freq.func().isValid() || m_prop.func().isValid()) {
		pyFunc func = m_freq.func().isValid() ? m_freq.func() : m_prop.func();
		//
		for (size_t fp = 0; fp < func.numArgs(); ++fp) {
			if (func.arg(fp) == "loc")
				by_locus = true;
			else if (func.arg(fp) != "subPop" && func.arg(fp) != "vsp")
				throw ValueError("Unacceptable parameter " + func.arg(fp));
		}
	}

	DBG_WARNIF(sz >= 1 && sz != loci.size() && sz != loci.size() * pop.ploidy(),
		(boost::format("Genotype [%1%] specified in operator InitGenotype has %2% alleles, which does not match"
			           " number of loci of individuals. This sequence will be truncated or repeated and may lead to erroneous results.")
		 % m_genotype % m_genotype.size()).str());
	for (size_t idx = 0; sp != sp_end; ++sp) {
		// will go through virtual subpopulation if sp is virtual
		pop.activateVirtualSubPop(*sp);
		if (!m_genotype.empty()) {
			// multi-thread write to compressed mutant is not allowed.
			IndIterator it = pop.indIterator(sp->subPop());
			for (size_t ii = 0; it.valid(); ++it, ++ii)
				for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
					for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc, ++idx)
						it->setAllele(TO_ALLELE(m_genotype[idx % sz]), *loc, static_cast<int>(*p));
		} else if (m_prop.valid()) {
			WeightedSampler ws;
			size_t sz = pop.subPopSize(*sp);
			// use haplotype
			if (!m_haplotypes.empty()) {
				IndIterator it = pop.indIterator(sp->subPop());
				if (m_prop.elems().empty()) {
					if (by_locus)
						throw ValueError("Return value of user provided function cannot have parameter loc when haplotypes are provided.");
					vectorf prop = getFreqOrProp(0, *sp);
					ws.set(prop.begin(), prop.end(), sz * ploidy.size());
				} else
					ws.set(m_prop.elems().begin(), m_prop.elems().end(), sz * ploidy.size());
				for (; it.valid(); ++it)
					for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p) {
						const vectori & haplotype = m_haplotypes[ws.draw()];
						size_t hapSz = haplotype.size();
						size_t j = 0;
						for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc, ++j)
							it->setAllele(TO_ALLELE(haplotype[j % hapSz]), *loc, static_cast<int>(*p));
					}
			} else if (!m_genotypes.empty()) {
				for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc) {
					if (m_prop.elems().empty()) {
						vectorf prop = getFreqOrProp(*loc, *sp);
						ws.set(prop.begin(), prop.end(), sz);
					} else
						ws.set(m_prop.elems().begin(), m_prop.elems().end(), sz);
					IndIterator it = pop.indIterator(sp->subPop());
					for (; it.valid(); ++it) {
						const vectori & genotype = m_genotypes[ws.draw()];
						for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
							it->setAllele(TO_ALLELE(genotype[*p]), *loc, static_cast<int>(*p));
					}
				}
			} else  {
				// initialize by allele. Gurantee proportion at each locus.
				for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc) {
					if (m_prop.elems().empty()) {
						vectorf prop = getFreqOrProp(*loc, *sp);
						ws.set(prop.begin(), prop.end(), sz * ploidy.size());
					} else
						ws.set(m_prop.elems().begin(), m_prop.elems().end(), sz * ploidy.size());
					IndIterator it = pop.indIterator(sp->subPop());
					for (; it.valid(); ++it)
						for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
							it->setAllele(TO_ALLELE(ws.draw()), *loc, static_cast<int>(*p));
				}
			}
		} else {
			if (!m_haplotypes.empty()) {
				// if user provides haplotypes
				WeightedSampler ws;
				IndIterator it = pop.indIterator(sp->subPop());
				// m_freq can be empty if ....
				if (!m_freq.valid()) {
					vectorf freq = vectorf(m_haplotypes.size(), 1. / m_haplotypes.size());
					ws.set(freq.begin(), freq.end());
				} else if (m_freq.elems().empty()) {
					if (by_locus)
						throw ValueError("Return value of user provided function cannot have parameter loc when haplotypes are provided.");
					vectorf freq = getFreqOrProp(0, *sp);
					ws.set(freq.begin(), freq.end());
				} else
					ws.set(m_freq.elems().begin(), m_freq.elems().end());
				//
				for (; it.valid(); ++it) {
					for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p) {
						const vectori & haplotype = m_haplotypes[ws.draw()];
						size_t hapSz = haplotype.size();
						size_t j = 0;
						for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc, ++j)
							it->setAllele(TO_ALLELE(haplotype[j % hapSz]), *loc, static_cast<int>(*p));
					}
				}
			} else if (!m_genotypes.empty()) {
				// if user provides haplotypes
				WeightedSampler ws;
				//
				for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc) {
					// m_freq can be empty if ....
					if (!m_freq.valid()) {
						vectorf freq = vectorf(m_genotypes.size(), 1. / m_genotypes.size());
						ws.set(freq.begin(), freq.end());
					} else if (m_freq.elems().empty()) {
						vectorf freq = getFreqOrProp(*loc, *sp);
						ws.set(freq.begin(), freq.end());
					} else
						ws.set(m_freq.elems().begin(), m_freq.elems().end());

					IndIterator it = pop.indIterator(sp->subPop());
					for (; it.valid(); ++it) {
						const vectori & genotype = m_genotypes[ws.draw()];
						for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
							it->setAllele(TO_ALLELE(genotype[*p]), *loc, static_cast<int>(*p));
					}
				}

			} else {
				// use allele frequency
				if (m_freq.func().isValid()) {
					for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc) {
						WeightedSampler ws(getFreqOrProp(*loc, *sp));
						//
						IndIterator it = pop.indIterator(sp->subPop());
						for (; it.valid(); ++it) {
							for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p) {
								it->setAllele(TO_ALLELE(ws.draw()), *loc, static_cast<int>(*p));
							}
						}
					}
				} else {
					// use specified frequencies
					IndIterator it = pop.indIterator(sp->subPop());
					WeightedSampler ws(m_freq.elems());
					for (; it.valid(); ++it) {
						for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p) {
							for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc)
								it->setAllele(TO_ALLELE(ws.draw()), *loc, static_cast<int>(*p));
						}
					}
				}
			}
		}
		pop.deactivateVirtualSubPop(sp->subPop());
	}
	return true;
}


InitLineage::InitLineage(const intList & lineage, int mode,
	const lociList & loci, const uintList & ploidy,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops,
	const stringList & infoFields)
	: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_lineage(lineage.elems()), m_loci(loci),
	m_ploidy(ploidy), m_mode(mode)
{
	PARAM_ASSERT(m_mode == PER_ALLELE || m_mode == PER_CHROMOSOME || m_mode == PER_PLOIDY ||
		m_mode == PER_INDIVIDUAL || m_mode == FROM_INFO || m_mode == FROM_INFO_SIGNED, ValueError,
		"Paramter mode of operator InitLineage can only be one of PER_ALLELE, PER_CHROMOSOME, "
		"PER_PLOIDY, PER_INDIVIDUAL, FROM_INFO, or FROM_INFO_SIGNED.");
	PARAM_FAILIF((m_mode == FROM_INFO || m_mode == FROM_INFO_SIGNED) && infoSize() != 1,
		ValueError, "An information field is needed for mode FROM_INFO or FROM_INFO_SIGNED");
	PARAM_FAILIF((m_mode == FROM_INFO || m_mode == FROM_INFO_SIGNED) && !m_lineage.empty(),
		ValueError, "A list of lineage values is specified for mode FROM_INFO or FROM_INFO_SIGNED");
}


string InitLineage::describe(bool /* format */) const
{
	string desc = "<simuPOP.InitLineage> initialize individual lineage ";

	if (!m_lineage.empty()) {
		desc += "using specified lineage at each ";
		if (m_mode == PER_ALLELE)
			desc += "loci";
		else if (m_mode == PER_CHROMOSOME)
			desc += "chromosome";
		else if (m_mode == PER_PLOIDY)
			desc += "ploidy";
		else
			desc += "individual";
	} else {
		desc += "using the field '" + infoField(0) + "'.";
	}
	return desc;
}


bool InitLineage::apply(Population & pop) const
{
#ifdef LINEAGE
	const subPopList subPops = applicableSubPops(pop);

	const vectoru & loci = m_loci.elems(&pop);

	vectoru ploidy = m_ploidy.elems();

	if (m_ploidy.allAvail())
		for (size_t i = 0 ; i < pop.ploidy(); ++i)
			ploidy.push_back(i);

	pop.syncIndPointers();

	subPopList::const_iterator sp = subPops.begin();
	subPopList::const_iterator sp_end = subPops.end();
	size_t sz = m_lineage.size();

	if (m_mode == PER_ALLELE) {
		DBG_WARNIF(sz >= 1 && sz != loci.size() && sz != loci.size() * pop.ploidy(),
			(boost::format("Lineage [%1%] specified in operator InitLineage has %2% "
				           " lineages, which does not match number of loci of individuals. "
				           "This sequence will be truncated or repeated and may lead to erroneous results.")
			 % m_lineage % m_lineage.size()).str());
	} else if (m_mode == PER_PLOIDY) {
		DBG_WARNIF(sz >= 1 && sz != pop.ploidy() && sz != pop.ploidy() * pop.popSize(),
			(boost::format("Lineage [%1%] specified in operator InitLineage has %2% lineages, "
				           "which does not match number of ploidy of individuals. "
				           "This sequence will be truncated or repeated and may lead to erroneous results.")
			 % m_lineage % m_lineage.size()).str());
	} else {
		DBG_WARNIF(sz > 1 && sz != pop.popSize(),
			(boost::format("Lineage [%1%] specified in operator InitLineage has %2%"
				           " lineages, which does not match number of individuals. "
				           "This sequence will be truncated or repeated and may lead to erroneous results.")
			 % m_lineage % m_lineage.size()).str());
	}
	size_t nCh = pop.numChrom();
	vectoru chromIndex;
	if (m_mode == PER_CHROMOSOME) {
		// find the chromosome index for each locus
		chromIndex.resize(loci.size());
		for (size_t i = 0; i < loci.size(); ++i)
			chromIndex[i] = pop.chromLocusPair(loci[i]).first;
	}
	for (size_t idx = 0; sp != sp_end; ++sp) {
		// will go through virtual subpopulation if sp is virtual
		pop.activateVirtualSubPop(*sp);
		if (m_mode != FROM_INFO && m_mode != FROM_INFO_SIGNED) {
#  pragma omp parallel firstprivate(idx) if(numThreads() > 1)
			{
#  ifdef _OPENMP
				size_t id = omp_get_thread_num();
				IndIterator it = pop.indIterator(sp->subPop(), id);
				if (m_mode == PER_ALLELE)
					idx += id * (pop.subPopSize(sp->subPop()) / numThreads()) *
					       (ploidy.end() - ploidy.begin()) * (loci.end() - loci.begin());
				else if (m_mode == PER_CHROMOSOME)
					idx += id * (pop.subPopSize(sp->subPop()) / numThreads()) *
					       (ploidy.end() - ploidy.begin()) * nCh;
				else if (m_mode == PER_PLOIDY)
					idx += id * (pop.subPopSize(sp->subPop()) / numThreads()) *
					       (ploidy.end() - ploidy.begin());
				else
					idx += id * (pop.subPopSize(sp->subPop()) / numThreads());
#  else
				IndIterator it = pop.indIterator(sp->subPop());
#  endif
				if (sz == 0) {
					for (; it.valid(); ++it) {
						if (m_mode == PER_ALLELE)
							for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
								for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc, ++idx)
									it->setAlleleLineage(idx + 1, *loc, static_cast<int>(*p));
						else if (m_mode == PER_CHROMOSOME)
							for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p, idx += nCh) {
								vectoru::const_iterator loc = loci.begin();
								vectoru::const_iterator locEnd = loci.end();
								vectoru::const_iterator chIdx = chromIndex.begin();
								for (; loc != locEnd; ++loc, ++chIdx)
									it->setAlleleLineage((idx + *chIdx + 1), *loc, static_cast<int>(*p));
							}
						else if (m_mode == PER_PLOIDY)
							for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p, ++idx)
								for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc)
									it->setAlleleLineage(idx + 1, *loc, static_cast<int>(*p));
						else {
							for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
								for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc)
									it->setAlleleLineage(idx + 1, *loc, static_cast<int>(*p));
							++idx;
						}
					}

				} else {
					for (; it.valid(); ++it) {
						if (m_mode == PER_ALLELE)
							for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
								for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc, ++idx)
									it->setAlleleLineage(m_lineage[idx % sz], *loc, static_cast<int>(*p));
						else if (m_mode == PER_CHROMOSOME)
							for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p, idx += nCh) {
								vectoru::const_iterator loc = loci.begin();
								vectoru::const_iterator locEnd = loci.end();
								vectoru::const_iterator chIdx = chromIndex.begin();
								for (; loc != locEnd; ++loc, ++chIdx)
									it->setAlleleLineage(m_lineage[(idx + *chIdx) % sz], *loc, static_cast<int>(*p));
							}
						else if (m_mode == PER_PLOIDY)
							for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p, ++idx)
								for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc)
									it->setAlleleLineage(m_lineage[idx % sz], *loc, static_cast<int>(*p));
						else {
							for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p)
								for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc)
									it->setAlleleLineage(m_lineage[idx % sz], *loc, static_cast<int>(*p));
							++idx;
						}
					}
				}
			}
#  ifdef _OPENMP
			if (m_mode == PER_ALLELE)
				idx += pop.subPopSize(sp->subPop()) *
				       (ploidy.end() - ploidy.begin()) * (loci.end() - loci.begin());
			else if (m_mode == PER_CHROMOSOME)
				idx += pop.subPopSize(sp->subPop()) * (ploidy.end() - ploidy.begin()) * nCh;
			else if (m_mode == PER_PLOIDY)
				idx += pop.subPopSize(sp->subPop()) * (ploidy.end() - ploidy.begin());
			else
				idx += pop.subPopSize(sp->subPop());
#  endif
		} else {
#  pragma omp parallel if(numThreads() > 1)
			{
#  ifdef _OPENMP
				IndIterator it = pop.indIterator(sp->subPop(), omp_get_thread_num());
#  else
				IndIterator it = pop.indIterator(sp->subPop());
#  endif
				size_t idIdx = pop.infoIdx(infoField(0));
				//UINT popPloidy = pop.ploidy();

				for (; it.valid(); ++it)
					for (vectoru::iterator p = ploidy.begin(); p != ploidy.end(); ++p) {
						int sign = m_mode == FROM_INFO ? 1 : (*p % 2 == 0 ? 1 : -1);
						long lineage = toLineage(it->info(idIdx) * sign);
						for (vectoru::const_iterator loc = loci.begin(); loc != loci.end(); ++loc)
							it->setAlleleLineage(lineage, *loc, static_cast<int>(*p));
					}
			}
		}
		pop.deactivateVirtualSubPop(sp->subPop());
	}
#endif  // LINEAGE
	(void)pop;            // avoid a compiler warning of unused variable
	return true;
}


}
