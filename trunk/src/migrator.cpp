/**
 *  $File: Migrator.cpp $
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

#include "virtualSubPop.h"
#include "migrator.h"

#include "boost_pch.hpp"
namespace simuPOP {

Migrator::Migrator(const floatMatrix & rate, int mode, const uintList & toSubPops,
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops, const stringList & infoFields)
	: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_rate(rate.elems()), m_mode(mode), m_to(toSubPops)
{
	DBG_FAILIF(mode != BY_IND_INFO && !subPops.empty() && subPops.size() != m_rate.size(),
		ValueError, "Length of param fromSubPop must match rows of rate matrix.");

	DBG_FAILIF(mode != BY_IND_INFO && !m_to.elems().empty() && m_to.elems().size() != m_rate[0].size(),
		ValueError, "Length of param toSubPop must match columns of rate matrix.");
}


string Migrator::describe(bool /* format */) const
{
	return "<simuPOP.Migrator>";
}


bool Migrator::apply(Population & pop) const
{
	// set info of individual
	size_t info = pop.infoIdx(infoField(0));

	subPopList fromSubPops = applicableSubPops(pop);

	DBG_FAILIF(m_mode != BY_IND_INFO && fromSubPops.size() != m_rate.size(),
		ValueError, "Number of 'from' subpopulations should match number of rows of migration rate matrix.");

	vectorf oldInfo;

	if (m_mode == BY_IND_INFO && !fromSubPops.empty()) {
		oldInfo.resize(pop.popSize());
		for (size_t i = 0; i < pop.popSize(); ++i)
			oldInfo[i] = pop.individual(static_cast<double>(i)).info(info);
	}

	if (m_mode != BY_IND_INFO || !fromSubPops.empty()) {
		for (size_t sp = 0; sp < pop.numSubPop(); ++sp) {
			RawIndIterator it = pop.rawIndBegin(sp);
			RawIndIterator it_end = pop.rawIndEnd(sp);
			if (numThreads() > 1) {
#ifdef _OPENMP
				size_t popSize = it_end - it;
#  pragma omp parallel firstprivate(it, it_end)
				{
					size_t id = omp_get_thread_num();
					it = it + id * (popSize / numThreads());
					it_end = id == numThreads() - 1 ? it_end : it + popSize / numThreads();
					for (; it != it_end; ++it)
						it->setInfo(static_cast<double>(sp), info);
				}
#endif
			} else {
				for (; it != it_end; ++it)
					it->setInfo(static_cast<double>(sp), info);
			}
		}

	}

	DBG_FAILIF(pop.hasActivatedVirtualSubPop(), RuntimeError,
		"Migration can not be applied to activated virtual subpopulations");

	// This one does not support ALL_AVAIL
	vectoru toSubPops = m_to.elems();
	if (m_to.allAvail())
		for (size_t i = 0; i < pop.numSubPop(); ++i)
			toSubPops.push_back(i);

	// real migration matrix might change from population to population because
	// of different number of subpopulations, and toSubPops can be ALL_AVAIL, and
	// then does not have to match subPops.
	matrixf migrationRate = m_rate;
	if (m_mode != BY_IND_INFO) {
		size_t szFrom = migrationRate.size();
		size_t szTo = migrationRate[0].size();

		// check parameters
		for (size_t i = 0; i < szFrom; ++i) {
			DBG_FAILIF(migrationRate[i].size() != szTo, ValueError,
				"Expecting a matrix of migration rate.");

			for (size_t j = 0; j < szTo; ++j) {
				DBG_FAILIF(fcmp_lt(migrationRate[i][j], 0.), ValueError,
					"Migration rate should be positive.");
				DBG_FAILIF(m_mode != BY_COUNTS && fcmp_gt(migrationRate[i][j], 1.), ValueError,
					"Migration rate should be in the range of [0,1]");
			}
		}

		// set r[i][i]--- may need to extend rate (to add i->i)
		if (m_mode == BY_PROBABILITY || m_mode == BY_PROPORTION) {
			for (size_t i = 0; i < szFrom; i++) {               // from
				// look for from=to cell.
				size_t spFrom = fromSubPops[i].subPop();
				double sum = accumulate(migrationRate[i].begin(), migrationRate[i].end(), 0.0);
				//
				vectoru::const_iterator spTo = find(toSubPops.begin(), toSubPops.end(), spFrom);
				if (spTo == toSubPops.end()) {                        // if no to, only check if sum <= 1
					if (fcmp_gt(sum, 1.0))
						throw ValueError("Sum of migrate rate from one subPop should <= 1");
					// adding i->i item
					migrationRate[i].push_back(1.0 - sum);
				} else {                                                          // if not, set r[i][i]
					double & self = migrationRate[i][ spTo - toSubPops.begin() ];
					sum -= self;
					if (fcmp_gt(sum, 1.0))
						throw ValueError("Sum of migrate rate from one subPop should <= 1");
					// reset to-my-self probability/proportion
					self = 1.0 - sum;
				}
			}
		}
	}

	for (size_t from = 0, fromEnd = fromSubPops.size(); from < fromEnd; ++from) {
		size_t spFrom = fromSubPops[from].subPop();
		// rateSize might be toSize + 1, the last one is from->from
		size_t toSize = toSubPops.size();
		size_t toIndex;

		// fromSubPops out of range....
		DBG_FAILIF(spFrom >= pop.numSubPop(), IndexError,
			(boost::format("Subpopulation index %1% out of range") % spFrom).str());

		size_t spSize = pop.subPopSize(fromSubPops[from]);

		if (fromSubPops[from].isVirtual())
			pop.activateVirtualSubPop(fromSubPops[from]);

		if (m_mode == BY_IND_INFO) {
			// restore information fields set by user so that other individuals
			// can stay at their original subpopulation.
			if (!oldInfo.empty()) {
				if (numThreads() > 1) {
#ifdef _OPENMP
#  pragma omp parallel
					{
						for (IndIterator ind = pop.indIterator(spFrom, omp_get_thread_num()); ind.valid(); ++ind)
							ind->setInfo(oldInfo[&*ind - &*pop.rawIndBegin()], info);
					}
#endif
				} else {
					for (IndIterator ind = pop.indIterator(spFrom); ind.valid(); ++ind)
						ind->setInfo(oldInfo[&*ind - &*pop.rawIndBegin()], info);
				}
			}
		} else if (m_mode == BY_PROBABILITY) {
			WeightedSampler ws(migrationRate[from]);

			// for each individual, migrate according to migration probability
			if (numThreads() > 1) {
#ifdef _OPENMP
#  pragma omp parallel private(toIndex)
				{
					for (IndIterator ind = pop.indIterator(spFrom, omp_get_thread_num()); ind.valid(); ++ind) {
						toIndex = ws.draw();
						DBG_ASSERT(toIndex < migrationRate[from].size(), ValueError,
							"Return index out of range.");
						if (toIndex < toSize && toSubPops[toIndex] != spFrom)
							ind->setInfo(static_cast<double>(toSubPops[toIndex]), info);
					}
				}
#endif
			} else {
				for (IndIterator ind = pop.indIterator(spFrom); ind.valid(); ++ind) {
					//toIndex = getRNG().randIntByFreq( rateSize, &migrationRate[from][0] ) ;
					toIndex = ws.draw();

					DBG_ASSERT(toIndex < migrationRate[from].size(), ValueError,
						"Return index out of range.");

					// rateSize == toSize (no i->i addition)
					//     toIndex < toSize
					// rateSize = toSize + 1, ignore i->1 (last one)
					//     toIndex < toSize
					if (toIndex < toSize && toSubPops[toIndex] != spFrom)
						ind->setInfo(static_cast<double>(toSubPops[toIndex]), info);
				}
			}
		} else {
			// 2nd, or 3rd method
			// first find out how many people will move to other subPop
			// then randomly assign individuals to move
			vectoru toNum(toSize);
			if (m_mode == BY_PROPORTION) {
				// in case that to sub is not in from sub, the last added
				// element is not used. sum of toNum is not spSize.
				for (size_t i = 0; i < toSize; ++i)
					toNum[i] = static_cast<ULONG>(spSize * migrationRate[from][i]);
			} else {                                                                      // by count
				for (size_t i = 0; i < toSize; ++i)
					toNum[i] = static_cast<ULONG>(migrationRate[from][i]);
			}
			// create a vector and assign indexes, then random shuffle
			// and assign info
			vectoru toIndices(spSize);
			size_t k = 0;
			for (size_t i = 0; i < toSize && k < spSize; ++i)
				for (size_t j = 0; j < toNum[i] && k < spSize; ++j)
					toIndices[k++] = toSubPops[i];

			// the rest of individuals will stay in their original subpopulation.
			while (k < spSize)
				toIndices[k++] = spFrom;

			getRNG().randomShuffle(toIndices.begin(), toIndices.end());
			IndIterator ind = pop.indIterator(spFrom);
			// set info
			for (size_t i = 0; ind.valid(); ++i, ++ind)
				// The previous migration_to value, if set by a previous vsp, will be overridden.
				ind->setInfo(static_cast<double>(toIndices[i]), info);
		}
		if (fromSubPops[from].isVirtual())
			pop.deactivateVirtualSubPop(spFrom);
	}   // for all subPop.

	// do migration.
	size_t oldNumSubPop = pop.numSubPop();
	pop.setSubPopByIndInfo(infoField(0));
	// try to keep number of subpopulations.
	if (pop.numSubPop() < oldNumSubPop) {
		vectorf split(oldNumSubPop - pop.numSubPop() + 1, 0);
		split[0] = static_cast<double>(pop.subPopSize(pop.numSubPop() - 1));
		pop.splitSubPop(pop.numSubPop() - 1, split);
	}
	DBG_ASSERT(pop.numSubPop() >= oldNumSubPop, RuntimeError,
		"Migrator should not decrease number of subpopulations.");

	return true;
}



BackwardMigrator::BackwardMigrator(const floatMatrix & rate, int mode, 
	int begin, int end, int step, const intList & at,
	const intList & reps, const subPopList & subPops, const stringList & infoFields)
	: BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
	m_rate(rate.elems()), m_inverse_rate(rate.elems().size(), rate.elems().size()), m_symmetric_matrix(true), m_mode(mode)
{
	DBG_FAILIF(!subPops.empty() && subPops.size() != m_rate.size(),
		ValueError, "Length of param subPop must match rows of rate matrix.");
	DBG_FAILIF(m_mode != BY_PROBABILITY && m_mode != BY_PROPORTION,
		ValueError, "Only BY_PROBABILITY and BY_PROPORTION is supported by BackMigrator");

	// find the inverse of B'
	size_t sz = m_rate.size();
	boost::numeric::ublas::matrix<double> Bt(sz, sz);
	for (size_t i = 0; i < sz; ++i) {
		if (m_rate[i].size() != sz) 
			throw ValueError("A m by m matrix is required for backward migration matrix.");
		//
		for (size_t j = 0; j < sz; ++j) {
			if (i == j) {
				Bt(i, i) = accumulate(m_rate[i].begin(), m_rate[i].end(), 0.0) - m_rate[i][i];
				Bt(i, i) = 1. - Bt(i, i);				
			} else
				Bt(i, j) = m_rate[j][i];
			if (fcmp_lt(Bt(i, j), 0.))
				throw ValueError((boost::format("Backward migration rate should be positive. %d provided.") % Bt(i,j)).str());
			if (fcmp_gt(Bt(i, j), 1.))
				throw ValueError((boost::format("Backward migration rate should be in the range of [0,1], %d provided.") % Bt(i,j)).str());
			if (j > i && m_rate[i][j] != m_rate[j][i])
				m_symmetric_matrix = false;
		}
	}
	// inverse
	boost::numeric::ublas::permutation_matrix<std::size_t> pm(sz);
	int res = lu_factorize(Bt, pm);
	if (res != 0)
		throw RuntimeError("Failed to convert backward matrix to forward migration matrix. (Matrix is not inversable).");
	m_inverse_rate.assign(boost::numeric::ublas::identity_matrix<double>(sz));
	// backsubstite to get the inverse		
	lu_substitute(Bt, pm, m_inverse_rate);

	/*
	DBG_DO(DBG_MIGRATOR, cerr << "Inverse of B' is " << m_inverse_rate << endl);
	*/
}


string BackwardMigrator::describe(bool /* format */) const
{
	return "<simuPOP.BackwardMigrator>";
}


bool BackwardMigrator::apply(Population & pop) const
{
	// set info of individual
	size_t info = pop.infoIdx(infoField(0));

	subPopList VSPs = applicableSubPops(pop);
	if (VSPs.size() <= 1)
		return true;
	
	DBG_FAILIF(VSPs.size() != m_rate.size(),
		ValueError, "Number of 'from' subpopulations should match number of rows of migration rate matrix.");
	
	vectoru subPops;
	for (size_t i = 0; i < VSPs.size(); ++i) {
		DBG_FAILIF(VSPs[i].isVirtual(), ValueError, 
			"BackwardMigrator does not support virtual subpupulations.")
		DBG_FAILIF(m_rate[i].size() != VSPs.size(), ValueError,
			"A square matrix is required for BackwardMigrator")
		subPops.push_back(VSPs[i].subPop());
	}

	// assign individuals their own subpopulation ID
	for (size_t sp = 0; sp < pop.numSubPop(); ++sp) {
		RawIndIterator it = pop.rawIndBegin(sp);
		RawIndIterator it_end = pop.rawIndEnd(sp);
		if (numThreads() > 1) {
#ifdef _OPENMP
			size_t popSize = it_end - it;
#  pragma omp parallel firstprivate(it, it_end)
			{
				size_t id = omp_get_thread_num();
				it = it + id * (popSize / numThreads());
				it_end = id == numThreads() - 1 ? it_end : it + popSize / numThreads();
				for (; it != it_end; ++it)
					it->setInfo(static_cast<double>(sp), info);
			}
#endif
		} else {
			for (; it != it_end; ++it)
				it->setInfo(static_cast<double>(sp), info);
		}
	}

	DBG_FAILIF(pop.hasActivatedVirtualSubPop(), RuntimeError,
		"Migration can not be applied to activated virtual subpopulations");

	// subpopulation size S (before migration)
	vectoru S;
	for (size_t i = 0; i < subPops.size(); ++i)
		S.push_back(pop.subPopSize(subPops[i]));

	// now, we need to calculate a forward migration matrix from the backward one
	// the formula is

	// S' = B'^-1 * S 
	// F = diag(S)^-1 B' diag(S')
	//
	// but for the special case of equal population size and symmetric matrix,
	// two matrices are the same
	bool simple_case = m_symmetric_matrix;
	// symmtrix matrix, check if equal population size
	if (simple_case) {
		for (size_t i = 1; i < subPops.size(); ++i)
			if (S[i] != S[i-1])
				simple_case = false;
	}
	
	size_t sz = m_rate.size();
	// if not the simple case, we need to calculate rate...
	matrixf migrationRate;
	if (simple_case)
		migrationRate = m_rate;
	else {
		// with Bt^-1, we can calculate expected population size
		vectorf Sp(sz);
		for (size_t i = 0; i < sz; ++i) {
			Sp[i] = 0;
			for (size_t j = 0; j < sz; ++j)
				Sp[i] += m_inverse_rate(i, j) * S[j];
		}
		DBG_DO(DBG_MIGRATOR, cerr << "Expected next population size is " << Sp << endl);
		for (size_t i = 0; i < sz; ++i) {
			if (Sp[i] <= 0)
				throw RuntimeError((boost::format("Failed to calculate forward migration matrix: negative expected population size %f for subpopulation %d")
					% Sp[i] % i).str());
		}
		// now, F = diag(S)^-1 * BT * diag (S')
		migrationRate = m_rate;
		for (size_t i = 0; i < sz; ++i)
			for (size_t j = 0; j < sz; ++j)
				// m_rate[i][i] might not be defined.
				migrationRate[i][j] = m_rate[j][i] * Sp[j] / S[i];
	}

	// check parameters
	for (size_t i = 0; i < sz; ++i) {
		for (size_t j = 0; j < sz; ++j) {
			if (fcmp_lt(migrationRate[i][j], 0.))
				throw ValueError("Converted forward migration rate should be positive.");
			if (fcmp_gt(migrationRate[i][j], 1.))
				throw ValueError("Converted forward migration rate should be in the range of [0,1]");
		}
		// look for from=to cell.
		double sum = accumulate(migrationRate[i].begin(), migrationRate[i].end(), 0.0);
		//
		double & self = migrationRate[i][i];
		sum -= self;
		if (fcmp_gt(sum, 1.0))
			throw ValueError("Sum of migrate rate from one subPop should <= 1");
		// reset to-my-self probability/proportion
		self = 1.0 - sum;

	}
	if (! simple_case) {
		DBG_DO(DBG_MIGRATOR, cerr << "Forward migration matrix is " << migrationRate << endl);
	}

	for (size_t from = 0, fromEnd = subPops.size(); from < fromEnd; ++from) {
		size_t spFrom = subPops[from];
		// rateSize might be toSize + 1, the last one is from->from
		size_t toIndex;

		size_t spSize = pop.subPopSize(spFrom);

		if (m_mode == BY_PROBABILITY) {
			WeightedSampler ws(migrationRate[from]);

			// for each individual, migrate according to migration probability
			for (IndIterator ind = pop.indIterator(spFrom); ind.valid(); ++ind) {
				//toIndex = getRNG().randIntByFreq( rateSize, &migrationRate[from][0] ) ;
				toIndex = ws.draw();

				DBG_ASSERT(toIndex < migrationRate[from].size(), ValueError,
					"Return index out of range.");

				if (toIndex < sz && subPops[toIndex] != spFrom)
					ind->setInfo(static_cast<double>(subPops[toIndex]), info);
			}
		} else {
			// 2nd, or 3rd method
			// first find out how many people will move to other subPop
			// then randomly assign individuals to move
			vectoru toNum(sz);
			// in case that to sub is not in from sub, the last added
			// element is not used. sum of toNum is not spSize.
			for (size_t i = 0; i < sz; ++i)
				toNum[i] = static_cast<ULONG>(spSize * migrationRate[from][i]);
			// create a vector and assign indexes, then random shuffle
			// and assign info
			vectoru toIndices(spSize);
			size_t k = 0;
			for (size_t i = 0; i < sz && k < spSize; ++i)
				for (size_t j = 0; j < toNum[i] && k < spSize; ++j)
					toIndices[k++] = subPops[i];

			// the rest of individuals will stay in their original subpopulation.
			while (k < spSize)
				toIndices[k++] = spFrom;

			getRNG().randomShuffle(toIndices.begin(), toIndices.end());
			IndIterator ind = pop.indIterator(spFrom);
			// set info
			for (size_t i = 0; ind.valid(); ++i, ++ind)
				// The previous migration_to value, if set by a previous vsp, will be overridden.
				ind->setInfo(static_cast<double>(toIndices[i]), info);
		}
	}   // for all subPop.

	// do migration.
	size_t oldNumSubPop = pop.numSubPop();
	pop.setSubPopByIndInfo(infoField(0));
	// try to keep number of subpopulations.
	if (pop.numSubPop() < oldNumSubPop) {
		vectorf split(oldNumSubPop - pop.numSubPop() + 1, 0);
		split[0] = static_cast<double>(pop.subPopSize(pop.numSubPop() - 1));
		pop.splitSubPop(pop.numSubPop() - 1, split);
	}
	DBG_ASSERT(pop.numSubPop() >= oldNumSubPop, RuntimeError,
		"Migrator should not decrease number of subpopulations.");

	return true;
}


struct compareVSP
{
	int operator()(const vspID & v1, const vspID & v2)
	{
		return v1.subPop() > v2.subPop();
	}


};

bool SplitSubPops::apply(Population & pop) const
{
	subPopList subPops = applicableSubPops(pop);

	// we have to split from top to bottom subpopulations
	// because split will change subpopulation index
	std::sort(subPops.begin(), subPops.end(), compareVSP());

	for (size_t i = 0; i < subPops.size(); ++i) {
		size_t sp = subPops[i].subPop();
		if (pop.subPopSize(sp) == 0)
			continue;

		// randomize indiviudlas
		if (m_randomize) {
			getRNG().randomShuffle(pop.rawIndBegin(sp), pop.rawIndEnd(sp));
			pop.setIndOrdered(false);
		}

		if (!m_subPopSizes.empty()) {
			vectorf sz(m_subPopSizes.size());
			for (size_t i = 0; i < sz.size(); ++i)
				sz[i] = static_cast<double>(m_subPopSizes[i]);
			pop.splitSubPop(sp, sz, m_names);
		} else if (!m_proportions.empty()) {
			pop.splitSubPop(sp, m_proportions, m_names);
		} else {
			// using an information field
			size_t idx = pop.infoIdx(infoField(0));
			std::sort(pop.rawIndBegin(sp), pop.rawIndEnd(sp), indCompare(idx));
			ConstRawIndIterator it = pop.rawIndBegin(sp);
			ConstRawIndIterator it_end = pop.rawIndEnd(sp);
			vectorf spSize(1, 1);
			double value = it->info(idx);
			for (++it; it != it_end; ++it) {
				if (fcmp_ne(it->info(idx), value)) {
					spSize.push_back(1);
					value = it->info(idx);
				} else
					++spSize.back();
			}
			pop.splitSubPop(sp, spSize, m_names);
		}
	}
	return true;
}


bool MergeSubPops::apply(Population & pop) const
{
	const subPopList sp = applicableSubPops(pop);

	vectoru subPops(sp.size());

	for (size_t i = 0; i < sp.size(); ++i)
		subPops[i] = sp[i].subPop();
	pop.mergeSubPops(subPops, m_name);
	return true;
}


bool ResizeSubPops::apply(Population & pop) const
{
	vectoru newSizes = pop.subPopSizes();

	const subPopList subPops = applicableSubPops(pop);

	DBG_FAILIF(subPops.empty() && m_sizes.size() != pop.numSubPop()
		&& m_proportions.size() != pop.numSubPop(), ValueError,
		"Please specify new subpopulation size for all subpopulations.");

	for (size_t i = 0; i < subPops.size(); ++i) {
		DBG_FAILIF(static_cast<size_t>(subPops[i].subPop()) >= pop.numSubPop(), IndexError,
			(boost::format("Subpopulation index %1% out of range of 0 ~ %2%") % subPops[i].subPop()
			 % (pop.numSubPop() - 1)).str());
		if (m_sizes.empty())
			newSizes[subPops[i].subPop()] = static_cast<ULONG>(newSizes[subPops[i].subPop()] * m_proportions[i]);
		else
			newSizes[subPops[i].subPop()] = m_sizes[i];
	}
	DBG_DO(DBG_MIGRATOR, cerr << "Resize subpopulations to size " << newSizes << endl);
	pop.resize(newSizes, m_propagate);
	return true;
}


}
