/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu                                                        *
*                                                                         *
*   $LastChangedDate$
*   $Rev$
*
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
*   This program is distributed in the hope that it will be useful,       *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*   GNU General Public License for more details.                          *
*                                                                         *
*   You should have received a copy of the GNU General Public License     *
*   along with this program; if not, write to the                         *
*   Free Software Foundation, Inc.,                                       *
*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/

#include "migrator.h"

namespace simuPOP {
void migrator::setRates(const matrix & rate, int mode)
{
	if (rate.empty() )
		return;

	UINT szFrom = rate.size();
	UINT szTo = rate[0].size();

	if (m_from.empty() )
		for (UINT i = 0; i < szFrom; ++i)
			m_from.push_back(i);

	if (m_to.empty() )
		for (UINT i = 0; i < szTo; ++i)
			m_to.push_back(i);

	m_mode = mode;

	if (m_mode != MigrByProbability && m_mode != MigrByProportion
	    && m_mode != MigrByCounts)
		throw ValueError("Migration mode can only be MigrByProbability), "
		                 " MigrByProportion or MigrByCounts");

	// check parameters
	for (UINT i = 0; i < szFrom; ++i) {
		if (rate[i].size() != szTo)
			throw ValueError("Expecting a matrix of migration rate.");

		for (size_t j = 0; j < szTo; ++j)
			if (fcmp_lt(rate[i][j], 0.)
			    || (m_mode != MigrByCounts && fcmp_gt(rate[i][j], 1.)))
				throw ValueError("Migration rate should be in the range of [0,1]");
	}

	m_rate = rate;

	// set r[i][i]--- may need to extend rate (to add i->i)
	if (m_mode == MigrByProbability || m_mode == MigrByProportion) {
		for (UINT i = 0; i < szFrom; i++) {               // from
			// look for from=to cell.
			UINT spFrom = m_from[i];
			double sum = accumulate(m_rate[i].begin(), m_rate[i].end(), 0.0);
			//
			vectoru::iterator spTo = find(m_to.begin(), m_to.end(), spFrom);
			if (spTo == m_to.end() ) {                        // if no to, only check if sum <= 1
				if (fcmp_gt(sum, 1.0) )
					throw ValueError("Sum of migrate rate from one subPop should <= 1");
				// adding i->i item
				m_rate[i].push_back(1.0 - sum);
			} else {                                                          // if not, set r[i][i]
				double & self = m_rate[i][ spTo - m_to.begin() ];
				sum -= self;
				if (fcmp_gt(sum, 1.0) )
					throw ValueError("Sum of migrate rate from one subPop should <= 1");
				// reset to-my-self probability/proportion
				self = 1.0 - sum;
			}
		}
	}
}


bool migrator::apply(population & pop)
{
	// set info of individual
	pop.setIndSubPopIDWithID();

	vectorlu toIndices(0);

	Weightedsampler ws(rng());

	for (UINT from = 0, fromEnd = m_from.size(); from < fromEnd; ++from) {
		UINT spFrom = m_from[from];
		// rateSize might be toSize + 1, the last one is from->from
		UINT toSize = m_to.size(), toIndex;

		// m_from out of range.... ignore.
		if (spFrom >= pop.numSubPop() )
			continue;

		if (m_mode == MigrByProbability) {                // migrate by probability
			ws.set(m_rate[from]);

			// for each individual, migrate according to migration probability
			for (IndIterator ind = pop.indBegin(spFrom); ind.valid();  ++ind) {
				//toIndex = rng().randIntByFreq( rateSize, &m_rate[from][0] ) ;
				toIndex = ws.get();

				DBG_ASSERT(toIndex < m_rate[from].size(), ValueError,
				    "Return index out of range.");

				// rateSize==toSize (no i->i addition)
				//   toIndex < toSize
				// rateSize = toSize + 1, ignore i->1 (last one)
				//  toIndex < toSize
				if (toIndex < toSize && m_to[toIndex] != spFrom)
					ind->setSubPopID(m_to[toIndex]);
			}
			continue;
		}

		// 2nd, or 3rd method
		// first find out how many people will move to other subPop
		// then randomly assign individuals to move
		vectorlu toNum(toSize);
		ULONG spSize = pop.subPopSize(spFrom);
		if (m_mode == MigrByProportion) {
			for (UINT i = 0; i < toSize; ++i)
				toNum[i] = static_cast<ULONG>(spSize * m_rate[from][i]);
		} else {                                                                  // by count
			for (UINT i = 0; i < toSize; ++i)
				toNum[i] = static_cast<ULONG>(m_rate[from][i]);
		}
		// create a vector and assign indexes, then random shuffle
		// and assign info
		toIndices.resize(spSize);
		UINT k = 0;
		for (UINT i = 0; i < toSize && k < spSize; ++i)
			for (UINT j = 0; j < toNum[i] && k < spSize; ++j)
				toIndices[k++] = m_to[i];

		while (k < spSize)
			toIndices[k++] = spFrom;

		random_shuffle(toIndices.begin(), toIndices.end());
		IndIterator ind = pop.indBegin(spFrom);
		// set info
		for (UINT i = 0; i < spSize; ++i)
			// SubPopID is signed short, to save a few bits
			(ind + i)->setSubPopID(static_cast<SubPopID>(toIndices[i]));
	}                                                                                         /// for all subPop.

	// do migration.
	// true: rearrange individuals
	pop.setSubPopByIndID();

	return true;
}


bool pyMigrator::apply(population & pop)
{
	if (m_rateFunc != NULL) {
		// get rate,
		matrix rate;
		PyObject * curSize = PyTuple_New(pop.numSubPop());

		for (size_t i = 0; i < pop.numSubPop(); ++i)
			PyTuple_SetItem(curSize, i, PyInt_FromLong(pop.subPopSize(i)));

		DBG_DO(DBG_MIGRATOR, cout << "Current population size " << pop.subPopSizes() << endl);

		PyCallFunc2(m_rateFunc, "(iO)", pop.gen(), curSize, rate, PyObj_As_Matrix);

		DBG_DO(DBG_MIGRATOR, cout << "Migration rate: " << rate << endl);
		Py_XDECREF(curSize);
		//
		setRates(rate, m_mode);
		// apply migrator
		return migrator::apply(pop);
	}
	// now, m_indFunc
	//
	// if loci is given
	vectora alleles;
	PyObject * numArray = NULL;
	UINT pld = pop.ploidy();
	//
	if (!m_loci.empty()) {
		alleles.resize(m_loci.size() * pld);
		numArray = Allele_Vec_As_NumArray(alleles.begin(), alleles.end());
	}
	// call the python function, pass the each individual to it.
	// get pop object
	for (IndIterator it = pop.indBegin(); it.valid(); ++it) {
		PyObject * indObj = pyIndObj(static_cast<void *>(& * it));
		// if pop is valid?
		if (indObj == NULL)
			throw SystemError("Could not pass population to the provided function. \n"
			                  "Compiled with the wrong version of SWIG?");

		// loci
		if (!m_loci.empty()) {
			for (size_t i = 0, iEnd = m_loci.size(), j = 0; i < iEnd; ++i)
				for (UINT p = 0; p < pld; ++p)
					alleles[j++] = it->allele(m_loci[i], p);
		}
		// hold the return subpopulation id.
		int resID = 0;
		// parenthesis is needed since PyCallFuncX are macros.
		if (m_param == NULL) {
			if (m_loci.empty()) {
				PyCallFunc(m_indFunc, "(O)", indObj, resID, PyObj_As_Int);
			} else {
				PyCallFunc2(m_indFunc, "(OO)", indObj, numArray, resID, PyObj_As_Int);
			}
		} else {
			if (m_loci.empty()) {
				PyCallFunc2(m_indFunc, "(OO)", indObj, m_param, resID, PyObj_As_Int);
			} else {
				PyCallFunc3(m_indFunc, "(OOO)", indObj, numArray, m_param, resID, PyObj_As_Int);
			}
		}
		it->setSubPopID(resID);
		Py_DECREF(indObj);
	}
	// do migration.
	// true: rearrange individuals
	pop.setSubPopByIndID();
	return true;
}


bool splitSubPop::apply(population & pop)
{
	// randomize indiviudlas
	if (m_randomize) {
		// random shuffle individuals
		for (ULONG it = 0; it < pop.subPopSize(m_which); ++it)
			pop.ind(it, m_which).setSubPopID(static_cast<SubPopID>(rng().randInt(MaxSubPopID)));
		std::sort(pop.indBegin(m_which), pop.indEnd(m_which));
		// not actully required since spliSubPop will do this.
		// this is to remind myself this step is important.
		pop.setShallowCopied(true);
	}
	if (!m_subPopSizes.empty())
		pop.splitSubPop(m_which, m_subPopSizes, m_subPopID);
	else
		pop.splitSubPopByProportion(m_which, m_proportions, m_subPopID);
	return true;
}


}
