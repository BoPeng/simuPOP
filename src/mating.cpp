/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu                                                        *
*                                                                         *
*   $LastChangedDate: 2006-02-21 15:27:25 -0600 (Tue, 21 Feb 2006)        *
*   $Rev: 191$                                                            *
*                                                                         *
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

#include "mating.h"

namespace simuPOP {


numOffspringGenerator::numOffspringGenerator(double numOffspring, PyObject * numOffspringFunc, UINT maxNumOffspring,
                                             UINT mode) :
	m_numOffspring(numOffspring), m_numOffspringFunc(NULL), m_maxNumOffspring(maxNumOffspring),
	m_mode(mode), m_firstOffspring(true)
{
	if (numOffspringFunc != NULL) {
		if (!PyCallable_Check(numOffspringFunc))
			throw ValueError("Passed variable is not a callable python function.");

		Py_INCREF(numOffspringFunc);
		m_numOffspringFunc = numOffspringFunc;
	}
	DBG_FAILIF(mode == MATE_BinomialDistribution && maxNumOffspring < 2,
	    ValueError, "If mode is MATE_BinomialDistribution, maxNumOffspring should be > 1");
	DBG_FAILIF(mode == MATE_UniformDistribution && maxNumOffspring < static_cast<UINT>(numOffspring),
	    ValueError, "If mode is MATE_UniformDistribution, maxNumOffspring should be greater than numOffspring");

}


bool numOffspringGenerator::fixedFamilySize()
{
	// not eachFamily...
	return m_mode == MATE_NumOffspring;
}


ULONG numOffspringGenerator::numOffspring(int gen)
{
	static double numOS = 0.;

	// use the same numOffspings each generation
	if (!m_firstOffspring && m_mode == MATE_NumOffspring)
		return static_cast<UINT>(numOS);

	if (m_numOffspringFunc == NULL)
		numOS = m_numOffspring;
	else if (m_mode == MATE_NumOffspring ||           // case 1, first time
	         // case 2, all the time; case 3,4,5, first time
	         m_mode == MATE_NumOffspringEachFamily ||
	         ( m_firstOffspring &&
	          ( m_mode == MATE_GeometricDistribution  ||
	           m_mode == MATE_PoissonDistribution   ||
	           m_mode == MATE_BinomialDistribution ))) {
		PyCallFunc(m_numOffspringFunc, "(i)", gen, numOS, PyObj_As_Double);
	}

	m_firstOffspring = false;

	if (m_mode == MATE_NumOffspring || m_mode == MATE_NumOffspringEachFamily) {
		DBG_FAILIF(numOS < 1, ValueError, "Need at least one offspring.");
		return static_cast<UINT>(numOS);
	} else if (m_mode == MATE_GeometricDistribution) {
		DBG_FAILIF(fcmp_lt(numOS, 0) || fcmp_gt(numOS, 1.), ValueError,
		    "P for a geometric distribution should be within [0,1], given " + toStr(numOS));
		UINT nos = rng().randGeometric(numOS);
		return nos;
	} else if (m_mode == MATE_PoissonDistribution) {
		// FIXME: numOS is no longer the average.
		UINT nos = rng().randPoisson(numOS) + 1;
		return nos;
	} else if (m_mode == MATE_BinomialDistribution) {
		DBG_FAILIF(fcmp_lt(numOS, 0) || fcmp_gt(numOS, 1.), ValueError,
		    "P for a Bionomial distribution should be within [0,1], given " + toStr(numOS));
		DBG_FAILIF(m_maxNumOffspring < 1, ValueError,
		    "Max number of offspring should be greater than 1. Given "
		    + toStr(m_maxNumOffspring));
		UINT nos = rng().randBinomial(m_maxNumOffspring - 1, numOS) + 1;
		return nos;
	} else if (m_mode == MATE_UniformDistribution) {
		// max: 5
		// num: 2
		// randint(4)  ==> 0, 1, 2, 3
		// + 2 ==> 2, 3, 4, 5
		UINT nos = rng().randInt(m_maxNumOffspring - static_cast<unsigned long>(m_numOffspring) + 1)
		           + static_cast<UINT>(m_numOffspring);
		return nos;
	} else
		throw ValueError("Wrong mating numoffspring mode. Should be one of \n"
		    "MATE_NumOffspring, MATE_NumOffspringEachFamily and MATE_GEometricDistribution");
}


offspringGenerator::offspringGenerator(const population & pop,
                                       vector<baseOperator *> & ops)
	: m_ops(ops)
{
#ifndef OPTIMIZED
	m_genoStruIdx = pop.genoStruIdx();
#endif
	m_formOffGenotype = checkFormOffspringGenotype();
}


bool offspringGenerator::checkFormOffspringGenotype()
{
	for (vector<baseOperator *>::const_iterator iop = m_ops.begin();
	     iop != m_ops.end(); ++iop) {
		if ((*iop)->formOffGenotype())
			return false;
	}
	return true;
}


void cloneOffspringGenerator::generateOffspring(population & pop, individual * parent, individual *, UINT numOff,
                                                population::IndIterator & it)
{
	// if population has changed.
	DBG_FAILIF(m_genoStruIdx != pop.genoStruIdx(), ValueError,
	    "Offspring generator is used for two different types of populations");

	// generate numOff offspring per mating, or until it  reaches offEnd
	UINT count = 0;
	bool accept = true;
	while (count < numOff) {
		if (m_formOffGenotype)
			// use deep copy!!!!!!!
			it->copyFrom(*parent);

		accept = true;
		// apply during mating operators
		vector<baseOperator *>::iterator iop = m_ops.begin();
		vector<baseOperator *>::iterator iopEnd = m_ops.end();
		for (; iop != iopEnd;  ++iop) {
			try {
				// During mating operator might reject this offspring.
				if (!(*iop)->applyDuringMating(pop, it, parent, NULL)) {
					accept = false;
					break;
				}
			} catch (...) {
				cout << "DuringMating operator " << (*iop)->__repr__()
				     << " throws an exception." << endl << endl;
				throw;
			}
		}                                                                         // all during-mating operators
		if (accept) {
			it++;
			count++;
		}
	}
}


mendelianOffspringGenerator::mendelianOffspringGenerator(const population & pop,
                                                         vector<baseOperator *> & ops)
	: offspringGenerator(pop, ops), m_bt(rng())
{
	m_hasSexChrom = pop.sexChrom();
	if (m_formOffGenotype) {
		vectorf prob(2 * pop.numChrom(), 0.5);
		m_bt.setParameter(prob, pop.popSize());
	}
	m_chIdx = pop.chromIndex();
}


void mendelianOffspringGenerator::formOffspringGenotype(individual * parent,
                                                        population::IndIterator & it, int ploidy, bool setSex)
{
	// current parental ploidy (copy from which chromosome copy)
	int parPloidy = 0;
	// pointer to parental, and offspring chromosome copies
	GenoIterator par[2];
	GenoIterator off;

	par[0] = parent->genoBegin(0);
	par[1] = parent->genoBegin(1);
	off = it->genoBegin(ploidy);

	UINT chEnd = parent->numChrom();
	int btShift = ploidy * chEnd;
#ifndef BINARYALLELE
	// the easy way to copy things.
	for (UINT ch = 0; ch < chEnd; ++ch) {
		parPloidy = m_bt.trialSucc(ch + btShift);
		for (size_t gt = m_chIdx[ch]; gt < m_chIdx[ch + 1]; ++gt)
			off[gt] = par[parPloidy][gt];
	}
#else
	//
	// 1. try to copy in blocks,
	// 2. if two chromosomes can be copied together, copy together
	// 3. if length is short, using the old method.
	//
	size_t parBegin = 0;
	size_t parEnd = 0;
	// first chromosome
	parPloidy = m_bt.trialSucc(btShift);
	//
	int nextParPloidy = 0;
	bool copyPar;
	for (UINT ch = 0; ch < chEnd; ++ch) {
		// if it is the last chromosome, copy anyway
		if (ch == chEnd - 1)
			copyPar = true;
		else {                                                                    // is there a different chromosome?
			nextParPloidy = m_bt.trialSucc(ch + 1 + btShift);
			copyPar = parPloidy != nextParPloidy;
		}
		if (copyPar) {
			// end of this chromosome, is the beginning of the next
			parEnd = m_chIdx[ch + 1];
			size_t length = parEnd - parBegin;
			//
			// the easiest case, try to get some speed up...
			if (length == 1)
				off[parBegin] = par[parPloidy][parBegin];
			else
				copyGenotype(par[parPloidy] + parBegin, off + parBegin, length);

			if (ch != chEnd - 1)
				parPloidy = nextParPloidy;
			parBegin = parEnd;
		}
	}
#endif

	if (setSex) {
		// last chromosome (sex chromosomes) determine sex
		if (m_hasSexChrom)
			it->setSex(parPloidy == 1 ? Male : Female);
		else
			it->setSex(rng().randInt(2) == 0 ? Male : Female);
	}
}


void mendelianOffspringGenerator::generateOffspring(population & pop, individual * dad, individual * mom,
                                                    UINT numOff, population::IndIterator & it)
{
	DBG_FAILIF(m_genoStruIdx != pop.genoStruIdx(), ValueError,
	    "Offspring generator is used for two different types of populations");

	// generate numOffspring offspring per mating
	UINT count = 0;
	bool accept = true;

	while (count < numOff) {
		if (m_formOffGenotype) {
			// m_bt 's width is 2*numChrom() and can be used for
			// the next two functions.
			m_bt.trial();
			formOffspringGenotype(mom, it, 0, false);
			formOffspringGenotype(dad, it, 1, true);
		}

		accept = true;
		// apply all during mating operators
		for (vector<baseOperator *>::iterator iop = m_ops.begin(), iopEnd = m_ops.end(); iop != iopEnd;  ++iop) {
			try {
				// During mating operator might reject this offspring.
				if (!(*iop)->applyDuringMating(pop, it, dad, mom)) {
					accept = false;
					break;
				}
			} catch (...) {
				cout << "DuringMating operator " << (*iop)->__repr__() << " throws an exception." << endl << endl;
				throw;
			}
		}
		if (accept) {
			it++;
			count++;
		}
	}                                                                                         // one offspring is successfully generated
}


void selfingOffspringGenerator::generateOffspring(population & pop, individual * parent,
                                                  individual *, UINT numOff, population::IndIterator & it)
{
	DBG_FAILIF(m_genoStruIdx != pop.genoStruIdx(), ValueError,
	    "Offspring generator is used for two different types of populations");

	// generate numOffspring offspring per mating
	UINT count = 0;
	bool accept = true;

	while (count < numOff) {
		if (m_formOffGenotype) {
			// m_bt 's width is 2*numChrom() and can be used for
			// the next two functions.
			m_bt.trial();
			// use the same parent to produce to copies of chromosomes
			formOffspringGenotype(parent, it, 0, false);
			formOffspringGenotype(parent, it, 1, true);
		}

		accept = true;
		// apply all during mating operators
		for (vector<baseOperator *>::iterator iop = m_ops.begin(), iopEnd = m_ops.end(); iop != iopEnd;  ++iop) {
			try {
				// During mating operator might reject this offspring.
				if (!(*iop)->applyDuringMating(pop, it, parent, NULL)) {
					accept = false;
					break;
				}
			} catch (...) {
				cout << "DuringMating operator " << (*iop)->__repr__() << " throws an exception." << endl << endl;
				throw;
			}
		}
		if (accept) {
			it++;
			count++;
		}
	}                                                                                         // one offspring is successfully generated
}


randomParentChooser::randomParentChooser(population & pop, size_t sp)
	: parentChooser(1), m_sampler(rng())
{
	m_selection = pop.selectionOn(sp);
	m_begin = pop.indBegin(sp);
	if (m_selection) {
		UINT fit_id = pop.infoIdx("fitness");
		// we need to order (true) the individual fitness.
		GappedInfoIterator fitness = pop.infoBegin(fit_id, true);
		// regardless of sex, get fitness for everyone.
		m_sampler.set(vectorf(fitness + pop.subPopBegin(sp),
		        fitness + pop.subPopEnd(sp) ) );
	} else
		m_size = pop.subPopSize(sp);
}


individual * randomParentChooser::chooseParent()
{
	// choose a parent
	if (m_selection)
		return & * (m_begin + m_sampler.get());
	else
		return & * (m_begin + rng().randInt(m_size));
}


randomParentsChooser::randomParentsChooser(population & pop, size_t sp)
	: parentChooser(2), m_maleIndex(0), m_femaleIndex(0),
	m_maleFitness(0), m_femaleFitness(0),
	m_malesampler(rng()), m_femalesampler(rng())
{
	ULONG spSize = pop.subPopSize(sp);

	// m_begin is relative to the beginning of the whole population.
	m_begin = pop.indBegin(sp);

	m_numMale = 0;
	m_numFemale = 0;
	population::IndIterator it = pop.indBegin(sp);
	population::IndIterator itEnd = pop.indEnd(sp);
	for (; it < itEnd; ++it)
		if (it->sex() == Male)
			m_numMale++;

	// allocate memory at first for performance reasons
	m_maleIndex.resize(m_numMale);
	m_femaleIndex.resize(spSize - m_numMale);

	m_selection = pop.selectionOn(sp);
	UINT fit_id = 0;
	if (m_selection) {
		fit_id = pop.infoIdx("fitness");
		m_maleFitness.resize(m_numMale);
		m_femaleFitness.resize(spSize - m_numMale);
	}

	m_numMale = 0;
	m_numFemale = 0;

	size_t idx = 0;
	population::IndIterator ind = pop.indBegin(sp);
	population::IndIterator indEnd = pop.indEnd(sp);
	for (; ind < indEnd; ind++) {
		if (ind->sex() == Male) {
			m_maleIndex[m_numMale] = idx;
			if (m_selection)
				m_maleFitness[m_numMale] = ind->info(fit_id);
			m_numMale++;
		} else {
			m_femaleIndex[m_numFemale] = idx;
			if (m_selection)
				m_femaleFitness[m_numFemale] = ind->info(fit_id);
			m_numFemale++;
		}
		idx++;
	}

	if (m_selection) {
		m_malesampler.set(m_maleFitness);
		m_femalesampler.set(m_femaleFitness);
		DBG_DO(DBG_DEVEL, cout << "Male fitness " << m_maleFitness << endl);
		DBG_DO(DBG_DEVEL, cout << "Female fitness " << m_femaleFitness << endl);
	}

	DBG_ASSERT(m_numFemale + m_numMale == spSize, SystemError,
	    "Wrong number of male/female.");
}


parentChooser::individualPair randomParentsChooser::chooseParents()
{
	individual * dad = NULL;
	individual * mom = NULL;

	if (m_selection) {                                        // with selection
		// using weidhted sampler.
		if (m_numMale != 0)
			dad = & * (m_begin + m_maleIndex[m_malesampler.get()]);
		else
			dad = & * (m_begin + m_femaleIndex[m_femalesampler.get()]);

		if (m_numFemale != 0)
			mom = & * (m_begin + m_femaleIndex[m_femalesampler.get()]);
		else
			mom = & * (m_begin + m_maleIndex[m_malesampler.get()]);
	} else {
		// using random sample.
		if (m_numMale != 0)
			dad = & * (m_begin + m_maleIndex[rng().randInt(m_numMale)]);
		else
			dad = & * (m_begin + m_femaleIndex[rng().randInt(m_numFemale)]);

		if (m_numFemale != 0)
			mom = & * (m_begin + m_femaleIndex[rng().randInt(m_numFemale)]);
		else
			mom = & * (m_begin + m_maleIndex[rng().randInt(m_numMale)]);
	}
	return std::make_pair(dad, mom);
}


pyParentsChooser::pyParentsChooser(population & pop, size_t sp, PyObject * pc)
	: parentChooser(0), m_generator(NULL), m_parIterator(NULL)
{
#if PY_VERSION_HEX < 0x02040000
	throw SystemError("Your Python version does not have good support for generator"
	    " so this python parent chooser can not be used.");
#else
 #ifndef OPTIMIZED
	m_size = pop.subPopSize(sp);
 #endif
	m_begin = pop.indBegin(sp);

	PyObject * popObj = pyPopObj(static_cast<void *>(&pop));
	// if pop is valid?
	DBG_FAILIF(popObj == NULL, SystemError,
	    "Could not pass population to the provided function. \n"
	    "Compiled with the wrong version of SWIG?");
	PyObject * arglist = Py_BuildValue("(Oi)", popObj, pop.gen());
	m_generator = PyEval_CallObject(pc, arglist);
	Py_XDECREF(arglist);

	// test if m_generator is a generator
	DBG_ASSERT(PyGen_Check(m_generator), ValueError,
	    "Passed function is not a python generator");

	m_parIterator = PyObject_GetIter(m_generator);

	// test if m_parIterator is iteratable.
	DBG_FAILIF(m_parIterator == NULL, ValueError,
	    "Can not iterate through parent generator");
#endif
}


parentChooser::individualPair pyParentsChooser::chooseParents()
{
	PyObject * item = PyIter_Next(m_parIterator);

	vectori parents;
	int parent;

	if (PySequence_Check(item)) {
		PyObj_As_IntArray(item, parents);
		DBG_ASSERT(parents.size() == 2, ValueError,
		    "Returned parents indexes should have size 2");
#ifndef OPTIMIZED
		DBG_ASSERT(static_cast<unsigned>(parents[0]) < m_size
		    && static_cast<unsigned>(parents[1]) < m_size,
		    ValueError, "Returned parent index (" + toStr(parents[0])
		    + ", and " + toStr(parents[1]) +
		    ") is greater than subpopulation size " + toStr(m_size));
#endif
		DBG_DO(DBG_MATING, cout << "choose parents " << parents[0]
		                        << " and " << parents[1] << endl;);
		Py_DECREF(item);
		return std::make_pair(& * (m_begin + parents[0]),
		           & * (m_begin + parents[1]));
	} else if (PyInt_Check(item) || PyLong_Check(item)) {
		PyObj_As_Int(item, parent);
#ifndef OPTIMIZED
		DBG_ASSERT(static_cast<unsigned>(parent) < m_size,
		    ValueError, "Returned index (" + toStr(parent) +
		    ") is greater than subpopulation size " + toStr(m_size));
#endif
		DBG_DO(DBG_MATING, cout << "choose parent " << parent
		                        << endl;);
		Py_DECREF(item);
		return parentChooser::individualPair(& * (m_begin + parent), NULL);
	} else
		DBG_ASSERT(false, ValueError,
		    "Invalid type of returned parent index(es)");
}


mating::mating(vectorlu newSubPopSize, string newSubPopSizeExpr, PyObject * newSubPopSizeFunc)
	: m_subPopSize(newSubPopSize),
	m_subPopSizeExpr(newSubPopSizeExpr, ""), m_subPopSizeFunc(NULL)
{
	DBG_FAILIF(!m_subPopSizeExpr.empty() && newSubPopSizeFunc != NULL,
	    ValueError, "Please only specify one of newSubPopSizeExpr and newSubPopSizeFunc.");


	if (newSubPopSizeFunc != NULL) {
		if (!PyCallable_Check(newSubPopSizeFunc))
			throw ValueError("Passed variable is not a callable python function.");

		Py_INCREF(newSubPopSizeFunc);
		m_subPopSizeFunc = newSubPopSizeFunc;
	}
}


void mating::prepareScratchPop(population & pop, population & scratch)
{
	// use population structure of pop
	if (m_subPopSize.empty() && m_subPopSizeExpr.empty() && m_subPopSizeFunc == NULL)
		scratch.setSubPopStru(pop.subPopSizes(), true);
	else if (!m_subPopSize.empty())  // set subPoplation size
		scratch.setSubPopStru(this->m_subPopSize, true);
	// evaluate from an expression
	else if (!m_subPopSizeExpr.empty()) {
		m_subPopSizeExpr.setLocalDict(pop.dict());
		vectorf sizef = m_subPopSizeExpr.valueAsArray();
		vectorlu sz(sizef.size());

		for (size_t i = 0, iEnd = sizef.size(); i < iEnd; i++)
			sz[i] = static_cast<ULONG>(sizef[i]);

		scratch.setSubPopStru(sz, true);
	} else {                                                                            // use m_subPopSizeFunc
		// get generation number
		int gen = pop.gen();
		// convert current pop size to a tuple
		PyObject * curSize = PyTuple_New(pop.numSubPop());

		DBG_ASSERT(curSize != NULL, SystemError, "Can not convert current pop size to a list");

		for (size_t i = 0; i < pop.numSubPop(); ++i)
			PyTuple_SetItem(curSize, i, PyInt_FromLong(pop.subPopSize(i)));

		vectorf res;
		PyCallFunc2(m_subPopSizeFunc, "(iO)", gen, curSize, res, PyObj_As_Array);
		Py_XDECREF(curSize);

		vectorlu sz(res.size());

		for (size_t i = 0; i < res.size(); i++)
			sz[i] = static_cast<ULONG>(res[i]);

		// allow change of pop size of scratch
		scratch.setSubPopStru(sz, true);
	}
	DBG_DO(DBG_SIMULATOR, cout << "New subpop size " << scratch.subPopSizes() << endl);

	DBG_FAILIF(scratch.numSubPop() != pop.numSubPop(),
	    ValueError, "number of subPopulaitons must agree.\n Pre: "
	    + toStr(pop.numSubPop()) + " now: " + toStr(scratch.numSubPop() ));
}


duplicateSplitter::duplicateSplitter(vectori const & offWeights)
	: virtualSplitter(offWeights)
{
}


infoSplitter::infoSplitter(vectori const & offWeights)
	: virtualSplitter(offWeights)
{
}


proportionSplitter::proportionSplitter(vectori const & offWeights)
	: virtualSplitter(offWeights)
{
}


rangeSplitter::rangeSplitter(vectori const & offWeights)
	: virtualSplitter(offWeights)
{
}


// nomating does nothing but applying during-mating operators
bool noMating::mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit)
{
	// apply during mating operators
	if (!ops.empty() ) {
		for (population::IndIterator it = pop.indBegin(), itEnd = pop.indEnd(); it != itEnd;  ++it) {
			for (vector<baseOperator *>::iterator iop = ops.begin(), iopEnd = ops.end(); iop != iopEnd;  ++iop) {
				(*iop)->applyDuringMating(pop, it, NULL, NULL);
			}                                                                         // all during-mating operators
		}
	}
	return true;
}


///
bool binomialSelection::mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit)
{
	m_numOffGen.resetNumOffspring();
	// scrtach will have the right structure.
	prepareScratchPop(pop, scratch);

	DBG_DO(DBG_MATING, m_famSize.clear());

	cloneOffspringGenerator og(pop, ops);

	// for each subpopulation
	for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
		UINT spSize = pop.subPopSize(sp);
		if (spSize == 0)
			continue;

		randomParentChooser pc(pop, sp);

		// choose a parent and genereate m_numOffspring offspring
		population::IndIterator it = scratch.indBegin(sp);
		population::IndIterator itEnd = scratch.indEnd(sp);
		UINT numOS;
		while (it != itEnd) {
			individual * parent = pc.chooseParent();

			numOS = m_numOffGen.numOffspring(pop.gen());
			// avoid warning assocaited with 'numOS > itEnd - it'
			if (it + numOS > itEnd)
				numOS = itEnd - it;
			// record family size, for debug reasons.
			DBG_DO(DBG_MATING, m_famSize.push_back(numOS));
			//
			og.generateOffspring(pop, parent, NULL, numOS, it);
		}                                                                                           // all offspring
	}                                                                                               // all subpopulation.

	if (submit)
		submitScratch(pop, scratch);
	return true;
}


bool randomMating::mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit)
{
	m_numOffGen.resetNumOffspring();
	// scrtach will have the right structure.
	prepareScratchPop(pop, scratch);

	DBG_DO(DBG_MATING, m_famSize.clear());

	mendelianOffspringGenerator og(pop, ops);

	/// random mating happens within each subpopulation
	for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
		ULONG spSize = pop.subPopSize(sp);
		if (spSize == 0)
			continue;

		randomParentsChooser pc(pop, sp);
		/// now, all individuals of needToFind sex is collected
		if ( (pc.numMale() == 0 || pc.numFemale() == 0 ) && !m_contWhenUniSex) {
			throw ValueError("Subpopulation becomes uni-sex. Can not continue. \n"
			    "You can use ignoreParentsSex (do not check parents' sex) or \ncontWhenUnixSex "
			    "(same sex mating if have to) options to get around this problem.");
		}

		// generate scratch.subPopSize(sp) individuals.
		population::IndIterator it = scratch.indBegin(sp);
		population::IndIterator itEnd = scratch.indEnd(sp);
		UINT numOS;
		while (it != itEnd) {
			parentChooser::individualPair const parents = pc.chooseParents();

			// record family size (this may be wrong for the last family)
			numOS = m_numOffGen.numOffspring(pop.gen());
			if (it + numOS > itEnd)
				numOS = itEnd - it;
			// record family size (this may be wrong for the last family)
			DBG_DO(DBG_MATING, m_famSize.push_back(numOS));
			//
			og.generateOffspring(pop, parents.first, parents.second, numOS, it);
		}
	}                                                                                         // each subPop

	if (submit)
		submitScratch(pop, scratch);
	return true;
}


/// CPPONLY
void countAlleles(population & pop, int subpop, const vectori & loci, const vectori & alleles,
                  vectorlu & alleleNum)
{
	int pldy = pop.ploidy();

	alleleNum = vectorlu(loci.size(), 0L);
	for (size_t l = 0; l < loci.size(); ++l) {
		int loc = loci[l];
		Allele ale = alleles[l];
		// go through all alleles
		if (subpop == -1) {
			if (pop.shallowCopied()) {
				for (population::IndIterator it = pop.indBegin(); it < pop.indEnd(); ++it)
					for (int p = 0; p < pldy; ++p)
						if (it->allele(loc, p) == ale)
							alleleNum[l]++;
			} else {
				for (GappedAlleleIterator a = pop.alleleBegin(loc, false),
				                          aEnd = pop.alleleEnd(loc, false); a != aEnd; ++a)
					if (AlleleUnsigned(*a) == ale)
						alleleNum[l]++;
			}
		} else {
			if (pop.shallowCopied()) {
				for (population::IndIterator it = pop.indBegin(subpop); it < pop.indEnd(subpop); ++it)
					for (int p = 0; p < pldy; ++p)
						if (it->allele(loc, p) == ale)
							alleleNum[l]++;
			} else {
				for (GappedAlleleIterator a = pop.alleleBegin(loc, subpop, false),
				                          aEnd = pop.alleleEnd(loc, subpop, false); a != aEnd; ++a)
					if (AlleleUnsigned(*a) == ale)
						alleleNum[l]++;
			}
		}
	}
}


bool controlledMating::mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit)
{
	// first call the function and get the range
	vectorf freqRange;

	PyCallFunc(m_freqFunc, "(i)", pop.gen(),
	    freqRange, PyObj_As_Array);

	DBG_ASSERT(freqRange.size() == m_loci.size() || freqRange.size() == 2 * m_loci.size(),
	    ValueError, "Length of returned range should have the same or double the number of loci.");

	vectorlu alleleNum;

#ifndef OPTIMIZED
	// calculate allele frequen at these loci
	// do not consider subpop
	countAlleles(pop, -1, m_loci, m_alleles, alleleNum);

	DBG_DO(DBG_MATING, cout << "Range of allele frequencies at generation "
	                        << pop.gen() << " is " << freqRange << endl);
#endif

	// compare integer is easier and more acurate, and we need to make sure
	// there is one allele when the frequency is greater than 0.
	vectorlu alleleRange(2 * m_loci.size(), 0);
	if (freqRange.size() == m_loci.size() ) {
		for (size_t i = 0; i < m_loci.size(); ++i) {
			if (freqRange[i] < 0.)
				freqRange[i] = 0.;
			if (freqRange[i] > 1.)
				freqRange[i] = 1.;

			alleleRange[2 * i] = static_cast<ULONG>(freqRange[i] * pop.popSize() * pop.ploidy());
			if (freqRange[i] > 0 && alleleRange[2 * i] == 0)
				alleleRange[2 * i] = 1;
			// upper bound using range.
			alleleRange[2 * i + 1] = static_cast<ULONG>((freqRange[i] + m_range) * pop.popSize() * pop.ploidy()) + 1;
#ifndef OPTIMIZED
			if (alleleNum[i] == 0 && alleleRange[2 * i] > 0)
				throw ValueError("No allele exists so there is no way to reach specified allele frequency.\n"
				    "Locus " + toStr(m_loci[i]) + " at generation " + toStr(pop.gen()) );
#endif
		}
	} else {
		// returned values are [low1, high1, low2, high2...]
		for (size_t i = 0; i < m_loci.size(); ++i) {
			if (freqRange[2 * i] < 0.)
				freqRange[2 * i] = 0.;
			if (freqRange[2 * i + 1] > 1.)
				freqRange[2 * i + 1] = 1.;

			alleleRange[2 * i] = static_cast<ULONG>(freqRange[2 * i] * pop.popSize() * pop.ploidy());

			DBG_FAILIF(freqRange[2 * i] > freqRange[2 * i + 1], ValueError,
			    "Incorrect frequency range: " + toStr(freqRange[2 * i]) +
			    " - " + toStr(freqRange[2 * i + 1]));

			if (freqRange[2 * i] > 0 && alleleRange[2 * i] == 0)
				alleleRange[2 * i] = 1;
			alleleRange[2 * i + 1] = static_cast<ULONG>(freqRange[2 * i + 1] * pop.popSize() * pop.ploidy()) + 1;
#ifndef OPTIMIZED
			if (alleleNum[i] == 0 && alleleRange[2 * i] > 0)
				throw ValueError("No allele exists so there is no way to reach specified allele frequency.\n"
				    "Locus " + toStr(m_loci[i]) + " at generation " + toStr(pop.gen()) );
#endif
		}
	}

	while (true) {
		// do the mating
		m_matingScheme->mate(pop, scratch, ops, false);

		// check allele frequency
		countAlleles(scratch, -1, m_loci, m_alleles, alleleNum);

		DBG_DO(DBG_MATING, cout << "mating finished, new count "
		                        << alleleNum << " range " << alleleRange << endl);
		//
		bool succ = true;
		for (size_t i = 0; i < m_loci.size(); ++i) {
			if (alleleNum[i] < alleleRange[2 * i] || alleleNum[i] > alleleRange[2 * i + 1]) {
				succ = false;
				break;
			}
		}
		// cout << "succ " << succ << "submit " << submit << endl;
		if (succ && submit) {
			// cout << "success" << endl;
			m_matingScheme->submitScratch(pop, scratch);
			break;
		}
	}
	return true;
}


/// give expected frequency for the whole population, or all subpopulations
/// return expected number of alleles at each subpopulations.
/// CPPONLY
void getExpectedAlleles(population & pop, vectorf & expFreq, const vectori & loci,
                        const vectori & alleles, vectoru & expAlleles)
{
	// determine expected number of alleles of each allele
	// at each subpopulation.
	UINT nLoci = loci.size();
	UINT numSP = pop.numSubPop();
	size_t pldy = pop.ploidy();
	size_t i, p;

	for (i = 0; i < nLoci; ++i) {
		if (expFreq[i] < 0.)
			expFreq[i] = 0.;
		if (expFreq[i] > 1.)
			expFreq[i] = 1.;
	}
	//
	// in the order of Loc1: sp1, sp2, sp3, sloc2: p1, sp2, sp3
	expAlleles.resize(nLoci * numSP);
	// if exp frequencies in subpopulation is not specified.
	//  I need to find the current allele frequencies
	// and use them as proportions for the next generation.
	if (numSP > 1 && expFreq.size() == nLoci) {
		vectorlu totalAlleles(nLoci);
		for (i = 0; i < nLoci; ++i) {
			int locus = loci[i];
			Allele allele = alleles[i];

			totalAlleles[i] = static_cast<ULONG>(expFreq[i] * pop.popSize() * pop.ploidy());
			// make sure one allele (our seed :-) exists
			if (expFreq[i] > 0. && totalAlleles[i] == 0)
				totalAlleles[i] = 1;

			// determine the number alleles at each subpopulation.
			vectorf curFreq(numSP);
			ULONG numOfAlleles = 0;
			for (size_t sp = 0; sp < numSP; ++sp) {
				ULONG n = 0;
				if (pop.shallowCopied()) {
					for (population::IndIterator it = pop.indBegin(sp); it < pop.indEnd(sp); ++it)
						for (p = 0; p < pldy; ++p)
							if (it->allele(locus, p) == allele)
								n++;
				} else {
					for (GappedAlleleIterator a = pop.alleleBegin(locus, sp, false),
					                          aEnd = pop.alleleEnd(locus, sp, false); a != aEnd; ++a)
						if (AlleleUnsigned(*a) == allele)
							n++;
				}
				numOfAlleles += n;
				curFreq[sp] = double (n) / (pop.subPopSize(sp) * pldy);
			}

			DBG_DO(DBG_MATING, cout << "Cur freq " << curFreq << endl);

			// if there is no alleles
			if (numOfAlleles == 0 && expFreq[i] > 0.)
				throw ValueError("No disease allele exists, but exp allele frequency is greater than 0.\n"
				    " Generation " + toStr(pop.gen()) );

			// calculate exp number of affected offspring in the next generation.
			//
			// step 1: totalsize*expFreq is the total number of disease alleles
			// step 2: assign these alleles to each subpopulation according to a multi-nomial
			// distribution with p_i beging allele frequency at each subpopulation.
			// assign these numbers to each subpopulation
			rng().randMultinomial(static_cast<unsigned int>(pop.popSize() * expFreq[i] * pldy),
			    curFreq, expAlleles.begin() + numSP * i);
		}
	}
	// simpler case, one subpopulation, or with gieven allele frequency
	else if (expFreq.size() == numSP * nLoci) {
		for (i = 0; i < nLoci; ++i) {
			for (size_t sp = 0; sp < numSP; ++sp) {
#ifndef OPTIMIZED
				int locus = loci[i];
				Allele allele = alleles[i];
				ULONG n = 0;
				// go through all alleles
				if (pop.shallowCopied()) {
					for (population::IndIterator it = pop.indBegin(sp); it < pop.indEnd(sp); ++it)
						for (p = 0; p < pldy; ++p)
							if (it->allele(locus, p) == allele)
								n++;
				} else {
					for (GappedAlleleIterator a = pop.alleleBegin(locus, sp, false),
					                          aEnd = pop.alleleEnd(locus, sp, false); a != aEnd; ++a) {
						if (AlleleUnsigned(*a) == allele)
							n++;
					}
				}

				// if there is no alleles
				if (n == 0 && expFreq[numSP * i + sp] > 0.)
					throw ValueError("No disease allele exists, but exp allele frequency is greater than 0.\n"
					    " Generation " + toStr(pop.gen()) );
#endif
				expAlleles[numSP * i + sp] = static_cast<UINT>(pop.subPopSize(sp) * pldy * expFreq[numSP * i + sp]);
				if (expFreq[numSP * i + sp] > 0. && expAlleles[numSP * i + sp] == 0)
					expAlleles[numSP * i + sp] = 1;
			}
		}                                                                                 // each locus
	} else {
		throw ValueError("Returned expected frequency has wrong length");
	}
}


//
// This mating scheme is very complicated, it is similar to randomMating, but it tries to control
// the mating process so that the total disease allele frequency controlled to pre-specified values.
//
bool controlledRandomMating::mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit)
{
	m_numOffGen.resetNumOffspring();
	// scrtach will have the right structure.
	prepareScratchPop(pop, scratch);

	size_t pldy = pop.ploidy(), nLoci = m_loci.size();
	size_t i;

	// expected frequency at each locus
	vectorf expFreq;
	PyCallFunc(m_freqFunc, "(i)", pop.gen(), expFreq, PyObj_As_Array);
	DBG_DO(DBG_MATING, cout << "expected freq " << expFreq << endl);

	// determine expected number of alleles of each allele
	// at each subpopulation.
	UINT numSP = pop.numSubPop();
	vectoru expAlleles(nLoci * numSP);
	getExpectedAlleles(pop, expFreq, m_loci, m_alleles, expAlleles);
	DBG_DO(DBG_MATING, cout << "expected alleles " << expAlleles << endl);

	/// whether or not use stack.
	bool useStack = m_numOffGen.fixedFamilySize();

	mendelianOffspringGenerator og(pop, ops);
	// use to go through offspring generation to count alleles
	UINT totNumLoci = pop.totNumLoci();

	/// controlled random mating happens within each subpopulation
	for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
		ULONG spSize = pop.subPopSize(sp);
		if (spSize == 0)
			continue;

		// reset stack each time
		if (useStack)
			m_stack = stack<population::IndIterator>();

		// total allowed disease alleles.
		vectoru totAllele(nLoci);
		// currently available disease allele. (in the offspring generation)
		vectoru curAllele(nLoci, 0);
		for (i = 0; i < nLoci; ++i) {
			totAllele[i] = expAlleles[sp + numSP * i];
			if (totAllele[i] > spSize * pldy) {
				cout << "Warning: number of planned affected alleles exceed population size.";
				totAllele[i] = spSize * pldy;
			}
		}

		randomParentsChooser pc(pop, sp);
		if (pc.numMale() == 0 || pc.numFemale() == 0) {
			if (m_contWhenUniSex)
				cout << "Warning: the subpopulation is uni-sex. Mating will continue with same-sex mate" << endl;
			else
				throw ValueError("Subpopulation becomes uni-sex. Can not continue. \n"
				    "You can use ignoreParentsSex (do not check parents' sex) or \ncontWhenUnixSex "
				    "(same sex mating if have to) options to get around this problem.");
		}

		// ploidy
		vectori na(nLoci, 0);
		// it is possible that no disease allele is required
		// so everyone is accepted
		bool freqRequMet = true;
		for (i = 0; i < nLoci; ++i) {
			if (totAllele[i] > 0) {
				freqRequMet = false;
				break;
			}
		}
		// if a family has disease allele.
		bool hasAff;
		bool stackStage = false;
		// start!
		population::IndIterator it = scratch.indBegin(sp);
		population::IndIterator itEnd = scratch.indEnd(sp);
		population::IndIterator itBegin = it;
		UINT numOS = 0;
		/// if mor ethan noAAattempt times, no pure homo found,
		/// accept non-homo cases.
		int AAattempt = 200;
		while (true) {
			// the logic is complicated here and I will try to be explicit
			// if already in the stackStage, it is taken from the stack.
			if (useStack && stackStage) {
				if (m_stack.empty()) {
					DBG_ASSERT(!freqRequMet, SystemError, "Empty stack should only happen when freq requirement is not met.");
					cout << "Warning: frequency requirement is not met, for subpop " << sp << " at generation "
					     << pop.gen() << ".\nThis is usually caused by multiple high disease allele frequency." << endl;
					break;
				}
				it = m_stack.top();
			}

			// randomly choose parents

			parentChooser::individualPair parents = pc.chooseParents();

			// generate m_numOffspring offspring per mating
			// record family size (this may be wrong for the last family)

			// this is the basic scheme,
			// use stage 1, stage 2 as described in the paper

			itBegin = it;
			numOS = m_numOffGen.numOffspring(pop.gen());
			// note that this is still valid in stack stage
			if (it + numOS > itEnd)
				numOS = itEnd - it;
			// generate numOffspring offspring per mating
			// it moves forward
			og.generateOffspring(pop, parents.first, parents.second, numOS, it);

			// count alleles in this family
			// count number of alleles in the family.
			for (i = 0; i < nLoci; ++i)
				na[i] = 0;
			hasAff = false;
			/*
			   for(population::IndIterator tmp=itBegin; tmp != it; ++tmp)
			   {
			   	// success count na, nu
			   	for(i=0; i<nLoci; ++i)
			   	{
			   		for(p=0; p<pldy; ++p)
			   		{
			   			if(tmp->allele(m_loci[i], p) == m_alleles[i])
			   			{
			   				na[i]++;
			   hasAff = true;
			   }
			   }
			   }
			   }								  // family generated.
			 */
			// we know that scratch population has ordered linear genotype

			// this, to my surprise, does not improve performance.
			for (i = 0; i < nLoci; ++i) {
				GenoIterator ptr = itBegin->genoBegin() + m_loci[i];
				for (size_t j = 0; j < numOS * pldy; ++j, ptr += totNumLoci) {
					if (*ptr == m_alleles[i]) {
						na[i]++;
						hasAff = true;
					}
				}
			}

			if (useStack) {
				// now check if this family is usable.
				bool accept = false;
				// all disease alleles have been satidfied.
				// only accept unaffected families
				// otherwise, accept any family that can help.
				if (freqRequMet) {
					// and not it == itEnd, that is to say, we only accept unaffected
					if (!hasAff) {
						// has AA, so we have less reason to compromise
						// but it might still be very difficult to find one,
						// so we still allow accepting of "wrong" individuals.
						AAattempt = 10000;
						accept = true;
					}
					// tried 100 times, no AA is found.
					else if (AAattempt == 0) {
						AAattempt = 100;
						accept = true;
					}
					AAattempt--;
				} else {
					// we accept affected helpful ones, and purely unsffected ones,
					if (hasAff) {
						// has the right kind of mutant?
						for (i = 0; i < nLoci; ++i) {
							// accept the whole family, if we need this allele
							if (curAllele[i] < totAllele[i] && na[i] > 0) {
								accept = true;
								break;
							}
						}
					} else {
						// in the stack stage, we only accept affected.
						if (!stackStage)
							accept = true;
					}
				}

				// reject this family
				if (!accept) {
					// it relocate to its begin point
					// DBG_DO(DBG_MATING, cout << "Reject " << na << endl);
					it = itBegin;
					continue;
				}
				DBG_DO(DBG_DEVEL, cout << "Accept " << na << " CUR " << curAllele << " TOT  " << totAllele << endl);

				// accpet this family, see if all done.
				if (!freqRequMet) {
					freqRequMet = true;
					for (i = 0; i < nLoci; ++i) {
						curAllele[i] += na[i];
						if (curAllele[i] < totAllele[i])
							freqRequMet = false;
					}
				}
				if (!stackStage && !freqRequMet && !hasAff) {
					// this family is in stack, might be
					m_stack.push(itBegin);
					DBG_DO(DBG_DEVEL, cout << "Push in stack " << m_stack.size() << endl);
				}
				// accepted,
				if (stackStage) {
					m_stack.pop();
					DBG_DO(DBG_DEVEL, cout << "Pop index " << m_stack.size() << endl);
				}
				// see if break
				if (it == itEnd) {
					stackStage = true;
					DBG_DO(DBG_MATING, cout << "Stack stage " << m_stack.size() << endl);
				}
				if (freqRequMet && stackStage) {
					DBG_DO(DBG_MATING, cout << "Finish generating offspring" << endl);
					break;
				}
			} else {
				// now check if this family is usable.
				bool accept = false;
				// all disease alleles have been satidfied.
				// only accept unaffected families
				// otherwise, accept any family that can help.
				if (freqRequMet) {
					if (!hasAff) {
						// has AA, so no need to compromise
						AAattempt = 10000;
						accept = true;
					}
					// tried 100 times, no AA is found.
					else if (AAattempt == 0) {
						AAattempt = 100;
						accept = true;
					}
					AAattempt--;
				} else {                                                       // do not use stack
					for (i = 0; i < nLoci; ++i) {
						// accept the whole family, if we need this allele
						if (curAllele[i] < totAllele[i] && na[i] > 0) {
							accept = true;
							break;
						}
					}
				}

				// reject this family
				if (!accept) {
					// it relocate to its begin point
					// DBG_DO(DBG_MATING, cout << "Reject " << na << endl);
					it = itBegin;
					continue;
				}
				DBG_DO(DBG_MATING, if (it - scratch.indBegin(sp) < 100) cout << "Accept " << na << endl;);

				// accpet this family, see if all done.
				if (!freqRequMet) {
					freqRequMet = true;
					for (i = 0; i < nLoci; ++i) {
						curAllele[i] += na[i];
						if (curAllele[i] < totAllele[i])
							freqRequMet = false;
					}
				}
				// see if break
				if (it == itEnd) {
					if (!freqRequMet)
						cout << "Warning: frequency requirement is not met, for subpop " << sp << " at generation "
						     << pop.gen() << ".\nThis is usually caused by multiple high disease allele frequency." << endl;
					break;
				}
			}
		}                                                                                           // nostack scheme
	}                                                                                               // each subPop

	if (submit)
		submitScratch(pop, scratch);
	return true;
}


pyMating::pyMating(vectori const & parentChoosers,
                   vectorobj const & pyChoosers,
                   vectori const & offspringGenerators,
                   vector<virtualSplitter *> const & splitters,
                   double numOffspring,
                   PyObject * numOffspringFunc,
                   UINT maxNumOffspring,
                   UINT mode,
                   vectorlu newSubPopSize,
                   string newSubPopSizeExpr,
                   PyObject * newSubPopSizeFunc
                   )
	: mating(newSubPopSize, newSubPopSizeExpr, newSubPopSizeFunc),
	m_numOffGen(numOffspring,
	            numOffspringFunc, maxNumOffspring, mode),
	m_parentChoosers(parentChoosers),
	m_pyChoosers(),
	m_offspringGenerators(offspringGenerators),
	m_splitters()
{
	vectorobj::const_iterator it = pyChoosers.begin();
	vectorobj::const_iterator it_end = pyChoosers.end();

	for (; it != it_end; ++it) {
		m_pyChoosers.push_back(*it);
		Py_XINCREF(*it);
	}
	vector<virtualSplitter *>::const_iterator sit = splitters.begin();
	vector<virtualSplitter *>::const_iterator sit_end = splitters.end();
	for (; sit != sit_end; ++sit)
		m_splitters.push_back((*sit)->clone());
}


bool pyMating::mate(population & pop, population & scratch, vector<baseOperator *> & ops, bool submit)
{
	m_numOffGen.resetNumOffspring();
	// scratach will have the right structure.
	prepareScratchPop(pop, scratch);

	DBG_DO(DBG_MATING, m_famSize.clear());

	DBG_ASSERT(m_parentChoosers.size() == 1 || m_parentChoosers.size() == pop.numSubPop(),
	    ValueError, "You should specify one parentChoosers, or one parentChoosers for each subpopulation. ("
	    + toStr(m_parentChoosers.size()) + ")");

	DBG_ASSERT(m_offspringGenerators.size() == 1 || m_offspringGenerators.size() == pop.numSubPop(),
	    ValueError, "You should specify one offspringGenerator, or one offspringGenerator for each subpopulation");

	size_t pcIdx = 0;
	/// random mating happens within each subpopulation
	for (UINT sp = 0; sp < pop.numSubPop(); ++sp) {
		ULONG spSize = pop.subPopSize(sp);
		if (spSize == 0)
			continue;

		int pcType = m_parentChoosers.size() > 1 ? m_parentChoosers[sp] : m_parentChoosers[0];
		parentChooser * pc = NULL;

		switch (pcType) {
		case MATE_RandomParentChooser:
			pc = new randomParentChooser(pop, sp);
			break;
		case MATE_RandomParentsChooser:
			pc = new randomParentsChooser(pop, sp);
			break;
		case MATE_PyParentsChooser:
			pc = new pyParentsChooser(pop, sp,
			         m_pyChoosers.size() > 1 ? m_pyChoosers[pcIdx++] : m_pyChoosers[0]);
			break;
		default:
			throw ValueError("Unknown parentChoosers type");
		}

		int ogType = m_offspringGenerators.size() > 1 ? m_offspringGenerators[sp] : m_offspringGenerators[0];
		offspringGenerator * og = NULL;
		switch (ogType) {
		case MATE_CloneOffspringGenerator:
			og = new cloneOffspringGenerator(pop, ops);
			// numParents can be 0 and its parents are
			// determined dynamically
			DBG_FAILIF(pc->numParents() == 2, ValueError,
			    "Imcompatible parent chooser and offspring generator");
			break;
		case MATE_MendelianOffspringGenerator:
			og = new mendelianOffspringGenerator(pop, ops);
			// numParents can be 0 and its parents are
			// determined dynamically
			DBG_FAILIF(pc->numParents() == 1, ValueError,
			    "Imcompatible parent chooser and offspring generator");
			break;
		case MATE_SelfingOffspringGenerator:
			og = new selfingOffspringGenerator(pop, ops);
			// numParents can be 0 and its parents are
			// determined dynamically
			DBG_FAILIF(pc->numParents() == 2, ValueError,
			    "Imcompatible parent chooser and offspring generator");
			break;
		default:
			throw ValueError("Unknown offspring generator type");
		}


		// generate scratch.subPopSize(sp) individuals.
		population::IndIterator it = scratch.indBegin(sp);
		population::IndIterator itEnd = scratch.indEnd(sp);
		UINT numOS;
		while (it != itEnd) {
			individual * dad = NULL;
			individual * mom = NULL;
			if (pc->numParents() == 1)
				dad = pc->chooseParent();
			// 0 or 2, that is to say, dad or mom can be NULL
			else {
				parentChooser::individualPair const parents = pc->chooseParents();
				dad = parents.first;
				mom = parents.second;
			}

			// record family size (this may be wrong for the last family)
			numOS = m_numOffGen.numOffspring(pop.gen());
			if (it + numOS > itEnd)
				numOS = itEnd - it;
			// record family size (this may be wrong for the last family)
			DBG_DO(DBG_MATING, m_famSize.push_back(numOS));
			//
			og->generateOffspring(pop, dad, mom, numOS, it);
		}
		delete og;
		delete pc;
	}                                                                                         // each subPop

	if (submit)
		submitScratch(pop, scratch);
	return true;
}


}
