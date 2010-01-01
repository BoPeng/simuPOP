/**
 *  $File: simulator.cpp $
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

#include "simulator.h"

// for file compression
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

using std::ostringstream;

namespace simuPOP {

population & pyPopIterator::next()
{
	if (m_index == m_end)
		throw StopIteration("");
	else
		return **(m_index++);
}


simulator::simulator(const population & pop, UINT rep)
	: m_numRep(rep)
{
	DBG_ASSERT(m_numRep >= 1, ValueError,
		"Number of replicates should be greater than or equal one.");

	DBG_DO(DBG_SIMULATOR, cerr << "Creating simulator " << endl);


	// create replicates of given population
	m_ptrRep = vector<population *>(m_numRep);
	m_scratchPop = new population(pop);

	DBG_FAILIF(m_scratchPop == NULL,
		SystemError, "Fail to create scratch population");

	try {
		for (UINT i = 0; i < m_numRep; ++i) {
			m_ptrRep[i] = new population(pop);

			DBG_FAILIF(m_ptrRep[i] == NULL,
				SystemError, "Fail to create new replicate");

			// set replication number
			m_ptrRep[i]->setRep(i);
		}
	} catch (...) {
		cerr << "Can not create " << m_numRep << " populations" << endl;
		throw RuntimeError("Failed to create a population.");
	}

	// set generation number for all replicates
	DBG_DO(DBG_SIMULATOR, cerr << "simulator created" << endl);
}


simulator::~simulator()
{
	// call the destructor of each replicates
	delete m_scratchPop;

	for (UINT i = 0; i < m_numRep; ++i)
		delete m_ptrRep[i];
}


simulator::simulator(const simulator & rhs) :
	m_numRep(rhs.m_numRep),
	m_ptrRep(0),
	m_scratchPop(NULL)
{
	m_scratchPop = rhs.m_scratchPop->clone();
	m_ptrRep = vector<population *>(m_numRep);
	for (size_t i = 0; i < m_numRep; ++i) {
		m_ptrRep[i] = rhs.m_ptrRep[i]->clone();
		m_ptrRep[i]->setRep(i);
	}
}


simulator * simulator::clone() const
{
	return new simulator(*this);
}


population & simulator::pop(UINT rep) const
{
	DBG_FAILIF(rep >= m_numRep, IndexError,
		"replicate index out of range. From 0 to numRep()-1 ");

	return *m_ptrRep[rep];
}


population & simulator::extract(UINT rep)
{
	DBG_FAILIF(rep >= m_numRep, IndexError,
		"replicate index out of range. From 0 to numRep()-1 ");

	population * pop = m_ptrRep[rep];
	m_ptrRep.erase(m_ptrRep.begin() + rep);
	--m_numRep;
	return *pop;
}


void simulator::add(const population & pop)
{
	++m_numRep;

	m_ptrRep.push_back(new population(pop));
	DBG_FAILIF(m_ptrRep.back() == NULL,
		RuntimeError, "Fail to add new population.");
}


string simulator::describe(const opList & initOps,
                           const opList & preOps,
                           const mating & matingScheme,
                           const opList & postOps,
                           const opList & finalOps,
                           int gen)
{
	if (initOps.empty() && preOps.empty() && postOps.empty() && finalOps.empty() && gen == -1)
		return "<simuPOP.simulator> a simulator with " + toStr(m_numRep) + " population" + (m_numRep == 1 ? "." : "s.");

	vectorstr allDesc(m_numRep, "");

	// assuming all active replicates.
	vector<bool> activeReps(m_numRep);

	for (UINT curRep = 0; curRep < m_numRep; curRep++) {
		ostringstream desc;

		if (initOps.empty())
			desc << "No operator is used to initialize population (initOps).\n";
		else {
			desc << "Apply pre-evolution operators to the initial population (initOps).\n<ul>\n";
			for (size_t it = 0; it < initOps.size(); ++it)
				desc << "<li>" << initOps[it]->describe(false) << " " << initOps[it]->applicability(true, false) << endl;
			desc << "</ul>\n";
		}
		if (gen < 0)
			desc << "\nEvolve a population indefinitely until an operator determines it." << endl;
		else
			desc	<< "\nEvolve a population for " << gen << " generations" << endl;
		desc << "<ul>\n";
		if (preOps.empty())
			desc << "<li>No operator is applied to the parental generation (preOps)." << endl;
		else {
			desc << "<li>Apply pre-mating operators to the parental generation (preOps)\n<ul>\n";
			for (size_t it = 0; it < preOps.size(); ++it)
				if (preOps[it]->isActive(curRep, 0, 0, activeReps, true))
					desc << "<li>" << preOps[it]->describe(false) << " " << preOps[it]->applicability() << endl;
			desc << "</ul>\n";
		}
		desc	<< "\n<li>Populate an offspring populaton from the parental population using mating scheme "
		        << matingScheme.describe(false) << endl;
		//
		if (postOps.empty())
			desc << "\n<li>No operator is applied to the offspring population (postOps)." << endl;
		else {
			desc << "\n<li>Apply post-mating operators to the offspring population (postOps).\n<ul>\n";
			for (size_t it = 0; it < postOps.size(); ++it)
				if (postOps[it]->isActive(curRep, 0, 0, activeReps, true))
					desc << "<li>" << postOps[it]->describe(false) << " " << postOps[it]->applicability() << endl;
			desc << "</ul>\n";
		}
		desc << "</ul>\n\n";
		if (finalOps.empty() )
			desc << "No operator is applied to the final population (finalOps)." << endl;
		else {
			desc << "Apply post-evolution operators (finalOps)\n<ul>\n";
			for (size_t it = 0; it < finalOps.size(); ++it)
				desc << "<li>" << finalOps[it]->describe(false) << " " << finalOps[it]->applicability(true, false) << endl;
			desc << "</ul>\n";
		}
		allDesc[curRep] = desc.str();
	}
	ostringstream desc;
	vectoru reps;
	for (UINT curRep = 0; curRep < m_numRep; curRep++) {
		if (reps.empty())
			reps.push_back(curRep);
		else {
			if (allDesc[curRep] == allDesc[curRep - 1])
				reps.push_back(curRep);
			else {
				desc << "Replicate";
				for (size_t i = 0; i < reps.size(); ++i)
					desc << " " << i;
				desc << ":\n" << allDesc[curRep - 1] << "\n";
				reps.clear();
				reps.push_back(curRep);
			}
		}
	}
	// reps should not be empty
	desc << "Replicate";
	for (size_t i = 0; i < reps.size(); ++i)
		desc << " " << i;
	desc << ":\n" << allDesc.back();
	return formatText(desc.str());
}


vectoru simulator::evolve(
                          const opList & initOps,
                          const opList & preOps,
                          const mating & matingScheme,
                          const opList & postOps,
                          const opList & finalOps,
                          int gens)
{
	if (numRep() == 0)
		return vectoru();

	// check compatibility of operators
	for (size_t i = 0; i < preOps.size(); ++i) {
		DBG_ASSERT(preOps[i]->isCompatible(*m_ptrRep[0]), ValueError,
			"Operator " + preOps[i]->describe() + " is not compatible.");
	}
	for (size_t i = 0; i < postOps.size(); ++i) {
		DBG_ASSERT(postOps[i]->isCompatible(*m_ptrRep[0]), ValueError,
			"Operator " + postOps[i]->describe() + " is not compatible.");
	}
	if (!matingScheme.isCompatible(pop(0)))
		throw ValueError("mating type is not compatible with current population settings.");

	vector<bool> activeReps(m_numRep);
	fill(activeReps.begin(), activeReps.end(), true);
	UINT numStopped = 0;

	// evolved generations, which will be returned.
	vectoru evolvedGens(m_numRep, 0U);

	// does not evolve.
	if (gens == 0)
		return evolvedGens;


	InitClock();

	// appy pre-op, most likely initializer. Do not check if they are active
	// or if they are successful
	if (!initOps.empty())
		apply(initOps);

	ElapsedTime("PreopDone");

    // make sure rep and gen exists in pop
    for (UINT curRep = 0; curRep < m_numRep; curRep++) {
        if (!m_ptrRep[curRep]->getVars().hasVar("gen"))
            m_ptrRep[curRep]->setGen(0);
        m_ptrRep[curRep]->setRep(curRep);
    }

	while (1) {
		// save refcount at the beginning
#ifdef Py_REF_DEBUG
		saveRefCount();
#endif

		for (UINT curRep = 0; curRep < m_numRep; curRep++) {
			population & curPop = *m_ptrRep[curRep];
            int curGen = curPop.gen();
	        int end = -1;
        	if (gens > 0)
		        end = curGen + gens - 1;
	        DBG_FAILIF(end < 0 && preOps.empty() && postOps.empty(), ValueError,
        		"Evolve with unspecified ending generation should have at least one terminator (operator)");

			DBG_ASSERT(static_cast<int>(curRep) == curPop.rep(), SystemError,
				"Replicate number does not match");

			if (!activeReps[curRep])
				continue;

			size_t it = 0;                                            // asign a value to reduce compiler warning

			if (PyErr_CheckSignals()) {
				cerr << "Evolution stopped due to keyboard interruption." << endl;
				fill(activeReps.begin(), activeReps.end(), false);
				numStopped = activeReps.size();
			}
			// apply pre-mating ops to current gen()
			if (!preOps.empty()) {
				for (it = 0; it < preOps.size(); ++it) {
					if (!preOps[it]->isActive(curRep, curGen, end, activeReps))
						continue;

					try {
						if (!preOps[it]->apply(curPop)) {
							DBG_DO(DBG_SIMULATOR, cerr << "Pre-mating Operator " + preOps[it]->describe() +
								" stops at replicate " + toStr(curRep) << endl);

							if (activeReps[curRep]) {
								numStopped++;
								activeReps[curRep] = false;
								break;
							}
						}
						if (PyErr_CheckSignals())
							throw StopEvolution("Evolution stopped due to keyboard interruption.");
					} catch (StopEvolution e) {
						DBG_DO(DBG_SIMULATOR, cerr	<< "All replicates are stopped due to a StopEvolution exception raised by "
							                        << "Pre-mating Operator " + preOps[it]->describe() +
							" stops at replicate " + toStr(curRep) << endl);
						if (e.message()[0] != '\0')
							cerr << e.message() << endl;
						fill(activeReps.begin(), activeReps.end(), false);
						numStopped = activeReps.size();
						break;
					}
					ElapsedTime("PreMatingOp: " + preOps[it]->describe());
				}
			}

			if (!activeReps[curRep])
				continue;
			// start mating:
			try {
				if (!const_cast<mating&>(matingScheme).mate(curPop, scratchpopulation())) {
					DBG_DO(DBG_SIMULATOR, cerr << "Mating stops at replicate " + toStr(curRep) << endl);

					numStopped++;
					activeReps[curRep] = false;
					// does not execute post-mating operator
					continue;
				}
				if (PyErr_CheckSignals())
					throw StopEvolution("Evolution stopped due to keyboard interruption.");
			} catch (StopEvolution e) {
				DBG_DO(DBG_SIMULATOR, cerr	<< "All replicates are stopped due to a StopEvolution exception raised by "
					                        << "During-mating Operator at replicate " + toStr(curRep) << endl);
				if (e.message()[0] != '\0')
					cerr << e.message() << endl;
				fill(activeReps.begin(), activeReps.end(), false);
				numStopped = activeReps.size();
				// does not execute post mating operator
				continue;
			}

			ElapsedTime("matingDone");

			// apply post-mating ops to next gen()
			if (!postOps.empty()) {
				for (it = 0; it < postOps.size(); ++it) {
					if (!postOps[it]->isActive(curRep, curGen, end, activeReps))
						continue;

					try {
						if (!postOps[it]->apply(curPop)) {
							DBG_DO(DBG_SIMULATOR, cerr << "Post-mating Operator " + postOps[it]->describe() +
								" stops at replicate " + toStr(curRep) << endl);
							numStopped++;
							activeReps[curRep] = false;
							// does not run the rest of the post-mating operators.
							break;
						}
						if (PyErr_CheckSignals())
							throw StopEvolution("Evolution stopped due to keyboard interruption.");
					} catch (StopEvolution e) {
						DBG_DO(DBG_SIMULATOR, cerr	<< "All replicates are stopped due to a StopEvolution exception raised by "
							                        << "Post-mating Operator " + postOps[it]->describe() +
							" stops at replicate " + toStr(curRep) << endl);
						if (e.message()[0] != '\0')
							cerr << e.message() << endl;
						fill(activeReps.begin(), activeReps.end(), false);
						numStopped = activeReps.size();
						// does not run the rest of the post-mating operators.
						break;
					}
					ElapsedTime("PostMatingOp: " + postOps[it]->describe());
				}
			}
			// if a replicate stops at a post mating operator, consider one evolved generation.
			++evolvedGens[curRep];
            curPop.setGen(curGen+1);
		}                                                                                       // each replicates

#ifdef Py_REF_DEBUG
		checkRefCount();
#endif

		--gens;
		//
		//   start 0, gen = 2
		//   0 -> 1 -> 2 stop (two generations)
		//
		//   step:
		//    cur, end = cur +1
		//    will go two generations.
		//  therefore, step should:
		if (numStopped == m_numRep || gens == 0)
			break;
	}                                                                                         // the big loop

	if (!finalOps.empty())
		apply(finalOps);

	// close every opened file (including append-cross-evolution ones)
	ostreamManager().closeAll();
	return evolvedGens;
}


bool simulator::apply(const opList & ops)
{
	for (size_t i = 0; i < ops.size(); ++i) {
		// check compatibility of operators
		DBG_ASSERT(ops[i]->isCompatible(*m_ptrRep[0]), ValueError,
			"Operator " + ops[i]->describe() + " is not compatible.");
	}

	// really apply
	for (UINT curRep = 0; curRep < m_numRep; curRep++) {
		population & curPop = *m_ptrRep[curRep];
		size_t it;

		// apply pre-mating ops to current gen
		for (it = 0; it < ops.size(); ++it) {
			vector<bool> activeReps(m_numRep);
			fill(activeReps.begin(), activeReps.end(), true);
			if (!ops[it]->isActive(curRep, 0, 0, activeReps, true))
				continue;

			ops[it]->apply(curPop);

			ElapsedTime("PrePost-preMatingop" + toStr(it));
		}
	}
	return true;
}


int simulator::__cmp__(const simulator & rhs) const
{
	if (numRep() != rhs.numRep())
		return 1;

	for (size_t i = 0; i < numRep(); ++i)
		if (pop(i).__cmp__(rhs.pop(i)) != 0)
			return 1;

	return 0;
}


}
