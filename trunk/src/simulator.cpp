/**
 *  $File: Simulator.cpp $
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

#include "simulator.h"

#include <sstream>
using std::ostringstream;

namespace simuPOP {

Population & pyPopIterator::next()
{
	if (m_index == m_end)
		throw StopIteration("");
	else
		return **(m_index++);
}


Simulator::Simulator(PyObject * pops, UINT rep, bool steal)
{
	PARAM_ASSERT(rep >= 1, ValueError,
		"Number of replicates should be greater than or equal one.");

	DBG_DO(DBG_SIMULATOR, cerr << "Creating Simulator " << endl);
	m_pops = vector<Population *>();
	// create replicates of given Population
	m_scratchPop = new Population();

	if (pops == NULL) {
		return;
	} else if (PySequence_Check(pops)) {
		size_t size = PySequence_Size(pops);
		for (size_t i = 0; i < size; ++i) {
			PyObject * item = PySequence_GetItem(pops, i);
			void * pop = pyPopPointer(item);
			PARAM_ASSERT(pop, ValueError, "Parameter pops should be a single population or a list populations.");
			if (steal) {
				Population * tmp = new Population();
				tmp->swap(*reinterpret_cast<Population *>(pop));
				m_pops.push_back(tmp);
			} else {
				try {
					m_pops.push_back(reinterpret_cast<Population *>(pop)->clone());
					DBG_FAILIF(m_pops.back() == NULL,
						SystemError, "Fail to create new replicate");
				} catch (...) {
					throw RuntimeError("Failed to create a population.");
				}
			}
			Py_DECREF(item);
		}
	} else {
		void * pop = pyPopPointer(pops);
		PARAM_ASSERT(pop, ValueError, "Parameter pops should be a single population or a list populations.");
		if (steal) {
			Population * tmp = new Population();
			tmp->swap(*reinterpret_cast<Population *>(pop));
			m_pops.push_back(tmp);
		} else {
			try {
				m_pops.push_back(reinterpret_cast<Population *>(pop)->clone());
				DBG_FAILIF(m_pops.back() == NULL,
					SystemError, "Fail to create new replicate");
			} catch (...) {
				throw RuntimeError("Failed to create a population.");
			}
		}
	}
	// parameter rep
	size_t numRep = m_pops.size();
	for (UINT i = 1; i < rep; ++i) {
		for (UINT j = 0; j < numRep; ++j) {
			try {
				m_pops.push_back(m_pops[j]->clone());
				DBG_FAILIF(m_pops.back() == NULL,
					SystemError, "Fail to create new replicate.");
			} catch (...) {
				throw RuntimeError("Failed to create a population.");
			}
		}
	}
	// set var "rep"
	for (UINT i = 0; i < m_pops.size(); ++i)
		m_pops[i]->setRep(i);

	DBG_FAILIF(m_scratchPop == NULL,
		SystemError, "Fail to create scratch population");

	// set generation number for all replicates
	DBG_DO(DBG_SIMULATOR, cerr << "Simulator created" << endl);
}


Simulator::~Simulator()
{
	// call the destructor of each replicates
	delete m_scratchPop;

	for (UINT i = 0; i < m_pops.size(); ++i)
		delete m_pops[i];
}


Simulator::Simulator(const Simulator & rhs) :
	m_pops(0),
	m_scratchPop(NULL)
{
	m_scratchPop = rhs.m_scratchPop->clone();
	m_pops = vector<Population *>(rhs.m_pops.size());
	for (size_t i = 0; i < m_pops.size(); ++i) {
		m_pops[i] = rhs.m_pops[i]->clone();
		m_pops[i]->setRep(i);
	}
}


Simulator * Simulator::clone() const
{
	return new Simulator(*this);
}


Population & Simulator::population(size_t rep) const
{
	DBG_WARNIF(true, "The returned object of function Simulator.population is a temporary reference "
		             "to a population inside a Simulator. It will become invalid once the simulator "
		             "changes. Please use Simulator.extract or Simulator.population(rep).clone() if"
		             "would like to get an independent population.");

	DBG_FAILIF(rep >= m_pops.size(), IndexError,
		"replicate index out of range. From 0 to numRep()-1 ");

	return *m_pops[rep];
}


Population & Simulator::extract(UINT rep)
{
	DBG_FAILIF(rep >= m_pops.size(), IndexError,
		"replicate index out of range. From 0 to numRep()-1 ");

	Population * pop = m_pops[rep];
	m_pops.erase(m_pops.begin() + rep);
	return *pop;
}


void Simulator::add(const Population & pop, bool steal)
{

	if (steal) {
		Population * tmp = new Population();
		const_cast<Population &>(pop).swap(*tmp);
		m_pops.push_back(tmp);
	} else
		m_pops.push_back(new Population(pop));
	PARAM_FAILIF(m_pops.back() == NULL,
		RuntimeError, "Fail to add new Population.");
}


string Simulator::describe(bool /* format */) const
{
	return (boost::format("<simuPOP.Simulator> a simulator with %1% populations") % m_pops.size()).str();
}


vectoru Simulator::evolve(
                          const opList & initOps,
                          const opList & preOps,
                          const MatingScheme & matingScheme,
                          const opList & postOps,
                          const opList & finalOps,
                          int gens, bool dryrun)
{
	if (dryrun) {
		cerr << describeEvolProcess(initOps, preOps, matingScheme, postOps, finalOps, gens, numRep()) << endl;
		return vectoru(numRep());
	}

	if (numRep() == 0)
		return vectoru();

	vector<bool> activeReps(m_pops.size());
	fill(activeReps.begin(), activeReps.end(), true);
	size_t numStopped = 0;

	// evolved generations, which will be returned.
	vectoru evolvedGens(m_pops.size(), 0U);

	// does not evolve.
	if (gens == 0)
		return evolvedGens;

	// make sure rep and gen exists in pop
	for (UINT curRep = 0; curRep < m_pops.size(); curRep++) {
		if (!m_pops[curRep]->getVars().hasVar("gen"))
			m_pops[curRep]->setGen(0);
		m_pops[curRep]->setRep(curRep);
	}

	initClock();

	// appy pre-op, most likely initializer. Do not check if they are active
	// or if they are successful
	if (!initOps.empty())
		apply(initOps);

	elapsedTime("Start evolution.");

	while (1) {
		// save refcount at the beginning
#ifdef Py_REF_DEBUG
		saveRefCount();
#endif

		for (size_t curRep = 0; curRep < m_pops.size(); curRep++) {
			Population & curPop = *m_pops[curRep];
			// sync population variable gen with gen(). This allows
			// users to set population variable to change generation number.
			long curGen = curPop.getVars().getVarAsInt("gen");
			if (curGen != static_cast<long>(curPop.gen()))
				curPop.setGen(curGen);

			ssize_t end = -1;
			if (gens > 0)
				end = curGen + gens - 1;
			//PARAM_FAILIF(end < 0 && preOps.empty() && postOps.empty(), ValueError,
			//	"Evolve with unspecified ending generation should have at least one terminator (operator)");

			DBG_ASSERT(curRep == curPop.rep(), SystemError,
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
							DBG_DO(DBG_SIMULATOR, cerr << "Pre-mating Operator " << preOps[it]->describe() <<
								" stops at replicate " << curRep << endl);

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
							                        << "Pre-mating Operator " << preOps[it]->describe() <<
							" stops at replicate " << curRep << endl);
						if (e.message()[0] != '\0')
							cerr << e.message() << endl;
						fill(activeReps.begin(), activeReps.end(), false);
						numStopped = activeReps.size();
						break;
					} catch (RevertEvolution e) {
						long newCurGen = curPop.getVars().getVarAsInt("gen");
						if (newCurGen != static_cast<long>(curPop.gen()))
							curPop.setGen(newCurGen);
						if (gens > 0)
							gens += curGen - newCurGen;
						curGen = newCurGen;
						DBG_DO(DBG_SIMULATOR, cerr << "Revert to generation " << curGen << endl);
					}

					elapsedTime("Applied " + preOps[it]->describe());
				}
			}

			if (!activeReps[curRep])
				continue;
			elapsedTime((boost::format("Start mating at generation %1%") % curGen).str());
			// start mating:
			try {
				if (!const_cast<MatingScheme &>(matingScheme).mate(curPop, scratchPopulation())) {
					DBG_DO(DBG_SIMULATOR, cerr << "Mating stops at replicate " << curRep << endl);

					numStopped++;
					activeReps[curRep] = false;
					// does not execute post-mating operator
					continue;
				}
				if (PyErr_CheckSignals())
					throw StopEvolution("Evolution stopped due to keyboard interruption.");
			} catch (StopEvolution e) {
				DBG_DO(DBG_SIMULATOR, cerr	<< "All replicates are stopped due to a StopEvolution exception raised by "
					                        << "During-mating Operator at replicate " << curRep << endl);

				fill(activeReps.begin(), activeReps.end(), false);
				numStopped = activeReps.size();
				// does not execute post mating operator
				break;
			} catch (RevertEvolution e) {
				long newCurGen = curPop.getVars().getVarAsInt("gen");
				if (newCurGen != static_cast<long>(curPop.gen()))
					curPop.setGen(newCurGen);
				if (gens > 0)
					gens += curGen - newCurGen;
				curGen = newCurGen;
				DBG_DO(DBG_SIMULATOR, cerr << "Revert to generation " << curGen << endl);
			}

			elapsedTime("Mating finished.");

			// apply post-mating ops to next gen()
			if (!postOps.empty()) {
				for (it = 0; it < postOps.size(); ++it) {
					if (!postOps[it]->isActive(curRep, curGen, end, activeReps))
						continue;

					try {
						if (!postOps[it]->apply(curPop)) {
							DBG_DO(DBG_SIMULATOR, cerr << "Post-mating Operator " + postOps[it]->describe() +
								" stops at replicate " << curRep << endl);
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
							" stops at replicate " << curRep << endl);
						if (e.message()[0] != '\0')
							cerr << e.message() << endl;
						fill(activeReps.begin(), activeReps.end(), false);
						numStopped = activeReps.size();
						// does not run the rest of the post-mating operators.
						break;
					} catch (RevertEvolution e) {
						long newCurGen = curPop.getVars().getVarAsInt("gen");
						if (newCurGen != static_cast<long>(curPop.gen()))
							curPop.setGen(newCurGen);
						if (gens > 0)
							gens += curGen - newCurGen;
						curGen = newCurGen;
						DBG_DO(DBG_SIMULATOR, cerr << "Revert to generation " << curGen << endl);
					}
					elapsedTime("Applied " + postOps[it]->describe());
				}
			}
			// if a replicate stops at a post mating operator, consider one evolved generation.
			++evolvedGens[curRep];
			curPop.setGen(curGen + 1);
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
		if (numStopped == m_pops.size() || gens == 0)
			break;
	}                                                                                         // the big loop

	if (!finalOps.empty())
		apply(finalOps);

	// close every opened file (including append-cross-evolution ones)
	ostreamManager().closeAll();
	cleanupCircularRefs();
	return evolvedGens;
}


bool Simulator::apply(const opList & ops)
{
	// really apply
	for (UINT curRep = 0; curRep < m_pops.size(); curRep++) {
		Population & curPop = *m_pops[curRep];
		size_t it;
		// apply pre-mating ops to current gen
		for (it = 0; it < ops.size(); ++it) {
			vector<bool> activeReps(m_pops.size());
			fill(activeReps.begin(), activeReps.end(), true);
			if (!ops[it]->isActive(curRep, 0, 0, activeReps, true))
				continue;

			try {
				ops[it]->apply(curPop);
			} catch (RevertEvolution e) {
				//
			}
			elapsedTime("Applied " + ops[it]->describe());
		}
	}
	return true;
}


int Simulator::__cmp__(const Simulator & rhs) const
{
	if (numRep() != rhs.numRep())
		return 1;

	for (size_t i = 0; i < numRep(); ++i)
		if (population(i).__cmp__(rhs.population(i)) != 0)
			return 1;

	return 0;
}


string describeEvolProcess(const opList & initOps,
                           const opList & preOps,
                           const MatingScheme & matingScheme,
                           const opList & postOps,
                           const opList & finalOps,
                           int gen, size_t numRep)
{
	vectorstr allDesc(numRep, "");

	// assuming all active replicates.
	vector<bool> activeReps(numRep);

	for (UINT curRep = 0; curRep < numRep; curRep++) {
		ostringstream desc;

		if (initOps.empty())
			desc << "No operator is used to initialize Population (initOps).\n";
		else {
			desc << "Apply pre-evolution operators to the initial population (initOps).\n<ul>\n";
			for (size_t it = 0; it < initOps.size(); ++it)
				desc << "<li>" << initOps[it]->describe(false) << " " << initOps[it]->applicability(true, false) << endl;
			desc << "</ul>\n";
		}
		if (gen < 0)
			desc << "\nEvolve a population indefinitely until an operator determines it." << endl;
		else
			desc << "\nEvolve a population for " << gen << " generations" << endl;
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
		if (finalOps.empty())
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
	for (UINT curRep = 0; curRep < numRep; curRep++) {
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
	return formatDescription(desc.str());
}


}
