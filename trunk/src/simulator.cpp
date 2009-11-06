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


simulator::simulator(const population & pop, mating & matingScheme, UINT rep)
	: m_gen(0), m_numRep(rep)
{
	DBG_ASSERT(m_numRep >= 1, ValueError,
		"Number of replicates should be greater than or equal one.");

	DBG_DO(DBG_SIMULATOR, cerr << "Creating simulator " << endl);

	m_matingScheme = matingScheme.clone();
	DBG_DO(DBG_SIMULATOR, cerr << "mating scheme copied" << endl);

	if (!m_matingScheme->isCompatible(pop))
		throw ValueError
			    ("mating type is not compatible with current population settings.");

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
	setGen(0);

	DBG_DO(DBG_SIMULATOR, cerr << "simulator created" << endl);
}


simulator::~simulator()
{
	// call the destructor of each replicates
	delete m_scratchPop;

	for (UINT i = 0; i < m_numRep; ++i)
		delete m_ptrRep[i];

	delete m_matingScheme;
}


simulator::simulator(const simulator & rhs) :
	m_gen(rhs.m_gen),
	m_matingScheme(NULL),
	m_numRep(rhs.m_numRep),
	m_ptrRep(0),
	m_scratchPop(NULL)
{
	m_matingScheme = rhs.m_matingScheme->clone();
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


void simulator::setMatingScheme(const mating & matingScheme)
{
	delete m_matingScheme;
	m_matingScheme = matingScheme.clone();
}


void simulator::add(const population & pop)
{
	++m_numRep;

	m_ptrRep.push_back(new population(pop));
	DBG_FAILIF(m_ptrRep.back() == NULL,
		RuntimeError, "Fail to add new population.");
}


void simulator::setGen(ULONG gen)
{
	m_gen = gen;
	// set gen for all replicates
	for (UINT i = 0; i < m_numRep; ++i)
		m_ptrRep[i]->setGen(gen, true);
}


string simulator::describe(const opList & initOps,
                           const opList & preOps,
                           const opList & duringOps,
                           const opList & postOps,
                           const opList & finalOps,
                           int gen)
{
	vectorstr allDesc(m_numRep, "");

	// assuming all active replicates.
	vector<bool> activeReps(m_numRep);

	for (UINT curRep = 0; curRep < m_numRep; curRep++) {
		ostringstream desc;

		if (initOps.empty())
			desc << "No operator is used to initialize population (initOps)." << endl;
		else {
			desc << "Apply pre-evolution operators to the initial population (initOps)." << endl;
			for (size_t it = 0; it < initOps.size(); ++it)
				desc << "    - " << initOps[it]->describe() << " " << initOps[it]->applicability(true, false) << endl;
		}
		if (gen < 0)
			desc << "\nEvolve a population indefinitely until an operator determines it." << endl;
		else
			desc << "\nEvolve a population for " << gen << " generations. "
			     << "(generations " << m_gen << " - " << m_gen + gen - 1
			     << ", stops at generation " << m_gen + gen << ")" << endl;
		if (preOps.empty())
			desc << "    No operator is applied to the parental generation (preOps)." << endl;
		else {
			desc << "    Apply pre-mating operators to the parental generation (preOps)" << endl;
			for (size_t it = 0; it < preOps.size(); ++it)
				if (preOps[it]->isActive(curRep, 0, 0, activeReps, true))
					desc << "    - " << preOps[it]->describe() << " " << preOps[it]->applicability() << endl;
		}
		desc << "\n    Populate an offspring populaton from the parental population." << endl
		     << "    " << m_matingScheme->describe() << endl;
		if (!duringOps.empty()) {
			desc << "    with additional during mating operators" << endl;
			for (size_t it = 0; it < duringOps.size(); ++it)
				desc << "        - " << duringOps[it]->describe() << " " << duringOps[it]->applicability() << endl;
		}
		//
		if (postOps.empty())
			desc << "\n    No operator is applied to the offspring population (postOps)." << endl;
		else {
			desc << "\n    Apply post-mating operators to the offspring population (postOps)." << endl;
			for (size_t it = 0; it < postOps.size(); ++it)
				if (postOps[it]->isActive(curRep, 0, 0, activeReps, true))
					desc << "    - " << postOps[it]->describe() << " " << postOps[it]->applicability() << endl;
		}
		if (finalOps.empty() )
			desc << "\nNo operator is applied to the final population (finalOps)." << endl;
		else {
			desc << "\nApply post-evolution operators (finalOps)" << endl;
			for (size_t it = 0; it < finalOps.size(); ++it)
				desc << "      - " << finalOps[it]->describe() << " " << finalOps[it]->applicability(true, false) << endl;
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
	// wrap the text .....
	string text = desc.str();
	string newtext;
	size_t newline = 0;
	size_t j = 0;
	int start = 0;
	int blank = 0;
	int lastblank = 0;
	bool continuation = false;
	while (true) {
nextline:
		lastblank = 0;
		if (!continuation) {
			for (j = 0; newline + j < text.size() && text[newline + j] == ' '; ++j) ;
			start = j;
			for (++j; newline + j < text.size() && text[newline + j] != ' '; ++j) ;
			blank = j + 1;
		}
		for (j = 0; newline + j < text.size() && text[newline + j] != '\n'; ++j) {
			if (j == 78) {
				newtext += text.substr(newline, lastblank) + "\n";
				if (start != 0)
					for (int i = 0; i < blank; ++i)
						newtext += " ";
				continuation = true;
				newline += lastblank + 1;
				goto nextline;
			}
			if (text[newline + j] == ' ')
				lastblank = j;
		}
		if (text[newline + j] == '\n' || newline + j >= text.size()) {
			newtext += text.substr(newline, j + 1);
			continuation = false;
			if (newline + j >= text.size())
				break;
			else {
				newline += j + 1;
				goto nextline;
			}
		}
	}
	return newtext;
}


vectoru simulator::evolve(
                          const opList & initOps,
                          const opList & preOps,
                          const opList & duringOps,
                          const opList & postOps,
                          const opList & finalOps,
                          int gens)
{
	// check compatibility of operators
	for (size_t i = 0; i < preOps.size(); ++i) {
		DBG_ASSERT(preOps[i]->isCompatible(*m_ptrRep[0]), ValueError,
			"Operator " + preOps[i]->describe() + " is not compatible.");
	}
	for (size_t i = 0; i < duringOps.size(); ++i) {
		DBG_ASSERT(duringOps[i]->isCompatible(*m_ptrRep[0]), ValueError,
			"Operator " + duringOps[i]->describe() + " is not compatible.");
	}
	for (size_t i = 0; i < postOps.size(); ++i) {
		DBG_ASSERT(postOps[i]->isCompatible(*m_ptrRep[0]), ValueError,
			"Operator " + postOps[i]->describe() + " is not compatible.");
	}

	vector<bool> activeReps(m_numRep);
	fill(activeReps.begin(), activeReps.end(), true);
	UINT numStopped = 0;

	// evolved generations, which will be returned.
	vectoru evolvedGens(m_numRep, 0U);

	// does not evolve.
	if (gens == 0)
		return evolvedGens;

	int end = -1;
	if (gens > 0)
		end = gen() + gens - 1;

	DBG_FAILIF(end < 0 && preOps.empty() && postOps.empty(), ValueError,
		"Evolve with unspecified ending generation should have at least one terminator (operator)");

	InitClock();

	// appy pre-op, most likely initializer. Do not check if they are active
	// or if they are successful
	if (!initOps.empty())
		apply(initOps);

	ElapsedTime("PreopDone");

	while (1) {
		// starting a new gen
		setGen(m_gen);

		// save refcount at the beginning
#ifdef Py_REF_DEBUG
		saveRefCount();
#endif

		for (UINT curRep = 0; curRep < m_numRep; curRep++) {
			population & curPop = *m_ptrRep[curRep];

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
					if (!preOps[it]->isActive(curRep, m_gen, end, activeReps))
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
						DBG_DO(DBG_SIMULATOR, cerr << "All replicates are stopped due to a StopEvolution exception raised by "
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
			// find out active during-mating operators
			vectorop activeDuringMatingOps;
			for (vectorop::const_iterator op = duringOps.begin(), opEnd = duringOps.end();
			     op != opEnd; ++op) {
				if ( (*op)->isActive(curRep, m_gen, end, activeReps))
					activeDuringMatingOps.push_back(*op);
			}

			try {
				if (!m_matingScheme->mate(curPop, scratchpopulation(), activeDuringMatingOps)) {
					DBG_DO(DBG_SIMULATOR, cerr << "During-mating Operator stops at replicate "
						+ toStr(curRep) << endl);

					numStopped++;
					activeReps[curRep] = false;
					// does not execute post-mating operator
					continue;
				}
				if (PyErr_CheckSignals())
					throw StopEvolution("Evolution stopped due to keyboard interruption.");
			} catch (StopEvolution e) {
				DBG_DO(DBG_SIMULATOR, cerr << "All replicates are stopped due to a StopEvolution exception raised by "
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
					if (!postOps[it]->isActive(curRep, m_gen, end, activeReps))
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
						DBG_DO(DBG_SIMULATOR, cerr << "All replicates are stopped due to a StopEvolution exception raised by "
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
		}                                                                                       // each replicates

#ifdef Py_REF_DEBUG
		checkRefCount();
#endif

		++m_gen;                                                                  // increase generation!
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


void simulator::save(string filename) const
{
	boost::iostreams::filtering_ostream ofs;

	ofs.push(boost::iostreams::gzip_compressor());
	ofs.push(boost::iostreams::file_sink(filename, std::ios::binary));

	if (!ofs)
		throw RuntimeError("Can not open file " + filename);

	boost::archive::text_oarchive oa(ofs);
	oa << *this;
}


void simulator::load(string filename)
{
	boost::iostreams::filtering_istream ifs;

	ifs.push(boost::iostreams::gzip_decompressor());
	ifs.push(boost::iostreams::file_source(filename, std::ios::binary));

	// do not need to test again
	if (!ifs)
		throw RuntimeError("Can not open file " + filename);

	try {
		boost::archive::text_iarchive ia(ifs);
		ia >> *this;
	} catch (...) {
		throw RuntimeError("Failed to load simulator. Your file may be corrupted.");
	}
}


simulator & LoadSimulator(const string & file, mating & matingScheme)
{
	population p;
	simulator * a = new simulator(p, matingScheme);

	a->load(file);
	return *a;
}


}
