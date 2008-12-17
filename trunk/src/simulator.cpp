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

#include "simulator.h"

// for file compression
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>


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

	DBG_DO(DBG_SIMULATOR, cout << "Creating simulator " << endl);

	m_matingScheme = matingScheme.clone();
	DBG_DO(DBG_SIMULATOR, cout << "mating scheme copied" << endl);

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
		cout << "Can not create " << m_numRep << " populations" << endl;
		throw RuntimeError("Failed to create a population.");
	}

	// use pop's geno structure
	this->setGenoStruIdx(m_ptrRep[0]->genoStruIdx());

	// set generation number for all replicates
	setGen(0);

	DBG_DO(DBG_SIMULATOR, cout << "simulator created" << endl);
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
	GenoStruTrait(rhs),
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


void simulator::addInfoField(const string & field, double init)
{
	vectorstr newfields;

	try {
		infoIdx(field);
	} catch (IndexError &) {
		newfields.push_back(field);
	}

	if (!newfields.empty())
		setGenoStructure(struAddInfoFields(newfields));
	// all replicate
	for (UINT i = 0; i < m_numRep; ++i) {
		m_ptrRep[i]->addInfoField(field, init);
		DBG_ASSERT(genoStruIdx() == m_ptrRep[i]->genoStruIdx(),
			ValueError, "Genotypic structure of one of the "
			            "replicates does not agree with the structure of the simulator");
	}
	// and the scratch pop
	m_scratchPop->addInfoField(field, init);
	DBG_ASSERT(genoStruIdx() == m_scratchPop->genoStruIdx(),
		ValueError, "Genotypic structure of one of the "
		            "replicates does not agree with the structure of the simulator");
}


void simulator::addInfoFields(const vectorstr & fields, double init)
{
	vectorstr newfields;

	for (vectorstr::const_iterator it = fields.begin(); it != fields.end(); ++it) {
		try {
			infoIdx(*it);
		} catch (IndexError &) {
			newfields.push_back(*it);
		}
	}

	if (!newfields.empty())
		setGenoStructure(struAddInfoFields(newfields));

	for (UINT i = 0; i < m_numRep; ++i) {
		m_ptrRep[i]->addInfoFields(fields, init);
		DBG_ASSERT(genoStruIdx() == m_ptrRep[i]->genoStruIdx(),
			ValueError, "Genotypic structure of one of the "
			            "replicates does not agree with the structure of the simulator");
	}
	m_scratchPop->addInfoFields(fields, init);
	DBG_ASSERT(genoStruIdx() == m_scratchPop->genoStruIdx(),
		ValueError, "Genotypic structure of one of the "
		            "replicates does not agree with the structure of the simulator");
}


void simulator::setAncestralDepth(UINT depth)
{
	for (UINT i = 0; i < m_numRep; ++i)
		m_ptrRep[i]->setAncestralDepth(depth);
	m_scratchPop->setAncestralDepth(depth);
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


void simulator::setPopulation(population & pop, UINT rep)
{
	DBG_FAILIF(rep >= m_numRep, IndexError,
		"replicate index out of range. From 0 to numRep()-1 ");
	// get a copy of population
	delete m_ptrRep[rep];
	m_ptrRep[rep] = new population(pop);

	if (pop.genoStru() != this->genoStru() ) {
		DBG_DO(DBG_SIMULATOR,  cout << "Warning: added population has different genotypic structure." << endl);
		setGenoStruIdx(pop.genoStruIdx());
	}
	if (pop.genoStru() != m_scratchPop->genoStru() ) {
		delete m_scratchPop;
		m_scratchPop = new population(pop);
	}
	m_ptrRep[rep]->setRep(rep);
}


void simulator::setGen(ULONG gen)
{
	m_gen = gen;
	// set gen for all replicates
	for (UINT i = 0; i < m_numRep; ++i)
		m_ptrRep[i]->setGen(gen, true);
}


vectoru simulator::evolve(const vectorop & ops,
                          const vectorop & preOps,
                          const vectorop & postOps,
                          int gens, bool dryrun)
{
	vectorop preMatingOps, durmatingOps, postMatingOps, activeDurmatingOps;

	for (size_t i = 0; i < ops.size(); ++i) {
		if (ops[i]->canApplyPreMating())
			preMatingOps.push_back(ops[i]);
		if (ops[i]->canApplyPostMating())
			postMatingOps.push_back(ops[i]);
	}
	// first put in DuringMating operators that can form genotype
	for (size_t i = 0; i < ops.size(); ++i)
		if (ops[i]->canApplyDuringMating() && ops[i]->formOffGenotype())
			durmatingOps.push_back(ops[i]);

	for (size_t i = 0; i < ops.size(); ++i)
		if (ops[i]->canApplyDuringMating() && !ops[i]->formOffGenotype())
			durmatingOps.push_back(ops[i]);

	// check compatibility of operators
	for (size_t i = 0; i < ops.size(); ++i) {
		DBG_ASSERT(ops[i]->isCompatible(*m_ptrRep[0]), ValueError,
			"Operator " + ops[i]->__repr__() + " is not compatible.");
	}

	if (dryrun) {
		cout << "Dryrun mode: display calling sequence" << endl;
		if (!preOps.empty() ) {
			cout << "Apply pre-evolution operators" << endl;
			apply(preOps, true);
		}
		cout << "Start evolution" << endl;

		for (UINT curRep = 0; curRep < m_numRep; curRep++) {
			cout << "  Replicate " << curRep << endl;
			// apply pre-mating ops to current gen()
			if (!preMatingOps.empty()) {
				cout << "    Pre-mating operators" << endl;
				for (size_t it = 0; it < preMatingOps.size(); ++it)
					if (preMatingOps[it]->isActive(curRep, m_numRep, 0, 0, true) )
						cout << "      - " << preMatingOps[it]->__repr__() << preMatingOps[it]->atRepr() << endl;
			}
			cout << "    Start mating" << endl;
			for (vectorop::iterator op = durmatingOps.begin(), opEnd = durmatingOps.end();
			     op != opEnd; ++op)
				cout << "      - " << (*op)->__repr__() << (*op)->atRepr() << endl;
			// apply post-mating ops to next gen()
			if (!postMatingOps.empty()) {
				cout << "    Apply post-mating operators" << endl;
				for (size_t it = 0; it < postMatingOps.size(); ++it)
					if (postMatingOps[it]->isActive(curRep, m_numRep, 0, 0, true) )
						cout << "      - " << postMatingOps[it]->__repr__() << postMatingOps[it]->atRepr() << endl;
			}
		}
		if (!postOps.empty() ) {
			cout << "Apply post-evolution operators: " << endl;
			apply(postOps, true);
		}
		return vectoru(m_numRep, 0);
	}

	// it is possible that a user changes the internal population's
	// genotype strucutre. It is therefore necessary to check if
	// all populations have the same structure.
#ifndef OPTIMIZED
	for (size_t i = 0; i < m_numRep; ++i) {
		DBG_FAILIF(genoStruIdx() != m_ptrRep[i]->genoStruIdx(),
			ValueError, "Genotypic structure of one of the \n"
			            "replicates does not agree with the simulator. It is likely that you\n"
			            "have changed the genotypic structure of a population obtained from \n"
			            "simu::population(rep). This is not allowed.\n");
	}
	DBG_FAILIF(genoStruIdx() != m_scratchPop->genoStruIdx(),
		ValueError, "Genotypic structure of one of the \n"
		            "replicates does not agree with the simulator. It is likely that you\n"
		            "have changed the genotypic structure of a population obtained from \n"
		            "simu::population(rep). This is not allowed.\n");
#endif
	// evolved generations, which will be returned.
	vectoru evolvedGens(m_numRep, 0U);

	// does not evolve.
	if (gens == 0)
		return evolvedGens;

	int end = -1;
	if (gens > 0)
		end = gen() + gens - 1;

	DBG_FAILIF(end < 0 && ops.empty(), ValueError,
		"Evolve with unspecified ending generation should have at least one terminator (operator)");

	InitClock();

	// appy pre-op, most likely initializer. Do not check if they are active
	// or if they are successful
	if (!preOps.empty())
		apply(preOps, false);

	ElapsedTime("PreopDone");

	vector<bool> stopped(m_numRep);
	fill(stopped.begin(), stopped.end(), false);
	UINT numStopped = 0;

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

			if (stopped[curRep])
				continue;

			// set selection off so that all selector has to be preMating
			// that is to say, if some one set selection=True in a post mating opertor
			// it will have no effect
			DBG_FAILIF(curPop.hasVar("selection") && curPop.getVarAsBool("selection"),
				ValueError, "Selection is on from previous generation. Did you use PostMating selector?");

			curPop.setBoolVar("selection", false);

			size_t it = 0;                                            // asign a value to reduce compiler warning

			// apply pre-mating ops to current gen()
			if (!preMatingOps.empty()) {
				for (it = 0; it < preMatingOps.size(); ++it) {
					if (!preMatingOps[it]->isActive(curRep, m_numRep, m_gen, end))
						continue;

					try {
						if (!preMatingOps[it]->apply(curPop)) {
							DBG_DO(DBG_SIMULATOR, cout << "Pre-mating Operator " + preMatingOps[it]->__repr__() +
								" stops at replicate " + toStr(curRep) << endl);

							if (!stopped[curRep]) {
								numStopped++;
								stopped[curRep] = true;
								break;
							}
						}
					} catch (StopEvolution e) {
						DBG_DO(DBG_SIMULATOR, cout << "All replicates are stopped due to a StopEvolution exception raised by "
							                       << "Pre-mating Operator " + preMatingOps[it]->__repr__() +
							" stops at replicate " + toStr(curRep) << endl);
						if (e.message())
							cout << e.message() << endl;
						fill(stopped.begin(), stopped.end(), true);
						numStopped = stopped.size();
						break;
					} catch (...) {
						cout << "PreMating operator " << preMatingOps[it]->__repr__() << " throws an exception." << endl << endl;
						throw;
					}
					ElapsedTime("PreMatingOp: " + preMatingOps[it]->__repr__());
				}
			}

			if (stopped[curRep])
				continue;
			// start mating:
			// find out active during-mating operators
			activeDurmatingOps.clear();
			for (vectorop::iterator op = durmatingOps.begin(), opEnd = durmatingOps.end();
			     op != opEnd; ++op) {
				if ( (*op)->isActive(curRep, m_numRep, m_gen, end))
					activeDurmatingOps.push_back(*op);
			}

			try {
				if (!dryrun && !m_matingScheme->mate(curPop, scratchpopulation(), activeDurmatingOps, true)) {
					DBG_DO(DBG_SIMULATOR, cout << "During-mating Operator stops at replicate "
						+ toStr(curRep) << endl);

					numStopped++;
					stopped[curRep] = true;
					// does not execute post-mating operator
					continue;
				}
			} catch (StopEvolution e) {
				DBG_DO(DBG_SIMULATOR, cout << "All replicates are stopped due to a StopEvolution exception raised by "
					                       << "During-mating Operator at replicate " + toStr(curRep) << endl);
				if (e.message())
					cout << e.message() << endl;
				fill(stopped.begin(), stopped.end(), true);
				numStopped = stopped.size();
				// does not execute post mating operator
				continue;
			} catch (...) {
				cout << "mating or one of the during mating operator throws an exception." << endl << endl;
				throw;
			}

			ElapsedTime("matingDone");

			// apply post-mating ops to next gen()
			if (!postMatingOps.empty()) {
				for (it = 0; it < postMatingOps.size(); ++it) {
					if (!postMatingOps[it]->isActive(curRep, m_numRep, m_gen, end))
						continue;

					try {
						if (!postMatingOps[it]->apply(curPop)) {
							DBG_DO(DBG_SIMULATOR, cout << "Post-mating Operator " + postMatingOps[it]->__repr__() +
								" stops at replicate " + toStr(curRep) << endl);
							numStopped++;
							stopped[curRep] = true;
							// does not run the rest of the post-mating operators.
							break;
						}
					} catch (StopEvolution e) {
						DBG_DO(DBG_SIMULATOR, cout << "All replicates are stopped due to a StopEvolution exception raised by "
							                       << "Post-mating Operator " + postMatingOps[it]->__repr__() +
							" stops at replicate " + toStr(curRep) << endl);
						if (e.message())
							cout << e.message() << endl;
						fill(stopped.begin(), stopped.end(), true);
						numStopped = stopped.size();
						// does not run the rest of the post-mating operators.
						break;
					} catch (...) {
						cout << "PostMating operator " << postMatingOps[it]->__repr__() << " throws an exception." << endl << endl;
						throw;
					}
					ElapsedTime("PostMatingOp: " + postMatingOps[it]->__repr__());
				}
			}                                                                                   // post mating ops
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

	if (!postOps.empty())
		apply(postOps, false);

	// close every opened file (including append-cross-evolution ones)
	ostreamManager().closeAll();
	return evolvedGens;
}


bool simulator::apply(const vectorop ops, bool dryrun)
{
	for (size_t i = 0; i < ops.size(); ++i) {
		if (ops[i]->canApplyDuringMating())
			throw ValueError("During-mating operator has to be called by a simulator.");
		// check compatibility of operators
		DBG_ASSERT(ops[i]->isCompatible(*m_ptrRep[0]), ValueError,
			"Operator " + ops[i]->__repr__() + " is not compatible.");
	}

	// really apply
	for (UINT curRep = 0; curRep < m_numRep; curRep++) {
		if (dryrun)
			cout << "  Replicate " << curRep << endl;

		population & curPop = *m_ptrRep[curRep];
		size_t it;

		// apply pre-mating ops to current gen
		for (it = 0; it < ops.size(); ++it) {

			if (dryrun) {
				cout << "      - " << ops[it]->__repr__() << ops[it]->atRepr() << endl;
				continue;
			}

			if (!ops[it]->isActive(curRep, m_numRep, 0, 0, true))
				continue;

			ops[it]->apply(curPop);

			ElapsedTime("PrePost-preMatingop" + toStr(it));
		}
	}
	return true;
}


void simulator::save(string filename) const
{
	boost::iostreams::filtering_ostream ofs;

	ofs.push(boost::iostreams::gzip_compressor());
	ofs.push(boost::iostreams::file_sink(filename));

	if (!ofs)
		throw ValueError("Can not open file " + filename);

	boost::archive::text_oarchive oa(ofs);
	oa << *this;
}


void simulator::load(string filename)
{
	boost::iostreams::filtering_istream ifs;

	if (isGzipped(filename))
		ifs.push(boost::iostreams::gzip_decompressor());
	ifs.push(boost::iostreams::file_source(filename));

	// do not need to test again
	if (!ifs)
		throw ValueError("Can not open file " + filename);

	try {
		boost::archive::text_iarchive ia(ifs);
		ia >> *this;
	} catch (...) {
		throw ValueError("Failed to load simulator. Your file may be corrupted, "
			             "or being a copy of non-transferrable file (.bin)");
	}                                                                                               // try bin
}


simulator & LoadSimulator(const string & file, mating & matingScheme)
{
	population p;
	simulator * a = new simulator(p, matingScheme);

#ifndef _NO_SERIALIZATION_
	a->load(file);
	return *a;
#else
	cout << "This feature is not supported in this platform" << endl;
#endif
	return *a;
}


}
