/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu                                                        *
 *                                                                         *
 *   $LastChangedDate: 2006-02-21 15:27:25 -0600 (Tue, 21 Feb 2006) $
 *   $Rev: 191 $
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

namespace io = boost::iostreams;

namespace simuPOP
{

	simulator::simulator(const population & pop,
		mating&  matingScheme,
		int rep , vectori grp)
		: m_gen(0), m_curRep(0), m_numRep(rep), m_groups(0),
		m_stopIfOneRepStop(false), m_applyOpToStoppedReps(false)
	{
		DBG_ASSERT(m_numRep >= 1, ValueError,
			"Number of replicates should be greater than or equal one.");

		DBG_DO(DBG_SIMULATOR, cout << "Creating simulator " << endl);

		m_matingScheme = matingScheme.clone();
		DBG_DO(DBG_SIMULATOR, cout << "mating scheme copied" << endl);

		if (!m_matingScheme->isCompatible(pop))
			throw TypeError
				("mating type is not compatible with current population settings.");

		// create replicates of given population
		m_ptrRep = new population *[m_numRep];
		m_scratchPop = new population(pop);

		DBG_FAILIF( m_scratchPop == NULL,
			SystemError, "Fail to create scratch population");

		try
		{
			for (UINT i = 0; i < m_numRep; ++i)
			{
				m_ptrRep[i] = new population(pop);

				DBG_FAILIF( m_ptrRep[i] == NULL,
					SystemError, "Fail to create new replicate");

				// set replication number
				m_ptrRep[i]->setRep(i);
			}
		}
		catch(OutOfMemory&)
		{
			cout << "Can not create " << m_numRep << " populations" << endl;
			throw OutOfMemory("Out of memory");
		}

		// use pop's geno structure
		this->setGenoStruIdx(m_ptrRep[0]->genoStruIdx());

		m_curRep = 0;

		setGroup(grp);

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

		delete[] m_ptrRep;

		delete m_matingScheme;

		DBG_DO(DBG_SIMULATOR, cout << "populations have been erased. " << endl
			<< "If you have ever used population() function to access some of the replicates, " << endl
			<< "these referenced population will not be working now." << endl);
	}

	simulator * simulator::clone() const
	{
		simulator * simu = new simulator(population(0),
			*m_matingScheme, m_numRep, m_groups);
		for(size_t i=0; i < m_numRep; ++i)
			simu->setPopulation(*m_ptrRep[i], i);
		return simu;
	}

	void simulator::setMatingScheme(const mating& matingScheme)
	{
		delete m_matingScheme;
		m_matingScheme = matingScheme.clone();
	}

	void simulator::setGroup(const vectori& grp)
	{
		m_groups = grp;

		if( m_groups.empty() )
		{
			m_groups.resize( m_numRep);
			for(UINT i=0; i < m_numRep; ++i)
				m_groups[i] = i;
		}

		DBG_ASSERT( m_groups.size() == m_numRep, ValueError,
			"Specified group index should have equal length to numRep()");

		for(UINT i=0; i < m_numRep; ++i)
			m_ptrRep[i]->setGrp(m_groups[i]);
	}

	bool simulator::evolve(const vectorop& ops,
		const vectorop& preOps,
		const vectorop& postOps,
		int end, bool dryrun)
	{
		DBG_DO(DBG_SIMULATOR, cout << "Starting generation: " << gen()
			<< " with ending generation " << end << endl);

		DBG_FAILIF(end > 0 && gen() > static_cast<UINT>(end), ValueError,
			"population gen is already greater than ending generation.");

		DBG_FAILIF( end < 0 && ops.empty(), ValueError,
			"Evolve with unspecified ending generation should have at least one terminator (operator)");

		vectorop preMatingOps, durmatingOps, postMatingOps, activeDurmatingOps;

		// an operator can belong to more than one groups.
		for (size_t i = 0; i < ops.size(); ++i)
		{
			if (ops[i]->canApplyPreMating())
				preMatingOps.push_back(ops[i]);
			if (ops[i]->canApplyPostMating())
				postMatingOps.push_back(ops[i]);
		}
		// first put in DuringMating operators that can form genotype
		for (size_t i = 0; i < ops.size(); ++i)
			if (ops[i]->canApplyDuringMating() && ops[i]->formOffGenotype())
				durmatingOps.push_back(ops[i]);

		if(durmatingOps.size() > 1)
			cout << "Warning: More than one during mating operators. Make sure they act on different replicates." << endl;

		for (size_t i = 0; i < ops.size(); ++i)
			if (ops[i]->canApplyDuringMating() && !ops[i]->formOffGenotype())
				durmatingOps.push_back(ops[i]);

		// check compatibility of operators
		for(size_t i = 0; i < ops.size(); ++i)
		{
			DBG_ASSERT(ops[i]->isCompatible(curpopulation()), ValueError,
				"Operator " + ops[i]->__repr__() + " is not compatible.");
		}

		InitClock();

		if( dryrun)
			cout << "Dryrun mode: display calling sequence" << endl;

		// appy pre-op, most likely initializer. Do not check if they are active
		// or if they are successful
		if( ! preOps.empty() )
		{
			if( dryrun )
				cout << "Apply pre-evolution operators" << endl;
			apply(preOps, dryrun);
		}

		ElapsedTime("PreopDone");

		vector < bool > stop(m_numRep);
		fill(stop.begin(), stop.end(), false);
		UINT numStopped = 0;

		// start the evolution loop
		if(dryrun)
			cout << "Start evolution:" << endl;

		while (1)
		{
			// starting a new gen
			if(! dryrun)
				setGen(m_gen);

			// save refcount at the beginning
#ifdef Py_REF_DEBUG
			saveRefCount();
#endif

			for (m_curRep = 0; m_curRep < m_numRep; m_curRep++)
			{
				DBG_ASSERT( static_cast<int>(m_curRep) == curpopulation().rep(), SystemError,
					"Replicate number does not match");

				if(dryrun)
					cout << "  Replicate " << m_curRep << endl;

				if (stop[m_curRep])
				{
					/// if apply to stopped reps, do it.
					if( !m_applyOpToStoppedReps)
						continue;
				}
				
				// set selection off so that all selector has to be preMating
				// that is to say, if some one set selection=True in a post mating opertor
				// it will have no effect
				DBG_FAILIF(curpopulation().hasVar("selection") && curpopulation().getVarAsBool("selection"),
					ValueError, "Selection is on from previous generation. Did you use PostMating selector?");

				curpopulation().setBoolVar("selection", false);

				size_t it = 0;					  // asign a value to reduce compiler warning

				// apply pre-mating ops to current gen()
				if (!preMatingOps.empty())
				{
					if(dryrun)
						cout << "    Pre-mating operators: " << endl;

					for (it = 0; it < preMatingOps.size(); ++it)
					{
						if( dryrun )
						{
							if(preMatingOps[it]->isActive(m_curRep, m_numRep, 0, 0, grp(), true) )
								cout << "      - " << preMatingOps[it]->__repr__() << preMatingOps[it]->atRepr() << endl;
							continue;
						}

						if( ! preMatingOps[it]->isActive(m_curRep, m_numRep, m_gen, end, grp()))
							continue;

						try
						{
							if (! preMatingOps[it]->applyWithScratch(curpopulation(), scratchpopulation(), PreMating))
							{
								DBG_DO(DBG_SIMULATOR, cout << "Pre-mating Operator " + preMatingOps[it]->__repr__() +
									" stops at replicate " + toStr(curRep()) << endl);

								if(! stop[m_curRep])
								{
									numStopped++;
									stop[m_curRep] = true;
								}
							}
						}
						catch(...)
						{
							cout << "PreMating operator " << preMatingOps[it]->__repr__() << " throws an exception." << endl << endl;
							throw;
						}
						ElapsedTime("PreMatingOp: " + preMatingOps[it]->__repr__());
					}
				}

				// start mating:
				// find out active during-mating operators
				if(dryrun)
					cout << "    Start mating." << endl;

				activeDurmatingOps.clear();
				for(vectorop::iterator op=durmatingOps.begin(), opEnd = durmatingOps.end();
					op != opEnd; ++op)
				{
					if(dryrun)
						cout << "      - " << (*op)->__repr__() << (*op)->atRepr() << endl;

					if( (*op)->isActive(m_curRep, m_numRep, m_gen, end, grp()))
						activeDurmatingOps.push_back( *op);
				}

				try
				{
					if (!dryrun && !m_matingScheme->mate(curpopulation(), scratchpopulation(), activeDurmatingOps, true))
					{
						DBG_DO(DBG_SIMULATOR, cout << "During-mating Operator stops at replicate "
							+ toStr(curRep()) << endl);

						if(! stop[m_curRep])
						{
							numStopped ++;
							stop[m_curRep] = true;
						}
					}
				}
				catch(...)
				{
					cout << "mating or one of the during mating operator throws an exception." << endl << endl;
					throw;
				}

				ElapsedTime("matingDone");

				// apply post-mating ops to next gen()
				if (!postMatingOps.empty())
				{
					if(dryrun )
						cout << "    Apply post-mating operators " << endl;

					for (it = 0; it < postMatingOps.size(); ++it)
					{

						if(dryrun)
						{
							if( postMatingOps[it]->isActive(m_curRep, m_numRep, 0, 0, grp(), true) )
								cout << "      - " << postMatingOps[it]->__repr__() << postMatingOps[it]->atRepr() << endl;
							continue;
						}

						if( ! postMatingOps[it]->isActive(m_curRep, m_numRep, m_gen, end, grp()))
							continue;

						try
						{
							if (! postMatingOps[it]-> applyWithScratch(curpopulation(), scratchpopulation(), PostMating))
							{
								DBG_DO(DBG_SIMULATOR, cout << "Post-mating Operator " + postMatingOps[it]->__repr__() +
									" stops at replicate " + toStr(curRep()) << endl);

								if(! stop[m_curRep])
								{
									numStopped ++;
									stop[m_curRep] = true;;
								}
							}
						}
						catch(...)
						{
							cout << "PostMating operator " << postMatingOps[it]->__repr__() << " throws an exception." << endl << endl;
							throw;
						}
						ElapsedTime("PostMatingOp: "+ postMatingOps[it]->__repr__());
					}
				}								  // post mating ops
			}									  // each replicates

			if(dryrun)
				break;

			// if one replicate stop and stopIfOneRepStop is set,
			// or if all replicates stop, or reach ending gen, stop iteration.
			// count the number of stopped replicates
			DBG_DO(DBG_SIMULATOR, cout << endl << "Number of stopped replicates: " << numStopped << endl);

#ifdef Py_REF_DEBUG
			checkRefCount();
#endif

			m_gen ++;							  // increase generation!
			// only if m_gen > end will simulation stop
			// that is to say, the ending generating will be executed
			//
			//   start 0, end = 1
			//   0 -> 1 -> 2 stop (two generations)
			//
			//   step:
			//    cur, end = cur +1
			//    will go two generations.
			//  therefore, step should:
			if ( (numStopped >= 1 && m_stopIfOneRepStop)
				|| numStopped >= m_numRep || (end >= 0
				&& gen() > static_cast<UINT>(end)))
				break;
		}										  // the big loop

		if( ! postOps.empty() )
		{
			if(dryrun)
				cout << "Apply post-evolution operators: " << endl;

			// finishing up, apply post-op
			apply(postOps, dryrun);
		}

		// close every opened file (including append-cross-evolution ones)
		ostreamManager().closeAll(true);
		return true;
	}

	bool simulator::apply(const vectorop ops, bool dryrun)
	{
		vectorop preMatingOps, postMatingOps;

		// an operator can belong to more than one groups.
		for (size_t i = 0; i < ops.size(); ++i)
		{
			if (ops[i]->canApplyPreMating())
				preMatingOps.push_back(ops[i]);
			if (ops[i]->canApplyDuringMating())
				throw TypeError("During-mating operator has to be called by a simulator.");
			if (ops[i]->canApplyPostMating())
				postMatingOps.push_back(ops[i]);
			// check compatibility of operators
			DBG_ASSERT(ops[i]->isCompatible(curpopulation()), ValueError,
				"Operator " + ops[i]->__repr__() + " is not compatible.");
		}

		// really apply
		for (m_curRep = 0; m_curRep < m_numRep; m_curRep++)
		{

			if(dryrun)
				cout << "  Replicate " << m_curRep << endl;

			size_t it;

			// apply pre-mating ops to current gen
			if (!preMatingOps.empty())
			{
				if(dryrun)
					cout << "    Apply pre-mating ops " << endl;

				for (it = 0; it < preMatingOps.size(); ++it)
				{

					if(dryrun)
					{
						cout << "      - " << preMatingOps[it]->__repr__() << preMatingOps[it]->atRepr() << endl;
						continue;
					}

					if( ! preMatingOps[it]->isActive(m_curRep, m_numRep, 0, 0, grp(), true))
						continue;

					preMatingOps[it] -> applyWithScratch(curpopulation(),
						scratchpopulation(), PreMating);

					ElapsedTime("PrePost-preMatingop"+toStr(it));

				}
			}

			// apply post-mating ops to next gen
			if (!postMatingOps.empty())
			{
				if(dryrun)
					cout << "    Applying post-mating ops " << endl;

				for (it = 0; it < postMatingOps.size(); ++it)
				{
					if(dryrun)
					{
						cout << "      - " << preMatingOps[it]->__repr__() << preMatingOps[it]->atRepr() << endl;
						continue;
					}

					if( ! postMatingOps[it]->isActive(m_curRep, m_numRep, 0, 0, grp(), true))
						continue;

					postMatingOps[it] -> applyWithScratch(curpopulation(), scratchpopulation(), PostMating);

					ElapsedTime("PrePost-postMatingop"+toStr(it));
				}
			}
		}
		return true;
	}

	void simulator::saveSimulator(string filename, string format, bool compress) const
	{
		io::filtering_ostream ofs;
		if(compress)
			ofs.push(io::gzip_compressor());
		ofs.push(io::file_sink(filename));

		if(!ofs)
			throw ValueError("Can not open file " + filename );

		if( format == "text"|| (format == "auto" && filename.substr(filename.size()-4,4) == ".txt" ))
		{
			boost::archive::text_oarchive oa(ofs);
			oa << *this;
		}
		else if (format == "xml" || (format == "auto" && filename.substr(filename.size()-4,4) == ".xml" ) )
		{
			boost::archive::xml_oarchive oa(ofs);
			oa << boost::serialization::make_nvp("simulator",*this);
		}
		else if (format == "bin"  ||  (format == "auto" && filename.substr(filename.size()-4,4) == ".bin" ))
		{
			boost::archive::binary_oarchive oa(ofs);
			oa << *this;
		}
		else
			throw ValueError("Wrong format type. Use one of text, xml, bin.");
	}

	void simulator::loadSimulator(string filename, string format)
	{
		io::filtering_istream ifs;
		if(isGzipped(filename))
			ifs.push(io::gzip_decompressor());
		ifs.push(io::file_source(filename));

		// do not need to test again
		if(!ifs)
			throw ValueError("Can not open file " + filename );

		try
		{
			if( format == "text" || (format == "auto" && filename.substr(filename.size()-4,4) == ".txt" ))
			{
				boost::archive::text_iarchive ia(ifs);
				ia >> *this;
			}
#ifndef __NO_XML_SUPPORT__
			else if( format == "xml"  ||  (format == "auto" && filename.substr(filename.size()-4,4) == ".xml" ))
			{
				boost::archive::xml_iarchive ia(ifs);
				ia >> boost::serialization::make_nvp("simulator", *this);
			}
#endif
			else if (format == "bin" || (format == "auto" && filename.substr(filename.size()-4,4) == ".bin" ) )
			{
				boost::archive::binary_iarchive ia(ifs);
				ia >> *this;
			}
			else
				throw;
		}
		catch(...)
		{
			// first close the file handle.

			DBG_DO(DBG_POPULATION,
				cout << "Can not determine file type, or file type is wrong. Trying different ways." << endl);

			// open a fresh ifstream
			io::filtering_istream ifbin;
			if(isGzipped(filename))
				ifbin.push(io::gzip_decompressor());
			ifbin.push(io::file_source(filename));

			// try to load the file using different iarchives.
			try									  // binary?
			{
				boost::archive::binary_iarchive ia(ifbin);
				ia >> *this;
			}
			catch(...)							  // not binary, text?
			{
				io::filtering_istream iftxt;
				if(isGzipped(filename))
					iftxt.push(io::gzip_decompressor());
				iftxt.push(io::file_source(filename));

				try
				{
					boost::archive::text_iarchive ia(iftxt);
					ia >> *this;
				}
				catch(...)						  // then xml?
				{
#ifndef __NO_XML_SUPPORT__
					io::filtering_istream ifxml;
					if(isGzipped(filename))
						ifxml.push(io::gzip_decompressor());
					ifxml.push(io::file_source(filename));
					try
					{
						boost::archive::xml_iarchive ia(ifxml);
						ia >> boost::serialization::make_nvp("simulator", *this);
					}
					catch(...)
					{
						throw ValueError("Failed to load simulator. Your file may be corrupted, "
							"or being a copy of non-transferrable file (.bin)");
					}
#else
					throw ValueError("Failed to load simulator. Your file may be corrupted, "
						"or being a copy of non-transferrable file (.bin)");
#endif
				}								  // try xml
			}									  // try text
		}										  // try bin
	}

	simulator& LoadSimulator(const string& file,
		mating& mate, string format)
	{
		population p;
		simulator * a = new simulator(
			p, mate );
#ifndef _NO_SERIALIZATION_
		a->loadSimulator(file, format);
		return *a;
#else
		cout << "This feature is not supported in this platform" << endl;
#endif
		return *a;
	}

}
