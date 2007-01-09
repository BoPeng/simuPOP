/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu
 *                                                                         *
 *   $LastChangedDate$
 *   $Rev$                                                            *
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

#ifndef _SIMULATOR_H
#define _SIMULATOR_H

/**
\file
\brief head file of class simulator and some other utility functions
*/

#include "utility.h"
#include "simuPOP_cfg.h"
#include "mating.h"
#include "operator.h"

#include <vector>
#include <algorithm>

using std::vector;
using std::fill;
using std::swap;

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
using boost::serialization::make_nvp;

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>

#ifndef OPTIMIZED
#include <time.h>								  // for clock() function

#define InitClock(); \
if(debug(DBG_PROFILE))  m_clock = clock();

#define ElapsedTime(name); \
if(debug(DBG_PROFILE)) \
{ \
	cout << name << ": " << static_cast<double>( clock() - m_clock )/CLOCKS_PER_SEC << "\n"; \
	m_clock = clock(); \
}
#else
#define InitClock();
#define ElapsedTime(name);
#endif

/** \brief all classes in simuPOP is defined in this namespace
 */
namespace simuPOP
{

	/**
	\brief simulator manage several replicates of a population,
	  evolve them using given mating scheme and operators.

	simulators combine three important components of simuPOP:
	population, mating scheme and operators together. A
	simulator is created with an instance of population, a
	replicate number and a mating scheme. It makes 'rep'
	replicates of this population and control the evolution
	process of these populations.

	The most important functions of a simulator is of course
	evolve(). It accepts arrays of operators as its parameters,
	among which, 'preOps' and 'postOps' will be applied to the
	populations at the begining and end of evolution, whereas
	'operators'will be applied at every generation.

	simulators separates operators into pre-, during- and post-
	mating operators. During evolution, simulator first
	apply all pre-mating operators and then call the mate()
	function of the given mating scheme, which will call
	during-mating operators during the birth of each offsrping.
	After the mating is finished, post-mating operators are
	applied in the order they apprear in the operator list.

	Since operators can apply to specific replicate or replicates
	group, and might not be active at all time, the isActive(m_curRep, m_numRep, m_gen, end, grp())
	function of each operator is called before it is applied
	to the populations.

	simulators can evolve a given number of generations (the
	'end' parameter of evolve), or evolve indefinitely using
	a certain type of operators called terminators. In this
	case, one or more terminators will check the status of
	evolution and determine if the simulation should be stopped.
	An obvious example of such a terminator is a fixation-checker.

	Finally, a simulator can be saved to a file in the format
	of 'text', 'bin', or 'xml'. This enables us to stop a
	simulation and ressume it at another time or on another
	machine. It is also a good idea to save a snapshot of a
	simulation every several generations.

	DEVONLY{This is a template class that manage the whole simulation process. It glues
	three major components of simuPOP together: population, mating, Operation.
	More specifically, a simuPOP<population> object manages several copies of
	population (or one of its subclasses) plus a 'scratch population'; it has
	a mating object that knows how to generate next gen; it controls
	the evolution process by applying pre- during- and post- mating Operators
	during evolution. }

	*/
	class simulator : public GenoStruTrait
	{
		public:
			/// vector of operator pointers.
			typedef std::vector< Operator * > vectorop;

		public:

			/// create a simulator
			/**
			\param population a population created by population()
			   function. This population will be copied to the simulator
			   so its content will not be changed.
			\param mate a mating scheme
			\param rep number of replicates. default to 1
			\param grp grp number for each replicate. For example,
			  we can seprate all replicates into two groups and
			  give them differnt initial values before evolution.
			\return none
			\sa population, mating

			DEVONLY{ m_curRep, gen  are reference to
			glocal shared variables. }
			*/
			simulator(const population & pop, mating&  matingScheme,
				int rep = 1, vectori grp=vectori());

			/// destroy a simulator along with all its populations
			/**
			\note pop = simulator::population()
			returns temporary reference to an internal population.
			After a simulator evolves another genertion or
			after the simulator is destroyed, this referenced
			population should *not* be used.
			*/
			~simulator();

			simulator * clone() const;

			/// the 'rep' replicate of this simulator
			/**
			\param rep number of replicate.
			\return reference to population rep.
			\note  The returned reference is temporary in the sense that
			  the refered population will be invalid after another round
			  of evolution. Therefore, the use of this function should
			  be limited to 'immediate after retrival'.
			  If you would like to get a persistent population, use getPopulation(rep)
			*/
			population& pop(UINT rep) const
			{
				DBG_FAILIF( rep >= m_numRep, IndexError,
					"replicate index out of range. From 0 to numRep()-1 ");

				return *m_ptrRep[rep];
			}

			/**
			this function returns a copy of population rep
			\param rep number of replicate.
			\return reference to a population.
			*/
			population& getPopulation(UINT rep)
			{
				return *new population( *m_ptrRep[rep]);
			}

			// set mating scheme
			void setMatingScheme(const mating& matingScheme);

			// set population... unsafe!
			void setPopulation(population& pop, UINT rep)
			{
				DBG_FAILIF( rep >= m_numRep, IndexError,
					"replicate index out of range. From 0 to numRep()-1 ");
				// get a copy of population
				delete m_ptrRep[rep];
				m_ptrRep[rep] = new population(pop);

				if( pop.genoStru() != this->genoStru() )
				{
					DBG_DO(DBG_SIMULATOR,  cout << "Warning: added population has different genotypic structure." << endl);
					setGenoStruIdx( pop.genoStruIdx());
				}
				if( pop.genoStru() != m_scratchPop->genoStru() )
				{
					delete m_scratchPop;
					m_scratchPop = new population(pop);
				}
				m_ptrRep[rep]->setRep(rep);
			}

			/// CPPONLY
			UINT curRep() const
			{
				return m_curRep;
			}

			UINT numRep() const
			{
				return m_numRep;
			}

			ULONG gen() const
			{
				return m_gen;
			}

			/// CPPONLY
			int grp()
			{
				return m_groups[m_curRep];
			}

			/// return group indices
			vectori group()
			{
				return m_groups;
			}

			/// set groups for replicates
			void setGroup(const vectori& grp);

			/// set generation number
			/**
			\param gen new generation number
			\note this will also set shared variable %gen
			*/
			void setGen(ULONG gen)
			{
				m_gen = gen;
				// set gen for all replicates
				for(UINT i=0; i < m_numRep; ++i)
					m_ptrRep[i]->setGen(gen, true);
			}

			/// evolve one step
			/**
			\sa simulator::evolve()
			*/
			bool step(const vectorop& ops = vectorop(),
				const vectorop& preOps = vectorop(),
				const vectorop& postOps = vectorop(), UINT steps=1,
				bool dryrun=false)
			{
				// since end gen itself will be executed
				// cur = 5,
				// end = 6
				// will go two steps.
				return evolve( ops, preOps, postOps, gen() + steps - 1, dryrun);
			}

			/// evolve till 'end' generation subject to given operators
			/**
			\param ops operators that will be applied at all generations.
			Of course they might not be active at all generations.
			\param preOps operators that will be applied before evolution.
			evolve() function will *not* check if they are active.
			\param postOps operators that will be applied after the last
			generation, or after a replicate is terminated.
			\param end ending generation. Should be -1 (no ending generation)
			or a number greater than current generation number. When end=-1,
			simulator can only be stopped by terminators.
			\result True if evolution finishs successfully.

			'operators' will be applied in the order of:

			*  all pre-mating opertors

			*  during-mating operators called by the mating scheme at the
			birth of each offspring

			*  all post-mating operators

			If any pre or post-mating opertor fails to apply, that
			replicate will be stopped. This is exactly how terminators
			work.

			\sa simulator::step()

			\note When end=-1, you can not specify negative generation
			parameters to operators. How would an operator know which
			genertion is the -1 genertion if no ending genertion is
			given?
			*/
			bool evolve(const vectorop& ops,
				const vectorop& preOps = vectorop(),
				const vectorop& postOps = vectorop(),
				int end =  -1, bool dryrun = false);

			/// apply some ops, geneartion of population does not change
			/// No mating is allowed.
			/**
			\param ops operators that will be applied at all generations.
			  Of course they might not be active at all generations.
			\result True if evolution finishs successfully.

			pre-mating oeprators are applied before post-mating operators.
			no during-mating operators are allowed.
			*/
			bool apply(const vectorop ops, bool dryrun=false);

			/// stop if one replicate stops or not
			/**
			\param on turn on or off stopIfOneRepStop

			if set, the simulator will stop evolution if one replicate
			stops. This is sometimes useful.
			*/
			bool setStopIfOneRepStop(bool on=true)
			{
				m_stopIfOneRepStop = on;
				return(on);
			}

			bool stopIfOneRepStop()
			{
				return(m_stopIfOneRepStop);
			}

			/// apply ops even if rep stops
			/**
			\param on turn on or off applyOpToStoppedReps flag

			if set, the simulator will continue to apply operators
			to all stopped repicates until all replicates are marked
			stopped. This is sometimes useful.
			*/
			bool setApplyOpToStoppedReps(bool on = true)
			{
				m_applyOpToStoppedReps = on;
				return(on);
			}

			bool applyOpToStoppedReps()
			{
				return(m_applyOpToStoppedReps);
			}

			/// get simulator namespace, if rep > 0 is given,
			/// return replicate rep namespace
			PyObject* vars(UINT rep, int subPop=-1)
			{
				if( static_cast<UINT>(rep) >= m_numRep )
					throw ValueError("Replicate index out of range.");

				return m_ptrRep[rep]->vars(subPop);
			}

			/// save simulator in 'text','bin' or 'xml' format
			/**
			\param filename    save to filename
			\param format format to save. Can be one of 'text', 'bin', 'xml'

			The default format is 'text' but the output is not suppored to be read.
			'bin' has smaller size and should be used for large populations.
			'xml' format is most readable and should be used when you would
			like to convert simuPOP populations to other formats.

			\sa global function loadsimulator
			*/
			void saveSimulator(string filename, string format="auto", bool compress=true) const;

			/// CPPONLY load simulator from a file
			/**
			\param filename load from filename
			\param format format to load. Can be one of "text", "bin", "xml". It should match
			  the format used to save the simulator.

			\sa saveSimulator
			*/
			void loadSimulator(string filename, string format="auto");

			// allow str(population) to get something better looking
			string __repr__()
			{
				return "<simulator with " + toStr(numRep()) + " replicates>";
			}

		private:

			friend class boost::serialization::access;

			template<class Archive>
				void save(Archive &ar, const UINT version) const
			{
				ULONG l_gen = gen();
				ar & make_nvp("generation", l_gen);
				ar & make_nvp("num_replicates", m_numRep);
				// ignore scratch population
				for( UINT i=0; i< m_numRep; i++)
					ar & make_nvp("populations", *m_ptrRep[i]);
				ar & make_nvp("groups", m_groups );
			}

			template<class Archive>
				void load(Archive &ar, const UINT version)
			{
				m_curRep = 1;
				ULONG l_gen;
				ar & make_nvp("generation", l_gen);

				ar & make_nvp("num_replicates", m_numRep);

				m_ptrRep = new population *[m_numRep];

				for (UINT i = 0; i < m_numRep; ++i)
				{
					m_ptrRep[i] = new population();
					ar & make_nvp("populations", *m_ptrRep[i]);
					m_ptrRep[i]->setRep(i);
				}
				m_scratchPop = new population(*m_ptrRep[0]);
				setGenoStruIdx(m_ptrRep[0]->genoStruIdx());

				ar & make_nvp("groups", m_groups);

				m_stopIfOneRepStop = false;
				m_applyOpToStoppedReps = false;

				// only after all replicates are ready do we set gen
				setGen(l_gen);
			}

			BOOST_SERIALIZATION_SPLIT_MEMBER();

		private:

			/// access current population
			population& curpopulation()
			{
				return *m_ptrRep[m_curRep];
			}

			/// access scratch population
			population&  scratchpopulation()
			{
				return *m_scratchPop;
			}

		private:

			/// current generation
			ULONG m_gen;

			/// index to current population, refernce to shared "%rep"
			UINT m_curRep;

			/// mating functor that will do the mating.
			mating* m_matingScheme;

			/// number of replicates of population
			UINT m_numRep;

			/// groups
			vectori m_groups;

			/// replicate pointers
			population **m_ptrRep;

			/// the scratch pop
			population *m_scratchPop;

			/// determine various behaviors of simulator
			/// for example: if stop after one replicate stops
			///   or stop after all stop
			///   or if apply operators if stopped.
			bool m_stopIfOneRepStop;

			bool m_applyOpToStoppedReps;

#ifndef OPTIMIZED
			clock_t m_clock;
#endif

	};

	simulator& LoadSimulator(const string& file,
		mating& mate,
		string format="auto");

}
#endif
