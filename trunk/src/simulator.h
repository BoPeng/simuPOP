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
\brief head file of class Simulator and some other utility functions
*/

#include "simupop_cfg.h"
#include "utility.h"
#include "mating.h"
#include "operator.h"

#include <vector>
#include <algorithm>

using std::vector;
using std::fill;
using std::swap;

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#ifndef __NO_XML_SUPPORT__
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#endif
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
using boost::serialization::make_nvp;

#ifndef OPTIMIZED
#include <time.h>                                 // for clock() function

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

  Simulators combine three important components of simuPOP:
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

  Simulators separates operators into pre-, during- and post-
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

  Simulators can evolve a given number of generations (the
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
  three major components of simuPOP together: Population, Mating, Operation.
  More specifically, a simuPOP<Population> object manages several copies of
  Population (or one of its subclasses) plus a 'scratch population'; it has
  a Mating object that knows how to generate next gen; it controls
  the evolution process by applying pre- during- and post- mating Operators
  during evolution. }

  */
  template < class Pop >
    class Simulator : public GenoStruTrait
  {
    public:
      /// expose population type.
      typedef Pop PopType;

      /// expose individual type
      typedef typename PopType::IndType IndType;

      /// vector of operator pointers.
      typedef std::vector< Operator< Pop> * > vectorop;

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
      Simulator(Pop & population, Mating<Pop>&  matingScheme, string varName = "simuVars",
        int rep = 1, vectori grp=vectori())
        : m_gen(0), m_curRep(0), m_numRep(rep), m_groups(0),
        m_stopIfOneRepStop(false), m_applyOpToStoppedReps(false)
      {
        DBG_ASSERT(m_numRep >= 1, ValueError,
          "Number of replicates should be greater than or equal one.");

        DBG_DO(DBG_SIMULATOR, cout << "Creating simulator " << endl);

        m_matingScheme = matingScheme.clone();

        DBG_DO(DBG_SIMULATOR, cout << "Mating scheme copied" << endl);

        if (!m_matingScheme->isCompatible(population))
          throw TypeError
            ("Mating type is not compatible with current population settings.");

        // create replicates of given population
        m_ptrRep = new Pop *[m_numRep];
        m_scratchPop = new Pop(population);

        DBG_FAILIF( m_scratchPop == NULL,
          SystemError, "Fail to create scratch population");

        try
        {
          for (UINT i = 0; i < m_numRep; ++i)
          {
            m_ptrRep[i] = new Pop(population);

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
        // mainVars().setIntVar("gen", 0);

        DBG_DO(DBG_SIMULATOR, cout << "Simulator created" << endl);

      }

      /// destroy a simulator along with all its populations
      /**
      \note pop = simulator::population()
      returns temporary reference to an internal population.
      After a simulator evolves another genertion or
      after the simulator is destroyed, this referenced
      population should *not* be used.
      */
      ~Simulator()
      {
        // call the destructor of each replicates
        delete m_scratchPop;

        for (UINT i = 0; i < m_numRep; ++i)
          delete m_ptrRep[i];

        delete[] m_ptrRep;

        delete m_matingScheme;

        DBG_DO(DBG_SIMULATOR, cout << "Populations have been erased. " << endl
          << "If you have ever used population() function to access some of the replicates, " << endl
          << "these referenced population will not be working now." << endl);
      };

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
      PopType& population(UINT rep)
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
      PopType& getPopulation(UINT rep)
      {
        return *new PopType( *m_ptrRep[rep]);
      }

      // set population... unsafe!
      void setPopulation(Pop& pop, UINT rep)
      {
        DBG_FAILIF( rep >= m_numRep, IndexError,
          "replicate index out of range. From 0 to numRep()-1 ");
        // get a copy of population
        delete m_ptrRep[rep];
        m_ptrRep[rep] = new Pop(pop);

        if( pop.genoStru() != this->genoStru() )
        {
          DBG_DO(DBG_SIMULATOR,  cout << "Warning: added population has different genotypic structure." << endl);
          setGenoStruIdx( pop.genoStruIdx());
        }
        if( pop.genoStru() != m_scratchPop->genoStru() )
        {
          delete m_scratchPop;
          m_scratchPop = new Pop(pop);
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
      void setGroup(const vectori& grp)
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
        //  mainVars().setIntVar("gen", gen);
      }

      /// evolve one step
      /**
      \sa simulator::evolve()
      */
      bool step(const vectorop& ops = vectorop(),
        const vectorop& preOps = vectorop(),
        const vectorop& postOps = vectorop(), UINT steps=1)
      {
        // since end gen itself will be executed
        // cur = 5,
        // end = 6
        // will go two steps.
        return evolve( ops, preOps, postOps, gen() + steps - 1);
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
        int end =  -1, bool dryrun = false)
      {

        DBG_DO(DBG_SIMULATOR, cout << "Starting generation: " << gen()
          << " with ending generation " << end << endl);

        DBG_FAILIF(end > 0 && gen() > static_cast<UINT>(end), ValueError,
          "Population gen is already greater than ending generation.");

        DBG_FAILIF( end < 0 && ops.empty(), ValueError,
          "Evolve with unspecified ending generation should have at least one terminator (operator)");

        vectorop preMatingOps, durMatingOps, postMatingOps, activeDurMatingOps;

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
            durMatingOps.push_back(ops[i]);

        if(durMatingOps.size() > 1)
          cout << "Warning: More than one during mating operators. Make sure they act on different replicates." << endl;

        for (size_t i = 0; i < ops.size(); ++i)
          if (ops[i]->canApplyDuringMating() && !ops[i]->formOffGenotype())
            durMatingOps.push_back(ops[i]);

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

          for (m_curRep = 0; m_curRep < m_numRep; m_curRep++)
          {
            DBG_ASSERT( static_cast<int>(m_curRep) == curPop().rep(), SystemError,
              "Replicate number does not match");

            if(dryrun)
              cout << "  Replicate " << m_curRep << endl;

            if (stop[m_curRep])
            {
              /// if apply to stopped reps, do it.
              if( !m_applyOpToStoppedReps)
                continue;
            }

            size_t it = 0;                        // asign a value to reduce compiler warning

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
                  if (! preMatingOps[it]->applyWithScratch(curPop(), scratchPop(), PreMating))
                  {
                    DBG_DO(DBG_SIMULATOR, cout << "Pre-Mating Operator " + preMatingOps[it]->__repr__() +
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
                  cout << "Premating operator " << preMatingOps[it]->__repr__() << " throws an exception." << endl << endl;
                  throw;
                }
                ElapsedTime("PreMatingOp: " + preMatingOps[it]->__repr__());
              }
            }

            // start mating:
            // find out active during-mating operators
            if(dryrun)
              cout << "    Start mating." << endl;

            activeDurMatingOps.clear();
            for(typename vectorop::iterator op=durMatingOps.begin(), opEnd = durMatingOps.end();
              op != opEnd; ++op)
            {
              if(dryrun)
                cout << "      - " << (*op)->__repr__() << (*op)->atRepr() << endl;

              if( (*op)->isActive(m_curRep, m_numRep, m_gen, end, grp()))
                activeDurMatingOps.push_back( *op);
            }

            try
            {
              if (!dryrun && !m_matingScheme->mate(curPop(), scratchPop(), activeDurMatingOps, true))
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
              cout << "Mating or one of the during mating operator throws an exception." << endl << endl;
              throw;
            }

            ElapsedTime("MatingDone");

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
                  if (! postMatingOps[it]-> applyWithScratch(curPop(), scratchPop(), PostMating))
                  {
                    DBG_DO(DBG_SIMULATOR, cout << "Post-Mating Operator " + postMatingOps[it]->__repr__() +
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
                  cout << "Postmating operator " << postMatingOps[it]->__repr__() << " throws an exception." << endl << endl;
                  throw;
                }
                ElapsedTime("PostMatingOp: "+ postMatingOps[it]->__repr__());
              }
            }                                     // post mating ops
          }                                       // each replicates

          if(dryrun)
            break;

          // if one replicate stop and stopIfOneRepStop is set,
          // or if all replicates stop, or reach ending gen, stop iteration.
          // count the number of stopped replicates
          DBG_DO(DBG_SIMULATOR, cout << endl << "Number of stopped replicates: " << numStopped << endl);

          m_gen ++;                               // increase generation!
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

        }                                         // the big loop

        if( ! postOps.empty() )
        {
          if(dryrun)
            cout << "Apply post-evolution operators: " << endl;

          // finishing up, apply post-op
          apply(postOps, dryrun);
        }

        // close every opened file (but not append-cross-evolution ones)
        ostreamManager().closeAll(false);
        return true;
      }

      /// apply some ops, geneartion of Population does not change
      /// No mating is allowed.
      /**
      \param ops operators that will be applied at all generations.
        Of course they might not be active at all generations.
      \result True if evolution finishs successfully.

      pre-mating oeprators are applied before post-mating operators.
      no during-mating operators are allowed.
      */
      bool apply(const vectorop ops, bool dryrun=false)
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

              preMatingOps[it] -> applyWithScratch(curPop(),
                scratchPop(), PreMating);

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

              postMatingOps[it] -> applyWithScratch(curPop(), scratchPop(), PostMating);

              ElapsedTime("PrePost-postMatingop"+toStr(it));
            }
          }
        }
        return true;
      }

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

      \sa global function loadSimulator
      */
      void saveSimulator(string filename, string format="auto")
      {

        ofstream ofs(filename.c_str());
        if(!ofs)
          throw ValueError("Can not open file " + filename );

        if( format == "text"|| (format == "auto" && filename.substr(filename.size()-4,4) == ".txt" ))
        {
          boost::archive::text_oarchive oa(ofs);
          oa << *this;
        }
#ifndef __NO_XML_SUPPORT__
        else if (format == "xml" || (format == "auto" && filename.substr(filename.size()-4,4) == ".xml" ) )
        {
          boost::archive::xml_oarchive oa(ofs);
          oa << boost::serialization::make_nvp("simulator",*this);
        }
#endif
        else if (format == "bin"  ||  (format == "auto" && filename.substr(filename.size()-4,4) == ".bin" ))
        {
          boost::archive::binary_oarchive oa(ofs);
          oa << *this;
        }
        else
        {
          ofs.close();
          throw ValueError("Wrong format type. Use one of text, xml, bin.");
        }
        ofs.close();
      }

      /// CPPONLY load simulator from a file
      /**
      \param filename load from filename
      \param format format to load. Can be one of "text", "bin", "xml". It should match
        the format used to save the simulator.

      \sa saveSimulator
      */
      void loadSimulator(string filename, string format="auto")
      {
        ifstream ifs(filename.c_str());
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
          ifs.close();

          DBG_DO(DBG_POPULATION,
            cout << "Can not determine file type, or file type is wrong. Trying different ways." << endl);

          // open a fresh ifstream
          ifstream ifbin(filename.c_str());
          // try to load the file using different iarchives.
          try                                     // binary?
          {
            boost::archive::binary_iarchive ia(ifbin);
            ia >> *this;
          }
          catch(...)                              // not binary, text?
          {
            ifbin.close();
            ifstream iftxt(filename.c_str());
            try
            {
              boost::archive::text_iarchive ia(iftxt);
              ia >> *this;
            }
            catch(...)                            // then xml?
            {
              iftxt.close();
#ifndef __NO_XML_SUPPORT__
              ifstream ifxml(filename.c_str());
              try
              {
                boost::archive::xml_iarchive ia(ifxml);
                ia >> boost::serialization::make_nvp("simulator", *this);
              }
              catch(...)
              {
                ifxml.close();
                throw ValueError("Failed to load simulator. Your file may be corrupted, "
                  "or being a copy of non-transferrable file (.bin)");
              }
#else
              throw ValueError("Failed to load simulator. Your file may be corrupted, "
                "or being a copy of non-transferrable file (.bin)");
#endif
            }                                     // try xml
          }                                       // try text
        }                                         // try bin
      }

      // allow str(population) to get something better looking
      string __repr__()
      {
        return "<Simulator with " + toStr(numRep()) + " replicates>";
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

        m_ptrRep = new Pop *[m_numRep];

        for (UINT i = 0; i < m_numRep; ++i)
        {
          m_ptrRep[i] = new PopType();
          ar & make_nvp("populations", *m_ptrRep[i]);
          m_ptrRep[i]->setRep(i);
        }
        m_scratchPop = new PopType(*m_ptrRep[0]);
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
      PopType& curPop()
      {
        return *m_ptrRep[m_curRep];
      }

      /// access scratch population
      PopType&  scratchPop()
      {
        return *m_scratchPop;
      }

    private:

      /// current generation
      ULONG m_gen;

      /// index to current population, refernce to shared "%rep"
      UINT m_curRep;

      /// mating functor that will do the mating.
      Mating<Pop>* m_matingScheme;

      /// number of replicates of population
      UINT m_numRep;

      /// groups
      vectori m_groups;

      /// replicate pointers
      Pop **m_ptrRep;

      /// the scratch pop
      Pop *m_scratchPop;

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

}
#endif
