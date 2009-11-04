/**
 *  $File: simulator.h $
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

#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>

#ifndef OPTIMIZED
#  include <time.h>

#  define InitClock(); \
    if (debug(DBG_PROFILE)) m_clock = clock();

#  define ElapsedTime(name); \
    if (debug(DBG_PROFILE)) \
	{ \
		cerr << name << ": " << static_cast<double>(clock() - m_clock) / CLOCKS_PER_SEC << "\n"; \
		m_clock = clock(); \
	}
#else
#  define InitClock();
#  define ElapsedTime(name);
#endif

/** \brief all classes in simuPOP is defined in this namespace
 */
namespace simuPOP {

/** This class implements a Python itertor class that can be used to iterate
    through populations in a population.
 */
class pyPopIterator
{
public:
	pyPopIterator(vector<population *>::iterator const begin,
		vector<population *>::iterator const end) :
		m_index(begin),
		m_end(end)
	{
	}


	~pyPopIterator()
	{
	}


	pyPopIterator __iter__()
	{
		return *this;
	}


	population & next();

private:
	// current (initial population)
	vector<population *>::iterator m_index;

	// ending idx
	vector<population *>::iterator m_end;
};


/** A simuPOP simulator is responsible for evolving one or more replicates
 *  of a \e population forward in time, subject to various \e operators.
 *  Populations in a simulator are created as identical copies of a population
 *  and will become different after evolution. A <em>mating scheme</em> needs
 *  to be specified, which will be used to generate offspring generations during
 *  evolution. A number of functions are provided to access simulator
 *  properties, access populations and their variables, copy, save and load a
 *  simulator.
 *
 *  The most important member function of a simulator is \c evolve, which
 *  evolves populations forward in time, subject to various \e operators. For
 *  convenience, member functions are provided to set virtual splitter, add
 *  information field and set ancestral depth to all populations in a
 *  simulator.
 */
class simulator
{
public:
	/** Create a simulator with \e rep replicates of population \e pop.
	 *  Population \e pop will be copied \e rep times (default to \c 1), while
	 *  keeping the passed population intact. A mating scheme \e matingScheme
	 *  will be used to evolve these populations.
	 *
	 *  \note Population \e pop is copied to a simulator so the input
	 *  population will be kept untouched.
	 */
	simulator(const population & pop, mating & matingScheme, UINT rep = 1);

	// destroy a simulator along with all its populations
	~simulator();

	/// CPPONLY Copy constructor
	simulator(const simulator & rhs);

	/** Clone a simulator, along with all its populations. Note that Python
	 *  assign statement <tt>simu1 = simu</tt> only creates a symbolic link to
	 *  an existing simulator.
	 *  <group>0-stru</group>
	 */
	simulator * clone() const;


	/** Return the current generation number, which is the initial generation
	 *  number (\c 0, or some value set by \c setGen(gen)) plus the total
	 *  number of generations evolved.
	 *  <group>1-gen</group>
	 */
	ULONG gen() const
	{
		return m_gen;
	}


	/** Set the current generation number of a simulator to \e gen.
	 *  <group>1-gen</group>
	 */
	void setGen(ULONG gen);


	/** Return the number of replicates.
	 *  <group>3-pop</group>
	 */
	UINT numRep() const
	{
		return m_numRep;
	}


	/** Return a reference to the \e rep-th population of a simulator. The
	 *  reference will become invalid once the simulator starts evolving or
	 *  becomes invalid (removed). Modifying the returned object is discouraged
	 *  because it will change the population within the simulator. If an
	 *  independent copy of the population is needed, use
	 *  <tt>simu.population(rep).clone()</tt>.
	 *  <group>3-pop</group>
	 */
	population & pop(UINT rep) const;

	/** Add a population \e pop to the end of an existing simulator. This
	 *  creates an cloned copy of \e pop in the simulator so the evolution
	 *  of the simulator will not change \e pop.
	 *  <group>4-modify</group>
	 */
	void add(const population & pop);

	/** Extract the \e rep-th population from a simulator. This will reduce
	 *  the number of populations in this simulator by one.
	 *  <group>3-pop</group>
	 */
	population & extract(UINT rep);

	/** Return a Python iterator that can be used to iterate through all
	 *  populations in a simulator.
	 *  <group>4-modify</group>
	 */
	pyPopIterator populations()
	{
		return pyPopIterator(m_ptrRep.begin(), m_ptrRep.end());
	}


	/** Evolve all populations \e gen generations, subject to several lists of
	 *  operators which are applied at different stages of an evolutionary
	 *  process. Operators \e initOps are applied to all populations (subject
	 *  to applicability restrictions of the operators, imposed by the \e rep
	 *  parameter of these operators) before evolution. They are used to
	 *  initialize populations before evolution. Operators \e finalOps are
	 *  applied to all populations after the evolution.
	 *
	 *  Operators \e preOps, \e duringOps and \e postOps are applied during the
	 *  life cycle of each generation. These operators can be applied at all or
	 *  some of the generations, to all or some of the evolving populations,
	 *  depending the \e begin, \e end, \e step, \e at and \e reps parameters
	 *  of these operators. These operators are applied in the order at which
	 *  they are specified. Populations in a simulator are evolved one by one.
	 *  At each generation, operators \e preOps are applied to the parental
	 *  generations. A mating scheme is then used to populate an offspring
	 *  generation. For each offspring, his or her sex is determined before
	 *  during-mating operators of the mating scheme are used to transmit
	 *  parental genotypes. During-mating operators specified in parameters
	 *  \e duringOps are applied afterwards. An offspring will be discarded if
	 *  any of the during-mating operator fails (return \c False). After an
	 *  offspring generation is successfully generated and becomes the current
	 *  generation, operators \e postOps are applied to the offspring
	 *  generation. If any of the \e preOps and \e postOps fails (return
	 *  \c False), the evolution of a population will be stopped. The
	 *  generation number of a population is increased by one if an offspring
	 *  generation has been successfully populated even if a post-during
	 *  operator fails.
	 *
	 *  Parameter \e gen can be set to a positive number, which is the number
	 *  of generations to evolve. Because a simulator always starts at the
	 *  beginning of a generation \c g (e.g. 0), a simulator will stop at the
	 *  beginning (instead of the end) of generation <tt>g + gen</tt> (e.g.
	 *  gen). If \e gen is negative (default), the evolution will continue
	 *  indefinitely, until all replicates are stopped by operators that return
	 *  \c False at some point (these operators are called \e terminators). At
	 *  the end of the evolution, the generations that each replicates have
	 *  evolved are returned. Note that \e finalOps are applied to all applicable
	 *  population, including those that have stopped before others.
	 *
	 *  The last parameter \e dryrun, if set to \c True, will print a
	 *  description of this evolutionary process. It can help you understand
	 *  what exactly will happen at each generation during an evolutionary
	 *  process.
	 *
	 *  <group>2-evolve</group>
	 */
	vectoru evolve(
		const opList & initOps = opList(),
		const opList & preOps = opList(),
		const opList & duringOps = opList(),
		const opList & postOps = opList(),
		const opList & finalOps = opList(),
		int gen = -1, bool dryrun = false);

	/// CPPONLY apply a list of operators to all populations
	bool apply(const opList & ops, bool dryrun = false);


	/** Set a new mating scheme \e matingScheme to a simulator.
	 *  <group>7-change</group>
	 */
	void setMatingScheme(const mating & matingScheme);

	/** Return the local namespace of the \e rep-th population, equivalent to
	 *  <tt>x.population(rep).vars(subPop)</tt>.
	 *  <group>9-var</group>
	 */
	PyObject * vars(UINT rep, vspID subPop = vspID())
	{
		if (static_cast<UINT>(rep) >= m_numRep)
			throw ValueError("Replicate index out of range.");
		return m_ptrRep[rep]->vars(subPop);
	}


	/// a Pyton function used to compare the simulator objects
	/// Note that mating schemes are not tested.
	int __cmp__(const simulator & rhs) const;

	/** Save a simulator to file \c filename, which can be loaded by a global
	 *  function \c LoadSimulator.
	 *  <group>0-stru</group>
	 */
	void save(string filename) const;

	/// CPPONLY load simulator from a file
	void load(string filename);

private:
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const UINT version) const
	{
		ULONG l_gen = gen();
		ar & l_gen;
		ar & m_numRep;

		DBG_DO(DBG_SIMULATOR, cerr << "Saving a simulator with "
			                       << m_numRep << " populations." << endl);

		// ignore scratch population
		for (UINT i = 0; i < m_numRep; i++)
			ar & (*m_ptrRep[i]);
	}


	template<class Archive>
	void load(Archive & ar, const UINT version)
	{
		ULONG l_gen;
		ar & l_gen;

		ar & m_numRep;

		DBG_DO(DBG_SIMULATOR, cerr << "Loading a simulator with "
			                       << m_numRep << " populations." << endl);

		m_ptrRep = vector<population *>(m_numRep);

		for (UINT i = 0; i < m_numRep; ++i) {
			m_ptrRep[i] = new population();
			ar & (*m_ptrRep[i]);
			m_ptrRep[i]->setRep(i);
		}
		m_scratchPop = new population(*m_ptrRep[0]);

		// only after all replicates are ready do we set gen
		setGen(l_gen);
	}


	BOOST_SERIALIZATION_SPLIT_MEMBER();

private:
	/// access scratch population
	population & scratchpopulation()
	{
		return *m_scratchPop;
	}


private:
	/// current generation
	ULONG m_gen;

	/// mating functor that will do the mating.
	mating * m_matingScheme;

	/// number of replicates of population
	UINT m_numRep;

	/// replicate pointers
	vector<population *> m_ptrRep;

	/// the scratch pop
	population * m_scratchPop;

#ifndef OPTIMIZED
	clock_t m_clock;
#endif

};

/// load a simulator from a file with the specified mating scheme. The file format is by default determined by file extension (<tt>format="auto"</tt>). Otherwise, \c format can be one of \c txt, \c bin, or \c xml.
simulator & LoadSimulator(const string & file,
	mating & matingScheme);

}

#ifndef SWIG
#  ifndef _NO_SERIALIZATION_
// version 0: base (reset for version 1.0)
BOOST_CLASS_VERSION(simuPOP::simulator, 0)
#  endif
#endif

#endif
