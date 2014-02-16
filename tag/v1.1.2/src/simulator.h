/**
 *  $File: simulator.h $
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

#ifndef _SIMULATOR_H
#define _SIMULATOR_H

/**
   \file
   \brief head file of class Simulator and some other utility functions
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

/** \brief all classes in simuPOP is defined in this namespace
 */
namespace simuPOP {

/** This class implements a Python itertor class that can be used to iterate
    through populations in a population.
 */
class pyPopIterator
{
public:
	pyPopIterator(vector<Population *>::iterator const begin,
		vector<Population *>::iterator const end) :
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


	Population & next();

	// Python 3 uses __next__ instead of next.
	Population & __next__()
	{
		return next();
	}


private:
	// current (initial population)
	vector<Population *>::iterator m_index;

	// ending idx
	vector<Population *>::iterator m_end;
};


/** A simuPOP simulator is responsible for evolving one or more populations
 *  forward in time, subject to various \e operators. Populations in a
 *  simulator are created from one or more replicates of specified populations.
 *  A number of functions are provided to access and manipulate populations,
 *  and most importantly, to evolve them.
 */
class Simulator
{
public:
	/** Create a simulator with \e rep (default to \c 1) replicates of
	 *  populations \e pops, which is a list of populations although a
	 *  single population object is also acceptable. Contents of passed
	 *  populations are by default moved to the simulator to avoid duplication
	 *  of potentially large population objects, leaving empty populations
	 *  behind. This behavior can be changed by setting \e stealPops to \c False,
	 *  in which case populations are copied to the simulator.
	 */
	Simulator(PyObject * pops, UINT rep = 1, bool stealPops = true);

	// destroy a simulator along with all its populations
	~Simulator();

	/// CPPONLY Copy constructor
	Simulator(const Simulator & rhs);

	/** Clone a simulator, along with all its populations. Note that Python
	 *  assign statement <tt>simu1 = simu</tt> only creates a symbolic link to
	 *  an existing simulator.
	 *  <group>0-stru</group>
	 */
	Simulator * clone() const;


	/** Return the number of replicates.
	 *  <group>3-pop</group>
	 */
	size_t numRep() const
	{
		return m_pops.size();
	}


	/** Return a reference to the \e rep-th population of a simulator. The
	 *  reference will become invalid once the simulator starts evolving or
	 *  becomes invalid (removed). If an independent copy of the population is
	 *  needed, you can use \c population.clone() to create a cloned copy or
	 *  \c simulator.extract() to remove the population from the simulator.
	 *  <group>3-pop</group>
	 */
	Population & population(size_t rep) const;

	/** Add a population \e pop to the end of an existing simulator. This
	 *  function by default moves \e pop to the simulator, leaving an empty
	 *  population for passed population object. If \e steal is set to
	 *  \c False, the population will be copied to the simulator, and thus
	 *  unchanged.
	 *  <group>4-modify</group>
	 */
	void add(const Population & pop, bool stealPop = true);

	/** Extract the \e rep-th population from a simulator. This will reduce
	 *  the number of populations in this simulator by one.
	 *  <group>3-pop</group>
	 */
	Population & extract(UINT rep);

	/** Return a Python iterator that can be used to iterate through all
	 *  populations in a simulator.
	 *  <group>4-modify</group>
	 */
	pyPopIterator populations()
	{
		return pyPopIterator(m_pops.begin(), m_pops.end());
	}


	/** CPPONLY describe a simulator.
	 */
	string describe(bool format = true) const;

	/** Evolve all populations \e gen generations, subject to several lists of
	 *  operators which are applied at different stages of an evolutionary
	 *  process. Operators \e initOps are applied to all populations (subject
	 *  to applicability restrictions of the operators, imposed by the \e rep
	 *  parameter of these operators) before evolution. They are used to
	 *  initialize populations before evolution. Operators \e finalOps are
	 *  applied to all populations after the evolution.
	 *
	 *  Operators \e preOps, and \e postOps are applied during the life cycle
	 *  of each generation. These operators can be applied at all or some of
	 *  the generations, to all or some of the evolving populations, depending
	 *  the \e begin, \e end, \e step, \e at and \e reps parameters of these
	 *  operators. These operators are applied in the order at which they are
	 *  specified. populations in a simulator are evolved one by one. At each
	 *  generation, operators \e preOps are applied to the parental
	 *  generations. A mating scheme is then used to populate an offspring
	 *  generation. For each offspring, his or her sex is determined before
	 *  during-mating operators of the mating scheme are used to transmit
	 *  parental genotypes. After an offspring generation is successfully
	 *  generated and becomes the current generation, operators \e postOps are
	 *  applied to the offspring generation. If any of the \e preOps and
	 *  \e postOps fails (return \c False), the evolution of a population will
	 *  be stopped. The generation number of a population, which is the
	 *  variable \c "gen" in each populations local namespace, is increased by
	 *  one if an offspring generation has been successfully populated even if
	 *  a post-mating operator fails. Another variable \c "rep" will also be
	 *  set to indicate the index of each population in the simulator. Note
	 *  that populations in a simulator does not have to have the same
	 *  generation number. You could reset a population's generation number
	 *  by changing this variable.
	 *
	 *  Parameter \e gen can be set to a non-negative number, which is the number
	 *  of generations to evolve. If a simulator starts at the beginning of a
	 *  generation \c g (for example 0), a simulator will stop at the beginning
	 *  (instead of the end) of generation <tt>g + gen</tt> (for example gen).
	 *  If \e gen is negative (default), the evolution will continue
	 *  indefinitely, until all replicates are stopped by operators that return
	 *  \c False at some point (these operators are called \e terminators). At
	 *  the end of the evolution, the generations that each replicates have
	 *  evolved are returned. Note that \e finalOps are applied to all applicable
	 *  population, including those that have stopped before others.
	 *
	 *  If parameter \e dryrun is set to \c True, this function will print a
	 *  description of the evolutionary process generated by function
	 *  \c describeEvolProcess() and exits.
	 *  <group>2-evolve</group>
	 */
	vectoru evolve(
		const opList & initOps = opList(),
		const opList & preOps = opList(),
		const MatingScheme & matingScheme = MatingScheme(),
		const opList & postOps = opList(),
		const opList & finalOps = opList(),
		int gen = -1, bool dryrun = false);


	/// CPPONLY apply a list of operators to all populations
	bool apply(const opList & ops);


	/** Return the local namespace of the \e rep-th population, equivalent to
	 *  <tt>x.Population(rep).vars(subPop)</tt>.
	 *  <group>9-var</group>
	 */
	PyObject * vars(UINT rep, vspID subPop = vspID())
	{
		if (static_cast<UINT>(rep) >= m_pops.size())
			throw ValueError("Replicate index out of range.");
		return m_pops[rep]->vars(subPop);
	}


	/// a Pyton function used to compare the simulator objects
	/// Note that mating schemes are not tested.
	int __cmp__(const Simulator & rhs) const;

private:
	/// access scratch population
	Population & scratchPopulation()
	{
		return *m_scratchPop;
	}


private:
	/// replicate pointers
	vector<Population *> m_pops;

	/// the scratch pop
	Population * m_scratchPop;

};


/** This function takes the same parameters as \c Simulator.evolve and
 *  output a description of how an evolutionary process will be executed.
 *  It is recommended that you call this function if you have any doubt how
 *  your simulation will proceed.
 */
string describeEvolProcess(
	const opList & initOps = opList(),
	const opList & preOps = opList(),
	const MatingScheme & matingScheme = MatingScheme(),
	const opList & postOps = opList(),
	const opList & finalOps = opList(),
	int gen = -1,
	size_t numRep = 1);

}

#endif
