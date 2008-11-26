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
		cout << name << ": " << static_cast<double>(clock() - m_clock) / CLOCKS_PER_SEC << "\n"; \
		m_clock = clock(); \
	}
#else
#  define InitClock();
#  define ElapsedTime(name);
#endif

/** \brief all classes in simuPOP is defined in this namespace
 */
namespace simuPOP {

/** CPPONLY
    this class implements a Python itertor class that can be used to iterate
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
 *  simulator.\n
 *
 *  The most important member function of a simulator is \c evolve, which
 *  evolves populations forward in time, subject to various \e operators.
 *  Because populations in a simulator have to keep the same genotypic
 *  structure, several functions are provided to change ancestral depth and
 *  information fields of all populations. These functions cannot be replaced
 *  by similar calls to all populations in a simulator because the genotypic
 *  structure of the simulator itself needs to be updated.
 */
class simulator : public GenoStruTrait
{
public:
	/** Create a simulator with \e rep replicates of population \e pop.
	 *  Population \e pop will be copied \e rep times (default to \c 1), while
	 *  keeping the passed population intact. A mating scheme \e matingScheme
	 *  will be used to evolve these populations.
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


	/** Return the \e rep-th population of a simulator, in the form of a
	 *  reference (<tt>byRef=True</tt> (default)), or a cloned copy. In the
	 *  first case, a temporary reference is returned, which will become
	 *  invalid once the simulator starts evolving or becomes invalid (removed).
	 *  Modifying the returned object is discouraged because it will change
	 *  the population within the simulator. In the second case, an independent,
	 *  cloned copy of the internal population is returned. Modifying the
	 *  returned population will not change the simulator.
	 *  <group>3-pop</group>
	 */
	population & pop(UINT rep, bool byRef = true) const;

	/** Extract the \e rep-th population from a simulator. This will reduce
	 *  the number of populations in this simulator by one.
	 *  <group>3-pop</group>
	 */
	population & extract(UINT rep);

	/** Return a Python iterator that can be used to iterate through all
	 *  populations in a simulator.
	 *  <group>3-pop</group>
	 */
	pyPopIterator populations()
	{
		return pyPopIterator(m_ptrRep.begin(), m_ptrRep.end());
	}


	/// CPPONLY set population... unsafe!
	void setPopulation(population & pop, UINT rep);

	/** Evolve all populations \e gen generations, subject to operators \e ops
	 *  \e preOps and \e postOps. Operators \e preOps are applied to all
	 *  populations (subject to applicability restrictions of the operators,
	 *  imposed by the \e rep parameter of these operators) before evolution.
	 *  They are usually used to initialize populations. Operators \e postOps
	 *  are applied to all populations after the evolution.\n
	 *
	 *  Operators \e ops are applied during the life cycle of each generation.
	 *  Depending on the stage of these operators, they can be applied before-,
	 *  during-, and/or post-mating. These operators can be applied at all or
	 *  some of the generations, depending the \e begin, \e end, \e step, and
	 *  \e at parameters of these operators. Populations in a simulator are
	 *  evolved one by one. At each generation, the applicability of these
	 *  operators are determined. Pre-mating operators are applied to a
	 *  population first. A mating scheme is then used to populate an offspring
	 *  generation, using applicable during-mating operators. After an
	 *  offspring generation is successfully generated and becomes the current
	 *  generation, applicable post-mating operators are applied to it. Because
	 *  the order at which operators are applied can be important, and
	 *  the stage(s) at which operators are applied are not always clear,
	 *  a parameter \e dryRun can be used. If set to \c True, this function
	 *  will print out the order at which all operators are applied, without
	 *  actually evolving the populations.\n
	 *
	 *  Parameter \e gen can be set to a positive number, which is the number
	 *  of generations to evolve. If \e gen is negative (default), the evolution
	 *  will continue indefinitely, until all replicates are stopped by a
	 *  special kind of operators called \e terminators. At the end of the
	 *  evolution, the generations that each replicates have evolved are
	 *  returned.
	 *  <group>2-evolve</group>
	 */
	vectoru evolve(const vectorop & ops,
		const vectorop & preOps = vectorop(),
		const vectorop & postOps = vectorop(),
		int gen = -1, bool dryrun = false);

	/// CPPONLY apply a list of operators to all populations
	bool apply(const vectorop ops, bool dryrun = false);


	/** Add an information field \e field to all populations in a simulator,
	 *  and update the genotypic structure of the simulator itself. The
	 *  information field will be initialized by value \e init.
	 *  <group>7-change</group>
	 */
	void addInfoField(const string & field, double init = 0);

	/** Add information fields \e fields to all populations in a simulator,
	 *  and update the genotypic structure of the simulator itself. The
	 *  information field will be initialized by value \e init.
	 *  <group>7-change</group>
	 */
	void addInfoFields(const vectorstr & fields, double init = 0);

	/** Set ancestral depth of all populations in a simulator.
	 *  <group>7-change</group>
	 */
	void setAncestralDepth(UINT depth);


	/** Set a new mating scheme \e matingScheme to a simulator.
	 *  <group>7-change</group>
	 */
	void setMatingScheme(const mating & matingScheme);

	/** Return the local namespace of the \e rep-th population, equivalent to
	 *  <tt>x.population(rep).vars()</tt>.
	 *  <group>9-var</group>
	 */
	PyObject * vars(UINT rep)
	{
		if (static_cast<UINT>(rep) >= m_numRep)
			throw ValueError("Replicate index out of range.");

		return m_ptrRep[rep]->vars();
	}

	/** Return a dictionary of subpopulation variables in a local namespace of
     *  the \e rep-th population, equivalent to
     *  <tt>x.population(rep).vars(subPop)</tt>.
	 *  <group>9-var</group>
	 */
	PyObject * vars(UINT rep, vspID subPop)
	{
		if (static_cast<UINT>(rep) >= m_numRep)
			throw ValueError("Replicate index out of range.");

		return m_ptrRep[rep]->vars(subPop);
	}

	/** Save a simulator to file \c filename, which can be loaded by a global
	 *  function \c LoadSimulator.
	 *  <group>0-stru</group>
	 */
	void save(string filename) const;

	/// CPPONLY load simulator from a file
	void load(string filename);

	/// used by Python print function to print out the general information
	/// of the simulator
	string __repr__()
	{
		return "<simulator with " + toStr(numRep()) + " replicates>";
	}


private:
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const UINT version) const
	{
		ULONG l_gen = gen();
		ar & l_gen;
		ar & m_numRep;

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

		m_ptrRep = vector<population *>(m_numRep);

		for (UINT i = 0; i < m_numRep; ++i) {
			m_ptrRep[i] = new population();
			ar & (*m_ptrRep[i]);
			m_ptrRep[i]->setRep(i);
		}
		m_scratchPop = new population(*m_ptrRep[0]);
		setGenoStruIdx(m_ptrRep[0]->genoStruIdx());

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
