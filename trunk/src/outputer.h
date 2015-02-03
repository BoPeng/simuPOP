/**
 *  $File: outputer.h $
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

#ifndef _OUTPUTER_H
#define _OUTPUTER_H
/**
   \file
   \brief head file of class outputer: public BaseOperator
 */
#include "utility.h"
#include "operator.h"
#include <iostream>

#include <iomanip>
using std::setw;
using std::hex;
using std::dec;

namespace simuPOP {

/** This operator outputs a given string when it is applied to a population.
 */
class PyOutput : public BaseOperator
{

public:
	/** Creates a \c PyOutput operator that outputs a string \e msg to
	 *  \e output (default to standard terminal output) when it is applied
	 *  to a population. Please refer to class \c BaseOperator for a detailed
	 *  description of common operator parameters such as \e stage, \e begin
	 *  and \e output.
	 */
	PyOutput(const string & msg = string(), const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = vectorstr()) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_string(msg)
	{
	}


	/// HIDDEN simply output some info
	virtual bool apply(Population & pop) const;

	/// destructor
	virtual ~PyOutput()
	{
	}


	/// HIDDEN Deep copy of a \e PyOutput operator.
	virtual BaseOperator * clone() const
	{
		return new PyOutput(*this);
	}


	/// HIDDEN
	string describe(bool format = true) const;

private:
	const string m_string;
};


/** This operator dumps the content of a population in a human readable format.
 *  Because this output format is not structured and can not be imported back
 *  to simuPOP, this operator is usually used to dump a small population to a
 *  terminal for demonstration and debugging purposes.
 */
class Dumper : public BaseOperator
{
public:
	/** Create a operator that dumps the genotype structure (if \e structure is
	 *  \c True) and genotype (if \e genotype is \c True) to an \e output (
	 *  default to standard terminal output). Because a population can be large,
	 *  this operator will only output the first 100 (parameter \e max)
	 *  individuals of the present generation (parameter \e ancGens). All loci
	 *  will be outputed unless parameter \e loci are used to specify a subset
	 *  of loci. This operator by default output values of all information fields
	 *  unless parameter \e infoFields is used to specify a subset of info fields
	 *  to display. If a list of (virtual) subpopulations are specified, this
	 *  operator will only output individuals in these outputs. Please refer to
	 *  class \c BaseOperator for a detailed explanation for common parameters
	 *  such as \e output and \e stage.
	 */
	Dumper(bool genotype = true, bool structure = true, const uintList & ancGens = uintList(NULL),
		int width = 1, UINT max = 100, const uintList & loci = vectoru(), const stringFunc & output = ">",
		int begin = 0, int end = -1, int step = 1, const intList & at = vectori(),
		const intList & reps = intList(), const subPopList & subPops = subPopList(),
		const stringList & infoFields = stringList()) :
		BaseOperator(output, begin, end, step, at, reps, subPops, infoFields),
		m_showGenotype(genotype), m_showStructure(structure), m_ancGens(ancGens), m_width(width),
		m_loci(loci.elems()), m_max(max)
	{
	}


	/// HIDDEN Deep copy of a Dumper operator.
	virtual BaseOperator * clone() const
	{
		return new Dumper(*this);
	}


	/// HIDDEN Apply a Dumper operator to population \e pop.
	virtual bool apply(Population & pop) const;

	/// destructor.
	virtual ~Dumper()
	{
	};

	/// HIDDEN
	string describe(bool format = true) const
	{
		(void)format;  // avoid warning about unused parameter
		return "<simuPOP.Dumper>" ;
	}


private:
	void displayStructure(const Population & pop, ostream & out) const;

	UINT displayGenotype(const Population & pop, const subPopList & subPops, ostream & out) const;

private:
	///
	const bool m_showGenotype;

	///
	const bool m_showStructure;

	///
	const uintList m_ancGens;

	/// disp width when outputing alleles
	const int m_width;

	///
	const vectoru m_loci;

	/// only output first ... individuals. Good for large population
	const UINT m_max;
};


/** An operator that save populations to specified files.
 */
class SavePopulation : public BaseOperator
{
public:
	/** Create an operator that saves a population to \e output when it is
	 *  applied to the population. This operator supports all output
	 *  specifications (\c '', \c 'filename', \c 'filename' prefixed by one
	 *  or more '>' characters, and \c '!expr') but output from different
	 *  operators will always replace existing files (effectively ignore
	 *  '>' specification). Parameter \e subPops is ignored. Please refer to
	 *  class \c BaseOperator for a detailed description about common operator
	 *  parameters such as \e stage and \e begin.
	 */
	SavePopulation(const stringFunc & output = "", int begin = 0, int end = -1,
		int step = 1, const intList & at = vectori(), const intList & reps = intList(),
		const subPopList & subPops = subPopList(), const stringList & infoFields = vectorstr()) :
		BaseOperator("", begin, end, step, at, reps, subPops, infoFields),
		m_filename(output.value())
	{
		DBG_WARNIF(output.empty(), "An empty output string is passed to operator SavePopulation. No file will be saved.");
	}


	/// destructor.
	~SavePopulation()
	{
	}


	/// HIDDEN Deep copy of a SavePopulation operator.
	virtual BaseOperator * clone() const
	{
		return new SavePopulation(*this);
	}


	/// HIDDEN apply operator to population \e pop.
	virtual bool apply(Population & pop) const;

	/// HIDDEN
	string describe(bool format = true) const;

private:
	/// filename,
	const string m_filename;
};

}
#endif
