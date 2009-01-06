/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu
*                                                                         *
*   $LastChangedDate: 2007-11-07 23:33:59 -0600 (Wed, 07 Nov 2007) $
*   $Rev: 1293 $                                                       *
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

#ifndef _VIRTUALSUBPOP_H
#define _VIRTUALSUBPOP_H
/**
   \file
   \brief head file of class mating and its subclasses
 */
#include "utility.h"
#include "simuPOP_cfg.h"

#include <cmath>

namespace simuPOP {

class individual;
class population;

/**
 * A class to specify virtual subpopulation, which is composed of a subPopulation
 * ID and a virtual subpopulation ID.
 *
 */
class vspID
{
public:
	vspID(const vectori & subPop)
	{
		DBG_FAILIF(subPop.size() > 2, ValueError,
			"VSP should be specified as a subPop and virtualSubPop ID pair");
		m_subPop = subPop.size() > 0 ? subPop[0] : InvalidSubPopID;
		m_virtualSubPop = subPop.size() > 1 ? subPop[1] : InvalidSubPopID;
	}


	vspID(SubPopID subPop = InvalidSubPopID, SubPopID virtualSubPop = InvalidSubPopID)
		: m_subPop(subPop), m_virtualSubPop(virtualSubPop)
	{
	}


	SubPopID subPop() const { return m_subPop; }
	SubPopID virtualSubPop() const { return m_virtualSubPop; }
	bool isVirtual() const { return m_virtualSubPop != InvalidSubPopID; }

private:
	SubPopID m_subPop;
	SubPopID m_virtualSubPop;
};


/** This class is the base class of all virtual subpopulation (VSP) splitters,
 *  which provide ways to define groups of individuals in a subpopulation who
 *  share certain properties. A splitter defines a fixed number of named VSPs.
 *  They do not have to add up to the whole subpopulation, nor do they have to
 *  be distinct. After a splitter is assigned to a population, many functions
 *  and operators can be applied to individuals within specified VSPs.
 *
 *  Only one VSP splitter can be assigned to a population, which defined VSPs
 *  for all its subpopulations. It different splitters are needed for different
 *  subpopulations, a \c combinedSplitter should be.
 */
class vspSplitter
{

public:
	/** This is a virtual class that cannot be instantiated.
	 */
	vspSplitter() : m_activated(InvalidSubPopID)
	{
	}


	/** All VSP splitter defines a \c clone() function to create an identical
	 *  copy of itself.
	 */
	virtual vspSplitter * clone() const = 0;

	virtual ~vspSplitter()
	{
	}


	/// Which subpopulation is activated.
	/// CPPONLY
	SubPopID activatedSubPop() const
	{
		return m_activated;
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	virtual ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const = 0;

	/** Return the number of VSPs defined by this splitter.
	 */
	virtual UINT numVirtualSubPop() = 0;

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	virtual void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		IterationType type) = 0;

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	virtual void deactivate(population & pop, SubPopID subPop) = 0;

	/** Return the name of VSP \e vsp (an index between \c 0 and
	 *  <tt>numVirtualSubPop()</tt>).
	 */
	virtual string name(SubPopID vsp) = 0;

protected:
	ULONG countVisibleInds(const population & pop, SubPopID sp) const;

	SubPopID m_activated;

	void resetSubPop(population & pop, SubPopID subPop);

};

typedef std::vector<vspSplitter *> vectorsplitter;

/** This splitter takes several splitters and stacks their VSPs together. For
 *  example, if the first splitter defines \c 3 VSPs and the second splitter
 *  defines \c 2, the two VSPs from the second splitter becomes the fourth
 *  (index \c 3) and the fifth (index \c 4) VSPs of the combined splitter.
 *  This splitter is usually used to define different types of VSPs to a
 *  population.
 */
class combinedSplitter : public vspSplitter
{
public:
	/** Create a combined splitter using a list of \e splitters. For example,
	 *  <tt>combinedSplitter([sexSplitter(), affectionSplitter()])</tt> defines
	 *  a combined splitter with four VSPs.
	 */
	combinedSplitter(const vectorsplitter & splitters = vectorsplitter());

	~combinedSplitter();

	/// HIDDEN
	vspSplitter * clone() const;

	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/** Return the number of VSPs defined by this splitter, which is the sum of
	 *  the number of VSPs of all combined splitters.
	 */
	UINT numVirtualSubPop()
	{
		return m_numVSP;
	}


	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		IterationType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/** Return the name of a VSP \e vsp, which is the name a VSP defined by one
	 *  of the combined splitters.
	 */
	string name(SubPopID vsp);

private:
	/// the splitters
	vector<vspSplitter *> m_splitters;

	/// total number of vsp
	int m_numVSP;
	/// the splitter correspond to a vsp
	vectori m_splitter;
	/// the real vsp of a splitter correspond to a vsp
	vectori m_vsp;
	/// currently activated splitter
	int m_curSplitter;
};


/** This splitter defines two VSPs by individual sex. The first VSP consists of
 *  all male individuals and the second VSP consists of all females in a
 *  subpopulation.
 */
class sexSplitter : public vspSplitter
{
public:
	/// Create a sex splitter that defines male and female VSPs.
	sexSplitter() : vspSplitter()
	{
	}


	/// HIDDEN
	vspSplitter * clone() const
	{
		return new sexSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/// Return \c 2.
	UINT numVirtualSubPop()
	{
		return 2;
	}


	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		IterationType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/// Return \c "Male" if \e vsp=0 and \c "Female" otherwise.
	string name(SubPopID vsp)
	{
		DBG_FAILIF(vsp > 1, IndexError, "Virtual subpopulation index out of range");
		return vsp == 0 ? "Male" : "Female";
	}


};


/** This class defines two VSPs according individual affection status. The
 *  first VSP consists of unaffected invidiauls and the second VSP consists
 *  of affected ones.
 */
class affectionSplitter : public vspSplitter
{
public:
	/// Create a splitter that defined two VSPs by affection status.
	affectionSplitter() : vspSplitter()
	{
	}


	/// HIDDEN
	vspSplitter * clone() const
	{
		return new affectionSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/// Return 2.
	UINT numVirtualSubPop()
	{
		return 2;
	}


	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		IterationType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/// Return \c "Unaffected" if \e vsp=0 and \c "Affected" if \e vsp=1.
	string name(SubPopID vsp)
	{
		DBG_FAILIF(vsp > 1, IndexError, "VSP index out of range");
		return vsp == 0 ? "Unaffected" : "Affected";
	}


};


/** This splitter defines VSPs according to the value of an information field
 *  of each indivdiual. A VSP is defined either by a value or a range of values.
 */
class infoSplitter : public vspSplitter
{
public:
	/** Create an infomration splitter using information field \e field. If
	 *  parameter \e values is specified, each item in this list defines a VSP
	 *  in which all individuals have this value at information field \e field.
	 *  If a set of cutoff values are defined in parameter \e cutoff,
	 *  individuals are grouped by intervals defined by these cutoff values.
	 *  For example, <tt>cutoff=[1,2]</tt> defines three VSPs with
	 *  <tt>v < 1</tt>, <tt>1 <= v < 2</tt> and <tt>v >=2</tt> where \c v is the
	 *  value of an individual at information field \e field. Of course, only
	 *  one of the parameters \e values and \e cutoff should be defined,
	 *  values in \e cutoff should be distinct, and in an increasing order.
	 */
	infoSplitter(string field, vectorinfo const & values = vectorinfo(),
		vectorf const & cutoff = vectorf());

	/// HIDDEN
	vspSplitter * clone() const
	{
		return new infoSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/** Return the number of VSPs defined by this splitter, which is the length
	 *  parameter \e values or the length of \e cutoff plus one, depending on
	 *  which parameter is specified.
	 */
	UINT numVirtualSubPop();

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		IterationType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/** Return the name of a VSP \e vsp, which is <tt>field = value</tt> if VSPs
	 *  are defined by values in parameter \e values, or <tt>field < value</tt>
	 *  (the first VSP), <tt>v1 <= field < v2</tt> and <tt>field >= v</tt> (the
	 *  last VSP) if VSPs are defined by cutoff values.
	 */
	string name(SubPopID vsp);

private:
	string m_info;
	//
	vectorinfo m_values;
	//
	vectorf m_cutoff;
};


/** This splitter divides subpopulations into several VSPs by proportion.
 */
class proportionSplitter : public vspSplitter
{
public:
	/** Create a splitter that divides subpopulations by \e proportions, which
	 *  should be a list of float numbers (between \c 0 and \c 1) that add up
	 *  to \c 1.
	 */
	proportionSplitter(vectorf const & proportions = vectorf());

	/// HIDDEN
	vspSplitter * clone() const
	{
		return new proportionSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/** Return the number of VSPs defined by this splitter, which is the length
	 *  of parameter \e proportions.
	 */
	UINT numVirtualSubPop();

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		IterationType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/** Return the name of VSP \e vsp, which is <tt>"Prop p"</tt> where
	 *  <tt>p=propotions[vsp]</tt>.
	 */
	string name(SubPopID vsp);

private:
	vectorf m_proportions;
};


/** This class defines a splitter that groups individuals in certain ranges
 *  into VSPs.
 */
class rangeSplitter : public vspSplitter
{
public:
	/** Create a splitter according to a number of individual ranges defined
	 *  in \e ranges. For example,
	 *  <tt>rangeSplitter(ranges=[[0, 20], [40, 50]])</tt> defines two VSPs.
	 *  The first VSP consists of individuals \c 0, \c 1, ..., \c 19, and the
	 *  sceond VSP consists of individuals \c 40, \c 41, ..., \c 49. Note that
	 *  a nested list has to be used even if only one range is defined.
	 */
	rangeSplitter(const intMatrix & ranges);

	/// HIDDEN
	vspSplitter * clone() const
	{
		return new rangeSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/** Return the number of VSPs, which is the number of ranges defined in
	 *  parameter \e ranges.
	 */
	UINT numVirtualSubPop();

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		IterationType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/** Return the name of VSP \e vsp, which is <tt>"Range [a, b]"</tt> where
	 *  <tt>[a, b]</tt> is range <tt>ranges[vsp]</tt>.
	 */
	string name(SubPopID vsp);

private:
	intMatrix m_ranges;
};


/** This class defines a VSP splitter that defines VSPs according to individual
 *  genotype at specified loci.
 */
class genotypeSplitter : public vspSplitter
{
public:
	/** Create a splitter that defines VSPs by individual genotype at \e loci
	 *  (a locus index or a list of loci indexes). Each list in a list
	 *  \e allele defines a VSP, which is a list of allowed alleles at these
	 *  \e loci. If only one VSP is defined, the outer list of the nested list
	 *  can be ignored. If phase if true, the order of alleles in each list is
	 *  significant. If more than one set of alleles are given, individuals
	 *  having either of them is qualified.\n
	 *
	 *  For example, in a haploid population, <tt>loci=1, alleles=[0, 1]</tt>
	 *  defines a VSP with individuals having allele \c 0 or \c 1 at locus \c 1,
	 *  <tt>alleles=[[0, 1], [2]]</tt> defines two VSPs with indivdiuals in the
	 *  second VSP having allele \c 2 at locus \c 1. If multiple loci are
	 *  involved, alleles at each locus need to be defined. For example,
	 *  VSP defined by <tt>loci=[0, 1], alleles=[0, 1, 1, 1]</tt> consists of
	 *  individuals having alleles <tt>[0, 1]</tt> or <tt>[1, 1]</tt> at loci
	 *  <tt>[0, 1]</tt>.\n
	 *
	 *  In a haploid population, <tt>loci=1, alleles=[0, 1]</tt> defines a VSP
	 *  with individuals having genotype <tt>[0, 1]</tt> or <tt>[1, 0]</tt> at
	 *  locus \c 1. <tt>alleles[[0, 1], [2, 2]]</tt> defines two VSPs with
	 *  indivdiuals in the second VSP having genotype <tt>[2, 2]</tt> at locus
	 *  \c 1. If \e phase is set to \c True, the first VSP will only has
	 *  individuals with genotype <tt>[0, 1]</tt>. In the multiple loci case,
	 *  alleles should be arranged by haplotypes, for example,
	 *  <tt>loci=[0, 1], alleles=[0, 0, 1, 1], phase=True</tt> defines a VSP with
	 *  individuals having genotype <tt>-0-0-, -1-1-</tt> at loci \c 0 and \c 1.
	 *  If <tt>phase=False</tt> (default), genotypes <tt>-1-1-, -0-0-</tt>,
	 *  <tt>-0-1-</tt> and <tt>-1-0-</tt> are all allowed.
	 */
	genotypeSplitter(const uintList & loci,
		const intMatrix & alleles, bool phase = false);

	/// HIDDEN
	vspSplitter * clone() const
	{
		return new genotypeSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/// number of virtual subpops of subpopulation sp
	UINT numVirtualSubPop();

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		IterationType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/** Return name of VSP \e vsp, which is <tt>"Genotype loc1,loc2:genotype"</tt>
	 *  as defined by parameters \e loci and \e alleles.
	 */
	string name(SubPopID vsp);

private:
	bool match(const individual * ind, const vectori & alleles) const;

	bool matchSingle(const individual * ind, const vectori & alleles) const;

private:
	uintList m_loci;
	intMatrix m_alleles;
	bool m_phase;
};

}
#endif
