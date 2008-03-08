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
/** virtual subpopss split a subpopulation into virtual sub-subpopulations.
   The virtual subpopulations do not have to add up to the whole
   subpopulation, nor they have to be distinct. For example,
   a virtual subpopulation may be all individuals in a population
   that is over the age of 30. Or, two virtual populations may
   overlap so some of the inviduals may belong to more than one virtual
   subpopulations.
 */
class vspSplitter
{
public:
	// iteratable and visible are two different concepts.
	// When a population is activated by setting the visible flag
	// All operations that work on this subpopulation will be limited
	// to visible individuals. To be exact, IndIterator would skip
	// invisible individuals.
	//
	// When a population is activated by setting the iteratable flag,
	// only operations that respect this flag would check it and
	// respond to it.
	//
	enum activateType {
		Iteratable,
		Visible,
	};

public:
	vspSplitter() : m_activated(false)
	{
	}


	virtual vspSplitter * clone() const = 0;

	virtual ~vspSplitter()
	{
	}


	/// if the virtual subpopulation is activated.
	/// CPPONLY
	bool activated() const
	{
		return m_activated;
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	virtual ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const = 0;

	/// number of virtual subpops of subpopulation sp
	virtual UINT numVirtualSubPop() = 0;

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	virtual void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		activateType type) = 0;

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	virtual void deactivate(population & pop, SubPopID subPop) = 0;

	/// name of a virtual subpopulation
	virtual string name(SubPopID sp) = 0;

protected:
	ULONG countVisibleInds(const population & pop, SubPopID sp) const;

	bool m_activated;

	void resetSubPop(population & pop, SubPopID subPop);

};

typedef std::vector<vspSplitter *> vectorvsp;

/** This plitter takes several splitters, and stacks their virtual
subpopulations together. For example, if the first splitter has
three vsp, the second has two. The two vsp from the second splitter
will be the fouth (index 3) and fifth (index 4) of the combined 
splitter. 
*/
class combinedSplitter : public vspSplitter
{
public:
	combinedSplitter(const vectorvsp & splitters = vectorvsp());
	
	~combinedSplitter();
	
	vspSplitter * clone() const;
	
	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/// number of virtual subpops of subpopulation sp
	UINT numVirtualSubPop()
	{
		return m_numVSP;
	}

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		activateType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/// name of a virtual subpopulation
	string name(SubPopID sp);

private:
	/// the splitters
	vector<vspSplitter*> m_splitters;
	
	/// total number of vsp
	int m_numVSP;
	/// currently activated splitter
	int m_curSplitter;
	/// the splitter correspond to a vsp
	vectori m_splitter;
	/// the real vsp of a splitter correspond to a vsp
	vectori m_vsp;
};


/** split the population into Male and Female virtual subpopulations
 */
class sexSplitter : public vspSplitter
{
public:
	sexSplitter() : vspSplitter()
	{
	}


	vspSplitter * clone() const
	{
		return new sexSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/// number of virtual subpops of subpopulation sp
	UINT numVirtualSubPop()
	{
		return 2;
	}


	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		activateType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/// name of a virtual subpopulation
	string name(SubPopID sp)
	{
		return sp == 0 ? "Male" : "Female";
	}


};


/** split a subpopulation into unaffected and affected virtual subpopulations.
 */
class affectionSplitter : public vspSplitter
{
public:
	affectionSplitter() : vspSplitter()
	{
	}


	vspSplitter * clone() const
	{
		return new affectionSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/// number of virtual subpops of subpopulation sp
	UINT numVirtualSubPop()
	{
		return 2;
	}


	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		activateType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/// name of a virtual subpopulation
	string name(SubPopID sp)
	{
		return sp == 0 ? "Unaffected" : "Affected";
	}


};


/** Split the population according to the value of an information field.
   A list of distinct values, or a cutoff vector can be given
   to determine how the virtual subpopulations are divided.
   Note that in the first case, an individual does not have to belong to
   any virtual subpopulation.
 */
class infoSplitter : public vspSplitter
{
public:
	/**
	\param info name of the information field
	\param values a list of values, each defines a virtual subpopulation
	\param cutoff a list of cutoff values. For example, cutoff=[1, 2]
		defines three virtual subpopulations with v < 1, 1 <= v < 2,
		and v >= 2.
	*/
	infoSplitter(string info, vectorinfo const & values = vectorinfo(),
	             vectorf const & cutoff = vectorf());

	/// CPPONLY
	vspSplitter * clone() const
	{
		return new infoSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/// number of virtual subpops of subpopulation sp
	UINT numVirtualSubPop();

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		activateType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/// name of a virtual subpopulation
	string name(SubPopID sp);

private:
	string m_info;
	//
	vectorinfo m_values;
	//
	vectorf m_cutoff;
};


/** Split the population according to a proportion */
class proportionSplitter : public vspSplitter
{
public:
	/** \param proportions A list of float numbers (between 0 and 1) that
		defines the proportion of individuals in each virtual subpopulation.
		These numbers should add up to one.
	*/
	proportionSplitter(vectorf const & proportions = vectorf());

	vspSplitter * clone() const
	{
		return new proportionSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/// number of virtual subpops of subpopulation sp
	UINT numVirtualSubPop();

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		activateType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/// name of a virtual subpopulation
	string name(SubPopID sp);

private:
	vectorf m_proportions;
};


/** split the population according to individual range. The ranges
can overlap and does not have to add up to the whole subpopulation.
 */
class rangeSplitter : public vspSplitter
{
public:
	/**
	 \param range a shortcut for ranges=[range]
	 \param ranges a list of ranges
	 */
	rangeSplitter(const intMatrix & ranges);

	vspSplitter * clone() const
	{
		return new rangeSplitter(*this);
	}


	/// the size of a given virtual subpopulation.
	/// CPPONLY
	ULONG size(const population & pop, SubPopID subPop, SubPopID virtualSubPop) const;

	/// number of virtual subpops of subpopulation sp
	UINT numVirtualSubPop();

	/// mark individuals in the given vsp as visible, and others invisible.
	/// CPPONLY
	void activate(population & pop, SubPopID subPop, SubPopID virtualSubPop,
		activateType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/// name of a virtual subpopulation
	string name(SubPopID sp);

private:
	intMatrix m_ranges;
};


/** split the population according to given genotype
 */
class genotypeSplitter : public vspSplitter
{
public:
	/**
	 \param locus a shortcut to loci=[locus]
	 \param loci A list of locus at which alleles are used to classify individuals
	 \param alleles a list (for each virtual subpopulation), of a list
		of alleles at each locus. If phase if true, the order of alleles
		is significant. If more than one set of alleles are given, 
		individuals having either of them is qualified.
	\param phase whether or not phase is respected.

	For example,
	Genotype Aa or aa at locus 1:
		locus = 1, alleles = [0, 1]
	Genotype Aa at locus 1 (assuming A is 1):
		locus = 1, alleles = [1, 0], phase = True
	Genotype AaBb at loci 1 and 2:
		loci = [1, 2], alleles = [1, 0, 1, 0], phase = True
	Two virtual subpopulations with Aa and aa
		locus = 1, alleles = [[1, 0], [0, 0]], phase = True
	A virtual subpopulation with Aa or aa
		locus = 1, alleles = [1, 0, 0, 0]
	Two virtual subpopulation with genotype AA and the rest
		locus = 1, alleles = [[1, 1], [1, 0, 0, 0]], phase = False
	 */
	genotypeSplitter(const vectori & loci,
	                 const intMatrix & alleles, bool phase = false);

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
		activateType type);

	/// deactivate. Namely make all individuals visible again.
	/// CPPONLY
	void deactivate(population & pop, SubPopID sp);

	/// name of a virtual subpopulation
	string name(SubPopID sp);

private:
	bool match(const individual * ind, const vectori & alleles) const;
	bool matchSingle(const individual * ind, const vectori & alleles) const;

private:
	vectori m_loci;
	intMatrix m_alleles;
	bool m_phase;
};

}
#endif
