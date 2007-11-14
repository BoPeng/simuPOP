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

/* this is a class that store subpopulation, or subpopulation+virtual subpop id.
   Basically, a integer number reprents a subpop, and a float number such as
   1.1 represents the second virtual subpop in the second subpopulation.
 */
class virtualSubPopID
{
public:
	virtualSubPopID() : m_subPop(MaxSubPopID), m_virtualSubPop(MaxSubPopID)
	{
	}
	
	virtualSubPopID(int id) : m_subPop(id), m_virtualSubPop(MaxSubPopID)
	{
	}


	virtualSubPopID(int id, int vid) : m_subPop(id),
		m_virtualSubPop(vid)
	{
	}


	virtualSubPopID(double id)
	{
		m_subPop = static_cast<SubPopID>(floor(id));
		m_virtualSubPop = static_cast<SubPopID>(floor((id - m_subPop) * 10));
		DBG_ASSERT(fcmp_eq(id, m_subPop + 0.1 * m_virtualSubPop), ValueError,
		    "Wrong virtual subpopulation id " + toStr(id) + " ("
		    + toStr(m_subPop) + ", " + toStr(m_virtualSubPop) + ")\n"
		                                                        "Please use virtualSubPopID(id, vid) instead.");
	}


	SubPopID id() const
	{
		return m_subPop;
	}


	SubPopID vid() const
	{
		return m_virtualSubPop;
	}


	bool isValid() const
	{
		return static_cast<unsigned long>(m_subPop) != MaxSubPopID;
	}

	bool isVirtual() const
	{
		return static_cast<unsigned long>(m_virtualSubPop) != MaxSubPopID;
	}


private:
	SubPopID m_subPop;
	SubPopID m_virtualSubPop;
};

/** virtual subpopss split a subpopulation into virtual sub-subpopulations.
   The virtual subpopulations do not have to add up to the whole
   subpopulation, nor they have to be distinct. For example,
   a virtual subpopulation may be all individuals in a population
   that is over the age of 30. Or, two virtual populations may
   overlap so some of the inviduals will go through more than one
   mating schemes.
 */
class vspSplitter
{
public:
	vspSplitter() : m_activated(false)
	{
	}


	virtual vspSplitter * clone() const = 0;

	virtual ~vspSplitter()
	{
	}


	bool activated() const
	{
		return m_activated;
	}


	virtual ULONG size(const population & pop, virtualSubPopID subPop) const = 0;

	// number of virtual subpops of subpopulation sp
	virtual UINT numVirtualSubPop() = 0;

	// prepare a subpopulation ssp, and return a real SP
	// id of the virtual subpopulation. The population
	// may be manipulated.
	virtual void activate(population & pop, virtualSubPopID subPop) = 0;

	// because prepareVirtualSubPop may manipulate subpopulation
	// this function removes such manipulations. At least, there
	// should be no subpopulation exist in this subpopulation.
	virtual void reset(population & pop, SubPopID subPop) = 0;

	virtual string name(SubPopID sp) = 0;

protected:
	bool m_activated;

	void resetSubPop(population & pop, SubPopID subPop);

};


/** This subpops does nothing. It treats the whole
   subpopulation as the only virtual subpopulation.
 */
class nullSplitter : public vspSplitter
{
public:
	nullSplitter() : vspSplitter()
	{
	}


	vspSplitter * clone() const
	{
		return new nullSplitter(*this);
	}


	ULONG size(const population & pop, virtualSubPopID subPop) const;


	UINT numVirtualSubPop()
	{
		return 1;
	}


	void activate(population & pop, virtualSubPopID subPop)
	{
		m_activated = true;
	}


	void reset(population & pop, SubPopID sp)
	{
		m_activated = false;
	}


	string name(SubPopID sp) { return "orignal pop"; }
};


/** duplicateSplitter does not split the parental subpopulation in any way.
   It presents the specified subpopulation as several subpopulations to the
   mating system, thus allow several different mating schemes to be appllied to the
   same subpopulation. For example, a selfing-mating scheme can be applied to
   the subpopulation and populate half of the offspring subpopulation, and a
   random-mating scheme can be used to generate the rest of the offspring
   subpopulation.
 */
class duplicateSplitter : public nullSplitter
{
public:
	duplicateSplitter(int num) : nullSplitter(), m_num(num)
	{
	}


	vspSplitter * clone() const
	{
		return new duplicateSplitter(*this);
	}


	UINT numVirtualSubPop()
	{
		return m_num;
	}


	ULONG size(const population & pop, virtualSubPopID subPop) const;

	string name(SubPopID sp)
	{
		return "replicate" + toStr(sp);
	}


private:
	int m_num;
};


/** split the population according to sex
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


	ULONG size(const population & pop, virtualSubPopID subPop) const;

	UINT numVirtualSubPop()
	{
		return 2;
	}


	void activate(population & pop, virtualSubPopID subPop);

	void reset(population & pop, SubPopID sp);

	string name(SubPopID sp)
	{
		return sp == 0 ? "Male" : "Female";
	}


};


/** split the population according to affection status
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


	ULONG size(const population & pop, virtualSubPopID subPop) const;

	UINT numVirtualSubPop()
	{
		return 2;
	}


	void activate(population & pop, virtualSubPopID subPop);

	void reset(population & pop, SubPopID sp);

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
	infoSplitter(string info, vectorinfo const & values = vectorinfo(),
	             vectorf const & cutoff = vectorf());

	vspSplitter * clone() const
	{
		return new infoSplitter(*this);
	}


	ULONG size(const population & pop, virtualSubPopID subPop) const;

	UINT numVirtualSubPop();

	void activate(population & pop, virtualSubPopID subPop);

	void reset(population & pop, SubPopID sp);

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
	proportionSplitter(vectorf const & proportions = vectorf());

	vspSplitter * clone() const
	{
		return new proportionSplitter(*this);
	}


	ULONG size(const population & pop, virtualSubPopID subPop) const;

	UINT numVirtualSubPop();

	void activate(population & pop, virtualSubPopID subPop);

	void reset(population & pop, SubPopID sp);

	string name(SubPopID sp);

private:
	vectorf m_proportions;
};


/** split the population according to individual range
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


	ULONG size(const population & pop, virtualSubPopID subPop) const;

	UINT numVirtualSubPop();

	void activate(population & pop, virtualSubPopID subPop);

	void reset(population & pop, SubPopID sp);

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
	 \param genotype a shortcut for genotypes=[genotype]
	 \param genotypes a list of genotypes
	 */
	genotypeSplitter(const vectori & loci,
	                 const intMatrix & alleles, bool phase = false);

	vspSplitter * clone() const
	{
		return new genotypeSplitter(*this);
	}


	ULONG size(const population & pop, virtualSubPopID subPop) const;

	UINT numVirtualSubPop();

	void activate(population & pop, virtualSubPopID subPop);

	void reset(population & pop, SubPopID sp);

	string name(SubPopID sp);

private:
	bool match(const individual * ind, const vectori & alleles) const;

private:
	vectori m_loci;
	intMatrix m_alleles;
	bool m_phase;
};

typedef vector<vspSplitter *> vectorvsp;
}
#endif
