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
#include "population.h"

#include <cmath>

namespace simuPOP {

/* this is a class that store subpopulation, or subpopulation+virtual subpop id.
Basically, a integer number reprents a subpop, and a float number such as
1.1 represents the second virtual subpop in the second subpopulation.
*/
class virtualSubPopID
{
	virtualSubPopID(int id) : m_subPop(id), m_virtualSubPop(MaxSubPopID)
	{
	}

	virtualSubPopID(double id)
	{
		m_subPop = static_cast<SubPopID>(floor(id));
		m_virtualSubPop = static_cast<SubPopID>(floor((id-m_subPop)*10));
		DBG_ASSERT(fcmp_eq(id - m_subPop - 0.1*m_virtualSubPop, 0.), ValueError,
			"Wrong virtual subpopulation id " + toStr(id));
	}

	SubPopID id() const 
	{
		return m_subPop;
	}

	SubPopID vid() const
	{
		return m_virtualSubPop;
	}

	bool isVirtual() const
	{
		return static_cast<unsigned long>(m_virtualSubPop) != MaxSubPopID;
	}

private:
	SubPopID m_subPop;
	SubPopID m_virtualSubPop;
};

/** virtual splitters split a subpopulation into virtual sub-subpopulations.
   The virtual subpopulations does not have to add up to the whole
   subpopulation, nor they have to be distinct. For example,
   a virtual subpopulation may be all individuals in a population
   that is over the age of 30. Or, two virtual populations may
   overlap so some of the inviduals will go through more than one
   mating schemes.

   When a subpopulation is splitted by a virtual splitter, only its
   virtual subpopulations will go through mating. The virtual subpopulations
   are identified as virtual-subPopulations, relative to their parental
   subpopulation.

   virtual splitter will split a parental subpopulation, and the corresponding
   offspring subpopulation ino the same number of virtual subpopulations.
   Whereas parental virtual subpopulations can be virtual, offspring virtual
   subpopulations has to occupy the whole offspring subpopulation, and
   can not overlap. A wight-system is used to determine the offspring
   sub-subpopulation sizes. For each sub-subpopulation, weight
 \li 0 no sub-subpopulation. The corresponding parental virtual subpopulation
   	is ignored.
 \li -1 The same size as the corresponding parental virtual subpopulation
 \li n a weight. After 0, -1 are handled, the rest of the weights are
   	added and a proportion is calculated for each offspring sub-subpopulation.

   Note that if any parental virtual subpopulation is empty, its corresponding
   entry in the offspring sub-subpopulations is removed.
 */
class virtualSubPops
{
public:
	virtualSubPops(vectori const & offWeights) : m_offWeights(offWeights)
	{
	}


	virtual virtualSubPops * clone() const = 0;

	virtual ~virtualSubPops()
	{
	}


	// number of virtual subpops of subpopulation sp
	virtual UINT numVirtualSubPops(UINT sp) = 0;

	// prepare a subpopulation ssp, and return a real SP
	// id of the virtual subpopulation. The population
	// may be manipulated.
	virtual UINT prepareVirtualSubPop(population & pop, population & scratch,
	                                  UINT sp, UINT ssp) = 0;

	// because prepareVirtualSubPop may manipulate subpopulation
	// this function removes such manipulations. At least, there
	// should be no subpopulation exist in this subpopulation.
	virtual void restoreSubPop(population & pop, population & scratch,
	                           UINT sp) = 0;

private:
	vectori m_offWeights;
};


/** This splitter does nothing. It treats the whole
   subpopulation as the only virtual subpopulation.
 */
class nullVirtualSubPops : public virtualSubPops
{
public:
	nullVirtualSubPops() : virtualSubPops(vectori(1, 1))
	{
	}

	virtualSubPops * clone() const
	{
		return new nullVirtualSubPops(*this);
	}

	UINT numVirtualSubPops(UINT sp)
	{
		return 1;
	}


	// no need to prepare anything.
	UINT prepareVirtualSubPop(population & pop, population & scratch, UINT sp, UINT ssp)
	{
		return sp;
	}

	void restoreSubPop(population & pop, population & scratch,
	                           UINT sp)
	{
	}
};


/** duplicateVirtualSubPops does not split the parental subpopulation in any way.
   It presents the specified subpopulation as several subpopulations to the
   mating system, thus allow several different mating schemes to be appllied to the
   same subpopulation. For example, a selfing-mating scheme can be applied to
   the subpopulation and populate half of the offspring subpopulation, and a
   random-mating scheme can be used to generate the rest of the offspring
   subpopulation.
 */
class duplicateVirtualSubPops : public virtualSubPops
{
public:
	duplicateVirtualSubPops(vectori const & offWeights);

	virtualSubPops * clone() const
	{
		return new duplicateVirtualSubPops(*this);
	}

	UINT numVirtualSubPops(UINT sp)
	{
		return 1;
	}

	// no need to prepare anything.
	UINT prepareVirtualSubPop(population & pop, population & scratch, UINT sp, UINT ssp)
	{
		return sp;
	}

	void restoreSubPop(population & pop, population & scratch,
	                           UINT sp)
	{
	}


};


/** Split the population according to the value of an information field.
   A cutoff vector of length \c n is given to split the subpopulation into
 \c n+1 distinct virtual subpopulations.
 */
class infoVirtualSubPops : public virtualSubPops
{
public:
	infoVirtualSubPops(vectori const & offWeights);

	virtualSubPops * clone() const
	{
		return new infoVirtualSubPops(*this);
	}


	UINT numVirtualSubPops(UINT sp)
	{
		return 1;
	}

	// no need to prepare anything.
	UINT prepareVirtualSubPop(population & pop, population & scratch, UINT sp, UINT ssp)
	{
		return sp;
	}

	void restoreSubPop(population & pop, population & scratch,
	                           UINT sp)
	{
	}


};

/** Split the population according to a proportion */
class proportionVirtualSubPops : public virtualSubPops
{
public:
	proportionVirtualSubPops(vectori const & offWeights);

	virtualSubPops * clone() const
	{
		return new proportionVirtualSubPops(*this);
	}


	UINT numVirtualSubPops(UINT sp)
	{
		return 1;
	}

	// no need to prepare anything.
	UINT prepareVirtualSubPop(population & pop, population & scratch, UINT sp, UINT ssp)
	{
		return sp;
	}

	void restoreSubPop(population & pop, population & scratch,
	                           UINT sp)
	{
	}
};

/** Split the population using given ranges. The duplicateVirtualSubPops
   can be a special case of this virtual splitter.
 */
class rangeVirtualSubPops : public virtualSubPops
{
public:
	rangeVirtualSubPops(vectori const & offWeights);

	virtualSubPops * clone() const
	{
		return new rangeVirtualSubPops(*this);
	}


	UINT numVirtualSubPops(UINT sp)
	{
		return 1;
	}

	// no need to prepare anything.
	UINT prepareVirtualSubPop(population & pop, population & scratch, UINT sp, UINT ssp)
	{
		return sp;
	}

	void restoreSubPop(population & pop, population & scratch,
	                           UINT sp)
	{
	}
};

}
#endif
