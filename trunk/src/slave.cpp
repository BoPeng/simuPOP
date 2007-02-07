/***************************************************************************
 *   Copyright (C) 2004 by Bo Peng                                         *
 *   bpeng@rice.edu                                                        *
 *                                                                         *
 *   $LastChangedDate: 2006-02-21 15:27:25 -0600 (Tue, 21 Feb 2006) $
 *   $Rev: 191 $
 *
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

#ifdef SIMUMPI
#include "slave.h"
#include "population.h"

#include <map>
using std::map;

#include <boost/parallel/mpi.hpp>
namespace mpi = boost::parallel::mpi;

namespace simuPOP
{

	// NOTE:
	// populations are saved in a map ('ID', population pointer) in slave's local
	// node. These population keeps minimal information to perform specified operation.
	//
	// For each operation requested, the master node usually send:
	//
	// 1. operation code
	// 2. population address in the master node (this) which is the ID in the slave code.
	// 3. paarmeters.
	//

	typedef map<int, population *> popMap;
	popMap g_popMap;
	int g_curID = 0;
	population * g_curPop = NULL;

#define setCurPop(id) \
	if (g_curID != id) \
	{ \
		g_curID = id; \
		g_curPop = g_popMap[id]; \
	}

	population * slavePopulationCreate()
	{
		//
		ULONG size;
		UINT ploidy;
		vectoru loci;
		bool sexChrom;
		vectorf lociPos;
		vectorlu subPop;
		int ancestralDepth;
		vectorstr alleleNames;
		vectorstr lociNames;
		UINT maxAllele;
		vectorstr infoFields;
		vectori chromMap;
		// get parameters
		mpiComm().recv(0, 2, size);
		mpiComm().recv(0, 3, ploidy);
		mpiComm().recv(0, 4, loci);
		mpiComm().recv(0, 5, sexChrom);
		mpiComm().recv(0, 6, lociPos);
		mpiComm().recv(0, 7, subPop);
		mpiComm().recv(0, 8, ancestralDepth);
		mpiComm().recv(0, 9, alleleNames);
		mpiComm().recv(0, 10, lociNames);
		mpiComm().recv(0, 11, maxAllele);
		mpiComm().recv(0, 12, infoFields);
		mpiComm().recv(0, 13, chromMap);
		// create population
		return new population(size, ploidy, loci, sexChrom,
			lociPos, subPop, ancestralDepth, alleleNames,
			lociNames, maxAllele, infoFields, chromMap);
	}

	Allele slavePopulationGetAllele(population * pop)
	{
		re
		mpiComm().recv(0, 2, idx);
	}
	
	bool slaveExecutionLoop()
	{
		// get population ID
		while(true)
		{
			int code = 0;
			ULONG id = 0;
			// operation code
			mpiComm().recv(0, 0, code);
			DBG_DO(DBG_MPI, cout << "Get code " << code << endl);
			if (code == SLAVE_TERMINATE)
				break;
			// object id
			mpiComm().recv(0, 1, id);
			DBG_DO(DBG_MPI, cout << "Object ID " << id << endl);
			//
			// operations:
			switch(code)
			{
				case SLAVE_POPULATION_CREATE:
					g_popMap[id] = slavePopulationCreate();
					setCurPop(id);
					break;
				case SLAVE_POPULATION_GET_ALLELE:
					setCurPop(id);
					slavePopulationGetAllele(g_curPop);
					break;
				default:
					throw ValueError("Wrong slave operation code");
			}
		}
		return true;
	}
}
#endif
