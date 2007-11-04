/***************************************************************************
*   Copyright (C) 2004 by Bo Peng                                         *
*   bpeng@rice.edu
*                                                                         *
*   $LastChangedDate$
*   $Rev$                                                       *
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

#ifndef _MISC_H
#define _MISC_H
/**
 \file
 \brief head file of class mating and its subclasses
 */
#include "utility.h"

namespace simuPOP {

//
// simulate trajectories of disease susceptibility loci using an extension of
// the backward method described in Slatkin 2001.
//
// Parameters:
//
//    curGen current generation number
//    N      constant population size N
//    NtFunc a python function that returns population size at each generation.
//           It should return an array of subpop sizes.
//           gen is defined in forward order. NtFunc(curGen) should be current
//           generation number.
//    freq   expected allele frequency of allele a.
//    fitness  constant fitness for [AA, Aa, aa, BB, Bb, bb ...]
//    fitnessFunc  a python function that returns selection pressure at each generation
//           the function expects parameters gen and freq. gen is current generation
//           number and freq is the allele frequency at all loci. This allows
//           frequency dependent selection. gen is defined in forward order.
//    minMutAge minimal generation number. The process will restart if the trajectory
//           is less than minGen. Default to 0.
//    maxMutAge maximum generation number. The process will terminate or restart if it
//           can not reach allele zero after T generations. Default to 100,000,
//           roughly 2,000,000 years which is longer than human history.
//    restartIfFail  If the process can not finish after T generations, restart if
//           restartIfFail=true, otherwise return. Default to false.
//    maxAttempts    How many times to try to get a valid path? Default 1000
//    allowFixation  return the trajectory, if fixation instead of absorption is
//           achieve. Default to false. This option has a side effect to consider
//           0, 2, .. as valid trajectory, and should not be used for simulations.
//
// Of course, you should specify only one of N/NtFunc and one of s/sFunc
//
// Tracking the allele frequency of allele a.
//
//
vectorf FreqTrajectoryStoch(ULONG curGen = 0,
			    double freq = 0,
			    long N = 0,
			    PyObject * NtFunc = NULL,
			    vectorf fitness = vectorf(),
			    PyObject * fitnessFunc = NULL,
			    ULONG minMutAge = 0,
			    ULONG maxMutAge = 100000,
			    int ploidy = 2,
			    bool restartIfFail = false,
			    long maxAttempts = 1000,
			    bool allowFixation = false);

/// CPPONLY a utility function to get marginal fitness given interaction and allele freq
vectorf MarginalFitness(unsigned nLoci, const vectorf & fitness, const vectorf & freq);

//
// simulate trajectories of disease susceptibility loci using an extension of
// the backward method described in Slatkin 2001.
//
// Parameters:
//
//    curGen current generation number
//    N      constant population size N
//    NtFunc a python function that returns population size at each generation.
//           gen is defined in forwrd order. NtFunc(curGen) should be current
//           generation number. NtFunc should return an array of size of subpops.
//    freq   expected allele frequencies of alleles of multiple unlinked loci
//    fitness  constant fitness for [AA, Aa, aa, BB, Bb, bb ...]
//    fitnessFunc  a python function that returns selection pressure at each generation
//           the function expects parameters gen and freq. gen is current generation
//           number and freq is the allele frequency at all loci. This allows
//           frequency dependent selection. gen is defined in forward order.
//    minMutAge minimal generation number. The process will restart if the trajectory
//           is less than minGen. Default to 0.
//    maxMutAge maximum generation number. The process will terminate or restart if it
//           can not reach allele zero after T generations. Default to 100,000,
//           roughly 2,000,000 years which is longer than human history.
//    restartIfFail  If the process can not finish after T generations, restart if
//           restartIfFail=true, otherwise return. Default to false.
//    ploidy    Number of chromosomes will be N*ploidy
//    maxAttempts    How many times to try to get a valid path? Default 10000
//
// Of course, you should specify only one of N/NtFunc and one of s/sFunc
//
// Tracking the allele frequency of allele a.
//
//
matrix FreqTrajectoryMultiStoch(
				ULONG curGen = 0,
				vectorf freq = vectorf(),
				long N = 0,
				PyObject * NtFunc = NULL,
				vectorf fitness = vectorf(),
				PyObject * fitnessFunc = NULL,
				ULONG minMutAge = 0,
				ULONG maxMutAge = 100000,
				int ploidy = 2,
				bool restartIfFail = false,
				long maxAttempts = 1000);

#ifndef OPTIMIZED
// simulate trajectory
vectorf FreqTrajectorySelSim(
			     double sel,                                                        // strength of selection coef  ::8
			     long Ne,                                                           // effective population size ::9
			     double freq,                                                       // initial freq ::10
			     double dom_h,                                                      // strength of dominance ::27
			     int selection                                                      // selection ::5
			     );

/*simulate the sample path of the frequency of disease allele,
   conditional on non-extinction and non-fixation
   array of the mutant allele frequency is backword, start from present-day,
   end at the founding time of the disease
   But, simulation is forward, start from the past when only one copy
   of disease allele, until the
   present-day (if #disease allele <5, poisson dis. conditional
   on #disease >0, else, normal dis.)*/
vectorf FreqTrajectoryForward(double lowbound, double highbound,
			      int disAge, double grate, long N0, double seleCo);
#endif
}
#endif
