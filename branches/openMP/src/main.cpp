/**
 *  $File: main.cpp $
 *  $LastChangedDate: 2010-12-03 16:49:34 -0600 (Fri, 03 Dec 2010) $
 *  $Rev: 3932 $
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

/**
 * \file main.cpp
 *
 * Because it is easier to trace or debug simuPOP as a regular single
 * executable, simuPOP can be compiled with the main() function provided in
 * this file. However, because a Python engine is unavailable, not all
 * features of simuPOP can be used in this way and relevant source code
 * should be commented out with a
 *    #ifndef STANDALINE_EXECUTABLE
 * macro definition.
 *
 * To compile single executable, please use
 *     scons debug
 * to generate build/MOD/simuPOP_MOD, or
 *     scons debug MOD
 * to generate executable for selected module.
 *
 * This file is NOT used when simuPOP is compiled as python modules.
 */

#ifdef STANDALONE_EXECUTABLE

//#include "config.h"
//#include "simuPOP_cfg.h"
#include "genoStru.h"
#include "individual.h"
#include "population.h"
#include "pedigree.h"
#include "virtualSubPop.h"
#include "operator.h"
#include "simulator.h"

#include "utility.h"
#include "pedigree.h"
#include "initializer.h"
#include "outputer.h"
#include "mating.h"
#include "tagger.h"
#include "stator.h"
#include "migrator.h"
#include "mutator.h"
#include "transmitter.h"
#include "selector.h"
#include "qtrait.h"
#include "penetrance.h"
#include "sandbox.h"

using namespace simuPOP;
using namespace std;

void banner()
{
	cout	<< "Standalone "
#ifdef OPTIMIZED
	        << "optimized "
#else
	        << "standard "
#endif
	// AlleleType
#ifdef LONGALLELE
	<< "long"
#else
#  ifdef BINARYALLELE
	<< "binary"
#  else
	<< "short"
#  endif
#endif
	<< " version of simuPOP for debugging purposes" << endl;
}


// a simple version of random mating
MatingScheme * RandomMating(int numOffspring = 1, int sexMode = RANDOM_SEX,
                            const vectoru & subPopSize = vectoru())
{
	floatListFunc no(numOffspring);
	floatListFunc sm(sexMode);
	BaseOperator * mendelian = new MendelianGenoTransmitter();
	opList duringOps(vectorop(1, mendelian));
	RandomParentsChooser pc(true, "fitness");
	OffspringGenerator og(duringOps, no, sm);
	uintListFunc sp(subPopSize);
	delete mendelian;
	return new HomoMating(pc, og, sp);
}


bool basicRandomMating(UINT size = 1000, UINT loci = 10, UINT gen = 10)
{
	/*
	   pop = Population(1000, loci=10)
	   pop.evolve(
	    initOps=InitSex(),
	    matingScheme=RandomMating(),
	    gen=10
	   )
	 */
	// create the population and simulator
	uintList sz(vectoru(1, size));
	uintList lo(vectoru(1, loci));
	Population pop(sz, 2, lo);
	Simulator sim(NULL);

	sim.add(pop);
	// operators
	BaseOperator * iSex = new InitSex();
	opList initOps(vectorop(1, iSex));
	opList preOps;
	opList postOps;
	opList finalOps;
	// mating scheme
	MatingScheme * ms = RandomMating();
	// evolve
	cout << "Begin evolving basicRandomMating with size " << size << endl;
	try {
		sim.evolve(initOps, preOps, *ms, postOps, finalOps, gen);
	} catch (Exception e) {
		cerr << e.message() << endl;
	}
	cout << "Done" << endl;

	delete iSex;
	delete ms;
	return true;
}


int main(int argc, char ** argv)
{

	if (argc != 5) {
		printf("usage: %s <num_threads> <size> <loci> <gen> \n", argv[0]);
		exit(0);
	}
	int nth = atoi(argv[1]);
	int size = atoi(argv[2]);
	int loci = atoi(argv[3]);
	int gen = atoi(argv[4]);

#ifdef _OPENMP
	omp_set_num_threads(nth);
	initialize();
#endif

	banner();
	basicRandomMating(size, loci, gen);
	return 0;
}


#endif
