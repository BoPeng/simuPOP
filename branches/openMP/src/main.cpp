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
 * This file is used for debugging only. Basically, it compiles with
 * the C++ version of simuPOP implement a test script in C/C++. To use this
 * script, please execute 'make test_simuPOP' under linux or create a 
 * VC project file in a similar fashion.
 */

#include "config.h"
#include "simuPOP_cfg.h"
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

int main()
{
	uintList sz(vectoru(1, 100));
	Population pop(sz, 2);
	Simulator sim(NULL);
	sim.add(pop);
	BaseOperator * iSex = new InitSex();
	BaseOperator * mendelian = new MendelianGenoTransmitter();
	opList initOps(vectorop(1, iSex));
	opList preOps;
	opList postOps;
	opList duringOps(vectorop(1, mendelian));
	opList finalOps;
	RandomParentsChooser pc;
	floatListFunc numOff(1);
	floatListFunc sexMode(RANDOM_SEX);
	OffspringGenerator og(duringOps, numOff, sexMode);
	MatingScheme ms = HomoMating(pc, og);
	sim.evolve(initOps, preOps, ms, postOps, finalOps, 1);
	return 0;
}


