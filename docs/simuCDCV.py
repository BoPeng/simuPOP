#!/usr/bin/env python

#
# $File: simuCDCV.py $
#
# This file is part of simuPOP, a forward-time population genetics
# simulation environment. Please visit http://simupop.sourceforge.net
# for details.
#
# Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

# This script is an example in the simuPOP user's guide. Please refer to
# the user's guide (http://simupop.sourceforge.net/manual) for a detailed
# description of this example.
#

#!/usr/bin/env python
#
# Author:  Bo Peng
# Purpose: A real world example for simuPOP user's guide.
#
'''
Simulation the evolution of allelic spectra (number and frequencies
of alleles at a locus), under the influence of sim.population expansion,
mutation, and natural selection.
'''
import simuOpt
simuOpt.setOptions(quiet=True, alleleType='long')
import simuPOP as sim
import sys, types, os, math
options = [
    {'name': 'demo',
     'default': 'instant',
     'label': 'Population growth model',
     'description': 'How does a sim.Population grow from N0 to N1.',
     'type': ('chooseOneOf', ['instant', 'exponential']),
    },
    {'name': 'N0',
     'default': 10000,
     'label': 'Initial sim.population size',
     'type': 'integer',
     'description': '''Initial sim.population size. This size will be maintained
                till the end of burnin stage''',
     'validator': simuOpt.valueGT(0)
    },
    {'name': 'N1',
     'default': 100000,
     'label': 'Final sim.population size',
     'type': 'integer',
     'description': 'Ending sim.population size (after sim.population expansion)',
     'validator': simuOpt.valueGT(0)
    }, 
    {'name': 'G0',
     'default': 500,
     'label': 'Length of burn-in stage',
     'type': 'integer',
     'description': 'Number of generations of the burn in stage.',
     'validator': simuOpt.valueGT(0)
    },
    {'name': 'G1',
     'default': 1000,
     'label': 'Length of expansion stage',
     'type': 'integer',
     'description': 'Number of geneartions of the sim.population expansion stage',
     'validator': simuOpt.valueGT(0)
    },
    {'name': 'spec',
     'default': [0.9] + [0.02]*5,
     'label': 'Initial allelic spectrum',
     'type': 'numbers',
     'description': '''Initial allelic spectrum, should be a list of allele
            frequencies, for allele 0, 1, 2, ... respectively.''',
     'validator': simuOpt.valueListOf(simuOpt.valueBetween(0, 1)),
    },
    {'name': 's',
     'default': 0.01,
     'label': 'Selection pressure',
     'type': 'number',
     'description': '''Selection coefficient for homozygtes (aa) genotype.
            A recessive selection model is used so the fitness values of
            genotypes AA, Aa and aa are 1, 1 and 1-s respectively.''',
     'validator': simuOpt.valueGT(-1),
    },
    {'name': 'mu',
     'default': 1e-4,
     'label': 'Mutation rate',
     'type': 'number',
     'description': 'Mutation rate of a k-allele mutation model',
     'validator': simuOpt.valueBetween(0, 1),
    },
    {'name': 'k',
     'default': 200,
     'label': 'Maximum allelic state',
     'type': 'integer',
     'description': 'Maximum allelic state for a k-allele mutation model',
     'validator': simuOpt.valueGT(1),
    },
]

def demo_model(type, N0=1000, N1=100000, G0=500, G1=500):
    '''Return a demographic function 
    type: linear or exponential
    N0:   Initial sim.population size.
    N1:   Ending sim.population size.
    G0:   Length of burn-in stage.
    G1:   Length of sim.population expansion stage.
    '''
    rate = (math.log(N1) - math.log(N0))/G1
    def ins_expansion(gen):
        if gen < G0:
            return N0
        else:
            return N1
    
    def exp_expansion(gen):
        if gen < G0:
            return N0
        else:            
            return int(N0 * math.exp((gen - G0) * rate))
    
    if type == 'instant':
        return ins_expansion
    elif type == 'exponential':
        return exp_expansion

class ne(sim.PyOperator):
    '''Define an operator that calculates effective number of
    alleles at given loci. The result is saved in a population
    variable ne.
    '''
    def __init__(self, loci, *args, **kwargs):
        self.loci = loci
        sim.PyOperator.__init__(self, func=self.calcNe, *args, **kwargs)
    
    def calcNe(self, pop):
        sim.stat(pop, alleleFreq=self.loci)
        ne = {}
        for loc in self.loci:
            freq = pop.dvars().alleleFreq[loc]
            sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
            if sumFreq == 0:
                ne[loc] = 0
            else:
                ne[loc] = 1. / sum([(freq[x]/sumFreq)**2 for x in list(freq.keys()) if x != 0])
        # save the result to the sim.Population.
        pop.dvars().ne = ne
        return True

def simuCDCV(model, N0, N1, G0, G1, spec, s, mu, k):
    '''Evolve a sim.Population using given demographic model
    and observe the evolution of its allelic spectrum.
    model: type of demographic model.
    N0, N1, G0, G1: parameters of demographic model.
    spec: initial allelic spectrum, should be a list of allele
        frequencies for each allele.
    s: selection pressure.
    mu: mutation rate.
    k: maximum allele for the k-allele model
    '''
    demo_func = demo_model(model, N0, N1, G0, G1)
    print(demo_func(0))
    pop = sim.Population(size=demo_func(0), loci=1, infoFields='fitness')
    pop.evolve(
        initOps=[
            sim.InitSex(),
            sim.InitGenotype(freq=spec, loci=0)
        ],
        matingScheme=sim.RandomMating(subPopSize=demo_func),
        postOps=[
            sim.KAlleleMutator(rates=mu, k=k),
            sim.MaSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
            ne(loci=(0,), step=100),
            sim.PyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
                step=100),
        ],
        gen = G0 + G1
    )
    return pop

if __name__ == '__main__':
    # get parameters
    par = simuOpt.Params(options, __doc__)
    if not par.getParam():
        sys.exit(1)
    
    if not sum(par.spec) == 1:
        print('Initial allelic spectrum should add up to 1.')
        sys.exit(1)
    # save user input to a configuration file
    par.saveConfig('simuCDCV.cfg')
    #
    simuCDCV(*par.asList())


