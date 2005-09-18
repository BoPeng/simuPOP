#!/usr/bin/env python
#
# Purpose:
#   selection.
# 
# Author:
#   Bo Peng (bpeng@rice.edu)
#
# load module.
from simuPOP import *

#
# 1. basic selection model
#
#    specify relative fitness: w11, w12/w21, w22
#
#

# random mating with subPopulations, multiple offspring
setArrayVar("popSize", [1000,1000])
simu = simulator(
    population(size=2000, ploidy=2, loci=[1], subPop=[1000,1000]),
    randomMating(numOffsprings=4, newSubPopSizeExpr="%popSize=%popSize*1.02"))

# initial allele frequency
count =  alleleCounter(
      alleles=[[0,1,2]],
      genotypes=[[0,1,1,1,2,2,2]]) 
sel = basicSelector(locus=0, fitness={101:0.6, 102:1, 202:.8},
                    begin=100)
visual = matlabPlotter(
    varDataSource("%gen%genotypeFreq0_sp0[102]%genotypeFereq0_sp1[102]"),
    win=300, rep=-1,
    update=5)
visual.addPostPlotCmd("ylim([0,1.3])")

simu.setGen(0)
simu.evolve([
    count,
    sel,
    visual,
    ],
    preop= [initByFreq(alleleFreq=[.2])],
    end=10)


