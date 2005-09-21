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
simu = simulator(
    population(size=1000, ploidy=2, loci=[1]),
    randomMating() )

# initial allele frequency
count =  alleleCounter(
      alleles=[[0,1,2]],
      genotypes=[[0,1,1,1,2,2,2]]) 
simu.apply([
    initByFreq(alleleFreq=[.2,.8]),
    count
    ])
#
listVars(1)

# genotype frequency does not change during evolution
#
visual = matlabPlotter(
    varDataSource("%gen%genotypeFreq0"),
    win=300,
    update=5)
visual.addPostPlotCmd("ylim([0,1.3])")
simu.setGen(0)
simu.evolve([
    count,
    visual
    ], end=100)


# now, we add selection
#
sel = basicSelector(locus=0, fitness={101:0.6, 102:1, 202:.8})

# now, mating will do the selection.
simu.evolve([
    count,
    sel,
    visual    
    ], end=300)


#
# now, let us study the change of allele frequency
# under some different selection model
#
#
# 1. directional selection
#   w11 > w12 > w22
#  p -> 1
#
sel = basicSelector(locus=0, fitness={101:1, 102:0.9, 202:.8})
visual = matlabPlotter(
    varDataSource("%gen%alleleFreq0[1]"),
    win=100,
    update=5)
visual.addPostPlotCmd("ylim([0,1.3])")
simu.setGen(0)
simu.evolve([
    count,
    sel,
    visual
    ],
    preop=[  initByFreq(alleleFreq=[.2,.8])],
    end=100)

#
# 2. heterozygote superiority
#   w11 < w12, w12 > w22
#  stable. 
# let
#    s1 = w12-  w11
#    s2 = w12 - w22
#  p_ = s2/ (s1+s2)
#
s1 = .1
s2 = .2
p = .2/ (.1+.2)
sel = basicSelector(locus=0, fitness={101:(1-s1), 102:1, 202:(1-s2)})
visual = matlabPlotter(
    varDataSource("%gen%alleleFreq0[1]"),
    win=500,
    update=50)
visual.addPostPlotCmd("""
ylim([0,1.3])
hold on
x=[1,500]
y=[%f %f]
plot(x,y)
hold off
""" % (p,  p) )
simu.setGen(0)
simu.evolve([
    count,
    sel,
    visual
    ],
    preop=[  initByFreq(alleleFreq=[.2,.8])],
    end=500)

#
# 2. heterozygote inferiority
#   w11 > w12, w12 < w22
#  p unstable, become fix
#
sel = basicSelector(locus=0, fitness={101:1, 102:.9, 202:1})
visual = matlabPlotter(
    varDataSource("%gen%alleleFreq0[1]"),
    win=100,
    update=10)
visual.addPostPlotCmd("ylim([0,1.3])")
simu.setGen(0)
simu.evolve([
    count,
    sel,
    visual
    ],
    preop=[  initByFreq(alleleFreq=[.2,.8])],
    end=100)



#
# 3. dominance 
#   1, 1-hs, 1-s
#  stable. 
# let
#  p_ = (h-1)/(2h-1)
#
#  h < 0, underdominance
#  h > 1, overdominance
#  0<h<1, degree of dominance
#
h = 1.2
s = 0.2
p = (h-1)/ (2*h-1)  
sel = basicSelector(locus=0, fitness={101:1, 102:(1-h*s), 202:(1-s)})
visual = matlabPlotter(
    varDataSource("%gen%alleleFreq0[1]"),
    win=500,
    update=50)
visual.addPostPlotCmd("""
ylim([0,1.3])
hold on
x=[1,500]
y=[%f %f]
plot(x,y)
hold off
""" % (p,  p) )
simu.setGen(0)
simu.evolve([
    count,
    sel,
    visual
    ],
    preop=[  initByFreq(alleleFreq=[.3,.7])],
    end=500)


# testing other mating type and parameter
#
#  
simu = simulator(
    population(size=1000, ploidy=2, loci=[1]),
    binomialSelection())

# initial allele frequency
count =  alleleCounter(
      alleles=[[0,1,2]],
      genotypes=[[0,1,1,1,2,2,2]]) 
sel = basicSelector(locus=0, fitness={101:0.6, 102:1, 202:.8},
                    begin=100)
visual = matlabPlotter(
    varDataSource("%gen%genotypeFreq0"),
    win=300,
    update=5)
visual.addPostPlotCmd("ylim([0,1.3])")

simu.setGen(0)
simu.evolve([
    count,
    visual,
    sel
    ],
    preop= [initByFreq(alleleFreq=[.2,.8])],
    end=300)

# because there is no exchange of chromosome,
# no 11 will be generated if lost,
# finally 12/21 will become fix
#

#
#
# random mating with subPopulations, multiple offspring
#
simu = simulator(
    population(size=2000, ploidy=2, loci=[1], subPop=[1000]*2),
    randomMating(numOffsprings=4))

# initial allele frequency
count =  alleleCounter(
      alleles=[[0,1,2]],
      genotypes=[[0,1,1,1,2,2,2]]) 
sel = basicSelector(locus=0, fitness={101:0.6, 102:1, 202:.8},
                    begin=100)
visual = matlabPlotter(
    varDataSource("%gen%genotypeFreq0_sp0[102]%genotypeFreq0_sp1[102]"),
    win=300, rep=-1,
    update=5)
visual.addPostPlotCmd("ylim([0,1.3])")

simu.setGen(0)
simu.evolve([
    count,
    sel,
    visual,
    ],
    preop= [initByFreq(alleleFreq=[.2,.8])],
    end=10)


