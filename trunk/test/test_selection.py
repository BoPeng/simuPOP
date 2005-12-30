## #!/usr/bin/env python
## #
## # Purpose:
## #   selection.
## # 
## # Author:
## #   Bo Peng (bpeng@rice.edu)
## #
## # load module.
## from simuPOP import *
## from simuUtil import *
## from simuRPy import *
## 
## #
## # 1. map selection model
## #
## #    specify relative fitness: w11, w12/w21, w22
## #
## #  
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[1]),
##     randomMating() )
## 
## # initial allele frequency
## count =  stat(
##       alleleFreq=[0],
##       genoFreq=[0]) 
## simu.apply([
##     initByFreq(alleleFreq=[.2,.8]),
##     count
##     ])
## #
## listVars(simu.dvars(0))
## visual = varPlotter('alleleFreq[0]', ylim=[0,1], varDim=3)  
## # genotype frequency does not change during evolution
## #
## simu.setGen(0)
## simu.evolve([
##     count,
##     visual
##     ], end=100)
## 
## 
## # now, we add selection
## #
## sel = mapSelector(locus=0, fitness={'1-1':0.6, '1-2':1, '2-2':.8})
## 
## # now, mating will do the selection.
## simu.evolve([
##     count,
##     sel,
##     visual    
##     ], end=300)
## 
## 
## #
## # now, let us study the change of allele frequency
## # under some different selection model
## #
## #
## # 1. directional selection
## #   w11 > w12 > w22
## #  p -> 1
## #
## 
## visual = varPlotter('alleleFreq[0]', ylim=[0,1], varDim=3)  
## sel = mapSelector(locus=0, fitness={'1-1':1, '1-2':0.9, '2-2':.8})
## simu.setGen(0)
## simu.evolve([
##     count,
##     sel,
##     visual
##     ],
##     preOps=[  initByFreq(alleleFreq=[.2,.8])],
##     end=100)
## 
## #
## # 2. heterozygote superiority
## #   w11 < w12, w12 > w22
## #  stable. 
## # let
## #    s1 = w12-  w11
## #    s2 = w12 - w22
## #  p_ = s2/ (s1+s2)
## #
## s1 = .1
## s2 = .2
## p = .2/ (.1+.2)
## visual = varPlotter('alleleFreq[0]', ylim=[0,1], varDim=3, update=100)  
## sel = mapSelector(locus=0,
##   fitness={'1-1':(1-s1), '1-2':1, '2-2':(1-s2)})
## simu.setGen(0)
## simu.evolve([
##     count,
##     sel,
##     visual
##     ],
##     preOps=[  initByFreq(alleleFreq=[.2,.8])],
##     end=500)
## 
## #
## # 2. heterozygote inferiority
## #   w11 > w12, w12 < w22
## #  p unstable, become fix
## #
## visual = varPlotter('alleleFreq[0]', ylim=[0,1], varDim=3, update=10)  
## sel = mapSelector(locus=0, fitness={'1-1':1, '1-2':.9, '2-2':1})
## simu.setGen(0)
## simu.evolve([
##     count,
##     sel,
##     visual
##     ],
##     preOps=[  initByFreq(alleleFreq=[.2,.8])],
##     end=100)
## 
## 
## 
## #
## # 3. dominance 
## #   1, 1-hs, 1-s
## #  stable. 
## # let
## #  p_ = (h-1)/(2h-1)
## #
## #  h < 0, underdominance
## #  h > 1, overdominance
## #  0<h<1, degree of dominance
## #
## h = 1.2
## s = 0.2
## p = (h-1)/ (2*h-1)
## visual = varPlotter('alleleFreq[0]', ylim=[0,1], varDim=3, update=100)  
## sel = mapSelector(locus=0, fitness={'1-1':1, '1-2':(1-h*s), '2-2':(1-s)})
## simu.setGen(0)
## simu.evolve([
##     count,
##     sel,
##     visual
##     ],
##     preOps=[  initByFreq(alleleFreq=[.3,.7])],
##     end=500)
## 
## 
## # testing other mating type and parameter
## #
## #  
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[1]),
##     binomialSelection())
## 
## # initial allele frequency
## count =  stat(
##       alleleFreq=[0],
##       genoFreq=[0])
## visual = varPlotter('alleleFreq[0]', ylim=[0,1], varDim=3, update=100)  
## 
## sel = mapSelector(locus=0, fitness={'1-1':0.6, '1-2':1, '2-2':.8},
##                     begin=100)
## simu.setGen(0)
## simu.evolve([
##     count,
##     visual,
##     sel
##     ],
##     preOps= [initByFreq(alleleFreq=[.2,.8])],
##     end=300)
## 
## # because there is no exchange of chromosome,
## # no 11 will be generated if lost,
## # finally 12/21 will become fix
## #
## 
## #
## #
## # random mating with subPopulations, multiple offspring
## #
## simu = simulator(
##     population(size=2000, ploidy=2, loci=[1], subPop=[1000]*2),
##     randomMating(numOffsprings=4))
## visual = varPlotter('alleleFreq[0]', ylim=[0,1], varDim=3, update=100)  
## # initial allele frequency
## count =  stat(
##       alleleFreq=[0],
##       genoFreq=[0]) 
## sel = mapSelector(locus=0, fitness={'1-1':0.6, '1-2':1, '2-2':.8},
##                     begin=100)
## simu.setGen(0)
## simu.evolve([
##     count,
##     sel,
##     visual,
##     ],
##     preOps= [initByFreq(alleleFreq=[.2,.8])],
##     end=200)
## 
## 
## 
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[1]),
##     randomMating())
## s1 = .1
## s2 = .2
## simu.evolve([
##     stat( alleleFreq=[0], genoFreq=[0]),
##     mapSelector(locus=0, fitness={'1-1':(1-s1), '1-2':1, '2-2':(1-s2)}),
##     pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=10)
##     ],
##     preOps=[  initByFreq(alleleFreq=[.2,.8])],
##     end=300
## )
## 
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[3]),
##     randomMating() )
## 
## s1 = .2
## s2 = .3
## def sel(arr):
##   if arr[0] == 1 and arr[1] == 1:
##     return 1 - s1
##   elif arr[0] == 1 and arr[1] == 2:
##     return 1
##   elif arr[0] == 2 and arr[1] == 1:
##     return 1
##   else:
##     return 1 - s2
## # test func
## print sel(carray('B',[1,1]))
## 
## simu.evolve([
##     stat( alleleFreq=[0], genoFreq=[0]),
##     pySelector(loci=[0,1],func=sel),
##     pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=10)
##     ],
##     preOps=[  initByFreq(alleleFreq=[.2,.8])],
##     end=100)
## 
## 
## # mlSelector
## 
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[2]),
##     randomMating())
## simu.evolve([
##     stat( alleleFreq=[0,1], genoFreq=[0,1]),
##     mlSelector([
##       mapSelector(locus=0, fitness={'1-1':1,'1-2':1,'2-2':.8}),
##       mapSelector(locus=1, fitness={'1-1':1,'1-2':1,'2-2':.8})
##       ], mode=SEL_Additive),
##     pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=10)
##     ],
##     preOps=[  initByFreq(alleleFreq=[.2,.8])],
##  end=100
## )
## 
## 
## 
## simu = simulator(
##     population(size=1000, ploidy=2, loci=[2]),
##     randomMating())
## simu.evolve([
##     stat( alleleFreq=[0,1], genoFreq=[0,1]),
##     mlSelector([
##       mapSelector(locus=0, fitness={'1-1':1,'1-2':1,'2-2':.8}),
##       maSelector(locus=1, wildtype=[1], fitness=[1,1,.8])
##       ], mode=SEL_Additive),
##     pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=10)
##     ],
##     preOps=[  initByFreq(alleleFreq=[.2,.8])],
##  end=100
## )
## 
## 
## 
