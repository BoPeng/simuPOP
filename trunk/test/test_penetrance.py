# test for penetrance
#
# create population, 400 loci on 20 chromosomes,
## from simuPOP import *
## from simuUtil import *
## #from simuRPy import *
## 
## pop = population(size=10, ploidy=2, loci=[20]*20,
##   alleleNames=['_','A','T','C','G'], maxAllele=4)
## simu = simulator(pop, randomMating(), rep=5)
## # evolve. 
## simu.evolve( 
##   preOps = [ initByFreq([.25 ]*4) ],
##   ops = [
##     # Juke-cantor model (4-allels model)
##     # independent on all loci
##     kamMutator(rate=0.001),
##     # recombination, uniform
##     recombinator(rate=0.0001),
##     # selection on a recessive disease at allele 30
##     # assume that allele 1 (A) is deleterious
##     maSelector(locus=30, wildtype=[2,3,4], fitness=[1,1,.9]),
##     # penetrance
##     maPenetrance(locus=30, wildtype=[2,3,4], penetrance=[0,0,.8]),
##     # count diseased individuals
##     stat(numOfAffected=1),
##     # output
##     pyEval(r'"%d\t" % numOfAffected'),
##     endl(rep=REP_LAST),
## ##     # sample at the last generation
##     caseControlSample(cases=50, controls=100, at=[-1])
##     ],
##   end=20
##   )
##     
