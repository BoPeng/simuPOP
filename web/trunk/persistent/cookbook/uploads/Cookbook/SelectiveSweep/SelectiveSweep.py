#
# Selection Sweep Model
#

"""
This program attempts to model selective sweep.
"""
from simuOpt import setOptions
setOptions(alleleType='long')

from simuPOP import *

try:
    from simuRPy import *
except:
    print "simuRPy import failed. Please check your rpy installation."
    useRPy = False
else:
    useRPy = True
        
def I2k(pop, loci):
    #Stat(pop, alleleFreq=loci)
    I2k ={}
    for loc in loci:
        #print pop.dvars().alleleFreq[loc].values()[0]
        I2k[loc] = sum([x*x for x in pop.dvars().alleleFreq[loc].values()])
        #print I2k
    pop.dvars().I2k = I2k
    return True            

def simulate(size=30000, numLoci=2, gen=500):
    pop = population(size=size, loci=[numLoci], infoFields='fitness')

    for i,ind in enumerate(pop.individuals()):
        if i < 5:
            ind.setAllele(0, 0, ploidy=0)
            ind.setAllele(0, 0, ploidy=1)
            ind.setAllele(1, 1, ploidy=0)
            ind.setAllele(1, 1, ploidy=1)
        else:
            ind.setAllele(2*i, 0, ploidy=0)
            ind.setAllele(2*i+1, 0, ploidy=1)
            ind.setAllele(0, 1, ploidy=0)
            ind.setAllele(0, 1, ploidy=1)

    simu = simulator(pop, randomMating())

    selLocus = 1

    if useRPy:
        plotter = varPlotter('alleleFreq[%d][1]' % selLocus,
            ylim=[0, 1], update=10, legend=['Allele freq'],
	    xlab="genenerations", ylab="alleleFrequency", saveAs='selectionsweep.png')
    else:
        plotter = noneOp()
        
    a=[]

    g = simu.evolve(
    preOps = [
        initSex()
        ],
    ops = [
        stat (alleleFreq=[selLocus], stage=PreMating),
        ifElse('alleleFreq[%d][1] == 0' % selLocus, stage=PreMating,
            ifOps=[
            pyEval (r'"introduce at gen %d\n"% gen'),
            pointMutator(inds=0, loci=selLocus, allele=1, ploidy=0, at=0, stage=PreMating),
            ]),
        
        recombinator(rates=r),

        #dumper(stage=PreMating),
        stat (alleleFreq=[selLocus], stage=PreMating),
        maSelector (loci = selLocus, fitness = [1, 1.+2*eta*0.05, 1+2*0.05]),
        #dumper(),
        stat(alleleFreq=[0,1]),
        pyOperator(func = I2k, param=[0]),
        #pyEval(r"'Gen %d: allele frequency: %.3e I2k: %.2e\n' % (gen, alleleFreq[1][1], I2k[0])"),
        terminateIf('len(alleleFreq[1]) == 1'),
        #plotter,
        #pause(at=gen-1)
        ],

    #dryrun=True,
    gen = gen
    )
    
    
    return simu.dvars(0).I2k[0], simu.dvars(0).alleleFreq[1][1], g[0]

resFile = open('result.txt', 'w')

for r, eta in zip((0.5, 0.01, 0.0033, 0.001, 0.001, 0.001, 0.001),
                 (0.5, 0.5, 0.5, 0.02, 0.1, 0.5, 0.9)):
  for rep in range(10):
      res = simulate(size=10000, numLoci=2, gen=2000)
      print >> resFile, r, eta, res, rep