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
    from simuPOP.plotter import VarPlotter
except:
    print "simuRPy import failed. Please check your rpy installation."
    useRPy = False
else:
    useRPy = True
        
def I2k(pop, param):
    #stat(pop, alleleFreq=loci)
    I2k ={}
    for loc in param:
        #print pop.dvars().alleleFreq[loc].values()[0]
        I2k[loc] = sum([x*x for x in pop.dvars().alleleFreq[loc].values()])
        #print I2k
    pop.dvars().I2k = I2k
    return True            

def simulate(r, eta, size=30000, numLoci=2, gen=500):
    pop = Population(size=size, loci=[numLoci], infoFields='fitness')

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


    selLocus = 1

    if useRPy:
        plotter = VarPlotter('alleleFreq[%d][1]' % selLocus,
            ylim=[0, 1], update=10, legend=['Allele freq'],
	    xlab="genenerations", ylab="alleleFrequency", saveAs='selectionsweep.png')
    else:
        plotter = NoneOp()
        
    a=[]

    g = pop.evolve(
        initOps = InitSex(),
        preOps = [
            Stat (alleleFreq=[selLocus]),
            IfElse('alleleFreq[%d][1] == 0' % selLocus, 
                ifOps=[
                    PyEval (r'"introduce at gen %d\n"% gen'),
                    PointMutator(inds=0, loci=selLocus, allele=1, ploidy=0, at=0),
                ]),
            MaSelector (loci = selLocus, fitness = [1, 1.+2*eta*0.05, 1+2*0.05]),
            Stat (alleleFreq=[selLocus]),
        ],
        matingScheme = RandomMating(ops=Recombinator(rates=r)),
        postOps = [
            #Dumper(),
            Stat(alleleFreq=[0,1]),
            PyOperator(func = I2k, param=[0]),
            #pyEval(r"'Gen %d: allele frequency: %.3e I2k: %.2e\n' % (gen, alleleFreq[1][1], I2k[0])"),
            TerminateIf('len(alleleFreq[1]) == 1'),
            #plotter,
            #Pause(at=gen-1)
        ],
        gen = gen
    )
    return pop.dvars().I2k[0], pop.dvars().alleleFreq[1][1], g

if __name__ == '__main__':
    #resFile = open('result.txt', 'w')
    #for r, eta in zip((0.5, 0.01, 0.0033, 0.001, 0.001, 0.001, 0.001),
    #             (0.5, 0.5, 0.5, 0.02, 0.1, 0.5, 0.9)):
    #    for rep in range(10):
    #        res = simulate(r, eta, size=10000, numLoci=2, gen=2000)
    #        print >> resFile, r, eta, res, rep
    simulate(0.5, 0.5, size=10000, numLoci=2, gen=200)
