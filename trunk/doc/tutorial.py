#!/usr/bin/env python
#
# python scripts for tutorial.lyx
#
from simuOpt import setOptions
setOptions(quiet=True)
#file log/tutorial_example1.log
from simuPOP import *
simu = simulator(
    population(size=1000, ploidy=2, loci=[2]),
    randomMating(),
    rep = 3)
simu.evolve(
    preOps = [initByValue([1,2,2,1])],  
    ops = [
        recombinator(rate=0.1),
        stat(LD=[0,1]),
        pyEval(r"'%3d   ' % gen", rep=0, step=10),
        pyEval(r"'%f    ' % LD[0][1]", step=10),
        pyEval(r"'\n'", rep=REP_LAST, step=10)
    ],
    end=100
)
#end
# #file log/tutorial_example2.log
# from simuPOP import *
# from simuRPy import *
# simu = simulator(
#     population(size=1000, ploidy=2, loci=[2]),
#     randomMating(),
#     rep = 3)
# simu.evolve(
#     preOps = [initByValue([1,2,2,1])],  
#     ops = [
#         recombinator(rate=0.1),
#         stat(LD=[0,1]),
#         varPlotter('LD[0][1]', numRep=3, step=10, 
#             saveAs='ld', ylim=[0,.25], 
#             lty=range(1, 4), col=range(2, 5),
#             xlab='generation', ylab='D', 
#             title='LD Decay'),
#     ],
#     end=100
# )
# #end
# #PS epstopdf ld10.eps
# #PS epstopdf ld20.eps
# #PS epstopdf ld30.eps
# #PS epstopdf ld40.eps
# #PS epstopdf ld50.eps
# #PS epstopdf ld60.eps
# #PS epstopdf ld70.eps
# #PS epstopdf ld80.eps
# #PS epstopdf ld90.eps
# #PS epstopdf ld100.eps
#file log/tutorial_popsize.log
def lin_inc(gen, oldsize=[]):
    return [10+gen]*5

simu = simulator(
    population(subPop=[5]*5, loci=[1]),
    randomMating(newSubPopSizeFunc=lin_inc)
)
simu.evolve(
    ops = [
        stat(popSize=True),
        pyEval(r'"%d %d\n"%(gen, subPop[0]["popSize"])'),
    ],
    end=5
)
#end
#file log/tutorial_stat.log
simu = simulator(
    population(subPop=[10000]*2, loci=[10]),
    randomMating()
)
simu.evolve(
    preOps = [
        initByFreq([0.2, 0.8], subPop=[0]),
        initByFreq([0.8, 0.2], subPop=[1]),
    ],
    ops = [
        stat(LD=[[0,1], [5,6]], Fst=range(10), step=100),
        migrator(rate=[[0, 0.01], [0, 0.02]]),
        pyEval(r'"Gen: %4d LD: %.3f R2: %.3f Fst: %.3f\n"' \
            ' % (gen, LD[0][1], R2[0][1], AvgFst)', 
            step=100)
    ],
    end=1000
)
#end
#file log/tutorial_pyPenetrance.log
def myPene(geno):
    'geno is the genotype at the two given loci'
    loc1 = geno[0] + geno[1]
    loc2 = geno[2] + geno[3]
    if (loc1 == 2 and loc2 < 2) or \
        (loc1 < 2 and loc2 == 2):
        return 0.1
    else:
        return 0.5

pop = population(subPop=[1000], loci=[6])
# initialize the population
InitByFreq(pop, [0.1, 0.9])
# apply penetrance and obtain affection status
PyPenetrance(pop, loci=[3, 5], func=myPene)
# draw case control sample
(sample,) = CaseControlSample(pop, cases=3, controls=3)
# save sample in Merlin QTDT format
from simuUtil import SaveQTDT
SaveQTDT(sample, output='sample', affectionCode=['U', 'A'], 
    fields=['affection'])
# have a look at the sample in Merlin-QTDT Format
print open('sample.map').read()
print open('sample.dat').read()
print open('sample.ped').read()
#end
#file log/tutorial_eff_alleles.log
def Ne(pop, loci):
    'Calculate effective number of alleles'
    Stat(pop, alleleFreq=loci)
    pop.dvars().Ne = {}
    v = pop.dvars().alleleFreq
    for locus in loci:
        f0 = 1 - v[locus][0]
        Ne = 1./f0*f0/sum([x*x for x in v[locus][1:]])
        pop.dvars().Ne[locus] = Ne
    return True

pop = population(1000, loci=[2,3])
      
#end
#file log/tutorial_hapmap.log
genes = [
    "rs1042522",
    "rs1625895",
    "rs1799793",
]
pops = []
for i in range(1, 23):
    print "Loading hapmap chromosome %d..." % i
    pop = LoadPopulation('hapmap_%d.bin' % i)
    markers = []
    for name in genes:
        try:
            idx = pop.locusByName(name)
            markers.append(idx)
        except:
            pass
    if len(markers) > 0:
        markers.sort()
        pop.removeLoci(keep=markers)
        pops.append(pop)

all = MergePopulationsByLoci(pops)
#end
