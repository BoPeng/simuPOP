# This script will genreate code/result pieces
# that will be inserted into simuPOP user's guide
# and reference manual
# 
# #file_begin file
#
# #file_end file
#
# will be used in this file so the running result can be
# separated into files specified by 'file'
#
# #PS after commands that will be executed.
#
# create directory log if not exist
import os, sys

if not os.path.isdir('log'):
  try:
    os.mkdir('log')
  except:
    print "Failed to make output directory log"
    sys.exit(1)

#

#file log/standard.log
from simuPOP import *
pop = population(10, loci=[2])
pop.locusPos(10)
pop.individual(20).setAllele(1, 0)
#end

# make sure each run generates the same output to avoid unnecessary
# documentation changes.
rng().setSeed(12345)

#file log/importSimuPOP.log
import simuOpt
simuOpt.setOptions(optimized=False, alleleType='long', quiet=False)
from simuPOP import *
#end

#file log/simpleExample.log
from simuPOP import *
pop = population(size=1000, loci=[2])
simu = simulator(pop, randomMating(), rep=3)
simu.evolve(
    preOps = initByValue([1, 2, 2, 1]),  
    ops = [
        recombinator(rates=0.01),
        stat(LD=[0, 1]),
        pyEval(r"'%.2f\t' % LD[0][1]", step=10),
        pyOutput('\n', rep=-1, step=10)
    ],
    gen=100
)
#end

#file log/help.log
help(population.addInfoFields)
#end

#file log/absIndex.log
pop = population(size=[10, 20], loci=[5, 7])
print pop.chromLocusPair(7)
print pop.absLocusIndex(1, 1)
print pop.absIndIndex(2, 1)
print pop.subPopIndPair(25)
#end

#file log/iterator.log
pop = population(size=2, loci=[5, 6])
InitByFreq(pop, [0.2, 0.3, 0.5])
for ind in pop.individuals():
    for loc in range(pop.chromBegin(1), pop.chromEnd(1)):
        print ind.allele(loc),
    print 

#end

#file log/carray.log
pop = population(size=2, loci=[3, 4])
InitByFreq(pop, [0.3, 0.5, 0.2])
ind = pop.individual(0)
arr = ind.genotype()    # a carray to the underlying genotype
geno = list(arr)        # a list of alleles
print arr
print geno
arr.count(1)           # count
arr.index(2)           # index 
ind.setAllele(5, 3)    # change underlying genotype using setAllele
print arr              # arr is change
print geno             # but not geno
arr[2:5] = 4           # can use regular Python slice operation
print ind.genotype()
#end

#file log/genoStru.log
pop = population(size=[2, 3], ploidy=2, loci=[5, 10],
    lociPos=[range(0, 5), range(0, 20, 2)], chromNames=['Chr1', 'Chr2'],
    alleleNames=['A', 'C', 'T', 'G'])
# access genotypic information from the population
pop.ploidy()
pop.ploidyName()
pop.numChrom()
pop.locusPos(2)
pop.alleleName(1)
# access from an individual
ind = pop.individual(2)
ind.numLoci(1)
ind.chromName(0)
ind.locusName(1)
# utility functions
ind.chromBegin(1)
ind.chromByName('Chr2')
# loci pos can be unordered within each chromosome
pop = population(loci=[2, 3], lociPos=[3, 1, 1, 3, 2],
    lociNames=['loc%d' % x for x in range(5)])
pop.lociPos()
pop.lociNames()
#end


#file log/haplodiploid.log
pop = population(size=[2,5], ploidy=Haplodiploid, loci=[3, 5])
InitByFreq(pop, [0.3, 0.7])
Dump(pop)
#end

#file log/chromType.log
pop = population(size=6, ploidy=2, loci=[3, 3, 6, 4, 4, 4],
    chromTypes=[Autosome]*2 + [ChromosomeX, ChromosomeY] + [Customized]*2)
InitByFreq(pop, [0.3, 0.7])
Dump(pop, structure=False) # does not display genotypic structure information
#end

#file log/infoField.log
pop = population(10, loci=[20], ancGen=1,
    infoFields=['father_idx', 'mother_idx'])
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = initByValue([0]*20+[1]*20),
    ops = [
        parentsTagger(),
        recombinator(rates=0.01)
    ],
    gen = 1
)
pop = simu.extract(0)
pop.indInfo('mother_idx')  # mother of all offspring
ind = pop.individual(0)
mom = pop.ancestor(ind.intInfo('mother_idx'), 1)
print ind.genotype(0)
print mom.genotype(0)
print mom.genotype(1)
#end


#file log/individual.log
pop = population([5, 4], loci=[2, 5], infoFields='x')
# get an individual
ind = pop.individual(3)
ind.ploidy()            # access to genotypic structure
ind.numChrom()
ind.affected()
ind.setAffected(True)   # access affection status,
ind.sex()               # sex,
ind.setInfo(4, 'x')     # and information fields
ind.info('x')
#end

#file log/individual_genotype.log
pop = population([2, 1], loci=[2, 5])
for ind in pop.individuals(1):
    for marker in range(pop.totNumLoci()):
        ind.setAllele(marker % 2, marker, 0)
        ind.setAllele(marker % 2, marker, 1)
        print '%d %d ' % (ind.allele(marker, 0), ind.allele(marker, 1))

ind = pop.individual(1)
geno = ind.genotype(1)      # the second homologous copy
geno
geno[2] = 3
ind.genotype(1)
geno[2:4] = [3, 4]          # direct modification of the underlying genotype
ind.genotype(1)
# set genotype (genotype, ploidy, chrom)
ind.setGenotype([2, 1], 1, 1)
geno
#end


#file log/subPop.log
pop = population(size=[3, 4, 5], ploidy=1, loci=[1], infoFields='x')
# individual 0, 1, 2, ... will have an allele 0, 1, 2, ...
pop.setGenotype(range(pop.popSize()))
#
pop.subPopSize(1)
# merge subpopulations
pop.mergeSubPops([1, 2])
# split subpopulations
pop.splitSubPop(1, [2, 7])
pop.subPopSizes()
# set information field to each individual's new subpopulation ID
pop.setIndInfo([0, 1, 2, -1, 0, 1, 2, -1, -1, 0, 1, 2], 'x')
# this manually triggers an migration, individuals with negative values
# at this information field are removed.
pop.setSubPopByIndInfo('x')
Dump(pop, width=2, structure=False)
#end


#file log/subPopName.log
pop = population(size=[3, 4, 5], subPopNames=['x', 'y', 'z'])
pop.removeSubPops([1])
pop.subPopNames()
pop.subPopByName('z')
pop.splitSubPop(1, [2, 3])
pop.subPopNames()
pop.setSubPopName('z-1', 1)
pop.subPopNames()
pop.subPopByName('z')
#end

#file log/virtualSplitter.log
import random
pop = population(size=[200, 400], loci=[30], infoFields='x')
# assign random information fields
pop.setIndInfo([random.randint(0, 3) for x in range(pop.popSize())], 'x')
# define a virtual splitter by information field 'x'
pop.setVirtualSplitter(infoSplitter(field='x', values=[0, 1, 2, 3]))
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 0])    # Each VSP has a name
pop.subPopSize([0, 0])    # Size of VSP 0 in subpopulation 0
pop.subPopSize([1, 0])    # Size of VSP 0 in subpopulation 1
# use a combined splitter that defines additional VSPs by sex
InitSex(pop)
pop.setSubPopName('subPop 1', 0)
pop.setVirtualSplitter(combinedSplitter([
    infoSplitter(field='x', values=[0, 1, 2, 3]),
    sexSplitter()])
)
pop.numVirtualSubPop()    # Number of defined VSPs
pop.subPopName([0, 4])    # VSP 4 is the first VSP defined by the sex splitter
pop.subPopSize([0, 4])    # Number of male individuals
#end


#file log/virtualSubPop.log
import random
pop = population(10, loci=[2, 3], infoFields='Sex')
InitSex(pop)
pop.setVirtualSplitter(sexSplitter())
# initialize male and females with different genotypes. Set initSex
# to False because this operator will by default also initialize sex.
InitByValue(pop, [[0]*5, [1]*5], subPops=([0, 0], [0, 1]), initSex=False)
# set Sex information field to 0 for all males, and 1 for all females
pop.setIndInfo([Male], 'Sex', [0, 0])
pop.setIndInfo([Female], 'Sex', [0, 1])
# Print individual genotypes, followed by values at information field Sex
Dump(pop, structure=False)
#end


#file log/accessIndividual.log
# create a population with two generations. The current generation has values
# 0-9 at information field x, the parental generation has values 10-19.
pop = population(size=[5, 5], loci=[2, 3], infoFields='x', ancGen=1)
pop.setIndInfo(range(11, 20), 'x')
pop1 = pop.clone()
pop1.setIndInfo(range(10), 'x')
pop.push(pop1)
#
ind = pop.individual(5)       # using absolute index
ind.info('x')
# use a for loop, and relative index
for idx in range(pop.subPopSize(1)):
    print pop.individual(idx, 1).info('x'),

# It is usually easier to use an iterator
for ind in pop.individuals(1):
    print ind.info('x'),

# Access individuals in VSPs
pop.setVirtualSplitter(infoSplitter(cutoff=[3, 7], field='x'))
for ind in pop.individuals([1, 1]):
    print ind.info('x'),

# Access individuals in ancetral generations
pop.ancestor(5, 1).info('x')        # absolute index
pop.ancestor(0, 1, 1).info('x')     # relative index
# Or make ancestral generation the current generation and use 'individual'
pop.useAncestralGen(1)
pop.individual(5).info('x')         # absolute index
pop.individual(0, 1).info('x')      # relative index
# 'ancestor' can still access the 'present' (generation 0) generation
pop.ancestor(5, 0).info('x')
#end


#file log/batchAccess.log
import random
pop = population(size=[4, 6], loci=[2], infoFields='x')
pop.setIndInfo([random.randint(0, 10) for x in range(10)], 'x')
pop.indInfo('x')
pop.setGenotype([0, 1, 2, 3], 0)
pop.genotype(0)
pop.setVirtualSplitter(infoSplitter(cutoff=[3], field='x'))
pop.setGenotype([0])    # clear all values
pop.setGenotype([5, 6, 7], [1, 1])
pop.indInfo('x', 1)
pop.genotype(1)
#end


#file log/popInfo.log
pop = population(10)
pop.setInfoFields(['a', 'b'])
pop.addInfoFields('c')
pop.addInfoFields(['d', 'e'])
pop.infoFields()
#
cIdx = pop.infoIdx('c')
eIdx = pop.infoIdx('e')
pop.setIndInfo([1], cIdx)
for ind in pop.individuals():
    ind.setInfo(ind.info(cIdx) + 1, eIdx)

print pop.indInfo(eIdx)
#end          


#file log/ancestralPop.log
simu = simulator(population(500, loci=[1]), randomMating())
simu.evolve(
    preOps = initByFreq([0.5, 0.5]),
    ops = [
        # start recording ancestral generations at generation 18
        setAncestralDepth(2, at=[-2]),
        stat(alleleFreq=[0], begin=-3),
        pyEval(r"'%.3f\n' % alleleFreq[0][0]", begin=-3)
    ],
    gen = 20
)
pop = simu.population(0)
# start from current generation
for i in range(pop.ancestralGens(), -1, -1):
  pop.useAncestralGen(i)
  Stat(pop, alleleFreq=[0])
  print '%d   %.3f' % (i, pop.dvars().alleleFreq[0][0])

# restore to the current generation  
pop.useAncestralGen(0)  
#end

#file log/addRemoveLoci.log
pop = population(10, loci=[3], chromNames=['chr1'])
# 1 1 1, 
pop.setGenotype([1])
# 1 1 1, 0 0 0
pop.addChrom(lociPos=[0.5, 1, 2], lociNames=['rs1', 'rs2', 'rs3'],
    chromName='chr2')
pop1 = population(10, loci=[3], chromNames=['chr3'],
    lociNames=['rs4', 'rs5', 'rs6'])
# 2 2 2,
pop1.setGenotype([2])
# 1 1 1, 0 0 0, 2 2 2
pop.addChromFrom(pop1)
# 1 1 1, 0 0 0, 2 0 2 2 0
pop.addLoci(chrom=[2, 2], pos=[1.5, 3.5], names=['rs7', 'rs8'])
# 1 1 1, 0 0 0, 2 0 2 0
pop.removeLoci([8])
Dump(pop)
#end

#file log/extract.log
import random
pop = population(size=[10, 10], loci=[5, 5],
    infoFields=['x', 'y'])
InitByValue(pop, range(10))
pop.setIndInfo([-1]*4 + [0]*3 + [-1]*3 + [2]*4 + [-1]*3 + [1]*4, 'x')
pop1 = pop.extract(field='x', loci=[1, 2, 3, 6, 7], infoFields='x')
Dump(pop1, structure=False)
#end

#file log/popVars.log
from pprint import pprint
pop = population(100, loci=[2])
InitByFreq(pop, [0.3, 0.7])
print pop.vars()    # No variable now
pop.dvars().myVar = 21
print pop.vars()
Stat(pop, popSize=1, alleleFreq=[0])
# pprint prints in a less messy format
pprint(pop.vars())
# print number of allele 1 at locus 0
print pop.vars()['alleleNum'][0][1]
# use the dvars() function to access dictionary keys as attributes
print pop.dvars().alleleNum[0][1]
print pop.dvars().alleleFreq[0]
#end


#file log/expression.log
simu = simulator(population(100, loci=[1]),
    randomMating(), 5)
simu.evolve(
    preOps = initByFreq([0.5, 0.5]),
    ops = [
        stat(alleleFreq=[0]),
        terminateIf('alleleFreq[0][0] == 0. or alleleFreq[0][0] == 1.')
    ]
)
#end


#file log/savePop.log
pop = population(100, loci=[5], chromNames=['chrom1'])
pop.dvars().name = 'my population'
pop.save('sample.pop')
pop1 = LoadPopulation('sample.pop')
pop1.chromName(0)
pop1.dvars().name
#end

os.remove('sample.pop')

#file log/stageAndGen.log
simu = simulator(population(100, loci=[20]), randomMating())
simu.evolve(
    preOps = initByFreq([0.2, 0.8]),
    ops = [
        stat(alleleFreq=[0], begin=80, step=10),
        pyEval(r"'After gen %d: allele freq: %.2f\n' % (gen, alleleFreq[0][0])",
            begin=80, step=10),
        pyEval(r"'Around gen %d: allele Freq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1], stage=PrePostMating)
    ],
    postOps = [savePopulation(output='sample.pop')],
    gen=100
)
#end

os.remove('sample.pop')

#file log/dryrun.log
simu = simulator(population(100, loci=[20]), randomMating())
simu.evolve(
    preOps = initByFreq([0.2, 0.8]),
    ops = [
        stat(alleleFreq=[0], begin=80, step=10),
        pyEval(r"'After gen %d: allele freq: %.2f\n' % (gen, alleleFreq[0][0])",
            begin=80, step=10),
        pyEval(r"'Around gen %d: alleleFreq: %.2f\n' % (gen, alleleFreq[0][0])",
            at = [-10, -1], stage=PrePostMating)
    ],
    postOps = [savePopulation(output='sample.pop')],
    gen=100,
    dryrun = True
)
#end


#file log/replicate.log
simu = simulator(population(100, loci=[20]), randomMating(), 5)
simu.evolve(
    preOps = initByFreq([0.2, 0.8]),
    ops = [
        stat(alleleFreq=[0], step=10),
        pyEval('gen', step=10, rep=0),
        pyEval(r"'\t%.2f' % alleleFreq[0][0]", step=10, rep=(0, 2, -1)),
        pyOutput('\n', step=10, rep=-1)
    ],
    gen=30,
)
#end


#file log/output.log
simu = simulator(population(size=1000, loci=[2]), randomMating(), rep=3)
simu.evolve(
    preOps = initByValue([1, 2, 2, 1]),  
    ops = [
        recombinator(rates=0.01),
        stat(LD=[0, 1]),
        pyEval(r"'%.2f\t' % LD[0][1]", step=20, output='>>LD.txt'),
        pyOutput('\n', rep=-1, step=20, output='>>LD.txt'),
        pyEval(r"'%.2f\t' % R2[0][1]", output='R2.txt'),
        pyEval(r"'%.2f\t' % LD[0][1]", step=20, output="!'>>LD_%d.txt' % rep"),
    ],
    gen=100
)
print open('LD.txt').read()
print open('R2.txt').read()    # Only the last write operation succeed.
print open('LD_2.txt').read()  # Each replicate writes to a different file.
#end

#file log/outputFunc.log
import logging
# logging to a file simulation.log, with detailed debug information
logging.basicConfig(
    filename='simulation.log',
    level=logging.DEBUG,
    format='%(levelname)s: %(message)s',
    filemode='w'
)
# logging to standard output with less information
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(message)s')
console.setFormatter(formatter)
logger = logging.getLogger('')
logger.addHandler(console)
#
simu = simulator(population(size=1000, loci=[2]), randomMating())
simu.evolve(
    preOps = initByValue([1, 2, 2, 1]),  
    ops = [
        recombinator(rates=0.01),
        stat(LD=[0, 1]),
        pyEval(r"'LD: %d, %.2f' % (gen, LD[0][1])", step=20,
            output=logger.info),   # send LD to console and a logfile
        pyEval(r"'R2: %d, %.2f' % (gen, R2[0][1])", step=20,
            output=logger.debug),  # send R2 only to a logfile
    ],
    gen=100
)
print open('simulation.log').read()
#end

logging.shutdown()
for file in ['LD.txt', 'LD_0.txt', 'LD_1.txt', 'LD_2.txt', 'R2.txt', 'LD_2.txt', 'simulation.log']:
    os.remove(file)

#file log/transmitter.log
simu = simulator(population(size=10000, loci=[2]), randomMating())
simu.evolve(
    preOps = initByValue([1, 2, 2, 1]),
    ops = [
        # Recombination only happens after generation 30. A
        # mendelianGenoTransmitter defined in randomMating is responsible
        # for genotype transmission before that.
        recombinator(rates=0.01, begin=30),
        stat(LD=[0, 1]),
        pyEval(r"'gen %d, LD: %.2f\n' % (gen, LD[0][1])", step=20)
    ],
    gen=100
)
#end

#file log/hybrid.log
def myPenetrance(geno):
    'A three-locus heterogeneity penetrance model'
    if sum(geno) < 2:
        return 0
    else:
        return sum(geno)*0.1

simu = simulator(population(1000, loci=[20]*3), randomMating())
simu.evolve(
    preOps = initByFreq([0.8, 0.2]),
    ops = [
        pyPenetrance(func=myPenetrance, loci=[10, 30, 50]),
        stat(numOfAffected=True),
        pyEval(r"'%d: %d\n' % (gen, numOfAffected)")
    ],
    gen = 5
)
#end


#file log/pyOperator.log
def dynaMutator(pop, param):
    '''This mutator mutates commom loci with low mutation rate and rare
    loci with high mutation rate, as an attempt to raise allele frequency
    of rare loci to an higher level.'''
    # unpack parameter
    (cutoff, mu1, mu2) = param;
    Stat(pop, alleleFreq=range(pop.totNumLoci()))
    for i in range(pop.totNumLoci()):
        # Get the frequency of allele 1 (disease allele)
        if pop.dvars().alleleFreq[i][1] < cutoff:
            KamMutate(pop, k=2, rates=mu1, loci=[i])
        else:
            KamMutate(pop, k=2, rates=mu2, loci=[i])
    return True

simu = simulator(population(size=10000, loci=[2, 3]),
    randomMating())
simu.evolve(
    preOps = [ 
        initByFreq([.99, .01], loci=[0, 2, 4]),
        initByFreq([.8, .2], loci=[1, 3])],
    ops = [ 
        pyOperator(func=dynaMutator, param=(.2, 1e-2, 1e-5), stage=PreMating),
        stat(alleleFreq=range(5), step=10),
        pyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
            step=10),
    ],
    gen = 31
)                
#end

#file log/pyDuringMatingOperator.log
def rejectInd(off):
    'reject an individual if it off.allele(0) == 1'
    return off.allele(0) == 0

simu = simulator(population(size=100, loci=[1]),
    randomMating())
simu.evolve(
    preOps = initByFreq([0.5, 0.5]),
    ops = [ 
        pyOperator(func=rejectInd, stage=DuringMating, offspringOnly=True),
    ],
    gen = 1
)
# You should see no individual with allele 1 at locus 0, ploidy 0.
simu.population(0).genotype()[:20]
#end


#file log/newOperator.log
class dynaMutator(pyOperator):
    '''This mutator mutates commom loci with low mutation rate and rare
    loci with high mutation rate, as an attempt to raise allele frequency
    of rare loci to an higher level.'''
    def __init__(self, cutoff, mu1, mu2, *args, **kwargs):
        self.cutoff = cutoff
        self.mu1 = mu1
        self.mu2 = mu2
        pyOperator.__init__(self, func=self.mutate, *args, **kwargs)
    #
    def mutate(self, pop):
        Stat(pop, alleleFreq=range(pop.totNumLoci()))
        for i in range(pop.totNumLoci()):
            # Get the frequency of allele 1 (disease allele)
            if pop.dvars().alleleFreq[i][1] < self.cutoff:
                KamMutate(pop, k=2, rates=self.mu1, loci=[i])
            else:
                KamMutate(pop, k=2, rates=self.mu2, loci=[i])
        return True

simu = simulator(population(size=10000, loci=[2, 3]),
    randomMating())
simu.evolve(
    preOps = [ 
        initByFreq([.99, .01], loci=[0, 2, 4]),
        initByFreq([.8, .2], loci=[1, 3])],
    ops = [ 
        dynaMutator(cutoff=.2, mu1=1e-2, mu2=1e-5, stage=PreMating),
        stat(alleleFreq=range(5), step=10),
        pyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
            step=10),
    ],
    gen = 31
)          
#end


#file log/InitByFreq.log
def InitByFreq(pop, *args, **kwargs):
    initByFreq(*args, **kwargs).apply(pop)

InitByFreq(pop, [.2, .3, .5])
#end

#file log/migrSize.log
simu = simulator(
    population(size=[500, 1000], infoFields='migrate_to'),
    randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize')
    ],
    gen = 3
)
#end

#file log/migrFixedSize.log
simu = simulator(
    population(size=[500, 1000], infoFields='migrate_to'),
    randomMating(subPopSize=[500, 1000]))
simu.evolve(
    preOps = initSex(),
    ops = [
        migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
        stat(popSize=True, stage=PrePostMating),
        pyEval(r'"%s\n" % subPopSize', stage=PrePostMating)
    ],
    gen = 3
)
#end

#file log/demoFunc.log
def demo(gen, oldSize=[]):
    return [500 + gen*10, 1000 + gen*10]

simu = simulator(
    population(size=[500, 1000], infoFields='migrate_to'),
    randomMating(subPopSize=demo))
simu.evolve(
    preOps = initSex(),
    ops = [
        migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize')
    ],
    gen = 3
)
#end

#file log/numOff.log
def checkNumOffspring(ms):
    '''Check the number of offspring for each family using
       information field father_idx
    '''
    simu = simulator(
        population(size=[30], infoFields=['father_idx', 'mother_idx']),
        matingScheme=ms)
    simu.evolve(
        preOps = initSex(),
        ops=[parentsTagger()],
        gen=1)
    # get the parents of each offspring
    parents = [(x, y) for x, y in zip(simu.population(0).indInfo('mother_idx'),
        simu.population(0).indInfo('father_idx'))]
    # Individuals with identical parents are considered as siblings.
    famSize = []
    lastParent = (-1, -1)
    for parent in parents:
        if parent == lastParent:
            famSize[-1] += 1
        else:
            lastParent = parent
            famSize.append(1)
    return famSize

# Case 1: produce the given number of offspring
checkNumOffspring(randomMating(numOffspring=2))
# Case 2: Use a Python function
import random
def func(gen):
    return random.randint(5, 8)

checkNumOffspring(randomMating(numOffspring=func))
# Case 3: A geometric distribution
checkNumOffspring(randomMating(numOffspring=(GeometricDistribution, 0.3)))
# Case 4: A Possition distribution
checkNumOffspring(randomMating(numOffspring=(PoissonDistribution, 3)))
# Case 5: A Binomial distribution
checkNumOffspring(randomMating(numOffspring=(BinomialDistribution, 0.1, 10)))
# Case 6: A uniform distribution
checkNumOffspring(randomMating(numOffspring=(UniformDistribution, 2, 6)))
#end


#file log/sexMode.log
def checkSexMode(ms):
    '''Check the assignment of sex to offspring'''
    simu = simulator(
        population(size=[40]),
        matingScheme=ms)
    simu.evolve(preOps = initSex(), ops=[], gen=1)
    # return individual sex as a string
    return ''.join([ind.sexChar() for ind in simu.population(0).individuals()])

# Case 1: NoSex (all male, randomMating will not continue)
checkSexMode(randomMating(sexMode=NoSex))
# Case 2: RandomSex (Male/Female with probability 0.5)
checkSexMode(randomMating(sexMode=RandomSex))
# Case 3: ProbOfMale (Specify probability of male)
checkSexMode(randomMating(sexMode=(ProbOfMale, 0.8)))
# Case 4: NumOfMale (Specify number of male in each family)
checkSexMode(randomMating(numOffspring=3, sexMode=(NumOfMale, 1)))
# Case 5: NumOfFemale (Specify number of female in each family)
checkSexMode(randomMating(
    numOffspring=(UniformDistribution, 4, 6),
    sexMode=(NumOfFemale, 2))
)
#end

#file log/monogamous.log
simu = simulator(population(20, infoFields=['father_idx', 'mother_idx']),
    monogamousMating(numOffspring=2, sexMode=(NumOfMale, 1)))
simu.evolve(
    preOps = [initSex(sex=(Male, Female))],
    ops = [parentsTagger()],
    gen = 5
)
pop = simu.extract(0)
[ind.sex() for ind in pop.individuals()]
[ind.intInfo('father_idx') for ind in pop.individuals()]
[ind.intInfo('mother_idx') for ind in pop.individuals()]
# count the number of distinct parents
len(set(pop.indInfo('father_idx')))
len(set(pop.indInfo('mother_idx')))
#end

#file log/polygamous.log
simu = simulator(population(100, infoFields=['father_idx', 'mother_idx']),
    polygamousMating(polySex=Male, polyNum=2))
simu.evolve(
    preOps = initSex(),
    ops = [parentsTagger()],
    gen = 5
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals()][:20]
[ind.intInfo('mother_idx') for ind in pop.individuals()][:20]
#end

#file log/randomSelection.log
simu = simulator(population(100, ploidy=1, loci=[5, 5], ancGen=1,
    infoFields='parent_idx'),
    randomSelection())
simu.evolve(
    preOps = [initByFreq([0.3, 0.7])],
    ops = [parentTagger()],
    gen = 5
)
pop = simu.extract(0)
ind = pop.individual(0)
par = pop.ancestor(ind.intInfo('parent_idx'), 1)
print ind.sex(), ind.genotype()
print par.sex(), par.genotype()
#end

#file log/alphaMating.log
simu = simulator(population(1000, loci=[5], 
    infoFields=['father_idx', 'mother_idx', 'fitness']),
    alphaMating(alphaSex=Male, alphaNum=2))
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [parentsTagger(),
        maSelector(loci=[0], fitness=[0.8, 0.8, 1]),
        stat(alleleFreq=[0]),
        pyEval(r'"%.2f\n" % alleleFreq[0][1]', step=5)
    ],
    gen = 20,
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals()][:10]
[ind.intInfo('mother_idx') for ind in pop.individuals()][:10]
#end

#file log/haplodiploidMating.log
pop = population(10, ploidy=Haplodiploid, loci=[5, 5],
    infoFields=['father_idx', 'mother_idx'])
pop.setVirtualSplitter(sexSplitter())
simu = simulator(pop, haplodiploidMating())
simu.evolve(
    preOps = [initSex(),
        initByValue([0]*10, subPops=[(0, 0)], initSex=False),
        initByValue([1]*10+[2]*10, subPops=[(0, 1)], initSex=False)],
    ops = [parentsTagger(),
        dumper(structure=False, stage=PrePostMating)],
    gen = 1
)
#end

#file log/selfMating.log
pop = population(20, loci=[8])
# every chromosomes are different. :-)
for idx, ind in enumerate(pop.individuals()):
    ind.setGenotype([idx*2], 0)
    ind.setGenotype([idx*2+1], 1)

simu = simulator(pop, selfMating())
simu.evolve(
    ops = [recombinator(rates=0.1)],
    gen = 1
)
Dump(simu.population(0), width=3, structure=False, max=10)
#end

#file log/heteroMatingSP.log
pop = population(size=[1000, 1000], loci=[2],
    infoFields=['father_idx', 'mother_idx'])
simu = simulator(pop, heteroMating(
    [randomMating(numOffspring=2, subPop=0),
     randomMating(numOffspring=4, subPop=1)
    ])
)
simu.evolve(
    preOps = initSex(),
    ops= [parentsTagger()],
    gen=10
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals(0)][:10]
[ind.intInfo('father_idx') for ind in pop.individuals(1)][:10]
#end

#file log/heteroMatingVSP.log
pop = population(size=[1000], loci=[2],
    infoFields=['father_idx', 'mother_idx'])
pop.setVirtualSplitter(proportionSplitter([0.2, 0.8]))
simu = simulator(pop, heteroMating(
    matingSchemes = [
        selfMating(subPop=(0, 0)),
        randomMating(subPop=(0, 1))
    ])
)
simu.evolve(
    preOps = initSex(),
    ops= [parentsTagger()],
    gen = 10
)
pop = simu.extract(0)
[ind.intInfo('father_idx') for ind in pop.individuals(0)][:15]
[ind.intInfo('mother_idx') for ind in pop.individuals(0)][:15]
#end

#file log/heteroMatingWeight.log
pop = population(size=[1000], loci=[2],
    infoFields='mark')
pop.setVirtualSplitter(rangeSplitter([[0, 500], [200, 1000]]))
def markOff(param):
    '''define a Python during mating operator that marks
       individual information field 'mark'
    '''
    def func(off, param):
        off.setInfo(param, 'mark')
        return True
    return pyOperator(func=func, param=param, stage=DuringMating,
        offspringOnly=True)

simu = simulator(pop, heteroMating(
    matingSchemes = [
        randomMating(subPop=0, weight=-0.5, ops=[markOff(0)]),
        randomMating(subPop=(0, 0), weight=2, ops=[markOff(1)]),
        randomMating(subPop=(0, 1), weight=3, ops=[markOff(2)])
    ])
)
simu.evolve(
    preOps = initSex(),
    ops= [],
    gen = 10
)
marks = list(simu.extract(0).indInfo('mark'))
marks.count(0.)
marks.count(1.)
marks.count(2.)
#end


#file log/randomMating.log
def mendelianOffspringGenerator(ops=[], *args, **kwargs):
    'An offspring generator that uses mendelianGenoTransmitter()'
    return  offspringGenerator([mendelianGenoTransmitter()] + ops, *args, **kwargs)

def randomMating(numOffspring = 1., sexMode = RandomSex, ops = [], subPopSize = [],
        subPop = (), weight = 0, selectionField = 'fitness'):
    'A basic diploid sexual random mating scheme.'
    return homoMating(
        chooser = randomParentsChooser(True, selectionField),
        generator = mendelianOffspringGenerator(ops, numOffspring, sexMode),
        subPopSize = subPopSize,
        subPop = subPop,
        weight = weight)
#end


#file log/sequentialSelfing.log
simu = simulator(population(100, loci=[5]*3, infoFields='parent_idx'),
    homoMating(sequentialParentChooser(), selfingOffspringGenerator()))
simu.evolve(
    preOps = [initByFreq([0.2]*5)],
    ops = [
        parentTagger(),
        dumper(structure=False, stage=PrePostMating, max=5)],
    gen = 1
)
#end

#file log/controlledOffGenerator.log
def traj(gen):
    return [0.5 + gen * 0.01]

simu = simulator(population(1000, loci=[10]*2),
    homoMating(randomParentChooser(),
        controlledOffspringGenerator(loci=[5],
            alleles=[0], freqFunc=traj,
            ops = [selfingGenoTransmitter()]))
)

# evolve the population while keeping allele frequency 0.5
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [stat(alleleFreq=[5, 15]),
        pyEval(r'"%.2f\t%.2f\n" % (alleleFreq[5][0], alleleFreq[15][0])')],
    gen = 5
)

#end

#file log/mitochondrial.log
pop = population(10, loci=[5]*5,
    # one autosome, two sex chromosomes, and two mitochondrial chromosomes
    chromTypes=[Autosome, ChromosomeX, ChromosomeY] + [Customized]*2,
    infoFields=['father_idx', 'mother_idx'])

simu = simulator(pop, randomMating(ops=[mitochondrialGenoTransmitter()]))

simu.evolve(
    preOps=[initByFreq([0.4] + [0.2]*3)],
    ops=[
        recombinator(rates=0.1),
        parentsTagger(),
        dumper(structure=False),
    ],
    gen = 2
)
#end


#file log/sexSpecificRec.log
class sexSpecificRecombinator(pyOperator):
    def __init__(self, intensity=0, rates=0, loci=[], convMode=NoConversion,
            maleIntensity=0, maleRates=0, maleLoci=[], maleConvMode=NoConversion,
            *args, **kwargs):
        # This operator is used to recombine maternal chromosomes
        self.recombinator = recombinator(rates, intensity, loci, convMode)
        # This operator is used to recombine paternal chromosomes
        self.maleRecombinator = recombinator(maleRates, maleIntensity,
            maleLoci, maleConvMode)
        #
        self.initialized = False
        # Note the use of parameter isTransmitter
        pyOperator.__init__(self, func=self.transmitGenotype,
            stage=DuringMating, isTransmitter=True, *args, **kwargs)
    #
    def transmitGenotype(self, pop, off, dad, mom):
        # Recombinators need to be initialized. Basically, they cache some
        # population properties to speed up genotype transmission.
        if not self.initialized:
            self.recombinator.initialize(pop)
            self.maleRecombinator.initialize(pop)
            self.initialized = True
        # Form the first homologous copy of offspring.
        self.recombinator.transmitGenotype(mom, off, 0)
        self.maleRecombinator.transmitGenotype(dad, off, 1)
        return True

pop = population(10, loci=[15]*2, infoFields=['father_idx', 'mother_idx'])
simu = simulator(pop, randomMating())
simu.evolve(
    preOps=[initByFreq([0.4] + [0.2]*3)],
    ops=[
        parentsTagger(),
        sexSpecificRecombinator(rates=0.1, maleRates=0),
        dumper(structure=False),
    ],
    gen = 2
)
#end


#file log/infoChooser.log
pop = population(100, loci=[10],
    infoFields=['father_idx', 'mother_idx', 'sibling'])

pop.setVirtualSplitter(sexSplitter())

def locate_sibling(pop):
    '''The population is arranged as MFMFMFMF... where MF are siblings, so the
    sibling of males are 1, 3, 5, .. and the slibling of females are 0, 2, 4, ...
    '''
    pop.setIndInfo([2*x+1 for x in range(pop.popSize()/2)], 'sibling', (0, 0))
    pop.setIndInfo([2*x for x in range(pop.popSize()/2)], 'sibling', (0, 1))

simu = simulator(pop, consanguineousMating(func=locate_sibling, infoFields='sibling',
    numOffspring=2, sexMode=(NumOfMale, 1)))
simu.evolve(
    preOps = [initByFreq([0.2, 0.8], sex=[Male, Female])],
    ops = [
        parentsTagger(),
        dumper(structure=False, max=6, at=[-1])
    ],
    gen = 2
)
#end


#file log/generator.log
def func():
    i = 1
    all = 0
    while i <= 5:
        all += 1./i
        i += 1
        yield all 

for i in func():
    print '%.3f' % i,

#end

#file log/pyParentsChooser.log
from random import randint
def randomChooser(pop, sp):
    males = []
    females = []
    # identify males and females in each social rank
    for rank in range(3):
        males.append([x for x in pop.individuals(sp) \
            if x.sex() == Male and x.info('rank') == rank])
        females.append([x for x in pop.individuals(sp) \
            if x.sex() == Female and x.info('rank') == rank])
    while True:
        # choose a rank randomly
        rank = pop.individual(randint(0, pop.subPopSize(sp) - 1), sp).intInfo('rank')
        yield males[rank][randint(0, len(males[rank]) - 1)], \
            females[rank][randint(0, len(females[rank]) - 1)]

def setRank(pop, dad, mom, off):
    'The rank of offspring can increase or drop to zero randomly'
    off.setInfo((dad.info('rank') + randint(-1, 1)) % 3, 'rank')

pop = population(size=[1000, 2000], loci=[1], infoFields='rank')
pop.setIndInfo([randint(0, 2) for x in range(pop.popSize())], 'rank')

simu = simulator(pop, homoMating(
    pyParentsChooser(randomChooser),
    mendelianOffspringGenerator()))
simu.evolve(
    preOps = initSex(),
    ops = [],
    gen = 5
)    

#end

classFile = open('log/myParentsChooser.h', 'w')
classFile.write('''#include <stdlib.h>
#include <vector>
#include <utility>
using std::pair;
using std::vector;

class myParentsChooser
{
public:
	// A constructor takes all locations of male and female.
	myParentsChooser(const std::vector<int> & m, const std::vector<int> & f)
		: male_idx(m), female_idx(f)
	{
		srand(time(0));
	}
	
	pair<unsigned long, unsigned long> chooseParents()
	{
		unsigned long male = rand() % male_idx.size();
		unsigned long female = rand() % male_idx.size();
		return std::make_pair(male, female);
	}

private:
	vector<int> male_idx;
	vector<int> female_idx;
};
''')
classFile.close()

interFile = open('log/myParentsChooser.i', 'w')
interFile.write('''%module myParentsChooser
%{
#include "myParentsChooser.h"
%}

// std_vector.i for std::vector
%include "std_vector.i"
%template() std::vector<int>;

// stl.i for std::pair
%include "stl.i"
%template() std::pair<unsigned long, unsigned long>;

%include "myParentsChooser.h"
''')
interFile.close()

setupFile = open('log/setup.py', 'w')
setupFile.write('''from distutils.core import setup, Extension
import sys
# Under linux/gcc, lib stdc++ is needed for C++ based extension.
if sys.platform == 'linux2':
    libs = ['stdc++']
else:
    libs = []

setup(name = "myParentsChooser",
    description = "A sample parent chooser",
    py_modules = ['myParentsChooser'],  # will be generated by SWIG
    ext_modules = [
        Extension('_myParentsChooser',
            sources = ['myParentsChooser.i'],
            swig_opts = ['-O', '-templatereduce', '-shadow',
                '-python', '-c++', '-keyword', '-nodefaultctor'],
            include_dirs = ["."],
    )
  ]
)
''')
setupFile.close()

sys.path.append('log')

try:
    import myParentsChooser
except:
    os.chdir('log')
    os.system('python setup.py build_ext --swig-opts="-O -templatereduce -shadow -c++ -keyword -nodefaultctor" install --install-purelib="." --install-platlib="."')
    os.chdir('..')

#file log/cppParentChooser.log
# The class myParentsChooser is defined in module myParentsChooser
from myParentsChooser import myParentsChooser

def parentsChooser(pop, sp):
    'How to call a C++ level parents chooser.'
    # create an object with needed information (such as x, y) ...
    pc = myParentsChooser(
        [x for x in range(pop.popSize()) if pop.individual(x).sex() == Male],
        [x for x in range(pop.popSize()) if pop.individual(x).sex() == Female])
    while True:
        # return indexes of parents repeatedly
        yield pc.chooseParents()

pop = population(100, loci=[1])
simu = simulator(pop,
    homoMating(pyParentsChooser(parentsChooser), mendelianOffspringGenerator())
)
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [],
    gen = 100
)

#end


#file log/simuGen.log
simu = simulator(population(50, loci=[10], ploidy=1),
    randomSelection(), rep=3)
simu.evolve(ops = [], gen = 5)
simu.gen()
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [
        stat(alleleFreq=[5]),
        ifElse('alleleNum[5][0] == 0',
            pyEval(r"'Allele 0 is lost in rep %d at gen %d\n' % (rep, gen)")),
        ifElse('alleleNum[5][0] == 50',
            pyEval(r"'Allele 0 is fixed in rep %d at gen %d\n' % (rep, gen)")),
        terminateIf('alleleNum[5][0] == 0 or alleleNum[5][0] == 50'),
    ],
)
simu.gen()
#end


#file log/twoStage.log
# First stage: use the standard random mating scheme, do not use any
# information field for efficiency considerations.
simu = simulator(population(500, loci=[10]), randomMating())
simu.evolve(preOps = [initByFreq([0.5, 0.5])],
    ops = [], gen = 50)
# Second stage: track parents and produce more offspring per mating
# event. In preparation for pedigree ascertainment.
for pop in simu.populations():
    pop.addInfoFields(['father_idx', 'mother_idx'])
    pop.setAncestralDepth(1)

simu.setMatingScheme(randomMating(numOffspring=2))
simu.evolve(
    ops = [
        parentsTagger(),
        maPenetrance(loci=0, penetrance=(0.2, 0.4, 0.5))
    ],
    gen = 5
)
# Sample affected sibpairs
pop = simu.extract(0)
sample = AffectedSibpairSample(pop, size=5)[0]
[ind.intInfo('father_idx') for ind in sample.individuals()]
#end

#file log/changeStru.log
import random
def mutator(pop, param):
    'Parameter has a length of region and a mutation rate at each basepair'
    region, rate = param
    # there are certainly more efficient algorithm, but this 
    # example just mutate each basepair one by one....
    for i in range(region):
        if random.random() < rate:
            try:
                idx = pop.addLoci(chrom=0, pos=i)[0]
            except:
                # position might duplicate
                continue
            # choose someone to mutate
            ind = pop.individual(random.randint(0, pop.popSize() - 1))
            ind.setAllele(1, idx)
    return True

# The populations start with no loci at all.
simu = simulator(population(1000, loci=[]), randomMating(), rep=3)
simu.evolve(
    preOps = initSex(),
    ops = [pyOperator(func=mutator, param=(10000, 2e-6))],
    gen = 200
)
for pop in simu.populations():
    print pop.totNumLoci(), pop.lociPos()

#end


#file log/simuFunc.log
simu = simulator(population(100, loci=[5, 10], infoFields='x'),
    randomMating(), rep=5)
simu.evolve(preOps=[initByFreq([0.4, 0.6])],
    ops=[], gen=10)
# clone
cloned = simu.clone()
# save and load, using a different mating scheme
simu.save("sample.sim")
loaded = LoadSimulator("sample.sim", randomMating(numOffspring=2))
# 
simu.numRep()
loaded.numRep()
for pop1,pop2 in zip(cloned.populations(), loaded.populations()):
    assert pop1 == pop2

# continue to evolve
simu.evolve(ops=[], gen=10)
simu.gen()
#end

os.remove('sample.sim')


#file log/initSex.log
pop = population(size=[1000, 1000])
InitSex(pop, maleFreq=0.3, subPops=0)
InitSex(pop, sex=[Male, Female, Female], subPops=1)
Stat(pop, numOfMale=True)
print pop.dvars(0).numOfMale
print pop.dvars(1).numOfMale
#end

#file log/initByFreq.log
pop = population(size=[2, 3], loci=[5, 7])
InitByFreq(pop, alleleFreq=[[.2, .8], [.8, .2]])
Dump(pop, structure=False)
#end

#file log/initByFreqIdenticalInds.log
pop = population(size=[2, 3], loci=[5, 7])
InitByFreq(pop, alleleFreq=[.2, .8], identicalInds=True)
Dump(pop, structure=False)
#end

#file log/initByValue.log
pop = population(size=[2, 3], loci=[5, 7])
InitByValue(pop, [1]*5 + [2]*7 + [3]*5 +[4]*7)
Dump(pop, structure=False)
#end

#file log/initByValueProp.log
pop = population(size=[6, 8], loci=[5, 7])
pop.setVirtualSplitter(sexSplitter())
# initialize sex and the first two loci
InitByValue(pop, loci=range(5), value=range(10))
# initialize all males
InitByValue(pop, loci=range(5, 12), value=[2]*7,
    subPops=[(0, 0), (1, 0)], initSex=False)
# initialize females by proportion
InitByValue(pop, loci=range(5, 12), ploidy=1, value=[[3]*7, [4]*7],
    initSex=False, subPops=[(0, 1), (1, 1)], proportions=[0.4, 0.6])
Dump(pop, structure=False)
#end

#file log/dumper.log
pop = population(size=[10, 10], loci=[20, 30], infoFields='gen',
    ancGen=-1)
pop.setVirtualSplitter(sexSplitter())
pop1 = pop.clone()
InitByFreq(pop, [0]*20 + [0.1]*10)
pop.setIndInfo(1, 'gen')
InitByFreq(pop1, [0]*50 + [0.1]*10)
pop1.setIndInfo(2, 'gen')
pop.push(pop1)
Dump(pop, width=3, loci=[5, 6, 30], subPops=([0, 0], [1, 1]),
    max=10, structure=False, ancGen=-1)

#end

#file log/savePopulation.log
simu = simulator(population(100, loci=[2]),
    randomMating(), rep=5)
simu.evolve(
    preOps = [initByFreq([0.2, 0.8])],
    ops = [
        savePopulation(output="!'snapshot_%d_%d.pop' % (rep, gen)",
            step = 10),
        ],
    gen = 50
)
#end

for rep in range(5):
    for gen in range(0, 50, 10):
        os.remove('snapshot_%d_%d.pop' % (rep, gen))

#file log/setAncDepth.log
simu = simulator(population(100, infoFields=['father_idx', 'mother_idx']),
    randomMating(), rep=5)
simu.evolve(
    preOps = [initByFreq([0.3, 0.7])],
    ops = [
        setAncestralDepth(2, at=-2),
        parentsTagger(begin=-2)
    ],
    gen = 100
)
pop = simu.population(3)
print pop.ancestralGens()
print pop.ancestor(10, 1).info('father_idx')
#end


#file log/ifElse.log
simu = simulator(
    population(size=1000, loci=[1]),
    randomMating(), rep=4)
simu.evolve(
    preOps = [
        initByFreq([0.5, 0.5]),
        pyExec('below40, above60 = 0, 0')
    ],
    ops = [
        stat(alleleFreq=[0]),
        ifElse('alleleFreq[0][1] < 0.4',
            pyExec('below40 += 1')),
        ifElse('alleleFreq[0][1] > 0.6',
            pyExec('above60 += 1')),
        ifElse('alleleFreq[0][1] == 0 or alleleFreq[0][1] == 1',
            pyExec('stoppedAt = gen')),
        terminateIf('alleleFreq[0][1] == 0 or alleleFreq[0][1] == 1')
    ]
)
for pop in simu.populations():
    print 'Overall: %4d, below 40%%: %4d, above 60%%: %4d' % \
        (pop.dvars().stoppedAt, pop.dvars().below40, pop.dvars().above60)

#end

#file log/debug.log
simu = simulator(population(100, loci=[1]), randomMating(), rep=5)
simu.evolve(
    preOps = [initByFreq([0.1, 0.9])],
    ops = [
        stat(alleleFreq=[0]),
        ifElse('alleleNum[0][0] == 0',
            turnOnDebug(DBG_MUTATOR),
            turnOffDebug(DBG_MUTATOR)),
        ifElse('alleleNum[0][0] == 0',
            pointMutator(loci=0, allele=0, inds=0)),
    ],
    gen = 100
)
#end

#file log/pause.log
simu = simulator(population(100), randomMating(), rep=10)
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [pause(stopOnKeyStroke=str(x), rep=x) for x in range(10)],
    gen = 100
)
#end

#file log/ticToc.log
simu = simulator(population(10000, loci=[100]*5), randomMating(), rep=2)
simu.evolve(
    preOps = [initByFreq([0.1, 0.9])],
    ops = [
        stat(alleleFreq=[0]),
        ticToc(step=50, rep=-1),
    ],
    gen = 101
)
#end

#file log/pyExec.log
simu = simulator(population(100, loci=[1]),
    randomMating(), rep=2)
simu.evolve(
    preOps = [
        initByFreq([0.2, 0.8]),
        pyExec('traj=[]')
    ],
    ops = [
        stat(alleleFreq=[0]),
        pyExec('traj.append(alleleFreq[0][1])'),
    ],
    gen=5
)
# print trajectory
print ', '.join(['%.3f' % x for x in simu.dvars(0).traj])
#end

#file log/pyEval.log
simu = simulator(population(1000, loci=[1],
    infoFields=['mother_idx', 'father_idx']),
    randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        stat(alleleFreq=[0]),
        parentsTagger(),
        pyEval(r'"gen %d, #father %d, #mother %d\n"' \
            ' % (gen, numFather, numMother)',
            stmts="numFather = len(set(pop.indInfo('father_idx')))\n"
                "numMother = len(set(pop.indInfo('mother_idx')))",
            exposePop='pop')
    ],
    gen=3
)
#end

#file log/infoEval.log
import random
pop = population(20, loci=[1], infoFields='a')
pop.setVirtualSplitter(infoSplitter('a', cutoff=[3]))
InitByFreq(pop, [0.2, 0.8])
pop.setIndInfo([random.randint(2, 5) for x in range(20)], 'a')
InfoEval(pop, 'a', subPops=[(0, 0)]);print
InfoEval(pop, 'ind.allele(0, 0)', exposeInd='ind');print
# use population variables
pop.dvars().b = 5
InfoEval(pop, '"%d " % (a+b)', usePopVars=True);print
#end

#file log/infoExec.log
pop = population(100, loci=[1], infoFields=['a', 'b', 'c'])
InitByFreq(pop, [0.2, 0.8])
InfoExec(pop, 'a=1')
print pop.indInfo('a')[:10]
InfoExec(pop, 'b=ind.sex()', exposeInd='ind')
print pop.indInfo('b')[:10]
InfoExec(pop, 'c=a+b')
print pop.indInfo('c')[:10]
pop.dvars().d = 5
InfoExec(pop, 'c+=d', usePopVars=True)
print pop.indInfo('c')[:10]
#end

#file log/migrateByProb.log
simu = simulator(
    population(size=[1000]*3, infoFields='migrate_to'),
    randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        migrator(rate=[
            [0, 0.1, 0.1],
            [0, 0, 0.1],
            [0, 0.1, 0]
        ]),
        stat(popSize=True),
        pyEval('subPopSize'),
        pyOutput('\n')
    ],
    gen = 5
)        
#end

#file log/migrateByPropAndCount.log
simu = simulator(
    population(size=[1000]*3, infoFields='migrate_to'),
    randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        migrator(rate=[[0.1], [0.2]],
            mode=ByProportion,
            subPops=[1, 2],
            toSubPops=[3]),
        stat(popSize=True),
        pyEval('subPopSize'),
        pyOutput('\n')
    ],
    gen = 5
)        
#
simu.evolve(
    ops = [
        migrator(rate=[[50, 50], [100, 50]],
            mode=ByCounts,
            subPops=[3, 2],
            toSubPops=[2, 1]),
        stat(popSize=True),
        pyEval('subPopSize'),
        pyOutput('\n')
    ],
    gen = 5
)        
#end

#file log/migrateVSP.log
pop = population(size=[1000]*2, infoFields='migrate_to')
pop.setVirtualSplitter(sexSplitter())
simu = simulator(pop, randomMating())
simu.evolve(
    # 500 males and 500 females
    preOps = [initSex(sex=[Male, Female])],
    ops = [
        migrator(rate=[
            [0, 0.10],
            [0, 0.05],
            ],
            mode = ByProportion,
            subPops=[(0, 0), (0, 1)]),
        stat(popSize=True, numOfMale=True, stage=PrePostMating),
        pyEval(r"'%d/%d\t%d/%d\n' % (subPop[0]['numOfMale'], subPopSize[0], "
            "subPop[1]['numOfMale'], subPopSize[1])", stage=PrePostMating),
    ],
    gen = 2
)   
#end

#file log/manualMigration.log
pop = population([10]*2, infoFields='migrate_to')
pop.setIndInfo([0, 1, 2, 3]*5, 'migrate_to')
Migrate(pop, mode=ByIndInfo)
pop.subPopSizes()
#end

#file log/splitBySize.log
simu = simulator(population(1000), randomSelection())
simu.evolve(
    ops = [
        splitSubPops(subPops=0, sizes=[300, 300, 400], at=2),
        stat(popSize=True),
        pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 4
)
#end

#file log/splitByProp.log
def demo(gen, oldSize=[]):
    if gen < 2:
        return 1000 + 100 * gen
    else:
        return [x + 50 * gen for x in oldSize]

simu = simulator(population(1000),
    randomSelection(subPopSize=demo))
simu.evolve(
    ops = [
        splitSubPops(subPops=0, proportions=[.5]*2, at=2),
        stat(popSize=True),
        pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 4
)
#end


#file log/splitByInfo.log
import random
pop = population([1000]*3, subPopNames=['a', 'b', 'c'], infoFields='x')
pop.setIndInfo([random.randint(0, 3) for x in range(1000)], 'x')
print pop.subPopSizes()
print pop.subPopNames()
SplitSubPops(pop, subPops=[0, 2], infoFields=['x'])
print pop.subPopSizes()
print pop.subPopNames()
#end

#file log/mergeSubPops.log
simu = simulator(population([500]*2),
    randomSelection())
simu.evolve(
    ops = [
        mergeSubPops(subPops=[0, 1], at=3),
        stat(popSize=True),
        pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 5
)
#end


#file log/resizeSubPops.log
simu = simulator(population([500]*2),
    randomSelection())
simu.evolve(
    ops = [
        resizeSubPops(proportions=(1.5, 2), at=3),
        stat(popSize=True),
        pyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
    ],
    gen = 5
)

#end

#file log/recRate.log
simu = simulator(population(size=[1000], loci=[100]),
    randomMating(), rep=2)
simu.evolve(
    preOps = [initByValue([0]*100 + [1]*100)],
    ops = [
        recombinator(rates=0.01, rep=0),
        recombinator(rates=[0.01]*10, loci=range(50, 60), rep=1),
        stat(LD=[[40, 55], [60, 70]]),
        pyEval(r'"%d:\t%.3f\t%.3f\t" % (rep, LD_prime[40][55], LD_prime[60][70])'),
        pyOutput('\n', rep=-1)
    ],
    gen = 5
)
#end

#file log/recIntensity.log
simu = simulator(population(size=[1000], loci=[3], lociPos=[0, 1, 1.1]),
    randomMating())
simu.evolve(
    preOps = [initByValue([0]*3 + [1]*3)],
    ops = [
        recombinator(intensity=0.01),
        stat(LD=[[0, 1], [1, 2]]),
        pyEval(r'"%.3f\t%.3f\n" % (LD_prime[0][1], LD_prime[1][2])', step=10)
    ],
    gen = 50
)
#end

#file log/conversion.log
simu = simulator(population(size=[1000], loci=[100]),
    randomMating(), rep=2)
simu.evolve(
    preOps = [initByValue([0]*100 + [1]*100)],
    ops = [
        recombinator(rates=0.01, loci=50, rep=0),
        recombinator(rates=0.01, loci=50, rep=1,
            convMode=(NumMarkers, 1, 10)),
        stat(LD=[[40, 55], [40, 70]]),
        pyEval(r'"%d:\t%.3f\t%.3f\t" % (rep, LD_prime[40][55], LD_prime[40][70])'),
        pyOutput('\n', rep=-1)
    ],
    gen = 5
)
#end

#file log/matrixMutator.log
simu = simulator(population(size=[2000], loci=[1]),
    randomMating())
simu.evolve(
    preOps = [initByFreq([0.2, 0.3, 0.5])],
    ops = [
        matrixMutator(rate = [
            [0, 1e-5, 1e-5],
            [1e-4, 0, 1e-4],
            [1e-3, 1e-3, 0]
        ]),
        stat(alleleFreq=[0], step=100),
        pyEval(r"', '.join(['%.3f' % alleleFreq[0][x] for x in range(3)]) + '\n'",
            step=100),
    ],
    gen=1000
)
#end


#file log/kamMutator.log
pop = population(size=[2000], loci=[1]*3)
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        kamMutator(k=5, rates=[1e-2, 1e-3], loci=[0, 1]),
        stat(alleleFreq=range(3), step=100),
        pyEval(r"', '.join(['%.3f' % alleleFreq[x][0] for x in range(3)]) + '\n'",
            step=100),
    ],
    gen=500
)
#end

#file log/snpMutator.log
pop = population(size=[2000], loci=[1, 1], infoFields='fitness')
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        snpMutator(u=0.001),
        maSelector(loci=0, fitness=[1, 0.99, 0.98]),
        stat(alleleFreq=[0, 1], step=100),
        pyEval(r"'%.3f\t%.3f\n' % (alleleFreq[0][1], alleleFreq[1][1])",
            step=100),
    ],
    gen=500
)
#end

#file log/acgtMutator.log
pop = population(size=[2000], loci=[1],
    alleleNames=['A', 'C', 'G', 'T'])
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = [initByFreq([.1, .1, .1, .7])],
    ops = [
        acgtMutator(rate=[1e-4, 0.5], model='K80'),
        stat(alleleFreq=[0], step=100),
        pyEval(r"', '.join(['%.3f' % alleleFreq[0][x] for x in range(4)]) + '\n'",
            step=100),
    ],
    gen=500
)
#end

#file log/smmMutator.log
simu = simulator(population(size=1000, loci=[1, 1]), randomMating())
simu.evolve(
    # all start from allele 50
    preOps = [initByFreq( [0]*50 + [1])],
    ops = [
        smmMutator(rates=1e-3, loci=0),
        smmMutator(rates=1e-3, incProb=0.6, loci=1,
            mutStep=(GeometricDistribution, 0.2)),
    ],
    gen=100
)
# count the average number tandem repeats at both loci
cnt0 = cnt1 = 0
for ind in simu.population(0).individuals():
    cnt0 += ind.allele(0, 0) + ind.allele(0, 1)
    cnt1 += ind.allele(1, 0) + ind.allele(1, 1)

print 'Average number of repeats at two loci are %.2f and %.2f.' % \
    (cnt0/2000., cnt1/2000.)

#end

#file log/pyMutator.log
import random
def incAllele(allele):
    return allele + random.randint(1, 5)

simu = simulator(population(size=1000, loci=[20]),
    randomMating())
simu.evolve(
    preOps = initSex(),
    ops = [
        pyMutator(func=incAllele, rates=[1e-4, 1e-3],
            loci=[2, 10])
    ],
    gen = 1000
)
# count the average number tandem repeats at both loci
def avgAllele(pop, loc):
    ret = 0
    for ind in pop.individuals():
        ret += ind.allele(loc, 0) + ind.allele(loc, 1)
    return ret / (pop.popSize() * 2.)

pop = simu.population(0)
print 'Average number of repeats at two loci are %.2f and %.2f.' % \
    (avgAllele(pop, 2), avgAllele(pop, 10))
#end

#file log/mixedMutator.log
simu = simulator(population(5000, loci=[1, 1]),
    randomMating())
simu.evolve(
    preOps = initByValue([50, 50]),
    ops = [
        # the first locus uses a pure stepwise mutation model
        smmMutator(rates=0.001, loci=0),
        # the second locus uses a mixed model
        mixedMutator(rates=0.001, loci=1, mutators=[        
            kamMutator(rates=1, k=100),
            smmMutator(rates=1)
        ], prob=[0.1, 0.9])],
    gen = 20
)
# what alleles are there?
geno0 = []
geno1 = []
for ind in simu.population(0).individuals():
    geno0.extend([ind.allele(0, 0), ind.allele(0, 1)])
    geno1.extend([ind.allele(1, 0), ind.allele(1, 1)])

print 'Locus 0 has alleles', ', '.join([str(x) for x in set(geno0)])
print 'Locus 1 has alleles', ', '.join([str(x) for x in set(geno1)])
#end

#file log/contextMutator.log
simu = simulator(population(5000, loci=[3, 3]),
    randomMating())
simu.evolve(
    # initialize locus by 0, 0, 0, 1, 0, 1
    preOps = initByValue([1, 1], loci=[3, 5]),
    ops = [
        contextMutator(mutators=[
            snpMutator(u=0.1),
            snpMutator(u=1),
            ],
            contexts=[(0, 0), (1, 1)],
            loci=[1, 4],
            rates=0.01
        ),
        stat(alleleFreq=[1, 4], step=5),
        pyEval(r"'Gen: %2d freq1: %.3f, freq2: %.3f\n'" + 
            " % (gen, alleleFreq[1][1], alleleFreq[4][1])", step=5)
    ], 
    gen = 20
)
#end

#file log/pyContextMutator.log
import random

simu = simulator(population(5000, loci=[3, 3]),
    randomMating())

def contextMut(allele, context):
    if context == [0, 0]:
        if allele == 0 and random.random() < 0.1:
            return 1
    elif context == [1, 1]:
        if allele == 0:
            return 1
    # do not mutate
    return allele

simu.evolve(
    # initialize locus by 0, 0, 0, 1, 0, 1
    preOps = initByValue([1, 1], loci=[3, 5]),
    ops = [
        pyMutator(func=contextMut, context=1,
            loci=[1, 4],  rates=0.01
        ),
        #snpMutator(u=0.01, v= 0.01, loci=[1, 4]),
        stat(alleleFreq=[1, 4], step=5),
        pyEval(r"'Gen: %2d freq1: %.3f, freq2: %.3f\n'" + 
            " % (gen, alleleFreq[1][1], alleleFreq[4][1])", step=5)
    ], 
    gen = 20
)
#end

#file log/pointMutator.log
pop = population(1000, loci=[1], infoFields='fitness')
simu = simulator(pop, randomSelection())
simu.evolve(
    preOps = pyOutput('Introducing alleles at generation'),
    ops = [
        stat(alleleFreq=[0]),
        ifElse('alleleNum[0][1] == 0', ifOps=[
            pyEval(r"' %d' % gen"),
            pointMutator(inds=0, loci=0, allele=1),
        ]),
        maSelector(loci=0, wildtype=0, fitness=[1, 1.05, 1.1]),
        ifElse('alleleFreq[0][1] > 0.05', ifOps=[
            pyEval(r"'.\nTerminate at generation %d at allele freq %.3f.\n'" +
                " % (gen, alleleFreq[0][1])"),
            terminateIf('True'),
        ])
    ],
)
#end

#file log/mutatorVSP.log
def fragileX(geno):
    '''A disease model where an individual has increased risk of 
    affected if the number of tandem repeats exceed 75.
    '''
    # Alleles A1, A2.
    maxRep = max(geno)
    if maxRep < 50:
        return 0
    else:
        # individuals with allele >= 70 will surely be affected
        return min(1, (maxRep - 50)*0.05)

def avgAllele(pop):
    'count number of alleles by affection status.'
    Stat(pop, numOfAffected=True)
    aff = pop.dvars().numOfAffected
    unaff = pop.dvars().numOfUnaffected
    #
    # This can be simplied after stator supports VSP.
    cnt = [0, 0, 0]
    for ind in pop.individuals():
        reps  = ind.allele(0, 0) + ind.allele(0, 1)
        if ind.affected():
            cnt[1] += reps
        else: # X, X
            cnt[0] += reps
        cnt[2] += ind.allele(1, 0) + ind.allele(1, 1)
    if unaff != 0:
        cnt[0] /= 2. * unaff
    if aff != 0:
        cnt[1] /= 2. * aff
    cnt[2] /= 2. * (aff + unaff)
    # male, female, loc2
    pop.dvars().avgAllele = cnt
    return True

pop = population(10000, loci=[1, 1])
pop.setVirtualSplitter(affectionSplitter())
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = initByValue([50, 50]),
    ops = [
        # determine affection status for each offspring (duringMating)
        pyPenetrance(func=fragileX, loci=0),
        # unaffected offspring, mutation rate is high to save some time
        smmMutator(rates=1e-3, loci=1),
        # unaffected offspring, mutation rate is high to save some time
        smmMutator(rates=1e-3, loci=0, subPops=[(0, 0)]),
        # affected offspring have high probability of mutating upward
        smmMutator(rates=1e-2, loci=0, subPops=[(0, 1)],
           incProb=0.7, mutStep=3),
        # number of affected
        pyOperator(func=avgAllele, step=20),
        pyEval(r"'Gen: %3d #Aff: %d AvgRepeat: %.2f (unaff), %.2f (aff), %.2f (unrelated)\n'"
            + " % (gen, numOfAffected, avgAllele[0], avgAllele[1], avgAllele[2])",
            step=20),
    ],
    gen = 101
)

#end

#file log/alleleMapping.log
pop = population(size=[2000], loci=[1])
simu = simulator(pop, randomMating())
simu.evolve(
    preOps = initByFreq([0]*4 + [0.1, 0.2, 0.3, 0.4]),
    ops = [
        kamMutator(k=4, rates=1e-4, mapIn=[0]*4 + range(4),
            mapOut=[4, 5, 6, 7]),
        stat(alleleFreq=[0], step=100),
        pyEval(r"', '.join(['%.2f' % alleleFreq[0][x] for x in range(8)]) + '\n'",
            step=100),
    ],
    gen=500
)
#end

#file log/rpy.log
from simuPOP import *
from simuRPy import varPlotter
pop = population(size=1000, loci=[2])
simu = simulator(pop, randomMating(), rep=3)
simu.evolve(
    preOps = initByValue([1, 2, 2, 1]),
    ops = [
        recombinator(rates=0.01),
        stat(LD=[0, 1]),
        varPlotter('LD[0][1]', step=5, update=40, saveAs='log/rpy.png',
            legend=['Replicate %d' % x for x in range(3)],
            ylab='LD between marker 1 and 2',
            ylim=[0, 0.25], main='LD decay', lty_rep=[1, 2, 3],
        ),
    ],
    gen=100
)
#end

#file log/rpyByRep.log
from simuPOP import *
from simuRPy import varPlotter
pop = population(size=1000, loci=[1]*4)
simu = simulator(pop, randomMating(), rep=3)
simu.evolve(
    preOps = [initByFreq([0.1*(x+1), 1-0.1*(x+1)], loci=x) for x in range(4)],
    ops = [
        stat(alleleFreq=range(4)),
        varPlotter('[alleleFreq[x][0] for x in range(4)]', byRep=True,
            update=10, saveAs='log/rpy_byRep.png',
            legend=['Locus %d' % x for x in range(4)],
            ylab='Allele frequency',
            ylim=[0, 1],
            main_rep=['Genetic drift, replicate %d' % x for x in range(3)],
        ),
    ],
    gen=100
)
#end

#file log/rpyByDim.log
from simuPOP import *
from simuRPy import varPlotter
pop = population(size=1000, loci=[1]*4)
simu = simulator(pop, randomMating(), rep=3)

def drawFrame(r, dim=None, **kwargs):
    '''Draw a frame around subplot dim. Parameter r is defined in the rpy
    module and is used for calling R functions. Parameter dim is the dimension
    index. Other parameters are ignored.
    '''
    r.axis(1)
    r.axis(2)
    r.grid()
    r.mtext({0:'A', 1:'B', 2:'C', 3:'D'}[dim], adj=1)

simu.evolve(
    preOps = [initByFreq([0.1*(x+1), 1-0.1*(x+1)], loci=x) for x in range(4)],
    ops = [
        stat(alleleFreq=range(4)),
        varPlotter('[alleleFreq[x][0] for x in range(4)]', byDim=True,
            update=10, saveAs='log/rpy_byDim.png',
            legend=['Replicate %d' % x for x in range(3)],
            ylab='Allele frequency',
            ylim=[0, 1],
            main_dim=['Genetic drift, freq=%.1f' % ((x+1)*0.10) for x in range(4)],
            col_rep=['red', 'blue', 'black'],
            lty_rep=[1, 2, 3],
            # the default png dimension is 800x600
            dev_print_width=600, dev_print_height=500,
            # do not draw axes in r.plot, leaving the job to drawFrame
            plot_axes=False,
            # plot frame, grid etc after each r.plot call
            plotHook = drawFrame,
        ),
    ],
    gen=100
)
#end

#file log/scatterPlotter.log
from simuPOP import *
from simuRPy import scatterPlotter
import random
pop = population([500], infoFields=['x', 'y', 'anc'])
# random sex
InitSex(pop)
# random geographic location
pop.setIndInfo([random.random() for i in range(500)], 'x')
pop.setIndInfo([random.random() for i in range(500)], 'y')
# anc is 0 or 1
pop.setIndInfo([random.randint(0, 1) for i in range(500)], 'anc')
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(infoSplitter('anc', cutoff=[0.2, 0.4, 0.6, 0.8]))
#
def passInfo(fields):
    'Parental fields will be passed as x1, y1, anc1, x2, y2, anc2'
    x1, y1, anc1, x2, y2, anc2 = fields
    anc = (anc1 + anc2)/2.
    x = (x1 + x2)/2 + random.normalvariate(anc - 0.5, 0.1)
    y = (y1 + y2)/2 + random.normalvariate(0, 0.1)
    return x, y, anc

simu = simulator(pop, randomMating())
simu.evolve(
    ops = [
        pyTagger(passInfo, infoFields=['x', 'y', 'anc']),
        scatterPlotter(['x', 'y'], stage=PreMating,            
            saveAs = 'log/scatterPlotter.png',
            subPops = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4)],
            ylim = [0, 1.2],
            main = "!'Ancestry distribution of individuals at generation %d' % gen",
            legend = ['anc < 0.2', '0.2 <= anc < 0.4', '0.4 <= anc < 0.6',
                '0.6 <= anc < 0.8', '0.8 <= anc'],
            plot_axes = False,
            par_mar = [0, 0, 2, 0],
        ),
    ],
    gen = 5,
)
#end

#file log/histPlotter.log
from simuPOP import *
from simuRPy import histPlotter, qqPlotter, boxPlotter
import random
pop = population([500], infoFields=['x', 'y', 'anc'])
# random sex
InitSex(pop)
# random geographic location
pop.setIndInfo([random.random() for i in range(500)], 'x')
pop.setIndInfo([random.random() for i in range(500)], 'y')
# anc is 0 or 1
pop.setIndInfo([random.randint(0, 1) for i in range(500)], 'anc')
# Defines VSP 0, 1, 2, 3, 4 by anc.
pop.setVirtualSplitter(sexSplitter())
#
def passInfo(fields):
    'Parental fields will be passed as x1, y1, anc1, x2, y2, anc2'
    x1, y1, anc1, x2, y2, anc2 = fields
    anc = (anc1 + anc2)/2.
    x = (x1 + x2)/2 + random.normalvariate(anc - 0.5, 0.1)
    y = (y1 + y2)/2 + random.normalvariate(0, 0.1)
    return x, y, anc

simu = simulator(pop, randomMating())
simu.evolve(
    ops = [
        pyTagger(passInfo, infoFields=['x', 'y', 'anc']),
        histPlotter(infoFields='anc', stage=PreMating,            
            subPops=[(0,0), (0,1)], col_sp=['blue', 'red'],
            saveAs='log/histPlotter.png',
            main="!'Histogram of ancestry values at generation %d' % gen",
        ),
        qqPlotter(infoFields='anc', stage=PreMating,            
            subPops=[(0,0), (0,1)], col_sp=['blue', 'red'],
            saveAs='log/qqPlotter.png',
            main="!'QQ plot of ancestry values at generation %d' % gen",
        ),
        boxPlotter(infoFields='anc', stage=PreMating,
            subPops=[(0,0), (0,1)],
            saveAs='whatever',
            dev_print_file='!"log/Gen%d.png" % gen',
            main="!'Boxplots of ancestry values at generation %d' % gen",
        ),
    ],
    gen = 5,
)
#end

#file log/getParam.log
import types, simuOpt

options = [
    {'arg': 'r:',
     'longarg': 'rate=',
     'default': [0.01],
     'useDefault': True,
     'label': 'Recombination rate',
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Recombination rate for each replicate. If a single value
            is given, it will be used for all replicates.''',
     'validate': simuOpt.valueListOf(simuOpt.valueBetween(0, 1))
    },
    {'longarg': 'rep=',
     'default': 5,
     'label': 'Number of replicates',
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'Number of replicates to simulate.',
     'validate': simuOpt.valueGT(0)
    }, 
    {'longarg': 'pop=',
     'default': 'CEU',
     'label': 'Initial population',
     'allowedTypes': [types.StringType],
     'description': '''Use one of the HapMap populations as the initial
            population for this simulation. You can choose from:
            |YRI: 33 trios from the Yoruba people in Nigeria (Africa)
            |CEU: 30 trios from Utah with European ancestry (European)
            |CHB+JPT: 90 unrelated individuals from China and Japan (Asia)
            ''',
     'chooseOneOf': ['CEU', 'YRI', 'CHB+JPT'],
     'validate': simuOpt.valueOneOf(['CEU', 'YRI', 'CHB+JPT'])
    }
]

pars = simuOpt.simuOpt(options, 'A demo simulation')
print pars.usage()
# You can manually feed parameters...
pars.processArgs(['--rep=10'])
pars.rep
#beginignore
oldArg = sys.argv[-1]
sys.argv.pop()
if not os.path.isfile('getParam.png'):
    print 'Run a GUI if getParam has not been runned'
else:
    sys.argv.extend(['--rate=[0.25]', '--rep=5', '--pop="CEU"'])
    simuOpt.setOptions(gui=False)
    pars.processArgs(sys.argv)

#endignore
# but simuOpt.getParam is the easiest to used
if not pars.getParam():
    sys.exit(1)

#beginignore
sys.argv[1] = oldArg
#endignore
# save prameters to a configuration file
pars.saveConfig('sample.cfg')
# post-process parameters
pars.rate
pars.rep
pars.rate = pars.rate * pars.rep
# extract parameters as a dictionary or a list
pars.asDict()
pars.asList()
# Default value of parameter rep is changed
# additional attribute is added.
par1 = simuOpt.simuOpt(options, # all parameters with default values
    rep=50,                     # default value of rep is changed
    additional=10               # derived parameters are added
)
# print all parameters except for derived ones.
print par1.asDict()
# All parameters are derived ...
par2 = simuOpt.simuOpt(rep=50, pop='CEU', rate=[0.5])
print par2.asDict()
print par2.rep, par2.pop
#end


#file log/reichDemo.log
import math
def demo_model(model, N0=1000, N1=100000, G0=500, G1=500):
    '''Return a demographic function 
    model: linear or exponential
    N0:   Initial population size.
    N1:   Ending population size.
    G0:   Length of burn-in stage.
    G1:   Length of population expansion stage.
    '''
    def ins_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:
            return N1
    rate = (math.log(N1) - math.log(N0))/G1
    def exp_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:            
            return int(N0 * math.exp((gen - G0) * rate))
    if model == 'instant':
        return ins_expansion
    elif model == 'exponential':
        return exp_expansion

# when needed, create a demographic function as follows
demo_func = demo_model('exponential', 1000, 100000, 500, 500)
# population size at generation 700
print demo_func(700)
#end


#file log/reichStat.log
class ne(pyOperator):
    '''Define an operator that calculates effective number of
    alleles at given loci. The result is saved in a population
    variable ne.
    '''
    def __init__(self, loci, *args, **kwargs):
        self.loci = loci
        pyOperator.__init__(self, func=self.calcNe, *args, **kwargs)
    #
    def calcNe(self, pop):
        Stat(pop, alleleFreq=self.loci)
        ne = {}
        for loc in self.loci:
            freq = pop.dvars().alleleFreq[loc][1:]
            sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
            if sumFreq == 0:
                ne[loc] = 0
            else:
                ne[loc] = 1. / sum([(x/sumFreq)**2 for x in freq])
        # save the result to the population.
        pop.dvars().ne = ne
        return True

def Ne(pop, loci):
    '''Function form of operator ne'''
    ne(loci).apply(pop)
    return pop.dvars().ne

pop = population(100, loci=[10])
InitByFreq(pop, [.2] * 5)
print Ne(pop, loci=[2, 4])
#end

#file log/reichEvolve.log
#beginignore
import math
def demo_model(model, N0=1000, N1=100000, G0=500, G1=500):
    '''Return a demographic function 
    model: linear or exponential
    N0:   Initial population size.
    N1:   Ending population size.
    G0:   Length of burn-in stage.
    G1:   Length of population expansion stage.
    '''
    def ins_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:
            return N1
    rate = (math.log(N1) - math.log(N0))/G1
    def exp_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:            
            return int(N0 * math.exp((gen - G0) * rate))
    if model == 'instant':
        return ins_expansion
    elif model == 'exponential':
        return exp_expansion

class ne(pyOperator):
    '''Define an operator that calculates effective number of
    alleles at given loci. The result is saved in a population
    variable ne.
    '''
    def __init__(self, loci, *args, **kwargs):
        self.loci = loci
        pyOperator.__init__(self, func=self.calcNe, *args, **kwargs)
    #
    def calcNe(self, pop):
        Stat(pop, alleleFreq=self.loci)
        ne = {}
        for loc in self.loci:
            freq = pop.dvars().alleleFreq[loc][1:]
            sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
            if sumFreq == 0:
                ne[loc] = 0
            else:
                ne[loc] = 1. / sum([(x/sumFreq)**2 for x in freq])
        # save the result to the population.
        pop.dvars().ne = ne
        return True

#endignore

def simulate(model, N0, N1, G0, G1, spec, s, mu, k):
    '''Evolve a population using given demographic model
    and observe the evolution of its allelic spectrum.
    model: type of demographic model.
    N0, N1, G0, G1: parameters of demographic model.
    spec: initial allelic spectrum, should be a list of allele
        frequencies for each allele.
    s: selection pressure.
    mu: mutation rate.
    k: k for the k-allele model
    '''
    demo_func = demo_model(model, N0, N1, G0, G1)
    simu = simulator(
        population(size=demo_func(0), loci=[1], infoFields='fitness'),
        randomMating(subPopSize=demo_func)
    )
    simu.evolve(
        preOps = initByFreq(loci=[0], alleleFreq=spec),
        ops=[
            kamMutator(k=k, rates=mu),
            maSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
            ne(loci=[0], step=100),
            pyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
                step=100),
        ],
        gen = G0 + G1
    )

simulate('instant', 1000, 10000, 500, 500, [0.9]+[0.02]*5, 0.01, 1e-4, 200)
#end


#file log/simuCDCV.log
#!/usr/bin/env python
#
# Author:  Bo Peng
# Purpose: A real world example for simuPOP user's guide.
#
'''
Simulation the evolution of allelic spectra (number and frequencies
of alleles at a locus), under the influence of population expansion,
mutation, and natural selection.
'''

import simuOpt
simuOpt.setOptions(quiet=True, alleleType='long')
from simuPOP import *

import sys, types, os, math

options = [
    {'longarg': 'demo=',
     'default': 'instant',
     'label': 'Population growth model',
     'description': 'How does a population grow from N0 to N1.',
     'chooseOneOf': ['instant', 'exponential'],
    },
    {'longarg': 'N0=',
     'default': 10000,
     'label': 'Initial population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': '''Initial population size. This size will be maintained
                till the end of burnin stage''',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'N1=',
     'default': 100000,
     'label': 'Final population size',
     'allowedTypes': [types.IntType, types.LongType],
     'description': 'Ending population size (after population expansion)',
     'validate': simuOpt.valueGT(0)
    }, 
    {'longarg': 'G0=',
     'default': 500,
     'label': 'Length of burn-in stage',
     'allowedTypes': [types.IntType],
     'description': 'Number of generations of the burn in stage.',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'G1=',
     'default': 1000,
     'label': 'Length of expansion stage',
     'allowedTypes': [types.IntType],
     'description': 'Number of geneartions of the population expansion stage',
     'validate': simuOpt.valueGT(0)
    },
    {'longarg': 'spec=',
     'default': [0.9] + [0.02]*5,
     'label': 'Initial allelic spectrum',
     'allowedTypes': [types.TupleType, types.ListType],
     'description': '''Initial allelic spectrum, should be a list of allele
            frequencies, for allele 0, 1, 2, ... respectively.''',
     'validate': simuOpt.valueListOf(simuOpt.valueBetween(0, 1)),
    },
    {'longarg': 's=',
     'default': 0.01,
     'label': 'Selection pressure',
     'allowedTypes': [types.IntType, types.FloatType],
     'description': '''Selection coefficient for homozygtes (aa) genotype.
            A recessive selection model is used so the fitness values of
            genotypes AA, Aa and aa are 1, 1 and 1-s respectively.''',
     'validate': simuOpt.valueGT(-1),
    },
    {'longarg': 'mu=',
     'default': 1e-4,
     'label': 'Mutation rate',
     'allowedTypes': [types.IntType, types.FloatType],
     'description': 'Mutation rate of a k-allele mutation model',
     'validate': simuOpt.valueBetween(0, 1),
    },
    {'longarg': 'k=',
     'default': 200,
     'label': 'Maximum allelic state',
     'allowedTypes': [types.IntType],
     'description': 'Maximum allelic state for a k-allele mutation model',
     'validate': simuOpt.valueGT(1),
    },
]

def demo_model(type, N0=1000, N1=100000, G0=500, G1=500):
    '''Return a demographic function 
    type: linear or exponential
    N0:   Initial population size.
    N1:   Ending population size.
    G0:   Length of burn-in stage.
    G1:   Length of population expansion stage.
    '''
    def ins_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:
            return N1
    rate = (math.log(N1) - math.log(N0))/G1
    def exp_expansion(gen, oldsize=[]):
        if gen < G0:
            return N0
        else:            
            return int(N0 * math.exp((gen - G0) * rate))
    if type == 'instant':
        return ins_expansion
    elif type == 'exponential':
        return exp_expansion

class ne(pyOperator):
    '''Define an operator that calculates effective number of
    alleles at given loci. The result is saved in a population
    variable ne.
    '''
    def __init__(self, loci, *args, **kwargs):
        self.loci = loci
        pyOperator.__init__(self, func=self.calcNe, *args, **kwargs)
    #
    def calcNe(self, pop):
        Stat(pop, alleleFreq=self.loci)
        ne = {}
        for loc in self.loci:
            freq = pop.dvars().alleleFreq[loc][1:]
            sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
            if sumFreq == 0:
                ne[loc] = 0
            else:
                ne[loc] = 1. / sum([(x/sumFreq)**2 for x in freq])
        # save the result to the population.
        pop.dvars().ne = ne
        return True

def simuCDCV(model, N0, N1, G0, G1, spec, s, mu, k):
    '''Evolve a population using given demographic model
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
    print demo_func(0)
    simu = simulator(
        population(size=demo_func(0), loci=[1], infoFields='fitness'),
        randomMating(subPopSize=demo_func)
    )
    simu.evolve(
        preOps=initByFreq(loci=[0], alleleFreq=spec),
        ops=[
            kamMutator(rate=mu, maxAllele=k),
            maSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
            ne(loci=[0], step=100),
            pyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
                step=100),
        ],
        gen = G0 + G1
    )
    return simu.extract(0)

if __name__ == '__main__':
    # get parameters
    par = simuOpt.simuOpt(options, __doc__)
    if not par.getParam():
        sys.exit(1)
    #
    if not sum(par.spec) == 1:
        print 'Initial allelic spectrum should add up to 1.'
        sys.exit(1)
    # save user input to a configuration file
    par.saveConfig('simuCDCV.cfg')
    #
    simuCDCV(*par.asList())

#end

out = os.popen('python log/simuCDCV.py -h')
hlp = open('log/simuCDCV.hlp', 'w')
print >> hlp, out.read()
hlp.close()

################################################

 

#### #file log/recombinator.log
simu = simulator(population(4, loci=[4, 5, 6], 
    infoFields=['father_idx', 'mother_idx']),
    randomMating())
simu.evolve(
    preOps = [initByFreq([.2, .2, .4, .2]), dumper(structure=False) ],
    ops = [parentsTagger()],
    postOps = [ dumper(structure=False)],
    gen=1
)
simu.evolve(
    ops = [
        parentsTagger(),
        recombinator(rate=[1, 1, 1], loci=[2, 6, 10])
    ],
    postOps = [dumper(structure=False)],
    gen=1
)
#end




#file log/mapSelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[1], infoFields='fitness'),
    randomMating())
s1 = .1
s2 = .2
simu.evolve(
    preOps=initByFreq(alleleFreq=[.2, .8]),
    ops = [
        stat(alleleFreq=[0], genoFreq=[0]),
        mapSelector(loci=0, fitness={'0-0':(1-s1), '0-1':1, '1-1':(1-s2)}),
        pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
    ],
    gen=300
)
#end

 
#file log/maSelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[1], infoFields='fitness'),
    randomMating())
s1 = .1
s2 = .2
simu.evolve(
    preOps=initByFreq(alleleFreq=[.2, .8]),
    ops = [
        stat(alleleFreq=[0], genoFreq=[0]),
        maSelector(loci=0, fitness=[1-s1, 1, 1-s2]),
        pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
    ],
    gen = 300)
#end

 
#file log/mlSelector.log
simu = simulator(
    population(size=10, ploidy=2, loci=[2], 
    infoFields=['fitness', 'spare']),
    randomMating())
simu.evolve(
    [ mlSelector([
         mapSelector(loci=0, fitness={'0-0':1, '0-1':1, '1-1':.8}),
         mapSelector(loci=1, fitness={'0-0':1, '0-1':1, '1-1':.8}),
         ], mode = Additive),
    ],
    preOps = initByFreq(alleleFreq=[.2, .8]),
    gen = 2
)
#end



#file log/pySelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[3], infoFields='fitness'),
    randomMating()
)

s1 = .2
s2 = .3
# the second parameter gen can be used for varying selection pressure
def sel(arr, gen=0):
  if arr[0] == 1 and arr[1] == 1:
    return 1 - s1
  elif arr[0] == 1 and arr[1] == 2:
    return 1
  elif arr[0] == 2 and arr[1] == 1:
    return 1
  else:
    return 1 - s2

# test func
print sel([1, 1])

simu.evolve(
    preOps=initByFreq(alleleFreq=[.2, .8]),
    ops = [
        stat(alleleFreq=[0], genoFreq=[0]),
        pySelector(loci=[0, 1], func=sel),
        pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=25)
    ],
    gen=100
)
#end


#file log/mapPenetrance.log
pop = population(size=[2,8], ploidy=2, loci=[2] )
InitByFreq(pop, [.2, .8])
MapPenetrance(pop, loci=0, 
    penetrance={'0-0':0, '0-1':1, '1-1':1})
Stat(pop, numOfAffected=1)
#end


#file log/mlPenetrance.log
pop = population(1000, loci=[3])
InitByFreq(pop, [0.3, 0.7])
pen = []
for loc in (0, 1, 2):
    pen.append(maPenetrance(loci=loc, wildtype=[1],
        penetrance=[0, 0.3, 0.6] ) )

# the multi-loci penetrance
MlPenetrance(pop, mode=Multiplicative, peneOps=pen)
Stat(pop, numOfAffected=True)
print pop.dvars().numOfAffected
#end

#file log/pyPenetrance.log
pop = population(1000, loci=[3])
InitByFreq(pop, [0.3, 0.7])
def peneFunc(geno):
    p = 1
    for l in range(len(geno)/2):
        p *= (geno[l*2]+geno[l*2+1])*0.3
    return p

PyPenetrance(pop, func=peneFunc, loci=(0, 1, 2))
Stat(pop, numOfAffected=True)
print pop.dvars().numOfAffected
#
# You can also define a function, that returns a penetrance
# function using given parameters
def peneFunc(table):
    def func(geno):
      return table[geno[0]][geno[1]]
    return func

# then, given a table, you can do
PyPenetrance(pop, loci=(0, 1, 2),
    func=peneFunc( ((0, 0.5), (0.3, 0.8)) ) )
#end

