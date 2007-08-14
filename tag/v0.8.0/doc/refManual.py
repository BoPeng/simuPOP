#
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

#file log/ref_importSimuPOPOpt.log
import simuOpt
simuOpt.setOptions(optimized=False, alleleType='long', quiet=True)
from simuPOP import *
print alleleType()
print optimized()
#end

#file log/ref_help.log
help(population.addInfoField)
#end

#file log/ref_standard.log
pop = population(10, loci=[2])
pop.locusPos(10)
pop.individual(20).setAllele(1, 0)
#end

#file log/ref_absIndex.log
pop = population(subPop=[20, 30], loci=[5, 6])
print pop.chromLocusPair(7)
print pop.absLocusIndex(1,1)
print pop.absIndIndex(10, 1)
print pop.subPopIndPair(40)
#end


#file log/ref_popAndOperator.log
simu = simulator(pop, randomMating(), rep=3)
simu.evolve(
    preOps = [ initByFreq([.8, .2])],
    ops = [
        stat(alleleFreq=[0,1], Fst=[1], step=10),
        kamMutator(rate=0.001, rep=1),
        kamMutator(rate=0.0001, rep=2)
    ],
    end=10
)
#end


#file log/ref_InitByFreq.log
def InitByFreq(pop, *args, **kwargs):
    initByFreq(*args, **kwargs).apply(pop)

InitByFreq(pop, [.2, .3, .4, .1])
#end


#file log/ref_carray.log
# obtain an object using one of the arrXXX functions
pop = population(loci=[3,4], lociPos=[1,2,3,4,5,6,7])
arr = pop.arrLociPos()
# print and expression (just like list)
print arr
str(arr)
# count
arr.count(2)
# index 
arr.index(2)
# can read write
arr[0] = 0.5
# the underlying locus position is also changed
print pop.lociPos()
# convert to list
arr.tolist()
# or simply
list(arr)
# compare to list directly
arr == [0.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
# you can also convert and compare
list(arr) == [0.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
# slice
arr[:] = [1,2,3,4,5,6,7]
# assign from another part
arr[1:3] = arr[3:5]
# arr1 is 1,2,3
arr1 = arr[:3]
# assign slice from a number
# arr will also be affected since arr1 point to a part of arr
arr1[:] = 10
# assign vector of the same length
len(arr1)
arr1[:] = [30,40, 50]
#end


#file log/ref_dumpPop.log
pop = population(size=5, ploidy=2, loci=[5,10],
        lociPos=[range(0,5),range(0,20,2)],
        alleleNames=['A','C','T','G'],
        subPop=[2,3], maxAllele=3)
# .apply form
initByFreq([.2, .3, .4, .1]).apply(pop)
# function form
Dump(pop)
#end


#file log/ref_popVars.log
from simuUtil import ListVars
pop = population(subPop=[1000, 2000], loci=[1])
InitByFreq(pop, [0.2, 0.8])
ListVars(pop.vars(), useWxPython=False)
Stat(pop, popSize=1, alleleFreq=[0])
# subPop is True by default, use name to limit the variables to display
ListVars(pop.vars(), useWxPython=False, subPop=False, name='alleleFreq')
# print number of allele 1 at locus 0
print pop.vars()['alleleNum'][0][1]
print pop.dvars().alleleNum[0][1]
print pop.dvars().alleleFreq[0]
print pop.dvars(1).alleleNum[0][1]
#end


#file log/ref_localNamespace.log
pop = population(subPop=[1000, 2000], loci=[1])
InitByFreq(pop, [0.2, 0.8])
Stat(pop, popSize=1, alleleFreq=[0])
print pop.evaluate('alleleNum[0][0] + alleleNum[0][1]')
pop.execute('newPopSize=int(popSize*1.5)')
ListVars(pop.vars(), level=1, useWxPython=False)
# this variable is 'local' to the population and is
# not available in the main namespace
newPopSize
#
simu = simulator(population(10),noMating(), rep=2)
# evaluate an expression in different areas
print simu.vars(1)
print simu.population(0).evaluate("grp*2")
print simu.population(1).evaluate("grp*2")
# a statement (no return value)
simu.population(0).execute("myRep=2+rep*rep")
simu.population(1).execute("myRep=2*rep")
print simu.vars(0)
#end


#file log/ref_pyEval.log
simu = simulator(population(100, loci=[1]),
    randomMating(), rep=2)
simu.evolve(
    preOps = [initByFreq([0.2, 0.8])],
    ops = [ stat(alleleFreq=[0]),
        pyExec('myNum = alleleNum[0][0] * 2'),
        pyEval(r'"gen %d, rep %d, num %d, myNum %d\n"' \
            ' % (gen, rep, alleleNum[0][0], myNum)')
        ],
    end=2
)
#end

#file log/ref_info1.log
pop = population(10, infoFields=['a', 'b'])
aIdx = pop.infoIdx('a')
bIdx = pop.infoIdx('b')
for ind in pop.individuals():
    a = ind.info(aIdx)
    ind.setInfo(a+1, bIdx)

print pop.indInfo(bIdx, False)
#end


#file log/ref_info2.log
pop = population(10, infoFields=['a', 'b'])
aIdx = pop.infoIdx('a')
bIdx = pop.infoIdx('b')
info = pop.arrIndInfo(False)
sz = pop.infoSize()
for idx in range(pop.popSize()):
    info[sz*idx + bIdx] = info[sz*idx + aIdx] + 1

print pop.indInfo(bIdx, False)
#end

#file log/ref_ancestralPop.log
simu = simulator(population(10000, loci=[2]), randomMating())
simu.evolve(
    ops = [
        setAncestralDepth(5, at=[-5]),
        kamMutator(rate=0.01, loci=[0], maxAllele=1),
        kamMutator(rate=0.001, loci=[1], maxAllele=1)
    ],
    end = 20
)
pop = simu.population(0)
# start from current generation
for i in range(pop.ancestralDepth()+1):
    pop.useAncestralPop(i)
    Stat(pop, alleleFreq=[0,1])
    print '%d     %5f     %5f' % \
        (i, pop.dvars().alleleFreq[0][1], pop.dvars().alleleFreq[1][1])

# restore to the current generation    
pop.useAncestralPop(0)    
#end



#turnOnDebug(DBG_SIMULATOR)
#turnOnDebug(DBG_UTILITY)


# Note that I can not use pop now, since it is
# obtained from simu.population(0) which is invalid now.


#file log/ref_saveFstat.log
from simuUtil import *
SaveFstat(pop, "pop.dat", maxAllele=9)
print open("pop.dat").read()
pop2 = LoadFstat("pop.dat")
#end

import os
os.remove('pop.xml')
os.remove('pop.bin')
os.remove('pop.dat')
os.remove('pop.txt')


#file log/ref_randomSample.log
# random sample
# [0]: RandomSample already return
#    a list of samples even if times=1 (default)
Dump( RandomSample(pop, 3)[0])
#end


#file log/ref_calcStat.log
Stat(pop, popSize=1, alleleFreq=range(0, pop.totNumLoci()),
    heteroFreq=range(0,pop.totNumLoci()), Fst=[0])
#end

#file log/ref_simulatorsaveload.log
simu.saveSimulator("s.txt")
simu.saveSimulator("s.xml", format="xml")
simu.saveSimulator("s.bin", format="bin")
simu1 = LoadSimulator("s.txt", randomMating())
simu2 = LoadSimulator("s.xml", randomMating(), format="xml")
simu3 = LoadSimulator("s.bin", randomMating(), format="bin")
#end

# remove these files
os.remove('s.txt')
os.remove('s.xml')
os.remove('s.bin')
os.remove('a0.txt')
os.remove('a1.txt')

#file log/ref_pyOperator.log
def dynaMutator(pop, param):
    ''' this mutator mutate common loci with low mutation rate
    and rare loci with high mutation rate, as an attempt to
    bring allele frequency of these loci at an equal level.'''
    # unpack parameter
    (cutoff, mu1, mu2) = param;
    Stat(pop, alleleFreq=range( pop.totNumLoci() ) )
    for i in range( pop.totNumLoci() ):
        # 1-freq of wild type = total disease allele frequency
        if 1-pop.dvars().alleleFreq[i][1] < cutoff:
            KamMutate(pop, maxAllele=2, rate=mu1, loci=[i])
        else:
            KamMutate(pop, maxAllele=2, rate=mu2, loci=[i])
    return True
#end

#file log/ref_pyOperatorUse.log
pop = population(size=10000, ploidy=2, loci=[2, 3])

simu = simulator(pop, randomMating())

simu.evolve(
    preOps = [ 
        initByFreq( [.6, .4], loci=[0,2,4]),
        initByFreq( [.8, .2], loci=[1,3]) ],
    ops = [ 
        pyOperator( func=dynaMutator, param=(.5, .1, 0) ),
        stat(alleleFreq=range(5)),
        pyEval(r'"%f\t%f\n"%(alleleFreq[0][1],alleleFreq[1][1])', step=10)
        ],
    end = 30
)                
#end



#file log/ref_pySubset.log
simu = simulator(population(subPop=[2,3], loci=[3,4], infoFields=['fitness']),
        randomMating())
simu.step([
        initByFreq([.3,.5,.2]),
        pySubset( [1,-1,-1,1,-1] ),
        dumper(alleleOnly=True, stage=PrePostMating)
     ])
#end



#turnOnDebug(DBG_ALL)
#turnOnDebug(DBG_SIMULATOR)

#file log/ref_varPlotter.log
from simuUtil import *
from simuRPy import *

simu = simulator( population(size=200, ploidy=2, loci=[3,4],
    subPop=[50,50,100]), randomMating(), rep=4)

# migrate
migr = migrator([[0,.2,.1],[.25,0,.1],[.1,.2,0]],
    mode=MigrByProbability)
# and count the size of subpopulations
stat = stat(popSize=1, stage=PreMating)
# plot subPopSize. 
simu.evolve([
     migr, 
     stat,
     varPlotter('subPopSize', numRep=4, byRep=1, 
         varDim=3, win=10, title='subPop size', saveAs='log/ref_simuDemo')
     ],
     end=30
)

#end
#PS /usr/bin/convert log/ref_simuDemo16.eps log/simuDemo16.png
#PS /bin/rm -f log/ref_simuDemo*.eps

#file log/ref_ifElse.log
from simuRPy import *
from simuUtil import *
numRep=4
popSize=100
endGen=50

simu = simulator(population(size=popSize, loci=[1]),
    randomMating(), rep=numRep)
simu.evolve(
    preOps = [ initByValue([1,1])],
    ops = [
        # penetrance, additve penetrance
        maPenetrance(locus=0, wildtype=[1], penetrance=[0,0.5,1]),
        # count number of affected
        stat(numOfAffected=True),
        # introduce disease if no one is affected
        ifElse(cond='numOfAffected==0',
            ifOp=kamMutator(rate=0.01, maxAllele=2)),
        # expose affected status
        pyExec('pop.exposeAffectedness()', exposePop=True),
        # plot affected status
        varPlotter(expr='affected',plotType="image", byRep=1, update=endGen, 
            varDim=popSize, win=endGen, numRep=numRep,
            title='affected status', saveAs="ifElse")
    ],
    end=endGen,
    dryrun=False
)
#end

#file log/ref_rng.log
print ListAllRNG()
# get system RNG
print rng().name()
# set system RNG
SetRNG("taus2", seed=10)
print rng().name()
# create a seprate RNG instance
r=RNG()
for n in range(1,5):
    print r.randBinomial(10, .7),

#end




#file log/ref_extgenostru.log
pop = population(1, loci=[2,3,4])
print pop.numLoci(1)
print pop.locusPos(2)
dis = pop.arrLociPos()
print dis
dis[2] = 0.5
print pop.locusPos(2)
print pop.arrLociPos()
#end

#file log/ref_extgenotype.log
InitByFreq(pop, [.2,.8])
Dump(pop, alleleOnly=1)
ind = pop.individual(0)
print ind.allele(1,1)
ind.setAllele(3,1,1)
Dump(pop, alleleOnly=1)
a = ind.arrGenotype()
print a
a = ind.arrGenotype(1)
print a
a = ind.arrGenotype(1,2)
print a
a[2]=4
# the allele on the third chromosome has been changed
Dump(pop, alleleOnly=1)
#end

#file log/ref_extother.log
print ind.sex()
print ind.sexChar()
ind.setSex(Female)
ind.setAffected(True)
print ind.tag()
ind.setTag([1,2])
Dump(pop)
#end

#file log/ref_extsimu.log
simu = simulator(pop, randomMating(), rep=3)
pop1 = simu.population(1)
ind1 = pop1.individual(0)
ind1.setAllele(3,0)
Dump(pop1)
#end

#file log/ref_extoperator.log
simu.step([ pyEval(stmts="pop=simu.population(rep)")])
#end


#file log/ref_tab.log
def tab(**kwargs):
    parm = ''    
    for (k,v) in kwargs.items():
        parm += ' , ' + str(k) + '=' + str(v)
    cmd = r'output( """\t""" ' + parm + ')'
    # print cmd
    return eval(cmd)
#end

#file log/ref_saveFstat.tmp
def saveInFstatFormat(pop, output='', outputExpr='', maxAllele=0):
    if output != '':
        file = output
    elif outputExpr != '':
        file = eval(outputExpr, globals(), pop.vars() )
    else:
        raise exceptions.ValueError, "Please specify output or outputExpr"
    # open file
    try:
        f = open(file, "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + " to write."
    #    
    # file is opened.
    np = pop.numSubPop()
    if np > 200:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 200 samples"
    nl = pop.totNumLoci()
    if nl > 100:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 100 loci"
    if maxAllele != 0:
        nu = maxAllele
    else:
        nu = pop.maxAllele()
    if nu > 999:
        print "Warning: Current version (2.93) of FSTAT can not handle more than 999 alleles at each locus"
        print "If you used simuPOP_la library, you can specify maxAllele in population constructure"
    if nu < 10:
        nd = 1
    elif nu < 100:
        nd = 2
    elif nu < 1000:
        nd = 3
    else: # FSTAT can not handle this now. how many digits?
        nd = len(str(nu))
    # write the first line
    f.write( '%d %d %d %d\n' % (np, nl, nu, nd) )
    # following lines with loci name.
    for ch in range(0, pop.numChrom()):
        for al in range(0, pop.numLoci(ch)):
            f.write( "loc_%d_%d\n" % (ch, al))
    # genoSize=totNumLoci()*ploidy()
    gs = pop.totNumLoci()
    for sp in range(0, pop.numSubPop()):
        # genotype of subpopulation sp, individuals are
        # rearranged in perfect order
        gt = pop.arrGenotype(sp)
        for ind in range(0, pop.subPopSize(sp)):
            f.write("%d " % (sp+1))
            p1 = 2*gs*ind                # begining of first hemo copy
            p2 = 2*gs*ind + gs     # second
            for al in range(0, gs): # allele
                ale1 = gt[p1+al]
                ale2 = gt[p2+al]
                if ale1 == 0 or ale2 == 0:
                    f.write('%%%dd' % (2*nd) % 0 )
                else:
                    f.write('%%0%dd%%0%dd ' % (nd, nd) % (ale1, ale2))
            f.write( "\n")
    f.close()
#end
#PS head -15 log/ref_saveFstat.tmp > log/saveInFstatFormat.log
#PS rm -f log/ref_saveFstat.tmp

#file log/ref_saveFstat.log
def saveFstat(output='', outputExpr='', **kwargs):
    # deal with additional arguments
    parm = ''
    for (k,v) in kwargs.items():
        parm += str(k) + '=' + str(v) + ', '
    # pyEval( exposePop=1, param?, stmts="""
    # saveInFSTATFormat( pop, rep=rep?, output=output?, outputExpr=outputExpr?)
    # """)
    opt = '''pyEval(exposePop=1, %s
        stmts=r\'\'\'saveInFstatFormat(pop, rep=rep, output=r"""%s""", 
        outputExpr=r"""%s""" )\'\'\')''' % ( parm, output, outputExpr) 
    # print opt
    return eval(opt)
#end


#file log/ref_expLD.log
#
# this is an example of observing decay of LD
from simuUtil import *
from simuRPy import *

simu = simulator(
        population(size=1000, ploidy=2, loci=[2]),
        randomMating(), rep=4 )

# see the change of allele/genotype/heplotype numbers as
# the result of genetic drift.
init = initByValue([1,2,2,1])
count = stat(LD=[0,1])
recombine = recombinator( rate=0.1 )
simu.evolve([
    recombine, count,
    pyEval(r'"%.4f\t" % LD[0][1]'),
    endl(rep=REP_LAST),
    #varPlotter(expr='LD[0][1]', title='Linkage disequilibrium',
    #    numRep = 4, ytitle='LD', saveAs='LD')
    ], preOps=[init],
    end=10
)

#end



#file log/ref_expcomplex.log
#

numSubPop = 100         # number of archipelagos
numFamilies = 10        # real simulation uses 1000
numOffspring = 4     # kind of family size
numReplicate = 1
loci = [20]*20            # 400 loci on 20 chromosomes
endGen = 10                 # should be at leat 1000
maxAllele = 30
mutationRate = 0.001
recombinationRate = 0.02

popSize = numFamilies*numOffspring*numSubPop
subPopSize = [numFamilies*numOffspring]*numSubPop

# intializer
init = initByFreq( alleleFreq=[1./maxAllele]*maxAllele )

# migration: island model
#     by proportion, .1 to all others
#
migrRate = .1
# rate[i->i] will be ignored so we can do the following
migrRates = [[migrRate/(numSubPop-1)]*numSubPop]*numSubPop 
migrMode = MigrByProbability
#
migrate = migrator(migrRates, mode=migrMode)

# mutation
mutate = kamMutator(rate=mutationRate, maxAllele=maxAllele)

# recombination
recombine = recombinator( rate = recombinationRate )

# create a simulator 
simu = simulator(
    population(size=popSize, ploidy=2, loci=loci,
        subPop=subPopSize),
    randomMating(numOffspring = numOffspring,
                         newSubPopSize=subPopSize) )
#
# evolve
simu.evolve([
    migrate, 
    recombine, 
    mutate,
    pyEval(r"gen", rep=0),    # report progress
    endl(rep=REP_LAST)
    ],
    preOps=[init],
    end=endGen)


#end

#file log/ref_expmigration.log
# this is an example of complex population size change.
# for endl and tab
from simuUtil import *

#number of cities
nc = 6

# how to change subpop size?
def changeSPSize(gen, oldSize):
    size = [0]*len(oldSize)
    for i in range(0, len(size)):
        size[i] = oldSize[i]*1.2
    if size[i] > 1000:
        size[i] /= 2
    return size

# migration between subpopulaitons
rates = []
for i in range(nc):
    rates.append([0.]*nc)
#
for i in range(1,nc-1):
    rates[i][i+1]=0.05
    rates[i][i-1]=0.05

#
rates[0][1] = 0.1
rates[nc-1][nc-2] = 0.1

# print rates
print rates
migr = migrator(rate=rates, mode=MigrByProbability)

# initially, we need to set everyone to middle subpop
initMigr = migrator(rate=[[1]], mode=MigrByProportion,
    fromSubPop=[0], toSubPop=[nc/2])

pop = population(size=500)

# the new popsize relies on a variable newSPSize
# which is calculated from subPopSize bu newSize operator
simu = simulator(pop,
    randomMating(newSubPopSizeFunc=changeSPSize) )

# evolve!
simu.evolve( [
    migr, stat(popSize=True),
    pyEval('list(subPopSize)'), endl()],
    preOps = [ initMigr ], end=10
)

#end

# need reich.py
#PS /bin/cp -f ../examples/Reich2002/reich.py log/ref_reich.py
