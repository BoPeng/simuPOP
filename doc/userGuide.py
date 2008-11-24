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

#
# make sure each run generates the same output to avoid unnecessary
# documentation changes.
rng().setSeed(12345)

#file log/standard.log
from simuPOP import *
pop = population(10, loci=[2])
pop.locusPos(10)
pop.individual(20).setAllele(1, 0)
#end

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
    preOps = [initByValue([1,2,2,1])],  
    ops = [
        recombinator(rate=0.01),
        stat(LD=[0,1]),
        pyEval(r"'%.2f\t' % LD[0][1]", step=10),
        pyOutput('\n', rep=-1, step=10)
    ],
    gen=100
)
#end

#file log/help.log
help(population.addInfoField)
#end

#file log/absIndex.log
print pop.chromLocusPair(7)
print pop.absLocusIndex(1,1)
#end

#file log/genoStru.log
# create a population, most parameters have default values
pop = population(size=[2,3], ploidy=2, loci=[5,10],
    lociPos=[range(0,5),range(0,20,2)],
    alleleNames=['A','C','T','G'], maxAllele=3)
print pop.popSize()
print pop.ploidy()
print pop.ploidyName()
print pop.numChrom()
print pop.locusPos(2)
print pop.alleleName(1)
#end

#file log/indGenoStru.log
# get the fourth individual of the population
ind = pop.individual(3)
# access genotypic structure info
print ind.ploidy()
print ind.numChrom()
print ind.numLoci(0)
print ind.genoSize()
# and from simulator level
simu = simulator(pop, randomMating(), rep=3)
print simu.numChrom()
#end



#file log/popInit.log
# a Wright-Fisher population
WF = population(size=100, ploidy=1, loci=[1])

# a diploid population of size 10
# there are two chromosomes with 5 and 7 loci respectively
pop = population(size=[2, 8], ploidy=2, loci=[5, 7])

# a population with SNP markers (with names A,C,T,G
# range() are python functions
pop = population(size=[2,3], ploidy=2, loci=[5,10],
    lociPos=[range(0,5),range(0,20,2)],
    alleleNames=['A','C','T','G'], maxAllele=3)
#end

#file log/popAndOperator.log
simu = simulator(pop, randomMating(), rep=3)
simu.evolve(
  preOps = [ initByFreq([.8, .2])],
  ops = [
    stat(alleleFreq=[0,1], Fst=[1], step=10),
    kamMutator(rate=0.001, rep=1),
    kamMutator(rate=0.0001, rep=2)
  ],
  gen=10
)
#end


#file log/InitByFreq.log
def InitByFreq(pop, *args, **kwargs):
  initByFreq(*args, **kwargs).apply(pop)

InitByFreq(pop, [.2, .3, .4, .1])
#end



#file log/dumpPop.log
pop = population(size=[2,3], ploidy=2, loci=[5,10],
    lociPos=[range(0,5),range(0,20,2)],
    alleleNames=['A','C','T','G'], maxAllele=3)
# .apply form
initByFreq([.2, .3, .4, .1]).apply(pop)
# function form
Dump(pop)
#end



#file log/popStru.log
print pop.popSize()
print pop.numSubPop()
print pop.subPopSize(0)
print pop.subPopSizes()
print pop.subPopBegin(1)
print pop.subPopEnd(1)
print pop.subPopIndPair(3)
print pop.absIndIndex(1,1)
#end

#file log/popStruManip.log
pop.setIndSubPopID([1,2,2,3,1])
pop.setSubPopByIndID()
pop.removeLoci(keep=range(2,7))
Dump(pop)
#end



#file log/ind.log
# get an individual
ind = pop.individual(9)
# oops, wrong index
ind = pop.individual(3)
# you can access genotypic structure info
print ind.ploidy()
print ind.numChrom()
# ...
# as well as genotype
print ind.allele(1) 
ind.setAllele(1,5)
print ind.allele(1)
# you can also use an overloaded function
# with a second parameter being the ploidy index
print ind.allele(1,1) # second locus at the second copy of chromosome
# other information
print ind.affected()
print ind.affectedChar()
ind.setAffected(1)
print ind.affectedChar()
print ind.sexChar()
#end



#file log/popVars.log
from simuUtil import ListVars
ListVars(pop.vars(), useWxPython=False)
Stat(pop, popSize=1, alleleFreq=[0])
# subPop is True by default, use name to limit the variables to display
ListVars(pop.vars(), useWxPython=False, subPop=False, name='alleleFreq')
# print number of allele 1 at locus 0
print pop.vars()['alleleNum'][0][1]
print pop.dvars().alleleNum[0][1]
print pop.dvars().alleleFreq[0]
#end


#file log/localNamespace.log
print pop.evaluate('alleleNum[0][1] + alleleNum[0][2]')
pop.execute('newPopSize=int(popSize*1.5)')
ListVars(pop.vars(), level=1, useWxPython=False)
# this variable is 'local' to the population and is
# not available in the main namespace
newPopSize
#end


#file log/ancestralPop.log
simu = simulator(population(10000, loci=[2]), randomMating())
simu.evolve(
  ops = [
    setAncestralDepth(5, at=[-5]),
    kamMutator(rate=0.01, atLoci=[0], maxAllele=1),
    kamMutator(rate=0.001, atLoci=[1], maxAllele=1)
  ],
  end = 20
)
pop = simu.population(0)
# start from current generation
for i in range(pop.ancestralDepth()+1):
  pop.useAncestralPop(i)
  Stat(pop, alleleFreq=[0,1])
  print '%d   %5f   %5f' % (i, pop.dvars().alleleFreq[0][1], pop.dvars().alleleFreq[1][1])

# restore to the current generation  
pop.useAncestralPop(0)  
#end



#turnOnDebug(DBG_SIMULATOR)
#turnOnDebug(DBG_UTILITY)

#file log/expr.log
simu = simulator(population(10),noMating(), rep=2)
# evaluate an expression in different areas
print simu.vars(0)
print simu.population(0).evaluate("gen+1")
# a statement (no return value)
simu.population(0).execute("myRep=2+rep*rep")
simu.population(1).execute("myRep=2*rep")
print simu.vars(0)
#end
#file log/expreval.log
simu.step([ pyExec("myRep=2+rep*rep") ])
print simu.vars(0)
#end


# Note that I can not use pop now, since it is
# obtained from simu.population(0) which is invalid now.

#file log/popSaveLoad.log
# save it in various formats, default format is "txt"
pop = population(1000, loci=[2, 5, 10])
pop.savePopulation("pop.txt")
pop.savePopulation("pop.txt", compress=False)
pop.savePopulation("pop.xml", format="xml")
pop.savePopulation("pop.bin", format="bin")

# load it in another population
pop1 = LoadPopulation("pop.xml", format="xml")
#end

#file log/saveFstat.log
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


#file log/randomSample.log
# random sample
# [0]: RandomSample already return
#  a list of samples even if times=1 (default)
Dump( RandomSample(pop, 3)[0])
#end

#file log/calcStat.log
Stat(pop, popSize=1, alleleFreq=range(0, pop.totNumLoci()),
  heteroFreq=range(0,pop.totNumLoci()), Fst=[0])
#end


#file log/generator.log
def func():
    i = 1
    all = 0
    while i < 10:
        all += 1./i
        i += 1
        yield all 

a = func()
a.next()
a.next()
for i in a:
    print '%.3f' % i,

#end

#file log/generator_random.log
def randomChooser(pop, sp):
    males = [x for x in range(pop.subPopSize(sp)) \
        if pop.individual(x, sp).sex() == Male]
    females = [x for x in range(pop.subPopSize(sp)) \
        if pop.individual(x, sp).sex() == Female]
    nm = len(males)
    nf = len(females)
    while True:
        yield males[rng().randInt(nm)], females[rng().randInt(nf)]

pop = population(size=[100, 20], loci=[1])
# this will initialize sex randomly
InitByFreq(pop, [0.2, 0.8])
rc1 = randomChooser(pop, 0)
for i in range(5):
    print rc1.next(),

rc2 = randomChooser(pop, 1)
for i in range(5):
    print rc2.next(),

#end

#file log/pyMating.log
simu = simulator(pop,
    pyMating(pyParentsChooser(randomChooser), 
        mendelianOffspringGenerator()))
simu.step()
#end

#file log/heteroMating.log
TurnOnDebug(DBG_MATING)
pop = population(100, loci=[2])
pop.setVirtualSplitter(proportionSplitter([0.2, 0.8]))
simu = simulator(pop, heteroMating(
    [selfMating(numOffspring=5, subPop=0, virtualSubPop=0),
    randomMating(numOffspring=20, subPop=0, virtualSubPop=1)]))

simu.step()
print simu.dvars(0).famSizes
TurnOffDebug(DBG_MATING)
#end


#file log/operatorstages.log
d = dumper()
print d.canApplyPreMating()
print d.canApplyDuringMating()
# so dumper is a post mating operator
print d.canApplyPostMating()
#end

#file log/operatorgen.log
simu = simulator(population(1),binomialSelection(), rep=3)
op1 = output("a", begin=5, end=20, step=3)
op2 = output("a", begin=-5, end=-1, step=2)
op3 = output("a", at=[2,5,10])
op4 = output("a", at=[-10,-5,-1])
simu.evolve( [ pyEval(r"str(gen)+'\n'", begin=5, end=-1, step=2)],
               gen=10)
#end

#file log/operatoroutput.log
simu = simulator(population(100), randomMating(), rep=2)
simu.step(
  preOps=[
    initByFreq([0.2, 0.8], rep=0),
    initByFreq([0.5, 0.5], rep=1) ],
  ops = [
    stat(alleleFreq=[0]),
    pyEval('alleleFreq[0][0]', output='a.txt')
  ]
)
# only from rep 1
print open('a.txt').read()

simu.step(
  ops = [
    stat(alleleFreq=[0]),
    pyEval('alleleFreq[0][0]', output='>>a.txt')
  ])
# from both rep0 and rep1
print open("a.txt").read()

outfile='>>>a.txt'
simu.step(
  ops = [
    stat(alleleFreq=[0]),
    pyEval('alleleFreq[0][0]', output=outfile),
    output("\t", output=outfile),
    output("\n", output=outfile, rep=0)
  ],
)
print open("a.txt").read()
#end

os.remove('a.txt')

#file log/operatoroutputexpr.log
outfile="'>>a'+str(rep)+'.txt'"
simu.step(
  ops = [
    stat(alleleFreq=[0]),
    pyEval('alleleFreq[0][0]', outputExpr=outfile)
  ]
)
print open("a0.txt").read()
print open("a1.txt").read()
#end

#file log/simulatorsaveload.log
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

#file log/pyOperator.log
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
      KamMutate(pop, maxAllele=2, rate=mu1, atLoci=[i])
    else:
      KamMutate(pop, maxAllele=2, rate=mu2, atLoci=[i])
  return True
#end

#file log/pyOperatorUse.log
pop = population(size=10000, ploidy=2, loci=[2, 3])

simu = simulator(pop, randomMating())

simu.evolve(
  preOps = [ 
    initByFreq( [.6, .4], atLoci=[0,2,4]),
    initByFreq( [.8, .2], atLoci=[1,3]) ],
  ops = [ 
    pyOperator( func=dynaMutator, param=(.5, .1, 0) ),
    stat(alleleFreq=range(5)),
    pyEval(r'"%f\t%f\n"%(alleleFreq[0][1],alleleFreq[1][1])', step=10)
    ],
  end = 30
)        
#end


#file log/initByFreq.log
simu = simulator( population(size=[2,3], loci=[5,7]),
    randomMating(), rep=1)
simu.step([
    initByFreq(alleleFreq=[ [.2,.8],[.8,.2]]),
    dumper(alleleOnly=True)
  ])
#end

#file log/initByValue.log
simu.step([
    initByValue([1]*5 + [2]*7 + [3]*5 +[4]*7),
    dumper(alleleOnly=True)])
#end

#file log/pyInit.log
def initAllele(ind, p, sp):
  return sp + ind + p

simu.step([
    pyInit(func=initAllele),
    dumper(alleleOnly=True, dispWidth=2)])
#end


#file log/pyMigrator.log
simu = simulator(population(size=[2,3], loci=[2,5]),
    randomMating())
# an Numeric array, force to Int type
spID = [2,2,1,1,0]
simu.step( [
    initByFreq([.2,.4,.4]), 
    dumper(alleleOnly=True, stage=PrePostMating),
    pyMigrator(subPopID=spID)
    ])
#end

#file log/kamMutator.log
simu = simulator(population(size=5, loci=[3,5]), noMating())
simu.step([
    kamMutator( rate=[.2,.6,.5], atLoci=[0,2,6], maxAllele=9),
    dumper(alleleOnly=True)])
#end

#file log/smmMutator.log
simu = simulator(population(size=3, loci=[3,5]), noMating())
simu.step([
    initByFreq( [.2,.3,.5]),
    smmMutator(rate=1,  incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)])
#end

#file log/gsmMutator.log
simu.step([
    initByFreq( [.2,.3,.5]),
    gsmMutator(rate=1, p=.8, incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)])

import random
def rndInt():
  return random.randrange(3,6)

simu.step([
    initByFreq( [.2,.3,.5]),
    gsmMutator(rate=1, func=rndInt, incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)])

#end

#file log/pyMutator.log
def mut(x):
  return 8

simu.step([
  pyMutator(rate=.5, atLoci=[3,4,5], func=mut),
  dumper(alleleOnly=True)])
#end


#file log/recombinator.log
simu = simulator(population(4, loci=[4,5,6], 
    infoFields=['father_idx', 'mother_idx']),
    randomMating())
simu.step([
  parentsTagger(),
  ],
  preOps = [initByFreq([.2,.2,.4,.2]), dumper(alleleOnly=True) ],
  postOps = [ dumper(alleleOnly=True)]
)
simu.step([
  parentsTagger(),
  recombinator(rate=[1,1,1], afterLoci=[2,6,10])
  ],
  postOps = [ dumper(alleleOnly=True)]
)
#end


#file log/basicSelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[1], infoFields=['fitness']),
    randomMating())
s1 = .1
s2 = .2
simu.evolve([
    stat( alleleFreq=[0], genoFreq=[0]),
    mapSelector(locus=0, fitness={'0-0':(1-s1), '0-1':1, '1-1':(1-s2)}),
    pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
    ],
    preOps=[  initByFreq(alleleFreq=[.2,.8])],
    gen=300)
#end

#file log/pySelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[3], infoFields=['fitness'] ),
    randomMating())

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
print sel([1,1])

simu.evolve([
    stat( alleleFreq=[0], genoFreq=[0]),
    pySelector(loci=[0,1],func=sel),
    pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=25)
    ],
    preOps=[  initByFreq(alleleFreq=[.2,.8])],
    gen=100)
#end

#file log/pySubset.log
simu = simulator(population(size=[2,3], loci=[3,4], infoFields=['fitness']),
    randomMating())
simu.step([
    initByFreq([.3,.5,.2]),
    pySubset( [1,-1,-1,1,-1] ),
    dumper(alleleOnly=True, stage=PrePostMating)
   ])
#end



#turnOnDebug(DBG_ALL)
#turnOnDebug(DBG_SIMULATOR)

#file log/varPlotter.log
from simuUtil import *
from simuRPy import *

simu = simulator( population(size=[50,50,100], ploidy=2, loci=[3,4]),
    randomMating(), rep=4)

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
     varDim=3, win=10, title='subPop size', saveAs='log/simuDemo')
   ],
   gen=30
)

#end
#PS /usr/bin/convert log/simuDemo16.eps log/simuDemo16.png
#PS /bin/rm -f log/simuDemo*.eps


#file log/rng.log
print ListAllRNG()
print rng().name()
SetRNG("taus2", seed=10)
print rng().name()
#end

#file log/rngrand.log
r=rng()
#help(RNG)
for n in range(1,10):
  print r.randBinomial(10, .7),
#end




#file log/extgenostru.log
pop = population(1, loci=[2,3,4])
print pop.numLoci(1)
print pop.locusPos(2)
dis = pop.arrLociPos()
print dis
dis[2] = 0.5
print pop.locusPos(2)
print pop.arrLociPos()
#end

#file log/extgenotype.log
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

#file log/extother.log
print ind.sex()
print ind.sexChar()
ind.setSex(Female)
ind.setAffected(True)
print ind.tag()
ind.setTag([1,2])
Dump(pop)
#end

#file log/extsimu.log
simu = simulator(pop, randomMating(), rep=3)
pop1 = simu.population(1)
ind1 = pop1.individual(0)
ind1.setAllele(3,0)
Dump(pop1)
#end

#file log/extoperator.log
simu.step([ pyEval(stmts="pop=simu.population(rep)")])
#end


#file log/tab.log
def tab(**kwargs):
  parm = ''  
  for (k,v) in kwargs.items():
    parm += ' , ' + str(k) + '=' + str(v)
  cmd = r'output( """\t""" ' + parm + ')'
  # print cmd
  return eval(cmd)
#end

#file log/saveFstat.tmp
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
      p1 = 2*gs*ind        # begining of first hemo copy
      p2 = 2*gs*ind + gs   # second
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
#PS head -15 log/saveFstat.tmp > log/saveInFstatFormat.log
#PS rm -f log/saveFstat.tmp

#file log/saveFstat.log
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


#file log/expLD.log
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
   endl(rep=-1),
   #varPlotter(expr='LD[0][1]', title='Linkage disequilibrium',
   #  numRep = 4, ytitle='LD', saveAs='LD')
   ], preOps=[init],
   gen=10
)

#end



#file log/expcomplex.log
#

numSubPop = 100     # number of archipelagos
numFamilies = 10    # real simulation uses 1000
numOffspring = 4   # kind of family size
numReplicate = 1
loci = [20]*20      # 400 loci on 20 chromosomes
endGen = 10         # should be at leat 1000
maxAllele = 30
mutationRate = 0.001
recombinationRate = 0.02

popSize = numFamilies*numOffspring*numSubPop
subPopSize = [numFamilies*numOffspring]*numSubPop

# intializer
init = initByFreq( alleleFreq=[1./maxAllele]*maxAllele )

# migration: island model
#   by proportion, .1 to all others
#
migrRate = .1
# rate[i->i] will be ignored so we can do the following
migrRates = [[migrRate/(numSubPop-1)]*numSubPop]*numSubPop 
migrMode  = MigrByProbability
#
migrate =  migrator(migrRates, mode=migrMode)

# mutation
mutate = kamMutator(rate=mutationRate, maxAllele=maxAllele)

# recombination
recombine = recombinator( rate = recombinationRate )

# create a simulator 
simu = simulator(
    population(size=subPopSize, ploidy=2, loci=loci),
    randomMating(numOffspring = numOffspring,
                       newSubPopSize=subPopSize) )
#
# evolve
simu.evolve([
    migrate, 
    recombine, 
    mutate,
    pyEval(r"gen", rep=0),  # report progress
    endl(rep=-1)
    ],
    preOps=[init],
    gen=endGen)


#end

#file log/expmigration.log
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
simu.evolve(
  [migr, stat(popSize=True),
   pyEval('list(subPopSize)'), endl()],
  preOps = [ initMigr ], gen=10
  )

#end

# need reich.py
#PS /bin/cp -f ../examples/Reich2002/reich.py log/reich.py
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

#file log/importSimuPOPOpt.log
import simuOpt
simuOpt.setOptions(optimized=False, alleleType='long', quiet=True)
from simuPOP import *
# make sure each run generates the same output to avoid unnecessary
# documentation changes.
rng().setSeed(12345)
# remember that global functions start with captical letters
print AlleleType()
print Optimized()
#end


#file log/absIndex.log
pop = population(size=[20, 30], loci=[5, 6])
print pop.chromLocusPair(7)
print pop.absLocusIndex(1,1)
print pop.absIndIndex(10, 1)
print pop.subPopIndPair(40)
#end


#file log/popAndOperator.log
simu = simulator(pop, randomMating(), rep=3)
simu.evolve(
    preOps = [ initByFreq([.8, .2])],
    ops = [
        stat(alleleFreq=[0,1], Fst=[1], step=10),
        kamMutator(rate=0.001, rep=1),
        kamMutator(rate=0.0001, rep=2)
    ],
    gen=10
)
#end


#file log/InitByFreq.log
def InitByFreq(pop, *args, **kwargs):
    initByFreq(*args, **kwargs).apply(pop)

InitByFreq(pop, [.2, .3, .4, .1])
#end


#file log/carray.log
# obtain an object using a genotype function
pop = population(size=2, loci=[1,2])
InitByValue(pop, [1,2,3])
arr = pop.genotype()
# print and expression (just like list)
print arr
str(arr)
# count
arr.count(2)
# index 
arr.index(2)
# can read write
arr[0] = 0
# the underlying locus position is also changed
print pop.genotype()
# convert to list
arr.tolist()
# or simply
list(arr)
# compare to list directly
arr == [0, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]
# you can also convert and compare
list(arr) == [0, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]
# slice
arr[:] = [3,2,1]*4
# assign from another part
arr[1:3] = arr[3:5]
#end

#file log/genotype.log
pop = population(size=[3, 2], loci=[2])
# single allele access
for ind in pop.individuals(1):
    for marker in range(pop.totNumLoci()):
        ind.setAllele(marker % 2, marker, 0)
        ind.setAllele(marker % 2, marker, 1)
        print '%d %d ' % (ind.allele(marker, 0), ind.allele(marker, 1))

# batch access
ind = pop.individual(4)
geno = ind.genotype()
print geno
geno[2] = 3
print ind.genotype()
# direct modification of the underlying genotype
geno[2:4] = [3, 4]
print ind.genotype()
# set genotype
ind.setGenotype([2, 1])
print geno
# print genotypes of all individuals in the second subpopulation.
print pop.genotype(1)
#end

#file log/dumpPop.log
pop = population(size=[2,3], ploidy=2, loci=[5,10],
        lociPos=[range(0,5),range(0,20,2)],
        alleleNames=['A','C','T','G'], maxAllele=3)
# .apply form
initByFreq([.2, .3, .4, .1]).apply(pop)
# function form
Dump(pop)
#end


#file log/popVars.log
from simuUtil import ListVars
pop = population(size=[1000, 2000], loci=[1])
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


#file log/localNamespace.log
pop = population(size=[1000, 2000], loci=[1])
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
# a statement (no return value)
simu.population(0).execute("myRep=2+rep*rep")
simu.population(1).execute("myRep=2*rep")
print simu.vars(0)
#end


#file log/pyEval.log
simu = simulator(population(100, loci=[1]),
    randomMating(), rep=2)
simu.evolve(
    preOps = [initByFreq([0.2, 0.8])],
    ops = [ stat(alleleFreq=[0]),
        pyExec('myNum = alleleNum[0][0] * 2'),
        pyEval(r'"gen %d, rep %d, num %d, myNum %d\n"' \
            ' % (gen, rep, alleleNum[0][0], myNum)')
        ],
    gen=3
)
#end

#file log/info1.log
pop = population(10, infoFields=['a', 'b'])
aIdx = pop.infoIdx('a')
bIdx = pop.infoIdx('b')
for ind in pop.individuals():
    a = ind.info(aIdx)
    ind.setInfo(a+1, bIdx)

print pop.indInfo(bIdx)
#end

#file log/info2.log
pop = population(5, infoFields=['a', 'b'])
InfoExec(pop, 'import random\na=random.randint(2,10)')
InfoExec(pop, 'b=a+a*2')
InfoEval(pop, r"'(%.0f, %.0f) ' % (a, b)")

# this is wrong because 'c' is not available
InfoExec(pop, 'b=c+a')
# but we can also make use of population variables.
pop.vars()['c'] = 6
InfoExec(pop, 'b=c+a', usePopVars=True)
print pop.indInfo('b')
#end

#file log/ancestralPop.log
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

#file log/virtualSubPop.log
import random
pop = population(1000, loci=[2, 3], infoFields=['age'])
InitByFreq(pop, [0.2, 0.8])
for ind in pop.individuals():
    ind.setInfo(random.randint(0,5), 'age')

# split by age
pop.setVirtualSplitter(infoSplitter('age', values=[2,4]))
pop.virtualSubPopSize(0, 0)
pop.virtualSubPopName(0, 1)

# split by genotype
a = pop.setVirtualSplitter(
    genotypeSplitter(locus=2, alleles=[[0,1], [1,1]], phase=True))
pop.virtualSubPopSize(0, 0)
pop.virtualSubPopSize(0, 1)

for ind in pop.individuals(0, 0):
    assert ind.allele(2, 0) == 0 and ind.allele(2, 1) == 1

#end



#file log/saveFstat.log
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


#file log/randomSample.log
# random sample
# [0]: RandomSample already return
#    a list of samples even if times=1 (default)
Dump( RandomSample(pop, 3)[0])
#end


#file log/calcStat.log
Stat(pop, popSize=1, alleleFreq=range(0, pop.totNumLoci()),
    heteroFreq=range(0,pop.totNumLoci()), Fst=[0])
#end

#file log/simulatorsaveload.log
simu.saveSimulator("s.sim")
simu1 = LoadSimulator("s.sim", randomMating())
#end

# remove these files
os.remove('s.sim')
os.remove('a0.txt')
os.remove('a1.txt')

#file log/pyOperator.log
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

#file log/pyOperatorUse.log
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



#file log/pySubset.log
simu = simulator(population(size=[2,3], loci=[3,4], infoFields=['fitness']),
        randomMating())
simu.step([
        initByFreq([.3,.5,.2]),
        pySubset( [1,-1,-1,1,-1] ),
        dumper(alleleOnly=True, stage=PrePostMating)
     ])
#end

#file log/generator_random.log
from random import randint

def randomChooser(pop, sp):
    males = [x for x in range(pop.subPopSize(sp)) \
        if pop.individual(x, sp).sex() == Male \
            and pop.individual(x, sp).info('age') > 30]
    females = [x for x in range(pop.subPopSize(sp)) \
        if pop.individual(x, sp).sex() == Female \
            and pop.individual(x, sp).info('age') > 30]
    nm = len(males)
    nf = len(females)
    while True:
        yield males[randint(0, nm-1)], females[randint(0, nf-1)]

pop = population(size=[1000, 200], loci=[1], infoFields=['age'])
# this will initialize sex randomly
InitByFreq(pop, [0.2, 0.8])
for ind in pop.individuals():
    ind.setInfo(randint(0, 60), 'age')

rc1 = randomChooser(pop, 0)
for i in range(5):
    print rc1.next(),

rc2 = randomChooser(pop, 1)
for i in range(5):
    print rc2.next(),

#end



#turnOnDebug(DBG_ALL)
#turnOnDebug(DBG_SIMULATOR)

#file log/varPlotter.log
from simuUtil import *
from simuRPy import *

simu = simulator( population(size=[50,50,100], ploidy=2, loci=[3,4]),
    randomMating(), rep=4)

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
         varDim=3, win=10, title='subPop size', saveAs='log/simuDemo')
     ],
     gen=30
)

#end
#PS /usr/bin/convert log/simuDemo16.eps log/simuDemo16.png
#PS /bin/rm -f log/simuDemo*.eps

#file log/ifElse.log
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
    gen=endGen,
    dryrun=False
)
#end

#file log/rng.log
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




#file log/extgenostru.log
pop = population(1, loci=[2,3,4])
print pop.numLoci(1)
print pop.locusPos(2)
dis = pop.arrLociPos()
print dis
dis[2] = 0.5
print pop.locusPos(2)
print pop.arrLociPos()
#end

#file log/extgenotype.log
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

#file log/extother.log
print ind.sex()
print ind.sexChar()
ind.setSex(Female)
ind.setAffected(True)
Dump(pop)
#end

#file log/extsimu.log
simu = simulator(pop, randomMating(), rep=3)
pop1 = simu.population(1)
ind1 = pop1.individual(0)
ind1.setAllele(3,0)
Dump(pop1)
#end

#file log/extoperator.log
simu.step([ pyEval(stmts="pop=simu.population(rep)")])
#end


#file log/tab.log
def tab(**kwargs):
    parm = ''    
    for (k,v) in kwargs.items():
        parm += ' , ' + str(k) + '=' + str(v)
    cmd = r'output( """\t""" ' + parm + ')'
    # print cmd
    return eval(cmd)
#end

#file log/saveFstat.tmp
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
#PS head -15 log/saveFstat.tmp > log/saveInFstatFormat.log
#PS rm -f log/saveFstat.tmp

#file log/saveFstat.log
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


#file log/expLD.log
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
    endl(rep=-1),
    #varPlotter(expr='LD[0][1]', title='Linkage disequilibrium',
    #    numRep = 4, ytitle='LD', saveAs='LD')
    ], preOps=[init],
    gen=10
)

#end



#file log/expcomplex.log
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
    population(size=subPopSize, ploidy=2, loci=loci),
    randomMating(numOffspring = numOffspring,
                         newSubPopSize=subPopSize) )
#
# evolve
simu.evolve([
    migrate, 
    recombine, 
    mutate,
    pyEval(r"gen", rep=0),    # report progress
    endl(rep=-1)
    ],
    preOps=[init],
    gen=endGen)


#end

#file log/expmigration.log
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
    preOps = [ initMigr ], gen=10
)

#end

# need reich.py
#PS /bin/cp -f ../examples/Reich2002/reich.py log/reich.py
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

#
from simuPOP import *


#file log/src_genoStruTrait.log
# create a population, most parameters have default values
pop = population(size=5, ploidy=2, loci=[5,10],
    lociPos=[range(0,5),range(0,20,2)],
    alleleNames=['A','C','T','G'],
    subPop=[2,3], maxAllele=3)
print pop.popSize()
print pop.ploidy()
print pop.ploidyName()
print pop.numChrom()
print pop.locusPos(2)
print pop.alleleName(1)
# get the fourth individual of the population
ind = pop.individual(3)
# access genotypic structure info
print ind.ploidy()
print ind.numChrom()
print ind.numLoci(0)
print ind.genoSize()
# and from simulator level
simu = simulator(pop, randomMating(), rep=3)
print simu.numChrom()
#end

#file log/src_population.log
# use of population function
# a Wright-Fisher population
WF = population(size=100, ploidy=1, loci=[1])

# a diploid population of size 10
# there are two chromosomes with 5 and 7 loci respectively
pop = population(size=10, ploidy=2, loci=[5, 7], subPop=[2, 8])

# a population with SNP markers (with names A,C,T,G)
# range() are python functions
pop = population(size=5, ploidy=2, loci=[5,10],
    lociPos=[range(0,5),range(0,20,2)],
    alleleNames=['A','C','T','G'],
    subPop=[2,3], maxAllele=3)

#
# population structure functions
print pop.popSize()
print pop.numSubPop()
print pop.subPopSize(0)
print pop.subPopSizes()
print pop.subPopBegin(1)
print pop.subPopEnd(1)
print pop.subPopIndPair(3)
print pop.absIndIndex(1,1)

#
# functions of setting population structure
pop.setIndSubPopID([1,2,2,3,1])
pop.setSubPopByIndID()
pop.removeLoci(keep=range(2,7))
Dump(pop)

#
# save and load population
# save it in various formats, default format is "txt"
pop = population(1000, loci=[2, 5, 10])
pop.savePopulation("pop.txt")
pop.savePopulation("pop.txt", compress=False)
pop.savePopulation("pop.xml", format="xml")
pop.savePopulation("pop.bin", format="bin")

# load it in another population
pop1 = LoadPopulation("pop.xml", format="xml")
#end


#file log/src_individual.log
pop = population(500, loci=[2, 5, 10])
# get an individual
ind = pop.individual(9)
# oops, wrong index
ind = pop.individual(3)
# you can access genotypic structure info
print ind.ploidy()
print ind.numChrom()
# ...
# as well as genotype
print ind.allele(1) 
ind.setAllele(1,5)
print ind.allele(1)
# you can also use an overloaded function
# with a second parameter being the ploidy index
print ind.allele(1,1) # second locus at the second copy of chromosome
# other information
print ind.affected()
print ind.affectedChar()
ind.setAffected(1)
print ind.affectedChar()
print ind.sexChar()
#end


#file log/src_operator.log
simu = simulator(population(1),binomialSelection(), rep=2)
op1 = pyOutput("a", begin=5, end=20, step=3)
op2 = pyOutput("a", begin=-5, end=-1, step=2)
op3 = pyOutput("a", at=[2,5,10])
op4 = pyOutput("a", at=[-10,-5,-1])
simu.evolve( [ pyEval(r"str(gen)+'\n'", begin=5, end=-1, step=2)],
                             end=10)
#
#
# operator group
from simuUtil import *
simu = simulator(population(1),binomialSelection(), rep=4,
    grp=[1,2,1,2])
simu.step(
    ops = [ 
        pyEval(r"grp+3", grp=1),
        pyEval(r"grp+6", grp=2),
        pyOutput('\n', rep=-1)
    ]
)

#
# parameter output 
simu = simulator(population(100), randomMating(), rep=2)
simu.step(
    preOps=[
        initByFreq([0.2, 0.8], rep=0),
        initByFreq([0.5, 0.5], rep=1) ],
    ops = [
        stat(alleleFreq=[0]),
        pyEval('alleleFreq[0][0]', output='a.txt')
    ]
)
# only from rep 1
print open('a.txt').read()

simu.step(
    ops = [
        stat(alleleFreq=[0]),
        pyEval('alleleFreq[0][0]', output='>>a.txt')
    ])
# from both rep0 and rep1
print open("a.txt").read()

outfile='>>>a.txt'
simu.step(
    ops = [
        stat(alleleFreq=[0]),
        pyEval('alleleFreq[0][0]', output=outfile),
        pyOutput("\t", output=outfile),
        pyOutput("\n", output=outfile, rep=0)
    ],
)
print open("a.txt").read()
#
# Output expression
outfile="'>>a'+str(rep)+'.txt'"
simu.step(
    ops = [
        stat(alleleFreq=[0]),
        pyEval('alleleFreq[0][0]', outputExpr=outfile)
    ]
)
print open("a0.txt").read()
print open("a1.txt").read()
#end
os.remove('a.txt')



#file log/src_initByFreq.log
simu = simulator( 
    population(subPop=[2,3], loci=[5,7], maxAllele=1),
    randomMating(), rep=1)
simu.step([
    initByFreq(alleleFreq=[ [.2,.8],[.8,.2]]),
    dumper(alleleOnly=True)
  ])
#end

#file log/src_initByValue.log
simu = simulator(
    population(subPop=[2,3], loci=[5,7], maxAllele=9),
    randomMating(), rep=1)
simu.step([
    initByValue([1]*5 + [2]*7 + [3]*5 +[4]*7),
    dumper(alleleOnly=True)])
#end

#file log/src_pyInit.log
def initAllele(ind, p, sp):
  return sp + ind + p

simu = simulator( 
    population(subPop=[2,3], loci=[5,7]),
    randomMating(), rep=1)
simu.step([
    pyInit(func=initAllele),
    dumper(alleleOnly=True, dispWidth=2)])
#end


#file log/src_kamMutator.log
simu = simulator(population(size=5, loci=[3,5]), noMating())
simu.step([
    kamMutator( rate=[.2,.6,.5], loci=[0,2,6], maxAllele=9),
    dumper(alleleOnly=True)])
#end

#file log/src_smmMutator.log
simu = simulator(population(size=3, loci=[3,5]), noMating())
simu.step([
    initByFreq( [.2,.3,.5]),
    smmMutator(rate=1,  incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)])
#end

#file log/src_gsmMutator.log
simu = simulator(population(size=3, loci=[3,5]), noMating())
simu.step([
    initByFreq( [.2,.3,.5]),
    gsmMutator(rate=1, p=.8, incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)])

import random
def rndInt():
  return random.randrange(3,6)

simu.step([
    initByFreq( [.2,.3,.5]),
    gsmMutator(rate=1, func=rndInt, incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)
    ]
)
#end

#file log/src_pyMutator.log
def mut(x):
  return 8

simu = simulator(population(size=3, loci=[3,5]), noMating())
simu.step([
  pyMutator(rate=.5, loci=[3,4,5], func=mut),
  dumper(alleleOnly=True)])
#end


#file log/src_recombinator.log
simu = simulator(
    population(4, loci=[4,5,6], maxAllele=9,
    infoFields=['father_idx', 'mother_idx']),
    randomMating())
simu.step([
  parentsTagger(),
  ],
  preOps = [initByFreq([.2,.2,.4,.2]), dumper(alleleOnly=True) ],
  postOps = [ dumper(alleleOnly=True)]
)
simu.step([
  parentsTagger(),
  recombinator(rate=[1,1,1], afterLoci=[2,6,10])
  ],
  postOps = [ dumper(alleleOnly=True)]
)
#end


#file log/src_pyOperator.log
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



#file log/src_mating.log
# arbitrary demographic model
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

#
# control the number of offspring per mating event
# famSizes is only defined when DBG_MATING is defined
TurnOnDebug(DBG_MATING)
simu = simulator(population(50, loci=[1]),
    randomMating(numOffspring=2, 
        maxNumOffspring=5,
        mode=MATE_UniformDistribution))
simu.step(ops=[])
print simu.population(0).dvars().famSizes
TurnOffDebug(DBG_MATING)
#end


#file log/src_mapSelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[1], infoFields=['fitness']),
    randomMating())
s1 = .1
s2 = .2
simu.evolve([
    stat( alleleFreq=[0], genoFreq=[0]),
    mapSelector(locus=0, fitness={'0-0':(1-s1), '0-1':1, '1-1':(1-s2)}),
    pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
    ],
    preOps=[  initByFreq(alleleFreq=[.2,.8])],
    end=300)
#end


#file log/src_maSelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[1], infoFields=['fitness']),
    randomMating())
s1 = .1
s2 = .2
simu.evolve(
    preOps=[initByFreq(alleleFreq=[.2,.8])],
    ops = [
        stat( alleleFreq=[0], genoFreq=[0]),
        maSelector(locus=0, fitness=[1-s1, 1, 1-s2]),
        pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=100)
    ],
    end=300)
#end

#file log/src_pySelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[3], infoFields=['fitness'] ),
    randomMating())

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
print sel([1,1])

simu.evolve([
    stat( alleleFreq=[0], genoFreq=[0]),
    pySelector(loci=[0,1],func=sel),
    pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=25)
    ],
    preOps=[  initByFreq(alleleFreq=[.2,.8])],
    end=100)
#end

#file log/src_pySubset.log
simu = simulator(
    population(subPop=[2,3], loci=[3,4]),
    randomMating())
simu.step([
    initByFreq([.3, .5, .2]),
    pySubset([0, -1, -1, 1, -1] ),
    dumper(alleleOnly=True, stage=PrePostMating)
   ])
#end


#file log/src_mlPenetrance.log
pop = population(1000, loci=[3])
InitByFreq(pop, [0.3, 0.7])
pen = []
for loc in (0, 1, 2):
    pen.append( maPenetrance(locus=loc, wildtype=[1],
        penetrance=[0, 0.3, 0.6] ) )

# the multi-loci penetrance
MlPenetrance(pop, mode=PEN_Multiplicative, peneOps=pen)
Stat(pop, numOfAffected=True)
print pop.dvars().numOfAffected
#end


#file log/src_pyPenetrance.log
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
    func=peneFunc( ((0,0.5),(0.3,0.8)) ) )
#end



#turnOnDebug(DBG_ALL)
#turnOnDebug(DBG_SIMULATOR)

#file log/src_ifElse.log
simu = simulator(
    population(size=1000, loci=[1]),
    randomMating(), rep=4)
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
    ifElse(cond='numOfAffected==0',
        ifOp=pyEval(r'"No affected at gen %d\n" % gen'))
  ],
  end=50
)
#end

#file log/src_noneOp.log
# this may be set from command line option
savePop = False
# then, saveOp is defined accordingly
if savePop:
    saveOp = savePopulation(output='a.txt')
else:
    saveOp = noneOp()

simu = simulator(population(10), randomMating())
simu.step([saveOp])
#end


#file log/src_rng.log
print ListAllRNG()
print rng().name()
SetRNG("taus2", seed=10)
print rng().name()
#end

#file log/src_rngrand.log
r=rng()
#help(RNG)
for n in range(1,10):
    print r.randBinomial(10, .7),
#end


#file log/src_pyIndOperator.log
def indFunc(ind, param):
    ind.setInfo(ind.info('info2')+param[1], 'info2')
    print ind.info('info2')
    return True

simu = simulator(
    population(5, loci=[2,5], infoFields=['info1', 'info2']),
    randomMating(), rep=1)
simu.step(
    ops = [
        pyIndOperator(func=indFunc, param=(3,2)), 
    ],
)
#end


#file log/src_spread.log
pop = population(ploidy=2, loci=[2], subPop=[4, 4])
InitByFreq(pop, [.5, .5])
Spread(pop, 2, subPop=[1]),
Dump(pop)
#end


#file log/src_splitSubPop.log
pop = population(size=10, loci=[2,6])
InitByFreq(pop, [.2,.4,.4])
SplitSubPop(pop, which=0, sizes=[2,8], randomize=False)
print pop.subPopSizes()
#end


#file log/src_mlSelector.log
simu = simulator(
    population(size=10, ploidy=2, loci=[2], 
    infoFields=['fitness', 'spare']),
    randomMating())
simu.evolve(
    [ mlSelector([
         mapSelector(locus=0, fitness={'0-0':1,'0-1':1,'1-1':.8}),
         mapSelector(locus=1, fitness={'0-0':1,'0-1':1,'1-1':.8}),
         ], mode=SEL_Additive),
    ],
    preOps = [
        initByFreq(alleleFreq=[.2,.8])
    ],
    end=2
)
#end


#file log/src_mapPenetrance.log
pop = population(size=10, ploidy=2, loci=[2], subPop=[2, 8])
InitByFreq(pop, [.2, .8])
MapPenetrance(pop, locus=0, 
    penetrance={'0-0':0, '0-1':1, '1-1':1})
Stat(pop, numOfAffected=1)
#end


#file log/src_maPenetrance.log
pop = population(size=10, ploidy=2, loci=[2], subPop=[2, 8])
InitByFreq(pop, [.2, .8])
MaPenetrance(pop, locus=0, wildtype=0, penetrance=[0, 1, 1])
Stat(pop, numOfAffected=1)
#end
