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
#file log/importSimuPOP.log
from simuPOP import *
#end

#file log/testRPy.log
from simuRPy import *
r.seq(1,10,2)
r.plot([1,2,3], type='l')
#end


# we assume that ../scripts/ is the path to simuLDDecay.py
import sys

#PS ../scripts/simuLDDecay.py --help > log/LDDecayHelp.log
#PS ../scripts/simuLDDecay.py --noDialog --size=1000 --endGen=50 --recRate=0.01 --replicate=5 --saveFigure=LDDecay 
#PS mv LDDecay50.eps log/LDDecay.eps
#PS convert -rotate 90 log/LDDecay.eps log/LDDecay.png
#PS /bin/rm -f log/LDDecay.eps
#PS ../scripts/simuComplexDisease.py --help > log/simuComplexDisease.log


#file log/importSimuPOPOpt.log
import simuOpt
simuOpt.setOptions(optimized=False, longAllele=True)
from simuPOP import *
#end

#file log/addSysPath.log
import sys
sys.path.append('/path/to/simuPOP')
from simuPOP import *
#end

#file log/simpleExample.log
from simuPOP import *
from simuUtil import *
simu = simulator(
    population(size=1000, ploidy=2, loci=[2]),
    randomMating(),
    rep = 3)
simu.evolve(
    preOps = [initByValue([1,2,2,1]) ],  
    ops = [
        recombinator( rate=0.1),
        stat( haploFreq=[[0,1]]),
        pyEval(r"'%5.0f\t' % haploNum['0-1']['1-2']"),
        endl(rep=REP_LAST)
        ],
    end=5
)

#end



#file log/genoStru.log
# create a population, most parameters have default values
pop = population(size=5, ploidy=2, loci=[5,10],
    lociDist=[range(0,5),range(0,20,2)],
    alleleNames=['_','A','C','T','G'],
    subPop=[2,3], maxAllele=4)
print pop.popSize()
print pop.ploidy()
print pop.ploidyName()
print pop.numChrom()
print pop.locusDist(2)
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

#file log/absIndex.log
print pop.chromLocusPair(7)
print pop.absLocusIndex(1,1)
#end

#file log/helpPopInit.log
help(population.__init__)
#end

#file log/popInit.log
# a Wright-Fisher population
WF = population(size=100, ploidy=1, loci=[1])

# a diploid population of size 10
# there are two chromosomes with 5 and 7 loci respectively
pop = population(size=10, ploidy=2, loci=[5, 7], subPop=[2, 8])

# a population with SNP markers (with names A,C,T,G
#  range() are python functions
pop = population(size=5, ploidy=2, loci=[5,10],
    lociDist=[range(0,5),range(0,20,2)],
    alleleNames=['_','A','C','T','G'],
    subPop=[2,3], maxAllele=4)
InitByFreq(pop, [.25]*4)
#end


#file log/popSaveLoad.log
# save it in various formats, default format is "txt"
pop.savePopulation("pop.txt")
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


#file log/InitByFreq.log
def InitByFreq(pop, *args, **kwargs):
  initByFreq(*args, **kwargs).apply(pop)

InitByFreq(pop, [.2, .3, .4, .1])
#end


#file log/dumpPop.log
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

#file log/randomSample.log
# random sample
# [0]: RandomSample already return
#  a list of samples even if times=1 (default)
Dump( RandomSample(pop, 3)[0])
#end

#file log/popVars.log
from simuUtil import listVars
listVars(pop.vars(), useWxPython=False)
Stat(pop, popSize=1, alleleFreq=[0])
listVars(pop.vars(), useWxPython=False)
# print number of allele 1 at locus 0
print pop.vars()['alleleNum'][0][1]
print pop.dvars().alleleNum[0][1]
#end

#file log/localNamespace.log
print pop.evaluate('alleleNum[0][1] + alleleNum[0][2]')
pop.execute('newPopSize=int(popSize*1.5)')
listVars(pop.vars(), level=1, useWxPython=False)
# this variable is 'local' to the population and is
# not available in the main namespace
newPopSize
#end

 
#turnOnDebug(DBG_SIMULATOR)
#turnOnDebug(DBG_UTILITY)

#file log/expr.log
simu = simulator(population(10),noMating(), rep=2)
# evaluate an expression in different areas
print simu.vars(0)
print simu.population(0).evaluate("grp*2")
print simu.population(1).evaluate("grp*2")
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

#file log/calcStat.log
Stat(pop, popSize=1, alleleFreq=range(0, pop.totNumLoci()),
  heteroFreq=range(0,pop.totNumLoci()), Fst=[0])
#end

#file log/calcFstH.log
def calc_Fst_H(pop, alleles):
  """ calculate expected heterozygosities at given loci
    Formula etc please refer to user's manual
  """
  s = pop.dvars()
  if len(alleles) == 0:
    raise exceptions.ValueError("Please specify alleles on which to calculate Fst_H")

  # FIXME: need better interface.
  for l in alleles:   # of form [locus, allele, allele ...]
    if (type(l) != type([]) and type(l) != type(())) or len(l) <= 1:
      raise exceptions.ValueError("Format [ [ locus, allele,...]. [...] ]");
    
    s.Fst_H = {}
    s.Fis_H = {}
    s.Fit_H = {}
    loc = l[0]
    for ale in l[1:]:
      # calculate Fst_H for each loc, ale pair.
      # H_I based on observed heterozygosities in individuals in subpopulations
      H_I = 0.0
      H_S = 0.0
      for sp in range(0, s.numSubPop):
        H_I = H_I + s.subPopSize[sp]*s.subPop[sp]['heteroFreq'][loc][ale]
        H_S = H_S + s.subPopSize[sp]*s.subPop[sp]['heteroFreq'][loc][0]
      H_I = H_I / s.popSize
      H_S = H_S / s.popSize
      H_T = s.heteroFreq[loc][0]
      s.Fst_H['%d-%d' % (loc,ale)] = (H_T - H_S)/H_T
      s.Fis_H['%d-%d' % (loc,ale)] = (H_S - H_I)/H_S
      s.Fit_H['%d-%d' % (loc,ale)] = (H_T - H_I)/H_T
#end

#file log/wrapFstH.log
def Fst_H(alleles,**kwargs):
  parm = ''  
  for (k,v) in kwargs.items():
    parm += ' , ' + str(k) + '=' + str(v)
  #  calc_Fst_H(loci= loci?, rep=rep)
  cmd = r'pyExec( exposePop=1, stmts=r"""calc_Fst_H(pop=pop, alleles= ' + \
    str(alleles) + ')""", %s)' % parm
  # print cmd
  return eval( cmd )
#end

#file log/useFstH.log
simu.apply([ initByFreq([.3,.5,.2]), 
  stat(popSize=1, heteroFreq=[0]), 
  Fst_H([[0,1,2]]) ] )
listVars(simu.vars(0), level=1, useWxPython=False)
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
               end=10)
#end


#file log/operatorgrp.log
from simuUtil import *
simu = simulator(population(1),binomialSelection(), rep=4,
                 grp=[1,2,1,2])
simu.apply([ pyEval(r"grp+3", grp=1),
             pyEval(r"grp+6", grp=2), tab(), endl() ])
#end

#file log/operatoroutput.log
simu = simulator(population(100),randomMating(), rep=2)
simu.step([ stat(alleleFreq=[0], output=">"), 
    tab(), endl() ],
    preOps=[initByFreq([0.2,.8])])

# output to another file (I use outfile since this is used by tab and endl
outfile="a.txt"
simu.step([ stat(alleleFreq=[0], output=outfile), 
    output("\t", output=outfile),
    output("\n", output=outfile, rep=REP_LAST) ],
    preOps=[initByFreq([0.2,.8])])
print open("a.txt").read()
outfile=">>a.txt"
simu.step([ stat(alleleFreq=[0], output=outfile), 
    output("\t", output=outfile),
    output("\n", output=outfile, rep=REP_LAST) ],
    preOps=[initByFreq([0.2,.8])])
print open("a.txt").read()
#end

os.remove('a.txt')

#file log/operatoroutputexpr.log
outfile="'>>a'+str(rep)+'.txt'"
simu.step([ stat(alleleFreq=[0], outputExpr=outfile), 
    output("\n", outputExpr=outfile) ],
    preOps=[initByFreq([0.2,.8])])
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
from simuRPy import *
import time

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

pop = population(size=10000, ploidy=2, loci=[2, 3])    
simu = simulator(pop, randomMating())
simu.evolve(
  preOps = [ 
    initByFreq( [.6, .4], atLoci=[0,2,4]),
    initByFreq( [.8, .2], atLoci=[1,3]) ],
  ops = [ 
    pyOperator( func=dynaMutator, param=(.5, .1, 0) ),
    stat(alleleFreq=range(5)),
    varPlotter('(alleleFreq[0][2],alleleFreq[1][2])', 
      title='allele frequency at loci 0 and 1',
      update=5, varDim=2),
    ],
  end = 30
)        
print "The R window will be closed after five seconds..."
time.sleep(5)


#end


#file log/initByFreq.log
help(initByFreq.__init__)
simu = simulator( population(subPop=[2,3], loci=[5,7]),
    randomMating(), rep=1)
simu.apply([
    initByFreq(alleleFreq=[ [.2,.8],[.8,.2]]),
    dumper(alleleOnly=True)
  ])
#end

#file log/initByValue.log
simu.apply([
    initByValue([1]*5 + [2]*7 + [3]*5 +[4]*7),
    dumper(alleleOnly=True)])
#end

#file log/pyInit.log
help(pyInit.__init__)
def initAllele(ind, p, sp):
  return sp + ind + p

simu.apply([
    pyInit(func=initAllele),
    dumper(alleleOnly=True, dispWidth=2)])
#end

#file log/migratorhelp.log
help(migrator.__init__)
#end

#file log/pyMigrator.log
help(pyMigrator.__init__)

simu = simulator(population(subPop=[2,3], loci=[2,5]),
    randomMating())
# an Numeric array, force to Int type
spID = carray('i',[2,2,1,1,0])
simu.apply( [
    initByFreq([.2,.4,.4]), 
    dumper(alleleOnly=True, stage=PrePostMating),
    pyMigrator(subPopID=spID)
    ])

#end

#file log/mutatorhelp.log
help(mutator.__init__)
#end

#file log/kamMutator.log
simu = simulator(population(size=5, loci=[3,5]), noMating())
simu.apply([
    kamMutator( rate=[.2,.6,.5], atLoci=[0,2,6], maxAllele=9),
    dumper(alleleOnly=True)])
#end

#file log/smmMutator.log
simu = simulator(population(size=3, loci=[3,5]), noMating())
simu.apply([
    initByFreq( [.2,.3,.5]),
    smmMutator(rate=1,  incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)])
#end

#file log/gsmMutator.log
simu.apply([
    initByFreq( [.2,.3,.5]),
    gsmMutator(rate=1, p=.8, incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)])

import random
def rndInt():
  return random.randrange(3,6)

simu.apply([
    initByFreq( [.2,.3,.5]),
    gsmMutator(rate=1, func=rndInt, incProb=.8),
    dumper(alleleOnly=True, stage=PrePostMating)])

#end

#file log/pyMutator.log
def mut(x):
  return 8

simu.apply([
  pyMutator(rate=.5, atLoci=[3,4,5], func=mut),
  dumper(alleleOnly=True)])
#end

#file log/recombinatorhelp.log
help(recombinator.__init__)
#end

#file log/recombinator.log
simu = simulator(population(4, loci=[4,5,6]),
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

#file log/selectorhelp.log
help(basicSelector.__init__)
#end

#file log/basicSelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[1]),
    randomMating())
s1 = .1
s2 = .2
simu.evolve([
    stat( alleleFreq=[0], genoFreq=[0]),
    basicSelector(locus=0, fitness={'1-1':(1-s1), '1-2':1, '2-2':(1-s2)}),
    pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=10)
    ],
    preOps=[  initByFreq(alleleFreq=[.2,.8])],
    end=300)
#end

#file log/pySelector.log
simu = simulator(
    population(size=1000, ploidy=2, loci=[3]),
    randomMating() )

s1 = .2
s2 = .3
def sel(arr):
  if arr[0] == 1 and arr[1] == 1:
    return 1 - s1
  elif arr[0] == 1 and arr[1] == 2:
    return 1
  elif arr[0] == 2 and arr[1] == 1:
    return 1
  else:
    return 1 - s2


# test func
print sel(carray('B',[1,1]))

simu.evolve([
    stat( alleleFreq=[0], genoFreq=[0]),
    pySelector(loci=[0,1],func=sel),
    pyEval(r"'%.4f\n' % alleleFreq[0][1]", step=10)
    ],
    preOps=[  initByFreq(alleleFreq=[.2,.8])],
    end=100)
#end

#file log/pySubset.log
simu = simulator(population(subPop=[2,3], loci=[3,4]),
    randomMating())
simu.apply([
    initByFreq([.3,.5,.2]),
    pySubset( carray('i',[1,0,0,1,0]) ),
    dumper(alleleOnly=True, stage=PrePostMating)
   ])
#end


#file log/expressionhelp.log
help(pyEval.__init__)
#end

#turnOnDebug(DBG_ALL)
#turnOnDebug(DBG_SIMULATOR)

#file log/varPlotter.log
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
     varDim=3, win=10, title='subPop size', saveAs='log/simuDemo')
   ],
   end=30
)

#end
#PS convert log/simuDemo16.eps log/simuDemo16.png
#PS /bin/rm -f log/simuDemo*.eps

#file log/scipy.log
from simuUtil import *
from simuSciPy import *
pop = population(size=200, ploidy=2, loci=[3,4],
                 subPop=[50,50,100])
simu = simulator(pop, randomMating())
migr = migrator([[0,.2,.1],[.25,0,.1],[.1,.2,0]],
                mode=MigrByProbability)
pSize = stat(popSize=1, stage=PreMating)
show = pyEval(r"str(subPopSize)+'\n'")

pPlot = varPlotter(expr='subPopSize', win=10, update=5, title='subPop size',
    ytitle='size', legend=["rep"], saveAs="log/subpop.png")

simu.setGen(0)
simu.evolve([
   pSize,
   migr,
   show,
   pPlot,
   ],
   end=20
)

#end

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
  end=endGen,
  dryrun=False
)
#end

#file log/rng.log
print listAllRNG()
print rng().name()
setRNG("ranlux389")
print rng().name()
#end

#file log/rngrand.log
r=rng()
#help(RNG)
for n in range(1,10):
  print r.randBinomial(10, .7),
#end




#file log/extgenostru.log
pop=population(1, loci=[2,3,4])
print pop.numLoci(1)
print pop.locusDist(2)
dis = pop.arrLociDist()
print dis
dis[2] = 0.5
print pop.locusDist(2)
print pop.arrLociDist()
#end

#file log/extgenotype.log
InitByFreq(pop, [.2,.8])
Dump(pop, alleleOnly=1)
ind = pop.individual(0)
print ind.allele(1,1)
ind.setAllele(3,1,1)
Dump(pop, alleleOnly=1)
a = ind.arrAlleles()
print a
a = ind.arrAlleles(1)
print a
a = ind.arrAlleles(1,2)
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
simu.apply([ pyEval(stmts="pop=simu.population(rep)")])
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
   endl(rep=REP_LAST),
   #varPlotter(expr='LD[0][1]', title='Linkage disequilibrium',
   #  numRep = 4, ytitle='LD', saveAs='LD')
   ], preOps=[init],
   end=10
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
    pyEval(r"gen", rep=0),  # report progress
    endl(rep=REP_LAST)
    ],
    preOps=[init],
    end=endGen)


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
import Numeric
rates = Numeric.array([[0.]*nc]*nc)
for i in range(1,nc-1):
  rates[i,i+1]=0.05
  rates[i,i-1]=0.05
rates[0,1] = 0.1
rates[nc-1,nc-2] = 0.1

# print rates
print rates
migr = migrator( rate=rates, mode=MigrByProbability)

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
  preOps = [ initMigr ], end=10
  )

#end

# need reich.py
#PS /bin/cp -f ../examples/Reich2002/reich.py log/reich.py
