"""
This simulation tries to repeat results from

Balloux & Goudet ( 2002)
  Statistical properties of population differentiation
  estimators under stepwise mutation in a finite island model

It originally use easyPOP to get samples and use Fstat to calculate
Fst. Since Fst is builtin in simuPOP, we do not have to go through
the trouble.
"""

import simuOpt
# need to use long allele library.
# set optimized=True will gain around 50% speed up.
simuOpt.setOptions(optimized=True, longAllele=True)
from simuPOP  import *
from simuUtil import *
from simuRPy  import *

import exceptions

# parameter
nLoci = 12
ma = 999

# do simulation
# return simu for further analysis.
def simulate(subPop, migrRate, mutaRates, nRep=5, endGen=10, visual=[],
             savePop=0, name='simu'):
  """
  subPop: subpopulation sizes
  migrRate: migration rate
  mutaRates: mutation rates, an array of length nLoci
  nRep: number of replicates
  endGen: generations to run
  visual: visualizations (choose from Fst, AvgFst, alleleFreq, heteroFreq
  savePop: whether or not save population. If 0, do not save, otherwise,
    savepopulation every savePop generations.
  name: prefix for saved results and populations
  """
  
  numSP = len(subPop)
  if len(mutaRates) != nLoci:
    raise exceptions.ValueError("Please specify mutation rate for each locus.")
  
  ## initialization (uniform, all possible alleles)
  init =  initByFreq( [1./ma]*ma )
  
  ## mutation
  ##  different mutation rate at each locus
  ##  step-wise mutation model
  mutate = smmMutator( rates = mutaRates, atLoci = range(0, nLoci),
     maxAllele=ma)
  
  ## migration, r_ii will be automatically set correctly
  ## DO NOT USE MigrByProportion. (0.1 average migrants will be intepreted as migrate 0)
  migrate = migrator(
    rates = [[migrRate/(numSP-1)]* numSP] * numSP, 
    mode=MigrByProbability)
  
  ## visualizers
  v1,v2,v3,v4 = noneOp(),noneOp(),noneOp(),noneOp()
  ## visualize all Fst's
  if "Fst" in visual:
    v1 = varPlotter('Fst[:4]', byRep=True, title='Fst', numRep=nRep,
        varDim=4, update=100, win=2000, step=10, saveAs="Fst")
  ## plot average Fst (using all loci)
  if "AvgFst" in visual:
    v2 = varPlotter('AvgFst', byVal=1, title='Fst', numRep=nRep,
        varDim=1, update=10, win=2000, step=10, saveAs="avgFst")
  ## allele frequency
  if "alleleFreq" in visual:
    v3 = varPlotter('[alleleFreq[0][x] for x in range(1,11)]',
        byVal=1, title='expected heterozygosity', numRep=nRep,
        varDim=10, step=10, win=2000,update=100, saveAs="alleleFreq")
  ## heterozygote frequency
  if "heteroFreq" in visual:
    v3 = varPlotter('[heteroFreq[x][0] for x in range(0,4)]',
        byVal=1, title='expected heterozygosity', numRep=nRep,
        varDim=4, step=10, win=3000, update=100, saveAs="heter")
  
  ## save population for post-morten analysis
  if savePop > 0:
    saveAt = range(5000, endGen, savePop)
  else:
    saveAt = []

  ## statistics:
  ## calculate Fst and save result (all replicate in one file)
  ## calculate every 50 steps to save time
  stats = stat(Fst=range(0,nLoci), step=50)
  saveFst =  pyEval(r"'%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t' % (gen,Fst[0],Fst[1],Fst[2],Fst[3],AvgFst)",
    step = 100, output=">>"+name+"_Fst.txt", name="save Fst values to file")
  ## output \n at the end of last replicate
  saveFst1 = endl(rep=REP_LAST, step=100, output=">>"+name+"_Fst.txt")

  ## report progress (and average Fst
  report1 = pyEval(r"'%.6f\t' % AvgFst", step=100)
  report2 = endl(rep=REP_LAST, step=100)
    
  # create population and simulator
  simu = simulator(
    population(subPop=subPop, ploidy=2, loci=[1]*nLoci,
      maxAllele=ma),
    ## keep population size
    randomMating(newSubPopSize=subPop),
    rep=nRep)
  
  ## start simulation
  simu.evolve(
    preOps = [ init ],
    ops = [
      pyEval(r"'%d\n' % gen", step=100, rep=0, name="output generation"),
      mutate, migrate, stats,
      saveFst, saveFst1, 
      v1, v2, v3, v4,
      ## report running time
      ticToc(step=1000,rep=0),
    ],
    end = endGen,
    dryrun = False,
    saveAs = name,
    saveAt = saveAt,
    format = 'bin',
  )
  return simu

## First run, try to get some figures
if False:
  simu = simulate(subPop=[1000]*2, migrRate=0.0001,
    mutaRates=[.01]*nLoci, endGen=10000, nRep=4,
    visual=["Fst","AvgFst"], savePop=0, name='pop')

## verify Fst results, output to fstat format and compare
## The results are identical (on Fst)
if False:
  pop = LoadPopulation("pop1476_2.bin")
  Stat(pop, Fst=range(0,nLoci))
  print pop.dvars().Fst
  print pop.dvars().AvgFst
  SaveFstat(pop, "pop1476_2.dat")


## run 99 replicates, without visualization
if True:
  subPops = [ [1000]*2, [400]*5, [100]*20]
  Nm = [0.1, 1., 10.]
  mutaRates = [0.01, 0.001, 0.0001]
  index = 0
  rec = open("simulation.log", "w")
  rec.write("scenario\tsubPopSize\tmigrRate\tmutaRate\n")
  for subPop in subPops:
    for nm in Nm:
      for mutaRate in mutaRates:
        migrRate = nm/subPop[0]
        name = "sce" + str(index) + '_'
        print "Simulating scenario ", index, "N ", subPop[0], "m ", migrRate, "mu ", mutaRate
        rec.write("%s\t%d\t%.5f\t%.5f\n" % (name, subPop[0], migrRate, mutaRate))
        if index == 1 :  # starting from # scenario
          simu = simulate(subPop=subPop, migrRate=migrRate,
            mutaRates=[mutaRate]*nLoci, endGen=10000, nRep=100,
            visual=[], savePop=2500, name=name)  
        index += 1
  rec.close()


## run the first senario 40000 generations

if False:
## we are not sure if Fst reaches equilibrium at 1000 gen,
## case 0
  nRep=100
  nLoci = 12
  ma = 999
  subPop = [1000]*2
  nm = 1.
  mutaRate = 0.0001
  migrRate = nm/subPop[0]
  name = "sce0_ext_"
  simu = simulator(
    population(subPop=subPop, ploidy=2, loci=[1]*nLoci,
      maxAllele=ma),
    randomMating(newSubPopSize=subPop),
    rep=nRep)
  numSP = len(subPop)
  mutate = smmMutator( rates = [mutaRate]*nLoci, atLoci = range(0, nLoci),
     maxAllele=ma)
  migrate = migrator(
    rates = [[migrRate/(numSP-1)]* numSP] * numSP, 
    mode=MigrByProbability)
  stats = stat(Fst=range(0,nLoci), step=50)
  # really append to file
  saveFst =  pyEval(r"'%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t' % (gen,Fst[0],Fst[1],Fst[2],Fst[3],AvgFst)",
    step = 100, output=">>"+name+"_Fst.txt", name="save Fst values to file")
  saveFst1 = endl(rep=REP_LAST, step=100, output=">>"+name+"_Fst.txt")
  # start simulation
  simu.evolve(
    ops = [
      # report progress
      pyEval(r"'%d\n' % gen", step=100, rep=REP_LAST, name="output generation"),
      mutate, migrate, stats,
      saveFst, saveFst1,
      ticToc(step=1000,rep=0),
      ],
      end = 40000
    )

if False:
  nLoci = 12
  ma = 999
  subPop = [1000]*2
  nm = 1.
  mutaRate = 0.0001
  migrRate = nm/subPop[0]
  name = "sce5_10000_"
  names=[]
  ## load populaitons
  for i in range(0,100):
    names.append( name + str(i) + '.bin')
  
  simu = LoadSimulatorFromFiles( names,
       randomMating(newSubPopSize=subPop))
  
  numSP = len(subPop)
  mutate = smmMutator( rates = [mutaRate]*nLoci, atLoci = range(0, nLoci),
     maxAllele=ma)
  migrate = migrator(
    rates = [[migrRate/(numSP-1)]* numSP] * numSP, 
    mode=MigrByProbability)
  stats = stat(Fst=range(0,nLoci), step=50)
  # really append to file
  saveFst =  pyEval(r"'%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t' % (gen,Fst[0],Fst[1],Fst[2],Fst[3],AvgFst)",
    step = 100, output=">>"+name+"_Fst.txt", name="save Fst values to file")
  saveFst1 = endl(rep=REP_LAST, step=100, output=">>"+name+"_Fst.txt")
  simu.setGen(10001)
  # start simulation
  simu.evolve(
    ops = [
      # report progress
      pyEval(r"'%d\n' % gen", step=100, rep=REP_LAST, name="output generation"),
      mutate, migrate, stats,
      saveFst, saveFst1,
      ticToc(step=1000,rep=0),
      ],
      end = 20000
    )

