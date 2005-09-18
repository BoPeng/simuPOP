# test plotting with Rpy package

from simuUtil import *
from simuRPy import *

# test aggregator

a = Aggregator(win=3, width=4)
a.push(0, [1,2,3,4])
print a
a.push(1, [2,2,3,4])
print a
a.push(1, [2,2,3,4])
print a
# the first column should be removed.
a.push(2, [2,3,4,5])
print a

print a.flatData()
#
a = Aggregator(win=3, width=4)
a.push(0,1,0)
print a
a.push(0,1,1)
print a
a.push(1,[1,2,3,4])
print a
a.push(1,23,3)
print a
a.push(2,22,2)
print a


# test plotter
p = VarPlotter_His()
p.plot(0,0,[2])
p.plot(1,0,[2])
p.plot(2,0,[3])
p.clear()

p = VarPlotter_His(varDim=3)
p.plot(0,0,[2,3,4])
p.plot(1,0,[4,5,6])
p.plot(2,0,[4,4,5])

p = VarPlotter_His(varDim=3, numRep=2, byRep=1)
p.plot(0,0,[2,3,4])
p.plot(1,0,[4,5,6])
p.plot(2,0,[4,4,5])
p.plot(0,1,[1,2,3])
p.plot(1,1,[2,3,4])
p.plot(2,1,[2,2,2])

# test the operator
simu = simulator(
  population(size=200, ploidy=2, loci=[3,4],
    subPop=[50,50,100]),
  randomMating(), rep=3)

migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
                mode=MigrByProbability)
stat = stat(popSize=1, stage=PreMating)

simu.setGen(0)
simu.evolve([
   migr, 
   stat,
   pyEval("list(subPopSize)"),
   tab(),
   endl(rep=REP_LAST),
   varPlotter('subPopSize', numRep=3, byRep=1, 
     varDim=3, win=10, title='subPop size')
   ],
   end=30
)



simu.setGen(0)
simu.evolve([
   migr, 
   stat,
   pyEval("list(subPopSize)"),
   tab(),
   endl(rep=REP_LAST),
   varPlotter('[x**2 for x in subPopSize]', numRep=3, byVal=1, 
     varDim=3, win=10, title='subPop size')
   ],
   end=30
)

# separate
simu.setGen(0)
simu.evolve([
   migr, 
   stat,
   pyEval("list(subPopSize)"),
   tab(),
   endl(rep=REP_LAST),
   varPlotter('[x**2 for x in subPopSize]', numRep=3, byVal=1, 
     varDim=3, win=10, title='subPop size', separate=1)
   ],
   end=30
)
# separate
simu.setGen(0)
simu.evolve([
   migr, 
   stat,
   pyEval("list(subPopSize)"),
   tab(),
   endl(rep=REP_LAST),
   varPlotter('subPopSize', numRep=3, byVal=1, 
     varDim=3, win=10, title='subPop size', separate=1)
   ],
   end=30
)
# each subplot has three lines
simu.setGen(0)
simu.evolve([
   migr, 
   stat,
   pyEval("list(subPopSize)"),
   tab(),
   endl(rep=REP_LAST),
   varPlotter('[x**2 for x in subPopSize]', numRep=3, byRep=1, 
     varDim=3, win=10, update=5, title='subPop size',
     saveAs='demo')
   ],
   end=30
)

#seprate lines in a subplot
simu.setGen(0)
simu.evolve([
   migr, 
   stat,
   pyEval("list(subPopSize)"),
   tab(),
   endl(rep=REP_LAST),
   varPlotter('subPopSize', numRep=3, byRep=1, ylim=[0,100],
     separate=True, varDim=3, win=10, update=5, title='subPop size', saveAs='demo')
   ],
   end=30
)

#

simu.setGen(0)
simu.evolve([
   migr, 
   stat,
   pyEval("list(subPopSize)"),
   tab(),
   endl(rep=REP_LAST),
   varPlotter('[x**2 for x in subPopSize]', history=False, numRep=3, byRep=1, 
     win=10, update=5, title='subPop size')
   ],
   end=30
)

# use separate panel.

simu.setGen(0)
simu.evolve([
   migr, 
   stat,
   pyEval("list(subPopSize)"),
   tab(),
   endl(rep=REP_LAST),
   varPlotter('[x**2 for x in subPopSize]', history=False, numRep=3,  
     win=10, update=5, separate=1, title='subPop size')
   ],
   end=30
)

# use image?

# test the operator
from simuRPy import *
from simuUtil import *
nr = 1
simu = simulator(
  population(size=200, ploidy=2, loci=[3,4],
    subPop=[50,50,100]),
  randomMating(), rep=nr)
migr = migrator([[0,.2,.1],[.25,0,.1],[.1,.2,0]],
                mode=MigrByProbability)
stat = stat(popSize=1, stage=PreMating)
simu.setGen(0)
simu.evolve([
   migr, 
   stat,
   pyEval("list(subPopSize)"),
   tab(),
   endl(rep=REP_LAST),
   varPlotter('subPopSize', history=True, numRep=nr, byRep=1,
     win=10, update=5, plotType='image', varDim=3, title='subPop size')
   ],
   end=10
)

# of course this is not the best image to draw.
# now, consider plotting penetrance?
# or affected status?

# check the spreading of disease
from simuRPy import *
from simuUtil import *
numRep=4
popSize=500
#turnOnDebug(DBG_MUTATOR)
#turnOnDebug(DBG_SIMULATOR)
simu = simulator(population(size=popSize, loci=[1]),
  randomMating(), rep=numRep)
simu.evolve(
  preOps = [ initByValue([1,1])],
  ops = [
    # penetrance, additve penetrance
    maPenetrance(locus=0, wildtype=[1], penetrance=[0,0.5,1],
       stage=PreMating, exposePenetrance=True),
    # count number of affected
    stat(numOfAffected=True),
    # introduce disease if no one is affected
    ifElse(cond='numOfAffected==0',
      ifOp=kamMutator(rate=0.01, maxAllele=2)),
    # expose affected status
    pyExec('pop.exposeAffectedness()', exposePop=True),
    # plot affected status
    varPlotter(expr='affected',plotType='image', byRep=1, update=50, 
      varDim=popSize, win=100, numRep=numRep, title='affectedness')
  ],
  end=100,
  dryrun=False
)

# testing histroy = false and plotType='image'
from simuRPy import *
from simuUtil import *
numRep=4
popSize=500
#turnOnDebug(DBG_MUTATOR)
#turnOnDebug(DBG_SIMULATOR)
simu = simulator(population(size=popSize, loci=[10]),
  randomMating(), rep=numRep)
simu.evolve(
  preOps = [ initByFreq([.2,.3,.5])],
  ops = [
    recombinator(rate=0.01),
    # count number of affected
    stat(LD=[ [x, y] for x in range(0,10) for y in range(0,10)]),
    # plot affected status
    varPlotter(expr='LD', plotType='image', numRep=4, byRep=1, update=5, 
      title='pairwise LD', history=False),
    # pause(rep=REP_LAST, step=5),
  ],
  end=20,
  dryrun=False
)

