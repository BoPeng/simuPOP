# test for plotting using SciPy


from simuUtil import *
from simuSciPy import *

pop = population(size=200, ploidy=2, loci=[3,4],
                 subPop=[50,50,100])
simu = simulator(pop, randomMating())
migr = migrator([[0,.2,.1],[.25,0,.1],[.1,.2,0]],
                mode=MigrByProbability)
stat = stat(popSize=1, stage=PreMating)
show = pyEval(r"str(list(subPopSize))+'\n'")
pInit = """
p = varPlotter(win=10, update=5, title='subPop size',
    ytitle='size', legend=["rep"])
"""
pPlot = "p.plot(gen, subPopSize)"
pSave = "p.save('subpop.png','png')"

simu.setGen(0)
simu.evolve([
   stat,
   migr,
   show,
   pyEval(pPlot, output='')
   ],
   preOps=[ pyEval(stmts=pInit, output='')],
   postOps=[ pyEval(stmts=pSave, output='')],
   end=20
)
