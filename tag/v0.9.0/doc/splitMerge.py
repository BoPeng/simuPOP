from simuOpt import *
setOptions(quiet=True)
#file splitAndMerge.out
from simuPOP import *
pop = population(1000, loci=[1], infoFields=['migrate_to'])
simu = simulator(pop, binomialSelection())
simu.evolve(
    ops=[
        splitSubPop(0, proportions=[0.2, 0.8], at = [3]),
        splitSubPop(1, proportions=[0.4, 0.6], at = [5]),
        mergeSubPops([0,2], at = [7]),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize'),
    ],
    gen = 10
)
#end
#file splitMigration.out
from simuPOP import *
pop = population(1000, loci=[1], infoFields=['migrate_to'])
simu = simulator(pop, binomialSelection())
simu.evolve(
    ops=[
        splitSubPop(0, proportions=[0.2, 0.3, 0.5], at = [3]),
        migrator(rate = [0.2], fromSubPop=[0], toSubPop=[1], 
            begin = 3, end = 4),
        migrator(rate = [
            [0, 0.2, 0.4],
            [0, 0,   0.1],
            [0.1, 0.1, 0]],
            begin = 4),
        stat(popSize=True),
        pyEval(r'"%s\n" % subPopSize'),
    ],
    gen = 10
)
#end
#file splitMigration2.out
from simuPOP import *
pop = population(1000, loci=[1], infoFields=['migrate_to'])
def popSize(gen, oldSize=[]):
    if gen < 3:
        return [1000]
    elif gen < 5:
        return [400, 500]
    else:
        return [300, 400, 600]

simu = simulator(pop, binomialSelection(newSubPopSizeFunc=popSize))
simu.evolve(
    ops=[
        splitSubPop(0, proportions=[0.3, 0.7], at = [3]),
        migrator(rate = [0.2], fromSubPop=[0], toSubPop=[1], 
            begin = 3, end = 4),
        splitSubPop(0, proportions=[0.3, 0.7], at = [5]),
        migrator(rate = [
            [0, 0.2, 0.4],
            [0, 0,   0.1],
            [0.1, 0.1, 0]],
            begin = 5),
        stat(popSize=True, stage=PreMating),
        pyEval(r'"From %s\t" % subPopSize', stage=PreMating),
        stat(popSize=True),
        pyEval(r'"to: %s\n" % subPopSize'),
    ],
    gen = 10
)
#end
#file splitMigration3.out
from simuPOP import *
pop = population(1000, loci=[1], infoFields=['migrate_to'])
def popSize(gen, oldSize=[]):
    return [x*2 for x in oldSize]

simu = simulator(pop, binomialSelection(newSubPopSizeFunc=popSize))
simu.evolve(
    ops=[
        splitSubPop(0, proportions=[0.3, 0.7], at = [3]),
        migrator(rate = [0.2], fromSubPop=[0], toSubPop=[1], 
            begin = 3, end = 4),
        splitSubPop(0, proportions=[0.3, 0.7], at = [5]),
        migrator(rate = [
            [0, 0.2, 0.4],
            [0, 0,   0.1],
            [0.1, 0.1, 0]],
            begin = 5),
        stat(popSize=True, stage=PrePostMating),
        pyEval(r'"From %s\t" % subPopSize', stage=PreMating),
        pyEval(r'"to: %s\n" % subPopSize'),
    ],
    gen = 10
)
#end



