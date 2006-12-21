from simuOpt import setOptions
setOptions(alleleType='long')
from simuPOP import *

# instantaneous population growth
def ins_exp(gen, oldSize=[]):
    if gen < burnin:
        return [initSize]
    else:
        return [finalSize]

def simulate(incSenario):
    simu = simulator(                                        # create a simulator
        population(subPop=incSenario(0), loci=[1,1],
            infoFields=['fitness']),                         # inital population
        randomMating(newSubPopSizeFunc=incSenario)           # random mating
    )

simulate(ins_exp)
