#!/usr/bin/env python
#
# Demonstrate the changes of allele frequencies among subpopulations due to migration
#

"""
This program demonstrates the effect of migration by rate which results in
changes of allele frequencies among subpopulations.
"""

import simuOpt, os, sys, types, time
from simuPOP import *

try:
    from simuPOP.plotter import VarPlotter
except:
    print "simuRPy import failed. Please check your rpy installation."
    print "allele frquencies in subpopulations will not be plotted"
    useRPy = False
else:
    useRPy = True

options = [
    {
     'name':'subPopSize',
     'default':5000,
     'label':'SubPopulation Size',
     'type':[int, long],
     'validator':simuOpt.valueGT(0),
     },
    {
     'name':'numOfSubPops',
     'default':5,
     'type':[int, long],
     'label':'Number of Subpopulations',
     'validator':simuOpt.valueGT(0)
     },
    {
     'name':'m',
     'default':0.05,
     'label':'Migration Rate',
     'type':[float],
     'validator':simuOpt.valueBetween(0., 1.),
     },
    {
     'name':'generations',
     'default':200,
     'label':'Generations to evolve',
     'description':'Length of evolution',
     'type':[int, long],
     'validator':simuOpt.valueGT(0)
     },
]

     
def simuMigration(subPopSize, numOfSubPops, m, generations):
    '''Simulate the change of allele frequencies among subpopulations as a result of migration.'''
    # diploid population, one chromosome with 1 locus
    # random mating with sex
    pop = Population(size=[subPopSize]*numOfSubPops, loci=[1], infoFields=['migrate_to'])
    # set initial allele frequencies to each subpopulation
    # Rule of initialization is to use a list of numbers, which have counts equal to number of
    # subpopulations and been set evenly distributed between 0 and 1.
    # Therefore, the average allele frequency named theoretical value among all subpopulations will be 0.5
    for i in range(numOfSubPops):
        initGenotype(pop, freq=[i*1./(numOfSubPops-1), 1 - i*1./(numOfSubPops-1)], subPops=[i])

    # check if plot
    if useRPy:
        plotter = VarPlotter('[0.5] + [subPop[i]["alleleFreq"][0][0] for i in range(%d)]' % numOfSubPops,
            ylim=[0, 1], update=generations-1, legend=['Theoretical'] + ['Allele freq at subpop %d' % x for x in range(numOfSubPops)],
	    xlab="generation", ylab="alleleFrequency", saveAs='migration.png')
    else:
        plotter = NoneOp()
        
    a=[]
    # set migration rate matrix a[]. Within each row i of matrix a[], all elements have the same value which is
    # m divided by number of subpopulations minus one, except the diagonal element, whose value is set to be 0.
    # This setting ensures that every individual in subpopulation i has probability m to become the Migrator,
    # and such a migrator.has equal possibility of entering any other subpopulation. 
    for i in range(numOfSubPops):
        b=[];
        for j in range(numOfSubPops):
            if j==i:
                b.append(0)
            else:
                b.append(m/(numOfSubPops-1))
        a.append(b)
    # if number of generations is smaller than 200, step will be set to 10 generations, otherwise it will be 20.
    if generations <= 200:
        s = 10
    else:
        s = 20
    pop.evolve(
        initOps = InitSex(),
        preOps = [
            Stat(alleleFreq=[0], vars=['alleleFreq_sp']), 
            PyEval(r'"Frequency at generation %d: %s\n" % (gen, ", ".join(["%.2f" % x for x in freq]))',
                stmts = 'freq = [subPop[x]["alleleFreq"][0][0] for x in range(%i)]' % numOfSubPops, step = s),
        ],
        matingScheme = RandomMating(),
        postOps = [
            Migrator(rate = a),
            plotter,
            ],
        gen = generations,
        )

if __name__ == '__main__':
    # get all parameters
    pars = simuOpt.Params(options, __doc__)
    if not pars.getParam():
        sys.exit(0)

    simuMigration(pars.subPopSize, pars.numOfSubPops, pars.m, pars.generations)
    #wait five seconds before exit
    if useRPy:
        print "Figure will be closed after five seconds."
        time.sleep(5)
