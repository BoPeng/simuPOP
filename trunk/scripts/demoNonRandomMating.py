#!/usr/bin/env python

# this script demonstrate the use of non-random mating of simuPOP
#
# Example 1:
#
# 

from simuPOP import *
import random

try:
    from rpy import *
    has_rpy = True
except:
    has_rpy = False

# try to load a more efficient version of parentsChooser
# parentsChooser_cpp is defined in demoNonRandomMating.i, and it must
# be compiled to be imported into this script.
try:
    from demoNonRandomMating import parentsChooser_cpp
    has_cpp = True
except:
    has_cpp = False


def locOfOffspring(loc):
    '''set offspring loc from parental locs'''
    # loc = (dad_x, dad_y, mom_x, mom_y)
    #
    # move to (dad_x + mom_x)/2 + N(0, 1), (dad_y + mom_y)/2 + N(0,1)
    new_x = (loc[0]+loc[2])/2. + random.normalvariate(0, 1)
    new_y = (loc[1]+loc[3])/2. + random.normalvariate(0, 1)
    # limit to region [-5, 5] x [-5, 5]
    new_x = min(new_x, 5)
    new_x = max(new_x, -5)
    new_y = min(new_y, 5)
    new_y = max(new_y, -5)
    #print '(%.2f, %.2f) + (%.2f, %.2f) ==> (%.2f, %.2f)' % \
    #    (loc[0], loc[2], loc[1], loc[3], new_x, new_y)
    return (new_x, new_y)


def plotInds(pop):
    '''plot the location of individuals. This requires R and rpy. '''
    if not has_rpy:
        return True
    r.postscript('loc_%d.eps' % pop.gen())
    r.plot(0, 0, xlim=[-5, 5], ylim=[-5, 5], type='n', 
        xlab='x', ylab='y',
        main='Locations of individuals at generation %d' % pop.gen())
    for ind in pop.individuals():
        r.points(ind.info(0), ind.info(1))
    r.dev_off()
    return True


def parentsChooser(pop, sp):
    '''Choose parents according to their locations. Because this is only a
       demonstration, performance is not under consideration.
    '''
    ################### The C++ version ############
    if has_cpp_chooser:
        pc = parentsChooser_cpp(
            [x.info('x') for x in pop.individuals() if x.sex() == Male],
            [x.info('y') for x in pop.individuals() if x.sex() == Male],
            [x.info('x') for x in pop.individuals() if x.sex() == Feale],
            [x.info('y') for x in pop.individuals() if x.sex() == Female],
            random.random())
        while True:
            print pc.chooseParents()
            yield pc.chooseParents()
    ##################### The Python version ###############
    males = [x for x in range(pop.popSize()) if pop.individual(x).sex() == Male]
    females = [x for x in range(pop.popSize()) if pop.individual(x).sex() == Female]
    if len(males) == 0 or len(females) == 0:
        print 'Lacking male or female. Existing'
        yield (None, None)
    while True:
        # randomly choose a male
        male = males[random.randint(0, len(males)-1)]
        # choose its closest female
        diff_x = [pop.individual(x).info(0) - pop.individual(male).info(0) for x in females]
        diff_y = [pop.individual(x).info(1) - pop.individual(male).info(1) for x in females]
        dist = [diff_x[i]**2 + diff_y[i]**2 for i in range(len(females))]
        female = females[dist.index(min(dist))]
        #print male, female
        yield (male, female)

    

#
pop = population(20, loci=[1], infoFields=['x', 'y'])
simu = simulator(pop,
    pyMating(pyParentsChooser(parentsChooser), mendelianOffspringGenerator())
)
simu.evolve(
    preOps = [initByFreq([0.5, 0.5])],
    ops = [
        pyEval(r'"%s\n" % gen'),
        pyTagger(func=locOfOffspring, infoFields=['x', 'y']),
        pyOperator(func=plotInds),
    ],
    end = 3
)
