#!/usr/bin/env python
#
# testing for simulator
# 
# Bo Peng (bpeng@rice.edu)
# 
# $LastChangedRevision$
# $LastChangedDate$
#
import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, os, sys

# record active generations in pop.dvars().hist
def genRecorder(pop):
    try:
        pop.dvars().hist.append(pop.gen())
    except:
        pop.dvars().hist = [pop.gen()]
    return True
    
def opRecorder(*args, **kwargs):
    return pyOperator(func=genRecorder, *args, **kwargs)
    
    
class TestSimulator(unittest.TestCase):

    def testPopulations(self):
        'Testing set/get populations'
        pop = population(size=2, loci=[1])
        simu = simulator(pop, randomMating(), rep=2)
        # pop is not affected if simu changes
        simu.population(0).individual(0).setAllele(1,0)
        self.assertEqual( pop.individual(0).allele(0), 0)
        # if we get a reference,...
        pop = simu.population(1)
        pop.individual(0).setAllele(1,0)
        self.assertEqual( simu.population(1).individual(0).allele(0), 1)
        # if we get a real copy out, they are independent
        pop1 = simu.getPopulation(1)
        simu.population(0).individual(1).setAllele(0,0)
        self.assertEqual( pop1.individual(0).allele(0), 1)
        # if we use simu.getPopulation(1, true), the population
        # within the simulator will be destroyed
        pop2 = simu.getPopulation(0, True)
        self.assertEqual(pop1, pop2)
        self.assertEqual(simu.population(0).popSize(), 0)
    
    def testProperties(self):
        'Testing simulator properties'
        pop = population(size=1, loci=[1])
        simu = simulator(pop, randomMating(), rep=3, stopIfOneRepStops=True)
        self.assertEqual( simu.numRep(), 3)
        self.assertEqual( simu.gen(), 0)
        self.assertEqual( simu.group(), (0,1,2) )
        simu.setGroup([1,1,2])
        self.assertEqual( simu.group(), (1,1,2) )
        simu.setGen(10)
        self.assertEqual( simu.gen(), 10)
        
    def testEvolve(self):
        'Testing function evolve and step'
        pop = population(size=1, loci=[1])
        simu = simulator(pop, randomMating(), rep=3)
        self.assertEqual( simu.gen(), 0)
        simu.step( ops=[] )
        # no terminator, no ending generation is specified
        self.assertRaises( exceptions.ValueError,
            simu.evolve, ops=[] )


    def testClone(self):
        'Testing cloning of simulator'
        pop = population(subPop=[10,4], loci=[2,5])
        simu = simulator(pop, randomMating(), rep = 3)
        simu.evolve(
            preOps = [initByFreq([0.3, .7])],
            ops = [stat(alleleFreq=range(pop.totNumLoci()))],
            end = 10
        )
        simu1 = simu.clone()
        for i in range(3):
            self.assertEqual(simu.population(i), simu1.population(i))
        self.assertEqual(simu1.gen(), simu.gen())
        # test if a cloned simulator can evolve again
        simu1.evolve(
            ops = [stat(alleleFreq=range(pop.totNumLoci()))],
            end = 20
        )
        
        

    def testGenoStru(self):
        'Testing genotypic structure related functions'
        # genetic structure can also be accessed from simulator
        if AlleleType() != 'binary':
            pop = population(size=100, ploidy=2, loci=[5, 7], 
                subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
                maxAllele=4, alleleNames=['_','A','C','T','G']) 
        else:
            pop = population(size=100, ploidy=2, loci=[5, 7], 
                subPop=[20, 80], lociPos=[ [2,3,4,5,6],[2,4,6,8,10,12,14]], 
                alleleNames=['1','2']) 
        #
        simu = simulator(pop, noMating() )
        #
        self.assertEqual(simu.ploidy(), 2)
        self.assertEqual(simu.ploidyName(), 'diploid')
        self.assertEqual(simu.numChrom(), 2)
        self.assertEqual(simu.numLoci(0), 5)
        self.assertEqual(simu.numLoci(1), 7)
        self.assertRaises(exceptions.IndexError, simu.numLoci, 2 )
        self.assertEqual(simu.locusPos(10), 12)
        self.assertRaises(exceptions.IndexError, simu.locusPos, 20 )
        self.assertRaises((exceptions.TypeError, exceptions.OverflowError), simu.locusPos, -1 )
        self.assertEqual(len(simu.arrLociPos()), 12)
        self.assertEqual(simu.arrLociPos().tolist(), [2,3,4,5,6,2,4,6,8,10,12,14])
        self.assertEqual(simu.chromBegin(0), 0)
        self.assertEqual(simu.chromBegin(1), 5)
        self.assertEqual(simu.chromEnd(0), 5)
        self.assertEqual(simu.chromEnd(1), 12)
        self.assertRaises(exceptions.IndexError, simu.chromBegin, 2 )
        self.assertRaises(exceptions.IndexError, simu.chromEnd, 2 )
        self.assertEqual(simu.absLocusIndex(1,5), 10)
        self.assertEqual(simu.locusPos(simu.absLocusIndex(1,2) ), 6)
        self.assertRaises(exceptions.IndexError, simu.absLocusIndex, 2, 5 )
        self.assertEqual(simu.chromLocusPair(10), (1,5) )
        self.assertRaises(exceptions.IndexError, simu.chromLocusPair, 50 )
        self.assertEqual(simu.totNumLoci(), 12)
        self.assertEqual(simu.genoSize(), simu.totNumLoci()*simu.ploidy() )
        if AlleleType() == 'binary':
            self.assertEqual(simu.alleleNames(), ('1','2') )
            self.assertEqual(simu.alleleName(0), '1')
            self.assertEqual(simu.alleleName(1), '2')
            # 5 is passed to be function as bool
            self.assertEqual(simu.alleleName(5), '2')
        else:
            self.assertEqual(simu.alleleName(0), '_')
            self.assertEqual(simu.alleleName(1), 'A')
            self.assertEqual(simu.alleleName(2), 'C')
            self.assertEqual(simu.alleleName(3), 'T')
            self.assertEqual(simu.alleleName(4), 'G')
            self.assertRaises(exceptions.IndexError, simu.alleleName, 5)
        # loci name, default, the name will be used by other programs
        # or file format, so we set it to be one based.
        self.assertEqual(simu.locusName(0), 'loc1-1')
        self.assertEqual(simu.locusName(1), 'loc1-2')
        self.assertEqual(simu.locusName(5), 'loc2-1')
        self.assertRaises(exceptions.IndexError, simu.locusName, 12)
        
    def testTerminator(self):
        'Testing terminator'
        simu = simulator(population(1), noMating() )
        simu.evolve( ops=[ terminateIf( 'gen==10' ) ] )
        # always point to the enxt gen
        self.assertEqual(simu.gen(), 11 )

    def testIntegrity(self):
        'Testing checking of simulator inegrity'
        simu = simulator(population(1), noMating())
        simu.population(0).addInfoField('something')
        # one can not change the genotype structure of the populations
        # in a simulator
        self.assertRaises(exceptions.ValueError, simu.evolve, ops=[])

    def testMultiRep(self):
        'Testing multi-replicates related functions'
        simu = simulator(population(1), noMating(), rep=3 )
        simu.evolve( 
            ops=[ 
                opRecorder(), 
                terminateIf( 'gen==10', rep=0 ),
                terminateIf( 'gen==15', rep=1 ),
                terminateIf( 'gen==20', rep=2 ) 
            ]
        )
        # by default, run until the last resplicate die
        self.assertEqual(simu.gen(), 21 )
        self.assertEqual(simu.population(0).dvars().hist, range(11))
        self.assertEqual(simu.population(1).dvars().hist, range(16))
        self.assertEqual(simu.population(2).dvars().hist, range(21))
        #
        # if    set stopIfOneRepStop
        simu = simulator(population(1), noMating(), rep=3, stopIfOneRepStops=True )
        simu.evolve( 
            ops=[ 
                opRecorder(), 
                terminateIf( 'gen==10', rep=0 ),
                terminateIf( 'gen==15', rep=1 ),
                terminateIf( 'gen==20', rep=2 ) 
            ]
        )
        # by default, run until the last resplicate die
        self.assertEqual(simu.gen(), 11 )
        self.assertEqual(simu.population(0).dvars().hist, range(11))
        self.assertEqual(simu.population(1).dvars().hist, range(11))
        self.assertEqual(simu.population(2).dvars().hist, range(11))
        #
        # if set applyOpToStoppedReps
        simu = simulator(population(1), noMating(), rep=3, applyOpToStoppedReps=True )
        simu.evolve( 
            ops=[ 
                opRecorder(), 
                terminateIf( 'gen==10', rep=0 ),
                terminateIf( 'gen==15', rep=1 ),
                terminateIf( 'gen==20', rep=2 ) 
            ]
        )
        # by default, run until the last resplicate die
        self.assertEqual(simu.gen(), 21 )
        self.assertEqual(simu.population(0).dvars().hist, range(21))
        self.assertEqual(simu.population(1).dvars().hist, range(21))
        self.assertEqual(simu.population(2).dvars().hist, range(21))


    def testAddInfoField(self):
        'Testing adding information fields to a simulator'
        simu = simulator(population(10, infoFields=['a']), randomMating(), rep=3)
        simu.addInfoField('b')
        self.assertEqual(simu.infoFields(), ('a', 'b'))
        simu.addInfoFields(['c', 'd'])
        simu.setAncestralDepth(2)
        simu.evolve(ops=[], end=10)
        self.assertEqual(simu.population(0).ancestralDepth(), 2)
        self.assertEqual(simu.infoFields(), ('a', 'b', 'c', 'd'))
        self.assertEqual(simu.population(0).infoFields(), ('a', 'b', 'c', 'd'))
        simu.population(0).addInfoField('l')
        self.assertRaises(exceptions.ValueError, simu.addInfoField, 'j')
        self.assertRaises(exceptions.ValueError, simu.evolve, ops=[])
        

        
if __name__ == '__main__':
    unittest.main()
