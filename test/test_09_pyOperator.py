#!/usr/bin/env python
#
# Testing python operator
#
# Bo Peng (bpeng@rice.edu)
# 
# $LastChangedRevision$
# $LastChangedDate$
# 

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, time, exceptions
        
class TestPyOperator(unittest.TestCase):
    
    def setUp(self):
        self.pop = population(size=10000, ploidy=2, 
            loci=[2, 3])

    # define a function
    def myFunc(self, pop):
        Stat(pop, alleleFreq=[0])
        # allelic frequency was assigned to be 0.2
        assert abs(pop.dvars().alleleFreq[0][0] - 0.2) < 0.05
        return True
        
    def testSimpleFunc(self):
        'Testing a simple example'
        InitByFreq(self.pop, [.2, .3, .5])
        simu = simulator(self.pop, randomMating())
        simu.evolve( ops=[pyOperator(self.myFunc)], 
            end=20)
        
    def testCopyClone(self):
        'Testing copy of python operator'
        op = pyOperator(self.myFunc)
        op1 = op
        op2 = op.clone()
        InitByFreq(self.pop, [.2, .3, .5])
        # all copied version are working fine.
        op.apply(self.pop)
        op1.apply(self.pop)
        op2.apply(self.pop)
        
    def myFuncWithParam(self, pop, para):
        ' para is (allele, freq) pair '
        Stat(pop, alleleFreq=[0])
        assert abs(pop.dvars().alleleFreq[0][para[0]] - para[1]) < 0.05
        return True

    def testFuncWithParam(self):
        'Testing python operator with parameters'
        InitByFreq(self.pop, [.2, .8])
        simu = simulator(self.pop, randomMating())
        simu.evolve( ops=[
            pyOperator(func=self.myFuncWithParam, param=(0,.2)),
            pyOperator(func=self.myFuncWithParam, param=(1,.8)),
            ], 
            end=2
        )

    def myFuncAsTerminator(self, pop):
        # print pop.gen()
        if pop.gen() == 3:
            return False
        else:
            return True

    def testTerminator(self):
        'Testing a python terminator'
        simu = simulator(self.pop, randomMating())
        simu.evolve( ops=[ pyOperator(self.myFuncAsTerminator) ],
            end = 10 )
        assert simu.gen() == 4
        
    def dynaMutator(self, pop, param):
        ''' this mutator mutate common loci with low mutation rate
        and rare loci with high mutation rate, as an attempt to
        bring allele frequency of these loci at an equal level.'''
        # unpack parameter
        (cutoff, mu1, mu2) = param;
        Stat(pop, alleleFreq=range( pop.totNumLoci() ) )
        if alleleType() == 'binary':
            for i in range( pop.totNumLoci() ):
                # 1-freq of wild type = total disease allele frequency
                if 1-pop.dvars().alleleFreq[i][0] < cutoff:
                    KamMutate(pop, rate=mu1, loci=[i])
                else:
                    KamMutate(pop, rate=mu2, loci=[i])
        else:
            for i in range( pop.totNumLoci() ):
                # 1-freq of wild type = total disease allele frequency
                if 1-pop.dvars().alleleFreq[i][1] < cutoff:
                    KamMutate(pop, maxAllele=2, rate=mu1, loci=[i])
                else:
                    KamMutate(pop, maxAllele=2, rate=mu2, loci=[i])
        return True

    def testDynaMutator(self):
        'Testing dynamic mutator (an example)'
        simu = simulator(self.pop, randomMating())
        simu.evolve(
            preOps = [ 
                initByFreq( [.6, .4], loci=[0,2,4]),
                initByFreq( [.8, .2], loci=[1,3]) ],
            ops = [ 
                pyOperator( func=self.dynaMutator, param=(.5, .1, 0) ),
                stat(alleleFreq=range(5)),
                terminateIf( 'alleleFreq[0][1] < 0.2' )
                ],
            end = 30
        )                
        # will not terminate when allele frequency get out of control
        self.assertEqual( simu.gen(), 31)

    def testPyIndOperator(self):
        'Testing pyIndOperator'
        def indFunc(ind):
            ind.setInfo(ind.info('info1')+1, 'info1')
            return True
        def testFunc(ind):
            self.assertEqual(ind.info('info1'), 1)
            return True
        def indFunc1(ind, param):
            ind.setInfo(ind.info('info2')+param[1], 'info2')
            return True
        def testFunc1(ind):
            self.assertEqual(ind.info('info2'), 2)
            return True
        def indFunc2(ind, genotype):
            self.assertEqual(len(genotype), 4)
            if alleleType() == 'binary':
                self.assertEqual(genotype, (1, 1, 1, 1))
            else:
                self.assertEqual(genotype, (1, 1, 2, 2))
            ind.setInfo(sum(genotype), 'info1')
            return True
        def testFunc2(ind):
            if alleleType() == 'binary':
                self.assertEqual(ind.info('info1'), 4)
            else:
                self.assertEqual(ind.info('info1'), 6)
            return True
        def indFunc3(ind, genotype, param):
            self.assertEqual(len(genotype), 4)
            if alleleType() == 'binary':
                self.assertEqual(genotype, (1, 1, 1, 1))
            else:
                self.assertEqual(genotype, (1, 1, 2, 2))
            ind.setInfo(sum(genotype) + param[0], 'info1')
            return True
        def testFunc3(ind):
            if alleleType() == 'binary':
                self.assertEqual(ind.info('info1'), 2)
            else:
                self.assertEqual(ind.info('info1'), 4)
            return True
        def indFunc4(ind):
            return [1.5]
        def testFunc4(ind):
            self.assertEqual(ind.info('info1'), 1.5)
            return True            
        def indFunc5(ind, param):
            return param
        def testFunc5(ind):
            self.assertEqual(ind.info('info1'), 3)
            self.assertEqual(ind.info('info2'), 4)
            return True            
        def indFunc6(ind, genotype):
            self.assertEqual(len(genotype), 4)
            if alleleType() == 'binary':
                self.assertEqual(genotype, (1, 1, 1, 1))
            else:
                self.assertEqual(genotype, (1, 1, 2, 2))
            return (sum(genotype),)
        def testFunc6(ind):
            if alleleType() == 'binary':
                self.assertEqual(ind.info('info1'), 4)
            else:
                self.assertEqual(ind.info('info1'), 6)
            return True            
        def indFunc7(ind, genotype, param):
            self.assertEqual(len(genotype), 4)
            if alleleType() == 'binary':
                self.assertEqual(genotype, (1, 1, 1, 1))
            else:
                self.assertEqual(genotype, (1, 1, 2, 2))
            return (sum(genotype) + param[0],)
        def testFunc7(ind):
            if alleleType() == 'binary':
                self.assertEqual(ind.info('info2'), 6)
            else:
                self.assertEqual(ind.info('info2'), 8)
            return True            
        # use noMating, otherwise info will *not* be inherited
        simu = simulator(population(10, loci=[2,5], infoFields=['info1', 'info2']), 
            noMating(), rep=1)
        simu.step(
            preOps = [
                initByValue(value=[1], loci=[2]),
                initByValue(value=[2], loci=[4]),
                ],
            ops = [
                pyIndOperator(func=indFunc),
                pyIndOperator(func=testFunc),
                pyIndOperator(func=indFunc1, param=(3,2)),
                pyIndOperator(func=testFunc1),
                pyIndOperator(func=indFunc2, loci=[2,4]),
                pyIndOperator(func=testFunc2),
                pyIndOperator(func=indFunc3, loci=[2,4], param=(-2,)),
                pyIndOperator(func=testFunc3),
                pyIndOperator(func=indFunc4, infoFields=['info1']),
                pyIndOperator(func=testFunc4),
                pyIndOperator(func=indFunc5, param=(3,4), infoFields=['info1', 'info2']),
                pyIndOperator(func=testFunc5),
                pyIndOperator(func=indFunc6, loci=[2,4], infoFields=['info1']),
                pyIndOperator(func=testFunc6),
                pyIndOperator(func=indFunc7, loci=[2,4], param=(2,), infoFields=['info2']),
                pyIndOperator(func=testFunc7 ),
            ],
        )
            
        
if __name__ == '__main__':
    unittest.main()
