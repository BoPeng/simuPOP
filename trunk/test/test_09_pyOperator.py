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
        InitSex(self.pop)

    # define a function
    def myFunc(self, pop):
        Stat(pop, alleleFreq=[0])
        # allelic frequency was assigned to be 0.2
        assert abs(pop.dvars().alleleFreq[0][0] - 0.2) < 0.05
        return True

    def testSimpleFunc(self):
        'Testing Python operator'
        InitByFreq(self.pop, [.2, .3, .5])
        self.pop.evolve( postOps = pyOperator(self.myFunc),
            matingScheme=randomMating(),
            gen=20)

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

    def myFuncWithParam(self, pop, param):
        ' para is (allele, freq) pair '
        Stat(pop, alleleFreq=[0])
        assert abs(pop.dvars().alleleFreq[0][param[0]] - param[1]) < 0.05
        return True

    def testFuncWithParam(self):
        'Testing python operator with parameters'
        InitByFreq(self.pop, [.2, .8])
        self.pop.evolve( postOps=[
            pyOperator(func=self.myFuncWithParam, param=(0,.2)),
            pyOperator(func=self.myFuncWithParam, param=(1,.8)),
            ],
            matingScheme=randomMating(),
            gen=2
        )

    def myFuncAsTerminator(self, pop):
        if pop.dvars().gen == 3:
            return False
        else:
            return True

    def testTerminator(self):
        'Testing hybrid terminator'
        simu = simulator(self.pop)
        simu.evolve(initOps = [initSex()],
            postOps = pyOperator(self.myFuncAsTerminator),
            matingScheme=randomMating(),
            gen = 10 )
        assert simu.dvars(0).gen == 4

    def dynaMutator(self, pop, param):
        ''' this mutator mutate common loci with low mutation rate
        and rare loci with high mutation rate, as an attempt to
        bring allele frequency of these loci at an equal level.'''
        # unpack parameter
        (cutoff, mu1, mu2) = param;
        Stat(pop, alleleFreq=range( pop.totNumLoci() ) )
        if ModuleInfo()['alleleType'] == 'binary':
            for i in range( pop.totNumLoci() ):
                # 1-freq of wild type = total disease allele frequency
                if 1-pop.dvars().alleleFreq[i][0] < cutoff:
                    KamMutate(pop, k=2, rates=mu1, loci=[i])
                else:
                    KamMutate(pop, k=2, rates=mu2, loci=[i])
        else:
            for i in range( pop.totNumLoci() ):
                # 1-freq of wild type = total disease allele frequency
                if 1-pop.dvars().alleleFreq[i][1] < cutoff:
                    KamMutate(pop, k=2, rates=mu1, loci=[i])
                else:
                    KamMutate(pop, k=2, rates=mu2, loci=[i])
        return True

    def testDynaMutator(self):
        'Testing dynamic mutator (an example)'
        simu = simulator(self.pop)
        simu.evolve(
            initOps = [
                initByFreq( [.6, .4], loci=[0,2,4]),
                initByFreq( [.8, .2], loci=[1,3]) ],
            matingScheme=randomMating(),
            postOps = [
                pyOperator( func=self.dynaMutator, param=(.5, .1, 0) ),
                stat(alleleFreq=range(5)),
                terminateIf( 'alleleFreq[0][1] < 0.2' )
                ],
            gen = 30
        )
        # will not terminate when allele frequency get out of control
        self.assertEqual(simu.dvars(0).gen, 30)


if __name__ == '__main__':
    unittest.main()
