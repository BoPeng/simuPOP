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
        self.pop = Population(size=10000, ploidy=2,
            loci=[2, 3])
        initSex(self.pop)

    # define a function
    def myFunc(self, pop):
        stat(pop, alleleFreq=[0])
        # allelic frequency was assigned to be 0.2
        self.assertTrue(abs(pop.dvars().alleleFreq[0][0] - 0.2) < 0.05, 
            "abs(pop.dvars().alleleFreq[0][0] - 0.2) is supposed to be less than 0.05. This test may occasionally fail due to the randomness of outcome.")
        return True

    def testSimpleFunc(self):
        'Testing Python operator'
        initGenotype(self.pop, freq= [.2, .3, .5])
        self.pop.evolve( postOps = PyOperator(self.myFunc),
            matingScheme=RandomMating(),
            gen=20)

    def testCopyClone(self):
        'Testing copy of python operator'
        op = PyOperator(self.myFunc)
        op1 = op
        op2 = op.clone()
        initGenotype(self.pop, freq= [.2, .3, .5])
        # all copied version are working fine.
        op.apply(self.pop)
        op1.apply(self.pop)
        op2.apply(self.pop)

    def myFuncWithParam(self, pop, param):
        ' para is (allele, freq) pair '
        stat(pop, alleleFreq=[0])
        self.assertTrue(abs(pop.dvars().alleleFreq[0][param[0]] - param[1]) < 0.05, 
            "abs(pop.dvars().alleleFreq[0][param[0]] - param[1]) is supposed to be less than 0.05. This test may occasionally fail due to the randomness of outcome.")
        return True

    def testFuncWithParam(self):
        'Testing python operator with parameters'
        initGenotype(self.pop, freq= [.2, .8])
        self.pop.evolve( postOps=[
            PyOperator(func=self.myFuncWithParam, param=(0,.2)),
            PyOperator(func=self.myFuncWithParam, param=(1,.8)),
            ],
            matingScheme=RandomMating(),
            gen=2
        )

    def myFuncAsTerminator(self, pop):
        if pop.dvars().gen == 3:
            return False
        else:
            return True

    def testTerminator(self):
        'Testing hybrid terminator'
        simu = Simulator(self.pop)
        simu.evolve(initOps = [InitSex()],
            postOps = PyOperator(self.myFuncAsTerminator),
            matingScheme=RandomMating(),
            gen = 10 )
        self.assertEqual(simu.dvars(0).gen, 4)

    def dynaMutator(self, pop, param):
        ''' this mutator mutate common loci with low mutation rate
        and rare loci with high mutation rate, as an attempt to
        bring allele frequency of these loci at an equal level.'''
        # unpack parameter
        (cutoff, mu1, mu2) = param;
        stat(pop, alleleFreq=range( pop.totNumLoci() ) )
        if moduleInfo()['alleleType'] == 'binary':
            for i in range( pop.totNumLoci() ):
                # 1-freq of wild type = total disease allele frequency
                if 1-pop.dvars().alleleFreq[i][0] < cutoff:
                    kAlleleMutate(pop, k=2, rates=mu1, loci=[i])
                else:
                    kAlleleMutate(pop, k=2, rates=mu2, loci=[i])
        else:
            for i in range( pop.totNumLoci() ):
                # 1-freq of wild type = total disease allele frequency
                if 1-pop.dvars().alleleFreq[i][1] < cutoff:
                    kAlleleMutate(pop, k=2, rates=mu1, loci=[i])
                else:
                    kAlleleMutate(pop, k=2, rates=mu2, loci=[i])
        return True

    def testDynaMutator(self):
        'Testing dynamic mutator (an example)'
        simu = Simulator(self.pop)
        simu.evolve(
            initOps = [
                InitGenotype(freq= [.6, .4], loci=[0,2,4]),
                InitGenotype(freq= [.8, .2], loci=[1,3]) ],
            matingScheme=RandomMating(),
            postOps = [
                PyOperator( func=self.dynaMutator, param=(.5, .1, 0) ),
                Stat(alleleFreq=range(5)),
                TerminateIf( 'alleleFreq[0][1] < 0.2' )
                ],
            gen = 30
        )
        # will not terminate when allele frequency get out of control
        self.assertEqual(simu.dvars(0).gen, 30)


if __name__ == '__main__':
    unittest.main()
