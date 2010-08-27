#!/usr/bin/env python
#
# testing operator behaviors for simupoop and general operators such as PyOperator and IfElse
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
        pop.dvars().hist.append(pop.dvars().gen)
    except:
        pop.dvars().hist = [pop.dvars().gen]
    return True

def opRecorder(*args, **kwargs):
    return PyOperator(func=genRecorder, *args, **kwargs)

class TestOperator(unittest.TestCase):

    def setUp(self):
        self.pop = Population(size=10000, ploidy=2,
            loci=[2, 3])
        initSex(self.pop)

    def testActiveGen(self):
        'Testing active generation specifications'
        def getActiveGens(endGen=20, *args, **kwargs):
            d = opRecorder(*args, **kwargs)
            simu = Simulator(Population())
            simu.evolve(postOps=d, gen=endGen)
            return simu.population(0).dvars().hist
        self.assertEqual(getActiveGens(begin=2, end=10),
            range(2,11))
        self.assertEqual(getActiveGens(begin=2, end=10, step=2),
            range(2,11,2))
        self.assertEqual(getActiveGens(begin=2, step=2),
            range(2,20,2))
        self.assertEqual(getActiveGens(step=2), range(0,19,2))
        self.assertEqual(getActiveGens(), range(0,20))
        self.assertEqual(getActiveGens(at=[2,5,9]), [2,5,9])
        self.assertEqual(getActiveGens(at=[2,5,-1]), [2,5,19])
        self.assertEqual(getActiveGens(begin=-10), range(10,20))
        # 20=-1, 16=-5
        self.assertEqual(getActiveGens(begin=-10, end=-5), range(10,16))
        #
        self.assertEqual(getActiveGens(begin=-10, step=2, end=-5), range(10,16,2))
        self.assertRaises( ValueError,
            getActiveGens, begin=-10, step=-3, end=-5 )

    def testReplicate(self):
        'Testing replicate related functions'
        simu = Simulator(Population(), rep=3)
        simu.evolve(
            postOps = opRecorder(reps=-1),
            matingScheme=CloneMating(),
            gen=10
        )
        try:
            simu.population(0).dvars().hist
        except AttributeError:
            pass
        try:
            simu.population(1).dvars().hist
        except AttributeError:
            pass
        self.assertEqual(simu.population(2).dvars().hist, range(10))

    def assertFileContent(self, file, text):
        f = open(file)
        t = f.read()
        f.close()
        self.assertEqual(t, text)

    def testOutput(self):
        'Testing output specifications'
        simu = Simulator( Population(), rep=5)
        simu.evolve(postOps = PyOutput("a", output=">a.txt"),
            matingScheme=CloneMating(),
            gen=10)
        # although everyone have written to this file,
        # only the last one will be kept
        self.assertFileContent("a.txt", 'a')
        os.remove('a.txt')
        #
        # you can ignore >
        for pop in simu.populations():
            pop.dvars().gen = 0
        simu.evolve(postOps = PyOutput("a", output="a.txt"),
            gen=10)
        # although everyone have written to this file,
        # only the last one will be kept
        self.assertFileContent("a.txt", 'a')
        os.remove('a.txt')
        #
        # >>
        for pop in simu.populations():
            pop.dvars().gen = 0
        simu.evolve(postOps = PyOutput("a", output=">>a.txt"),
            gen=10)
        # a is appended 5 rep * 11 generations
        self.assertFileContent("a.txt", 'a'*50)
        os.remove('a.txt')
        #
        # rep = ...
        for pop in simu.populations():
            pop.dvars().gen = 0
        simu.evolve(postOps = PyOutput("a", output=">>a.txt", reps=-1),
            gen=10)
        # a is appended 5 rep * 11 generations
        self.assertFileContent("a.txt", 'a'*10)
        os.remove('a.txt')

    def testOutputExpr(self):
        'Testing the usage of output expression'
        simu = Simulator( Population(), rep=5)
        # each replicate
        simu.evolve(postOps = PyOutput("a", output="!'rep%d.txt'%rep"),
            matingScheme=CloneMating(),
            gen=10)
        # although everyone have written to this file,
        # only the last one will be kept
        for i in range(5):
            self.assertFileContent("rep%d.txt"%i, 'a')
            os.remove('rep%d.txt'%i)
        #
        # you can ignore >
        for pop in simu.populations():
            pop.dvars().gen = 0
        simu.evolve(postOps = PyOutput("a", output="!'>rep%d.txt'%rep"),
            matingScheme=CloneMating(),
            gen=10)
        # although everyone have written to this file,
        # only the last one will be kept
        for i in range(5):
            self.assertFileContent("rep%d.txt"%i, 'a')
            os.remove('rep%d.txt'%i)
        #
        # >>
        for pop in simu.populations():
            pop.dvars().gen = 0
        simu.evolve(postOps = PyOutput("a", output="!'>>rep%d.txt'%rep"),
            matingScheme=CloneMating(),
            gen=10)
        # a is appended 1 rep * 11 generations
        for i in range(5):
            self.assertFileContent("rep%d.txt"%i, 'a'*10)
            os.remove('rep%d.txt'%i)
        # each generation?
        for pop in simu.populations():
            pop.dvars().gen = 0
        simu.evolve(postOps = PyOutput("a", output="!'>>gen%d.txt'%gen"),
            matingScheme=CloneMating(),
            gen=10)
        # a is appended 1 rep * 11 generations
        for i in range(10):
            self.assertFileContent("gen%d.txt"%i, 'a'*5)
            os.remove('gen%d.txt'%i)

    def testOutputFunc(self):
        '''Testing output to a function'''
        simu = Simulator(Population(), rep=5)
        def func1(msg):
            self.assertEqual(msg, 'func1')
        def func2(msg):
            self.assertEqual(msg, 'func2')
        # each replicate
        simu.evolve(
            matingScheme=CloneMating(),
            postOps = [
            PyOutput("func1", output=func1),
            PyOutput("func2", output=func2),
        ], gen=10)

    def testInfoEval(self):
        '''Testing operator InfoEval'''
        pop = Population(10, infoFields=['a', 'b'])
        infoEval(pop, expr='b', stmts='b=a+1', output='')
        # information field b is updated
        self.assertEqual(pop.indInfo('b'), tuple([1.]*10))
        #
        # use population variable
        pop.vars()['c'] = 5
        # usePopVars is needed
        infoEval(pop, 'c+4', output='')


    def testInfoExec(self):
        '''Testing operator InfoExec'''
        pop = Population(10, infoFields=['a', 'b'])
        infoExec(pop, 'b=a+1')
        self.assertEqual(pop.indInfo('b'), tuple([1]*10))
        infoExec(pop, 'a+=1')
        self.assertEqual(pop.indInfo('a'), tuple([1]*10))
        # this will not do anything because there is no c to be updated.
        infoExec(pop, 'c=a+b')
        #
        # use population variable
        pop.vars()['c'] = 5
        infoExec(pop, 'b=c+4')
        self.assertEqual(pop.indInfo('b'), tuple([9]*10))
        #
        # as an operator
        pop.evolve(
            initOps = [InfoExec('b=0')],
            matingScheme=CloneMating(),
            postOps = [
                InfoEval(r"'\t%.1f' % b", output=''),
                InfoExec('b+=1', output=''),
                PyOutput('\n', output=''),
            ],
            gen = 4
        )


    def testcloseOutput(self):
        '''Testing global function closeOutput'''
        pop = Population(100, loci=[2])
        dump(pop, output='a.pop')
        size = len(open('a.pop').read())
        self.assertRaises(RuntimeError, closeOutput, 'a.pop')
        dump(pop, output='>>a.pop')
        closeOutput('a.pop')
        self.assertEqual(len(open('a.pop').read()), size)
        self.assertRaises(RuntimeError, closeOutput, 'a.pop')
        self.assertRaises(RuntimeError, closeOutput, 'b.pop')
        dump(pop, output='>>a.pop')
        dump(pop, output='>>a.pop')
        self.assertEqual(len(open('a.pop').read()), size * 2)
        #
        dump(pop, output='>>>a.pop')
        dump(pop, output='>>>a.pop')
        self.assertEqual(len(open('a.pop').read()), size * 4)
        closeOutput('a.pop')
        self.assertRaises(RuntimeError, closeOutput, 'a.pop')
        os.remove('a.pop')


    # define a function
    def myFunc(self, pop):
        stat(pop, alleleFreq=[0])
        # allelic frequency was assigned to be 0.2
        self.assertTrue(abs(pop.dvars().alleleFreq[0][0] - 0.2) < 0.05, 
            "Expression abs(pop.dvars().alleleFreq[0][0] - 0.2) (test value %f) be less than 0.05. This test may occasionally fail due to the randomness of outcome." % (abs(pop.dvars().alleleFreq[0][0] - 0.2)))
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
            "Expression abs(pop.dvars().alleleFreq[0][param[0]] - param[1]) (test value %f) be less than 0.05. This test may occasionally fail due to the randomness of outcome." % (abs(pop.dvars().alleleFreq[0][param[0]] - param[1])))
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


    def testIfElseOperator(self):
        'Testing opeartor IfElse'
        pop = Population(1000, loci=[2])
        # now if we want to flip a lot of alleles whenever it reaches 0.2
        pop.evolve(
            # init to 1,2 (or 0, 1 in the binary case)
            initOps = [ InitSex(), InitGenotype(freq=[.5,.5]) ],
            matingScheme = RandomMating(),
            postOps = [
                # count number of allels at this locus
                Stat(alleleFreq=[0]),
                # inject 50% of allele 2 if this allele get low freq
                IfElse('alleleFreq[0][0]<0.2',
                    KAlleleMutator(rates=.6, k=2, loci=[0]) ),
                # the other way around?
                IfElse('alleleFreq[0][0]>0.8',
                    KAlleleMutator(rates=.6, k=2, loci=[0]) ),
                # terminate if
                TerminateIf('alleleFreq[0][0]<0.1 or alleleFreq[0][0]>0.9')
            ],
            gen = 1000
        )
        self.assertEqual(pop.dvars().gen, 1000)

    def testIfElseConstantExpre(self):
        'Testing operator IfElse with constant expression'
        def evolveGen(cond):
            pop = Population(100, loci=[2])
            return pop.evolve(
                postOps = IfElse(cond, NoneOp(), TerminateIf("True")),
                gen = 10)
        self.assertEqual(evolveGen(''), 1)
        self.assertEqual(evolveGen([]), 1)
        self.assertEqual(evolveGen(False), 1)
        self.assertEqual(evolveGen(True), 10)
        self.assertEqual(evolveGen('False'), 1)
        self.assertEqual(evolveGen('True'), 10)
        self.assertRaises(RuntimeError, evolveGen, 'nonsense')

    def testIfElseOperators(self):
        'Testing opeartor IfElse with multiple operators'
        simu = Simulator(Population(1000, loci=[2]))
        # now if we want to flip a lot of alleles whenever it reaches 0.2
        simu.evolve(
            # init to 1,2 (or 0, 1 in the binary case)
            initOps = [ InitSex(), InitGenotype(freq=[.5,.5]) ],
            matingScheme = RandomMating(),
            postOps = [
                # count number of allels at this locus
                Stat(alleleFreq=[0]),
                # inject 50% of allele 2 if this allele get low freq
                IfElse('alleleFreq[0][0]<0.2', ifOps=[
                    KAlleleMutator(rates=.6, k=2, loci=[0]),
                    NoneOp(),
                    ]),
                # the other way around?
                IfElse('alleleFreq[0][0]>0.8', ifOps = [
                    KAlleleMutator(rates=.6, k=2, loci=[0]),
                    NoneOp()
                    ]),
                # terminate if
                TerminateIf('alleleFreq[0][0]<0.1 or alleleFreq[0][0]>0.9')
            ],
            gen = 1000
        )
        self.assertEqual(simu.dvars(0).gen, 1000)


if __name__ == '__main__':
    unittest.main()



