#!/usr/bin/env python
#
# testing operator behaviors for simupoop and general operators such as PyOperator and IfElse
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision$
# $LastChangedDate$
#

import unittest, os, sys
from simuOpt import setOptions
setOptions(quiet=True)
new_argv = []
for arg in sys.argv:
    if arg in ['short', 'long', 'binary', 'mutant', 'lineage']:
        setOptions(alleleType = arg)
    elif arg.startswith('-j'):
        setOptions(numThreads = int(arg[2:]))
    else:
        new_argv.append(arg) 

sys.argv=new_argv

from simuPOP import *

# record active generations in pop.dvars().hist
def genRecorder(pop):
    try:
        pop.dvars().hist.append(pop.dvars().gen)
    except:
        pop.dvars().hist = [pop.dvars().gen]
    return True

def opRecorder(*args, **kwargs):
    return PyOperator(func=genRecorder, *args, **kwargs)

class func:
    def __init__(self):
        pass

    def __call__(self, *args):
        return sum(args)

class func1:
    def __init__(self):
        pass

    def __call__(self, x1, x2):
        return sum([x1, x2])

class OutputCollector:
    def __init__(self):
        self.messages = []

    def func1(self, msg):
        self.messages.append(msg)

    def func2(self, msg):
        self.messages.append(msg)

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
            list(range(2,11)))
        self.assertEqual(getActiveGens(begin=2, end=10, step=2),
            list(range(2,11,2)))
        self.assertEqual(getActiveGens(begin=2, step=2),
            list(range(2,20,2)))
        self.assertEqual(getActiveGens(step=2), list(range(0,19,2)))
        self.assertEqual(getActiveGens(), list(range(0,20)))
        self.assertEqual(getActiveGens(at=[2,5,9]), [2,5,9])
        self.assertEqual(getActiveGens(at=[2,5,-1]), [2,5,19])
        self.assertEqual(getActiveGens(begin=-10), list(range(10,20)))
        # 20=-1, 16=-5
        self.assertEqual(getActiveGens(begin=-10, end=-5), list(range(10,16)))
        #
        self.assertEqual(getActiveGens(begin=-10, step=2, end=-5), list(range(10,16,2)))
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
        self.assertEqual(simu.population(2).dvars().hist, list(range(10)))

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

    def testWithArgs(self):
        '''Testing class WithArgs'''
        def blah(*args):
            return sum(args)
        pop = Population(100, infoFields=['trait'] + ['x%d' % x for x in range(20)])
        pyQuanTrait(pop, WithArgs(blah, ['x%d' % x for x in range(20)]),
            infoFields='trait')
        pyQuanTrait(pop, WithArgs(func(), ['x%d' % x for x in range(20)]),
            infoFields='trait')
        pop = Population(100, infoFields=['trait'] + ['x%d' % x for x in range(20)])
        pyQuanTrait(pop, func1(), infoFields='trait')

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


    def testOutputMemberFunc(self):
        '''Testing output to a function'''

        simu = Simulator(Population(), rep=5)
        o = OutputCollector()
        # each replicate
        simu.evolve(
            matingScheme=CloneMating(),
            postOps = [
            PyOutput("func1", output=o.func1),
            PyOutput("func2", output=o.func2),
        ], gen=10)
        self.assertEqual(o.messages, ['func1', 'func2'] * 50)

    def testOutputFileHandler(self):
        '''Testing output to a file handler'''
        simu = Simulator(Population(), rep=5)
        with open('test1.txt', 'w') as out:
            # each replicate
            simu.evolve(
                matingScheme=CloneMating(),
                postOps = [
                PyOutput("func1", output=out.write),
            ], gen=10)
        with open('test1.txt') as test1txt:
            self.assertEqual(test1txt.read(), 'func1'*50)
        #
        simu = Simulator(Population(), rep=5)
        with open('test1.txt', 'w') as out:
            # each replicate
            simu.evolve(
                matingScheme=CloneMating(),
                postOps = [
                PyOutput("func2", output=out),
            ], gen=10)
        with open('test1.txt') as test1txt:
            self.assertEqual(test1txt.read(), 'func2'*50)
        # runtime error only raise in python 3
        simu = Simulator(Population(), rep=5)
        with open('test1.txt', 'wb') as out:
            # each replicate
            self.assertRaises(RuntimeError, simu.evolve,
                matingScheme=CloneMating(),
                postOps = [
                PyOutput("func2", output=out),
            ], gen=10)
        #
        simu = Simulator(Population(), rep=5)
        with open('test1.txt', 'wb') as out:
            # each replicate
            simu.evolve(
                matingScheme=CloneMating(),
                postOps = [
                PyOutput('func2', output=WithMode(out, 'b')),
            ], gen=10)
        with open('test1.txt') as test1txt:
            self.assertEqual(test1txt.read(), 'func2'*50)

    def testInfoEval(self):
        '''Testing operator InfoEval'''

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

    def testIdTagger(self):
        '''Testing operator IdTagger'''
        pop = Population(10000, infoFields='ind_id', ancGen=-1)
        pop.evolve(
            initOps=[
                InitSex(),
                IdTagger(),
            ],
            matingScheme=RandomMating(
                ops=IdTagger()
            ),
            gen=10
        )
        IDs = [ind.ind_id for ind in pop.allIndividuals()]
        self.assertEqual(len(set(IDs)), 11*10000)
    
    def testTicToc(self):
        '''Testing operator TicToc'''
        pop = Population(10000)
        gen = pop.evolve(
            initOps=InitSex(),
            preOps=TicToc(stopAfter=0.5, output=''),
            matingScheme=RandomMating(),
            gen=10000
        )
        self.assertTrue(gen < 10000)
        gen = pop.evolve(
            initOps=InitSex(),
            matingScheme=RandomMating(
                ops=[
                    MendelianGenoTransmitter(),
                    TicToc(stopAfter=1, output='')
                ],
            ),
            gen=10000
        )
        self.assertTrue(gen < 20000)
    

    def testPyEval(self):
        '''Testing operator PyEval'''
        pop = Population([100,2000])
        stat(pop, popSize=True)
        self.assertEqual(pyEval(pop, 'popSize'), '2100')
        stat(pop, popSize=True, vars=['popSize_sp'])
        self.assertEqual(pyEval(pop, 'popSize', subPops=ALL_AVAIL), '1002000')
        self.assertEqual(pyEval(pop, r'"%d\t" % popSize', subPops=ALL_AVAIL),
            '100\t2000\t')
        pyExec(pop, 'p=popSize+200')
        self.assertRaises(RuntimeError, pyEval, pop, 'p', subPops=ALL_AVAIL)
        self.assertRaises(ValueError, PyEval, expr=" 1.2")        
        self.assertRaises(ValueError, PyEval, expr="a", stmts=' a=1')
        self.assertRaises(ValueError, PyEval, expr="\ta", stmts=' a=1')
        pyExec(pop, 'q=popSize+200', subPops=ALL_AVAIL)
        self.assertEqual(pyEval(pop, r'"%d\t" % q', subPops=ALL_AVAIL),
            '300\t2200\t')

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


    def testDumper(self):
        '''Testing operator Dumper'''
        pop = Population(100, loci=10, infoFields=('a', 'b'))
        dump(pop, genotype=False, output='a.dump')
        dump(pop, structure = False, output='a.dump')
        dump(pop, width=2, output='a.dump')
        dump(pop, infoFields='a', output='a.dump')
        dump(pop, infoFields=('b', 'a', 'b'), output='a.dump')

    def testcloseOutput(self):
        '''Testing global function closeOutput'''
        pop = Population(100, loci=[2])
        dump(pop, output='a.pop')
        with open('a.pop') as apop:
            size = len(apop.read())
        self.assertRaises(RuntimeError, closeOutput, 'a.pop')
        dump(pop, output='>>a.pop')
        closeOutput('a.pop')
        with open('a.pop') as apop:
            self.assertEqual(len(apop.read()), size)
        self.assertRaises(RuntimeError, closeOutput, 'a.pop')
        self.assertRaises(RuntimeError, closeOutput, 'b.pop')
        dump(pop, output='>>a.pop')
        dump(pop, output='>>a.pop')
        with open('a.pop') as apop:
            self.assertEqual(len(apop.read()), size * 2)
        #
        dump(pop, output='>>>a.pop')
        dump(pop, output='>>>a.pop')
        with open('a.pop') as apop:
            self.assertEqual(len(apop.read()), size * 4)
        closeOutput('a.pop')
        self.assertRaises(RuntimeError, closeOutput, 'a.pop')
        os.remove('a.pop')


    # define a function
    def myFunc(self, pop):
        stat(pop, alleleFreq=[0])
        # allelic frequency was assigned to be 0.2
        # self.assertTrue(abs(pop.dvars().alleleFreq[0][0] - 0.2) < 0.05, 
        #     "Expression abs(pop.dvars().alleleFreq[0][0] - 0.2) (test value %f) be less than 0.05. This test may occasionally fail due to the randomness of outcome." % (abs(pop.dvars().alleleFreq[0][0] - 0.2)))
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
        # self.assertTrue(abs(pop.dvars().alleleFreq[0][param[0]] - param[1]) < 0.05, 
        #     "Expression abs(pop.dvars().alleleFreq[0][param[0]] - param[1]) (test value %f) be less than 0.05. This test may occasionally fail due to the randomness of outcome." % (abs(pop.dvars().alleleFreq[0][param[0]] - param[1])))
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
        stat(pop, alleleFreq=list(range( pop.totNumLoci())) )
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
                Stat(alleleFreq=list(range(5))),
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

    def testSubPopsOfDuringMatingOperator(self):
        'Testing subPops parameter of during mating operators'
        pop = Population([100]*2, infoFields='a')
        pop.evolve(
            initOps=InitSex(),
            matingScheme=RandomMating(ops=
                [MendelianGenoTransmitter()] + 
                 [InfoExec('a=%d' % i, subPops=i) for i in range(10)]),
            gen=1
        )
        self.assertEqual(pop.indInfo('a'), tuple([0.0]*100 + [1.0]*100))
        # virtual subpopulation
        pop.setVirtualSplitter(SexSplitter())
        pop.evolve(
            initOps=InfoExec('a=0'),
            matingScheme=RandomMating(ops=
                [MendelianGenoTransmitter(),
                 InfoExec('a=1', subPops=[(ALL_AVAIL,0)]),
                 InfoExec('a=2', subPops=[(ALL_AVAIL,1)])]),
            gen=1
        )
        for ind in pop.individuals():
            if ind.sex() == MALE:
                self.assertEqual(ind.a, 1.0)
            else:
                self.assertEqual(ind.a, 2.0)

    def testRevertFixedSites(self):
        'Testing operator TestRevertFixedSites'
        #
        # for diploid
        pop = Population([1000]*2, loci=10)
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop)
        stat(pop, alleleFreq=ALL_AVAIL)
        for loc in range(10):
            if loc in [2, 4]:
                self.assertEqual(pop.dvars().alleleFreq[loc][0], 1)
            else:
                self.assertNotEqual(pop.dvars().alleleFreq[loc][0], 0)
        # for one subpopulation
        pop = Population([1000]*2, loci=10)
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop, subPops=[1])
        stat(pop, alleleFreq=ALL_AVAIL)
        for loc in range(10):
            if loc in [2, 4]:
                self.assertEqual(pop.dvars().alleleFreq[loc][0], 0.5)
            else:
                self.assertNotEqual(pop.dvars().alleleFreq[loc][0], 1)
        # parameter loci
        pop = Population([1000]*2, loci=10)
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop, loci=range(4))
        stat(pop, alleleFreq=ALL_AVAIL)
        for loc in range(10):
            if loc == 2:
                self.assertEqual(pop.dvars().alleleFreq[loc][0], 1)
            else:
                self.assertNotEqual(pop.dvars().alleleFreq[loc][0], 1)
        # for one virtual subpopulation
        #
        pop = Population([1000]*2, loci=10)
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop, subPops=[(0,0), (1,1)])
        stat(pop, alleleFreq=ALL_AVAIL, vars='alleleFreq_sp')
        for loc in range(10):
            if loc == 2:
                self.assertEqual(pop.dvars(0).alleleFreq[loc][0], 0.5)
                self.assertEqual(pop.dvars(1).alleleFreq[loc][0], 0.5)
            else:
                self.assertNotEqual(pop.dvars(0).alleleFreq[loc][0], 0)
                self.assertNotEqual(pop.dvars(1).alleleFreq[loc][0], 0)
        # 
        # haploid, sex chromosome etc
        #
        pop = Population([1000]*2, ploidy=1, loci=10)
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop)
        stat(pop, alleleFreq=ALL_AVAIL)
        for loc in range(10):
            if loc in [2, 4]:
                self.assertEqual(pop.dvars().alleleFreq[loc][0], 1)
            else:
                self.assertNotEqual(pop.dvars().alleleFreq[loc][0], 0)
        # for one subpopulation
        pop = Population([1000]*2,ploidy=1,  loci=10)
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop, subPops=[1])
        stat(pop, alleleFreq=ALL_AVAIL)
        for loc in range(10):
            if loc in [2, 4]:
                self.assertEqual(pop.dvars().alleleFreq[loc][0], 0.5)
            else:
                self.assertNotEqual(pop.dvars().alleleFreq[loc][0], 1)
        # parameter loci
        pop = Population([1000]*2, ploidy=1, loci=10)
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop, loci=range(4))
        stat(pop, alleleFreq=ALL_AVAIL)
        for loc in range(10):
            if loc == 2:
                self.assertEqual(pop.dvars().alleleFreq[loc][0], 1)
            else:
                self.assertNotEqual(pop.dvars().alleleFreq[loc][0], 1)
        # for one virtual subpopulation
        #
        pop = Population([1000]*2, ploidy=1, loci=10)
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop, subPops=[(0,0), (1,1)])
        stat(pop, alleleFreq=ALL_AVAIL, vars='alleleFreq_sp')
        for loc in range(10):
            if loc == 2:
                self.assertEqual(pop.dvars(0).alleleFreq[loc][0], 0.5)
                self.assertEqual(pop.dvars(1).alleleFreq[loc][0], 0.5)
            else:
                self.assertNotEqual(pop.dvars(0).alleleFreq[loc][0], 0)
                self.assertNotEqual(pop.dvars(1).alleleFreq[loc][0], 0)
        # 
        # sex chromosome
        #
        pop = Population([1000]*2, loci=[3, 7], chromTypes=[CHROMOSOME_X, CHROMOSOME_Y])
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop)
        stat(pop, alleleFreq=ALL_AVAIL)
        for loc in range(10):
            if loc in [2, 4]:
                self.assertEqual(pop.dvars().alleleFreq[loc][0], 1)
            else:
                self.assertNotEqual(pop.dvars().alleleFreq[loc][0], 0)
        # for one subpopulation
        pop = Population([1000]*2, loci=[3, 7], chromTypes=[CHROMOSOME_X, CHROMOSOME_Y])
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop, subPops=[1])
        stat(pop, alleleFreq=ALL_AVAIL)
        for loc in range(10):
            if loc in [2, 4]:
                self.assertEqual(pop.dvars().alleleFreq[loc][0], 0.5)
            else:
                self.assertNotEqual(pop.dvars().alleleFreq[loc][0], 1)
        # parameter loci
        pop = Population([1000]*2, loci=[3, 7], chromTypes=[CHROMOSOME_X, CHROMOSOME_Y])
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop, loci=range(4))
        stat(pop, alleleFreq=ALL_AVAIL)
        for loc in range(10):
            if loc == 2:
                self.assertEqual(pop.dvars().alleleFreq[loc][0], 1)
            else:
                self.assertNotEqual(pop.dvars().alleleFreq[loc][0], 1)
        # for one virtual subpopulation
        #
        pop = Population([1000]*2, loci=[3, 7], chromTypes=[CHROMOSOME_X, CHROMOSOME_Y])
        initSex(pop, sex=[MALE, FEMALE])
        pop.setVirtualSplitter(SexSplitter())
        initGenotype(pop, freq=[0.4, 0.6])
        if moduleInfo()['alleleType'] != 'binary':
            initGenotype(pop, loci=[2, 4], freq=[0, 1])
        else:
            initGenotype(pop, loci=[2, 4], freq=[0, 0, 1])
        revertFixedSites(pop, subPops=[(0,0), (1,0)])
        stat(pop, alleleFreq=ALL_AVAIL, subPops=[(0,0), (1,0)], vars='alleleFreq_sp')
        for loc in range(10):
            if loc in [2, 4]:
                self.assertEqual(pop.dvars((0,0)).alleleFreq[loc][0], 1)
                self.assertEqual(pop.dvars((1,0)).alleleFreq[loc][0], 1)
            else:
                self.assertNotEqual(pop.dvars((0,0)).alleleFreq[loc][0], 0)
                self.assertNotEqual(pop.dvars((1,0)).alleleFreq[loc][0], 0)

    def locateLoci(self):
        return [1,2]

    def locateEmptyLoci(self, pop):
        stat(pop, alleleFreq=ALL_AVAIL)
        return [x for x in range(pop.totNumLoci()) if pop.dvars().alleleFreq[x][1] == 0]   

    def testDynamicLoci(self):
        pop = Population(100, loci=10)
        initGenotype(pop, freq=[0.2, 0.8], loci=self.locateLoci)
        stat(pop, alleleFreq=ALL_AVAIL)
        self.assertGreater(pop.dvars().alleleFreq[1][1], 0)
        self.assertGreater(pop.dvars().alleleFreq[2][1], 0)
        self.assertEqual(pop.dvars().alleleFreq[3][1], 0)
        #
        initGenotype(pop, freq=[0.2, 0.8], loci=self.locateEmptyLoci)
        stat(pop, alleleFreq=ALL_AVAIL)
        for i in range(10):
            self.assertGreater(pop.dvars().alleleFreq[i][1], 0)
        #
        pop = Population(100, loci=10)
        initGenotype(pop, freq=[0.2, 0.8], loci=lociFunc(1))
        stat(pop, alleleFreq=ALL_AVAIL)
        for i in range(10):
            if i == 1:
                self.assertGreater(pop.dvars().alleleFreq[i][1], 0)
            else:
                self.assertEqual(pop.dvars().alleleFreq[i][1], 0)


    def testMemberFunc(self):
        for i in range(100):
            pop = Population(10)
            pop.evolve(
                matingScheme=RandomSelection(),
                postOps=PyOperator(func=MemberFunc().my_func, param=self),
                gen=10
            )
        self.assertEqual(MemberFunc.instance_count, 0)
        #
        me = MemberFunc()
        for i in range(100):
            pop = Population(10)
            pop.evolve(
                matingScheme=RandomSelection(),
                postOps=PyOperator(func=me.my_func, param=self),
                gen=10
            )
        self.assertEqual(MemberFunc.instance_count, 1)


    def testFunctor(self):
        Functor.instance_count = 0
        for i in range(100):
            pop = Population(10)
            pop.evolve(
                matingScheme=RandomSelection(),
                postOps=PyOperator(func=Functor(), param=self),
                gen=10
            )
        self.assertLessEqual(Functor.instance_count, 1)
        

    def testNonListOpList(self):
        '''Test if simuPOP accepts list of operators such as tuple, set, and DictKey'''
        pop = Population(10, loci=2)
        p1 = InitSex()
        p2 = InitGenotype(freq=[0.6, 0.4])
        for ops in [
            # single
            p1,
            # sequence
            [p1, p2],
            (p1, p2),
            # iterator
            (x for x in [p1, p1, p2]),
            # set
            {p1, p2},
            {'1': p1, '2': p2} .values(),
        ]:
            pop.evolve( initOps=ops, gen=1)

    def testWrapperOp(self):
        # test memory leak of Function or class that involves self
        # simuPOP automatically clears memory but the first instance
        # is not cleared (due to the fact that we do not know the 
        # life of the object at the C++ level.
        for i in range(100):
            pop = Population(10)
            pop.evolve(
                matingScheme=RandomSelection(),
                postOps=WrapperOpFunc(param=self),
                gen=10
            )
        self.assertLessEqual(WrapperOpFunc.instance_count, 1)
        #
        for i in range(100):
            pop = Population(10)
            pop.evolve(
                matingScheme=RandomSelection(),
                postOps=WrapperOpCall(param=self),
                gen=10
            )
        self.assertLessEqual(WrapperOpCall.instance_count, 1)
        #
        for i in range(100):
            pop = Population(10)
            pop.evolve(
                matingScheme=RandomSelection(),
                postOps=WrapperOpSelf(param=self),
                gen=10
            )
        self.assertLessEqual(WrapperOpSelf.instance_count, 1)
        #

class lociFunc:
    def __init__(self, a):
        self.a = a

    def __call__(self, pop):
        return self.a

class WrapperOpFunc(PyOperator):
    instance_count = 0
    def __init__(self, param):
        self.my_member = 23434
        WrapperOpFunc.instance_count += 1
        PyOperator.__init__(self, func=self.my_func, param=param)

    def my_func(self, pop, param):
        # can access to self from within PyOperator
        param.assertEqual(self.my_member, 23434)
        return True

    def __del__(self):
        WrapperOpFunc.instance_count -= 1

class WrapperOpCall(PyOperator):
    instance_count = 0
    def __init__(self, param):
        self.my_member = 23434
        WrapperOpCall.instance_count += 1
        PyOperator.__init__(self, func=self.__call__, param=param)

    def __call__(self, pop, param):
        # can access to self from within PyOperator
        param.assertEqual(self.my_member, 23434)
        return True

    def __del__(self):
        WrapperOpCall.instance_count -= 1

class WrapperOpSelf(PyOperator):
    instance_count = 0
    def __init__(self, param):
        self.my_member = 23434
        WrapperOpSelf.instance_count += 1
        PyOperator.__init__(self, func=self, param=param)

    def __call__(self, pop, param):
        # can access to self from within PyOperator
        param.assertEqual(self.my_member, 23434)
        return True

    def __del__(self):
        WrapperOpSelf.instance_count -= 1


class MemberFunc:
    instance_count = 0
    def __init__(self):
        MemberFunc.instance_count += 1
        self.my_member = 23434

    def my_func(self, pop, param):
        # can access to self from within PyOperator
        param.assertEqual(self.my_member, 23434)
        return True

    def __del__(self):
        MemberFunc.instance_count -= 1

class Functor:
    instance_count = 0
    def __init__(self):
        Functor.instance_count += 1
        self.my_member = 23434

    def __call__(self, pop, param):
        # can access to self from within PyOperator
        param.assertEqual(self.my_member, 23434)
        return True

    def __del__(self):
        Functor.instance_count -= 1

if __name__ == '__main__':
    unittest.main()



