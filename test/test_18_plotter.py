#!/usr/bin/env python
#
# Testing plotting with R
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision: 149 $
# $LastChangedDate: 2006-02-03 15:51:04 -0600 (Fri, 03 Feb 2006) $
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
from time import sleep

hasRPy = True
try:
    from simuPOP.plotter import *
except (ImportError, RuntimeError, AttributeError):
    print("simuRPy can not be imported. Either rpy or r is not installed properly.")
    hasRPy = False

class TestPlotter(unittest.TestCase):
    def testDerivedArgs(self):
        'Testing class derivedARgs'
        if not hasRPy:
            return True
        pop = Population(0)
        pop.dvars().gen = 100
        args = DerivedArgs(
            defaultFuncs = ['plot', 'lines'],
            allFuncs = ['plot', 'lines', 'par'],
            suffixes = ['rep', 'dim'],
            lty = 1,
            cex_rep=[1, 2, 3],
            pch_dim=[5,6,7,8],
            main = 'g',
            par_val=[1,3],
            par_blah_dim=[1,2],
            some_var=4,
            expr = '!gen',
        )
        self.assertEqual(args.getArgs('plot', pop),
            {'expr': 100, 'lty': 1, 'some_var':4, 'main':'g'})
        self.assertEqual(args.getArgs('plot', pop, dim=1),
            {'expr': 100, 'pch': 6, 'lty': 1, 'some_var':4, 'main':'g'})
        self.assertEqual(args.getArgs('plot', pop, dim=4),
            {'expr': 100, 'pch': 5, 'lty': 1, 'some_var':4, 'main':'g'})
        self.assertEqual(args.getArgs('plot', pop, rep=4, dim=5),
            {'expr': 100, 'cex':2, 'pch': 6, 'lty': 1, 'some_var':4, 'main':'g'})
        self.assertEqual(args.getArgs('par', pop, rep=4, dim=5),
            {'val':[1,3], 'blah': 2})
        #
        self.assertRaises(ValueError, args.getArgs, 'par1', pop)
        self.assertEqual(args.getArgs('par', pop, blah=1), {'val': [1, 3], 'blah': 1})
        #
        self.assertEqual(args.getArgs('lines', pop, rep=1),
            {'expr': 100, 'some_var':4, 'main': 'g', 'lty':1, 'cex':2})
        #
        self.assertEqual(args.getLegendArgs('lines', pop, ['cex', 'pch'], 'rep', list(range(5))),
            {'cex': [1, 2, 3, 1, 2]})
        self.assertEqual(args.getLegendArgs('lines', pop, ['cex', 'pch', 'lty'], 'rep', list(range(5))),
            {'cex': [1, 2, 3, 1, 2], 'lty':[1, 1, 1, 1, 1]})
        self.assertEqual(args.getLegendArgs('lines', pop, ['cex', 'pch', 'lty'],
            ['rep', 'dim'], [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]),
            {'cex': [1, 1, 2, 2, 3, 3], 'lty':[1, 1, 1, 1, 1, 1], 'pch': [5, 6, 5, 6, 5, 6]})
        #
        args = DerivedArgs(
            defaultFuncs = ['plot', 'lines'],
            allFuncs = ['plot', 'lines', 'par'],
            suffixes = ['rep', 'dim'],
            defaultParams = {'plot_lty': 1, 'par_blah': 2, 'plot_blah': 'blah'},
            lty = 1,
            cex_rep=[1, 2, 3],
            main = 'g',
            par_val=[1,3],
            par_blah_dim=[1,2],
        )
        self.assertEqual(args.params, {'lty':1, 'cex_rep': [1,2,3], 'lty':1,
            'plot_blah': 'blah', 'par_blah_dim': [1,2], 'par_val':[1,3], 'main': 'g'})


    def testVarPlotterBase(self):
        'Testing byRep parameter of VarPlotter'
        if not hasRPy:
            return True
        simu = Simulator(
            Population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            rep=3)
        migr = Migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]], mode=BY_PROBABILITY)
        stator = Stat(popSize=1)
        simu.evolve(
            initOps = InitSex(),
            preOps = [stator, migr],
            matingScheme=RandomMating(),
            postOps = [VarPlotter('subPopSize[0]', update=10,
                win=10, main="!'Three colorful lines, no legend, win=10, gen=%d' % gen")],
            gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterByRep(self):
        'Testing byRep parameter of VarPlotter'
        if not hasRPy:
            return True
        simu = Simulator(
            Population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            rep=3)
        migr = Migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = Stat(popSize=1)
        simu.evolve(
            initOps = [InitSex()],
            preOps = [stator, migr],
            matingScheme=RandomMating(),
            postOps = [
             VarPlotter('subPopSize', byRep=True, lty_dim=[1, 2, 3],
                main='3 rep, 3 colorful thick lines, ylabs differ',
                lwd=2, update=10,
                ylab_rep=['subPopSize (rep %d)' % x for x in range(3)])
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterByDim(self):
        'Testing byDim paramter of VarPlotter'
        if not hasRPy:
            return True
        simu = Simulator(
            Population(size=[200, 100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            rep=3)
        migr = Migrator(rate=[[0,.2],[.25,0]], mode=BY_PROBABILITY)
        stator = Stat(popSize=1)
        simu.evolve(
            initOps = [InitSex()],
            preOps = [stator, migr], 
            matingScheme=RandomMating(),
            postOps = [
             VarPlotter('[x**2 for x in subPopSize]', ylab='sp', 
                 col_rep=['red', 'green'], update=10,
                 byDim=True, win=10, main='win=10, 2 dim, 2 colorful lines, legend',
                 legend=['a', 'b'])
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterSelectedRep(self):
        'Testing byDim paramter of VarPlotter using selected replicates'
        if not hasRPy:
            return True
        simu = Simulator(
            Population(size=[50, 50, 100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            rep=5)
        migr = Migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = Stat(popSize=1)
        simu.evolve(
            initOps = [InitSex()],
            matingScheme=RandomMating(),
            preOps = [stator, migr],
            postOps = VarPlotter('[x**2 for x in subPopSize]', ylab='sp', type='l',
                 col_rep=['red', 'green', 'blue'], byDim=True, update=10,
                 main='3 out of 5 reps, 3 dim plots', reps=[0, 2, 3]),
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterTogether(self):
        'Testing plotting all lines together using VarPlotter'
        if not hasRPy:
            return True
        simu = Simulator(
            Population(size=[50, 50, 100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            rep=3)
        migr = Migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = Stat(popSize=1)
        simu.evolve(
            initOps = [InitSex()],
            preOps = [stator, migr],
            matingScheme=RandomMating(),
            postOps = VarPlotter('[x**2 for x in subPopSize]', ylab='sp',
                 col_rep=['red', 'green', 'blue'], update=10,
                 main='9 lines, col rep, lty dim',
                 lty_dim=list(range(1, 4)), ylim=[0, 10000],
                 legend=['rep1-sp1', 'rep1-sp2', 'rep1-sp3', 'rep2-sp1',
                    'rep2-sp2', 'rep2-sp3', 'rep3-sp1', 'rep3-sp2', 'rep3-sp3']),
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterByRepDim(self):
        'Testing byDim paramter of VarPlotter'
        if not hasRPy:
            return True
        simu = Simulator(
            Population(size=[200, 100], ploidy=2, loci=[3, 4], infoFields = ['migrate_to']),
            rep=2)
        migr = Migrator(rate=[[0, .2],[.25, 0]],
            mode=BY_PROBABILITY)
        stator = Stat(popSize=1)
        simu.evolve(
            initOps = [InitSex()],
            preOps = [stator, migr],
            matingScheme=RandomMating(),
            postOps = [
             VarPlotter('[x**2 for x in subPopSize]', ylab='sp', 
                 byDim=True, byRep=True, update=10,
                 win=10, main_repdim=['rep dim %d' % x for x in range(4)],
                 col_rep=['red', 'green'], lty_dim=list(range(1, 3)))
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterSaveAs(self):
        'Testing saveAs parameter of VarPlotter'
        if not hasRPy:
            return True
        simu = Simulator(
            Population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            rep=3)
        migr = Migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = Stat(popSize=1)
        simu.evolve(
            initOps = [InitSex()],
            preOps = [stator, migr],
            matingScheme=RandomMating(),
            postOps = [
             VarPlotter('[x**2 for x in subPopSize]', byRep=True,
                 win=20, update=10, ylim=[0, 10000],
                 main='Save as, 3 rep, 3 colors', ylab='sp',
                 saveAs='demo.eps')
             ],
             gen = 31
        )
        sleep(1)
        for x in range(3):
            simu.population(x).dvars().gen = 0
        # automatically create directory
        simu.evolve(
            initOps = [InitSex()],
            preOps = [stator, migr],
            matingScheme=RandomMating(),
            postOps = [
             VarPlotter('[x**2 for x in subPopSize]', byRep=True,
                 win=20, update=10, ylim=[0, 10000],
                 main='Save as, 3 rep, 3 colors', ylab='sp',
                 saveAs='NonExistDir/demo.eps')
             ],
             gen = 31
        )
        # if failed to create directory
        self.assertRaises(ValueError, simu.evolve,
            initOps = [InitSex()],
            preOps = [stator, migr],
            matingScheme=RandomMating(),
            postOps = [
             VarPlotter('[x**2 for x in subPopSize]', byRep=True,
                 win=20, update=10, ylim=[0, 10000],
                 main='Save as, 3 rep, 3 colors', ylab='sp',
                 saveAs='NonExistDir/demo_30.eps/abs.eps')
             ],
             gen = 31
        )
        for f in range(10,31,10):
            self.assertEqual(os.path.isfile('demo_%d.eps' % f), True)
            self.assertEqual(os.path.isfile(os.path.join('NonExistDir', 'demo_%d.eps' % f)), True)
            os.remove('demo_%d.eps' % f)
            os.remove(os.path.join('NonExistDir', 'demo_%d.eps' % f))
        os.rmdir('NonExistDir')
        r.dev_off()

    def testVarPlotterPar(self):
        'Testing parameter passing of VarPlotter'
        if not hasRPy:
            return True
        simu = Simulator(
            Population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields=['migrate_to']),
            rep=2)
        migr = Migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = Stat(popSize=1)
        simu.evolve(
            initOps = [InitSex()],
            preOps = [stator, migr],
            matingScheme=RandomMating(),
            postOps = [
             VarPlotter('subPopSize', byRep=True, ylim=[0,100],
                 ylab='sp', win=10, update=10,
                 par_mfrow=[1, 2], main='mfrow=[1, 3]')
             ],
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testVarPlotterHook(self):
        'Testing ylim parameter of VarPlotter'
        if not hasRPy:
            return True
        simu = Simulator(
            Population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields=['migrate_to']),
            rep=3)
        migr = Migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = Stat(popSize=1)
        def setPar(r):
            r.par(mar=[2]*4)
        def mtext(r):
            r.mtext('A marginal text', outer=True, side=3)
        def drawFrame(r, **kwargs):
            r.axis(1)
            r.axis(2)
            r.grid()
        #
        simu.evolve(
            initOps = [InitSex()],
            preOps = [stator, migr],
            matingScheme=RandomMating(),
            postOps = VarPlotter('subPopSize', byRep=True, ylim=[0,100], plot_axes=False,
                 ylab='sp', col_dim=['red', 'blue', 'black'],
                 update=5, main_rep=['dimension %d' % x for x in range(3)],
                 preHook=setPar,
                 postHook=mtext,
                 plotHook=drawFrame),
             gen = 30
        )
        sleep(1)
        r.dev_off()

    def testScatterPlotter(self):
        'Testing ScatterPlotter'
        if not hasRPy:
            return True
        import random
        pop = Population([100, 200], infoFields=['x', 'y'])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(SexSplitter())
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x'),
                InheritTagger(MATERNAL, infoFields='y'),
            ]),
            postOps = ScatterPlotter(['x', 'y'], main='B/W, 300 points', step=2),
            gen = 5,
        )
        sleep(1)
        r.dev_off()

    def testScatterPlotterSP(self):
        'Testing ScatterPlotter with multiple virtual subpopulations'
        if not hasRPy:
            return True
        import random
        pop = Population([100, 200], infoFields=['x', 'y'])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(SexSplitter())
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x'),
                InheritTagger(MATERNAL, infoFields='y'),
            ]),
            postOps = [
                ScatterPlotter(['x', 'y'],
                    subPops = [(0, 0), (0, 1), (1, 0), (1, 1)],
                    main='Color, 300 points, left right does not mix'),
            ],
            gen = 5,
        )
        sleep(1)
        r.dev_off()

    def testScatterPlotterSubSet(self):
        'Testing ScatterPlotter with partial individuals'
        if not hasRPy:
            return True
        import random
        pop = Population([100, 200], infoFields=['x', 'y'])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(SexSplitter())
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x'),
                InheritTagger(MATERNAL, infoFields='y'),
            ]),
            postOps = [
                ScatterPlotter(['x', 'y'],
                    subPops = [(0, 0), (0, 1)],
                    xlim = [0, 1],
                    main='Twoo colors, 100 points xlim=[0, 1]',
                    legend = ['Male', 'Female']),
                #Pause(),
            ],
            gen = 5,
        )
        sleep(1)
        r.dev_off()

    def testInfoPlotterBase(self):
        'Testing basic histogram using InfoPlotter'
        if not hasRPy:
            return True
        import random
        pop = Population([500, 1000], infoFields=['x'])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x')
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x')]),
            postOps = HistPlotter(infoFields='x', main='histogram of x in green',
                    angle=60, col='green', step=2),
                #Pause(),
            gen = 5,
        )
        sleep(1)
        r.dev_off()

    def testInfoPlotterFields(self):
        'Testing Stat plotter with multiple fields and subpopulations'
        if not hasRPy:
            return True
        import random
        pop = Population([500, 1000], infoFields=['x', 'y'])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(SexSplitter())
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x'),
                InheritTagger(MATERNAL, infoFields='y'),
            ]),
            postOps = [
                HistPlotter(infoFields=['x', 'y'],
                    subPops = [(0, 0), (0, 1)],
                    col_fld=['red', 'green'],
                    density_sp=[5, 20], step=2,
                    main_spfld=['Field x, Male', 'Field y, Male',
                        'Field x, Female', 'Field y, Female'],
                ),
                #Pause(),
            ],
            gen = 5,
        )
        sleep(1)
        r.dev_off()

    def testInfoPlotterQQplot(self):
        'Testing barplotter with multiple fields and subpopulations'
        if not hasRPy:
            return True
        import random
        pop = Population([500, 100], infoFields=['x', 'y'])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(SexSplitter())
        def qqline(r, data=None, **kwargs):
            r.qqline(data)
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x'),
                InheritTagger(MATERNAL, infoFields='y'),
            ]),
            postOps = [
                QQPlotter(infoFields=['x', 'y'],
                    subPops = [(0, 0), (0, 1)],
                    col_fld=['red', 'green'],
                    pch_sp=[1, 2],
                    main_spfld=['Field x, Male', 'Field y, Male',
                        'Field x, Female', 'Field y, Female'],
                    plotHook=qqline, step=2
                ),
                #Pause(),
            ],
            gen = 5,
        )
        sleep(1)
        r.dev_off()

    def testInfoPlotterNoFunc(self):
        'Testing the Stat plotter when no function is specified'
        if not hasRPy:
            return True
        import random
        pop = Population([500, 100], infoFields=['x', 'y'])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(SexSplitter())
        def qqplot(r, data=None, field=None, subPop=None, **kwargs):
            r.qqnorm(data, main='QQ drawn in plotHook fld: %s sp: %s' % (field, subPop))
            r.qqline(data)
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x'),
                InheritTagger(MATERNAL, infoFields='y'),
            ]),
            postOps = [
                InfoPlotter(infoFields=['x', 'y'],
                    subPops = [(0, 0), (0, 1)],
                    plotHook=qqplot, step=2
                ),
                #Pause(),
            ],
            gen = 5,
        )
        sleep(1)
        r.dev_off()


    def testBoxPlotterBase(self):
        'Testing the base boxplotter'
        if not hasRPy:
            return True
        import random
        pop = Population([500, 100], infoFields=['x'])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x')
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x')]),
            postOps = [
                BoxPlotter(infoFields='x', xlab='x', main='boxplot for field x',
                    step=2), 
                #Pause(),
            ],
            gen = 5,
        )
        sleep(1)
        r.dev_off()

    def testBoxPlotterFields(self):
        'Testing barplotter with multiple fields and subpopulations'
        if not hasRPy:
            return True
        import random
        pop = Population([500, 100], infoFields=['x', 'y'])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(SexSplitter())
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x'),
                InheritTagger(MATERNAL, infoFields='y'),
            ]),
            postOps = BoxPlotter(infoFields=['x', 'y'], step=2),
                #Pause(),
            gen = 5,
        )
        sleep(1)
        r.dev_off()

    def testBoxPlotterFieldsAndSubPop(self):
        'Testing BoxPlotter with both fields and subpopulation'
        if not hasRPy:
            return True
        import random
        pop = Population([500, 100], infoFields=['x', 'y'],
            subPopNames=['sp1', 'sp2'])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(SexSplitter())
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x'),
                InheritTagger(MATERNAL, infoFields='y'),
            ]),
            postOps = [
                BoxPlotter(infoFields=['x', 'y'], 
                    subPops=[(0, 0), (0, 1)], style='quantile',
                    step=2, main='Boxplot with 2 fields x 2 subpops',
                    col='green', horizontal=True,
                ),
                #Pause(),
            ],
            gen = 5,
        )
        sleep(1)
        r.dev_off()

    def testBoxPlotterByField(self):
        'Testing BoxPlotter separated by information field'
        if not hasRPy:
            return True
        import random
        pop = Population([500, 100], infoFields=['x', 'y'],
            subPopNames=['sp1', 'sp2'])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(SexSplitter())
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x'),
                InheritTagger(MATERNAL, infoFields='y'),
            ]),
            postOps = BoxPlotter(infoFields=['x', 'y'], byField=True,
                    subPops=[(0,0), (0,1)],
                    step=2, col_fld=['yellow', 'green'],
                ),
                #Pause(),
            gen = 5,
        )
        sleep(1)
        r.dev_off()

    def testBoxPlotterBySubPop(self):
        'Testing BoxPlotter separated by subpopulation'
        if not hasRPy:
            return True
        import random
        pop = Population([500, 100], infoFields=['x', 'y'], subPopNames=['', ''])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(SexSplitter())
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x'),
                InheritTagger(MATERNAL, infoFields='y'),
            ]),
            postOps = BoxPlotter(infoFields=['x', 'y'], bySubPop=True,
                    subPops=[(0,0), (0,1)],
                    step=2, horizontal_sp=[True, False],
                ),
                #Pause(),
            gen = 5,
        )
        sleep(1)
        r.dev_off()

    def testBoxPlotterByFieldSubPop(self):
        'Testing BoxPlotter separated by both field and subpopulation'
        if not hasRPy:
            return True
        import random
        pop = Population([500, 100], infoFields=['x', 'y'], subPopNames=['', ''])
        initSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(SexSplitter())
        pop.evolve(
            matingScheme=RandomMating(ops=[MendelianGenoTransmitter(),
                InheritTagger(PATERNAL, infoFields='x'),
                InheritTagger(MATERNAL, infoFields='y'),
            ]),
            postOps = BoxPlotter(infoFields=['x', 'y'], bySubPop=True,
                    subPops=[(0,0), (0,1)], byField=True,
                    step=2, col_spfld = ['red', 'green', 'yellow', 'blue'],
                ),
                #Pause(),
            gen = 5,
        )
        sleep(1)
        r.dev_off()

if __name__ == '__main__':
    unittest.main()
