#!/usr/bin/env python
#
# Testing plotting with R
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision: 149 $
# $LastChangedDate: 2006-02-03 15:51:04 -0600 (Fri, 03 Feb 2006) $
#

import simuOpt
simuOpt.setOptions(quiet=True)

from simuPOP import *
import unittest, sys, os, exceptions

from time import sleep

hasRPy = True
try:
    from simuPOP.plotter import *
except exceptions.ImportError:
    print "simuRPy can not be imported. Either rpy or r is not installed properly."
    hasRPy = False

class TestPlotter(unittest.TestCase):
    def testDerivedArgs(self):
        'Testing class derivedARgs'
        pop = population(0)
        pop.dvars().gen = 100
        args = derivedArgs(
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
        self.assertRaises(exceptions.ValueError, args.getArgs, 'par1', pop)
        self.assertEqual(args.getArgs('par', pop, blah=1), {'val': [1, 3], 'blah': 1})
        #
        self.assertEqual(args.getArgs('lines', pop, rep=1),
            {'expr': 100, 'some_var':4, 'main': 'g', 'lty':1, 'cex':2})
        #
        self.assertEqual(args.getLegendArgs('lines', pop, ['cex', 'pch'], 'rep', range(5)),
            {'cex': [1, 2, 3, 1, 2]})
        self.assertEqual(args.getLegendArgs('lines', pop, ['cex', 'pch', 'lty'], 'rep', range(5)),
            {'cex': [1, 2, 3, 1, 2], 'lty':[1, 1, 1, 1, 1]})
        self.assertEqual(args.getLegendArgs('lines', pop, ['cex', 'pch', 'lty'],
            ['rep', 'dim'], [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]),
            {'cex': [1, 1, 2, 2, 3, 3], 'lty':[1, 1, 1, 1, 1, 1], 'pch': [5, 6, 5, 6, 5, 6]})
        #
        args = derivedArgs(
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
        'Testing byRep parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]], mode=BY_PROBABILITY)
        stator = stat(popSize=1)
        simu.evolve(
            initOps = initSex(),
            preOps = [stator, migr],
            postOps = [varPlotter('subPopSize[0]', update=10,
                win=10, main="!'Three colorful lines, no legend, win=10, gen=%d' % gen")],
            gen = 30
        )
        sleep(1)
        plotter.r.dev_off()

    def testVarPlotterByRep(self):
        'Testing byRep parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = stat(popSize=1)
        simu.evolve(
            initOps = [initSex()],
            preOps = [stator, migr],
            postOps = [
             varPlotter('subPopSize', byRep=True, lty_dim=[1, 2, 3],
                main='3 rep, 3 colorful thick lines, ylabs differ',
                lwd=2, update=10,
                ylab_rep=['subPopSize (rep %d)' % x for x in range(3)])
             ],
             gen = 30
        )
        sleep(1)
        plotter.r.dev_off()

    def testVarPlotterByDim(self):
        'Testing byDim paramter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[200, 100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2],[.25,0]], mode=BY_PROBABILITY)
        stator = stat(popSize=1)
        simu.evolve(
            initOps = [initSex()],
            preOps = [stator, migr],
            postOps = [
             varPlotter('[x**2 for x in subPopSize]', ylab='sp', 
                 col_rep=['red', 'green'], update=10,
                 byDim=True, win=10, main='win=10, 2 dim, 2 colorful lines, legend',
                 legend=['a', 'b'])
             ],
             gen = 30
        )
        sleep(1)
        plotter.r.dev_off()

    def testVarPlotterSelectedRep(self):
        'Testing byDim paramter of varPlotter using selected replicates'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50, 50, 100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=5)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = stat(popSize=1)
        simu.evolve(
            initOps = [initSex()],
            preOps = [stator, migr],
            postOps = varPlotter('[x**2 for x in subPopSize]', ylab='sp', type='l',
                 col_rep=['red', 'green', 'blue'], byDim=True, update=10,
                 main='3 out of 5 reps, 3 dim plots', reps=[0, 2, 3]),
             gen = 30
        )
        sleep(1)
        plotter.r.dev_off()

    def testVarPlotterTogether(self):
        'Testing plotting all lines together using varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50, 50, 100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = stat(popSize=1)
        simu.evolve(
            initOps = [initSex()],
            preOps = [stator, migr],
            postOps = varPlotter('[x**2 for x in subPopSize]', ylab='sp',
                 col_rep=['red', 'green', 'blue'], update=10,
                 main='9 lines, col rep, lty dim',
                 lty_dim=range(1, 4), ylim=[0, 10000],
                 legend=['rep1-sp1', 'rep1-sp2', 'rep1-sp3', 'rep2-sp1',
                    'rep2-sp2', 'rep2-sp3', 'rep3-sp1', 'rep3-sp2', 'rep3-sp3']),
             gen = 30
        )
        sleep(1)
        plotter.r.dev_off()

    def testVarPlotterByRepDim(self):
        'Testing byDim paramter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[200, 100], ploidy=2, loci=[3, 4], infoFields = ['migrate_to']),
            randomMating(), rep=2)
        migr = migrator(rate=[[0, .2],[.25, 0]],
            mode=BY_PROBABILITY)
        stator = stat(popSize=1)
        simu.evolve(
            initOps = [initSex()],
            preOps = [stator, migr],
            postOps = [
             varPlotter('[x**2 for x in subPopSize]', ylab='sp', 
                 byDim=True, byRep=True, update=10,
                 win=10, main_repdim=['rep dim %d' % x for x in range(4)],
                 col_rep=['red', 'green'], lty_dim=range(1, 3))
             ],
             gen = 30
        )
        sleep(1)
        plotter.r.dev_off()

    def testVarPlotterSaveAs(self):
        'Testing saveAs parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields = ['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = stat(popSize=1)
        simu.evolve(
            initOps = [initSex()],
            preOps = [stator, migr],
            postOps = [
             varPlotter('[x**2 for x in subPopSize]', byRep=True,
                 win=20, update=10, ylim=[0, 10000],
                 main='Save as, 3 rep, 3 colors', ylab='sp',
                 saveAs='demo.eps')
             ],
             gen = 31
        )
        for f in range(10,31,10):
            self.assertEqual(os.path.isfile('demo_%d.eps' % f), True)
            os.remove('demo_%d.eps' % f)
        sleep(1)
        plotter.r.dev_off()

    def testVarPlotterPar(self):
        'Testing parameter passing of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields=['migrate_to']),
            randomMating(), rep=2)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = stat(popSize=1)
        simu.evolve(
            initOps = [initSex()],
            preOps = [stator, migr],
            postOps = [
             varPlotter('subPopSize', byRep=True, ylim=[0,100],
                 ylab='sp', win=10, update=10,
                 par_mfrow=[1, 2], main='mfrow=[1, 3]')
             ],
             gen = 30
        )
        sleep(1)
        plotter.r.dev_off()

    def testVarPlotterHook(self):
        'Testing ylim parameter of varPlotter'
        if not hasRPy:
            return True
        simu = simulator(
            population(size=[50,50,100], ploidy=2, loci=[3,4], infoFields=['migrate_to']),
            randomMating(), rep=3)
        migr = migrator(rate=[[0,.2,.1],[.25,0,.1],[.1,.2,0]],
            mode=BY_PROBABILITY)
        stator = stat(popSize=1)
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
            initOps = [initSex()],
            preOps = [stator, migr],
            postOps = varPlotter('subPopSize', byRep=True, ylim=[0,100], plot_axes=False,
                 ylab='sp', col_dim=['red', 'blue', 'black'],
                 update=5, main_rep=['dimension %d' % x for x in range(3)],
                 preHook=setPar,
                 postHook=mtext,
                 plotHook=drawFrame),
             gen = 30
        )
        sleep(1)
        plotter.r.dev_off()

    def testScatterPlotter(self):
        'Testing scatterPlotter'
        import random
        pop = population([100, 200], infoFields=['x', 'y'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps = [
                inheritTagger(PATERNAL, infoFields='x'),
                inheritTagger(MATERNAL, infoFields='y'),
            ],
            postOps = scatterPlotter(['x', 'y'], main='B/W, 300 points', step=2),
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

    def testScatterPlotterSP(self):
        'Testing scatterPlotter with multiple virtual subpopulations'
        import random
        pop = population([100, 200], infoFields=['x', 'y'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps = [
                inheritTagger(PATERNAL, infoFields='x'),
                inheritTagger(MATERNAL, infoFields='y'),
            ],
            postOps = [
                scatterPlotter(['x', 'y'],
                    subPops = [(0, 0), (0, 1), (1, 0), (1, 1)],
                    main='Color, 300 points, left right does not mix'),
            ],
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

    def testScatterPlotterSubSet(self):
        'Testing scatterPlotter with partial individuals'
        import random
        pop = population([100, 200], infoFields=['x', 'y'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps = [
                inheritTagger(PATERNAL, infoFields='x'),
                inheritTagger(MATERNAL, infoFields='y'),
            ],
            postOps = [
                scatterPlotter(['x', 'y'],
                    subPops = [(0, 0), (0, 1)],
                    xlim = [0, 1],
                    main='Twoo colors, 100 points xlim=[0, 1]',
                    legend = ['MALE', 'FEMALE']),
                #pause(),
            ],
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

    def testInfoPlotterBase(self):
        'Testing basic histogram using infoPlotter'
        import random
        pop = population([500, 1000], infoFields=['x'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x')
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps = inheritTagger(PATERNAL, infoFields='x'),
            postOps = histPlotter(infoFields='x', main='histogram of x in green',
                    angle=60, col='green', step=2),
                #pause(),
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

    def testInfoPlotterFields(self):
        'Testing stat plotter with multiple fields and subpopulations'
        import random
        pop = population([500, 1000], infoFields=['x', 'y'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps = [
                inheritTagger(PATERNAL, infoFields='x'),
                inheritTagger(MATERNAL, infoFields='y'),
            ],
            postOps = [
                histPlotter(infoFields=['x', 'y'],
                    subPops = [(0, 0), (0, 1)],
                    col_fld=['red', 'green'],
                    density_sp=[5, 20], step=2,
                    main_spfld=['Field x, MALE', 'Field y, MALE',
                        'Field x, FEMALE', 'Field y, FEMALE'],
                ),
                #pause(),
            ],
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

    def testInfoPlotterQQplot(self):
        'Testing barplotter with multiple fields and subpopulations'
        import random
        pop = population([500, 100], infoFields=['x', 'y'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        def qqline(r, data=None, **kwargs):
            r.qqline(data)
        simu.evolve(
            duringOps = [
                inheritTagger(PATERNAL, infoFields='x'),
                inheritTagger(MATERNAL, infoFields='y'),
            ],
            postOps = [
                qqPlotter(infoFields=['x', 'y'],
                    subPops = [(0, 0), (0, 1)],
                    col_fld=['red', 'green'],
                    pch_sp=[1, 2],
                    main_spfld=['Field x, MALE', 'Field y, MALE',
                        'Field x, FEMALE', 'Field y, FEMALE'],
                    plotHook=qqline, step=2
                ),
                #pause(),
            ],
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

    def testInfoPlotterNoFunc(self):
        'Testing the stat plotter when no function is specified'
        import random
        pop = population([500, 100], infoFields=['x', 'y'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        def qqplot(r, data=None, field=None, subPop=None, **kwargs):
            r.qqnorm(data, main='QQ drawn in plotHook fld: %s sp: %s' % (field, subPop))
            r.qqline(data)
        simu.evolve(
            duringOps = [
                inheritTagger(PATERNAL, infoFields='x'),
                inheritTagger(MATERNAL, infoFields='y'),
            ],
            postOps = [
                infoPlotter(infoFields=['x', 'y'],
                    subPops = [(0, 0), (0, 1)],
                    plotHook=qqplot, step=2
                ),
                #pause(),
            ],
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()


    def testBoxPlotterBase(self):
        'Testing the base boxplotter'
        import random
        pop = population([500, 100], infoFields=['x'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x')
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps = inheritTagger(PATERNAL, infoFields='x'),
            postOps = [
                boxPlotter(infoFields='x', xlab='x', main='boxplot for field x',
                    step=2), 
                #pause(),
            ],
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

    def testBoxPlotterFields(self):
        'Testing barplotter with multiple fields and subpopulations'
        import random
        pop = population([500, 100], infoFields=['x', 'y'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps = [
                inheritTagger(PATERNAL, infoFields='x'),
                inheritTagger(MATERNAL, infoFields='y'),
            ],
            postOps = boxPlotter(infoFields=['x', 'y'], step=2),
                #pause(),
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

    def testBoxPlotterFieldsAndSubPop(self):
        'Testing boxPlotter with both fields and subpopulation'
        import random
        pop = population([500, 100], infoFields=['x', 'y'],
            subPopNames=['sp1', 'sp2'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps = [
                inheritTagger(PATERNAL, infoFields='x'),
                inheritTagger(MATERNAL, infoFields='y'),
            ],
            postOps = [
                boxPlotter(infoFields=['x', 'y'], 
                    subPops=[(0, 0), (0, 1)], style='quantile',
                    step=2, main='Boxplot with 2 fields x 2 subpops',
                    col='green', horizontal=True,
                ),
                #pause(),
            ],
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

    def testBoxPlotterByField(self):
        'Testing boxPlotter separated by information field'
        import random
        pop = population([500, 100], infoFields=['x', 'y'],
            subPopNames=['sp1', 'sp2'])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps = [
                inheritTagger(PATERNAL, infoFields='x'),
                inheritTagger(MATERNAL, infoFields='y'),
            ],
            postOps = boxPlotter(infoFields=['x', 'y'], byField=True,
                    subPops=[(0,0), (0,1)],
                    step=2, col_fld=['yellow', 'green'],
                ),
                #pause(),
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

    def testBoxPlotterBySubPop(self):
        'Testing boxPlotter separated by subpopulation'
        import random
        pop = population([500, 100], infoFields=['x', 'y'], subPopNames=['', ''])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps = [
                inheritTagger(PATERNAL, infoFields='x'),
                inheritTagger(MATERNAL, infoFields='y'),
            ],
            postOps = boxPlotter(infoFields=['x', 'y'], bySubPop=True,
                    subPops=[(0,0), (0,1)],
                    step=2, horizontal_sp=[True, False],
                ),
                #pause(),
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

    def testBoxPlotterByFieldSubPop(self):
        'Testing boxPlotter separated by both field and subpopulation'
        import random
        pop = population([500, 100], infoFields=['x', 'y'], subPopNames=['', ''])
        InitSex(pop)
        pop.setIndInfo([random.random() for i in range(100)], 'x', 0)
        pop.setIndInfo([1 + random.random() for i in range(100)], 'x', 1)
        pop.setIndInfo([random.random() for i in range(300)], 'y')
        pop.setVirtualSplitter(sexSplitter())
        simu = simulator(pop, randomMating())
        simu.evolve(
            duringOps = [
                inheritTagger(PATERNAL, infoFields='x'),
                inheritTagger(MATERNAL, infoFields='y'),
            ],
            postOps = boxPlotter(infoFields=['x', 'y'], bySubPop=True,
                    subPops=[(0,0), (0,1)], byField=True,
                    step=2, col_spfld = ['red', 'green', 'yellow', 'blue'],
                ),
                #pause(),
            gen = 5,
        )
        sleep(1)
        plotter.r.dev_off()

if __name__ == '__main__':
    unittest.main()
