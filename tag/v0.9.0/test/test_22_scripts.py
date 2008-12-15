#!/usr/bin/env python
#
# Purpose:
#    testing the bundled scripts
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedRevision: 475 $
# $LastChangedDate: 2006-10-04 01:00:38 -0500 (Wed, 04 Oct 2006) $
#
#
import simuOpt
simuOpt.setOptions(quiet=False)

from simuPOP import *
import unittest, os, sys, exceptions

sys.path.append('../scripts')

class TestScripts(unittest.TestCase):

    def TestLoadHapMap(self):
        if os.path.isfile('hapmap_21.bin'):
            os.remove('hapmap_21.bin')
        if os.path.isfile('hapmap_22.bin'):
            os.remove('hapmap_22.bin')
        import loadHapMap
        loadHapMap.loadHapMap([21, 22], '.')
        assert os.path.isfile('hapmap_21.bin')
        assert os.path.isfile('hapmap_22.bin')
        pass

    def testSimuCDCV(self):
        if AlleleType() != 'short':
             return
        import simuCDCV
        if not os.path.isdir('cdcv'):
            os.mkdir('cdcv')
        simuCDCV.simuCDCV(
            numDSL=2,
            initSpec=[[0.9] + [0.02]*5]*2,
            selModel='recessive',
            selModelAllDSL='multiplicative',
            selCoef=[0.02, 0.01],
            mutaModel='k-allele',
            maxAllele=255,
            mutaRate=[0.01]*2,
            initSize=1000,
            finalSize=2000,
            burnin=20,
            noMigrGen=20,
            mixingGen=20,
            growth='exponential',
            numSubPop=2,
            migrModel='island',
            migrRate=0.05,
            update=10,
            dispPlot=False,
            saveAt=[],
            savePop='cdcv.bin',
            resume='',
            resumeAtGen=[],
            name='cdcv',
            dryrun=False)

    def TestSimuForward(self):
        'Testing simuForward.py'
        #FIXME: this script currently does not work.
        import simuForward
        simuForward.simuForward(
            numChrom=2,
            numLoci=20,
            markerType='microsatellite',
            DSLafter=[2, 30],
            DSLdist=[0.5, 0.5],
            initSize=1000,
            endingSize=2000,
            growthModel='exponential',
            burninGen=20,
            introLen=20,
            splitGen=50,
            mixingGen=70,
            endingGen=80,
            introSel=0.01,
            minAlleleFreq=[0.001]*2,
            maxAlleleFreq=[0.20]*2,
            numSubPop=2,
            migrModel='stepping stone',
            migrRate=0.0001,
            fitness=[1, 1.01, 1.02],
            mlSelModel='none',
            mutaRate=0.0001,
            recRate=[0.0005],
            savedGen=2,
            numOffspring=2,
            numOffMode='constant',
            dryrun=False,
            savePop=[],
            filename='forward.bin',
            format='bin')

    def testDemoNonRandomMating_cpp(self):
        pass

if __name__ == '__main__':
    unittest.main()




