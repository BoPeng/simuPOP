#!/usr/bin/env python
#
# Purpose:
#   Test the python binding of coaSim (http://www.birc.dk/Software/CoaSim/)
#
# Author:
#   Bo Peng
#
# $LastChangedRevision: 154 $
# $LastChangedDate: 2006-02-06 01:05:52 -0600 (Mon, 06 Feb 2006) $
# 

from coaSim import *
import unittest, os, sys, exceptions

class TestCoaSim(unittest.TestCase):
  
  def assertMarkerProperties(self, m, pos, name, default, run_first):
    'Testing behavior to all marker types'
    self.assertEqual(m.position(), pos)
    m.position(0.5)
    self.assertEqual(m.position(), 0.5)
    self.assertEqual(m.type(), name)
    self.assertEqual(m.default_value(), default)
    # wrong position
    self.assertRaises(exceptions.ValueError,
      m.position, 2)
    # run_first (??)
    self.assertEqual(m.run_first(), run_first)
    
  def testSNPMarker(self):
    'Testing SNP markers'
    a = SNPMarker(0.11, 0.01, 0.90)
    # 
    self.assertMarkerProperties(a, pos=0.11, 
      name='snp', default=0, run_first=False)
    #
    self.assertEqual(a.low_freq(), 0.01)
    self.assertEqual(a.high_freq(), 0.90)
    # can not add_value
    self.assertRaises(exceptions.ValueError,
      a.add_value, 2)

  def testMicroSatelliteMarker(self):
    'Testing microsatellite markers'
    a = MicroSatelliteMarker(0.11, theta=1e-6, K=50)
    #
    self.assertMarkerProperties(a, pos=0.11, 
      name='ms', default=0, run_first=False)
    #
    self.assertEqual( a.theta(), 1e-6)
    self.assertEqual( a.K(), 50)
    # add value should be allowed FIXME: why add_value?

  def testTraitMarker(self):
    'Testing trait markers'
    a = TraitMarker(0.11, low_freq=0.01, high_freq=0.90)
    #
    self.assertMarkerProperties(a, pos=0.11, 
      name='trait', default=0, run_first=True)
    #
    self.assertEqual(a.low_freq(), 0.01)
    self.assertEqual(a.high_freq(), 0.90)
    # can not add value 
    self.assertRaises(exceptions.ValueError,
      a.add_value, 2)

  def testDescender(self):
    'Testing descender'
    #dec = Descender(conf)
    #dec.evolve()

  def assertEpochProperties(self, epo, st, et):
    'Assert properties for all epochs'
    self.assertEqual(epo.start_time(), st)
    self.assertEqual(epo.end_time(), et)
    self.assertEqual(epo.earliest_event(), 0)
    
  def testMigrationEvent(self):
    'Testing migration events (epochs.hh, configuration.hh)'
    a = Migration(source=0, destination=1,
      migration_rate=0.001, start_time=0,
      end_time= 0.5)
    #
    self.assertEpochProperties(a, st=0, et=0.5)
    self.assertEqual(a.source(), 0)
    self.assertEqual(a.destination(), 1)
    self.assertEqual(a.migration_rate(), 0.001)
  
  def testBottleNeckEpoch(self):
    'Testing BottleNeckEpoch (epochs.hh)'
    a = BottleNeckEpoch(population=0,  
      scale_fraction=1.2, 
      start_time=0, end_time= 0.5)
    #
    self.assertEpochProperties(a, st=0, et=0.5)
    self.assertEqual( a.scale_fraction(), 1.2)
    
  def testGrowthEpoch(self):
    'Testing GrowthEpoch (epochs.hh)'
    a = GrowthEpoch(population=0,  
      beta=1.2, start_time=0, end_time= 0.5)
    #
    self.assertEpochProperties(a, st=0, et=0.5)
    self.assertEqual( a.beta(), 1.2)    
   
  def testScheduler(self):
    'Testing scheduler (configuration.hh)'
    #a = Scheduler()
    #a.add_event(e)
    #a.remove_event(e)
    #print next_event
    
  def testConfiguration(self):
    'Testing configure object (configuration.hh)'
    conf = Configuration(
      popSizes=[1000,2000,3000], 
      markers = [
        SNPMarker(0.1, 0.1, 0.9), 
        SNPMarker(0.2, 0.1, 0.9)
        ],
      epochs = [ GrowthEpoch(0, 1.2, 0, 0.5)], 
      rho=0.1, Q=0.01, gamma=0.001, growth=0.0001)
    # 
    self.assertEqual( conf.pop_sizes(), (1000,2000,3000) )
    # FIXME: is no_leaves always total pop size?
    self.assertEqual( conf.no_leaves(), 6000)
    self.assertEqual( conf.no_markers(), 2)
    self.assertEqual( conf.position(0), 0.1)
    self.assertEqual( conf.position(1), 0.2)
    self.assertRaises( exceptions.IndexError,
      conf.position, 2)
    self.assertEqual( conf.rho(), 0.1)
    self.assertEqual( conf.Q(), 0.01)
    self.assertEqual( conf.gamma(), 0.001)
    self.assertEqual( conf.growth(), 0.0001)
    #
    self.assertEqual( conf.marker(0).type(), 'snp')
    self.assertRaises( exceptions.IndexError,
      conf.marker, 2)
    # FIXME: what are first_marker and plain_marker
    assert not conf.is_first_marker(0)
    assert not conf.is_first_marker(1)
    assert conf.is_plain_marker(0)
    assert conf.is_plain_marker(1)
    self.assertRaises( exceptions.IndexError,
      conf.first_marker, 0)
    self.assertEqual(conf.plain_marker(0).position(), 0.1)
    self.assertEqual(conf.plain_marker(1).position(), 0.2)

  def testDistributionFunctions(self):
    'Testing distribution functions'
    # Should we simply trust these functions?
    expdev(2.)
    expdev(2, .2)
    expdist(.2, 3.)
    uniform()
    uniform(1,2)
    # what is this guy?
    uniform(1,2,3)  
    random_sign()
    irand(10)
    # can not test two_int_rand


  def testPopulation(self):
    'Testing population'
    # need ARG...

  def testSimulate(self):
    'Test simulation'
    # need to translate the following into python
## 
## (define markers
##   (list (snpMarker   0.1 0.1 0.9)
## 	(snpMarker   0.2 0.1 0.9)
## 	(traitMarker 0.3 0.18 0.22)
## 	(msMarker    0.4 0.5 20)
## 	(snpMarker   0.6 0.1 0.9)))
## 
## (define arg (simulate markers 10 :rho 50 :randomSeed 100))
## 
## (define seqs (sequences arg))
## (display seqs)(newline)
## (newline)
## 
## (define retiredIntervals (intervals arg))
## (define startPositions (map intervalStart retiredIntervals))
## (define endPositions (map intervalEnd retiredIntervals))
## (define branchLengths (map totalBranchLength retiredIntervals))
## (display branchLengths)(newline)
## (newline)
## 
## 
## (define (print x) (display x)(display " "))
## (define (combine f g) (lambda (x) (f (g x))))
## (map (combine print totalBranchLength) (intervals arg))(newline)
## (map (combine print totalBranchLength) (localTrees arg))(newline)
## (newline)
## 
## (display (simulateSequences markers 10 :rho 50 :randomSeed 100))(newline)
## 
## (newline)
## (display (simulateSequences markers 
## 			     '(population 1 (merge 0.2
## 						   (population 1 (sample 10))
## 						   (population 1 (sample 10))))
##                              :randomSeed 100))
## (newline)
## 
## (newline)
## (newline)
## (display (simulateSequences markers 
## 			     '(population 1 (merge 0.2
## 						   (population .2 (sample 10))
## 						   (population .9 (sample 10))))
##                              :randomSeed 100))
## (newline)
## (newline)
## (display (simulateSequences markers 
## 			     '(population 1 (merge 0.2
## 						   (population p1 .2 (sample 10))
## 						   (population p2 .9 (sample 10))))
## 			     :migration '((migration p1 p2 0.1 0 0.2)
## 					  (migration p2 p1 0.2 0 0.2))
##                              :randomSeed 100))
## (newline)
## 
## 
## (newline)
## (display "now testing validation of incorrect input\n")
## (display "this will print some error messages\n")
## (newline)
## 
## (catch 'wrongNumberOfArgs
##        (lambda () (simulate))
##        (lambda (key . args) (display key)(newline)))
## 
## (catch 'wrongNumberOfArgs
##        (lambda () (simulate markers))
##        (lambda (key . args) (display key)(newline)))
## 
## (catch #t
##        (lambda () (simulate 0 10))
##        (lambda (key . args) (display key)(newline)))
## 
## (define overlappingMarkers
##   (list (snpMarker 0 0.1 0.9) (snpMarker 0 0.1 0.9)))
## (catch 'outOfSequence
##        (lambda () (simulate overlappingMarkers 10))
##        (lambda (key . args) (display key)(display " ")(display args)(newline)))
## 
## (define outOfSequenceMarkers
##   (list (snpMarker 0.1 0.1 0.9) (snpMarker 0.0 0.1 0.9)))
## (catch 'outOfSequence
##        (lambda () (simulate outOfSequenceMarkers 10))
##        (lambda (key . args) (display key)(display " ")(display args)(newline)))
## 
## (catch 'nonPositiveSampleSize
##        (lambda () (simulate '() 10))
##        (lambda (key . args) (display key)(display " ")(display args)(newline)))
    
  def testPopulationStructure(self):
    'Testing population structure'
    # need to translate the following into python
## 
## (load "simulate.scm")
## (useModules (ice9 receive))
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; Small test framework ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## (defmacro unless (pred . body)
##   `(if ,pred
##        '()
##        (begin
## 	 ,@body)))
## 
## (define (runTest test)
##   (catch #t
## 	 (lambda ()
## 	   (values (apply test '()) #t))
## 	 (lambda (key . args)
## 	   (values (list* key args) #f))))
## 
## (define (deftestFun name test testPrintedVer expectedRes)
##   (receive (res goodTestRunP) (runTest test)
## 	   (unless (and goodTestRunP 
## 			(equal? res expectedRes))
## 		   (forceOutput)
## 		   (format #t "===== Test FAILED: ~a =====~%" name)
## 		   (format #t "Expected:~%~y~%Got:~%~a~%"
## 			   expectedRes res))))
## 
## (defmacro deftest (name test expectedRes)
##   `(deftestFun ',name (lambda () ,test) ',test ,expectedRes))
## 
## 
## ;;(deftest TestSuccess 42 42)
## ;;(deftest TestFail 42 43)
## ;;(deftest TestFailExeption (throw 'error "hehe") 42)
## 
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; Test data and tests ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## (deftest compile1
##   (compile 
##    '(population p1 1 :beta 77 :epochs ((bottleneck 111 1 2)) (sample 1))
##    '())
##   '(:sampleSizes (1) :epochs ((bottleneck 0 1 0 1) (growth 0 77 0 1) (bottleneck 0 111 1 2))))
## 
## (deftest compile2
##   (compile 
##    '(population p1 5 (merge 11 
## 			    (population p2 2 :epochs (bottleneck 7 0 5) (sample 2))
## 			    (population p3 3 (sample 3))))
##    '((migration p2 p3 0.2 3 5) (migration p2 p3 0.7)))
##   '(:sampleSizes (2 3)
## 		  :epochs ((populationMerge 11 0 1)
## 			   (bottleneck 0 5 11 1)
## 			   (bottleneck 0 2 0 11)
## 			   (bottleneck 0 7 0 5)
## 			   (bottleneck 1 3 0 11)
## 			   (migration 0 1 0.2 3 5)
## 			   (migration 0 1 0.7 0 11))))
## 
## (deftest compile3
##   (compile 
##    '(population p1 7 
## 		 (merge 22
## 			(population p2 2 (sample 2))
## 			(population p3 5 
## 				    (merge 11
## 					   (population p4 2 (sample 2))
## 					   (population p5 3 (sample 3))))))
##    '())
##   '(:sampleSizes (2 2 3)
## 		  :epochs ((populationMerge 22 0 1)
## 			   (populationMerge 11 1 2)
## 			   (bottleneck 0 7 22 1)
## 			   (bottleneck 0 2 0 22)
## 			   (bottleneck 1 5 11 22)
## 			   (bottleneck 1 2 0 11) 
## 			   (bottleneck 2 3 0 11))))
## 
## (deftest compile4
##   (compile 
##    '(population 1 (sample 1))
##    '())
##   '(:sampleSizes (1) :epochs ((bottleneck 0 1 0 1))))
## 
## (deftest compile5
##   (compile 
##    '(population 1 :epochs () (sample 1))
##    '())
##   '(:sampleSizes (1) :epochs ((bottleneck 0 1 0 1))))
## 
## (deftest compile6
##   (compile 
##    '(population 1 (sample 1))
##    '())
##   '(:sampleSizes (1) :epochs ((bottleneck 0 1 0 1))))
## 
## (deftest compile7
##   (compile 
##    '(population 1 :beta 55 (sample 1))
##    '())
##   '(:sampleSizes (1) :epochs ((bottleneck 0 1 0 1) (growth 0 55 0 1))))
## 
## 
## (deftest implicitEnd
##   (compile
##    '(population 1 :epochs (bottleneck 2 0.5) (sample 1))
##    '())
##   '(:sampleSizes (1) :epochs ((bottleneck 0 1 0 1) (bottleneck 0 2 0.5 1))))
## 
## (deftest implicitEnd2
##   (compile
##    '(population 1 (merge 1
## 			 (population 2 :epochs (growth 2 0.5) (sample 1))
## 		         (population 4 (sample 1))))
##    '())
##   '(:sampleSizes (1 1)
##     :epochs ((populationMerge 1 0 1)
## 	     (bottleneck 0 1 1 1)
## 	     (bottleneck 0 2 0  1)
## 	     (growth 0 2 0.5 1)
## 	     (bottleneck 1 4 0 1))))
## 
## 
## (deftest implicitEnd3
##   (compile
##    '(population 1 :epochs (growth 2 0.5) (sample 1))
##    '())
##   '(:sampleSizes (1) :epochs ((bottleneck 0 1 0 1) (growth 0 2 0.5 1))))
## 
## 
## (deftest doubleBottleneck
##   (compile
##    '(population 1 :epochs ((bottleneck 2 0.004 0.04) (bottleneck .5 0.04)) (sample 10))
##    '())
##   '(:sampleSizes (10) :epochs ((bottleneck 0 1   0     1)
## 				(bottleneck 0 2   0.004  0.04)
## 				(bottleneck 0 0.5 0.04  1))))
##   
## 
## ;;;;(deftest compile4
## ;;;;  (compile 
## ;;;;   '(population p1 7 (merge 11
## ;;;;			    (population p2 2 (sample 2))
## ;;;;			    (population p3 5 (merge 2
## ;;;;						    (population p4 2 (sample 2))
## ;;;;						    (population p5 3 (sample 3)))))))
## ;;;;  '())
## ;;;;
## ;;;;(deftest compile5
## ;;;;  (compile 
## ;;;;   '(population p1 7 (merge 11
## ;;;;			    (population p2 2 (sample 2))
## ;;;;			    (merge 2
## ;;;;				   (population p4 2 (sample 2))
## ;;;;				   (population p5 3 (sample 3))))))
## ;;;;  '())
## ;;;;
## ;;;;(deftest compile6
## ;;;;  (compile 
## ;;;;   '(population p1 8 (merge 11
## ;;;;			    (population p2 2 (sample 2)))))
## ;;;;  '())
## ;;;;
## 
## 
## (define *epochs0*
##   '(population :name p1 :size 2 :epochs ((bottleneck 2 0 2)) :subtree (sample 2)))
## 
## (define *epochs1* 
##   '(population 
##     :name p0
##     :size 4 
##     :subtree (merge 2
## 		    (population :name p1 :size 2 :subtree (sample 2))
## 		    (population :name p2 :size 2 :subtree (sample 2)))))
## 
## (define *epochs2*
##   '(population :name p1
## 	       :size 6
## 	       :subtree (merge 1.5
## 			       (population :name p2 :size 2 :subtree (sample 2))
## 			       (population :name p3 :size 4
## 					   :subtree (merge 1
## 							   (population :name p4 :size 2
## 								       :subtree (sample 2))
## 							   (population :name p5 :size 2
## 								       :subtree (sample 2)))))))
## 
## 
## (deftest populationTimeFramesTest1
##   (populationTimeFrames 7 *epochs0*)
##   '((p1 (0 7))))
## 
## (deftest populationTimeFramesTest2
##   (populationTimeFrames 7 *epochs1*)
##   '((p0 (2 7)) (p1 (0 2)) (p2 (0 2))))
## 
## (deftest populationTimeFramesTest3
##   (populationTimeFrames 7 *epochs2*)
##   '((p1 (1.5 7)) (p2 (0 1.5)) (p3 (1 1.5)) (p4 (0 1)) (p5 (0 1))))
## 
## 
## 
## 
## (deftest getPopSizesTest
##   (getPopSizes (getPopulations 1 *epochs2*))
##   '(6 2 4 2 2))
## 
## (deftest findPopIndexTest
##   (findPopIndex 'p4 (getPopulations 1 (addIndexes *epochs2*)))
##   1)
## 
## 
## (deftest getMergeTimes0
##   (getMergeTimes *epochs0*)
##   '())
## 
## (deftest getMergeTimes1 
##   (getMergeTimes *epochs1*)
##   '(2))
## 
## (deftest getMergeTimes2
##   (getMergeTimes *epochs2*)
##   '(1.5 1))
## 
## 
## 
## (deftest getPopulations1
##   (getPopulations 1 *epochs0*)
##   '((population :name p1 :index #f :size 2 :epochs ((bottleneck 2 0 1) (bottleneck 2 0 2))
## 		:timeFrame (0 1))))
## 
## (deftest getPopulations2
##   (getPopulations 1 *epochs1*)
##   '((population :name p0 :index #f :size 4 :epochs ((bottleneck 4 2 1)) :timeFrame (2 1))
##     (population :name p1 :index #f :size 2 :epochs ((bottleneck 2 0 2)) :timeFrame (0 2))
##     (population :name p2 :index #f :size 2 :epochs ((bottleneck 2 0 2)) :timeFrame (0 2))))
## 
## (deftest getPopulations3
##   (getPopulations 1 *epochs2*)
##   '((population :name p1 :index #f :size 6 :epochs ((bottleneck 6 1.5 1)) :timeFrame (1.5 1))
##     (population :name p2 :index #f :size 2 :epochs ((bottleneck 2 0 1.5)) :timeFrame (0 1.5))
##     (population :name p3 :index #f :size 4 :epochs ((bottleneck 4 1 1.5)) :timeFrame (1 1.5))
##     (population :name p4 :index #f :size 2 :epochs ((bottleneck 2 0 1)) :timeFrame (0 1))
##     (population :name p5 :index #f :size 2 :epochs ((bottleneck 2 0 1)) :timeFrame (0 1))))
## 
## 
## 
## 
## (deftest collectMerges0
##   (collectMerges (addIndexes *epochs0*))
##   '())
## 
## (deftest collectMerges1
##   (collectMerges (addIndexes *epochs1*))
##   '((populationMerge 2 0 1)))
## 
## (deftest collectMerges2
##   (collectMerges (addIndexes *epochs2*))
##   '((populationMerge 1.5 0 1) (populationMerge 1 1 2)))
## 
## 
## 
## (deftest checkMergeTimes0
##   (checkMergeTimes 1 *epochs0* )
##   #t)
## 
## (deftest checkMergeTimes1
##   (checkMergeTimes 1 *epochs1*)
##   #t)
## 
## (deftest checkMergeTimes2
##   (checkMergeTimes 1 *epochs2* )
##   #t)
## 
## (deftest checkMergeTimesError1
##   (checkMergeTimes 1
##    '(population :name p1 :size 1 
## 		:subtree (merge 2
## 				(population :name p2 :size 2 :subtree (sample 2))
## 				(population :name p3 :size 4
## 					    :subtree (merge 2
## 							    (population :name p4 :size 2
## 									:subtree(sample 2))
## 							    (population :name p5 :size 2
## 									:subtree (sample 2)))))))
##   #f)
## 
## (deftest checkMergeTimesError2
##   (checkMergeTimes 1
##    '(population :name p1 :size 2
## 		:subtree (merge 1
## 				(population :name p2 :size 2 :subtree (sample 2))
## 				(population :name p3 :size 4 
## 					    :subtree (merge 2
## 							    (population :name p4 :size 2 
## 									:subtree (sample 2))
## 							    (population :name p5 :size 2
## 									:subtree (sample 2)))))))
##   #f)
## 
## 
## 
## (deftest noMergeMerge0
##   (noMergeMerge *epochs0*)
##   #t)
## 
## (deftest noMergeMerge1
##   (noMergeMerge *epochs1*)
##   #t)
## 
## (deftest noMergeMerge2
##   (noMergeMerge *epochs2*)
##   #t)
## 
## (deftest noMergeMergeError1
##   (noMergeMerge
##    '(population :name p1 :size 1
## 		:subtree(merge 2
## 			       (population :name p2 :size 2 :subtree (sample 2))
## 			       (merge 2
## 				      (population :name p4 :size 2 :subtree (sample 2))
## 				      (population :name p5 :size 2 :subtree (sample 2))))))
##   #f)
## 
## 
    
  def testMarkers(self):
    'Testing markers'
    # need to translate the following into python
## 
## (define tm (traitMarker 0.1 0.18 0.22))
## (define sm (snpMarker   0.2 0.10 0.90))
## (define mm (msMarker    0.3 0.4 4))
## 
## (display tm)(newline)
## (display sm)(newline)
## (display mm)(newline)
## 
## (display (position tm))(newline)
## (display (position sm))(newline)
## (display (position mm))(newline)
## 
## 
## (display (traitMarker? tm))(display (snpMarker? tm))(display (msMarker? tm))
## (newline)
## (display (traitMarker? sm))(display (snpMarker? sm))(display (msMarker? sm))
## (newline)
## (display (traitMarker? mm))(display (snpMarker? mm))(display (msMarker? mm))
## (newline)
    
  def testEpochs(self):
    'Testing epochs'
## 
## (display (bottleneck 0 0.1 1 2))(newline)
## (display (growth 0 10 1 2))(newline)
## (display (migration 0 1 .98 0 1))(newline)
## (display (populationMerge  .98 0 1))(newline)
## 
## (catch 'illegalEpoch
##        (lambda () (bottleneck 0 0.1 2 1))
##        (lambda (key . args) (display key)(newline)))
## 
## (catch 'illegalEpoch
##        (lambda () (bottleneck 0 0.1 1 2))
##        (lambda (key . args) (display key)(newline)))
## 
## 
## (display (bottleneck 0 0.1 1))(newline)
## (display (growth 0 10 1))(newline)
## (display (migration 0 1 .98 0 1))(newline)
## (display (populationMerge  .98 0 1))(newline)
## (display (populationMerge  .98 0 1 2))(newline)
    

if __name__ == '__main__':
  unittest.main()
