%module coaSim

//
//interface file for the python/swig binding of coasim
//
//
/////////////// VECTORS AND MATRICES FOR PARAMETER INPUT/////////////

%include "std_vector.i"
%include "std_string.i"
%include "stl.i"

namespace std
{
  %template()     vector< unsigned int >;
}



///////////////////////// EXCEPTION HANDLING ////////////////////////

%include exception.i

// exceptions will be passed to and correctly handled by python
%exception
{
  try
  {
    $function
  }
  catch(core::illegal_position e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(core::illegal_value e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(core::empty_interval e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(core::interval_out_of_range e)
  {
    SWIG_exception(SWIG_IndexError, e.what());
  }
  catch(core::non_pos_pop_size e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(core::out_of_sequence e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch(std::out_of_range e)
  {
    SWIG_exception(SWIG_IndexError, e.what());
  }
  catch(core::uninitialized_marker e)
  {
    SWIG_exception(SWIG_IndexError, e.what());
  }
  catch(...)
  {
    SWIG_exception(SWIG_RuntimeError, "Unknown runtime error happened.");
  }
}


////////////////////////// INCLUDE FILES ////////////////////////////

%{
#include "all_markers.hh"
#include "builder.hh"
#include "builder_events.hh"
#include "compile_options.hh"
#include "configuration.hh"
#include "descender.hh"
#include "dist_funcs.hh"
#include "interval.hh"
#include "marker.hh"
#include "epochs.hh"
#include "micro_satellite_marker.hh"
#include "monitor.hh"
#include "node.hh"
#include "retired_interval.hh"
#include "simulator.hh"
#include "snp_marker.hh"
#include "trait_marker.hh" 
%}

namespace std
{
  %template(MarkerVec) vector< core::Marker* >;
  %template(EpochVec)  vector< core::Event* >;
  %template(timeEvent) pair<double, core::Event*>;
}

////////////////////////// SWIG_INIT FUNCTION ///////////////////////
%inline
%{
  // required, but we have nothing to initialize yet.
  bool initialize()
  {
    return true;
  }  
%}

%init
%{
  initialize();
%}

////////////////////////// WRAP THESE CLASSES ///////////////////////
%ignore core::Configuration::Configuration(PopSizeItr, PopSizeItr, MakerItr, MarkerItrm EpochItr, EpochItr, double, double, double, double);

%include "all_markers.hh"
%include "builder.hh"
%include "marker.hh"
%include "configuration.hh"
%include "builder_events.hh"
%include "compile_options.hh"
%include "descender.hh"
%include "dist_funcs.hh"
%include "interval.hh"
%include "epochs.hh"
%include "micro_satellite_marker.hh"
%include "monitor.hh"
%include "node.hh"
%include "retired_interval.hh"
%include "simulator.hh"
%include "snp_marker.hh"
%include "trait_marker.hh"

//////////////////////// PYTHON UTILITY FUNCTION /////////////////////

%pythoncode %{

## Code that need to be translated to python
## function names have been translated 
##  by s/-\(.\)/\U\1/g in vim
##

## Simulate.scm ##########################################
## 
## (useModules (ice9 optargs)
## 	     (ice9 receive)
## 	     (ice9 format)
## 	     (ice9 prettyPrint)
## 	     (srfi srfi1))
## 
## (readSet! keywords 'prefix)
## 
## (defmacro unless (pred . body)
##   `(if ,pred
##        '()
##        (begin
## 	 ,@body)))
## 
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; dispatch framework ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## (define (defdispatchHelper name ast args)
##   (lambda (m)
##     (let ((tag (car  m))
## 	  (fun (cadr m)))
##       (if args
## 	  `(',tag (apply ,fun (list* ,@args (cdr ,ast))))
## 	  `(',tag (apply ,fun (cdr ,ast)))))))
## 
## (defmacro* defdispatch (name ast :key args . maps)
##   (let ((maps (if (equal? :args (car maps))
## 		  (listTail maps 2)
## 		  maps)))
##     `(define (,name ,@(if args args '()) ,ast)
##        (if (pair? ,ast)
## 	   (case (car ,ast)
## 	     ,@(map (defdispatchHelper name ast args) maps)
## 	     (else (throw 'unknownNode ,ast)))
## 	   (throw 'error "No atoms in this DSL" ,ast)))))
## 
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; getSampleSizes ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## (defdispatch getSampleSizes treeDescr
##   (merge      getSampleSizesMerge)
##   (sample     getSampleSizesSample)
##   (population getSampleSizesPop))
## 
## (define (getSampleSizesMerge time . subtrees) 
##   (apply append (map getSampleSizes subtrees)))
## 
## (define (getSampleSizesSample size) (list size))
## 
## (define* (getSampleSizesPop :key name index size epochs subtree)
##   (getSampleSizes subtree))
## 
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; getMergeTimes ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## (defdispatch getMergeTimes treeDescr
##   (merge      getMergeTimesMerge)
##   (sample     getMergeTimesSample)
##   (population getMergeTimesPop))
## 
## (define (getMergeTimesMerge time . subtrees)
##   (cons time (apply append (map getMergeTimes subtrees))))
## 
## (define (getMergeTimesSample size) '())
## 
## (define* (getMergeTimesPop :key name index size epochs subtree)
##   (getMergeTimes subtree))
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; getPopulations ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## (define (getPopKeyValue key pop)
##   (let ((res (findTail (lambda (e) 
## 			  (equal? e key))
## 			pop)))
##     (if res
## 	(second res)
## 	#f)))
## 
## (define (getPopSizes pops) 
##   (map (lambda (pop)
## 	 (getPopKeyValue :size pop))
##        pops))
## 
## (define (getPopNames pops) 
##   (map (lambda (pop)
## 	 (getPopKeyValue :name pop))
##        pops))
## 
## ;; Moej kode!!!!!!
## (define (getPopEpochs pops) 
##   (let ((tmp (map (lambda (pop)
## 		    (cons (getPopKeyValue :index pop)
## 			  (getPopKeyValue :epochs pop)))
## 		  pops)))
##     
##       (apply append 
## 	     (map (lambda (e)
## 		    (let ((index (first  e))
## 			  (epochs (cdr e)))
## 		      (map (lambda (epoch)
## 			     (list* (car epoch)
## 				    index
## 				    (cdr epoch)))
## 			   epochs)))
## 		  tmp))))
## 
## (define (findPopIndex name pops)
##   (getPopKeyValue :index
## 		     (find (lambda (pop) 
## 			     (equal? name (getPopKeyValue :name pop)))
## 			   pops)))
## 
## (define (findPopTimeFrame name pops)
##   (getPopKeyValue :timeFrame
## 		     (find (lambda (pop) 
## 			     (equal? name (getPopKeyValue :name pop)))
## 			   pops)))
## 
## (defdispatch getPopulationsHelper treeDescr
##   (merge      getPopulationsHelperMerge)
##   (sample     getPopulationsHelperSample)
##   (population getPopulationsHelperPop))
## 
## (define (getPopulationsHelperMerge time . subtrees)
##   (apply append (map getPopulationsHelper subtrees)))
## 
## (define (getPopulationsHelperSample size) '())
## 
## (define* (getPopulationsHelperPop :key name index size (epochs '()) subtree)
##   (cons `(population :name ,name :index ,index :size ,size :epochs ,epochs) 
## 	(getPopulationsHelper subtree)))
## 
## 
## (define (addTimeFrames pops timeFrames)  
##   (map (lambda (pop timeFrame)
## 	 (if (equal? (third pop) (car timeFrame))
## 	     (append pop (list :timeFrame) (cdr timeFrame))
## 	     (throw 'error "This should not happen 7962297043!!!")))
##        pops
##        timeFrames))
## 
## (define (checkForDuplicateNames pops)
##   (let ((names (filter (lambda (x) x)
## 		       (map (lambda (pop)
## 			      (getPopKeyValue :name pop))
## 			    pops))))
##     (unless (equal? (length names) (length (deleteDuplicates names)))
## 	    (throw 'syntaxError "Populations must have unique names" names))))
## 
## (define (inf<= . rest)
##   (let ((splitIndex (listIndex (lambda (x)
## 				   (equal? 1 x))
## 				 rest)))
##     (receive (xs infs) (if splitIndex
## 			   (splitAt rest splitIndex)
## 			   (values rest '()))
## 	     (if (every (lambda (x) (equal? 1 x)) infs)
## 		 (apply <= xs)
## 		 #f))))
## 
## (define (firstInsideSecond? tf1 tf2)
##   (and (inf<= (first tf2) (first  tf1) (second tf2))
##        (inf<= (first tf2) (second tf1) (second tf2))))
## 
## (define (firstLater=ThanSecond? tf1 tf2)
##   (and (inf<= (first tf2) (second tf2) (first  tf1) (second tf1))))
## 
## (define (crossBoundaries? tf1 tf2)
##   (not (or (firstInsideSecond? tf1 tf2)
## 	   (firstInsideSecond? tf2 tf1)
## 	   (firstLater=ThanSecond? tf1 tf2)
## 	   (firstLater=ThanSecond? tf2 tf1))))
## 
## (define (getEpochsTimeFrames pop)
##   (let ((epochs (getPopKeyValue :epochs pop)))
##     (map (lambda (epoch)
## 	   (takeRight epoch 2))
## 	 epochs)))
## 
## (define (checkForEpochsOverlap pop)
##   (let ((timeFrames (deleteDuplicates (getEpochsTimeFrames pop))))
##     (forEach (lambda (timeFrame)
## 		(if (find (lambda (timeFrame2)
## 			    (crossBoundaries? timeFrame timeFrame2))
## 			  timeFrames)
## 		    (throw 'syntaxError "Overlapping epochs exists" pop)))
## 	      timeFrames)))
## 
## (define (normaliseEpochs pops)
##   (map (lambda (pop)
## 	 (apply (lambda* (tag :key name index size epochs timeFrame)
## 			 (let* ((implicitGrowth
## 				 (map (lambda (epoch)
## 					(if (equal? 'beta (car epoch))
## 					    `(growth ,(second epoch) 
## 						     ,(first timeFrame) ,(second timeFrame))
## 					    epoch))
## 				      epochs))
## 				(implicitEnd 
## 				 (map (lambda (epoch)
## 					(cond ((and (equal? 'bottleneck (car epoch))
## 						    (= (length epoch) 3))
## 					       `(bottleneck ,(second epoch) ,(third epoch)
## 							    ,(second timeFrame)))
## 					      ((and (equal? 'growth (car epoch))
## 						    (= (length epoch) 3))
## 					       `(growth ,(second epoch) ,(third epoch)
## 							,(second timeFrame)))
## 					      (else epoch)))
## 				      implicitGrowth))
## 				(newEpochs implicitEnd))
## 			   `(population
## 			     :name ,name
## 			     :index ,index
## 			     :size ,size
## 			     :epochs ,(list* `(bottleneck ,size ,(first timeFrame) ,(second timeFrame))
## 					     (if newEpochs newEpochs '()))
## 			     :timeFrame ,timeFrame)))
## 		pop))
##        pops))
## 
## (define (getPopulations time treeDescr)
##   (resetIndex)
## 
##   (let* ((popsNoTime           (getPopulationsHelper treeDescr))
## 	 (timeFrames            (populationTimeFrames time treeDescr))
## 	 (popsNoDefaultEpochs (addTimeFrames popsNoTime timeFrames))
## 	 (pops                   (normaliseEpochs popsNoDefaultEpochs)))
## 
##     (checkForDuplicateNames pops)
##     (map checkForEpochsOverlap pops)
## 
##     pops))
## 
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; addIndexes ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## (define *index* 0)
## 
## (define (resetIndex) (set! *index* 0) *index*)
## 
## (define (incIndex)
##   (let ((i *index*))
##     (set! *index* (+ *index* 1))
##     i))
## 
## (defdispatch addIndexesHelper treeDescr
##   (merge      addIndexesHelperMerge)
##   (sample     addIndexesHelperSample)
##   (population addIndexesHelperPop))
## 
## (define (addIndexesHelperMerge time . subtrees)
##   (let* ((res   (map addIndexesHelper subtrees))
## 	 (index (caar res))
## 	 (trees (map (lambda (i/tree) (cdr i/tree)) res)))
##     (cons index `(merge ,time ,@trees))))
## 
## (define (addIndexesHelperSample size)
##   (cons (incIndex) `(sample ,size)))
## 
## (define* (addIndexesHelperPop :key name size (epochs '()) subtree)
##   (let ((i/tree (addIndexesHelper subtree)))
##   (cons (car i/tree)
## 	`(population :name ,name 
## 		     :index ,(car i/tree)
## 		     :size ,size
## 		     :epochs ,epochs
## 		     :subtree ,(cdr i/tree)))))
## 
## (define (addIndexes treeDescr)
##   (resetIndex)
##   (cdr (addIndexesHelper treeDescr)))
## 
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; collectMerges ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## (defdispatch collectMerges treeDescr
##   (merge      collectMergesMerge)
##   (sample     collectMergesSample)
##   (population collectMergesPop))
## 
## (define (collectMergesMerge time . subtrees)
##   (cons (list* 'populationMerge 
## 	       time
## 	       (map (lambda (pop)
## 		      (getPopulationIndex pop))
## 		    subtrees))
## 	(apply append (map collectMerges subtrees))))
## 
## (define (collectMergesSample size) '())
## 
## (define* (collectMergesPop :key name index size epochs subtree)
##   (collectMerges subtree))
## 
## (define (getPopulationIndex pop)
##   (apply (lambda* (tag :key index :allowOtherKeys) index) pop))
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; checkMergeTimes ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## ;; Checks if nested merges nodes are descending in time
## 
## (defdispatch checkMergeTimes treeDescr :args (currentTime)
##   (merge      checkMergeTimesMerge)
##   (sample     checkMergeTimesSample)
##   (population checkMergeTimesPop))
## 
## (define (checkMergeTimesMerge currentTime time . subtrees)
##   (if (or (= 1 currentTime) ;;infinity
## 	  (> currentTime time))
##       (every (lambda (x) x)
## 	     (map (lambda (subtree)
## 		    (checkMergeTimes time subtree))
## 		  subtrees))
##       #f))
## 
## (define (checkMergeTimesSample currentTime size) #t)
## 
## (define* (checkMergeTimesPop currentTime :key name index size epochs subtree)
##   (checkMergeTimes currentTime subtree))
## 
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; noMergeMerge ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## ;; Checks that a merge node cannot be a direct child of another merge node
## 
## (defdispatch noMergeMergeHelper treeDescr :args (mergep)
##   (merge      noMergeMergeMerge)
##   (sample     noMergeMergeSample)
##   (population noMergeMergePop))
## 
## (define (noMergeMergeMerge mergep time . subtrees)
##   (if mergep
##       #f
##       (every (lambda (x) x)
## 	     (map (lambda (subtree)
## 		    (noMergeMergeHelper #t subtree))
## 		  subtrees))))
## 
## (define (noMergeMergeSample mergep size) #t)
## 
## (define* (noMergeMergePop mergep :key name index size epochs subtree)
##   (noMergeMergeHelper #f subtree))
## 
## (define (noMergeMerge treeDescr)
##   (noMergeMergeHelper #f treeDescr))
## 
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; populationTimeFrames ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## ;; Checks that population sizes are consistent
## 
## (defdispatch populationTimeFrames tree :args (startTime)
##   (merge      populationTimeFramesMerge)
##   (sample     populationTimeFramesSample)
##   (population populationTimeFramesPop))
## 
## (define (populationTimeFramesMerge startTime time . subtrees)
##   (cons time (apply append (map (lambda (subtree) (populationTimeFrames time subtree)) subtrees))))
## 
## (define (populationTimeFramesSample startTime size)
##   (list 0))
## 
## (define* (populationTimeFramesPop startTime :key name index size epochs subtree)
##   (let ((res (populationTimeFrames #f subtree)))
##     (cons (list name (list (car res) startTime)) (cdr res))))
## 
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; systaxCheck ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## (define *list* '())
## (define (initListIterator lst) (set! *list* lst))
## (define (peek) (car *list*))
## (define (pop) 
##   (let ((e (peek))) 
##     (if (null? *list*)
## 	(throw 'error "Trying to pop the empty list")
## 	(set! *list* (cdr *list*))) 
##     e))
## 
## (define* (optionalArg test :key default)
##   (if (test (peek))
##       (pop)
##       default))
## 
## (define* (mandatoryArg test :key error)
##   (if (test (peek))
##       (pop)
##       (apply throw error)))
## 
## (define* (optionalKeywordArg key :key default)
##   (if (equal? key (peek))
##       (begin (pop) (pop))
##       default))
## 
## (defdispatch syntaxCheckAndNormalization treeDescr
##   (merge      syntaxCheckAndNormalizationMerge)
##   (sample     syntaxCheckAndNormalizationSample)
##   (population syntaxCheckAndNormalizationPop))
## 
## (define* (syntaxCheckAndNormalizationMerge . in)
##   (unless (>= (length in) 2)
## 	  (throw 'syntaxError "Merge takes 3 or more arguments" `(merge ,@in)))
##   (let ((time     (first in))
## 	(subtrees (cdr   in)))
##     (unless (number? time)
## 	    (throw 'systaxError "In merge: time has to be a number" time))
##     (unless (every pair? subtrees)
## 	    (throw 'systaxError "In merge: all subtrees are not nodes" `(merge ,@in)))
##     `(merge ,time ,@(map syntaxCheckAndNormalization subtrees))))
## 
## (define* (syntaxCheckAndNormalizationSample . in)
##   (unless (= (length in) 1)
## 	  (throw 'syntaxError "Sample takes one argument" `(sample ,@in)))
##   (let ((size (first  in)))
##     (unless (number? size)
## 	    (throw 'systaxError "In sample: size has to be a number" size))
##     `(sample ,size)))
## 
## (define (makeLocalEpochs epochs beta)
##   (if (and (pair? epochs)
## 	   (not (pair? (car epochs))))
##       (list* epochs (if beta
## 			`((beta .beta))
## 			'()))
##       (if beta
## 	  (list* `(beta ,beta) epochs)
## 	  epochs)))
## 
## (define* (syntaxCheckAndNormalizationPop . in)
##   (initListIterator in)
##   (let ((name (optionalArg symbol? :default #f))
## 	(size (mandatoryArg number?
## 			     :error `(syntaxError
## 				      "In population: size must be a number" ,`(population ,@in))))
## 	(beta   (optionalKeywordArg :beta :default #f))
## 	(epochs (optionalKeywordArg :epochs :default '()))
## 	(subtree (mandatoryArg pair?
## 				:error `(syntaxError
## 					 "In population: error in subtree definition"
## 					 ,`(population ,@in)))))
##     `(population :name ,name 
## 		 :size ,size 
## 		 :epochs ,(makeLocalEpochs epochs beta)
## 		 :subtree ,(syntaxCheckAndNormalization subtree))))
## 
## 
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;; The compiler ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## 
## (define (checkThatPopsHaveNames migrations pops)
##   (let ((mNames (deleteDuplicates 
## 		  (apply append 
## 			 (map (lambda (m) (list (second m) (third m)))
## 			      migrations))))
## 	(pNames (deleteDuplicates (getPopNames pops))))
##     (unless (or (null? mNames)
## 		(null? (lsetDifference equal? mNames pNames)))
## 	    (throw 'syntaxError "Populations have to have names in order to specify migrations"
## 		   (lsetDifference equal? mNames pNames)))))
## 
## (define (normalizeMigrations migrations pops)
##   (map (lambda (migration) 
## 	 (if (= 4 (length migration))
## 	     (let ((tfFrom (findPopTimeFrame (second migration) pops))
## 		   (tfTo   (findPopTimeFrame (third  migration) pops)))
## 	       (unless (and (not (firstLater=ThanSecond? tfFrom tfTo))
## 			    (not (firstLater=ThanSecond? tfTo   tfFrom)))
## 		       (throw 'syntaxError "Populations not valid in the same time interval" migration))
## 	       (append migration (list (max (first  tfFrom) (first  tfTo))
## 				       (min (second tfFrom) (second tfTo)))))
## 	     migration))
##        migrations))
## 
## (define (addMigrations migrations pops)
##   (checkThatPopsHaveNames migrations pops) 
##   (map (lambda (migration)
## 	 (list* 'migration 
## 		(findPopIndex (second migration) pops)
## 		(findPopIndex (third migration) pops)
## 		(listTail migration 3)))
##        (normalizeMigrations migrations pops)))
## 
## (define (compile tree migrations)
##   (let ((tree2 (addIndexes (syntaxCheckAndNormalization tree))))
##     
##     (unless (checkMergeTimes 1 tree2)
## 	    (throw 'astCheck "MergeTimes not descending"))
##     (unless (noMergeMerge tree2)
## 	    (throw 'astCheck "A mergeNode cannot be followed by a mergeNode"))
## 
##     (let* ((pops         (getPopulations 1 tree2))
## 	   (sampleSizes (getSampleSizes tree2))
## 	   (popSizes    (getPopSizes    pops))
## 	   (merges       (collectMerges   tree2)))
## 
##       `(:sampleSizes ,sampleSizes
## 		      :epochs ,(append merges
## 				       (getPopEpochs pops)
## 				       (addMigrations migrations pops))))))
## 
## 
## 
## 
## (define (simulate ms program . args)
##    (letKeywords args #f ((rho   0)
##                           (gamma 0)
##                           (Q     0)
##                           (beta  0)
## 			  (migration '())
##                           (coalescenceCallback '())
##                           (recombinationCallback '())
##                           (geneconversionCallback '())
##                           (keepEmptyIntervals #f)
##                           (randomSeed 0))
##      (if (number? program)
##          ;; simple simulation run...
##          (cSimulate ms (list program)
##                       (list rho gamma Q beta)
##                       (list coalescenceCallback
##                             recombinationCallback
##                             geneconversionCallback)
##                       '() ;; no epochs
##                       keepEmptyIntervals
##                       randomSeed)
## 	 ;; otherwise, compile the program to get epochs and samples
## 	 (let ((popStructure (compile program migration)))
## 	   (letKeywords popStructure #f ((sampleSizes '()) (epochs '()))
## 			 (let ((realEpochs (map (lambda (e) (primitiveEval e)) epochs)))
## 			   (cSimulate ms sampleSizes
## 				       (list rho gamma Q beta)
## 				       (list coalescenceCallback
## 					     recombinationCallback
## 					     geneconversionCallback)
## 				       realEpochs
## 				       keepEmptyIntervals
## 				       randomSeed)))))))
## 
## 
## 
%}
