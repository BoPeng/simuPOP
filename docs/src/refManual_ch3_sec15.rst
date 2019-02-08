Function form of operators
==========================


Function acgtMutate
-------------------


.. function:: acgtMutate(pop, *args, **kwargs)

   Function form of operator :class:`AcgtMutator`


Function contextMutate
----------------------


.. function:: contextMutate(pop, *args, **kwargs)

   Function form of operator :class:`ContextMutator`


Function discardIf
------------------


.. function:: discardIf(pop, *args, **kwargs)

   Apply operator :class:`DiscardIf` to population *pop* to remove individuals according
   to an expression or a Python function.


Function dump
-------------


.. function:: dump(pop, *args, **kwargs)

   Apply operator :class:`Dumper` to population *pop*.


Function infoEval
-----------------


.. function:: infoEval(pop, *args, **kwargs)

   Evaluate *expr* for each individual, using information fields as variables.
   Please refer to operator :class:`InfoEval` for details.


Function infoExec
-----------------


.. function:: infoExec(pop, *args, **kwargs)

   Execute *stmts* for each individual, using information fields as variables.
   Please refer to operator :class:`InfoExec` for details.


Function initGenotype
---------------------


.. function:: initGenotype(pop, *args, **kwargs)

   Apply operator :class:`InitGenotype` to population *pop*.


Function initInfo
-----------------


.. function:: initInfo(pop, *args, **kwargs)

   Apply operator :class:`InitInfo` to population *pop*.


Function initSex
----------------


.. function:: initSex(pop, *args, **kwargs)

   Apply operator :class:`InitSex` to population *pop*.


Function kAlleleMutate
----------------------


.. function:: kAlleleMutate(pop, *args, **kwargs)

   Function form of operator :class:`KAlleleMutator`


Function maPenetrance
---------------------


.. function:: maPenetrance(pop, loci, penetrance, wildtype=0, ancGens=True, *args, **kwargs)

   Apply opertor :class:`MaPenetrance` to population *pop*. Unlike the
   operator form of this operator that only handles the current generation,
   this function by default assign affection status to all generations.


Function mapPenetrance
----------------------


.. function:: mapPenetrance(pop, loci, penetrance, ancGens=True, *args, **kwargs)

   Apply opertor :class:`MapPenetrance` to population *pop*. Unlike the
   operator form of this operator that only handles the current generation,
   this function by default assign affection status to all generations.


Function matrixMutate
---------------------


.. function:: matrixMutate(pop, *args, **kwargs)

   Function form of operator :class:`MatrixMutator`


Function mergeSubPops
---------------------


.. function:: mergeSubPops(pop, *args, **kwargs)

   Merge subpopulations *subPops* of population *pop* into a single
   subpopulation. Please refer to the operator form of this funciton
   (:class:`MergeSubPops`) for details


Function migrate
----------------


.. function:: migrate(pop, *args, **kwargs)

   Function form of operator :class:`Migrator`.


Function backwardMigrate
------------------------


.. function:: backwardMigrate(pop, *args, **kwargs)

   Function form of operator :class:`BackwardMigrator`.


Function mixedMutate
--------------------


.. function:: mixedMutate(pop, *args, **kwargs)

   Function form of operator :class:`MixedMutator`


Function mlPenetrance
---------------------


.. function:: mlPenetrance(pop, ops, mode, ancGens=True, *args, **kwargs)

   Apply opertor :class:`MapPenetrance` to population *pop*. Unlike the
   operator form of this operator that only handles the current generation,
   this function by default assign affection status to all generations.


Function pointMutate
--------------------


.. function:: pointMutate(pop, *args, **kwargs)

   Function form of operator :class:`PointMutator`


Function pyEval
---------------


.. function:: pyEval(pop, *args, **kwargs)

   Evaluate statements *stmts* (optional) and expression *expr* in
   population *pop*\ 's local namespace and return the result of *expr*.
   If *exposePop* is given, population *pop* will be exposed in its local
   namespace as a variable with a name specified by *exposePop*. Unlike its
   operator counterpart, this function returns the result of *expr* rather
   than writting it to an output.


Function pyExec
---------------


.. function:: pyExec(pop, *args, **kwargs)

   Execute *stmts* in population *pop*\ 's local namespace.


Function pyMutate
-----------------


.. function:: pyMutate(pop, *args, **kwargs)

   Function form of operator :class:`PyMutator`


Function pyPenetrance
---------------------


.. function:: pyPenetrance(pop, func, loci=[], ancGens=True, *args, **kwargs)

   Apply opertor :class:`PyPenetrance` to population *pop*. Unlike the
   operator form of this operator that only handles the current generation,
   this function by default assign affection status to all generations.


Function pyMlPenetrance
-----------------------


.. function:: pyMlPenetrance(pop, func, mode, loci=[], ancGens=True, *args, **kwargs)

   Apply opertor :class:`PyMlPenetrance` to population *pop*. Unlike the
   operator form of this operator that only handles the current generation,
   this function by default assign affection status to all generations.


Function pyQuanTrait
--------------------


.. function:: pyQuanTrait(pop, func, loci=[], ancGens=True, *args, **kwargs)

   Apply opertor :class:`PyQuanTrait` to population *pop*. Unlike the
   operator form of this operator that only handles the current generation,
   this function by default assign affection status to all generations.


Function resizeSubPops
----------------------


.. function:: resizeSubPops(pop, *args, **kwargs)

   Resize subpopulations *subPops* of population *pop* into new sizes
   *size*. Individuals will be added or removed accordingly. Please refer to
   the operator form of this funciton (:class:`ResizeSubPops`) for details


Function snpMutate
------------------


.. function:: snpMutate(pop, *args, **kwargs)

   Function form of operator :class:`SNPMutator`


Function splitSubPops
---------------------


.. function:: splitSubPops(pop, *args, **kwargs)

   Split subpopulations (*subPops*) of population *pop* according to either
   *sizes* or *proportions* of the resulting subpopulations, or an information
   field. Please refer to the operator form of this function (``splitSubPop``)
   for details.


Function stat
-------------


.. function:: stat(pop, *args, **kwargs)

   Apply operator :class:`Stat` with specified parameters to population *pop*.
   Resulting statistics could be accessed from the local namespace of ``pop``
   using functions ``pop.vars()`` or ``pop.dvars()``


Function stepwiseMutate
-----------------------


.. function:: stepwiseMutate(pop, *args, **kwargs)

   Function form of operator :class:`StepwiseMutator`


Function tagID
--------------


.. function:: tagID(pop, reset=False, *args, **kwargs)

   Apply operator :class:`IdTagger` to population *pop* to assign a unique ID
   to all individuals in the population. Individuals ID will starts from a
   system wide index. You can reset this start ID using parameter ``reset``
   which can be ``True`` (reset to 1) or a non-negative number (start from
   this number).


