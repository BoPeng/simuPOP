Global functions
================


Function closeOutput
--------------------


.. function:: closeOutput(output="")

   Output files specified by ``'>'`` are closed immediately after they
   are written. Those specified by ``'>>'`` and ``'>>>'`` are closed
   by a simulator after ``Simulator.evolve()``. However, these files
   will be kept open if the operators are applied directly to a
   population using the operators' function form. In this case,
   function :func:`closeOutput` can be used to close a specific file
   *output*, and close all unclosed files if *output* is unspecified.
   An exception will be raised if *output* does not exist or it has
   already been closed.


Function describeEvolProcess
----------------------------


.. function:: describeEvolProcess(initOps=[], preOps=[], matingScheme=MatingScheme, postOps=[], finalOps=[], gen=-1, numRep=1)

   This function takes the same parameters as :meth:`Simulator.evolve`
   and output a description of how an evolutionary process will be
   executed. It is recommended that you call this function if you have
   any doubt how your simulation will proceed.


Function loadPopulation
-----------------------


.. function:: loadPopulation(file)

   load a population from a file saved by ``Population::save()``.


Function loadPedigree
---------------------


.. function:: loadPedigree(file, idField="ind_id", fatherField="father_id", motherField="mother_id", ploidy=2, loci=[], chromTypes=[], lociPos=[], chromNames=[], alleleNames=[], lociNames=[], subPopNames=[], infoFields=[])

   Load a pedigree from a file saved by operator
   :class:`PedigreeTagger` or function :meth:`Pedigree.save`. This
   file contains the ID of each offspring and their parent(s) and
   optionally sex ('M' or 'F'), affection status ('A' or 'U'), values
   of information fields and genotype at some loci. IDs of each
   individual and their parents are loaded to information fields
   *idField*, *fatherField* and *motherField*. Only numeric IDs are
   allowed, and individual IDs must be unique across all generations.

   Because this file does not contain generation information,
   generations to which offspring belong are determined by the parent-
   offspring relationships. Individuals without parents are assumed to
   be in the top-most ancestral generation. This is the case for
   individuals in the top-most ancestral generation if the file is
   saved by function ``Pedigree.save()``, and for individuals who only
   appear as another individual's parent, if the file is saved by
   operator :class:`PedigreeTagger`. The order at which offsprng is
   specified is not important because this function essentially
   creates a top-most ancestral generation using IDs without parents,
   and creates the next generation using offspring of these parents,
   and so on until all generations are recreated. That is to say, if
   you have a mixture of pedigrees with different generations, they
   will be lined up from the top most ancestral generation.

   If individual sex is not specified, sex of of parents are
   determined by their parental roles (father or mother) but the sex
   of individuals in the last generation can not be determined so they
   will all be males. If additional information fields are given,
   their names have to be specified using parameter *infoFields*. The
   rest of the columns are assued to be alleles, arranged *ploidy*
   consecutive columns for each locus. If paraemter *loci* is not
   specified, the number of loci is calculated by number of columns
   divided by *ploidy* (default to 2). All loci are assumed to be on
   one chromosome unless parameter *loci* is used to specified number
   of loci on each chromosome. Additional parameters such as *ploidy*,
   *chromTypes*, *lociPos*, *chromNames*, *alleleNames*, *lociNames*
   could be used to specified the genotype structured of the loaded
   pedigree. Please refer to class :class:`Population` for details
   about these parameters.


Function moduleInfo
-------------------


.. function:: moduleInfo()

   Return a dictionary with information regarding the currently loaded
   simuPOP module. This dictionary has the following keys:

   + ``revision:`` revision number.

   + ``version:`` simuPOP version string.

   + ``optimized:`` Is this module optimized (``True`` or ``False``).

   + ``alleleType:`` Allele type of the module (``short``, ``long`` or
     ``binary``).

   + ``maxAllele:`` the maximum allowed allele state, which is ``1``
     for binary modules, ``255`` for short modules and ``65535`` for
     long modules.

   + ``compiler:`` the compiler that compiles this module.

   + ``date:`` date on which this module is compiled.

   + ``python:`` version of python.

   + ``platform:`` platform of the module.

   + ``wordsize:`` size of word, can be either 32 or 64.

   + ``alleleBits:`` the number of bits used to store an allele

   + ``maxNumSubPop:`` maximum number of subpopulations.

   + ``maxIndex:`` maximum index size (limits population size * total
     number of marker).

   + ``debug:`` A dictionary with debugging codes as keys and the
     status of each debugging code (``True`` or ``False``) as their
     values.


Function getRNG
---------------


.. function:: getRNG()

   return the currently used random number generator


Function setRNG
---------------


.. function:: setRNG(name='', seed=0)

   Set random number generator. This function is obsolete but is provided
   for compatibility purposes. Please use setOptions instead


Function setOptions
-------------------


.. function:: setOptions(numThreads=-1, name=None, seed=0)

   First argument is to set number of thread in openMP. The number of
   threads can be be positive, integer (number of threads) or 0, which
   implies all available cores, or a number set by environmental
   variable ``OMP_NUM_THREADS``. Second and third argument is to set
   the type or seed of existing random number generator using
   RNG*name* with *seed*. If using openMP, it sets the type or seed of
   random number generator of each thread.


Function turnOnDebug
--------------------


.. function:: turnOnDebug(code="")

   Set debug code *code*. More than one code could be specified using
   a comma separated string. Name of available codes are available
   from ``moduleInfo()['debug'].keys()``.


Function turnOffDebug
---------------------


.. function:: turnOffDebug(code="DBG_ALL")

   Turn off debug code *code*. More than one code could be specified
   using a comma separated string. Default to turn off all debug
   codes.


