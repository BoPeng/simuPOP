Debug-related functions and operators \*
========================================

Debug information can be useful when something looks suspicious. By turnning on
certain debug code, simuPOP will print out some internal information before and
during evolution. Functions :func:`turnOnDebug`\ (``code``) and
:func:`turnOffDebug`\ (``code``) could be used to turn on and off some debug
information.

For example, the following code might crash simuPOP:

::

   >>> Population(1, loci=[100]).individual(0).genotype()

It is unclear why this simple command causes us trouble, instead of outputting
the genotype of the only Individual of this population. However, the reason is
clear if you turn on debug information:

**Example**: *Turn on/off debug information*

::

   >>> turnOnDebug(DBG_POPULATION)
   >>> Population(1, loci=100).individual(0).genotype()
   Constructor of population is called
   Destructor of population is called 
   Segmentation fault (core dumped)

:class:`Population`\ (``1, loci=[100]``) creates a temporary object that is
destroyed right after the execution of the command. When Python tries to display
the genotype, it will refer to an invalid location. The correct method to print
the genotype is to create a persistent population object:

::

   >>> pop = Population(1, loci=[100])
   >>> pop.individual(0).genotype()

Another useful debug code is ``DBG_WARNING``. When this code is set, it will
output warning messages for some common misuse of simuPOP. For example, it will
warn you that population object returned by function
:meth:`Simulator.population`\ () is a temporary object that will become invalid
once a simulator is changed. If you are new to simuPOP, it is recommended that
you use   ::

   import simuOpt
   simuOpt.setOptions(optimized=False, debug='DBG_WARNING')

when you develop your script.

Besides functions :func:`turnOnDebug`\ ``(code)``\ and :func:`turnOffDebug`\
(``code``), you can set environmental variable ``SIMUDEBUG=code`` where ``code``
is a comma separated debug codes.\ ````A list of valid debug code could be found
in function :func:`moduleInfo``\ `['debug']`\ `. Note that debug information is
only available in standard (non-optimized) modules.

The amount of output can be overwhelming in some cases which makes it necessary
to limit the debug information to certain generations, or triggered by certain
conditions. In addition, debugging information may interfere with your regular
output so you may want to direct such output to another destination, such as a
dedicated file.

Example :ref:`debug <debug>` demonstrates how to turn on debug information
conditionally and turn it off afterwards, using operator :class:`PyOperator`. It
also demonstrates how to redirect debug output to a file but redefining system
standard error output. Note that "``is None``" is used to make sure the lamdba
functions return ``True`` so that the evolutionary process can continue after
the python operator.

.. _debug:

**Example**: *Turn on and off debug information during evolution.*

::

   >>> import simuPOP as sim
   >>> # redirect system stderr
   >>> import sys
   >>> debugOutput = open('debug.txt', 'w')
   >>> old_stderr = sys.stderr
   >>> sys.stderr = debugOutput
   >>> # start simulation
   >>> simu = sim.Simulator(sim.Population(100, loci=1), rep=5)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.1, 0.9])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.IfElse('alleleNum[0][0] == 0',
   ...             ifOps=[
   ...                 # the is None part makes the function return True
   ...                 sim.PyOperator(lambda : sim.turnOnDebug("DBG_MUTATOR") is None),
   ...                 sim.PointMutator(loci=0, allele=0, inds=0),
   ...             ],
   ...             elseOps=sim.PyOperator(lambda : sim.turnOffDebug("DBG_MUTATOR") is None)),
   ...     ],
   ...     gen = 100
   ... )
   (100, 100, 100, 100, 100)
   >>> # replace standard stdandard error
   >>> sys.stderr = old_stderr
   >>> debugOutput.close()
   >>> print(''.join(open('debug.txt').readlines()[:5]))
   Mutate locus 0 at ploidy 0 to allele 0 at generation 12
   Mutate locus 0 at ploidy 0 to allele 0 at generation 13
   Mutate locus 0 at ploidy 0 to allele 0 at generation 15
   Mutate locus 0 at ploidy 0 to allele 0 at generation 16
   Mutate locus 0 at ploidy 0 to allele 0 at generation 21


   now exiting runScriptInteractively...

`Download debug.py <debug.py>`_


