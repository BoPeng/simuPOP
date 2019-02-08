Miscellaneous operators
=======================


An operator that does nothing (operator :class:`NoneOp`)
--------------------------------------------------------

Operator :class:`NoneOp` does nothing when it is applied to a population. It
provides a placeholder when an operator is needed but no action is required.
Example :ref:`NoneOp <NoneOp>` demonstrates a typical usage of this operator

.. _NoneOp:

::

   if hasSelection:
       sel = MapSelector(loci=[0], fitness=[1, 0.99, 0.98])
   else:
       sel = NoneOp()
   #
   simu.evolve(
       preOps=[sel], # and other operators
       matingScheme=RandomMating(),
       gen=10
   )


dump the content of a population (operator :class:`Dumper`)
-----------------------------------------------------------

Operator :class:`Dumper` and its function form :func:`dump` has been used
extensively in this guide. They are prefect for demonstration and debugging
purposes because they display all properties of a population in a human readable
format. They are, however, rarely used in realistic settings because outputting
a large population to your terminal can be disastrous.

Even with modestly-sized populations, it is a good idea to dump only parts of
the population that you are interested. For example, you can use parameter
``genotype=False`` to stop outputting individual genotype, ``structure=False``
to stop outtputing genotypic and population structure information,
``loci=range(5)`` to output genotype only at the first five loci, ``max=N`` to
output only the first ``N`` individuals (default to ``100``), ``subPops=[(0,
0)]`` to output, for example, only the first virtual subpopulation in
subpopulation 0. Multiple virtual subpopulations are allowed and you can even
use ``subPops=[(ALL_AVAIL, 0)]`` to go through a specific virtual subpopulation
of all subpopulations. This operator by default only dump the present generation
but you can set ``ancGens`` to a list of generation numbers or ``ALL_AVAIL`` to
dump part or all ancestral generations. Finally, if there are more than 10
alleles, you can set the ``width`` at which each allele will be printed. The
following example (Example :ref:`Dumper <Dumper>`) presents a rather complicated
usage of this operator.

.. _Dumper:

**Example**: *dump the content of a population*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[10, 10], loci=[20, 30], infoFields='gen',
   ...     ancGen=-1)
   >>> sim.initSex(pop)
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> pop1 = pop.clone()
   >>> sim.initGenotype(pop, freq=[0]*20 + [0.1]*10)
   >>> pop.setIndInfo(1, 'gen')
   >>> sim.initGenotype(pop1, freq=[0]*50 + [0.1]*10)
   >>> pop1.setIndInfo(2, 'gen')
   >>> pop.push(pop1)
   >>> sim.dump(pop, width=3, loci=[5, 6, 30], subPops=([0, 0], [1, 1]),
   ...     max=10, structure=False)
   SubPopulation 0,0 (Male), 5 Individuals:
      2: MU  56 54 52 |  58 54 51 |  2
      3: MU  52 50 51 |  56 51 50 |  2
      4: MU  50 53 52 |  52 59 56 |  2
      5: MU  57 54 56 |  57 57 53 |  2
      6: MU  59 54 54 |  57 51 50 |  2
   SubPopulation 1,1 (Female), 7 Individuals:
     10: FU  54 53 57 |  59 59 59 |  2
     11: FU  55 59 51 |  59 51 58 |  2
     12: FU  55 58 58 |  57 54 58 |  2
     14: FU  53 57 52 |  51 54 58 |  2
     15: FU  51 58 59 |  54 52 54 |  2

   >>> # list all male individuals in all subpopulations
   >>> sim.dump(pop, width=3, loci=[5, 6, 30], subPops=[(sim.ALL_AVAIL, 0)],
   ...     max=10, structure=False)
   SubPopulation 0,0 (Male), 5 Individuals:
      2: MU  56 54 52 |  58 54 51 |  2
      3: MU  52 50 51 |  56 51 50 |  2
      4: MU  50 53 52 |  52 59 56 |  2
      5: MU  57 54 56 |  57 57 53 |  2
      6: MU  59 54 54 |  57 51 50 |  2
   SubPopulation 1,0 (Male), 3 Individuals:
     13: MU  55 52 53 |  57 56 52 |  2
     17: MU  55 51 51 |  57 55 51 |  2
     19: MU  56 54 53 |  58 58 56 |  2


   now exiting runScriptInteractively...

`Download Dumper.py <Dumper.py>`_


Save a population during evolution (operator :class:`SavePopulation`)
---------------------------------------------------------------------

Because it is usually not feasible to store all parental generations of an
evolving population, it is a common practise to save snapshots of a population
during an evolutionary process for further analysis. Operator
:class:`SavePopulation` is designed for this purpose. When it is applied to a
population, it will save the population to a file specified by parameter
``output``.

The tricky part is that populations at different generations need to be saved to
different filenames so the expression version of parameter ``output`` needs to
be used (see operator :class:`BaseOperator` for details). For example,
expression ``'snapshot_%d_%d.pop' % (rep, gen)`` is used in Example
:ref:`SavePopulation <SavePopulation>` to save population to files such as
``snapshot_5_20.pop`` during the evolution.

.. _SavePopulation:

**Example**: *Save snapshots of an evolving population*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(sim.Population(100, loci=2),
   ...     rep=5)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.2, 0.8])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=sim.SavePopulation(output="!'snapshot_%d_%d.pop' % (rep, gen)",
   ...             step = 10),
   ...     gen = 50
   ... )
   (50, 50, 50, 50, 50)

   now exiting runScriptInteractively...

`Download SavePopulation.py <SavePopulation.py>`_


Pause and resume an evolutionary process (operator :class:`Pause`) \*
---------------------------------------------------------------------

If you are presenting an evolutinary process in public, you might want to
temporarily stop the evolution so that your audience can have a better look at
intermediate results or figures. If you have an exceptionally long evolutionary
process, you might want to examine the status of the evolution process from time
to time. These can be done using a :class:`Pause` operator.

The :class:`Pause` operator can stop the evolution at specified generations, or
when you press a key. In the first case, you usually specify the generations to
Pause (e.g. :class:`Pause`\ (``step=1000``)) so that you can examine the status
of a simulation from time to time. In the second case, you can apply the
operator at each generation and Pause the simulation when you press a key (e.g.
:class:`Pause`\ (``stopOnKeyStroke=True``)). A specific key can be specified so
that you can use different keys to stop different populations, as shown in
Example :ref:`Pause <Pause>`.

.. _Pause:

**Example**: *Pause the evolution of a simulation*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(sim.Population(100), rep=10)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[sim.Pause(stopOnKeyStroke=str(x), reps=x) for x in range(10)],
   ...     gen = 100
   ... )
   (100, 100, 100, 100, 100, 100, 100, 100, 100, 100)

   now exiting runScriptInteractively...

`Download Pause.py <Pause.py>`_

When a simulation is Paused, you are given the options to resume evolution, stop
the evolution of the Paused population or all populations, or enter an
interactive Python shell to examine the status of a population, which will be
available in the Python shell as ``pop_X_Y`` where ``X`` and ``Y`` are
generation and replicate number of the population, respectively. The evolution
will resume after you exit the Python shell.


Measuring execution time of operators (operator :class:`TicToc`) \*
-------------------------------------------------------------------

The :class:`TicToc` operator can be used to measure the time between two events
during an evolutionary process. It outputs the elapsed time since the last time
it is called, and the overall time since the operator is created. It is very
flexible in that you can measure the time spent for mating in an evolutionary
cycle if you apply it before and after mating, and you can measure time spent
for several evolutionary cycles using generation applicability parameters such
as ``step`` and ``at``. The latter usage is demonstrated in Example :ref:`TicToc
<TicToc>`.

.. _TicToc:

**Example**: *Monitor the performance of operators*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(sim.Population(10000, loci=[100]*5), rep=2)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.1, 0.9])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.TicToc(step=50, reps=-1),
   ...     ],
   ...     gen = 101
   ... )
   Start stopwatch.
   Elapsed time: 5.00s	 Overall time: 5.00s
   Elapsed time: 4.00s	 Overall time: 9.00s
   (101, 101)

   now exiting runScriptInteractively...

`Download TicToc.py <TicToc.py>`_


