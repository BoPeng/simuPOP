Conditional operators
=====================


Conditional operator (operator :class:`IfElse`) \*
--------------------------------------------------

Operator :class:`IfElse` provides a simple way to conditionally apply an
operator. The condition can be a fixed condition, a expression (a string) that
will be evaluated in a population's local namespace or a user-defined function
when it is applied to the population.

The first case is used to control the execution of certain operators depending
on user input. For example, Example :ref:`IfElseFixed <IfElseFixed>` determines
whether or not some outputs should be given depending on a variable ``verbose``.
Note that the applicability of the conditional operators are determined by the
:class:`IfElse` operator and individual opearators. That is to say, the
parameters ``begin``, ``step``, ``end``, ``at``, and ``reps`` of operators in
``ifOps`` and ``elseOps`` are only honored when operator :class:`IfElse` is
applied.

.. _IfElseFixed:

**Example**: *A conditional opeartor with fixed condition*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=1000, loci=1)
   >>> verbose = True
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=sim.IfElse(verbose,
   ...         ifOps=[
   ...             sim.Stat(alleleFreq=0),
   ...             sim.PyEval(r"'Gen: %3d, allele freq: %.3f\n' % (gen, alleleFreq[0][1])",
   ...                 step=5)
   ...         ],
   ...         begin=10),
   ...     gen = 30
   ... )
   Gen:  10, allele freq: 0.483
   Gen:  15, allele freq: 0.455
   Gen:  20, allele freq: 0.481
   Gen:  25, allele freq: 0.481
   30

   now exiting runScriptInteractively...

`Download IfElseFixed.py <IfElseFixed.py>`_

When a string is specified, it will be considered as an expression and be
evaluated in a population's namespace. The return value will be used to
determine if an operator should be executed. For example, you can re-introduce a
mutant if it gets lost in the population, output a warning when certain
condition is met, or record the occurance of certain events in a population. For
example, Example :ref:`IfElse <IfElse>` records the number of generations the
frequency of an allele goes below 0.4 and beyong 0.6 before it gets lost or
fixed in the population. Note that a list of else-operators can also be executed
when the condition is not met.

.. _IfElse:

**Example**: *A conditional opeartor with dynamic condition*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(
   ...     sim.Population(size=1000, loci=1),
   ...     rep=4)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...         sim.PyExec('below40, above60 = 0, 0')
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.IfElse('alleleFreq[0][1] < 0.4',
   ...             sim.PyExec('below40 += 1')),
   ...         sim.IfElse('alleleFreq[0][1] > 0.6',
   ...             sim.PyExec('above60 += 1')),
   ...         sim.IfElse('len(alleleFreq[0]) == 1',
   ...             sim.PyExec('stoppedAt = gen')),
   ...         sim.TerminateIf('len(alleleFreq[0]) == 1')
   ...     ]
   ... )
   (892, 1898, 4001, 2946)
   >>> for pop in simu.populations():
   ...     print('Overall: %4d, below 40%%: %4d, above 60%%: %4d' % \
   ...         (pop.dvars().stoppedAt, pop.dvars().below40, pop.dvars().above60))
   ... 
   Overall:  891, below 40%:   20, above 60%:  515
   Overall: 1897, below 40%: 1039, above 60%:   51
   Overall: 4000, below 40%: 2878, above 60%:    0
   Overall: 2945, below 40%:  198, above 60%: 1731

   now exiting runScriptInteractively...

`Download IfElse.py <IfElse.py>`_

In the last case, a user-defined function can be specified. This function should
accept parameter ``pop`` when the operator is applied to a population, and one
or more parameters ``pop``, ``off``, ``dad`` and ``mom`` when it is applied
during-mating. The later could be used to apply different during-mating
operators for different types of parents or offspring. For example, Example
:ref:`pedigreeMatingAgeStructured <pedigreeMatingAgeStructured>` in Chapter 6
uses a :class:`CloneGenoTransmitter` when only one parent is available (when
parameter ``mom`` is ``None``), and a :class:`MendelianGenoTransmitter` when two
parents are available.


Conditionally terminate an evolutionary process (operator :class:`TerminateIf`)
-------------------------------------------------------------------------------

Operator :class:`TerminateIf` has been described and used in several examples
such as Example :ref:`simuGen <simuGen>`, :ref:`expression <expression>` and
:ref:`IfElse <IfElse>`. This operator accept an Python expression and terminate
the evolution of the population being applied if the expression is evaluated to
be ``True``. This operator is well suited for situations where the number of
generations to evolve cannot be determined in advance.

If a :class:`TerminateIf` operator is applied to the offspring generation, the
evolutionary cycle is considered to be completed. If the evolution is terminated
before mating, the evolutionary cycle is condered to be incomplete. Such a
difference can be important if the number of generations that have been involved
is important for your analysis.

A less-known feature of operator :class:`TerminateIf` is its ability to
terminate the evolution of all replicates, using parameter ``stopAll=True``. For
example, Example :ref:`TerminateIf <TerminateIf>` terminates the evolution of
all populations when one of the populations gets fixed. The return value of
``simu.evolve`` shows that some populations have evolved one generation less
than the population being fixed.

.. _TerminateIf:

**Example**: *Terminate the evolution of all populations in a simulator*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(
   ...     sim.Population(size=100, loci=1),
   ...     rep=10)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.TerminateIf('len(alleleFreq[0]) == 1', stopAll=True)
   ...     ]
   ... )
   (88, 88, 88, 88, 87, 87, 87, 87, 87, 87)
   >>> 

   now exiting runScriptInteractively...

`Download TerminateIf.py <TerminateIf.py>`_


Conditional removal of individuals (operator :class:`DiscardIf`)
----------------------------------------------------------------

Operator :class:`DiscardIf` accepts a fixed condition or probability, or a
condition or a Python function that returns either ``True``/``False`` or a
probability to remove an individual. When it is applied during mating, it will
evaluate the condition or call the function for each offspring, and discard the
offspring if the return value of the expression or function is True, or remove
at a probability if the return value is a number between 0 and 1. The python
expression accepts information fields as variables so operator
:class:`DiscardIf`\ (``'age > 80'``) will discard all individuals with age > 80,
and :class:`DiscardIf`\ (``'1-fitness'``) will remove individuals according to 1
minus their fitness. Optionally, the offspring itself can be used in the
expression if parameter exposeInd is used to set the variable name of the
offspring.

Alternatively, a Python function can be passed to this operator. This function
should be defined with parameters ``pop``, ``off``, ``mom``, ``dad`` or names of
information fields. For example, :class:`DiscardIf`\ (``lambda age: age > 80``)
will remove individuals with age > 80.

A constant expression is also allowed in this operator. A fixed condition or
number is acceptable so :class:`DiscardIf`\ (``0.1``) will randomly remove 10%
of all individuals. Although it does not make sense to use :class:`DiscardIf`\
(``True``) because all offspring will be discarded, it is quite useful to use
this operator in the context of :class:`DiscardIf`\ (``True, subPops=[(0, 0)]``)
to remove all individuals in a virtual subpopulation. If virtual subpopulation
``(0, 0)`` is defined as all individuals with age > 80, the last method achieves
the same effect as the first two methods.

Example :ref:`DiscardIf <DiscardIf>` demonstrates an interesting application of
this operator. This example evolves a population for one generation. Instead of
keeping all offspring, it keeps only 500 affected and 500 unaffected offspring.
This is achieved by defining virtual subpopulations by affection status and
range, and discard the first 500 offspring if they are unaffected, and the last
500 offspring if they are affected.

.. _DiscardIf:

**Example**: *Use operator DiscardIf to generate case control samples*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=500, loci=1)
   >>> pop.setVirtualSplitter(sim.ProductSplitter([
   ...     sim.AffectionSplitter(),
   ...     sim.RangeSplitter([[0,500], [500, 1000]]),
   ...     ])
   ... )
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...     ],
   ...     matingScheme=sim.RandomMating(
   ...         ops=[
   ...             sim.MendelianGenoTransmitter(),
   ...             sim.MaPenetrance(loci=0, penetrance=[0, 0.01, 0.1]),
   ...             sim.DiscardIf(True, subPops=[
   ...                 (0, 'Unaffected, Range [0, 500)'),
   ...                 (0, 'Affected, Range [500, 1000)')])
   ...         ],
   ...         subPopSize=1000,
   ...     ),
   ...     gen = 1
   ... )
   1
   >>> sim.stat(pop, numOfAffected=True)
   >>> print(pop.dvars().numOfAffected, pop.dvars().numOfUnaffected)
   500 500

   now exiting runScriptInteractively...

`Download DiscardIf.py <DiscardIf.py>`_


