Simulator
=========

A simuPOP simulator evolves one or more copies of a population forward in time,
subject to various operators. Although a population could evolve by itself using
function :meth:`Population.evolve`, a simulator with one replicate is actually
used.


Add, access and remove populations from a simulator
---------------------------------------------------

A simulator could be created by one or more replicates of a list of populations.
For example, you could create a simulator from five replicates of a population
using   ::

   Simulator(pop, rep=5)

or from a list of populations using   ::

   Simulator([pop, pop1, pop2])

. ``pop``, ``pop1`` and ``pop2`` do not have to have the same genotypic
structure. In order to avoid duplication of potentially large populations, a
population is by default *stolen* after it is used to create a simulator. If you
would like to keep the populations, you could set parameter ``stealPops`` to
``False`` so that the populations will be copied to the simulator. Populations
in a simulator can be added or removed using functions :meth:`Simulator.add`\ ()
and :meth:`Simulator.extract`\ (``idx``).

When a simulator is created, you can access populations in this simulator using
function :meth:`Simulator.population`\ (``idx``) or iterate through all
populations using function :meth:`Simulator.populations`\ (). These functions
return references to the populations so that you can access populations.
Modifying these references will change the corresponding populations within the
simulator. The references will become invalid once the simulator object is
destoryed.

Example :ref:`Simulator <Simulator>` demonstrates different ways to create a
simulator and how to access populations within it.

.. _Simulator:

**Example**: *Create a simulator and access populations*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(100, loci=10)
   >>> # five copies of the same population
   >>> simu = sim.Simulator(pop, rep=5)
   >>> simu.numRep()
   5
   >>> # evolve for ten generations and save the populations
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.3, 0.7])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     finalOps=sim.SavePopulation('!"pop%d.pop"%rep'),
   ...     gen=10
   ... )
   (10, 10, 10, 10, 10)
   >>> # load the population and create another Simulator
   >>> simu = sim.Simulator([sim.loadPopulation('pop%d.pop' % x) for x in range(5)])
   >>> # continue to evolve
   >>> simu.evolve(
   ...     matingScheme=sim.RandomMating(),
   ...     gen=10
   ... )
   (10, 10, 10, 10, 10)
   >>> # print out allele frequency
   >>> for pop in simu.populations():
   ...     sim.stat(pop, alleleFreq=0)
   ...     print('%.2f' % pop.dvars().alleleFreq[0][0])
   ... 
   0.36
   0.30
   0.28
   0.01
   0.11
   >>> # get a population
   >>> pop = simu.extract(0)
   >>> simu.numRep()
   4

   now exiting runScriptInteractively...

`Download Simulator.py <Simulator.py>`_


Number of generations to evolve
-------------------------------

A simulator usually evolves a specific number of generations according to
parameter ``gen`` of the ``evolve`` function. A generation number is used to
track the number of generations a simulator has evolved. Because a new
population has generation number 0, a population would be at the beginning of
generation :math:`n` after it evolves :math:`n` generations. The generation
number would increase if the simulator continues to evolve. During evoluting,
variables ``rep`` (replicate number) and ``gen`` (current generation number) are
set to each population's local namespace.

It is not always possible to know in advance the number of generations to
evolve. For example, you may want to evolve a population until a specific allele
gets fixed or lost in the population. In this case, you can let the simulator
run indefinitely (do not set the ``gen`` parameter) and depend on a *terminator
*to terminate the evolution of a population. The easiest method to do this is to
use population variables to track the status of a population, and use a
:class:`TerminateIf` operator to terminate the evolution according to the value
of an expression. Example :ref:`simuGen <simuGen>` demonstrates the use of such
a terminator, which terminates the evolution of a population if allele 0 at
locus 5 is fixed or lost. It also shows the application of an interesting
operator :class:`IfElse`, which applies an operator, in this case
:class:`PyEval`, only when an expression returns ``True``. Note that this
example calls the ``simulator.evolve`` function twice. The first call does not
specify a mating scheme so a default empty mating scheme (:class:`MatingScheme`)
that does not transmit genotype is used. Populations start from the beginning of
the fifth generation when the second ``simulator.evole`` function is called.

The generation number is stored in each Population using population variable
``gen``.You can access these numbers from a simulator using function
:meth:`Simulator.dvars`\ (``idx``) or from a population using function
:meth:`Population.dvars`\ (). If needed, **you can reset generation numbers by
changing these variables.**

.. _simuGen:

**Example**: *Generation number of a simulator*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(sim.Population(50, loci=[10], ploidy=1),
   ...     rep=3)
   >>> simu.evolve(gen = 5)
   (5, 5, 5)
   >>> simu.dvars(0).gen
   5
   >>> simu.evolve(
   ...     initOps=[sim.InitGenotype(freq=[0.5, 0.5])],
   ...     matingScheme=sim.RandomSelection(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=5),
   ...         sim.IfElse('alleleNum[5][0] == 0',
   ...             sim.PyEval(r"'Allele 0 is lost in rep %d at gen %d\n' % (rep, gen)")),
   ...         sim.IfElse('alleleNum[5][0] == 50',
   ...             sim.PyEval(r"'Allele 0 is fixed in rep %d at gen %d\n' % (rep, gen)")),
   ...         sim.TerminateIf('len(alleleNum[5]) == 1'),
   ...     ],
   ... )
   Allele 0 is fixed in rep 2 at gen 29
   Allele 0 is fixed in rep 1 at gen 74
   Allele 0 is lost in rep 0 at gen 120
   (116, 70, 25)
   >>> [simu.dvars(x).gen for x in range(3)]
   [121, 75, 30]

   now exiting runScriptInteractively...

`Download simuGen.py <simuGen.py>`_


Evolve populations in a simulator
---------------------------------

There are a number of rules about when and how operators are applied during the
evolution of a population. In summary, in the order at which operators are
processed and applied,

* Operators specified in parameter ``initOps`` of function
  :meth:`Simulator.evolve` will be applied to the initial population before
  evolution, subject to replicate applicability restraint specified by parameter
  ``reps``.

* Operators specified in parameter ``preOps`` of function
  :meth:`Simulator.evolve` will be applied to the parental population at each
  generation, subject to replicate and generation applicability restraint
  specified by parameters ``begin``, ``end``, ``step``, ``at``, and ``reps``.

* During-mating operators specified in the ``ops`` parameter of a mating scheme
  will be called during mating to transmit genotype (and possibly information
  fields etc) from parental to offspring, subject to replicate and generation
  applicability restraint specified by parameters ``begin``, ``end``, ``step``,
  ``at``, and ``reps``.

* Operators specified in parameter postOps of function :meth:`Simulator.evolve`
  will be applied to the offspring population at each generation, subject to
  replicate and generation applicability restraint specified by parameters
  ``begin``, ``end``, ``step``, ``at``, and ``reps``.

* Operators specified in parameter ``finalOps`` of function
  :meth:`Simulator.evolve` will be applied to the final population after
  evolution, subject to replicate applicability restraint specified by parameter
  ``reps``.

Figure :ref:`fig_operator_orders <fig_operator_orders>` illustrated how
operators are applied to an evolutionary process. It worth noting that a default
during-mating operator is defined for each mating scheme. User-specfied
operators will **replace** the default operator so you need to explicitly
specify the default operator if you intent to add another one.

**Figure**: *Orders at which operators are applied during an evolutionary process*

.. _fig_operator_orders:

.. figure:: /Users/bpeng1/simuPOP/simuPOP/doc/figures/operators.png
   :width: 680


If you suspect that your simulation is not running as expected, you can have a
close look at your evolutionary process by setting the ``dryrun`` parameter of
an ``evolve`` function to ``True``, or by calling function
:func:`describeEvolProcess`\ (). This function takes the same set of parameters
as :meth:`Simulator.evolve`\ () and returns a description of the evolution
process, which might help you identify misuse of operators.

.. _describe:

**Example**: *describe an evolutionary process*

::

   >>> import simuPOP as sim
   >>> 
   >>> def outputstat(pop):
   ...     'Calculate and output statistics, ignored'
   ...     return True
   ... 
   >>> # describe this evolutionary process
   >>> print(sim.describeEvolProcess(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitInfo(lambda: random.randint(0, 75), infoFields='age'),
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...         sim.IdTagger(),
   ...         sim.PyOutput('Prevalence of disease in each age group:\n'),
   ...     ],
   ...     preOps=sim.InfoExec('age += 1'),
   ...     matingScheme=sim.HeteroMating([
   ...         sim.CloneMating(subPops=[(0,0), (0,1), (0,2)], weight=-1),
   ...         sim.RandomMating(ops=[
   ...             sim.IdTagger(),
   ...             sim.Recombinator(intensity=1e-4)
   ...         ], subPops=[(0,1)]),
   ...     ]),
   ...     postOps=[
   ...         sim.MaPenetrance(loci=0, penetrance=[0.01, 0.1, 0.3]),
   ...         sim.PyOperator(func=outputstat)
   ...     ],
   ...     gen = 100,
   ...     numRep = 3
   ... ))     
   Replicate 0 1 2:
   Apply pre-evolution operators to the initial population (initOps).
      * <simuPOP.InitSex> initialize sex randomly 
      * <simuPOP.InitInfo> initialize information field age using a Python
        function <lambda> 
      * <simuPOP.InitGenotype> initialize individual genotype acccording to 
        allele frequencies. 
      * <simuPOP.IdTagger> assign an unique ID to individuals 
      * <simuPOP.PyOutput> write 'Prevalence of disease in each age group:... '
        to output 

   Evolve a population for 100 generations
      * Apply pre-mating operators to the parental generation (preOps)
         # <simuPOP.InfoExec> execute statement age += 1 using information fields
           as variables. 

      * Populate an offspring populaton from the parental population using mating
        scheme <simuPOP.HeteroMating> a heterogeneous mating scheme with 2
        homogeneous mating schemes:
         # <simuPOP.HomoMating> a homogeneous mating scheme that uses
            - <simuPOP.SequentialParentChooser> chooses a parent sequentially
            - <simuPOP.OffspringGenerator> produces offspring using operators
               . <simuPOP.CloneGenoTransmitter> clone genotype, sex and
                 information fields of parent to offspring 
           in subpopulations (0, 0), (0, 1), (0, 2).
         # <simuPOP.HomoMating> a homogeneous mating scheme that uses
            - <simuPOP.RandomParentsChooser> chooses two parents randomly
            - <simuPOP.OffspringGenerator> produces offspring using operators
               . <simuPOP.IdTagger> assign an unique ID to individuals 
               . <simuPOP.Recombinator> genetic recombination. 
           in subpopulations (0, 1).


      * Apply post-mating operators to the offspring population (postOps).
         # <simuPOP.MaPenetrance> multiple-alleles penetrance 
         # <simuPOP.PyOperator> calling a Python function outputstat 

   No operator is applied to the final population (finalOps).


   now exiting runScriptInteractively...

`Download describe.py <describe.py>`_


