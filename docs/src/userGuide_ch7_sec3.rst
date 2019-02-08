Module :mod:`simuPOP.demography`
================================


.. _subsec_Predefined_migration_models:

Predefined migration models
---------------------------

The following functions are defined to generate migration matrixes for popular
migration models.

* ``migrIslandRates(r, n)`` returns a :math:`n\times n` migration matrix

  .. math::

      \left(\begin{array}{ccccc}
      1-r & \frac{r}{n-1} & ... & ... & \frac{r}{n-1}\\
      \frac{r}{n-1} & 1-r & ... & ... & \frac{r}{n-1}\\
       &  & ...\\
      \frac{r}{n-1} & ... & ... & 1-r & \frac{r}{n-1}\\
      \frac{r}{n-1} & ... & ... & \frac{r}{n-1} & 1-r
      \end{array}\right)

  for a traditional **island model** where individuals have equal probability of
  migrating to any other subpopulations. This model is also called a **migrant-
  pool island model**.

* ``migrHierarchicalIslandRates(r1, r2, n)`` models a **hierarchical island
  model** in which local populations are grouped into neighborhoods within which
  there is considerable gene flow and between which there is less gene flow.
  :math:`n` should be a list of group size. :math:`r_{1}` is the within-group
  migration rate and :math:`r_{2}` is the cross-group migration rate. That is to
  say, an individual in an island has probability :math:`1-r_{1}-r_{2}` to stay,
  :math:`r_{1}` to be a migratant to other islands in the group (migration rate
  depending on the size of group), and :math:`r_{2}` to be a migrant to other
  islands in another group (migration rate depending on the number of islands in
  other groups). Both :math:`r_{1}` and :math:`r_{2}` can vary across groups of
  islands. For example, ``migrHierarchicalIslandRates([r11, r12], r2, [3, 2])``
  returns a :math:`5\times5` migration matrix

  .. math::

      \left(\begin{array}{ccccc}
      1-r_{11}-r_{2} & \frac{r_{11}}{2} & \frac{r_{11}}{2} & \frac{r_{2}}{2} & \frac{r_{2}}{2}\\
      \frac{r_{11}}{2} & 1-r_{11}-r_{2} & \frac{r_{11}}{2} & \frac{r_{2}}{2} & \frac{r_{2}}{2}\\
      \frac{r_{11}}{2} & \frac{r_{11}}{2} & 1-r_{11}-r_{2} & \frac{r_{2}}{2} & \frac{r_{2}}{2}\\
      \frac{r_{2}}{3} & \frac{r_{2}}{3} & \frac{r_{2}}{3} & 1-r_{12}-r_{2} & r_{12}\\
      \frac{r_{2}}{3} & \frac{r_{2}}{3} & \frac{r_{2}}{3} & r_{12} & 1-r_{12}-r_{2}
      \end{array}\right)

* ``migrSteppingStoneRates(r, n, circular=False)`` returns a :math:`n\times n`
  migration matrix

  .. math::

      \left(\begin{array}{ccccc}
      1-r & r\\
      r/2 & 1-r & r/2\\
       &  & ...\\
       &  & r/2 & 1-r & r/2\\
       &  &  & r & 1-r
      \end{array}\right)

  and if ``circular=True``, returns

  .. math::

      \left(\begin{array}{ccccc}
      1-r & r/2 &  &  & r/2\\
      r/2 & 1-r & r/2\\
       &  & ...\\
       &  & r/2 & 1-r & r/2\\
      r/2 &  &  & r/2 & 1-r
      \end{array}\right)

* ``migr2DSteppingStoneRates(r, m, n, diagonal=False, circular=False)``\ models
  a 2D stepping stone model in which local populations are arranged into a lattice
  of :math:`m\times n` (:math:`m` rows, :math:`n` columns) patches. The population
  thus needs to have :math:`m\times n` subpopulations with subpopulation indexes
  counted by row. In this model, an individual in a center patch has a probability
  of :math:`1-r` to stay, and :math:`r/4` to migrate to its neighbor patches if
  ``diagonal`` is set to False, or :math:`r/8` to migrate to 8 neighbors
  (including diagnal ones) if ``range`` is set to 8. If ``circular`` is set to
  ``False``, the corner patch has a probability of :math:`r/2` or :math:`r/3` (if
  range=8) to migrate, and a side patch has a probability :math:`r/3` or
  :math:`r/5` to migrate. If ``circular`` is set to ``True``, the lattice will be
  conceptually connected to a ball so that there is no boundary effect. For
  example, for a 3 by 2 lattice

  .. math::

      \left(\begin{array}{cc}
      0 & 1\\
      2 & 3\\
      4 & 5
      \end{array}\right)

  with ``diagonal=False`` and ``circular=False``, the migration matrix will be

  .. math::

      \left(\begin{array}{cccccc}
      1-r & \frac{r}{2} & \frac{r}{2}\\
      \frac{r}{2} & 1-r &  & \frac{r}{2}\\
      \frac{r}{3} &  & 1-r & \frac{r}{3} & \frac{r}{3}\\
       & \frac{r}{3} & \frac{r}{3} & 1-r &  & \frac{r}{3}\\
       &  & \frac{r_{2}}{2} &  & 1-r & \frac{r}{2}\\
       &  &  & \frac{r}{2} & \frac{r}{2} & 1-r
      \end{array}\right)

Many more migration models have been proposed and studied, sometimes under
different names with slightly different definitions. If you cannot find your
model there, it should not be too difficult to construct a migration rate matrix
for it. I will be glad to add such functions to this module if you could provide
a reference and your implementation of the model.


Uniform interface of demographic models
---------------------------------------

A realistic demographic models can be very complex that involves population
growth, population bottleneck, subdivided populations, migration, population
split and admixture for a typical demographic model for human populations, and
carrying capacity, fecunity, sex distribution and many more factors for more
complex ones (e.g. models for animal populations under continuous habitat). The
goal of this module is to provide a common interface for demographic models,
classes for frequently used demographic models, and several pre-defined
demographic models for human populations. More complex demographic models will
be added if needed.

A demographic model usually consists of the following components:

* An initial population size that is used to initialize a population (the
  ``size`` parameter of ``sim.Population``)

* One or more operators to split and merge populations (e.g. Operators
  :class:`SplitSubPops`)

* One or more operators to migrate individuals across subpopulations (e.g.
  operator :class:`Migrator`)

* Determine sizes of subpopulations before mating (parameter ``subPopSize`` of a
  mating scheme)

* Number of generations to evolve (parameter ``gen`` of the ``evolve`` function)
  or operators to terminate the evolution conditionally (e.g. operator
  :class:`TerminateIf`)

Using an object-oriented approach, a demographic model defined in this module
encapsulates all these in a single object. More specifically, a demographic
object ``model`` is a callable Python object that

* has attribute ``model.init_size`` and ``model.info_fields`` to determine the
  initial population size and required information fields to construct an initial
  population (e.g., ``sim.Population(size=model.init_size,
  infoFields=model.info_fields + ['my_fields'])``)

* handles population split, merge, migration etc internally before mating when
  it is passed to parameter ``subPopSize`` of a mating scheme. (e.g.
  :class:`RandomMating`\ (``subPopSize=model``))

* has attribute ``model.num_gens`` to determine the number of generations to
  evolve (e.g. ``pop.evolve(..., gen=model.num_gens)``). The model can optionally
  terminate the evolution by returnning an empty offspring population size before
  mating.

* provides a function ``model.plot(filename='', title='')`` to plot the
  demographic function. It by default prints out population sizes whenever
  population size changes. If a ``filename`` is specified and if module
  ``matplotlib`` is available, it will plot the demographic model and save it to
  filename. A ``title`` can be specified for the figure. This function actually
  use the demographic model to evolve a haploid population using
  :class:`RandomSelection` mating scheme, which is a good way to test if your
  demographic model works properly.

* saves population sizes of evolved generations, which makes it possible to
  revert an evolutionary process to an previous state using operator
  :class:`RevertIf`.

A demographic model can be defined in two ways. The first approach is to specify
the size of subpopulations at each generation, and the second approach is to
specify the events that change population sizes. The :mod:`simuPOP.demography`
module provides functions and classes to define demographic models using both
approaches and you can use the one that is most convenient for your model.


Demographic models defined by outcomes
--------------------------------------

The :mod:`simuPOP.demography` module defines a number of widely used demographic
models, including linear and exponential population growth with carrying
capacity, shrink, split and merge, and bottleneck.

For example,

*  ::

     InstantChangeModel(T=1000, N0=1000, G=500, NG=2000)

  defines an instant population growth model that expands a population of size
  from 1000 to 2000 instantly at generation 500

*  ::

     InstantChangeModel(T=1000, N0=1000, G=[500, 600], NG=[100, 1000])

  defines a bottleneck model that introduces a bottleneck of size 100 between
  generation 500 and 600 to a population of size 1000

*  ::

     InstantChangeModel(T=1000, N0=1000, G=500, NG=[[400, 600]])

  defines a bottleneck model that split a population of size into two
  subpopulations of sizes 400 and 600 at generation 500

*  ::

     ExponentialGrowthModel(T=100, N0=1000, NT=10000)

  expands a population of size 1000 to 10000 in 100 generations

*  ::

     ExponentialGrowthModel(T=100, N0=[200, 800], r=[0.02, 0.01],
         ops=Migrator(rate=[[0, 0.1], [0.1, 0]])

  expands a population of two subpopulation sizes at rate ``0.02`` and ``0.01``
  for ``100`` generations, with migration between these two subpopulations. The
  initial population will be resized (split if necessary) to two populations of
  sizes 200 and 800.

*  ::

     LinearGrowthModel(N0=(200, 'A'), r=0.02, NT=1000)

  expands a population of size ``200`` at a rate 0f ``0.02`` (add 4 individuals at
  each generation) until it reaches size ``1000``. Here the initial size is
  expressed as a size name tuple, which directs the demographic model to assign
  the name ``A`` to the initial population. Such named size is acceptable for all
  places where population size is needed.

Here we specify only two of the three parameters for linear and exponential
growth models and allow simuPOP to figure out the rest. If all three parameters
are specified, the ending population size will be interpretted as carraying
capacity, namely population growth (or decline of negative rates are specified)
will stop after it reaches the specified size.

A demographic model does not have to have a fixed initial population size. If an
initial population size is not provided, its size will be determined from the
population when it is first applied to. For example

*  ::

     InstantChangeModel(T=100, G=50, NT=[0.5, 0.5])

  split a population into two equally sized subpopulations at generation 50. The
  ending population size is set to ``[0.5, 0.5]``, which means 50% of the size at
  time ``G``.

*  ::

     InstantChangeModel(T=100, G=50, NT=[None, 100])

  forks a population of size 100 from the main population at generation 50.
  ``NT=[None, 100]`` is equivalent to ``NT=[1.0, 100]`` in this case.

*  ::

     InstantChangeModel(T=0, removEmptySubPops=True)

  removes all empty subpopulations from the existing subpopulation. Here we do not
  specify an input population size because the the size of the input population
  will be kept.

*  ::

     InstantChangeMoel(T=0, N0=[None, 0, None], removEmptySubPops=True)

  removes the second of the three subpopulations while keep other two
  subpopulations intact. The input population of this demographic model must have
  three subpopulations.

*  ::

     ExponentialGrowthModel(T=100, NT=[10000, 20000])

  expands a population of two subpopulations to sizes ``10000`` and ``20000`` in
  ``100`` generations. An error will be raised if the population does not have two
  subpopulations.

*  ::

     ExponentialGrowthModel(T=100, N0=[1., 400], NT=[10000, 20000], 
         ops=Migrator(rate=[[0, 0.1], [0.1, 0]])

  split a population into two subpopulations. The first one keeps all individuals
  (100%), the second one with 400 individuals, and then expands them, with
  migration, to sizes ``10000`` and ``20000`` in ``100`` generations.

The demography model also defines two models for population admxture. The HI
model (Hybrid Isolation) model creates a separate subpopulation with :math:`\mu`
and :math:`1-\mu` individuals from two specified subpopulations. The CGF
(Continuous Gene Flow) model replaces :math:`1-\mu` individuals from the doner
population at each generation, thus keep both the recipient and doner population
constant in size. For example,

*  ::

     AdmixtureModel(model=('HI', 1, 3, 0.5, 'Admixed'), T=10)

  Creates a separate population with 50% of individuals from subpopulation 1 and
  50% of individuals from subpopulation 3, regardless if population sizes 1 and 3
  have the same number of individuals. An optional name Admixed is assigned to the
  new subpopulation. The admixed population will evolve independently for 10
  generations.

*  ::

     AdmixtureModel(model=('CGF', 1, 3, 0.9), T=10)

  Replaces 10% of individuals in subpopulation 1 with individuals from
  subpopulation 3 for 10 generations.

As you can imagine, these models do not provide a valid ``init_size`` to
initialize a population. As a matter of fact, they are mostly stacked to other
demographic models to form more complex demographic models, in model
``MultiStageModel``. For example,

*  ::

     MultiStageModel([
         InstantChangeModel(T=1000, N0=1000, G=[500, 600], NG=[100, 1000]),  
         ExponentialGrowthModel(T=100, NT=10000)
     ])

  defines a demographic model with a bottleneck followed by exponential population
  growth. ``N0`` of the second stage is not specified because it is determined
  from its previous stage.

*  ::

     MultiStageModel([
         LinearGrowthModel(T=100, N0=1000, r=0.01),  
         ExponentialGrowthModel(T=100, N0=[0.4, 0.6], r=0.001),
         ExponentialGrowthModel(r=0.01, NT=[2000, 4000]),
         AdmixtureModel(model=('HI', 0, 1, 0.8, 'admixed'), T=10)
     ])

  defines a demographic model that expands a single population linearly for 100
  generations, split into two subpopulations and grow exponentially at a rate of
  0.001, and growth at a higher rate of 0.01 until they reaches sizes 2000 and
  4000 respectively. This stage is tricky because one of the subpopulations will
  reach its carrying capacity sooner and keep a contant population size
  afterwards. As the last step, the two populations admixed and formed a new
  subpopulation called ``admixed``. The model is depicted in figure
  :ref:`fig_multi_stage <fig_multi_stage>`

  **Figure**: *A linear and two stage exponential population growth model, followed by population admixture*

  .. _fig_multi_stage:

  .. figure:: log/MultiStage.png
     :width: 680
  

Example :ref:`demoModel <demoModel>` defines a demographic model use it to
evolve a population. The demographic model is depicted in Figure
:ref:`fig_demoModel_example <fig_demoModel_example>`.

.. _demoModel:

**Example**: *A demographic model for human population*

::

   >>> import simuPOP as sim
   >>> from simuPOP.demography import *
   >>> model = MultiStageModel([
   ...     InstantChangeModel(T=200, 
   ...         # start with an ancestral population of size 1000
   ...         N0=(1000, 'Ancestral'),
   ...         # change population size at 50 and 60
   ...         G=[50, 60], 
   ...         # change to population size 200 and back to 1000
   ...         NG=[(200, 'bottleneck'), (1000, 'Post-Bottleneck')]),
   ...     ExponentialGrowthModel(
   ...         T=50, 
   ...         # split the population into two subpopulations
   ...         N0=[(400, 'P1'), (600, 'P2')],
   ...         # expand to size 4000 and 5000 respectively
   ...         NT=[4000, 5000])]
   ...     )
   >>> #
   >>> # model.init_size returns the initial population size
   >>> # migrate_to is required for migration
   >>> pop = sim.Population(size=model.init_size, loci=1,
   ...     infoFields=model.info_fields)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     matingScheme=sim.RandomMating(subPopSize=model),
   ...     finalOps=
   ...         sim.Stat(alleleFreq=0, vars=['alleleFreq_sp']),
   ...     gen=model.num_gens
   ... )
   250
   >>> # print out population size and frequency
   >>> for idx, name in enumerate(pop.subPopNames()):
   ...     print('%s (%d): %.4f' % (name, pop.subPopSize(name), 
   ...         pop.dvars(idx).alleleFreq[0][0]))
   ... 
   P1 (4000): 0.6185
   P2 (5000): 0.7218
   >>> # get a visual presentation of the demographic model
   >>> model.plot('log/demoModel.png',
   ...     title='A bottleneck + exponential growth demographic model')
   A bottleneck + exponential growth demographic model
   0: 1000 (Ancestral)
   50: 200 (bottleneck)
   60: 1000 (Post-Bottleneck)
   200: 419 (P1), 626 (P2)
   201: 439 (P1), 653 (P2)
   202: 459 (P1), 681 (P2)
   203: 481 (P1), 711 (P2)
   204: 504 (P1), 742 (P2)
   205: 527 (P1), 774 (P2)
   206: 552 (P1), 807 (P2)
   207: 578 (P1), 842 (P2)
   208: 605 (P1), 879 (P2)
   209: 634 (P1), 917 (P2)
   210: 664 (P1), 957 (P2)
   211: 695 (P1), 998 (P2)
   212: 728 (P1), 1041 (P2)
   213: 762 (P1), 1086 (P2)
   214: 798 (P1), 1133 (P2)
   215: 836 (P1), 1183 (P2)
   216: 875 (P1), 1234 (P2)
   217: 916 (P1), 1287 (P2)
   218: 960 (P1), 1343 (P2)
   219: 1005 (P1), 1401 (P2)
   220: 1052 (P1), 1462 (P2)
   221: 1102 (P1), 1525 (P2)
   222: 1154 (P1), 1591 (P2)
   223: 1208 (P1), 1660 (P2)
   224: 1265 (P1), 1732 (P2)
   225: 1325 (P1), 1807 (P2)
   226: 1387 (P1), 1885 (P2)
   227: 1452 (P1), 1967 (P2)
   228: 1521 (P1), 2052 (P2)
   229: 1592 (P1), 2141 (P2)
   230: 1667 (P1), 2234 (P2)
   231: 1746 (P1), 2331 (P2)
   232: 1828 (P1), 2432 (P2)
   233: 1915 (P1), 2537 (P2)
   234: 2005 (P1), 2647 (P2)
   235: 2099 (P1), 2761 (P2)
   236: 2198 (P1), 2881 (P2)
   237: 2302 (P1), 3006 (P2)
   238: 2410 (P1), 3136 (P2)
   239: 2524 (P1), 3272 (P2)
   240: 2643 (P1), 3414 (P2)
   241: 2767 (P1), 3562 (P2)
   242: 2898 (P1), 3716 (P2)
   243: 3034 (P1), 3877 (P2)
   244: 3177 (P1), 4045 (P2)
   245: 3327 (P1), 4220 (P2)
   246: 3484 (P1), 4403 (P2)
   247: 3648 (P1), 4593 (P2)
   248: 3820 (P1), 4792 (P2)
   249: 4000 (P1), 5000 (P2)
   Traceback (most recent call last):
     File "/var/folders/ys/gnzk0qbx5wbdgm531v82xxljv5yqy8/T/tmpdvg5jvxd", line 2, in <module>
       #begin_ignore
     File "/Users/bpeng1/anaconda3/envs/sos/lib/python3.6/site-packages/simuPOP/demography.py", line 446, in plot
       region = region.reshape(region.size / 4, 4)
   TypeError: 'float' object cannot be interpreted as an integer

   now exiting runScriptInteractively...

`Download demoModel.py <demoModel.py>`_

**Figure**: *A exponential population growth followed by bottleneck demographic model*

.. _fig_demoModel_example:

.. figure:: log/demoModel.png
   :width: 680



Demographic models defined by population changes (events)
---------------------------------------------------------

Another way to define a demographic model is to specify the events that changes
population sizes. This approach can be easier to use because it conforms with
the way many demographic models are specified, also because the events can be
specified for a subset of subpopulations so you can, for example, split one
subpopulation without worrying about its impact on other subpopulations.

A event-based demographic model is defined using  ::

   EventBasedModel(events=[], T=None, N0=None, ops=[], infoFields=[])

where\ ``T`` and ``N0`` are the duration and initial size of the demographic
model, respectively, and ops is the operators that will be applied to the
population (without checking applicability). Parameter ``events`` acepts one or
more of ``DemographicEvent`` and its derived classes. For example,

::

   ExpansionEvent(rates=0.05, begin=500)

expands all subpopulations exponentially at a rate of 0.05, and

::

   ExpansionEvent(rates=[0.05, 0.01], capacity=10000, subPops=[0, 2], begin=500)

expands two subpopulations at rates 0.05 and 0.01 respectively, until they reach
10000 individuals in each subpopulation.

::

   ExpansionEvent(slopes=500, subPops=[0, 2], begin=500)

expands the populations linearly by adding 500 individuals to each subpopulation
at each generation. These events happen at each generation starting from
generation 500.

Simiarly, you can split, merge, and resize subpopulations using events
``SplitEvent``, ``MergeEvent``, and ``ResizeEvent``. For example,

::

   SplitEvent(subPops='AF', sizes=[500, 500], names=['AF', 'EU'], at=-4000)

splits an ancestral population named AF to two populations AF and EU at 4000
generations before the end of the demographic model. The AF population will be
expanded automatically if it does not have 1000 individuals.

Finally, an ``AdmixtureEvent`` mix two or more subpopulations by certain
proportions, and either create a new subpopulation or replace an existing
subpopulation. In particular,

::

   AdmixtureEvent(subPops=['MX', 'EU'], at=-10, sizes=[0.4, 0.6], name='MXL')

creates a new admixed population called MXL with 40% of individuals from the MX
population, and the rest from the EU population. The admixture process happens
once and follows an Hybrid Isolation model. Alternatively,

::

   AdmixtureEvent(subPops=['MX', 'EU'], begin=-10, sizes=[0.8, 0.2], toSubPop='MX')

will create an admixed population with 80% MX and 20% EU individuals for 10
generations. Because 20% of the admixed population will be replaced by
individuals from the EU population, this models a continuous gene flow model of
admixture. If you would like to control the exact size of the admixed
population, you can specify the number of individuals as integer numbers instead
of proportions:

::

   AdmixtureEvent(subPops=['MX', 'EU'], begin=-10, sizes=[int(1400*0.8), int(1400*0.2)], toSubPop='MX')

Note that the type of elements in parameter ``sizes`` is important, ``1.``
stands for all subpopulation and ``1`` stands for one individual from it.

Example\ :ref:`demoEventModel <demoEventModel>` defines the same model as
:ref:`demoModel <demoModel>` using an event based demographic model. The result
is depicted in Figure :ref:`fig_demoEventModel_example
<fig_demoEventModel_example>`. These two models look similar but the event-based
model does not have the same final population sizes as the previous model. This
is because the population size of the previous model was calculated by
:math:`N(t)=N(0)\exp(rt)` whereas the event based model was calculated using
:math:`N(t)=\mbox{round}(N(t-1)\*(1+r))` for each generation, and the integer
rounding error accumulates over time.

.. _demoEventModel:

**Example**: *A event-based demographic model*

::

   >>> import simuPOP as sim
   >>> from simuPOP.demography import *
   >>> import math
   >>> model = EventBasedModel(
   ...     N0=(1000, 'Ancestral'),
   ...     T=250,
   ...     events=[
   ...         ResizeEvent(at=50, sizes=200),
   ...         ResizeEvent(at=60, sizes=1000),
   ...         SplitEvent(sizes=[0.4, 0.6], names=['P1', 'P2'], at=200),
   ...         ExpansionEvent(rates=[math.log(4000/400)/50, math.log(5000/600)/50], begin=200)
   ...     ]
   ... )
   >>> #
   >>> # model.init_size returns the initial population size
   >>> # migrate_to is required for migration
   >>> pop = sim.Population(size=model.init_size, loci=1,
   ...     infoFields=model.info_fields)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     matingScheme=sim.RandomMating(subPopSize=model),
   ...     finalOps=
   ...         sim.Stat(alleleFreq=0, vars=['alleleFreq_sp']),
   ...     gen=model.num_gens
   ... )
   250
   >>> # print out population size and frequency
   >>> for idx, name in enumerate(pop.subPopNames()):
   ...     print('%s (%d): %.4f' % (name, pop.subPopSize(name), 
   ...         pop.dvars(idx).alleleFreq[0][0]))
   ... 
   P1 (4000): 0.6185
   P2 (5000): 0.7218
   >>> # get a visual presentation of the demographic model
   >>> model.plot('log/demoEventModel.png',
   ...     title='A event-based bottleneck + exponential growth demographic model')
   A event-based bottleneck + exponential growth demographic model
   0: 1000 (Ancestral)
   50: 200 (Ancestral)
   60: 1000 (Ancestral)
   200: 419 (P1), 626 (P2)
   201: 439 (P1), 653 (P2)
   202: 459 (P1), 681 (P2)
   203: 481 (P1), 711 (P2)
   204: 504 (P1), 742 (P2)
   205: 527 (P1), 774 (P2)
   206: 552 (P1), 807 (P2)
   207: 578 (P1), 842 (P2)
   208: 605 (P1), 879 (P2)
   209: 634 (P1), 917 (P2)
   210: 664 (P1), 957 (P2)
   211: 695 (P1), 998 (P2)
   212: 728 (P1), 1041 (P2)
   213: 762 (P1), 1086 (P2)
   214: 798 (P1), 1133 (P2)
   215: 836 (P1), 1183 (P2)
   216: 875 (P1), 1234 (P2)
   217: 916 (P1), 1287 (P2)
   218: 960 (P1), 1343 (P2)
   219: 1005 (P1), 1401 (P2)
   220: 1052 (P1), 1462 (P2)
   221: 1102 (P1), 1525 (P2)
   222: 1154 (P1), 1591 (P2)
   223: 1208 (P1), 1660 (P2)
   224: 1265 (P1), 1732 (P2)
   225: 1325 (P1), 1807 (P2)
   226: 1387 (P1), 1885 (P2)
   227: 1452 (P1), 1967 (P2)
   228: 1521 (P1), 2052 (P2)
   229: 1592 (P1), 2141 (P2)
   230: 1667 (P1), 2234 (P2)
   231: 1746 (P1), 2331 (P2)
   232: 1828 (P1), 2432 (P2)
   233: 1915 (P1), 2537 (P2)
   234: 2005 (P1), 2647 (P2)
   235: 2099 (P1), 2761 (P2)
   236: 2198 (P1), 2881 (P2)
   237: 2302 (P1), 3006 (P2)
   238: 2410 (P1), 3136 (P2)
   239: 2524 (P1), 3272 (P2)
   240: 2643 (P1), 3414 (P2)
   241: 2767 (P1), 3562 (P2)
   242: 2898 (P1), 3716 (P2)
   243: 3034 (P1), 3877 (P2)
   244: 3177 (P1), 4045 (P2)
   245: 3327 (P1), 4220 (P2)
   246: 3484 (P1), 4403 (P2)
   247: 3648 (P1), 4593 (P2)
   248: 3820 (P1), 4792 (P2)
   249: 4000 (P1), 5000 (P2)
   >>> 

   now exiting runScriptInteractively...

`Download demoEventModel.py <demoEventModel.py>`_

**Figure**: *A event-based demographic model*

.. _fig_demoEventModel_example:

.. figure:: /Users/bpeng1/simuPOP/simuPOP/doc/log/demoEventModel.png
   :width: 680



Predefined demographic models for human populations
---------------------------------------------------

The :mod:`simuPOP.demography` module currently defines the following models

* Out of Africa model for YRI, CEU and CHB populations (:ref:`fig_Out_of_Africa
  <fig_Out_of_Africa>`),

  ::

     OutOfAfricaModel(10000).plot('OutOfAfrica.png')

  **Figure**: *Out of Africa model for YRI, CEU, and CHB populations*

  .. _fig_Out_of_Africa:

  .. figure:: log/OutOfAfrica.png
     :width: 680
  

* The settlement of new world model for Mexican American
  (:ref:`fig_Settlement_of_New <fig_Settlement_of_New>`) ( Gutenkunst, 2009, PLoS
  Genetics). In this model, the simulated CHB and MX populations are mixed to
  produce an admixed population at the last generation.

  ::

     SettlementOfNewWorldModel(10000).plot('SettlementOfNewWorld.png')

  **Figure**: *Settlement of New World model for Mexican America population*

  .. _fig_Settlement_of_New:

  .. figure:: log/SettlementOfNewWorld.png
     :width: 680
  

* The demographic model developed by cosi (Schaffner, 2005, genome research).

  ::

     CosiModel(20000).plot('Cosi.png')

  **Figure**: *Demographic models for African, Asian and European populations (cosi)*

  .. _fig_cosi:

  .. figure:: log/Cosi.png
     :width: 680
  

These functions all accept a parameter scale. If specified, it will scale all
population sizes and generation numbers by the specified scaling factor. For
example

::

   CosiModel(20000, scale=10)

will result in a demographic model that evolves 2000 instead of 20000
generations, with all population sizes reduced by a factor of 10. Note that the
burn-in period of the examples above are relatively short and you might need to
use a longer burn-in period (e.g. T=100,000 generations for a burn-in period of
about 80,000 generations).


Demographic model without predefined generations to evolve \*
-------------------------------------------------------------

All migration models accept one or more operators that will be applied to the
population before population population changes are applied. The most frequently
application of this operator is to pass a migrator to the model, but we can also
pass an operator to terminate a demographic model under certain conditions. For
example, Example :ref:`demoTerminate <demoTerminate>` defines a demographic
model that starts with a burn-in stage with indefinite size and will stop if the
average allele frequency at segregating sites exceeds 0.1. It splits to two
equally sized subpopulations and expand rate a rate of 0.01 to size 2000 and
5000 respectively.

.. _demoTerminate:

**Example**: *A demographic model with a terminator*

::

   >>> import simuPOP as sim
   simuPOP Version 1.1.9 : Copyright (c) 2004-2016 Bo Peng
   Revision 4583 (Oct 10 2018) for Python 3.6.6 (64bit, 0thread)
   Random Number Generator is set to mt19937 with random seed 0x81aae4a664e115de.
   This is the standard short allele version with 256 maximum allelic states.
   For more information, please visit http://simupop.sourceforge.net,
   or email simupop-list@lists.sourceforge.net (subscription required).
   >>> import simuPOP.demography as demo
   >>> 
   >>> model = demo.MultiStageModel([
   ...     demo.InstantChangeModel(N0=1000, 
   ...         ops=[
   ...             sim.Stat(alleleFreq=sim.ALL_AVAIL, numOfSegSites=sim.ALL_AVAIL),
   ...             # terminate if the average allele frequency of segregating sites
   ...             # are more than 0.1 
   ...             sim.TerminateIf('sum([x[1] for x in alleleFreq.values() if '
   ...                 'x[1] != 0])/(1 if numOfSegSites==0 else numOfSegSites) > 0.1')
   ...         ]
   ...     ),
   ...     demo.ExponentialGrowthModel(N0=[0.5, 0.5], r=0.01, NT=[2000, 5000])
   ...     ]
   ... )
   >>> 
   >>> pop = sim.Population(size=model.init_size, loci=100)
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     preOps=sim.SNPMutator(u=0.001, v=0.001),
   ...     matingScheme=sim.RandomMating(subPopSize=model),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=sim.ALL_AVAIL, numOfSegSites=sim.ALL_AVAIL,
   ...             popSize=True, step=50),
   ...         sim.PyEval(r'"%d: %s, %.3f\n" % (gen, subPopSize, sum([x[1] for x '
   ...             'in alleleFreq.values() if x[1] != 0])/(1 if numOfSegSites == 0 '
   ...             'else numOfSegSites))', step=50)
   ...     ],
   ... )
   0: [1000], 0.001
   50: [1000], 0.047
   100: [1000], 0.089
   150: [738, 738], 0.128
   200: [1218, 1218], 0.166
   250: [2000, 2007], 0.199
   300: [2000, 3310], 0.230
   343
   >>> 

   now exiting runScriptInteractively...

`Download demoTerminate.py <demoTerminate.py>`_


