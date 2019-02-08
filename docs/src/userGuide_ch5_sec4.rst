Demographic changes
===================

A mating scheme controls the size of an offspring generation using parameter
``subPopSize``. This parameter has been described in detail in section
:ref:`subsec_offspring_size <subsec_offspring_size>`. In summary,

* The subpopulation sizes of the offspring generation will be the same as the
  parental generation if subPopSize is not set.

* The offspring generation will have a fixed size if ``subPopSize`` is set to a
  number (no subpopulation) or a list of subpopulation sizes.

* The subpopulation sizes of an offspring generation will be determined by the
  return value of a demographic function if ``subPopSize`` is set to such a
  function (a function that returns subpopulation sizes at each generation).

.. note::

   Parameter ``subPopSize`` only controls subpopulation sizes of an offspring
   generation immediately after it is generated. population or subpopulation sizes
   could be changed by other operators.

During mating, a mating scheme goes through each parental subpopulation and
populates its corresponding offspring subpopulation. This implies that

* Parental and offspring populations should have the same number of
  subpopulations.

* Mating happens strictly within each subpopulation.

This section will introduce several operators that allow you to move dndividuals
across the boundary of subpopulations (migration), and change the number of
subpopulations during evolution (split and merge). Please refer to
:ref:`subsec_offspring_size <subsec_offspring_size>` (control the size of the
offspring generation section of chapter mating scheme) for more details. For
more advanced demographic models, please refer to the :mod:`simuPOP.demography`
module.


Migration (operator :class:`Migrator`)
--------------------------------------


Migration by probability
^^^^^^^^^^^^^^^^^^^^^^^^

Operator :class:`Migrator` (and its function form ``migrate``) migrates
individuals from one subpopulation to another. The key parameters are

* *from* subpopulations (parameter ``subPops``). A list of subpopulations from
  which individuals migrate. Default to all subpopulations.

* *to* subpopulations (parameter ``toSubPops``). A list of subpopulations to
  which individuals migrate. Default to all subpopulations. **A new subpopulation
  ID can be specified to create a new subpopulation from migrants.**

* A migration rate matrix (parameter ``rate``). A :math:`m` by :math:`n` matrix
  ( a nested list in Python) that specifies migration rate from each source to
  each destination subpopulation. That is to say, :math:`\mbox{rate}{}_{i,j}`
  specifies migration rate from :math:`\mbox{subPops}_{i}` to
  :math:`\mbox{toSubPops}_{j}`. Needless to say, :math:`m` and :math:`n` are
  determined by the number of *from* and *to* subpopulations.

Example :ref:`migrateByProb <migrateByProb>` demonstrate the use of a
:class:`Migrator` to migrate individuals between three subpopulations. Note that

* Operator :class:`Migrator` relies on an information field ``migrate_to``
  (configurable) to record destination subpopulation of each individual so this
  information field needs to be added to a population befor migration.

* Migration rates to subpopulation themselves are determined automatically so
  they can be left unspecified.

.. _migrateByProb:

**Example**: *Migration by probability*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[1000]*3, infoFields='migrate_to')
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     preOps=sim.Migrator(rate=[
   ...             [0, 0.1, 0.1],
   ...             [0, 0, 0.1],
   ...             [0, 0.1, 0]
   ...         ]), 
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval('subPopSize'),
   ...         sim.PyOutput('\n')
   ...     ],
   ...     gen = 5
   ... )        
   [762, 1108, 1130]
   [601, 1175, 1224]
   [490, 1233, 1277]
   [395, 1282, 1323]
   [320, 1300, 1380]
   5

   now exiting runScriptInteractively...

`Download migrateByProb.py <migrateByProb.py>`_


Migration by proportion and counts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Migration rate specified in the rate parameter in Example :ref:`migrateByProb
<migrateByProb>` is intepreted as probabilities. That is to say, a migration
rate :math:`r_{m,n}` is interpreted as the probability at which any individual
in subpopulation :math:`m` migrates to subpopulation :math:`n`. The exact number
of migrants are randomly distributed.

If you would like to specify exactly how many migrants migrate from a
subpopulation to another, you can specify parameter ``mode`` of operator
:class:`Migrator` to ``BY_PROPORTION`` or ``BY_COUNTS``. The ``BY_PROPORTION``
mode interpret :math:`r_{m,n}` as proportion of individuals who will migrate
from subpopulation :math:`m` to :math:`n` so the number of :math:`m\rightarrow
n` migrant will be exactly :math:`r_{m,n}\times`\ subPopSize(m). In the
``BY_COUNTS`` mode, :math:`r_{m,n}` is interpretted as number of migrants,
regardless the size of subpopulation :math:`m`. Example
:ref:`migrateByPropAndCount <migrateByPropAndCount>` demonstrates these two
migration modes, as well as the use of parameters ``subPops`` and ``toSubPops.``

.. _migrateByPropAndCount:

**Example**: *Migration by proportion and count*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[1000]*3, infoFields='migrate_to')
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     preOps=sim.Migrator(rate=[[0.1], [0.2]],
   ...             mode=sim.BY_PROPORTION,
   ...             subPops=[1, 2],
   ...             toSubPops=[3]),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval('subPopSize'),
   ...         sim.PyOutput('\n')
   ...     ],
   ...     gen = 5
   ... )        
   [1000, 900, 800, 300]
   [1000, 810, 640, 550]
   [1000, 729, 512, 759]
   [1000, 657, 410, 933]
   [1000, 592, 328, 1080]
   5
   >>> #
   >>> pop.evolve(
   ...     preOps=sim.Migrator(rate=[[50, 50], [100, 50]],
   ...             mode=sim.BY_COUNTS,
   ...             subPops=[3, 2],
   ...             toSubPops=[2, 1]),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval('subPopSize'),
   ...         sim.PyOutput('\n')
   ...     ],
   ...     gen = 5
   ... )        
   [1000, 692, 328, 980]
   [1000, 792, 328, 880]
   [1000, 892, 328, 780]
   [1000, 992, 328, 680]
   [1000, 1092, 328, 580]
   5

   now exiting runScriptInteractively...

`Download migrateByPropAndCount.py <migrateByPropAndCount.py>`_


Theoretical migration models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To facilitate the use of widely used theoretical migration models, a few
functions are defined in module :mod:`simuPOP.demography`
:ref:`subsec_Predefined_migration_models <subsec_Predefined_migration_models>`.
These functions generate migration matrixes that can be plugged in to the
:class:`Migrator` operator.


migrate from virtual subpopulations \*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under a realistic eco-social settings, individuals in a subpopulation rarely
have the same probability to migrate. Genetic evidence has shown that female has
a higher migrate rate than male in humans, perhaps due to migration patterns
related to inter-population marriages. Such sex-biased migration also happens in
other large migration events such as slave trade.

It is easy to simulate most of such complex migration models by migrating from
virtual subpopulations. For example, if you define virtual subpopulations by
sex, you can specify different migration rates for males and females and control
the proportion of males among migrants, by specifying virtual subpopulations in
parameter ``subPops``. Parameter ``toSubPops`` does not accept virtual
subpopulations because you cannot, for example, migrate to females in a
subpopulation.

Example :ref:`migrateVSP <migrateVSP>` demonstrate a sex-biased migration model
where males dominate migrants from subpopulation 0. To avoid confusing, this
example uses the proportion migration mode. At the beginning of the first
generation, there are 500 males and 500 females in each subpopulation. A 10%
male migration rate and 5% female migration rate leads to 50 male migrants and
25 female migrants. Subpopulation sizes and number of males in each
subpopulation before mating are therefore:

* Subpopulation 0: male 500-50, female 500-25, total 925

* Subpopulation 1: male 500+50, female 500+25, total 1075

Note that the unspecified *to* subpopulations are subpopulation 0 and 1, which
cannot be virtual.

.. _migrateVSP:

**Example**: *Migration from virtual subpopulations*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[1000]*2, infoFields='migrate_to')
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> pop.evolve(
   ...     # 500 males and 500 females
   ...     initOps=sim.InitSex(sex=[sim.MALE, sim.FEMALE]),
   ...     preOps=[
   ...         sim.Migrator(rate=[
   ...             [0, 0.10],
   ...             [0, 0.05],
   ...             ],
   ...             mode = sim.BY_PROPORTION,
   ...             subPops=[(0, 0), (0, 1)]),
   ...         sim.Stat(popSize=True, numOfMales=True, vars='numOfMales_sp'),
   ...         sim.PyEval(r"'%d/%d\t%d/%d\n' % (subPop[0]['numOfMales'], subPopSize[0], "
   ...             "subPop[1]['numOfMales'], subPopSize[1])"),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(popSize=True, numOfMales=True, vars='numOfMales_sp'),
   ...         sim.PyEval(r"'%d/%d\t%d/%d\n' % (subPop[0]['numOfMales'], subPopSize[0], "
   ...             "subPop[1]['numOfMales'], subPopSize[1])"),
   ...     ],
   ...     gen = 2
   ... )   
   450/925	550/1075
   426/925	520/1075
   384/859	562/1141
   425/859	582/1141
   2

   now exiting runScriptInteractively...

`Download migrateVSP.py <migrateVSP.py>`_


Arbitrary migration models \*\*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If none of the described migration mothods fits your need, you can always resort
to manual migration. One such example is when you need to mimick an existing
evolutionary scenario so you know exactly which subpopulation each individual
will migrate to.

Manual migration is actually very easy. All you need to do is specifying the
destination subpopulation of all individuals in the *from* subpopulations
(parameter ``subPops``), using an information field (usually ``migrate_to``).
You can then call the :class:`Migrator` using ``mode=BY_IND_INFO``. Example
:ref:`manualMigration <manualMigration>` shows how to manually move individuals
around. This example uses the function form of :class:`Migrator`. You usually
need to use a Python operator to set destination subpopulations if you would
like to manually migrate individuals during an evolutionary process.

.. _manualMigration:

**Example**: *Manual migration*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([10]*2, infoFields='migrate_to')
   >>> pop.setIndInfo([0, 1, 2, 3]*5, 'migrate_to')
   >>> sim.migrate(pop, mode=sim.BY_IND_INFO)
   >>> pop.subPopSizes()
   (5, 5, 5, 5)

   now exiting runScriptInteractively...

`Download manualMigration.py <manualMigration.py>`_

.. note::

   individuals with an invalid destination subpopulation ID (e.g. an negative
   number) will be discarded silently. Although not recommended, this feature can
   be used to remove individuals from a subpopulation.


Migration using backward migration matrix (operator :class:`BackwardMigrator`)
------------------------------------------------------------------------------

Backward migration matrices are widely used in theoretical population genetics
and coalescent based simulations. Instead of specifying the probability of
migrating from one subpopulation to another (namely how migration happens), such
matrices specify the probability that individuals in a subpopulation originate
from others (namely the result of migration). simuPOP simulates such models by
converting backward migration matrices to foward ones using the theory described
below. Due to the limit of such models, simuPOP cannot simulate migration
from/to virtual subpopulatons, creation of new subpopulation, different source
and destination subpopulations, and will generate an error if the conversion
process fails.

To explain the differences between forward and backward migration matrices, let
us assume that there are :math:`d` subpopulations with population sizes
:math:`S=\left[S_{1},S_{2},...,S_{d}\right]`, and a forward migration matrix

.. math::

      F=\left[\begin{array}{cccc}
      f_{11} & f_{12} & \cdots & f_{1d}\\
      f_{21} & f_{22} & \cdots & f_{2d}\\
      \vdots &  &  & \vdots\\
      f_{d1} & f_{d2} & \cdots & f_{dd}
      \end{array}\right]

where :math:`f_{ij}` is the probability that an individual will migration from
subpopulation :math:`i` to :math:`j`. After migration happens, subppulation
sizes are changed to :math:`S'=\left[S'_{1},S'_{2},...,S'_{d}\right]`, and the
origin of individuals in each subpopulation can be described by the backward
migration matrix

.. math::

      B=\left[\begin{array}{cccc}
      b_{11} & b_{12} & \cdots & b_{1d}\\
      b_{21} & b_{22} & \cdots & b_{2d}\\
      \vdots &  &  & \vdots\\
      b_{d1} & b_{d2} & \cdots & b_{dd}
      \end{array}\right]

where :math:`b_{ij}` is the probability that an individual in subpopulation
:math:`i` originates from subpopulation :math:`j`.

These qualities can be derived from original population sizes and the forward
migration matrix. That is to say, the size of new subpopulation :math:`k` is the
sum of all migrants to this subpopulation

.. math::

      S'_{k}={\displaystyle \sum_{i=1}^{d}}S_{i}f_{ik}

and the size of the original population :math:`k` is the sum of all migrants
from this subpopulation

.. math::

      S_{k}={\displaystyle \sum_{i=1}^{d}}S'_{i}b_{ik}

and the composition of subpopulation :math:`k` (e.g. individuals originate from
subpopulation :math:`j`) is

.. math::

      b_{kj}=\dfrac{S_{j}f_{jk}}{S'_{k}}

In matrix form, these formulas can be written as

.. math::

      S'=F^{T}S

.. math::

      S=B^{T}S'

and

.. math::

      B=diag(S')^{-1}F^{T}diag(S)

Therefore, given a backward migration matrix :math:`B` and current population
size :math:`S`, we can derive a forward migration matrix using

.. math::

      S'=\left(B^{T}\right)^{-1}S

and

.. math::

      F=diag(S)^{-1}B^{T}diag(S')

Note that :math:`F=B` is always true if :math:`B` is symmetric and
:math:`S_{i}=S_{j}` (equal subpopulation size) so simuPOP will use :math:`B`
directly in this case. Also note that :math:`B` might not be inversable and
:math:`S'` and :math:`F` might be invalid (e.g. negative population size or
forward migration rate) for given :math:`B` and :math:`S`. simuPOP will
terminate with an error message in these cases.

The following example :ref:`backwardMigration <backwardMigration>` demonstrates
how to use a backward migration matrix to perform migration. It initializes all
individuals with indexes of subpopulations they belong to before migration and
calculates the percent of individuals from each source population using a
PyOperator with function originOfInds. The so-called overseved backward
migration matrix is similar to specified migration matrix despite of stochastic
effects. This example also uses turnOnDebug function to let the operator print
the expected subpopulation size (:math:`S'`) and calculate forward migration
matrix (:math:`F`) at each generation, which, as expected, vary from generation
to generation.

.. _backwardMigration:

**Example**: *Migration using a backward migration matrix*

::

   >>> import simuPOP as sim
   >>> sim.turnOnDebug('DBG_MIGRATOR')
   >>> pop = sim.Population(size=[10000, 5000, 8000], infoFields=['migrate_to', 'migrate_from'])
   >>> def originOfInds(pop):
   ...     print('Observed backward migration matrix at generation {}'.format(pop.dvars().gen))
   ...     for sp in range(pop.numSubPop()): 
   ...         # get source subpop for all individuals in subpopulation i
   ...         origins = pop.indInfo('migrate_from', sp)
   ...         spSize = pop.subPopSize(sp)
   ...         B_sp = [origins.count(j) * 1.0 /spSize for j in range(pop.numSubPop())]
   ...         print('    ' + ', '.join(['{:.3f}'.format(x) for x in B_sp]))
   ...     return True
   ... 
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     preOps=
   ...         # mark the source subpopulation of each individual
   ...         [sim.InitInfo(i, subPops=i, infoFields='migrate_from') for i in range(3)] + [
   ...         # perform migration
   ...         sim.BackwardMigrator(rate=[
   ...             [0, 0.04, 0.02],
   ...             [0.05, 0, 0.02],
   ...             [0.02, 0.01, 0]
   ...         ]),
   ...         # calculate and print observed backward migration matrix 
   ...         sim.PyOperator(func=originOfInds),
   ...         # calculate population size
   ...         sim.Stat(popSize=True),
   ...         # and print it
   ...         sim.PyEval(r'"Pop size after migration: {}\n".format(", ".join([str(x) for x in subPopSize]))'),
   ...         ], 
   ...     matingScheme=sim.RandomMating(),
   ...     gen = 5
   ... )        
   Expected next population size is 10211.4, 4851.8, 7936.84
   Forward migration matrix is 0.959867, 0.024259, 0.0158737, 0.0816908, 0.902435, 0.0158737, 0.0255284, 0.0121295, 0.962342
   Observed backward migration matrix at generation 0
       0.939, 0.040, 0.021
       0.051, 0.927, 0.022
       0.020, 0.010, 0.969
   Pop size after migration: 10218, 4859, 7923
   Expected next population size is 10453.6, 4690.64, 7855.79
   Forward migration matrix is 0.961671, 0.0229529, 0.0153764, 0.0860553, 0.897777, 0.0161675, 0.0263879, 0.0118406, 0.961772
   Observed backward migration matrix at generation 1
       0.942, 0.038, 0.020
       0.049, 0.932, 0.020
       0.023, 0.010, 0.968
   Pop size after migration: 10417, 4706, 7877
   Expected next population size is 10675.5, 4517.1, 7807.37
   Forward migration matrix is 0.963329, 0.0216814, 0.0149897, 0.0907397, 0.89267, 0.0165902, 0.0271056, 0.0114691, 0.961425
   Observed backward migration matrix at generation 2
       0.942, 0.039, 0.020
       0.048, 0.930, 0.022
       0.020, 0.010, 0.970
   Pop size after migration: 10660, 4536, 7804
   Expected next population size is 10946, 4323.5, 7730.53
   Forward migration matrix is 0.965217, 0.0202791, 0.0145038, 0.0965253, 0.886432, 0.0170426, 0.0280522, 0.0110802, 0.960868
   Observed backward migration matrix at generation 3
       0.940, 0.040, 0.020
       0.050, 0.930, 0.020
       0.020, 0.011, 0.969
   Pop size after migration: 10942, 4321, 7737
   Expected next population size is 11260.4, 4079.55, 7660
   Forward migration matrix is 0.967357, 0.0186417, 0.0140011, 0.104239, 0.878033, 0.0177274, 0.0291081, 0.0105456, 0.960346
   Observed backward migration matrix at generation 4
       0.937, 0.043, 0.021
       0.046, 0.933, 0.021
       0.019, 0.009, 0.972
   Pop size after migration: 11331, 4042, 7627
   5

   now exiting runScriptInteractively...

`Download backwardMigrate.py <backwardMigrate.py>`_


Split subpopulations (operators :class:`SplitSubPops`)
------------------------------------------------------

.. index:: single: SplitSubPops

Operator ``SplitSubPops`` splits one or more subpopulations into finer
subpopulations. It can be used to simulate populations that originate from the
same founder population. For example, a population of size 1000 in Example
:ref:`splitBySize <splitBySize>` is split into three subpopulations of sizes
300, 300 and 400 respectively, after evolving as a single population for two
generations.

.. _splitBySize:

**Example**: *Split subpopulations by size*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000)
   >>> pop.evolve(
   ...     preOps=[
   ...         sim.SplitSubPops(subPops=0, sizes=[300, 300, 400], at=2),
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
   ...     ],
   ...     matingScheme=sim.RandomSelection(),
   ...     gen = 4
   ... )
   Gen 0:	[1000]
   Gen 1:	[1000]
   Gen 2:	[300, 300, 400]
   Gen 3:	[300, 300, 400]
   4

   now exiting runScriptInteractively...

`Download splitBySize.py <splitBySize.py>`_

Operator :class:`SplitSubPops` splits a subpopulation by sizes of the resulting
subpopulations. It is often easier to do so with proportions. In addition, if a
demographic function is used, you should make sure that the number of
subpopulations will be the same before and after mating at any generation. One
way of doing this is to apply a :class:`SplitSubPops` operator at the right
generation. Example :ref:`splitByProp <splitByProp>` demonstrates such an
evolutionary scenario. However, it is often easier to split the population in
the demographic function in such case (see section
:ref:`subsec_Advanced_demo_func <subsec_Advanced_demo_func>` for details).

.. _splitByProp:

**Example**: *Split subpopulations by proportion*

::

   >>> import simuPOP as sim
   >>> def demo(gen, pop):
   ...     if gen < 2:
   ...         return 1000 + 100 * gen
   ...     else:
   ...         return [x + 50 * gen for x in pop.subPopSizes()]
   ... 
   >>> pop = sim.Population(1000)
   >>> pop.evolve(
   ...     preOps=[
   ...         sim.SplitSubPops(subPops=0, proportions=[.5]*2, at=2),
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
   ...     ],
   ...     matingScheme=sim.RandomSelection(subPopSize=demo),
   ...     gen = 4
   ... )
   Gen 0:	[1000]
   Gen 1:	[1000]
   Gen 2:	[550, 550]
   Gen 3:	[650, 650]
   4

   now exiting runScriptInteractively...

`Download splitByProp.py <splitByProp.py>`_

Either by *sizes* or by *proportions*, individuals in a subpopulation are
divided randomly. It is, however, also possible to split subpopulations
according to individual information fields. In this case, individuals with
different values at a given information field will be split into different
subpopulations. This is demonstrated in Example :ref:`splitByInfo <splitByInfo>`
where the function form of operator :class:`SplitSubPops` is used.

.. _splitByInfo:

**Example**: *Split subpopulations by individual information field*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population([1000]*3, subPopNames=['a', 'b', 'c'], infoFields='x')
   >>> pop.setIndInfo([random.randint(0, 3) for x in range(1000)], 'x')
   >>> print(pop.subPopSizes())
   (1000, 1000, 1000)
   >>> print(pop.subPopNames())
   ('a', 'b', 'c')
   >>> sim.splitSubPops(pop, subPops=[0, 2], infoFields=['x'])
   >>> print(pop.subPopSizes())
   (243, 244, 262, 251, 1000, 243, 244, 262, 251)
   >>> print(pop.subPopNames())
   ('a', 'a', 'a', 'a', 'b', 'c', 'c', 'c', 'c')

   now exiting runScriptInteractively...

`Download splitByInfo.py <splitByInfo.py>`_


Merge subpopulations (operator :class:`MergeSubPops`)
-----------------------------------------------------

Operator :class:`MergeSubPops` merges specified subpopulations into a single
subpopulation. This operator can be used to simulate admixed populations where
two or more subpopulations merged into one subpopulation and continue to evolve
for a few generations. Example :ref:`MergeSubPops <MergeSubPops>` simulates such
an evolutionary scenario. A demographic model could be added similar to Example
:ref:`splitByProp <splitByProp>`.

.. _MergeSubPops:

**Example**: *Merge multiple subpopulations into a single subpopulation*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([500]*2)
   >>> pop.evolve(
   ...     preOps=[
   ...         sim.MergeSubPops(subPops=[0, 1], at=3),
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
   ...     ],
   ...     matingScheme=sim.RandomSelection(),
   ...     gen = 5
   ... )
   Gen 0:	[500, 500]
   Gen 1:	[500, 500]
   Gen 2:	[500, 500]
   Gen 3:	[1000]
   Gen 4:	[1000]
   5

   now exiting runScriptInteractively...

`Download MergeSubPops.py <MergeSubPops.py>`_


Resize subpopulations (operator :class:`ResizeSubPops`)
-------------------------------------------------------

Whenever possible, it is recommended that subpopulation sizes are changed
naturally, namely through the population of an offspring generation. However, it
is sometimes desired to change the size of a population forcefully. Examples of
such applications include immediate expansion of a small population before
evolution, and the simulation of sudden population size change caused by natural
disaster. By default, new individuals created by such sudden population
expansion get their genotype from existing individuals. Example
:ref:`ResizeSubPops <ResizeSubPops>` shows a scenario where two subpopulations
expand instantly at generation 3.

.. _ResizeSubPops:

**Example**: *Resize subpopulation sizes*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([500]*2)
   >>> pop.evolve(
   ...     preOps=[
   ...         sim.ResizeSubPops(proportions=(1.5, 2), at=3),
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
   ...     ],
   ...     matingScheme=sim.RandomSelection(),
   ...     gen = 5
   ... )
   Gen 0:	[500, 500]
   Gen 1:	[500, 500]
   Gen 2:	[500, 500]
   Gen 3:	[750, 1000]
   Gen 4:	[750, 1000]
   5

   now exiting runScriptInteractively...

`Download ResizeSubPops.py <ResizeSubPops.py>`_


Time-dependent migration rate
-----------------------------

In evolutionary scenarios with complex demographic models, number of
subpopulations and migration rate might change from generation to generation.
For example, if one of the subpopulations is split into two, the migration
matrix has to be changed to accommendate increased number of subpopulations.

If there are a limited number of demographic changes and a few number of pre-
determined migration matrices. You can use a number of ``Migrators`` that are
applied at different generations. For example, you can use the following
operators to apply the first migration scheme during first ten generations (0,
..., 9), and the second migration scheme during the rest of the evolutionary
process::

   preOps=[
       Migrator(rate=M1, end=9),
       Migrator(rate=M2, begin=10),
   ]

If changes of demographies are frequent or stochastic so that migration matrices
can only be determined programmatically, it is easier to use a
:class:`PyOperator` to migrate populations using the function form of a
:class:`Migrator`. This is demonstrated in Example :ref:`varyingMigr
<varyingMigr>` where migration matrixes are computed dynamically due to random
split of subpopulations.

.. _varyingMigr:

**Example**: *Varying migration rate*

::

   >>> import simuPOP as sim
   >>> 
   >>> from simuPOP.utils import migrIslandRates
   >>> import random
   >>> 
   >>> def demo(pop):
   ...   # this function randomly split populations
   ...   numSP = pop.numSubPop()
   ...   if random.random() > 0.3:
   ...       pop.splitSubPop(random.randint(0, numSP-1), [0.5, 0.5])
   ...   return pop.subPopSizes()
   ... 
   >>> def migr(pop):
   ...   numSP = pop.numSubPop()
   ...   sim.migrate(pop, migrIslandRates(0.01, numSP))
   ...   return True
   ... 
   >>> pop = sim.Population(10000, infoFields='migrate_to')
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     preOps=[
   ...         sim.PyOperator(func=migr),
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')
   ...     ],
   ...     matingScheme=sim.RandomMating(subPopSize=demo),
   ...     gen = 5
   ... )
   Gen 0:	[10000]
   Gen 1:	[4982, 5018]
   Gen 2:	[2495, 2505, 5000]
   Gen 3:	[2509, 2517, 4974]
   Gen 4:	[2512, 2512, 4976]
   5

   now exiting runScriptInteractively...

`Download VaryingMigr.py <VaryingMigr.py>`_


