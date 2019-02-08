.. _sec_Mating_Schemes:

Mating Schemes
==============

.. index:: single: mating scheme

Mating schemes are responsible for populating an offspring generation from the
parental generation. There are currently two types of mating schemes

* A **homogeneous mating scheme** is the most flexible and most frequently used
  mating scheme and is the center topic of this section. A homogeneous mating is
  composed of a *parent chooser* that is responsible for choosing parent(s) from a
  (virtual) subpopulation and an *offspring generator* that is used to populate
  all or part of the offspring generation. During-mating operators are used to
  transmit genotypes from parents to offspring. Figure
  :ref:`fig_homogeneous_mating_scheme <fig_homogeneous_mating_scheme>`
  demonstrates this process.

* A **heterogeneous mating scheme** applies several homogeneous mating scheme to
  different (virtual) subpopulations. Because the division of virtual
  subpopulations can be arbitrary, this mating scheme can be used to simulate
  mating in heterogeneous populations such as populations with age structure.

* A **pedigree mating scheme** evolves a population by following the pedigree
  structure of a pedigree. This mating scheme is used to a replay a recorded or
  manually created evolutionary process.

This section describes some standard features of mating schemes and most pre-
defined mating schemes. The next section will demonstrate how to build complex
nonrandom mating schemes from scratch.

**Figure**: *A homogeneous mating scheme*

.. _fig_homogeneous_mating_scheme:

.. figure:: /Users/bpeng1/simuPOP/simuPOP/doc/figures/HomoMatingScheme.png
   :width: 680


A homogeneous mating scheme is responsible to choose parent(s) from a
subpopulation or a virtual subpopulation, and population part or all of the
corresponding offspring subpopulation. A parent chooser is used to choose one or
two parents from the parental generation, and pass it to an offspring generator,
which produces one or more offspring. During mating operators such as taggers
and Recombinator can be applied when offspring is generated.


.. _subsec_offspring_size:

Control the size of the offspring generation
--------------------------------------------

A mating scheme goes through each subpopulation and populates the subpopulations
of an offspring generation sequentially. The number of offspring in each
subpopulation is determined by the mating scheme, following the following rules:

* A simuPOP mating scheme, by default, produces an offspring generation that has
  the same subpopulation sizes as the parental generation. This does not guarantee
  a constant population size because some operators, such as a :class:`Migrator`
  and :class:`DiscardIf` can change population or subpopulation sizes.

* If fixed subpopulation sizes are given to parameter ``subPopSize``. A mating
  scheme will generate an offspring generation with specified sizes even if an
  operator has changed parental population sizes.

* A **demographic function** can be specified to parameter ``subPopSize``. This
  function should take one of the two forms ``func(gen)`` or ``func(gen, pop)``
  where ``gen`` is the current generation number and ``pop`` is the parental
  population just before mating. This function should return an array of new
  subpopulation sizes. A single number can be returned if there is only one
  subpopulation. The :mod:`simuPOP.demography` module provides a number of
  demography-related functions for complex evolutionary secenarios. **Please
  consider contributing to this module if you have implemented demographic models
  for particular populations.**

The following examples demonstrate these cases. Example :ref:`migrSize
<migrSize>` uses a default :class:`RandomMating`\ () scheme that keeps parental
subpopulation sizes. Because migration between two subpopulations are
asymmetric, the size of the first subpopulation increases at each generation,
although the overall population size keeps constant.

.. _migrSize:

**Example**: *Free change of subpopulation sizes*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[500, 1000], infoFields='migrate_to')
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     preOps=sim.Migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"%s\n" % subPopSize')
   ...     ],
   ...     gen = 3
   ... )
   [843, 657]
   [948, 552]
   [1010, 490]
   3

   now exiting runScriptInteractively...

`Download migrSize.py <migrSize.py>`_

Example :ref:`migrFixedSize <migrFixedSize>` uses the same Migrator to move
individuals between two subpopulations. Because a constant subpopulation size is
specified, the offspring generation always has 500 and 1000 individuals in its
two subpopulations. Note that operators :class:`Stat` and :class:`PyEval` are
applied both before and after mating. It is clear that subpopulation sizes
changes before mating as a result of migration, although the pre-mating
population sizes vary because of uncertainties of migration.

.. _migrFixedSize:

**Example**: *Force constant subpopulation sizes*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[500, 1000], infoFields='migrate_to')
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     preOps=[
   ...         sim.Migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"%s\n" % subPopSize')
   ...     ],
   ...     matingScheme=sim.RandomMating(subPopSize=[500, 1000]),
   ...     postOps=[
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"%s\n" % subPopSize')
   ...     ],
   ...     gen = 3
   ... )
   [843, 657]
   [500, 1000]
   [795, 705]
   [500, 1000]
   [821, 679]
   [500, 1000]
   3

   now exiting runScriptInteractively...

`Download migrFixedSize.py <migrFixedSize.py>`_

Example :ref:`demoFunc <demoFunc>` uses a demographic function to control the
subpopulation size of the offspring generation. This example implements a linear
population expansion model but arbitrarily complex demographic model can be
implemented similarly.

.. _demoFunc:

**Example**: *Use a demographic function to control population size*

::

   >>> import simuPOP as sim
   >>> def demo(gen):
   ...     return [500 + gen*10, 1000 + gen*10]
   ... 
   >>> pop = sim.Population(size=[500, 1000], infoFields='migrate_to')
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     preOps=sim.Migrator(rate=[[0.8, 0.2], [0.4, 0.6]]),
   ...     matingScheme=sim.RandomMating(subPopSize=demo),
   ...     postOps=[
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"%s\n" % subPopSize')
   ...     ],
   ...     gen = 3
   ... )
   [500, 1000]
   [510, 1010]
   [520, 1020]
   3

   now exiting runScriptInteractively...

`Download demoFunc.py <demoFunc.py>`_

If the size of the offspring generation can not be determined directly from
generation number, you can pass the parental population as parameter ``pop`` to
the demographic function. For example, Example :ref:`demoFunc1 <demoFunc1>`
implements a demographic model where a population expand at random numbers at
each generation.

.. _demoFunc1:

**Example**: *Use parental population to determine the size of offspring population*

::

   >>> import simuPOP as sim
   >>> import random
   >>> def demo(pop):
   ...     return [x + random.randint(50, 100) for x in pop.subPopSizes()]
   ... 
   >>> pop = sim.Population(size=[500, 1000], infoFields='migrate_to')
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     matingScheme=sim.RandomMating(subPopSize=demo),
   ...     postOps=[
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"%s\n" % subPopSize')
   ...     ],
   ...     gen = 3
   ... )
   [586, 1075]
   [649, 1128]
   [742, 1214]
   3

   now exiting runScriptInteractively...

`Download demoFunc1.py <demoFunc1.py>`_

In all the above examples, migration and demographic changes are introduced
manually to influence the evolution of populations. However, the demographic
changes might be driven by other factors such as natural selection so that it is
difficult to predict the size of offspring generations in advance. In this case,
you can manually remove individuals from parental (or offspring) populations
using appropriate operators.

For example, a population in Example :ref:`demoBySelection <demoBySelection>`
suffers from a sudden reduction of population size (due to perhaps a famine) at
generation 3, and a gradual reduction of population size (due to perhaps an
outburst of an infectious disease) after generation 5. The first event is
implemented using a :class:`ResizeSubPops` operator that directly shrink the
population size in half. The second event is implemented using a
:class:`MaPenetrance` and a :class:`DiscardIf` operator. The first operator
assigns affection status of each individual using a disease model that involves
individual genotype. The second operator discard all individuals that are
affected with the disease. Despite of these unfortunate events, the population
tries to expand exponentially with offspring population sizes set to 105% of
their parental populations.

.. _demoBySelection:

**Example**: *Change of  population size caused by natural selection*

::

   >>> import simuPOP as sim
   >>> def demo(pop):
   ...     return int(pop.popSize() * 1.05)
   ... 
   >>> pop = sim.Population(size=10000, loci=1)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.7, 0.3])
   ...     ],
   ...     preOps=[
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"%d %s --> " % (gen, subPopSize)'),
   ...         sim.ResizeSubPops(0, proportions=[0.5], at=2),
   ...         sim.MaPenetrance(loci=0, penetrance=[0.01, 0.2, 0.6], begin=4),
   ...         sim.DiscardIf('ind.affected()', exposeInd='ind', begin=4),
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"%s --> " % subPopSize'),
   ...     ],
   ...     matingScheme=sim.RandomMating(subPopSize=demo),
   ...     postOps=[
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"%s\n" % subPopSize')
   ...     ],
   ...     gen = 6
   ... )
   0 [10000] --> [10000] --> [10500]
   1 [10500] --> [10500] --> [11025]
   2 [11025] --> [5512] --> [5787]
   3 [5787] --> [5787] --> [6076]
   4 [6076] --> [5188] --> [5447]
   5 [5447] --> [4845] --> [5087]
   6

   now exiting runScriptInteractively...

`Download demoBySelection.py <demoBySelection.py>`_


.. _subsec_Advanced_demo_func:

Advanced use of demographic functions \*
----------------------------------------

The parental population passed to a demographic function is usually used to
determine offspring population size from parental population size. However,
because this function is called immediately before mating happens, it provides a
good opportunity for you to prepare the parental generation for mating. Such
activities could generally be done by operators, but operations related to
demographic changes could be done here. For example, Example
:ref:`advancedDemoFunc <advancedDemoFunc>` uses a demographic function to split
populations at certain generation. The advantage of this method over the use of
a :class:`SplitSubPops` operator (for example as in Example :ref:`splitByProp
<splitByProp>`) is that all demographic information presents in the same
function so you do not have to worry about changing an operator when your
demographic model changes.

.. _advancedDemoFunc:

**Example**: *Use a demographic function to split parental population*

::

   >>> import simuPOP as sim
   >>> def demo(gen, pop):
   ...     if gen < 2:
   ...         return 1000 + 100 * gen
   ...     if gen == 2:
   ...         # this happens right before mating at generation 2
   ...         size = pop.popSize()
   ...         pop.splitSubPop(0, [size // 2, size - size//2]) 
   ...     # for generation two and later
   ...     return [x + 50 * gen for x in pop.subPopSizes()]
   ... 
   >>> pop = sim.Population(1000)
   >>> pop.evolve(
   ...     preOps=[
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"Gen %d:\t%s (before mating)\t" % (gen, subPopSize)')
   ...     ],
   ...     matingScheme=sim.RandomSelection(subPopSize=demo),
   ...     postOps=[
   ...         sim.Stat(popSize=True),
   ...         sim.PyEval(r'"%s (after mating)\n" % subPopSize')
   ...     ],
   ...     gen = 5
   ... )
   Gen 0:	[1000] (before mating)	[1000] (after mating)
   Gen 1:	[1000] (before mating)	[1100] (after mating)
   Gen 2:	[1100] (before mating)	[650, 650] (after mating)
   Gen 3:	[650, 650] (before mating)	[800, 800] (after mating)
   Gen 4:	[800, 800] (before mating)	[1000, 1000] (after mating)
   5

   now exiting runScriptInteractively...

`Download advancedDemoFunc.py <advancedDemoFunc.py>`_


.. _subsec_number_of_offspring:

Determine the number of offspring during mating
-----------------------------------------------

simuPOP by default produces only one offspring per mating event. Because more
parents are involved in the production of offspring, this setting leads to
larger effective population sizes than mating schemes that produce more
offspring at each mating event. However, various situations require a larger
family size or even varying family sizes. In these cases, parameter
``numOffspring`` can be used to control the number of offspring that are
produced at each mating event. This parameter takes the following types of
inputs

* If a single number is given, ``numOffspring`` offspring are produced at each
  mating event.

* If a Python function is given, this function will be called each time when a
  mating event happens. Generation number can be passed to this function as
  parameter ``gen`` to allow different numbers of offspring at different
  generations. A python generator function can also be passed to provide an
  iterator interface to yield number of offspring for all mating events.

* If a tuple (or list) with more than one numbers is given, the first number
  must be one of ``GEOMETRIC_DISTRIBUTION``, ``POISSON_DISTRIBUTION``,
  ``BINOMIAL_DISTRIBUTION`` and ``UNIFORM_DISTRIBUTION``, with one or two
  additional parameters.

The number of offspring in the last case will then follow a specific statistical
distribution. More specifically,

* ``numOffspring=(GEOMETRIC_DISTRIBUTION, p)``: The number of offspring for each
  mating event follows a geometric distribution with mean :math:`1/p` and variance
  :math:`\left(1-p\right)/p^{2}`:

  .. math::

      \mbox{Pr}\left(k\right)=p\left(1-p\right)^{k-1}\;\textrm{ for }k\geq1

* ``numOffspring=(POISSON_DISTRIBUTION, p)``: The number of offspring for each
  mating event follows a Poisson distribution with mean :math:`p` and variance
  :math:`p`. The distribution is

  .. math::

      \mbox{Pr}\left(k\right)=\frac{p^{k}e^{-p}}{k!}\;\textrm{ for }k\geq0

  Note that, however, because families with zero offspring are ignored, the
  distribution of the observed number of offspring (excluding zero) follows a
  zero-truncated Poission distribution with probability

  .. math::

      \mbox{Pr}\left(k\right)=\frac{p^{k}e^{-p}}{k!\left(1-e^{-p}\right)}\;\textrm{ for }k\geq1

  The mean number of offspring is therefore :math:`\frac{1}{1-e^{-p}}p`, which is
  2.31 for :math:`p=2`.

* ``numOffspring=(BINOMIAL_DISTRIBUTION, p, n):``\ The number of offspring for
  each mating event follows a Binomial distribution with mean :math:`np` and
  variance :math:`np\left(1-p\right)`.

  .. math::

      \mbox{Pr}\left(k\right)=\frac{n!}{k!\left(n-k\right)!}p^{k}\left(1-p\right)^{n-k}\;\textrm{ for }n\geq k\geq0

  Because families with zero offspring are ignored, the distribution of the
  observed number of offspring (excluding zero) follows a zero-truncated Bionimial
  distribution, with mean number of offspring being
  :math:`\frac{np}{\left(1-p\right)^{n}}`.

* ``numOffspring=(UNIFORM_DISTRIBUTION, a, b):`` The number of offspring for
  each mating event follows a discrete uniform distribution with lower bound
  :math:`a` and upper bound :math:`b`.

  .. math::

      \mbox{Pr}\left(k\right)=\frac{1}{b-a+1}\;\textrm{ for }b\geq k\geq a

  The lower bound of this distribution can be ``0`` but is identical to the case
  with :math:`a=1`.

Example :ref:`numOff <numOff>` demonstrates how to use parameter
``numOffspring``. In this example, a function ``checkNumOffspring`` is defined.
It takes a mating scheme as its input parameter and use it to evolve a
population with 30 individuals. After evolving a population for one generation,
parental indexes are used to identify siblings, and then the number of offspring
per mating event.

.. _numOff:

**Example**: *Control the number of offspring per mating event.*

::

   >>> import simuPOP as sim
   >>> def checkNumOffspring(numOffspring, ops=[]):
   ...     '''Check the number of offspring for each family using
   ...        information field father_idx
   ...     '''
   ...     pop = sim.Population(size=[30], loci=1, infoFields=['father_idx', 'mother_idx'])
   ...     pop.evolve(
   ...         initOps=[
   ...             sim.InitSex(),
   ...             sim.InitGenotype(freq=[0.5, 0.5]),
   ...         ],
   ...         matingScheme=sim.RandomMating(ops=[
   ...             sim.MendelianGenoTransmitter(),
   ...             sim.ParentsTagger(),
   ...             ] + ops,
   ...             numOffspring=numOffspring),
   ...         gen=1)
   ...     # get the parents of each offspring
   ...     parents = [(x, y) for x, y in zip(pop.indInfo('mother_idx'),
   ...         pop.indInfo('father_idx'))]
   ...     # Individuals with identical parents are considered as siblings.
   ...     famSize = []
   ...     lastParent = (-1, -1)
   ...     for parent in parents:
   ...         if parent == lastParent:
   ...             famSize[-1] += 1
   ...         else:
   ...             lastParent = parent
   ...             famSize.append(1)
   ...     return famSize
   ... 
   >>> # Case 1: produce the given number of offspring
   >>> checkNumOffspring(numOffspring=2)
   [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
   >>> # Case 2: Use a Python function
   >>> import random
   >>> def func(gen):
   ...     return random.randint(5, 8)
   ... 
   >>> checkNumOffspring(numOffspring=func)
   [5, 7, 5, 5, 6, 2]
   >>> # Case 3: A geometric distribution
   >>> checkNumOffspring(numOffspring=(sim.GEOMETRIC_DISTRIBUTION, 0.3))
   [3, 1, 1, 3, 4, 1, 1, 1, 2, 1, 1, 4, 6, 1]
   >>> # Case 4: A Possition distribution
   >>> checkNumOffspring(numOffspring=(sim.POISSON_DISTRIBUTION, 1.6))
   [2, 2, 1, 5, 3, 3, 1, 1, 2, 3, 3, 2, 2]
   >>> # Case 5: A Binomial distribution
   >>> checkNumOffspring(numOffspring=(sim.BINOMIAL_DISTRIBUTION, 0.1, 10))
   [1, 4, 1, 1, 2, 1, 1, 3, 1, 1, 1, 3, 2, 2, 1, 1, 1, 2, 1]
   >>> # Case 6: A uniform distribution
   >>> checkNumOffspring(numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 6))
   [4, 4, 2, 6, 6, 2, 2, 2, 2]
   >>> # Case 7: With selection on offspring
   >>> checkNumOffspring(numOffspring=8,
   ...     ops=[sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):0.8, (1,1):0.5})])
   [8, 5, 7, 6, 4]

   now exiting runScriptInteractively...

`Download numOff.py <numOff.py>`_

However, **the actual number of offspring can be less than specified because
offspring can be discarded during mating.** More specifically, if any during-
mating generator, such as a during-mating selector, returns ``False`` during the
production of offspring, the offspring will be discarded so the total number of
offspring will be reduced. This is the case in the seventh case of Example
:ref:`numOff <numOff>` where offspring with certain genotypes have lower
probabilities to survive. If you would like to control size of families in the
presence of natural selection, you could set a larger ``numOffspring`` use a
:class:`OffspringTagger` to mark the index of offspring, and discard offspring
conditionally using operator :class:`DiscardIf` . Please refer to example
:ref:`OffspringTagger <OffspringTagger>` for details.


Dynamic population size determined by number of offspring \*
------------------------------------------------------------

What we have described so far requires you to determine the size of offspring
population in advance. Each mating event produces a number of offspring that is
determined by parameter ``NumOffspring``. The mating process stops when the
offspring population is filled. This works for most scenarios but there are
cases where the offspring population size is determined dynamically from a fixed
number of mating events with random number of offspring. For example, you might
design a mating scheme where all males in a population mate only once and
produce random number of offspring.

These kind of mating schemes can be simulated using a demographic model that
calculates offspring population size from pre-simulated number of offspring for
each family. More specifically, we

* Define a demogrphic function (model) that will be called before mating
  happens.

* This function determines and save the number of offspring for each mating
  event, and return the total number of offspring as offspring population size.

* Pass a function or generator to parameter numOffspring to pass pre-determined
  number of offspring. This function will be called each time when number of
  offspring is needed.

The number of offspring could be saved and retrieved as global variable but a
more clever method is to store the numbers of offspring in a demographic model
(class). Example :ref:`dynamicNumOff <dynamicNumOff>` demonstrates this method
by implementing a demographic model that simulate, save, and return the number
of offspring. Note that although we determine the number of mating events from
number of males in the parental population, a random mating scheme will choose
parents with replacement so it is likely that some parents will be chosen
multiple times while some others are not chosen at all. Please refer to section
"Non-random and customized mating schemes" to learn how to define a mating
scheme that picks parents without replacement.

.. _dynamicNumOff:

**Example**: *Dynamic population size determined by number of offspring*

::

   >>> import simuPOP as sim
   >>> 
   >>> import random
   >>> 
   >>> class RandomNumOff:
   ...     # a demographic model
   ...     def __init__(self):
   ...         self.numOff = []
   ...     
   ...     def getNumOff(self):
   ...         # return the pre-simulated number of offspring as a generator function
   ...         for item in self.numOff:
   ...             yield item
   ...     
   ...     def __call__(self, pop):
   ...         # define __call__ so that a RandomNumOff object is callable.
   ...         #
   ...         # Each male produce from 1 to 3 offspring. For large population, get the
   ...         # number of males instead of checking the sex of each individual
   ...         self.numOff = [random.randint(1, 3) for ind in pop.individuals() if ind.sex() == sim.MALE]
   ...         # return the total population size
   ...         print('{} mating events with number of offspring {}'.format(len(self.numOff), self.numOff))
   ...         return sum(self.numOff)
   ... 
   >>> 
   >>> pop = sim.Population(10)
   >>> 
   >>> # create a demogranic model
   >>> numOffModel = RandomNumOff()
   >>> 
   >>> pop.evolve(
   ...     preOps=sim.InitSex(),
   ...     matingScheme=sim.RandomMating(
   ...         # the model will be called before mating to deteremine
   ...         # family and population size
   ...         subPopSize=numOffModel,
   ...         # the getNumOff function (generator) returns number of offspring
   ...         # for each mating event
   ...         numOffspring=numOffModel.getNumOff
   ...     ),
   ...     gen=3
   ... )
   5 mating events with number of offspring [3, 2, 2, 3, 3]
   6 mating events with number of offspring [3, 2, 3, 1, 2, 3]
   6 mating events with number of offspring [2, 1, 1, 2, 3, 2]
   3
   >>> 

   now exiting runScriptInteractively...

`Download dynamicNumOff.py <dynamicNumOff.py>`_


.. _subsec_offspring_sex:

Determine sex of offspring
--------------------------

Because sex can influence how genotypes are transmitted (e.g. sex chromosomes,
haplodiploid population), simuPOP determines offspring sex before it passes an
offspring to a *genotype transmitter* (during-mating operator) to transmit
genotype from parents to offspring. The default ``sexMode`` in almost all mating
schemes is ``RandomSex``, in which case simuPOP assign ``Male`` or ``Female`` to
offspring with equal probability.

Other sex determination methods are also available:

* ``sexMode=RANDOM_SEX``: Sex is determined randomly, with equal probability for
  ``MALE`` and ``FEMALE``. This is the default mode for sexual mating schemes such
  as random mating.

* ``sexMode=NO_SEX``: Sex is not simulated so everyone is ``MALE``. This is the
  default mode for asexual mating schemes.

* ``sexMode=(PROB_OF_MALES, prob)``: Produce males with given probability.

* ``sexMode=(NUM_OF_MALES, n)``: The first ``n`` offspring in each family will
  be ``Male``. If the number of offspring at a mating event is less than or equal
  to ``n``, all offspring will be male.

* ``sexMode=(NUM_OF_FEMALES, n)``: The first ``n`` offspring in each family will
  be ``Female``.

* ``sexMode=(SEQUENCE_OF_SEX, s1, s2, ...)``: Use sequence ``s1``, ``s2``, ...
  for offspring in each mating event.

* ``sexMode=(GLOBAL_SEQUENCE_OF_SEX, s1, s2, ...)``: Use sequence ``s1``,
  ``s2``, ... for all offspring in a subpopulation. Because other mode of sex
  determination works within each mating event, this is the only way to ensure
  proportion of sex in a subpopulation. For example, ``(GLOBAL_SEQUENCE_OF_SEX,
  MALE, FEMALE)`` will gives ``MALE`` and ``FEMALE`` iteratively to all offspring,
  making sure there are equal number of males and females (if there are even
  number of offspring).

* ``sexMode=func`` or ``sexMode=generator_func``: In this last case, a Python
  function or a Python generator function can be specified to provide sex to each
  offspring. The function is called whenever an offspring is created. The
  generator function is called for each subpopulation, and provides an iterator
  that provides sex for all offspring in a subpopulation.

``NumOfMales`` and ``NumOfFemales`` are useful in theoretical studies where the
sex ratio of a population needs to be controlled strictly, or in special mating
schemes, usually for animal populations, where only a certain number of male or
female Individuals are allowed in a family. It worth noting that a genotype
transmitter can override specified offspring sex. This is the case for
:class:`CloneGenoTransmitter` where an offspring inherits both genotype and sex
from his/her parent.

Example :ref:`sexMode <sexMode>` demonstrates how to use parameter ``sexMode``.
In this example, a function ``checkSexMode`` is defined. It takes a mating
scheme as its input parameter and use it to evolve a population with 40
individuals. After evolving a population for one generation, sexes of all
offspring are returned as a string.

.. _sexMode:

**Example**: *Determine the sex of offspring*

::

   >>> import simuPOP as sim
   >>> def checkSexMode(ms):
   ...     '''Check the assignment of sex to offspring'''
   ...     pop = sim.Population(size=[40])
   ...     pop.evolve(initOps=sim.InitSex(), matingScheme=ms, gen=1)
   ...     # return individual sex as a string
   ...     return ''.join(['M' if ind.sex() == sim.MALE else 'F' for ind in pop.individuals()])
   ... 
   >>> # Case 1: sim.NO_SEX (all male, sim.RandomMating will not continue)
   >>> checkSexMode(sim.RandomMating(sexMode=sim.NO_SEX))
   'MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM'
   >>> # Case 2: sim.RANDOM_SEX (sim.Male/Female with probability 0.5)
   >>> checkSexMode(sim.RandomMating(sexMode=sim.RANDOM_SEX))
   'MFFFFFFMFFFMFFFMMFFMFFMMMFFMMFMFFFFFFFMF'
   >>> # Case 3: sim.PROB_OF_MALES (Specify probability of male)
   >>> checkSexMode(sim.RandomMating(sexMode=(sim.PROB_OF_MALES, 0.8)))
   'MMFMMFFMMFFMMMMMMMMMFMMFFMMMMMMMMMMMMMMM'
   >>> # Case 4: sim.NUM_OF_MALES (Specify number of male in each family)
   >>> checkSexMode(sim.RandomMating(numOffspring=3, sexMode=(sim.NUM_OF_MALES, 1)))
   'MFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFMFFM'
   >>> # Case 5: sim.NUM_OF_FEMALES (Specify number of female in each family)
   >>> checkSexMode(sim.RandomMating(
   ...     numOffspring=(sim.UNIFORM_DISTRIBUTION, 4, 6),
   ...     sexMode=(sim.NUM_OF_FEMALES, 2))
   ... )
   'FFMMFFMMMFFMMFFMMMMFFMMFFMMMMFFMMFFMMFFM'
   >>> # Case 6: sim.SEQUENCE_OF_SEX
   >>> checkSexMode(sim.RandomMating(
   ...     numOffspring=4, sexMode=(sim.SEQUENCE_OF_SEX, sim.MALE, sim.FEMALE))
   ... )
   'MFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMF'
   >>> # Case 7: sim.GLOBAL_SEQUENCE_OF_SEX
   >>> checkSexMode(sim.RandomMating(
   ...     numOffspring=3, sexMode=(sim.GLOBAL_SEQUENCE_OF_SEX, sim.MALE, sim.FEMALE))
   ... )
   'MFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMF'
   >>> # Case 8: A generator function
   >>> def sexFunc():
   ...     i = 0
   ...     while True:
   ...         i += 1
   ...         if i % 2 == 0:
   ...             yield sim.MALE
   ...         else:
   ...             yield sim.FEMALE
   ... 
   >>> checkSexMode(sim.RandomMating(numOffspring=3, sexMode=sexFunc))
   'FMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFMFM'

   now exiting runScriptInteractively...

`Download sexMode.py <sexMode.py>`_


Monogamous mating
-----------------

Monogamous mating (monogamy) in simuPOP refers to mating schemes in which each
parent mates only once. In an asexual setting, this implies parents are chosen
without replacement. In sexual mating schemes, this means that parents are
chosen without replacement, they have only one spouse during their life time so
that all siblings have the same parents (no half-sibling).

simuPOP provides a diploid sexual monogamous mating scheme
:class:`MonogamousMating`. However, without careful planning, this mating scheme
can easily stop working due to the lack of parents. For example, if a population
has 40 males and 55 females, only 40 successful mating events can happen and
result in 40 offspring in the offspring generation. :class:`MonogamousMating`
will exit if the offspring generation is larger than 40.

Example :ref:`monogamous <monogamous>` demonstrates one scenario of using a
monogamous mating scheme where sex of parents and offspring are strictly
specified so that parents will not be exhausted. The sex initializer
:class:`InitSex` assigns exactly 10 males and 10 females to the initial
population. Because of the use of ``numOffspring=2, sexMode=(NUM_OF_MALES, 1)``,
each mating event will produce exactly one male and one female. Unlike a random
mating scheme that only about 80% of parents are involved in the production of
an offspring population with the same size, this mating scheme makes use of all
parents.

.. _monogamous:

**Example**: *Sexual monogamous mating*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(20, infoFields=['father_idx', 'mother_idx'])
   >>> pop.evolve(
   ...     initOps=sim.InitSex(sex=(sim.MALE, sim.FEMALE)),
   ...     matingScheme=sim.MonogamousMating(
   ...         numOffspring=2,
   ...         sexMode=(sim.NUM_OF_MALES, 1),
   ...         ops=[
   ...             sim.MendelianGenoTransmitter(),
   ...             sim.ParentsTagger(),
   ...         ],
   ...     ),
   ...     gen = 5
   ... )
   5
   >>> [ind.sex() for ind in pop.individuals()]
   [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
   >>> [int(ind.father_idx) for ind in pop.individuals()]
   [16, 16, 2, 2, 4, 4, 8, 8, 0, 0, 14, 14, 10, 10, 12, 12, 18, 18, 6, 6]
   >>> [int(ind.mother_idx) for ind in pop.individuals()]
   [13, 13, 17, 17, 1, 1, 15, 15, 19, 19, 9, 9, 3, 3, 5, 5, 7, 7, 11, 11]
   >>> # count the number of distinct parents
   >>> len(set(pop.indInfo('father_idx')))
   10
   >>> len(set(pop.indInfo('mother_idx')))
   10

   now exiting runScriptInteractively...

`Download monogamous.py <monogamous.py>`_


Polygamous mating
-----------------

In comparison to monogamous mating, parents in a polygamous mate with more than
one spouse during their life-cycle. Both *polygany* (one man has more than one
wife) and ``polyandry`` (one woman has more than one husband) are supported.

Other than regular parameters such as ``numOffspring``, mating scheme
``PolygamousMating`` accepts parameters ``polySex`` (default to ``Male``) and
``polyNum`` (default to 1). During mating, an individual with ``polySex`` is
selected and then mate with ``polyNum`` randomly selected spouse. Example
:ref:`polygamous <polygamous>` demonstrates the use of this mating schemes. Note
that this mating scheme support natural selection, but does not yet handle
varying ``polyNum`` and selection of parents without replacement.

.. _polygamous:

**Example**: *Sexual polygamous mating*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(100, infoFields=['father_idx', 'mother_idx'])
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     matingScheme=sim.PolygamousMating(polySex=sim.MALE, polyNum=2,
   ...         ops=[sim.ParentsTagger(),
   ...             sim.MendelianGenoTransmitter()],
   ...     ),
   ...     gen = 5
   ... )
   5
   >>> [int(ind.father_idx) for ind in pop.individuals()][:20]
   [67, 67, 42, 42, 91, 91, 25, 25, 65, 65, 47, 47, 18, 18, 16, 16, 96, 96, 57, 57]
   >>> [int(ind.mother_idx) for ind in pop.individuals()][:20]
   [58, 58, 58, 0, 68, 32, 37, 89, 6, 85, 12, 58, 36, 12, 66, 44, 51, 85, 60, 29]

   now exiting runScriptInteractively...

`Download polygamous.py <polygamous.py>`_


Asexual random mating
---------------------

Mating scheme :class:`RandomSelection` implements an asexual random mating
scheme. It randomly select parents from a parental population (with replacement)
and copy them to an offspring generation. Both genotypes and sex of the parents
are copied because genotype and sex are sometimes related. This mating scheme
can be used to simulate the evolution of haploid sequences in a standard haploid
Wright-Fisher model.

Example :ref:`RandomSelection <RandomSelection>` applies a
:class:`RandomSelection` mating scheme to a haploid population with 100
sequences. A ``parentTagger`` is used to track the parent of each individual.
Although sex information is not used in this mating scheme, Individual sexes are
initialized and passed to offspring.

.. _RandomSelection:

**Example**: *Asexual random mating*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(100, ploidy=1, loci=[5, 5], ancGen=1,
   ...     infoFields='parent_idx')
   >>> pop.evolve(
   ...     initOps=sim.InitGenotype(freq=[0.3, 0.7]),
   ...     matingScheme=sim.RandomSelection(ops=[
   ...         sim.ParentsTagger(infoFields='parent_idx'),
   ...         sim.CloneGenoTransmitter(),
   ...     ]),
   ...     gen = 5
   ... )
   5
   >>> ind = pop.individual(0)
   >>> par = pop.ancestor(ind.parent_idx, 1)
   >>> print(ind.sex(), ind.genotype())
   1 [1, 1, 0, 1, 1, 0, 1, 0, 0, 0]
   >>> print(par.sex(), par.genotype())
   1 [1, 1, 0, 0, 1, 1, 1, 1, 0, 1]

   now exiting runScriptInteractively...

`Download RandomSelection.py <RandomSelection.py>`_


Mating in haplodiploid populations
----------------------------------

Male individuals in a haplodiploid population are derived from unfertilized eggs
and thus have only one set of chromosomes. Mating in such a population is
handled by a special mating scheme called ``haplodiplodMating``. This mating
scheme chooses a pair of parents randomly and produces some offspring. It
transmit maternal chromosomes and paternal chromosomes (the only copy) to female
offspring, and only maternal chromosomes to male offspring. Example
:ref:`HaplodiploidMating <HaplodiploidMating>` demonstrates how to use this
mating scheme. It uses three initializers because sex has to be initialized
before two other intializers can initialize genotype by sex.

.. _HaplodiploidMating:

**Example**: *Random mating in haplodiploid populations*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(10, ploidy=sim.HAPLODIPLOID, loci=[5, 5],
   ...     infoFields=['father_idx', 'mother_idx'])
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[0]*10, subPops=[(0, 'Male')]),
   ...         sim.InitGenotype(genotype=[1]*10+[2]*10, subPops=[(0, 'Female')])
   ...     ],
   ...     preOps=sim.Dumper(structure=False),
   ...     matingScheme=sim.HaplodiploidMating(
   ...         ops=[sim.HaplodiploidGenoTransmitter(), sim.ParentsTagger()]),
   ...     postOps=sim.Dumper(structure=False),
   ...     gen = 1
   ... )
   SubPopulation 0 (), 10 Individuals:
      0: FU 11111 11111 | 22222 22222 |  0 0
      1: FU 11111 11111 | 22222 22222 |  0 0
      2: MU 00000 00000 | _____ _____ |  0 0
      3: MU 00000 00000 | _____ _____ |  0 0
      4: MU 00000 00000 | _____ _____ |  0 0
      5: MU 00000 00000 | _____ _____ |  0 0
      6: MU 00000 00000 | _____ _____ |  0 0
      7: FU 11111 11111 | 22222 22222 |  0 0
      8: FU 11111 11111 | 22222 22222 |  0 0
      9: FU 11111 11111 | 22222 22222 |  0 0

   SubPopulation 0 (), 10 Individuals:
      0: MU 11111 11111 | _____ _____ |  4 9
      1: MU 11111 22222 | _____ _____ |  4 8
      2: MU 22222 11111 | _____ _____ |  6 8
      3: MU 22222 11111 | _____ _____ |  3 8
      4: MU 22222 22222 | _____ _____ |  2 8
      5: MU 22222 22222 | _____ _____ |  6 9
      6: FU 22222 22222 | 00000 00000 |  2 1
      7: FU 22222 22222 | 00000 00000 |  2 1
      8: FU 22222 22222 | 00000 00000 |  3 9
      9: FU 11111 11111 | 00000 00000 |  5 8

   1

   now exiting runScriptInteractively...

`Download HaplodiploidMating.py <HaplodiploidMating.py>`_

Note that this mating scheme does not support recombination and the standard
Recombinator does not work with haplodiploid populations. Please refer to the
next Chapter for how to define a customized genotype transmitter to handle such
a situation.


Self-fertilization
------------------

Some plant populations evolve through self-fertilization. That is to say, a
parent fertilizes with itself during the production of offspring (seeds). In a
:class:`SelfMating` mating scheme, parents are chosen randomly (one at a time),
and are used twice to produce two homologous sets of offspring chromosomes. The
standard Recombinator can be used with this mating scheme. Example
:ref:`SelfMating <SelfMating>` initializes each chromosome with different
alleles to demonstrate how these alleles are transmitted in this population.

.. _SelfMating:

**Example**: *Selfing mating scheme*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(20, loci=8)
   >>> # every chromosomes are different. :-)
   >>> for idx, ind in enumerate(pop.individuals()):
   ...     ind.setGenotype([idx*2], 0)
   ...     ind.setGenotype([idx*2+1], 1)
   ... 
   >>> pop.evolve(
   ...     matingScheme=sim.SelfMating(ops=sim.Recombinator(rates=0.01)),
   ...     gen = 1
   ... )
   1
   >>> sim.dump(pop, width=3, structure=False, max=10)
   SubPopulation 0 (), 20 Individuals:
      0: FU  36 36 36 36 36 36 36 36 |  36 36 36 36 36 36 36 36 
      1: FU   6  6  6  6  6  6  6  6 |   7  7  7  7  7  7  7  7 
      2: MU  33 33 33 33 33 33 33 33 |  33 33 33 33 33 33 33 33 
      3: MU  22 22 22 22 22 23 23 23 |  22 22 22 22 22 22 22 22 
      4: FU  27 27 27 27 27 27 27 27 |  27 27 27 27 27 27 27 27 
      5: MU  15 15 15 15 15 15 15 15 |  15 15 15 15 15 15 15 15 
      6: MU  35 35 35 35 34 34 34 34 |  34 34 34 34 34 34 34 34 
      7: FU  11 11 11 11 11 11 11 11 |  10 10 10 10 10 10 10 10 
      8: MU  11 11 11 11 11 11 11 11 |  11 11 11 11 11 11 11 11 
      9: FU  24 24 24 24 24 24 24 24 |  25 25 25 25 25 25 25 25 


   now exiting runScriptInteractively...

`Download SelfMating.py <SelfMating.py>`_


Heterogeneous mating schemes \*
-------------------------------

Different groups of individuals in a population may have different mating
patterns. For example, individuals with different properties can have varying
fecundity, represented by different numbers of offspring generated per mating
event. This can be extended to aged populations in which only adults (may be
defined by age > 20 and age < 40) can produce offspring, where other individuals
will either be copied to the offspring generation or die.

A heterogeneous mating scheme (:class:`HeteroMating`) accepts a list of mating
schemes that are applied to different subpopulation or virtual subpopulations.
If multiple mating schemes are applied to the same subpopulation, each of them
only populate part of the offspring subpopulation. This is illustrated in Figure
:ref:`fig_heterogenous_mating <fig_heterogenous_mating>`.

**Figure**: *Illustration of a heterogeneous mating scheme*

.. _fig_heterogenous_mating:

.. figure:: /Users/bpeng1/simuPOP/simuPOP/doc/figures/MatingScheme.png
   :width: 680


A heterogeneous mating scheme that applies homogeneous mating schemes MS0,
MS0.0, MS0.1, MS1, MS2.0 and MS2.1 to subpopulation 0, the first and second
virtual subpopulation in subpopulation 0, subpopulation 1, the first and second
virtual subpopulation in subpopulation 2, respectively. Note that VSP 0 and 1 in
subpopulation 0 overlap, and do not add up to subpopulation 0.

For example, Example :ref:`hateroMatingSP <hateroMatingSP>` applies two random
mating schemes to two subpopulations. The first mating scheme produces two
offspring per mating event, and the second mating scheme produces four.

.. _hateroMatingSP:

**Example**: *Applying different mating schemes to different subpopulations*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[1000, 1000], loci=2,
   ...     infoFields=['father_idx', 'mother_idx'])
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     matingScheme=sim.HeteroMating([
   ...         sim.RandomMating(numOffspring=2, subPops=0,
   ...             ops=[sim.MendelianGenoTransmitter(), sim.ParentsTagger()]
   ...         ),
   ...         sim.RandomMating(numOffspring=4, subPops=1,
   ...             ops=[sim.MendelianGenoTransmitter(), sim.ParentsTagger()]
   ...         )
   ...     ]),
   ...     gen=10
   ... )
   10
   >>> [int(ind.father_idx) for ind in pop.individuals(0)][:10]
   [134, 134, 451, 451, 780, 780, 443, 443, 457, 457]
   >>> [int(ind.father_idx) for ind in pop.individuals(1)][:10]
   [1978, 1978, 1978, 1978, 1582, 1582, 1582, 1582, 1322, 1322]

   now exiting runScriptInteractively...

`Download HeteroMatingSP.py <HeteroMatingSP.py>`_

The real power of heterogeneous mating schemes lies on their ability to apply
different mating schemes to different virtual subpopulations. For example, due
to different micro-environmental factors, plants in the same population may
exercise both self and cross-fertilization. Because of the randomness of such
environmental factors, it is difficult to divide a population into self and
cross-mating subpopulations. Applying different mating schemes to groups of
individuals in the same subpopulation is more appropriate.

Example :ref:`hateroMatingVSP <hateroMatingVSP>` applies two mating schemes to
two VSPs defined by proportions of individuals. In this mating scheme, 20% of
individuals go through self-mating and 80% of individuals go through random
mating. This can be seen from the parental indexes of individuals in the
offspring generation: individuals whose ``mother_idx`` are ``-1`` are
genetically only derived from their fathers.

It might be surprising that offspring resulted from two mating schemes mix with
each other so the same VSPs in the next generation include both selfed and
cross-fertilized offspring. If this not desired, you can set parameter
``shuffleOffspring=False`` in :class:`HeteroMating`\ (). Because the number of
offspring that are produced by each mating scheme is proportional to the size of
parental (virtual) subpopulation, the first 20% of individuals that are produced
by self-fertilization will continue to self-fertilize.

.. _hateroMatingVSP:

**Example**: *Applying different mating schemes to different virtual subpopulations*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[1000], loci=2,
   ...     infoFields=['father_idx', 'mother_idx'])
   >>> pop.setVirtualSplitter(sim.ProportionSplitter([0.2, 0.8]))
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     matingScheme=sim.HeteroMating(matingSchemes=[
   ...         sim.SelfMating(subPops=[(0, 0)],
   ...             ops=[sim.SelfingGenoTransmitter(), sim.ParentsTagger()]
   ...         ),
   ...         sim.RandomMating(subPops=[(0, 1)],
   ...             ops=[sim.SelfingGenoTransmitter(), sim.ParentsTagger()]
   ...         )
   ...     ]),
   ...     gen = 10
   ... )
   10
   >>> [int(ind.father_idx) for ind in pop.individuals(0)][:15]
   [789, 666, 145, 125, 681, 183, 727, 308, 392, 11, 183, 223, 208, 29, 309]
   >>> [int(ind.mother_idx) for ind in pop.individuals(0)][:15]
   [370, 272, -1, 520, 121, 91, 220, 519, 101, 271, -1, 263, 663, -1, 286]

   now exiting runScriptInteractively...

`Download HeteroMatingVSP.py <HeteroMatingVSP.py>`_

Because there is no restriction on the choice of VSPs, mating schemes can be
applied to overlapped (virtual) subpopulations. For example,

::

   HeteroMating(
       matingSchemes = [
           SelfMating(subPops=[(0, 0)]),
           RandomMating(subPops=0)
           ]
   )

will apply SelfMating to the first 20% individuals, and RandomMating will be
applied to all individuals. Similarly,  ::

   HeteroMating(
       matingSchemes = [
           SelfMating(subPops=0),
           RandomMating(subPops=0)
           ]
   )

will allow all individuals to be involved in both :class:`SelfMating` and
:class:`RandomMating`.

This raises the question of how many offspring each mating scheme will produce.
By default, the number of offspring produced will be proportional to the size of
parental (virtual) subpopulations. In the last example, because both mating
schemes are applied to the same subpopulation, half of all offspring will be
produced by selfing and the other half will be produced by random mating.

This behavior can be changed by a weighting scheme controlled by parameter
``weight`` of each homogeneous mating scheme. Briefly speaking, a positive
weight will be compared against other mating schemes. a negative weight is
considered proportional to the existing (virtual) subpopulation size. Negative
weights are considered before positive or zero weights.

This weighting scheme is best explained by an example. Assuming that there are
three mating schemes working on the same parental subpopulation

* Mating scheme A works on the whole subpopulation of size 1000

* Mating scheme B works on a virtual subpopulation of size 500

* Mating scheme C works on another virtual subpopulation of size 800

Assuming the corresponding offspring subpopulation has :math:`N` individuals,

* If all weights are 0, the offspring subpopulation is divided in proportion to
  parental (virtual) subpopulation sizes. In this example, the mating schemes will
  produce :math:`\frac{10}{23}N`, :math:`\frac{5}{23}N`, :math:`\frac{8}{23}N`
  individuals respectively.

* If all weights are negative, they are multiplied to their parental (virtual)
  subpopulation sizes. For example, weight (-1, -2, -0.5) will lead to sizes
  (1000, 1000, 400) in the offspring subpopulation. If :math:`N\ne2400` in this
  case, an error will be raised.

* If all weights are positive, the number of offspring produced from each mating
  scheme is proportional to these weights. For example, weights (1, 2, 3) will
  lead to :math:`\frac{1}{6}N`, :math:`\frac{2}{6}N`, :math:`\frac{1}{3}N`
  individuals respectively. In this case, 0 weights will produce no offspring.

* If there are mixed positive and negative weights, the negative weights are
  processed first, and the rest of the individuals are divided using non-negative
  weights. For example, three mating schemes with weights (-0.5, 2, 3) will
  produce 500, :math:`\frac{2}{5}\left(N-500\right)`,
  :math:`\frac{3}{5}\left(N-500\right)` individuals respectively.

The last case is demonstrated in Example :ref:`HeteroMatingWeight
<HeteroMatingWeight>` where three random mating schemes are applied to
subpopulation ``0``, virtual subpopulation\ ``(0, 0)`` and virtual subpopulation
``(0, 1)``, with weights ``-``\ 0.5, ``2``, and ``3`` respectively. This example
uses an advanced features that will be described in the next section. Namely,
three during-mating Python operators are passed to each mating scheme to mark
their offspring with different numbers.

.. _HeteroMatingWeight:

**Example**: *A weighting scheme used by heterogeneous mating schemes.*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[1000], loci=2,
   ...     infoFields='mark')
   >>> pop.setVirtualSplitter(sim.RangeSplitter([[0, 500], [200, 1000]]))
   >>> 
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     matingScheme=sim.HeteroMating([
   ...         sim.RandomMating(subPops=0, weight=-0.5,
   ...             ops=[sim.InfoExec('mark=0'), sim.MendelianGenoTransmitter()]),
   ...         sim.RandomMating(subPops=[(0, 0)], weight=2,
   ...             ops=[sim.InfoExec('mark=1'), sim.MendelianGenoTransmitter()]),
   ...         sim.RandomMating(subPops=[(0, 1)], weight=3,
   ...             ops=[sim.InfoExec('mark=2'), sim.MendelianGenoTransmitter()])
   ...     ]),
   ...     gen = 10
   ... )
   10
   >>> marks = list(pop.indInfo('mark'))
   >>> marks.count(0.)
   500
   >>> marks.count(1.)
   200
   >>> marks.count(2.)
   300

   now exiting runScriptInteractively...

`Download HeteroMatingWeight.py <HeteroMatingWeight.py>`_

As a special case that can be quite annoying during the simulation of small
populations, a (virtual) subpopulation can have no male and/or female. If the
parental (virtual) subpopulation is empty, it will produce no offspring
regardless of its weight. However, if the parental (virtual) subpopulation is
not empty, it will be expected to produce some offspring, which is not possible
if a sexual mating scheme is used. In this case, you can use a parameter
``weightBy`` to specify how parental (virtual) population sizes are calculated.
This parameter accepts values ``ANY_SEX`` (default), ``MALE_ONLY``,
``FEMALE_ONLY``, ``PAIR_ONLY``, and use all individuals, number of male
individuals, number of female individuals, and number of male/female pairs
(basically the less of numbers of males and females) as the size of parental
(virtual) subpopulation, respectively. When ``weightBy=PAIR_ONLY`` is used,
parental (virtual) subpopulations with only males or females will appear to be
empty and produce no offspring. Note that in this mode (also ``MALE_ONLY``,
``FEMALE_ONLY``), the perceived parental population sizes are no longer the
actual parental population sizes so you might need to adjust parameter
``weight`` (e.g. ``weight=-2``) to produce correct number of offspring.


Conditional mating schemes
--------------------------

A :class:`ConditionalMating` mating scheme allows you to apply different mating
schemes to populations with different properties. The condition can be a
constant (True or False), an expression that will be evaluated in the local
namspace of the parental population, or a function that can take parental
population as its input paramter (with parameter name ``pop``).

Using variable ``rep`` and ``gen`` in the local namespace of the parental
population, we can use this mating scheme to apply different mating schemes to
different replicates and/or at different generations. For example,
:ref:`matingSchemeByRepAndGen <matingSchemeByRepAndGen>` simulates the evolution
of three replicates. The first replicate uses regular mating scheme, the third
replicate uses a mating scheme that produces 70% of males, and the second
replicate do this only for the first 5 generations. Because there are three
cases, a nested :class:`ConditionalMating` is used.

.. _matingSchemeByRepAndGen:

**Example**: *Apply different mating schemes for different replicates at different generations*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(sim.Population(1000, loci=[10]), rep=3)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     matingScheme=sim.ConditionalMating('rep == 0', 
   ...         # the first replicate use standard random mating
   ...         sim.RandomMating(),
   ...         sim.ConditionalMating('rep == 1 and gen >= 5',
   ...             # the second replicate produces more males for the first 5 generations
   ...             sim.RandomMating(),
   ...             # the last replicate produces more males all the time
   ...             sim.RandomMating(sexMode=(sim.PROB_OF_MALES, 0.7))
   ...             )
   ...         ),
   ...     postOps=[
   ...         sim.Stat(numOfMales=True),
   ...         sim.PyEval("'gen=%d' % gen", reps=0),
   ...         sim.PyEval(r"'\t%d' % numOfMales"),
   ...         sim.PyOutput('\n', reps=-1)
   ...     ],        
   ...     gen=10
   ... )
   gen=0	477	686	718
   gen=1	477	689	698
   gen=2	519	692	713
   gen=3	479	709	704
   gen=4	539	710	688
   gen=5	496	482	698
   gen=6	489	488	701
   gen=7	495	508	715
   gen=8	497	488	688
   gen=9	528	498	698
   (10, 10, 10)

   now exiting runScriptInteractively...

`Download matingSchemeByRepAndGen.py <matingSchemeByRepAndGen.py>`_

A function can be passed as the condition of a :class:`ConditionalMating` mating
scheme. This allows you to apply operators such as :class:`Stat` to examine the
condition of populations more closely and determine which mating scheme to use.


