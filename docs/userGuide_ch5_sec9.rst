Natural Selection
=================


Natural selection through the selection of parents
--------------------------------------------------

In the simplest scenario, natural selection is implemented in two steps:

* Before mating happens, an operator (called a **selector**) goes through a
  population and assign each individual a fitness value. The fitness values are
  stored in an information field called ``fitness``.

* When mating happens, parents are chosen with probabilities that are
  proportional to their fitness values. For example, assuming that a parental
  population consists of four Individuals with fitness values 1, 2, 3, and 4,
  respectively, the probability that they are picked to produce offspring are
  :math:`1/\left(1+2+3+4\right)=0.1`, :math:`0.2`, :math:`0.3`, and :math:`0.4`
  respectively. As you can image, if the offspring population has 10 individuals,
  the four parents will on average parent 1, 2, 3 and 4 offspring.

Because parents with lower fitness values have less chance to be produce
offspring, their genotypes have less chance to be passed to an offspring
generation. If the decreased fitness is caused by the presence of certain mutant
(e.g. a mutant causing a serious disease), individuals with that mutant will
have less change to survive and effecitively reduce or eleminate that mutant
from the population.

Example :ref:`selectParents <selectParents>` gives an example of natural
selection. In this example, a :class:`MapSelector` is used to explicitly assign
fitness value to genotypes at the first locus. The fitness values are ``1``,
``0.98``, ``0.97`` for genotypes ``00``, ``01`` and ``11`` respectively. The
selector set individual fitness values to information field ``fitness`` before
mating happens. The :class:`RandomMating` mating scheme then selects parents
according to parental fitness values.

.. _selectParents:

**Example**: *Natural selection through the selection of parents*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(4000, loci=1, infoFields='fitness')
   >>> simu = sim.Simulator(pop, rep=3)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     preOps=sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):0.98, (1,1):0.97}),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0, step=10),
   ...         sim.PyEval("'Gen:%3d ' % gen", reps=0, step=10),
   ...         sim.PyEval(r"'%.3f\t' % alleleFreq[0][1]", step=10),
   ...         sim.PyOutput('\n', reps=-1, step=10)
   ...     ],
   ...     gen = 50
   ... )
   Gen:  0 0.490	0.492	0.487	
   Gen: 10 0.433	0.430	0.431	
   Gen: 20 0.403	0.390	0.419	
   Gen: 30 0.343	0.325	0.383	
   Gen: 40 0.303	0.297	0.334	
   (50, 50, 50)

   now exiting runScriptInteractively...

`Download selectParents.py <selectParents.py>`_

.. note::

   The selection algorithm used in simuPOP is called *fitness proportionate
   selection*, or *roulette-wheel selection*. simuPOP does not use the more
   efficient *stochastic universal sampling* algorithm because the number of needed
   offspring is unknown in advance.


Natural selection through the selection of offspring \*
-------------------------------------------------------

Natural selection can also be implemented as selection of offspring. Remember
that an individual will be discarded if one of the during-mating operators fails
(return ``False``), **a** **during-mating selector** **discards offspring
according to fitness values of offspring**. Instead of relative fitness that
will be compared against other individuals during the selection of parents,
**fitness values of a during-mating selector are considered as absolute fitness
which are probabilities to survive** and have to be between 0 and 1.

A during-mating selector works as follows:

#. During evolution, parents are chosen randomly to produce one or more
   offspring. (Nothing prevents you from choosing parents according to their
   fitness values, but it is rarely justifiable to apply natural selection to both
   parents and offspring.)

#. A selection operator is applied to each offspring during mating and
   determines his or her fitness value. The fitness value is considered as
   probability to survive so an offspring will be discarded (operator returns
   ``False``) if the fitnessvalue is larger than an uniform random number.

#. Repeat steps 1 and 2 until the offspring generation is populated.

Because many offspring will be generated and discarded, especially when
offspring fitness values are low, selection through offspring is less efficient
than selection through parents. In addition, absolute fitness is usually more
difficult to estimate than relative fitness. So, unless there are compelling
reasons (e.g. simulating realistic scenarios of survival competition among
offspring), selection through parents are recommended.

Example :ref:`selectOffspring <selectOffspring>` gives an example of natural
selection through the selection of offspring. This example looks almost
identical to Example :ref:`selectParents <selectParents>` but the underlying
selection mechanism is quite different. Note that selection through offspring
does not save fitness values to an information field so you do not need to add
information field fitness to the population.

.. _selectOffspring:

**Example**: *Natural selection through the selection of offspring*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(10000, loci=1)
   >>> simu = sim.Simulator(pop, rep=3)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     matingScheme=sim.RandomMating(ops=[
   ...         sim.MendelianGenoTransmitter(),
   ...         sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):0.98, (1,1):0.97}),
   ...     ]),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0, step=10),
   ...         sim.PyEval("'Gen:%3d ' % gen", reps=0, step=10),
   ...         sim.PyEval(r"'%.3f\t' % alleleFreq[0][1]", step=10),
   ...         sim.PyOutput('\n', reps=-1, step=10)
   ...     ],
   ...     gen = 50
   ... )
   Gen:  0 0.493	0.493	0.496	
   Gen: 10 0.461	0.464	0.465	
   Gen: 20 0.436	0.445	0.442	
   Gen: 30 0.389	0.386	0.385	
   Gen: 40 0.370	0.345	0.348	
   (50, 50, 50)

   now exiting runScriptInteractively...

`Download selectOffspring.py <selectOffspring.py>`_


Are two selection scenarios equivalent? \*\*
--------------------------------------------

If you look closely at Examples :ref:`selectParents <selectParents>` and
:ref:`selectOffspring <selectOffspring>`, you will notice that their results are
quite similar. This is actually what you should expect in most cases. Let us
look at the theoretical consequence of selection through parents or offspring in
a simple case with asexual mating.

Assuming a diallelic marker with three genotypes :math:`g_{AA}`, :math:`g_{Aa}`
and :math:`g_{aa}`, with frequencies :math:`P_{AA}`, :math:`P_{Aa}` and
:math:`P_{aa}`, and relative fitness values :math:`w_{AA}`, :math:`w_{Aa}`, and
:math:`w_{22}` respectively. If we select through offspring, the proportion of
genotype :math:`g_{AA}` etc., should be

.. math::

      P_{AA}'=\frac{P_{AA}w_{AA}}{P_{AA}w_{AA}+P_{Aa}w_{Aa}+P_{aa}w_{aa}}

.. math::

      P_{Aa}'=\frac{P_{Aa}w_{Aa}}{P_{AA}w_{AA}+P_{Aa}w_{Aa}+P_{aa}w_{aa}}

.. math::

      P_{aa}'=\frac{P_{aa}w_{aa}}{P_{AA}w_{AA}+P_{Aa}w_{Aa}+P_{aa}w_{aa}}

because offspring genotypes are randomly drawn from the parental generation, and
each offspring has certain probability to survive.

Now, if we select through parents, the proportion of parents with genotype
:math:`AA` will be the number of :math:`AA` individuals times its probability to
be chosen:

.. math::

      n_{AA}\frac{w_{AA}}{\sum_{n=1}^{N}w_{n}}

This is, however, exactly

.. math::

      n_{AA}\frac{w_{AA}}{\sum_{n=1}^{N}w_{n}}=\frac{n_{AA}w_{AA}}{n_{AA}w_{AA}+n_{Aa}w_{Aa}+n_{aa}w_{aa}}=\frac{P_{AA}w_{AA}}{P_{AA}w_{AA}+P_{Aa}w_{Aa}+P_{aa}w_{aa}}=P_{AA}'

which corresponds to the proportion of offspring with such genotype. That is to
say, **in this simple case, two types of selection scenarios yield identical
results**.

These two types of selection scenarios do not have to always yield identical
results. Exceptions exist in cases with more than one offspring or sexual mating
with sex-specific survival rate. simuPOP provides both selection implementations
and you should choose one of them for your particular simulation.


Map selector (operator :class:`MapSelector`)
--------------------------------------------

A map selector uses a Python dictionary to provide fitness values for each type
of genotype. For example, Example :ref:`MapSelector <MapSelector>` uses a
dictionary with keys ``(0,0)``, ``(0,1)`` and ``(1,1)`` to specify fitness
values for individuals with these genotypes at locus 0. This example is a
typical example of heterozygote advantage. When :math:`w_{11}<w_{12}>w_{22},`
the genotype frequencies will go to an equilibrium state. Theoretically, if
:math:`s_{1}=w_{12}-w_{11}` and :math:`s_{2}=w_{12}-w_{22}`, the stable allele
frequency of allele 0 is

.. math::

      p=\frac{s_{2}}{s_{1}+s_{2}}

which is :math:`\frac{2}{3}` in the example (:math:`s_{1}=.1`,
:math:`s_{2}=.2`).

.. _MapSelector:

**Example**: *A selector that uses pre-defined fitness value*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=1000, loci=1, infoFields='fitness')
   >>> s1 = .1
   >>> s2 = .2
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.2, .8])
   ...     ],
   ...     preOps=sim.MapSelector(loci=0, fitness={(0,0):1-s1, (0,1):1, (1,1):1-s2}),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.PyEval(r"'%.4f\n' % alleleFreq[0][0]", step=100)
   ...     ],
   ...     gen=301
   ... )
   0.2250
   0.6605
   0.6530
   0.6870
   301
   >>> 

   now exiting runScriptInteractively...

`Download MapSelector.py <MapSelector.py>`_

The above example assumes that the fitness value for individuals with genotypes
``(0,1)`` and ``(1,0)`` are the same. This assumption is usually valid but can
be vialoated with impriting. In that case, you can specify fitness for both
types of genotypes. The underlying mechanism is that the :class:`MapSelector`
looks up a genotype in the dictionary first directly, and then without phase
information if a genotype is not found.

This operator supports haplodiploid populations and sex chromosomes. In these
cases, only valid alleles should be listed which can lead to dictionary keys
with different lengths. In addition, although less used because of potentially a
large number of keys, this operator can act on multiple loci. Please refer to
:class:`MapPenetrance` for details.


Multi-allele selector (operator :class:`MaSelector`)
----------------------------------------------------

A multi-allele selector divides alleles into two groups, wildtype *A* and
mutants *a*, and treat alleles within each group as the same. The fitness model
is therefore simplified to

* Two fitness values for genotype :math:`A`, :math:`a` in the haploid case

* Three fitness values for genotype *AA*, *Aa* and *aa* in the diploid single
  locus case. Genotype Aa and aA are assumed to have the same impact on fitness.

The default wildtype group contains allele 0 so the two allele groups are zero
and non-zero alleles. Example :ref:`MaSelector <MaSelector>` demonstrates the
use of this operator. This example is identical to Example :ref:`MapSelector
<MapSelector>` except that there are five alleles at locus 0 and alleles 1, 2,
3, 4 are treated as a single non-widetype group.

.. _MaSelector:

**Example**: *A multi-allele selector*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=1000, loci=1, infoFields='fitness')
   >>> s1 = .1
   >>> s2 = .2
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.2] * 5)
   ...     ],
   ...     preOps=sim.MaSelector(loci=0, fitness=[1-s1, 1, 1-s2]),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.PyEval(r"'%.4f\n' % alleleFreq[0][0]", step=100)
   ...     ],
   ...     gen = 301)
   0.2250
   0.6605
   0.6530
   0.6870
   301

   now exiting runScriptInteractively...

`Download MaSelector.py <MaSelector.py>`_

Operator :class:`MaSelector` also supports multiple loci by specifying fitness
values for all combination of genotype at specified loci. In the case of two
loci, this operator requires

* Four fitness values for genotype ``AB``, ``Ab``, ``aB`` and ``ab`` in the
  haploid case,

* Nine fitness values for genotype ``AABB``, ``AABb``, ``AAbb``, ``AaBB``,
  ``AaBb``, ``Aabb``, ``aaBB``, ``aaBb``, and ``aabb`` in the haploid case.

In general, :math:`2^{n}` values are needed for haploid populations and
:math:`3^{n}` values are needed for diploid populations where :math:`n` is the
number of loci. This operator does not yet support haplodiploid populations and
sex chromosomes. Example :ref:`MaSelectorHaploid <MaSelectorHaploid>`
demonstrates the use of a multi-locus model in a haploid population.

.. _MaSelectorHaploid:

**Example**: *A multi-locus multi-allele selection model in a haploid population*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=10000, ploidy=1, loci=[1,1], infoFields='fitness')
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.5, .5])
   ...     ],
   ...     # fitness values for AB, Ab, aB and ab
   ...     preOps=sim.MaSelector(loci=[0,1], fitness=[1, 1, 1, 0.95]),
   ...     matingScheme=sim.RandomSelection(),
   ...     postOps=[
   ...         sim.Stat(haploFreq=[0, 1], step=25),
   ...         sim.PyEval(r"'%.3f\t%.3f\t%.3f\t%.3f\n' % (haploFreq[(0,1)][(0,0)],"
   ...                 "haploFreq[(0,1)][(0,1)], haploFreq[(0,1)][(1,0)],"
   ...                 "haploFreq[(0,1)][(1,1)])", step=25)
   ...     ],
   ...     gen = 100
   ... )
   0.264	0.243	0.252	0.240
   0.292	0.294	0.321	0.093
   0.339	0.330	0.303	0.027
   0.310	0.383	0.297	0.009
   100

   now exiting runScriptInteractively...

`Download MaSelectorHaploid.py <MaSelectorHaploid.py>`_


Multi-locus selection models (operator :class:`MlSelector`)
-----------------------------------------------------------

Although an individual's fitness can be affected by several factors, each of
which can be modeled individually, **only one fitness value is used to determine
a person's ability to pass all these factors to his or her offspring**. Although
in theory we sometimes assume independent evolution of disease predisposing loci
(mostly for mathematical reasons), in practise we have to use a multi-locus
selection model to combine single-locus models.

This multi-loci selector applies several selectors to each individual and
computes an overall fitness value from the fitness values provided by these
selectors. Although this selector is designed to obtain multi-loci fitness
values from several single-locus fitness models, any selector, including those
obtain their fitness values from multiple disease predisposing loci, can be used
in this selector. This selector uses parameter ``mode`` to control how
individual fitness values are combined. More specifically, if :math:`f_{i}` are
fitness values obtained from individual selectors, this selector returns

* :math:`\Pi_{i}f_{i}` if ``mode=MULTIPLICATIVE``, and

* :math:`1-\sum_{i}\left(1-f_{i}\right)` if ``mode=ADDITIVE``, and

* :math:`1-\Pi_{i}\left(1-f_{i}\right)` if ``mode=HETEROGENEITY``

0 will be returned if the returned fitness value is less than 0.

This operator simply combines individual fitness values and it is your
responsibility to apply and interpret these models. For example, if relative
fitness values are greater than one, the heterogeneity model hardly makes sense.
Example :ref:`MlSelector <MlSelector>` demonstrates the use of this operator
using an additive multi-locus model over an additive and a recessive single-
locus model at two diesease predisposing loci. For comparison, we simulate two
additional replicates with selection only applying to one of the two loci. It
would be interesting to see if these two loci evolve more or less independently
by comparing allele freqency trajectories of these two replicates to those in
the first replicate.

.. _MlSelector:

**Example**: *A multi-loci selector*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=10000, loci=2, infoFields='fitness')
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.5, .5])
   ...     ],
   ...     preOps=[
   ...         sim.MlSelector([
   ...             sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):1, (1,1):.8}),
   ...             sim.MapSelector(loci=1, fitness={(0,0):1, (0,1):0.9, (1,1):.8}),
   ...             ], mode = sim.ADDITIVE, reps=0),
   ...         sim.MapSelector(loci=0, fitness={(0,0):1, (0,1):1, (1,1):.8}, reps=1),
   ...         sim.MapSelector(loci=1, fitness={(0,0):1, (0,1):0.9, (1,1):.8}, reps=2)
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...          sim.Stat(alleleFreq=[0,1]),
   ...          sim.PyEval(r"'REP %d:\t%.3f\t%.3f\t' % (rep, alleleFreq[0][1], alleleFreq[1][1])"),
   ...          sim.PyOutput('\n', reps=-1),
   ...     ],
   ...     gen = 5
   ... )
   REP 0:	0.472	0.465	
   REP 0:	0.452	0.429	
   REP 0:	0.429	0.397	
   REP 0:	0.405	0.378	
   REP 0:	0.382	0.355	
   5

   now exiting runScriptInteractively...

`Download MlSelector.py <MlSelector.py>`_


A hybrid selector (operator :class:`PySelector`)
------------------------------------------------

When your selection model involves multiple interacting genetic and
environmental factors, it might be easier to calculate a fitness value
explicitly using a Python function. A hybrid selector can be used for this
purpose. If your selection model depends solely on genotype, you can define a
function such as

::

   def fitness_func(geno):
       # calculate fitness according to genotype at specified loci
       # genotypes are arrange locus by locus, namely A1,A2,B1,B2 for loci A and B
       return val

and use this function in an operator :class:`PySelector`\ (``func=fitness_func,
loci=loci``). If your selection model depends on genotype as well as some
information fields, you can define a function in the form of

::

   def fitness_func(geno, field1, field2):
       # calculate fitness according to genotype at specified loci
       # and values at specified informaton fields.
       return val

where ``field1``, ``field2`` are names of information fields. simuPOP will pass
genotype and value of specified fields according to name of the passed function.
Note that genotypes are arrange locus by locus, namely in the order of
``A1``,``A2``,``B1``,``B2`` for loci ``A`` and ``B``. Other parameters such as
``gen``, ``ind``, and ``pop`` are also allowed. Please check the reference
manual for details.

When a :class:`PySelector` is used to calculate fitness for an individual
(parents if applied pre-mating, offspring if applied during-mating), it will
collect his or her genotype at specified loci, optional values at specified
information fields, generation number, or individual to a user-specified Python
function, and take its return value as fitness. As you can imagine, the
incorporation of information fields and generation number allow the
implementation of very complex selection scenarios such as gene environment
interaction and varying selection pressures.

Example :ref:`PySelector <PySelector>` demonstrates how to use a
:class:`PySelector` to specify fitness values according to a fitness table and
the smoking status of each individual.

.. _PySelector:

**Example**: *A hybrid selector*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(size=2000, loci=[1]*2, infoFields=['fitness', 'smoking'])
   >>> s1 = .02
   >>> s2 = .03
   >>> # the second parameter gen can be used for varying selection pressure
   >>> def sel(geno, smoking):
   ...     #     BB  Bb   bb
   ...     # AA  1   1    1
   ...     # Aa  1   1-s1 1-s2
   ...     # aa  1   1    1-s2
   ...     #
   ...     # geno is (A1 A2 B1 B2)
   ...     if geno[0] + geno[1] == 1 and geno[2] + geno[3] == 1:
   ...         v = 1 - s1  # case of AaBb
   ...     elif geno[2] + geno[3] == 2:
   ...         v = 1 - s2  # case of ??bb
   ...     else:                
   ...         v = 1       # other cases
   ...     if smoking:
   ...         return v * 0.9
   ...     else:
   ...         return v
   ... 
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.5, .5])
   ...     ],
   ...     preOps=sim.PySelector(loci=[0, 1], func=sel),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         # set smoking status randomly
   ...         sim.InitInfo(lambda : random.randint(0,1), infoFields='smoking'),
   ...         sim.Stat(alleleFreq=[0, 1], step=20),
   ...         sim.PyEval(r"'%.4f\t%.4f\n' % (alleleFreq[0][1], alleleFreq[1][1])", step=20)
   ...     ],
   ...     gen = 50
   ... )
   0.4943	0.4890
   0.4880	0.4285
   0.4898	0.4073
   50

   now exiting runScriptInteractively...

`Download PySelector.py <PySelector.py>`_


Multi-locus random fitness effects (operator :class:`PyMlSelector`)
-------------------------------------------------------------------

If the fitness of individuals is determined by fitness effects over a large
number of loci, both :class:`MlSelector` and :class:`PySelector` are difficult
to use because the former requires a large number of single-locus selectors, and
the latter requires the processing long genome sequences. If the overall fitness
can be determined by fitness effects of mutants, a :class:`PyMlSelector` can be
used. This operator

* Calls a user-provided call-back function for each locus with at least a mutant
  (non-zero allele). The function can accept location and genotype so the fitness
  can be location and genotype dependent. The return value is cached so the
  function will be called only once for each locus-genotype pair.

* The fitness of each individual is determined by fitness values of loci with at
  least one mutant, using the same methods as operator :class:`MlSelector`. This
  implicitly assumes that loci without any mutant have fitness value 1 and will
  not contribute to the final fitness value.

Example :ref:`PySelector <PySelector>` demonstrates how to use a
:class:`PyMlSelector` to implement a fitness model where each mutant has a
random fitness drawn from a Gamma distribution. An additive model is used so a
homozygote will have a fitness penalty that doubles that of a heterozygote.
Because the fitness values of heterozygote and homozygote at each locus are
requested separately, a class is used to store locus-specific s values.

The fitness value of each locus-genotype pair is outputted to a file, and it
should be interesting to plot the distribution of allele frequency at each locus
against the fitness values, because mutants that suffer from stronger negative
natural selection are supposed to be rarer.

.. _PyMlSelector:

**Example**: *Random fitness effect*

::

   >>> import simuOpt
   >>> simuOpt.setOptions(quiet=True, alleleType='mutant')
   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(size=2000, loci=[10000], infoFields=['fitness'])
   >>> 
   >>> class GammaDistributedFitness:
   ...     def __init__(self, alpha, beta):
   ...         self.coefMap = {}
   ...         self.alpha = alpha
   ...         self.beta = beta
   ...      
   ...     def __call__(self, loc, alleles):
   ...         # because s is assigned for each locus, we need to make sure the
   ...         # same s is used for fitness of genotypes 01 (1-s) and 11 (1-2s)
   ...         # at each locus
   ...         if loc in self.coefMap:
   ...             s = self.coefMap[loc]
   ...         else:
   ...             s = random.gammavariate(self.alpha, self.beta)
   ...             self.coefMap[loc] = s
   ...         #
   ...         if 0 in alleles:
   ...             return 1. - s
   ...         else:
   ...             return 1. - 2.*s
   ... 
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     preOps=[
   ...         sim.AcgtMutator(rate=[0.00001], model='JC69'),
   ...         sim.PyMlSelector(GammaDistributedFitness(0.23, 0.185),
   ...             output='>>sel.txt'),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(numOfSegSites=sim.ALL_AVAIL, step=50),
   ...         sim.PyEval(r"'Gen: %2d #seg sites: %d\n' % (gen, numOfSegSites)",
   ...             step=50)
   ...     ],
   ...     gen = 201
   ... )
   Gen:  0 #seg sites: 180
   Gen: 50 #seg sites: 1310
   Gen: 100 #seg sites: 1479
   Gen: 150 #seg sites: 1511
   Gen: 200 #seg sites: 1579
   201
   >>> print(''.join(open('sel.txt').readlines()[:5]))
   5855	1	0	0.978125
   1085	2	0	0.340724
   2907	0	1	0.998146
   7773	0	1	0.927273
   1835	0	2	0.999976


   now exiting runScriptInteractively...

`Download PyMlSelector.py <PyMlSelector.py>`_


Alternative implementations of natural selection
------------------------------------------------

If you know how natural selection works in simuPOP, you do not have to use a
selector to perform natural selection. For example,

* If you choose to use fitness values of parents to perform probabilistic
  natural selection during mating, you just need to set individual fitness in some
  way before mating. (You do not even have to use information field ``fitness``
  because you can specify which information field to use in a mating scheme using
  parameter ``selectionField``). This can be done through a penetrance model (as
  shown in the following example) where affected individuals are selected against
  during mating, a quantitative trait model (where a trait is defined to control
  individual fitness), or by setting information field fitness manually through a
  Python operator.

* If you would like to perform deterministic selection on certain phenotype, you
  can explicitly remove individuals before or during mating. More explicitly, you
  can use an operator :class:`DiscardIf` to remove parents before mating or remove
  offspring during mating according to certain status (disease status or
  quantitative trait), provided that the trait status is defined before this
  operator is applied.

Example :ref:`peneSelector <peneSelector>` demonstrates a commonly used case
where parents who are affected with certain disease are excluded from producing
offspring. In this example, a penetrance model (operator :class:`MaPenetrance`)
is applied to the parental generation to determine who will be affected. An
:class:`InfoExec` operator is used to set individual fitness to 1 if he or she
is unaffected, and 0 if he or she is affected. Due to the way parents are
selected, affected parents will not be able to produce offspring as long as
there is any unaffected individual. Because individual affection status is
determined by his or her genotype, this genotype - affection status - fitness
relationship could be implemented using an equivalent :class:`MaSelector`. This
method could be extended to :class:`InfoExec`\ (``'fitness = 1 -
0.01*ind.affected()', exposeInd='ind'``) to select against, but not remove,
affected parents, and similarly :class:`InfoExec`\ (``'fitness = 1 - 0.01*(LDL >
250)'``) to select against individuals according to a quantitative trait. For
this particular example, a :class:`DiscardIf` operator could be used, although
it can be slower because of the explicit removal of parents.

.. _peneSelector:

**Example**: *Natural selection according to individual affection status*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=2000, loci=1, infoFields='fitness')
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.5, .5])
   ...     ],
   ...     preOps=[
   ...         sim.MaPenetrance(loci=0, penetrance=[0.01, 0.1, 0.2]),
   ...         sim.Stat(numOfAffected=True, step=25, vars='propOfAffected'),
   ...         sim.PyEval(r"'Percent of affected: %.3f\t' % propOfAffected", step=50),
   ...         sim.InfoExec('fitness = not ind.affected()', exposeInd='ind')
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.PyEval(r"'%.4f\n' % alleleFreq[0][1]", step=50)
   ...     ],
   ...     gen=151
   ... )
   Percent of affected: 0.110	0.4713
   Percent of affected: 0.009	0.0095
   Percent of affected: 0.013	0.0000
   Percent of affected: 0.008	0.0000
   151

   now exiting runScriptInteractively...

`Download peneSelector.py <peneSelector.py>`_


Frequency dependent or dynamic selection pressure \*
----------------------------------------------------

If individual fitness depends on individual information fields and/or population
variables, you will have to calculate individual fitness using expressions or
functions. In order to access individual information fields and population
variable and calculate individual fitness, you have the option to

* Use a :class:`PySelector` and pass genotype, values of information fields,
  references to individual and population to a user-provided function, which
  returns fitness value for each individual.

* Use of :class:`PyOperator` to obtain information of the population (e.g.
  variables) and all individuals. Determine individual fitness and set information
  field ``fitness`` of all individuals.

* Use an operator :class:`InfoExec` to calculate individual fitness using
  expressions. This method can be more efficient than others because simuPOP does
  not have to call a user-provided function.

Example :ref:`freqDependentSelection <freqDependentSelection>` demonstrates an
example where the fitness values of individuals are calculated from allele
frequencies calculated using a :class:`Stat` operator. Because the fitness
values of individuals are 1, :math:`1-(p-0.5)\*0.1`, :math:`1-(p-0.5)\*0.2` for
genotype 00, 01 and 11 where :math:`p` is the frequency of allele 1, this allele
will be under purifying selection if its frequency is over 0.5, and positive
selection if its frequency is less than 0.5. Consequently, the frequency of this
allele will oscillate around 0.5 during evolution, as shown in the result of
this example.

.. _freqDependentSelection:

**Example**: *Frequency dependent selection*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=2000, loci=1, infoFields='fitness')
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.5, .5])
   ...     ],
   ...     preOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.InfoExec('''fitness = {
   ...             0: 1,
   ...             1: 1 - (alleleFreq[0][1] - 0.5)*0.1, 
   ...             2: 1 - (alleleFreq[0][1] - 0.5)*0.2}[ind.allele(0,0)+ind.allele(0,1)]''',
   ...             exposeInd='ind'),
   ...         sim.Stat(meanOfInfo='fitness'),
   ...         sim.PyEval(r"'alleleFreq=%.3f, mean fitness=%.5f\n' % (alleleFreq[0][1], meanOfInfo['fitness'])",
   ...             step=25),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     gen=151
   ... )
   alleleFreq=0.495, mean fitness=1.00045
   alleleFreq=0.504, mean fitness=0.99955
   alleleFreq=0.484, mean fitness=1.00150
   alleleFreq=0.492, mean fitness=1.00076
   alleleFreq=0.499, mean fitness=1.00005
   alleleFreq=0.526, mean fitness=0.99726
   alleleFreq=0.514, mean fitness=0.99856
   151

   now exiting runScriptInteractively...

`Download freqDependentSelector.py <freqDependentSelector.py>`_


.. _subsec_vspSelection:

Support for virtual subpopulations \*
-------------------------------------

Support for virtual subpopulations allows you to use different selectors for
different (virtual) subpopulations. Because virtual subpopulations may overlap,
and they do not have to cover all individuals in a subpopulation, it is
important to remember that

* If virtual subpopulations overlap, the fitness value set by the last selector
  will be used.

* If an individual is not included in any of the virtual subpopulation, its
  fitness value will be zero which will prevent them from producing any offspring.

Example :ref:`vspSelector <vspSelector>` demonstrates how to apply selectors to
virtual subpopulations. This example has two subpopulations, each having two
virtual subpopulations defined by sex. Natural selection is applied to male
individuals in the first subpopulation, and female individuals in the second
subpopulation. However, because the sex of offspring is randomly determined, the
selection actually decreases the disease allele frequency for all inviduals.

.. _vspSelector:

**Example**: *Selector in virtual subpopulations*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[5000, 5000], loci=1, infoFields='fitness')
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.5, .5])
   ...     ],
   ...     preOps=[
   ...         sim.MaSelector(loci=0, fitness=[1, 1, 0.98], subPops=[(0,0), (1,1)]),
   ...         sim.MaSelector(loci=0, fitness=[1, 0.99, 0.98], subPops=[(0,1), (1,0)]),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=[0], subPops=[(sim.ALL_AVAIL, sim.ALL_AVAIL)],
   ...             vars='alleleFreq_sp', step=50),
   ...         sim.PyEval(r"'%.4f\t%.4f\t%.4f\t%.4f\n' % "
   ...             "tuple([subPop[x]['alleleFreq'][0][1] for x in ((0,0),(0,1),(1,0),(1,1))])",
   ...             step=50)
   ...     ],
   ...     gen=151
   ... )
   0.5022	0.5083	0.4970	0.5020
   0.4086	0.4054	0.3849	0.3817
   0.3275	0.3259	0.2435	0.2532
   0.2715	0.2662	0.1305	0.1338
   151

   now exiting runScriptInteractively...

`Download vspSelector.py <vspSelector.py>`_

Selecting through offspring can also be applied to virtual subpopulations. For
example, Example :ref:`vspDuringMatingSelector <vspDuringMatingSelector>` moves
the selectors to the ``ops`` parameter of :class:`RandomMating`. In this way,
male and female offspring will have different survival probabilities according
to their genotype.

.. _vspDuringMatingSelector:

**Example**: *Selection against offspring in virtual subpopulations*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[5000, 5000], loci=1, infoFields='fitness')
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.5, .5])
   ...     ],
   ...     matingScheme=sim.RandomMating(ops=[
   ...         sim.MendelianGenoTransmitter(),
   ...         sim.MaSelector(loci=0, fitness=[1, 1, 0.98], subPops=[(0,0), (1,1)]),
   ...         sim.MaSelector(loci=0, fitness=[1, 0.99, 0.98], subPops=[(0,1), (1,0)]),
   ...         ]),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=[0], subPops=[(sim.ALL_AVAIL, sim.ALL_AVAIL)],
   ...             vars='alleleFreq_sp', step=50),
   ...         sim.PyEval(r"'%.4f\t%.4f\t%.4f\t%.4f\n' % "
   ...             "tuple([subPop[x]['alleleFreq'][0][1] for x in ((0,0),(0,1),(1,0),(1,1))])",
   ...             step=50)
   ...     ],
   ...     gen=151
   ... )
   0.5018	0.5034	0.4941	0.4853
   0.3652	0.3728	0.3820	0.3766
   0.2882	0.2920	0.2590	0.2667
   0.2083	0.1994	0.2378	0.2356
   151

   now exiting runScriptInteractively...

`Download vspDuringMatingSelector.py <vspDuringMatingSelector.py>`_


Natural selection in heterogeneous mating schemes \*\*
------------------------------------------------------

Multiple mating schemes could be applied to the same subpopulation in a
heterogeneous mating scheme (:class:`HeteroMating`). These mating schemes may or
may not support natural selection, may be applied to different virtual
subpopulations of population, and they may see Individuals differently in terms
of individual fitness. Parameter ``fitnessField`` of a mating scheme could be
used to handle such cases. More specifically,

* You can turn off the natural selection support of a mating scheme by setting
  ``fitnessField=''``.

* If a mating scheme uses a different set of fitness values, you can add an
  information field (e.g. ``fitness1``), setting individual fitness to this
  information field using a selector (with parameter ``infoFields='fitness1'``)
  and tells a mating scheme to look in this information field for fitness values
  (using parameter ``fitnessField='fitness1'``).


