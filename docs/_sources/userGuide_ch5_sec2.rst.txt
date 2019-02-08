Initialization
==============

simuPOP provides three operators to initialize individual sex, information
fields and genotype at the population level. A number of parameter are provided
to cover most commonly used initialization scenarios. A Python operator can be
used to intialize a population explicitly if none of the operators fits your
need.


Initialize individual sex (operator :class:`InitSex`)
-----------------------------------------------------

Operator :class:`InitSex`\ () and function ``initSex()`` initialize individual
sex either randomly or using a given sequence. In the first case, individuals
are assigned ``MALE`` or ``FEMALE`` with equal probability unless parameter
*maleFreq* is used to specify the probability of having a male Individual.
Alternatively, parameter *maleProp* can be used to specify exact proportions of
male individuals so that you will have exactly 1000 males and 1000 females if
you apply :class:`InitSex`\ (``maleProp=0.5``) to a population of 2000
individuals.

Both parameters ``maleFreq`` and ``maleProp`` assigns individual sex randomly.
If for some reason you need to specify individual sex explicitly, you could use
a sequence of sex (``MALE`` or ``FEMALE``) to assign sex to individuals
succesively. The list will be reused if needed. If a list of (virtual)
subpopulations are given, this operator will only initialize individuals in
these (virtual) subpopulations. Example :ref:`InitSex <InitSex>` demonstrates
how to use two :class:`InitSex` operators to initialize two subpopulations.

.. _InitSex:

**Example**: *Initialize individual sex*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[1000, 1000])
   >>> sim.initSex(pop, maleFreq=0.3, subPops=0)
   >>> sim.initSex(pop, sex=[sim.MALE, sim.FEMALE, sim.FEMALE], subPops=1)
   >>> sim.stat(pop, numOfMales=True, vars='numOfMales_sp')
   >>> print(pop.dvars(0).numOfMales)
   290
   >>> print(pop.dvars(1).numOfMales)
   334

   now exiting runScriptInteractively...

`Download InitSex.py <InitSex.py>`_


Initialize genotype (operator :class:`InitGenotype`)
----------------------------------------------------

Operator :class:`InitGenotype` (and its function form ``initGenotype``)
initializes individual genotype by allele frequency, allele proportion,
haplotype frequency, haplotype proportions or a list of genotypes:

* By frequency of alleles. For example, :class:`InitGenotype`\ (``freq=(0, 0.2,
  0.4, 0.2)``) will assign allele 0, 1, 2, and 3 with probability 0, 0.2, 0.4 and
  0.2 respectively.

* By proportions of alleles. For example, :class:`InitGenotype`\ (``prop=(0,
  0.2, 0.4, 0.2)``) will assign 400 allele 1, 800 allele 2 and 400 allele 3 to a
  diploid population with 800 individuals.

* By frequency of haplotypes. For example, :class:`InitGenotype`\
  (``haplotypes=[[0, 0], [1,1], [0,1],[1,1]]``) will assign four haplotypes with
  equal probabilities. :class:`InitGenotype`\ (``haplotypes=[[0, 0], [1,1],
  [0,1],[1,1]], freq=[0.2, 0.2, 0.3, 0.3]``) will assign these haplotypes with
  different frequencies. If there are more than two loci, the haplotypes will be
  repeated.

* By frequency of haplotypes. For example, :class:`InitGenotype`\
  (``haplotypes=[[0, 0], [1,1], [0,1],[1,1]], prop=[0.2, 0.2, 0.3, 0.3]``) will
  assign four haplotypes with exact proportions.

* By a list of genotype. For example, :class:`InitGenotype`\ (``genotype=[1, 2,
  2, 1]``) will assign genotype ``1``, ``2``, ``2``, ``1`` repeatedly to a
  population. If individuals in this population has two homologous copies of a
  chromosome with two loci, this operator will assign haplotype ``1``, ``2`` to
  the first homologous copy of the chromosome, and ``2``, ``1`` to the second
  copy.

* By multiple allele frequencies or proportions returned by a function passed to
  parameter ``freq`` or ``prop`` (new in version 1.1.7). This function can accept
  parameters ``loc``, ``subPop`` or ``vsp`` and returns locus, subpopopulation or
  virtual subpopulation specific allele frequencies. For example, if you would
  like to initialize genotypes with random allele frequency, you can set
  ``freq=lambda : random.random()`` so that a new frequency is drawn from an
  uniform distribution for each new locus. Note that simuPOP expects the return
  value of this function to be a list of frequencies for alleles 0, 1, ..., but
  treats a single return value *x* as [*x, 1-x*] for simplicity.

Parameter ``loci`` and ``ploidy`` can be used to specify a subset of loci and
homologous sets of chromosomes to initialize, and parameter ``subPops`` can be
used to specify subsets of individuals to initialize. Example :ref:`InitGenotype
<InitGenotype>` demonstrates how to use these the :class:`InitGenotype`
operator, including examples on how to define and use virtual subpopulations to
initialize individual genotype by sex or by proportion.

.. _InitGenotype:

**Example**: *Initialize individual genotype*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[2000, 3000], loci=[5, 7])
   >>> # by allele frequency
   >>> def printFreq(pop, loci):
   ...     sim.stat(pop, alleleFreq=loci)
   ...     print(', '.join(['{:.3f}'.format(pop.dvars().alleleFreq[x][0]) for x in loci]))
   ... 
   >>> sim.initGenotype(pop, freq=[.4, .6])
   >>> sim.dump(pop, max=6, structure=False)
   SubPopulation 0 (), 2000 Individuals:
      0: MU 11000 0011111 | 11111 0101110 
      1: MU 00000 1111111 | 11101 1111001 
      2: MU 10111 0111100 | 01111 1011111 
      3: MU 11011 1101010 | 11010 1011111 
      4: MU 11011 0011010 | 10011 1001110 
      5: MU 00001 1010011 | 11111 1111110 
   SubPopulation 1 (), 3000 Individuals:
   2000: MU 10011 0010100 | 01001 0011010 

   >>> printFreq(pop, range(5))
   0.397, 0.404, 0.400, 0.402, 0.406
   >>> # by proportion
   >>> sim.initGenotype(pop, prop=[0.4, 0.6])
   >>> printFreq(pop, range(5))
   0.400, 0.400, 0.400, 0.400, 0.400
   >>> # by haplotype frequency
   >>> sim.initGenotype(pop, freq=[.4, .6], haplotypes=[[1, 1, 0, 1], [0, 0, 1]])
   >>> sim.dump(pop, max=6, structure=False)
   SubPopulation 0 (), 2000 Individuals:
      0: MU 11011 1011101 | 00100 1001001 
      1: MU 11011 1011101 | 11011 1011101 
      2: MU 00100 1001001 | 00100 1001001 
      3: MU 00100 1001001 | 00100 1001001 
      4: MU 11011 1011101 | 11011 1011101 
      5: MU 00100 1001001 | 11011 1011101 
   SubPopulation 1 (), 3000 Individuals:
   2000: MU 00100 1001001 | 00100 1001001 

   >>> printFreq(pop, range(5))
   0.597, 0.597, 0.403, 0.597, 0.597
   >>> # by haplotype proportion
   >>> sim.initGenotype(pop, prop=[0.4, 0.6], haplotypes=[[1, 1, 0], [0, 0, 1, 1]])
   >>> printFreq(pop, range(5))
   0.600, 0.600, 0.400, 0.000, 0.600
   >>> # by genotype
   >>> pop = sim.Population(size=[2, 3], loci=[5, 7])
   >>> sim.initGenotype(pop, genotype=[1]*5 + [2]*7 + [3]*5 +[4]*7)
   >>> sim.dump(pop, structure=False)
   SubPopulation 0 (), 2 Individuals:
      0: MU 11111 2222222 | 33333 4444444 
      1: MU 11111 2222222 | 33333 4444444 
   SubPopulation 1 (), 3 Individuals:
      2: MU 11111 2222222 | 33333 4444444 
      3: MU 11111 2222222 | 33333 4444444 
      4: MU 11111 2222222 | 33333 4444444 

   >>> # 
   >>> # use virtual subpopulation
   >>> pop = sim.Population(size=[2000, 3000], loci=[5, 7])
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> sim.initSex(pop)
   >>> sim.initGenotype(pop, genotype=range(10), loci=range(5))
   >>> # initialize all males
   >>> sim.initGenotype(pop, genotype=[2]*7, loci=range(5, 12),
   ...     subPops=[(0, 0), (1, 0)])
   >>> # assign genotype by proportions
   >>> pop.setVirtualSplitter(sim.ProportionSplitter([0.4, 0.6]))
   >>> sim.initGenotype(pop, freq=[0.2, 0.8], subPops=[(0,0)])
   >>> sim.initGenotype(pop, freq=[0.5, 0.5], subPops=[(0,1)])
   >>> #
   >>> # initialize by random allele frequency
   >>> import random
   >>> sim.initGenotype(pop, freq=lambda : random.random())
   >>> printFreq(pop, range(5))
   0.580, 0.239, 0.100, 0.576, 0.674
   >>> # initialize with loci specific frequency. here
   >>> # lambda loc: 0.01*loc is equivalent to 
   >>> # lambda loc: [0.01*loc, 1-0.01*loc]
   >>> sim.initGenotype(pop,
   ...     freq=lambda loc: 0.01*loc)
   >>> printFreq(pop, range(5))
   0.000, 0.009, 0.018, 0.029, 0.041
   >>> # initialize with VSP-specific frequency
   >>> sim.initGenotype(pop,
   ...     freq=lambda vsp: [[0.2, 0.8], [0.5, 0.5]][vsp[1]],
   ...     subPops=[(0, 0), (0, 1)])
   >>> 

   now exiting runScriptInteractively...

`Download InitGenotype.py <InitGenotype.py>`_


Initialize information fields (operator :class:`InitInfo`)
----------------------------------------------------------

Operator :class:`InitInfo` and its function form ``initInfo`` initialize one or
more information fields of all individuals or Individuals in selected (virtual)
subpopulations using either a list of values or a Python function. If a value or
a list of value is given, it will be used repeatedly to assign values of
specified information fields of all applicable individuals. For example,
``initInfo(pop, values=1, infoFields='x')`` will assign value ``1`` to
information field ``x`` of all individuals, and   ::

   initInfo(pop, values=[1, 2, 3], infoFields='x', subPops=[(0,1)])

will assign values ``1``, ``2``, ``3``, ``1``, ``2``, ``3``... to information
field ``x`` of individuals in the second virtual subpopulation of subpopulation
0.

The ``values`` parameter also accepts a Python function. This feature is usually
used to assign random values to an information field. For example,
``values=random.random`` would assign a random value between 0 and 1. If a
function takes parameters, a lambda function can be used. For example,   ::

   initInfo(pop, lambda : random.randint(2, 5), infoFields=['x', 'y'])

assigns random integers between 2 and 5 to information fields ``x`` and ``y`` of
all individuals in *pop*. Example :ref:`InitInfo <InitInfo>` demonstrates these
usages.

.. _InitInfo:

**Example**: *initialize information fields*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(size=[5], loci=[2], infoFields=['sex', 'age'])
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> sim.initSex(pop)
   >>> sim.initInfo(pop, 0, subPops=[(0,0)], infoFields='sex')
   >>> sim.initInfo(pop, 1, subPops=[(0,1)], infoFields='sex')
   >>> sim.initInfo(pop, lambda: random.randint(20, 70), infoFields='age')
   >>> sim.dump(pop, structure=False)
   SubPopulation 0 (), 5 Individuals:
      0: FU 00 | 00 |  1 39
      1: FU 00 | 00 |  1 29
      2: MU 00 | 00 |  0 68
      3: MU 00 | 00 |  0 50
      4: MU 00 | 00 |  0 21


   now exiting runScriptInteractively...

`Download InitInfo.py <InitInfo.py>`_


