An overview of simuPOP concepts
===============================

A simuPOP **population** consists of **individuals** of the same **genotype
structure**, which includes properties such as number of homologous sets of
chromosomes (ploidy), number of chromosomes, and names and locations of markers
on each chromosome. In addition to basic information such as genotypes and sex,
individuals can have arbitray auxillary values as **information fields**.
Individuals in a population can be divided into **subpopulations** that can be
further grouped into **virtual subpopulations** according to individual
properties such as sex, affection status, or arbitrary auxiliary information
such as age. Whereas subpopulations define boundaries of individuals that
restrict the flow of individuals and their genotypes (mating happens within
subpopulations), virtual subpopulations are groups of individuals who share the
same properties, with membership of individuals change easily with change of
individual properties.

**Figure**: *A life cycle of an evolutionary process*

.. _fig_life_cycle:

.. figure:: /Users/bpeng1/simuPOP/simuPOP/doc/figures/evolve.png
   :width: 680


Illustration of the discrete-generation evolutionary model used by simuPOP.

**Operators** are Python objects that act on a population. They can be applied
to a population before or after mating during a life cycle of an evolutionary
process (Figure :ref:`fig_life_cycle <fig_life_cycle>`), or to parents and
offspring during the production of each offspring. Arbitrary numbers of
operators can be applied to an evolving population.

A simuPOP **mating scheme** is responsible for choosing parent or parents from a
parental (virtual) subpopulation and for populating an offspring subpopulation.
simuPOP provides a number of pre-defined **homogeneous mating schemes**, such as
random, monogamous or polygamous mating, selfing, and haplodiploid mating in
hymenoptera. More complicated nonrandom mating schemes such as mating in age-
structured populations can be constructed using **heterogeneous mating
schemes**, which applies multiple homogeneous mating schemes to different
(virtual) subpopulations.

simuPOP evolves a population generation by generation, following the
evolutionary cycle depicted in Figure :ref:`fig_life_cycle <fig_life_cycle>`.
Briefly speaking, a number of **operators** such as a :class:`KAlleleMutator`
are applied to a population before a mating scheme repeatedly chooses a parent
or parents to produce offspring. **During-mating operators** such as
:class:`Recombinator` can be applied by a mating scheme to transmit parental
genotype to offspring. After an offspring population is populated, other
**operators** can be applied, for example, to calculate and output population
statistics. The offspring population will then become the parental population of
the next evolutionary cycle. Many simuPOP operators can be applied in different
stages so the type of an operator is determined by the stage at which it is
applied. Several populations, or replicates of a single population, could form a
**simulator** and evolve together.

.. _simple_example:

**Example**: *A simple example*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=1000, loci=2)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[1, 2, 2, 1])
   ...     ],
   ...     matingScheme=sim.RandomMating(ops=sim.Recombinator(rates=0.01)),
   ...     postOps=[
   ...         sim.Stat(LD=[0, 1], step=10),
   ...         sim.PyEval(r"'%.2f\n' % LD[0][1]", step=10),
   ...     ],
   ...     gen=100
   ... )
   0.25
   0.23
   0.20
   0.20
   0.18
   0.15
   0.12
   0.10
   0.10
   0.11
   100

   now exiting runScriptInteractively...

`Download simpleExample.py <simpleExample.py>`_

Some of these concepts are demonstrated in Example :ref:`simple_example
<simple_example>`, where a standard diploid Wright-Fisher model with
recombination is simulated. The first line imports the standard simuPOP module.
The second line creates a diploid population with 1000 individuals, each having
one chromosome with two loci. The ``evolve()`` function evolves the population
using a random mating scheme and four operators.

Operators :class:`InitSex` and :class:`InitGenotype` are applied at the
beginning of the evolutionary process. Operator :class:`InitSex` initializes
individual sex randomly and :class:`InitGenotype` initializes all individuals
with the same genotype ``12/21``. The populations are then evolved for 100
generations. A random mating scheme is used to generate offspring. Instead of
using the default Mendelian genotype transmitter, a :class:`Recombinator`
(during-mating operator) is used to recombine parental chromosomes with the
given recombination rate ``0.01`` during the generation of offspring. The other
operators are applied to the offspring generation (post-mating) at every 10
generations (parameter ``step``). Operator :class:`Stat` calculates linkage
disequilibrium between the first and second loci. The results of this operator
are stored in a local variable space of the Population. The last operator
:class:`PyEval` outputs calculated linkage disequilibrium values with a trailing
new line. The result represents the decay of linkage disequilibrium of this
population at 10 generation intervals. The return value of the ``evolve``
function, which is the number of evolved generations, is also printed.


