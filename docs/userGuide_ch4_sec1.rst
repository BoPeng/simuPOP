.. _sec_Genotypic_structure:

Genotypic structure
===================

.. index:: single: genotypic structure

.. index::
   single: GenoStruTrait; ploidy
   single: GenoStruTrait; ploidyName
   single: GenoStruTrait; numChrom
   single: GenoStruTrait; chromType
   single: GenoStruTrait; chromName
   single: GenoStruTrait; numLoci
   single: GenoStruTrait; locusPos
   single: GenoStruTrait; infoField
   single: GenoStruTrait; infoFields

Genotypic structure refers to structural information shared by all individuals
in a population, including number of homologous copies of chromosomes (c.f.
``ploidy(), ploidyName()``), chromosome types and names (c.f. ``numChrom(),
chromType(), chromName()``), position and name of each locus (c.f.
``numLoci(ch),`` ``locusPos(loc),`` ``locusName(loc)``), and axillary
information attached to each individual (c.f. ``infoField(idx), infoFields()``).
In addition to property access functions, a number of utility functions are
provided to, for example, look up the index of a locus by its name (c.f.
``locusByName()``, ``chromBegin()``, ``chromLocusPair()``).

In simuPOP, locus is a (named) position and alleles are just different numbers
at that position. **A locus can be a gene, a nucleotide, or even a deletion,
depending on how you define alleles and mutations**. For example,

* A codon can be simulated as a locus with 64 allelic states, or three locus
  each with 4 allelic states. Alleles in the first case would be codons such as
  ``AAC`` and a mutation event would mutate one codon to another (e.g. ``AAC`` ->
  ``ACC``). Alleles in the second case would be ``A``, ``C``, ``T`` or ``G``, and
  a mutation event would mutate one nucleotide to another (e.g. ``A`` -> ``G``).

* You can use 0 and 1 (and the binary module of simuPOP) to simulate SNP
  (single-nucleotide polymorphism) markers and ignore the exact meaning of 0 and
  1, or use 0, 1, 2, 3 to simulate different nucleotide (A, C, T, or G) in these
  markers. The mutation model in the second case would be more complex.

* For microsatellite markers, alleles are usually interpreted as the number of
  tandem repeats. It would be difficult (though doable) to simulate the expansion
  and contraction of genome caused by the mutation of microsatellite markers.

* The infinite site and infinite allele mutation models could be simulated using
  either a continuous sequence of nucleotides with a simple 2-allele mutation
  model, or a locus with a large number of possible allelic states. It is also
  possible to simulate an empty region (without any locus) with loci introduced by
  mutation events.

* If you consider deletion as a special allelic state, you can simulate gene
  deletions without shrinking a chromosome. For example, a deletion mutation event
  can set the allelic state of one or more loci to 0, which can no longer be
  mutated.

* Alleles in different individuals could be interpretted differently. For
  example, if you would like to simulate major chromosomal mutations such as
  inversion, you could use a super set of markers for different types of
  chromosomes and use an indicator (information field) to mark the type of
  chromosome and which markers are valid. Using virtual subpopulations, these
  individuals could be handled differently during mating.

* In an implementation of an infinite-sites model, **Individual loci are used to
  store mutation events**. In this example (Example :ref:`infiniteSites
  <infiniteSites>`), 100 loci are allocated for each individual and they are used
  to store mutation events (location of a mutation) that happens in a 10Mb region.
  Whenever a mutation event happens, its location is stored as an allele of an
  individual. At the end of the evolution, each individual has a list of mutation
  events which can be readily translated to real alleles. Similar ideas could be
  used to simulate the accumulation of recombination events.

In summary, the exact meaning of loci and their alleles are user defined. With
appropriate mutation model and mating scheme, it is even possible to simulate
phenotypic traits using this mechanism, although it is more natual to use
information fields for quatitative traits.

A genotypic structure can be retrieved from *Individual* and *Population*
objects. **Because a population consists of individuals of the same type,
genotypic information can only be changed for all individuals at the population
level**. populations in a simulator usually have the same genotypic structure
because they are created by as replicates, but their structure may change during
evolution. Example :ref:`genostructure <genostructure>` demonstrates how to
access genotypic structure functions at the population and individual levels.
Note that ``lociPos`` determines the order at which loci are arranged on a
chromosome. Loci positions and names will be rearranged if given ``lociPos`` is
unordered.

.. _genostructure:

**Example**: *Genotypic structure functions*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[2, 3], ploidy=2, loci=[5, 10],
   ...     lociPos=list(range(0, 5)) + list(range(0, 20, 2)), chromNames=['Chr1', 'Chr2'],
   ...     alleleNames=['A', 'C', 'T', 'G'])
   >>> # access genotypic information from the sim.Population
   >>> pop.ploidy()
   2
   >>> pop.ploidyName()
   'diploid'
   >>> pop.numChrom()
   2
   >>> pop.locusPos(2)
   2.0
   >>> pop.alleleName(1)
   'C'
   >>> # access from an individual
   >>> ind = pop.individual(2)
   >>> ind.numLoci(1)
   10
   >>> ind.chromName(0)
   'Chr1'
   >>> ind.locusName(1)
   ''
   >>> # utility functions
   >>> ind.chromBegin(1)
   5
   >>> ind.chromByName('Chr2')
   1
   >>> # loci pos can be unordered within each chromosome
   >>> pop = sim.Population(loci=[2, 3], lociPos=[3, 1, 1, 3, 2],
   ...     lociNames=['loc%d' % x for x in range(5)])
   >>> pop.lociPos()
   (1.0, 3.0, 1.0, 2.0, 3.0)
   >>> pop.lociNames()
   ('loc1', 'loc0', 'loc2', 'loc4', 'loc3')

   now exiting runScriptInteractively...

`Download genoStru.py <genoStru.py>`_

.. note::

   simuPOP does not assume any unit for loci positions. Depending on your
   application, it can be basepair (bp), kilo-basepair (kb), mega base pair (mb) or
   even using genetic-map distance such as centiMorgan. It is your responsibility
   to interpret and use loci positions properly. For example, recombination rate
   between two adjacent markers can be specified as the product between their
   physical distance and a recombination intensity. The scale of this intensity
   will vary by the unit assumed.

.. note::

   Names of loci, alleles and subpopulations are optional. Empty names will be used
   if they are not specified. Whereas ``locusName``, ``subPopName`` and
   ``alleleName`` always return a value (empty string or specified value) for any
   locus, subpopulation or allele, respectively, ``lociNames``, ``subPopNames`` and
   ``alleleNames`` only return specified values, which can be empty lists.


Haploid, diploid and haplodiploid populations
---------------------------------------------

simuPOP is most widely used to study human (diploid) populations. A large number
of mating schemes, operators and population statistics are designed around the
evolution of such a population. simuPOP also supports haploid and haplodiploid
populations although there are fewer choices of mating schemes and operators.
simuPOP can also support other types of populations such as triploid and
tetraploid populations, but these features are largely untested due to their
limited usage. It is expected that supports for these populations would be
enhanced over time with additional dedicated operators and functions.

For efficiency considerations, simuPOP saves the same numbers of homologous sets
of chromosomes even if some individuals have different numbers of homologous
sets in a population. For example, in a haplodiploid population, because male
individuals have only one set of chromosomes, their second homologous set of
chromosomes are *unused*, which are labeled as ``'_'``, as shown in Example
:ref:`haplodiploid <haplodiploid>`.

.. _haplodiploid:

**Example**: *An example of haplodiploid population*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[2,5], ploidy=sim.HAPLODIPLOID, loci=[3, 5])
   >>> sim.initGenotype(pop, freq=[0.3, 0.7])
   >>> sim.dump(pop)
   Ploidy: 2 (haplodiploid)
   Chromosomes:
   1:  (AUTOSOME, 3 loci)
      (1),  (2),  (3)
   2:  (AUTOSOME, 5 loci)
      (1),  (2),  (3),  (4),  (5)
   population size: 7 (2 subpopulations with 2, 5 Individuals)
   Number of ancestral populations: 0

   SubPopulation 0 (), 2 Individuals:
      0: MU 111 00001 | ___ _____ 
      1: MU 111 01110 | ___ _____ 
   SubPopulation 1 (), 5 Individuals:
      2: MU 111 11110 | ___ _____ 
      3: MU 101 11111 | ___ _____ 
      4: MU 110 11111 | ___ _____ 
      5: MU 101 11101 | ___ _____ 
      6: MU 110 11001 | ___ _____ 


   now exiting runScriptInteractively...

`Download haplodiploid.py <haplodiploid.py>`_


Autosomes, sex chromosomes, mitochondrial, and other types of chromosomes \*
----------------------------------------------------------------------------

The default chromosome type is autosome, which is the *normal* chromosomes in
diploid, and in haploid populations. simuPOP supports four other types of
chromosomes, namely *chromosome X*, *chromosome Y, mitochondrial,* and*
customized* chromosome types. Sex chromosomes are only valid in haploid
populations where chromosomes X and Y are used to determine the sex of an
offspring. Mitochondrial DNAs can exists in haploid or diploid populations, and
are inherited maternally. Customized chromosomes rely on user defined functions
and operators to be passed from parents to offspring.

Example :ref:`subPopName <subPopName>` shows how to specify different chromosome
types, and how genotypes of these special chromosomes are arranged.

.. _chromtypes:

**Example**: *Different chromosome types*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=6, ploidy=2, loci=[3, 3, 3, 2, 2, 4, 4],
   ...     chromTypes=[sim.AUTOSOME]*2 + [sim.CHROMOSOME_X, sim.CHROMOSOME_Y, sim.MITOCHONDRIAL]
   ...         + [sim.CUSTOMIZED]*2)
   >>> sim.initGenotype(pop, freq=[0.3, 0.7])
   >>> sim.dump(pop, structure=False) # does not display genotypic structure information
   SubPopulation 0 (), 6 Individuals:
      0: MU 111 000 011 __ 11 1111 1101 | 110 000 ___ 11 __ 1111 1011 
      1: MU 111 111 101 __ 11 1110 1011 | 111 011 ___ 11 __ 1110 1011 
      2: MU 110 101 011 __ 11 1011 0011 | 110 100 ___ 11 __ 1010 1111 
      3: MU 010 011 111 __ 11 1111 1111 | 110 010 ___ 11 __ 1111 0111 
      4: MU 101 000 111 __ 01 0111 0100 | 110 111 ___ 00 __ 0111 0001 
      5: MU 111 010 111 __ 10 0111 1011 | 111 111 ___ 11 __ 0111 1011 


   now exiting runScriptInteractively...

`Download chromType.py <chromType.py>`_

The evolution of sex chromosomes follow the following rules

* There can be only one X chromosome and one Y chromosome. It is not allowed to
  have only one kind of sex chromosome.

* The Y chromosome of female individuals are ignored. The second homologous copy
  of the X chromosome and the first copy of the Y chromosome are ignored for male
  individuals.

* During mating, female parent pass one of her X chromosome to her offspring,
  male parent pass chromosome X or Y to his offspring. Recombination is allowed
  for the X chromosomes of females, but not allowed for males.

* The sex of offspring is determined by the types of sex chromosomes he/she
  inherits, XX for female, and XY for male.

The evolution of mitochonrial DNAs follow the following rules

* There can be only one copy of mitochondrial DNA, exists for both males and
  females.

* In a non-haploid population where all chromosomes have multiple homologous
  copies, only the first copy is used for mitochondrial DNA.

* mtDNAs are inherited maternally

Customized chromosomes are used to model more complex types of chromosomes. They
rely on customized operators for inheritence. For example, if you would like to
model multiple copies of mitochondrial DNAs (cytohets with multiple organellar
chromosomes) in a cell, and the process of genetic drift of somatic cytoplasmic
segregation of mtDNAs, you can use multiple customized chromosomes to model
multiple cytohets (see section :ref:`subsec_Pre_defined_genotype_transmitters
<subsec_Pre_defined_genotype_transmitters>` for an Example). Figure
:ref:`fig_chromTypes <fig_chromTypes>` depicts the possible chromosome structure
of two diploid parents, and how offspring chromosomes are formed. It uses two
customized chromosomes to model multiple copies of mitochondrial chromosomes
that are passed randomly from mother to offspring. The second homologous copy of
customized chromosomes are unused in this example.

**Figure**: *Inheritance of different types of chromosomes in a diploid population*

.. _fig_chromTypes:

.. figure:: /Users/bpeng1/simuPOP/simuPOP/doc/figures/chromType.png
   :width: 680


individuals in this population have five chromosomes, one autosome (A), one X
chromosome (X), one Y chromosome (Y) and two customized chromosomes (C). The
customized chromosomes model multiple copies of mitochondrial chromosomes that
are passed randomly from mother to offspring. Y chromosomes for the female
parent, the second copy of chromosome X and the first copy of chromosome Y for
the male parent, and the second copy of customized chromosomes are unused (gray
chromosome regions). A male offspring inherits one copy of autosome from his
mother (with recombination), one copy of autosome from his father (with
recombination), an X chromosome from his mother (with recombination), a Y
chromosome from his father (without recombination), and two copies of the first
customized chromosome.


.. _subsec_stru_infoFields:

Information fields
------------------

Different kinds of simulations require different kinds of individuals.
individuals with only genotype information are sufficient to simulate the basic
Wright-Fisher model. Sex is needed to simulate such a model in diploid
populations with sex. individual fitness may be needed if selection is induced,
and age may be needed if the population is age-structured. In addition,
different types of quantitative traits or affection status may be needed to
study the impact of genotype on Individual phenotype. Because it is infeasible
to provide all such information to an individual, simuPOP keeps genotype, sex
(``MALE`` or ``FEMALE``) and affection status as *built-in properties* of an
individual, and all others as optional *information fields* (float numbers)
attached to each individual.

Information fields can be specified when a population is created, or added later
using population member functions. They are essential for proper operation of
many simuPOP operators. For example, all selection operators require information
field ``fitness`` to store evaluated fitness values for each individual.
Operator :class:`Migrator` uses information field ``migrate_to`` to store the ID
of subpopulation an individual will migrate to. An error will be raised if these
operators are applied to a population without needed information fields.

.. _basicInfoFields:

**Example**: *Basic usage of information fields*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(10, loci=[20], ancGen=1,
   ...     infoFields=['father_idx', 'mother_idx'])
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[0]*20+[1]*20)
   ...     ],
   ...     matingScheme=sim.RandomMating(
   ...         ops=[
   ...             sim.Recombinator(rates=0.01),
   ...             sim.ParentsTagger()
   ...         ]
   ...     ),
   ...     gen = 1
   ... )
   1
   >>> pop.indInfo('mother_idx')  # mother of all offspring
   (9.0, 8.0, 8.0, 0.0, 8.0, 9.0, 8.0, 7.0, 7.0, 9.0)
   >>> ind = pop.individual(0)
   >>> mom = pop.ancestor(ind.mother_idx, 1)
   >>> print(ind.genotype(0))
   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
   >>> print(mom.genotype(0))
   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
   >>> print(mom.genotype(1))
   [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

   now exiting runScriptInteractively...

`Download infoField.py <infoField.py>`_

Example :ref:`basicInfoFields <basicInfoFields>` demonstrates the basic usage of
information fields. In this example, a population with two information fields
``mother_idx`` and ``father_idx`` are created. Besides the present generation,
this population keeps one ancestral generations (``ancGen=1``, see Section
:ref:`subsec_Ancestral_populations <subsec_Ancestral_populations>` for details).
After initializing each individual with two chromosomes with all zero and all
one alleles respectively, the population evolves one generation, subject to
recombination at rate 0.01. Parents of each individual are recorded, by operator
:class:`ParentsTagger`, to information fields ``mother_idx`` and ``father_idx``
of each offspring\ ``.``

After evolution, the population is extracted from the simulator, and the values
of information field ``mother_idx`` of all individuals are printed. The next
several statements get the first Individual from the population, and his mother
from the parental generation using the indexes stored in this individual's
information fields. Genotypes at the first homologous copy of this individual's
chromosome is printed, along with two parental chromosomes.

**Information fields can only be added or removed at the population level**
because all individuals need to have the same set of fields. Values of
information fields could be accessed at Individual or population levels, using
functions such as :meth:`Individual.info`, :meth:`Individual.setInfo`,
``population.indInfo``, :meth:`Population.setIndInfo`. These functions will be
introduced in their respective classes.

.. note::

   Information fields can be located both by names and by indexes**,** the former
   provides better readability at a slight cost of performance because these names
   have to be translated into indexes each time. However, use of names are
   recommended in most cases for readability considerations.


