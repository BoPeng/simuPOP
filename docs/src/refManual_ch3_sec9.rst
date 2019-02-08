Natural selection
=================

.. index:: single: selection


class BaseSelector
------------------

.. class:: BaseSelector

   This class is the base class to all selectors, namely operators
   that perform natural selection. It defines a common interface for
   all selectors.

   A selector can be applied before mating or during mating. If a
   selector is applied to one or more (virtual) subpopulations of a
   parental population before mating, it sets individual fitness
   values to all involved parents to an information field (default to
   *fitness*). When a mating scheme that supports natural selection is
   applied to the parental population, it will select parents with
   probabilities that are proportional to individual fitness stored in
   an information field (default to *fitness*).  Individual fitness is
   considered **relative** fitness and can be any non-negative number.
   This simple process has some implications that can lead to advanced
   usages of natural selection in  simuPOP:

   + It is up to the mating scheme how to handle individual fitness.
     Some mating schemes do not support natural selection at all.

   + A mating scheme performs natural selection according to fitness
     values stored in an information field. It does not care how these
     values are set. For example, fitness values can be inherited from
     a parent using a tagging operator, or set directly using a Python
     operator.

   + A mating scheme can treat any information field as fitness field.
     If an specified information field does not exist, or if all
     individuals have the same fitness values (e.g. 0), the mating
     scheme selects parents randomly.

   + Multiple selectors can be applied to the same parental
     generation. individual fitness is determined by the last fitness
     value it is assigned.

   + A selection operator can be applied to virtual subpopulations and
     set fitness values only to part of the individuals.

   + individuals with zero fitness in a subpopulation with anyone
     having a positive fitness value will not be selected to produce
     offspring. This can sometimes lead to unexpected behaviors. For
     example, if you only assign fitness value to part of the
     individuals in a subpopulation, the rest of them will be
     effectively discarded. If you migrate individuals with valid
     fitness values to a subpopulation with all individuals having
     zero fitness, the migrants will be the only mating parents.

   + It is possible to assign multiple fitness values to different
     information fields so that different homogeneous mating schemes
     can react to different fitness schemes when they are used in a
     heterogeneous mating scheme.

   + You can apply a selector to the offspring generation using the
     *postOps* parameter of :meth:`Simulator.evolve`, these fitness
     values will be used when the offspring generation becomes
     parental generation in the next generation.

   Alternatively, a selector can be used as a during mating operator.
   In this case, it caculates fitness value for each offspring which
   will be treated as **absolute** fitness, namely the probability for
   each offspring to survive. This process uses the fact that an
   individual will be discarded when any of the during mating
   operators returns *False*. It is important to remember that:

   + individual fitness needs to be between 0 and 1 in this case.

   + Fitness values are not stored so the population does not need an
     information field *fitness*.

   + This method applies natural selection to offspring instead of
     parents. These two implementation can be identical or different
     depending on the mating scheme used.

   + Seleting offspring is less efficient than the selecting parents,
     especially when fitness values are low.

   + Parameter *subPops* are applied to the offspring population and
     is used to judge if an operator should be applied. It thus does
     not make sense to apply a selector to a virtual subpopulation
     with affected individuals.


   .. method:: BaseSelector(output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)


      Create a base selector object. This operator should not be
      created directly.



class MapSelector
-----------------

.. class:: MapSelector

   This selector assigns individual fitness values using a user-
   specified dictionary. This operator can be applied to populations
   with arbitrary number of homologous chromosomes.


   .. method:: MapSelector(loci, fitness, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)


      Create a selector that assigns individual fitness values using a
      dictionary *fitness* with genotype at *loci* as keys, and
      fitness as values. Parameter *loci* can be a list of indexes,
      loci names, list of chromosome position pairs, ``ALL_AVAIL``, or
      a function with optional parameter ``pop`` that will be called
      at each ganeeration to determine indexes of loci. For each
      individual (parents if this operator is applied before mating,
      and offspring if this operator is applied during mating),
      genotypes at *loci* are collected one by one (e.g. p0_loc0,
      p1_loc0, p0_loc1, p1_loc1... for a diploid individual, with
      number of alleles varying for sex and mitochondrial DNAs) and
      are looked up in the dictionary. If a genotype cannot be found,
      it will be looked up again without phase information (e.g.
      ``(1,0)`` will match key ``(0,1)``). If the genotype still can
      not be found, a :class:`ValueError` will be raised. This
      operator supports sex chromosomes and haplodiploid populations.
      In these cases, only valid genotypes should be used to generator
      the dictionary keys.



class MaSelector
----------------

.. class:: MaSelector

   This operator is called a 'multi-allele' selector because it groups
   multiple alleles into two groups: wildtype and non-wildtype
   alleles. Alleles in each allele group are assumed to have the same
   effect on individual fitness. If we denote all wildtype alleles as
   ``A``, and all non-wildtype alleles ``a``, this operator assign
   individual fitness according to genotype ``AA``, ``Aa``, ``aa`` in
   the diploid case, and ``A`` and ``a`` in the haploid case.


   .. method:: MaSelector(loci, fitness, wildtype=0, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)


      Creates a multi-allele selector that groups multiple alleles
      into a wildtype group (with alleles *wildtype*, default to
      ``[0]``), and a non-wildtype group. A list of fitness values is
      specified through parameter *fitness*, for genotypes at one or
      more *loci*. Parameter *loci* can be a list of indexes, loci
      names , list of chromosome position pairs, ``ALL_AVAIL``, or a
      function with optional parameter ``pop`` that will be called at
      each ganeeration to determine indexes of loci. If we denote
      wildtype alleles using capital letters ``A``, ``B`` ... and non-
      wildtype alleles using small letters ``a``, ``b`` ..., the
      fitness values should be for

      + genotypes ``A`` and ``a`` for the haploid single-locus case,

      + genotypes ``AB``, ``Ab``, ``aB`` and ``bb`` for haploid
        two=locus cases,

      + genotypes ``AA``, ``Aa`` and ``aa`` for diploid single-locus
        cases,

      + genotypes ``AABB``, ``AABb``, ``AAbb``, ``AaBB``, ``AaBb``,
        ``Aabb``, ``aaBB``, ``aaBb``, and ``aabb`` for diploid two-
        locus cases,

      + and in general 2**n for diploid and 3**n for haploid cases if
        there are ``n`` loci.

      This operator does not support haplodiploid populations, sex and
      mitochondrial chromosomes.



class MlSelector
----------------

.. class:: MlSelector

   This selector is created by a list of selectors. When it is applied
   to an individual, it applies these selectors to the individual,
   obtain a list of fitness values, and compute a combined fitness
   value from them. ADDITIVE, multiplicative, and a heterogeneour
   multi-locus model are supported.


   .. method:: MlSelector(ops, mode=MULTIPLICATIVE, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)


      Create a multiple-locus selector from a list selection operator
      *selectors*. When this operator is applied to an individual
      (parents when used before mating and offspring when used during
      mating), it applies these operators to the individual and obtain
      a list of (usually single-locus) fitness values. These fitness
      values are combined to a single fitness value using

      + *Prod(f_i)*, namely the product of individual fitness if
        *mode* = ``MULTIPLICATIVE``,

      + *1-sum(1 - f_i)* if *mode* = ``ADDITIVE``,

      + *1-Prod(1 - f_i)* if *mode* = ``HETEROGENEITY``, and

      + *exp(- sum(1 - f_i))* if *mode* = ``EXPONENTIAL``,

      zero will be returned if the combined fitness value is less than
      zero.

      Applicability parameters (begin, end, step, at, reps, subPops)
      could be used in both :class:`MlSelector` and selectors in
      parameter *ops*, but parameters in :class:`MlSelector` will be
      interpreted first.



class PySelector
----------------

.. class:: PySelector

   This selector assigns fitness values by calling a user provided
   function. It accepts a list of loci (parameter *loci*) and a Python
   function ``func`` which should be defined with one or more of
   parameters ``geno``, ``mut``, ``gen``, ``ind``, ``pop`` or names of
   information fields. Parameter *loci* can be a list of loci indexes,
   names, list of chromosome position pairs, ``ALL_AVAIL``, or a
   function with optional parameter ``pop`` that will be called at
   each ganeeration to determine indexes of loci. When this operator
   is applied to a population, it passes genotypes or mutants at
   specified loci, generation number, a reference to an individual, a
   reference to the current population (usually used to retrieve
   population variable), and values at specified information fields to
   respective parameters of this function. Genotypes are passed as a
   tuple of alleles arranged locus by locus (in the order of
   A1,A2,B1,B2 for loci A and B). Mutants are passed as a default
   dictionary of loci index (with respect to all genotype of
   individuals, not just the first ploidy) and alleles. The returned
   value will be used to determine the fitness of each individual.


   .. method:: PySelector(func, loci=[], begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, output="", subPops=ALL_AVAIL, infoFields=ALL_AVAIL)


      Create a Python hybrid selector that passes genotype at
      specified *loci*, values at specified information fields (if
      requested) and a generation number to a user-defined function
      *func*. The return value will be treated as individual fitness.



class PyMlSelector
------------------

.. class:: PyMlSelector

   This selector is a multi-locus Python selector that assigns fitness
   to individuals by combining locus and genotype specific fitness
   values. It differs from a :class:`PySelector` in that the python
   function is responsible for assigning fitness values for each
   gentoype type at each locus, which can potentially be random, and
   locus or gentoype-specific.


   .. method:: PyMlSelector(func, mode=EXPONENTIAL, loci=ALL_AVAIL, output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)


      Create a selector that assigns individual fitness values by
      combining locus-specific fitness values that are determined by a
      Python call-back function. The callback function accepts
      parameter *loc*, *alleles* (both optional) and returns location-
      or genotype-specific fitness values that can be constant or
      random. The fitness values for each genotype will be cached so
      the same fitness values will be assigned to genotypes with
      previously assigned values. Note that a function that does not
      examine the genotype naturally assumes a dominant model where
      genotypes with one or two mutants have the same fitness effect.
      Because genotypes at a locus are passed separately and in no
      particular order, this function is also responsible for
      assigning consistent fitness values for genotypes at the same
      locus (a class is usually used). This operator currently ignores
      chromosome types so unused alleles will be passed for loci on
      sex or mitochondrial chromosomes. It also ignores phase of
      genotype so it will use the same fitness value for genotype
      (a,b) and (b,a).

      Individual fitness will be combined in ``ADDITIVE``,
      ``MULTIPLICATIVE``, ``HETEROGENEITY``, or ``EXPONENTIAL`` mode
      from fitness values of loci with at least one non-zero allele
      (See :class:`MlSelector` for details). If an output is given,
      location, genotype, fitness and generation at which the new
      genotype is assgined the value will be written to the output, in
      the format of 'loc a1 a2 fitness gen' for loci on autosomes of
      diploid populations.



