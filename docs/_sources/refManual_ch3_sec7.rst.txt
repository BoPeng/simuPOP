Penetrance
==========

.. index:: single: penetrance


class BasePenetrance
--------------------

.. class:: BasePenetrance

   A penetrance model models the probability that an individual has a
   certain disease provided that he or she has certain genetic
   (genotype) and environmental (information field) riske factors. A
   penetrance operator calculates this probability according to
   provided information and set his or her affection status randomly.
   For example, an individual will have probability 0.8 to be affected
   if the penetrance is 0.8. This class is the base class to all
   penetrance operators and defines a common interface for all
   penetrance operators.

   A penetrance operator can be applied at any stage of an
   evolutionary cycle. If it is applied before or after mating, it
   will set affection status of all parents and offspring,
   respectively. If it is applied during mating, it will set the
   affection status of each offspring. You can also apply a penetrance
   operator to an individual using its ``applyToIndividual`` member
   function.

   By default, a penetrance operator assigns affection status of
   individuals but does not save the actual penetrance value. However,
   if an information field is specified, penetrance values will be
   saved to this field for future analysis.

   When a penetrance operator is applied to a population, it is only
   applied to the current generation. You can, however, use parameter
   *ancGens* to set affection status for all ancestral generations
   (``ALL_AVAIL``), or individuals in specified generations if a list
   of ancestral generations is specified. Note that this parameter is
   ignored if the operator is applied during mating.


   .. method:: BasePenetrance(ancGens=UNSPECIFIED, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a base penetrance operator. This operator assign
      individual affection status in the present generation (default).
      If ``ALL_AVAIL`` or a list of ancestral generations are spcified
      in parameter *ancGens*, individuals in specified ancestral
      generations will be processed. A penetrance operator can be
      applied to specified (virtual) subpopulations (parameter
      *subPops*) and replicates (parameter *reps*). If an informatio
      field is given, penetrance value will be stored in this
      information field of each individual.


   .. method:: BasePenetrance.apply(pop)

      set penetrance to all individuals and record penetrance if
      requested

   .. method:: BasePenetrance.applyToIndividual(ind, pop=None)

      Apply the penetrance operator to a single individual *ind* and
      set his or her affection status. A population reference can be
      passed if the penetrance model depends on population properties
      such as generation number. This function returns the affection
      status.


class MapPenetrance
-------------------

.. class:: MapPenetrance

   This penetrance operator assigns individual affection status using
   a user-specified penetrance dictionary.


   .. method:: MapPenetrance(loci, penetrance, ancGens=UNSPECIFIED, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a penetrance operator that get penetrance value from a
      dictionary *penetrance* with genotype at *loci* as keys, and
      *penetrance* as values. For each individual, genotypes at *loci*
      are collected one by one (e.g. p0_loc0, p1_loc0, p0_loc1,
      p1_loc1... for a diploid individual) and are looked up in the
      dictionary. Parameter *loci* can be a list of loci indexes,
      names, list of chromosome position pairs, ``ALL_AVAIL``, or a
      function with optional parameter ``pop`` that will be called at
      each ganeeration to determine indexes of loci. If a genotype
      cannot be found, it will be looked up again without phase
      information (e.g. ``(1,0)`` will match key ``(0,1)``). If the
      genotype still can not be found, a :class:`ValueError` will be
      raised. This operator supports sex chromosomes and haplodiploid
      populations. In these cases, only valid genotypes should be used
      to generator the dictionary keys.



class MaPenetrance
------------------

.. class:: MaPenetrance

   This operator is called a 'multi-allele' penetrance operator
   because it groups multiple alleles into two groups: wildtype and
   non-wildtype alleles. Alleles in each allele group are assumed to
   have the same effect on individual penetrance. If we denote all
   wildtype alleles as ``A``, and all non-wildtype alleles ``a``, this
   operator assign  Individual penetrance according to genotype
   ``AA``, ``Aa``, ``aa`` in the diploid case, and ``A`` and ``a`` in
   the haploid case.


   .. method:: MaPenetrance(loci, penetrance, wildtype=0, ancGens=UNSPECIFIED, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Creates a multi-allele penetrance operator that groups multiple
      alleles into a wildtype group (with alleles *wildtype*, default
      to ``[0]``), and a non-wildtype group. A list of penetrance
      values is specified through parameter *penetrance*, for
      genotypes at one or more *loci*. Parameter *loci* can be a list
      of loci indexes, names, list of chromosome position pairs,
      ``ALL_AVAIL``, or a function with optional parameter ``pop``
      that will be called at each ganeeration to determine indexes of
      loci. If we denote wildtype alleles using capital letters ``A``,
      ``B`` ... and non-wildtype alleles using small letters ``a``,
      ``b`` ..., the penetrance values should be for

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

      This operator does not support haplodiploid populations and sex
      chromosomes.



class MlPenetrance
------------------

.. class:: MlPenetrance

   This penetrance operator is created by a list of penetrance
   operators. When it is applied to an individual, it applies these
   penetrance operators to the individual, obtain a list of penetrance
   values, and compute a combined penetrance value from them and
   assign affection status accordingly. ADDITIVE, multiplicative, and
   a heterogeneour multi-locus model are supported. Please refer to
   Neil Rish (1989) "Linkage Strategies for

   Genetically Complex Traits" for some analysis of these models.


   .. method:: MlPenetrance(ops, mode=MULTIPLICATIVE, ancGens=UNSPECIFIED, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a multiple-locus penetrance operator from a list
      penetrance operator *ops*. When this operator is applied to an
      individual (parents when used before mating and offspring when
      used during mating), it applies these operators to the
      individual and obtain a list of (usually single-locus)
      penetrance values. These penetrance values are combined to a
      single penetrance value using

      + *Prod(f_i)*, namely the product of individual penetrance if
        *mode* = ``MULTIPLICATIVE``,

      + *sum(f_i)* if *mode* = ``ADDITIVE``, and

      + *1-Prod(1 - f_i)* if *mode* = ``HETEROGENEITY``

      0 or 1 will be returned if the combined penetrance value is less
      than zero or greater than 1.

      Applicability parameters (begin, end, step, at, reps, subPops)
      could be used in both :class:`MlSelector` and selectors in
      parameter *ops*, but parameters in :class:`MlSelector` will be
      interpreted first.



class PyPenetrance
------------------

.. class:: PyPenetrance

   This penetrance operator assigns penetrance values by calling a
   user provided function. It accepts a list of loci (parameter
   ``loci``), and a Python function ``func`` which should be defined
   with one or more of parameters ``geno``, ``mut``, ``gen``, ``ind``,
   ``pop``, or names of information fields. When this operator is
   applied to a population, it passes genotypes or mutants (non-zero
   alleles) at specified loci at specified loci, generation number, a
   reference to an individual, a reference to the current population
   (usually used to retrieve population variables) and values at
   specified information fields to respective parameters of this
   function. Genotypes of each individual are passed as a tuple of
   alleles arranged locus by locus (in the order of A1,A2,B1,B2 for
   loci A and B). Mutants are passed as a default dictionary of loci
   index (with respect to all genotype of individuals, not just the
   first ploidy) and alleles. The returned penetrance values will be
   used to determine the affection status of each individual.


   .. method:: PyPenetrance(func, loci=[], ancGens=UNSPECIFIED, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a Python hybrid penetrance operator that passes genotype
      at specified *loci*, values at specified information fields (if
      requested), and a generation number to a user-defined function
      *func*. Parameter *loci* can be a list of loci indexes, names,
      list of chromosome position pairs, ``ALL_AVAIL``, or a function
      with optional parameter ``pop`` that will be called at each
      ganeeration to determine indexes of loci. The return value will
      be treated as  Individual penetrance.



class PyMlPenetrance
--------------------

.. class:: PyMlPenetrance

   This penetrance operator is a multi-locus Python penetrance
   operator that assigns penetrance values by combining locus and
   genotype specific penetrance values. It differs from a
   :class:`PyPenetrance` in that the python function is responsible
   for penetrance values values for each gentoype type at each locus,
   which can potentially be random, and locus or gentoype-specific.


   .. method:: PyMlPenetrance(func, mode=MULTIPLICATIVE, loci=ALL_AVAIL, ancGens=UNSPECIFIED, output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a penetrance operator that assigns individual affection
      status according to penetrance values combined from locus-
      specific penetrance values that are determined by a Python call-
      back function. The callback function accepts parameter *loc*,
      *alleles* (both optional) and returns location- or genotype-
      specific penetrance values that can be constant or random. The
      penetrance values for each genotype will be cached so the same
      penetrance values will be assigned to genotypes with previously
      assigned values. Note that a function that does not examine the
      genotype naturally assumes a dominant model where genotypes with
      one or two mutants have the same penetrance value. Because
      genotypes at a locus are passed separately and in no particular
      order, this function is also responsible for assigning
      consistent fitness values for genotypes at the same locus (a
      class is usually used). This operator currently ignores
      chromosome types so unused alleles will be passed for loci on
      sex or mitochondrial chromosomes. This operator also ignores the
      phase of genotype so genotypes (a,b) and (b,a) are assumed to
      have the same fitness effect.

      Individual penetrance will be combined in ``ADDITIVE``,
      ``MULTIPLICATIVE``, or ``HETEROGENEITY`` mode from penetrance
      values of loci with at least one non-zero allele (See
      :class:`MlPenetrance` for details).



