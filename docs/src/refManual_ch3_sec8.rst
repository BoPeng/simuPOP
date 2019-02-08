Quantitative Trait
==================

.. index:: single: quantitative trait


class BaseQuanTrait
-------------------

.. class:: BaseQuanTrait

   A quantitative trait in  simuPOP is simply an information field. A
   quantitative trait model simply assigns values to one or more
   information fields (called trait fields) of each individual
   according to its genetic (genotype) and environmental (information
   field) factors. It can be applied at any stage of an evolutionary
   cycle. If a quantitative trait operator is applied before or after
   mating, it will set the trait fields of all parents and offspring.
   If it is applied during mating, it will set the trait fields of
   each offspring.

   When a quantitative trait operator is applied to a population, it
   is only applied to the current generation. You can, however, use
   parameter *ancGen=-1* to set the trait field of all ancestral
   generations, or a generation index to apply to only ancestral
   generation younger than *ancGen*. Note that this parameter is
   ignored if the operator is applied during mating.


   .. method:: BaseQuanTrait(ancGens=UNSPECIFIED, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a base quantitative trait operator. This operator assigns
      one or more quantitative traits to trait fields in the present
      generation (default). If ``ALL_AVAIL`` or a list of ancestral
      generations are specified, this operator will be applied to
      individuals in these generations as well. A quantitative trait
      operator can be applied to specified (virtual) subpopulations
      (parameter *subPops*) and replicates (parameter *reps*).


   .. method:: BaseQuanTrait.apply(pop)

      set ``qtrait`` to all individual


class PyQuanTrait
-----------------

.. class:: PyQuanTrait

   This quantitative trait operator assigns a trait field by calling a
   user provided function. It accepts a list of loci (parameter
   *loci*), and a Python function ``func`` which should be defined
   with one or more of parameters ``geno``, ``mut``, ``gen``, ``ind``,
   or names of information fields. When this operator is applied to a
   population, it passes genotypes or mutants (non-zero alleles) of
   each individual at specified loci, generation number, a reference
   to an individual, and values at specified information fields to
   respective parameters of this function. Genotypes of each
   individual are passed as a tuple of alleles arranged locus by locus
   (in the order of A1,A2,B1,B2 for loci A and B). Mutants are passed
   as a default dictionary of loci index (with respect to all genotype
   of individuals, not just the first ploidy) and alleles. The return
   values will be assigned to specified trait fields.


   .. method:: PyQuanTrait(func, loci=[], ancGens=UNSPECIFIED, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a Python hybrid quantitative trait operator that passes
      genotype at specified *loci*, optional values at specified
      information fields (if requested), and an optional generation
      number to a user-defined function *func*. Parameter *loci* can
      be a list of loci indexes, names, or ``ALL_AVAIL``. The return
      value will be assigned to specified trait fields (*infoField*).
      If only one trait field is specified, a number or a sequence of
      one element is acceptable. Otherwise, a sequence of values will
      be accepted and be assigned to each trait field.



