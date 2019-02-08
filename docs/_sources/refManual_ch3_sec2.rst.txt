Initialization
==============

.. index:: single: initializer


class InitSex
-------------

.. class:: InitSex

   This operator initializes sex of individuals, either randomly or
   use a list of sexes.


   .. method:: InitSex(maleFreq=0.5, maleProp=-1, sex=[], begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create an operator that initializes individual sex to ``MALE``
      or ``FEMALE``. By default, it assigns sex to individuals
      randomly, with equal probability of having a male or a female.
      This probabability can be adjusted through parameter *maleFreq*
      or be made to exact proportions by specifying parameter
      *maleProp*. Alternatively, a fixed sequence of sexes can be
      assigned. For example, if ``sex=[MALE, FEMALE]``, individuals
      will be assigned ``MALE`` and ``FEMALE`` successively. Parameter
      *maleFreq* or *maleProp* are ignored if *sex* is given. If a
      list of (virtual) subpopulation is specified in parameter
      *subPop*, only individuals in these subpopulations will be
      initialized. Note that the *sex* sequence, if used, is assigned
      repeatedly regardless of (virtual) subpopulation boundaries so
      that you can assign *sex* to all individuals in a population.



class InitInfo
--------------

.. class:: InitInfo

   This operator initializes given information fields with a sequence
   of values, or a user-provided function such as ``random.random``.


   .. method:: InitInfo(values, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create an operator that initialize individual information fields
      *infoFields* using a sequence of values or a user-defined
      function. If a list of values are given, it will be used
      sequentially for all individuals. The values will be reused if
      its length is less than the number of individuals. The values
      will be assigned repeatedly regardless of subpopulation
      boundaries. If a Python function is given, it will be called,
      without any argument, whenever a value is needed. If a list of
      (virtual) subpopulation is specified in parameter *subPop*, only
      individuals in these subpopulations will be initialized.



class InitGenotype
------------------

.. class:: InitGenotype

   This operator assigns alleles at all or part of loci with given
   allele frequencies, proportions or values. This operator
   initializes all chromosomes, including unused genotype locations
   and customized chromosomes.


   .. method:: InitGenotype(freq=[], genotype=[], prop=[], haplotypes=[], genotypes=[], loci=ALL_AVAIL, ploidy=ALL_AVAIL, begin=0, end=1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      This function creates an initializer that initializes individual
      genotypes with random alleles, genotypes, or haplotypes with
      specified frequencies (parameter *freq*) or proportions
      (parameter *prop*). If parameter *genotypes* or *haplotypes* is
      not specified, *freq* specifies the allele frequencies of
      alleles ``0``, ``1``, ``2``... respectively. Alternatively, you
      can use parameter *prop* to specified the exact proportions of
      alleles ``0``, ``1``, ..., although alleles with small
      proportions might not be assigned at all.

      Values of parameter *prob* or *prop* should add up to 1. In
      addition to a vector, parameter *prob* and *prop* can also be a
      function that accepts optional parameters *loc*, *subPop* or
      *vsp* and returns a list of requencies for alleles ``0``, ``1``,
      etc, or a number for frequency of allele ``0`` as a speciail
      case for each locus, subpopulation (parameter *subPop*), or
      virtual subpopulations (parameter *vsp*, pass as a tuple).

      If parameter *genotypes* is specified, it should contain a list
      of genotypes (alleles on different strand of chromosomes) with
      length equal to population ploidy. Parameter *prob* and *prop*
      then specifies frequencies or proportions of each genotype,
      which can vary for each subpopulation but not each locus if the
      function form of parameters is used.

      If parameter *haplotypes* is specified, it should contain a list
      of haplotypes (alleles on the same strand of chromosome) and
      parameter *prob* or *prop* specifies frequencies or proportions
      of each haplotype.

      If *loci*, *ploidy* and/or *subPop* are specified, only
      specified loci, ploidy, and individuals in these (virtual)
      subpopulations will be initialized. Parameter *loci* can be a
      list of loci indexes, names or ``ALL_AVAIL``. If the length of a
      haplotype is not enough to fill all loci, the haplotype will be
      reused. If a list (or a single) haplotypes are specified without
      *freq* or *prop*, they are used with equal probability.

      In the last case, if a sequence of genotype is specified through
      parameter *genotype* (not *genotypes*), it will be used
      repeatedly to initialize all alleles sequentially. This works
      similar to function ``Population.setGenotype()`` except that you
      can limit the initialization to certain *loci* and *ploidy*.



class InitLineage
-----------------

.. class:: InitLineage

   This operator assigns lineages at all or part of loci with given
   values. This operator initializes all chromosomes, including unused
   lineage locations and customized chromosomes.


   .. method:: InitLineage(lineage=[], mode=PER_ALLELE, loci=ALL_AVAIL, ploidy=ALL_AVAIL, begin=0, end=1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=["ind_id"])


      This function creates an initializer that initializes lineages
      with either a specified set of values or from the field
      *infoFields* (default to ``ind_id``), whose value will be saved
      as the lineage of modified alleles. If a list of values is
      specified in parameter *lineage*, each value in this list is
      applied to one or more alleles so that each allele
      (``PER_ALLELE``, default mode), alleles on each chromosome
      (``PER_CHROMOSOME``), on chromosomes of each ploidy
      (``PER_PLOIDY``), or for each individual (``PER_INDIVIDUAL``)
      have the same lineage. A single value is allowed and values in
      *lineage* will be re-used if not enough values are provided. If
      an empty list is provided, values 1, 2, 3, .. will be used to
      provide an unique identify for each allele, genotype,
      chromosome, etc. If a valid field is specified (default to
      ``ind_id``), the value of this field will be used for all
      alleles of each individual if *mode* is set to ``FROM_INFO``, or
      be adjusted to produce positive values for alleles on the frist
      ploidy, and negative values for the second ploidy (and so on) if
      *mode* equals to ``FROM_INFO_SIGNED``. If *loci*, *ploidy*
      and/or *subPops* are specified, only specified loci, ploidy, and
      individuals in these (virtual) subpopulations will be
      initialized.



