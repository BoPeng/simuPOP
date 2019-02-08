Individual, Population, pedigree and Simulator
==============================================

.. index::
   single: Individual
   single: Population
   single: Pedigree
   single: Simulator


class GenoStruTrait
-------------------

.. class:: GenoStruTrait

   All individuals in a population share the same genotypic properties
   such as number of chromosomes, number and position of loci, names
   of markers, chromosomes, and information fields. These properties
   are stored in this :class:`GenoStruTrait` class and are accessible
   from both :class:`Individual` and :class:`Population` classes.
   Currently, a genotypic structure consists of

   + Ploidy, namely the number of homologous sets of chromosomes, of a
     population. Haplodiploid population is also supported.

   + Number of chromosomes and number of loci on each chromosome.

   + Positions of loci, which determine the relative distance between
     loci on the same chromosome. No unit is assumed so these
     positions can be ordinal (``1``, ``2``, ``3``, ..., the default),
     in physical distance (``bp``, ``kb`` or ``mb``), or in map
     distance (e.g. ``centiMorgan``) depending on applications.

   + Names of alleles, which can either be shared by all loci or be
     specified for each locus.

   + Names of loci and chromosomes.

   + Names of information fields attached to each individual.

   In addition to basic property access functions, this class provides
   some utility functions such as ``locusByName``, which looks up a
   locus by its name.


   .. method:: GenoStruTrait()


      A :class:`GenoStruTrait` object is created with the construction
      of a :class:`Population` object and cannot be initialized
      directly.


   .. method:: GenoStruTrait.absLocusIndex(chrom, locus)

      return the absolute index of locus *locus* on chromosome
      *chrom*. c.f. ``chromLocusPair``.

   .. method:: GenoStruTrait.alleleName(allele, locus=0)

      return the name of allele *allele* at *lcous* specified by the
      *alleleNames* parameter of the :class:`Population` function.
      *locus* could be ignored if alleles at all loci share the same
      names. If the name of an allele is unspecified, its value
      (``'0'``, ``'1'``, ``'2'``, etc) is returned.

   .. method:: GenoStruTrait.alleleNames(locus=0)

      return a list of allele names at \locus given by the
      *alleleNames* parameter of the :class:`Population` function.
      *locus* could be ignored if alleles at all loci share the same
      names. This list does not have to cover all possible allele
      states of a population so ``alleleNames()[``*allele*``]`` might
      fail (use ``alleleNames(``*allele*``)`` instead).

   .. method:: GenoStruTrait.chromBegin(chrom)

      return the index of the first locus on chromosome *chrom*.

   .. method:: GenoStruTrait.chromByName(name)

      return the index of a chromosome by its *name*.

   .. method:: GenoStruTrait.chromEnd(chrom)

      return the index of the last locus on chromosome *chrom* plus 1.

   .. method:: GenoStruTrait.chromLocusPair(locus)

      return the chromosome and relative index of a locus using its
      absolute index *locus*. c.f. ``absLocusIndex``.

   .. method:: GenoStruTrait.chromName(chrom)

      return the name of a chromosome *chrom*.

   .. method:: GenoStruTrait.chromNames()

      return a list of the names of all chromosomes.

   .. method:: GenoStruTrait.chromType(chrom)

      return the type of a chromosome *chrom* (``CUSTOMIZED``,
      ``AUTOSOME``, ``CHROMOSOME_X``, ``CHROMOSOME_Y`` or
      ``MITOCHONDRIAL``.

   .. method:: GenoStruTrait.chromTypes()

      return the type of all chromosomes (``CUSTOMIZED``,
      ``AUTOSOME``, ``CHROMOSOME_X``, ``CHROMOSOME_Y``, or
      ``MITOCHONDRIAL``).

   .. method:: GenoStruTrait.indexesOfLoci(loci=ALL_AVAIL)

      return the indexes of loci with positions *positions* (list of
      (chr, pos) pairs). Raise a :class:`ValueError` if any of the
      loci cannot be found.

   .. method:: GenoStruTrait.infoField(idx)

      return the name of information field *idx*.

   .. method:: GenoStruTrait.infoFields()

      return a list of the names of all information fields of the
      population.

   .. method:: GenoStruTrait.infoIdx(name)

      return the index of information field *name*. Raise an
      :class:`IndexError` if *name* is not one of the information
      fields.

   .. method:: GenoStruTrait.lociByNames(names)

      return the indexes of loci with names *names*. Raise a
      :class:`ValueError` if any of the loci cannot be found.

   .. method:: GenoStruTrait.lociDist(locus1, locus2)

      Return the distance between loci *locus1* and *locus2* on the
      same chromosome. A negative value will be returned if *locus1*
      is after *locus2*.

   .. method:: GenoStruTrait.lociNames()

      return the names of all loci specified by the *lociNames*
      parameter of the :class:`Population` function. An empty list
      will be returned if *lociNames* was not specified.

   .. method:: GenoStruTrait.lociPos()

      return the positions of all loci, specified by the *lociPos*
      prameter of the :class:`Population` function. The default
      positions are 1, 2, 3, 4, ... on each chromosome.

   .. method:: GenoStruTrait.locusByName(name)

      return the index of a locus with name *name*. Raise a
      :class:`ValueError` if no locus is found. Note that empty
      strings are used for loci without name but you cannot lookup
      such loci using this function.

   .. method:: GenoStruTrait.locusName(locus)

      return the name of locus *locus* specified by the *lociNames*
      parameter of the :class:`Population` function. An empty string
      will be returned if no name has been given to locus *locus*.

   .. method:: GenoStruTrait.locusPos(locus)

      return the position of locus *locus* specified by the *lociPos*
      parameter of the :class:`Population` function.

   .. method:: GenoStruTrait.numChrom()

      return the number of chromosomes.

   .. method:: GenoStruTrait.numLoci(chrom)

      return the number of loci on chromosome *chrom*.

   .. method:: GenoStruTrait.numLoci()

      return a list of the number of loci on all chromosomes.

   .. method:: GenoStruTrait.ploidy()

      return the number of homologous sets of chromosomes, specified
      by the *ploidy* parameter of the :class:`Population` function.
      Return 2 for a haplodiploid population because two sets of
      chromosomes are stored for both males and females in such a
      population.

   .. method:: GenoStruTrait.ploidyName()

      return the ploidy name of this population, can be one of
      ``haploid``, ``diploid``, ``haplodiploid``, ``triploid``,
      ``tetraploid`` or ``#-ploid`` where ``#`` is the ploidy number.

   .. method:: GenoStruTrait.totNumLoci()

      return the total number of loci on all chromosomes.


class Individual
----------------

.. class:: Individual

   A :class:`Population` consists of individuals with the same
   genotypic structure. An :class:`Individual` object cannot be
   created independently, but refences to inidividuals can be
   retrieved using member functions of a :class:`Population` object.
   In addition to structural information shared by all individuals in
   a population (provided by class :class:`GenoStruTrait`), the
   :class:`Individual` class provides member functions to get and set
   *genotype*, *sex*, *affection status* and *information fields* of
   an individual.

   Genotypes of an individual are stored sequentially and can be
   accessed locus by locus, or in batch. The alleles are arranged by
   position, chromosome and ploidy. That is to say, the first allele
   on the first chromosome of the first homologous set is followed by
   alleles at other loci on the same chromsome, then markers on the
   second and later chromosomes, followed by alleles on the second
   homologous set of the chromosomes for a diploid individual. A
   consequence of this memory layout is that alleles at the same locus
   of a non-haploid individual are separated by
   ``Individual::totNumLoci()`` loci. It is worth noting that access
   to invalid chromosomes, such as the Y chromosomes of female
   individuals, is not restricted.


   .. method:: Individual()


      An :class:`Individual` object cannot be created directly. It has
      to be accessed from a :class:`Population` object using functions
      such as ``Population::Individual(idx)``.


   .. method:: Individual.affected()

      Return ``True`` if this individual is affected.

   .. method:: Individual.allele(idx, ploidy=-1, chrom=-1)

      return the current allele at a locus, using its absolute index
      *idx*. If a ploidy *ploidy* and/or a chromosome indexes is
      given, *idx* is relative to the beginning of specified
      homologous copy of chromosomes (if *chrom=-1*) or the beginning
      of the specified homologous copy of specified chromosome (if
      *chrom* >= 0).

   .. method:: Individual.alleleChar(idx, ploidy=-1, chrom=-1)

      return the name of ``allele(idx, ploidy, chrom)``. If idx is
      invalid (e.g. second homologus copy of chromosome Y), '_' is
      returned.

   .. method:: Individual.alleleLineage(idx, ploidy=-1, chrom=-1)

      return the lineage of the allele at a locus, using its absolute
      index *idx*. If a ploidy *ploidy* and/or a chromosome indexes is
      given, *idx* is relative to the beginning of specified
      homologous copy of chromosomes (if *chrom=-1*) or the beginning
      of the specified homologous copy of specified chromosome (if
      *chrom* >= 0). **This function returns 0 for modules without
      lineage information.**

   .. method:: Individual.__cmp__(rhs)

      a python function used to compare the individual objects

   .. method:: Individual.genotype(ploidy=ALL_AVAIL, chroms=ALL_AVAIL)

      return an editable array (a ``carray`` object) that represents
      all alleles of an individual. If *ploidy* or *chroms* is given,
      only alleles on the specified chromosomes and homologous copy of
      chromosomes will be returned. If multiple chromosomes are
      specified, there should not be gaps between chromosomes. This
      function ignores type of chromosomes so it will return unused
      alleles for sex and mitochondrial chromosomes.

   .. method:: Individual.info(field)

      Return the value of an information field *filed* (by index or
      name). ``ind.info(name)`` is equivalent to ``ind.name`` although
      the function form allows the use of indexes of information
      fieldes.

   .. method:: Individual.lineage(ploidy=ALL_AVAIL, chroms=ALL_AVAIL)

      return an editable array (a ``carray_lineage`` object) that
      represents the lineages of all alleles of an individual. If
      *ploidy* or *chroms* is given, only lineages on the specified
      chromosomes and homologous copy of chromosomes will be returned.
      If multiple chromosomes are specified, there should not be gaps
      between chromosomes. This function ignores type of chromosomes
      so it will return lineage of unused alleles for sex and
      mitochondrial chromosomes. A ``None`` object will be returned
      for modules without lineage information.

   .. method:: Individual.mutants(ploidy=ALL_AVAIL, chroms=ALL_AVAIL)

      return an itertor that iterate through all mutants (non-zero
      alleles) of an individual. Each mutant is presented as a tuple
      of (index, value) where index is the index of mutant ranging
      from zero to  totNumLoci() *  ploidy() - 1, so you will have to
      adjust indexes to check multiple alleles at a locus. If *ploidy*
      or *chroms* is given, only alleles on the specified chromosomes
      and homologous copy of chromosomes will be iterated. If multiple
      chromosomes are specified, there should not be gaps between
      chromosomes. This function ignores type of chromosomes so it
      will return unused alleles for sex and mitochondrial
      chromosomes.

   .. method:: Individual.setAffected(affected)

      set affection status to *affected* (``True`` or ``False``).

   .. method:: Individual.setAllele(allele, idx, ploidy=-1, chrom=-1)

      set allele *allele* to a locus, using its absolute index *idx*.
      If a ploidy *ploidy* and/or a chromosome indexes are given,
      *idx* is relative to the beginning of specified homologous copy
      of chromosomes (if *chrom=-1*) or the beginning of the specified
      homologous copy of specified chromosome (if *chrom* >= 0).

   .. method:: Individual.setAlleleLineage(lineage, idx, ploidy=-1, chrom=-1)

      set lineage *lineage* to an allele, using its absolute index
      *idx*. If a ploidy *ploidy* and/or a chromosome indexes are
      given, *idx* is relative to the beginning of specified
      homologous copy of chromosomes (if *chrom=-1*) or the beginning
      of the specified homologous copy of specified chromosome (if
      *chrom* >= 0). This function does nothing for modules without
      lineage information.

   .. method:: Individual.setGenotype(geno, ploidy=ALL_AVAIL, chroms=ALL_AVAIL)

      Fill the genotype of an individual using a list of alleles
      *geno*. If parameters *ploidy* and/or *chroms* are specified,
      alleles will be copied to only all or specified chromosomes on
      selected homologous copies of chromosomes. ``geno`` will be
      reused if its length is less than number of alleles to be
      filled. This function ignores type of chromosomes so it will set
      genotype for unused alleles for sex and mitochondrial
      chromosomes.

   .. method:: Individual.setInfo(value, field)

      set the value of an information field *field* (by index or name)
      to *value*. ``ind.setInfo(value, field)`` is equivalent to
      ``ind.field = value`` although the function form allows the use
      of indexes of information fieldes.

   .. method:: Individual.setLineage(lineage, ploidy=ALL_AVAIL, chroms=ALL_AVAIL)

      Fill the lineage of an individual using a list of IDs *lineage*.
      If parameters *ploidy* and/or *chroms* are specified, lineages
      will be copied to only all or specified chromosomes on selected
      homologous copies of chromosomes. ``lineage`` will be reused if
      its length is less than number of allelic lineage to be filled.
      This function ignores type of chromosomes so it will set lineage
      to unused alleles for sex and mitochondrial chromosomes. It does
      nothing for modules without lineage information.

   .. method:: Individual.setSex(sex)

      set individual sex to ``MALE`` or ``FEMALE``.

   .. method:: Individual.sex()

      return the sex of an individual, ``1`` for male and ``2`` for
      female.


class Population
----------------

.. class:: Population

   A  simuPOP population consists of individuals of the same genotypic
   structure, organized by generations, subpopulations and virtual
   subpopulations. It also contains a Python dictionary that is used
   to store arbitrary population variables.

   In addition to genotypic structured related functions provided by
   the :class:`GenoStruTrait` class, the population class provides a
   large number of member functions that can be used to

   + Create, copy and compare populations.

   + Manipulate subpopulations. A population can be divided into
     several subpopulations. Because individuals only mate with
     individuals within the same subpopulation, exchange of genetic
     information across subpopulations can only be done through
     migration. A number of functions are provided to access
     subpopulation structure information, and to merge and split
     subpopulations.

   + Define and access virtual subpopulations. A *virtual
     subpopulation splitter* can be assigned to a population, which
     defines groups of individuals called *virtual subpopulations*
     (VSP) within each subpopulation.

   + Access individuals individually, or through iterators that
     iterate through individuals in (virtual) subpopulations.

   + Access genotype and information fields of individuals at the
     population level. From a population point of view, all genotypes
     are arranged sequentially individual by individual. Please refer
     to class :class:`Individual` for an introduction to genotype
     arragement of each individual.

   + Store and access *ancestral generations*. A population can save
     arbitrary number of ancestral generations. It is possible to
     directly access an ancestor, or make an ancestral generation the
     current generation for more efficient access.

   + Insert or remove loci, resize (shrink or expand) a population,
     sample from a population, or merge with other populations.

   + Manipulate population variables and evaluate expressions in this
     *local namespace*.

   + Save and load a population.


   .. method:: Population(size=[], ploidy=2, loci=[], chromTypes=[], lociPos=[], ancGen=0, chromNames=[], alleleNames=[], lociNames=[], subPopNames=[], infoFields=[])


      The following parameters are used to create a population object:

      size
         A list of subpopulation sizes. The length of this list determines the
         number of subpopulations of this population. If there is no
         subpopulation, ``size=[popSize]`` can be written as
         ``size=popSize``.

      ploidy
         Number of homologous sets of chromosomes. Default to ``2`` (diploid).
         For efficiency considerations, all chromosomes have the same
         number of homologous sets, even if some customized
         chromosomes or some individuals (e.g. males in a haplodiploid
         population) have different numbers of homologous sets. The
         first case is handled by setting *chromTypes* of each
         chromosome. Only the haplodiploid populations are handled for
         the second case, for which ``ploidy=HAPLODIPLOID`` should be
         used.

      loci
         A list of numbers of loci on each chromosome. The length of this
         parameter determines the number of chromosomes. If there is
         only one chromosome, ``numLoci`` instead of ``[numLoci]`` can
         be used.

      chromTypes
         A list that specifies the type of each chromosome, which can be
         ``AUTOSOME``, ``CHROMOSOME_X``, ``CHROMOSOME_Y``, or
         ``CUSTOMIZED``. All chromosomes are assumed to be autosomes
         if this parameter is ignored. Sex chromosome can only be
         specified in a diploid population where the sex of an
         individual is determined by the existence of these
         chromosomes using the ``XX`` (``FEMALE``) and ``XY``
         (``MALE``) convention. Both sex chromosomes have to be
         available and be specified only once. Because chromosomes
         ``X`` and ``Y`` are treated as two chromosomes, recombination
         on the pseudo-autosomal regions of the sex chromsomes is not
         supported. ``CUSTOMIZED`` chromosomes are special chromosomes
         whose inheritance patterns are undefined. They rely on user-
         defined functions and operators to be passed from parents to
         offspring. Multiple customized chromosomes have to be
         arranged consecutively.

      lociPos
         Positions of all loci on all chromosome, as a list of float numbers.
         Default to ``1``, ``2``, ... etc on each chromosome.
         *lociPos* should be arranged chromosome by chromosome. If
         ``lociPos`` are not in order within a chromosome, they will
         be re-arranged along with corresponding *lociNames* (if
         specified).

      ancGen
         Number of the most recent ancestral generations to keep during
         evolution. Default to ``0``, which means only the current
         generation will be kept. If it is set to ``-1``, all
         ancestral generations will be kept in this population (and
         exhaust your computer RAM quickly).

      chromNames
         A list of chromosome names. Default to ``''`` for all chromosomes.

      alleleNames
         A list or a nested list of allele names. If a list of alleles is
         given, it will be used for all loci in this population. For
         example, ``alleleNames=('A','C','T','G')`` gives names ``A``,
         ``C``, ``T``, and ``G`` to alleles ``0``, ``1``, ``2``, and
         ``3`` respectively. If a nested list of names is given, it
         should specify alleles names for all loci.

      lociNames
         A list of names for each locus. It can be empty or a list of unique
         names for each locus. If loci are not specified in order,
         loci names will be rearranged according to their position on
         the chromosome.

      subPopNames
         A list of subpopulation names. All subpopulations will have name
         ``''`` if this parameter is not specified.

      infoFields
         Names of information fields (named float number) that will be attached
         to each individual.


   .. method:: Population.absIndIndex(idx, subPop)

      return the absolute index of an individual *idx* in
      subpopulation *subPop*.

   .. method:: Population.addChrom(lociPos, lociNames=[], chromName="", alleleNames=[], chromType=AUTOSOME)

      Add chromosome *chromName* with given type *chromType* to a
      population, with loci *lociNames* inserted at position
      *lociPos*. *lociPos* should be ordered. *lociNames* and
      *chromName* should not exist in the current population. Allele
      names could be specified for all loci (a list of names) or
      differently for each locus (a nested list of names), using
      parameter *alleleNames*. Empty loci names will be used if
      *lociNames* is not specified. The newly added alleles will have
      zero lineage in modules wiht lineage information.

   .. method:: Population.addChromFrom(pop)

      Add chromosomes in population *pop* to the current population.
      population *pop* should have the same number of individuals as
      the current population in the current and all ancestral
      generations. Chromosomes of *pop*, if named, should not conflict
      with names of existing chromosome. This function merges
      genotypes on the new chromosomes from population ``pop``
      individual by individual.

   .. method:: Population.addIndFrom(pop)

      Add all individuals, including ancestors, in *pop* to the
      current population. Two populations should have the same
      genotypic structures and number of ancestral generations.
      Subpopulations in population *pop* are kept.

   .. method:: Population.addInfoFields(fields, init=0)

      Add a list of information fields *fields* to a population and
      initialize their values to *init*. If an information field
      alreay exists, it will be re-initialized.

   .. method:: Population.addLoci(chrom, pos, lociNames=[], alleleNames=[])

      Insert loci *lociNames* at positions *pos* on chromosome
      *chrom*. These parameters should be lists of the same length,
      although *names* may be ignored, in which case empty strings
      will be assumed. Single-value input is allowed for parameter
      *chrom* and *pos* if only one locus is added. Alleles at
      inserted loci are initialized with zero alleles. Note that loci
      have to be added to existing chromosomes. If loci on a new
      chromosome need to be added, function ``addChrom`` should be
      used. Optionally, allele names could be specified either for all
      loci (a single list) or each loci (a nested list). This function
      returns indexes of the inserted loci. Newly inserted alleles
      will have zero lineage in modules with lineage information.

   .. method:: Population.addLociFrom(pop, byName=False)

      Add loci from population *pop*. By default, chromosomes are
      merged by index and names of merged chromosomes of population
      *pop* will be ignored (merge of two chromosomes with different
      names will yield a warning). If *byName* is set to ``True``,
      chromosomes in *pop* will be merged to chromosomes with
      identical names. Added loci will be inserted according to their
      position. Their position and names should not overlap with any
      locus in the current population. population *pop* should have
      the same number of individuals as the current population in the
      current and all ancestral generations. Allele lineages are also
      copied from *pop* in modules with lineage information.

   .. method:: Population.ancestor(idx, gen, subPop=[])

      Return a reference to individual ``idx`` in ancestral generation
      ``gen``. The correct individual will be returned even if the
      current generation is not the present one (see also
      ``useAncestralGen``). If a valid *subPop* is specified, *index*
      is relative to that *subPop*. Virtual subpopulation is not
      supported. Note that a float *idx* is acceptable as long as it
      rounds closely to an integer.

   .. method:: Population.ancestralGens()

      Return the actual number of ancestral generations stored in a
      population, which does not necessarily equal to the number set
      by ``setAncestralDepth()``.

   .. method:: Population.clone()

      Create a cloned copy of a population. Note that Python statement
      ``pop1 = pop`` only creates a reference to an existing
      population ``pop``.

   .. method:: Population.__cmp__(rhs)

      a python function used to compare the population objects

   .. method:: Population.dvars(subPop=[])

      Return a wrapper of Python dictionary returned by
      ``vars(subPop)`` so that dictionary keys can be accessed as
      attributes.

   .. method:: Population.extractIndividuals(indexes=[], IDs=[], idField="ind_id", filter=None)

      Extract individuals with given absolute indexes (parameter
      *indexes*), IDs (parameter *IDs*, stored in information field
      *idField*, default to ``ind_id``), or a filter function
      (parameter *filter*). If a list of absolute indexes are
      specified, the present generation will be extracted and form a
      one-generational population. If a list of IDs are specified,
      this function will look through all ancestral generations and
      extract individuals with given ID. Individuals with shared IDs
      are allowed. In the last case, a user-defined Python function
      should be provided. This function should accept parameter
      ``"ind"`` or one or more of the information fields. All
      individuals, including ancestors if there are multiple ancestral
      generations, will be passed to this function. Individuals that
      returns ``True`` will be extracted. Extracted individuals will
      be in their original ancestral generations and subpopulations,
      even if some subpopulations or generations are empty. An
      :class:`IndexError` will be raised if an index is out of bound
      but no error will be given if an invalid ID is encountered.

   .. method:: Population.extractSubPops(subPops=ALL_AVAIL, rearrange=False)

      Extract a list of (virtual) subpopulations from a population and
      create a new population. If *rearrange* is ``False`` (default),
      structure and names of extracted subpopulations are kept
      although extracted subpopulations can have fewer individuals if
      they are created from extracted virtual subpopulations. (e.g. it
      is possible to extract all male individuals from a subpopulation
      using a ``SexSplitter()``). If *rearrange* is ``True``, each
      (virtual) subpopulation in *subPops* becomes a new subpopulation
      in the extracted population in the order at which they are
      specified. Because each virtual subpopulation becomes a
      subpopulation, this function could be used, for example, to
      separate male and female individuals to two subpopulations (
      ``subPops=[(0,0), (0,1)]``). If overlapping (virtual)
      subpopulations are specified, individuals will be copied
      multiple times. This function only extract individuals from the
      present generation.

   .. method:: Population.genotype(subPop=[])

      Return an editable array of the genotype of all individuals in a
      population (if ``subPop=[]``, default), or individuals in a
      subpopulation *subPop*. Virtual subpopulation is unsupported.

   .. method:: Population.indByID(id, ancGens=ALL_AVAIL, idField="ind_id")

      Return a reference to individual with *id* stored in information
      field *idField* (default to ``ind_id``). This function by
      default search the present and all ancestral generations
      (``ancGen=ALL_AVAIL``), but you can limit the search in specific
      generations if you know which generations to search
      (``ancGens=[0,1]`` for present and parental generations) or
      ``UNSPECIFIED`` to search only the current generation. If no
      individual with *id* is found, an :class:`IndexError` will be
      raised. A float *id* is acceptable as long as it rounds closely
      to an integer. Note that this function uses a dynamic searching
      algorithm which tends to be slow. If you need to look for
      multiple individuals from a static population, you might want to
      convert a population object to a pedigree object and use
      function :meth:`Pedigree.indByID`.

   .. method:: Population.indInfo(field, subPop=[])

      Return the values (as a list) of information field ``field`` (by
      index or name) of all individuals (if ``subPop=[]``, default),
      or individuals in a (virtual) subpopulation (if ``subPop=sp`` or
      ``(sp, vsp)``).

   .. method:: Population.individual(idx, subPop=[])

      Return a refernce to individual *idx* in the population (if
      ``subPop=[]``, default) or a subpopulation (if ``subPop=sp``).
      Virtual subpopulation is not supported. Note that a float *idx*
      is acceptable as long as it rounds closely to an integer.

   .. method:: Population.individuals(subPop=[])

      Return an iterator that can be used to iterate through all
      individuals in a population (if ``subPop=[]``, default), or a
      (virtual) subpopulation (``subPop=spID`` or ``(spID,  vspID)``).
      If you would like to iterate through multiple subpopulations in
      multiple ancestral generations, please use function
      ``Population.allIndividuals()``.

   .. method:: Population.lineage(subPop=[])

      Return an editable array of the lineage of alleles for all
      individuals in a population (if ``subPop=[]``, default), or
      individuals in a subpopulation *subPop*. Virtual subpopulation
      is unsupported. **This function returns ``None`` for modules
      without lineage information.**

   .. method:: Population.mergeSubPops(subPops=ALL_AVAIL, name="", toSubPop=-1)

      Merge subpopulations *subPops*. If *subPops* is ``ALL_AVAIL``
      (default), all subpopulations will be merged. *subPops* do not
      have to be adjacent to each other. They will all be merged to
      the subpopulation with the smallest subpopulation ID, unless a
      subpopulation ID is specified using parameter ``toSubPop``.
      Indexes of the rest of the subpopulation may be changed. A new
      name can be assigned to the merged subpopulation through
      parameter *name* (an empty *name* will be ignored). This
      function returns the ID of the merged subpopulation.

   .. method:: Population.mutants(subPop=[])

      Return an iterator that iterate through mutants of all
      individuals in a population (if ``subPop=[]``, default), or
      individuals in a subpopulation *subPop*. Virtual subpopulation
      is unsupported. Each mutant is presented as a tuple of (index,
      value) where index is the index of mutant (from 0 to
      totNumLoci()*ploidy()) so you will have to adjust its value to
      check multiple alleles at a locus. This function ignores type of
      chromosomes so non-zero alleles in unused alleles of sex and
      mitochondrial chromosomes are also iterated.

   .. method:: Population.numSubPop()

      Return the number of subpopulations in a population. Return 1 if
      there is no subpopulation structure.

   .. method:: Population.numVirtualSubPop()

      Return the number of virtual subpopulations (VSP) defined by a
      VSP splitter. Return ``0`` if no VSP is defined.

   .. method:: Population.popSize(ancGen=-1, sex=ANY_SEX)

      Return the total number of individuals in all subpopulations of
      the current generation (default) or the an ancestral generation
      *ancGen*. This function by default returns number of all
      individuals (``sex=ANY_SEX``), but it will return number of
      males (if ``sex=MALE_ONLY``), number of females (if
      ``sex=MALE_ONLY``), and number of male/female pairs (if
      ``sex=PAIR_ONLY``) which is essentially less of the number of
      males and females.

   .. method:: Population.push(pop)

      Push population *pop* into the current population. Both
      populations should have the same genotypic structure. The
      current population is discarded if *ancestralDepth* (maximum
      number of ancestral generations to hold) is zero so no ancestral
      generation can be kept. Otherise, the current population will
      become the parental generation of *pop*. If *ancGen* of a
      population is positive and there are already *ancGen* ancestral
      generations (c.f. ``ancestralGens()``), the greatest ancestral
      generation will be discarded. In any case,  Population*pop*
      becomes invalid as all its individuals are absorbed by the
      current population.

   .. method:: Population.recodeAlleles(alleles, loci=ALL_AVAIL, alleleNames=[])

      Recode alleles at *loci* (can be a list of loci indexes or
      names, or all loci in a population (``ALL_AVAIL``)) to other
      values according to parameter *alleles*. This parameter can a
      list of new allele numbers for alleles ``0``, ``1``, ``2``, ...
      (allele ``x`` will be recoded to ``newAlleles[x]``, ``x``
      outside of the range of *newAlleles* will not be recoded,
      although a warning will be given if ``DBG_WARNING`` is defined)
      or a Python function, which should accept one or both parameters
      ``allele`` (existing allele) and ``locus`` (index of locus). The
      return value will become the new allele. This function is
      intended to recode some alleles without listing all alleles in a
      list. It will be called once for each existing allele so it is
      not possible to recode an allele to different alleles. A new
      list of allele names could be specified for these *loci*.
      Different sets of names could be specified for each locus if a
      nested list of names are given. This function recode alleles for
      all subpopulations in all ancestral generations.

   .. method:: Population.removeIndividuals(indexes=[], IDs=[], idField="ind_id", filter=None)

      remove individual(s) by absolute indexes (parameter *index*) or
      their IDs (parameter *IDs*), or using a filter function
      (paramter *filter*). If indexes are used, only individuals at
      the current generation will be removed. If IDs are used, all
      individuals with one of the IDs at information field *idField*
      (default to ``"ind_id"``) will be removed. Although ``"ind_id"``
      usually stores unique IDs of individuals, this function is
      frequently used to remove groups of individuals with the same
      value at an information field. An :class:`IndexError` will be
      raised if an index is out of bound, but no error will be given
      if an invalid ID is specified. In the last case, a user-defined
      function should be provided. This function should accept
      parameter ``"ind"`` or one or more of the information fields.
      All individuals, including ancestors if there are multiple
      ancestral generations, will be passed to this function.
      Individuals that returns ``True`` will be removed. This function
      does not affect subpopulation structure in the sense that a
      subpopulation will be kept even if all individuals from it are
      removed.

   .. method:: Population.removeInfoFields(fields)

      Remove information fields *fields* from a population.

   .. method:: Population.removeLoci(loci=UNSPECIFIED, keep=UNSPECIFIED)

      Remove *loci* (absolute indexes or names) and genotypes at these
      loci from the current population. Alternatively, a parameter
      *keep* can be used to specify loci that will not be removed.

   .. method:: Population.removeSubPops(subPops)

      Remove (virtual) subpopulation(s) *subPops* and all their
      individuals. This function can be used to remove complete
      subpopulations (with shifted subpopulation indexes) or
      individuals belonging to virtual subpopulations of a
      subpopulation. In the latter case, the subpopulations are kept
      even if all individuals have been removed. This function only
      handles the present generation.

   .. method:: Population.resize(sizes, propagate=False)

      Resize population by giving new subpopulation sizes *sizes*.
      individuals at the end of some subpopulations will be removed if
      the new subpopulation size is smaller than the old one. New
      individuals will be appended to a subpopulation if the new size
      is larger. Their genotypes will be set to zero (default), or be
      copied from existing individuals if *propagate* is set to
      ``True``. More specifically, if a subpopulation with ``3``
      individuals is expanded to ``7``, the added individuals will
      copy genotypes from individual ``1``, ``2``, ``3``, and ``1``
      respectively. Note that this function only resizes the current
      generation.

   .. method:: Population.save(filename)

      Save population to a file *filename*, which can be loaded by a
      global function ``loadPopulation(filename)``.

   .. method:: Population.setAncestralDepth(depth)

      set the intended ancestral depth of a population to *depth*,
      which can be ``0`` (does not store any ancestral generation),
      ``-1`` (store all ancestral generations), and a positive number
      (store *depth* ancestral generations. If there exists more than
      *depth* ancestral generations (if *depth* > 0), extra ancestral
      generations are removed.

   .. method:: Population.setGenotype(geno, subPop=[])

      Fill the genotype of all individuals in a population (if
      ``subPop=[]``) or in a (virtual) subpopulation *subPop* (if
      ``subPop=sp`` or ``(sp, vsp)``) using a list of alleles *geno*.
      *geno* will be reused if its length is less than
      ``subPopSize(subPop)*totNumLoci()*ploidy()``.

   .. method:: Population.setIndInfo(values, field, subPop=[])

      Set information field ``field`` (specified by index or name) of
      all individuals (if ``subPop=[]``, default), or individuals in a
      (virtual) subpopulation (``subPop=sp`` or ``(sp, vsp)``) to
      *values*. *values* will be reused if its length is smaller than
      the size of the population or (virtual) subpopulation.

   .. method:: Population.setInfoFields(fields, init=0)

      Set information fields *fields* to a population and initialize
      them with value *init*. All existing information fields will be
      removed.

   .. method:: Population.setLineage(geno, subPop=[])

      Fill the lineage of all individuals in a population (if
      ``subPop=[]``) or in a (virtual) subpopulation *subPop* (if
      ``subPop=sp`` or ``(sp, vsp)``) using a list of IDs *lineage*.
      *lineage* will be reused if its length is less than
      ``subPopSize(subPop)*totNumLoci()*ploidy()``. This function
      returns directly for modules without lineage information.

   .. method:: Population.setSubPopByIndInfo(field)

      Rearrange individuals to their new subpopulations according to
      their integer values at information field *field* (value
      returned by ``Individual::info(field)``). individuals with
      negative values at this *field* will be removed. Existing
      subpopulation names are kept. New subpopulations will have empty
      names.

   .. method:: Population.setSubPopName(name, subPop)

      Assign a name *name* to subpopulation *subPop*. Note that
      subpopulation names do not have to be unique.

   .. method:: Population.setVirtualSplitter(splitter)

      Set a VSP *splitter* to the population, which defines the same
      VSPs for all subpopulations. If different VSPs are needed for
      different subpopulations, a :class:`CombinedSplitter` can be
      used to make these VSPs available to all subpopulations.

   .. method:: Population.sortIndividuals(infoFields, reverse=False)

      Sort individuals according to values at specified information
      fields (*infoFields*). Individuals will be sorted at an
      increasing order unless *reverse* is set to ``true``.

   .. method:: Population.splitSubPop(subPop, sizes, names=[])

      Split subpopulation *subPop* into subpopulations of given
      *sizes*, which should add up to the size of subpopulation
      *subPop* or *1*, in which case *sizes* are treated as
      proportions. If *subPop* is not the last subpopulation, indexes
      of subpopulations after *subPop* are shifted. If *subPop* is
      named, the same name will be given to all new subpopulations
      unless a new set of *names* are specified for these
      subpopulations. This function returns the IDs of split
      subpopulations.

   .. method:: Population.subPopBegin(subPop)

      Return the index of the first individual in subpopulation
      *subPop*.

   .. method:: Population.subPopByName(name)

      Return the index of the first subpopulation with name *name*. An
      :class:`IndexError` will be raised if subpopulations are not
      named, or if no subpopulation with name *name* is found. Virtual
      subpopulation name is not supported.

   .. method:: Population.subPopEnd(subPop)

      Return the index of the last individual in subpopulation
      *subPop* plus ``1``, so that ``range(subPopBegin(subPop)``,
      ``subPopEnd(subPop)`` can iterate through the index of all
      individuals in subpopulation *subPop*.

   .. method:: Population.subPopIndPair(idx)

      Return the subpopulation ID and relative index of an individual,
      given its absolute index ``idx``.

   .. method:: Population.subPopName(subPop)

      Return the "spName - vspName" (virtual named subpopulation), ""
      (unnamed non-virtual subpopulation), "spName" (named
      subpopulation) or "vspName" (unnamed virtual subpopulation),
      depending on whether subpopulation is named or if *subPop* is
      virtual.

   .. method:: Population.subPopNames()

      Return the names of all subpopulations (excluding virtual
      subpopulations). An empty string will be returned for unnamed
      subpopulations.

   .. method:: Population.subPopSizes(ancGen=-1)

      Return the sizes of all subpopulations at the current generation
      (default) or specified ancestral generation *ancGen*. Virtual
      subpopulations are not considered.

   .. method:: Population.swap(rhs)

      Swap the content of two population objects, which can be handy
      in some particular circumstances. For example, you could swap
      out a population in a simulator.

   .. method:: Population.updateInfoFieldsFrom(fields, pop, fromFields=[], ancGens=ALL_AVAIL)

      Update information fields *fields* from *fromFields* of another
      population (or  Pedigree) *pop*. Two populations should have the
      same number of individuals. If *fromFields* is not specified, it
      is assumed to be the same as *fields*. If *ancGens* is not
      ``ALL_AVAIL``, only the specified ancestral generations are
      updated.

   .. method:: Population.useAncestralGen(idx)

      Making ancestral generation *idx* (``0`` for current generation,
      ``1`` for parental generation, ``2`` for grand-parental
      generation, etc) the current generation. This is an efficient
      way to access  Population properties of an ancestral generation.
      ``useAncestralGen(0)`` should always be called afterward to
      restore the correct order of ancestral generations.

   .. method:: Population.vars(subPop=[])

      return variables of a population as a Python dictionary. If a
      valid subpopulation *subPop* is specified, a dictionary
      ``vars()["subPop"][subPop]`` is returned. A :class:`ValueError`
      will be raised if key *subPop* does not exist in ``vars()``, or
      if key *subPop* does not exist in ``vars()["subPop"]``.

   .. method:: Population.virtualSplitter()

      Return the virtual splitter associated with the population,
      ``None`` will be returned if there is no splitter.

   .. method:: Population.asPedigree(idField='ind_id', fatherField='father_id', motherField='mother_id')

      Convert the existing population object to a pedigree. After this
      function pedigree function should magically be usable for this
      function.

   .. method:: Population.subPopSize(subPop=[], ancGen=-1, sex=ANY_SEX)

      Return the size of a subpopulation (``subPopSize(sp)``) or a
      virtual subpopulation (``subPopSize([sp, vsp])``) in the current
      generation (default) or a specified ancestral generation
      *ancGen*. If no *subpop* is given, it is the same as
      ``popSize(ancGen, sex)``.  Population and virtual subpopulation
      names can be used. This function by default returns number of
      all individuals (``sex=ANY_SEX``), but it will return number of
      males (if ``sex=MALE_ONLY``), number of females (if
      ``sex=MALE_ONLY``), and number of male/female pairs (if
      ``sex=PAIR_ONLY``) which is essentially less of the number of
      males and females. <group>2-subpopsize</grouplociList()>

   .. method:: Population.allIndividuals(subPops=ALL_AVAIL, ancGens=True)

      Return an iterator that iterat through all (virtual)
      subpopulations in all ancestral generations. A list of (virtual)
      subpopulations (*subPops*) and a list of ancestral generations
      (*ancGens*, can be a single number) could be specified to
      iterate through only selected subpopulation and generations.
      Value ``ALL_AVAIL`` is acceptable in the specification of ``sp``
      and/or ``vsp`` in specifying a virtual subpopulation ``(sp,
      vsp)`` for the iteration through all or specific virtual
      subpopulation in all or specific subpopulations.

   .. method:: Population.evolve(initOps=[], preOps=[], matingScheme=MatingScheme(), postOps=[], finalOps=[], gen=-1, dryrun=False)

      Evolve the current population *gen* generations using mating
      scheme *matingScheme* and operators *initOps* (applied before
      evolution), *preOps* (applied to the parental population at the
      beginning of each life cycle), *postOps* (applied to the
      offspring population at the end of each life cycle) and
      *finalOps* (applied at the end of evolution). More specifically,
      this function creates a *Simulator* using the current
      population, call its *evolve* function using passed parameters
      and then replace the current population with the evolved
      population. Please refer to function  ``Simulator.evolve`` for
      more details about each parameter.


class Pedigree
--------------

.. class:: Pedigree

   The pedigree class is derived from the population class. Unlike a
   population class that emphasizes on individual properties, the
   pedigree class emphasizes on relationship between individuals. An
   unique ID for all individuals is needed to create a pedigree object
   from a population object. Compared to the :class:`Population`
   class, a :class:`Pedigree` object is optimized for access
   individuals by their IDs, regardless of population structure and
   ancestral generations. Note that the implementation of some
   algorithms rely on the fact that parental IDs are smaller than
   their offspring because individual IDs are assigned sequentially
   during evolution. Pedigrees with manually assigned IDs should try
   to obey such a rule.


   .. method:: Pedigree(pop, loci=[], infoFields=[], ancGens=ALL_AVAIL, idField="ind_id", fatherField="father_id", motherField="mother_id", stealPop=False)


      Create a pedigree object from a population, using a subset of
      loci (parameter *loci*, can be a list of loci indexes, names, or
      ``ALL_AVAIL``, default to no locus), information fields
      (parameter *infoFields*, default to no information field besides
      *idField*, *fatherField* and *motherField*), and ancestral
      generations (parameter *ancGens*, default to all ancestral
      generations). By default, information field ``father_id``
      (parameter *fatherField*) and ``mother_id`` (parameter
      *motherField*) are used to locate parents identified by
      ``ind_id`` (parameter *idField*), which should store an unique
      ID for all individuals. Multiple individuls with the same ID are
      allowed and will be considered as the same individual, but a
      warning will be given if they actually differ in genotype or
      information fields. Operators :class:`IdTagger` and
      :class:`PedigreeTagger` are usually used to assign such IDs,
      although function ``sampling.indexToID`` could be used to assign
      unique IDs and construct parental IDs from index based
      relationship recorded by operator :class:`ParentsTagger`. A
      pedigree object could be constructed with one or no parent but
      certain functions such as relative tracking will not be
      available for such pedigrees. In case that your are no longer
      using your population object, you could steal the content from
      the population by setting *stealPop* to ``True``.


   .. method:: Pedigree.clone()

      Create a cloned copy of a  Pedigree.

   .. method:: Pedigree.identifyAncestors(IDs=ALL_AVAIL, subPops=ALL_AVAIL, ancGens=ALL_AVAIL)

      If a list of individuals (*IDs*) is given, this function traces
      backward in time and find all ancestors of these individuals. If
      *IDs* is ``ALL_AVAIL``, ancestors of all individuals in the
      present generation will be located. If a list of (virtual)
      subpopulations (*subPops*) or ancestral geneartions (*ancGens*)
      is given, the search will be limited to individuals in these
      subpopulations and generations. This could be used to, for
      example, find all fathers of *IDs*. This function returns a list
      of IDs, which includes valid specified IDs. Invalid IDs will be
      silently ignored. Note that parameters *subPops* and *ancGens*
      will limit starting IDs if ``IDs`` is set to ``ALL_AVAIL``, but
      specified IDs will not be trimmed according to these parameters.

   .. method:: Pedigree.identifyFamilies(pedField="", subPops=ALL_AVAIL, ancGens=ALL_AVAIL)

      This function goes through all individuals in a pedigree and
      group related individuals into families. If an information field
      *pedField* is given, indexes of families will be assigned to
      this field of each family member. The return value is a list of
      family sizes corresponding to families 0, 1, 2, ... etc. If a
      list of (virtual) subpopulations (parameter *subPops*) or
      ancestral generations are specified (parameter *ancGens*), the
      search will be limited to individuals in these subpopulations
      and generations.

   .. method:: Pedigree.identifyOffspring(IDs=[], subPops=ALL_AVAIL, ancGens=ALL_AVAIL)

      This function traces forward in time and find all offspring of
      individuals specified in parameter *IDs*. If a list of (virtual)
      subpopulations (*subPops*) or ancestral geneartions (*ancGens*)
      is given, the search will be limited to individuals in these
      subpopulations and generations. This could be used to, for
      example, find all male offspring of *IDs*. This function returns
      a list of IDs, which includes valid starting *IDs*. Invalid IDs
      are silently ignored. Note that parameters *subPops* and
      *ancGens* will limit search result but will not be used to trim
      specified *IDs*.

   .. method:: Pedigree.indByID(id)

      Return a reference to individual with *id*. An
      :class:`IndexError` will be raised if no individual with *id* is
      found. An float *id* is acceptable as long as it rounds closely
      to an integer.

   .. method:: Pedigree.individualsWithRelatives(infoFields, sex=[], affectionStatus=[], subPops=ALL_AVAIL, ancGens=ALL_AVAIL)

      Return a list of IDs of individuals who have non-negative values
      at information fields *infoFields*. Additional requirements
      could be specified by parameters *sex* and *affectionStatus*.
      *sex* can be ``ANY_SEX`` (default), ``MALE_ONLY``,
      ``FEMALE_ONLY``, ``SAME_SEX`` or ``OPPOSITE_SEX``, and
      *affectionStatus* can be ``AFFECTED``, ``UNAFFECTED`` or
      ``ANY_AFFECTION_STATUS`` (default). This function by default
      check all individuals in all ancestral generations, but you
      could limit the search using parameter *subPops* (a list of
      (virtual) subpopulations) and ancestral generations *ancGens*.
      Relatives fall out of specified subpopulations and ancestral
      generaions will be considered invalid.

   .. method:: Pedigree.locateRelatives(relType, resultFields=[], sex=ANY_SEX, affectionStatus=ANY_AFFECTION_STATUS, ancGens=ALL_AVAIL)

      This function locates relatives (of type *relType*) of each
      individual and store their IDs in information fields
      *relFields*. The length of *relFields* determines how many
      relatives an individual can have.

      Parameter *relType* specifies what type of relative to locate,
      which can be

      + ``SPOUSE`` locate spouses with whom an individual has at least
        one common offspring.

      + ``OUTBRED_SPOUSE`` locate non-slibling spouses, namely spouses
        with no shared parent.

      + ``OFFSPRING`` all offspring of each individual.

      + ``COMMON_OFFSPRING`` common offspring between each individual
        and its spouse (located by ``SPOUSE`` or ``OUTBRED_SPOUSE``).
        *relFields* should consist of an information field for spouse
        and ``m-1`` fields for offspring where ``m`` is the number of
        fields.

      + ``FULLSIBLING`` siblings with common father and mother,

      + ``SIBLING`` siblings with at least one common parent.

      Optionally, you can specify the sex and affection status of
      relatives you would like to locate, using parameters *sex* and
      *affectionStatus*. *sex* can be ``ANY_SEX`` (default),
      ``MALE_ONLY``, ``FEMALE_ONLY``, ``SAME_SEX`` or
      ``OPPOSITE_SEX``, and *affectionStatus* can be ``AFFECTED``,
      ``UNAFFECTED`` or ``ANY_AFFECTION_STATUS`` (default). Only
      relatives with specified properties will be located.

      This function will by default go through all ancestral
      generations and locate relatives for all individuals. This can
      be changed by setting parameter *ancGens* to certain ancestral
      generations you would like to process.

   .. method:: Pedigree.save(filename, infoFields=[], loci=[])

      Save a pedigree to file *filename*. This function goes through
      all individuals of a pedigree and outputs in each line the ID of
      individual, IDs of his or her parents, sex (``'M'`` or ``'F'``),
      affection status (``'A'`` or ``'U'``), values of specified
      information fields *infoFields* and genotypes at specified loci
      (parameter ``loci``, which can be a list of loci indexes, names,
      or ``ALL_AVAIL``). Allele numbers, instead of their names are
      outputed. Two columns are used for each locus if the population
      is diploid. This file can be loaded using function
      :func:`loadPedigree` although additional information such as
      names of information fields need to be specified. This format
      differs from a ````.ped file used in some genetic analysis
      software in that there is no family ID and IDs of all
      individuals have to be unique. Note that parental IDs will be
      set to zero if the parent is not in the pedigree object.
      Therefore, the parents of individuals in the top-most ancestral
      generation will always be zero.

   .. method:: Pedigree.traceRelatives(fieldPath, sex=[], affectionStatus=[], resultFields=[], ancGens=ALL_AVAIL)

      Trace a relative path in a population and record the result in
      the given information fields *resultFields*. This function is
      used to locate more distant relatives based on the relatives
      located by function ``locateRelatives``. For example, after
      siblings and offspring of all individuals are located, you can
      locate mother's sibling's offspring using a *relative path*, and
      save their indexes in each individuals information fields
      *resultFields*.

      A *relative path* consits of a *fieldPath* that specifies which
      information fields to look for at each step, a *sex* specifies
      sex choices at each generation, and a *affectionStatus* that
      specifies affection status at each generation. *fieldPath*
      should be a list of information fields, *sex* and
      *affectionStatus* are optional. If specified, they should be a
      list of ``ANY_SEX``, ``MALE_ONLY``, ``FEMALE_ONLY``,
      ``SAME_SEX`` and ``OppsiteSex`` for parameter *sex*, and a list
      of ``UNAFFECTED``, ``AFFECTED`` and ``ANY_AFFECTION_STATUS`` for
      parameter *affectionStatus*.

      For example, if ``fieldPath = [['father_id', 'mother_id'],
      ['sib1', 'sib2'], ['off1', 'off2']]``, and ``sex = [ANY_SEX,
      MALE_ONLY, FEMALE_ONLY]``, this function will locate
      ``father_id`` and ``mother_id`` for each individual, find all
      individuals referred by ``father_id`` and ``mother_id``, find
      informaton fields ``sib1`` and ``sib2`` from these parents and
      locate male individuals referred by these two information
      fields. Finally, the information fields ``off1`` and ``off2``
      from these siblings are located and are used to locate their
      female offspring. The results are father or mother's brother's
      daughters. Their indexes will be saved in each individuals
      information fields *resultFields*. If a list of ancestral
      generations is given in parameter *ancGens* is given, only
      individuals in these ancestral generations will be processed.

   .. method:: Pedigree.asPopulation()

      Convert the existing pedigree object to a population. This
      function will behave like a regular population after this
      function call.


class Simulator
---------------

.. class:: Simulator

   A  simuPOP simulator is responsible for evolving one or more
   populations forward in time, subject to various *operators*.
   Populations in a simulator are created from one or more replicates
   of specified populations. A number of functions are provided to
   access and manipulate populations, and most importantly, to evolve
   them.


   .. method:: Simulator(pops, rep=1, stealPops=True)


      Create a simulator with *rep* (default to ``1``) replicates of
      populations *pops*, which is a list of populations although a
      single population object is also acceptable. Contents of passed
      populations are by default moved to the simulator to avoid
      duplication of potentially large population objects, leaving
      empty populations behind. This behavior can be changed by
      setting *stealPops* to ``False``, in which case populations are
      copied to the simulator.


   .. method:: Simulator.add(pop, stealPop=True)

      Add a population *pop* to the end of an existing simulator. This
      function by default moves *pop* to the simulator, leaving an
      empty population for passed population object. If *steal* is set
      to ``False``, the population will be copied to the simulator,
      and thus unchanged.

   .. method:: Simulator.clone()

      Clone a simulator, along with all its populations. Note that
      Python assign statement ``simu1 = simu`` only creates a symbolic
      link to an existing simulator.

   .. method:: Simulator.__cmp__(rhs)

      a Pyton function used to compare the simulator objects Note that
      mating schemes are not tested.

   .. method:: Simulator.dvars(rep, subPop=[])

      Return a wrapper of Python dictionary returned by ``vars(rep,
      subPop)`` so that dictionary keys can be accessed as attributes.

   .. method:: Simulator.evolve(initOps=[], preOps=[], matingScheme=MatingScheme, postOps=[], finalOps=[], gen=-1, dryrun=False)

      Evolve all populations *gen* generations, subject to several
      lists of operators which are applied at different stages of an
      evolutionary process. Operators *initOps* are applied to all
      populations (subject to applicability restrictions of the
      operators, imposed by the *rep* parameter of these operators)
      before evolution. They are used to initialize populations before
      evolution. Operators *finalOps* are applied to all populations
      after the evolution.

      Operators *preOps*, and *postOps* are applied during the life
      cycle of each generation. These operators can be applied at all
      or some of the generations, to all or some of the evolving
      populations, depending the *begin*, *end*, *step*, *at* and
      *reps* parameters of these operators. These operators are
      applied in the order at which they are specified. populations in
      a simulator are evolved one by one. At each generation,
      operators *preOps* are applied to the parental generations. A
      mating scheme is then used to populate an offspring generation.
      For each offspring, his or her sex is determined before during-
      mating operators of the mating scheme are used to transmit
      parental genotypes. After an offspring generation is
      successfully generated and becomes the current generation,
      operators *postOps* are applied to the offspring generation. If
      any of the *preOps* and *postOps* fails (return ``False``), the
      evolution of a population will be stopped. The generation number
      of a population, which is the variable ``"gen"`` in each
      populations local namespace, is increased by one if an offspring
      generation has been successfully populated even if a post-mating
      operator fails. Another variable ``"rep"`` will also be set to
      indicate the index of each population in the simulator. Note
      that populations in a simulator does not have to have the same
      generation number. You could reset a population's generation
      number by changing this variable.

      Parameter *gen* can be set to a non-negative number, which is
      the number of generations to evolve. If a simulator starts at
      the beginning of a generation ``g`` (for example 0), a simulator
      will stop at the beginning (instead of the end) of generation
      ``g + gen`` (for example gen). If *gen* is negative (default),
      the evolution will continue indefinitely, until all replicates
      are stopped by operators that return ``False`` at some point
      (these operators are called *terminators*). At the end of the
      evolution, the generations that each replicates have evolved are
      returned. Note that *finalOps* are applied to all applicable
      population, including those that have stopped before others.

      If parameter *dryrun* is set to ``True``, this function will
      print a description of the evolutionary process generated by
      function ``describeEvolProcess()`` and exits.

   .. method:: Simulator.extract(rep)

      Extract the *rep-th* population from a simulator. This will
      reduce the number of populations in this simulator by one.

   .. method:: Simulator.numRep()

      Return the number of replicates.

   .. method:: Simulator.population(rep)

      Return a reference to the *rep-th* population of a simulator.
      The reference will become invalid once the simulator starts
      evolving or becomes invalid (removed). If an independent copy of
      the population is needed, you can use ``population.clone()`` to
      create a cloned copy or ``simulator.extract()`` to remove the
      population from the simulator.

   .. method:: Simulator.populations()

      Return a Python iterator that can be used to iterate through all
      populations in a simulator.

   .. method:: Simulator.vars(rep, subPop=[])

      Return the local namespace of the *rep-th* population,
      equivalent to ``x.Population(rep).vars(subPop)``.


