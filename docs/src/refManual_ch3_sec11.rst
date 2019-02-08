Statistics Calculation
======================


class Stat
----------

.. class:: Stat

   Operator :class:`Stat` calculates various statistics of the
   population being applied and sets variables in its local namespace.
   Other operators or functions can retrieve results from or evalulate
   expressions in this local namespace after :class:`Stat` is applied.


   .. method:: Stat(popSize=False, numOfMales=False, numOfAffected=False, numOfSegSites=[], numOfMutants=[], alleleFreq=[], heteroFreq=[], homoFreq=[], genoFreq=[], haploFreq=[], haploHeteroFreq=[], haploHomoFreq=[], sumOfInfo=[], meanOfInfo=[], varOfInfo=[], maxOfInfo=[], minOfInfo=[], LD=[], association=[], neutrality=[], structure=[], HWE=[], inbreeding=[], effectiveSize=[], vars=ALL_AVAIL, suffix="", output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a :class:`Stat` operator that calculates specified
      statistics of a population when it is applied to this
      population. This operator can be applied to specified replicates
      (parameter *rep*) at specified generations (parameter *begin*,
      *end*, *step*, and *at*). This operator does not produce any
      output (ignore parameter *output*) after statistics are
      calculated. Instead, it stores results in the local namespace of
      the population being applied. Other operators can retrieve these
      variables or evalulate expression directly in this local
      namespace. Please refer to operator :class:`BaseOperator` for a
      detailed explanation of these common operator parameters.

      :class:`Stat` supports parameter *subPops*. It usually calculate
      the same set of statistics for all subpopulations
      (``subPops=subPopList()``). If a list of (virtual)
      subpopulations are specified, statistics for only specified
      subpopulations will be calculated. However, different statistics
      treat this parameter differently and it is very important to
      check its reference before you use *subPops* for any statistics.

      Calculated statistics are saved as variables in a population's
      local namespace. These variables can be numbers, lists or
      dictionaries and can be retrieved using functions
      ``Population.vars()`` or ``Population.dvars()``. A special
      default dictionary (:class:`defdict`) is used for dictionaries
      whose keys are determined dynamically. Accessing elements of
      such a dictionary with an invalid key will yield value 0 instead
      of a ``KeyError``. If the same variables are calculated for one
      or more (virtual) subpopulation, the variables are stored in
      ``vars()['subPop'][sp]['var']`` where sp is a subpopulation ID
      (``sp``) or a tuple of virtual subpopulation ID (``(sp, vsp)``).
      ``Population.vars(sp)`` and ``Population.dvars(sp)`` provide
      shortcuts to these variables.

      Operator *Stat* outputs a number of most useful variables for
      each type of statistic. For example, ``alleleFreq`` calculates
      both allele counts and allele frequencies and it by default sets
      variable ``alleleFreq`` (``dvars().alleleFreq``) for all or
      specified subpopulations. If this does not fit your need, you
      can use parameter *vars* to output additional parameters, or
      limit the output of existing parameters. More specifically, for
      this particular statistic, the available variables are
      ``'alleleFreq'``, ``'alleleNum'``, ``'alleleFreq_sp'``
      (``'alleleFreq'`` in each subpopulation), and ``'alleleNum_sp'``
      (``'alleleNum'`` in each subpopulation). You can set
      ``vars=['alleleNum_sp']`` to output only subpopulation specific
      allele count. An optional suffix (parameter *suffix*) can be
      used to append a suffix to default parameter names. This
      parameter can be used, for example, to calculate and store the
      same statistics for different subpopulations (e.g. pairwise
      ``Fst``).

      Operator :class:`Stat` supports the following statistics:

      **popSize**: If *popSize=True*, number of individuals in all or
      specified subpopulations (parameter *subPops*) will be set to
      the following variables:

      + ``popSize`` (default): Number of individuals in all or
        specified subpopulations. Because *subPops* does not have to
        cover all individuals, it may not be the actual population
        size.

      + ``popSize_sp:`` Size of (virtual) subpopulation ``sp``.

      + ``subPopSize`` (default): A list of (virtual) subpopulation
        sizes. This variable is easier to use than accessing popSize
        from each (virtual) subpopulation.

      **numOfMales**: If *numOfMales=True*, number of male individuals
      in all or specified (virtual) subpopulations will be set to the
      following variables:

      + ``numOfMales`` (default): Total number of male individuals in
        all or specified (virtual) subpopulations.

      + ``numOfFemales`` (default): Total number of female individuals
        in all or specified (virtual) subpopulations.

      + ``propOfMales:`` Proportion of male individuals.

      + ``propOfFemales:`` Proportion of female individuals.

      + ``numOfMales_sp:`` Number of male individuals in each
        (virtual) subpopulation.

      + ``numOfFemales_sp:`` Number of female individuals in each
        (virtual) subpopulation.

      + ``propOfMales_sp:`` Proportion of male individuals in each
        (virtual) subpopulation.

      + ``propOfFemales_sp:`` Proportion of female individuals in each
        (virtual) subpopulation.

      **numOfAffected**: If *numOfAffected=True*, number of affected
      individuals in all or specified (virtual) subpopulations will be
      set to the following variables:

      + ``numOfAffected`` (default): Total number of affected
        individuals in all or specified (virtual) subpopulations.

      + ``numOfUnaffected`` (default): Total number of unaffected
        individuals in all or specified (virtual) subpopulations.

      + ``propOfAffected:`` Proportion of affected individuals.

      + ``propOfUnaffected:`` Proportion of unaffected individuals.

      + ``numOfAffected_sp:`` Number of affected individuals in each
        (virtual) subpopulation.

      + ``numOfUnaffected_sp:`` Number of unaffected individuals in
        each (virtual) subpopulation.

      + ``propOfAffected_sp:`` Proportion of affected individuals in
        each (virtual) subpopulation.

      + ``propOfUnaffected_sp:`` Proportion of unaffected individuals
        in each (virtual) subpopulation.

      **numOfSegSites**: Parameter *numOfSegSites* accepts a list of
      loci (loci indexes, names, or ``ALL_AVAIL``) and count the
      number of loci with at least two different alleles (segregating
      sites) or loci with only one non-zero allele (no zero allele,
      not segragating) for individuals in all or specified (virtual)
      subpopulations. This parameter sets variables

      + ``numOfSegSites`` (default): Number of segregating sites in
        all or specified (virtual) subpopulations.

      + ``numOfSegSites_sp:`` Number of segregating sites in each
        (virtual) subpopulation.

      + ``numOfFixedSites:`` Number of sites with one non-zero allele
        in all or specified (virtual) subpopulations.

      + ``numOfFixedSites_sp:`` Number of sites with one non-zero
        allele in in each (virtual) subpopulations.

      + ``segSites:`` A list of segregating sites in all or specified
        (virtual) subpopulations.

      + ``segSites_sp:`` A list of segregating sites in each (virtual)
        subpopulation.

      + ``fixedSites:`` A list of sites with one non-zero allele in
        all or specified (virtual) subpopulations.

      + ``fixedSites_sp:`` A list of sites with one non-zero allele in
        in each (virtual) subpopulations.

      **numOfMutants**: Parameter *numOfMutants* accepts a list of
      loci (loci indexes, names, or ``ALL_AVAIL``) and count the
      number of mutants (non-zero alleles) for individuals in all or
      specified (virtual) subpopulations. It sets variables

      + ``numOfMutants`` (default): Number of mutants in all or
        specified (virtual) subpopulations.

      + ``numOfMutants_sp:`` Number of mutants in each (virtual)
        subpopulations.

      **alleleFreq**: This parameter accepts a list of loci (loci
      indexes, names, or ``ALL_AVAIL``), at which allele frequencies
      will be calculated. This statistic outputs the following
      variables, all of which are dictionary (with loci indexes as
      keys) of default dictionaries (with alleles as keys). For
      example, ``alleleFreq[loc][a]`` returns 0 if allele ``a`` does
      not exist.

      + ``alleleFreq`` (default): ``alleleFreq[loc][a]`` is the
        frequency of allele ``a`` at locus \loc for all or specified
        (virtual) subpopulations.

      + ``alleleNum`` (default): ``alleleNum[loc][a]`` is the number
        of allele ``a`` at locus \loc for all or specified (virtual)
        subpopulations.

      + ``alleleFreq_sp:`` Allele frequency in each (virtual)
        subpopulation.

      + ``alleleNum_sp:`` Allele count in each (virtual)
        subpopulation.

      **heteroFreq** and **homoFreq**: These parameters accept a list
      of loci (by indexes or names), at which the number and frequency
      of homozygotes and/or heterozygotes will be calculated. These
      statistics are only available for diploid populations. The
      following variables will be outputted:

      + ``heteroFreq`` (default for parameter *heteroFreq*): A
        dictionary of proportion of heterozygotes in all or specified
        (virtual) subpopulations, with loci indexes as dictionary
        keys.

      + ``homoFreq`` (default for parameter *homoFreq*): A dictionary
        of proportion of homozygotes in all or specified (virtual)
        subpopulations.

      + ``heteroNum:`` A dictionary of number of heterozygotes in all
        or specified (virtual) subpopulations.

      + ``homoNum:`` A dictionary of number of homozygotes in all or
        specified (virtual) subpopulations.

      + ``heteroFreq_sp:`` A dictionary of proportion of heterozygotes
        in each (virtual) subpopulation.

      + ``homoFreq_sp:`` A dictionary of proportion of homozygotes in
        each (virtual) subpopulation.

      + ``heteroNum_sp:`` A dictionary of number of heterozygotes in
        each (virtual) subpopulation.

      + ``homoNum_sp:`` A dictionary of number of homozygotes in each
        (virtual) subpopulation.

      **genoFreq**: This parameter accept a list of loci (by indexes
      or names) at which number and frequency of all genotypes are
      outputed as a dictionary (indexed by loci indexes) of default
      dictionaries (indexed by tuples of possible indexes). This
      statistic is available for all population types with genotype
      defined as ordered alleles at a locus. The length of genotype
      equals the number of homologous copies of chromosomes (ploidy)
      of a population. Genotypes for males or females on sex
      chromosomes or in haplodiploid populations will have different
      length. Because genotypes are ordered, ``(1, 0)`` and ``(0, 1)``
      (two possible genotypes in a diploid population) are considered
      as different genotypes. This statistic outputs the following
      variables:

      + ``genoFreq`` (default): A dictionary (by loci indexes) of
        default dictionaries (by genotype) of genotype frequencies.
        For example, ``genoFreq[1][(1, 0)]`` is the frequency of
        genotype (1, 0) at locus 1.

      + ``genoNum`` (default): A dictionary of default dictionaries of
        genotype counts of all or specified (virtual) subpopulations.

      + ``genoFreq_sp:`` genotype frequency in each specified
        (virtual) subpopulation.

      + ``genoFreq_sp:`` genotype count in each specified (virtual)
        subpopulation.

      **haploFreq**: This parameter accepts one or more lists of loci
      (by index) at which number and frequency of haplotypes are
      outputted as default dictionaries. ``[(1,2)]`` can be
      abbreviated to ``(1,2)``. For example, using parameter
      ``haploFreq=(1,2,4)``, all haplotypes at loci ``1``, ``2`` and
      ``4`` are counted. This statistic saves results to dictionary
      (with loci index as keys) of default dictionaries (with
      haplotypes as keys) such as ``haploFreq[(1,2,4)][(1,1,0)]``
      (frequency of haplotype ``(1,1,0)`` at loci ``(1,2,3)``). This
      statistic works for all population types. Number of haplotypes
      for each individual equals to his/her ploidy number.
      Haplodiploid populations are supported in the sense that the
      second homologous copy of the haplotype is not counted for male
      individuals. This statistic outputs the following variables:

      + ``haploFreq`` (default): A dictionary (with tuples of loci
        indexes as keys) of default dictionaries of haplotype
        frequencies. For example, ``haploFreq[(0, 1)][(1,1)]`` records
        the frequency of haplotype ``(1,1)`` at loci ``(0, 1)`` in all
        or specified (virtual) subpopulations.

      + ``haploNum`` (default): A dictionary of default dictionaries
        of haplotype counts in all or specified (virtual)
        subpopulations.

      + ``haploFreq_sp:`` Halptype frequencies in each (virtual)
        subpopulation.

      + ``haploNum_sp:`` Halptype count in each (virtual)
        subpopulation.

      **haploHeteroFreq** and **haploHomoFreq**: These parameters
      accept a list of haplotypes (list of loci), at which the number
      and frequency of haplotype homozygotes and/or heterozygotes will
      be calculated. Note that these statistics are **observed** count
      of haplotype heterozygote. The following variables will be
      outputted:

      + ``haploHeteroFreq`` (default for parameter *haploHeteroFreq*):
        A dictionary of proportion of haplotype heterozygotes in all
        or specified (virtual) subpopulations, with haplotype indexes
        as dictionary keys.

      + ``haploHomoFreq`` (default for parameter *haploHomoFreq*): A
        dictionary of proportion of homozygotes in all or specified
        (virtual) subpopulations.

      + ``haploHeteroNum:`` A dictionary of number of heterozygotes in
        all or specified (virtual) subpopulations.

      + ``haploHomoNum:`` A dictionary of number of homozygotes in all
        or specified (virtual) subpopulations.

      + ``haploHeteroFreq_sp:`` A dictionary of proportion of
        heterozygotes in each (virtual) subpopulation.

      + ``haploHomoFreq_sp:`` A dictionary of proportion of
        homozygotes in each (virtual) subpopulation.

      + ``haploHeteroNum_sp:`` A dictionary of number of heterozygotes
        in each (virtual) subpopulation.

      + ``haploHomoNum_sp:`` A dictionary of number of homozygotes in
        each (virtual) subpopulation.

      **sumOfinfo**, **meanOfInfo**, **varOfInfo**, **maxOfInfo** and
      **minOfInfo**: Each of these five parameters accepts a list of
      information fields. For each information field, the sum, mean,
      variance, maximum or minimal (depending on the specified
      parameter(s)) of this information field at iddividuals in all or
      specified (virtual) subpopulations will be calculated. The
      results will be put into the following population variables:

      + ``sumOfInfo`` (default for *sumOfInfo*): A dictionary of the
        sum of specified information fields of individuals in all or
        specified (virtual) subpopulations. This dictionary is indexed
        by names of information fields.

      + ``meanOfInfo`` (default for *meanOfInfo*): A dictionary of the
        mean of information fields of all individuals.

      + ``varOfInfo`` (default for *varOfInfo*): A dictionary of the
        sample variance of information fields of all individuals.

      + ``maxOfInfo`` (default for *maxOfInfo*): A dictionary of the
        maximum value of information fields of all individuals.

      + ``minOfInfo`` (default for *minOfInfo*): A dictionary of the
        minimal value of information fields of all individuals.

      + ``sumOfInfo_sp:`` A dictionary of the sum of information
        fields of individuals in each subpopulation.

      + ``meanOfInfo_sp:`` A dictionary of the mean of information
        fields of individuals in each subpopulation.

      + ``varOfInfo_sp:`` A dictionary of the sample variance of
        information fields of individuals in each subpopulation.

      + ``maxOfInfo_sp:`` A dictionary of the maximum value of
        information fields of individuals in each subpopulation.

      + ``minOfInfo_sp:`` A dictionary of the minimal value of
        information fields of individuals in each subpopulation.

      **LD**: Parameter ``LD`` accepts one or a list of loci pairs
      (e.g. ``LD=[[0,1], [2,3]]``) with optional primary alleles at
      both loci (e.g. ``LD=[0,1,0,0]``). For each pair of loci, this
      operator calculates linkage disequilibrium and optional
      association statistics between two loci. When primary alleles
      are specified, signed linkage disequilibrium values are
      calculated with non-primary alleles are combined. Otherwise,
      absolute values of diallelic measures are combined to yield
      positive measure of LD. Association measures are calculated from
      a ``m`` by ``n`` contigency of haplotype counts (``m=n=2`` if
      primary alleles are specified). Please refer to the  simuPOP
      user's guide for detailed information. This statistic sets the
      following variables:

      + ``LD`` (default) Basic LD measure for haplotypes in all or
        specified (virtual) subpopulations. Signed if primary alleles
        are specified.

      + ``LD_prime`` (default) Lewontin's D' measure for haplotypes in
        all or specified (virtual) subpopulations. Signed if primary
        alleles are specified.

      + ``R2`` (default) Correlation LD measure for haplotypes in all
        or specified (virtual) subpopulations.

      + ``LD_ChiSq`` ChiSq statistics for a contigency table with
        frequencies of haplotypes in all or specified (virtual)
        subpopulations.

      + ``LD_ChiSq_p`` Single side p-value for the ChiSq statistic.
        Degrees of freedom is determined by number of alleles at both
        loci and the specification of primary alleles.

      + ``CramerV`` Normalized ChiSq statistics.

      + ``LD_sp`` Basic LD measure for haplotypes in each (virtual)
        subpopulation.

      + ``LD_prime_sp`` Lewontin's D' measure for haplotypes in each
        (virtual) subpopulation.

      + ``R2_sp`` R2 measure for haplotypes in each (virtual)
        subpopulation.

      + ``LD_ChiSq_sp`` ChiSq statistics for each (virtual)
        subpopulation.

      + ``LD_ChiSq_p_sp`` p value for the ChiSq statistics for each
        (virtual) subpopulation.

      + ``CramerV_sp`` Cramer V statistics for each (virtual)
        subpopulation.

      **association**: Parameter ``association`` accepts a list of
      loci, which can be a list of indexes, names, or ``ALL_AVAIL``.
      At each locus, one or more statistical tests will be performed
      to test association between this locus and individual affection
      status. Currently,  simuPOP provides the following tests:

      + An allele-based Chi-square test using alleles counts. This
        test can be applied to loci with more than two alleles, and to
        haploid populations.

      + A genotype-based Chi-square test using genotype counts. This
        test can be applied to loci with more than two alleles (more
        than 3 genotypes) in diploid populations. ``aA`` and ``Aa``
        are considered to be the same genotype.

      + A genotype-based Cochran-Armitage trend test. This test can
        only be applied to diallelic loci in diploid populations. A
        codominant model is assumed.

      This statistic sets the following variables:

      + ``Allele_ChiSq`` A dictionary of allele-based Chi-Square
        statistics for each locus, using cases and controls in all or
        specified (virtual) subpopulations.

      + ``Allele_ChiSq_p`` (default) A dictionary of *p-values* of the
        corresponding Chi-square statistics.

      + ``Geno_ChiSq`` A dictionary of genotype-based Chi-Square
        statistics for each locus, using cases and controls in all or
        specified (virtual) subpopulations.

      + ``Geno_ChiSq_p`` A dictionary of *p-values* of the
        corresponding genotype-based Chi-square test.

      + ``Armitage_p`` A dictionary of *p-values* of the Cochran-
        Armitage tests, using cases and controls in all or specified
        (virtual) subpopulations.

      + ``Allele_ChiSq_sp`` A dictionary of allele-based Chi-Square
        statistics for each locus, using cases and controls from each
        subpopulation.

      + ``Allele_ChiSq_p_sp`` A dictionary of p-values of allele-based
        Chi-square tests, using cases and controls from each (virtual)
        subpopulation.

      + ``Geno_ChiSq_sp`` A dictionary of genotype-based Chi-Square
        tests for each locus, using cases and controls from each
        subpopulation.

      + ``Geno_ChiSq_p_sp`` A dictionary of p-values of genotype-based
        Chi-Square tests, using cases and controls from each
        subpopulation.

      + ``Armitage_p_sp`` A dictionary of *p-values* of the Cochran-
        Armitage tests, using cases and controls from each
        subpopulation.

      **neutrality**: This parameter performs neutrality tests
      (detection of natural selection) on specified loci, which can be
      a list of loci indexes, names or ``ALL_AVAIL``. It currently
      only outputs *Pi*, which is the average number of pairwise
      difference between loci. This statistic outputs the following
      variables:

      + ``Pi`` Mean pairwise difference between all sequences from all
        or specified (virtual) subpopulations.

      + ``Pi_sp`` Mean paiewise difference between all sequences in
        each (virtual) subpopulation.

      **structure**: Parameter ``structure`` accepts a list of loci at
      which statistics that measure population structure are
      calculated. *structure* accepts a list of loci indexes, names or
      ``ALL_AVAIL``. This parameter currently supports the following
      statistics:

      + Weir and Cockerham's Fst (1984). This is the most widely used
        estimator of Wright's fixation index and can be used to
        measure  Population differentiation. However, this method is
        designed to estimate Fst from samples of larger populations
        and might not be appropriate for the calculation of Fst of
        large populations.

      + Nei's Gst (1973). The Gst estimator is another estimator for
        Wright's fixation index but it is extended for multi-allele
        (more than two alleles) and multi-loci cases. This statistics
        should be used if you would like to obtain a *true* Fst value
        of a large  Population. Nei's Gst uses only allele frequency
        information so it is available for all population type
        (haploid, diploid etc). Weir and Cockerham's Fst uses
        heterozygosity frequency so it is best for autosome of diploid
        populations. For non-diploid population, sex, and
        mitochondrial DNAs,  simuPOP uses expected heterozygosity (1 -
        sum p_i^2) when heterozygosity is needed. These statistics
        output the following variables:

      + ``F_st`` (default) The WC84 *Fst* statistic estimated for all
        * specified loci.

      + ``F_is`` The WC84 *Fis* statistic estimated for all specified
        loci.

      + ``F_it`` The WC84 *Fit* statistic estimated for all specified
        loci.

      + ``f_st`` A dictionary of locus level WC84 *Fst* values.

      + ``f_is`` A dictionary of locus level WC84 *Fis* values.

      + ``f_it`` A dictionary of locus level WC84 *Fit* values.

      + ``G_st`` Nei's Gst statistic estimated for all specified loci.

      + ``g_st`` A dictionary of Nei's Gst statistic estimated for
        each locus.

      **HWE**: Parameter ``HWE`` accepts a list of loci at which exact
      two-side tests for Hardy-Weinberg equilibrium will be performed.
      This statistic is only available for diallelic loci in diploid
      populations. *HWE* can be a list of loci indexes, names or
      ``ALL_AVAIL``. This statistic outputs the following variables:

      + ``HWE`` (default) A dictionary of p-values of HWE tests using
        genotypes in all or specified (virtual) subpopulations.

      + ``HWE_sp`` A dictionary of p-values of HWS tests using
        genotypes in each (virtual) subpopulation.

      **inbreeding**: Inbreeding measured by Identitcal by Decent (and
      by State). This statistics go through all loci of individuals in
      a diploid population and calculate the number and proportions of
      alleles that are identitcal by decent and by state. Because
      ancestral information is only available in lineage module,
      variables IBD_freq are always set to zero in other modules. Loci
      on sex and mitochondrial chromosomes, and non-diploid
      populations are currently not supported. This statistic outputs
      the following variables:

      + ``IBD_freq`` (default) The frequency of IBD pairs among all
        allele pairs. To use this statistic, the population must be
        initialized by operator InitLineage() to assign each ancestral
        allele an unique identify.

      + ``IBS_freq`` (default) The proportion of IBS pairs among all
        allele pairs.

      + ``IBD_freq_sp`` frequency of IBD in each (virtual)
        subpopulations.

      + ``IBS_freq_sp`` frequency of IBS in each (virtual)
        subpopulations.

      **effectiveSize**: Parameter ``effectiveSize`` accepts a list of
      loci at which the effective population size for the whole or
      specified (virtual) subpopulations is calculated.
      *effectiveSize* can be a list of loci indexes, names or
      ``ALL_AVAIL``. Parameter *subPops* is usually used to define
      samples from which effective sizes are estimated. This statistic
      allows the calculation of true effective size based on number of
      gametes each parents transmit to the offspring population (per-
      locus before and after mating), and estimated effective size
      based on sample genotypes. Due to the temporal natural of some
      methods, more than one  Stat operators might be needed to
      calculate effective size. The *vars* parameter specified which
      method to use and which variable to set. Acceptable values
      include:

      + ``Ne_demo_base`` When this variable is set before mating, it
        stores IDs of breeding parents and, more importantly, assign
        an unique lineage value to alleles at specified loci of each
        individual. **This feature is only available for lineage
        modules and will change lineage values at specified loci of
        all individuals**.

      + ``Ne_demo_base_sp`` Pre-mating information for each (virtual)
        subpopulation, used by variable ``Ne_demo_sp``.

      + ``Ne_demo`` A dictionary of locus-specific demographic
        effective population size, calculated using number of gemetes
        each parent transmits to the offspring population. The method
        is vased on Crow & Denniston 1988 (Ne = KN-1/k-1+Vk/k) and
        need variable ``Ne_demo_base`` set before mating. **Effective
        size estimated from this formula is model dependent and might
        not be applicable to your mating schemes.**

      + ``Ne_demo_sp`` Calculate subpopulation-specific effective
        size.

      + ``Ne_temporal_base`` When this variable is set in parameter
        *vars*, the  Stat operator saves baseline allele frequencies
        and other information in this variable, which are used by
        temporary methods to estimate effective population size
        according to changes in allele frequency between the baseline
        and present generations. This variable could be set repeatedly
        to change baselines.

      + ``Ne_temporal_base_sp`` Set baseline information for each
        (virtual) subpopulation specified.

      + ``Ne_tempoFS_P1`` Effective population size, 2.5% and 97.5%
        confidence interval for sampling plan 1 as a list of size 3,
        estimated using a temporal method as described in Jorde &
        Ryman (2007), and as implemented by software tempoFS
        (http://www.zoologi.su.se/~ryman/). This variable is set to
        census population size if no baseline has been set, and to the
        temporal effective size between the present and the baseline
        generation otherwise. This method uses population size or sum
        of subpopulation sizes of specified (virtual) subpopulations
        as census population size for the calculation based on plan 1.

      + ``Ne_tempoFS_P2`` Effective population size, 2.5% and 97.5%
        confidence interval for sampling plan 2 as a list of size 6,
        estimated using a temporal method as described in Jorde &
        Ryman (2007). This variable is set to census population size
        no baseline has been set, and to the temporal effective size
        between the present and the baseline generation otherwise.
        This method assumes that the sample is drawn from an
        infinitely-sized population.

      + ``Ne_tempoFS`` deprecated, use ``Ne_tempoFS_P2`` instead.

      + ``Ne_tempoFS_P1_sp`` Estimate effective size of each (virtual)
        subpopulation using method Jorde & Ryman 2007, assuming
        sampling plan 1. The census population sizes for sampling plan
        1 are the sizes for each subpopulation that contain the
        specified (virtual) subpopulations.

      + ``Ne_tempoFS_P2_sp`` Estimate effective size of each (virtual)
        subpopulation using method Jorde & Ryman 2007, assuming
        sampling plan 2.

      + ``Ne_tempoFS_sp`` deprecated, use ``Ne_tempoFS_P2_sp``
        instead.

      + ``Ne_waples89_P1`` Effective population size, 2.5% and 97.5%
        confidence interval for sampling plan 1 as a list of size 6,
        estimated using a temporal method as described in Waples 1989,
        Genetics. Because this is a temporal method, Ne_waples89
        estimates effective size between the present and the baseline
        generation set by variable ``Ne_temporal_base``. Census
        population size will be resutned if no baseline has been set.
        This method uses population size or sum of subpopulation sizes
        of specified (virtual) subpopulations as census population
        size for the calculation based on plan 1.

      + ``Ne_waples89_P2`` Effective population size, 2.5% and 97.5%
        confidence interval for sampling plan 2 as a list of size 6,
        estimated using a temporal method as described in Waples 1989,
        Genetics. Because this is a temporal method, Ne_waples89
        estimates effective size between the present and the baseline
        generation set by variable ``Ne_temporal_base``. Census
        population size will be returned if no baseline has been set.

      + ``Ne_waples89_P1_sp`` Estimate effective size for each
        (virtual) subpopulation using method Waples 89, assuming
        sampling plan 1. The census population sizes are the sizes for
        each subpopulation that contain the specified (virtual)
        subpopulation.

      + ``Ne_waples89_P2_sp`` Estimate effective size for each
        (virtual) subpopulation using method Waples 89, assuming
        sampling plan 2.

      + ``Ne_waples89_sp`` deprecated, use ``Ne_waples89_P2_sp``
        instead.

      + ``Ne_LD`` Lists of length three for effective population size,
        2.5% and 97.% confidence interval for cutoff allele frequency
        0., 0.01, 0.02 and 0.05 (as dictionary keys), using a
        parametric method, estimated from linkage disequilibrim
        information of one sample, using LD method developed by Waples
        & Do 2006 (LDNe). This method assumes unlinked loci and uses
        LD measured from genotypes at loci. Because this is a sample
        based method, it should better be applied to a random sample
        of the population. 95% CI is calculated using a Jackknife
        estimated effective number of independent alleles. Please
        refer to relevant papers and the LDNe user's guide for
        details.

      + ``Ne_LD_sp`` Estimate LD-based effective population size for
        each specified (virtual) subpopulation.

      + ``Ne_LD_mono`` A version of Ne_LD that assumes monogamy (see
        Waples 2006 for details.

      + ``Ne_LD_mono_sp`` Ne_LD_mono calculated for each (virtual)
        subpopulation.



