Genotype transmitters
=====================


class GenoTransmitter
---------------------

.. class:: GenoTransmitter

   This during mating operator is the base class of all genotype
   transmitters. It is made available to users because it provides a
   few member functions that can be used by derived transmitters, and
   by customized Python during mating operators.


   .. method:: GenoTransmitter(output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a base genotype transmitter.


   .. method:: GenoTransmitter.clearChromosome(ind, ploidy, chrom)

      Clear (set alleles to zero) chromosome *chrom* on the *ploidy-
      th* homologous set of chromosomes of individual *ind*. It is
      equivalent to ``ind.setGenotype([0], ploidy, chrom)``, except
      that it also clears allele lineage if it is executed in a module
      with lineage allele type.

   .. method:: GenoTransmitter.copyChromosome(parent, parPloidy, offspring, ploidy, chrom)

      Transmit chromosome *chrom* on the *parPloidy* set of homologous
      chromosomes from *parent* to the *ploidy* set of homologous
      chromosomes of *offspring*. It is equivalent to
      ``offspring.setGenotype(parent.genotype(parPloidy, chrom),
      polidy, chrom)``, except that it also copies allelic lineage
      when it is executed in a module with lineage allele type.

   .. method:: GenoTransmitter.copyChromosomes(parent, parPloidy, offspring, ploidy)

      Transmit the *parPloidy* set of homologous chromosomes from
      *parent* to the *ploidy* set of homologous chromosomes of
      *offspring*. Customized chromosomes are not copied. It is
      equivalent to
      ``offspring.setGenotype(parent.genotype(parPloidy), ploidy)``,
      except that it also copies allelic lineage when it is executed
      in a module with lineage allele type.


class CloneGenoTransmitter
--------------------------

.. class:: CloneGenoTransmitter

   This during mating operator copies parental genotype directly to
   offspring. This operator works for all mating schemes when one or
   two parents are involved. If both parents are passed, maternal
   genotype are copied. In addition to genotypes on all non-customized
   or specified chromosomes, sex and information fields are by default
   also coped copied from parent to offspring.


   .. method:: CloneGenoTransmitter(output="", chroms=ALL_AVAIL, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)


      Create a clone genotype transmitter (a during-mating operator)
      that copies genotypes from parents to offspring. If two parents
      are specified, genotypes are copied maternally. After genotype
      transmission, offspring sex and affection status is copied from
      the parent even if sex has been determined by an offspring
      generator. All or specified information fields (parameter
      *infoFields*, default to ``ALL_AVAIL``) will also be copied from
      parent to offspring. Parameters *subPops* is ignored. This
      operator by default copies genotypes on all autosome and sex
      chromosomes (excluding customized chromosomes), unless a
      parameter *chroms* is used to specify which chromosomes to copy.
      This operator also copies allelic lineage when it is executed in
      a module with lineage allele type.



class MendelianGenoTransmitter
------------------------------

.. class:: MendelianGenoTransmitter

   This Mendelian offspring generator accepts two parents and pass
   their genotypes to an offspring following Mendel's laws. Sex
   chromosomes are handled according to the sex of the offspring,
   which is usually determined in advance by an offspring generator.
   Customized chromosomes are not handled.


   .. method:: MendelianGenoTransmitter(output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a Mendelian genotype transmitter (a during-mating
      operator) that transmits genotypes from parents to offspring
      following Mendel's laws. Autosomes and sex chromosomes are
      handled but customized chromosomes are ignored. Parameters
      *subPops* and *infoFields* are ignored. This operator also
      copies allelic lineage when it is executed in a module with
      lineage allele type.


   .. method:: MendelianGenoTransmitter.transmitGenotype(parent, offspring, ploidy)

      Transmit genotype from parent to offspring, and fill the
      *ploidy* homologous set of chromosomes. This function does not
      set genotypes of customized chromosomes and handles sex
      chromosomes properly, according to offspring sex and ``ploidy``.


class SelfingGenoTransmitter
----------------------------

.. class:: SelfingGenoTransmitter

   A genotype transmitter (during-mating operator) that transmits
   parental genotype of a parent through self-fertilization. That is
   to say, the offspring genotype is formed according to Mendel's
   laws, only that a parent serves as both maternal and paternal
   parents.


   .. method:: SelfingGenoTransmitter(output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a self-fertilization genotype transmitter that transmits
      genotypes of a parent to an offspring through self-
      fertilization. Cutsomized chromosomes are not handled.
      Parameters *subPops* and *infoFields* are ignored. This operator
      also copies allelic lineage when it is executed in a module with
      lineage allele type.



class HaplodiploidGenoTransmitter
---------------------------------

.. class:: HaplodiploidGenoTransmitter

   A genotype transmitter (during-mating operator) for haplodiploid
   populations. The female parent is considered as diploid and the
   male parent is considered as haploid (only the first homologous
   copy is valid). If the offspring is ``FEMALE``, she will get a
   random copy of two homologous chromosomes of her mother, and get
   the only paternal copy from her father. If the offspring is
   ``MALE``, he will only get a set of chromosomes from his mother.


   .. method:: HaplodiploidGenoTransmitter(output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a haplodiploid genotype transmitter (during-mating
      operator) that transmit parental genotypes from parents to
      offspring in a haplodiploid population. Parameters *subPops* and
      *infoFields* are ignored. This operator also copies allelic
      lineage when it is executed in a module with lineage allele
      type.



class MitochondrialGenoTransmitter
----------------------------------

.. class:: MitochondrialGenoTransmitter

   This geno transmitter transmits the first homologous copy of a
   ``Mitochondrial`` chromosome. If no mitochondrial chromosome is
   present, it assumes that the first homologous copy of several (or
   all) ``Customized`` chromosomes are copies of mitochondrial
   chromosomes. This operator transmits the mitochondrial chromosome
   from the female parent to offspring for sexsual reproduction, and
   any parent to offspring for asexual reproduction. If there are
   multiple chromosomes, the organelles are selected randomly. If this
   transmitter is applied to populations with more than one homologous
   copies of chromosomes, it transmits the first homologous copy of
   chromosomes and clears alleles (set to zero) on other homologous
   copies.


   .. method:: MitochondrialGenoTransmitter(output="", chroms=ALL_AVAIL, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Createa a mitochondrial genotype transmitter that treats the
      Mitochondiral chromosome, or Customized chromosomes if no
      Mitochondrial chromosome is specified, or a list of chromosomes
      specified by *chroms*, as human mitochondrial chromosomes. These
      chromosomes should have the same length and the same number of
      loci. This operator transmits these chromosomes randomly from
      the female parent to offspring of both sexes. It also copies
      allelic lineage when it is executed in a module with lineage
      allele type.



class Recombinator
------------------

.. class:: Recombinator

   A genotype transmitter (during-mating operator) that transmits
   parental chromosomes to offspring, subject to recombination and
   gene conversion. This can be used to replace
   :class:`MendelianGenoTransmitter` and
   :class:`SelfingGenoTransmitter`. It does not work in haplodiploid
   populations, although a customized genotype transmitter that makes
   uses this operator could be defined. Please refer to the  simuPOP
   user's guide or online cookbook for details.

   Recombination could be applied to all adjacent markers or after
   specified loci. Recombination rate between two adjacent markers
   could be specified directly, or calculated using physical distance
   between them. In the latter case, a recombination intensity is
   multiplied by physical distance between markers.

   Gene conversion is interpreted as double-recombination events. That
   is to say, if a recombination event happens, it has a certain
   probability (can be 1) to become a conversion event, namely
   triggering another recombination event down the chromosome. The
   length of the converted chromosome can be controlled in a number of
   ways.

   .. note::

      simuPOP does not assume any unit to loci positions so
      recombination intensity could be explained differntly (e.g.
      cM/Mb, Morgan/Mb) depending on your intepretation of loci
      positions. For example, if basepair is used for loci position,
      ``intensity=10^-8`` indicates ``10^-8`` per basepair, which is
      equivalent to ``10^-2`` per Mb or 1 cM/Mb. If ``Mb`` is used for
      physical positions, the same recombination intensity could be
      achieved by ``intensity=0.01``.


   .. method:: Recombinator(rates=[], intensity=-1, loci=ALL_AVAIL, convMode=NO_CONVERSION, output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a  Recombinator (a mendelian genotype transmitter with
      recombination and gene conversion) that passes genotypes from
      parents (or a parent in case of self-fertilization) to
      offspring.

      Recombination happens by default between all adjacent markers
      but can be limited to a given set of *loci*, which can be a list
      of loci indexes, names, list of chromosome position pairs,
      ``ALL_AVAIL``, or a function with optional parameter ``pop``
      that will be called at each ganeeration to determine indexes of
      loci. Each locus in this list specifies a recombination point
      between the locus and the locus immediately **after** it. Loci
      that are the last locus on each chromosome are ignored.

      If a single recombination rate (parameter *rates*) is specified,
      it will used for all loci (all loci or loci specified by
      parameter *loci*), regardless of physical distances between
      adjacent loci.

      If a list of recombination rates are specified in *rates*,
      different recombination rates could be applied after a list of
      specified loci (between loci and their immediate neighbor to the
      right). The loci should be specified by parameter *loci* as a
      list with the same length as *rates*, or *ALL_AVAIL* (default)
      in which case the length of *rates* should equal to the total
      number of loci. Note that recombination rates specified for the
      last locus on each chromosome are ignored because  simuPOP
      assumes free recombination between chromosomes.

      A recombination intensity (*intensity*) can be used to specify
      recombination rates that are proportional to physical distances
      between adjacent markers. If the physical distance between two
      markers is ``d``, the recombination rate between them will be
      ``intensity * d``. No unit is assume for loci position and
      recombination intensity.

      Gene conversion is controlled using parameter *convMode*, which
      can be

      + ``NoConversion``: no gene conversion (default).

      + ``(NUM_MARKERS, prob, n)``: With probability *prob*, convert a
        fixed number (*n*) of markers if a recombination event
        happens.

      + ``(GEOMETRIC_DISTRIBUTION, prob, p)``: With probability
        *prob*, convert a random number of markers if a recombination
        event happens. The number of markes converted follows a
        geometric distribution with probability *p*.

      + ``(TRACT_LENGTH, prob, n)``: With probability *prob*, convert
        a region of fixed tract length (*n*) if a recombination event
        happens. The actual number of markers converted depends on
        loci positions of surrounding loci. The starting position of
        this tract is the middle of two adjacent markers. For example,
        if four loci are located at ``0, 1, 2, 3`` respectively, a
        conversion event happens between ``0`` and ``1``, with a tract
        length 2 will start at 0.5 and end at 2.5, covering the second
        and third loci.

      + ``(EXPONENTIAL_DISTRIBUTION, prob, p)``: With probability
        *prob*, convert a region of random tract length if a
        recombination event happens. The distribution of tract length
        follows a exponential distribution with probability ``p``. The
        actual number of markers converted depends on loci positions
        of surrounding loci.

      simuPOP uses this probabilistic model of gene conversion because
      when a recombination event happens, it may become a
      recombination event if the if the Holliday junction is
      resolved/repaired successfully, or a conversion event if the
      junction is not resolved/repaired. The probability, however, is
      more commonly denoted by the ratio of conversion to
      recombination events in the literature. This ratio varies
      greatly from study to study, ranging from 0.1 to 15 (Chen et al,
      Nature Review Genetics, 2007). This translate to 0.1/0.9~0.1 to
      15/16~0.94 of the gene conversion probability.

      A :class:`Recombinator` usually does not send any output.
      However, if an information field is given (parameter
      *infoFields*), this operator will treat this information field
      as an unique ID of parents and offspring and output all
      recombination events in the format of ``offspring_id parent_id
      starting_ploidy loc1 loc2 ... `` where ``starting_ploidy``
      indicates which homologous copy genotype replication starts from
      (``0`` or ``1``), ``loc1``, ``loc2`` etc are loci after which
      recombination events happens. If there are multiple chromosomes
      on the genome, you will see a lot of (fake) recombination events
      because of independent segregation of chromosomes. Such a record
      will be generated for each set of homologous chromosomes so an
      diploid offspring will have two lines of output. Note that
      individual IDs need to be set (using a :class:`IdTagger`
      operator) before this  Recombinator is applied.

      In addition to genotypes, this operator also copies alleleic
      lineage if it is executed in a module with lineage allele type.

      .. note::

         There is no recombination between sex chromosomes
         (Chromosomes X and Y), although recombination is possible
         between pesudoautosomal regions on these chromosomes. If such
         a feature is required, you will have to simulate the
         pesudoautosomal regions as separate chromosomes.

   .. method:: Recombinator.transmitGenotype(parent, offspring, ploidy)

      This function transmits genotypes from a *parent* to the
      *ploidy-th* homologous set of chromosomes of an *offspring*. It
      can be used, for example, by a customized genotype transmitter
      to use sex-specific recombination rates to transmit parental
      genotypes to offspring.


