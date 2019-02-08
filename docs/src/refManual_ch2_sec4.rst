Pre-defined mating schemes
==========================


class CloneMating
-----------------

.. class:: CloneMating

   A homogeneous mating scheme that uses a sequential parent chooser and
   a clone offspring generator.

   .. method:: CloneMating.CloneMating(numOffspring=1, sexMode=None, ops=CloneGenoTransmitter(), subPopSize=[], subPops=ALL_AVAIL, weight=0, selectionField=None)

      Create a clonal mating scheme that clones parents to offspring using
      a :class:`CloneGenoTransmitter`. Please refer to class :class:`OffspringGenerator`
      for parameters *ops* and *numOffspring*, and to class :class:`HomoMating` for
      parameters  *subPopSize*, *subPops* and *weight*. Parameters *sexMode* and
      *selectionField* are ignored because this mating scheme does not support
      natural selection, and :class:`CloneGenoTransmitter` copies sex from parents
      to offspring. Note that :class:`CloneGenoTransmitter` by default also copies
      all parental information fields to offspring.



class RandomSelection
---------------------

.. class:: RandomSelection

   A homogeneous mating scheme that uses a random single-parent parent
   chooser with replacement, and a clone offspring generator. This mating
   scheme is usually used to simulate the basic haploid Wright-Fisher model
   but it can also be applied to diploid populations.

   .. method:: RandomSelection.RandomSelection(numOffspring=1, sexMode=None, ops=CloneGenoTransmitter(), subPopSize=[], subPops=ALL_AVAIL, weight=0, selectionField='fitness')

      Create a mating scheme that select a parent randomly and copy him or
      her to the offspring population. Please refer to class 
      :class:`RandomParentChooser` for parameter *selectionField*, to class
      :class:`OffspringGenerator` for parameters *ops* and *numOffspring*, and to
      class :class:`HomoMating` for parameters *subPopSize*, *subPops* and *weight*.
      Parameter *sexMode* is ignored because ``cloneOffspringGenerator`` copies
      sex from parents to offspring.



class RandomMating
------------------

.. class:: RandomMating

   A homogeneous mating scheme that uses a random parents chooser with
   replacement and a Mendelian offspring generator. This mating scheme is
   widely used to simulate diploid sexual Wright-Fisher random mating.

   .. method:: RandomMating.RandomMating(numOffspring=1, sexMode=RANDOM_SEX, ops=MendelianGenoTransmitter(), subPopSize=[], subPops=ALL_AVAIL, weight=0, selectionField='fitness')

      Creates a random mating ssheme that selects two parents randomly and
      transmit genotypes according to Mendelian laws. Please refer to class
      :class:`RandomParentsChooser` for parameter *selectionField*, to class
      :class:`OffspringGenerator` for parameters *ops*, *sexMode* and
      *numOffspring*, and to class :class:`HomoMating` for parameters
      *subPopSize*, *subPops* and *weight*.



class MonogamousMating
----------------------

.. class:: MonogamousMating

   A homogeneous mating scheme that uses a random parents chooser without
   replacement and a Mendelian offspring generator. It differs from the basic
   random mating scheme in that each parent can mate only once so there is no
   half-sibling in the population.

   .. method:: MonogamousMating.MonogamousMating(numOffspring=1, sexMode=RANDOM_SEX, ops=MendelianGenoTransmitter(), subPopSize=[], subPops=ALL_AVAIL, weight=0, selectionField=None)

      Creates a monogamous mating scheme that selects each parent only
      once. Please refer to class :class:`OffspringGenerator` for parameters
      *ops*, *sexMode* and *numOffspring*, and to class :class:`HomoMating` for
      parameters *subPopSize*, *subPops* and *weight*. Parameter
      *selectionField* is ignored because this mating scheme does not
      support natural selection.



class PolygamousMating
----------------------

.. class:: PolygamousMating

   A homogeneous mating scheme that uses a multi-spouse parents chooser
   and a Mendelian offspring generator. It differs from the basic random
   mating scheme in that each parent of sex *polySex* will have *polyNum*
   spouses.

   .. method:: PolygamousMating.PolygamousMating(polySex=MALE, polyNum=1, numOffspring=1, sexMode=RANDOM_SEX, ops=MendelianGenoTransmitter(), subPopSize=[], subPops=ALL_AVAIL, weight=0, selectionField='fitness')

      Creates a polygamous mating scheme that each parent mates with
      multiple spouses. Please refer to class :class:`PolyParentsChooser` for
      parameters *polySex*, *polyNum* and *selectionField*, to class
      :class:`OffspringGenerator` for parameters *ops*,  *sexMode* and
      *numOffspring*, and to class :class:`HomoMating` for parameters
      *subPopSize*, *subPops* and *weight*.



class HaplodiploidMating
------------------------

.. class:: HaplodiploidMating

   A homogeneous mating scheme that uses a random parents chooser with
   replacement and a haplodiploid offspring generator. It should be used
   in a haplodiploid population where male individuals only have one set
   of homologous chromosomes.

   .. method:: HaplodiploidMating.HaplodiploidMating(numOffspring=1.0, sexMode=RANDOM_SEX, ops=HaplodiploidGenoTransmitter(), subPopSize=[], subPops=ALL_AVAIL, weight=0, selectionField='fitness')

      Creates a mating scheme in haplodiploid populations. Please refer
      to class :class:`RandomParentsChooser` for parameter *selectionField*, to
      class :class:`OffspringGenerator` for parameters *ops*, *sexMode* and
      *numOffspring*, and to class :class:`HomoMating` for parameters
      *subPopSize*, *subPops* and *weight*.



class SelfMating
----------------

.. class:: SelfMating

   A homogeneous mating scheme that uses a random single-parent parent
   chooser with or without replacement (parameter *replacement*) and a
   selfing offspring generator. It is used to mimic self-fertilization
   in certain plant populations.

   .. method:: SelfMating.SelfMating(replacement=True, numOffspring=1, sexMode=RANDOM_SEX, ops=SelfingGenoTransmitter(), subPopSize=[], subPops=ALL_AVAIL, weight=0, selectionField='fitness')

      Creates a selfing mating scheme where two homologous copies of
      parental chromosomes are transmitted to offspring according to
      Mendelian laws. Please refer to class :class:`RandomParentChooser` for
      parameter *replacement* and  *selectionField*, to class
      :class:`OffspringGenerator` for parameters *ops*, *sexMode* and
      *numOffspring*, and to class :class:`HomoMating` for parameters
      *subPopSize*, *subPops* and *weight*.



class HermaphroditicMating
--------------------------

.. class:: HermaphroditicMating

   A hermaphroditic mating scheme that chooses two parents randomly
   from the population regardless of sex. The parents could be chosen
   with or without replacement (parameter *replacement*). Selfing (if
   the same parents are chosen) is allowed unless *allowSelfing* is 
   set to *False*

   .. method:: HermaphroditicMating.HermaphroditicMating(replacement=True, allowSelfing=True, numOffspring=1, sexMode=RANDOM_SEX, ops=MendelianGenoTransmitter(), subPopSize=[], subPops=ALL_AVAIL, weight=0, selectionField='fitness')

      Creates a hermaphroditic mating scheme where individuals can
      serve as father or mother, or both (self-fertilization). Please 
      refer to class :class:`CombinedParentsChooser` for parameter *allowSelfing``,
      to :class:`RandomParentChooser` for parameter *replacement* and
      *selectionField*, to class :class:`OffspringGenerator` for parameters *ops*,
      *sexMode* and *numOffspring*, and to class :class:`HomoMating` for parameters
      *subPopSize*, *subPops* and *weight*.



class ControlledRandomMating
----------------------------

.. class:: ControlledRandomMating

   A homogeneous mating scheme that uses a random sexual parents chooser
   with replacement and a controlled offspring generator using Mendelian
   genotype transmitter. It falls back to a regular random mating scheme
   if there is no locus to control or no trajectory is defined.

   .. method:: ControlledRandomMating.ControlledRandomMating(loci=[], alleles=[], freqFunc=None, numOffspring=1, sexMode=RANDOM_SEX, ops=MendelianGenoTransmitter(), subPopSize=[], subPops=ALL_AVAIL, weight=0, selectionField='fitness')

      Creates a random mating scheme that controls allele frequency at
      loci *loci*. At each generation, function *freqFunc* will be called to
      called to obtain intended frequencies of alleles *alleles* at loci
      *loci*. The controlled offspring generator will control the acceptance
      of offspring so that the generation reaches desired allele frequencies
      at these loci. If *loci* is empty or *freqFunc* is ``None``, this mating
      scheme works identically to a ``RandomMating scheme``. Rationals and
      applications of this mating scheme is described in details in a paper *Peng
      et al, 2007 (PLoS Genetics)*. Please refer to class :class:`RandomParentsChooser`
      for parameters *selectionField*, to class :class:`ControlledOffspringGenerator`
      for parameters *loci*, *alleles*, *freqFunc*, to class
      :class:`OffspringGenerator` for parameters *ops*, *sexMode* and *numOffspring*,
      and to class :class:`HomoMating` for parameters *subPopSize*, *subPops* and
      *weight*.



