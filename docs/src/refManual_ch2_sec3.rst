Mating Schemes
==============

.. index:: single: mating scheme


class MatingScheme
------------------

.. class:: MatingScheme

   This mating scheme is the base class of all mating schemes. It
   evolves a population generation by generation but does not actually
   transmit genotype.


   .. method:: MatingScheme(subPopSize=[])


      Create a base mating scheme that evolves a population without
      transmitting genotypes. At each generation, this mating scheme
      creates an offspring generation according to parameter
      *subPopSize*, which can be a list of subpopulation sizes (or a
      number if there is only one subpopulation) or a Python function
      which will be called at each generation, just before mating, to
      determine the subpopulation sizes of the offspring generation.
      The function should be defined with one or both parameters of
      ``gen`` and ``pop`` where ``gen`` is the current generation
      number and ``pop`` is the parental population just before
      mating. The return value of this function should be a list of
      subpopulation sizes for the offspring generation. A single
      number can be returned if there is only one subpopulation. The
      passed parental population is usually used to determine
      offspring population size from parental population size but you
      can also modify this population to prepare for mating. A common
      practice is to split and merge parental populations in this
      function so that you demographic related information and actions
      could be implemented in the same function.



class HomoMating
----------------

.. class:: HomoMating

   A homogeneous mating scheme that uses a parent chooser to choose
   parents from a prental generation, and an offspring generator to
   generate offspring from chosen parents. It can be either used
   directly, or within a heterogeneous mating scheme. In the latter
   case, it can be applied to a (virtual) subpopulation.


   .. method:: HomoMating(chooser, generator, subPopSize=[], subPops=ALL_AVAIL, weight=0)


      Create a homogeneous mating scheme using a parent chooser
      *chooser* and an offspring generator *generator*.

      If this mating scheme is used directly in a simulator, it will
      be responsible for creating an offspring population according to
      parameter *subPopSize*. This parameter can be a list of
      subpopulation sizes (or a number if there is only one
      subpopulation) or a Python function which will be called at each
      generation to determine the subpopulation sizes of the offspring
      generation. Please refer to class :class:`MatingScheme` for
      details about this parameter.

      If this mating shcme is used within a heterogeneous mating
      scheme. Parameters *subPops* and *weight* are used to determine
      which (virtual) subpopulations this mating scheme will be
      applied to, and how many offspring this mating scheme will
      produce. Please refer to mating scheme :class:`HeteroMating` for
      the use of these two parameters.



class HeteroMating
------------------

.. class:: HeteroMating

   A heterogeneous mating scheme that applies a list of homogeneous
   mating schemes to different (virtual) subpopulations.


   .. method:: HeteroMating(matingSchemes, subPopSize=[], shuffleOffspring=True, weightBy=ANY_SEX)


      Create a heterogeneous mating scheme that will apply a list of
      homogeneous mating schemes *matingSchemes* to different
      (virtual) subpopulations. The size of the offspring generation
      is determined by parameter *subPopSize*, which can be a list of
      subpopulation sizes or a Python function that returns a list of
      subpopulation sizes at each generation. Please refer to class
      :class:`MatingScheme` for a detailed explanation of this
      parameter.

      Each mating scheme defined in *matingSchemes* can be applied to
      one or more (virtual) subpopulation. If parameter *subPops* is
      not specified, a mating scheme will be applied to all
      subpopulations. If a list of (virtual) subpopulation is
      specified, the mating scheme will be applied to specific
      (virtual) subpopulations.

      If multiple mating schemes are applied to the same
      subpopulation, a weight (parameter *weight*) can be given to
      each mating scheme to determine how many offspring it will
      produce. The default \weight for all mating schemes are ``0``.
      In this case, the number of offspring each mating scheme
      produces is proportional to the number of individuals in its
      parental (virtual) subpopulation (default to all parents, but
      can be father for ``weightBy=MALE_ONLY``, mother for
      ``weightBy=FEMALE_ONLY``, or father mother pairs (less of number
      of father and mothers) for ``weightBy=PAIR_ONLY``). If all
      weights are negative, the numbers of offspring are determined by
      the multiplication of the absolute values of the weights and
      their respective parental (virtual) subpopulation sizes. If all
      weights are positive, the number of offspring produced by each
      mating scheme is proportional to these weights, except for
      mating schemes with zero parental population size (or no father,
      no mother, or no pairs, depending on value of parameter
      ``weightBy``). Mating schemes with zero weight in this case will
      produce no offspring. If both negative and positive weights are
      present, negative weights are processed before positive ones.

      A sexual mating scheme might fail if a parental (virtual)
      subpopulation has no father or mother. In this case, you can set
      ``weightBy`` to ``PAIR_ONLY`` so a (virtual) subpopulation will
      appear to have zero size, and will thus contribute no offspring
      to the offspring population. Note that the perceived parental
      (virtual) subpopulation size in this mode (and in modes of
      ``MALE_ONLY``, ``FEMALE_ONLY``) during the calculation of the
      size of the offspring subpopulation will be roughly half of the
      actual population size so you might have to use ``weight=-2`` if
      you would like to have an offspring subpopulation that is
      roughly the same size of the parental (virtual) subpopulation.

      If multiple mating schemes are applied to the same
      subpopulation, offspring produced by these mating schemes are
      shuffled randomly. If this is not desired, you can turn off
      offspring shuffling by setting parameter *shuffleOffspring* to
      ``False``.



class ConditionalMating
-----------------------

.. class:: ConditionalMating

   A conditional mating scheme that applies different mating schemes
   according to a condition (similar to operator  IfElse). The
   condition can be a fixed condition, an expression or a user-defined
   function, to determine which mating scheme to be used.


   .. method:: ConditionalMating(cond, ifMatingScheme, elseMatingScheme)


      Create a conditional mating scheme that applies mating scheme
      *ifMatingScheme* if the condition *cond* is ``True``, or
      *elseMatingScheme* if *cond* is ``False``. If a Python
      expression (a string) is given to parameter *cond*, the
      expression will be evalulated in parental population's local
      namespace. When a Python function is specified, it accepts
      parameter ``pop`` for the parental population. The return value
      of this function should be ``True`` or ``False``. Otherwise,
      parameter *cond* will be treated as a fixed condition (converted
      to ``True`` or ``False``) upon which *ifMatingScheme* or
      *elseMatingScheme* will alway be applied.



class PedigreeMating
--------------------

.. class:: PedigreeMating

   This mating scheme evolves a population following an existing
   pedigree structure. If the :class:`Pedigree` object has ``N``
   ancestral generations and a present generation, it can be used to
   evolve a population for ``N`` generations, starting from the
   topmost ancestral generation. At the *k-th* generation, this mating
   scheme produces an offspring generation according to subpopulation
   structure of the ``N-k-1`` ancestral generation in the pedigree
   object (e.g. producing the offspring population of generation 0
   according to the ``N-1`` ancestral generation of the pedigree
   object ). For each offspring, this mating scheme copies individual
   ID and sex from the corresponing individual in the pedigree object.
   It then locates the parents of each offspring using their IDs in
   the pedigree object. A list of during mating operators are then
   used to transmit parental genotype to the offspring. The population
   being evolved must have an information field ``'ind_id'``.


   .. method:: PedigreeMating(ped, ops, idField="ind_id")


      Creates a pedigree mating scheme that evolves a population
      according to :class:`Pedigree` object *ped*. The evolved
      population should contain individuals with ID (at information
      field *idField*, default to ``'ind_id'``) that match those
      individual in the topmost ancestral generation who have
      offspring. After parents of each individuals are determined from
      their IDs, a list of during-mating operators *ops* are applied
      to transmit genotypes. The return value of these operators are
      not checked.


   .. method:: PedigreeMating.parallelizable()

      FIXME: No document


class SequentialParentChooser
-----------------------------

.. class:: SequentialParentChooser

   This parent chooser chooses a parent from a parental (virtual)
   subpopulation sequentially. Natural selection is not considered. If
   the last parent is reached, this parent chooser will restart from
   the beginning of the (virtual) subpopulation.


   .. method:: SequentialParentChooser(sexChoice=ANY_SEX)


      Create a parent chooser that chooses a parent from a parental
      (virtual) subpopulation sequentially. Parameter *choice* can be
      ``ANY_SEX`` (default), ``MALE_ONLY`` and ``FEMALE_ONLY``. In the
      latter two cases, only male or female individuals are selected.
      A :class:`RuntimeError` will be raised if there is no male or
      female individual from the population.


   .. method:: SequentialParentChooser.chooseParents()

      Return chosen parents from a population if the parent chooser
      object is created with a population.

   .. method:: SequentialParentChooser.initialize(pop, subPop)

      Initialize a parent chooser for subpopulation *subPop* of
      *population* pop.


class SequentialParentsChooser
------------------------------

.. class:: SequentialParentsChooser

   This parent chooser chooses two parents (a father and a mother)
   sequentially from their respective sex groups. Selection is not considered.
   If all fathers (or mothers) are exhausted, this parent chooser will choose
   fathers (or mothers) from the beginning of the (virtual) subpopulation
   again.

   .. method:: SequentialParentsChooser.SequentialParentsChooser()

      Create a parent chooser that chooses two parents sequentially from a
      parental (virtual) subpopulation.



class RandomParentChooser
-------------------------

.. class:: RandomParentChooser

   This parent chooser chooses a parent randomly from a (virtual)
   parental subpopulation. Parents are chosen with or without
   replacement. If parents are chosen with replacement, a parent can
   be selected multiple times. If individual fitness values are
   assigned to individuals ( stored in an information field
   *selectionField* (default to ``"fitness"``), individuals will be
   chosen at a probability proportional to his or her fitness value.
   If parents are chosen without replacement, a parent can be chosen
   only once. An :class:`RuntimeError` will be raised if all parents
   are exhausted. Natural selection is disabled in the without-
   replacement case.


   .. method:: RandomParentChooser(replacement=True, selectionField="fitness", sexChoice=ANY_SEX)


      Create a random parent chooser that choose parents with or
      without replacement (parameter *replacement*, default to
      ``True``). If selection is enabled and information field
      *selectionField* exists in the passed population, the
      probability that a parent is chosen is proportional to his/her
      fitness value stored in *selectionField*. This parent chooser by
      default chooses parent from all individuals (``ANY_SEX``), but
      it can be made to select only male (``MALE_ONLY``) or female
      (``FEMALE_ONLY``) individuals by setting parameter *sexChoice*.


   .. method:: RandomParentChooser.chooseParents()

      Return chosen parents from a population if the parent chooser
      object is created with a population.

   .. method:: RandomParentChooser.initialize(pop, subPop)

      Initialize a parent chooser for subpopulation *subPop* of
      *population* pop.


class RandomParentsChooser
--------------------------

.. class:: RandomParentsChooser

   This parent chooser chooses two parents, a male and a female,
   randomly from a (virtual) parental subpopulation. Parents are
   chosen with or without replacement from their respective sex group.
   If parents are chosen with replacement, a parent can be selected
   multiple times. If individual fitness values are assigned (stored
   in information field *selectionField*, default to ``"fitness"``,
   the probability that an individual is chosen is proportional to
   his/her fitness value among all individuals with the same sex. If
   parents are chosen without replacement, a parent can be chosen only
   once. An :class:`RuntimeError` will be raised if all males or
   females are exhausted. Natural selection is disabled in the
   without-replacement case.


   .. method:: RandomParentsChooser(replacement=True, selectionField="fitness")


      Create a random parents chooser that choose two parents with or
      without replacement (parameter *replacement*, default to
      ``True``). If selection is enabled and information field
      *selectionField* exists in the passed population, the
      probability that a parent is chosen is proportional to his/her
      fitness value stored in *selectionField*.


   .. method:: RandomParentsChooser.chooseParents()

      Return chosen parents from a population if the parent chooser
      object is created with a population.

   .. method:: RandomParentsChooser.initialize(pop, subPop)

      Initialize a parent chooser for subpopulation *subPop* of
      *population* pop.


class PolyParentsChooser
------------------------

.. class:: PolyParentsChooser

   This parent chooser is similar to random parents chooser but
   instead of selecting a new pair of parents each time, one of the
   parents in this parent chooser will mate with several spouses
   before he/she is replaced. This mimicks multi-spouse mating schemes
   such as polygyny or polyandry in some populations. Natural
   selection is supported for both sexes.


   .. method:: PolyParentsChooser(polySex=MALE, polyNum=1, selectionField="fitness")


      Create a multi-spouse parents chooser where each father (if
      *polySex* is MALE) or mother (if *polySex* is FEMALE) has
      *polyNum* spouses. The parents are chosen with replacement. If
      individual fitness values are assigned (stored to information
      field ``selectionField``, default to ``"fitness"``), the
      probability that an individual is chosen is proportional to
      his/her fitness value among all individuals with the same sex.


   .. method:: PolyParentsChooser.chooseParents()

      Return chosen parents from a population if the parent chooser
      object is created with a population.

   .. method:: PolyParentsChooser.initialize(pop, subPop)

      Initialize a parent chooser for subpopulation *subPop* of
      *population* pop.


class CombinedParentsChooser
----------------------------

.. class:: CombinedParentsChooser

   This parent chooser accepts two parent choosers. It takes one
   parent from each parent chooser and return them as father and
   mother. Because two parent choosers do not have to choose parents
   from the same virtual subpopulation, this parent chooser allows you
   to choose parents from different subpopulations.


   .. method:: CombinedParentsChooser(fatherChooser, motherChooser, allowSelfing=True)


      Create a Python parent chooser using two parent choosers
      *fatherChooser* and *motherChooser*. It takes one parent from
      each parent chooser and return them as father and mother. If two
      valid parents are returned, the first valid parent (father) will
      be used for *fatherChooser*, the second valid parent (mother)
      will be used for *motherChooser*. Although these two parent
      choosers are supposed to return a father and a mother
      respectively, the sex of returned parents are not checked so it
      is possible to return parents with the same sex using this
      parents chooser. This choose by default allows the selection of
      the same parents as father and mother (self-fertilization),
      unless a parameter *allowSelfing* is used to disable it.


   .. method:: CombinedParentsChooser.chooseParents()

      Return chosen parents from a population if the parent chooser
      object is created with a population.

   .. method:: CombinedParentsChooser.initialize(pop, subPop)

      Initialize a parent chooser for subpopulation *subPop* of
      *population* pop.


class PyParentsChooser
----------------------

.. class:: PyParentsChooser

   This parent chooser accepts a Python generator function that
   repeatedly yields one or two parents, which can be references to
   individual objects or indexes relative to each subpopulation. The
   parent chooser calls the generator function with parental
   population and a subpopulation index for each subpopulation and
   retrieves parents repeatedly using the iterator interface of the
   generator function.

   This parent chooser does not support virtual subpopulation
   directly. However, because virtual subpopulations are defined in
   the passed parental population, it is easy to return parents from a
   particular virtual subpopulation using virtual subpopulation
   related functions.


   .. method:: PyParentsChooser(generator)


      Create a Python parent chooser using a Python generator function
      *parentsGenerator*. This function should accept one or both of
      parameters *pop* (the parental population) and *subPop* (index
      of subpopulation) and return the reference or index (relative to
      subpopulation) of a parent or a pair of parents repeatedly using
      the iterator interface of the generator function.


   .. method:: PyParentsChooser.chooseParents()

      Return chosen parents from a population if the parent chooser
      object is created with a population.

   .. method:: PyParentsChooser.initialize(pop, subPop)

      Initialize a parent chooser for subpopulation *subPop* of
      *population* pop.


class OffspringGenerator
------------------------

.. class:: OffspringGenerator

   An *offspring generator* generates offspring from parents chosen by
   a parent chooser. It is responsible for creating a certain number
   of offspring, determinning their sex, and transmitting genotypes
   from parents to offspring.


   .. method:: OffspringGenerator(ops, numOffspring=1, sexMode=RANDOM_SEX)


      Create a basic offspring generator. This offspring generator
      uses *ops* genotype transmitters to transmit genotypes from
      parents to offspring.

      A number of *during-mating operators* (parameter *ops*) can be
      used to, among other possible duties such as setting information
      fields of offspring, transmit genotype from parents to
      offspring. This general offspring generator does not have any
      default during-mating operator but all stock mating schemes use
      an offspring generator with a default operator. For example, a
      ``mendelianOffspringGenerator`` is used by :class:`RandomMating`
      to trasmit genotypes. Note that applicability parameters
      ``begin``, ``step``, ``end``, ``at`` and ``reps`` could be used
      in these operators but negative population and generation
      indexes are unsupported.

      Parameter *numOffspring* is used to control the number of
      offspring per mating event, or in another word the number of
      offspring in each family. It can be a number, a Python function
      or generator, or a mode parameter followed by some optional
      arguments. If a number is given, given number of offspring will
      be generated at each mating event. If a Python function is
      given, it will be called each time when a mating event happens.
      When a generator function is specified, it will be called for
      each subpopulation to provide number of offspring for all mating
      events during the populating of this subpopulation. Current
      generation number will be passed to this function or generator
      function if parameter "gen" is used in this function. In the
      last case, a tuple (or a list) in one of the following forms can
      be given:

      + ``(GEOMETRIC_DISTRIBUTION, p)``

      + ``(POISSON_DISTRIBUTION, p)``, p > 0

      + ``(BINOMIAL_DISTRIBUTION, p, N)``, 0 < p <=1, N > 0

      + ``(UNIFORM_DISTRIBUTION, a, b)``, 0 <= a <= b.

      In this case, the number of offspring will be determined
      randomly following the specified statistical distributions.
      Because families with zero offspring are silently ignored, the
      distribution of the observed number of offspring per mating
      event (excluding zero) follows zero-truncated versions of these
      distributions.

      Parameter *numOffspring* specifies the number of offspring per
      mating event but the actual surviving offspring can be less than
      specified. More spefically, if any during-mating operator
      returns ``False``, an offspring will be discarded so the
      actually number of offspring of a mating event will be reduced.
      This is essentially how during-mating selector works.

      Parameter *sexMode* is used to control the sex of each
      offspring. Its default value is usually *RANDOM_SEX* which
      assign ``MALE`` or ``FEMALE`` to each individual randomly, with
      equal probabilities. If ``NO_SEX`` is given, offspring sex will
      not be changed. *sexMode* can also be one of

      + ``(PROB_OF_MALES, p)`` where ``p`` is the probability of male
        for each offspring,

      + ``(NUM_OF_MALES, n)`` where ``n`` is the number of males in a
        mating event. If ``n`` is greater than or equal to the number
        of offspring in this family, all offspring in this family will
        be ``MALE``.

      + ``(NUM_OF_FEMALES, n)`` where ``n`` is the number of females
        in a mating event,

      + ``(SEQUENCE_OF_SEX, s1, s2 ...)`` where ``s1``, ``s2`` etc are
        MALE or FEMALE. The sequence will be used for each mating
        event. It will be reused if the number of offspring in a
        mating event is greater than the length of sequence.

      + ``(GLOBAL_SEQUENCE_OF_SEX, s1, s2, ...)`` where ``s1``, ``s2``
        etc are MALE or FEMALE. The sequence will be used across
        mating events. It will be reused if the number of offspring in
        a subpopulation is greater than the length of sequence.

      Finally, parameter *sexMode* accepts a function or a generator
      function. A function will be called whenever an offspring is
      produced. A generator will be created at each subpopulation and
      will be used to produce sex for all offspring in this
      subpopulation. No parameter is accepted.



class ControlledOffspringGenerator
----------------------------------

.. class:: ControlledOffspringGenerator

   This offspring generator populates an offspring population and
   controls allele frequencies at specified loci. At each generation,
   expected allele frequencies at these loci are passed from a user
   defined allele frequency *trajectory* function. The offspring
   population is populated in two steps. At the first step, only
   families with disease alleles are accepted until until the expected
   number of disease alleles are met. At the second step, only
   families with wide type alleles are accepted to populate the rest
   of the offspring generation. This method is described in detail in
   "Peng et al, (2007) PLoS Genetics".


   .. method:: ControlledOffspringGenerator(loci, alleles, freqFunc, ops=[], numOffspring=1, sexMode=RANDOM_SEX)


      Create an offspring generator that selects offspring so that
      allele frequency at specified loci in the offspring generation
      reaches specified allele frequency. At the beginning of each
      generation, expected allele frequency of *alleles* at *loci* is
      returned from a user-defined trajectory function *freqFunc*.
      Parameter *loci* can be a list of loci indexes, names, or
      ALL_AVAIL. If there is no subpopulation, this function should
      return a list of frequencies for each locus. If there are
      multiple subpopulations, *freqFunc* can return a list of allele
      frequencies for all subpopulations or combined frequencies that
      ignore population structure. In the former case, allele
      frequencies should be arranged by loc0_sp0, loc1_sp0, ...
      loc0_sp1, loc1_sp1, ..., and so on. In the latter case, overall
      expected number of alleles are scattered to each subpopulation
      in proportion to existing number of alleles in each
      subpopulation, using a multinomial distribution.

      After the expected alleles are calculated, this offspring
      generator accept and reject families according to their genotype
      at *loci* until allele frequecies reach their expected values.
      The rest of the offspring generation is then filled with
      families without only wild type alleles at these *loci*.

      This offspring generator is derived from class
      *OffspringGenerator*. Please refer to class *OffspringGenerator*
      for a detailed description of parameters *ops*, *numOffspring*
      and *sexMode*.



