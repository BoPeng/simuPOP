Mutation
========

.. index:: single: Mutation


class BaseMutator
-----------------

.. class:: BaseMutator

   Class ``mutator`` is the base class of all mutators. It handles all
   the work of picking an allele at specified loci from certain
   (virtual) subpopulation with certain probability, and calling a
   derived mutator to mutate the allele. Alleles can be changed before
   and after mutation if existing allele numbers do not match those of
   a mutation model.


   .. method:: BaseMutator(rates=[], loci=ALL_AVAIL, mapIn=[], mapOut=[], context=0, output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=["ind_id"], lineageMode=FROM_INFO)


      A mutator mutates alleles from one state to another with given
      probability. This base mutator does not perform any mutation but
      it defines common behaviors of all mutators.

      By default, a mutator mutates all alleles in all populations of
      a simulator at all generations. A number of parameters can be
      used to restrict mutations to certain generations (parameters
      *begin*, *end*, *step* and *at*), replicate populations
      (parameter *rep*), (virtual) subpopulations (parameter
      *subPops*) and loci (parameter *loci*). Parameter *loci* can be
      a list of loci indexes, names, list of chromosome position
      pairs, ``ALL_AVAIL``, or a function with optional parameter
      ``pop`` that will be called at each ganeeration to determine
      indexes of loci. Please refer to class :class:`BaseOperator` for
      a detailed explanation of these parameters.

      Parameter *rate* or its equivalence specifies the probability
      that a a mutation event happens. The exact form and meaning of
      *rate* is mutator-specific. If a single rate is specified, it
      will be applied to all *loci*. If a list of mutation rates are
      given, they will be applied to each locus specified in parameter
      *loci*. Note that not all mutators allow specification of
      multiple mutation rate, especially when the mutation rate itself
      is a list or matrix.

      Alleles at a locus are non-negative numbers 0, 1, ... up to the
      maximum allowed allele for the loaded module (1 for binary, 255
      for short and 65535 for long modules). Whereas some general
      mutation models treat alleles as numbers, other models assume
      specific interpretation of alleles. For example, an
      :class:`AcgtMutator` assumes alleles ``0``, ``1``, ``2`` and
      ``3`` as nucleotides ``A``, ``C``, ``G`` and ``T``. Using a
      mutator that is incompatible with your simulation will certainly
      yield erroneous results.

      If your simulation assumes different alleles with a mutation
      model, you can map an allele to the allele used in the model and
      map the mutated allele back. This is achieved using a *mapIn*
      list with its ``i-th`` item being the corresponding allele of
      real allele ``i``, and a *mapOut* list with its *i-th* item
      being the real allele of allele ``i`` assumed in the model. For
      example ``mapIn=[0, 0, 1]`` and ``mapOut=[1, 2]`` would allow
      the use of a :class:`SNPMutator` to mutate between alleles 1 and
      2, instead of 0 and 1. Parameters *mapIn* and *mapOut* also
      accept a user-defined Python function that returns a
      corresponding allele for a given allele. This allows easier
      mapping between a large number of alleles and advanced models
      such as random emission of alleles.

      If a valid information field is specified for parameter
      *infoFields* (default to ``ind_id``) for modules with lineage
      allele type, the lineage of the mutated alleles will be the ID
      (stored in the first field of *infoFields*) of individuals that
      harbor the mutated alleles if *lineageMode* is set to
      ``FROM_INFO`` (default). If *lineageMode* is set to
      ``FROM_INFO_SIGNED``, the IDs will be assigned a sign depending
      on the ploidy the mutation happens (1 for ploidy 0, -1 for
      ploidy 1, etc). The lineage information will be transmitted
      along with the alleles so this feature allows you to track the
      source of mutants during evolution.A

      A mutator by default does not produce any output. However, if an
      non-empty output is specified, the operator will output
      generation number, locus, ploidy, original allele, mutant, and
      values of all information field specified by parameter
      ``infoFields`` (e.g. individual ID if ``ind_id`` is specified).

      Some mutation models are context dependent. Namely, how an
      allele mutates will depend on its adjecent alleles. Whereas most
      simuPOP mutators are context independent, some of them accept a
      parameter *context* which is the number of alleles to the left
      and right of the mutated allele. For example *context=1* will
      make the alleles to the immediate left and right to a mutated
      allele available to a mutator. These alleles will be mapped in
      if parameter *mapIn* is defined. How exactly a mutator makes use
      of these information is mutator dependent.



class MatrixMutator
-------------------

.. class:: MatrixMutator

   A matrix mutator mutates alleles ``0``, ``1``, ..., ``n-1`` using a
   ``n`` by ``n`` matrix, which specifies the probability at which
   each allele mutates to another. Conceptually speaking, this mutator
   goes through all mutable allele and mutate it to another state
   according to probabilities in the corresponding row of the rate
   matrix. Only one mutation rate matrix can be specified which will
   be used for all specified loci. #


   .. method:: MatrixMutator(rate, loci=ALL_AVAIL, mapIn=[], mapOut=[], output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=["ind_id"], lineageMode=FROM_INFO)


      Create a mutator that mutates alleles ``0``, ``1``, ..., ``n-1``
      using a ``n`` by ``n`` matrix ``rate``. Item ``(i,j)`` of this
      matrix specifies the probability at which allele *i* mutates to
      allele *j*. Diagnal items ``(i, i)`` are ignored because they
      are automatically determined by other probabilities. Only one
      mutation rate matrix can be specified which will be used for all
      loci in the applied population, or loci specified by parameter
      *loci*. If alleles other than ``0``, ``1``, ..., ``n-1`` exist
      in the population, they will not be mutated although a warning
      message will be given if debugging code ``DBG_WARNING`` is
      turned on. Please refer to classes ``mutator`` and
      :class:`BaseOperator` for detailed explanation of other
      parameters.



class KAlleleMutator
--------------------

.. class:: KAlleleMutator

   This mutator implements a *k-allele* mutation model that assumes
   *k* allelic states (alleles 0, 1, 2, ..., *k-1*) at each locus.
   When a mutation event happens, it mutates an allele to any other
   states with equal probability.


   .. method:: KAlleleMutator(k, rates=[], loci=ALL_AVAIL, mapIn=[], mapOut=[], output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=["ind_id"], lineageMode=FROM_INFO)


      Create a k-allele mutator that mutates alleles to one of the
      other ``k-1`` alleles with equal probability. This mutator by
      default applies to all loci unless parameter *loci* is
      specified. A single mutation rate will be used for all loci if a
      single value of parameter *rates* is given. Otherwise, a list of
      mutation rates can be specified for each locus in parameter
      *loci*. If the mutated allele is larger than or equal to ``k``,
      it will not be mutated. A warning message will be displayed if
      debugging code ``DBG_WARNING`` is turned on. Please refer to
      classes ``mutator`` and :class:`BaseOperator` for descriptions
      of other parameters.



class StepwiseMutator
---------------------

.. class:: StepwiseMutator

   A stepwise mutation model treats alleles at a locus as the number
   of tandem repeats of microsatellite or minisatellite markers. When
   a mutation event happens, the number of repeats (allele) either
   increase or decrease. A standard stepwise mutation model increases
   of decreases an allele by 1 with equal probability. More complex
   models (generalized stepwise mutation model) are also allowed. Note
   that an allele cannot be mutated beyond boundaries (0 and maximum
   allowed allele).


   .. method:: StepwiseMutator(rates=[], loci=ALL_AVAIL, incProb=0.5, maxAllele=0, mutStep=[], mapIn=[], mapOut=[], output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=["ind_id"], lineageMode=FROM_INFO)


      Create a stepwise mutation mutator that mutates an allele by
      increasing or decreasing it. This mutator by default applies to
      all loci unless parameter *loci* is specified. A single mutation
      rate will be used for all loci if a single value of parameter
      *rates* is given. Otherwise, a list of mutation rates can be
      specified for each locus in parameter *loci*.

      When a mutation event happens, this operator increases or
      decreases an allele by *mutStep* steps. Acceptable input of
      parameter *mutStep* include

      + A number: This is the default mode with default value 1.

      + ``(GEOMETRIC_DISTRIBUTION, p)``: The number of steps follows a
        a geometric distribution with parameter *p*.

      + A Python function: This user defined function accepts the
        allele being mutated and return the steps to mutate.

      The mutation process is usually neutral in the sense that
      mutating up and down is equally likely. You can adjust parameter
      *incProb* to change this behavior.

      If you need to use other generalized stepwise mutation models,
      you can implement them using a :class:`PyMutator`. If
      performance becomes a concern, I may add them to this operator
      if provided with a reliable reference.



class PyMutator
---------------

.. class:: PyMutator

   This hybrid mutator accepts a Python function that determines how
   to mutate an allele when an mutation event happens.


   .. method:: PyMutator(rates=[], loci=ALL_AVAIL, func=None, context=0, mapIn=[], mapOut=[], output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=["ind_id"], lineageMode=FROM_INFO)


      Create a hybrid mutator that uses a user-provided function to
      mutate an allele when a mutation event happens. This function
      (parameter *func*) accepts the allele to be mutated as parameter
      ``allele``, locus index ``locus``, and optional array of alleles
      as parameter ``context``, which are *context* alleles the left
      and right of the mutated allele. Invalid context alleles (e.g.
      left allele to the first locus of a chromosome) will be marked
      by -1. The return value of this function will be used to mutate
      the passed allele. The passed, returned and context alleles
      might be altered if parameter *mapIn* and *mapOut* are used.
      This mutator by default applies to all loci unless parameter
      *loci* is specified. A single mutation rate will be used for all
      loci if a single value of parameter *rates* is given. Otherwise,
      a list of mutation rates can be specified for each locus in
      parameter *loci*. Please refer to classes ``mutator`` and
      :class:`BaseOperator` for descriptions of other parameters.



class MixedMutator
------------------

.. class:: MixedMutator

   This mixed mutator accepts a list of mutators and use one of them
   to mutate an allele when an mutation event happens.


   .. method:: MixedMutator(rates=[], loci=ALL_AVAIL, mutators=[], prob=[], mapIn=[], mapOut=[], context=0, output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=["ind_id"], lineageMode=FROM_INFO)


      Create a mutator that randomly chooses one of the specified
      *mutators* to mutate an allele when a mutation event happens.
      The mutators are choosen according to a list of probabilities
      (\parameter *prob*) that should add up to ``1``. The passed and
      returned alleles might be changed if parameters *mapIn* and
      *mapOut* are used. Most parameters, including *loci*, *mapIn*,
      *mapOut*, *rep*, and *subPops* of mutators specified in
      parameter *mutators* are ignored. This mutator by default
      applies to all loci unless parameter *loci* is specified. Please
      refer to classes ``mutator`` and :class:`BaseOperator` for
      descriptions of other parameters.



class ContextMutator
--------------------

.. class:: ContextMutator

   This context-dependent mutator accepts a list of mutators and use
   one of them to mutate an allele depending on the context of the
   mutated allele.


   .. method:: ContextMutator(rates=[], loci=ALL_AVAIL, mutators=[], contexts=[], mapIn=[], mapOut=[], output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=["ind_id"], lineageMode=FROM_INFO)


      Create a mutator that choose one of the specified *mutators* to
      mutate an allele when a mutation event happens. The mutators are
      choosen according to the context of the mutated allele, which is
      specified as a list of alleles to the left and right of an
      allele (\parameter *contexts*). For example, ``contexts=[(0,0),
      (0,1), (1,1)]`` indicates which mutators should be used to
      mutate allele ``X`` in the context of ``0X0``, ``0X1``, and
      ``1X1``. A context can include more than one alleles at both
      left and right sides of a mutated allele but all contexts should
      have the same (even) number of alleles. If an allele does not
      have full context (e.g. when a locus is the first locus on a
      chromosome), unavailable alleles will be marked as -1. There
      should be a mutator for each context but an additional mutator
      can be specified as the default mutator for unmatched contexts.
      If parameters *mapIn* is specified, both mutated allele and its
      context alleles will be mapped. Most parameters, including
      *loci*, *mapIn*, *mapOut*, *rep*, and *subPops* of mutators
      specified in parameter *mutators* are ignored. This mutator by
      default applies to all loci unless parameter *loci* is
      specified. Please refer to classes ``mutator`` and
      :class:`BaseOperator` for descriptions of other parameters.



class PointMutator
------------------

.. class:: PointMutator

   A point mutator is different from all other mutators because
   mutations in this mutator do not happen randomly. Instead, it
   happens to specific loci and mutate an allele to a specific state,
   regardless of its original state. This mutator is usually used to
   introduce a mutant to a population.


   .. method:: PointMutator(loci, allele, ploidy=0, inds=[], output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=0, infoFields=["ind_id"], lineageMode=FROM_INFO)


      Create a point mutator that mutates alleles at specified *loci*
      to a given *allele* of individuals *inds*. If there are multiple
      alleles at a locus (e.g. individuals in a diploid population),
      only the first allele is mutated unless indexes of alleles are
      listed in parameter *ploidy*. This operator is by default
      applied to individuals in the first subpopulation but you can
      apply it to a different or more than one (virtual)
      subpopulations using parameter *subPops* (``AllAvail`` is also
      accepted). Please refer to class :class:`BaseOperator` for
      detailed descriptions of other parameters.



class SNPMutator
----------------

.. class:: SNPMutator

   A mutator model that assumes two alleles 0 and 1 and accepts mutation
   rate from 0 to 1, and from 1 to 0 alleles.

   .. method:: SNPMutator.SNPMutator(u=0, v=0, loci=True, mapIn=[], mapOut=[], output='', begin=0, end=-1, step=1, at=[], reps=True, subPops=ALL_AVAIL, infoFields=['ind_id'], lineageMode=115)

      Return a :class:`MatrixMutator` with proper mutate matrix for a two-allele
      mutation model using mutation rate from allele 0 to 1 (parameter ``u``)
      and from 1 to 0 (parameter ``v``)



class AcgtMutator
-----------------

.. class:: AcgtMutator

   This mutation operator assumes alleles 0, 1, 2, 3 as nucleotides ``A``,
   ``C``, ``G`` and ``T`` and use a 4 by 4 mutation rate matrix to mutate them.
   Although a general model needs 12 parameters, less parameters are needed
   for specific nucleotide mutation models (parameter ``model``). The length
   and meaning of parameter ``rate`` is model dependent.

   .. method:: AcgtMutator.AcgtMutator(rate=[], model='general', loci=True, mapIn=[], mapOut=[], output='', begin=0, end=-1, step=1, at=[], reps=True, subPops=ALL_AVAIL, infoFields=['ind_id'], lineageMode=115)

      Create a mutation model that mutates between nucleotides ``A``,
      ``C``, ``G``, and ``T`` (alleles are coded in that order as 0, 1, 2
      and 3). Currently supported models are Jukes and Cantor 1969 model
      (``JC69``), Kimura's 2-parameter model (``K80``), Felsenstein 1981
      model (``F81``), Hasgawa, Kishino and Yano 1985 model (``HKY85``),
      Tamura 1992 model (``T92``), Tamura and Nei 1993 model (``TN93``),
      Generalized time reversible model (``GTR``), and a general model
      (``general``) with 12 parameters. Please refer to the simuPOP user's
      guide for detailed information about each model.



