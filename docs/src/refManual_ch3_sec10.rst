Tagging operators
=================


class IdTagger
--------------

.. class:: IdTagger

   An  IdTagger gives a unique ID for each individual it is applies
   to. These ID can be used to uniquely identify an individual in a
   multi-generational population and be used to reliably reconstruct a
   Pedigree.

   To ensure uniqueness across populations, a single source of ID is
   used for this operator. individual IDs are assigned consecutively
   starting from 1. Value 1 instead of 0 is used because most software
   applications use 0 as missing values for parentship. If you would
   like to reset the sequence or start from a different number, you
   can call the ``reset(startID)`` function of any :class:`IdTagger`.

   An :class:`IdTagger` is usually used during-mating to assign ID to
   each offspring. However, if it is applied directly to a population,
   it will assign unique IDs to all individuals in this population.
   This property is usually used in the ``preOps`` parameter of
   function :meth:`Simulator.evolve` to assign initial ID to a
   population.


   .. method:: IdTagger(begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, output="", infoFields=["ind_id"])


      Create an :class:`IdTagger` that assign an unique ID for each
      individual it is applied to. The IDs are created sequentially
      and are stored in an information field specified in parameter
      *infoFields* (default to ``ind_id``). This operator is
      considered a during-mating operator but it can be used to set ID
      for all individuals of a population when it is directly applied
      to the population.


   .. method:: IdTagger.reset(startID=1)

      Reset the global individual ID number so that IdTaggers will
      start from id (default to 1) again.


class InheritTagger
-------------------

.. class:: InheritTagger

   An inheritance tagger passes values of parental information
   field(s) to the corresponding fields of offspring. If there are two
   parental values from parents of a sexual mating event, a parameter
   *mode* is used to specify how to assign offspring information
   fields.


   .. method:: InheritTagger(mode=PATERNAL, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, output="", infoFields=[])


      Creates an inheritance tagger that passes values of parental
      information fields (parameter *infoFields*) to the corresponding
      fields of offspring. If there is only one parent, values at the
      specified information fields are copied directly. If there are
      two parents, parameter *mode* specifies how to pass them to an
      offspring. More specifically,

      + ``mode=MATERNAL`` Passing the value from mother.

      + ``mode=PATERNAL`` Passing the value from father.

      + ``mode=MEAN`` Passing the average of two values.

      + ``mode=MAXIMUM`` Passing the maximum value of two values.

      + ``mode=MINIMUM`` Passing the minimum value of two values.

      + ``mode=SUMMATION`` Passing the summation of two values.

      + ``mode=MULTIPLICATION`` Passing the multiplication of two
        values.

      An :class:`RuntimeError` will be raised if any of the parents
      does not exist. This operator does not support parameter
      *subPops* and does not output any information.



class SummaryTagger
-------------------

.. class:: SummaryTagger

   A summary tagger summarize values of one or more parental
   information field to another information field of an offspring. If
   mating is sexual, two sets of parental values will be involved.


   .. method:: SummaryTagger(mode=MEAN, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, output="", infoFields=[])


      Creates a summary tagger that summarize values of one or more
      parental information field (*infoFields*[:-1]) to an offspring
      information field (*infoFields*[-1]). A parameter *mode*
      specifies how to pass summarize parental values. More
      specifically,

      + ``mode=MEAN`` Passing the average of values.

      + ``mode=MAXIMUM`` Passing the maximum value of values.

      + ``mode=Minumum`` Passing the minimum value of values.

      + ``mode=SUMMATION`` Passing the sum of values.

      + ``mode=MULTIPLICATION`` Passing the multiplication of values.

      This operator does not support parameter *subPops* and does not
      output any information.



class ParentsTagger
-------------------

.. class:: ParentsTagger

   This tagging operator records the indexes of parents (relative to
   the parental generation) of each offspring in specified information
   fields ( default to ``father_idx`` and ``mother_idx``). Only one
   information field should be specified if an asexsual mating scheme
   is used so there is one parent for each offspring. Information
   recorded by this operator is intended to be used to look up parents
   of each individual in multi-generational  Population.


   .. method:: ParentsTagger(begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, output="", infoFields=["father_idx", "mother_idx"])


      Create a parents tagger that records the indexes of parents of
      each offspring when it is applied to an offspring during-mating.
      If two information fields are specified (parameter *infoFields*,
      with default value ``['father_idx', 'mother_idx']``), they are
      used to record the indexes of each individual's father and
      mother. Value ``-1`` will be assigned if any of the parent is
      missing. If only one information field is given, it will be used
      to record the index of the first valid parent (father if both
      parents are valid). This operator ignores parameters *output*
      and *subPops*.



class OffspringTagger
---------------------

.. class:: OffspringTagger

   This tagging operator records the indexes of offspring within a
   family (sharing the same parent or parents) in specified
   information field (default to ``offspring_idx``). This tagger can
   be used to control the number of offspring during mating.


   .. method:: OffspringTagger(begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, output="", infoFields=ALL_AVAIL)


      Create an offspring tagger that records the indexes of offspring
      within a family. The index is determined by successful
      production of offspring during a mating events so the it does
      not increase the index if a previous offspring is discarded, and
      it resets index even if adjacent families share the same
      parents. This operator ignores parameters *stage*, *output*, and
      *subPops*.



class PedigreeTagger
--------------------

.. class:: PedigreeTagger

   This tagging operator records the ID of parents of each offspring
   in specified information fields (default to ``father_id`` and
   ``mother_id``). Only one information field should be specified if
   an asexsual mating scheme is used so there is one parent for each
   offspring. Information recorded by this operator is intended to be
   used to record full pedigree information of an evolutionary
   process.


   .. method:: PedigreeTagger(idField="ind_id", output="", outputFields=[], outputLoci=[], begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=["father_id", "mother_id"])


      Create a pedigree tagger that records the ID of parents of each
      offspring when it is applied to an offspring during-mating. If
      two information fields are specified (parameter *infoFields*,
      with default value ``['father_id', 'mother_id']``), they are
      used to record the ID of each individual's father and mother
      stored in the *idField* (default to ``ind_id``) field of the
      parents. Value ``-1`` will be assigned if any of the parent is
      missing. If only one information field is given, it will be used
      to record the ID of the first valid parent (father if both
      pedigree are valid).

      This operator by default does not send any output. If a valid
      output stream is given (should be in the form of
      ``'>>filename'`` so that output will be concatenated), this
      operator will output the ID of offspring, IDs of his or her
      parent(s), sex and affection status of offspring, and values at
      specified information fields (*outputFields*) and loci
      (*outputLoci*) in the format of ``off_id father_id mother_id M/F
      A/U fields genotype``. ``father_id`` or ``mother_id`` will be
      ignored if only one parent is involved. This file format can be
      loaded using function :func:`loadPedigree`.

      Because only offspring will be outputed, individuals in the top-
      most ancestral generation will not be outputed. This is usually
      not a problem because individuals who have offspring in the next
      generation will be constructed by function :func:`loadPedigree`,
      although their information fields and genotype will be missing.
      If you would like to create a file with complete pedigree
      information, you can apply this operator before evolution in the
      *initOps* parameter of functions :meth:`Population.evolve` or
      :meth:`Simulator.evolve`. This will output all individuals in
      the initial population (the top-most ancestral population after
      evolution) in the same format. Note that sex, affection status
      and genotype can be changed by other operators so this operator
      should usually be applied after all other operators are applied.



class PyTagger
--------------

.. class:: PyTagger

   A Python tagger takes some information fields from both parents,
   pass them to a user provided Python function and set the offspring
   individual fields with the return values.


   .. method:: PyTagger(func=None, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, output="", infoFields=[])


      Create a hybrid tagger that provides an user provided function
      *func* with values of specified information fields (determined
      by parameter names of this function) of parents and assign
      corresponding information fields of offspring with its return
      value. If more than one parent are available, maternal values
      are passed after paternal values. For example, if a function
      ``func(A, B)`` is passed, this operator will send two tuples
      with parental values of information fields ``'A'`` and ``'B'``
      to this function and assign its return values to fields ``'A'``
      and ``'B'`` of each offspring. The return value of this function
      should be a list, although a single value will be accepted if
      only one information field is specified. This operator ignores
      parameters *stage*, *output* and *subPops*.



