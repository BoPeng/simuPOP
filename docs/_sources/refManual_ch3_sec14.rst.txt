Miscellaneous operators
=======================


class NoneOp
------------

.. class:: NoneOp

   This operator does nothing when it is applied to a population. It
   is usually used as a placeholder when an operator is needed
   syntactically.


   .. method:: NoneOp(output=">", begin=0, end=0, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a :class:`NoneOp`.



class Dumper
------------

.. class:: Dumper

   This operator dumps the content of a population in a human readable
   format. Because this output format is not structured and can not be
   imported back to  simuPOP, this operator is usually used to dump a
   small population to a terminal for demonstration and debugging
   purposes.


   .. method:: Dumper(genotype=True, structure=True, ancGens=UNSPECIFIED, width=1, max=100, loci=[], output=">", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)


      Create a operator that dumps the genotype structure (if
      *structure* is ``True``) and genotype (if *genotype* is
      ``True``) to an *output* ( default to standard terminal output).
      Because a population can be large, this operator will only
      output the first 100 (parameter *max*) individuals of the
      present generation (parameter *ancGens*). All loci will be
      outputed unless parameter *loci* are used to specify a subset of
      loci. This operator by default output values of all information
      fields unless parameter *infoFields* is used to specify a subset
      of info fields to display. If a list of (virtual) subpopulations
      are specified, this operator will only output individuals in
      these outputs. Please refer to class :class:`BaseOperator` for a
      detailed explanation for common parameters such as *output* and
      *stage*.



class SavePopulation
--------------------

.. class:: SavePopulation

   An operator that save populations to specified files.


   .. method:: SavePopulation(output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create an operator that saves a population to *output* when it
      is applied to the population. This operator supports all output
      specifications (``''``, ``'filename'``, ``'filename'`` prefixed
      by one or more '>' characters, and ``'!expr'``) but output from
      different operators will always replace existing files
      (effectively ignore '>' specification). Parameter *subPops* is
      ignored. Please refer to class :class:`BaseOperator` for a
      detailed description about common operator parameters such as
      *stage* and *begin*.



class Pause
-----------

.. class:: Pause

   This operator pauses the evolution of a simulator at given
   generations or at a key stroke. When a simulator is stopped, you
   can go to a Python shell to examine the status of an evolutionary
   process, resume or stop the evolution.


   .. method:: Pause(stopOnKeyStroke=False, prompt=True, output=">", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create an operator that pause the evolution of a population when
      it is applied to this population. If *stopOnKeyStroke* is
      ``False`` (default), it will always pause a population when it
      is applied, if this parameter is set to ``True``, the operator
      will pause a population if *any* key has been pressed. If a
      specific character is set, the operator will stop when this key
      has been pressed. This allows, for example, the use of several
      pause operators to pause different populations.

      After a population has been paused, a message will be displayed
      (unless *prompt* is set to ``False``) and tells you how to
      proceed. You can press ``'s'`` to stop the evolution of this
      population, ``'S'`` to stop the evolution of all populations, or
      ``'p'`` to enter a Python shell. The current population will be
      available in this Python shell as ``"pop_X_Y"`` when ``X`` is
      generation number and ``Y`` is replicate number. The evolution
      will continue after you exit this interactive Python shell.

      .. note::

         Ctrl-C will be intercepted even if a specific character is
         specified in parameter *stopOnKeyStroke*.


class TicToc
------------

.. class:: TicToc

   This operator, when called, output the difference between current
   and the last called clock time. This can be used to estimate
   execution time of each generation. Similar information can also be
   obtained from ``turnOnDebug("DBG_PROFILE")``, but this operator has
   the advantage of measuring the duration between several generations
   by setting ``step`` parameter. As an advanced feature that mainly
   used for performance testing, this operator accepts a parameter
   *stopAfter* (seconds), and will stop the evolution of a population
   if the overall time exceeds *stopAfter*. Note that elapsed time is
   only checked when this operator is applied to a population so it
   might not be able to stop the evolution process right after
   *stopAfter* seconds. This operator can also be applied during
   mating. Note that to avoid excessive time checking, this operator
   does not always check system time accurately.


   .. method:: TicToc(output=">", stopAfter=0, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a :class:`TicToc` operator that outputs the elapsed since
      the last time it was applied, and the overall time since the
      first time this operator is applied.



