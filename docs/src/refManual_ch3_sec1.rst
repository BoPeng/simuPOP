Base class for all operators
============================


class BaseOperator
------------------

.. class:: BaseOperator

   Operators are objects that act on populations. They can be applied
   to populations directly using their function forms, but they are
   usually managed and applied by a simulator. In the latter case,
   operators are passed to the ``evolve`` function of a simulator, and
   are applied repeatedly during the evolution of the simulator.

   The *BaseOperator* class is the base class for all operators. It
   defines a common user interface that specifies at which
   generations, at which stage of a life cycle, to which populations
   and subpopulations an operator is applied. These are achieved by a
   common set of parameters such as ``begin``, ``end``, ``step``,
   ``at``, ``stage`` for all operators. Note that a specific operator
   does not have to honor all these parameters. For example, a
   Recombinator can only be applied during mating so it ignores the
   ``stage`` parameter.

   An operator can be applied to all or part of the generations during
   the evolution of a simulator. At the beginning of an evolution, a
   simulator is usually at the beginning of generation ``0``. If it
   evolves ``10`` generations, it evolves generations ``0``, ``1``,
   ,,,., and ``9`` (``10`` generations) and stops at the begging of
   generation ``10``. A negative generation number ``a`` has
   generation number ``10 + a``, with -1 referring to the last evolved
   generation ``9``. Note that the starting generation number of a
   simulator can be changed by its ``setGen()`` member function.

   Output from an operator is usually directed to the standard output
   (``sys.stdout``). This can be configured using a output
   specification string, which can be ``''`` for no output, ``'>'``
   standard terminal output (default), a filename prefixed by one or
   more ``'>'`` characters or a Python expression indicated by a
   leading exclamation mark (``'!expr'``). In the case of
   ``'>filename'`` (or equivalently ``'filename'``), the output from
   an operator is written to this file. However, if two operators
   write to the same file ``filename``, or if an operator writes to
   this file more than once, only the last write operation will
   succeed. In the case of ``'>>filename'``, file ``filename`` will be
   opened at the beginning of the evolution and closed at the end.
   Outputs from multiple operators are appended. ``>>>filename`` works
   similar to ``>>filename`` but ``filename``, if it already exists at
   the beginning of an evolutionary process, will not be cleared. If
   the output specification is prefixed by an exclamation mark, the
   string after the mark is considered as a Python expression. When an
   operator is applied to a population, this expression will be
   evaluated within the population's local namespace to obtain a
   population specific output specification. As an advanced feature, a
   Python function can be assigned to this parameter. Output strings
   will be sent to this function for processing. Lastly, if the output
   stream only accept a binary output (e.g. a gzip stream),
   ``WithMode(output, 'b')`` should be used to let  simuPOP convert
   string to bytes before writing to the output.


   .. method:: BaseOperator(output, begin, end, step, at, reps, subPops, infoFields)


      The following parameters can be specified by all operators.
      However, an operator can ignore some parameters and the exact
      meaning of a parameter can vary.

      output
         A string that specifies how output from an operator is written, which
         can be ``''`` (no output), ``'>'`` (standard output),
         ``'filename'`` prefixed by one or more '>', or an Python
         expression prefixed by an exclamation mark (``'!expr'``). If
         a ``file`` object, or any Python object with a ``write``
         function is provided, the output will be write to this file.
         Alternatively, a Python function or a file object (any Python
         object with a ``write`` function) can be given which will be
         called with a string of output content. A global function
         :class:`WithMode` can be used to let  simuPOP output bytes
         instead of string.

      begin
         The starting generation at which an operator will be applied. Default
         to ``0``. A negative number is interpreted as a generation
         counted from the end of an evolution (-1 being the last
         evolved generation).

      end
         The last generation at which an operator will be applied. Default to
         ``-1``, namely the last generation.

      step
         The number of generations between applicable generations. Default to
         ``1``.

      at
         A list of applicable generations. Parameters ``begin``, ``end``, and
         ``step`` will be ignored if this parameter is specified. A
         single generation number is also acceptable.

      reps
         A list of applicable replicates. A common default value ``ALL_AVAIL``
         is interpreted as all replicates in a simulator. Negative
         indexes such as ``-1`` (last replicate) is acceptable.
         ``rep=idx`` can be used as a shortcut for ``rep=[idx]``.

      subPops
         A list of applicable (virtual) subpopulations, such as ``subPops=[sp1,
         sp2, (sp2, vsp1)]``. ``subPops=[sp1]`` can be simplied as
         ``subPops=sp1``. Negative indexes are not supported. A common
         default value (``ALL_AVAIL``) of this parameter reprents all
         subpopulations of the population being aplied. Suport for
         this parameter vary from operator to operator and some
         operators do not support virtual subpopulations at all.
         Please refer to the reference manual of individual operators
         for their support for this parameter.

      infoFields
         A list of information fields that will be used by an operator. You
         usually do not need to specify this parameter because
         operators that use information fields usually have default
         values for this parameter.


   .. method:: BaseOperator.apply(pop)

      Apply an operator to population *pop* directly, without checking
      its applicability.

   .. method:: BaseOperator.clone()

      Return a cloned copy of an operator. This function is available
      to all operators.


