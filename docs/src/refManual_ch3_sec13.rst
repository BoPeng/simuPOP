The Python operator
===================


class PyOperator
----------------

.. class:: PyOperator

   An operator that calls a user-defined function when it is applied
   to a population (pre- or post-mating) or offsprings (during-
   mating). The function can have have parameters ``pop`` when the
   operator is applied pre- or post-mating, ``pop, off, dad, mom``
   when the operator is applied during-mating. An optional parameter
   can be passed if parameter *param* is given. In the during-mating
   case, parameters ``pop``, ``dad`` and ``mom`` can be ignored if
   *offspringOnly* is set to ``True``.


   .. method:: PyOperator(func, param=None, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a pure-Python operator that calls a user-defined function
      when it is applied. If this operator is applied before or after
      mating, your function should have form ``func(pop)`` or
      ``func(pop, param)`` where ``pop`` is the population to which
      the operator is applied, ``param`` is the value specified in
      parameter *param*. ``param`` will be ignored if your function
      only accepts one parameter. Althernatively, the function should
      have form ``func(ind)`` with optional parameters ``pop`` and
      ``param``. In this case, the function will be called for all
      individuals, or individuals in subpopulations *subPops*.
      Individuals for which the function returns ``False`` will be
      removed from the population. This operator can therefore perform
      similar functions as operator :class:`DiscardIf`.

      If this operator is applied during mating, your function should
      accept parameters ``pop``, ``off`` (or ``ind``), ``dad``,
      ``mom`` and ``param`` where ``pop`` is the parental population,
      and ``off`` or ``ind``, ``dad``, and ``mom`` are offspring and
      their parents for each mating event, and ``param`` is an
      optional parameter. If *subPops* are provided, only offspring in
      specified (virtual) subpopulations are acceptable.

      This operator does not support parameters *output*, and
      *infoFields*. If certain output is needed, it should be handled
      in the user defined function *func*. Because the status of files
      used by other operators through parameter *output* is
      undetermined during evolution, they should not be open or closed
      in this Python operator.



