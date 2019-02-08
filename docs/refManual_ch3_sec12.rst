Conditional operators
=====================


class IfElse
------------

.. class:: IfElse

   This operator uses a condition, which can be a fixed condition, an
   expression or a user-defined function, to determine which operators
   to be applied when this operator is applied. A list of if-operators
   will be applied when the condition is ``True``. Otherwise a list of
   else-operators will be applied.


   .. method:: IfElse(cond, ifOps=[], elseOps=[], output=">", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a conditional operator that will apply operators *ifOps*
      if condition *cond* is met and *elseOps* otherwise. If a Python
      expression (a string) is given to parameter *cond*, the
      expression will be evalulated in each population's local
      namespace when this operator is applied. When a Python function
      is specified, it accepts parameter ``pop`` when it is applied to
      a population, and one or more parameters ``pop``, ``off``,
      ``dad`` or ``mom`` when it is applied during mating. The return
      value of this function should be ``True`` or ``False``.
      Otherwise, parameter *cond* will be treated as a fixed condition
      (converted to ``True`` or ``False``) upon which one set of
      operators is always applied. The applicability of *ifOps* and
      *elseOps* are controlled by parameters *begin*, *end*, *step*,
      *at* and *rep* of both the :class:`IfElse` operator and
      individual operators but *ifOps* and *elseOps* opeartors does
      not support negative indexes for replicate and generation
      numbers.



class TerminateIf
-----------------

.. class:: TerminateIf

   This operator evaluates an expression in a population's local
   namespace and terminate the evolution of this population, or the
   whole simulator, if the return value of this expression is
   ``True``. Termination caused by an operator will stop the execution
   of all operators after it. The generation at which the population
   is terminated will be counted in the *evolved generations* (return
   value from ``Simulator::evolve``) if termination happens after
   mating.


   .. method:: TerminateIf(condition="", stopAll=False, message="", output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create a terminator with an expression *condition*, which will
      be evalulated in a population's local namespace when the
      operator is applied to this population. If the return value of
      *condition* is ``True``, the evolution of the population will be
      terminated. If *stopAll* is set to ``True``, the evolution of
      all replicates of the simulator will be terminated. If this
      operator is allowed to write to an *output* (default to ""), the
      generation number, proceeded with an optional *message*.



class DiscardIf
---------------

.. class:: DiscardIf

   This operator discards individuals according to either an
   expression that evaluates according to individual information
   field, or a Python function that accepts individual and its
   information fields.


   .. method:: DiscardIf(cond, exposeInd="", output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create an operator that discard individuals according to an
      expression or the return value of a Python function (parameter
      *cond*). This operator can be applied to a population before or
      after mating, or to offspring during mating. If an expression is
      passed to *cond*, it will be evalulated with each individual's
      information fields (see operator :class:`InfoEval` for details).
      If *exposeInd* is non-empty, individuals will be available for
      evaluation in the expression as an variable with name spacied by
      *exposeInd*. If the expression is evaluated to be ``True``,
      individuals (if applied before or after mating) or offspring (if
      applied during mating) will be removed or discard. Otherwise the
      return value should be either ``False`` (not discard), or a
      float number between ``0`` and ``1`` as the probability that the
      individual is removed. If a function is passed to *cond*, it
      should accept paramters *ind* and *pop* or names of information
      fields when it is applied to a population (pre or post mating),
      or parameters *off*, *dad*, *mom*, *pop* (parental population),
      or names of information fields if the operator is applied during
      mating. Individuals will be discarded if this function returns
      ``True`` or at a probability if a float number between 0 and 1
      is returned. A constant expression (e.g. ``True``, ``False``,
      ``0.4``) is also acceptable, with the last example
      (``cond=0.1``) that removes 10% of individuals at randomly. This
      operator supports parameter *subPops* and will remove only
      individuals belonging to specified (virtual) subpopulations.



