Expression and Statements
=========================


class PyOutput
--------------

.. class:: PyOutput

   This operator outputs a given string when it is applied to a
   population.


   .. method:: PyOutput(msg="", output=">", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Creates a :class:`PyOutput` operator that outputs a string *msg*
      to *output* (default to standard terminal output) when it is
      applied to a population. Please refer to class
      :class:`BaseOperator` for a detailed description of common
      operator parameters such as *stage*, *begin* and *output*.



class PyEval
------------

.. class:: PyEval

   A :class:`PyEval` operator evaluates a Python expression in a
   population's local namespace when it is applied to this population.
   The result is written to an output specified by parameter *output*.


   .. method:: PyEval(expr="", stmts="", exposePop="", output=">", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=Py_False, infoFields=[])


      Create a :class:`PyEval` operator that evaluates a Python
      expression *expr* in a population's local namespaces when it is
      applied to this population. This namespace can either be the
      population's local namespace (``pop.vars()``), or namespaces
      ``subPop[sp]`` for (virtual) subpop (``pop.vars(subpop)``) in
      specified *subPops*. If Python statements *stmts* is given (a
      single or multi-line string), the statement will be executed
      before *expr*. If *exposePop* is set to an non-empty string, the
      current population will be exposed in its own local namespace as
      a variable with this name. This allows the execution of
      expressions such as ``'pop.individual(0).allele(0)'``. The
      result of *expr* will be sent to an output stream specified by
      parameter ``output``. The exposed population variable will be
      removed after *expr* is evaluated. Please refer to class
      :class:`BaseOperator` for other parameters.

      .. note::

         Although the statements and expressions are evaluated in a
         population's local namespace, they have access to a global
         namespace which is the module global namespace. It is
         therefore possible to refer to any module variable in these
         expressions. Such mixed use of local and global variables is,
         however, strongly discouraged.


class PyExec
------------

.. class:: PyExec

   This operator executes given Python statements in a population's
   local namespace when it is applied to this population.


   .. method:: PyExec(stmts="", exposePop="", output=">", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=Py_False, infoFields=[])


      Create a :class:`PyExec` operator that executes statements
      *stmts* in a population's local namespace when it is applied to
      this population. This namespace can either be the population's
      local namespace (``pop.vars()``), or namespaces ``subPop[sp]``
      for each (virtual) subpop (``pop.vars(subpop)``) in specified
      *subPops*. If *exposePop* is given, current population will be
      exposed in its local namespace as a variable named by
      *exposePop*. Although multiple statements can be executed, it is
      recommended that you use this operator to execute short
      statements and use :class:`PyOperator` for more complex once.
      Note that exposed population variables will be removed after the
      statements are executed.



class InfoEval
--------------

.. class:: InfoEval

   Unlike operator :class:`PyEval` and :class:`PyExec` that work at
   the population level, in a population's local namespace, operator
   :class:`InfoEval` works at the individual level, working with
   individual information fields. When this operator is applied to a
   population, information fields of eligible individuals are put into
   the local namespace of the population. A Python expression is then
   evaluated for each individual. The result is written to an output.


   .. method:: InfoEval(expr="", stmts="", usePopVars=False, exposeInd="", output=">", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create an operator that evaluate a Python expression *expr*
      using individual information fields and population variables as
      variables. If *exposeInd* is not empty, the individual itself
      will be exposed in the population's local namespace as a
      variable with name specified by *exposeInd*.

      A Python expression (*expr*) is evaluated for each individual.
      The results are converted to strings and are written to an
      output specified by parameter *output*. Optionally, a statement
      (or several statements separated by newline) can be executed
      before *expr* is evaluated. The evaluation of this statement may
      change the value of information fields.

      Parameter *usePopVars* is obsolete because population variables
      are always usable in such expressions.



class InfoExec
--------------

.. class:: InfoExec

   Operator :class:`InfoExec` is similar to :class:`InfoEval` in that
   it works at the individual level, using individual information
   fields as variables. This is usually used to change the value of
   information fields. For example, ``"b=a*2"`` will set the value of
   information field ``b`` to ``a*a`` for all individuals.


   .. method:: InfoExec(stmts="", usePopVars=False, exposeInd="", output="", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])


      Create an operator that executes Python statements *stmts* using
      individual information fields and population variables as
      variables. If *exposeInd* is not empty, the individual itself
      will be exposed in the population's local namespace as a
      variable with name specified by *exposeInd*.

      One or more python statements (*stmts*) are executed for each
      individual. Information fields of these individuals are then
      updated from the corresponding variables. For example, ``a=1``
      will set information field *a* of all individuals to ``1``,
      ``a=b`` will set information field *a* of all individuals to
      information field ``b`` or a population variable ``b`` if ``b``
      is not an information field but a population variable, and
      ``a=ind.sex()`` will set information field *a* of all
      individuals to its sex (needs ``exposeInd='ind'``.

      Parameter *usePopVars* is obsolete because population variables
      will always be usable.



