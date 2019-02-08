Expressions and statements
==========================


Output a Python string (operator :class:`PyOutput`)
---------------------------------------------------

Operator :class:`PyOutput` is a simple operator that prints a Python string when
it is applied to a population. It is commonly used to print the progress of a
simulation (e.g. :class:`PyOutput`\ (``'start migration\\n', at=200``)) or
output separators to beautify outputs from :class:`PyEval` outputs (e.g.
:class:`PyOutput`\ (``'\\n', rep=-1``).


Execute Python statements (operator :class:`PyExec`)
----------------------------------------------------

Operator :class:`PyExec` executes Python statements in a population's local
namespace when it is applied to that population. This operator is designed to
execute short Python statements but multiple statements separated by newline
characters are allowed.

Example :ref:`PyExec <PyExec>` uses two :class:`PyExec` operators to create and
use a variable ``traj`` in each population's local namespace. The first operator
initialize this variable as an empty list. During evolution, the frequency of
allele 1 at locus 0 is calcuated (operator :class:`Stat`) and appended to this
variable (operator :class:`PyExec`). The result is a trajectory of allele
frequencies during evolution.

.. _PyExec:

**Example**: *Execute Python statements during evolution*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(sim.Population(100, loci=1),
   ...     rep=2)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.2, 0.8]),
   ...         sim.PyExec('traj=[]')
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.PyExec('traj.append(alleleFreq[0][1])'),
   ...     ],
   ...     gen=5
   ... )
   (5, 5)
   >>> # print Trajectory
   >>> print(', '.join(['%.3f' % x for x in simu.dvars(0).traj]))
   0.775, 0.790, 0.760, 0.750, 0.750

   now exiting runScriptInteractively...

`Download PyExec.py <PyExec.py>`_


Evaluate and output Python expressions (operator :class:`PyEval`)
-----------------------------------------------------------------

Operator :class:`PyEval` evaluate a given Python expression in a population's
local namespace and output its return value. This operator has been widely used
(e.g. Example :ref:`simple_example <simple_example>`, :ref:`ancestralPop
<ancestralPop>`, :ref:`applicableGen <applicableGen>` and :ref:`output
<output>`) to output statistics of populations and report progress.

Two additional features of this operator may become handy from time to time.
First, an optional Python statements (parameter *stmts*) can be specified which
will be executed before the expression is evaluated. Second, the population
being applied can be exposed in its own namespace as a variable (parameter
*exposePop*). This makes it possible to access properties of a population other
than its variables. Example :ref:`PyEval <PyEval>` demonstrates both features.
In this example, two statements are executed to count the number of unique
parents in an offspring population and save them as variables ``numFather`` and
``numMother``. The operator outputs these two variables alone with a generation
number.

.. _PyEval:

**Example**: *Evaluate a expression and statements in a population's local namespace.*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, loci=1,
   ...     infoFields=['mother_idx', 'father_idx'])
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     matingScheme=sim.RandomMating(ops=[
   ...         sim.MendelianGenoTransmitter(),
   ...         sim.ParentsTagger(),
   ...     ]),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.PyEval(r'"gen %d, #father %d, #mother %d\n"' \
   ...             ' % (gen, numFather, numMother)',
   ...             stmts="numFather = len(set(pop.indInfo('father_idx')))\n"
   ...                 "numMother = len(set(pop.indInfo('mother_idx')))",
   ...             exposePop='pop')
   ...     ],
   ...     gen=3
   ... )
   gen 0, #father 439, #mother 433
   gen 1, #father 433, #mother 432
   gen 2, #father 449, #mother 420
   3

   now exiting runScriptInteractively...

`Download PyEval.py <PyEval.py>`_

Note that the function form of this operator (:func:`pyEval`) returns the result
of the expression rather than writting it to an output.


Expression and statement involving individual information fields (operator :class:`InfoEval` and :class:`InfoExec`) \*
----------------------------------------------------------------------------------------------------------------------

Operators :class:`PyEval` and :class:`PyExec` work at the population level,
using the local namespace of populations. Operator :class:`InfoEval` and
:class:`InfoExec`, on the contraray, work at the individual level, using
individual information fields (and population variables) as variables. In this
case, individual information fields are copied to the population namespace one
by one before expression or statements are executed for each individual.
Optionally, the individual object can be exposed to these namespace using a
user-specified name (parameter *exposeInd*). Individual information fields will
be updated if the value of these fields are changed.

Operator :class:`InfoEval` evaluates an expression and outputs its value.
Operator :class:`InfoExec` executes one or more statements and does not produce
any output. Operator :class:`InfoEval` is usually used to output individual
information fields and properties in batch mode. It is faster and sometimes
easier to use than corresponding for loop plus individual level operations. For
example

* :class:`InfoEval`\ (``r'''%.2f\\t'' % a'``) outputs the value of information
  field a for all individuals, separated by tabs.

* :class:`InfoEval`\ (``'ind.sexChar()', exposeInd='ind'``) outputs the sex of
  all individuals using an exposed individual object ``ind``.

* :class:`InfoEval`\ (``'a+b**2'``) outputs :math:`a+b^{2}` for information
  fields :math:`a` and :math:`b` for all individuals.

Example :ref:`InfoEval <InfoEval>` demonstrates the use of this operator.

.. _InfoEval:

**Example**: *Evaluate expressions using individual information fields*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(20, loci=1, infoFields='a')
   >>> pop.setVirtualSplitter(sim.InfoSplitter('a', cutoff=[3]))
   >>> sim.initGenotype(pop, freq=[0.2, 0.8])
   >>> pop.setIndInfo([random.randint(2, 5) for x in range(20)], 'a')
   >>> sim.infoEval(pop, 'a', subPops=[(0, 0)]);print(' ')
   2.02.02.02.0 
   >>> sim.infoEval(pop, 'ind.allele(0, 0)', exposeInd='ind');print(' ')
   11011111111100111111 
   >>> # use sim.population variables
   >>> pop.dvars().b = 5
   >>> sim.infoEval(pop, '"%d " % (a+b)');print(' ')
   8 9 10 8 9 10 8 9 10 10 9 7 9 7 9 7 9 7 9 8  

   now exiting runScriptInteractively...

`Download InfoEval.py <InfoEval.py>`_

Operator :class:`InfoExec` is usually used to set individual information fields.
For example

* :class:`InfoExec`\ (``'age += 1'``) increases the age of all individuals by
  one.

* :class:`InfoExec`\ (``'risk = 2 if packPerYear > 10 else 1.5'``) sets
  information field ``risk`` to ``2`` if ``packPerYear`` is greater than ``10``,
  and ``1.5`` otherwise. Note that conditional expression is only available for
  Python version 2.5 or later.

* :class:`InfoExec`\ (``'a = b*c'``) sets the value of information field ``a``
  to the product of ``b`` and ``c``.

Example :ref:`InfoExec <InfoExec>` demonstrates the use of this operator, using
its function form ``infoExec``.

.. _InfoExec:

**Example**: *Execute statements using individual information fields*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(100, loci=1, infoFields=['a', 'b', 'c'])
   >>> sim.initSex(pop)
   >>> sim.initGenotype(pop, freq=[0.2, 0.8])
   >>> sim.infoExec(pop, 'a=1')
   >>> print(pop.indInfo('a')[:10])
   (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
   >>> sim.infoExec(pop, 'b=ind.sex()', exposeInd='ind')
   >>> print(pop.indInfo('b')[:10])
   (2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0)
   >>> sim.infoExec(pop, 'c=a+b')
   >>> print(pop.indInfo('c')[:10])
   (3.0, 3.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0)
   >>> pop.dvars().d = 5
   >>> sim.infoExec(pop, 'c+=d')
   >>> print(pop.indInfo('c')[:10])
   (8.0, 8.0, 7.0, 7.0, 7.0, 7.0, 7.0, 8.0, 8.0, 8.0)
   >>> # the operator can update population variable as well
   >>> sim.infoExec(pop, 'd+=c*c')
   >>> print(pop.dvars().d)
   5835.0

   now exiting runScriptInteractively...

`Download InfoExec.py <InfoExec.py>`_

Note that a statement can also be specified for operator :class:`InfoEval`,
which will be executed before an expression is evaluated.


Using functions in external modules in simuPOP expressions and statements
-------------------------------------------------------------------------

All simuPOP expressions and statements are evaluated in a population's local
namespace, which is a dictionary with no access to external modules. If you
would like to use external modules (e.g. functions from the ``random`` module),
you will have to import them to the namespace explicitly, using something like

::

   exec('import random', pop.vars(), pop.vars())

before you evolve the population.

Example :ref:`outputByInterval <outputByInterval>` demonstrates the application
of this technique. This example imports the ``time`` module in the population's
local namespace and set ``init_time`` and ``last_time`` before evolution. During
evolution, an\ ``IfElse`` operator is used to output the status of the
simulation for every 5 seconds using expression ``time.time() - last_time > 5``.
``last_time`` is reset using the :class:`PyExec` operator. The evolution will
last 20 seconds and be terminated by the Terminator with expression
``time.time() - init_time > 20.``

.. _outputByInterval:

**Example**: *Write the status of an evolutionary process every 10 seconds*

::

   >>> import simuPOP as sim
   >>> import time
   >>> pop = sim.Population(1000, loci=10)
   >>> pop.dvars().init_time = time.time()
   >>> pop.dvars().last_time = time.time()
   >>> exec('import time', pop.vars(), pop.vars())
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.IfElse('time.time() - last_time > 5', [
   ...             sim.PyEval(r'"Gen: %d\n" % gen'),
   ...             sim.PyExec('last_time = time.time()')
   ...             ]),
   ...         sim.TerminateIf('time.time() - init_time > 20')
   ...     ]
   ... )
   Gen: 5043
   Gen: 9971
   Gen: 14997
   19925
   >>>         

   now exiting runScriptInteractively...

`Download outputByInterval.py <outputByInterval.py>`_


