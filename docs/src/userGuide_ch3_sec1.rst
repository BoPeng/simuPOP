Pythonic issues
===============


``from simuPOP import`` \* v.s. ``import simuPOP``
--------------------------------------------------

Generally speaking, it is recommended to use ``import simuPOP`` rather than
``from simuPOP import *`` to import a simuPOP module. That is to say, instead of
using

::

   from simuPOP import *
   pop = Population(size=100, loci=[5])
   simu = Simulator(pop, RandomMating())

it is recommended that you use simuPOP like

::

   import simuPOP
   pop = simuPOP.Population(size=100, loci=[5])
   simu = simuPOP.Simulator(pop, simuPOP.RandomMating())

The major problem with ``from simuPOP import *`` is that it imports all simuPOP
symbols to the global namespace and increases the likelihood of name clashes.
For example, if you import a module ``myModule`` after simuPOP, which happens to
have a variable named ``MALE``, the following code might lead to a ``TypeError``
indicating your input for parameter sex is wrong.

::

   from simuPOP import *
   from myModule import *
   pop = Population(size=100, loci=[5])
   initSex(pop, sex=[MALE, FEMALE])

It can be even worse if the definition of ``MALE`` is changed to a different
value of the same type (e.g. to ``FEMALE``) and your simulation might produce
erroranous result without a hint.

For the sake of brevity, all examples in this user's guide use ``import simuPOP
as sim`` as an alternative form of the ``import simuPOP`` style. This saves some
keystrokes by referring simuPOP functions as ``sim.Population()`` instead of
``simuPOP.Population()``. Note that simuPOP has a number of submodules, which
are not imported by default. The recommended syntax to load these modules is:

::

   # import and use submodule simuPOP.utils
   from simuPOP import utils
   utils.simulateBackwardTrajectory(N=1000, endGen=100, endFreq=0.1)


References and the ``clone()``\ member function
-----------------------------------------------

Assignment in Python only creates a new reference to an existing object. For
example,

::

   pop = Population()
   pop1 = pop

creates a reference ``pop1`` to population ``pop``. Modifying ``pop1`` will
modify ``pop`` as well and the removal of ``pop`` will invalidate ``pop1``. For
example, a reference to the first Population in a simulator is returned from
function ``func()`` in Example :ref:`lst_Reference_to_Population
<lst_Reference_to_Population>`. The subsequent use of this ``pop`` object may
crash simuPOP because the simulator ``simu`` is destroyed, along with all its
internal populations, after ``func()`` is finished, leaving ``pop`` referring to
an invalid object.

.. _lst_Reference_to_Population:

**Example**: *Reference to a population in a
simulator*

::

   def func():
       simu = Simulator(Population(10), RandomMating(), rep=5)
       # return a reference to the first Population in the simulator
       return simu.population(0)

   pop = func()
   # simuPOP will crash because pop refers to an invalid Population.
   pop.popSize()

If you would like to have an independent copy of a population, you can use the
``clone()`` member function. Example :ref:`lst_Reference_to_Population
<lst_Reference_to_Population>` would behave properly if the ``return`` statement
is replaced by

::

   return simu.population(0).clone()

although in this specific case, extracting the first population from the
simulator using the ``extract`` function

::

   return simu.extract(0)

would be more efficient.

The ``clone()`` function exists for all simuPOP classes (objects) such as
*simulator*, *mating schemes* and *operators*. simuPOP also supports the
standard Python shallow and deep copy operations so you can also make a cloned
copy of ``pop`` using the ``deepcopy`` function defined in the Python ``copy``
module

::

   import copy
   pop1 = copy.deepcopy(pop)


Zero-based indexes, absolute and relative indexes
-------------------------------------------------

**All arrays in simuPOP start at index 0**. This conforms to Python and C++
indexes. To avoid confusion, I will refer the first locus as locus zero, the
second locus as locus one; the first individual in a population as Individual
zero, and so on.

.. index::
   single: index; absolute
   single: index; relative
   single: index; absolute
   single: index; relative

Another two important concepts are the *absolute index* and *relative index* of
a locus. The former index ignores chromosome structure. For example, if there
are 5 and 7 loci on the first two chromosomes, the absolute indexes of the two
chromosomes are (0, 1, 2, 3, 4), (5, 6, 7, 8, 9, 10, 11) and the relative
indexes are (0, 1, 2, 3, 4), (0, 1, 2, 3, 4, 5, 6). Absolute indexes are more
frequently used because they avoid the trouble of having to use two numbers
(chrom, index) to refer to a locus. Two functions ``chromLocusPair(idx)`` and
``absLocusIndex(chrom,index)`` are provided to convert between these two kinds
of indexes. An individual can also be referred by its *absolute index* and
*relative index* where *relative index* is the index in its subpopulation.
Related member functions are ``subPopIndPair(idx)`` and ``absIndIndex(idx,
subPop)``. Example :ref:`absIndex <absIndex>` demonstrates the use of these
functions.

.. _absIndex:

**Example**: *Conversion between absolute and relative indexes*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[10, 20], loci=[5, 7])
   >>> print(pop.chromLocusPair(7))
   (1, 2)
   >>> print(pop.absLocusIndex(1, 1))
   6
   >>> print(pop.absIndIndex(2, 1))
   12
   >>> print(pop.subPopIndPair(25))
   (1, 15)

   now exiting runScriptInteractively...

`Download absIndex.py <absIndex.py>`_


Ranges and iterators
--------------------

Ranges in simuPOP also conform to Python ranges. That is to say, a range has the
form of ``[a,b)``\ where ``a``\ belongs to the range, and ``b``\ does not. For
example, ``pop.chromBegin(1)``\ refers to the index of the first locus on
chromosome 1 (actually exists), and ``pop.chromEnd(1)``\ refers to the index of
the last locus on chromosome 1 **plus 1**, which might or might not be a valid
index.

A number of simuPOP functions return Python iterators that can be used to
iterate through an internal array of objects. For example,
``Population.Individuals([subPop])`` returns an iterator iterates through all
individuals, or all individuals in a (virtual) subpoulation.
:meth:`Simulator.populations`\ () can be used to iterate through all populations
in a simulator. Example :ref:`iterator <iterator>` demonstrates the use of
ranges and iterators in simuPOP.

.. _iterator:

**Example**: *Ranges and iterators*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=2, loci=[5, 6])
   >>> sim.initGenotype(pop, freq=[0.2, 0.3, 0.5])
   >>> for ind in pop.individuals():
   ...     for loc in range(pop.chromBegin(1), pop.chromEnd(1)):
   ...         print(ind.allele(loc))
   ... 
   0
   2
   2
   1
   1
   1
   1
   2
   2
   2
   1
   2

   now exiting runScriptInteractively...

`Download iterator.py <iterator.py>`_


Empty, ``ALL_AVAIL`` and dynamic values for parameters ``loci``, ``reps``, ``ancGen`` and ``subPops``
-----------------------------------------------------------------------------------------------------

Parameters ``loci``, ``reps`` and ``subPops`` are widely used in simuPOP to
specify which loci, replicates, ancestral generations, or (virtual) subpulations
a function or operator is applied to. These parameter accepts a list of indexes
such as ``[1, 2]``, names such as ``['a', 'b']``, and take single form inputs
(e.g. ``loci=1`` is equivalent to ``loci=[1]``). For example,

* :class:`Recombinator`\ (``loci=[]``) recombine at no locus, and

* :class:`Recombinator`\ (``loci=1``) recombine at locus 1

* :class:`Recombinator`\ (``loci=[1,2,4]``) recombine at loci 1, 2, and 4

* :class:`Recombinator`\ (``loci=[('1', 20), ('1', 25)]``) recombine at loci
  with position ``20`` and\ ``25`` on chromosome ``1``. This usage is only
  available for parameter ``loci``.

* :class:`Recombinator`\ (``loci=['a2', 'a4']``) recombine at loci ``'a2'`` and
  ``'a4'``.

The last method is easier to understand in some cases. Moreover, when you use
loci names instead of indexes in an operator, this operator can be applied to
populations with loci at different locations. For example  ::

   MaSelector(loci='a2', fitness=[1,1.01,1.02])

will be applied to locus ``a2`` regardless the actual location of this locus in
the population to which this operator is applied.

However, in the majority of the cases, these parameters take a default value
``ALL_AVAIL`` which applies the function or operator to all available loci,
replicates or subpopulations. That is to say, :class:`Recombinator`\ () or
:class:`Recombinator`\ (``loci=ALL_AVAIL``) will recombine at all applicable
loci, which will vary from population to population. Value ``UNSPECIFIED`` is
sometimes used as default parameter of these parameters, indicating that no
value has been specified. Similarly, ``subPops=[0, 'Male']`` can be used to
refer a virtual subpopulation with name ``'Male'``, regardless its virtual
subpopulation index.

Besides ``subPops=ALL_AVAIL``, which means ``subPops=[0,1,2,3]`` for a
population with 4 subpopulations, ``ALL_AVAIL`` could also be used as
``subPops=[(ALL_AVAIL, 1)]`` to specify a specific virtual subpopulation for all
subpopulations, or ``subPops=[(1, ALL_AVAIL)]`` or even ``subPops=[(ALL_AVAIL,
ALL_AVAIL)]`` to specify all virtual subpopulations in specified or all
subpopulations. This becomes handy when you, for example, would like to list all
male individuals in a population, regardless of number of subpopulations.


User-defined functions and class :class:`WithArgs` \*
-----------------------------------------------------

Some simuPOP objects call user-defined functions to perform customized
operations. For example, a penetrance operator can call a user-defined function
with genotype at specified loci and use its return value to determine the
affection status of an individual.

simuPOP uses parameter names to determine which information should be passed to
such a function. For example, a :class:`PyOperator` will pass a reference to
each offspring to a function defined with parameter ``off`` (e.g. ``func(off)``)
and references to offspring and his/her parents to a function defined with
parameters ``off``, ``dad``, and ``mom`` (e.g. ``func(off, dad, mom)``). For
example, Example :ref:`userFunc <userFunc>` defines a function ``func(geno,
smoking)`` using parameters ``geno`` and ``smoking`` so operator
:class:`PyPenetrance` will pass genotype at specified loci and value at
information field ``smoking`` to this function.

.. _userFunc:

**Example**: *Use of user-defined Python function in simuPOP*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(1000, loci=1, infoFields='smoking')
   >>> sim.initInfo(pop, lambda:random.randint(0,1), infoFields='smoking')
   >>> sim.initGenotype(pop, freq=[0.3, 0.7])
   >>> 
   >>> # a penetrance function that depends on smoking
   >>> def func(geno, smoking):
   ...     if smoking:
   ...         return (geno[0]+geno[1])*0.4
   ...     else:
   ...         return (geno[0]+geno[1])*0.1
   ... 
   >>> sim.pyPenetrance(pop, loci=0, func=func)
   >>> sim.stat(pop, numOfAffected=True)
   >>> print(pop.dvars().numOfAffected)
   352
   >>> 

   now exiting runScriptInteractively...

`Download userFunc.py <userFunc.py>`_

However, there are circumstances that you do not know the number or names of
parameters in advance so it is difficult to define such a function. For example,
your function may use an information field with programmed name
'``off``'``+str(numOffspring)`` where ``numOffspring`` is a parameter. In this
case, you can create a wrapper function object using :class:`WithArgs`\ (``func,
args``) and list passed arguments in ``args`` (e.g. :class:`WithArgs`\ ``(func,
args=['geno', 'off' + str(numOffspring)]``). As long as simuPOP knows which
arguments to pass, your function can be defined in any format you want (e.g. use
\*args parameters). Example :ref:`WithArgs <WithArgs>` provides such an example.

.. _WithArgs:

**Example**: *Specify arguments of user-provided function using function WithArgs*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(1000, loci=1, infoFields=('x', 'y'))
   >>> sim.initInfo(pop, lambda:random.randint(0,1), infoFields=('x', 'y'))
   >>> sim.initGenotype(pop, freq=[0.3, 0.7])
   >>> 
   >>> # a penetrance function that depends on unknown information fields
   >>> def func(*fields):
   ...     return 0.4*sum(fields)
   ... 
   >>> # function WithArgs tells PyPenetrance that func accepts fields x, y so that
   >>> # it will pass values at fields x and y to func.
   >>> sim.pyPenetrance(pop, loci=0, func=sim.WithArgs(func, pop.infoFields()))
   >>> sim.stat(pop, numOfAffected=True)
   >>> print(pop.dvars().numOfAffected)
   427

   now exiting runScriptInteractively...

`Download WithArgs.py <WithArgs.py>`_


Exception handling \*
---------------------

As shown in Examples :ref:`lst_Use_of_standard_module
<lst_Use_of_standard_module>` and :ref:`lst_Use_of_optimized_module
<lst_Use_of_optimized_module>`, optimized modules raise less exceptions than
standard modules. More specifically, the standard modules check for invalid
inputs frequently and raise exceptions (e.g. out of bound loci indexes). In
constrast, the optimized modules only raise exceptions where proper values could
not be pre-determined (e.g. looking for an individual in a population with an
ID). **Only exceptions that are raised in both types of modules are documented
in the simuPOP reference manual**.

Generally speaking, **you should avoid using exceptions to direct the logic of
your script** (e.g. use a ``try ... except ...`` statement around a function to
find a valid input value). Because the optimized modules might not raise these
exceptions, such a script may crash or yield invalid results when an optimized
module is used. If you have to use such a structure, please check the reference
manual and see whether or not an exception will be raised in optimized modules.


