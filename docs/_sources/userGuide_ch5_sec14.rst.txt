Hybrid and Python operators
===========================


Hybrid operators
----------------

Despite the large number of built-in operators, it is obviously not possible to
implement every genetics models available. For example, although simuPOP
provides several penetrance models, a user may want to try a customized one. In
this case, one can use a *hybrid operator*.

A *hybrid operator* is an operator that calls a user-defined function when its
applied to a population. The number and meaning of input parameters and return
values vary from operator to operator. For example, a hybrid mutator sends a to-
be-mutated allele to a user-defined function and use its return value as a
mutant allele. A hybrid selector uses the return value of a user defined
function as individual fitness. Such an operator handles the routine part of the
work (e.g. scan through a chromosome and determine which allele needs to be
mutated), and leave the creative part to users. Such a mutator can be used to
implement complicated genetic models such as an asymmetric stepwise mutation
model for microsatellite markers.

**simuPOP operators use parameter names to determine which information should be
passed to a user-defined function**. For example, a hybrid quantitative trait
operator recognizes parameters ``ind``, ``geno``, ``gen`` and names of
information fields such as ``smoking``. If your model depends on genotype, you
could provide a function with parameter geno (e.g. ``func(geno)``); if your
model depends on smoking and genotype, you could provide a function with
parameters geno and smoking (e.g. ``func(geno, smoking)``); if you model depends
on individual sex, you can use a function that passes the whole individual (e.g.
``func(ind)``) so that you could check individual sex. When a hybrid operator is
applied to a population, it will check the parameter names of provided Python
function and send requested information automatically.

For example, Example :ref:`hybridOperator <hybridOperator>` defines a three-
locus heterogeneity penetrance model Risch1990 that yields positive penetrance
only when at least two disease susceptibility alleles are available. The
underlying mechanism of this operator is that for each individual, simuPOP will
collect genotype at specified loci (parameter ``loci``) and send them to
function ``myPenetrance`` and evaluate. The return values are used as the
penetrance value of the individual, which is then interpreted as the probability
that this individual will become affected.

.. _hybridOperator:

**Example**: *Use a hybrid operator*

::

   >>> import simuPOP as sim
   >>> def myPenetrance(geno):
   ...     'A three-locus heterogeneity penetrance model'
   ...     if sum(geno) < 2:
   ...         return 0
   ...     else:
   ...         return sum(geno)*0.1
   ... 
   >>> pop = sim.Population(1000, loci=[20]*3)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.8, 0.2])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.PyPenetrance(func=myPenetrance, loci=[10, 30, 50]),
   ...         sim.Stat(numOfAffected=True),
   ...         sim.PyEval(r"'%d: %d\n' % (gen, numOfAffected)")
   ...     ],
   ...     gen = 5
   ... )
   0: 97
   1: 96
   2: 78
   3: 95
   4: 80
   5

   now exiting runScriptInteractively...

`Download hybrid.py <hybrid.py>`_


Python operator :class:`PyOperator` \*
--------------------------------------

If hybrid operators are still not flexible enough, you can always resort to a
pure-Python operator :class:`PyOperator`. This operator has full access to the
evolving population (or parents and offspring when aplied during-mating), and
can therefore perform arbitrary operations.

A :class:`PyOperator` that is applied pre- or post- mating expects a function
with one or both parameters ``pop`` and ``param``, where\ ``pop`` is the
population being applied, and ``param`` is optional, depending on whether or not
a parameter is passed to the :class:`PyOperator`\ () constructor. Function
``func`` can perform arbitrary action to ``pop`` and must return ``True`` or
``False``. **The evolution of pop will be stopped if this function returns
False.** This is essentially how operator :class:`TerminateIf` works.
Alternatively, this callback function can accept ``ind`` as one of the
parameters. In this case, the function will be called for all individuals or
individuals in specified (virtual) subpopulations. **Individuals will be removed
from the populaton if this function returns False**.

Example :ref:`PyOperator <PyOperator>` defines such a function. It accepts a
cutoff value and two mutation rates as parameters. It then calculate the
frequency of allele 1 at each locus and apply a two-allele model at high
mutation rate if the frequency is lower than the cutoff and a low mutation rate
otherwise. The :func:`kAlleleMutate` function is the function form of a mutator
:class:`KAlleleMutator` (see Section :ref:`subsec_Function_form
<subsec_Function_form>` for details).

.. _PyOperator:

**Example**: *A frequency dependent mutation operator*

::

   import simuPOP as sim
   def dynaMutator(pop, param):
       '''This mutator mutates commom loci with low mutation rate and rare
       loci with high mutation rate, as an attempt to raise allele frequency
       of rare loci to an higher level.'''
       # unpack parameter
       (cutoff, mu1, mu2) = param;
       sim.stat(pop, alleleFreq=range(pop.totNumLoci()))
       for i in range(pop.totNumLoci()):
           # Get the frequency of allele 1 (disease allele)
           if pop.dvars().alleleFreq[i][1] < cutoff:
               sim.kAlleleMutate(pop, k=2, rates=mu1, loci=[i])
           else:
               sim.kAlleleMutate(pop, k=2, rates=mu2, loci=[i])
       return True
`Download PyOperator.py <PyOperator.py>`_

Example :ref:`usePyOperator <usePyOperator>` demonstrates how to use this
operator. It first initializes the population using two :class:`InitGenotype`
operators that initialize loci with different allele frequencies. It applies a
``PyOperator``\ with function ``dynaMutator`` and a tuple of parameters. Allele
frequencies at all loci are printed at generation ``0``, ``10``, ``20``, and
``30``. Note that this :class:`PyOperator` is applied at to the parental
generation so allele frequencies have to be recalculated to be used by post-
mating operator :class:`PyEval`.

.. _usePyOperator:

**Example**: *Use a PyOperator during evolution*

::

   >>> pop = sim.Population(size=10000, loci=[2, 3])
   >>> pop.evolve(
   ...     initOps=[ 
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.99, .01], loci=[0, 2, 4]),
   ...         sim.InitGenotype(freq=[.8, .2], loci=[1, 3])
   ...     ],
   ...     preOps=sim.PyOperator(func=dynaMutator, param=(.2, 1e-2, 1e-5)),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=range(5), step=10),
   ...         sim.PyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
   ...             step=10),
   ...     ],
   ...     gen = 31
   ... )                
   0.02 0.20 0.02 0.20 0.02
   0.11 0.22 0.11 0.20 0.11
   0.19 0.21 0.20 0.20 0.18
   0.21 0.21 0.22 0.21 0.21
   31

   now exiting runScriptInteractively...

`Download PyOperator.py <PyOperator.py>`_


During-mating Python operator \*
--------------------------------

A :class:`PyOperator` can also be applied during-mating. They can be used to
filter out unwanted offspring (by returning ``False`` in a user-defined
function), modify offspring, calculate statistics, or pass additional
information from parents to offspring. Depending the names of parameters of your
function, the Python operator will pass offspring (parameter ``off``), his or
her parents (parameter ``dad`` and ``mom``), the whole population (parameter
``pop``) and an optional parameter (parameter ``param``) to this function. For
example, function ``func(off)`` will accept references to an offspring, and
``func(off, mom, dad)`` will accept references to both offspring and his or her
parents.

Example :ref:`duringMatingPyOperator <duringMatingPyOperator>` demonstrates the
use of a during-mating Python operator. This operator rejects an offspring if it
has allele 1 at the first locus of the first homologous chromosome, and results
in an offspring population without such individuals.

.. _duringMatingPyOperator:

**Example**: *Use a during-mating PyOperator*

::

   >>> import simuPOP as sim
   >>> def rejectInd(off):
   ...     'reject an individual if it off.allele(0) == 1'
   ...     return off.allele(0) == 0
   ... 
   >>> pop = sim.Population(size=100, loci=1)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     matingScheme=sim.RandomMating(
   ...         ops=[
   ...             sim.MendelianGenoTransmitter(),
   ...             sim.PyOperator(func=rejectInd)
   ...         ]),
   ...     gen = 1
   ... )
   1
   >>> # You should see no individual with allele 1 at locus 0, ploidy 0.
   >>> pop.genotype()[:20]
   [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]

   now exiting runScriptInteractively...

`Download pyDuringMatingOperator.py <pyDuringMatingOperator.py>`_

:class:`PyOperator` is the most powerful operator in simuPOP and has been widely
used, for example, to calculate statistics and is not supported by the
:class:`Stat`\ () operator, to examine population property during evolution, or
prepare populations for a special mating scheme. However, because
``PyOperator``\ works in the Python interpreter, it is expected that it runs
slower than operators that are implemented at the C/C++ level. If performance
becomes an issue, you can re-implement part or all the operator in C++. Section
:ref:`subsec_Using_C++ <subsec_Using_C++>` describes how to do this.


Define your own operators \*
----------------------------

:class:`PyOperator` is a Python class so you can derive your own operator from
this operator. The tricky part is that the constructor of the derived operator
needs to call the ``__init__`` function of :class:`PyOperator` will proper
functions. This technique has been used by simuPOP in a number of occasions. For
example, the ``VarPlotter`` operator defined in ``plotter.py`` is derived from
:class:`PyOperator`. This class encapsulates several different plot class that
uses ``rpy`` to plot python expressions. One of the plotters is passed to the
func parameter of ``PyOperator.__init__`` so that it can be called when this
operator is applied.

Example :ref:`sequentialSelfing <sequentialSelfing>` rewrites the
``dynaMutator`` defined in Example :ref:`PyOperator <PyOperator>` into a derived
operator. The parameters are now passed to the constructor of ``dynaMutator``
and are saved as member variables. A member function ``mutate`` is defined and
is passed to the constructor of :class:`PyOperator`. Other than making
``dynaMutator`` look like a real simuPOP operator, this example does not show a
lot of advantage over defining a function. However, when the operator gets
complicated (as in the case for ``VarPlotter``), the object oriented
implementation will prevail.

.. _newOperator:

**Example**: *Define a new Python operator*

::

   >>> import simuPOP as sim
   >>> class dynaMutator(sim.PyOperator):
   ...     '''This mutator mutates commom loci with low mutation rate and rare
   ...     loci with high mutation rate, as an attempt to raise allele frequency
   ...     of rare loci to an higher level.'''
   ...     def __init__(self, cutoff, mu1, mu2, *args, **kwargs):
   ...         self.cutoff = cutoff
   ...         self.mu1 = mu1
   ...         self.mu2 = mu2
   ...         sim.PyOperator.__init__(self, func=self.mutate, *args, **kwargs)
   ...     #
   ...     def mutate(self, pop):
   ...         sim.stat(pop, alleleFreq=range(pop.totNumLoci()))
   ...         for i in range(pop.totNumLoci()):
   ...             # Get the frequency of allele 1 (disease allele)
   ...             if pop.dvars().alleleFreq[i][1] < self.cutoff:
   ...                 sim.kAlleleMutate(pop, k=2, rates=self.mu1, loci=[i])
   ...             else:
   ...                 sim.kAlleleMutate(pop, k=2, rates=self.mu2, loci=[i])
   ...         return True
   ... 
   >>> pop = sim.Population(size=10000, loci=[2, 3])
   >>> pop.evolve(
   ...     initOps=[ 
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.99, .01], loci=[0, 2, 4]),
   ...         sim.InitGenotype(freq=[.8, .2], loci=[1, 3])
   ...     ],
   ...     preOps=dynaMutator(cutoff=.2, mu1=1e-2, mu2=1e-5),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=range(5), step=10),
   ...         sim.PyEval(r"' '.join(['%.2f' % alleleFreq[x][1] for x in range(5)]) + '\n'",
   ...             step=10),
   ...     ],
   ...     gen = 31
   ... )          
   0.02 0.20 0.02 0.20 0.02
   0.11 0.22 0.11 0.20 0.11
   0.19 0.21 0.20 0.20 0.18
   0.21 0.21 0.22 0.21 0.21
   31

   now exiting runScriptInteractively...

`Download newOperator.py <newOperator.py>`_

New during-mating operators can be defined similarly. They are usually used to
define customized genotype transmitters. Section
:ref:`subsec_Customized_genotype_transmitter
<subsec_Customized_genotype_transmitter>` will describe this feature in detail.


