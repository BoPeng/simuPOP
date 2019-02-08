Introduction to operators
=========================

Operators are objects that act on populations. There are two types of operators:

* **Operators that are applied to populations**. These operators are used in the
  ``initOps``, ``preOps``, ``postOps`` and ``finalOps`` parameters of the
  ``evolve`` function. The ``initOps`` operators are applied before an
  evolutionary process, the ``preOps`` operators are applied to the parental
  population at each generation before mating, the ``postOps`` operators are
  applied to the offspring population at each generation after mating, and the
  ``finalOps`` operators are applied after an evolutionary process. Examples of
  such operators include :class:`MergeSubPops` to merge subpopulations and
  :class:`StepwiseMutator` to mutate individuals using a stepwise mutation model.

* **Operators that are applied to individuals** (offspring) during mating. These
  operators are used in the ``ops`` parameter of a mating scheme. They are usually
  used to transmit genotype or other information from parents to offspring.
  Examples of such operators include :class:`MendelianGenoTransmitter` that
  transmit parental genotype to offspring according to Mendelian laws and
  :class:`ParentsTagger` that record the indexes of parents in the parental
  population to each offspring.

Some mutators could be applied both to populations and individuals. For example,
an :class:`IdTagger` could be applied to a whole population and assign an unique
ID to all individuals, or to offspring during mating.

The following sections will introduce common features of all operators. The next
chapter will explain all simuPOP operators in detail.


Apply operators to selected replicates and (virtual) subpopulations at selected generations
-------------------------------------------------------------------------------------------

Operators are, by default, applied to all generations during an evolutionary
process. This can be changed using the ``begin``, ``end``, ``step`` and ``at``
parameters. As their names indicate, these parameters control the starting
generation (``begin``), ending generation (``end``), generations between two
applicable generations (``step``), and an explicit list of applicable
generations (``at``, a single generation number is also acceptable). Other
parameters will be ignored if ``at`` is specified. It is worth noting that, if
an evolutionary process has a pre-sepcified ending generation, negative
generations numbers are allowed. They are counted backward from the ending
generation.

For example, if a simulator starts at generation ``0``, and the ``evolve``
function has parameter ``gen=10``, the simulator will stop at the *beginning* of
generation ``10``. Generation ``-1`` refers to generation ``9``, and generation
``-2`` refers to generation ``8``, and so on. Example :ref:`applicableGen
<applicableGen>` demonstrates how to set applicable generations of an operator.
In this example, a population is initialized before evolution using an
:class:`InitGenotype` operator. allele frequency at locus ``0`` is calculated at
generation ``80``, ``90``, but not ``100`` because the evolution stops at the
beginning of generation ``100``. A :class:`PyEval` operator outputs generation
number and allele frequency at the end of generation ``80`` and ``90``. Another
:class:`PyEval` operator outputs similar information at generation ``90`` and
``99``, before and after mating. Note, however, because allele frequencies are
only calculated twice, the pre-mating allele frequency at generation ``90`` is
actually calculated at generation ``80``, and the allele frequencies display for
generation ``99`` are calculated at generation ``90``. At the end of the
evolution, the population is saved to a file using a :class:`SavePopulation`
operator.

.. _applicableGen:

**Example**: *Applicable generations of an operator.*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, loci=[20])
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.8, 0.2])
   ...     ],
   ...     preOps=[
   ...         sim.PyEval(r"'At the beginning of gen %d: allele Freq: %.2f\n' % (gen, alleleFreq[0][0])",
   ...             at = [-10, -1])
   ...     ],
   ...     matingScheme = sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0, begin=80, step=10),
   ...         sim.PyEval(r"'At the end of gen %d: allele freq: %.2f\n' % (gen, alleleFreq[0][0])",
   ...             begin=80, step=10),
   ...         sim.PyEval(r"'At the end of gen %d: allele Freq: %.2f\n' % (gen, alleleFreq[0][0])",
   ...             at = [-10, -1])
   ...     ],
   ...     finalOps=sim.SavePopulation(output='sample.pop'),
   ...     gen=100
   ... )
   At the end of gen 80: allele freq: 0.92
   At the beginning of gen 90: allele Freq: 0.92
   At the end of gen 90: allele freq: 0.93
   At the end of gen 90: allele Freq: 0.93
   At the beginning of gen 99: allele Freq: 0.93
   At the end of gen 99: allele Freq: 0.93
   100

   now exiting runScriptInteractively...

`Download applicableGen.py <applicableGen.py>`_


Applicable populations and (virtual) subpopulations
---------------------------------------------------

A simulator can evolve multiple replicates of a population simultaneously.
Different operators can be applied to different replicates of this population.
This allows side by side comparison between simulations.

Parameter ``reps`` is used to control which replicate(s) an operator can be
applied to. This parameter can be a list of replicate numbers or a single
replicate number. Negative index is allowed where ``-1`` refers to the last
replicate. This technique has been widely used to produce table-like output
where a :class:`PyOutput` outputs a newline when it is applied to the last
replicate of a simulator. Example :ref:`hybridOperator <hybridOperator>`
demonstrates how to use this ``reps`` parameter. It is worth noting that
negative indexes are *dynamic* indexes relative to number of active populations.
For example, ``rep=-1`` will refer to a previous population if the last
population has stopped evolving. Use a non-negative replicate number if this is
not intended.

.. _replicate:

**Example**: *Apply operators to a subset of populations*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(sim.Population(100, loci=[20]), 5)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.2, 0.8])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0, step=10),
   ...         sim.PyEval('gen', step=10, reps=0),
   ...         sim.PyEval(r"'\t%.2f' % alleleFreq[0][0]", step=10, reps=(0, 2, -1)),
   ...         sim.PyOutput('\n', step=10, reps=-1)
   ...     ],
   ...     gen=30,
   ... )
   0	0.23	0.22	0.29
   10	0.15	0.23	0.21
   20	0.04	0.07	0.10
   (30, 30, 30, 30, 30)

   now exiting runScriptInteractively...

`Download replicate.py <replicate.py>`_

An operator can also be applied to specified (virtual) subpopulations. For
example, an ``initializer`` can be applied to male individuals in the first
subpopulation, and everyone in the second subpopulation using parameter
``subPops=[(0,0)``, 1], if a virtual subpopulation is defined by individual sex.
Generally speaking,

* ``subPops=[]`` applies the operator to all subpopulation. This is usually the
  default value of an operator.

* ``subPops=[vsp1, vsp2,...]`` applies the operator all specified (virtual)
  subpopulations. (e.g. ``subPops=[(0,0)``, 1]).

* ``subPops=sp`` is an abbreviation for ``subPops=[sp]``. If ``sp`` is virtual,
  it has to be written as ``[sp]`` because ``subPops=(0, 1)`` is intepreted as two
  non-virtual subpopulation.

However, not all operators support this parameter, and even if they do, their
interpretations of parameter input may vary. Please refer to documentation for
individual operators in *the simuPOP reference manual* for details.


Dynamically determined loci (parameter ``loci``) \*
---------------------------------------------------

Many operators accept a parameter ``loci`` to specify the applicable loci. This
parameter can be

* ``ALL_AVAIL``: all available loci of the population to which the operator is
  applied.

* [1, 2, 4, 5]: A list of loci indexes. When the operator is applied to a
  population, it will be applied to the specified loci.

* ``[('chr1', 5), ('chr1', 10), ('chr2', 5)]``: A list of chromosome position
  pairs. That is to say, when the operator is applied to a population, it will
  find loci at specified position of specified chromosome. Here chromosome names
  are names specified by parameter ``chromNames`` of the :class:`Population`
  constructor. That is to say, the operator can be applied to all population with
  such chromosomes and loci at specified locations.

* func: A function with an optional parameter ``pop``. When the operator is
  applied to a population, it will call this function, optionally pass the
  population to be applied to this function, and use its output as indexes of
  loci.

The last usage is very interesting because it allows the determination of loci
according to population property. For example, Example :ref:`dynamicLoci
<dynamicLoci>` shows an example with a :class:`MaSelector` that is applied to
the locus with highest frequency at each generation by calling function
``mostPopular``, which calculates allele frequency and pick the locus with
highest allele frequency, This example looks silly, but the technique is very
useful in simulating the introduction of disease loci by, for example, adding
positive selection pressure to one of the chosen loci.

.. _dynamicLoci:

**Example**: *Natural selection with dynamically determined loci*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(100, loci=[10], infoFields='fitness')
   >>> 
   >>> def mostPopular(pop):
   ...     sim.stat(pop, alleleFreq=sim.ALL_AVAIL)
   ...     freq = [pop.dvars().alleleFreq[x][1] for x in range(pop.totNumLoci())]
   ...     max_freq = max(freq)
   ...     pop.dvars().selLoci = (freq.index(max_freq), max_freq)
   ...     return [freq.index(max_freq)]
   ... 
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.6, 0.4]),
   ...     ],
   ...     preOps=[
   ...         sim.MaSelector(fitness=[1, 0.9, 0.8], loci=mostPopular),
   ...         sim.PyEval(r"'gen=%d, select against %d with frequency %.2f\n' % (gen, selLoci[0], selLoci[1])"),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     gen=10,
   ... )
   gen=0, select against 6 with frequency 0.45
   gen=1, select against 7 with frequency 0.46
   gen=2, select against 2 with frequency 0.51
   gen=3, select against 2 with frequency 0.48
   gen=4, select against 2 with frequency 0.45
   gen=5, select against 9 with frequency 0.45
   gen=6, select against 3 with frequency 0.46
   gen=7, select against 9 with frequency 0.44
   gen=8, select against 7 with frequency 0.47
   gen=9, select against 3 with frequency 0.44
   10

   now exiting runScriptInteractively...

`Download dynamicLoci.py <dynamicLoci.py>`_


Write output of operators to one or more files
----------------------------------------------

All operators we have seen, except for the :class:`SavePopulation` operator in
Example :ref:`applicableGen <applicableGen>`, write their output to the standard
output, namely your terminal window. However, it would be much easier for
bookkeeping and further analysis if these output can be redirected to disk
files. Parameter ``output`` is designed for this purpose.

Parameter ``output`` can take the following values:

* ``''`` (an empty string): No output.

* ``'>'``: Write to standard output.

* ``'filename'`` or ``'>filename'``: Write the output to a file named filename.
  If multiple operators write to the same file, or if the same operator writes to
  the file file several times, only the last write operation will succeed.

* ``'>>filename'``: Append the output to a file named filename. The file will be
  opened at the beginning of ``evolve`` function and closed at the end. An
  existing file will be cleared.

* ``'>>>filename'``: This is similar to the ``'>>'`` form but the file will not
  be cleared at the beginning of the ``evolve`` function.

* ``'!expr'``: ``expr`` is considered as a Python expression that will be
  evaluated at a population's local namespace whenever an output string is needed.
  For example, ``'!''%d.txt'' % gen'`` would return ``0.txt``, ``1.txt`` etc at
  generation ``0``, ``1``, ....

* File handle of an opened file. Actually any python object with a ``write``
  function.

* A Python function that can accept a string as its only parameter
  (``func(msg)``). When an operator outputs a message, this function will be
  called with this message.

* A :class:`WithMode`\ (``output, 'b'``) object with ``output`` being the any of
  the allowed output string or function. This object tells simuPOP that the output
  is opened in binary model so that it should output bytes instead of texts to it.
  This is mostly designed for Python 3 because file objects in Python 2 accepts
  string even if they are opened in binary mode.

Because a table output such as the one in Example :ref:`hybridOperator
<hybridOperator>` is written by several operators, it is clear that all of them
need to use the ``'>>'`` output format.

The :class:`SavePopulation` operator in Example :ref:`applicableGen
<applicableGen>` write to file ``sample.pop``. This works well if there is only
one replicate but not so when the operator is applied to multiple populations.
Only the last population will be saved successfully! In this case, the
expression form of parameter ``output`` should be used.

The expression form of this parameter accepts a Python expression. Whenever a
filename is needed, this expression is evaluated against the local namespace of
the population it is applied to. Because the ``evolve`` function automatically
sets variables ``gen`` and ``rep`` in a population's local namespace, such
information can be used to produce an output string. Of course, any variable in
this namespace can be used so you are not limited to these two variable.

Example :ref:`hybridOperator <hybridOperator>` demonstrates the use of these two
parameters. In this example, a table is written to file ``LD.txt`` using
``output='>>LD.txt'``. Similar operation to ``output='R2.txt'`` fails because
only the last :math:`R^{2}` value is written to this file. The last operator
writes output for each replicate to their respective output file such as
``LD_0.txt``, using an expression that involves variable ``rep``.

.. _output:

**Example**: *Use the output and outputExpr parameters*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(sim.Population(size=1000, loci=2), rep=3)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[1, 2, 2, 1])
   ...     ],
   ...     matingScheme = sim.RandomMating(ops=sim.Recombinator(rates=0.01)),
   ...     postOps=[
   ...         sim.Stat(LD=[0, 1]),
   ...         sim.PyEval(r"'%.2f\t' % LD[0][1]", step=20, output='>>LD.txt'),
   ...         sim.PyOutput('\n', reps=-1, step=20, output='>>LD.txt'),
   ...         sim.PyEval(r"'%.2f\t' % R2[0][1]", output='R2.txt'),
   ...         sim.PyEval(r"'%.2f\t' % LD[0][1]", step=20, output="!'>>LD_%d.txt' % rep"),
   ...     ],
   ...     gen=100
   ... )
   (100, 100, 100)
   >>> print(open('LD.txt').read())
   0.25	0.24	0.24	
   0.21	0.20	0.21	
   0.16	0.15	0.17	
   0.15	0.13	0.13	
   0.11	0.10	0.13	

   >>> print(open('R2.txt').read())    # Only the last write operation succeed.
   0.20	
   >>> print(open('LD_2.txt').read())  # Each replicate writes to a different file.
   0.24	0.21	0.17	0.13	0.13	

   now exiting runScriptInteractively...

`Download output.py <output.py>`_

Example :ref:`outputFunc <outputFunc>` demonstrates an advanced usage of the
``output`` parameter. In this example, a logging object is created to write to a
logfile as well as the standard output. The ``info`` and ``debug`` functions of
this object are assigned to two operators so that their outputs can be sent to
both a logfile and to the console window. One of the advantages of using a
logging mechanism is that debugging output could be suppressed easily by
adjusting the logging level of the logging object. Note that function
``logging.info()`` automatically adds a new line to its input messages before it
writes them to an output.

.. _outputFunc:

**Example**: *Output to a Python function*

::

   >>> import simuPOP as sim
   >>> import logging
   >>> # logging to a file simulation.log, with detailed debug information
   >>> logging.basicConfig(
   ...     filename='simulation.log',
   ...     level=logging.DEBUG,
   ...     format='%(levelname)s: %(message)s',
   ...     filemode='w'
   ... )
   >>> formatter = logging.Formatter('%(message)s')
   >>> logger = logging.getLogger('')
   >>> pop = sim.Population(size=1000, loci=2)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[1, 2, 2, 1])
   ...     ],
   ...     matingScheme = sim.RandomMating(ops=sim.Recombinator(rates=0.01)),
   ...     postOps=[
   ...         sim.Stat(LD=[0, 1]),
   ...         sim.PyEval(r"'LD: %d, %.2f' % (gen, LD[0][1])", step=20,
   ...             output=logger.info),   # send LD to console and a logfile
   ...         sim.PyEval(r"'R2: %d, %.2f' % (gen, R2[0][1])", step=20,
   ...             output=logger.debug),  # send R2 only to a logfile
   ...     ],
   ...     gen=100
   ... )
   100
   >>> print(open('simulation.log').read())
   INFO: LD: 0, 0.25
   DEBUG: R2: 0, 0.97
   INFO: LD: 20, 0.20
   DEBUG: R2: 20, 0.64
   INFO: LD: 40, 0.18
   DEBUG: R2: 40, 0.51
   INFO: LD: 60, 0.12
   DEBUG: R2: 60, 0.25
   INFO: LD: 80, 0.10
   DEBUG: R2: 80, 0.17


   now exiting runScriptInteractively...

`Download outputFunc.py <outputFunc.py>`_


During-mating operators
-----------------------

All operators in Examples :ref:`applicableGen <applicableGen>`, :ref:`replicate
<replicate>` and :ref:`output <output>` are applied before or after mating.
There is, however, a hidden during-mating operator that is called by
:class:`RandomMating`\ (). This operator is called
:class:`MendelianGenoTransmitter`\ () and is responsible for transmitting
genotype from parents to offspring according to Mendel's laws. All pre-defined
mating schemes (see Section :ref:`sec_Mating_Schemes <sec_Mating_Schemes>`) use
a special kind of during-mating operator to transmit genotypes. They are called
**genotype transmitters** just to show the kind of task they perform. More
during mating operators could be specified by replacing the default operator
used in the ``ops`` parameter of a mating scheme (or an offspring generator if
you are defining your own mating scheme).

Operators used in a mating scheme honor applicability parameters ``begin``,
``step``, ``end``, ``at`` and ``reps`` although they do not support negative
population and replicate indexes. It is therefore possible to apply different
during-mating operators at different generations. For example, a
:class:`Recombinator` is used in Example :ref:`transmitter <transmitter>` to
transmit parental genotypes to offspring after generation 30 while the
:class:`MendelianGenoTransmitter` is applied before that.

.. _transmitter:

**Example**: *Genotype transmitters*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=10000, loci=2)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[1, 2, 2, 1])
   ...     ],
   ...     matingScheme = sim.RandomMating(ops=[
   ...         sim.MendelianGenoTransmitter(end=29),
   ...         sim.Recombinator(rates=0.01, begin=30),
   ...     ]),
   ...     postOps=[
   ...         sim.Stat(LD=[0, 1]),
   ...         sim.PyEval(r"'gen %d, LD: %.2f\n' % (gen, LD[0][1])", step=20)
   ...     ],
   ...     gen=100
   ... )
   gen 0, LD: 0.25
   gen 20, LD: 0.25
   gen 40, LD: 0.23
   gen 60, LD: 0.19
   gen 80, LD: 0.15
   100

   now exiting runScriptInteractively...

`Download transmitter.py <transmitter.py>`_

During-mating operators can be applied to (virtual) subpopulations using
parameter ``subPops``, which **refers to (virtual) subpopulations in the
offspring population**. Section :ref:`subsec_Pre_defined_genotype_transmitters
<subsec_Pre_defined_genotype_transmitters>` and :ref:`sec_Genotype_transmitters
<sec_Genotype_transmitters>` list all genotype transmitters, Section
:ref:`subsec_Customized_genotype_transmitter
<subsec_Customized_genotype_transmitter>` demonstrates how to define your own
genotype transmitter, Section :ref:`subsec_vspSelection <subsec_vspSelection>`
demonstrates the use of during-mating operator in virtual subpopulations.


.. _subsec_Function_form:

Function form of an operator
----------------------------

Operators are usually applied to populations through a simulator but they can
also be applied to a population directly. For example, it is possible to create
an :class:`InitGenotype` operator and apply to a population as follows:

::

   InitGenotype(freq=[.3, .2, .5]).apply(pop)

Similarly, you can apply the hybrid penetrance model defined in Example
:ref:`hybridOperator <hybridOperator>` to a population by

::

   PyPenetrance(func=myPenetrance, loci=[10, 30, 50]).apply(pop)

This usage is used so often that it deserves some simplification. Equivalent
functions are defined for most operators. For example, function ``initGenotype``
is defined for operator :class:`InitGenotype` as follows

.. _funcform:

**Example**: *The function form of operator \texttt{InitGenotype*

::

   >>> from simuPOP import InitGenotype, Population
   >>> def initGenotype(pop, *args, **kwargs):
   ...     InitGenotype(*args, **kwargs).apply(pop)
   ... 
   >>> pop = Population(1000, loci=[2,3])
   >>> initGenotype(pop, freq=[.2, .3, .5])

   now exiting runScriptInteractively...

`Download funcform.py <funcform.py>`_

These functions are called function form of operators. Using these functions,
the above two example can be written as

::

   initGenotype(pop, freq=[.3, .2, .5])

and

::

   pyPenetrance(pop, func=myPenetrance, loci=[10, 30, 50])

respectively. Note that applicability parameters such as ``begin`` and ``end``
can still be passed, but they are ignored by these functions.

Finally, it is worth noting that, if you have a function that manipulates
population, you can make it an operator by wrapping it in a :class:`PyOperator`
so that it can be called repeatedly during evolution. For example, for a
function ``myFunc`` that works on a population, you can define a wrapper
function

::

   def Func(pop):
       # call myFunc
       myFunc(pop)
       return True

which can then use it in a :class:`PyOperator` as follows:

::

   PyOperator(func=Func)

The wrapper function is not needed if myFunc returns ``True`` by itself. It can
also be simplifed to a lambda function

::

   PyOperator(func=lambda pop: myFunc(pop) is None)

if you are certain that ``myFunc`` does not return any value (return ``None``).

.. note::

   Whereas output files specified by ``'>'`` are closed immediately after they are
   written, those specified by ``'>>'`` and ``'>>>'`` are not closed after the
   operator is applied to a population. This is not a problem when operators are
   used in a simulator because :meth:`Simulator.evolve` closes all files opened by
   operators, but can cause trouble when the operator is applied directly to a
   population. For example, multiple calls to :func:`dump`\ (``pop,
   output='>>file'``) will dump pop to ``file`` repeatedly but ``file`` will not be
   closed afterward. In this case, :func:`closeOutput`\ (``'file'``) should be used
   to explicitly close the file.


