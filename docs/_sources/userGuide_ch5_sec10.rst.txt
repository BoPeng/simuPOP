Tagging operators
=================

In simuPOP, tagging refers to the action of setting various information fields
of offspring, usually using various parental information during the production
of offspring. simuPOP provides a number of tagging operators (called taggers)
for various purposes. Because tagging operators are during-mating operators,
parameter ``subPops`` can be used to tag only offspring that belong to specified
virtual subpopulation. (e.g. all male offspring)


Inheritance tagger (operator :class:`InheritTagger`)
----------------------------------------------------

An inheritance tagger passes values of parental information field(s) to the
corresponding offspring information field(s). Depending on the parameters, an
InheritTagger can

* For asexual mating schemes, pass one or more information fields from parent to
  offspring.

* Pass one or more information fields from father to offspring
  (``mode=PATERNAL``).

* Pass one or more information fields from mother to offspring
  (``mode=MATERNAL``).

* Pass the maximal, minimal, sum, multiplcation or average of values of one or
  more information fields of both parents (``mode=MAXIMUM``, ``MINIMUM``,
  ``ADDITION, MULTIPLICATION`` or ``MEAN``).

This can be used to track the spread of certain information during evolution.
For example, Example\ :ref:`InheritTagger <InheritTagger>` tags the first
individuals of ten subpopulations of size 1000. individuals in the offspring
generation inherits the maximum value of field ``x`` from his/her parents so
``x`` is inherited regardless of the sex of parents. A Stat operator is used to
calculate the number of offspring having this tag in each subpopulation. The
results show that some tagged ancestors have many offspring, and some have none.
If you run this simulation long enough, you can see that all ancestors become
the ancestor of either none or all indiviudals in a population. Note that this
simulation only considers genealogical inheritance and ancestors do not have to
pass any genotype to the last generation.

.. _InheritTagger:

**Example**: *Use an inherit tagger to track offspring of individuals*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[1000]*10, loci=1, infoFields='x')
   >>> # tag the first individual of each subpopulation.
   >>> for sp in range(pop.numSubPop()):
   ...     pop.individual(0, sp).x = 1
   ... 
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     matingScheme=sim.RandomMating(ops=[
   ...         sim.MendelianGenoTransmitter(),
   ...         sim.InheritTagger(mode=sim.MAXIMUM, infoFields='x'),
   ...     ]),
   ...     postOps=[
   ...         sim.Stat(sumOfInfo='x', vars=['sumOfInfo_sp']),
   ...         sim.PyEval(r'", ".join(["%3d" % subPop[i]["sumOfInfo"]["x"] for i in range(10)])+"\n"'),
   ...     ],
   ...     gen = 5
   ... )
     2,   1,   0,   1,   1,   2,   3,   3,   1,   1
     5,   1,   0,   1,   1,   3,   3,   5,   3,   0
     9,   2,   0,   2,   2,   7,   9,   5,  13,   0
    21,   4,   0,   2,   5,  18,  11,   9,  27,   0
    39,   5,   0,   6,   8,  36,  23,  20,  67,   0
   5

   now exiting runScriptInteractively...

`Download InheritTagger.py <InheritTagger.py>`_


Summarize parental informatin fields (operator ``SummaryTagger)``
-----------------------------------------------------------------

A :class:`SummaryTagger` summarize values of one or more parental information
fields and place the result in an offspring information field. If mating is
sexual, two sets of values will be involved. Summarization methods include
``MEAN``, ``MINIMUM``, ``MAXIMUM``, ``SUMMATION`` and ``MULTIPLICATION``. The
operator is usually used to summarize certain characteristic of parents of each
offspring. For example, a :class:`SummaryTagger` is used in  Example
:ref:`SummaryTagger <SummaryTagger>` to calculate the mean fitness of parents
during each  mating event. The results are saved in the ``avgFitness`` field of
offspring. Because allele 1 at locus 0 is under purifying selection, the allele
frequency of this allele decreases. In the mean time, fitness of parents
increases because less and less parents have this allele.

.. _SummaryTagger:

**Example**: *Using a summary tagger to calculate mean fitness of parents.*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, loci=1, infoFields=['fitness', 'avgFitness'])
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...     ],
   ...     preOps=sim.MaSelector(loci=0, wildtype=0, fitness=[1, 0.99, 0.95]),
   ...     matingScheme=sim.RandomMating(ops=[
   ...         sim.MendelianGenoTransmitter(),
   ...         sim.SummaryTagger(mode=sim.MEAN, infoFields=['fitness', 'avgFitness']),
   ...     ]),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0, meanOfInfo='avgFitness', step=10),
   ...         sim.PyEval(r"'gen %d: allele freq: %.3f, average fitness of parents: %.3f\n' % "
   ...             "(gen, alleleFreq[0][1], meanOfInfo['avgFitness'])", step=10)
   ...     ],
   ...     gen = 50,
   ... )
   gen 0: allele freq: 0.473, average fitness of parents: 0.984
   gen 10: allele freq: 0.421, average fitness of parents: 0.986
   gen 20: allele freq: 0.388, average fitness of parents: 0.988
   gen 30: allele freq: 0.288, average fitness of parents: 0.991
   gen 40: allele freq: 0.256, average fitness of parents: 0.993
   50

   now exiting runScriptInteractively...

`Download SummaryTagger.py <SummaryTagger.py>`_


Tracking parents (operator :class:`ParentsTagger`)
--------------------------------------------------

A parents tagger is used to record the indexes of parents (in the parental
population) in the information fields (default to ``father_idx``,
``mother_idx``) of their offspring. These indexes provide a way to track down an
individuals parents, offspring and consequently all relatives in a multi-
generation population. Because this operator has been extensively used in this
guide, please refer to other sections for an Example (e.g. Example
:ref:`basicInfoFields <basicInfoFields>`).

As long as parental generations do not change after the offspring generation is
created, recorded parental indexes can be used to locate parents of an
individual. However, in certain applications when parental generations change
(e.g. to draw a pedigree from a large population), or when individuals can not
be looked up easily using indexes (e.g. after individuals are saved to a file),
giving every Individual an unique ID and refer to them using ID will be a better
choice.


Tracking index of offspring within families (operator :class:`OffspringTagger`)
-------------------------------------------------------------------------------

An offspring tagger is used to record the index of offspring within each family
in an information field (default to ``offspring_idx``) of offspring. Because the
index is reset for each mating event, the index will be reset even if two
adjacent families share the same parents. In addition, this operator records the
relative index of an offspring so the index will not change if an offspring is
re-generated when the previous offspring is discarded for any reason.

Because during-mating selection operator discards offspring according their
genotypes, a mating scheme can produce families with varying sizes even if
``numOffspring`` is set to a constant number. On the other hand, if we would
like to ensure equal family size *N* in the presence of natural selection, we
will have to produce more offspring so that there can be at least *N* offspring
in each family after selection. Once *N* offspring have been generated,
excessive offspring can be discarded according to ``offspring_idx``. The
following example demonstrates such a simulation scenario:

.. _OffspringTagger:

**Example**: *Keeping constant family size in the presence of natural selection against offspring*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, loci=1, infoFields='offspring_idx')
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...     ],
   ...     matingScheme=sim.RandomMating(ops=[
   ...         sim.MendelianGenoTransmitter(),
   ...         # lethal recessive alleles
   ...         sim.MaSelector(loci=0, wildtype=0, fitness=[1, 0.90, 0.5]),
   ...         sim.OffspringTagger(),
   ...         sim.DiscardIf('offspring_idx > 4'),
   ...     ], numOffspring=10),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0, step=10),
   ...         sim.PyEval(r"'gen %d: allele freq: %.3f\n' % "
   ...             "(gen, alleleFreq[0][1])", step=10)
   ...     ],
   ...     gen = 50,
   ... )
   gen 0: allele freq: 0.445
   gen 10: allele freq: 0.187
   gen 20: allele freq: 0.089
   gen 30: allele freq: 0.087
   gen 40: allele freq: 0.059
   50

   now exiting runScriptInteractively...

`Download OffspringTagger.py <OffspringTagger.py>`_

Because families with lethal alleles produce the same number of offspring as
families without such alleles, natural selection happens within each families
and is weaker than the case when natural selection is used to all offspring.
This phenomena is generally referred to as reproductive compensation.


Assign unique IDs to individuals (operator :class:`IdTagger`)
-------------------------------------------------------------

Although it is possible to use generation number and individual indexes to
locate individuals in an evolving population, an unique I D makes it much easier
to identify individuals when migration is involved, and to analyze an
evolutionary process outside of simuPOP. An operator :class:`IdTagger` (and its
function form :func:`tagID`) is provided by simuPOP to assign an unique ID to
all individuals during evolution.

The IDs of individuals are usually stored in an information field named
``ind_id``. To ensure uniqueness across populations, a single source of ID is
used for this operator. individual IDs are assigned consecutively starting from
0. If you would like to reset the sequence or start from a different number, you
can call the ``reset(startID)`` function of any :class:`IdTagger`.

An :class:`IdTagger` is usually used during-mating to assign ID to each
offspring. However, if it is applied directly to a population, it will assign
unique IDs to all individuals in this population. This property is usually used
in the ``preOps`` parameter of function :meth:`Simulator.evolve` to assign
initial ID to a population. For example, two :class:`IdTagger` operators are
used in Example :ref:`IdTagger <IdTagger>` to assign IDs to all individuals.
Although different operators are used, different IDs are assigned to
individuals.

.. _IdTagger:

**Example**: *Assign unique IDs to individuals*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(10, infoFields='ind_id', ancGen=1)
   >>> pop.evolve(
   ...     initOps=sim.IdTagger(),
   ...     matingScheme=sim.RandomSelection(ops=[
   ...         sim.CloneGenoTransmitter(),
   ...         sim.IdTagger(),
   ...     ]),
   ...     gen = 1
   ... )
   1
   >>> print([int(ind.ind_id) for ind in pop.individuals()])
   [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
   >>> pop.useAncestralGen(1)
   >>> print([int(ind.ind_id) for ind in pop.individuals()])
   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
   >>> sim.tagID(pop) # re-assign ID
   >>> print([int(ind.ind_id) for ind in pop.individuals()])
   [21, 22, 23, 24, 25, 26, 27, 28, 29, 30]

   now exiting runScriptInteractively...

`Download IdTagger.py <IdTagger.py>`_


Tracking Pedigrees (operator :class:`PedigreeTagger`)
-----------------------------------------------------

A :class:`PedigreeTagger` is similar to a :class:`ParentsTagger` in that it
records parental information in offspring's information fields. However, instead
of indexes of parents, this operator records an unique ID of each parent to make
it easier to study and reconstruct a complete pedigree of a whole evolutionary
process. The default information fields are ``father_id`` and ``mother_id``.

By default, the :class:`PedigreeTagger` does not produce any output. However, if
a valid output string (or function) is specified, it will output the ID of
offspring and their parents, sex and affection status of offspring, and
optionally values at specified information fields (parameter ``outputFields``)
and genotype at specified loci (parameter ``outputLoci``). Because this operator
only outputs offspring, the saved file does not have detailed information of
individuals in the top-most ancestral generation. If you would like to record
complete pedigree information, you can apply :class:`PedigreeTagger` in the
``initOps`` operator of function :meth:`Simulator.evolve` or
:meth:`Population.evolve` to output information of the initial population.
Although this operator is primarily used to output pedigree information, values
at specified information fields and genotypes at specified loci could also be
outputed.

Example :ref:`PedigreeTagger <PedigreeTagger>` demonstrates how to output the
complete pedigree of an evolutionary process. Note that :class:`IdTagger` has to
be applied before :class:`PedigreeTagger` so that IDs of offspring could be
assigned before they are outputted.

.. _PedigreeTagger:

**Example**: *Output a complete pedigree of an evolutionary process*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(100, infoFields=['ind_id', 'father_id', 'mother_id'])
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.IdTagger(),
   ...         sim.PedigreeTagger(output='>>pedigree.txt'),
   ...     ],
   ...     matingScheme=sim.RandomMating(ops=[
   ...         sim.IdTagger(),
   ...         sim.PedigreeTagger(output='>>pedigree.txt'),
   ...         sim.MendelianGenoTransmitter()]
   ...     ),
   ...     gen = 100
   ... )
   100
   >>> ped = open('pedigree.txt')
   >>> lines = ped.readlines()
   >>> ped.close()
   >>> # first few lines, saved by the first PedigreeTagger
   >>> print(''.join(lines[:3]))
   1 0 0 F U
   2 0 0 F U
   3 0 0 M U

   >>> # last several lines, saved by the second PedigreeTagger
   >>> print(''.join(lines[-3:]))
   10098 9974 9915 F U
   10099 9967 9997 M U
   10100 9945 9936 M U

   >>> # load this file
   >>> ped = sim.loadPedigree('pedigree.txt')
   >>> # should have 100 ancestral generations (plus one present generation)
   >>> ped.ancestralGens()
   100

   now exiting runScriptInteractively...

`Download PedigreeTagger.py <PedigreeTagger.py>`_


A hybrid tagger (operator :class:`PyTagger`)
--------------------------------------------

A :class:`PyTagger` uses a user-defined function to pass parental information
fields to offspring. When a mating event happens, this operator collect values
of specified information fields of parents, pass them to a user-provided
function, and use the return values to set corresponding offspring information
fields. A typical usage of this operator is to set random environmental factors
that are affected by parental values. Example :ref:`PyTagger <PyTagger>`
demonstrates such an example where the location of each offspring (x, y) is
randomly assigned around the middle position of his or her parents.

.. _PyTagger:

**Example**: *Use of a hybrid tagger to pass parental information to offspring*

::

   >>> import simuPOP as sim
   >>> import random
   >>> def randomMove(x, y):
   ...     '''Pass parental information fields to offspring'''
   ...     # shift right with high concentration of alleles... 
   ...     off_x = random.normalvariate((x[0]+x[1])/2., 0.1)
   ...     off_y = random.normalvariate((y[0]+y[1])/2., 0.1)
   ...     return off_x, off_y
   ... 
   >>> pop = sim.Population(1000, loci=[1], infoFields=['x', 'y'])
   >>> pop.setVirtualSplitter(sim.GenotypeSplitter(loci=0, alleles=[[0, 0], [0,1], [1, 1]]))
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...         sim.InitInfo(random.random, infoFields=['x', 'y'])
   ...     ],
   ...     matingScheme=sim.RandomMating(ops=[
   ...         sim.MendelianGenoTransmitter(),
   ...         sim.PyTagger(func=randomMove),
   ...     ]),
   ...     postOps=[
   ...         sim.Stat(minOfInfo='x', maxOfInfo='x'),
   ...         sim.PyEval(r"'Range of x: %.2f, %.2f\n' % (minOfInfo['x'], maxOfInfo['x'])")
   ...     ],
   ...     gen = 5
   ... )
   Range of x: -0.17, 1.12
   Range of x: -0.05, 1.14
   Range of x: 0.01, 1.01
   Range of x: 0.01, 1.04
   Range of x: 0.06, 0.95
   5
   >>> 

   now exiting runScriptInteractively...

`Download PyTagger.py <PyTagger.py>`_


Tagging that involves other parental information
------------------------------------------------

If the way how parental information fields pass to their offspring is affected
by parental genotype, sex, or affection status, you could use a Python operator
(:class:`PyOperator`) during mating to explicitly obtain parental information
and set offspring information fields.

Alternatively, you can add another information field, translate needed
information to this field and pass the genotype information in the form of
information field. Operator :class:`InfoExec` could be helpful in this case.
Example :ref:`otherTagging <otherTagging>` demonstrates such an example where
the number of affected parents are recorded in an information field. Before
mating happens, a penetrance operator is used to assign affection status to
parents. The affection status is then copied to an information field affected so
that operator :class:`SummaryTagger` could be used to count the number of
affected parents. Two :class:`MaPenetrance` operators are used both before and
after mating to assign affection status to both parental and offspring
generations. This helps dividing the offspring generation into affected and
unaffected virtual subpopulations. Not surprisingly, the average number of
affected parents is larger for affected individuals than unaffected individuals.

.. _otherTagging:

**Example**: *Tagging that involves other parental information*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, loci=[1], infoFields=['aff', 'numOfAff'])
   >>> # define virtual subpopulations by affection sim.status
   >>> pop.setVirtualSplitter(sim.AffectionSplitter())
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...     ],
   ...     preOps=[
   ...         # get affection sim.status for parents
   ...         sim.MaPenetrance(loci=0, wildtype=0, penetrance=[0.1, 0.2, 0.4]),
   ...         # set 'aff' of parents
   ...         sim.InfoExec('aff = ind.affected()', exposeInd='ind'),
   ...     ],
   ...         # get number of affected parents for each offspring and store in numOfAff
   ...     matingScheme=sim.RandomMating(ops=[
   ...         sim.MendelianGenoTransmitter(),
   ...         sim.SummaryTagger(mode=sim.SUMMATION, infoFields=['aff', 'numOfAff'])]),
   ...     postOps=[
   ...         # get affection sim.status for offspring
   ...         sim.MaPenetrance(loci=0, wildtype=0, penetrance=[0.1, 0.2, 0.4]),
   ...         # calculate mean 'numOfAff' of offspring, for unaffected and affected subpopulations.
   ...         sim.Stat(meanOfInfo='numOfAff', subPops=[(0,0), (0,1)], vars=['meanOfInfo_sp']),
   ...         # print mean number of affected parents for unaffected and affected offspring.
   ...         sim.PyEval(r"'Mean number of affected parents: %.2f (unaff), %.2f (aff)\n' % "
   ...             "(subPop[(0,0)]['meanOfInfo']['numOfAff'], subPop[(0,1)]['meanOfInfo']['numOfAff'])")
   ...     ],
   ...     gen = 5
   ... )
   Mean number of affected parents: 0.41 (unaff), 0.44 (aff)
   Mean number of affected parents: 0.41 (unaff), 0.54 (aff)
   Mean number of affected parents: 0.47 (unaff), 0.55 (aff)
   Mean number of affected parents: 0.47 (unaff), 0.55 (aff)
   Mean number of affected parents: 0.42 (unaff), 0.45 (aff)
   5
   >>> 

   now exiting runScriptInteractively...

`Download otherTagging.py <otherTagging.py>`_


