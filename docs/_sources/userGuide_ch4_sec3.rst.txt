Population
==========

.. index:: single: Population

The ``Population`` object is the most important object of simuPOP. It consists
of one or more generations of individuals, grouped by subpopulations, and a
local Python dictionary to hold arbitrary population information. This class
provides a large number of functions to access and modify population structure,
individuals and their genotypes and information fields. The following sections
explain these features in detail.


Access and change individual genotype
-------------------------------------

From a user's point of view, genotypes of all individuals in a population are
arranged sequentially. Similar to functions :meth:`Individual.genotype`\ () and
:meth:`Individual.setGenotype`\ (), genotypes of a population can be accessed in
batch using functions :meth:`Population.genotype`\ () and
:meth:`Population.setGenotype`\ (). However, because it is error prone to locate
an allele of a particular individual in this long array, these functions are
usually used to perform population-level genotype operations such as clearing
all alleles (e.g. ``pop.setGenotype(0)``) or counting the number of a particular
allele across all individuals (e.g. ``pop.genotype().count(1)``).

Another way to change alleles across the whole population is to recode existing
alleles to other numbers. This is sometimes needed if you need to change allele
states to conform with a particular mutation model, assumptions of other
software applications or genetic samples. For example, if your dataset uses 1,
2, 3, 4 for A, C, T, G alleles, and you would like to use alleles 0, 1, 2 and 3
for A, C, G, T (a convention for simuPOP when nucleotide mutation models are
involved), you can use  ::

   pop.recodeAlleles([0, 0, 1, 3, 2], alleleNames=['A', 'C', 'G', 'T'])

to convert and rename the alleles (1 allele to 0, 2 allele to 1, etc). This
operation will be applied to all subpopulations for all ancestral generations,
but can be restricted to selected loci.


Subpopulations
--------------

A simuPOP population consists of one or more subpopulations. **If a population
is not structured, it has one subpopulation that is the population itself.**
Subpopulations serve as barriers of individuals in the sense that mating only
happens between individuals in the same subpopulation. A number of functions are
provided to merge, remove, resize subpopulations, and move individuals between
subpopulations (migration).

Example :ref:`subPopName <subPopName>` demonstrates how to use some of the
subpopulation related functions.

.. _subPop:

**Example**: *Manipulation of subpopulations*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[3, 4, 5], ploidy=1, loci=1, infoFields='x')
   >>> # individual 0, 1, 2, ... will have an allele 0, 1, 2, ...
   >>> pop.setGenotype(range(pop.popSize()))
   >>> #
   >>> pop.subPopSize(1)
   4
   >>> # merge subpopulations
   >>> pop.mergeSubPops([1, 2])
   1
   >>> # split subpopulations
   >>> pop.splitSubPop(1, [2, 7])
   (1, 2)
   >>> pop.subPopSizes()
   (3, 2, 7)
   >>> # remove subpopulations
   >>> pop.removeSubPops(1)
   >>> pop.subPopSizes()
   (3, 7)

   now exiting runScriptInteractively...

`Download subPop.py <subPop.py>`_

Some population operations change the IDs of subpopulations. For example, if a
population has three subpopulations 0, 1, and 2, and subpopulation 1 is split
into two subpouplations, subpopulation 2 will become subpopulation 3. Tracking
the ID of a subpopulation can be problematic, especially when conditional or
random subpopulation operations are involved. In this case, you can specify
names to subpopulations. These names will follow their associated subpopulations
during population operations so you can identify the ID of a subpopulation by
its name. Note that simuPOP allows duplicate subpopulation names.

.. _subPopName:

**Example**: *Use of subpopulation names*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[3, 4, 5], subPopNames=['x', 'y', 'z'])
   >>> pop.removeSubPops([1])
   >>> pop.subPopNames()
   ('x', 'z')
   >>> pop.subPopByName('z')
   1
   >>> pop.splitSubPop(1, [2, 3])
   (1, 2)
   >>> pop.subPopNames()
   ('x', 'z', 'z')
   >>> pop.setSubPopName('z-1', 1)
   >>> pop.subPopNames()
   ('x', 'z-1', 'z')
   >>> pop.subPopByName('z')
   2

   now exiting runScriptInteractively...

`Download subPopName.py <subPopName.py>`_


Virtual subpopulations and virtual splitters \*
-----------------------------------------------

simuPOP subpopulations can be further divided into virtual subpopulations (VSP),
which are groups of individuals who share certain properties. For example, all
male individuals, all unaffected individuals, all individuals with information
field age > 20, all individuals with genotype 0, 0 at a given locus, can form
VSPs. VSPs do not have to add up to the whole subpopulation, nor do they have to
be non-overlapping. Unlike subpopulations that have strict boundaries, VSPs
change easily with the changes of individual properties.

VSPs are defined by virtual splitters. **It is a definition for groups of
individuals in each subpopulation.** A splitter defines the same number of VSPs
in all subpopulations, although sizes of these VSPs vary across subpopulations
due to subpopulation differences. For example, a :class:`SexSplitter`\ ()
defines two VSPs, the first with all male individuals and the second with all
female individuals, and a :class:`InfoSplitter`\ (``field='x', values=[1, 2,
4]``) defines three VSPs whose members have values ``1``, ``2`` and ``4`` at
information field ``x``, respectively. This splitter also allows the use of
cutoff values and ranges to define VSPs. If different types of VSPs are needed,
a combined splitter can be used to combine VSPs defined by several splitters.

A VSP is represented by a ``[sp, vsp]`` pair where ``sp`` and ``vsp`` can be
subpopulation indexes or names. Its name and size can be obtained using
functions ``subPopName()`` and ``subPopSize()``. Example :ref:`virtualSplitter
<virtualSplitter>` demonstrates how to apply virtual splitters to a population,
and how to check VSP names and sizes.

.. _virtualSplitter:

**Example**: *Define virtual subpopulations in a population*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(size=[200, 400], loci=[30], infoFields='x')
   >>> # assign random information fields
   >>> sim.initSex(pop)
   >>> sim.initInfo(pop, lambda: random.randint(0, 3), infoFields='x')
   >>> # define a virtual splitter by sex
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> pop.numVirtualSubPop()    # Number of defined VSPs
   2
   >>> pop.subPopName([0, 0])    # Each VSP has a name
   'Male'
   >>> pop.subPopSize([0, 1])    # Size of VSP 1 in subpopulation 0
   109
   >>> pop.subPopSize([0, 'Female'])    # Refer to vsp by its name
   109
   >>> # define a virtual splitter by information field 'x'
   >>> pop.setVirtualSplitter(sim.InfoSplitter(field='x', values=[0, 1, 2, 3]))
   >>> pop.numVirtualSubPop()    # Number of defined VSPs
   4
   >>> pop.subPopName([0, 0])    # Each VSP has a name
   'x = 0'
   >>> pop.subPopSize([0, 0])    # Size of VSP 0 in subpopulation 0
   46
   >>> pop.subPopSize([1, 0])    # Size of VSP 0 in subpopulation 1
   92

   now exiting runScriptInteractively...

`Download virtualSplitter.py <virtualSplitter.py>`_

VSP provides an easy way to access groups of individuals in a subpopulation and
allows finer control of an evolutionary process. For example, mating schemes can
be applied to VSPs which makes it possible to apply different mating schemes to,
for example, individuals with different ages. By applying migration, mutation
etc to VSPs, it is easy to implement advanced features such as sex-biased
migrations, different mutation rates for individuals at different stages of a
disease. Example :ref:`virtualSubPop <virtualSubPop>` demonstrates how to
initialize genotype and information fields to individuals in male and female
VSPs.

.. _virtualSubPop:

**Example**: *Applications of virtual subpopulations*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(10, loci=[2, 3], infoFields='Sex')
   >>> sim.initSex(pop)
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> # initialize male and females with different genotypes. 
   >>> sim.initGenotype(pop, genotype=[0]*5, subPops=[(0, 0)])
   >>> sim.initGenotype(pop, genotype=[1]*5, subPops=[(0, 1)])
   >>> # set Sex information field to 0 for all males, and 1 for all females
   >>> pop.setIndInfo([sim.MALE], 'Sex', [0, 0])
   >>> pop.setIndInfo([sim.FEMALE], 'Sex', [0, 1])
   >>> # Print individual genotypes, followed by values at information field Sex
   >>> sim.dump(pop, structure=False)
   SubPopulation 0 (), 10 Individuals:
      0: FU 11 111 | 11 111 |  2
      1: FU 11 111 | 11 111 |  2
      2: MU 00 000 | 00 000 |  1
      3: MU 00 000 | 00 000 |  1
      4: MU 00 000 | 00 000 |  1
      5: MU 00 000 | 00 000 |  1
      6: MU 00 000 | 00 000 |  1
      7: FU 11 111 | 11 111 |  2
      8: FU 11 111 | 11 111 |  2
      9: FU 11 111 | 11 111 |  2


   now exiting runScriptInteractively...

`Download virtualSubPop.py <virtualSubPop.py>`_


Advanced virtual subpopulation splitters \*\*
---------------------------------------------

simuPOP provides a number of virtual splitters that can define VSPs using
specified properties. For example, :class:`InfoSplitter`\ (``field='a',
values=[1,2,3]``) defines three VSPs whose individuals have values ``1``, ``2``,
and ``3`` at information field ``a``, respectively; :class:`SexSplitter`\ ()
defines two VSPs of male and female individuals, respectively; and
:class:`RangeSplitter`\ (``ranges=[[0, 2000], [2000, 5000]]``) defines two VSPs
using two blocks of individuals.

A :class:`CombinedSplitter` can be used if your simulation needs more than one
sets of VSPs. For example, you may want to split your subpopulations both by sex
and by affection status. In this case, you can define a combined splitter using
::

   CombinedSplitter(splitters=[SexSplitter(), AffectionSplitter()])

This splitter simply stacks VSPs defined in :class:`AffectionSplitter`\ () after
:class:`SexSplitter`\ () so that unaffected and affected VSPs are now VSPs 2 and
3 (0 and 1 are used for male and female VSPs).

There are also scenarios when you would like to define finer VSPs with
individuals belonging to more than one VSPs. For example, you may want to have a
look of frequencies of certain alleles in affected male vs affected females, or
count the number of males and females with certain value at an information
field. In this case, a :class:`ProductSplitter` can be used to define VSPs using
interactions of several VSPs. For example,  ::

   ProductSplitter(splitters=[SexSplitter(), AffectionSplitter()])

defines 4 subpopulations by splitting VSPs defined by :class:`SexSplitter`\ ()
with affection status. These four VSPs will then have unaffected male, affected
male, unaffected female and affected female individuals, respectively.

If you consider :class:`ProductSplitter` as an intersection splitter that
defines new VSPs as intersections of existing VSPs, you may wonder how to define
unions of VSPs. For example, you can make a case where you want to consider
Individuals with information field a < 0 or a > 100 together. A regular
:class:`InfoSplitter`\ (``field='a', cutoff=[0, 100]``) cannot do that because
it defines three VSPs with :math:`a<0`, :math:`0\leq a<100` and
:math:`a\geq100`, respectively. The trick here is to use parameter ``vspMap`` of
a :class:`CombinedSplitter`. If this parameter is defined, multiple VSPs could
be groups or reordered to define a new set of VSPs. For example,   ::

   CombinedSplitter(splitters=[InfoSplitter(field='a', cutoff=[0, 100])], vspMap=[[0,2], 1])

defines two VSPs using VSPs 0 and 2, and VSP 1 defined by the
:class:`InfoSplitter` so that the first VSP contains individuals with
:math:`a<0` or :math:`a\geq100`.

Example :ref:`advancedVSP <advancedVSP>` demonstrates some advanced usages of
virtual splitters.

.. _advancedVSP:

**Example**: *Advanced virtual subpopulation usages.*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(size=[2000, 4000], loci=[30], infoFields='x')
   >>> # assign random information fields
   >>> sim.initSex(pop)
   >>> sim.initInfo(pop, lambda: random.randint(0, 3), infoFields='x')
   >>> #
   >>> # 1, use a combined splitter
   >>> pop.setVirtualSplitter(sim.CombinedSplitter(splitters = [
   ...     sim.SexSplitter(),
   ...     sim.InfoSplitter(field='x', values=[0, 1, 2, 3])
   ... ]))
   >>> pop.numVirtualSubPop()    # Number of defined VSPs
   6
   >>> pop.subPopName([0, 0])    # Each VSP has a name
   'Male'
   >>> pop.subPopSize([0, 0])    # sim.MALE
   1011
   >>> pop.subPopSize([1, 4])    # individuals in sp 1 with value 2 at field x
   1048
   >>> #
   >>> # use a product splitter that defines additional VSPs by sex and info
   >>> pop.setVirtualSplitter(sim.ProductSplitter(splitters = [
   ...     sim.SexSplitter(names=['M', 'F']),  # give a new set of names
   ...     sim.InfoSplitter(field='x', values=[0, 1, 2, 3])
   ... ]))
   >>> pop.numVirtualSubPop()    # Number of defined VSPs
   8
   >>> pop.subPopName([0, 0])    # Each VSP has a name
   'M, x = 0'
   >>> pop.subPopSize([0, 0])    # sim.MALE with value 1 in sp 0
   240
   >>> pop.subPopSize([1, 5])    # sim.FEMALE with value 1 in sp 1
   453
   >>> #
   >>> # use a combined splitter to join VSPs defined by a
   >>> # product splitter
   >>> pop.setVirtualSplitter(sim.CombinedSplitter([
   ...     sim.ProductSplitter([
   ...         sim.SexSplitter(),
   ...         sim.InfoSplitter(field='x', values=[0, 1, 2, 3])])],
   ...     vspMap = [[0,1,2], [4,5,6], [7]],
   ...     names = ['Male x<=3', 'Female x<=3', 'Female x=4']))
   >>> pop.numVirtualSubPop()    # Number of defined VSPs
   3
   >>> pop.subPopName([0, 0])    # Each VSP has a name
   'Male x<=3'
   >>> pop.subPopSize([0, 0])    # sim.MALE with value 0, 1, 2 at field x
   770
   >>> pop.subPopSize([1, 1])    # sim.FEMALE with value 0, 1 or 2 at field x
   1493

   now exiting runScriptInteractively...

`Download advancedVSP.py <advancedVSP.py>`_


.. _subsec_Individuals:

Access individuals and their properties
---------------------------------------

There are many ways to access individuals of a population. For example, function
``Population.Individual(idx)`` returns a reference to the ``idx``\ -th
individual in a population. An optional parameter ``subPop`` can be specified to
return the ``idx``\ -th individual in the ``subPop``\ -th subpopulation.

If you would like to access a group of individuals, either from a whole
population, a subpopulation, or from a virtual subpopulation,
:meth:`Population.individuals`\ (``[subPop]``) is easier to use. This function
returns a Python iterator that can be used to iterate through individuals. An
advantage of this function is that ``subPop``\ can be a virtual subpopulation
which makes it easy to iterate through Individuals with certain properties (such
as all male Individuals). If you would like to iterate through multiple virtual
subpopulations in one or more ancestral generations, you can use another
function :meth:`Population.allIndividuals`\ (``subPops, ancGens``).

If more than one generations are stored in a population, function
``ancestor(idx, [subPop], gen)`` can be used to access Individual from an
ancestral generation (see Section :ref:`subsec_Ancestral_populations
<subsec_Ancestral_populations>` for details). Because there is no group access
function for ancestors, it may be more convenient to use ``useAncestralGen`` to
make an *ancestral* generation the *current* generation, and use
``Population.Individuals``. Note that ancestor() function can always access
individuals at a certain generation, regardless which generation the current
generation is. Example :ref:`batchAccess <batchAccess>` demonstrates how to use
all these Individual-access functions.

If an unique ID is assigned to all individuals in a population, you can look up
individuals from their IDs using function :meth:`Population.indByID`\ (). The
information field to save individual ID is usually ``ind_id`` and you can use
operator :class:`IdTagger` and its function form :func:`tagID` to set this
field. Note that this function can be used to look up individuals in the present
and all ancestral generations, although a parameter (``ancGen``) can be used to
limit search to a specific generation if you know in advance which generation
the individual locates.

.. _accessIndividual:

**Example**: *Access individuals of a population*

::

   >>> import simuPOP as sim
   >>> # create a sim.population with two generations. The current generation has values
   >>> # 0-9 at information field x, the parental generation has values 10-19.
   >>> pop = sim.Population(size=[5, 5], loci=[2, 3], infoFields='x', ancGen=1)
   >>> pop.setIndInfo(range(10, 20), 'x')
   >>> pop1 = pop.clone()
   >>> pop1.setIndInfo(range(10), 'x')
   >>> pop.push(pop1)
   >>> #
   >>> ind = pop.individual(5)       # using absolute index
   >>> ind.x
   5.0
   >>> ind.x       # the same as ind.x
   5.0
   >>> # use a for loop, and relative index
   >>> for idx in range(pop.subPopSize(1)):
   ...     print(pop.individual(idx, 1).x)
   ... 
   5.0
   6.0
   7.0
   8.0
   9.0
   >>> # It is usually easier to use an iterator
   >>> for ind in pop.individuals(1):
   ...     print(ind.x)
   ... 
   5.0
   6.0
   7.0
   8.0
   9.0
   >>> # Access individuals in VSPs
   >>> pop.setVirtualSplitter(sim.InfoSplitter(cutoff=[3, 7, 17], field='x'))
   >>> for ind in pop.individuals([1, 1]):
   ...     print(ind.x)
   ... 
   5.0
   6.0
   >>> # Access all individuals in all ancestral generations
   >>> print([ind.x for ind in pop.allIndividuals()])
   [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0]
   >>> # or only specified subpopulations or ancestral generations
   >>> print([ind.x for ind in pop.allIndividuals(subPops=[(0,2), (1,3)], ancGens=1)])
   [10.0, 11.0, 12.0, 13.0, 14.0, 17.0, 18.0, 19.0]
   >>> 
   >>> # Access individuals in ancetral generations
   >>> pop.ancestor(5, 1).x        # absolute index
   15.0
   >>> pop.ancestor(0, 1, 1).x     # relative index
   15.0
   >>> # Or make ancestral generation the current generation and use 'individual'
   >>> pop.useAncestralGen(1)
   >>> pop.individual(5).x         # absolute index
   15.0
   >>> pop.individual(0, 1).x      # relative index
   15.0
   >>> # 'ancestor' can still access the 'present' (generation 0) generation
   >>> pop.ancestor(5, 0).x
   5.0
   >>> # access individual by ID
   >>> pop.addInfoFields('ind_id')
   >>> sim.tagID(pop)
   >>> [int(ind.ind_id) for ind in pop.individuals()]
   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
   >>> # access individual by ID. Note that individual 12 is in the parental generation
   >>> pop.indByID(12).x
   1.0

   now exiting runScriptInteractively...

`Download accessIndividual.py <accessIndividual.py>`_

Although it is easy to access individuals in a population, it is often more
efficient to access genotypes and information fields in batch mode. For example,
functions ``genotype()`` and\ ``setGenotype()`` can read/write genotype of all
individuals in a population or (virtual) subpopulation, functions ``indInfo()``
and ``setIndInfo()`` can read/write certain information fields in a population
or (virtual) subpopulation. The write functions work in a circular manner in the
sense that provided values are reused if they are not enough to fill all
genotypes or information fields. Example :ref:`batchAccess <batchAccess>`
demonstrates the use of such functions.

.. _batchAccess:

**Example**: *Access Individual properties in batch mode*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(size=[4, 6], loci=2, infoFields='x')
   >>> pop.setIndInfo([random.randint(0, 10) for x in range(10)], 'x')
   >>> pop.indInfo('x')
   (7.0, 5.0, 8.0, 10.0, 7.0, 0.0, 8.0, 4.0, 4.0, 10.0)
   >>> pop.setGenotype([0, 1, 2, 3], 0)
   >>> pop.genotype(0)
   [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
   >>> pop.setVirtualSplitter(sim.InfoSplitter(cutoff=[3], field='x'))
   >>> pop.setGenotype([0])    # clear all values
   >>> pop.setGenotype([5, 6, 7], [1, 1])
   >>> pop.indInfo('x', 1)
   (7.0, 0.0, 8.0, 4.0, 4.0, 10.0)
   >>> pop.genotype(1)
   [5, 6, 7, 5, 0, 0, 0, 0, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6]

   now exiting runScriptInteractively...

`Download batchAccess.py <batchAccess.py>`_


.. _sec_Information_fields:

Attach arbitrary auxillary information using information fields
---------------------------------------------------------------

Information fields are usually set during population creation, using the
``infoFields`` parameter of the population constructor. It can also be set or
added using functions ``setInfoFields, addInfoField``\ and ``addInfoFields``.
Example :ref:`popInfo <popInfo>` demonstrates how to read and write information
fields from an individual, or from a population in batch mode. Note that
functions :meth:`Population.indInfo` and :meth:`Population.setIndInfo` can be
applied to (virtual) subpopulation using a optional parameter subPop.

.. _popInfo:

**Example**: *Add and use of information fields in a population*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(10)
   >>> pop.setInfoFields(['a', 'b'])
   >>> pop.addInfoFields('c')
   >>> pop.addInfoFields(['d', 'e'])
   >>> pop.infoFields()
   ('a', 'b', 'c', 'd', 'e')
   >>> #
   >>> # information fields can be accessed in batch mode
   >>> pop.setIndInfo([1], 'c')
   >>> # as well as individually.
   >>> for ind in pop.individuals():
   ...     ind.e = ind.c + 1
   ... 
   >>> print(pop.indInfo('e'))
   (2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)

   now exiting runScriptInteractively...

`Download popInfo.py <popInfo.py>`_


.. _subsec_Ancestral_populations:

Keep track of ancestral generations
-----------------------------------

A simuPOP population usually holds individuals in one generation. During
evolution, an offspring generation will replace the parental generation and
become the present generation (population), after it is populated from a
parental population. The parental generation is discarded.

This is usually enough when only the present generation is of interest. However,
parental generations can provide useful information on how genotype and other
information are passed from parental to offspring generations. simuPOP provides
a mechanism to store and access arbitrary number of ancestral generations in a
population object. Applications of this feature include pedigree tracking,
reconstruction, and pedigree ascertainments.

A parameter ``ancGen`` is used to specify how many generations a population
object *can* store (which is usually called the *ancestral depth* of a
population). This parameter is default to ``0``, meaning keeping no ancestral
population. You can specify a positive number ``n`` to store n most recent
generations; or -\ ``1`` to store all generations. Of course, storing all
generations during an evolutionary process is likely to exhaust the RAM of your
computer quickly.

Several member functions can be used to manipulate ancestral generations:

* ``ancestralGens()``\ returns the number of ancestral generations stored in a
  population.

* ``setAncestralDepth(depth)`` resets the number of generations a population can
  store.

* ``push(pop)`` will push population ``pop`` into the current population.
  ``pop`` will become the current generation, and the current generation will
  either be removed (if ancGen == 0), or become the parental generation of pop.
  The greatest ancestral generation may be removed. This function is rarely used
  because populations with ancestral generations are usually created during an
  evolutionary process.

* ``useAncestralGen(idx)`` set the present generation to ``idx`` generation.
  ``idx`` ``= 1`` for the parental generation, ``2`` for grand-parental, ..., and
  ``0`` for the present generation. This is useful because most population
  functions act on the *present* generation. You should always call
  ``setAncestralPop(0)`` after you examined the ancestral generations.

If a population has several ancestral generations, they are referred by their
indexes 0 (the latest generation), 1 (parental generation), ... and :math:`k`
(top-most ancestral generation) where :math:`k` equals to ``ancestralGens()``.
In many cases, you can retrieve the properties of ancestral generations
directly, using functions such as

* ``popSize(ancGen=-1), subPopSizes(ancGen=-1), subPopSize(subPop, ancGen=-1)``:
  population and subpopulation sizes of ancestral generation ``ancGen``.

* ``ancestor(index, ancGen)``: Get a reference to the ``index`` individual of
  ancestral generation ``ancGen``.

However, most population member functions work at the current generation so you
will need to switch to an ancestral generation using function
``useAncestralGen()`` if you would like to manipulate an ancestral generation.
For example, you can remove the second subpopulation of the parental generation
using functions::

   pop.useAncestralGen(1)
   pop.removeSubPops(1)

A typical use of ancestral generations is demonstrated in example :ref:`extract
<extract>`. In this example, a population is created and is initialized with
allele frequency 0.5. Its ancestral depth is set to 2 at the beginning of
generation 18 so that it can hold parental generations at generation 18 and 19.
The allele frequency at each generation is calculated and displayed, both during
evolution using a :class:`Stat` operator, and after evolution using the function
form this operator. Note that setting the ancestral depth at the end of an
evolutionary process is a common practice because we are usually only interested
in the last few generations.

.. _ancestralPop:

**Example**: *Ancestral populations*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(500, loci=1, ancGen=2)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     matingScheme = sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0, begin=-3),
   ...         sim.PyEval(r"'%.3f\n' % alleleFreq[0][0]", begin=-3)
   ...     ],
   ...     gen = 20
   ... )
   0.495
   0.510
   0.506
   20
   >>> # information
   >>> pop.ancestralGens()
   2
   >>> pop.popSize(ancGen=1)
   500
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> # number of males in the current and parental generation
   >>> pop.subPopSize((0,0)), pop.subPopSize((0,0), ancGen=1)
   (254, 249)
   >>> # start from current generation
   >>> for i in range(pop.ancestralGens(), -1, -1):
   ...   pop.useAncestralGen(i)
   ...   sim.stat(pop, alleleFreq=0)
   ...   print('%d   %.3f' % (i, pop.dvars().alleleFreq[0][0]))
   ... 
   2   0.495
   1   0.510
   0   0.506
   >>> # restore to the current generation  
   >>> pop.useAncestralGen(0)  

   now exiting runScriptInteractively...

`Download ancestralPop.py <ancestralPop.py>`_


Change genotypic structure of a population
------------------------------------------

Several functions are provided to remove, add empty loci or chromosomes, and to
merge loci or chromosomes from another population. They can be used to trim
unneeded loci, expand existing population or merge two populations. Example
:ref:`extract <extract>` demonstrates how to use these populations. Note that
function :meth:`Population.addLociFrom` by default merges chromosomes one by one
according to chromosome index. If ``byName`` is set to True, it will try to
match chromosomes by name and merge them. This example also demonstrates the use
of ``DBG_WARNING`` flag, which will trigger a warning message when chromosomes
with different names are merged.

.. _addRemoveLoci:

**Example**: *Add and remove loci and chromosomes*

::

   >>> import simuOpt
   >>> simuOpt.setOptions(debug='DBG_WARNING')
   >>> import simuPOP as sim
   Turn on debug 'DBG_WARNING'
   >>> pop = sim.Population(10, loci=3, chromNames=['chr1'])
   >>> # 1 1 1, 
   >>> pop.setGenotype([1])
   >>> # 1 1 1, 0 0 0
   >>> pop.addChrom(lociPos=[0.5, 1, 2], lociNames=['rs1', 'rs2', 'rs3'],
   ...     chromName='chr2')
   >>> pop1 = sim.Population(10, loci=3, chromNames=['chr3'],
   ...     lociNames=['rs4', 'rs5', 'rs6'])
   >>> # 2 2 2,
   >>> pop1.setGenotype([2])
   >>> # 1 1 1, 0 0 0, 2 2 2
   >>> pop.addChromFrom(pop1)
   >>> # 1 1 1, 0 0 0, 2 0 2 2 0
   >>> pop.addLoci(chrom=[2, 2], pos=[1.5, 3.5], lociNames=['rs7', 'rs8'])
   (7, 10)
   >>> # 1 1 1, 0 0 0, 2 0 2 0
   >>> pop.removeLoci(8)
   >>> # loci names can also be used.
   >>> pop.removeLoci(['rs1', 'rs7'])
   >>> sim.dump(pop)
   Ploidy: 2 (diploid)
   Chromosomes:
   1: chr1 (AUTOSOME, 3 loci)
      (1),  (2),  (3)
   2: chr2 (AUTOSOME, 2 loci)
     rs2 (1), rs3 (2)
   3: chr3 (AUTOSOME, 3 loci)
     rs4 (1), rs6 (3), rs8 (3.5)
   population size: 10 (1 subpopulations with 10 Individuals)
   Number of ancestral populations: 0

   SubPopulation 0 (), 10 Individuals:
      0: MU 111 00 220 | 111 00 220 
      1: MU 111 00 220 | 111 00 220 
      2: MU 111 00 220 | 111 00 220 
      3: MU 111 00 220 | 111 00 220 
      4: MU 111 00 220 | 111 00 220 
      5: MU 111 00 220 | 111 00 220 
      6: MU 111 00 220 | 111 00 220 
      7: MU 111 00 220 | 111 00 220 
      8: MU 111 00 220 | 111 00 220 
      9: MU 111 00 220 | 111 00 220 

   >>> # add loci from another population 
   >>> pop2 = sim.Population(10, loci=2, lociPos=[0.1, 2.2], chromNames='chr3')
   >>> pop.addLociFrom(pop2)
   WARNING: Chromosome 'chr3' is merged to chromosome 'chr1'.
   >>> pop.addLociFrom(pop2, byName=2)
   >>> sim.dump(pop, genotype=False)
   Ploidy: 2 (diploid)
   Chromosomes:
   1: chr1 (AUTOSOME, 5 loci)
      (0.1),  (1),  (2),  (2.2),  (3)
   2: chr2 (AUTOSOME, 2 loci)
     rs2 (1), rs3 (2)
   3: chr3 (AUTOSOME, 5 loci)
      (0.1), rs4 (1),  (2.2), rs6 (3), rs8 (3.5)
   population size: 10 (1 subpopulations with 10 Individuals)
   Number of ancestral populations: 0


   now exiting runScriptInteractively...

`Download addRemoveLoci.py <addRemoveLoci.py>`_


Remove or extract individuals and subpopulations from a population
------------------------------------------------------------------

Functions :meth:`Population.removeIndividuals` and
:meth:`Population.removeSubPops` remove selected individuals or groups of
individuals from a population. Functions :meth:`Population.extractIndividuals`
and :meth:`Population.extractSubPops` extract individuals and subpopulations
from an existing population and form a new one.

Functions ``removeIndividauls`` and ``extractIndividuals`` could be used to
remove or extract individuals from the present generation by indexes or from all
ancestral generations by IDs or a Python filter function. This function should
accept parameter ``ind`` or one or more information fields. simuPOP will pass
individual for parameter ``ind``, and values at specified information fields
(``age`` in this example) of each individual to this function. The present
population structure will be kept, even if some subpopulations are left empty.
For example, you could remove the first thirty individuals of a population using
::

   pop.removeIndividuals(indexes=range(30))

or remove all individuals at age 20 or 30 using   ::

   pop.removeIndividuals(IDs=(20, 30), idField='age')

or remove all individuals with age between 20 and 30 using   ::

   pop.removeIndividuals(filter=lambda age: age >=20 and age <=30)

. In the last example, a Python lambda function is defined to avoid the
definition of a named function.

Functions ``removeSubPops`` or ``extractSubPops`` could be used to remove or
extract subpopulations, or goups of individuals defined by virtual
subpopulations from a population. The latter case is very interesting because it
could be used to remove or extract individuals with similar properties, such as
all individuals between the ages 40 and 60, as demonstrated in Example
:ref:`extract <extract>`.

.. _extract:

**Example**: *Extract individuals*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(size=[200, 200], loci=[5, 5], infoFields='age')
   >>> sim.initGenotype(pop, genotype=range(10))
   >>> sim.initInfo(pop, lambda: random.randint(0,75), infoFields='age')
   >>> pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[20, 60]))
   >>> # remove individuals
   >>> pop.removeIndividuals(indexes=range(0, 300, 10))
   >>> print(pop.subPopSizes())
   (180, 190)
   >>> # remove individuals using IDs
   >>> pop.setIndInfo([1, 2, 3, 4], field='age')
   >>> pop.removeIndividuals(IDs=[2, 4], idField='age')
   >>> # remove indiviuals using a filter function
   >>> sim.initSex(pop)
   >>> pop.removeIndividuals(filter=lambda ind: ind.sex() == sim.MALE)
   >>> print([pop.individual(x).sex() for x in range(8)])
   [2, 2, 2, 2, 2, 2, 2, 2]
   >>> #
   >>> # remove subpopulation
   >>> pop.removeSubPops(1)
   >>> print(pop.subPopSizes())
   (56,)
   >>> # remove virtual subpopulation (people with age between 20 and 60)
   >>> pop.removeSubPops([(0, 1)])
   >>> print(pop.subPopSizes())
   (56,)
   >>> # extract another virtual subpopulation (people with age greater than 60)
   >>> pop1 = pop.extractSubPops([(0,2)])
   >>> sim.dump(pop1, structure=False, max=10)
   SubPopulation 0 (), 0 Individuals:


   now exiting runScriptInteractively...

`Download extract.py <extract.py>`_


.. _subsec_Population_Variables:

Store arbitrary population information as population variables
--------------------------------------------------------------

.. index:: single: operator; Stat

Each simuPOP population has a Python dictionary that can be used to store
arbitrary Python variables. These variables are usually used by various
operators to share information between them. For example, the ``Stat`` operator
calculates population statistics and stores the results in this Python
dictionary. Other operators such as the :class:`PyEval` and ``TerminateIf``\
read from this dictionary and act upon its information.

.. index::
   single: Population; vars
   single: Population; Population

simuPOP provides two functions, namely ``Population.vars()`` and
``Population.dvars()`` to access a population dictionary. These functions return
the same dictionary object but :func:`dvars`\ () returns a wrapper class so that
you can access this dictionary as attributes. For example,
``pop.vars()['alleleFreq'][0]`` is equivalent to ``pop.dvars().alleleFreq[0]``.
Because dictionary ``subPop[spID]`` is frequently used by operators to store
variables related to a particular (virtual) subpopulation, function
``pop.vars(subPop)`` is provided as a shortcut to
``pop.vars()['subPop'][spID]``. Example :ref:`popVars <popVars>` demonstrates
how to set and access population variables.

.. _popVars:

**Example**: *population variables*

::

   >>> import simuPOP as sim
   >>> from pprint import pprint
   >>> pop = sim.Population(100, loci=2)
   >>> sim.initGenotype(pop, freq=[0.3, 0.7])
   >>> print(pop.vars())    # No variable now
   {}
   >>> pop.dvars().myVar = 21
   >>> print(pop.vars())
   {'myVar': 21}
   >>> sim.stat(pop, popSize=1, alleleFreq=0)
   >>> # pprint prints in a less messy format
   >>> pprint(pop.vars())
   {'alleleFreq': {0: defdict({0: 0.275, 1: 0.725})},
    'alleleNum': {0: defdict({0: 55.0, 1: 145.0})},
    'myVar': 21,
    'popSize': 100,
    'subPopSize': [100]}
   >>> # print number of allele 1 at locus 0
   >>> print(pop.vars()['alleleNum'][0][1])
   145.0
   >>> # use the dvars() function to access dictionary keys as attributes
   >>> print(pop.dvars().alleleNum[0][1])
   145.0
   >>> print(pop.dvars().alleleFreq[0])
   defdict({0: 0.275, 1: 0.725})

   now exiting runScriptInteractively...

`Download popVars.py <popVars.py>`_

It is important to understand that this dictionary forms a **local namespace**
in which Python expressions can be evaluated. This is the basis of how
expression-based operators work. For example, the ``PyEval``\ operator in
example :ref:`simple_example <simple_example>` evaluates expression
````'%.2f\\t' % LD[0][1]''`` in each population's local namespace when it is
applied to that population. This yields different results for different
population because their LD values are different. In addition to Python
expressions, Python statements can also be executed in the local namespace of a
population, using the ``stmts`` parameter of the :class:`PyEval` or
:class:`PyExec` operator. Example :ref:`expression <expression>` demonstrates
the use of a simuPOP terminator, which terminates the evolution of a population
when its expression is evaluated as ``True``. Note that The ``evolve()``\
function of this example does not specify how many generations to evolve so it
will stop only after all replicates stop. The return value of this function
indicates how many generations each replicate has evolved. This example also
demonstrates how to run multiple replicates of an evolutionary process, which we
will discuss in detail latter.

.. _expression:

**Example**: *Expression evaluation in the local namespace of a population*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(sim.Population(100, loci=1), rep=5)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     matingScheme = sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.TerminateIf('len(alleleFreq[0]) == 1')
   ...     ]
   ... )
   (129, 1540, 180, 247, 242)

   now exiting runScriptInteractively...

`Download expression.py <expression.py>`_


.. _subsec_Save_and_Load:

Save and load a population
--------------------------

.. index::
   single: Population; save
   single: function; loadPopulation

simuPOP populations can be saved to and loaded from disk files using
``Population.save(file)`` member function and global function
:func:`loadPopulation`. **Virtual splitters are not saved** because they are
considered as runtime definitions. Although files in any extension can be used,
extension ``.pop`` is recommended. Example :ref:`savePop <savePop>` demonstrates
how to save and load a population in the native simuPOP format.

.. _savePop:

**Example**: *Save and load a population*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(100, loci=5, chromNames=['chrom1'])
   >>> pop.dvars().name = 'my sim.Population'
   >>> pop.save('sample.pop')
   >>> pop1 = sim.loadPopulation('sample.pop')
   >>> pop1.chromName(0)
   'chrom1'
   >>> pop1.dvars().name
   'my sim.Population'

   now exiting runScriptInteractively...

`Download savePop.py <savePop.py>`_

The native simuPOP format is portable across different platforms but is not
human readable and is not recognized by other applications. If you need to save
a simuPOP population in a format that is recognizable by a particular software,
you can use functions ``importPopulation``, ``export``, and operator
``Exporter`` if you would like to export populations during evolution. These
functions are defined in module :mod:`simuPOP.utils`.


Import and export datasets in unsupported formats \*
----------------------------------------------------

simuPOP provides a few utility functions to import and export populations in
common formats such as GENEPOP, Phylip, and STRUCTURE (see chapter utility
modules for details). If you need to import data from a file in a format that is
not currently supported, you generally need to first scan the file for
information such as number and names of chromosomes, loci, alleles,
subpopulation, and individuals. After you create a population without genotype
information from these parameters, you can scan the file for the second time and
fill the population with genotypes and other information. Example
:ref:`importData <importData>` demonstrates how to define a function to import
from a file that is saved by function :func:`~simuPOP.utils.saveCSV`.

.. _importData:

**Example**: *Import a population from another file format*

::

   >>> import simuPOP as sim
   >>> def importData(filename):
   ...     'Read data from ``filename`` and create a population'
   ...     data = open(filename)
   ...     header = data.readline()
   ...     fields = header.split(',')
   ...     # columns 1, 3, 5, ..., without trailing '_1'
   ...     names = [fields[x].strip()[:-2] for x in range(1, len(fields), 2)]
   ...     popSize = 0
   ...     alleleNames = set()
   ...     for line in data.readlines():
   ...         # get all allele names
   ...         alleleNames |= set([x.strip() for x in line.split(',')[1:]])
   ...         popSize += 1
   ...     # create a population
   ...     alleleNames = list(alleleNames)
   ...     pop = sim.Population(size=popSize, loci=len(names), lociNames=names,
   ...         alleleNames=alleleNames)
   ...     # start from beginning of the file again
   ...     data.seek(0)
   ...     # discard the first line
   ...     data.readline()
   ...     for ind, line in zip(pop.individuals(), data.readlines()):
   ...         fields = [x.strip() for x in line.split(',')]
   ...         sex = sim.MALE if fields[0] == '1' else sim.FEMALE
   ...         ploidy0 = [alleleNames.index(fields[x]) for x in range(1, len(fields), 2)]
   ...         ploidy1 = [alleleNames.index(fields[x]) for x in range(2, len(fields), 2)]
   ...         ind.setGenotype(ploidy0, 0)
   ...         ind.setGenotype(ploidy1, 1)
   ...         ind.setSex(sex)
   ...     # close the file
   ...     data.close()
   ...     return pop
   ... 
   >>> from simuPOP.utils import saveCSV
   >>> pop = sim.Population(size=[10], loci=[3, 2], lociNames=['rs1', 'rs2', 'rs3', 'rs4', 'rs5'],
   ...     alleleNames=['A', 'B'])
   >>> sim.initSex(pop)
   >>> sim.initGenotype(pop, freq=[0.5, 0.5])
   >>> # output sex but not affection status.
   >>> saveCSV(pop, filename='sample.csv', affectionFormatter=None,
   ...     sexFormatter={sim.MALE:1, sim.FEMALE:2})
   >>> # have a look at the file
   >>> print(open('sample.csv').read())
   sex, rs1_1, rs1_2, rs2_1, rs2_2, rs3_1, rs3_2, rs4_1, rs4_2, rs5_1, rs5_2
   2, B, B, B, B, B, A, A, B, B, A
   2, B, A, B, A, B, A, A, A, A, B
   1, B, B, B, B, B, B, B, B, B, A
   1, B, A, B, A, B, B, B, A, A, A
   1, B, B, B, B, B, B, A, A, B, A
   1, A, B, B, A, B, B, B, A, B, B
   1, B, B, B, B, B, B, B, B, A, A
   2, B, B, A, A, B, A, A, A, B, A
   2, A, B, B, B, A, B, B, A, A, B
   2, B, A, A, B, A, A, B, B, B, A

   >>> pop1 = importData('sample.csv')
   >>> sim.dump(pop1)
   Ploidy: 2 (diploid)
   Chromosomes:
   1:  (AUTOSOME, 5 loci)
     rs1 (1), rs2 (2), rs3 (3), rs4 (4), rs5 (5)
   population size: 10 (1 subpopulations with 10 Individuals)
   Number of ancestral populations: 0

   SubPopulation 0 (), 10 Individuals:
      0: FU BBBAB | BBABA 
      1: FU BBBAA | AAAAB 
      2: MU BBBBB | BBBBA 
      3: MU BBBBA | AABAA 
      4: MU BBBAB | BBBAA 
      5: MU ABBBB | BABAB 
      6: MU BBBBA | BBBBA 
      7: FU BABAB | BAAAA 
      8: FU ABABA | BBBAB 
      9: FU BAABB | ABABA 


   now exiting runScriptInteractively...

`Download importData.py <importData.py>`_

Unless there are specific requirements in the order and labeling of individuals,
exporting a simuPOP population is usually straightforward. Functions that are
useful in such occasions include structural functions
:meth:`Population.numSubPop`\ (), :meth:`Population.subPopName`\ ``,
Population.popSize()`` and :meth:`Population.subPopSizes`\ (), ``and``
individual access functions :meth:`Population.individual`\ () and
:meth:`Population.individuals`\ () and individual population access functions
such as :meth:`Individual.allele`\ () and :meth:`Individual.info`\ (). Function
``saveFSTAT`` in the cookbook module ``fstatUtil`` or ``saveCSV`` in module
:mod:`simuPOP.utils` are good examples you can follow.


