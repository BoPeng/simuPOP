Pedigrees
=========


Create a pedigree object
------------------------

A :class:`Pedigree` object is basically a static population object that is used
to track relationship between individuals. An unique ID is required for all
individuals so that individuals could be identified easily using their IDs.
Individuals in a pedigree usually have one or two information fields to record
the IDs of their parents. Operators :class:`IdTagger` and
:class:`PedigreeTagger` are usually used to maintain these information fields
which are, although customizable, almost always ``ind_id``, ``father_id`` and
``mother_id``. After pedigrees are identified, population operations could be
applied, for example, to extracted identified pedigrees from an existing
population. This is basically how module :mod:`simuPOP.sampling` works.

A new pedigree can be created from a population object with an ID field (default
to ``ind_id``), and two optional parental ID fields (default to ``father_id``
and ``mother_id``). For example,   ::

   ped = Pedigree(pop, infoFields=ALL_AVAIL)

will create a pedigree object from population ``pop`` with information fields
``ind_id``, ``father_id`` and ``mother_id``, copying all available information
fields. The ID field should have an unique ID for each individual and the
parental ID fields should record the ID of his or her parents. Genotype
information and additional information fields can be copied to a pedigree object
if needed. The population object is unchanged.

Another method is to directly convert a population object to a pedigree object,
using member function ``asPedigree`` of a population class. For example,  ::

   pop.asPedigree()

will convert the existing population to a pedigree object. Object pop can then
be able to call all pedigree member functions. Once your task is done, you can
convert the object back to a population using the :meth:`Pedigree.asPopulation`\
() member function of the object.

A pedigree object can also be created from a file saved by function
:meth:`Pedigree.save`\ () or operator :class:`PedigreeTagger` using function
:func:`loadPedigree`. Please refer to section *save and load pedigrees* in
details.


Locate close and remote relatives of each individual
----------------------------------------------------

A pedigree object provides several functions for you to identify spouse, sibling
and more distant relatives of each individual. The results are stored to
additional information fields of each individual. For example, if you would like
to know the offspring of all individuals, you can call function
:meth:`Pedigree.locateRelatives` as follows:

::

   offFields = ['off1', 'off2', 'off3']
   ped.addInfoFields(offFields)
   ped.locateRelatives(OFFSPRING, resultFields=offFields)

This function will locate up to 3 (determined by the length of ``resultFields``)
offspring of each individual and put their IDs in specified informaton fields.
This function allows you to identify spouses (it is common to have multiple
spouses when random mating is used), outbred spouse (exclude spouses who share
at least one of the parents), offspring (all offspring) and common offspring
with a specified spouse, siblings (share at least one parent) and full siblings
(share two parents). It also allows you to limit the result by sex and affection
status (e.g. find only affected female offspring).

More distant relationship can be derived from these relationship using function
:meth:`Pedigree.traceRelatives`. This function accepts a path of information
fields and follows the path to identify relatives. For example

::

   sibFields = ['sib1', 'sib2']
   offFields = ['off1', 'off2', 'off3']
   cousinFields = ['cousin1', 'cousin2', 'cousin3']
   ped.addInfoFields(sibFields + offFields + cousinFields)
   ped.locateRelatives(FULLSIBLING, resultFields=sibFields)
   ped.locateRelatives(OFFSPRING, resultFields=offFields)
   ped.traceRelatives([['father_id', 'mother_id'], sibFields, offFields],
       sex=[ANY_SEX, MALE_ONLY, FEMALE_ONLY],
       resultField=cousinFields)

would first identify full siblings and offspring of all individuals and then
locate father or mother's male sibling's daughters. As you can imagine, this
function can be used to track very complicated relationships.

This function also provides a function for you to identify individuals with
specified relatives. Example :ref:`locateRelative <locateRelative>` gives an
example how to locate a grandfather with at least five grandchildren. With such
information, functions such as :meth:`Population.extractIndividuals`\ () could
be used to extract Pedigrees from a population. This is basically how
:mod:`simuPOP.sampling` module works.

.. _locateRelative:

**Example**: *Locate close and distant relatives of individuals*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, ancGen=2, infoFields=['ind_id', 'father_id', 'mother_id'])
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.IdTagger(),
   ...     ],
   ...     matingScheme=sim.RandomMating(
   ...         numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 4),
   ...         ops=[
   ...             sim.MendelianGenoTransmitter(),
   ...             sim.IdTagger(),
   ...             sim.PedigreeTagger()
   ...         ],
   ...     ),
   ...     gen = 5
   ... )
   5
   >>> ped = sim.Pedigree(pop)
   >>> offFields = ['off%d' % x for x in range(4)]
   >>> grandOffFields = ['grandOff%d' % x for x in range(5)]
   >>> ped.addInfoFields(['spouse'] + offFields + grandOffFields)
   >>> # only look spouse for fathers...
   >>> ped.locateRelatives(sim.OUTBRED_SPOUSE, ['spouse'], sex=sim.FEMALE_ONLY)
   >>> ped.locateRelatives(sim.COMMON_OFFSPRING, ['spouse'] + offFields)
   >>> # trace offspring of offspring
   >>> ped.traceRelatives([offFields, offFields], resultFields=grandOffFields)
   True
   >>> # 
   >>> IDs = ped.individualsWithRelatives(grandOffFields)
   >>> # check on ID.
   >>> grandFather = IDs[0]
   >>> grandMother = ped.indByID(grandFather).spouse
   >>> # some ID might be invalid.
   >>> children = [ped.indByID(grandFather).info(x) for x in offFields]
   >>> childrenSpouse = [ped.indByID(x).spouse for x in children if x >= 1]
   >>> childrenParents = [ped.indByID(x).father_id for x in children if x >= 1] \
   ...     + [ped.indByID(x).mother_id for x in children if x >= 1]
   >>> grandChildren = [ped.indByID(grandFather).info(x) for x in grandOffFields]
   >>> grandChildrenParents = [ped.indByID(x).father_id for x in grandChildren if x >= 1] \
   ...     + [ped.indByID(x).mother_id for x in grandChildren if x >= 1]
   >>> 
   >>> def idString(IDs):
   ...     uniqueIDs = list(set(IDs))
   ...     uniqueIDs.sort()
   ...     return ', '.join(['%d' % x for x in uniqueIDs if x >= 1])
   ... 
   >>> print('''GrandParents: %d, %d
   ... Children: %s
   ... Spouses of children: %s
   ... Parents of children: %s
   ... GrandChildren: %s
   ... Parents of grandChildren: %s ''' % \
   ... (grandFather, grandMother, idString(children), idString(childrenSpouse),
   ...     idString(childrenParents), idString(grandChildren), idString(grandChildrenParents)))
   GrandParents: 3040, 3847
   Children: 4078, 4079, 4080
   Spouses of children: 4446, 4797
   Parents of children: 3040, 3847
   GrandChildren: 5188, 5189, 5879, 5880, 5881
   Parents of grandChildren: 4078, 4079, 4446, 4797 
   >>> 
   >>> # let us look at the structure of this complete pedigree using another method
   >>> famSz = ped.identifyFamilies()
   >>> # it is amazing that there is a huge family that connects almost everyone
   >>> len(famSz), max(famSz)
   (533, 2383)
   >>> # if we only look at the last two generations, things are much better
   >>> ped.addInfoFields('ped_id')
   >>> famSz = ped.identifyFamilies(pedField='ped_id', ancGens=[0,1])
   >>> len(famSz), max(famSz)
   (664, 114)

   now exiting runScriptInteractively...

`Download locateRelative.py <locateRelative.py>`_


Identify pedigrees (related individuals)
----------------------------------------

The :class:`Pedigree` class provides some other functions that allows you to
identify related individuals. For example,

* Function :meth:`Pedigree.identifyAncestors` identifies all ancestors of
  specified individuals or all individuals at the present generation. In a
  diaploid population when there is only one parent, you can see that only a small
  portion of ancestors have offspring in the last generation.

* Function :meth:`Pedigree.identifyOffspring` identifies all offspring of
  specified individuals across multiple generations.

* Function :meth:`Pedigree.identifyFamilies` groups all related individuals into
  families and assign a family ID to all family members. You might be surprised by
  how large this kind of family can be when parents are allowed to have multiple
  spouses.

All these functions support parameters ``subPops`` and ``ancGens`` so that you
can limit your search in specific subpopulations and ancestral generations. For
example, you can limit your search to all male individuals to find out someone's
male offspring. Example :ref:`locateFamilies <locateFamilies>` demonstrates how
to use these functions to analyze the structure of a complete pedigree.

.. _locateFamilies:

**Example**: *Identify all ancestors*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, ancGen=-1, infoFields=['ind_id', 'father_id', 'mother_id'])
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.IdTagger(),
   ...     ],
   ...     matingScheme=sim.RandomMating(
   ...         numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 4),
   ...         ops=[
   ...             sim.MendelianGenoTransmitter(),
   ...             sim.IdTagger(),
   ...             sim.PedigreeTagger()
   ...         ],
   ...     ),
   ...     gen = 19
   ... )
   19
   >>> # we now have the complete pedigree of 20 generations
   >>> pop.asPedigree()
   >>> # total number of individuals should be 20 * 1000
   >>> # how many families do we have?
   >>> fam = pop.identifyFamilies()
   >>> len(fam)
   525
   >>> # but how many families with more than 1 individual?
   >>> # The rest of them must be in the initial generation
   >>> len([x for x in fam if x > 1])
   18
   >>> # let us look backward. allAnc are the ancestors who have offspring in the
   >>> # last generation. You can see this is a small number compared the number of
   >>> # ancestors.
   >>> allAnc = pop.identifyAncestors()
   >>> len(allAnc)
   8614

   now exiting runScriptInteractively...

`Download locateFamilies.py <locateFamilies.py>`_


Save and load pedigrees
-----------------------

A complete pedigree, including ID, sex and affection status of each individual,
IDs of their parents, and optionally values of some information fields and
genotypes at some loci could be saved to a file, and be loaded using function
:func:`loadPedigree`. The loaded pedigree could be analyzed using pedigree
functions, or be used to direct the evolution of another evolutionary process
using a pedigree mating scheme.

A pedigree could be saved in two ways. In the first method, a pedigree could be
created using the methods described above and be saved using function
:meth:`Pedigree.save`\ (). However, if the population is large, recording all
ancestral generations may not be feasible. If this is the case, you can use a
:class:`PedigreeTagger` operator to save individual information during the
evolution. If you do not care about details of the top-most ancestral
generation, a PedigreeTagger used in a mating scheme should be enough to record
pedigree information of all offspring. Individual in the top-most generation who
have offspring in the next generation will be constructed in
:func:`loadPedigree`. If you would like to include detailed information about
all individuals in the top-most ancestral generation, you can use a
:class:`PedigreeTagger` in the ``initOps`` parameter of the
:meth:`Simulator.evolve`\ () or :meth:`Population.evolve`\ () function.

Example :ref:`saveLoadPedigree <saveLoadPedigree>` demonstrates how to use these
functions to analyze the structure of a complete pedigree.

.. _saveLoadPedigree:

**Example**: *Save and load a complete pedigree*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(4, loci=1, infoFields=['ind_id', 'father_id', 'mother_id'],
   ...     ancGen=-1)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.IdTagger(),
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...         sim.PedigreeTagger(output='>>pedigree.ped', outputLoci=0)
   ...     ],
   ...     matingScheme=sim.RandomMating(
   ...         ops=[
   ...             sim.MendelianGenoTransmitter(),
   ...             sim.IdTagger(),
   ...             sim.PedigreeTagger(output='>>pedigree.ped', outputLoci=0)
   ...         ],
   ...     ),
   ...     gen = 2
   ... )
   2
   >>> #
   >>> print(open('pedigree.ped').read())
   1 0 0 F U 0 0
   2 0 0 F U 0 1
   3 0 0 M U 1 1
   4 0 0 M U 1 1
   5 4 1 M U 0 1
   6 4 2 F U 1 1
   7 3 2 F U 0 1
   8 3 2 M U 1 1
   9 8 7 F U 1 1
   10 5 6 M U 1 1
   11 5 6 M U 1 1
   12 5 7 F U 0 1

   >>> pop.asPedigree()
   >>> pop.save('pedigree1.ped', loci=0)
   >>> print(open('pedigree1.ped').read())
   1 0 0 F U 0 0
   2 0 0 F U 0 1
   3 0 0 M U 1 1
   4 0 0 M U 1 1
   5 4 1 M U 0 1
   6 4 2 F U 1 1
   7 3 2 F U 0 1
   8 3 2 M U 1 1
   9 8 7 F U 1 1
   10 5 6 M U 1 1
   11 5 6 M U 1 1
   12 5 7 F U 0 1

   >>> # 
   >>> ped = sim.loadPedigree('pedigree1.ped')
   >>> sim.dump(ped, ancGens=range(3))
   Ploidy: 2 (diploid)
   Chromosomes:
   1:  (AUTOSOME, 1 loci)
      (1)
   Information fields: 
   ind_id father_id mother_id 
   population size: 4 (1 subpopulations with 4 Individuals)
   Number of ancestral populations: 2

   SubPopulation 0 (), 4 Individuals:
      0: FU 1 | 1 |  9 8 7
      1: MU 1 | 1 |  10 5 6
      2: MU 1 | 1 |  11 5 6
      3: FU 0 | 1 |  12 5 7

   Ancestral population 1
   SubPopulation 0 (), 4 Individuals:
      0: MU 0 | 1 |  5 4 1
      1: FU 1 | 1 |  6 4 2
      2: FU 0 | 1 |  7 3 2
      3: MU 1 | 1 |  8 3 2

   Ancestral population 2
   SubPopulation 0 (), 4 Individuals:
      0: FU 0 | 0 |  1 0 0
      1: FU 0 | 1 |  2 0 0
      2: MU 1 | 1 |  3 0 0
      3: MU 1 | 1 |  4 0 0


`Download saveLoadPedigree.py <saveLoadPedigree.py>`_


