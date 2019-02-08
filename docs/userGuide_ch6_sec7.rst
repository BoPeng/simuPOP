Evolve a population following a specified pedigree structure \*\*
=================================================================

There are some applications where you would like to repeat the same evolutionary
process repeatedly using the same pedigree structure. For example, a gene-
dropping simulation method basically initialize leaves of a pedigree with random
genotypes and pass the genotypes along the pedigree according to Mendelian laws.
This can be done in simuPOP using a pedigree mating scheme.

A pedigree mating scheme :class:`PedigreeMating` evolves a population following
an existing pedigree structure. If the  :class:`Pedigree` object has ``N``
ancestral generations and a present generation, it can be used to evolve a
population for ``N`` generations, starting from the topmost ancestral
generation. At the *k*-th generation, this mating scheme produces an offspring
generation according to subpopulation structure of the ``N-k-1`` ancestral
generation in the pedigree object (e.g. producing the offspring population of
generation 0 according to the ``N-1`` ancestral generation of the pedigree
object). For each offspring, this mating scheme copies individual ID and sex
from the corresponing individual in the pedigree object. It then locates the
parents of each offspring using their IDs in the pedigree object. A list of
during mating operators are then used to transmit parental genotype to the
offspring.

To use this mating scheme, you should

* Prepare a pedigree object with ``N`` ancestral generations (and a present
  generation). Parental information should be available at the present, parental,
  ..., and ``N-1`` ancestral generations. This object could be created by evolving
  a population with ``ancGen`` set to -1 with parental information tracked by
  operators ``idTagger()`` and ``pedigreeTagger()``.

* Prepare the population so that it contains individuals with IDs matching this
  generation, or at least individuals who have offspring in the next topmost
  ancestral generation. Because individuals in such a population will parent
  offsprings at the ``N-1`` ancestral generation of the pedigree object, it is a
  good idea to assign ``ind_id`` using ``ped.indInfo('father_id')`` and
  ``ped.infInfo('mother_id')`` of the ``N-1`` ancestral generation of ``ped``.

* Evolve the population using a :class:`PedigreeMating` mating scheme for ``N``
  or less generations. Because parents are chosen by their IDs, subpopulation
  structure is ignored and migration will have no effect on the evolutionary
  process. No :class:`IdTagger` should be used to assign IDs to offspring because
  re-labeling IDs will confuse this mating scheme. This mating scheme copies
  individual sex from pedigree individual to each offspring because individual sex
  may affect the way genotypes are transmitted (e.g. a
  :class:`MendelianGenoTransmitter`\ () with sex chromosomes).

Example :ref:`pedigreeMating <pedigreeMating>` demonstrates how to create a
complete pedigree by evolving a population without genotype, and then replay the
evolutionary process using another population.

.. _pedigreeMating:

**Example**: *Use a pedigree mating scheme to replay an evolutionary process.*

::

   >>> import simuPOP as sim
   >>> # create a population without any genotype
   >>> from simuPOP.utils import migrSteppingStoneRates
   >>> ped = sim.Population(size=[1000]*5, ancGen=-1, 
   ...     infoFields=['ind_id', 'father_id', 'mother_id', 'migrate_to'])
   >>> ped.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.IdTagger(),
   ...     ],
   ...     preOps=sim.Migrator(rate=migrSteppingStoneRates(0.1, 5)),
   ...     matingScheme=sim.RandomMating(
   ...         numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 4),
   ...         ops=[
   ...             # we do not even need a genotype transmitter...
   ...             sim.IdTagger(),
   ...             sim.PedigreeTagger(),
   ...         ]),
   ...     gen=100
   ... )
   100
   >>> # convert itself to a pedigree object
   >>> ped.asPedigree()
   >>> # we should have 100 ancestral generations
   >>> N = ped.ancestralGens()
   >>> # We should have 101 * 1000 * 5 individuals, but how many actually
   >>> # contribute genotype to the last generation?
   >>> anc = ped.identifyAncestors()
   >>> len(anc)
   205647
   >>> # remove individuals who do not contribute genotype to the last generation
   >>> allIDs = [x.ind_id for x in ped.allIndividuals()]
   >>> removedIDs = list(set(allIDs) - set(anc))
   >>> ped.removeIndividuals(IDs=removedIDs)
   >>> # now create a top most population, but we do not need all of them
   >>> # so we record only used individuals
   >>> IDs = [x.ind_id for x in ped.allIndividuals(ancGens=N)]
   >>> sex = [x.sex() for x in ped.allIndividuals(ancGens=N)]
   >>> # create a population, this time with genotype. Note that we do not need
   >>> # populaton structure because PedigreeMating disregard population structure.
   >>> pop = sim.Population(size=len(IDs), loci=1000, infoFields='ind_id')
   >>> # manually initialize ID and sex
   >>> sim.initInfo(pop, IDs, infoFields='ind_id')
   >>> sim.initSex(pop, sex=sex)
   >>> pop.evolve(
   ...     initOps=sim.InitGenotype(freq=[0.4, 0.6]),
   ...     # we do not need migration, or set number of offspring,
   ...     # or demographic model, but we do need a genotype transmitter
   ...     matingScheme=sim.PedigreeMating(ped, 
   ...         ops=sim.MendelianGenoTransmitter()),
   ...     gen=100
   ... )
   100
   >>> # let us compare the pedigree and the population object
   >>> print(ped.indInfo('ind_id')[:5])
   (500001.0, 500002.0, 500003.0, 500004.0, 500005.0)
   >>> print(pop.indInfo('ind_id')[:5])
   (500001.0, 500002.0, 500003.0, 500004.0, 500005.0)
   >>> print([ped.individual(x).sex() for x in range(5)])
   [1, 2, 1, 1, 2]
   >>> print([pop.individual(x).sex() for x in range(5)])
   [1, 2, 1, 1, 2]
   >>> print(ped.subPopSizes())
   (663, 1254, 1213, 1230, 640)
   >>> print(pop.subPopSizes())
   (663, 1254, 1213, 1230, 640)

   now exiting runScriptInteractively...

`Download pedigreeMating.py <pedigreeMating.py>`_

As long as unique IDs are used for individuals in different generations, the
same technique could be used for overlapping generations as well. Even if some
individuals are copied from generation to generation, separate IDs should be
assigned to these individuals so that a pedigree could be correctly constructed.
Because these individuals are copied from a single parent, the pedigree object
will have mixed number of parents (some individuals have one parent, some have
two). If :class:`PedigreeTagger` operators are used to record parental
information, such a pedigree could be loaded by function :func:`loadPedigree`.
Example :ref:`pedigreeMatingAgeStructured <pedigreeMatingAgeStructured>` evolves
an age-structured population. Instead of saving all ancestral generations to a
population object and convert it to a pedigree, this example saves the complete
pedigree to file ``structure.ped`` and load the pedigree using function
:func:`loadPedigree`.

.. _pedigreeMatingAgeStructured:

**Example**: *Replay an evolutionary process of an age-structured population*

::

   >>> import simuPOP as sim
   >>> 
   >>> import random
   >>> N = 10000
   >>> pop = sim.Population(N, infoFields=['age', 'ind_id', 'father_id', 'mother_id'])
   >>> # we simulate age 0, 1, 2, 3 
   >>> pop.setVirtualSplitter(sim.InfoSplitter(field='age', values=[0, 1, 2, 3]))
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         # random assign age
   ...         sim.InitInfo(lambda: random.randint(0, 3), infoFields='age'),
   ...         # random genotype
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...         # assign an unique ID to everyone.
   ...         sim.IdTagger(),
   ...     ],
   ...     # increase the age of everyone by 1 before mating.
   ...     preOps=sim.InfoExec('age += 1'),
   ...     matingScheme=sim.HeteroMating([
   ...         # age 1, 2 will be copied
   ...         sim.CloneMating(
   ...             ops=[
   ...                 # This will set offspring ID
   ...                 sim.CloneGenoTransmitter(),
   ...                 # new ID for offspring in order to track pedigree
   ...                 sim.IdTagger(),
   ...                 # both offspring and parental IDs will be the same
   ...                 sim.PedigreeTagger(output='>>structured.ped'),
   ...             ],
   ...             subPops=[(0,1), (0,2)],
   ...             weight=-1
   ...         ),
   ...         # age 2 produce offspring
   ...         sim.RandomMating(
   ...             ops=[
   ...                 # new ID for offspring
   ...                 sim.IdTagger(),
   ...                 # record complete pedigree
   ...                 sim.PedigreeTagger(output='>>structured.ped'),
   ...                 sim.MendelianGenoTransmitter(),   # transmit genotype
   ...             ],
   ...             subPops=[(0,2)]
   ...         )]
   ...     ),
   ...     gen=20
   ... )
   20
   >>> 
   >>> # use a pedigree object recovered from a file saved by operator PedigreeTagger
   >>> ped = sim.loadPedigree('structured.ped')
   >>> # create a top most population, but we do not need all of them
   >>> # so we record only used individuals
   >>> IDs = [x.ind_id for x in ped.allIndividuals(ancGens=ped.ancestralGens())]
   >>> sex = [x.sex() for x in ped.allIndividuals(ancGens=ped.ancestralGens())]
   >>> # create a population, this time with genotype. Note that we do not need
   >>> # populaton structure because PedigreeMating disregard population structure.
   >>> pop = sim.Population(size=len(IDs), loci=1000, infoFields='ind_id')
   >>> # manually initialize ID and sex
   >>> sim.initInfo(pop, IDs, infoFields='ind_id')
   >>> sim.initSex(pop, sex=sex)
   >>> pop.evolve(
   ...     initOps=sim.InitGenotype(freq=[0.4, 0.6]),
   ...     # we do not need migration, or set number of offspring,
   ...     # or demographic model, but we do need a genotype transmitter
   ...     matingScheme=sim.PedigreeMating(ped, 
   ...         ops=sim.IfElse(lambda mom: mom is None,
   ...                 sim.CloneGenoTransmitter(),
   ...                 sim.MendelianGenoTransmitter())
   ...     ),
   ...     gen=100
   ... )
   20
   >>> # 
   >>> print(pop.indInfo('ind_id')[:5])
   (200001.0, 200002.0, 200003.0, 200004.0, 200005.0)
   >>> print([pop.individual(x).sex() for x in range(5)])
   [1, 2, 2, 1, 1]
   >>> # The pedigree object does not have population structure
   >>> print(pop.subPopSizes())
   (10000,)

   now exiting runScriptInteractively...

`Download pedigreeMatingAgeStructured.py <pedigreeMatingAgeStructured.py>`_

The pedigree is then used to repeat the evolutionary process. However, because
some individuals were produced sexually using :class:`MendelianGenoTransmitter`
and some were copied using ``CloneGenoTransitter``, an :class:`IfElse` operator
has to be used to transmit genotypes correctly. This example uses the function
condition of the :class:`IfElse` operator and makes use of the fact that parent
``mom`` will be ``None`` if an individual is copied from his or her father.

plainnat simuPOP


