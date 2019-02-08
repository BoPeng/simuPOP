Age structured populations with overlapping generations \*\*
============================================================

Age is an important factor in many applications because it is related to many
genetic (most obviously mating) and environmental factors that influence the
evolution of a population. The evolution of age structured populations will lead
to overlapping generations because parents can co-exist with their offspring in
such a population. Although simuPOP is based on a discrete generation model, it
can be used to simulate age structured populations.

To evolve an age structured population, you will need to

* Define an information field ``age`` and use it to store age of all
  individuals. Age is usally assigned randomly at the beginning of a simulation.

* Define a virtual splitter that splits the parental population into several
  virtual subpopulation. The most important VSP consists of mating individuals
  (e.g. individuals with age between 20 and 40). Advanced features of virtual
  splitters can be used to define complex VSPs such as males between age 20 - 40
  and females between age 15-30 (use a :class:`ProductSplitter` to split
  subpopulations by sex and age, and then a :class:`CombinedSplitter` to join
  several smaller VSPs together).

* Use a heterogeneous mating scheme that clones most individuals to the next
  generation (year) and produce offspring from the mating VSP.

Example :ref:`ageStructured <ageStructured>` gives an example of the evolution
of age-structured population.

* Information fields ``ind_id``, ``father_id`` and ``mother_id`` and operators
  :class:`IdTagger` and :class:`PedigreeTagger` are used to track pedigree
  information during evolution.

* A :class:`CloneMating` mating scheme is used to copy surviving individuals and
  a :class:`RandomMating` mating scheme is used to produce offspring.

* :class:`IdTagger` and :class:`PedigreeTagger` are used in the ``ops``
  parameter of :class:`RandomMating` because only new offspring should have a new
  ID and record parental IDs. If you use these operators in the ``duringOps``
  parameter of the ``evolve`` function, individuals copied by :class:`CloneMating`
  will have a new ID, and a missing parental ID.

* The resulting population is age-structured so Pedigrees could be extracted
  from such a population.

* The penetrance function is age dependent. Because this penetrance function is
  applied to all individuals at each year and an individual will have the disease
  once he or she is affected, this penetrance function is more or less a hazard
  function.

.. _ageStructured:

**Example**: *Example of the evolution of age-structured population.*

::

   >>> import simuPOP as sim
   >>> import random
   >>> N = 10000
   >>> pop = sim.Population(N, loci=1, infoFields=['age', 'ind_id', 'father_id', 'mother_id'])
   >>> pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[20, 50, 75]))
   >>> def demoModel(gen, pop):
   ...     '''A demographic model that keep a constant supply of new individuals'''
   ...     # number of individuals that will die
   ...     sim.stat(pop, popSize=True, subPops=[(0,3)])
   ...     # individuals that will be kept, plus some new guys.
   ...     return pop.popSize() - pop.dvars().popSize + N // 75
   ... 
   >>> def pene(geno, age, ind):
   ...     'Define an age-dependent penetrance function'
   ...     # this disease does not occur in children
   ...     if age < 16:
   ...         return 0
   ...     # if an individual is already affected, keep so
   ...     if ind.affected():
   ...         return 1
   ...     # the probability of getting disease increases with age
   ...     return (0., 0.001*age, 0.001*age)[sum(geno)]
   ... 
   >>> def outputstat(pop):
   ...     'Calculate and output statistics'
   ...     sim.stat(pop, popSize=True, numOfAffected=True,
   ...         subPops=[(0, sim.ALL_AVAIL)],
   ...         vars=['popSize_sp', 'propOfAffected_sp'])
   ...     for sp in range(3):
   ...         print('%s: %.3f%% (size %d)' % (pop.subPopName((0,sp)),
   ...             pop.dvars((0,sp)).propOfAffected * 100.,
   ...             pop.dvars((0,sp)).popSize))
   ...     #
   ...     return True
   ... 
   >>> 
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         # random assign age
   ...         sim.InitInfo(lambda: random.randint(0, 75), infoFields='age'),
   ...         # random genotype
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...         # assign an unique ID to everyone.
   ...         sim.IdTagger(),
   ...         sim.PyOutput('Prevalence of disease in each age group:\n'),
   ...     ],
   ...     # increase the age of everyone by 1 before mating.
   ...     preOps=sim.InfoExec('age += 1'),
   ...     matingScheme=sim.HeteroMating([
   ...         # all individuals with age < 75 will be kept. Note that
   ...         # CloneMating will keep individual sex, affection status and all
   ...         # information fields (by default).
   ...         sim.CloneMating(subPops=[(0,0), (0,1), (0,2)], weight=-1),
   ...         # only individuals with age between 20 and 50 will mate and produce
   ...         # offspring. The age of offspring will be zero.
   ...         sim.RandomMating(ops=[
   ...             sim.IdTagger(),                   # give new born an ID
   ...             sim.PedigreeTagger(),             # track parents of each individual
   ...             sim.MendelianGenoTransmitter(),   # transmit genotype
   ...         ],
   ...         numOffspring=(sim.UNIFORM_DISTRIBUTION, 1, 3),
   ...         subPops=[(0,1)]),],
   ...         subPopSize=demoModel),
   ...     # number of individuals?
   ...     postOps=[
   ...         sim.PyPenetrance(func=pene, loci=0),
   ...         sim.PyOperator(func=outputstat, step=20)
   ...     ],
   ...     gen = 200
   ... )
   Prevalence of disease in each age group:
   age < 20: 0.578% (size 2596)
   20 <= age < 50: 2.649% (size 4002)
   50 <= age < 75: 4.217% (size 3249)
   age < 20: 0.526% (size 2660)
   20 <= age < 50: 27.627% (size 3931)
   50 <= age < 75: 50.317% (size 3313)
   age < 20: 0.489% (size 2660)
   20 <= age < 50: 28.470% (size 3927)
   50 <= age < 75: 61.757% (size 3347)
   age < 20: 0.639% (size 2660)
   20 <= age < 50: 29.449% (size 3990)
   50 <= age < 75: 62.384% (size 3246)
   age < 20: 0.526% (size 2660)
   20 <= age < 50: 27.694% (size 3990)
   50 <= age < 75: 64.030% (size 3325)
   age < 20: 0.865% (size 2660)
   20 <= age < 50: 28.070% (size 3990)
   50 <= age < 75: 60.782% (size 3325)
   age < 20: 0.489% (size 2660)
   20 <= age < 50: 29.624% (size 3990)
   50 <= age < 75: 60.812% (size 3325)
   age < 20: 0.526% (size 2660)
   20 <= age < 50: 29.273% (size 3990)
   50 <= age < 75: 61.714% (size 3325)
   age < 20: 0.789% (size 2660)
   20 <= age < 50: 27.769% (size 3990)
   50 <= age < 75: 61.233% (size 3325)
   age < 20: 0.639% (size 2660)
   20 <= age < 50: 29.073% (size 3990)
   50 <= age < 75: 59.669% (size 3325)
   200
   >>> 
   >>> # draw two Pedigrees from the last age-structured population
   >>> from simuPOP import sampling
   >>> sample = sampling.drawNuclearFamilySample(pop, families=2, numOffspring=(2,3),
   ...     affectedParents=(1,2), affectedOffspring=(1,3))
   >>> sim.dump(sample)
   Ploidy: 2 (diploid)
   Chromosomes:
   1:  (AUTOSOME, 1 loci)
      (1)
   Information fields: 
   age ind_id father_id mother_id 
   population size: 8 (1 subpopulations with 8 Individuals)
   Number of ancestral populations: 0

   SubPopulation 0 (), 8 Individuals:
      0: MA 1 | 0 |  37 31578 27047 27596
      1: MU 1 | 0 |  29 32638 29986 29012
      2: MA 1 | 0 |  37 31579 27047 27596
      3: FA 1 | 0 |  57 29012 25317 22955
      4: MU 0 | 0 |  49 29986 27087 25888
      5: FA 1 | 1 |  67 27596 24124 24202
      6: FA 1 | 0 |  29 32637 29986 29012
      7: MA 1 | 0 |  71 27047 23653 20932

   >>> 

   now exiting runScriptInteractively...

`Download ageStructured.py <ageStructured.py>`_


