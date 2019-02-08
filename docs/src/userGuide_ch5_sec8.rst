Quantitative trait
==================

Quantitative traits are naturally stored in information fields of each
individual. A quantitative trait operator assigns quantitative trait fields
according to individual genetic (genotype) and environmental (other information
fields) information. Although a large number of quantitative trait models have
been used in theoretical and empirical studies, no model is popular enough to
deserve a specialized operator. Therefore, only one hybrid operator is currently
provided in simuPOP.


A hybrid quantitative trait operator (operator :class:`PyQuanTrait`)
--------------------------------------------------------------------

Operator :class:`PyQuanTrait` accepts a user defined function that returns
quantitative trait values for specified information fields. This operator can
comunicate with functions in one of the forms of ``func(geno)``, ``func(geno,
field_name, ...)`` or ``func(geno, field_name, gen)`` where ``field_name``
should be name of existing fields. simuPOP will pass genotype and value of
specified fields according to name of the passed function. Note that geno are
arrange locus by locus, namely in the order of ``A1``,``A2``,``B1``,``B2`` for
loci ``A`` and ``B``.

A quantitative trait operator can be applied before or after mating and assign
values to the trait fields of all parents or offspring, respectively. It can
also be applied during mating to assign trait values to offspring. Example
:ref:`PyQuanTrait <PyQuanTrait>` demonstrates the use of this operator, using
two trait fields ``trait1`` and ``trait2`` which are determined by individual
genotype and age. This example also demonstrates how to calculate statistics
within virtual subpopulations (defined by age).

.. _PyQuanTrait:

**Example**: *A hybrid quantitative trait model*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(size=5000, loci=2, infoFields=['qtrait1', 'qtrait2', 'age'])
   >>> pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[40]))
   >>> def qtrait(geno, age):
   ...     'Return two traits that depends on genotype and age'
   ...     return random.normalvariate(age * sum(geno), 10), random.randint(0, 10*sum(geno))
   ... 
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.2, 0.8]),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         # use random age for simplicity
   ...         sim.InitInfo(lambda:random.randint(20, 75), infoFields='age'),
   ...         sim.PyQuanTrait(loci=(0,1), func=qtrait, infoFields=['qtrait1', 'qtrait2']),
   ...         sim.Stat(meanOfInfo=['qtrait1'], subPops=[(0, sim.ALL_AVAIL)],
   ...             vars='meanOfInfo_sp'),
   ...         sim.PyEval(r"'Mean of trait1: %.3f (age < 40), %.3f (age >=40)\n' % "
   ...             "(subPop[(0,0)]['meanOfInfo']['qtrait1'], subPop[(0,1)]['meanOfInfo']['qtrait1'])"),
   ...     ],
   ...     gen = 5
   ... )
   Mean of trait1: 92.876 (age < 40), 183.515 (age >=40)
   Mean of trait1: 94.041 (age < 40), 183.374 (age >=40)
   Mean of trait1: 95.447 (age < 40), 183.288 (age >=40)
   Mean of trait1: 95.017 (age < 40), 183.919 (age >=40)
   Mean of trait1: 94.769 (age < 40), 185.430 (age >=40)
   5
   >>> 

   now exiting runScriptInteractively...

`Download PyQuanTrait.py <PyQuanTrait.py>`_


