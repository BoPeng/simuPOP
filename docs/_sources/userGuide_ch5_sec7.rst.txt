Penetrance
==========

Penetrance is the probability for an individual to be affected with a disease
conditioning on his or her genotype and other risk factors. A penetrance model
calculates such a probability for an individual and assign affection status
randomly according to this probability. For example, if an individual with
genotype ``10`` has probability 0.2 to be affected according to a penetrance
model, he or she will be affected with probability 0.2. Note that simuPOP
supports only one affection status. If there are multiple affection outcomes
involved, you can treat them as binary quantitative traits and use information
fields to store them.

A penetrance operator can be applied before or after mating, to assign affection
status to all individuals in the parental or offspring generation, respectively.
It can also be applied during mating and assign affection status to each
offspring. The latter could be used to assit natural selection through the
selection of offspring. You can also assign affection status to all individuals
in a population using the function form of a penetrance operator (e.g. function
``mapPenetrance`` for operator :class:`MapPenetrance`). Compared the penetrance
operators that assign affection status to only the current generation, **these
functions by default assign affection status to all ancestral generations as
well**.

A penetrance operator usually do not store the penetrance values. However, if an
information field is given, penetrance values will be saved to this information
field before it is used to determine individual affection status.


Map penetrance model (operator :class:`MapPenetrance`)
------------------------------------------------------

A map penetrance opertor uses a Python dictionary to provide penetrance values
for each type of genotype. For example, Example :ref:`MapPenetrance
<MapPenetrance>` uses a dictionary with keys ``(0,0)``, ``(0,1)`` and ``(1,1)``
to specify penetrance for individuals with these genotypes at locus 0.

.. _MapPenetrance:

**Example**: *A penetrance model that uses pre-defined fitness value*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=2000, loci=2)
   >>> sim.initGenotype(pop, freq=[.2, .8])
   >>> sim.mapPenetrance(pop, loci=0,
   ...     penetrance={(0,0):0, (0,1):.2, (1,1):.3})
   >>> sim.stat(pop, genoFreq=0, numOfAffected=1, vars='genoNum')
   >>> # number of affected individuals
   >>> pop.dvars().numOfAffected
   531
   >>> # which should be roughly (#01 + #10) * 0.2 + #11 * 0.3
   >>> (pop.dvars().genoNum[0][(0,1)] + pop.dvars().genoNum[0][(1,0)]) * 0.2 \
   ... + pop.dvars().genoNum[0][(1,1)] * 0.3
   514.2

   now exiting runScriptInteractively...

`Download MapPenetrance.py <MapPenetrance.py>`_

The above example assumes that penetrance for individuals with genotypes
``(0,1)`` and ``(1,0)`` are the same. This assumption is usually valid but can
be vialoated with impriting. In that case, you can specify fitness for both
types of genotypes. The underlying mechanism is that the :class:`MapPenetrance`
looks up a genotype in the dictionary first directly, and then without phase
information if a genotype is not found.

This operator supports haplodiploid populations and sex chromosomes. In these
cases, only valid alleles should be listed which can lead to dictionary keys
with different lengths. In addition, although less used because of potentially a
large number of keys, this operator can act on multiple loci. For example,

* keys ``(a1,a2)`` and ``(a1,)`` can be used to specify fitness values for
  female and male individuals in a haplodiploid population, respectively

* keys ``(x1,x2)`` and ``(x1,)`` can be used to specify fitness for female and
  male individuals according to a locus on the X chromosome in a diploid
  population, respectively. Similarly, keys ``()`` and ``(y,)`` for a locus on
  chromosome Y.

* keys ``(a1,a2,b1,b2)`` can be used to specify fitness values according to
  genotype at two loci in a diploid population.


Multi-allele penetrance model (operator :class:`MaPenetrance`)
--------------------------------------------------------------

A multi-allele penetrance model divides alleles into two groups, wildtype *A*
and mutants *a*, and treat alleles within each group as the same. The penetrance
model is therefore simplified to

* Two fitness values for genotype :math:`A`, :math:`a` in the haploid case

* Three fitness values for genotype *AA*, *Aa* and *aa* in the diploid single
  locus case. Genotype *Aa* and *aA* are assumed to have the same impact on
  fitness.

The default wildtype group contains allele 0 so the two allele groups are zero
and non-zero alleles. Example :ref:`MaPenetrance <MaPenetrance>` demonstrates
the use of this operator.

.. _MaPenetrance:

**Example**: *A multi-allele penetrance model*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(5000, loci=3)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.9] + [0.02]*5)
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.MaPenetrance(loci=0, penetrance=(0.01, 0.2, 0.3)),
   ...         sim.Stat(numOfAffected=True, vars='propOfAffected'),
   ...         sim.PyEval(r"'Gen: %d Prevalence: %.1f%%\n' % (gen, propOfAffected*100)"),
   ...     ],
   ...     gen = 5
   ... )
   Gen: 0 Prevalence: 4.4%
   Gen: 1 Prevalence: 4.4%
   Gen: 2 Prevalence: 4.7%
   Gen: 3 Prevalence: 4.4%
   Gen: 4 Prevalence: 4.3%
   5

   now exiting runScriptInteractively...

`Download MaPenetrance.py <MaPenetrance.py>`_

Operator :class:`MaPenetrance` also supports multiple loci by specifying fitness
values for all combination of genotype at specified loci. In the case of two
loci, this operator requires

* Four fitness values for genotype ``AB``, ``Ab``, ``aB`` and ``ab`` in the
  haploid case,

* Nine fitness values for genotype ``AABB``, ``AABb``, ``AAbb``, ``AaBB``,
  ``AaBb``, ``Aabb``, ``aaBB``, ``aaBb``, and ``aabb`` in the haploid case.

In general, :math:`2^{n}` values are needed for haploid populations and
:math:`3^{n}` values are needed for diploid populations where :math:`n` is the
number of loci. This operator does not yet support haplodiploid populations and
sex chromosomes.


Multi-loci penetrance model (operator :class:`MlPenetrance`)
------------------------------------------------------------

Although an individual's affection status can be affected by several factors,
each of which can be modeled individually, **only one penetrance value is used
to determine a person's affection status** and we have to use a multi-locus
penetrance model to combine single-locus models.

This multi-loci penetrance model applies several penetrance models to each
Individual and computes an overall penetrance value from the penetrance values
provided by these operators. Although this selector is designed to obtain multi-
loci penetrance values from several single-locus penetrance models, any
penetrance operator, including those obtain their penetrance values from
multiple disease predisposing loci, can be used in this operator. This operator
uses parameter ``mode`` to control how Individual penetrance values are
combined. More specifically, if :math:`f_{i}` are penetrance values obtained
from individual selectors, this selector returns

* :math:`\Pi_{i}f_{i}` if ``mode=MULTIPLICATIVE``, and

* :math:`\sum_{i}f_{i}` if ``mode=ADDITIVE``, and

* :math:`1-\Pi_{i}\left(1-f_{i}\right)` if ``mode=HETEROGENEITY``

0 or 1 will be returned if the returned fitness value is out of range of
``[0,1]``.

Example :ref:`MlPenetrance <MlPenetrance>` demonstrates the use of this operator
using an multiplicative multi-locus model over three additive single-locus
models at three diesease predisposing loci.

.. _MlPenetrance:

**Example**: *A multi-loci penetrance model*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(5000, loci=3)
   >>> sim.initGenotype(pop, freq=[0.2]*5)
   >>> # the multi-loci penetrance
   >>> sim.mlPenetrance(pop, mode=sim.MULTIPLICATIVE,
   ...     ops = [sim.MaPenetrance(loci=loc,
   ...         penetrance=[0, 0.3, 0.6]) for loc in range(3)])
   >>> # count the number of affected individuals.
   >>> sim.stat(pop, numOfAffected=True)
   >>> pop.dvars().numOfAffected
   542

   now exiting runScriptInteractively...

`Download MlPenetrance.py <MlPenetrance.py>`_


Hybrid penetrance model (operator :class:`PyPenetrance`)
--------------------------------------------------------

When your selection model involves multiple interacting genetic and
environmental factors, it might be easier to calculate a penetrance value
explicitly using a Python function. A hybrid penetrance operator can be used for
this purpose. If your penetrance model depends solely on genotype, you can
define a function such as

::

   def pfunc(geno):
       # calculate penetrance according to genotype at specified loci
       # in the order of A1,A2,B1,B2,C1,C2 for loci A,B,C (for diploid)
       return val

and use this function in an operator :class:`PySelector`\ (``func=pfunc,
loci=loci``). If your penetrance model depends on genotype as well as some
information fields, you can define a function in the form of

::

   def pfunc(geno, fields):
       # calculate penetrance according to genotype at specified loci
       # and values at specified informaton fields.
       return val

and use this function in an operator :class:`PySelector`\ (``func=pfunc,
loci=loci, paramFields=fields``). If the function you provide accepts three
arguments, :class:`PyPenetrance` will pass generation number as the third
argument so that you could implement generation-specific penetrance models (e.g.
``pfunc(geno, fields, gen)``).

When a :class:`PyPenetrance` operator is used to calculate penetrance for an
individual, it will collect his or her genotype at specified loci, optional
values at specified information fields, and the generation number to a user-
specified Python function, and take its return value as penetrance. As you can
imagine, the incorporation of information fields and generation number allow the
implementation of very complex penetrance scenarios such as gene environment
interaction and varying selection pressures. Note that this operator does not
pass sex and affection status to the user-defined function. If your selection
model is sex-dependent, you can define an information field ``sex``, synchronize
its value with individual sex (e.g. using operator :class:`InfoExec`\
``('sex=ind.sex()', exposeInd='ind'``) and pass this information to the user-
defined function (:class:`PySelector`\ (``func=func, paramFields='sex'``)).

Example :ref:`PySelector <PySelector>` demonstrates how to use a
:class:`PyPenetrance` to specify penetrance values according to a fitness table
and the smoking status of each individual. In this example, Individual risk is
doubled when he or she smokes. The disease prevalence is therefore much higher
in smokers than in non-smokers.

.. _PyPenetrance:

**Example**: *A hybrid penetrance model*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(size=2000, loci=[1]*2, infoFields=['p', 'smoking'])
   >>> pop.setVirtualSplitter(sim.InfoSplitter(field='smoking', values=[0,1]))
   >>> # the second parameter gen can be used for varying selection pressure
   >>> def penet(geno, smoking):
   ...     #     BB     Bb      bb
   ...     # AA  0.01   0.01    0.01
   ...     # Aa  0.01   0.03    0.03
   ...     # aa  0.01   0.03    0.05
   ...     #
   ...     # geno is (A1 A2 B1 B2)
   ...     if geno[0] + geno[1] == 1 and geno[2] + geno[3] != 0:
   ...         v = 0.03   # case of AaBb
   ...     elif geno[0] + geno[1] == 2 and geno[2] + geno[3] == 1:
   ...         v = 0.03   # case of aaBb
   ...     elif geno[0] + geno[1] ==2 and geno[2] + geno[3] == 2:
   ...         v = 0.05   # case of aabb
   ...     else:                
   ...         v = 0.01   # other cases
   ...     if smoking:
   ...         return v * 2
   ...     else:
   ...         return v
   ... 
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.5, .5]),
   ...         sim.PyOutput('Calculate prevalence in smoker and non-smokers\n'),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         # set smoking status randomly
   ...         sim.InitInfo(lambda : random.randint(0,1), infoFields='smoking'),
   ...         # assign affection status
   ...         sim.PyPenetrance(loci=[0, 1], func=penet),
   ...         sim.Stat(numOfAffected=True, subPops=[(0, sim.ALL_AVAIL)], 
   ...             vars='propOfAffected_sp', step=20),
   ...         sim.PyEval(r"'Non-smoker: %.2f%%\tSmoker: %.2f%%\n' % "
   ...             "(subPop[(0,0)]['propOfAffected']*100, subPop[(0,1)]['propOfAffected']*100)",
   ...             step=20)
   ...     ],
   ...     gen = 50
   ... )
   Calculate prevalence in smoker and non-smokers
   Non-smoker: 2.24%	Smoker: 4.52%
   Non-smoker: 2.29%	Smoker: 3.61%
   Non-smoker: 1.85%	Smoker: 3.80%
   50
   >>> 

   now exiting runScriptInteractively...

`Download PyPenetrance.py <PyPenetrance.py>`_


