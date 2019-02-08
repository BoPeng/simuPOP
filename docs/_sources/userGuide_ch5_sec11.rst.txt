Statistics calculation (operator :class:`Stat`)
===============================================


How statistics calculation works
--------------------------------

A :class:`Stat` operator calculates specified statistics of a population when it
is applied to this population. This operator can be applied to specified
replicates (parameter *rep*) at specified generations (parameter *begin*, *end*,
*step*, and *at*). This operator does not produce any output (ignore parameter
*output*) after statistics are calculated. Instead, it stores results in the
local namespace of the population being applied. Other operators can retrieve
these variables or evalulate expression directly in this local namespace.

The :class:`Stat` operator is usually used in conjunction with a :class:`PyEval`
or :class:`PyExec` operator which execute Python statements and/or expressions
in a population's local namespace. For example, operators  ::

   ops = [
       Stat(alleleFreq=[0]),
       PyEval("'%.2f' % alleleFreq[0][0]")
   ]

in the ``ops`` parameter of the :meth:`Simulator.evolve` function will be
applied to populations during evolution. The first operator calculates allele
frequency at the first locus and store the results in each population's local
namespace. The second operator formats and outputs one of the variables. Because
of the flexiblity of the :class:`PyEval` operator, you can output statistics,
even simple derived statistics, in any format. For example, you can output
expected heterozygosity (:math:`1-\sum p_{i}^{2}`) using calculated allele
frequencies as follows::

   PyEval("'H_exp=%.2f' % (1-sum([x*x for x in alleleFreq[0].values()]))")

Note that ``alleleFreq[0]`` is a dictionary.

You can also retrieve variables in a population directly using functions
:meth:`Population.vars`\ () or :meth:`Population.dvars`\ (). The only difference
between these functions is that ``vars`` returns a dictionary and :func:`dvars`\
() returns a Python object that uses variable names as attributes
(``vars()['alleleFreq']`` is equivalent to :func:`dvars`\ ``.alleleFreq``). This
method is usually used when the function form of the :class:`Stat` operator is
used. For example,  ::

   stat(pop, alleleFreq=[0])
   H_exp = 1 - sum([x*x for x in pop.dvars().alleleFreq[0].values()])

uses the ``stat`` function (note the capital S) to count frequencies of alleles
for a given population and calculates expected heterozygosity using these
variables.


:class:`defdict` datatype
-------------------------

simuPOP uses dictionaries to save statistics such as allele frequencies. For
example, ``alleleFreq[5]`` can be ``{0:0.2, 3:0.8}`` meaning there are 20%
allele 0 and 80% allele 3 at locus 5 in a population. However, because it is
sometimes unclear whether or not a particular allele exists in a population,
``alleleFreq[5][allele]`` can fail with a ``KeyError`` exception if
``alleleFreq[5]`` does not have key ``allele``.

To address this problem, a special default dictionary type :class:`defdict` is
used for dictionaries with keys determined from a population. This derived
dictionary type works just like a regular dictionay, but it returns 0, instead
of raising a ``KeyError`` exception, when an invalid key is used. For example,
subpopulations in Example :ref:`defdictType <defdictType>` have different
alleles. Although ``pop.dvars(sp).alleleFreq[0]`` have only two keys for
``sp=0`` or ``1``, ``pop.dvars(sp).alleleFreq[0][x]`` are used to print
frequencies of alleles ``0``, ``1`` and ``2``.

.. _defdictType:

**Example**: *The defdict datatype*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([100]*2, loci=1)
   >>> sim.initGenotype(pop, freq=[0, 0.2, 0.8], subPops=0)
   >>> sim.initGenotype(pop, freq=[0.2, 0.8], subPops=1)
   >>> sim.stat(pop, alleleFreq=0, vars=['alleleFreq_sp'])
   >>> for sp in range(2):
   ...     print('Subpop %d (with %d alleles): ' % (sp, len(pop.dvars(sp).alleleFreq[0])))
   ...     for a in range(3):
   ...         print('%.2f ' % pop.dvars(sp).alleleFreq[0][a])
   ... 
   Subpop 0 (with 2 alleles): 
   0.00 
   0.21 
   0.79 
   Subpop 1 (with 2 alleles): 
   0.21 
   0.79 
   0.00 

   now exiting runScriptInteractively...

`Download defdict.py <defdict.py>`_

.. note::

   The standard ``collections`` module of Python has a ``defaultdict`` type that
   accepts a default factory function that will be used when an invalid key is
   encountered. The :class:`defdict` type is similar to ``defaultdict(int)`` but
   with an important difference: when an invalid key is encountered, ``d[key]``
   with a default value will be inserted to a ``defaultdict(int)``, but will not be
   inserted to a :class:`defdict`. That is to say, it is safe to use
   ``alleleFreq[loc].keys()`` to get available alleles after non-assignment
   ``alleleFreq[loc][allele]`` operations.


Support for virtual subpopulations
----------------------------------

The :class:`Stat` operator supports parameter *subPops* and can calculate
statistics in specified subpopulations. For example  ::

   Stat(alleleFreq=[0], subPops=[(0, 0), (1, 0)])

will calculate the frequencies of alleles at locus 0, among Individuals in two
virtual subpopulations. If the virtual subpopulation is defined by sex (using a
:class:`SexSplitter`), the above operator will calculate allele frequency among
all males in the first and second subpopulations (not separately!). If
``subPops`` is not specified, allele frequency of the whole population (all
subpopulations) will be calculated.

Although many statistics could be calculated and outputted, the :class:`Stat`
operator by default outputs a selected number of variables for each statisic
calculated. Other statistics could be calculated and outputted if their names
are specified in parameter ``vars``. Variable names ending with ``_sp`` is
interpreted as variables that will be calculated and outputted in all or
specified (virtual) subpopulations. For example, parameter ``vars`` in    ::

   Stat(alleleFreq=[0], subPops=[0, (1, 0)], vars=['alleleFreq_sp', 'alleleNum_sp'])

tells this operator to output numbers and frequencies of alleles at locus ``0``
in subpopulation ``0`` and virtual subpopulation ``(1,0)``. These variables will
be saved in dictionaries ``subPop[sp]`` of the local namespace. For example, the
above operator will write variables such as ``subPop[0]['alleleFreq'],
subPop[(1,0)]['alleleFreq'] and subPop[(1,0)]['alleleNum']``. Functions
:meth:`Population.vars`\ (``sp``) and :meth:`Population.dvars`\ (``sp``) are
provided as shortcuts to access these variables but the full variable names have
to be specified if these variables are used in expressions.

By default, the same variables will be set for a statistic, regardless of the
values of the ``loci`` and ``subPops`` parameter. This can be a problem if
multiple :class:`Stat` operators are used to calculate the same statistics for
different sets of loci (e.g. for each chromosome) or subpopulations. To avoid
name conflict, you can use parameter *suffix* to add a suffix to all variables
outputted by a Stat operator. For example, Example :ref:`statSuffix
<statSuffix>` uses 4 :class:`Stat` operators to calculate overall and pairwise
:math:`F_{ST}` values for three subpopulations. Different suffixes are used for
pairwise :math:`F_{ST}` estimators so that variables set by these operators will
not override each other.

.. _statSuffix:

**Example**: *Add suffixes to variables set by multiple Stat operators*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([5000]*3, loci=5)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(structure=range(5), subPops=(0, 1), suffix='_01', step=40),
   ...         sim.Stat(structure=range(5), subPops=(1, 2), suffix='_12', step=40),
   ...         sim.Stat(structure=range(5), subPops=(0, 2), suffix='_02', step=40),
   ...         sim.Stat(structure=range(5), step=40),
   ...         sim.PyEval(r"'Fst=%.3f (pairwise: %.3f %.3f %.3f)\n' % (F_st, F_st_01, F_st_12, F_st_02)",
   ...             step=40),
   ...     ],
   ...     gen = 200
   ... )
   Fst=0.000 (pairwise: 0.000 0.000 0.000)
   Fst=0.004 (pairwise: 0.006 0.003 0.004)
   Fst=0.012 (pairwise: 0.017 0.015 0.004)
   Fst=0.008 (pairwise: 0.012 0.010 0.001)
   Fst=0.008 (pairwise: 0.007 0.009 0.007)
   200

   now exiting runScriptInteractively...

`Download statSuffix.py <statSuffix.py>`_

.. note::

   The :class:`Stat` opeartor accepts overlapping or even duplicate virtual
   subpopulations. During the calculation of summary statistics, these
   subpopulations are treated as separate subpopulations so some individuals can be
   counted more than once. For example, individuals in virtual subpopulation (0, 1)
   will be counted twice during the calculation of allele frequency and population
   size in operator   ::

      Stat(alleleFreq=[0], popSize=True, subPops=[0, (0, 1)])


Counting individuals by sex and affection status
------------------------------------------------

Parameters *popSize*, *numOfMales* and *numOfAffected* provide basic Individual
counting statistics. They count the number of all, male/female,
affected/unaffected individuals in all or specified (virtual) subpopulations,
and set variables such as ``popSize``, ``numOfMales``, ``numOfFemales``,
``numOfAffected``, ``numOfUnaffected``. Proportions and statistics for
subpopulations are available if variables such as ``propOfMales``,
``numOfAffected_sp`` are specified in parameter vars. Another variable
``subPopSize`` is defined for parameter ``popSize=True``. It is a list of sizes
of all or specified subpopulations and is easier to use than referring to
variable ``popSize`` from individual subpopulations.

Example :ref:`statCount <statCount>` demonstrates how to use these parameters in
operator :class:`Stat`. It defines four VSPs by sex and affection status (using
a ``stackedSplitter``) and count individuals by sex and affection status. It is
worth noting that ``pop.dvars().popSize`` in the first example is the total
number of individuals in two virtual subpopulations ``(0,0)`` and ``(0,2)``,
which are all male indiviudals, and all unaffected individuals. Because these
two VSPs overlap, this variable can be larger than actual population size.

.. _statCount:

**Example**: *Count individuals by sex and/or affection status*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(10000, loci=1)
   >>> pop.setVirtualSplitter(sim.CombinedSplitter(
   ...     [sim.SexSplitter(), sim.AffectionSplitter()]))
   >>> sim.initSex(pop)
   >>> sim.initGenotype(pop, freq=[0.2, 0.8])
   >>> sim.maPenetrance(pop, loci=0, penetrance=[0.1, 0.2, 0.5])
   >>> # Count sim.population size
   >>> sim.stat(pop, popSize=True, subPops=[(0, 0), (0, 2)])
   >>> # popSize is the size of two VSPs, does not equal to total sim.population size.
   >>> # Because two VSPs overlap (all males and all unaffected), popSize can be
   >>> # greater than real sim.population size.
   >>> print(pop.dvars().subPopSize, pop.dvars().popSize)
   [5052, 6080] 11132
   >>> # print popSize of each virtual subpopulation.
   >>> sim.stat(pop, popSize=True, subPops=[(0, 0), (0, 2)], vars='popSize_sp')
   >>> # Note the two ways to access variable in (virtual) subpopulations.
   >>> print(pop.dvars((0,0)).popSize, pop.dvars().subPop[(0,2)]['popSize'])
   5052 6080
   >>> # Count number of male (should be the same as the size of VSP (0,0).
   >>> sim.stat(pop, numOfMales=True)
   >>> print(pop.dvars().numOfMales)
   5052
   >>> # Count the number of affected and unaffected male individual
   >>> sim.stat(pop, numOfMales=True, subPops=[(0, 2), (0, 3)], vars='numOfMales_sp')
   >>> print(pop.dvars((0,2)).numOfMales, pop.dvars((0,3)).numOfMales)
   3056 1996
   >>> # or number of affected male and females
   >>> sim.stat(pop, numOfAffected=True, subPops=[(0, 0), (0, 1)], vars='numOfAffected_sp')
   >>> print(pop.dvars((0,0)).numOfAffected, pop.dvars((0,1)).numOfAffected)
   1996 1924
   >>> # These can also be done using a sim.ProductSplitter...
   >>> pop.setVirtualSplitter(sim.ProductSplitter(
   ...     [sim.SexSplitter(), sim.AffectionSplitter()]))
   >>> sim.stat(pop, popSize=True, subPops=[(0, x) for x in range(4)])
   >>> # counts for male unaffected, male affected, female unaffected and female affected
   >>> print(pop.dvars().subPopSize)
   [3056, 1996, 3024, 1924]

   now exiting runScriptInteractively...

`Download statCount.py <statCount.py>`_


Number of segregating and fixed sites
-------------------------------------

Parameter *numOfSegSites* counts the number of segregating sites for specified
or all loci, for all individuals or individuals in specified (virtual)
subpopulations. It can also be used to count the number of fixed sites . This
parameter sets variables ``numOfSegSites`` and ``numOfFixedSites``. Here we
defined fixed sites as loci with only one non-zero allele (e.g. fixed to a non-
zero allele). Other numbers, such as all loci with only one allele (including
zero), or loci with all wildtype alleles (only zero), can be derived from these
two counts. Starting from version 1.1.3, variables ``segSites`` and
``fixedSites`` can be used to return a list of segregating and fixed sites.

For example, Example :ref:`numSegSites <numSegSites>` demonstrates how to use
this operator to calculate the number of segregating sites (sites with alleles 0
and 1), number of fixed sites (sites with only allele 1), and number of loci
with only wildtype alleles (loci with only allele 0). As you can see, the
population starts with 100 segregating sites. During evolution, alleles at some
loci get lost and some get fixed, and there should be no segregating site if we
evolve the population for long enough.

.. _numSegSites:

**Example**: *Count number of segregating and fixed sites*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(100, loci=[1]*100)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.3, 0.7]),
   ...         sim.PyOutput('#all 0\t#seg sites\t#all 1\n'),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(numOfSegSites=sim.ALL_AVAIL,
   ...             vars=['numOfSegSites', 'numOfFixedSites']),
   ...         sim.PyEval(r'"%d\t%d\t%d\n" % (100-numOfSegSites-numOfFixedSites,'
   ...             'numOfSegSites, numOfFixedSites)',
   ...             step=50)
   ...         ],
   ...     gen=500
   ... )
   #all 0	#seg sites	#all 1
   0	100	0
   0	93	7
   3	76	21
   7	55	38
   12	40	48
   17	31	52
   19	23	58
   22	19	59
   26	14	60
   28	10	62
   500
   >>> # output a list of segregating sites
   >>> sim.stat(pop, numOfSegSites=sim.ALL_AVAIL, vars='segSites')
   >>> print(pop.dvars().segSites)
   [11, 15, 20, 32, 39, 43, 44, 51, 86, 95]

   now exiting runScriptInteractively...

`Download statNumOfSegSites.py <statNumOfSegSites.py>`_


Allele count and frequency
--------------------------

Parameter *alleleFreq* accepts a list of markers at which allele frequencies in
all or specified (virtual) subpopulations will be calculated. This statistic
sets variables ``alleleFreq[loc][allele]`` and ``alleleNum[loc][allele]`` which
are frequencies and numbers of allele ``allele`` at locus ``loc``, respectively.
If variables ``alleleFreq_sp`` and ``alleleNum_sp`` are specified in parameter
*vars*, these variables will be set for all or specified (virtual)
subpopulations. **At the Python level, these variables are dictionaries of
default dictionaries.** That is to say, ``alleleFreq[loc]`` at a unspecified
locus will raise a ``KeyError`` exception, and ``alleleFreq[loc][allele]`` of an
invalid allele will return 0.

Example :ref:`statAlleleFreq <statAlleleFreq>` demonstrates an advanced usage of
allele counting statistic. In this example, two virtual subpopulations are
defined by individual affection status. During evolution, a multi-allele
penetrance operator is used to determine individual affection status and a
:class:`Stat` operator is used to calculate allele frequencies in these two
virtual subpopulations, and in the whole population. Because the simulated
disease is largely caused by the existence of allele 1 at the first locus, it is
expected that the frequency of allele 1 is higher in the case group than in the
control group. It is worth noting that ``alleleFreq[0][1]`` in this example is
the frequency of allele 1 in the whole population because these two virtual
subpopulations add up to the whole population.

.. _statAlleleFreq:

**Example**: *Calculate allele frequency in affected and unaffected individuals*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(10000, loci=1)
   >>> pop.setVirtualSplitter(sim.AffectionSplitter())
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(loci=0, freq=[0.8, 0.2])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.MaPenetrance(penetrance=[0.1, 0.4, 0.6], loci=0),
   ...         sim.Stat(alleleFreq=0, subPops=[(0, 0), (0, 1)],
   ...             vars=['alleleFreq', 'alleleFreq_sp']),
   ...         sim.PyEval(r"'Gen: %d, freq: %.2f, freq (aff): %.2f, freq (unaff): %.2f\n' % " + \
   ...             "(gen, alleleFreq[0][1], subPop[(0,1)]['alleleFreq'][0][1]," + \
   ...             "subPop[(0,0)]['alleleFreq'][0][1])"),
   ...     ],
   ...     gen = 5
   ... )
   Gen: 0, freq: 0.20, freq (aff): 0.41, freq (unaff): 0.14
   Gen: 1, freq: 0.20, freq (aff): 0.40, freq (unaff): 0.14
   Gen: 2, freq: 0.20, freq (aff): 0.41, freq (unaff): 0.14
   Gen: 3, freq: 0.20, freq (aff): 0.41, freq (unaff): 0.14
   Gen: 4, freq: 0.19, freq (aff): 0.41, freq (unaff): 0.14
   5

   now exiting runScriptInteractively...

`Download statAlleleFreq.py <statAlleleFreq.py>`_


Genotype count and frequency
----------------------------

Parameter *genoFreq* accepts a list of loci at which genotype counts and
frequencies are calculated and outputted. A genotype is represented as a tuple
of alleles at a locus. The length of the tupples** **is determined by the number
of homologous copy of chromosomes in a population. For example, genotypes in a
diploid population are ordered pairs such as ``(1, 2)`` where 1 and 2 are
alleles at a locus on, respectively, the first and second homologous copies of
chromosomes. ``(1, 2)`` and ``(2, 1)`` are different genotypes. This statistic
sets dictionaries (with locus indexes as keys) of default dictionaries (with
genotypes as keys) ``genoFreq`` and ``genoNum``.

Example :ref:`statGenoFreq <statGenoFreq>` creates a small population and
initializes a locus with rare alleles 0, 1 and a common allele 2. A function
``stat`` (the function form of operator :class:`Stat`) is used to count the
available genotypes. Note that ``pop.dvars().genoFreq[0][(i,j)]`` can be used to
print frequencies of all genotypes even when not all genotypes are available in
the population.

.. _statGenoFreq:

**Example**: *Counting genotypes in a population*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(100, loci=[1, 1, 1], lociNames=['A', 'X', 'Y'],
   ...     chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y])
   >>> sim.initGenotype(pop, freq=[0.01, 0.05, 0.94])
   >>> sim.stat(pop, genoFreq=['A', 'X']) # both loci indexes and names can be used.
   >>> print('Available genotypes on autosome:', list(pop.dvars().genoFreq[0].keys()))
   Available genotypes on autosome: [(0, 2), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]
   >>> for i in range(3):
   ...     for j in range(3):
   ...         print('%d-%d: %.3f' % (i, j, pop.dvars().genoFreq[0][(i,j)]))
   ... 
   0-0: 0.000
   0-1: 0.000
   0-2: 0.020
   1-0: 0.000
   1-1: 0.030
   1-2: 0.070
   2-0: 0.010
   2-1: 0.040
   2-2: 0.830
   >>> print('Genotype frequency on chromosome X:\n', \
   ...     '\n'.join(['%s: %.3f' % (x,y) for x,y in pop.dvars().genoFreq[1].items()]))
   Genotype frequency on chromosome X:
    (0,): 0.020
   (1,): 0.030
   (2,): 0.950

   now exiting runScriptInteractively...

`Download statGenoFreq.py <statGenoFreq.py>`_


Homozygote and heterozygote count and frequency
-----------------------------------------------

In a diploid population, a heterozygote is a genotype with two different alleles
and a homozygote is a genotype with two identical alleles. Parameter
``heteroFreq`` accepts a list of loci and outputs variables ``heteroFreq`` which
is a dictionary of heterozygote frequencies at specfied loci. Optional variables
``heteroNum``, ``homoFreq`` and ``homoNum`` can be outputted for all and each
(virtual) subpopulations. Example :ref:`statHeteroFreq <statHeteroFreq>`
demonstrates the decay of heterozygosity of a locus due to genetic drift.

.. _statHeteroFreq:

**Example**: *Counting homozygotes and heterozygotes in a population*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(100, loci=1)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(heteroFreq=0, step=10),
   ...         sim.PyEval(r"'Gen: %d, HeteroFreq: %.2f\n' % (gen, heteroFreq[0])", step=20)
   ...     ],
   ...     gen = 100
   ... )
   Gen: 0, HeteroFreq: 0.45
   Gen: 20, HeteroFreq: 0.44
   Gen: 40, HeteroFreq: 0.55
   Gen: 60, HeteroFreq: 0.46
   Gen: 80, HeteroFreq: 0.40
   100

   now exiting runScriptInteractively...

`Download statHeteroFreq.py <statHeteroFreq.py>`_


Haplotype count and frequency
-----------------------------

Haplotypes refer to alleles on the same homologous copy of a chromosome at
specified loci. For example, an diploid individual can have haplotypes ``(0, 2,
1)`` and ``(0, 1, 1)`` at loci ``(2, 3, 5)`` if he or she has genotype ``(0,
0)``, ``(2, 1)`` and ``(1,1)`` at loci 2, 3 and 5 respectively. Parameter
*haploFreq* accept one or more lists of loci specifying one or more haplotype
sites (e.g. ``haploFreq=[(0,1,2), (2,3)]`` specifies two haplotype sites). The
results are saved to dictionaries (with haplotype site as keys) of default
dictionaries (with haplotype as keys). For example,
``haploFreq[(0,1,2)][(0,1,1)]`` will be the frequency of haplotype ``(0, 1, 1)``
at loci ``(0, 1, 2)``. Example :ref:`statHaploFreq <statHaploFreq>` prints the
numbers of genotypes and haplotypes at loci 0, 1 and 2 of a small population.
Note that the ``viewVars`` function defined in module ``simuUtil`` can make use
of a wxPython window to view all variables if it is called in GUI mode.

.. _statHaploFreq:

**Example**: *Counting haplotypes in a population*

::

   >>> import simuPOP as sim
   >>> from simuPOP.utils import viewVars
   >>> pop = sim.Population(100, loci=3)
   >>> sim.initGenotype(pop, freq=[0.2, 0.4, 0.4], loci=0)
   >>> sim.initGenotype(pop, freq=[0.2, 0.8], loci=2)
   >>> sim.stat(pop, genoFreq=[0, 1, 2], haploFreq=[0, 1, 2],
   ...     vars=['genoNum', 'haploFreq'])
   >>> viewVars(pop.vars(), gui=False)
   {'genoNum': {0: {(0, 0): 3.0,
                    (0, 1): 7.0,
                    (0, 2): 5.0,
                    (1, 0): 9.0,
                    (1, 1): 14.0,
                    (1, 2): 16.0,
                    (2, 0): 8.0,
                    (2, 1): 14.0,
                    (2, 2): 24.0},
                1: defdict({(0, 0): 100.0}),
                2: {(0, 0): 4.0,
                    (0, 1): 19.0,
                    (1, 0): 15.0,
                    (1, 1): 62.0}},
    'haploFreq': {(0, 1, 2): {(0, 0, 0): 0.03,
                              (0, 0, 1): 0.145,
                              (1, 0, 0): 0.055,
                              (1, 0, 1): 0.315,
                              (2, 0, 0): 0.125,
                              (2, 0, 1): 0.33}}}

   now exiting runScriptInteractively...

`Download statHaploFreq.py <statHaploFreq.py>`_

.. note::

   *haploFreq* does not check if loci in a haplotype site belong to the same
   chromosome, or if loci are duplicated or in order. It faithfully assemble
   alleles at specified loci as haplotypes although these haplotypes might not be
   biologically meaningful.

.. note::

   Counting a large number of haplotypes on long haplotype sites may exhaust the
   RAM of your computer.


Summary statistics of information fields
----------------------------------------

Parameter ``sumOfInfo``, ``meanOfInfo``, ``varOfInfo``, ``maxOfInfo`` and
``minOfInfo`` are used to calculate the sum, mean, sample variance
(:math:`\frac{1}{n-1}\sum_{i=1}^{n}\left(x_{i}-\bar{x}\right)^{2}`), max and min
of specified information fields of individuals in all or specified (virtual)
subpopulations. The results are saved in dictionaries ``sumOfInfo``,
``meanOfInfo``, ``varOfInfo``, ``maxOfInfo`` and ``minOfInfo`` with information
fields as keys. For example, parameter  ``meanOfInfo='age'`` calculates the mean
age of all individuals and set variable ``meanOfInfo['age']``.

Example :ref:`statInfo <statInfo>` demonstrates a mixing process of two
populations. The population starts with two types of individuals with ancestry
values 0 or 1 (information field ``anc``). During the evolution, parents mate
randomly and the ancestry of offspring is the mean of parental ancestry values.
A :class:`Stat` operator is used to calculate the mean and variance of
individual ancestry values, and the number of individuals in five ancestry
groups. It is not surprising that whereas population mean ancestry does not
change, more and more people have about the same number of ancestors from each
group and have an ancestry value around 0.5. The variance of ancestry values
therefore decreases gradually.

.. _statInfo:

**Example**: *Calculate summary statistics of information fields*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population([500], infoFields='anc')
   >>> # Defines VSP 0, 1, 2, 3, 4 by anc.
   >>> pop.setVirtualSplitter(sim.InfoSplitter('anc', cutoff=[0.2, 0.4, 0.6, 0.8]))
   >>> #
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         # anc is 0 or 1
   ...         sim.InitInfo(lambda : random.randint(0, 1), infoFields='anc')
   ...     ],
   ...     matingScheme=sim.RandomMating(ops=[
   ...         sim.MendelianGenoTransmitter(),
   ...         sim.InheritTagger(mode=sim.MEAN, infoFields='anc')
   ...     ]),
   ...     postOps=[
   ...         sim.Stat(popSize=True, meanOfInfo='anc', varOfInfo='anc',
   ...             subPops=[(0, sim.ALL_AVAIL)]),
   ...         sim.PyEval(r"'Anc: %.2f (%.2f), #inds: %s\n' %" + \
   ...             "(meanOfInfo['anc'], varOfInfo['anc'], " + \
   ...             "', '.join(['%4d' % x for x in subPopSize]))")
   ...     ],
   ...     gen = 5,
   ... )
   Anc: 0.51 (0.12), #inds:  118,    0,  251,    0,  131
   Anc: 0.51 (0.06), #inds:   27,  121,  190,  137,   25
   Anc: 0.52 (0.03), #inds:   14,  143,  138,  181,   24
   Anc: 0.52 (0.02), #inds:    4,   85,  267,  137,    7
   Anc: 0.52 (0.01), #inds:    0,   40,  385,   75,    0
   5

   now exiting runScriptInteractively...

`Download statInfo.py <statInfo.py>`_


Linkage disequilibrium
----------------------

Parameter *LD* accepts a list of loci-pairs (e.g. ``LD=[(0,1),(2,3)]``) with
optional primary alleles at two loci (e.g. ``LD=[(0,1,0,0),(2,3)]``). For each
pair of loci, this operator calculates linkage disequilibrium and optional
association measures between them.

Assuming that two loci are both diallelic, one with alleles :math:`A` and
:math:`a`, and the other with alleles :math:`B` and :math:`b`. If we denote
:math:`P_{x}`, :math:`P_{xy}` as allele and haplotype frequencies for allele
:math:`x` and haplotype :math:`xy`, respectively, the linkage disequilibrium
measures **with respect to primaries alleles** *A* and *B* are

* Basic LD measure :math:`D`:

  .. math::

      D=P_{AB}-P_{A}P_{B}

  *D* ranges from -0.25 to 0.25. The sign depends on the choice of alleles (*A*
  and *B*) at two loci.

* Lewontin's :math:`D'=D/D_{max}` where

  .. math::

      D_{max}=\begin{cases}
      \min\left(P_{A}\left(1-P_{B}\right),\left(1-P_{A}\right)P_{B}\right) & \textrm{if }D>0\\
      \min\left(P_{A}P_{B},\left(1-P_{A}\right)\left(1-P_{B}\right)\right) & \textrm{if }D<0
      \end{cases}

  *D'* ranges from -1 to 1. The sign depends on the choice of alleles (*A* and
  *B*) at two loci.

* :math:`r^{2}` (:math:`\Delta^{2}` in Devlin1995)

  .. math::

      r^{2}=\frac{D^{2}}{P_{A}\left(1-P_{A}\right)P_{B}\left(1-P_{B}\right)}

If one or both loci have more than 2 alleles, or if no primary allele is
specified, the LD measures are calculated as follows:

* If primary alleles are specified, all other alleles are considered as minor
  alleles with combined frequency (e.g. :math:`1-P_{A}`). The same formulas apply
  which lead to signed :math:`D` and :math:`D'` measures.

* If primary alleles are not specified, these LD measures are calculated as the
  average of the absolute value of diallelic measures of all allele pairs. For
  example, the multi-allele version of :math:`r^{2}` is

  .. math::

      r^{2}=\sum_{i}\sum_{j}P_{i}P_{j}\left|r_{ij}^{2}\right|=\sum_{i}\sum_{j}\frac{D_{ij}^{2}}{\left(1-P_{i}\right)\left(1-P_{j}\right)}

  where :math:`i` and :math:`j` iterate through all alleles at the two loci. **In
  the diallelic case, LD measures will be the absolute value of the single
  measures** because :math:`D_{ij}` and :math:`D'_{ij}` only differ by signs.

In another word,

* ``LD=[loc1, loc2]`` will yield positive :math:`D` and :math:`D'` measures.

* ``LD=[loc1, loc2, allele1, allele2]`` will yield signed :math:`D` and
  :math:`D'` measures.

* In the diallelic case, both cases yield identical results except for signs of
  :math:`D` and :math:`D'`.

* In the multi-allelic case, the results can be different because ``LD=[loc1,
  loc2, allele1, allele2]`` combines non-primary alleles and gives a single
  diallelic measure.

.. note::

   A large number of linkage disequilibrium measures have been used in different
   disciplines but not all of them are well-accepted. Requests of adding a
   particular LD measure will be considered when a reliable reference is provided.

Association tests between specified loci could also be calculated using a
:math:`m` by :math:`n` table of haplotype frequencies. If primary alleles are
specified, non-primary alleles are combined to form a 2 by 2 table
(:math:`m=n=2`). Otherwise, :math:`m` and :math:`n` are respective numbers of
alleles at two loci.

* :math:`\chi^{2}` and its :math:`p`\ -value (variable ``LD_ChiSq`` and
  ``LD_ChiSq_p``, respectively). A one-side :math:`\chi^{2}` test with
  :math:`\left(m-1\right)\times\left(n-1\right)` degrees of freedom will be used.

* Cramer V statistic (variable ``CramerV``):

  .. math::

      V=\sqrt{\frac{\chi^{2}}{N\times\mbox{min}\left(m-1,n-1\right)}}

  where :math:`N` equals the total number of haplotypes
  (:math:`2\times\mbox{popSize}` for autosomes in diploid populations).

This statistic sets variables ``LD``, ``LD_prime``, ``R2``, and optionally
``ChiSq``, ``ChiSq_p`` and ``CramerV``. SubPopulation specific variables can be
calculated by specifying variables such as ``LD_sp`` and ``R2_sp``. Example
:ref:`statLD <statLD>` demonstrates how to calculate various LD measures and
output selected variables. Note that the significant overall LD between two loci
is an artifact of population structure because loci are in linkage equilibrium
in each subpopulation.

.. _statLD:

**Example**: *Linkage disequilibrium measures*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([1000]*2, loci=3)
   >>> sim.initGenotype(pop, freq=[0.2, 0.8], subPops=0)
   >>> sim.initGenotype(pop, freq=[0.8, 0.2], subPops=1)
   >>> sim.stat(pop, LD=[[0, 1, 0, 0], [1, 2]],
   ...     vars=['LD', 'LD_prime', 'R2', 'LD_ChiSq', 'LD_ChiSq_p', 'CramerV',
   ...         'LD_prime_sp', 'LD_ChiSq_p_sp'])
   >>> from pprint import pprint
   >>> pprint(pop.vars())
   {'CramerV': {0: defdict({1: 0.3355834766347789}),
                1: defdict({2: 0.39144946095755695})},
    'LD': {0: defdict({1: 0.08387987499999999}),
           1: defdict({2: 0.09783043749999992})},
    'LD_ChiSq': {0: defdict({1: 450.4650791611408}),
                 1: defdict({2: 612.9307219358476})},
    'LD_ChiSq_p': {0: defdict({1: 0.0}), 1: defdict({2: 0.0})},
    'LD_prime': {0: defdict({1: 0.3425347836362625}),
                 1: defdict({2: 0.4057999832524774})},
    'R2': {0: defdict({1: 0.1126162697902852}),
           1: defdict({2: 0.15323268048396166})},
    'subPop': {0: {'LD_ChiSq_p': {0: defdict({1: 0.03843990070970382}),
                                  1: defdict({2: 0.5110492462003573})},
                   'LD_prime': {0: defdict({1: -0.17661111690962444}),
                                1: defdict({2: 0.016760924318107204})}},
               1: {'LD_ChiSq_p': {0: defdict({1: 0.8024214035646771}),
                                  1: defdict({2: 0.11685510935577492})},
                   'LD_prime': {0: defdict({1: -0.02259456714902688}),
                                1: defdict({2: 0.035632559660018596})}}}}

   now exiting runScriptInteractively...

`Download statLD.py <statLD.py>`_


Genetic association
-------------------

Genetic association refers to association between individual genotype (alleles
or genotype) and phenotype (affection status). There are a large number of
statistics tests based on different study designs (e.g. case-control, Pedigree,
longitudinal) with different covariate variables. Although specialized software
applications should be used for sophisticated statistical analysis, simuPOP
provides a number of simple genetic association tests for convenience. These
tests

* Are single-locus tests that test specified loci separately.

* Are based on individual affection status. Associations between genotype and
  quantitative traits are currently unsupported.

* Apply to all individuals in specified (virtual) subpopulations. Because a
  population usually has much more unaffected individuals than affected ones, it
  is a common practice to draw certain types of samples (e.g. a case-control
  sample with the same number of cases and controls) before statistical tests are
  applied.

simuPOP currently supports the following tests:

* **Allele-based Chi-square test**: This is the basic allele-based
  :math:`\chi^{2}` test that can be applied to diploid as well as haploid
  populations. Basically, a 2 by :math:`n` contigency table is set up for each
  locus with :math:`n_{ij}` being the number of alleles :math:`j` in cases
  :math:`\left(i=0\right)` and controls :math:`\left(i=1\right)`. A
  :math:`\chi^{2}` test is applied to each locus and set variables
  ``Allele_ChiSq`` and ``Allele_ChiSq_p`` to the :math:`\chi^{2}` statistic and
  its two-sided :math:`p` value (with degrees freedom :math:`n-1`). Note that
  genotype information is not preserved in such a test.

* **Genotype-based Chi-square test**: This is the genotype-based
  :math:`\chi^{2}` test for diploid populations. Basically, a 2 by :math:`n`
  contigency table is set up for each locus with :math:`n_{ij}` being the number
  of genotype :math:`j` (unordered pairs of alleles) in cases
  :math:`\left(i=0\right)` and controls :math:`\left(i=1\right)`. A
  :math:`\chi^{2}` test is applied to each locus and set variables ``Geno_ChiSq``
  and ``Geno_ChiSq_p`` to the :math:`\chi^{2}` statistic and its two-sided
  :math:`p` value (with degrees freedom :math:`n-1`). This test is usually applied
  to diallelic loci with 3 genotypes (*AA*, *Aa* and *aa*) but it can be applied
  to loci with more than two alleles as well.

* **Genotype-based trend test**: This Cochran-Armitage test can only be applied
  to diallelic loci in diploid populations. For each locus, a 2 by 3 contigency
  table is set up with :math:`n_{ij}` being the number of genotype :math:`j`
  (*AA*, *Aa* and *aa* with *A* being the wildtype allele) in cases
  :math:`\left(i=0\right)` and controls :math:`\left(i=1\right)`. A Cochran-
  Armitage trend test is applied to each locus and set variables ``Armitage_p`` to
  its two-sided :math:`p` value.

Example :ref:`statAssociation <statAssociation>` demonstrates how to apply a
penetrance model, draw a case-control sample and apply genetic association tests
to an evolving population. In this example, a penetrance model is applied to a
locus (locus 3). A Python operator is then used to draw a case-control sample
from the population and test genetic association at two surrounding loci.
Because these two loci are tightly linked to the disease predisposing locus,
they are in strong association with the disease initially. However, because of
recombination, such association decays with time at rates depending on their
genetic distances to the disease predisposing locus.

.. _statAssociation:

**Example**: *Genetic association tests*

::

   >>> import simuPOP as sim
   >>> from simuPOP.utils import *
   >>> from simuPOP.sampling import drawCaseControlSample
   >>> def assoTest(pop):
   ...     'Draw case-control sample and apply association tests'
   ...     sample = drawCaseControlSample(pop, cases=500, controls=500)
   ...     sim.stat(sample, association=(0, 2), vars=['Allele_ChiSq_p', 'Geno_ChiSq_p', 'Armitage_p'])
   ...     print('Allele test: %.2e, %.2e, Geno test: %.2e, %.2e, Trend test: %.2e, %.2e' \
   ...         % (sample.dvars().Allele_ChiSq_p[0], sample.dvars().Allele_ChiSq_p[2],
   ...         sample.dvars().Geno_ChiSq_p[0], sample.dvars().Geno_ChiSq_p[2],
   ...         sample.dvars().Armitage_p[0], sample.dvars().Armitage_p[2]))
   ...     return True
   ... 
   >>> pop = sim.Population(size=100000, loci=3)
   >>> pop.setVirtualSplitter(sim.ProportionSplitter([0.5, 0.5]))
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[0]*3, subPops=[(0,0)]),
   ...         sim.InitGenotype(genotype=[1]*3, subPops=[(0,1)]),
   ...     ],
   ...     matingScheme=sim.RandomMating(ops=sim.Recombinator(loci=[0, 1], rates=[0.01, 0.005])),
   ...     postOps=[
   ...         sim.MaPenetrance(loci=1, penetrance=[0.1, 0.2, 0.4]),
   ...         sim.PyOperator(func=assoTest, step=20),
   ...     ],
   ...     gen = 100
   ... )
   Allele test: 0.00e+00, 0.00e+00, Geno test: 0.00e+00, 0.00e+00, Trend test: 0.00e+00, 0.00e+00
   Allele test: 1.14e-13, 4.44e-16, Geno test: 3.09e-13, 2.66e-15, Trend test: 7.66e-14, 2.22e-16
   Allele test: 1.71e-08, 8.55e-15, Geno test: 4.95e-08, 3.45e-13, Trend test: 1.62e-08, 7.36e-14
   Allele test: 8.57e-09, 7.99e-15, Geno test: 3.09e-08, 2.18e-14, Trend test: 7.05e-09, 2.66e-15
   Allele test: 3.12e-06, 9.05e-09, Geno test: 5.95e-06, 8.83e-08, Trend test: 2.12e-06, 1.26e-08
   100

   now exiting runScriptInteractively...

`Download statAssociation.py <statAssociation.py>`_


population structure
--------------------

Parameter ``structure`` measures the structure of a population using the
following statistics:

* The :math:`G_{ST}` statistic developed by Nei Nei1973. This statistic is
  equivalent to Wright's fixation index :math:`F_{ST}` in the diallelic case so it
  can be considered as the multi-allele and multi-locus extension of Wright's
  :math:`F_{ST}`. It assumes known genotype frequency so it can be used to
  calculate true :math:`F_{ST}` of a population when all genotype information is
  available. This statistic sets a dictionary of locus level :math:`G_{ST}`
  (variable ``g_st``) and a summary statistics for all loci (variable ``G_st``).

* Wright's fixation index :math:`F_{ST}` calculated using an algorithm developed
  by Weir1984. This statistic considers existing populations as random samples
  from an infinite pool of populations with the same ancestral population so it is
  best to be applied to random samples where true genotype frequencies are
  unknown. This statistic sets dictionaries of locus level :math:`F_{ST}`,
  :math:`F_{IT}` and :math:`F_{IS}` (variables ``f_st``, ``f_is`` and ``f_it``),
  and summary statistics for all loci (variables ``F_st``, ``F_is`` and ``F_it``)
  . When hetergozygote count is unavailable (non-diploid population, loci on sex
  chromosomes and mitochondrial chromosomes), simuPOP uses expected heterozygosity
  to estimate this quantity.

These statistics by default uses all existing subpopulations, but it can also be
applied to a subset of subpopulations, or even virtual subpopulations using
parameter *subPops*. That is to say, you can measure the genetic difference
between males and females using ``subPops=[(0,0), (0,1)]`` if a SexSplitter is
used to define two virtual subpopulations with male and female individuals
respectively.

Example :ref:`statStructure <statStructure>` demonstrate a simulation with two
replicates. In the first replicate, three subpopulations evolve separately
without migration and become more and more genetically distinct. In the second
replicate, a low level migration is applied between subpopulations so the
population structure is kept at a low level.

.. _statStructure:

**Example**: *Measure of population structure*

::

   >>> import simuPOP as sim
   >>> from simuPOP.utils import migrIslandRates
   >>> simu = sim.Simulator(sim.Population([5000]*3, loci=10, infoFields='migrate_to'),
   ...     rep=2)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     preOps=sim.Migrator(rate=migrIslandRates(0.01, 3), reps=1),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(structure=range(10), step=40),
   ...         sim.PyEval("'Fst=%.3f (rep=%d without migration) ' % (F_st, rep)", step=40, reps=0),
   ...         sim.PyEval("'Fst=%.3f (rep=%d with migration) ' % (F_st, rep)", step=40, reps=1),
   ...         sim.PyOutput('\n', reps=-1, step=40)
   ...     ],
   ...     gen = 200
   ... )
   Fst=0.000 (rep=0 without migration) Fst=0.000 (rep=1 with migration) 
   Fst=0.003 (rep=0 without migration) Fst=0.002 (rep=1 with migration) 
   Fst=0.006 (rep=0 without migration) Fst=0.002 (rep=1 with migration) 
   Fst=0.008 (rep=0 without migration) Fst=0.003 (rep=1 with migration) 
   Fst=0.010 (rep=0 without migration) Fst=0.001 (rep=1 with migration) 
   (200, 200)

   now exiting runScriptInteractively...

`Download statStructure.py <statStructure.py>`_


Hardy-Weinberg equilibrium test
-------------------------------

Parameter ``HWE`` accepts a list of loci at which exact Hardy Weinberg
equilibrium tests are applied. The *p*-values of the tests are assigned to a
dictionary ``HWE``. Example :ref:`statHWE <statHWE>` demonstrates how Hardy
Weinberg equilibrium is reached in one generation.

.. _statHWE:

**Example**: *Hardy Weinberg Equilibrium test*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([1000], loci=1)
   >>> pop.setVirtualSplitter(sim.ProportionSplitter([0.4, 0.4, 0.2]))
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[0,0], subPops=[(0,0)]),
   ...         sim.InitGenotype(genotype=[0,1], subPops=[(0,1)]),
   ...         sim.InitGenotype(genotype=[1,1], subPops=[(0,2)]),
   ...     ],
   ...     preOps=[
   ...         sim.Stat(HWE=0, genoFreq=0),
   ...         sim.PyEval(r'"HWE p-value: %.5f (AA: %.2f, Aa: %.2f, aa: %.2f)\n" % (HWE[0], '
   ...             'genoFreq[0][(0,0)], genoFreq[0][(0,1)] + genoFreq[0][(1,0)], genoFreq[0][(1,1)])'),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(HWE=0, genoFreq=0),
   ...         sim.PyEval(r'"HWE p-value: %.5f (AA: %.2f, Aa: %.2f, aa: %.2f)\n" % (HWE[0], '
   ...             'genoFreq[0][(0,0)], genoFreq[0][(0,1)] + genoFreq[0][(1,0)], genoFreq[0][(1,1)])'),
   ...     ],
   ...     gen = 1
   ... )
   HWE p-value: 0.00000 (AA: 0.40, Aa: 0.40, aa: 0.20)
   HWE p-value: 0.93636 (AA: 0.38, Aa: 0.48, aa: 0.15)
   1

   now exiting runScriptInteractively...

`Download statHWE.py <statHWE.py>`_


Measure of Inbreeding
---------------------

Inbreeding coefficient at a generation is defined as the probability that the
two alleles in a given individual are identical by decent (IBD). Although it is
usually very difficult to estimate this quantity, it is easy to observe it
directly during evolution if the ancestors of alleles are tracked. This can be
done using the lineage module of simuPOP where allelic lineage is tracked during
evolution. For example, Example :ref:`statIBD <statIBD>` output the frequency of
IBD loci in a population of size 500. It also outputs the frequency of IBS
(Identical by State), which should always be larger than IBD frequency, and
theoretical estimate of the decay of inbreeding coefficient.

.. _statIBD:

**Example**: *Frequency of IBD as a measure of inbreeding coefficient*

::

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='lineage')
   >>> import simuPOP as sim
   >>> pop = sim.Population([500], loci=[1]*100)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitLineage(),
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.2]*5),
   ...     ],
   ...     preOps=[
   ...         sim.Stat(inbreeding=sim.ALL_AVAIL, popSize=True, step=10),
   ...         sim.PyEval(r'"gen %d: IBD freq %.4f, IBS freq %.4f, est: %.4f\n" % '
   ...             '(gen, sum(IBD_freq.values()) /len(IBD_freq), '
   ...             ' sum(IBS_freq.values()) /len(IBS_freq), '
   ...             ' 1 - (1-1/(2.*popSize))**gen)', step=10)
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     gen = 100
   ... )
   gen 0: IBD freq 0.0000, IBS freq 0.1994, est: 0.0000
   gen 10: IBD freq 0.0084, IBS freq 0.2072, est: 0.0100
   gen 20: IBD freq 0.0167, IBS freq 0.2142, est: 0.0198
   gen 30: IBD freq 0.0266, IBS freq 0.2204, est: 0.0296
   gen 40: IBD freq 0.0380, IBS freq 0.2292, est: 0.0392
   gen 50: IBD freq 0.0486, IBS freq 0.2383, est: 0.0488
   gen 60: IBD freq 0.0577, IBS freq 0.2457, est: 0.0583
   gen 70: IBD freq 0.0689, IBS freq 0.2566, est: 0.0676
   gen 80: IBD freq 0.0782, IBS freq 0.2616, est: 0.0769
   gen 90: IBD freq 0.0887, IBS freq 0.2638, est: 0.0861
   100

   now exiting runScriptInteractively...

`Download statIBD.py <statIBD.py>`_


Effective population size
-------------------------

Effective population size is an important, yet complicated concept in population
genetics. Simply put, the effective population size is determined by a mating
scheme, namely how parents are selected and how offsprings are generated. In the
context of forward-time simulation, if we populate an offspring population from
a parental population, a true effective population size can be calculated, under
certain assumptions, as

.. math::

      Ne=\frac{kN-1}{k-1+V_{k}/k}

where :math:`k` and :math:`V_{k}` are the mean and variance of the number of
gametes each parent  transmits to the offspring generation. Naturally, the
number of sex chromosomes transmitted will be different for males and females.
This effective size is independent of genotypes and is called the demographic
effective size.

Because the calculation of demographic effective size needs to track which
alleles are transmitted from parental to offspring population, it has to collect
information from both parental and offspring populations, and can only be
calculated using the lineage modules of simuPOP. As shown in Example
:ref:`statNeDemographic <statNeDemographic>`, a :class:`Stat` operator is
applied before mating to mark lineage of alleles of each locus with an
individual index, and save the IDs of parents in a variable ``Ne_demo_base``.
After mating, another :class:`Stat` operator is used to count how many alleles
each parent has contributed to the offspring generation, and calculate
demographic effective size accordingly. This example uses three virtual
subpopulations, a whole subpopulation, all male individuals, and all female
individuals, and calculated effective size for loci on an autosome, an X
chromosome, and a Y chromosome. As we can imagine, the effective size is 0 at
the Y chromosome for all females, because no such chromsome is transmitted from
the parental population.

.. _statNeDemographic:

**Example**: *Demographic effective population size*

::

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='lineage', quiet=True)
   >>> import simuPOP as sim
   >>> pop = sim.Population([2000], loci=[1]*3,
   ...     chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y])
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.3, 0.7]),
   ...     ],
   ...     preOps=[
   ...         sim.Stat(effectiveSize=range(3), subPops=[0, (0,0), (0,1)],
   ...             vars='Ne_demo_base_sp'),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(effectiveSize=range(3), subPops=[0, (0,0), (0,1)],
   ...             vars='Ne_demo_sp'),
   ...         sim.PyEval(r'"Demographic Ne: %.1f (auto), %.1f (X), %.1f (Y), '
   ...             r'Males: %.1f, %.1f, %.1f, Females: %.1f, %.1f, %.1f\n"'
   ...             '% tuple([subPop[0]["Ne_demo"][x] for x in (0, 1, 2)] + '
   ...             '[subPop[(0,0)]["Ne_demo"][x] for x in (0, 1, 2)] + '
   ...             '[subPop[(0,1)]["Ne_demo"][x] for x in (0, 1, 2)])')
   ...     ],
   ...     gen = 5
   ... )
   Demographic Ne: 2021.2 (auto), 1808.8 (X), 1056.1 (Y), Males: 1038.4, 1049.4, 1056.1, Females: 983.8, 983.8, nan
   Demographic Ne: 2024.8 (auto), 1886.4 (X), 918.2 (Y), Males: 965.7, 1014.2, 918.2, Females: 1063.3, 1063.3, nan
   Demographic Ne: 2048.7 (auto), 1858.5 (X), 969.2 (Y), Males: 1023.0, 1037.4, 969.2, Females: 1025.1, 1025.1, nan
   Demographic Ne: 1955.0 (auto), 1790.6 (X), 956.8 (Y), Males: 958.8, 985.2, 956.8, Females: 996.5, 996.5, nan
   Demographic Ne: 2000.5 (auto), 1811.7 (X), 955.1 (Y), Males: 983.8, 966.2, 955.1, Females: 1016.8, 1016.8, nan
   5

   now exiting runScriptInteractively...

`Download statNeDemographic.py <statNeDemographic.py>`_

Effective population sizes could also be estimated from genotypes because
changes of genotypes reflects properties of the mating scheme. However, it is
important to realize that **evolving a population for one generation is only one
realization of many possible realizations of the same mating scheme** (effective
size). If we consider the demographic effective size as the average effective
size of all realizations, estimating effective size from genotypes will be
inaccurate unless a large number of unlinked loci are used. The temporal methods
essentially try to get better estimate by averaging such realizations across
multiple generations, although the demographic effective size might vary due to
change of population size.

simuPOP currently provides two temporal methods proposed by Waples (1989) and
Jorde & Ryman's (2007). Because these methods estimate effective population size
using changes of allele frequencies of samples at two generations, it is
necessary to set a baseline generation before any temporal method could be
applied.

The baseline information is saved to variable ``Ne_temporal_base`` when this
variable is specified in the ``vars`` parameter of the :class:`Stat` operator.
After the baseline is set, for example, at generation 0, if the operator
:class:`Stat` is applied at generations 0, 20, and 40, it will set variable
``Ne_waples89_P1``, ``Ne_waples89_P2``\ (for Waples 1989) and ``Ne_tempoFS_P1``,
``Ne_tempoFS_P``\ 2 (for Jorde & Ryman 2007, as implemented in a package
``TempoFS``) as the census population size at generation 0, estimated effective
population sizes between generation 0 and 20 at generation 20, and estimates
between 0 and 40 at generation 40. The variables are lists of three elements:
the estimated Ne and lower and upper boundaries of the 95% confidence interval.

Sampling plan 1 assumes that samples are drawn with replacement at the first
time point so that some of the individuals sampled in the first time period
could have contributed genes to subsequent generations (see Nei and Tajima, 1981
Genetics and other papers). simuPOP uses census population (or subpopulation if
the statistics are calcuated for each subpopulations) size as :math:`N` and
consider the sample being a subset of the population (or subpopulation), it
should be applied to a virtual subpopulation (e.g. a subset of individuals
defined by a :class:`RangeSplitter`) of the whole population. Sample plan 2
treats the sample as a sample from an infinitely-sized population, and should be
applied to a population (sample) that is actually extracted from a larger
population. Results under both assumptions are calculated and provided so you
should choose the ones that match your sampling plan.

Example :ref:`statNeTemporal <statNeTemporal>` demonstrates how to calculate
temporal effective population sizes at a 20 generation interval during
evolution, using a fixed baseline generation at generation 0. The statistics are
estimated from genotypes at 50 unlinked loci from 500 random samples from a
population of size 2000. Instead of drawing random samples explicitly, this
example defines a virtual subpopulation that consists of the first 500
individuals in the population. The Stat operator is applied at generations 0,
20, 40, ..., 100 to this virtual subpopulation, with the first output being the
census size (of the sample). Because a standard Wright-Fisher random mating
scheme is used, the true effective population size should be around 2000. It
would be interesting to adjust this evolutionary process (with population
expansion, with varying number of offspring etc) and the method of estimation
(sample size, generations between estimates) to see how well this statistic
estimate effective population size under different scenarios.

.. _statNeTemporal:

**Example**: *Temporal effective population size using a fixed baseline sample*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([2000], loci=[1]*50)
   >>> pop.setVirtualSplitter(sim.RangeSplitter([0, 500]))
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.3, 0.7]),
   ...         sim.Stat(effectiveSize=range(50), subPops=[(0,0)],
   ...             vars='Ne_temporal_base'),
   ...     ],
   ...     preOps=[
   ...         sim.Stat(effectiveSize=range(50), subPops=[(0,0)],
   ...             vars=['Ne_waples89_P1', 'Ne_tempoFS_P1'], step=20),
   ...         sim.PyEval(r'"Waples Ne: %.1f (%.1f - %.1f), TempoFS: '
   ...             r'%.1f (%.1f - %.1f), at generation %d\n" % '
   ...             'tuple(Ne_waples89_P1 + Ne_tempoFS_P1 + [gen])', step=20)
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     gen = 101
   ... )
   Waples Ne: 500.0 (500.0 - 500.0), TempoFS: 500.0 (500.0 - 500.0), at generation 0
   Waples Ne: 1853.1 (1155.2 - 3536.1), TempoFS: 1843.2 (1255.1 - 3467.7), at generation 20
   Waples Ne: 1537.9 (979.7 - 2452.6), TempoFS: 1565.7 (1117.0 - 2617.2), at generation 40
   Waples Ne: 1843.3 (1178.0 - 2872.4), TempoFS: 1963.4 (1332.2 - 3730.9), at generation 60
   Waples Ne: 1783.0 (1143.4 - 2710.7), TempoFS: 1807.2 (1291.5 - 3008.7), at generation 80
   Waples Ne: 1572.7 (1011.2 - 2346.6), TempoFS: 1639.5 (1205.1 - 2563.6), at generation 100
   101

   now exiting runScriptInteractively...

`Download statNeTemporal.py <statNeTemporal.py>`_

Instead of using a fixed baseline generation, it is also possible to reset
baseline generation during evolution. For example, Example :ref:`statNeInterval
<statNeInterval>` demonstrates how to calculate temporal effective population
sizes at a 20 generation interval during evolution. This example sets variable
``Ne_temporal_base`` with ``Ne_waples89_P1`` whenever the Stat operator is
applied. This effectively resets the baseline generation to the present
generation at generations 0, 20, 40, etc, so baseline generations 0, 20, 40, ...
are used at generations 20, 40, .... This example also demonstrates how to use
the suffix parameter to apply the same statistics with different parameters.

.. _statNeInterval:

**Example**: *Temporal effective population size between consecutive samples*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([2000], loci=[1]*50)
   >>> pop.setVirtualSplitter(sim.RangeSplitter([0, 500]))
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.3, 0.7]),
   ...         sim.Stat(effectiveSize=range(50), subPops=[(0,0)],
   ...             vars='Ne_temporal_base'),
   ...     ],
   ...     preOps=[
   ...         sim.Stat(effectiveSize=range(50), subPops=[(0,0)], 
   ...             vars='Ne_waples89_P1', step=20),
   ...         sim.Stat(effectiveSize=range(50), subPops=[(0,0)], step=20,
   ...             suffix='_i', vars=['Ne_temporal_base', 'Ne_waples89_P1']),
   ...         sim.PyEval(r'"Waples Ne (till %d): %.1f (%.1f - %.1f), '
   ...             r'(interval) %.1f (%.1f - %.1f)\n" % '
   ...             'tuple([gen] + Ne_waples89_P1 + Ne_waples89_P1_i)',
   ...             step=20)
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     gen = 101
   ... )
   Waples Ne (till 0): 500.0 (500.0 - 500.0), (interval) 500.0 (500.0 - 500.0)
   Waples Ne (till 20): 1853.1 (1155.2 - 3536.1), (interval) 1853.1 (1155.2 - 3536.1)
   Waples Ne (till 40): 1537.9 (979.7 - 2452.6), (interval) 2063.7 (1281.1 - 4094.1)
   Waples Ne (till 60): 1843.3 (1178.0 - 2872.4), (interval) 1681.9 (1052.1 - 3112.9)
   Waples Ne (till 80): 1783.0 (1143.4 - 2710.7), (interval) 1872.7 (1167.0 - 3586.3)
   Waples Ne (till 100): 1572.7 (1011.2 - 2346.6), (interval) 2056.1 (1276.6 - 4073.3)
   101

   now exiting runScriptInteractively...

`Download statNeInterval.py <statNeInterval.py>`_

Linkage disequilibrium method is another popular method to estimate effective
population size. Compared to temporal methods, it has the distinct advantage
that it requires only one sample. simuPOP provides a method that is developed by
Waples in his 2006 paper. To use this method, you will need to specify variable
``Ne_LD`` for a random mating scheme, or ``Ne_LD_mono`` for a monogamous mating
scheme. :ref:`statNeLD <statNeLD>` demonstrates this usage. Note that because
the LDNe mehod is sensitive to rare alleles (which can lead to inflated measure
of LD), simuPOP provides estimates that ignores alleles with frequencies less
than 0 (all alleles are kept), 0.01, 0.02 and 0.05. The results are saved in
variable ``Ne_LD`` as a dictionary with keys 0, 0.01, 0.02, 0.05, and values as
lists of estimated effective population sizes and their 95% confidence
intervals. Because of the existence of many rare alleles, the example gives
quite different estimates with and without rare alleles (using cutoff=0.02).

.. _statNeLD:

**Example**: *Effective population size estimated using a LD based method*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([2000], loci=[1]*50)
   >>> pop.setVirtualSplitter(sim.RangeSplitter([0, 500]))
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.005]*4 + [0.015]*2 + [0.25, 0.7]),
   ...     ],
   ...     preOps=[
   ...         sim.Stat(effectiveSize=sim.ALL_AVAIL, subPops=[(0,0)], 
   ...             vars='Ne_LD', step=20),
   ...         sim.PyEval(r'"LD Ne (gen %d): %.1f (%.1f - %.1f)'
   ...             r', %.1f (%.1f - %.1f, adjusted)\n" % '
   ...             'tuple([gen] + Ne_LD[0.] + Ne_LD[0.02])',
   ...             step=20)
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     gen = 101
   ... )
   LD Ne (gen 0): 30623.2 (5220.9 - inf), inf (8071.2 - inf, adjusted)
   LD Ne (gen 20): 6297.4 (2574.4 - inf), 1900.0 (1160.3 - 4647.8, adjusted)
   LD Ne (gen 40): 2187.6 (1554.1 - 3589.2), 2535.5 (1459.2 - 8173.8, adjusted)
   LD Ne (gen 60): 2757.8 (1799.2 - 5619.3), 3510.9 (1801.6 - 32066.7, adjusted)
   LD Ne (gen 80): 2574.0 (1729.7 - 4828.9), 1813.2 (1197.7 - 3501.7, adjusted)
   LD Ne (gen 100): 3234.6 (1819.5 - 12210.9), 2834.8 (1603.4 - 10168.4, adjusted)
   101

   now exiting runScriptInteractively...

`Download statNeLD.py <statNeLD.py>`_

simuPOP allows you to estimate effective population size using genotypes at
selected loci from selected individuals. It is up to you, however, to decide
when to apply the operator (pre- or post-mating), how to draw samples, and
select the right method for your data. For example, the temporal methods assume
discrete generations and no (or slight) selection, migration, and mutation. The
LD method assumes that markers are selectively neutral and independent;
population has discrete generations and is closed to immigration; and sampling
is random. In addition, to keep the interface simple, simuPOP does not provide
many options as dedicated programs do (e.g. TempoFS). Please export your samples
in other formats (e.g. use operator ``Export(format=''GENEPOP'')`` or function
``export(pop, format=''GENEPOP'')`` from module ``simuPOP.utiles``) and use
these programs if you need such flexibilities.


Other statistics
----------------

If you need other statistics, a popular approach is to define them using Python
operators. If your statistics is based on existing statistics such as allele
frequency, it is a good idea to calculate existing statistics using a ``stat``
function and derive your statistics from population variables. Please refer to
the last chapter of this guide on an example.

If you would like to calculate some summary statistics that involves individual
information fields but cannot be calculated using parameters such as minOfInfo,
you can try to use operators such as InfoExec to process individuals one by one
and collect result. For example, you can use operators   ::

   PyExec('s=0')
   InfoExec('s+=x*x')
   PyEval('s')

to calculate and report :math:`s=\sum x^{2}` where x is an information field
during evolution. This makes use of the fact that operator :class:`InfoExec`
goes through all individuals and evaluate the statement.

If performance becomes a problem, you might want to have a look at the source
code of simuPOP and implement your statistics at the C++ level. If you believe
that your statistics are popular enough, please send your implementation to the
simuPOP mailinglist for possible inclusion of your statistics into simuPOP.


Support for sex and customized chromosome types
-----------------------------------------------

simuPOP supports statistics calculation for loci on sex chromosomes. For
example, when pair-wise difference between haplotypes is calculated using
parameter ``neutrality``, it will pick the right haplotypes for X, and Y
chromosomes. However, because neutrality is calculated based on a group of
haplotypes of all specified loci, even if the loci are collected across
chromosomes, you can not use operator   ::

   Stat(neutrality=ALL_AVAIL)

if the loci are selected from chromosomes of different types, because different
numbers of haplotypes exists on these chromosomes. To calculate ``Pi`` for these
chromosomes, you would have to calculate them separately, using operators such
as  ::

   Stat(neutrality=range(30,40), suffix='_X')
   Stat(neutrality=range(40,50), suffix='_Y')

so that all specified loci are on the same type of chromosomes. Here we use
parameter ``suffix`` to avoid conflict of variable names because both operator
would produce the same variable ``Pi`` without this parameter.

The case with customized chromosomes are more complex because the meaning of
these chromosomes are defined by users. If these chromosomes are mitochondrial
DNAs, only chromosomes from the females are carrying useful information. If you
would like to calculate, for example, the ``Pi`` statistics for these
chromosomes, you will have to explicitly selected females for calculation. This
can be done by operator  ::

   Stat(neutrality=range(50,60), vsps=[(ALL_AVAIL, 'FEMALE')], suffix='_mt')

if VSPs have been created by a :class:`SexSplitter`.

Example :ref:`statChromTypes <statChromTypes>` demonstrates the use of these
operators. This example intentionally initializes all individuals with the same
haplotypes on all chromosomes (the :class:`InitGenotype` operator ignores
chromosome types). Because of different chromosome types, four :class:`Stat`
operators are used to get the ``Pi`` statistics for them. These operators return
different results because different sets of haplotypes are picked for the
calculation of this statistics.

.. _statChromTypes:

**Example**: *Statistics for sex and customized chromosome types*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, loci=[5]*4,
   ...     chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y, sim.MITOCHONDRIAL])
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(haplotypes=[ [0, 1, 2, 0, 1]*4, [2, 1, 0, 2, 3]*4 ],
   ...             prop=[0.4, 0.6]),
   ...     ],
   ...     matingScheme=sim.RandomMating(
   ...         ops=[
   ...             sim.MendelianGenoTransmitter(),
   ...             sim.MitochondrialGenoTransmitter()]),
   ...     preOps=[
   ...         sim.Stat(neutrality=range(5)),
   ...         sim.Stat(neutrality=range(5, 10), suffix='_X'),
   ...         sim.Stat(neutrality=range(10, 15), suffix='_Y'),
   ...         sim.Stat(neutrality=range(15, 20), suffix='_mt'),
   ...         sim.PyEval(r'"%.3f %.3f %.3f %.3f\n" % (Pi, Pi_X, Pi_Y, Pi_mt)'),
   ...     ],
   ...     gen = 2
   ... )
   1.921 1.900 1.973 1.914
   1.931 1.921 1.957 1.945
   2

   now exiting runScriptInteractively...

`Download statChromTypes.py <statChromTypes.py>`_


