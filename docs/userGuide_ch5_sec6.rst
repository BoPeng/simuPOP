Mutation
========

A mutator (a mutation operator) mutates alleles at certain loci from one allele
to another. Because alleles are simple non-nagative numbers that can be
intrepreted as nucleotides, codons, squences of nucleotides or even genetic
deletions, appropriate mutation models have to be chosen for different types of
loci. Please refer to Section :ref:`sec_Genotypic_structure
<sec_Genotypic_structure>` for a few examples.

A mutator will mutate alleles at all loci unless parameter ``loci`` is used to
specify a subset of loci. Different mutators have different concepts and forms
of mutation rates. If a mutator accepts only a single mutation rate (which can
be in the form of a list or a matrix), it uses parameter ``rate`` and applies
the same mutation rate to all loci. If a mutator accepts a list of mutation
rates (each of which is a single number), it uses parameter ``rates`` and
applies different mutation rates to different loci if multiple loci are
specified. Note that parameter ``rates`` also accepts single form inputs (e.g.
``rates=0.01``) in which case the same mutation rate will be applied to all
loci.


Mutation models specified by rate matrixes (:class:`MatrixMutator`)
-------------------------------------------------------------------

A mutation model can be defined as a **mutation rate matrix**
:math:`\left(p_{ij}\right)_{n\times n}` where :math:`p_{ij}` is the probability
that an allele :math:`i` mutates to :math:`j` per generation per locus. Although
mathematical formulation of :math:`p_{ij}` are sometimes unscaled, simuPOP
assumes :math:`\sum_{j=0}^{n-1}p_{ij}=1` for all :math:`i` and requires such
rate matrixes in the specification of a mutation model. :math:`p_{ii}` of such a
matrix are ignored because they are automatically calculated from
:math:`p_{ii}=1-\sum_{j\ne i}p_{ij}`.

A :class:`MatrixMutator` is defined to mutate between alleles 0, 1, ...,
:math:`n-1` according to a given rate matrix. Conceptually speaking, this
mutator goes through each mutable allele and mutates it to allele
:math:`0,1,..,n-1` according to probabilities :math:`p_{ij}`,
:math:`j=0,...,n-1`. Most alleles will be kept intact because mutations usually
happen at low probability (with :math:`p_{ii}` close to 1). For example, Example
:ref:`MatrixMutator <MatrixMutator>` simulates a locus with 3 alleles. Because
the rate at which allele 2 mutats to alleles 0 and 1 is higher than the rate
alleles 0 and 2 mutate to allele 2, the frequency of allele 2 decreases over
time.

.. _MatrixMutator:

**Example**: *General mutator specified by a mutation rate matrix*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[2000], loci=1)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.2, 0.3, 0.5])
   ...     ],
   ...     preOps=sim.MatrixMutator(rate = [
   ...             [0, 1e-5, 1e-5],
   ...             [1e-4, 0, 1e-4],
   ...             [1e-3, 1e-3, 0]
   ...         ]),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0, step=100),
   ...         sim.PyEval(r"', '.join(['%.3f' % alleleFreq[0][x] for x in range(3)]) + '\n'",
   ...             step=100),
   ...     ],
   ...     gen=1000
   ... )
   0.192, 0.302, 0.505
   0.241, 0.292, 0.467
   0.328, 0.273, 0.399
   0.270, 0.322, 0.408
   0.312, 0.412, 0.276
   0.330, 0.344, 0.327
   0.332, 0.424, 0.244
   0.426, 0.372, 0.201
   0.413, 0.384, 0.203
   0.395, 0.408, 0.198
   1000

   now exiting runScriptInteractively...

`Download MatrixMutator.py <MatrixMutator.py>`_

.. note::

   Alleles other than 0, 1, ..., $n-1$ will not be mutated because their mutation
   rates are undefined. A warning message will be displayed for this case when
   debugging code ``DBG_WARNING`` is turnned on.


k-allele mutation model (:class:`KAlleleMutator`)
-------------------------------------------------

A :math:`k`\ -allele model assumes :math:`k` alleles
:math:`\left(0,1,...,k-1\right)` at a locus and mutate between them using rate
matrix

.. math::

      p_{ij}=\left(\begin{array}{cccc}
      1-\mu & \frac{\mu}{k-1} & \cdots & \frac{\mu}{k-1}\\
      \frac{\mu}{k-1} & 1-\mu & \cdots & \frac{\mu}{k-1}\\
      \vdots & \vdots & \ddots & \vdots\\
      \frac{\mu}{k-1} & \frac{\mu}{k-1} & \cdots & 1-\mu
      \end{array}\right)

The only parameter :math:`\mu` is the mutation rate, which is the rate at which
an allele mutates to any other allele with equal probability.

This mutation model is a special case of the :class:`MatrixMutator` but a
specialized :class:`KAlleleMutator` is recommended because it provides better
performance, especially when :math:`k` is large. In addition, this operator
allows different mutation rates at different loci. When :math:`k` is not
specified, it is assumed to be the number of allowed alleles (e.g. 2 for binary
modules). Example :ref:`KAlleleMutator <KAlleleMutator>` desmonstrates the use
of this operator where parameters ``rate`` and ``loci`` are used to specify
different mutation rates for different loci. Because this operator treats all
alleles equally, all alleles will have the same allele frequency in the long
run.

.. _KAlleleMutator:

**Example**: *A k-allele mutation model*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[2000], loci=1*3)
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.KAlleleMutator(k=5, rates=[1e-2, 1e-3], loci=[0, 1]),
   ...         sim.Stat(alleleFreq=range(3), step=100),
   ...         sim.PyEval(r"', '.join(['%.3f' % alleleFreq[x][0] for x in range(3)]) + '\n'",
   ...             step=100),
   ...     ],
   ...     gen=500
   ... )
   0.991, 0.999, 1.000
   0.368, 0.918, 1.000
   0.300, 0.815, 1.000
   0.257, 0.639, 1.000
   0.209, 0.573, 1.000
   500

   now exiting runScriptInteractively...

`Download KAlleleMutator.py <KAlleleMutator.py>`_

.. note::

   If alleles $k$ and higher exist in the population, they will not be mutated
   because their mutation rates are undefined. A warning message will be displayed
   for this case when debugging code ``DBG_WARNING`` is turnned on.


Diallelic mutation models (:class:`SNPMutator`)
-----------------------------------------------

:class:`MatrixMutator` and :class:`KAlleleMutator` are general purpose mutators
in the sense that they do not assume a type for the mutated alleles. This and
the following sections describe mutation models for specific types of alleles.

If there are only two alleles at a locus, a diallelic mutation model should be
used. Because single nucleotide polymorphisms (SNPs) are the most widely
avaiable diallelic markers, a :class:`SNPMutator` is provided to mutate such
markers using a mutate rate matrix

.. math::

      R=\left(\begin{array}{cc}
      1-u & u\\
      v & 1-v
      \end{array}\right).

Despite of its name, this mutator can be used in many theoretical models
assuming :math:`\mbox{Pr}\left(A\rightarrow a\right)=u` and
:math:`\mbox{Pr}\left(a\rightarrow A\right)=v`. If :math:`v=0`, mutations will
be directional. Example :ref:`SNPMutator <SNPMutator>` applies such a
directional mutaton model to two loci, but with a purifying selection applied to
the first locus. Because of the selection pressure, the frequency of allele 1 at
the first locus does not increase indefinitely as allele 1 at the second locus.

.. _SNPMutator:

**Example**: *A diallelic directional mutation model*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[2000], loci=[1, 1], lociNames=['A', 'B'],
   ...     infoFields='fitness')
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     preOps=[
   ...         sim.SNPMutator(u=0.001),
   ...         sim.MaSelector(loci='A', fitness=[1, 0.99, 0.98]),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=['A', 'B'], step=100),
   ...         sim.PyEval(r"'%.3f\t%.3f\n' % (alleleFreq[0][1], alleleFreq[1][1])",
   ...             step=100),
   ...     ],
   ...     gen=500
   ... )
   0.001	0.001
   0.077	0.087
   0.099	0.192
   0.099	0.300
   0.085	0.400
   500

   now exiting runScriptInteractively...

`Download SNPMutator.py <SNPMutator.py>`_


Nucleotide mutation models (:class:`AcgtMutator`)
-------------------------------------------------

Mutations in these models assume alleles 0, 1, 2, 3 as nucleotides A, C, G, and
T. The operator is named :class:`AcgtMutator` to remind you the alphabetic order
of these nucleotides. This mutation model is specified by a rate matrix

.. math::

      \begin{array}{ccccc}
       & A & C & G & T\\
      A & - & x_{1} & x_{2} & x_{3}\\
      C & x_{4} & - & x_{5} & x_{6}\\
      G & x_{7} & x_{8} & - & x_{9}\\
      T & x_{10} & x_{11} & x_{12} & -
      \end{array}

which is determined by 12 parameters. However, several simpler models with fewer
parameters can be used. In addition to parameters shared by all mutation
operators, a nucleotide mutator is specified by a parameter list and a model
name. For example:

::

   AcgtMutator(rate=[1e-5, 0.5], model='K80')

specifies a nucleotide mutator using Kimura's 2-parameter model with
:math:`\mu=10^{-5}` and :math:`\kappa=0.5`. Because multiple parameters could be
involved for a particular mutation model, **the definition of a mutation rate
and other paramters are model dependent and may varying with different
mathematical representation of the models**.

The names and acceptable parameters of acceptable models are listed below:

#. Jukes and Cantor 1969 model: ``model='JC69'``, rate=[:math:`\mu` ]

   The Jukes and Cantor model is similar to a :math:`4`\ -allele model but its
   definition of :math:`\mu` is different. More specifically, when a mutation event
   happens at rate :math:`\mu`, an allele will have equal probability to mutate to
   any of the 4 allelic states.

   .. math::

      R=\left(\begin{array}{cccc}
      - & \frac{\mu}{4} & \frac{\mu}{4} & \frac{\mu}{4}\\
      \frac{\mu}{4} & - & \frac{\mu}{4} & \frac{\mu}{4}\\
      \frac{\mu}{4} & \frac{\mu}{4} & - & \frac{\mu}{4}\\
      \frac{\mu}{4} & \frac{\mu}{4} & \frac{\mu}{4} & -
      \end{array}\right)

#. Kimura's 2-parameter 1980 model: ``model='K80'``, rate=[:math:`\mu`,
   :math:`\kappa`]

   Kimura 's model distinguishes transitions (:math:`A\longleftrightarrow G`, and
   :math:`C\leftrightarrow T` namely :math:`0\longleftrightarrow2` and
   :math:`1\longleftrightarrow3` with probability :math:`\frac{\mu}{4}\kappa`) and
   transversions (others) with probability :math:`\frac{\mu}{4}`. It would be a
   Jukes and Cantor model if :math:`\kappa=1`.

   .. math::

      R=\left(\begin{array}{cccc}
      - & \frac{\mu}{4} & \frac{\mu}{4}\kappa & \frac{\mu}{4}\\
      \frac{\mu}{4} & - & \frac{\mu}{4} & \frac{\mu}{4}\kappa\\
      \frac{\mu}{4}\kappa & \frac{\mu}{4} & - & \frac{\mu}{4}\\
      \frac{\mu}{4} & \frac{\mu}{4}\kappa & \frac{\mu}{4} & -
      \end{array}\right)

#. Felsenstein 1981 model: ``model='F81'``, rate=[:math:`\mu`, :math:`\pi_{A}`,
   :math:`\pi_{C}`, :math:`\pi_{G}`].

   This model assumes different base frequencies but the same probabilities for
   transitions and transversions. :math:`\pi_{T}` is calculated from
   :math:`\pi_{A}`, :math:`\pi_{C}` and :math:`\pi_{G}`.

   .. math::

      R=\left(\begin{array}{cccc}
      - & \mu\pi_{C} & \mu\pi_{G} & \mu\pi_{T}\\
      \mu\pi_{A} & - & \mu\pi_{G} & \mu\pi_{T}\\
      \mu\pi_{A} & \mu\pi_{C} & - & \mu\pi_{T}\\
      \mu\pi_{A} & \mu\pi_{C} & \mu\pi_{G} & -
      \end{array}\right)

#. Hasegawa, Kishino and Yano 1985 model: ``model='HKY85'``, rate=[:math:`\mu`,
   :math:`\kappa`, :math:`\pi_{A}`, :math:`\pi_{C}`, :math:`\pi_{G}`]

   This model replaces 1/4 frequency used in the Kimura's 2-parameter model with
   nucleotide-specific frequencies.

   .. math::

      R=\left(\begin{array}{cccc}
      - & \mu\pi_{C} & \mu\kappa\pi_{G} & \mu\pi_{T}\\
      \mu\pi_{A} & - & \mu\pi_{G} & \mu\kappa\pi_{T}\\
      \mu\kappa\pi_{A} & \mu\pi_{C} & - & \mu\pi_{T}\\
      \mu\pi_{A} & \mu\kappa\pi_{C} & \mu\pi_{G} & -
      \end{array}\right)

#. Tamura 1992 model: ``model='T92'``, rate=[:math:`\mu`, :math:`\pi_{GC}`]

   This model is a HKY85 model with :math:`\pi_{G}=\pi_{C}=\pi_{GC}/2` and
   :math:`\pi_{A}=\pi_{T}=\pi_{AT}/2=\left(1-\pi_{GC}\right)/2`,

   .. math::

      R=\left(\begin{array}{cccc}
      - & \frac{1}{2}\mu\pi_{GC} & \frac{1}{2}\mu\nu\pi_{GC} & \frac{1}{2}\mu\pi_{AT}\\
      \frac{1}{2}\mu\pi_{AT} & - & \frac{1}{2}\mu\pi_{GC} & \frac{1}{2}\mu\nu\pi_{AT}\\
      \frac{1}{2}\mu\nu\pi_{AT} & \frac{1}{2}\mu\pi_{GC} & - & \frac{1}{2}\mu\pi_{AT}\\
      \frac{1}{2}\mu\pi_{AT} & \frac{1}{2}\mu\nu\pi_{GC} & \frac{1}{2}\mu\pi_{GC} & -
      \end{array}\right)

#. Tamura and Nei 1993 model: ``model='TN93'``, rate=[:math:`\mu`,
   :math:`\kappa_{1}`, :math:`\kappa_{2}`, :math:`\pi_{A}`, :math:`\pi_{C}`,
   :math:`\pi_{G}`]

   This model extends the HKY1985 model by distinguishing
   :math:`A\longleftrightarrow G` transitions (namely
   :math:`0\longleftrightarrow2`) and :math:`C\leftrightarrow T` transitions
   (:math:`1\longleftrightarrow3`) with different :math:`\kappa`.

   .. math::

      R=\left(\begin{array}{cccc}
      - & \mu\pi_{C} & \mu\kappa_{1}\pi_{G} & \mu\pi_{T}\\
      \mu\pi_{A} & - & \mu\pi_{G} & \mu\kappa_{2}\pi_{T}\\
      \mu\kappa_{1}\pi_{A} & \mu\pi_{C} & - & \mu\pi_{T}\\
      \mu\pi_{A} & \mu\kappa_{2}\pi_{C} & \mu\pi_{G} & -
      \end{array}\right)

#. Generalized time reversible model: ``model='GTR'``, rate=[:math:`x_{1}`,
   :math:`x_{2}`, :math:`x_{3}`, :math:`x_{4}`, :math:`x_{5}`, :math:`x_{6}`,
   :math:`\pi_{A}`, :math:`\pi_{C}`, :math:`\pi_{G}`]

   The generalized time reviersible model is the most general neutral,
   indepdendent, finite-sites, time-reversible model possible. It is specified by
   six parameters and base frequencies. Its rate matrix is defined as

   .. math::

      R=\left(\begin{array}{cccc}
      - & \frac{\pi_{A}x_{1}}{\pi_{C}} & \frac{\pi_{A}x_{2}}{\pi_{G}} & \frac{\pi_{A}x_{3}}{\pi_{T}}\\
      x_{1} & - & \frac{\pi_{C}x_{4}}{\pi_{G}} & \frac{\pi_{C}x_{5}}{\pi_{T}}\\
      x_{2} & x_{4} & - & \frac{\pi_{G}x_{6}}{\pi_{T}}\\
      x_{3} & x_{5} & x_{6} & -
      \end{array}\right)

#. General model: ``model='general'`` (default), rate=[:math:`x_{1}`,
   :math:`x_{2}`, :math:`x_{3}`, :math:`x_{4}`, :math:`x_{5}`, :math:`x_{6}`,
   :math:`x_{7}`, :math:`x_{8}`, :math:`x_{9}`, :math:`x_{10}`, :math:`x_{11}`,
   :math:`x_{12}`].

   This is the most general model with 12 parameters:

   .. math::

      R=\left(\begin{array}{cccc}
      - & x_{1} & x_{2} & x_{3}\\
      x_{4} & - & x_{5} & x_{6}\\
      x_{7} & x_{8} & - & x_{9}\\
      x_{10} & x_{11} & x_{12} & -
      \end{array}\right)

   It is not surprising that all other models are implemented as special cases of
   this model.

Example :ref:`AcgtMutator <AcgtMutator>` applies a Kimmura's 2-parameter
mutation model to a population with a single nucleotide marker.

.. _AcgtMutator:

**Example**: *A Kimura's 2 parameter mutation model*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[2000], loci=1,
   ...     alleleNames=['A', 'C', 'G', 'T'])
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[.1, .1, .1, .7])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     preOps=[
   ...         sim.AcgtMutator(rate=[1e-4, 0.5], model='K80'),
   ...         sim.Stat(alleleFreq=0, step=100),
   ...         sim.PyEval(r"', '.join(['%.3f' % alleleFreq[0][x] for x in range(4)]) + '\n'",
   ...             step=100),
   ...     ],
   ...     gen=500
   ... )
   0.093, 0.101, 0.094, 0.712
   0.142, 0.073, 0.084, 0.701
   0.135, 0.160, 0.083, 0.623
   0.230, 0.128, 0.013, 0.628
   0.293, 0.189, 0.008, 0.510
   500

   now exiting runScriptInteractively...

`Download AcgtMutator.py <AcgtMutator.py>`_


Mutation model for microsatellite markers (:class:`StepwiseMutator`)
--------------------------------------------------------------------

The **stepwise mutation model** (SMM) was proposed by Ohta1973 to model the
mutation of Variable Number Tandem Repeat (VNTR), which consists of tandem
repeat of sequences. VNTR markers consisting of short sequences (e.g. 5 basepair
or less) are also called microsatellite markers. A mutation event of a VNTR
marker either increase of decrease the number of repeats, as a result of
slipped-strand mispairing or unequal sister chromatid exchange and genetic
recombination.

A :class:`StepwiseMutator` assumes that alleles at a locus are the number of
tandem repeats and mutates them by increasing or decreasing the number of
repeats during a mutation event. By adjusting parameters ``incProb``,
``maxAllele`` and ``mutStep``, this operator can be used to simulate the
standard neutral stepwise mutation model and a number of **generalized stepwise
mutation models**. For example, Example :ref:`StepwiseMutator <StepwiseMutator>`
uses two :class:`StepwiseMutator` to mutate two microsatellite markers, using a
standard and a generalized model where a geometric distribution is used to
determine the number of steps.

.. _StepwiseMutator:

**Example**: *A standard and a generalized stepwise mutation model*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=1000, loci=[1, 1])
   >>> pop.evolve(
   ...     # all start from allele 50
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq= [0]*50 + [1])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     preOps=[
   ...         sim.StepwiseMutator(rates=1e-3, loci=0),
   ...         sim.StepwiseMutator(rates=1e-3, incProb=0.6, loci=1,
   ...             mutStep=(sim.GEOMETRIC_DISTRIBUTION, 0.2)),
   ...     ],
   ...     gen=100
   ... )
   100
   >>> # count the average number tandem repeats at both loci
   >>> cnt0 = cnt1 = 0
   >>> for ind in pop.individuals():
   ...     cnt0 += ind.allele(0, 0) + ind.allele(0, 1)
   ...     cnt1 += ind.allele(1, 0) + ind.allele(1, 1)
   ... 
   >>> print('Average number of repeats at two loci are %.2f and %.2f.' % \
   ...     (cnt0/2000., cnt1/2000.))
   Average number of repeats at two loci are 50.03 and 49.70.

   now exiting runScriptInteractively...

`Download StepwiseMutator.py <StepwiseMutator.py>`_


Simulating arbitrary mutation models using a hybrid mutator (:class:`PyMutator`)\*
----------------------------------------------------------------------------------

A hybrid mutator :class:`PyMutator` mutates random alleles at selected loci
(parameter ``loci``), replicates (parameter ``loci``), subpopulations (parameter
``subPop``) with specified mutation rate (parameter ``rate``). Instead of
mutating the alleles by itself, it passes the alleles to a user-defined function
and use it return values as the mutated alleles. Arbitrary mutation models could
be implemented using this operator.

Example :ref:`PyMutator <PyMutator>` applies a simple mutation model where an
allele is increased by a random number between 1 and 5 when it is mutated. Two
different mutation rates are used for two different loci so average alleles at
these two loci are different.

.. _PyMutator:

**Example**: *A hybrid mutation model*

::

   >>> import simuPOP as sim
   >>> import random
   >>> def incAllele(allele):
   ...     return allele + random.randint(1, 5)
   ... 
   >>> pop = sim.Population(size=1000, loci=[20])
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=sim.PyMutator(func=incAllele, rates=[1e-4, 1e-3],
   ...             loci=[2, 10]),
   ...     gen = 1000
   ... )
   1000
   >>> # count the average number tandem repeats at both loci
   >>> def avgAllele(pop, loc):
   ...     ret = 0
   ...     for ind in pop.individuals():
   ...         ret += ind.allele(loc, 0) + ind.allele(loc, 1)
   ...     return ret / (pop.popSize() * 2.)
   ... 
   >>> print('Average number of repeats at two loci are %.2f and %.2f.' % \
   ...     (avgAllele(pop, 2), avgAllele(pop, 10)))
   Average number of repeats at two loci are 0.01 and 2.19.

   now exiting runScriptInteractively...

`Download PyMutator.py <PyMutator.py>`_


Mixed mutation models (:class:`MixedMutator`) \*\*
--------------------------------------------------

Mixed mutation models are sometimes used to model real data. For example, a
:math:`k`\ -allele model can be used to explain extremely large or small number
of tandem repeats at a microsatellite marker which are hard to justify using a
standard stepwise mutation model. A mixed mutation model would apply two or more
mutation models at pre-specified probabilities.

A :class:`MixedMutator` is constructed by a list of mutators and their
respective probabilities. It accepts regular mutator parameters such as
``rates``, ``loci``, ``subPops``, ``mapIn and mapOut`` and mutates aleles at
specified rate. When a mutation event happens, it calls one of the mutators to
mutate the allele. For example, Example :ref:`MixedMutator <MixedMutator>`
applies a mixture of :math:`k`\ -allele model and stepwise model to mutate a
micosatellite model.

.. _MixedMutator:

**Example**: *A mixed k-allele and stepwise mutation model*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(5000, loci=[1, 1])
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[50, 50])
   ...     ],
   ...     preOps=[
   ...         # the first locus uses a pure stepwise mutation model
   ...         sim.StepwiseMutator(rates=0.001, loci=0),
   ...         # the second locus uses a mixed model
   ...         sim.MixedMutator(rates=0.001, loci=1, mutators=[        
   ...             sim.KAlleleMutator(rates=1, k=100),
   ...             sim.StepwiseMutator(rates=1)
   ...         ], prob=[0.1, 0.9])],
   ...     matingScheme=sim.RandomMating(),
   ...     gen = 20
   ... )
   20
   >>> # what alleles are there?
   >>> geno0 = []
   >>> geno1 = []
   >>> for ind in pop.individuals():
   ...     geno0.extend([ind.allele(0, 0), ind.allele(0, 1)])
   ...     geno1.extend([ind.allele(1, 0), ind.allele(1, 1)])
   ... 
   >>> print('Locus 0 has alleles', ', '.join([str(x) for x in set(geno0)]))
   Locus 0 has alleles 49, 50, 51
   >>> print('Locus 1 has alleles', ', '.join([str(x) for x in set(geno1)]))
   Locus 1 has alleles 67, 49, 50, 51, 88

   now exiting runScriptInteractively...

`Download MixedMutator.py <MixedMutator.py>`_

When a mutation event happens, mutators in Example :ref:`MixedMutator
<MixedMutator>` mutate the allele with probability (mutation rate) 1. If
different mutation rates are specified, the overall mutation rates would be the
product of mutation rate of :class:`MixedMutator` and the passed mutators.
However, it is extremely important to understand that although
:class:`MixedMutator`\ (``rates=mu``) with :class:`StepwiseMutator`\
(``rates=1``) and :class:`MixedMutator`\ ``(rates=1)``\ with
:class:`StepwiseMutator`\ (``rates=mu``) mutate alleles at the same mutation
rate, the former is much more efficient because it triggers far less mutation
events.


Context-dependent mutation models (:class:`ContextMutator`)\*\*
---------------------------------------------------------------

All mutation models we have seen till now are context independent. That is to
say, how an allele is mutated depends only on the allele itself. However, it is
understood that DNA and amino acid substitution rates are highly sequence
context-dependent, e.g., C :math:`\rightarrow` T substitutions in vertebrates
may occur much more frequently at CpG sites. To simulate such models, a mutator
must consider the context of a mutated allele, e.g. certain number of alleles to
the left and right of this allele, and mutate the allele accordingly.

A :class:`ContextMutator` can be used to mutate an allele depending on its
surrounding loci. This mutator is constructed by a list of mutators and their
respective contexts. It accepts regular mutator parameters such as ``rates``,
``loci``, ``subPops``, ``mapIn and mapOut`` and mutates aleles at specified
rate. When a mutation event happens, it checks the context of the mutaed allele
and choose a corresponding mutator to mutate the allele. An additional mutator
can be specified to mutate alleles with unknown context. Example
:ref:`ContextMutator <ContextMutator>` applies two :class:`SNPMutator` at
different rates under different contexts.

.. _ContextMutator:

**Example**: *A context-dependent mutation model*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(5000, loci=[3, 3])
   >>> pop.evolve(
   ...     # initialize locus by 0, 0, 0, 1, 0, 1
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[1, 1], loci=[3, 5])
   ...     ],
   ...     preOps=[
   ...         sim.ContextMutator(mutators=[
   ...             sim.SNPMutator(u=0.1),
   ...             sim.SNPMutator(u=1),
   ...             ],
   ...             contexts=[(0, 0), (1, 1)],
   ...             loci=[1, 4],
   ...             rates=0.01
   ...         ),
   ...         sim.Stat(alleleFreq=[1, 4], step=5),
   ...         sim.PyEval(r"'Gen: %2d freq1: %.3f, freq2: %.3f\n'" + 
   ...             " % (gen, alleleFreq[1][1], alleleFreq[4][1])", step=5)
   ...     ], 
   ...     matingScheme=sim.RandomMating(),
   ...     gen = 20
   ... )
   Gen:  0 freq1: 0.001, freq2: 0.010
   Gen:  5 freq1: 0.005, freq2: 0.059
   Gen: 10 freq1: 0.007, freq2: 0.108
   Gen: 15 freq1: 0.015, freq2: 0.142
   20

   now exiting runScriptInteractively...

`Download ContextMutator.py <ContextMutator.py>`_

Note that although   ::

   ContextMutator(mutators=[
       SNPMutator(u=0.1),
       SNPMutator(u=1)],
       contexts=[(0, 0), (1, 1)],
       rates=0.01
   )    

and   ::

   ContextMutator(mutators=[
       SNPMutator(u=0.001),
       SNPMutator(u=0.01)],
       contexts=[(0, 0), (1, 1)],
       rates=1
   )    

both apply two :class:`SNPMutator` at mutation rates ``0.001`` and ``0.01``, the
former is more efficient because it triggers less mutation events.

Context-dependent mutator can also be implemented by a :class:`PyMutator`. When
a non-zero parameter ``context`` is specified, this mutator will collect
``context`` number of alleles to the left and right of a mutated allele and pass
them as a second parameter of the user-provided mutation function. Example
:ref:`pyContextMutator <pyContextMutator>` applies the same mutation model as
Example :ref:`ContextMutator <ContextMutator>` using a :class:`PyMutator`.

.. _pyContextMutator:

**Example**: *A hybrid context-dependent mutation model*

::

   >>> import simuPOP as sim
   >>> import random
   >>> pop = sim.Population(5000, loci=[3, 3])
   >>> def contextMut(allele, context):
   ...     if context == [0, 0]:
   ...         if allele == 0 and random.random() < 0.1:
   ...             return 1
   ...     elif context == [1, 1]:
   ...         if allele == 0:
   ...             return 1
   ...     # do not mutate
   ...     return allele
   ... 
   >>> pop.evolve(
   ...     # initialize locus by 0, 0, 0, 1, 0, 1
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[1, 1], loci=[3, 5])
   ...     ],
   ...     preOps=[
   ...         sim.PyMutator(func=contextMut, context=1,
   ...             loci=[1, 4],  rates=0.01
   ...         ),
   ...         #sim.SNPMutator(u=0.01, v= 0.01, loci=[1, 4]),
   ...         sim.Stat(alleleFreq=[1, 4], step=5),
   ...         sim.PyEval(r"'Gen: %2d freq1: %.3f, freq2: %.3f\n'" + 
   ...             " % (gen, alleleFreq[1][1], alleleFreq[4][1])", step=5)
   ...     ], 
   ...     matingScheme=sim.RandomMating(),
   ...     gen = 20
   ... )
   Gen:  0 freq1: 0.000, freq2: 0.000
   Gen:  5 freq1: 0.000, freq2: 0.000
   Gen: 10 freq1: 0.000, freq2: 0.000
   Gen: 15 freq1: 0.000, freq2: 0.000
   20

   now exiting runScriptInteractively...

`Download pyContextMutator.py <pyContextMutator.py>`_


Manually-introduced mutations (:class:`PointMutator`)
-----------------------------------------------------

Operator :class:`PointMutator` is different from all other mutators in that it
mutates specified alleles of specified individuals. It is usually used to
manually introduce one or more mutants to a population. Although it is not a
recommended method to introduce a disease predisposing allele, the following
example (Example :ref:`PointMutator <PointMutator>`) demonstrates an
evolutionary process where mutants are repeatedly introduced and raised by
positive selection until it reaches an appreciable allele frequency. This
example uses two :class:`IfElse` operators. The first one introduces a mutant
when there is no mutant in the population, and the second one terminate the
evolution when the frequency of the mutant reaches 0.05.

.. _PointMutator:

**Example**: *Use a point mutator to introduce a disease predisposing allele*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, loci=1, infoFields='fitness')
   >>> pop.evolve(
   ...     initOps=sim.PyOutput('Introducing alleles at generation'),
   ...     preOps=sim.MaSelector(loci=0, wildtype=0, fitness=[1, 1.05, 1.1]),
   ...     matingScheme=sim.RandomSelection(),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=0),
   ...         sim.IfElse('alleleNum[0][1] == 0', ifOps=[
   ...             sim.PyEval(r"' %d' % gen"),
   ...             sim.PointMutator(inds=0, loci=0, allele=1),
   ...         ]),
   ...         sim.IfElse('alleleFreq[0][1] > 0.05', ifOps=[
   ...             sim.PyEval(r"'.\nTerminate at generation %d at allele freq %.3f.\n'" +
   ...                 " % (gen, alleleFreq[0][1])"),
   ...             sim.TerminateIf('True'),
   ...         ])
   ...     ],
   ... )
   Introducing alleles at generation 0 1 2 16 17 18 22 30 32 33 34 41 81 82 83.
   Terminate at generation 111 at allele freq 0.051.
   112

   now exiting runScriptInteractively...

`Download PointMutator.py <PointMutator.py>`_


Apply mutation to (virtual) subpopulations \*
---------------------------------------------

A mutator is usually applied to all individuals in a population. However, you
can restrict its use to specified subpopulations and/or virtual subpopulations
using parameter ``subPop``. For example, you can use ``subPop=[0, 2]`` to apply
the mutator only to individuals in subpopulations 0 and 2.

Virtual subpopulations can also be specified in this parameter. For example, you
can apply different mutation models to male and female individuals, to
unaffected or affected individuals, to patients at different stages of a cancer.
Example :ref:`mutatorVSP <mutatorVSP>` demonstrate a mutation model where
individuals with more tandem repeats at a disease predisposing locus are more
likely to develop a disease (e.g. fragile-X). Affected individuals are then
subject to a non-neutral mutation model at an accerlerated mutation rate.

.. _mutatorVSP:

**Example**: *Applying mutation to virtual subpopulations.*

::

   >>> import simuPOP as sim
   >>> def fragileX(geno):
   ...     '''A disease model where an individual has increased risk of 
   ...     affected if the number of tandem repeats exceed 75.
   ...     '''
   ...     # Alleles A1, A2.
   ...     maxRep = max(geno)
   ...     if maxRep < 50:
   ...         return 0
   ...     else:
   ...         # individuals with allele >= 70 will surely be affected
   ...         return min(1, (maxRep - 50)*0.05)
   ... 
   >>> def avgAllele(pop):
   ...     'Get average allele by affection sim.status.'
   ...     sim.stat(pop, alleleFreq=(0,1), subPops=[(0,0), (0,1)],
   ...         numOfAffected=True, vars=['alleleNum', 'alleleNum_sp'])
   ...     avg = []
   ...     for alleleNum in [\
   ...             pop.dvars((0,0)).alleleNum[0],  # first locus, unaffected
   ...             pop.dvars((0,1)).alleleNum[0],  # first locus, affected
   ...             pop.dvars().alleleNum[1],       # second locus, overall
   ...         ]:
   ...         alleleSum = numAllele = 0
   ...         for idx,cnt in enumerate(alleleNum):
   ...             alleleSum += idx * cnt
   ...             numAllele += cnt
   ...         if numAllele == 0:
   ...             avg.append(0)
   ...         else:
   ...             avg.append(alleleSum * 1.0 /numAllele)
   ...     # unaffected, affected, loc2
   ...     pop.dvars().avgAllele = avg
   ...     return True
   ... 
   >>> pop = sim.Population(10000, loci=[1, 1])
   >>> pop.setVirtualSplitter(sim.AffectionSplitter())
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[50, 50])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         # determine affection sim.status for each offspring (duringMating)
   ...         sim.PyPenetrance(func=fragileX, loci=0),
   ...         # unaffected offspring, mutation rate is high to save some time
   ...         sim.StepwiseMutator(rates=1e-3, loci=1),
   ...         # unaffected offspring, mutation rate is high to save some time
   ...         sim.StepwiseMutator(rates=1e-3, loci=0, subPops=[(0, 0)]),
   ...         # affected offspring have high probability of mutating upward
   ...         sim.StepwiseMutator(rates=1e-2, loci=0, subPops=[(0, 1)],
   ...            incProb=0.7, mutStep=3),
   ...         # number of affected
   ...         sim.PyOperator(func=avgAllele, step=20),
   ...         sim.PyEval(r"'Gen: %3d #Aff: %d AvgRepeat: %.2f (unaff), %.2f (aff), %.2f (unrelated)\n'"
   ...             + " % (gen, numOfAffected, avgAllele[0], avgAllele[1], avgAllele[2])",
   ...             step=20),
   ...     ],
   ...     gen = 101
   ... )
   Gen:   0 #Aff: 0 AvgRepeat: 1.01 (unaff), 0.00 (aff), 1.01 (unrelated)
   Gen:  20 #Aff: 6 AvgRepeat: 1.53 (unaff), 0.50 (aff), 1.52 (unrelated)
   Gen:  40 #Aff: 20 AvgRepeat: 2.56 (unaff), 2.04 (aff), 1.53 (unrelated)
   Gen:  60 #Aff: 46 AvgRepeat: 2.56 (unaff), 2.04 (aff), 2.04 (unrelated)
   Gen:  80 #Aff: 55 AvgRepeat: 3.08 (unaff), 1.53 (aff), 2.04 (unrelated)
   Gen: 100 #Aff: 48 AvgRepeat: 2.04 (unaff), 1.52 (aff), 2.04 (unrelated)
   101

   now exiting runScriptInteractively...

`Download mutatorVSP.py <mutatorVSP.py>`_

At the beginning of a simulation, all individuals have 50 copies of a tandem
repeat and the mutation follows a standard neutral stepwise mutation model.
individuals with more than 50 repeats will have an increasing probability to
develop a disease (:math:`\mbox{Pr}\left(\mbox{affected}\mid
n\right)=\left(n-50\right)\*0.05`) for :math:`50\le n\le70`). The averge repeat
number therefore increases for affected individuals. In contrast, the mean
number of repeats at locus 1 on a separate chromosome oscillate around 50.


Allele mapping \*\*
-------------------

If alleles in your simulation do not follow the convention of a mutation model,
you may want to use the ``pop.recodeAlleles()`` function to recode your alleles
so that appropriate mutation models could be applied. If this is not possible,
you can use a general mutation model with your own mutation matrix, or an
advanced feature called **allele mapping**.

Allele mapping is done through two parameters *mapIn* and *mapOut*, which map
alleles in your population to and from alleles assumed in a mutation model. For
example, an :class:`AcgtMutator` mutator assumes alleles ``A``, ``C``, ``G`` and
``T`` for alleles 0, 1, 2, and 3 respectively. If for any reason the alleles in
your application does not follow this order, you will need to map these alleles
to the alleles assumed in the mutator. For example, if you assumes ``C``, ``G``,
``A``, ``T`` for alleles 0, 1, 2, and 3 respectively, you can use parameters

::

   mapIn=[1, 2, 0, 3], mapOut=[2, 0, 1, 3]

to map your alleles (``C(0)->C(1)``, ``G(1)->G(2)``, ``A(2)->A(0)``,
``T(3)->T(3)``) to alleles :class:`AcgtMutator` assumes, and then map mutated
alleles (``A(0)->A(2)``, ``C(1)->C(0)``, ``G(2)->G(1)``, ``T(3)->T(3)``) back.
Example :ref:`alleleMapping <alleleMapping>` gives another example where alleles
4, 5, 6 and 7 are mutated using a 4-allele model.

.. _alleleMapping:

**Example**: *Allele mapping for mutation operators*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[2000], loci=1)
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0]*4 + [0.1, 0.2, 0.3, 0.4])
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     postOps=[
   ...         sim.KAlleleMutator(k=4, rates=1e-4, mapIn=[0]*4 + list(range(4)),
   ...             mapOut=[4, 5, 6, 7]),
   ...         sim.Stat(alleleFreq=0, step=100),
   ...         sim.PyEval(r"', '.join(['%.2f' % alleleFreq[0][x] for x in range(8)]) + '\n'",
   ...             step=100),
   ...     ],
   ...     gen=500
   ... )
   0.00, 0.00, 0.00, 0.00, 0.09, 0.20, 0.30, 0.41
   0.00, 0.00, 0.00, 0.00, 0.13, 0.20, 0.40, 0.26
   0.00, 0.00, 0.00, 0.00, 0.17, 0.20, 0.31, 0.31
   0.00, 0.00, 0.00, 0.00, 0.19, 0.18, 0.26, 0.37
   0.00, 0.00, 0.00, 0.00, 0.18, 0.24, 0.23, 0.34
   500

   now exiting runScriptInteractively...

`Download alleleMapping.py <alleleMapping.py>`_

These two parameters also accept Python functions which should return
corresponding mapped-in or out allele for a given allele. These two functions
can be used to explore very fancy mutation models. For example, you can
categorize a large number of alleles into alleles assumed in a mutation model,
and emit random alleles from a mutated allele.


Mutation rate and transition matrix of a ``MatrixMutator``\*\*
--------------------------------------------------------------

A :class:`MatrixMutator` is specified by a mutation rate matrix. Although
mutation rates of this mutator is typically allele-dependent, the
:class:`MatrixMutator` is implemented as a two-step process where mutation
events are triggered independent to allelic states. This section describes these
two steps which can be useful if you need to use a ``maxtrixMutator`` in a
:class:`MixedMutator` or :class:`ContextMutator`, and would like to factor out
an allele-independent mutation rate to the wrapper mutator.

Because alleles usually have different probabilities of mutating to other
alleles, **a mutation process is usually allele dependent**. Given a mutation
model :math:`\left(p_{ij}\right)`, it is obviously inefficient to go through all
mutable alleles and determine whether or not to mutate it using :math:`p_{ij},`
:math:`j=0,...,1-n`. simuPOP uses a two step procedure to mutate a large number
of alleles. More specifically, for each mutation model, we determine
:math:`\mu=\max_{i=0}^{n-1}\left(1-p_{ii}\right)` as the overall mutation rate,
and then

#. For each allele, trigger a mutation event with probability :math:`\mu`.
   Because :math:`\mu` is usually very small and is the same for all alleles, this
   step can be implemented efficiently.

#. When a mutation event happens, mutation allele :math:`i` to allele :math:`j`
   with probability

   .. math::

      \mbox{Pr}\left(i\rightarrow j\right)=\begin{cases}
      1-\frac{1}{\mu}\left(1-p_{ii}\right) & \mbox{if }i=j\\
      \frac{p_{ij}}{\mu} & \mbox{if }i\ne j
      \end{cases}

Because steps 1 and 2 are independent, it is easy to verify that

.. math::

      p_{ij}=\mu\mbox{Pr}\left(i\rightarrow j\right)

if :math:`i\ne j` and

.. math::

      p_{ii}=\left(1-\mu\right)+\mu\mbox{Pr}\left(i\rightarrow i\right)

where the first and second items are probabilities of no-mutation at steps 1 and
2. :math:`\mu` was chosen as the smallest :math:`\mu` that makes
:math:`0\leq\mbox{Pr}\left(i\rightarrow i\right)\leq1` for all :math:`i`.

For example, for a :math:`k`\ -allele model with

.. math::

      p_{ij}=\left(\begin{array}{cccc}
      1-\mu & \frac{\mu}{k-1} & \cdots & \frac{\mu}{k-1}\\
      \frac{\mu}{k-1} & 1-\mu & \cdots & \frac{\mu}{k-1}\\
      \vdots & \vdots & \ddots & \vdots\\
      \frac{\mu}{k-1} & \frac{\mu}{k-1} & \cdots & 1-\mu
      \end{array}\right)

:math:`\mu` is directly :math:`\mu` for the first step and

.. math::

      \mbox{Pr}\left(i\rightarrow j\right)=\begin{cases}
      0 & \mbox{if }i=j\\
      \frac{1}{k-1} & \mbox{if }i\ne j
      \end{cases}

for the second step. Therefore, mutation rate :math:`\mu` in a :math:`k`\
-allele model could be interpreted as the probability of mutation, and a
mutation event would mutate an allele to any other allele with equal
probability.

For a classical mutation model with :math:`P\left(A\rightarrow a\right)=u` and
:math:`P\left(a\rightarrow A\right)=v`,

.. math::

      p_{ij}=\left(\begin{array}{cc}
      1-u & u\\
      v & 1-v
      \end{array}\right)

if :math:`u=0.001` and :math:`v=0.0005`, :math:`\mu=\max\left(u,v\right)=0.001`,

.. math::

      \mbox{Pr}\left(i\rightarrow j\right)=\left(\begin{array}{cc}
      0 & 1\\
      \frac{v}{u}=0.5 & 1-\frac{v}{u}=0.5
      \end{array}\right)

That is to say, we would mutate at a mutation rate :math:`u=0.001`, mutate
allele :math:`A` to :math:`a` with probability 1 and mutate allele :math:`a` to
:math:`A` with probability 0.5.


Infinite-sites model and other simulation techniques \*\*
---------------------------------------------------------

Infinite-sites and infinite-alleles models have some similarities. If you assume
that mutation is the only force to create new mutants, you can treat a long
chromosomal region as a locus and use the infinite-alleles model, actually a
:math:`k`\ -allele model with large :math:`k`, to mimic the infinite-site model.
This assumption is certainly wrong with the infinite-site model when
recombination is involved, because recombination creates new haplotypes
(alleles) under the infinite-site model. However, for short regions where
recombination can be ignored, an :math:`k`\ -allele model can be an easy and
fast way to mimic an infinite-site model. That statement basically says that you
have a choice between two models if you would like to simulate the evolution of
this gene, namely considering the gene as a locus and simulating variants as
alleles, or considering the gene as a sequence and simulating haplotypes as
alleles.

For example, the CFTR gene (for cystic fibrosis) can have many alleles (thinking
in terms of infinite-allele model) which are nucleotide mutations on tens of
locations (infinite-site model). In order to simulate the evolution of this
gene, you have a choice between two models, namely considering the gene as a
locus and simulating variants as alleles, or considering the gene as a sequence
and simulating haplotypes as alleles. Because there is supposed to be only one
mutant at each site, you can assign a unique *location* for each allele of an
infinite-allele model and convert multi-allelic datasets simulated by an
infinite-allele model to sequences of diallelic markers. Note that mutation
rates are interpreted differently for these two models.

If specific location of such a mutation is needed, it is possible to record the
location of mutations during an evolution and minic an infinite-sites model. For
example, alleles in Example :ref:`infiniteSites <infiniteSites>` are used to
store location of a mutation event. When a mutation event happens, the location
of the new allele (rather the allele itself) is recorded on the chromosome
(actually list of mutation events) of an individual. The transmission of
chromosomes proceed normally and effectively transmit mutants from parents to
offspring. At the end of the simulation, each individual accumulates a number of
mutation events and they are essentially alleles at their respective locations.

.. _infiniteSites:

**Example**: *Mimicking an infinite-sites model using mutation events as alleles*

::

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='long')
   >>> import simuPOP as sim
   >>> 
   >>> def infSitesMutate(pop, param):
   ...     '''Apply an infinite mutation model'''
   ...     (startPos, endPos, rate) = param
   ...     # for each individual
   ...     for ind in pop.individuals():
   ...         # for each homologous copy of chromosomes
   ...         for p in range(2):
   ...             # using a geometric distribution to determine
   ...             # the first mutation location
   ...             loc = sim.getRNG().randGeometric(rate)
   ...             # if a mutation happens, record the mutated location
   ...             if startPos + loc < endPos:
   ...                 try:
   ...                     # find the first non-zero location
   ...                     idx = ind.genotype(p).index(0)
   ...                     # record mutation here
   ...                     ind.setAllele(startPos + loc, idx, ploidy=p)
   ...                 except:
   ...                     raise
   ...                     print('Warning: more than %d mutations have accumulated' % pop.totNumLoci())
   ...                     pass
   ...     return True
   ... 
   >>> pop = sim.Population(size=[2000], loci=[100])
   >>> pop.evolve(
   ...     initOps=sim.InitSex(),
   ...     preOps=[
   ...         # mutate in a 10Mb region at rate 1e-8
   ...         sim.PyOperator(func=infSitesMutate, param=(1, 10000000, 1e-8)),
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     gen = 100
   ... )
   100
   >>> # now, we get a sim.Population. Let us have a look at the 'alleles'.
   >>> # print the first five mutation locations
   >>> print(pop.individual(0).genotype()[:5])
   [1527502, 4774892, 7979220, 3671118, 395142]
   >>> # how many alleles are there (does not count 0)?
   >>> print(len(set(pop.genotype())) - 1)
   2700
   >>> # Allele count a simple count of alleles.
   >>> cnt = {}
   >>> for allele in pop.genotype():
   ...     if allele == 0:
   ...         continue
   ...     if allele in cnt:
   ...         cnt[allele] += 1
   ...     else:
   ...         cnt[allele] = 1
   ... 
   >>> # highest allele frequency?
   >>> print(max(cnt.values()) *0.5 / pop.popSize())
   0.05475

   now exiting runScriptInteractively...

`Download infiniteSites.py <infiniteSites.py>`_

All mutation models in simuPOP apply to existing alleles at pre-specified loci.
However, if the location of loci cannot be determined beforehand, it is
sometimes desired to create new loci as a result of mutation. A customized
operator can be used for this purpose (see Example :ref:`newOperator
<newOperator>`), but extra attention is needed to make sure that other operators
are applied to the correct loci because loci indexes will be changed with the
insertion of new loci. This technique could also be used to simulate mutations
over long sequences.


Recording and tracing individual mutants \*\*
---------------------------------------------

Mutation operators mutate alleles in place and by default do not generate any
output. If you are interested in knowing the source of each mutant, you can
specify an output stream and let the mutation operators dump details of each
mutation event, which consists of generation number, locus index, ploidy,
original allele, and mutated allele. If a list of information fields are
specified through parameter ``infoFields``, values at these information fields
will also be outputted (if they exist in the population. The default information
field is ``ind_id``, which allow you to record the ID of individuals harboring
the mutants.

Example :ref:`countMutants <countMutants>` demonstrates how to use this feature
to count the number of mutants at each locus. Instead of sending the output to a
file (e.g. ``output='>>mutants.txt'``), this example sends the output to a
Python function, which parses input string and counts the number of mutants at
each locus using a global dictionary variable. As we can see from the output,
because the KAlleleMutator uses a higher mutation rate (0.01) at locus 1 than
mutation rate (0.001) at locus 0, there are 10 times more mutants at the second
locus. There are about 3/4 mutations on the locus on chromosome X and 1/4
mutations on the locus on chromosome Y, for obvious reasons.

.. _countMutants:

**Example**: *Count number of mutants from mutator outputs*

::

   >>> import simuPOP as sim
   >>> from collections import defaultdict
   >>> # count number of mutants at each locus
   >>> counter = defaultdict(int)
   >>> def countMutants(mutants):
   ...     global counter
   ...     for line in mutants.split('\n'):
   ...         # a trailing \n will lead to an empty string
   ...         if not line:  
   ...             continue
   ...         (gen, loc, ploidy, a1, a2, id) = line.split('\t')
   ...         counter[int(loc)] += 1
   ... 
   >>> pop = sim.Population([5000]*3, loci=[2,1,1], infoFields='ind_id',
   ...     chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y])
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5]),
   ...         sim.IdTagger(),
   ...     ],
   ...     preOps=[
   ...         sim.KAlleleMutator(rates=[0.001] + [0.01]*3,
   ...             loci=range(4), k=100, output=countMutants),
   ...     ],
   ...     matingScheme=sim.RandomMating(
   ...         ops=[
   ...             sim.IdTagger(),
   ...             sim.MendelianGenoTransmitter()
   ...         ]),
   ...     gen = 10
   ... )
   10
   >>> print(counter.items())
   dict_items([(0, 308), (1, 2984), (2, 2319), (3, 768)])

   now exiting runScriptInteractively...

`Download countMutants.py <countMutants.py>`_


