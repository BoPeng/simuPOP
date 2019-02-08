.. _sec_Genotype_transmitters:

Genotype transmitters
=====================


Generic genotype transmitters (operators :class:`GenoTransmitter`, :class:`CloneGenoTransmitter`, :class:`MendelianGenoTransmitter`, :class:`SelfingGenoTransmitter`, :class:`HaplodiploidGenoTransmitter`, and :class:`MitochondrialGenoTransmitter`) \*
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

A number of during-mating operators are defined to transmit genotype from
parent(s) to offspring. They are rarely used or even seen directly because they
are used as genotype transmitters of mating schemes.

* :class:`GenoTransmitter`: This genotype transmitter is usually used by
  customized genotype transmitters because it provides some utility functions that
  are more efficient than their Pythonic counterparts.

* :class:`CloneGenoTransmitter`: Copy all genotype on non-customized chromosomes
  from a parent to an offspring. It also copies parental sex to the offspring
  because sex can be genotype determined. This genotype transmitter is used by
  mating scheme :class:`CloneMating`. This genotype transmitter can be applied to
  populations of **any ploidy** type. If you would like to copy part of the
  chromosomes, or customized chromosomes, a parameter chroms could be used to
  specify chromosomes to copy.

* :class:`MendelianGenoTransmitter`: Copy genotypes from two parents (a male and
  a female) to an offspring following Mendel's laws, used by mating scheme
  ``RandomMating.``\ This genotype transmitter can only be applied to **diploid**
  populations.

* :class:`SelfingGenoTransmitter`: Copy genotypes from one parent to an
  offspring using self-fertilization, used by mating scheme :class:`SelfMating`.
  This genotype transmitter can only be applied to **diploid** populations.

* :class:`HaplodiploidGenoTransmitter`: Set genotype to male and female
  offspring differently in a haplodiploid population, used by mating scheme
  :class:`HaplodiploidMating`. This genotype transmitter can only be applied to
  **haplodiploid** populations.

* :class:`MitochondrialGenoTransmitter`: Treat a single mitochondrial
  chromosome, or all customized chromosomes, or specified chromosomes as
  mitochondrial chromosomes and transmit maternal mitochondrial chromosomes
  randomly to an offspring. This genotype transmitter can be applied to
  populations of **any ploidy** type. It trasmits the first homologous copy of
  chromosomes maternally and clears alleles on other homologous copies of
  chromosomes of an offspring.


Recombination (Operator :class:`Recombinator`)
----------------------------------------------

The generic genotype transmitters do not handle genetic recombination. A
genotype transmitter :class:`Recombinator` is provided for such purposes, and
can be used with :class:`RandomMating` and :class:`SelfMating` (replace
:class:`MendelianGenoTransmitter` and :class:`SelfingGenoTransmitter` used in
these mating schemes).

Recombination rate is implemented **between adjacent markers**. There can be
only one recombination event between adjacent markers no matter how far apart
they are located on a chromosome. In practise, a :class:`Recombinator` goes
along chromosomes and determine, between each adjacent loci, whether or not a
recombination happens.

Recombination rates could be specified in the following ways:

#. If a single recombination rate is specified through paramter ``rate``\ s, it
   will be the recombination rate between all adjacent loci, regardless of loci
   position.

#. If recombination happens only after certain loci, you can specify these loci
   using parameter ``loci``. For example,   ::

      Recombinator(rates=0.1, loci=[2, 5])

   recombines a chromosome only **after** loci 2 (between 2 and 3) and 5 (between 5
   and 6).

#. If parameter ``loci`` is given with a list of loci, different recombination
   rate can be given to each of them. The two lists should have the same length.
   For example  ::

      Recombinator(rates=[0.1, 0.05], loci=[2, 5])

   uses two different recombination rates after loci 2 and 5.

#. If parameter ``loci`` is not given (default to ``loci=ALL_AVAIL``) but a list
   of recombination rates is assigned, the rates will be assigned to each locus.
   The length of prameter ``rates`` should equal to total number of loci but the
   recombiantion rates for the locus at the end of each chromosome will be ignored
   (assumed to be 0.5). For example  ::

      Recombinator(rates=[0.1]*5 + [0.2]*5)

   uses two different recombination rates for two chromosomes with 5 loci.

#. If recombination rates vary across your chromosomes, a long list of ``rate``
   and ``loci`` may be needed to specify recombination rates one by one. An
   alternative method is to specify a **recombination intensity**. Recombination
   rate between two adjacent loci is calculated as the product of this intensity
   and distance between them. For example, if you apply operator  ::

      Recombinator(intensity=0.1)

   to a population  ::

      Population(size=100, loci=[4], lociPos=[0.1, 0.2, 0.4, 0.8])

   The recombination rates between adjacent markers will be ``0.1*0.1``,
   ``0.1*0.2`` and ``0.1*0.4`` respectively.

.. _recRate:

**Example**: *Genetic recombination at all and selected loci*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(sim.Population(size=[1000], loci=[100]),
   ...     rep=2)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[0]*100 + [1]*100)
   ...     ],
   ...     matingScheme=sim.RandomMating(ops = [
   ...         sim.Recombinator(rates=0.01, reps=0),
   ...         sim.Recombinator(rates=[0.01]*10, loci=range(50, 60), reps=1),
   ...     ]),
   ...     postOps=[
   ...         sim.Stat(LD=[[40, 55], [60, 70]]),
   ...         sim.PyEval(r'"%d:\t%.3f\t%.3f\t" % (rep, LD_prime[40][55], LD_prime[60][70])'),
   ...         sim.PyOutput('\n', reps=-1)
   ...     ],
   ...     gen = 5
   ... )
   0:	0.741	0.806	1:	0.904	1.000	
   0:	0.658	0.715	1:	0.882	1.000	
   0:	0.491	0.668	1:	0.843	1.000	
   0:	0.435	0.610	1:	0.818	1.000	
   0:	0.383	0.567	1:	0.763	1.000	
   (5, 5)

   now exiting runScriptInteractively...

`Download recRate.py <recRate.py>`_

Example :ref:`recRate <recRate>` demonstrates how to specify recombination rates
for all loci or for specified loci. In this example, two replicates of a
population are evolved, subject to two different Recombinators. The first
Recombinator applies the same recombination rate between all adjacent loci, and
the second Recombinator recombines only after loci 50 - 59. Because there is no
recombination event between loci 60 and 70 for the second replicate, linkage
disequilibrium values between these two loci does not decrease as what happens
in the first replicate.

.. _recIntensity:

**Example**: *Genetic recombination rates specified by intensity*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(size=[1000], loci=3, lociPos=[0, 1, 1.1])
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[0]*3 + [1]*3)
   ...     ],
   ...     matingScheme=sim.RandomMating(ops=sim.Recombinator(intensity=0.01)),
   ...     postOps=[
   ...         sim.Stat(LD=[[0, 1], [1, 2]]),
   ...         sim.PyEval(r'"%.3f\t%.3f\n" % (LD_prime[0][1], LD_prime[1][2])', step=10)
   ...     ],
   ...     gen = 50
   ... )
   0.988	0.998
   0.912	0.996
   0.836	0.991
   0.896	0.982
   0.814	0.991
   50

   now exiting runScriptInteractively...

`Download recIntensity.py <recIntensity.py>`_

Example :ref:`recIntensity <recIntensity>` demonstrates the use of the
``intensity`` parameter. In this example, the distances between the first two
loci and the latter two loci are 1 and 0.1 respectively. This leads
recombination rates 0.01 and 0.001 respectively with a recombination intensity
0.01. Consequently, LD between the first two loci decay much faster than the
latter two.

If more advanced recombination model is desired, a customized genotype
transmitter can be used. For example, Example :ref:`sexSpecificRec
<sexSpecificRec>` uses two Recombinators to implement sex-specific
recombination.

.. note::

   Both loci positions and recombination intensity are unitless. You can assume
   different unit for loci position and recombination intensity as long as the
   resulting recombination rate makes sense.


Gene conversion (Operator :class:`Recombinator`) \*
---------------------------------------------------

simuPOP uses the Holliday junction model to simulate gene conversion. This model
treats recombination and conversion as a unified process. The key features of
this model is

* Two (out of four) chromatids pair and a single strand cut is made in each
  chromatid

* Strand exchange takes place between the chromatids

* Ligation occurs yielding two completely intact DNA molecules

* Branch migration occurs, giving regions of heteroduplex DNA

* Resolution of the Holliday junction gives two DNA molecules with heteroduplex
  DNA. Depending upon how the holliday junction is resolved, we either observe no
  exchange of flanking markers, or an exchange of flanking markers. The former
  forms a conversion event, which can be considered as a double recombination.

In practise, gene conversion can be considered as a double recombination event.
That is to say, when a recombination event happens, it has certain probability
to trigger a second recombination event along the chromosome. The distance
between the two locations where recombination events happen is the tract length
of this conversion event.

The probability at which gene conversion happens, and how tract length is
determined is specify using parameter ``convMode`` of a Recombinator. This
parameter can be

* ``NoConversion`` No gene conversion. (default)

* ``(NUM_MARKERS, prob, N)`` Convert a fixed number ``N`` of markers at
  probability ``prob``.

* ``(TRACT_LENGTH, prob, N)`` Convert a fixed length ``N`` of chromosome regions
  at probability ``prob``. This can be used when markers are not equally spaced on
  chromosomes.

* ``(GEOMETRIC_DISTRIBUTION, prob, p)`` When a conversion event happens at
  probability ``prob``, convert a random number of markers, with a geometric
  distribution with parameter ``p``.

* ``(EXPONENTIAL_DISTRIBUTION, prob, p)`` When a conversion event happens at
  probability ``prob``, convert a random length of chromosome region, using an
  exponential distribution with parameter ``p``.

Note that

* If tract length is determined by length (``TractLength`` or
  ``ExponentialDistribution``), the starting point of the flanking region is
  uniformly distributed between marker :math:`i-1` and :math:`i`, if the
  recombination happens at marker :math:`i`. That is to say, it is possible that
  no marker is converted with a positive tract length.

* A conversion event will act like a recombination event if its flanking region
  exceeds the end of a chromosome, or if another recombination event happens
  before the end of the flanking region.

Example :ref:`conversion <conversion>` compares two Recombinators. The first
Recombinator is a regular Recombinator that recombine between loci 50 and 51.
The second Recombinator is a conversion operator because every recombination
event will become a conversion event (prob=1). Because a second recombination
event will surely happen between loci 60 and 61, there will be either no or
double recombination events between loci 40, 70. LD between these two loci
therefore does not decrease, although LD between locus 55 and these two loci
will decay.

.. _conversion:

**Example**: *Gene conversion*

::

   >>> import simuPOP as sim
   >>> simu = sim.Simulator(sim.Population(size=[1000], loci=[100]),
   ...     rep=2)
   >>> simu.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(genotype=[0]*100 + [1]*100)
   ...     ],
   ...     matingScheme=sim.RandomMating(ops=[
   ...         sim.Recombinator(rates=0.01, loci=50, reps=0),
   ...         sim.Recombinator(rates=0.01, loci=50, reps=1, convMode=(sim.NUM_MARKERS, 1, 10)),
   ...     ]),
   ...     postOps=[
   ...         sim.Stat(LD=[[40, 55], [40, 70]]),
   ...         sim.PyEval(r'"%d:\t%.3f\t%.3f\t" % (rep, LD_prime[40][55], LD_prime[40][70])'),
   ...         sim.PyOutput('\n', reps=-1)
   ...     ],
   ...     gen = 5
   ... )
   0:	0.988	0.988	1:	0.980	1.000	
   0:	0.982	0.982	1:	0.982	1.000	
   0:	0.982	0.982	1:	0.974	1.000	
   0:	0.974	0.974	1:	0.954	1.000	
   0:	0.960	0.960	1:	0.940	1.000	
   (5, 5)

   now exiting runScriptInteractively...

`Download conversion.py <conversion.py>`_


Tracking all recombination events \*\*
--------------------------------------

To understand the evolutionary history of a simulated population, it is
sometimes needed to track down all ancestral recombination events. In order to
do that, you will first need to give an unique ID to each individual so that you
could make sense of the dumped recombination events. Although this is routinely
done using operator :class:`IdTagger` (see example :ref:`IdTagger <IdTagger>`
for details), it is a little tricky here because you need to place the during-
mating :class:`IdTagger` before a :class:`Recombinator` in the ``ops`` parameter
of a mating scheme so that offspring ID could be set and outputted correctly.

After setting the name of the ID field (usually ``ind_id``) to the ``infoField``
parameter of a :class:`Recombinator`, it can dump a list of recombinatin events
(loci after which recombinatin events happened) for each set of homologous
chromosomes of an offspring. Each line is in the format of

::

   offspringID parentID startingPloidy rec1 rec2 ....

Example :ref:`trackRec <trackRec>` gives an example how the output looks like.

.. _trackRec:

**Example**: *Tracking all recombination events*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, loci=[1000, 2000], infoFields='ind_id')
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.IdTagger(),
   ...     ],
   ...     matingScheme=sim.RandomMating(ops = [
   ...         sim.IdTagger(),
   ...         sim.Recombinator(rates=0.001, output='>>rec.log', infoFields='ind_id')]),
   ...     gen = 5
   ... )
   5
   >>> rec = open('rec.log')
   >>> # print the first three lines of the log file
   >>> print(''.join(rec.readlines()[:4]))
   1001 642 0 381 999 1490
   1001 250 1 908 999 1315 2134
   1002 847 1 999
   1002 91 0 975 999 1245 2546


   now exiting runScriptInteractively...

`Download trackRec.py <trackRec.py>`_


