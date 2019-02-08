Module :mod:`simuPOP.sampling`
==============================


Introduction
------------

Sampling, in simuPOP term, is the action of extracting individuals from a large,
potentially multi-generational, population according to certain criteria. the
:mod:`simuPOP.sampling` module provides several classes and functions and allows
you to define more complicated sampling schemes by deriving from its these
class. For example, you can use ``drawRandomSample(pop, size=100)`` to select
100 random individuals from a population, or use
``drawAffectedSibpairSample(pop, families=100)`` to select 100 pairs of affected
invididuals with their parents from a multi-generational population, or a age-
structured population with parents and offspring in the same generation.

The :mod:`simuPOP.sampling` module currently support random, case control,
affected sibpair, nuclear family and three-generation family sampling types, and
a combined sampling type that allows you to draw different types of samples. For
each sampling type ``X``, a sampler class and two functions ``DrawXSample`` and
``DrawXSamples`` are provided The first function returns a population with all
sampled individuals and the second function returns a list of sample
populations.

If you would like to define your own sampling type, you can derive your sampler
from one of the existing sampler classes. These sampler classes provide member
functions ``prepareSample``, ``drawSample`` and ``drawSamples`` and you
typically only need to extend ``prepareSample`` of an appropriate base class.


Sampling individuals randomly (class ``RandomSampler``, functions ``drawRandomSample`` and ``drawRandomSamples``)
-----------------------------------------------------------------------------------------------------------------

Functions ``drawRandomSample`` and ``drawRandomSamples`` draw random invidiauls
from a given population. If a simple number is given to parameter ``size``,
population structure will be ignored so individuals will be drawn from all
subpopulations. If a list of numbers are given, this function will draw
specified numbers of individuals from each subpopulation. This function does not
need parental information. If your population does not have an ID field, you
will not be able to locate extracted individuals in the original population.

Example :ref:`randomSample <randomSample>` demonstrates how to draw a random
sample from the whole population, and from each subpopulation. Because sample
populations keep the population structure of the source population (this might
change when parameter ``subPops`` is used, see a later section for details), we
can use ``sample.subPopSizes()`` to check how many individuals are sampled from
each subpopulation.

.. _randomSample:

**Example**: *Draw random samples from a structured population*

::

   >>> import simuPOP as sim
   >>> from simuPOP.sampling import drawRandomSample
   >>> pop = sim.Population([2000]*5, loci=1)
   >>> # sample from the whole population
   >>> sample = drawRandomSample(pop, sizes=500)
   >>> print(sample.subPopSizes())
   (104, 105, 110, 81, 100)
   >>> # sample from each subpopulation
   >>> sample = drawRandomSample(pop, sizes=[100]*5)
   >>> print(sample.subPopSizes())
   (100, 100, 100, 100, 100)

   now exiting runScriptInteractively...

`Download randomSample.py <randomSample.py>`_


Sampling cases and controls (class ``CaseControlSampler``, functions ``CaseControlSample`` and ``CaseControlSamples``)
----------------------------------------------------------------------------------------------------------------------

Functions ``drawCaseControlSample`` and ``drawCaseControlSamples`` draw cases
(affected individuals) and controls (unaffected invidiauls) from a given
population. If a simple number is given to parameter ``cases`` and ``controls``,
population structure will be ignored so individuals will be drawn from all
subpopulations. If a list of numbers are given, this function will draw
specified number of cases and controls from each subpopulation.

Example :ref:`caseControlSample <caseControlSample>` demonstrates how to draw
multiple case-control samples from a population, and perform case-control
assocition tests using the ``stat`` function.

.. _caseControlSample:

**Example**: *Draw case control samples from a population and perform association test*

::

   >>> import simuPOP as sim
   >>> from simuPOP.sampling import drawCaseControlSamples
   >>> pop = sim.Population([10000], loci=5)
   >>> sim.initGenotype(pop, freq=[0.2, 0.8])
   >>> sim.maPenetrance(pop, loci=2, penetrance=[0.11, 0.15, 0.20])
   >>> # draw multiple case control sample
   >>> samples = drawCaseControlSamples(pop, cases=500, controls=500, numOfSamples=5)
   >>> for sample in samples:
   ...     sim.stat(sample, association=range(5))
   ...     print(', '.join(['%.6f' % sample.dvars().Allele_ChiSq_p[x] for x in range(5)]))
   ... 
   0.694748, 0.333041, 0.001039, 0.078127, 0.774085
   0.261750, 0.954592, 0.031830, 0.737788, 0.865679
   0.954949, 0.371093, 0.092487, 0.622153, 0.075739
   0.654721, 0.433848, 0.002859, 0.696375, 0.956630
   0.439721, 1.000000, 0.069651, 0.471087, 0.238199

   now exiting runScriptInteractively...

`Download caseControlSample.py <caseControlSample.py>`_


Sampling Pedigrees (functions ``indexToID`` and ``plotPedigree``)
-----------------------------------------------------------------

If your sampling scheme involves parental information, you need to prepare your
population so that it has

* an ID field (usually ``'ind_id'``) that stores a unique ID for each
  individual.

* two information fields (usually ``'father_id'``, and ``'mother_id'``) that
  stores the ID of parents of each individual. Although simuPOP supports one-
  parent Pedigrees, this feature will not be discussed in this guide.

The preferred method to prepare such a population is to add information fields
``ind_id``, ``father_id`` and ``mother_id`` to a population and track ID based
Pedigrees during evolution. More specifically, you can use operators
:class:`IdTagger` and :class:`PedigreeTagger` to assign IDs and record parental
IDs of each offspring during mating. This method supports age-structured
population when parents and offspring can be stored in the same generation.

You can also use information fields ``father_idx`` and ``mother_idx`` and
operator :class:`ParentsTagger` to track indexes of parents in the parental
generations. Before sampling, you can use function
:func:`~simuPOP.sampling.indexToID` to add needed information fields and convert
index based parental relationship to ID based relationshop. Because parents have
to stay in ancestral generations, this method does not support age-structured
population.

If you have R and rpy installed on your system, you can install the ``kinship``
library of R and use it to analyze Pedigree. The :mod:`simuPOP.sampling` module
provides a function ``plotPedigree`` to use this library to plot Pedigrees.
Example :ref:`plotPedigree <plotPedigree>` demonstrates how to use function
sampling.indexToID to prepare a pedigree and how to use sampling.DrawPedigree to
plot it.

Figure :ref:`fig_Pedigree <fig_Pedigree>` plots a small three-generational
population with 15 individuals at each generation. It is pretty clear that
random mating produces bad pedigree structure because it is common that one
parent would have multiple spouses.


Sampling affected sibpairs (class ``AffectedSibpairSampler``, functions ``drawAffectedSibpairSample(s)``)
---------------------------------------------------------------------------------------------------------

An affected sibpair family consists of two parents and their affected offspring.
Such families are useful in linkage analysis because of high likelihood of
shared disease predisposing alleles between siblings. :mod:`simuPOP.sampling`
module provides functions ``drawAffectedSibpairSample`` and
``drawAffectedSibpairSamples`` to draw such families from a population. Example
:ref:`sampleAffectedSibpair <sampleAffectedSibpair>` draws two affected sibpair
from the pedigree created in Example :ref:`plotPedigree <plotPedigree>`, with
samples plotted in Figure :ref:`fig_affectedSibpair <fig_affectedSibpair>`.

.. _sampleAffectedSibpair:

**Example**: *Draw affected sibpairs from a population*

::

   >>> import simuPOP as sim
   >>> from simuPOP.sampling import indexToID
   >>> pop = sim.Population(size=15, loci=5, infoFields=['father_idx', 'mother_idx'], ancGen=2)
   >>> pop.evolve(
   ...     preOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.7, 0.3]),
   ...     ],
   ...     matingScheme=sim.RandomMating(numOffspring=(sim.UNIFORM_DISTRIBUTION, 2, 4),
   ...         ops=[sim.MendelianGenoTransmitter(), sim.ParentsTagger()]),
   ...     postOps=sim.MaPenetrance(loci=3, penetrance=(0.1, 0.4, 0.7)),
   ...     gen = 5
   ... )
   5
   >>> indexToID(pop, reset=True)
   >>> # three information fields were added
   >>> print(pop.infoFields())
   ('father_idx', 'mother_idx', 'ind_id', 'father_id', 'mother_id')
   >>> # save this population for future use
   >>> pop.save('log/pedigree.pop')
   >>> 
   >>> from simuPOP.sampling import drawAffectedSibpairSample
   >>> pop = sim.loadPopulation('log/pedigree.pop')
   >>> sample = drawAffectedSibpairSample(pop, families=2)
   Warning: number of requested Pedigrees 2 is greater than what exists (0).
   Warning: not enough non-overlapping Pedigrees are found (requested 2, found 0).

   now exiting runScriptInteractively...

`Download sampleAffectedSibpair.py <sampleAffectedSibpair.py>`_


Sampling nuclear families (class ``NuclearFamilySampler``, functions ``drawNuclearFamilySample`` and ``drawNuclearFamilySample``\ s)
------------------------------------------------------------------------------------------------------------------------------------

A nuclear family consists of two parents and their offspring. Functions
``drawNuclearFamilySample`` and ``drawNuclearFamilySamples`` to draw such
families from a population, with restrictions on number of offspring, number of
affected parents and number of affected offspring. Although fixed numbers could
be given, a range with minimal and maximal acceptable numbers are usually
provided. Example :ref:`sampleNuclearFamily <sampleNuclearFamily>` draws two
nuclear families from the pedigree created in Example :ref:`plotPedigree
<plotPedigree>`. The samples are plotted in Figure :ref:`fig_nuclearFamily
<fig_nuclearFamily>`.

.. _sampleNuclearFamily:

**Example**: *Draw nuclear families from a population*

::

   >>> import simuPOP as sim
   >>> from simuPOP.sampling import drawNuclearFamilySample
   >>> pop = sim.loadPopulation('log/pedigree.pop')
   >>> sample = drawNuclearFamilySample(pop, families=2, numOffspring=(2,4),
   ...     affectedParents=(1,2), affectedOffspring=(1, 3))
   Warning: number of requested Pedigrees 2 is greater than what exists (0).
   Warning: not enough non-overlapping Pedigrees are found (requested 2, found 0).
   >>> # try to separate two families?
   >>> sample.asPedigree()
   >>> #= sim.Pedigree(sample, loci=sim.ALL_AVAIL, infoFields=sim.ALL_AVAIL)
   >>> sample.addInfoFields('ped_id')
   >>> # return size of families
   >>> sz = sample.identifyFamilies(pedField='ped_id')
   >>> print(sz)
   ()
   >>> ped1 = sample.extractIndividuals(IDs=0, idField='ped_id')
   >>> # print the ID of all individuals in the first pedigree
   >>> print([ind.ind_id for ind in ped1.allIndividuals()])
   []

   now exiting runScriptInteractively...

`Download sampleNuclearFamily.py <sampleNuclearFamily.py>`_


Sampling three-generation families (class ``ThreeGenFamilySampler``, functions ``drawThreeGenFamilySample and drawThreeGenFamilySamples``)
------------------------------------------------------------------------------------------------------------------------------------------

A three-generation family consists of two parents, their common offspring,
offspring's spouses, and their common offspring (grandchidren). individuals in
sampled families have either no or two parents. Functions
``drawThreeGenFamilySample`` and ``drawThreeGenFamilySamples`` to draw such
families from a population, with restrictions on number of offspring, total
number of individuals and number of affected individuals in the Pedigree. These
parameters (``numOffspring``, ``pedSize`` and ``numAffected``) could be a fixed
number of a range with minimal and maximal acceptable numbers. Example
:ref:`sampleNuclearFamily <sampleNuclearFamily>` draws two three generation
families from the pedigree created in Example :ref:`plotPedigree
<plotPedigree>`. The samples are plotted in Figure :ref:`fig_nuclearFamily
<fig_nuclearFamily>`.

.. _sampleThreeGenFamily:

**Example**: *Draw three-generation families  from a population*

::

   >>> import simuPOP as sim
   >>> from simuPOP.sampling import drawThreeGenFamilySample
   >>> pop = sim.loadPopulation('log/pedigree.pop')
   >>> sample = drawThreeGenFamilySample(pop, families=2, numOffspring=(1, 3),
   ...     pedSize=(8, 15), numOfAffected=(2, 5))

   now exiting runScriptInteractively...

`Download sampleThreeGenFamily.py <sampleThreeGenFamily.py>`_


Sampling different types of samples (class ``CombinedSampler``, functions ``drawCombinedSample`` and ``drawCombinedSamples``)
-----------------------------------------------------------------------------------------------------------------------------

Samples in real world studies sometimes do not have uniform types so it is
useful to draw samples of different types from the same population. Although it
is possible to draw samples using different functions and combine them, handling
of overlapping individuals, namely individuals who are chosen by multiple
samplers, can be a headache. The combined sampler of :mod:`simuPOP.sampling` is
designed to overcome this problem. This sampler takes a list of sampler objects
and apply them to a population sequentially. The extracted sample will not have
overlapping individuals.

Example :ref:`combinedSampling <combinedSampling>` draws an affected sibpair
family and a nuclear family from the pedigree created in Example
:ref:`plotPedigree <plotPedigree>`. The samples are plotted in Figure
:ref:`combinedSampling <combinedSampling>`.

.. _combinedSampling:

**Example**: *Draw different types of samples from a population*

::

   >>> import simuPOP as sim
   >>> from simuPOP.sampling import drawCombinedSample, AffectedSibpairSampler, NuclearFamilySampler
   >>> pop = sim.loadPopulation('log/pedigree.pop')
   >>> sample = drawCombinedSample(pop, samplers = [
   ...     AffectedSibpairSampler(families=1),
   ...     NuclearFamilySampler(families=1, numOffspring=(2,4), affectedParents=(1,2), affectedOffspring=(1,3))
   ...     ])
   Warning: number of requested Pedigrees 1 is greater than what exists (0).
   Warning: not enough non-overlapping Pedigrees are found (requested 1, found 0).
   Warning: number of requested Pedigrees 1 is greater than what exists (0).
   Warning: not enough non-overlapping Pedigrees are found (requested 1, found 0).

   now exiting runScriptInteractively...

`Download combinedSampling.py <combinedSampling.py>`_


Sampling from subpopulations and virtual subpopulations \*
----------------------------------------------------------

Virtual subpopulations (VSPs) could be specified in the ``subPops`` parameter of
sampling classes and functions. This can be used to limit your samples to
individuals with certain properties. For example, you may want to match the age
of cases and controls in a case-control association study by selecting your
samples from a certain age group. For examples, Example :ref:`samplingVSP
<samplingVSP>` draws 500 cases and 500 controls from two a VSP with individual
ages between 40 and 60.

.. _samplingVSP:

**Example**: *Draw samples from a virtual subpopulation.*

::

   >>> import simuPOP as sim
   >>> # create an age-structured population with a disease
   >>> import random
   >>> pop = sim.Population(10000, loci=10, infoFields='age')
   >>> sim.initGenotype(pop, freq=[0.3, 0.7])
   >>> sim.initInfo(pop, lambda: random.randint(0, 70), infoFields='age')
   >>> pop.setVirtualSplitter(sim.InfoSplitter(cutoff=(40, 60), field='age'))
   >>> sim.maPenetrance(pop, loci=5, penetrance=(0.1, 0.2, 0.3))
   >>> #
   >>> from simuPOP.sampling import drawCaseControlSample
   >>> sample = drawCaseControlSample(pop, cases=500, controls=500, subPops=[(0,1)])
   >>> ageInSample = sample.indInfo('age')
   >>> print(min(ageInSample), max(ageInSample))
   40.0 59.0

   now exiting runScriptInteractively...

`Download samplingVSP.py <samplingVSP.py>`_

If a list of sample sizes is given, specified number of samples will be drawn
from each subpopulation. For example, if you have an age-structured population
when individuals with different ages have different risk to a disease, you might
want to draw affected individuals from different age groups and perform
association analyses. Function ``drawCaseControlSample`` cannot be used because
both groups are affected, but you can ``drawRandomSample`` from two VSPs defined
by age. Example :ref:`samplingSeparateVSPs <samplingSeparateVSPs>` demonstrates
how to use this method.

.. _samplingSeparateVSPs:

**Example**: *Sampling separately from different virtual subpopulations*

::

   >>> import simuPOP as sim
   >>> # create an age-structured population with a disease
   >>> import random
   >>> pop = sim.Population(10000, loci=10, infoFields='age')
   >>> sim.initGenotype(pop, freq=[0.3, 0.7])
   >>> sim.initInfo(pop, lambda: random.randint(0, 70), infoFields='age')
   >>> pop.setVirtualSplitter(sim.InfoSplitter(cutoff=(20, 40), field='age'))
   >>> # different age group has different penetrance
   >>> sim.maPenetrance(pop, loci=5, penetrance=(0.1, 0.2, 0.3), subPops=[(0,1)])
   >>> sim.maPenetrance(pop, loci=5, penetrance=(0.2, 0.4, 0.6), subPops=[(0,2)])
   >>> # count the number of affected individuals in each group
   >>> sim.stat(pop, numOfAffected=True, subPops=[(0,1), (0,2)], vars='numOfAffected_sp')
   >>> print(pop.dvars((0,1)).numOfAffected, pop.dvars((0,2)).numOfAffected)
   668 2215
   >>> #
   >>> from simuPOP.sampling import drawRandomSample
   >>> sample = drawRandomSample(pop, sizes=[500, 500], subPops=[(0,1), (0,2)])
   >>> # virtual subpopulations are rearranged to different subpopulations.
   >>> print(sample.subPopSizes())
   (500, 500)

   now exiting runScriptInteractively...

`Download samplingSeparateVSPs.py <samplingSeparateVSPs.py>`_


