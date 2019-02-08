Module :mod:`simuPOP.sampling`
==============================


.. module:: simuPOP.sampling

This module provides classes and functions that could be used to draw samples
from a simuPOP population. These functions accept a list of parameters such
as ``subPops`` ((virtual) subpopulations from which samples will be drawn) and
``numOfSamples`` (number of samples to draw) and return a list of populations. Both
independent individuals and dependent individuals (Pedigrees) are supported.

Independent individuals could be drawn from any Population. pedigree
information is not necessary and is usually ignored. Unique IDs are not needed
either although such IDs could help you identify samples in the parent
Population.

Pedigrees could be drawn from multi-generational populations or age-structured
populations. All individuals are required to have a unique ID (usually tracked
by operator ``IdTagger`` and are stored in information field ``ind_id``).
Parents of individuals are usually tracked by operator ``PedigreeTagger`` and
are stored in information fields ``father_id`` and ``mother_id``. If parental
information is tracked using operator ``ParentsTagger`` and information fields
``father_idx`` and ``mother_idx``, a function ``sampling.indexToID`` can be
used to convert index based pedigree to ID based Pedigree. Note that
``ParentsTagger`` can not be used to track Pedigrees in age-structured
populations because they require parents of each individual resides in a
parental generation.

All sampling functions support virtual subpopulations through parameter
``subPops``, although sample size specification might vary. This feature
allows you to draw samples with specified properties. For example, you
could select only female individuals for cases of a female-only disease,
or select individuals within certain age-range. If you specify a list
of (virtual) subpopulations, you are usually allowed to draw certain
number of individuals from each subpopulation.


class BaseSampler
-----------------

.. class:: BaseSampler

   A sampler extracts individuals from a simuPOP population and return them
   as separate populations. This base class defines the common interface of
   all sampling classes, including how samples prepared and returned.

   .. method:: BaseSampler.BaseSampler(subPops=ALL_AVAIL)

      Create a sampler with parameter ``subPops``, which will be used
      to prepare population for sampling. ``subPops`` should be a list of
      (virtual) subpopulations from which samples are drawn. The default
      value is ALL_AVAIL, which means all available subpopulations of a
      Population.

   .. method:: BaseSampler.drawSample(pop)

      Draw and return a sample.

   .. method:: BaseSampler.drawSamples(pop, numOfSamples)

      Draw multiple samples and return a list of populations.

   .. method:: BaseSampler.prepareSample(pop, rearrange)

      Prepare passed population object for sampling according to parameter
      ``subPops``. If samples are drawn from the whole population, a
      Population will be trimmed if only selected (virtual) subpopulations
      are used. If samples are drawn separately from specified subpopulations,
      Population ``pop`` will be rearranged (if ``rearrange==True``) so that
      each subpoulation corresponds to one element in parameter ``subPops``.



class RandomSampler
-------------------

.. class:: RandomSampler

   A sampler that draws individuals randomly.

   .. method:: RandomSampler.RandomSampler(sizes, subPops=ALL_AVAIL)

      Creates a random sampler with specified number of individuals.

   .. method:: RandomSampler.drawSample(input_pop)

      Draw a random sample from passed population.

   .. method:: RandomSampler.drawSamples(pop, numOfSamples)

      Draw multiple samples and return a list of populations.

   .. method:: RandomSampler.prepareSample(pop, rearrange)

      Prepare passed population object for sampling according to parameter
      ``subPops``. If samples are drawn from the whole population, a
      Population will be trimmed if only selected (virtual) subpopulations
      are used. If samples are drawn separately from specified subpopulations,
      Population ``pop`` will be rearranged (if ``rearrange==True``) so that
      each subpoulation corresponds to one element in parameter ``subPops``.



Function drawRandomSample
-------------------------


.. function:: drawRandomSample(pop, sizes, subPops=ALL_AVAIL)

   Draw ``sizes`` random individuals from a population. If a single ``sizes``
   is given, individuals are drawn randomly from the whole population or
   from specified (virtual) subpopulations (parameter ``subPops``). Otherwise,
   a list of numbers should be used to specify number of samples from each
   subpopulation, which can be all subpopulations if ``subPops=ALL_AVAIL``
   (default), or from each of the specified (virtual) subpopulations. This
   function returns a population with all extracted individuals.


Function drawRandomSamples
--------------------------


.. function:: drawRandomSamples(pop, sizes, numOfSamples=1, subPops=ALL_AVAIL)

   Draw ``numOfSamples`` random samples from a population and return a list of
   populations. Please refer to function ``drawRandomSample`` for more details
   about parameters ``sizes`` and ``subPops``.


class CaseControlSampler
------------------------

.. class:: CaseControlSampler

   A sampler that draws affected and unaffected individuals randomly.

   .. method:: CaseControlSampler.CaseControlSampler(cases, controls, subPops=ALL_AVAIL)

      Ceates a case-control sampler with specified number of cases
      and controls.

   .. method:: CaseControlSampler.drawSample(input_pop)

      Draw a case control sample

   .. method:: CaseControlSampler.drawSamples(pop, numOfSamples)

      Draw multiple samples and return a list of populations.

   .. method:: CaseControlSampler.prepareSample(input_pop)

      Find out indexes all affected and unaffected individuales.



Function drawCaseControlSample
------------------------------


.. function:: drawCaseControlSample(pop, cases, controls, subPops=ALL_AVAIL)

   Draw a case-control samples from a population with ``cases``
   affected and ``controls`` unaffected individuals. If single ``cases`` and
   ``controls`` are given, individuals are drawn randomly from the whole
   Population or from specified (virtual) subpopulations (parameter
   ``subPops``). Otherwise, a list of numbers should be used to specify
   number of cases and controls from each subpopulation, which can be all
   subpopulations if ``subPops=ALL_AVAIL`` (default), or from each of the
   specified (virtual) subpopulations. This function returns a population with
   all extracted individuals.


Function drawCaseControlSamples
-------------------------------


.. function:: drawCaseControlSamples(pop, cases, controls, numOfSamples=1, subPops=ALL_AVAIL)

   Draw ``numOfSamples`` case-control samples from a population with ``cases``
   affected and ``controls`` unaffected individuals and return a list of
   populations. Please refer to function ``drawCaseControlSample`` for a
   detailed descriptions of parameters.


class PedigreeSampler
---------------------

.. class:: PedigreeSampler

   The base class of all pedigree based sampler.

   .. method:: PedigreeSampler.PedigreeSampler(families, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id')

      Creates a pedigree sampler with parameters
      
      families
          number of families. This can be a number or a list of numbers. In
          the latter case, specified families are drawn from each
          subpopulation.
      
      subPops
          A list of (virtual) subpopulations from which samples are drawn.
          The default value is ALL_AVAIL, which means all available
          subpopulations of a population.

   .. method:: PedigreeSampler.drawSample(input_pop)

      Randomly select Pedigrees

   .. method:: PedigreeSampler.drawSamples(pop, numOfSamples)

      Draw multiple samples and return a list of populations.

   .. method:: PedigreeSampler.family(id)

      Get the family of individual with id.

   .. method:: PedigreeSampler.prepareSample(pop, loci=[], infoFields=[], ancGens=True)

      Prepare self.pedigree, some pedigree sampler might need additional loci and
      information fields for this sampler.



class AffectedSibpairSampler
----------------------------

.. class:: AffectedSibpairSampler

   A sampler that draws a nuclear family with two affected offspring.

   .. method:: AffectedSibpairSampler.AffectedSibpairSampler(families, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id')

      Initialize an affected sibpair sampler.

   .. method:: AffectedSibpairSampler.drawSample(input_pop)

      Randomly select Pedigrees

   .. method:: AffectedSibpairSampler.drawSamples(pop, numOfSamples)

      Draw multiple samples and return a list of populations.

   .. method:: AffectedSibpairSampler.family(id)

      Return id, its spouse and their children

   .. method:: AffectedSibpairSampler.prepareSample(input_pop)

      Find the father or all affected sibpair families



Function drawAffectedSibpairSample
----------------------------------


.. function:: drawAffectedSibpairSample(pop, families, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id')

   Draw affected sibpair samples from a population. If a single
   ``families`` is given, affected sibpairs and their parents are drawn
   randomly from the whole population or from specified (virtual)
   subpopulations (parameter ``subPops``). Otherwise, a list of numbers should
   be used to specify number of families from each subpopulation, which can be
   all subpopulations if ``subPops=ALL_AVAIL`` (default), or from each of the
   specified (virtual) subpopulations. This function returns a population that
   contains extracted individuals.


Function drawAffectedSibpairSamples
-----------------------------------


.. function:: drawAffectedSibpairSamples(pop, families, numOfSamples=1, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id')

   Draw ``numOfSamples`` affected sibpair samplesa from population ``pop`` and
   return a list of populations. Please refer to function
   ``drawAffectedSibpairSample`` for a description of other parameters.


class NuclearFamilySampler
--------------------------

.. class:: NuclearFamilySampler

   A sampler that draws nuclear families with specified number of affected
   parents and offspring.

   .. method:: NuclearFamilySampler.NuclearFamilySampler(families, numOffspring, affectedParents=0, affectedOffspring=0, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id')

      Creates a nuclear family sampler with parameters
      
      families
          number of families. This can be a number or a list of numbers. In the latter
          case, specified families are drawn from each subpopulation.
      
      numOffspring
          number of offspring. This can be a fixed number or a range [min, max].
      
      affectedParents
          number of affected parents. This can be a fixed number or a range [min, max].
      
      affectedOffspring
          number of affected offspring. This can be a fixed number of a range [min, max].
      
      subPops
          A list of (virtual) subpopulations from which samples are drawn.
          The default value is ALL_AVAIL, which means all available
          subpopulations of a population.

   .. method:: NuclearFamilySampler.drawSample(input_pop)

      Randomly select Pedigrees

   .. method:: NuclearFamilySampler.drawSamples(pop, numOfSamples)

      Draw multiple samples and return a list of populations.

   .. method:: NuclearFamilySampler.family(id)

      Return id, its spouse and their children

   .. method:: NuclearFamilySampler.prepareSample(input_pop)

      Prepare self.pedigree, some pedigree sampler might need additional loci and
      information fields for this sampler.



Function drawNuclearFamilySample
--------------------------------


.. function:: drawNuclearFamilySample(pop, families, numOffspring, affectedParents=0, affectedOffspring=0, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id')

   Draw nuclear families from a population. Number of offspring, number of
   affected parents and number of affected offspring should be specified using
   parameters ``numOffspring``, ``affectedParents`` and ``affectedOffspring``,
   which can all be a single number, or a range ``[a, b]`` (``b`` is incldued).
   If a single ``families`` is given, Pedigrees are drawn randomly from the
   whole population or from specified (virtual) subpopulations (parameter
   ``subPops``). Otherwise, a list of numbers should be used to specify
   numbers of families from each subpopulation, which can be all
   subpopulations if ``subPops=ALL_AVAIL`` (default), or from each of the
   specified (virtual) subpopulations. This function returns a population that
   contains extracted individuals.


Function drawNuclearFamilySamples
---------------------------------


.. function:: drawNuclearFamilySamples(pop, families, numOffspring, affectedParents=0, affectedOffspring=0, numOfSamples=1, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id')

   Draw ``numOfSamples`` affected sibpair samplesa from population ``pop`` and
   return a list of populations. Please refer to function
   ``drawNuclearFamilySample`` for a description of other parameters.


class ThreeGenFamilySampler
---------------------------

.. class:: ThreeGenFamilySampler

   A sampler that draws three-generation families with specified pedigree
   size and number of affected individuals.

   .. method:: ThreeGenFamilySampler.ThreeGenFamilySampler(families, numOffspring, pedSize, numOfAffected=0, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id')

      families
          number of families. This can be a number or a list of numbers. In the latter
          case, specified families are drawn from each subpopulation.
      
      numOffspring
          number of offspring. This can be a fixed number or a range [min, max].
      
      pedSize
          number of individuals in the Pedigree. This can be a fixed number or
          a range [min, max].
      
      numAfffected
          number of affected individuals in the Pedigree. This can be a fixed number
          or a range [min, max]
      
      subPops
          A list of (virtual) subpopulations from which samples are drawn.
          The default value is ALL_AVAIL, which means all available
          subpopulations of a population.

   .. method:: ThreeGenFamilySampler.drawSample(input_pop)

      Randomly select Pedigrees

   .. method:: ThreeGenFamilySampler.drawSamples(pop, numOfSamples)

      Draw multiple samples and return a list of populations.

   .. method:: ThreeGenFamilySampler.family(id)

      Return id, its spouse, their children, children's spouse and grandchildren

   .. method:: ThreeGenFamilySampler.prepareSample(input_pop)

      Prepare self.pedigree, some pedigree sampler might need additional loci and
      information fields for this sampler.



Function drawThreeGenFamilySample
---------------------------------


.. function:: drawThreeGenFamilySample(pop, families, numOffspring, pedSize, numOfAffected=0, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id')

   Draw three-generation families from a population. Such families consist
   of grant parents, their children, spouse of these children, and grand
   children. Number of offspring, total number of individuals, and total
   number of affected individuals in a pedigree should be specified using
   parameters ``numOffspring``, ``pedSize`` and ``numOfAffected``, which can all
   be a single number, or a range ``[a, b]`` (``b`` is incldued). If a single
   ``families`` is given, Pedigrees are drawn randomly from the whole
   Population or from specified (virtual) subpopulations (parameter
   ``subPops``). Otherwise, a list of numbers should be used to specify
   numbers of families from each subpopulation, which can be all
   subpopulations if ``subPops=ALL_AVAIL`` (default), or from each of the
   specified (virtual) subpopulations. This function returns a population that
   contains extracted individuals.


Function drawThreeGenFamilySamples
----------------------------------


.. function:: drawThreeGenFamilySamples(pop, families, numOffspring, pedSize, numOfAffected=0, numOfSamples=1, subPops=ALL_AVAIL, idField='ind_id', fatherField='father_id', motherField='mother_id')

   Draw ``numOfSamples`` three-generation pedigree samples from population ``pop``
   and return a list of populations. Please refer to function
   ``drawThreeGenFamilySample`` for a description of other parameters.


class CombinedSampler
---------------------

.. class:: CombinedSampler

   A combined sampler accepts a list of sampler objects, draw samples and
   combine the returned sample into a single population. An id field is
   required to use this sampler, which will be used to remove extra copies of
   individuals who have been drawn by different samplers.

   .. method:: CombinedSampler.CombinedSampler(samplers=[], idField='ind_id')

      samplers
          A list of samplers

   .. method:: CombinedSampler.drawSample(pop)

      Draw and return a sample.

   .. method:: CombinedSampler.drawSamples(pop, numOfSamples)

      Draw multiple samples and return a list of populations.

   .. method:: CombinedSampler.prepareSample(pop, rearrange)

      Prepare passed population object for sampling according to parameter
      ``subPops``. If samples are drawn from the whole population, a
      Population will be trimmed if only selected (virtual) subpopulations
      are used. If samples are drawn separately from specified subpopulations,
      Population ``pop`` will be rearranged (if ``rearrange==True``) so that
      each subpoulation corresponds to one element in parameter ``subPops``.



Function drawCombinedSample
---------------------------


.. function:: drawCombinedSample(pop, samplers, idField='ind_id')

   Draw different types of samples using a list of ``samplers``. A
   Population consists of all individuals from these samples will
   be returned. An ``idField`` that stores an unique ID for all individuals
   is needed to remove duplicated individuals who are drawn multiple
   numOfSamples from these samplers.


Function drawCombinedSamples
----------------------------


.. function:: drawCombinedSamples(pop, samplers, numOfSamples=1, idField='ind_id')

   Draw combined samples ``numOfSamples`` numOfSamples and return a list of populations.
   Please refer to function ``drawCombinedSample`` for details about
   parameters ``samplers`` and ``idField``.


