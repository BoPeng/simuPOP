Module :mod:`simuPOP.demography`
================================


.. module:: simuPOP.demography

This module provides some commonly used demographic models. In addition
to several migration rate generation functions, it provides models that
encapsulate complete demographic features of one or more populations (
population growth, split, bottleneck, admixture, migration). These models
provides:

1. The model itself can be passed to parameter subPopSize of a mating
   scheme to determine the size of the next generation. More importantly,
   it performs necessary actions of population size change when needed.

2. The model provides attribute num_gens, which can be passed to parameter
   ``gens`` of ``Simulator.evolve`` or ``Population.evolve`` function.
   A demographic model can also terminate an evolutionary process by
   returnning an empty list so ``gens=model.num_gens`` is no longer required.


Function migrIslandRates
------------------------


.. function:: migrIslandRates(r, n)

   migration rate matrix
   
   ::
   
            x m/(n-1) m/(n-1) ....
            m/(n-1) x ............
            .....
            .... m/(n-1) m/(n-1) x
   
       where x = 1-m


Function migrHierarchicalIslandRates
------------------------------------


.. function:: migrHierarchicalIslandRates(r1, r2, n)

   Return the migration rate matrix for a hierarchical island model
   where there are different migration rate within and across groups
   of islands.
   
   r1
       Within group migration rates. It can be a number or a list of numbers
       for each group of the islands.
   
   r2
       Across group migration rates which is the probability that someone will
       migrate to a subpopulation outside of his group. A list of r2 could be
       specified for each group of the islands.
   
   n
       Number of islands in each group. E.g. n=[5, 4] specifies two groups of
       islands with 5 and 4 islands each.
   
   For individuals in an island, the probability that it remains in the same
   island is 1-r1-r2 (r1, r2 might vary by island groups), that it migrates
   to another island in the same group is r1 and migrates to another island
   outside of the group is r2. migrate rate to a specific island depends on
   the size of group.


Function migrSteppingStoneRates
-------------------------------


.. function:: migrSteppingStoneRates(r, n, circular=False)

   migration rate matrix for circular stepping stone model (X=1-m)
   
   ::
   
              X   m/2               m/2
              m/2 X   m/2           0
              0   m/2 x   m/2 ......0
              ...
              m/2 0 ....       m/2  X
   
   or non-circular
   
   ::
   
              X   m/2               m/2
              m/2 X   m/2           0
              0   m/2 X   m/2 ......0
              ...
              ...              m   X
   
       This function returns [[1]] when there is only one subpopulation.


Function migrtwoDSteppingStoneRates
-----------------------------------


.. function:: migr2DSteppingStoneRates(r, m, n, diagonal=False, circular=False)

   migration rate matrix for 2D stepping stone model, with or without
   diagonal neighbors (4 or 8 neighbors for central patches). The boundaries
   are connected if circular is True. Otherwise individuals from corner and
   bounary patches will migrate to their neighbors with higher probability.


class EventBasedModel
---------------------

.. class:: EventBasedModel

   An event based demographic model in which the demographic changes are 
   triggered by demographic events such as population growth, split, join, and 
   admixture. The population size will be kept constant if no event is applied
   at a certain generation.

   .. method:: EventBasedModel.EventBasedModel(events=[], T=None, N0=[], ops=[], infoFields=[])

      A demographic model that is driven by a list of demographic events.
      The events should be subclasses of ``DemographicEvent``, which have similar
      interface as regular operators with the exception that applicable parameters
      ``begin``, ``end``, ``step``, ``at`` are relative to the demographic model,
      not the population.

   .. method:: EventBasedModel.plot(filename='', title='', initSize=[])

      Evolve a haploid population using a :class:`RandomSelection` mating scheme
      using the demographic model. Print population size changes duringe evolution.
      An initial population size could be specified using parameter ``initSize``
      for a demographic model with dynamic initial population size. If a filename
      is specified and if matplotlib is available, this function draws a figure
      to depict the demographic model and save it to ``filename``. An optional
      ``title`` could be specified to the figure. Note that this function can
      not be plot demographic models that works for particular mating schemes
      (e.g. genotype dependent).



class DemographicEvent
----------------------

.. class:: DemographicEvent

   A demographic events that will be applied to one or more populations at
   specified generations. The interface of a DemographicEvent is very similar to
   an simuPOP operator, but the applicable parameters are handled so that
   the generations are relative to the demographic model, not the populations
   to which the event is applied.

   .. method:: DemographicEvent.DemographicEvent(ops=[], output='', begin=0, end=-1, step=1, at=[], reps=True, subPops=ALL_AVAIL, infoFields=[])

      Create a demographic event that will be applied at specified
      generations according to applicability parameters ``reps``, ``begin``, 
      ``end``, ``step`` and ``at``. Parameter ``subPops`` is usually used
      to specify the subpopulations affected by the event. One or more simuPOP
      operators, if specified in ``ops``, will be applied when the event
      happens. Parameters ``output`` and ``infoFields`` are currently ignored.

   .. method:: DemographicEvent.apply(pop)




class ExpansionEvent
--------------------

.. class:: ExpansionEvent

   A demographic event that increase applicable population size by
   ``N*r`` (to size ``N*(1+r)``), or ``s`` (to size ``N+s``) at each applicable
   generation. The first model leads to an exponential population expansion
   model with rate ``r`` (``N(t)=N(0)*exp(r*t)``), where the second model leads to
   an linear population growth model (``N(t)=N(0)+s*t``) and this is why the
   parameter is called ``slopes``. Note that if both population
   size and ``r`` are small (e.g. ``N*r<1``),  the population might not expand
   as expected.

   .. method:: ExpansionEvent.ExpansionEvent(rates=[], slopes=[], capacity=[], name='', ops=[], output='', begin=0, end=-1, step=1, at=[], reps=True, subPops=ALL_AVAIL, infoFields=[])

      A demographic event that expands all or specified subpopulations
      (``subPops``) exponentially by a rate of ``rates``, or linearly by a slope
      of ``slopes``, unless carray capacity (``capacity``) of the population has
      been reached. Parameter ``rates`` can be a single number or a list of rates
      for all subpopulations. Parameter ``slopes`` should be a number, or a list
      of numbers for all subpopulations. ``subPops`` can be a ``ALL_AVAIL`` or a list
      of subpopulation index or names. ``capacity`` can be empty (no limit on
      carrying capacity), or one or more numbers for each of the subpopulations.

   .. method:: ExpansionEvent.apply(pop)




class ResizeEvent
-----------------

.. class:: ResizeEvent

   A demographic event that resize specified subpopulations

   .. method:: ResizeEvent.ResizeEvent(sizes=[], names=[], removeEmptySubPops=False, ops=[], output='', begin=0, end=-1, step=1, at=[], reps=True, subPops=ALL_AVAIL, infoFields=[])

      A demographic event that resizes given subpopulations ``subPops`` to new
      ``sizes`` (integer type), or sizes proportional to original sizes (if a float
      number is given). For example, ``sizes=[0.5, 500]`` will resize the first
      subpopulation to half of its original size, and the second subpopulation to
      size ``500``. If the new size is larger, existing individuals will be copied
      to sequentially, and repeatedly if needed. If the size of a subpopulation is
      0 and ``removeEmptySubPops`` is ``True``, empty subpopulations will be
      removed. A new set of names could be assigned to the population being resized.

   .. method:: ResizeEvent.apply(pop)




class SplitEvent
----------------

.. class:: SplitEvent

   A demographic event that splits a specified population into two or
   more subpopulations.

   .. method:: SplitEvent.SplitEvent(sizes=[], names=[], ops=[], output='', begin=0, end=-1, step=1, at=[], reps=True, subPops=ALL_AVAIL, infoFields=[])

      A demographic event that splits a subpopulation specified by
      ``subPops`` to two or more subpopulations, with specified ``sizes``
      and ``names``. ``sizes`` can be a list of numbers, proportions
      (e.g. ``[1., 500]`` keeps the original population and copies 500
      individuals to create a new subpupulation). Note that ``sizes``
      and ``names``, if specified, should include the source subpopulation
      as its first element.

   .. method:: SplitEvent.apply(pop)




class MergeEvent
----------------

.. class:: MergeEvent

   A demographic event that merges one or more subpopulation to
   a single one.

   .. method:: MergeEvent.MergeEvent(name='', ops=[], output='', begin=0, end=-1, step=1, at=[], reps=True, subPops=ALL_AVAIL, infoFields=[])

      A demographic event that merges subpopulations into a single subpopulation.
      The merged subpopulation will have the name of the first merged subpopulation
      unless a separate ``name`` is supported.

   .. method:: MergeEvent.apply(pop)




class AdmixtureEvent
--------------------

.. class:: AdmixtureEvent

   This event implements a population admixture event that mix
   individuals from specified subpopulations to either a new 
   subpopulation or an existing subpopulation.

   .. method:: AdmixtureEvent.AdmixtureEvent(sizes=[], toSubPop=None, name='', ops=[], output='', begin=0, end=-1, step=1, at=[], reps=True, subPops=ALL_AVAIL, infoFields=[])

      Create an admixed population by choosing individuals
      from all or specified subpopulations (``subPops``) and creating
      an admixed population ``toSubPop``. The admixed population will
      be appended to the population as a new subpopulation with name
      ``name`` if ``toSubPop`` is ``None`` (default), or replace an
      existing subpopulation with name or index ``toSubPop``. The admixed
      population consists of individuals from ``subPops`` according to
      specified ``sizes``. Its size is maximized to have the largest
      number of individuals from the source population when a new population
      is created, or equal to the size of the existing destination population.
      The parameter ``sizes`` should be a list of float numbers 
      between 0 and 1, and add up to 1 (e.g. ``[0.4, 0.4, 0.2]``, although
      this function ignores the last element and set it to 1 minus the 
      sum of the other numbers). Alternatively, parameter ``sizes`` can
      be a list of numbers used to explicitly specify the size of admixed
      population and number of individuals from each source subpopulation.
      In all cases, the size of source populations will be kept constant.

   .. method:: AdmixtureEvent.apply(pop)




class InstantChangeModel
------------------------

.. class:: InstantChangeModel

   A model for instant population change (growth, resize, merge, split).

   .. method:: InstantChangeModel.InstantChangeModel(T=None, N0=[], G=[], NG=[], ops=[], infoFields=[], removeEmptySubPops=False)

      An instant population growth model that evolves a population
      from size ``N0`` to ``NT`` for ``T`` generations with population
      size changes at generation ``G`` to ``NT``. If ``G`` is a list,
      multiple population size changes are allowed. In that case, a list
      (or a nested list) of population size should be provided to parameter
      ``NT``. Both ``N0`` and ``NT`` supports fixed (an integer), dynamic
      (keep passed poulation size) and proportional (an float number) population
      size. Optionally, one or more operators (e.g. a migrator) ``ops``
      can be applied to population. Required information fields by these
      operators should be passed to parameter ``infoFields``. If ``removeEmpty``
      option is set to ``True``, empty subpopulation will be removed. This
      option can be used to remove subpopulations.

   .. method:: InstantChangeModel.plot(filename='', title='', initSize=[])

      Evolve a haploid population using a :class:`RandomSelection` mating scheme
      using the demographic model. Print population size changes duringe evolution.
      An initial population size could be specified using parameter ``initSize``
      for a demographic model with dynamic initial population size. If a filename
      is specified and if matplotlib is available, this function draws a figure
      to depict the demographic model and save it to ``filename``. An optional
      ``title`` could be specified to the figure. Note that this function can
      not be plot demographic models that works for particular mating schemes
      (e.g. genotype dependent).



class ExponentialGrowthModel
----------------------------

.. class:: ExponentialGrowthModel

   A model for exponential population growth with carry capacity

   .. method:: ExponentialGrowthModel.ExponentialGrowthModel(T=None, N0=[], NT=None, r=None, ops=[], infoFields=[])

      An exponential population growth model that evolves a population from size
      ``N0`` to ``NT`` for ``T`` generations with ``r*N(t)`` individuals added
      at each generation. ``N0``, ``NT`` and ``r`` can be a list of population
      sizes or growth rates for multiple subpopulations. The initial population
      will be resized to ``N0`` (split if necessary). Zero or negative growth
      rates are allowed. The model will automatically determine ``T``, ``r``
      or ``NT`` if one of them is unspecified. If all of them are specified,
      ``NT`` is intepretted as carrying capacity of the model, namely the 
      population will keep contant after it reaches size ``NT``. Optionally,
      one or more operators (e.g. a migrator) ``ops`` can be applied to 
      population.

   .. method:: ExponentialGrowthModel.plot(filename='', title='', initSize=[])

      Evolve a haploid population using a :class:`RandomSelection` mating scheme
      using the demographic model. Print population size changes duringe evolution.
      An initial population size could be specified using parameter ``initSize``
      for a demographic model with dynamic initial population size. If a filename
      is specified and if matplotlib is available, this function draws a figure
      to depict the demographic model and save it to ``filename``. An optional
      ``title`` could be specified to the figure. Note that this function can
      not be plot demographic models that works for particular mating schemes
      (e.g. genotype dependent).



class LinearGrowthModel
-----------------------

.. class:: LinearGrowthModel

   A model for linear population growth with carry capacity.

   .. method:: LinearGrowthModel.LinearGrowthModel(T=None, N0=[], NT=None, r=None, ops=[], infoFields=[])

      An linear population growth model that evolves a population from size
      ``N0`` to ``NT`` for ``T`` generations with ``r*N0`` individuals added
      at each generation. ``N0``, ``NT`` and ``r`` can be a list of population
      sizes or growth rates for multiple subpopulations. The initial population
      will be resized to ``N0`` (split if necessary). Zero or negative growth
      rates are allowed. The model will automatically determine ``T``, ``r``
      or ``NT`` if one of them is unspecified. If all of them are specified,
      ``NT`` is intepretted as carrying capacity of the model, namely the 
      population will keep contant after it reaches size ``NT``. Optionally,
      one or more operators (e.g. a migrator) ``ops`` can be applied to 
      population.

   .. method:: LinearGrowthModel.plot(filename='', title='', initSize=[])

      Evolve a haploid population using a :class:`RandomSelection` mating scheme
      using the demographic model. Print population size changes duringe evolution.
      An initial population size could be specified using parameter ``initSize``
      for a demographic model with dynamic initial population size. If a filename
      is specified and if matplotlib is available, this function draws a figure
      to depict the demographic model and save it to ``filename``. An optional
      ``title`` could be specified to the figure. Note that this function can
      not be plot demographic models that works for particular mating schemes
      (e.g. genotype dependent).



class MultiStageModel
---------------------

.. class:: MultiStageModel

   A multi-stage demographic model that connects a number of demographic
   models.

   .. method:: MultiStageModel.MultiStageModel(models, ops=[], infoFields=[])

      An multi-stage demographic model that connects specified
      demographic models ``models``. It applies a model to the population
      until it reaches ``num_gens`` of the model, or if the model returns
      ``[]``. One or more operators could be specified, which will be applied
      before a demographic model is applied. Note that the last model will be
      ignored if it lasts 0 generation.

   .. method:: MultiStageModel.plot(filename='', title='', initSize=[])

      Evolve a haploid population using a :class:`RandomSelection` mating scheme
      using the demographic model. Print population size changes duringe evolution.
      An initial population size could be specified using parameter ``initSize``
      for a demographic model with dynamic initial population size. If a filename
      is specified and if matplotlib is available, this function draws a figure
      to depict the demographic model and save it to ``filename``. An optional
      ``title`` could be specified to the figure. Note that this function can
      not be plot demographic models that works for particular mating schemes
      (e.g. genotype dependent).



class OutOfAfricaModel
----------------------

.. class:: OutOfAfricaModel

   A dempgraphic model for the CHB, CEU, and YRI populations, as defined in
   Gutenkunst 2009, Plos Genetics. The model is depicted in Figure 2, and the 
   default parameters are listed in Table 1 of this paper.

   .. method:: OutOfAfricaModel.OutOfAfricaModel(T0, N_A=7300, N_AF=12300, N_B=2100, N_EU0=1000, r_EU=0.004, N_AS0=510, r_AS=0.0055, m_AF_B=0.00025, m_AF_EU=3e-05, m_AF_AS=1.9e-05, m_EU_AS=9.6e-05, T_AF=8800, T_B=5600, T_EU_AS=848, ops=[], infoFields=[], outcome=['AF', 'EU', 'AS'], scale=1)

      Counting **backward in time**, this model evolves a population for ``T0``
      generations (required parameter). The ancient population ``A`` started at
      size ``N_A`` and expanded at ``T_AF`` generations from now, to pop ``AF``
      with size ``N_AF``. Pop ``B`` split from pop ``AF`` at ``T_B`` generations
      from now, with size ``N_B``; Pop ``AF`` remains as ``N_AF`` individuals. 
      Pop ``EU`` and  ``AS`` split from pop ``B`` at ``T_EU_AS`` generations
      from now; with size ``N_EU0`` individuals and ``N_ASO`` individuals,
      respectively. Pop ``EU`` grew exponentially with rate ``r_EU``; Pop
      ``AS`` grew exponentially with rate ``r_AS``. The ``YRI``, ``CEU`` and
      ``CHB`` samples are drawn from ``AF``, ``EU`` and ``AS`` populations
      respectively. Additional operators could be added to ``ops``. Information
      fields required by these operators should be passed to ``infoFields``. If 
      a scaling factor ``scale`` is specified, all population sizes and
      generation numbers will be divided by a factor of ``scale``. This demographic
      model by default returns all populations (``AF``, ``EU``, ``AS``) but
      you can choose to keep only selected subpopulations using parameter
      ``outcome`` (e.g. ``outcome=['EU', 'AS']``).
      
      This model merges all subpopulations if it is applied to an initial 
      population with multiple subpopulation.

   .. method:: OutOfAfricaModel.plot(filename='', title='', initSize=[])

      Evolve a haploid population using a :class:`RandomSelection` mating scheme
      using the demographic model. Print population size changes duringe evolution.
      An initial population size could be specified using parameter ``initSize``
      for a demographic model with dynamic initial population size. If a filename
      is specified and if matplotlib is available, this function draws a figure
      to depict the demographic model and save it to ``filename``. An optional
      ``title`` could be specified to the figure. Note that this function can
      not be plot demographic models that works for particular mating schemes
      (e.g. genotype dependent).



class SettlementOfNewWorldModel
-------------------------------

.. class:: SettlementOfNewWorldModel

   A dempgraphic model for settlement of the new world of Americans, as defined
   in Gutenkunst 2009, Plos Genetics. The model is depicted in Figure 3, and the 
   default parameters are listed in Table 2 of this paper.

   .. method:: SettlementOfNewWorldModel.SettlementOfNewWorldModel(T0, N_A=7300, N_AF=12300, N_B=2100, N_EU0=1500, r_EU=0.0023, N_AS0=590, r_AS=0.0037, N_MX0=800, r_MX=0.005, m_AF_B=0.00025, m_AF_EU=3e-05, m_AF_AS=1.9e-05, m_EU_AS=1.35e-05, T_AF=8800, T_B=5600, T_EU_AS=1056, T_MX=864, f_MX=0.48, ops=[], infoFields=[], outcome='MXL', scale=1)

      Counting **backward in time**, this model evolves a population for ``T0``
      generations. The ancient population ``A`` started at size ``N_A`` and
      expanded at ``T_AF`` generations from now, to pop ``AF`` with size ``N_AF``.
      Pop ``B`` split from pop ``AF`` at ``T_B`` generations from now, with
      size ``N_B``; Pop ``AF`` remains as ``N_AF`` individuals. Pop ``EU`` and 
      ``AS`` split from pop ``B`` at ``T_EU_AS`` generations from now; with 
      size ``N_EU0`` individuals and ``N_ASO`` individuals, respectively. Pop
      ``EU`` grew exponentially with final population size ``N_EU``; Pop
      ``AS`` grew exponentially with final populaiton size ``N_AS``. Pop ``MX``
      split from pop ``AS`` at ``T_MX`` generations from now with size ``N_MX0``,
      grew exponentially to final size ``N_MX``. Migrations are allowed between
      populations with migration rates ``m_AF_B``, ``m_EU_AS``, ``m_AF_EU``,
      and ``m_AF_AS``. At the end of the evolution, the ``AF`` and ``CHB``
      populations are removed, and the ``EU`` and ``MX`` populations are merged
      with ``f_MX`` proportion for ``MX``. The Mexican American<F19> sample could
      be sampled from the last single population. Additional operators could
      be added to ``ops``. Information fields required by these operators 
      should be passed to ``infoFields``. If a scaling factor ``scale``
      is specified, all population sizes and generation numbers will be divided by
      a factor of ``scale``. This demographic model by default only returns the
      mixed Mexican America model (``outputcom='MXL'``) but you can specify any
      combination of ``AF``, ``EU``, ``AS``, ``MX`` and ``MXL``.
      
      This model merges all subpopulations if it is applied to an initial population
      with multiple subpopulation.

   .. method:: SettlementOfNewWorldModel.plot(filename='', title='', initSize=[])

      Evolve a haploid population using a :class:`RandomSelection` mating scheme
      using the demographic model. Print population size changes duringe evolution.
      An initial population size could be specified using parameter ``initSize``
      for a demographic model with dynamic initial population size. If a filename
      is specified and if matplotlib is available, this function draws a figure
      to depict the demographic model and save it to ``filename``. An optional
      ``title`` could be specified to the figure. Note that this function can
      not be plot demographic models that works for particular mating schemes
      (e.g. genotype dependent).



class CosiModel
---------------

.. class:: CosiModel

   A dempgraphic model for Africa, Asia and Europe, as described in 
   Schaffner et al, Genome Research, 2005, and implemented in the coalescent
   simulator cosi.

   .. method:: CosiModel.CosiModel(T0, N_A=12500, N_AF=24000, N_OoA=7700, N_AF1=100000, N_AS1=100000, N_EU1=100000, T_AF=17000, T_OoA=3500, T_EU_AS=2000, T_AS_exp=400, T_EU_exp=350, T_AF_exp=200, F_OoA=0.085, F_AS=0.067, F_EU=0.02, F_AF=0.02, m_AF_EU=3.2e-05, m_AF_AS=8e-06, ops=[], infoFields=[], scale=1)

      Counting **backward in time**, this model evolves a population for a
      total of ``T0`` generations. The ancient population ``Ancestral`` started
      at size ``N_Ancestral`` and expanded at ``T_AF`` generations from now,
      to pop ``AF`` with size ``N_AF``. The Out of Africa population split from
      the ``AF`` population at ``T_OoA`` generations ago. The ``OoA`` population
      split into two subpopulations ``AS`` and ``EU`` but keep the same size.
      At the generations of ``T_EU_exp``, ``T_AS_exp``, and ``T_AF_exp`` ago,
      three populations expanded to modern population sizes of ``N_AF1``, 
      ``N_AS1`` and ``N_EU1`` exponentially, respectively. Migrations are
      allowed between ``AF`` and ``EU`` populations
      with rate ``m_AF_EU``, and between ``AF`` and ``AS`` with rate ``m_AF_AS``.
      
      Four bottlenecks happens in the ``AF``, ``OoA``, ``EU`` and ``AS`` populations.
      They are supposed to happen 200 generations after population split and last
      for 200 generations. The intensity is parameterized in F, which is number
      of generations devided by twice the effective size during bottleneck.
      So the bottleneck size is 100/F. 
      
      This model merges all subpopulations if it is applied to a population with
      multiple subpopulation. Although parameters are configurable, we assume
      the order of events so dramatically changes of parameters might need
      to errors.  If a scaling factor ``scale`` is specified, all population
      sizes and generation numbers will be divided by, and migration rates
      will be multiplied by a factor of ``scale``.

   .. method:: CosiModel.plot(filename='', title='', initSize=[])

      Evolve a haploid population using a :class:`RandomSelection` mating scheme
      using the demographic model. Print population size changes duringe evolution.
      An initial population size could be specified using parameter ``initSize``
      for a demographic model with dynamic initial population size. If a filename
      is specified and if matplotlib is available, this function draws a figure
      to depict the demographic model and save it to ``filename``. An optional
      ``title`` could be specified to the figure. Note that this function can
      not be plot demographic models that works for particular mating schemes
      (e.g. genotype dependent).



