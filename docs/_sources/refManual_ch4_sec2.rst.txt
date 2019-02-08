Module :mod:`simuPOP.utils`
===========================


.. module:: simuPOP.utils

This module provides some commonly used operators
and format conversion utilities.


class Trajectory
----------------

.. class:: Trajectory

   A ``Trajectory`` object contains frequencies of one or more loci in one
   or more subpopulations over several generations. It is usually returned by
   member functions of class ``TrajectorySimulator`` or equivalent global
   functions ``simulateForwardTrajectory`` and ``simulateBackwardTrajectory``.
   
   The ``Trajectory`` object provides several member functions to facilitate
   the use of Trajectory-simulation techiniques. For example,
   ``Trajectory.func()`` returns a trajectory function that can be provided
   directly to a :class:`ControlledOffspringGenerator`; ``Trajectory.mutators()``
   provides a list of :class:`PointMutator` that insert mutants at the right
   generations to initialize a trajectory.
   
   For more information about Trajectory simulation techniques and related
   controlled random mating scheme, please refer to the simuPOP user's guide,
   and Peng et al (PLoS Genetics 3(3), 2007).

   .. method:: Trajectory.Trajectory(endGen, nLoci)

      Create a ``Trajectory`` object of alleles at *nLoci* loci with
      ending generation *endGen*. *endGen* is the generation when expected
      allele frequencies are reached after mating. Therefore, a trajectory
      for 1000 generations should have ``endGen=999``.

   .. method:: Trajectory.freq(gen, subPop)

      Return frequencies of all loci in subpopulation *subPop* at
      generation *gen* of the simulated Trajectory. Allele frequencies are
      assumed to be zero if *gen* is out of range of the simulated
      Trajectory.

   .. method:: Trajectory.func()

      Return a Python function that returns allele frequencies for each
      locus at specified loci. If there are multiple subpopulations, allele
      frequencies are arranged in the order of ``loc0_sp0``, ``loc1_sp0``,
      ..., ``loc0_sp1``, ``loc1_sp1``, ... and so on. The returned function
      can be supplied directly to the ``freqFunc`` parameter of a controlled
      random mating scheme (:class:`ControlledRandomMating`) or a homogeneous
      mating scheme that uses a controlled offspring generator
      (:class:`ControlledOffspringGenerator`).

   .. method:: Trajectory.mutants()

      Return a list of mutants in the form of (loc, gen, subPop)

   .. method:: Trajectory.mutators(loci, inds=0, allele=1, *args, **kwargs)

      Return a list of :class:`PointMutator` operators that introduce mutants
      at the beginning of simulated trajectories. These mutators should be
      added to the ``preOps`` parameter of :meth:`Simulator.evolve` function to
      introduce a mutant at the beginning of a generation with zero allele
      frequency before mating, and a positive allele frequency after mating.
      A parameter ``loci`` is needed to specify actual loci indexes in the
      real forward simulation. Other than default parameters ``inds=0`` and
      ``allele=1``, additional parameters could be passed to point mutator
      as keyward parameters.



class TrajectorySimulator
-------------------------

.. class:: TrajectorySimulator

   A Trajectory Simulator takes basic demographic and genetic (natural
   selection) information of an evolutionary process of a diploid population
   and allow the simulation of Trajectory of allele frequencies of one or
   more loci. Trajectories could be simulated in two ways: forward-time and
   backward-time. In a forward-time simulation, the simulation starts from
   certain allele frequency and simulate the frequency at the next generation
   using given demographic and genetic information. The simulation continues
   until an ending generation is reached. A Trajectory is successfully
   simulated if the allele frequency at the ending generation falls into a
   specified range. In a backward-time simulation, the simulation starts from
   the ending generation with a desired allele frequency and simulate the
   allele frequency at previous generations one by one until the allele gets
   lost (allele frequency equals zero).
   
   The result of a trajectory simulation is a trajectory object which can be
   used to direct the simulation of a special random mating process that
   controls the evolution of one or more disease alleles so that allele 
   frequencies are consistent across replicate simulations. For more
   information about Trajectory simulation techniques and related controlled
   random mating scheme, please refer to the simuPOP user's guide, and Peng et
   al (PLoS Genetics 3(3), 2007).

   .. method:: TrajectorySimulator.TrajectorySimulator(N, nLoci=1, fitness=None, logger=None)

      Create a trajectory Simulator using provided demographic and genetic
      (natural selection) parameters. Member functions *simuForward* and
      *simuBackward* can then be used to simulate trajectories within certain
      range of generations. This class accepts the following parameters
      
      N
          Parameter *N* accepts either a constant number for population size
          (e.g. N=1000), a list of subpopulation sizes (e.g. N=[1000, 2000]),
          or a demographic function that returns population or subpopulation
          sizes at each generation. During the evolution, multiple
          subpopulations can be merged into one, and one population can be
          split into several subpopulations. The number of subpopulation is
          determined by the return value of the demographic function. Note
          that *N* should be considered as the population size at the end of
          specified generation.
      
      nLoci
          Number of unlinked loci for which trajectories of allele
          frequencies are simulated. We assume a diploid population with
          diallelic loci. The Trajectory represents frequencies of a 
      
      fitness
          Parameter fitness can be ``None`` (no selection), a list of fitness
          values for genotype with 0, 1, and 2 disease alleles (*AA*, *Aa*,
          and *aa*) at one or more loci; or a function that returns fitness
          values at each generation. When multiple loci are involved
          (*nLoci*), *fitness* can be a list of 3 (the same fitness values
          for all loci), a list of 3*nLoci (different fitness values for each
          locus) or a list of 3**nLoci (fitness value for each combination of
          genotype). The fitness function should accept generation number and
          a subpopulation index. The latter parameter allows, and is the only
          way to specify different fitness in each subpopulation.
      
      logger
          A logging object (see Python module ``logging``) that can be used
          to output intermediate results with debug information.

   .. method:: TrajectorySimulator.simuBackward(endGen, endFreq, minMutAge=None, maxMutAge=None, maxAttempts=1000)

      Simulate trajectories of multiple disease susceptibility loci using
      a forward time approach. This function accepts allele frequencies of
      alleles of multiple unlinked loci (*endFreq*) at the end of generation
      *endGen*. Depending on the number of loci and subpopulations, parameter
      *beginFreq* can be a number (same frequency for all loci in all
      subpopulations), or a list of frequencies for each locus (same
      frequency in all subpopulations), or a list of frequencies for each
      locus in each subpopulation in the order of ``loc0_sp0``, ``loc1_sp0``,
      ..., ``loc0_sp1``, ``loc1_sp1``, ... and so on.
      
      This simulator will simulate a trajectory generation by generation and
      restart if the disease allele got fixed (instead of lost), or if the 
      length simulated Trajectory does not fall into *minMutAge* and
      *maxMutAge* (ignored if ``None`` is given). This simulator will return
      ``None`` if no valid Trajectory is found after ``maxAttempts`` attemps.

   .. method:: TrajectorySimulator.simuForward(beginGen, endGen, beginFreq, endFreq, maxAttempts=10000)

      Simulate trajectories of multiple disease susceptibility loci using a
      forward time approach. This function accepts allele frequencies of
      alleles of multiple unlinked loci at the beginning generation (``freq``)
      at generation ``beginGen``, and expected *range* of allele frequencies
      of these alleles (``endFreq``) at the end of generation ``endGen``.
      Depending on the number of loci and subpopulations, these parameters
      accept the following inputs:
      
      beginGen
          Starting generation. The initial frequecies are considered as
          frequencies at the *beginning* of this generation.
      
      endGen
          Ending generation. The ending frequencies are considerd as
          frequencies at the *end* of this generation.
          
      beginFreq
          The initial allele frequency of involved loci in all subpopulations.
          It can be a number (same frequency for all loci in all
          subpopulations), or a list of frequencies for each locus (same
          frequency in all subpopulations), or a list of frequencies for each
          locus in each subpopulation in the order of ``loc0_sp0``,
          ``loc1_sp0``, ..., ``loc0_sp1``, ``loc1_sp1``, ... and so on.
      
      endFreq
          The range of acceptable allele frequencies at the ending generation.
          The ranges can be specified for all loci in all subpopulations,
          for all loci (allele frequency in the whole population is
          considered), or for all loci in all subpopulations, in the order
          of ``loc0_sp0``, ``loc1_sp0``, .... ``loc0_sp1``, ... and so on.
      
      This simulator will simulate a trajectory generation by generation and
      restart if the resulting frequencies do not fall into specified range
      of frequencies. This simulator will return ``None`` if no valid
      Trajectory is found after ``maxAttempts`` attemps.



Function simulateForwardTrajectory
----------------------------------


.. function:: simulateForwardTrajectory(N, beginGen, endGen, beginFreq, endFreq, nLoci=1, fitness=None, maxAttempts=10000, logger=None)

   Given a demographic model (*N*) and the fitness of genotype at one or
   more loci (*fitness*), this function simulates a trajectory of one or more
   unlinked loci (*nLoci*) from allele frequency *freq* at generation
   *beginGen* forward in time, until it reaches generation *endGen*. A
   ``Trajectory`` object will be returned if the allele frequency falls
   into specified ranges (*endFreq*). ``None`` will be returned if no valid
   Trajectory is simulated after ``maxAttempts`` attempts. Please refer to
   class ``Trajectory``, ``TrajectorySimulator`` and their member functions
   for more details about allowed input for these parameters. If a *logger*
   object is given, it will send detailed debug information at ``DEBUG``
   level and ending allele frequencies at the ``INFO`` level. The latter
   can be used to adjust your fitness model and/or ending allele frequency
   if a trajectory is difficult to obtain because of parameter mismatch.


Function simulateBackwardTrajectory
-----------------------------------


.. function:: simulateBackwardTrajectory(N, endGen, endFreq, nLoci=1, fitness=None, minMutAge=None, maxMutAge=None, maxAttempts=1000, logger=None)

   Given a demographic model (*N*) and the fitness of genotype at one or
   more loci (*fitness*), this function simulates a trajectory of one or more
   unlinked loci (*nLoci*) from allele frequency *freq* at generation *endGen*
   backward in time, until all alleles get lost. A ``Trajectory`` object will
   be returned if the length of simulated Trajectory with ``minMutAge`` and
   ``maxMutAge`` (if specified). ``None`` will be returned if no valid
   Trajectory is simulated after ``maxAttempts`` attempts. Please refer to
   class ``Trajectory``, ``TrajectorySimulator`` and their member functions
   for more details about allowed input for these parameters. If a *logger*
   object is given, it will send detailed debug information at ``DEBUG``
   level and ending generation and frequency at the ``INFO`` level. The latter
   can be used to adjust your fitness model and/or ending allele frequency
   if a trajectory is difficult to obtain because of parameter mismatch.


class ProgressBar
-----------------

.. class:: ProgressBar

   The ``ProgressBar`` class defines a progress bar. This class will use a
   text-based progress bar that outputs progressing dots (.) with intermediate
   numbers (e.g. 5 for 50%) under a non-GUI mode (``gui=False``) or not displaying
   any progress bar if ``gui='batch'``. In the GUI mode, a Tkinter or wxPython 
   progress dialog will be used (``gui=Tkinter``  or ``gui=wxPython``). The default
   mode is determined by the global gui mode of simuPOP
   (see also :func:`simuOpt.setOptions`).
   
   This class is usually used as follows::
   
       progress = ProgressBar("Start simulation", 500)
       for i in range(500):
           # i+1 can be ignored if the progress bar is updated by 1 step
           progress.update(i+1)   
       # if you would like to make sure the done message is displayed.
       progress.done()

   .. method:: ProgressBar.ProgressBar(message, totalCount, progressChar='.', block=2, done=' Done.\n', gui=None)

      Create a progress bar with ``message``, which will be the title of
      a progress dialog or a message for textbased progress bar. Parameter
      ``totalCount`` specifies total expected steps. If a text-based progress
      bar is used, you could specified progress character and intervals at
      which progresses will be displayed using parameters ``progressChar``
      and ``block``. A ending message will also be displayed in text mode.

   .. method:: ProgressBar.done()

      Finish progressbar, print 'done' message if in text-mode.

   .. method:: ProgressBar.update(count=None)

      Update the progreebar with ``count`` steps done. The dialog or textbar
      may not be updated if it is updated by full percent(s). If ``count`` is
      ``None``, the progressbar increases by one step (not percent).



Function viewVars
-----------------


.. function:: viewVars(var, gui=None)

   list a variable in tree format, either in text format or in a
       wxPython window.
   
   var
       A dictionary variable to be viewed. Dictionary wrapper objects returned
       by ``Population.dvars()`` and ``Simulator.dvars()`` are also acceptable.
   
   gui
       If gui is ``False`` or ``'Tkinter'``, a text presentation (use the
       pprint module) of the variable will be printed to the screen. If gui is
       ``'wxPython'`` and wxPython is available, a wxPython windows will be
       used. The default mode is determined by the global gui mode (see also
       :func:`simuOpt.setOptions`).


Function saveCSV
----------------


.. function:: saveCSV(pop, filename='', infoFields=[], loci=True, header=True, subPops=ALL_AVAIL, genoFormatter=None, infoFormatter=None, sexFormatter={1: 'M', 2: 'F'}, affectionFormatter={True: 'A', False: 'U'}, sep=', ', **kwargs)

   This function is deprecated. Please use ``export(format='csv')`` instead.
   Save a simuPOP population ``pop`` in csv format. Columns of this
   file is arranged in the order of information fields (``infoFields``),
   sex (if ``sexFormatter`` is not ``None``), affection status (if
   ``affectionFormatter`` is not ``None``), and genotype (if ``genoFormatter`` is
   not ``None``). This function only output individuals in the present
   generation of population ``pop``. This function accepts the following
   parameters:
   
   pop
       A simuPOP population object.
   
   filename
       Output filename. Leading '>' characters are ignored. However, if the first
       character of this filename is '!', the rest of the name will be evalulated
       in the population's local namespace. If ``filename`` is empty, the content
       will be written to the standard output.
   
   infoFileds
       Information fields to be outputted. Default to none.
   
   loci
       If a list of loci is given, only genotype at these loci will be
       written. Default to ``ALL_AVAIL``, meaning all available loci. You can
       set this parameter to ``[]`` if you do not want to output any genotype.
   
   header
       Whether or not a header should be written. These headers will include
       information fields, sex (if ``sexFormatter`` is not ``None``), affection
       status (if ``affectionFormatter`` is not ``None``) and loci names. If
       genotype at a locus needs more than one column, ``_1``, ``_2`` etc will
       be appended to loci names. Alternatively, a complete header (a string)
       or a list of column names could be specified directly.
   
   subPops
       A list of (virtual) subpopulations. If specified, only individuals
       from these subpopulations will be outputed.
   
   infoFormatter
       A format string that is used to format all information fields. If
       unspecified, ``str(value)`` will be used for each information field.
   
   genoFormatter
       How to output genotype at specified loci. Acceptable values include
       ``None`` (output allele names), a dictionary with genotype as keys,
       (e.g. ``genoFormatter={(0,0):1, (0,1):2, (1,0):2, (1,1):3}``, or a function
       with genotype (as a tuple of integers) as inputs. The dictionary value
       or the return value of this function can be a single or a list of
       number or strings.
   
   sexFormatter
       How to output individual sex. Acceptable values include ``None`` (no
       output) or a dictionary with keys ``MALE`` and ``FEMALE``.
   
   affectionFormatter
       How to output individual affection status. Acceptable values include
       ``None`` (no output) or a dictionary with keys ``True`` and ``False``.
   
   Parameters ``genoCode``, ``sexCode``, and ``affectionCode`` from version
   1.0.0 have been renamed to ``genoFormatter``, ``sexFormatter`` and 
   ``affectionFormatter`` but can still be used.


class Exporter
--------------

.. class:: Exporter

   An operator to export the current population in specified format.
   Currently supported file formats include:
   
   STRUCTURE (http://pritch.bsd.uchicago.edu/structure.html). This format
   accepts the following parameters:
   
   markerNames
       If set to True (default), output names of loci that are specified by parameter
       *lociNames* of the :class:`Population` class. No names will be outputted if loci are
       anonymous. A list of loci names are acceptable which will be outputted directly.
   
   recessiveAlleles
       If specified, value of this parameter will be outputted after the marker names
       line.
   
   interMarkerDistances
       If set to True (default), output distances between markers. The first marker
       of each chromosome has distance -1, as required by this format.
   
   phaseInformation
       If specified, output the value (0 or 1) of this parameter after the inter marker
       distances line. Note that simuPOP populations always have phase information.
   
   label
       Output 1-based indexes of individuals if this parameter is true (default)
   
   popData
       Output 1-based index of subpopulation if this parameter is set to true (default).
   
   popFlag
       Output value of this parameter (0 or 1) after popData if this parameter specified.
   
   locData
       Name of an information field with location information of each individual. Default
       to None (no location data)
   
   phenotype
       Name of an information field with phenotype information of each individual. Default
       to None (no phenotype)
   
       
   Genotype information are always outputted. Alleles are coded the same way (0, 1, 2, etc)
   as they are stored in simuPOP.
   
   GENEPOP (http://genepop.curtin.edu.au/). The genepop format accepts the following
   parameters:
   
   title
       The tile line. If unspecified, a line similar to 'produced by simuPOP on XXX'
       will be outputted.
   
   adjust
       Adjust values of alleles by specified value (1 as default). This adjustment is
       necessary in many cases because GENEPOP treats allele 0 as missing values, and 
       simuPOP treats allele 0 as a valid allele. Exporting alleles 0 and 1 as 1 and 2
       will allow GENEPOP to analyze simuPOP-exported files correctly.
   
   Because 0 is reserved as missing data in this format, allele A is outputted as A+adjust.
   simuPOP will use subpopulation names (if available) and 1-based individual index
   to output individual label (e.g. SubPop2-3). If parameter subPops is used to output
   selected individuals, each subpop will be outputted as a separate subpopulation even 
   if there are multiple virtual subpopulations from the same subpopulation. simuPOP 
   currently only export diploid populations to this format.
   
   FSTAT (http://www2.unil.ch/popgen/softwares/fstat.htm). The fstat format accepts
   the following parameters:
   
   lociNames
       Names of loci that will be outputted. If unspecified, simuPOP will try to use
       names of loci that are specified by parameter *lociNames* of the :class:`Population`
       class, or names in the form of chrX-Y.
   
   adjust
       Adjust values of alleles by specified value (1 as default). This adjustment is
       necessary in many cases because FSTAT treats allele 0 as missing values, and 
       simuPOP treats allele 0 as a valid allele. Exporting alleles 0 and 1 as 1 and 2
       will allow FSTAT to analyze simuPOP-exported files correctly.
       
   MAP (marker information format) output information about each loci. Each line of
   the map file describes a single marker and contains chromosome name, locus name,
   and position. Chromosome and loci names will be the names specified by parameters
   ``chromNames`` and ``lociNames`` of the :class:`Population` object, and will be
   chromosome index + 1, and '.' if these parameters are not specified. This
   format output loci position to the third column. If the unit assumed in your
   population does not match the intended unit in the MAP file, (e.g. you would like
   to output position in basepair while the population uses Mbp), you can use parameter
   ``posMultiplier`` to adjust it. This format accepts the following parameters:
   
   posMultiplier
       A number that will be multiplied to loci positions (default to 1). The result
       will be outputted in the third column of the output.
   
   
   PED (Linkage Pedigree pre MAKEPED format), with columns of family, individual,
   father mother, gender, affection status and genotypes. The output should be 
   acceptable by HaploView or plink, which provides more details of this format in
   their documentation. If a population does not have ``ind_id``, ``father_id`` or 
   ``mother_id``, this format will output individuals in specified (virtual) 
   subpopulations in the current generation (parental generations are ignored) 
   as unrelated individuals with 0, 0 as parent IDs. An incremental family
   ID will be assigned for each individual. If a population have ``ind_id``,
   ``father_id`` and ``mother_id``, parents will be recursively traced to separate
   all individuals in a (multigenerational) population into families of related
   individuals. father and mother id will be set to zero if one of them does not
   exist. This format uses 1 for MALE, 2 for FEMALE. If phenoField is ``None``,
   individual affection status will be outputted with 1 for Unaffected and 2
   for affected. Otherwise, values of an information field will be outputted as
   phenotype. Because 0 value indicates missing value, values of alleles will
   be adjusted by 1 by default, which should be avoided if you are using non-zero
   alleles to model ACTG alleles in simuPOP. This format will ignore subpopulation
   structure because parents might belong to different subpopulations. This format
   accepts the following parameters:
   
   idField
       A field for individual id, default to ``ind_id``. Value at this field will be
       individual ID inside a pedigree.
   
   fatherField
       A field for father id, default to ``father_id``. Value at this field will be
       used to output father of an individual, if an individual with this ID exists
       in the population.
   
   motherField
       A field for mother id, default to ``mother_id``. Value at this field will be
       used to output mother of an individual, if an individual with this ID exists
       in the population.
   
   phenoField
       A field for individual phenotype that will be outputted as the sixth column of
       the PED file. If ``None`` is specified (default), individual affection status
       will be outputted (1 for unaffected and 2 for affected).
   
   adjust
       Adjust values of alleles by specified value (1 as default). This adjustment
       is necessary in many cases because LINKAGE/PED format treats allele 0 as
       missing values, and simuPOP treats allele 0 as a valid allele. You should set
       this paremter to zero if you have already used alleles 1, 2, 3, 4 to model 
       A, C, T, and G alleles.
       
   Phylip (Joseph Felsenstein's Phylip format). Phylip is generally used for nuclotide
   sequences and protein sequences. This makes this format suitable for simulations
   of haploid populations (ploidy=1) with nucleotide or protein sequences (number of
   alleles = 4 or 24 with alleleNames as nucleotide or amino acid names). If your
   population does satisfy these conditions, you can still export it, with homologous
   chromosomes in a diploid population as two sequences, and with specified allele
   names for allele 0, 1, 2, .... This function outputs sequence name as SXXX where
   XXX is the 1-based index of individual and SXXX_Y (Y=1 or 2) for diploid individuals,
   unless names of sequences are provided by parameter seqNames. This format supports
   the following parameters:
   
   alleleNames
       Names of alleles 0, 1, 2, ... as a single string (e.g. 'ACTG') or a list of 
       single-character strings (e.g. ['A', 'C', 'T', 'G']). If this parameter is
       unspecified (default), this program will try to use names of alleles
       specified in alleleNames parameter of a Population, and raise an error if no
       name could be found.
   
   seqNames
       Names of each sequence outputted, for each individual, or for each sequences
       for non-haploid population. If unspecified, default names such as SXXX or
       SXXX_Y will be used.
   
   style
       Output style, can be 'sequential' (default) or 'interleaved'. For sequential
       output, each sequence consists of for the first line a name and 90 symbols
       starting from column 11, and subsequent lines of 100 symbols. The interleaved
       style have subsequent lines as separate blocks.
       
   MS (output from Richard R. Hudson's MS or msHOT program). This format records
   genotypes of SNP markers at segregating site so all non-zero genotypes are
   recorded as 1. simuPOP by default outputs a single block of genotypes at
   all loci on the first chromosome, and for all individuals, unless parameter
   ``splitBy`` is specified to separate genotypes by chromosome or subpopulations.
   
   splitBy:
       simuPOP by default output segregating sites at all loci on the first 
       chromosome for all individuals. If ``splitBy`` is set to ``'subPop'``,
       genotypes for individuals in all or specified (parameter ``subPops``) 
       subpopulations are outputted in separate blocks. The subpopulations should
       have the same number of individuals to produce blocks of the same number
       of sequences. Alternatively, ``splitBy`` can be set to ``chrom``, for
       which genotypes on different chromosomes will be outputted separately.
   
   
   CSV (comma separated values). This is a general format that output genotypes in
   comma (or tab etc) separated formats. The function form of this operator 
   ``export(format='csv')`` is similar to the now-deprecated ``saveCSV`` function,
   but its interface has been adjusted to match other formats supported by this
   operator. This format outputs a header (optiona), and one line for each individual
   with values of specified information fields, sex, affection status, and genotypes.
   All fields except for genotypes are optional. The output format is controlled by the
   following parameters:
   
   infoFileds
       Information fields to be outputted. Default to none.
   
   header
       Whether or not a header should be written. These headers will include
       information fields, sex (if ``sexFormatter`` is not ``None``), affection
       status (if ``affectionFormatter`` is not ``None``) and loci names. If
       genotype at a locus needs more than one column, ``_1``, ``_2`` etc will
       be appended to loci names. Alternatively, a complete header (a string)
       or a list of column names could be specified directly.
   
   infoFormatter
       A format string that is used to format all information fields. If
       unspecified, ``str(value)`` will be used for each information field.
   
   genoFormatter
       How to output genotype at specified loci. Acceptable values include
       ``None`` (output allele values), a dictionary with genotype as keys,
       (e.g. ``genoFormatter={(0,0):1, (0,1):2, (1,0):2, (1,1):3}``, or a function
       with genotype (as a tuple of integers) as inputs. The dictionary value
       or the return value of this function can be a single or a list of
       number or strings.
   
   sexFormatter
       How to output individual sex. Acceptable values include ``None`` (no
       output) or a dictionary with keys ``MALE`` and ``FEMALE``.
   
   affectionFormatter
       How to output individual affection status. Acceptable values include
       ``None`` (no output) or a dictionary with keys ``True`` and ``False``.
   
   delimiter
       Delimiter used to separate values, default to ','.
   
   subPopFormatter
       How to output population membership. Acceptable values include
       ``None`` (no output), a string that will be used for the column name, or
       ``True`` which uses 'pop' as the column name. If present, the column is
       written with the string represenation of the (virtual) subpopulation.
   
   This operator supports the usual applicability parameters such as begin,
   end, step, at, reps, and subPops. If subPops are specified, only
   individuals from specified (virtual) subPops are exported. Similar to
   other operators, parameter ``output`` can be an output specification string
   (``filename``, ``>>filename``, ``!expr``), filehandle (or any Python object
   with a ``write`` function), any python function. Unless explicitly stated for
   a particular format, this operator exports individuals from the current
   generation if there are multiple ancestral generations in the population.
   
   The Exporter class will make use of a progress bar to show the progress. The
   interface of the progress bar is by default determined by the global GUI status
   but you can also set it to, for example, ``gui=False`` to forcefully use a 
   text-based progress bar, or ``gui='batch'`` to suppress the progress bar.

   .. method:: Exporter.Exporter(format, output, begin=0, end=-1, step=1, at=[], reps=True, subPops=ALL_AVAIL, infoFields=[], gui=None, *args, **kwargs)

      Usage:
      
          PyOperator(func, param=None, begin=0, end=-1, step=1, at=[],
            reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])
      
      Details:
      
          Create a pure-Python operator that calls a user-defined function
          when it is applied. If this operator is applied before or after
          mating, your function should have form func(pop) or func(pop,
          param) where pop is the population to which the operator is
          applied, param is the value specified in parameter param. param
          will be ignored if your function only accepts one parameter.
          Althernatively, the function should have form func(ind) with
          optional parameters pop and param. In this case, the function will
          be called for all individuals, or individuals in subpopulations
          subPops. Individuals for which the function returns False will be
          removed from the population. This operator can therefore perform
          similar functions as operator DiscardIf.  If this operator is
          applied during mating, your function should accept parameters pop,
          off (or ind), dad, mom and param where pop is the parental
          population, and off or ind, dad, and mom are offspring and their
          parents for each mating event, and param is an optional parameter.
          If subPops are provided, only offspring in specified (virtual)
          subpopulations are acceptable.  This operator does not support
          parameters output, and infoFields. If certain output is needed, it
          should be handled in the user defined function func. Because the
          status of files used by other operators through parameter output
          is undetermined during evolution, they should not be open or
          closed in this Python operator.



Function importPopulation
-------------------------


.. function:: importPopulation(format, filename, *args, **kwargs)

   This function import and return a population from a file *filename* in
   specified *format*. Format-specific parameters can be used to define how the
   input should be interpreted and imported. This function supports the following
   file format.
   
   GENEPOP (http://genepop.curtin.edu.au/). For input file of this format, this
   function ignores the first title line, load the second line as loci names,
   and import genotypes of different POP sections as different subpopulations.
   This format accepts the following parameters:
   
   adjust
       Adjust alleles by specified value (default to 0 for no adjustment). This
       parameter is mostly used to convert alleles 1 and 2 in a GenePop file to
       alleles 0 and 1 (with adjust=-1) in simuPOP. Negative allele (e.g. missing
       value 0) will be imported as regular allele with module-dependent values
       (e.g. -1 imported as 255 for standard module).
   
   
   FSTAT (http://www2.unil.ch/popgen/softwares/fstat.htm). This format accepts
   the following parameters:
   
   adjust
       Adjust alleles by specified value (default to 0 for no adjustment). This
       parameter is mostly used to convert alleles 1 and 2 in a GenePop file to
       alleles 0 and 1 (with adjust=-1) in simuPOP. Negative allele (e.g. missing
       value 0) will be imported as regular allele with module-dependent values
       (e.g. -1 imported as 255 for standard module).
   
   Phylip (Joseph Felsenstein's Phylip format). This function ignores sequence
   names and import sequences in a haploid (default) or diploid population (if
   there are even number of sequences). An list of allele names are required to
   translate symbols to allele names. This format accepts the following
   parameters:
   
   alleleNames
       Names of alleles 0, 1, 2, ... as a single string (e.g. 'ACTG') or a list of 
       single-character strings (e.g. ['A', 'C', 'T', 'G']). This will be used to
       translate symbols into numeric alleles in simuPOP. Allele names will continue
       to be used as allele names of the returned population.
   
   ploidy
       Ploidy of the returned population, default to 1 (haploid). There should be
       even number of sequences if ploidy=2 (haploid) is specified.
   
   MS (output from Richard R. Hudson's MS or msHOT program). The ms program generates
   npop blocks of nseq haploid chromosomes for command starting with 
   ``ms nsample nrepeat``. By default, the result is imported as a haploid
   population of size nsample. The population will have nrepeat subpopulations
   each with the same number of loci but different number of segregating sites.
   This behavior could be changed by the following parameters:
   
   ploidy
       If ``ploidy`` is set to 2, the sequenences will be paired so the population
       will have ``nseq/2`` individuals. An error will be raised if an odd number
       of sequences are simulated.
       
   mergeBy
       By default, replicate samples will be presented as subpopulations. All
       individuals have the same number of loci but individuals in different
       subpopulations have different segregating sites. If ``mergeBy`` is set
       to ``"chrom"``, the replicates will be presented as separate chromosomes,
       each with a different set of loci determined by segregating sites.


Function export
---------------


.. function:: export(pop, format, *args, **kwargs)

   Apply operator ``Exporter`` to population *pop* in format *format*.


