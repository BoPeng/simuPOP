Module :mod:`simuPOP.utils`
===========================

The :mod:`simuPOP.utils` module provides a few utility functions and classes.
They do not belong to the simuPOP core but are distributed with simuPOP because
they are frequently used and play an important role in some specialized
simulation techniques. Please refer to the simuPOP online cookbook
(``http://simupop.sourceforge.net/cookbook``) for more utility modules and
functions.


Trajectory simulation (classes ``Trajectory`` and ``TrajectorySimulator``)
--------------------------------------------------------------------------

A forward-time simulation, by its nature, is directly influenced by random
genetic drift. Starting from the same parental generation, allele frequencies in
the offspring generation would vary from simulation to simulation, with perhaps
a predictable mean frequency which is determined by factors such as parental
allele frequency, natural selection, mutation and migration.

Genetic drift is unavoidable and is in many cases the target of theoretical and
simulation studies. However, in certain types of studies, there is often a need
to control the frequencies of certain alleles in the present generation. For
example, if we are studying a particular penetrance model with pre-specified
frequencies of disease predisposing alleles, the simulated populations would
better have consistent allele frequencies at the disease predisposing loci, and
consequently consistent disease prevalence.

simuPOP provides a special offspring generator
:class:`ControlledOffspringGenerator` and an associated mating scheme called
:class:`ControlledRandomMating` that can be used to generate offspring
generations conditioning on frequencies of one or more alleles. This offspring
generator essentially uses a reject-sampling algorithm to select (or reject)
offspring according to their genotypes at specified loci. A detailed description
of this algorithm is given in Peng2007a.

The controlled random mating scheme accepts a user-defined trajectory function
that tells the mating scheme the desired allele frequencies at each generation.
Example :ref:`controlledOffGenerator <controlledOffGenerator>` uses a manually
defined function that raises the frequency of an allele steadily. However, given
known demographic and genetic factors, **a trajectory should be simulated
randomly so that it represents a random sample from all possible trajectories
that match the allele frequency requirement**. If such a condition is met, the
controlled evolutionary process can be considered as a random process
conditioning on allele frequencies at the present generation. Please refer to
Peng2007a for a detailed discussion about the theoretical requirements of a
valid trajectory simulator.

The ``simuUtil`` module provides functions and classes that implement two
trajectory simulation methods that can be used in different situations. The
first class is ``TrajectorySimulator`` which takes a demographic model and a
selection model as its input and simulates allele frequency trajectories using a
forward or backward algorithm. The demographic model is given by parameter
``N``, which can be a constant (e.g. ``N=1000``) for constant population size, a
list of subpopulation sizes (e.g. ``N=[1000, 2000]``) for a structured
population with constant size, or a demographic function that returns population
or subpopulation sizes at each generation. In the last case, subpopulations can
be split or merged with the constrait that subpopulations can be merged into
one, from split from one population.

A fitness model specifies the fitness of genotypes at one or more loci using
parameter ``fitness``. It can be a list of three numbers (e.g. ``fitness=[1,
1.001, 1.003]``), repsenting the fitness of genotype ``AA``, ``Aa`` and ``aa``
at one or more loci; or different fitness for genotypes at each locus (e.g.
``fitness=[1, 1.001, 1.003, 1, 1, 1.002]``), or for each combination or genotype
(interaction). In the last case, :math:`3^{n}` values are needed for each
genotype if there are :math:`n` loci. This trajectory simulator also accepts
generation-specific fitness values by accepting a function that returns fitness
values at each generation.

The simulator then simulates trajectories of allele frequencies and return them
as objects of class ``Trajectory``. This object can be used provide a trajectory
function that can be used directly in a :class:`ControlledRandomMating` mating
scheme (function :meth:`~simuPOP.utils.Trajectory.func`\ ()) or provide a list
of :class:`PointMutator` to introduce mutants at appropriate generations
(function :meth:`~simuPOP.utils.Trajectory.mutators`\ ()). If a simulation
failed after specified number of attempts, a ``None`` object will be returned.


Forward-time trajectory simulations (function ``simulateForwardTrajectory``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A forward simulation starts from a specified generation with specified allele
frequencies at one or more loci. The simulator simulates allele frequencies
forward-in-time, until it reaches a specified ending generation. A trajectory
object will be returned if the simulated allele frequencies fall into specified
ranges. Example :ref:`forwardTrajectory <forwardTrajectory>` demonstrates how to
use this simulation method to obtain and use a simulated trajectory, for two
unlinked loci under different selection pressure.

.. _forwardTrajectory:

**Example**: *Simulation and use of forward-time simulated trajectories.*

::

   >>> import simuOpt
   >>> simuOpt.setOptions(quiet=True)
   >>> import simuPOP as sim
   >>> from simuPOP.utils import Trajectory, simulateForwardTrajectory
   >>> 
   >>> traj = simulateForwardTrajectory(N=[2000, 4000], fitness=[1, 0.99, 0.98],
   ...     beginGen=0, endGen=100, beginFreq=[0.2, 0.3],
   ...     endFreq=[[0.1, 0.11], [0.2, 0.21]])
   >>> # 
   >>> #traj.plot('log/forwardTrajectory.png', set_ylim_top=0.5,
   >>> #    plot_c_sp=['r', 'b'], set_title_label='Simulated Trajectory (forward-time)')
   >>> pop = sim.Population(size=[2000, 4000], loci=10, infoFields='fitness')
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.8, 0.2], subPops=0),
   ...         sim.InitGenotype(freq=[0.7, 0.3], subPops=1),
   ...         sim.PyOutput('Sp0: loc2\tloc5\tSp1: loc2\tloc5\n'),
   ...     ],
   ...     matingScheme=sim.ControlledRandomMating(
   ...         ops=[sim.Recombinator(rates=0.01)],
   ...         loci=5, alleles=1, freqFunc=traj.func()),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=[2, 5], vars=['alleleFreq_sp'], step=20),
   ...         sim.PyEval(r"'%.2f\t%.2f\t%.2f\t%.2f\n' % (subPop[0]['alleleFreq'][2][1],"
   ...             "subPop[0]['alleleFreq'][5][1], subPop[1]['alleleFreq'][2][1],"
   ...             "subPop[1]['alleleFreq'][5][1])", step=20)
   ...     ],
   ...     gen = 101
   ... )
   Sp0: loc2	loc5	Sp1: loc2	loc5
   0.19	0.20	0.30	0.29
   0.20	0.20	0.29	0.27
   0.20	0.14	0.28	0.27
   0.17	0.13	0.27	0.26
   0.14	0.13	0.31	0.23
   0.13	0.10	0.27	0.20
   101

   now exiting runScriptInteractively...

`Download forwardTrajectory.py <forwardTrajectory.py>`_

Figure :ref:`fig_forwardTrajectory <fig_forwardTrajectory>` plots simulated
trajectories of one locus in two subpopulations. The plot function uses either
rpy or matplotlib as the underlying plotting library.

**Figure**: *Simulated trajectories of one locus in two subpopulations*

.. _fig_forwardTrajectory:

.. figure:: log/forwardTrajectory.png
   :width: 680



Backward-time trajectory simulations (function ``simulateBackwardTrajectory``).
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A backward simulation starts from specified frequencies at the present
generation. In the single-allele case, the simulations goes backward-in-time
until an allele gets lost. The length of such a trajectory is random, which is
usually a desired property because the age of a mutant in the present generation
is usually unknown and is assumed to be random.

This trajectory simulation technique is usually used as follows:

#. Determine a demographic and a natural selection model using which a forward-
   time simulation will be performed.

#. Given current disease allele frequencies, simulate trajectories of allele
   frequencies at each DSL using a backward approach.

#. Evolve a population forward-in-time, using designed demographic and selection
   models. A :class:`ControlledRandomMating` scheme instead of the usual
   :class:`RandomMating` scheme should be used.

Figure :ref:`fig_backTrajectory <fig_backTrajectory>` plots simulated
trajectories of two unlinked loci.

**Figure**: *Simulated trajectories of two unlinked loci*

.. _fig_backTrajectory:

.. figure:: log/backTrajectory.png
   :width: 680


The trajectory is used in a :class:`ControlledRandomMating` scheme in the
following evolutionary scenario:

.. _backTrajectory:

**Example**: *Simulation and use of backward-time simulated trajectories.*

::

   >>> import simuPOP as sim
   >>> from simuPOP.utils import Trajectory, simulateBackwardTrajectory
   >>> from math import exp
   >>> def Nt(gen):
   ...     'An exponential sim.Population growth demographic model.'
   ...     return int((5000) * exp(.00115 * gen))
   ... 
   >>> def fitness(gen, sp):
   ...     'Constant positive selection pressure.'
   ...     return [1, 1.01, 1.02]
   ... 
   >>> # simulate a trajectory backward in time, from generation 1000
   >>> traj = simulateBackwardTrajectory(N=Nt, fitness=fitness, nLoci=2,
   ...      endGen=1000, endFreq=[0.1, 0.2])
   >>> # matplotlib syntax
   >>> #traj.plot('log/backTrajectory.png', set_ylim_top=0.3, set_ylim_bottom=0,
   >>> #        plot_c_loc=['r', 'b'], set_title_label='Simulated Trajectory (backward-time)')
   >>> 
   >>> print('Trajectory simulated with length %s ' % len(traj.traj))
   Trajectory simulated with length 834 
   >>> pop = sim.Population(size=Nt(0), loci=[1]*2)
   >>> # save Trajectory function in the sim.population's local namespace
   >>> # so that the sim.PyEval operator can access it.
   >>> pop.dvars().traj = traj.func()
   >>> pop.evolve(
   ...     initOps=[sim.InitSex()],
   ...     preOps=traj.mutators(loci=[0, 1]),
   ...     matingScheme=sim.ControlledRandomMating(loci=[0, 1], alleles=[1, 1],
   ...         subPopSize=Nt, freqFunc=traj.func()),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=[0, 1], begin=500, step=100),
   ...         sim.PyEval(r"'%4d: %.3f (exp: %.3f), %.3f (exp: %.3f)\n' % (gen, alleleFreq[0][1],"
   ...             "traj(gen)[0], alleleFreq[1][1], traj(gen)[1])",
   ...             begin=500, step=100)
   ...     ],
   ...     gen=1001  # evolve 1001 generations to reach the end of generation 1000
   ... )
    500: 0.013 (exp: 0.013), 0.000 (exp: 0.000)
    600: 0.005 (exp: 0.005), 0.003 (exp: 0.003)
    700: 0.011 (exp: 0.011), 0.008 (exp: 0.008)
    800: 0.012 (exp: 0.013), 0.031 (exp: 0.031)
    900: 0.037 (exp: 0.037), 0.092 (exp: 0.092)
   1000: 0.101 (exp: 0.100), 0.200 (exp: 0.200)
   1001

   now exiting runScriptInteractively...

`Download backTrajectory.py <backTrajectory.py>`_


Graphical or text-based progress bar (class ``ProgressBar``)
------------------------------------------------------------

If your simulation takes a while to finish, you could use a progress bar to
indicate its progress. The ``ProgressBar`` class is provided for such a purpose.
Basically, you create a ``ProgressBar`` project with intended total steps, and
calls its ``update`` member function with each progress. Depending on available
graphical toolkit and the global or local GUI settings, a ``wxPython`` based
dialog, a ``Tkinter`` based dialog, or a text-based dialog will be used. Example
:ref:`ProgressBar <ProgressBar>` demonstrates how to use a text-based progress
bar. If the progress bar is updated at each step (such as in this example),
function ``update()`` can be called without parameter because it updates the
progress bar at an increment of 1 in this case.

.. _ProgressBar:

**Example**: *Using a text-based progress bar*

::

   >>> import simuPOP as sim
   >>> from simuPOP.utils import ProgressBar
   >>> pop = sim.Population(10000, loci=[10], infoFields='index')
   >>> prog = ProgressBar('Setting individual genotype...\n', pop.popSize(), gui=False)
   Setting individual genotype...
   >>> for idx in range(pop.popSize()):
   ...     # do something to each individaul
   ...     pop.individual(idx).index = idx
   ...     # idx + 1 can be ignored in this case.
   ...     prog.update(idx + 1)
   ... 
   ....1....2....3....4....5....6....7....8....9.... Done.

   now exiting runScriptInteractively...

`Download ProgressBar.py <ProgressBar.py>`_


Display population variables (function ``viewVars``)
----------------------------------------------------

If a population has a large number of variables, or if you are not sure which
variable to output, you could use function ``viewVars`` to view the population
variables in a tree form. If wxPython is available, a dialog could be used to
view the variables interactively. Example :ref:`viewVars <viewVars>`
demonstrates how to use this function. The wxPython-based dialog is displayed in
Figure :ref:`viewVars <viewVars>`.

.. _viewVars:

**Example**: *Using function viewVars to display population variables*

::

   import simuPOP as sim
   from simuPOP.utils import viewVars
   pop = sim.Population([1000, 2000], loci=3)
   sim.initGenotype(pop, freq=[0.2, 0.4, 0.4], loci=0)
   sim.initGenotype(pop, freq=[0.2, 0.8], loci=2)
   sim.stat(pop, genoFreq=[0, 1, 2], haploFreq=[0, 1, 2],
       alleleFreq=range(3),
       vars=['genoFreq', 'genoNum', 'haploFreq', 'alleleNum_sp'])
   viewVars(pop.vars())

`Download viewVars.py <viewVars.py>`_

**Figure**: *Using wxPython to display population variables*

.. _fig_viewvars:

.. figure:: /Users/bpeng1/simuPOP/simuPOP/doc/figures/viewVars.png
   :width: 680



Import simuPOP population from files in ``GENEPOP, PHYLIP`` and ``FSTAT`` formats (function ``importPopulation``)
-----------------------------------------------------------------------------------------------------------------

A function importPopulation is provided in the :mod:`simuPOP.utils` module to
import populations from files in ``GENEPOP, PHYLIP`` and ``FSTAT`` formats.
Because these formats do not support many of the features of a simuPOP
population, this function can only import genotype and basic information of a
population. Because formats GENEPOP and FSTAT formats uses allele 0 to indicate
missing value, true alleles in these formats start at value 1. If you would like
to import alleles with starting value 0, you can use parameter adjust=-1 to
adjust imported values, if you data do not have any missing value.


Export simuPOP population to files in ``STRUCTURE, GENEPOP``, ``FSTAT, Phylip, PED, MAP, MS,`` and ``CSV`` formats (function ``export`` and operator ``Exporter``)
------------------------------------------------------------------------------------------------------------------------------------------------------------------

simuPOP uses a program-specific binary format to save and load populations but
you can use the ``export`` function to export a simuPOP population in other
formats if you would like to use other programs to analyze simulated
populations. An operator Exporter is also provided so that you could export
populations during evolution. Operator arameters such as output, begin, end,
step, at, reps, and subPops are supported so that you could export subsets of
individuals at multiple generations using different file names (e.g.
``output='!''%d.ped'' % gen'`` to output to different files at different
generations).

Commonly used population genetics file formats such as GENEPOP, FSTAT, Phylip,
MS, and STRUCTURE are supported. Because these formats cannot store all
information in a simuPOP population, export and import operations can lose
information. Also, because the processing application have different
assumptions, some conversion of genotypes might be needed. For example, because
GENEPOP uses allele 0 as missing genotype, function ``export(format='genepop')``
accepts a parameter ``adjust`` with default value ``1`` to export alleles 0, 1
etc to 1, 2, .... The same applies to function ``importPopulation`` where some
file formats accepts a parameter ``adjust`` (with default value 1) to adjust
allele values. Please refer to the simuPOP reference manual for a detailed list
of acceptable parameters for each format.

Example :ref:`importExport <importExport>` demonstrates how to import and export
a population in formats FSTAT and STRUCTURE. For the FSTAT format, because the
population is exported with allele values shifted by 1, the imported population
has different alleles than the original population. This can be fixed by adding
parameter ``adjust=-1`` to the ``importPopulation`` function.

.. _importExport:

**Example**: *Save and load a population*

::

   >>> import simuPOP as sim
   >>> from simuPOP.utils import importPopulation, export
   >>> pop = sim.Population([2,4], loci=5, lociNames=['a1', 'a2', 'a3', 'a4', 'a5'],
   ...     infoFields='BMI')
   >>> sim.initGenotype(pop, freq=[0.3, 0.5, 0.2])
   >>> sim.initSex(pop)
   >>> sim.initInfo(pop, [20, 30, 40, 50, 30, 25], infoFields='BMI')
   >>> export(pop, format='fstat', output='fstat.txt')
   Exporting....1....2....3....4....5....6....7....8....9.... Done.
   >>> print(open('fstat.txt').read())
   2 5 3 1
   a1
   a2
   a3
   a4
   a5
   1 21 21 23 12 12
   1 22 23 22 22 21
   2 31 21 22 11 13
   2 22 22 33 23 21
   2 22 32 33 22 21
   2 33 33 22 21 32

   >>> export(pop, format='structure', phenotype='BMI', output='stru.txt')
   Exporting....1....2....3....4....5....6....7....8....9.... Done.
   >>> print(open('stru.txt').read())
   a1	a2	a3	a4	a5
   -1	1.0	1.0	1.0	1.0
   1	1	20	1	1	1	0	0
   1	1	20	0	0	2	1	1
   2	1	30	1	1	1	1	1
   2	1	30	1	2	1	1	0
   1	2	40	2	1	1	0	0
   1	2	40	0	0	1	0	2
   2	2	50	1	1	2	1	1
   2	2	50	1	1	2	2	0
   3	2	30	1	2	2	1	1
   3	2	30	1	1	2	1	0
   4	2	25	2	2	1	1	2
   4	2	25	2	2	1	0	1

   >>> pop1 = importPopulation(format='fstat', filename='fstat.txt')
   >>> sim.dump(pop1)
   Ploidy: 2 (diploid)
   Chromosomes:
   1:  (AUTOSOME, 5 loci)
     a1 (1), a2 (2), a3 (3), a4 (4), a5 (5)
   population size: 6 (2 subpopulations with 2 (1), 4 (2) Individuals)
   Number of ancestral populations: 0

   SubPopulation 0 (1), 2 Individuals:
      0: MU 22211 | 11322 
      1: MU 22222 | 23221 
   SubPopulation 1 (2), 4 Individuals:
      2: MU 32211 | 11213 
      3: MU 22322 | 22331 
      4: MU 23322 | 22321 
      5: MU 33223 | 33212 


   now exiting runScriptInteractively...

`Download importExport.py <importExport.py>`_

Because coalescent simulations are increasingly used to generate initial
populations in equilibrium stats, importing data in MS format is very useful.
Because MS only simulates haploid sequences with genotype only at segregating
sites, you might have to simulate an even number of sequences and use option
ploidy=2 to import the simulated data as a haploid population. In addition, a
parameter mergeBy is provided to import multiple replicates as multiple
subpopulations or chromosomes. This corresponds to the splitBy parameter when
you export your data in MS format. Example :ref:`importMS <importMS>`
demonstrates how to use these parameters.

.. _importMS:

**Example**: *Export and import in MS format*

::

   >>> import simuPOP as sim
   >>> from simuPOP.utils import importPopulation, export
   >>> pop = sim.Population([20,20], loci=[10, 10])
   >>> # simulate a population but mutate only a subset of loci
   >>> pop.evolve(
   ...     preOps=[
   ...         sim.InitSex(),
   ...         sim.SNPMutator(u=0.1, v=0.01, loci=range(5, 17))
   ...     ],
   ...     matingScheme=sim.RandomMating(),
   ...     gen=100
   ... )
   100
   >>> # export first chromosome, all individuals
   >>> export(pop, format='ms', output='ms.txt')
   Exporting....1....2....3....4....5....6....7....8....9.... Done.
   >>> # export first chromosome, subpops as replicates
   >>> export(pop, format='ms', output='ms_subPop.txt', splitBy='subPop')
   Exporting....1....2....3....4....5....6....7....8....9.... Done.
   >>> # export all chromosomes, but limit to all males in subPop 1
   >>> pop.setVirtualSplitter(sim.SexSplitter())
   >>> export(pop, format='ms', output='ms_chrom.txt', splitBy='chrom', subPops=[(1,0)])
   Exporting....1....2....3....4....5....6....7....8....9.... Done.
   >>> # 
   >>> print(open('ms_chrom.txt').read())
   simuPOP_export 20 2
   30164 48394 29292

   //
   segsites: 5
   positions: 6.0 7.0 8.0 9.0 10.0
   11110
   11111
   11110
   11111
   11011
   11111
   01111
   10111
   11111
   11111
   01111
   01111
   11011
   11111
   01111
   11011
   11101
   10111
   11111
   11111

   //
   segsites: 7
   positions: 1.0 2.0 3.0 4.0 5.0 6.0 7.0
   1101111
   1110011
   1101110
   1111111
   0111110
   1111111
   1110001
   1111111
   0111110
   1111111
   1111111
   1111111
   1111111
   1011111
   1111111
   1111111
   1011111
   1111111
   1111111
   1011111

   >>> # import as haploid sequence
   >>> pop = importPopulation(format='ms', filename='ms.txt')
   >>> # import as diploid 
   >>> pop = importPopulation(format='ms', filename='ms.txt', ploidy=2)
   >>> # import as a single chromosome
   >>> pop = importPopulation(format='ms', filename='ms_subPop.txt', mergeBy='subPop')

   now exiting runScriptInteractively...

`Download importMS.py <importMS.py>`_

If the file format you are interested in is not supported, you can export data
in csv format and convert the file by yourself. You can also try to write your
own import or export functions as described in the advanced topics section of
this guide.


Export simuPOP population in csv format (function ``saveCSV``, deprecated)
--------------------------------------------------------------------------

Function ``saveCSV`` is provided in the :mod:`simuPOP.utils` module to save (the
present generation of) a simuPOP population in comma separated formats. It
allows you to save individual information fields, sex, affection status and
genotype (in that order). Because this function allows you to output these
information in different formats using parameters ``infoFormatter``,
``sexFormatter``, ``affectionFormatter``, and ``genoFormatter``, this function
can already be used to export a simuPOP population to formats that are
recognizable by some populat software applications. Example :ref:`saveCSV
<saveCSV>` creates a small population and demonstrates how to save it in
different formats.

.. _saveCSV:

**Example**: *Using function saveCSV to save a simuPOP population in different formats*

::

   >>> import simuPOP as sim
   >>> from simuPOP.utils import saveCSV
   >>> pop = sim.Population(size=[10], loci=[2, 3],
   ...     lociNames=['r11', 'r12', 'r21', 'r22', 'r23'],
   ...     alleleNames=['A', 'B'], infoFields='age')
   >>> sim.initSex(pop)
   >>> sim.initInfo(pop, [2, 3, 4], infoFields='age')
   >>> sim.initGenotype(pop, freq=[0.4, 0.6])
   >>> sim.maPenetrance(pop, loci=0, penetrance=(0.2, 0.2, 0.4))
   >>> # no filename so output to standard output
   >>> saveCSV(pop, infoFields='age')
   age, sex, aff, r11_1, r11_2, r12_1, r12_2, r21_1, r21_2, r22_1, r22_2, r23_1, r23_2
   2.0, F, A, B, B, B, B, B, A, B, B, B, A
   3.0, F, U, B, A, B, A, B, A, A, A, A, B
   4.0, M, U, B, B, B, B, B, B, B, B, B, A
   2.0, M, U, B, A, B, A, B, B, B, B, B, A
   3.0, M, A, B, B, B, B, B, B, A, A, B, A
   4.0, M, U, A, B, B, A, B, B, B, B, B, B
   2.0, M, U, B, B, B, B, B, B, B, B, A, A
   3.0, F, U, B, B, A, A, B, B, A, A, B, B
   4.0, F, U, A, B, B, B, B, B, B, A, B, B
   2.0, F, A, B, A, A, B, A, A, B, B, B, A
   >>> # change affection code and how to output genotype
   >>> saveCSV(pop, infoFields='age', affectionFormatter={True: 1, False: 2},
   ...     genoFormatter={(0,0):'AA', (0,1):'AB', (1,0):'AB', (1,1):'BB'})
   age, sex, aff, r11, r12, r21, r22, r23
   2.0, F, 1, BB, BB, AB, BB, AB
   3.0, F, 2, AB, AB, AB, AA, AB
   4.0, M, 2, BB, BB, BB, BB, AB
   2.0, M, 2, AB, AB, BB, BB, AB
   3.0, M, 1, BB, BB, BB, AA, AB
   4.0, M, 2, AB, AB, BB, BB, BB
   2.0, M, 2, BB, BB, BB, BB, AA
   3.0, F, 2, BB, AA, BB, AA, BB
   4.0, F, 2, AB, BB, BB, AB, BB
   2.0, F, 1, AB, AB, AA, BB, AB
   >>> # save to a file
   >>> saveCSV(pop, filename='pop.csv', infoFields='age', affectionFormatter={True: 1, False: 2},
   ...     genoFormatter=lambda geno: (geno[0] + 1, geno[1] + 1), sep=' ')
   >>> print(open('pop.csv').read())
   age sex aff r11_1 r11_2 r12_1 r12_2 r21_1 r21_2 r22_1 r22_2 r23_1 r23_2
   2.0 F 1 2 2 2 2 2 1 2 2 2 1
   3.0 F 2 2 1 2 1 2 1 1 1 1 2
   4.0 M 2 2 2 2 2 2 2 2 2 2 1
   2.0 M 2 2 1 2 1 2 2 2 2 2 1
   3.0 M 1 2 2 2 2 2 2 1 1 2 1
   4.0 M 2 1 2 2 1 2 2 2 2 2 2
   2.0 M 2 2 2 2 2 2 2 2 2 1 1
   3.0 F 2 2 2 1 1 2 2 1 1 2 2
   4.0 F 2 1 2 2 2 2 2 2 1 2 2
   2.0 F 1 2 1 1 2 1 1 2 2 2 1


   now exiting runScriptInteractively...

`Download saveCSV.py <saveCSV.py>`_

**This function is now deprecated with the introduction of function
**``export``** and operator **``Exporter``**.**


