Non-random and customized mating schemes \*
===========================================


The structure of a homogeneous mating scheme \*
-----------------------------------------------

A *homogeneous mating scheme* (:class:`HomoMating`) populates an offspring
generation as follows:

#. Create an empty offspring population (generation) with appropriate size.
   Parental and offspring generation can differ in size but they must have the same
   number of subpopulations.

#. For each subpopulation, repeatedly choose a parent or a pair of parents from
   the parental generation. This is done by a simuPOP object called a **parent
   chooser**.

#. One or more offspring are produced from the chosen parent(s) and are placed
   in the offspring population. This is done by a simuPOP **offspring generator**.
   An offspring generator uses one or more during-mating operators to transmit
   parental genotype to offspring. These operators are called **genotype
   transmitters**.

#. After the offspring generation is populated, it will replace the parental
   generation and becomes the present generation of a population.

To define a homogeneous mating scheme, you will need to provide a ``chooser`` (a
*parent chooser* that is responsible for choosing one or two parents from the
parental generation) and a ``generator`` (an *offspring generator* that is
responsible for generating a number of offspring from the chosen parents). For
example, a ``selfingMating`` mating scheme uses a :class:`RandomParentChooser`
to choose a parent randomly from a population, possibly according to individual
fitness, it uses a standard :class:`OffspringGenerator` that uses a
``selfingOffspringGenerator`` to transmit genotype. The constructor of
:class:`HomoMating` also accepts parameters ``subPopSize`` (parameter to control
offspring subpopulation sizes), ``subPops`` (applicable subpopulatiosn or
virtual subpopulations), and ``weight`` (weighting parameter when used in a
heterogeneous mating scheme). When this mating scheme is applied to the whole
population, ``subPopSize`` is used to determine the subpopulation sizes of the
offspring generation (see Section :ref:`subsec_offspring_size
<subsec_offspring_size>` for details), parameters ``subPops`` and ``weight`` are
ignored. Otherwise, the number of offspring this mating scheme will produce is
determined by the heterogeneous mating scheme.

Example :ref:`RandomMating <RandomMating>` demonstrates how the most commonly
used mating scheme, the diploid sexual :class:`RandomMating` mating scheme is
defined in ``simuPOP.py``. This mating scheme uses a
:class:`RandomParentsChooser` with replacement, and a standard
:class:`OffspringGenerator` using a default :class:`MendelianGenoTransmitter`.

.. _RandomMating:

**Example**: *Define a random mating scheme*

::

   def RandomMating(numOffspring=1., sexMode=RANDOM_SEX,
           ops=MendelianGenoTransmitter(), subPopSize=[],
           subPops=ALL_AVAIL, weight=0, selectionField='fitness'):
       'A basic diploid sexual random mating scheme.'
       return HomoMating(
           chooser=RandomParentsChooser(True, selectionField),
           generator=OffspringGenerator(ops, numOffspring, sexMode),
           subPopSize=subPopSize,
           subPops=subPops,
           weight=weight)


`Download RandomMating.py <RandomMating.py>`_

Different parent choosers and offspring generators can be combined to define a
large number of homogeneous mating schemes. Some of the parent choosers return
one parent so they work with offspring generators that need one parent (e.g.
selfing or clone offspring generator); some of the parent choosers return two
parents so they work with offspring generators that need two parents (e.g.
Mendelian offspring generator). For example, the standard :class:`SelfMating`
mating scheme uses a :class:`RandomParentChooser` but you can easily use a
:class:`SequentialParentChooser` to choose parents sequentially and self-
fertilize parents one by one. This is demonstrated in Example
:ref:`sequentialSelfing <sequentialSelfing>`.

.. _sequentialSelfing:

**Example**: *Define a sequential selfing mating scheme*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(100, loci=5*3, infoFields='parent_idx')
   >>> pop.evolve(
   ...     initOps=[sim.InitGenotype(freq=[0.2]*5)],
   ...     preOps=sim.Dumper(structure=False, max=5),
   ...     matingScheme=sim.HomoMating(
   ...         sim.SequentialParentChooser(),
   ...         sim.OffspringGenerator(ops=[
   ...             sim.SelfingGenoTransmitter(),
   ...             sim.ParentsTagger(infoFields='parent_idx'),
   ...         ])
   ...     ),
   ...     postOps=sim.Dumper(structure=False, max=5),
   ...     gen = 1
   ... )
   SubPopulation 0 (), 100 Individuals:
      0: MU 441000142224423 | 431303440010114 |  0
      1: MU 334442443034342 | 113203441333201 |  0
      2: MU 034344042424240 | 344304121430212 |  0
      3: MU 132322330420043 | 141300223114240 |  0
      4: MU 111123040033342 | 344344221133120 |  0

   SubPopulation 0 (), 100 Individuals:
      0: MU 441000142224423 | 431303440010114 |  0
      1: FU 334442443034342 | 113203441333201 |  1
      2: MU 344304121430212 | 034344042424240 |  2
      3: FU 141300223114240 | 132322330420043 |  3
      4: FU 344344221133120 | 111123040033342 |  4

   1

   now exiting runScriptInteractively...

`Download sequentialSelfing.py <sequentialSelfing.py>`_


Offspring generators \*
-----------------------

An :class:`OffspringGenerator` accepts a parameters ``ops`` (a list of during-
mating operators), ``numOffspring`` (control number of offspring per mating
event) and ``sexMode`` (control offspring sex). We have examined the last two
parameters in detail in sections :ref:`subsec_number_of_offspring
<subsec_number_of_offspring>` and :ref:`subsec_offspring_sex
<subsec_offspring_sex>`.

The most tricky parameter is the ``ops`` parameter. It accepts a list of during
mating operators that are used to transmit genotypes from parent(s) to offspring
and/or set individual information fields. The standard
:class:`OffspringGenerator` does not have any default operator so no genotype
will be transmitted by default. All stock mating schemes use a default genotype
transmitter. (e,g, a :class:`MendelianGenoTransmitter` in Example
:ref:`RandomMating <RandomMating>` is passed to the offspring generator used in
:class:`RandomMating`). Note that you need to specify all needed operators if
you use parameter ``ops`` to change the operators used in a mating scheme (see
Example :ref:`HeteroMatingWeight <HeteroMatingWeight>`). That is to say, you can
use ``ops=Recombinator()`` to replace a default
:class:`MendelianGenoTransmitter`\ (), but you have to use ``ops=[IdTagger(),
MendelianGenoTransmitter()]`` if you would like to add a during-mating operator
to the default one.

Another offspring generator is provided in simuPOP. This
``ControlledOffspringGenerator``\ is used to control an evolutionary process so
that the allele frequencies at certain loci follows some pre-simulated
*frequency trajectories*. Please refer to Peng2007a for rationals behind such an
offspring generator and its applications in the simulation of complex human
diseases.

Example :ref:`controlledOffGenerator <controlledOffGenerator>` demonstrates the
use of such a controlled offspring generator. Instead of using a realistic
frequency trajectory function, it forces allele frequency at locus 5 to increase
linearly. In contrast, the allele frequency at locus 15 on the second chromosome
oscillates as a result of genetic drift. Note that the random mating version of
this mating scheme is defined in simuPOP as :class:`ControlledRandomMating`.

.. _controlledOffGenerator:

**Example**: *A controlled random mating scheme*

::

   >>> import simuPOP as sim
   >>> def traj(gen):
   ...     return [0.5 + gen * 0.01]
   ... 
   >>> pop = sim.Population(1000, loci=[10]*2)
   >>> # evolve the sim.Population while keeping allele frequency 0.5
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.5, 0.5])
   ...     ],
   ...     matingScheme=sim.HomoMating(sim.RandomParentChooser(),
   ...         sim.ControlledOffspringGenerator(loci=5,
   ...             alleles=[0], freqFunc=traj,
   ...             ops = sim.SelfingGenoTransmitter())),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=[5, 15]),
   ...         sim.PyEval(r'"%.2f\t%.2f\n" % (alleleFreq[5][0], alleleFreq[15][0])')
   ...     ],
   ...     gen = 5
   ... )
   0.50	0.51
   0.51	0.51
   0.52	0.51
   0.53	0.52
   0.54	0.54
   5

   now exiting runScriptInteractively...

`Download controlledOffGenerator.py <controlledOffGenerator.py>`_


.. _subsec_Pre_defined_genotype_transmitters:

Genotype transmitters \*
------------------------

Although any during mating operators can be used in parameter ``ops``\ of an
offspring generator, those that transmit genotype from parents to offspring are
customarily called **genotype transmitters**. simuPOP provides a number of
genotype transmitters including clone, Mendelian, selfing, haplodiploid,
genotype transmitter, and a Recombinator. They are usually used implicitly in a
mating scheme, but they can also be used explicitly.

Although simuPOP provides a number of genotype transmitters, they may still be
cases where customized genotype transmitter is needed. For example, a
Recombinator can be used to recombine parental chromosomes but it is well known
that male and female individuals differ in recombination rates. How can you
apply two different Recombinators to male and female Individuals separately?

An immediate thought can be the use of virtual subpopulations. If you apply two
random mating schemes to two virtual subpopulations defined by sex,
:class:`RandomParentsChooser` will not work because no opposite sex can be found
in each virtual subpopulation. In this case, a customized genotype transmitter
can be used.

A customized genotype transmitter is only a Python during-mating operator.
Although it is possible to define a function and use a PyOperator directly
(Example :ref:`PyOperator <PyOperator>`), it is much better to derive an
operator from PyOperator, as the case in Example :ref:`newOperator
<newOperator>`.

Example :ref:`sexSpecificRec <sexSpecificRec>` defines a
``sexSpecificRecombinator`` that uses, internally, two different Recombinators
to recombine male and female parents. The key statement is the
``PyOperator.__init__`` line which initializes a Python operator with given
function ``self.transmitGenotype``. Example :ref:`sexSpecificRec
<sexSpecificRec>` outputs the population in two generations. You should notice
that paternal chromosome are not recombined when they are transmitted to
offspring.

.. _sexSpecificRec:

**Example**: *A customized genotype transmitter for sex-specific recombination*

::

   >>> from simuPOP import *
   >>> class sexSpecificRecombinator(PyOperator):
   ...     def __init__(self, intensity=0, rates=0, loci=[], convMode=NO_CONVERSION,
   ...             maleIntensity=0, maleRates=0, maleLoci=[], maleConvMode=NO_CONVERSION,
   ...             *args, **kwargs):
   ...         # This operator is used to recombine maternal chromosomes
   ...         self.Recombinator = Recombinator(rates, intensity, loci, convMode)
   ...         # This operator is used to recombine paternal chromosomes
   ...         self.maleRecombinator = Recombinator(maleRates, maleIntensity,
   ...             maleLoci, maleConvMode)
   ...         #
   ...         PyOperator.__init__(self, func=self.transmitGenotype, *args, **kwargs)
   ...     #
   ...     def transmitGenotype(self, pop, off, dad, mom):
   ...         # Form the first homologous copy of offspring.
   ...         self.Recombinator.transmitGenotype(mom, off, 0)
   ...         # Form the second homologous copy of offspring.
   ...         self.maleRecombinator.transmitGenotype(dad, off, 1)
   ...         return True
   ... 
   >>> pop = Population(10, loci=[15]*2, infoFields=['father_idx', 'mother_idx'])
   >>> pop.evolve(
   ...     initOps=[
   ...         InitSex(),
   ...         InitGenotype(freq=[0.4] + [0.2]*3)
   ...     ],
   ...     matingScheme=RandomMating(ops=[
   ...         sexSpecificRecombinator(rates=0.1, maleRates=0),
   ...         ParentsTagger()
   ...     ]),
   ...     postOps=Dumper(structure=False),
   ...     gen = 2
   ... )
   SubPopulation 0 (), 10 Individuals:
      0: FU 230000130212000 130110020112120 | 310300000030330 000113003202000 |  6 7
      1: FU 110100000002000 223313300111002 | 331311301000220 002330110020020 |  6 7
      2: MU 230301121003012 032010332330303 | 303303022100031 310232031321031 |  5 0
      3: MU 103001320130222 031300110100023 | 303303022100031 003000012020002 |  5 9
      4: FU 210230113000000 231111000121000 | 303303022100031 003000012020002 |  5 8
      5: MU 322030133101023 110323303020211 | 322111021000001 301200303300133 |  2 8
      6: MU 210230113000000 231111000121000 | 331303300011323 310232031321031 |  5 8
      7: FU 200331312001001 200011203020203 | 031032120003212 101032020302120 |  3 1
      8: FU 230000130212000 223313300111002 | 303303022100031 003000012020002 |  5 7
      9: FU 200331312001001 130301011230300 | 322111021000001 320103032303101 |  2 1

   SubPopulation 0 (), 10 Individuals:
      0: MU 230000130212000 223313300111002 | 322030133101023 301200303300133 |  5 8
      1: MU 230000130212000 130110020112120 | 303303022100031 310232031321031 |  2 0
      2: FU 303303022100031 003000012020002 | 322111021000001 301200303300133 |  5 4
      3: FU 331311301000220 223313300111002 | 322111021000001 110323303020211 |  5 1
      4: MU 200331312001001 101032020302120 | 230301121003012 032010332330303 |  2 7
      5: FU 031032120003212 200011203020203 | 103001320130222 031300110100023 |  3 7
      6: FU 200331312001001 320103032303101 | 303303022100031 032010332330303 |  2 9
      7: FU 200331312001001 320103032303101 | 303303022100031 310232031321031 |  2 9
      8: FU 200331312001001 130301011230300 | 303303022100031 031300110100023 |  3 9
      9: MU 303303022100031 003000012020002 | 210230113000000 231111000121000 |  6 4

   2

   now exiting runScriptInteractively...

`Download sexSpecificRec.py <sexSpecificRec.py>`_


A Python parent chooser \*
--------------------------

Parent choosers are responsible for choosing one or two parents from a parental
(virtual) subpopulation. simuPOP defines a few parent choosers that choose
parent(s) sequentially, randomly (with or without replacement), or with
additional conditions. Some of these parent choosers support natual selection.
We have seen sequential and random parent choosers in Examples
:ref:`sequentialSelfing <sequentialSelfing>` and :ref:`controlledOffGenerator
<controlledOffGenerator>`. Please refer to the simuPOP reference manual for
details about these objects.

A parent choosing scheme can be quite complicated in reality. For example,
salamanders along a river may mate with their neighbors and form several
subspecies. This behavior cannot be readily simulated using any pre-define
parent choosers so a hybrid parent chooser :class:`PyParentsChooser`\ () should
be used.

A :class:`PyParentsChooser` accepts a user-defined Python generator function,
instead of a normal python function, that returns a parent, or a pair of parents
repeatedly. Briefly speaking, when a generator function is called, it returns a
*generator* object that provides an iterator interface. Each time when this
iterator iterates, this function resumes where it was stopped last time,
executes and returns what the next *yield* statement returns. For example,
example :ref:`generator <generator>` defines a function that calculate
:math:`f\left(k\right)=\sum_{i=1}^{k}\frac{1}{i}` for :math:`k=1,...,5`. It does
not calculate each :math:`f\left(k\right)` repeatedly but returns
:math:`f\left(1\right)`, :math:`f\left(2\right)`, ... sequentially.

.. _generator:

**Example**: *A sample generator function*

::

   >>> import simuPOP as sim
   >>> def func():
   ...     i = 1
   ...     all = 0
   ...     while i <= 5:
   ...         all += 1./i
   ...         i += 1
   ...         yield all 
   ... 
   >>> for i in func():
   ...     print('%.3f' % i)
   ... 
   1.000
   1.500
   1.833
   2.083
   2.283

   now exiting runScriptInteractively...

`Download generator.py <generator.py>`_

A :class:`PyParentsChooser` accepts a parent generator function, which takes a
population and a subpopulation index as parameters. When this parent chooser is
applied to a subpopulation, it will call this generator function and ask the
generated generator object repeated for either a parent, or a pair of parents
(*references to individual objects or indexes relative to a subpopulation*).
Note that :class:`PyParentsChooser` does not support virtual subpopulation but
you can mimic the effect by returning only parents from certain virtual
subpopulations.

Example :ref:`PyParentsChooser <PyParentsChooser>` implements a hybrid parent
chooser that chooses parents with equal social status (``rank``). In this parent
chooser, all males and females are categorized by their sex and social status. A
parent is chosen randomly, and then his/her spouse is chosen from females/males
with the same social status. The rank of their offspring can increase or
decrease randomly. It becomes obvious now that whereas a python function can
return random male/female pair, the generator interface is much more efficient
because the identification of sex/status groups is done only once.

.. _PyParentsChooser:

**Example**: *A hybrid parent chooser that chooses parents by their social status*

::

   >>> import simuPOP as sim
   >>> from random import randint
   >>> def randomChooser(pop, subPop):
   ...     males = []
   ...     females = []
   ...     # identify males and females in each social rank
   ...     for rank in range(3):
   ...         males.append([x for x in pop.individuals(subPop) \
   ...             if x.sex() == sim.MALE and x.rank == rank])
   ...         females.append([x for x in pop.individuals(subPop) \
   ...             if x.sex() == sim.FEMALE and x.rank == rank])
   ...     #
   ...     while True:
   ...         # choose a rank randomly
   ...         rank = int(pop.individual(randint(0, pop.subPopSize(subPop) - 1), subPop).rank)
   ...         yield males[rank][randint(0, len(males[rank]) - 1)], \
   ...             females[rank][randint(0, len(females[rank]) - 1)]
   ... 
   >>> def setRank(rank):
   ...     'The rank of offspring can increase or drop to zero randomly'
   ...     # only use rank of the father
   ...     return (rank[0] + randint(-1, 1)) % 3
   ... 
   >>> pop = sim.Population(size=[1000, 2000], loci=1, infoFields='rank')
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitInfo(lambda : randint(0, 2), infoFields='rank')
   ...     ],
   ...     matingScheme=sim.HomoMating(
   ...         sim.PyParentsChooser(randomChooser),
   ...         sim.OffspringGenerator(ops=[
   ...             sim.MendelianGenoTransmitter(),
   ...             sim.PyTagger(setRank),
   ...             ])
   ...     ),
   ...     gen = 5
   ... )    
   5

   now exiting runScriptInteractively...

`Download PyParentsChooser.py <PyParentsChooser.py>`_

Built-in parent choosers could be used in a :class:`PyParentsChooser` to choose
parents. The parent chooser needs to be initialized with the parental population
and subpopulation index. Calling the ``chooseParents`` function repeatedly will
return pairs of individuals from the population (``None`` will be returned for
one of the parents if the parent chooser only returns one parent). The use of
built-in parent choosers can improve the performance of your
:class:`PyParentsChooser`, especially for complex selection patterns (e.g. with
natural selection). For example, :ref:`BuiltInParentsChooser
<BuiltInParentsChooser>` implements a similar mating scheme as Example
:ref:`PyParentsChooser <PyParentsChooser>` but uses a
:class:`RandomParentChooser` to choose males randomly.

.. _BuiltInParentsChooser:

**Example**: *Use built-in parent choosers in a Python parent chooser*

::

   >>> import simuPOP as sim
   >>> from random import randint
   >>> 
   >>> def randomChooser(pop, subPop):
   ...     maleChooser = sim.RandomParentChooser(sexChoice=sim.MALE_ONLY)
   ...     maleChooser.initialize(pop, subPop)
   ...     females = []
   ...     # identify females in each social rank
   ...     for rank in range(3):
   ...         females.append([x for x in pop.individuals(subPop) \
   ...             if x.sex() == sim.FEMALE and x.rank == rank])
   ...     #
   ...     while True:
   ...         # choose a random male
   ...         m = maleChooser.chooseParents()[0]
   ...         rank = int(m.rank)
   ...         # find a female in the same rank
   ...         yield m, females[rank][randint(0, len(females[rank]) - 1)]
   ... 
   >>> def setRank(rank):
   ...     'The rank of offspring can increase or drop to zero randomly'
   ...     # only use rank of the father
   ...     return (rank[0] + randint(-1, 1)) % 3
   ... 
   >>> pop = sim.Population(size=[1000, 2000], loci=1, infoFields='rank')
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitInfo(lambda : randint(0, 2), infoFields='rank')
   ...     ],
   ...     matingScheme=sim.HomoMating(
   ...         sim.PyParentsChooser(randomChooser),
   ...         sim.OffspringGenerator(ops=[
   ...             sim.MendelianGenoTransmitter(),
   ...             sim.PyTagger(setRank),
   ...             ])
   ...     ),
   ...     gen = 5
   ... )    
   5

   now exiting runScriptInteractively...

`Download BuiltInParentsChooser.py <BuiltInParentsChooser.py>`_


.. _subsec_Using_C++:

Using C++ to implement a parent chooser \*\*
--------------------------------------------

A user defined parent chooser can be fairly complex and computationally
intensive. For example, if a parent tends to find a spouse in his/her vincinity,
geometric distances between all qualified individuals and a chosen parent need
to be calculated for each mating event. If the optimization of the parent
chooser can speed up the simulation significantly, it may be worthwhile to write
the parent chooser in C++.

Although it is feasible, and sometimes easier to derive a class from class
:class:`ParentChooser` in mating.h (.cpp), modifying simuPOP source code is not
recommended because you would have to modify a new version of simuPOP whenever
you upgrade your simuPOP distribution. Implementing your parent choosing
algorithm in another Python module is preferred.

The first step is to write your own parent chooser in C/C++. Basically, you will
need to pass all necessary information to the C++ level and implement an
algorithm to choose parents randomly. Although simple function based solutions
are possible, a C++ level class such as the ``myParentsChooser``\ class defined
in Example :ref:`parentChooseHeader <parentChooseHeader>` is recommended. This
class is initialized with indexes of male and female individuals and use a
function ``chooseParents`` to return a pair of parents randomly. This parent
chooser is very simple but more complicated parent selection scenarios can be
implemented similarly.

.. _parentChooseHeader:

**Example**: *Implement a parent chooser in C++*

::

   #include <stdlib.h>
   #include <vector>
   #include <utility>
   using std::pair;
   using std::vector;
   class myParentsChooser
   {
   public:
           // A constructor takes all locations of male and female.
           myParentsChooser(const std::vector<int> & m, const std::vector<int> & f)
                   : male_idx(m), female_idx(f)
           {
                   srand(time(0));
           }

           pair<unsigned long, unsigned long> chooseParents()
           {
                   unsigned long male = rand() % male_idx.size();
                   unsigned long female = rand() % male_idx.size();
                   return std::make_pair(male, female);
           }
   private:
           vector<int> male_idx;
           vector<int> female_idx;
   };

`Download myParentsChooser.h <myParentsChooser.h>`_

The second step is to wrap your C++ functions and classes to a Python module.
There are many tools available but SWIG (``www.swig.org``) is arguably the most
convenient and powerful one. To use SWIG, you will need to prepare an interface
file, which basically tells SWIG which functions and classes you would like to
expose and how to pass parameters between Python and C++. Example
:ref:`parentsChooserInterface <parentsChooserInterface>` lists an interface file
for the C++ class defined in Example :ref:`parentChooseHeader
<parentChooseHeader>`. Please refer to the SWIG reference manual for details.

.. _parentsChooserInterface:

**Example**: *An interface file for the myParentsChooser class*

::

   %module myParentsChooser
   %{
   #include "myParentsChooser.h"
   %}
   // std_vector.i for std::vector
   %include "std_vector.i"
   %template() std::vector<int>;
   // stl.i for std::pair
   %include "stl.i"
   %template() std::pair<unsigned long, unsigned long>;
   %include "myParentsChooser.h"

`Download myParentsChooser.i <myParentsChooser.i>`_

The exact procedure to generate and compile a wrapper file varies from system to
system, and from compiler to compiler. Fortunately, the standard Python module
setup process supports SWIG. All you need to do is to write a Python
``setup.py`` file and let the ``distutil`` module of Python handle all the
details for you. A typical ``setup.py`` file is demonstrated in Example
:ref:`parentsChooserSetup <parentsChooserSetup>`.

.. _parentsChooserSetup:

**Example**: *Building and installing the myParentsChooser module*

::

   from distutils.core import setup, Extension
   import sys
   # Under linux/gcc, lib stdc++ is needed for C++ based extension.
   if sys.platform == 'linux2':
       libs = ['stdc++']
   else:
       libs = []
   setup(name = "myParentsChooser",
       description = "A sample parent chooser",
       py_modules = ['myParentsChooser'],  # will be generated by SWIG
       ext_modules = [
           Extension('_myParentsChooser',
               sources = ['myParentsChooser.i'],
               swig_opts = ['-O', '-shadow', '-c++', '-keyword',],
               include_dirs = ["."],
       )
     ]
   )

`Download setup.py <setup.py>`_

You parent chooser can now be compiled and installed using the standard Python
``setup.py`` commands such as

::

   python setup.py install

Please refer to the Python reference manual for other building and installation
options. Note that Python 2.4 and earlier do not support option swig_opts well
so you might have to pass these options using command

::

   python setup.py build_ext --swig-opts=-O -templatereduce \
       -shadow -c++ -keyword -nodefaultctor install

Example :ref:`parentChooseHeader <parentChooseHeader>` demonstrates how to use
such a C++ parents chooser in your simuPOP script. It uses the same Python
parent chooser interface as in :ref:`PyParentsChooser <PyParentsChooser>`, but
leaves all the (potentially) computationally intensive parts to the C++ level
``myParentsChooser`` object.

.. _cppParentChooser:

**Example**: *Implement a parent chooser in C++*

::

   import simuPOP as sim

   # The class myParentsChooser is defined in module myParentsChooser
   try:
       from myParentsChooser import myParentsChooser
   except ImportError:
       # if failed to import the C++ version, use a Python version
       import random
       class myParentsChooser:
           def __init__(self, maleIndexes, femaleIndexes):
               self.maleIndexes = maleIndexes
               self.femaleIndexes = femaleIndexes
           def chooseParents(self):
               return self.maleIndexes[random.randint(0, len(self.maleIndexes)-1)],\
                   self.femaleIndexes[random.randint(0, len(self.femaleIndexes)-1)]

   def parentsChooser(pop, sp):
       'How to call a C++ level parents chooser.'
       # create an object with needed information (such as x, y) ...
       pc = myParentsChooser(
           [x for x in range(pop.popSize()) if pop.individual(x).sex() == sim.MALE],
           [x for x in range(pop.popSize()) if pop.individual(x).sex() == sim.FEMALE])
       while True:
           # return indexes of parents repeatedly
           yield pc.chooseParents()

   pop = sim.Population(100, loci=1)
   simu.evolve(
       initOps=[
           sim.InitSex(),
           sim.InitGenotype(freq=[0.5, 0.5])
       ],
       matingScheme=sim.HomoMating(sim.PyParentsChooser(parentsChooser),
           sim.OffspringGenerator(ops=sim.MendelianGenoTransmitter())),
       gen = 100
   )

`Download cppParentChooser.py <cppParentChooser.py>`_


