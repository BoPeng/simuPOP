Tracing allelic lineage \*
==========================

Lineage of alleles consists of information such as the distribution of alleles
(how many people carry this allele, and the relationship between carriers) and
age of alleles (when the alleles were introduced to the population). These
information are important for the study of evolutionary history of mutants. They
are not readily available for normal simulations, and even if you can track the
generations when mutants are introduced, alleles in the present generation that
are of the same type (Identity by Stat, IBS) do not necessarily have the same
ancestral origin (Identity by Decent, IBD).

The lineage modules of simuPOP provides facilities to track allelic lineage.
More specifically,

* Each allele is associated with an integer number (an allelic lineage) that
  identifies the origin, or the source of the allele.

* The lineage of each allele is transmitted along with the allele during
  evolution. New alleles will be introduced with their own lineage, even if they
  share the same states with existing alleles.

* Origin of alleles can be accessed using member functions of the
  :class:`Individual` and :class:`Population` classes.

Example :ref:`geneticContribution <geneticContribution>` demonstrates how to
determine the contribution of genetic information from each ancestor. For this
simulation, the alleles of each ancestor are associated with individual-specific
numbers. During evolution, some alleles might get lost, some are copied, and
pieces of chromosomes are mixed due to genetic recombination. At the end of
simulation, the average number of 'contributors' of genetic information to each
individual is calculated, as well as the percent of genetic information from
each ancestor. Although this particular simulation can be mimicked using pure-
genotype simulations by using special alleles for each ancestor, the combined
information regarding the state and origin of each allele will be very useful
for genetic studies that involve IBD and IBS.

.. _geneticContribution:

**Example**: *Contribution of genetic information from ancestors*

::

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='lineage', quiet=True)
   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, loci=[10]*4)
   >>> 
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.25]*4),
   ...         sim.InitLineage(range(1000), mode=sim.PER_INDIVIDUAL),
   ...     ],
   ...     matingScheme=sim.RandomMating(ops=sim.Recombinator(rates=0.001)),
   ...     gen = 100
   ... )
   100
   >>> # average number of 'contributors'
   >>> num_contributors = [len(set(ind.lineage())) for ind in pop.individuals()]
   >>> print('Average number of contributors is %.2f' % (sum(num_contributors) / float(pop.popSize())))
   Average number of contributors is 13.98
   >>> # percent of genetic information from each ancestor (baseline is 1/1000)
   >>> lineage = pop.lineage()
   >>> lin_perc = [lineage.count(x)/float(len(lineage)) for x in range(1000)]
   >>> # how many of ancestors do not have any allele left?
   >>> print('Number of ancestors with no allele left: %d' % lin_perc.count(0.))
   Number of ancestors with no allele left: 817
   >>> # top five contributors
   >>> lin_perc.sort()
   >>> lin_perc.reverse()
   >>> print('Top contributors (started with 0.001): %.5f %.5f %.5f' % (lin_perc[0], lin_perc[1], lin_perc[2]))
   Top contributors (started with 0.001): 0.03474 0.03058 0.02475

   now exiting runScriptInteractively...

`Download geneticContribution.py <geneticContribution.py>`_

Example :ref:`geneticContribution <geneticContribution>` uses operator
:class:`InitLineage` to explictly assign lineage to alleles of each individual.
You can also track the fate of finer genetic pieces by assigning different
lineage values to chromosomes, or each loci using different ``mode``. This
operator can also assign lineage of alleles to an ID stored in an information
field, which is usually ``ind_id``, a field used by operators such as
:class:`IdTagger` and :class:`PedigreeTagger` to assign and trace the pedigree
(parentship) information during evolution. More interesting, when such a field
is present, mutation operators will assign the IDs of recipients of mutants as
the lineage of these mutants. This makes it possible to track the origin of
mutants. Moreover, when a mode ``FROM_INFO_SIGNED`` is used, additional ploidy
information will be tagged to lineage values (negative values for mutants on the
second homologous copy of chromosomes) so that you can track the inheritance of
haplotypes.

To make use of these features, it is important to assign IDs to individuals
before these operators are applied. Example :ref:`ageOfMutants <ageOfMutants>`
demonstrates how to use the lineage information to determine the age of mutants.
This example evolves a constant population of size 10,000. An :class:`IdTagger`
is used before :class:`InitGenotype` so individual IDs will be assigned as
allelic lineages. Because all offspring get their own IDs during evolution, the
IDs of individuals are assigned to mutants as their lineages, and can be used to
determine the age of these mutants. This is pretty easy to do in this example
because of constant population size. For more complex demographic models, you
might have to record the minimal and maximum IDs of each generation in order to
determine the age of mutants.

.. _ageOfMutants:

**Example**: *Distribution of age of mutants*

::

   >>> import simuOpt
   >>> simuOpt.setOptions(alleleType='lineage', quiet=True)
   >>> import simuPOP as sim
   >>> pop = sim.Population(size=10000, loci=[10]*10, infoFields='ind_id')
   >>> # just to make sure IDs starts from 1
   >>> sim.IdTagger().reset(1)
   >>> pop.evolve(
   ...     initOps = [
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.2, 0.3, 0.4, 0.1]),
   ...         sim.IdTagger(),
   ...         sim.InitLineage(mode=sim.FROM_INFO),
   ...     ],
   ...     # an extremely high mutation rate, just for demonstration
   ...     preOps = sim.AcgtMutator(rate=0.01, model='JC69'),
   ...     matingScheme=sim.RandomMating(
   ...         ops=[
   ...             sim.IdTagger(),
   ...             sim.MendelianGenoTransmitter(),
   ...         ]
   ...     ),
   ...     gen = 10
   ... )
   10
   >>> lin = pop.lineage()
   >>> # Number of alleles from each generation
   >>> for gen in range(10):
   ...     id_start = gen*10000 + 1
   ...     id_end = (gen+1)*10000
   ...     num_mut = len([x for x in lin if x >= id_start and x <= id_end])
   ...     print('Gen %d: %5.2f %%' % (gen, num_mut / (2*10000*100.) * 100))
   ... 
   Gen 0: 93.40 %
   Gen 1:  0.72 %
   Gen 2:  0.71 %
   Gen 3:  0.70 %
   Gen 4:  0.74 %
   Gen 5:  0.76 %
   Gen 6:  0.73 %
   Gen 7:  0.74 %
   Gen 8:  0.75 %
   Gen 9:  0.75 %

   now exiting runScriptInteractively...

`Download ageOfMutants.py <ageOfMutants.py>`_


