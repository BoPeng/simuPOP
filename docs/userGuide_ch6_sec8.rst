Simulation of mitochondrial DNAs (mtDNAs) \*
============================================

Mitochondrial DNAs resides in human mitochondrion. A zygote inherits its
organelles from the cytoplasm of the egg, and thus organelle inheritance is
generally maternal. Whereas there is only one copy of a nuclear chromosome per
gamete, there are man copies of an organellar chromosome, forming a population
of identical organelle chromosomes that is transmitted to the offspring through
the egg. Because these organellar chromosomes are identical, they are modelled
in simuPOP as a single chromosome with type ``MITOCHONDRIAL``. In order to
simulate mitochondrial DNAs, it is important to remember:

* :class:`MendelianGenoTransmitter` and :class:`Recombinator` do not handle
  mitochondrial DNAs so you will have to explicitly use
  :class:`MitochondrialGenoTransmitter` to transmit the mitochondrial DNAs from
  mother to offspring. Note that :class:`CloneGenoTransmitter` is a special
  transmitter that will copy everything including sex, information fields to
  offspring.

* The :class:`Stat` operator recognizes this chromosome type and will report
  allele, haplotype, and genotype counts, and other statistics correctly, although
  some diploid-specific statistics are not applicable.

* Natural selections on mtDNAs is usually performed using operator
  :class:`MapSelector` where single alleles are assigned a fitness value. Operator
  :class:`MaSelector` assumes two alleles and is not applicable.

Example :ref:`mitochondrial <mitochondrial>` demonstrates the use of a
:class:`Recombinator` to recombine an autosome and two sex chromosomes, and a
:class:`MitochondrialGenoTransmitter` to transmit mitochondrial chromosomes.
Natural selection is applied to allele 3 at the 3rd locus on the mitochondrial
DNA, whose frequency in the population decreases as a result.

.. _mitochondrial:

**Example**: *Transmission of mitochondrial chromosomes*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(1000, loci=[5]*4,
   ...     # one autosome, two sex chromosomes, and one mitochondrial chromosomes
   ...     chromTypes=[sim.AUTOSOME, sim.CHROMOSOME_X, sim.CHROMOSOME_Y, sim.MITOCHONDRIAL],
   ...     infoFields=['fitness'])
   >>> pop.evolve(
   ...     initOps=[
   ...         sim.InitSex(),
   ...         sim.InitGenotype(freq=[0.25]*4)
   ...     ],
   ...     preOps=[
   ...         sim.MapSelector(loci=17, fitness={(0,): 1, (1,): 1, (2,): 1, (3,): 0.4})
   ...     ],
   ...     matingScheme=sim.RandomMating(ops= [
   ...         sim.Recombinator(rates=0.1),
   ...         sim.MitochondrialGenoTransmitter(),
   ...     ]),
   ...     postOps=[
   ...         sim.Stat(alleleFreq=17, step=10),
   ...         sim.PyEval(r'"%.2f %.2f %.2f %.2f\n" % (alleleNum[17][0],'
   ...             'alleleNum[17][1], alleleNum[17][2], alleleNum[17][3])', step=10),
   ...     ],
   ...     gen = 100
   ... )
   1288.00 273.00 325.00 114.00
   1384.00 245.00 371.00 0.00
   1492.00 138.00 370.00 0.00
   1461.00 69.00 470.00 0.00
   1449.00 65.00 486.00 0.00
   1536.00 17.00 447.00 0.00
   1624.00 7.00 369.00 0.00
   1538.00 0.00 462.00 0.00
   1619.00 0.00 381.00 0.00
   1623.00 0.00 377.00 0.00
   100

   now exiting runScriptInteractively...

`Download mitochondrial.py <mitochondrial.py>`_

You might wonder how a mutation can change the allele of all organelles in the
mitochondrion. This is generally believed to be done through natural drift
during cytoplasmic segreagation, which is not a mitotic process because it takes
place in dividing asexual cells. Because only one mitochondrial chromosome is
allowed in simuPOP, you will have to use customized chromosome types if you
would like to simulate this process. Fortunately, operator
:class:`MitochondrialGenoTransmitter` can select random organelles from multiple
customized chromosomes, if no chromosome of type ``MITOCHONDRIAL`` is present.

Example :ref:`mtDNA_evolve <mtDNA_evolve>` demonstrates the fixation of mutant
in cells with multiple organelles. Althogh mutations are introduced to only one
of the organelles, after a number of cell divisions, the majority of the cells
now have only one type of allele. This example uses a :class:`RandomSelection`
mating scheme to select cells randomly from the parental population. Because no
sexual reproduction is involved, :class:`MitochondrialGenoTransmitter` passes
the parental genotype to offspring regardless of sex of parent. This example
also demonstrates a disadvantage of using customized chromosomes in that you
will have to calculate statistics by yourself because only you know the meaning
of these chromosomes. In this example, a function is written to count the number
of mutants in each cell (individual), and summarize the number of cells with 0,
1, 2, 3, 4, and 5 copies of the mutant.

.. _mtDNA_evolve:

**Example**: *Evolution of multiple organelles in mitochondrion*

::

   >>> import simuPOP as sim
   >>> 
   >>> def alleleCount(pop):
   ...     summary = [0]* 6
   ...     for ind in pop.individuals():
   ...         geno = ind.genotype(ploidy=0)
   ...         summary[geno[0] + geno[2] + geno[4] + geno[6] + geno[8]] += 1
   ...     print('%d %s' % (pop.dvars().gen, summary))
   ...     return True
   ... 
   >>> pop = sim.Population(1000, loci=[2]*5, chromTypes=[sim.CUSTOMIZED]*5)
   >>> pop.evolve(
   ...     # every one has miDNAs 10, 00, 00, 00, 00
   ...     initOps=[
   ...         sim.InitGenotype(haplotypes=[[1]+[0]*9]),
   ...     ],
   ...     # random select cells for cytoplasmic segregation
   ...     matingScheme=sim.RandomSelection(ops= [
   ...         sim.MitochondrialGenoTransmitter(),
   ...     ]),
   ...     postOps=sim.PyOperator(func=alleleCount, step=10),
   ...     gen = 51
   ... )
   0 [333, 408, 219, 38, 2, 0]
   10 [806, 16, 14, 16, 11, 137]
   20 [816, 1, 1, 3, 0, 179]
   30 [833, 0, 0, 0, 0, 167]
   40 [805, 0, 0, 0, 0, 195]
   50 [849, 0, 0, 0, 0, 151]
   51

   now exiting runScriptInteractively...

`Download mtDNA_evolve.py <mtDNA_evolve.py>`_


