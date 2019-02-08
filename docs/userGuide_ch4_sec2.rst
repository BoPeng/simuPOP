Individual
==========

individuals are building blocks of a population. An individual object cannot be
created independently, but references to inidividuals can be retrieved using
member functions of a population object.


Access individual genotype
--------------------------

From a user's point of view, genotypes of an individual are stored sequentially
and can be accessed locus by locus, or in batch. The alleles are arranged by
position, chromosome and ploidy. That is to say, the first allele on the first
chromosome of the first homologous set is followed by alleles at other loci on
the same chromosome, then markers on the second and later chromosomes, followed
by alleles on the second homologous set of the chromosomes for a diploid
individual. A consequence of this memory layout is that alleles at the same
locus of a non-haploid individual are separated by ``Individual.totNumLoci()``
loci. The memory layout of a diploid individual with two chromosomes is
illustrated in Figure :ref:`fig_genotype_layout <fig_genotype_layout>`.

**Figure**: *Memory layout of individual genotype*

.. _fig_genotype_layout:

.. figure:: /Users/bpeng1/simuPOP/simuPOP/doc/figures/genotype.png
   :width: 680


simuPOP provides several functions to read/write individual genotype. For
example, :meth:`Individual.allele`\ () and :meth:`Individual.setAllele`\ () can
be used to read and write single alleles. You could also access alleles in batch
mode using functions :meth:`Individual.genotype`\ () and
:meth:`Individual.setGenotype`\ (). It is worth noting that, instead of copying
genotypes of an individual to a Python tuple or list, the return value of
function ``genotype([p, [ch]])`` is a special python carray object that reflects
the underlying genotypes. This object behaves like a regular Python list except
that the underlying genotype will be changed if elements of this object are
changed. Only ``count(x)`` and\ ``index(x, [start, [stop]])`` member functions
can be used, but all comparison, assignment and slice operations are allowed. If
you would like to copy the content of this ``carray`` to a Python list, use the
``list`` function. Example :ref:`individualGenotype <individualGenotype>`
demonstrates the use of these functions.

.. _individualGenotype:

**Example**: *Access individual genotype*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([2, 1], loci=[2, 5])
   >>> for ind in pop.individuals(1):
   ...     for marker in range(pop.totNumLoci()):
   ...         ind.setAllele(marker % 2, marker, 0)
   ...         ind.setAllele(marker % 2, marker, 1)
   ...         print('%d %d ' % (ind.allele(marker, 0), ind.allele(marker, 1)))
   ... 
   0 0 
   1 1 
   0 0 
   1 1 
   0 0 
   1 1 
   0 0 
   >>> ind = pop.individual(1)
   >>> geno = ind.genotype(1)      # the second homologous copy
   >>> geno
   [0, 0, 0, 0, 0, 0, 0]
   >>> geno[2] = 3
   >>> ind.genotype(1)
   [0, 0, 3, 0, 0, 0, 0]
   >>> geno[2:4] = [3, 4]          # direct modification of the underlying genotype
   >>> ind.genotype(1)
   [0, 0, 3, 4, 0, 0, 0]
   >>> # set genotype (genotype, ploidy, chrom)
   >>> ind.setGenotype([2, 1], 1, 1)
   >>> geno
   [0, 0, 2, 1, 2, 1, 2]
   >>> #
   >>> geno.count(1)           # count
   2
   >>> geno.index(2)           # index 
   2
   >>> ind.setAllele(5, 3)    # change underlying genotype using setAllele
   >>> print(geno)              # geno is change
   [0, 0, 2, 1, 2, 1, 2]
   >>> print(geno)             # but not geno
   [0, 0, 2, 1, 2, 1, 2]
   >>> geno[2:5] = 4           # can use regular Python slice operation
   >>> print(ind.genotype())
   [0, 0, 0, 5, 0, 0, 0, 0, 0, 4, 4, 4, 1, 2]

   now exiting runScriptInteractively...

`Download individualGenotype.py <individualGenotype.py>`_

The same object will also be returned by function :meth:`Population.genotype`\
().


individual sex, affection status and information fields
-------------------------------------------------------

In addition to structural information shared by all individuals in a population,
the individual class provides member functions to get and set *genotype*, *sex*,
*affection status* and *information fields* of an individual. Example
:ref:`individuals <individuals>` demonstrates how to access and modify
individual sex, affection status and information fields. Note that **information
fields can be accessed as attributes of individuals**. For example,
``ind.info('father_idx')`` is equivalent to ``ind.father_idx`` and
``ind.setInfo(35, 'age')`` is equivalent to ``ind.age = 35``.

.. _individuals:

**Example**: *Access Individual properties*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population([5, 4], loci=[2, 5], infoFields='x')
   >>> # get an individual
   >>> ind = pop.individual(3)
   >>> ind.ploidy()            # access to genotypic structure
   2
   >>> ind.numChrom()
   2
   >>> ind.affected()
   False
   >>> ind.setAffected(True)   # access affection sim.status,
   >>> ind.sex()               # sex,
   1
   >>> ind.setInfo(4, 'x')     # and information fields
   >>> ind.x = 5               # the same as ind.setInfo(4, 'x')
   >>> ind.info('x')           # get information field x
   5.0
   >>> ind.x                   # the same as ind.info('x')
   5.0

   now exiting runScriptInteractively...

`Download individual.py <individual.py>`_


