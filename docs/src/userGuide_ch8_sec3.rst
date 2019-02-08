Mutation and selection models
=============================

The thoretical model empolyees an infinite allele model where there is a single
wild type allele and an infinite number of disease alleles. Each mutation would
introduce a new disease allele and there is no back mutation (mutation from
disease allele to wild type allele).

This mutation model can be mimicked by a :math:`k`\ -allele model with
resaonably large :math:`k`. We initialize all alleles to 0 which is the wild
type (:math:`A`) and all other alleles are considered as disease alleles
(:math:`a`). Because an allele in a :math:`k-`\ allele mutation model can mutate
to any other allele with equal probability, :math:`P\left(A\rightarrow
a\right)\gg P\left(a\rightarrow A\right)` since there are many more disease
alleles than the wild type allele. If we choose a smaller :math:`k` (e.g.
:math:`k=20`), recurrent and back mutations can on longer be ignored but it
would be interesting to simulate such cases because they are more realistic than
the infinite allele model in some cases.

A :math:`k`\ -allele model can be simulated using the :class:`KAlleleMutator`
operator which accepts a mutation rate and a maximum allelic state as
parameters.  ::

   KAlleleMutator(k=k, rates=mu)

Because there are many possible disease alleles, a multi-allelic selector
(:class:`MaSelector`) could be used to select against the disease alleles. This
operator accept a single or a list of wild type alleles (``[0]`` in this case)
and treat all other alleles as disease alleles. A penetrance table is needed
which specified the fitness of each individual when they have 0, 1 or 2 disease
alleles respectively. In this example, we assume a recessive model in which only
genotype :math:`aa` causes genetic disadvantages. If we assume a selection
pressure parameter :math:`s`, the operator to use is  ::

   MaSelector(loci=0, wildtype=0, penetrance=[1, 1, 1-s])

Note that the use of this selector requires a population information field
``fitness``.

This example uses a single-locus selection model but the complete script allows
the use of different kinds of multi-locus selection model. If we assume a
multiplicative multi-locus selection model where fitness values at different
loci are combined (multiplied), a multi-locus selection model
(:class:`MlSelector`) could be used as follows:

::

   MlSelector([
       MaSelector(loci=loc1, fitness=[1,1,1-s1], wildtype=0),
       MaSelector(loci=loc2, fitness=[1,1,1-s2], wildtype=0)],
       mode=MULTIPLICATIVE
   )

These multi-locus model treat disease alleles at different loci more or less
independently. If more complex multi-locus models (e.g. models involve gene -
gene and/or gene - interaction) are involved, a multi-locus selector that uses a
multi-locus penetrance table could be used.


