Output statistics
=================

We first want to output total disease allele frequency of each locus. This is
easy because :class:`Stat`\ () operator can calculate allele frequency for us.
What we need to do is use a :class:`Stat`\ () operator to calculate allele
frequency and get the result from population variable ``alleleFreq``. Because
allele frequcies add up to one, we can get the total disease allele frequency
using the allele frequency of the wild type allele 0
(:math:`\sum_{i=1}^{\infty}f_{i}=1-f_{0}`). The actual code would look more or
less like this:

::

   Stat(alleleFreq=[0,1]),
   PyEval(r'"%.2f" % (1-alleleFreq[0][0])')

We are also interested in the effective number of alleles Reich2001a at a locus.
Because simuPOP does not provide an operator or function to calculate this
statistic, we will have to calculate it manually. Fortunately, this is not
difficult because effective number of alleles can be calculated from existing
allele frequencies, using formula

.. math::

      n_{e}=\left(\sum_{i=1}^{\infty}\left(\frac{f_{i}}{1-f_{0}}\right)^{2}\right)^{-1}

where :math:`f_{i}` is the allele frequency of disease allele :math:`i`.

A quick-and-dirty way to output :math:`n_{e}` at a locus (e.g. locus 0) can be:

::

   PyEval('1./sum([(alleleFreq[0][x]/(1-alleleFreq[0][0]))**2 for x in alleleFreq[0].keys() if x != 0])')

but this expression looks complicated and does not handle the case when
:math:`f_{0}=1`. A more robust method would involve the ``stmts`` parameter of
:class:`PyEval`, which will be evaluated before parameter ``expr``:

::

   PyEval(stmts='''if alleleFreq[0][0] == 1:
       ne = 0
   else:
       freq = [freq[0][x] for x in alleleFreq[0].keys() if x != 0]
       ne = 1./sum([(f/(1-alleleFreq[0][0])**2 for x in freq])
   ''', expr=r'"%.3f" % ne')

However, this piece of code does not look nice with the multi-line string, and
the operator is not really reusable (only valid for locus o). It makes sense to
define a function to calculate :math:`n_{e}` generally:

::

   def ne(pop, loci):
       ' calculate effective number of alleles at given loci'
       stat(pop, alleleFreq=loci)
       ne = {}
       for loc in loci:
           freq = [y for x,y in pop.dvars().alleleFreq[loc].iteritems() if x != 0]
           sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
           if sumFreq == 0:
               ne[loc] = 0
           else:
               ne[loc] = 1. / sum([(x/sumFreq)**2 for x in freq])
       # save the result to the population.
       pop.dvars().ne = ne
       return True

When it is needed to calculate effective number of alleles, a Python operator
that uses this function can be used. For example, operator

::

   PyOperator(func=ne, param=[0], step=5)
   PyEval(r'"%.3f" % ne[0]', step=5)

would calculate effective number of alleles at locus 0 and output it.

The biggest difference between :class:`PyEval` and :class:`PyOperator` is that
:class:`PyOperator` is no longer evaluated in the population's local namespace.
You will have to get the variables explicitly using the ``pop.dvars()``
function, and the results have to be explicitly saved to the population's local
namespace.

The final implementation, as a way to demonstrate how to define a new statistics
that hides all the details, defines a new operator by inheriting a class from
:class:`PyOperator`. The resulting operator could be used as a regular operator
(e.g., ``ne(loci=[0])``). A function ``Ne`` is also defined as the function form
of this operator. The code is listed in Example :ref:`reichstat <reichstat>`

.. _reichstat:

**Example**: *A customized operator to calculate effective number of alleles*

::

   >>> import simuPOP as sim
   >>> class ne(sim.PyOperator):
   ...     '''Define an operator that calculates effective number of
   ...     alleles at given loci. The result is saved in a population
   ...     variable ne.
   ...     '''
   ...     def __init__(self, loci, *args, **kwargs):
   ...         self.loci = loci
   ...         sim.PyOperator.__init__(self, func=self.calcNe, *args, **kwargs)
   ...     #
   ...     def calcNe(self, pop):
   ...         sim.stat(pop, alleleFreq=self.loci)
   ...         ne = {}
   ...         for loc in self.loci:
   ...             freq = pop.dvars().alleleFreq[loc]
   ...             sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
   ...             if sumFreq == 0:
   ...                 ne[loc] = 0
   ...             else:
   ...                 ne[loc] = 1. / sum([(freq[x]/sumFreq)**2 for x in list(freq.keys()) if x != 0])
   ...         # save the result to the sim.Population.
   ...         pop.dvars().ne = ne
   ...         return True
   ... 
   >>> def Ne(pop, loci):
   ...     '''Function form of operator ne'''
   ...     ne(loci).apply(pop)
   ...     return pop.dvars().ne
   ... 
   >>> pop = sim.Population(100, loci=[10])
   >>> sim.initGenotype(pop, freq=[.2] * 5)
   >>> print(Ne(pop, loci=[2, 4]))
   {2: 3.9565470135154768, 4: 3.948841408365935}

   now exiting runScriptInteractively...

`Download reichstat.py <reichstat.py>`_


