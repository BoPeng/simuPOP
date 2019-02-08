Initialize and evolve the population
====================================

With appropriate operators to perform mutation, selection and output statistics,
it is relatively easy to write a simulator to perform a simulation. This
simulator would create a population, initialize alleles with an initial allic
spectrum, and then evolve it according to specified demographic model. During
the evolution, mutation and selection will be applied, statistics will be
calculated and outputed.

.. _reichEvolve:

**Example**: *Evolve a population subject to mutation and selection*

::

   >>> import simuPOP as sim
   >>> 
   >>> 
   >>> def simulate(model, N0, N1, G0, G1, spec, s, mu, k):
   ...     '''Evolve a sim.Population using given demographic model
   ...     and observe the evolution of its allelic spectrum.
   ...     model: type of demographic model.
   ...     N0, N1, G0, G1: parameters of demographic model.
   ...     spec: initial allelic spectrum, should be a list of allele
   ...         frequencies for each allele.
   ...     s: selection pressure.
   ...     mu: mutation rate.
   ...     k: k for the k-allele model
   ...     '''
   ...     demo_func = demo_model(model, N0, N1, G0, G1)
   ...     pop = sim.Population(size=demo_func(0), loci=1, infoFields='fitness')
   ...     pop.evolve(
   ...         initOps=[
   ...             sim.InitSex(),
   ...             sim.InitGenotype(freq=spec, loci=0)
   ...         ],
   ...         matingScheme=sim.RandomMating(subPopSize=demo_func),
   ...         postOps=[
   ...             sim.KAlleleMutator(k=k, rates=mu),
   ...             sim.MaSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
   ...             ne(loci=[0], step=100),
   ...             sim.PyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
   ...                 step=100),
   ...         ],
   ...         gen = G0 + G1
   ...     )
   ... 
   >>> simulate('instant', 1000, 10000, 500, 500, [0.9]+[0.02]*5, 0.01, 1e-4, 200)
   0: 0.09	4.91
   100: 0.12	2.63
   200: 0.09	1.22
   300: 0.02	2.85
   400: 0.02	2.12
   500: 0.05	1.02
   600: 0.06	1.51
   700: 0.08	1.58
   800: 0.09	1.80
   900: 0.08	1.79

   now exiting runScriptInteractively...

`Download reichEvolve.py <reichEvolve.py>`_


