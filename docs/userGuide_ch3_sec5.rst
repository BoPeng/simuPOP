Random number generator \*
==========================

.. index::
   single: moduleInfo
   single: setRNG

When simuPOP is loaded, it creates a default random number generator
(:class:`RNG`) of type ``mt19937`` for each thread. It uses a random seed for
the first RNG and uses seeds derived from the first seed to initialize RNGs for
other threads. The seed is drawn from a system random number generator that
guarantees random seeds for all instances of simuPOP even if they are
initialized at the same time. After simuPOP is loaded, you can reset this system
RNG with a different random number generator (c.f. :func:`moduleInfo`\
``['availableRNGs']``) or use a specified seed using function , ``setRNG(name,
seed)``.

:func:`getRNG`\ ``.seed()`` returns the seed of the simuPOP random number
generator. It can be used to replay your simulation if :func:`getRNG`\ () is
your only source of random number generator. If you also use the Python
``random`` module, it is a good practise to set its seed using
``random.seed(getRNG().seed())``. Example :ref:`randomSeed <randomSeed>`
demonstrates how to use these functions to replay an evolutionary process.
simuPOP uses a single seed to initialize multiple random number generators used
for different threads (seeds for other threads are determined from the first
seed) so you only need to save the head seed (:func:`getRNG`\ ``.seed()``)

.. _randomSeed:

**Example**: *Use saved random seed to replay an evolutionary process*

::

   >>> import simuPOP as sim
   >>> import random
   >>> def simulate():
   ...     pop = sim.Population(1000, loci=10, infoFields='age')
   ...     pop.evolve(
   ...         initOps=[
   ...             sim.InitSex(),
   ...             sim.InitGenotype(freq=[0.5, 0.5]),
   ...             sim.InitInfo(lambda: random.randint(0, 10), infoFields='age')
   ...         ],
   ...         matingScheme=sim.RandomMating(),
   ...         finalOps=sim.Stat(alleleFreq=0),
   ...         gen=100
   ...     )
   ...     return pop.dvars().alleleFreq[0][0]
   ... 
   >>> seed = sim.getRNG().seed()
   >>> random.seed(seed)
   >>> print('%.4f' % simulate())
   0.5780
   >>> # will yield different result
   >>> print('%.4f' % simulate())
   0.6355
   >>> sim.setRNG(seed=seed)
   >>> random.seed(seed)
   >>> # will yield identical result because the same seeds are used
   >>> print('%.4f' % simulate())
   0.5780

   now exiting runScriptInteractively...

`Download randomSeed.py <randomSeed.py>`_


