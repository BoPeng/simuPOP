Option handling
===============

Everything seems to be perfect until you need to

#. Run more simulations with different parameters such as initial population
   size and mutaion rate. This requires the script to get its parameters from
   command line (or a configuration file) and executes in batch mode, perhaps on a
   cluster system.

#. Allow users who are not familiar with the script to run it. This would better
   be achieved by a graphical user interface.

#. Allow other Python scripts to import your script and run the simulation
   function directly.

Although a number of Python modules such as ``getopt`` are available, the
simuPOP :mod:`simuOpt` module is especially designed to allow a simuPOP script
to be run both in batch and in GUI mode, in standard and optimized mode. Example
:ref:`reich <reich>` makes use of this module.

.. _reich:

**Example**: *A complete simulation script*

::

   #!/usr/bin/env python
   #
   # Author:  Bo Peng
   # Purpose: A real world example for simuPOP user's guide.
   #
   '''
   Simulation the evolution of allelic spectra (number and frequencies
   of alleles at a locus), under the influence of sim.population expansion,
   mutation, and natural selection.
   '''
   import simuOpt
   simuOpt.setOptions(quiet=True, alleleType='long')
   import simuPOP as sim
   import sys, types, os, math
   options = [
       {'name': 'demo',
        'default': 'instant',
        'label': 'Population growth model',
        'description': 'How does a sim.Population grow from N0 to N1.',
        'type': ('chooseOneOf', ['instant', 'exponential']),
       },
       {'name': 'N0',
        'default': 10000,
        'label': 'Initial sim.population size',
        'type': 'integer',
        'description': '''Initial sim.population size. This size will be maintained
                   till the end of burnin stage''',
        'validator': simuOpt.valueGT(0)
       },
       {'name': 'N1',
        'default': 100000,
        'label': 'Final sim.population size',
        'type': 'integer',
        'description': 'Ending sim.population size (after sim.population expansion)',
        'validator': simuOpt.valueGT(0)
       }, 
       {'name': 'G0',
        'default': 500,
        'label': 'Length of burn-in stage',
        'type': 'integer',
        'description': 'Number of generations of the burn in stage.',
        'validator': simuOpt.valueGT(0)
       },
       {'name': 'G1',
        'default': 1000,
        'label': 'Length of expansion stage',
        'type': 'integer',
        'description': 'Number of geneartions of the sim.population expansion stage',
        'validator': simuOpt.valueGT(0)
       },
       {'name': 'spec',
        'default': [0.9] + [0.02]*5,
        'label': 'Initial allelic spectrum',
        'type': 'numbers',
        'description': '''Initial allelic spectrum, should be a list of allele
               frequencies, for allele 0, 1, 2, ... respectively.''',
        'validator': simuOpt.valueListOf(simuOpt.valueBetween(0, 1)),
       },
       {'name': 's',
        'default': 0.01,
        'label': 'Selection pressure',
        'type': 'number',
        'description': '''Selection coefficient for homozygtes (aa) genotype.
               A recessive selection model is used so the fitness values of
               genotypes AA, Aa and aa are 1, 1 and 1-s respectively.''',
        'validator': simuOpt.valueGT(-1),
       },
       {'name': 'mu',
        'default': 1e-4,
        'label': 'Mutation rate',
        'type': 'number',
        'description': 'Mutation rate of a k-allele mutation model',
        'validator': simuOpt.valueBetween(0, 1),
       },
       {'name': 'k',
        'default': 200,
        'label': 'Maximum allelic state',
        'type': 'integer',
        'description': 'Maximum allelic state for a k-allele mutation model',
        'validator': simuOpt.valueGT(1),
       },
   ]

   def demo_model(type, N0=1000, N1=100000, G0=500, G1=500):
       '''Return a demographic function 
       type: linear or exponential
       N0:   Initial sim.population size.
       N1:   Ending sim.population size.
       G0:   Length of burn-in stage.
       G1:   Length of sim.population expansion stage.
       '''
       rate = (math.log(N1) - math.log(N0))/G1
       def ins_expansion(gen):
           if gen < G0:
               return N0
           else:
               return N1

       def exp_expansion(gen):
           if gen < G0:
               return N0
           else:            
               return int(N0 * math.exp((gen - G0) * rate))

       if type == 'instant':
           return ins_expansion
       elif type == 'exponential':
           return exp_expansion

   class ne(sim.PyOperator):
       '''Define an operator that calculates effective number of
       alleles at given loci. The result is saved in a population
       variable ne.
       '''
       def __init__(self, loci, *args, **kwargs):
           self.loci = loci
           sim.PyOperator.__init__(self, func=self.calcNe, *args, **kwargs)

       def calcNe(self, pop):
           sim.stat(pop, alleleFreq=self.loci)
           ne = {}
           for loc in self.loci:
               freq = pop.dvars().alleleFreq[loc]
               sumFreq = 1 - pop.dvars().alleleFreq[loc][0]
               if sumFreq == 0:
                   ne[loc] = 0
               else:
                   ne[loc] = 1. / sum([(freq[x]/sumFreq)**2 for x in list(freq.keys()) if x != 0])
           # save the result to the sim.Population.
           pop.dvars().ne = ne
           return True

   def simuCDCV(model, N0, N1, G0, G1, spec, s, mu, k):
       '''Evolve a sim.Population using given demographic model
       and observe the evolution of its allelic spectrum.
       model: type of demographic model.
       N0, N1, G0, G1: parameters of demographic model.
       spec: initial allelic spectrum, should be a list of allele
           frequencies for each allele.
       s: selection pressure.
       mu: mutation rate.
       k: maximum allele for the k-allele model
       '''
       demo_func = demo_model(model, N0, N1, G0, G1)
       print(demo_func(0))
       pop = sim.Population(size=demo_func(0), loci=1, infoFields='fitness')
       pop.evolve(
           initOps=[
               sim.InitSex(),
               sim.InitGenotype(freq=spec, loci=0)
           ],
           matingScheme=sim.RandomMating(subPopSize=demo_func),
           postOps=[
               sim.KAlleleMutator(rates=mu, k=k),
               sim.MaSelector(loci=0, fitness=[1, 1, 1 - s], wildtype=0),
               ne(loci=(0,), step=100),
               sim.PyEval(r'"%d: %.2f\t%.2f\n" % (gen, 1 - alleleFreq[0][0], ne[0])',
                   step=100),
           ],
           gen = G0 + G1
       )
       return pop

   if __name__ == '__main__':
       # get parameters
       par = simuOpt.Params(options, __doc__)
       if not par.getParam():
           sys.exit(1)

       if not sum(par.spec) == 1:
           print('Initial allelic spectrum should add up to 1.')
           sys.exit(1)
       # save user input to a configuration file
       par.saveConfig('simuCDCV.cfg')
       #
       simuCDCV(*par.asList())


`Download simuCDCV.py <simuCDCV.py>`_

Example :ref:`reich <reich>` uses a programming style that is used by almost all
simuPOP scripts. I highly recommend this style because it makes your script
seld-documentary and work well under a variety of environments. A script written
in this style follows the following order:

#. First comment block

   The first line of the script should always be   ::

      #!/usr/bin/env python

   This line tells a Unix shell which program should be used to process this script
   if the script to set to be executable. This line is ignored under windows. It is
   customary to put author and date information at the top of a script as Python
   comments.

#. Module doc string

   The first string in a script is the module docstring, which can be referred by
   variable ``__doc__`` in the script. It is a good idea to describe what this
   script does in detail here. As you will see later, this docstring will be used
   in the ``simuOpt.getParam()`` function and be outputed in the usage information
   of the script.

#. Loading simuPOP and other Python modules

   simuPOP and other modules are usually imported after module docstring. This is
   where you specify which simuPOP module to use. Although a number of parameters
   could be used, usually only ``alleleType`` is specified because other parameters
   such as ``gui`` and ``optimized`` should better be controlled from command line.

#. Parameter description list

   A list of parameter description dictionaries are given here. This list specifies
   what parameters will be used in this script and describes the type, default
   value, name of command line option, label of the parameter in the parameter
   input dialog in detail. Although some directionary items can be ignored, it is a
   good practice to give detailed information about each parameter here.

#. Helper functions and classes

   Helper functions and classes are given before the main simulation function.

#. Main simulation function

   The main simulation function preforms the main functionality of the whole
   script. It is written as a function so that it can be imported and executed by
   another script. The parameter processing part of the script would be ignored in
   this case.

#. Script execution part conditioned by ``__name__ == '__main__'``

   The execution part of a script should always be inside of a ``if __name__ ==
   '__main__'`` block so that the script will not be executed when it is imported
   by another script. The first few lines of this execution block are almost always

   ::

      par = simuOpt.Params(options, __doc__)
      if not par.getParam():
          sys.exit(1)

   which creates a simuOpt object and tries to get parameters from command line
   option, a configuration file, a parameter input dialog or interactive user
   input, depending on how this script is executed. Optionally, you can use

   ::

      par.saveConfig('file.cfg')

   to save the current configuration to a file so that the same parameters could be
   retrieved later using parameter ``--config file.cfg``.

   After simply parameter validation, the main simulation function can be called.
   This example uses ``simuCDCV(*par.asList())`` because the parameter list in the
   ``par`` object match the parameter list of function ``simuCDCV`` exactly. If
   there are a large number of parameters, it may be better to pass the
   :mod:`simuOpt` object directly in the main simulation function.

The script written in this style could be executed in a number of ways.

#. If a user executes the script directly, a Tkinter or wxPython dialog will be
   displayed for users to input parameters. This parameter is shown in Figure
   :ref:`fig_simuCDCV_dialog <fig_simuCDCV_dialog>`.

   **Figure**: *Parameter input dialog of the simuCDCV script*

   .. _fig_simuCDCV_dialog:

   .. figure:: /Users/bpeng1/simuPOP/simuPOP/doc/figures/simuCDCV.png
      :width: 5in
   

#. The help message of this script could be displayed using the Help button of
   the parameter input dialog, or using command ``simuCDCV.py -h``.

#. Using parameter ``--gui=False``, the script will be run in batch mode. You
   can specify parameters using  ::

      simuCDCV.py --gui=False --config file.cfg

   if a parameter file is available, or use command line options such as  ::

      simuCDCV.py --gui=False --demo='instant' --N0=10000 --N1=100000 \
          --G0=500 --G1=500 --spec='[0.9]+[0.02]*5' --s=0.01 \
          --mu='1e-4' --k=200

   Note that parameters with ``useDefault`` set to ``True`` can be ignored if the
   default parameter is used. In addition, parameter ``--optimized`` could be used
   to load the optimized version of a simuPOP module. For this particular
   configuration, the optimized module is 30% faster (62s vs. 40s) than the
   standard module.

#. The simulation function could be imported to another script as follows  ::

      from simuCDCV import simuCDCV
      simuCDCV(model='instant', N0=10000, N1=10000, G0=500, G1=500,
          spec=[0.9]+[0.02]*5, s=0.01, mu=1e-4, k=200)

document

