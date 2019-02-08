Loading simuPOP modules
=======================


Short, long, binary, mutant and lineage modules and their optimized versions
----------------------------------------------------------------------------

There are ten flavors of the core simuPOP module: short, long, binary, mutant,
and lineage allele modules, and their optimized versions.

* The short allele modules use *8 bits* to store each allele which limits the
  possible allele states to 256. This is enough most of the times so this is the
  default module of simuPOP.

* If you need to a large number of allele states to simulate, for example the
  infinite allele model, you should use the long allele version of the modules,
  which use *32 or 64 bits* for each allele and can have :math:`2^{32}` or
  :math:`2^{64}` possible allele states depending on your platform.

* If you would like to simulate a large number of binary (SNP) markers, binary
  libraries can save you a lot of RAM because they use *1 bit* for each allele.

* If you are simulating long sequence regions with rare variants, you can use
  the mutant module. This module uses compression technology that ignores wildtype
  alleles and is not efficient if you need to traverse all alleles frequently. The
  maximum allele state is 255 for this module. Because this module stores location
  and value of each allele, it uses at least 64 + 8 bits for each allele on a 64
  bit system. The complexity of the storage also prevents simultaneous write
  access to genotypes so this module does not benefit much from running in multi-
  thread mode.

* If you are interested in tracing the lineage of each allele (e.g. the ID of
  individuals to whom the allele was introduced), you can use the lineage module
  for which each allele is attached with information about its origin. The maximum
  allele state is 255 for this module, and the cost of storing each allele is 8
  (value) + 32 (lineage) bits.

Despite of differences in internal memory layout, all these modules have the
same interface, although some functions behave differently in terms of
functionality and performance.

Standard libraries have detailed debug and run-time validation mechanism to make
sure a simulation executes correctly. Whenever something unusual is detected,
simuPOP would terminate with detailed error messages. The cost of such run-time
validation varies from case to case but can be high under some extreme
circumstances. Because of this, optimized versions for all modules are provided.
They bypass most parameter checking and run-time validations and will simply
crash if things go wrong. It is recommended that you use standard libraries
whenever possible and only use the optimized version when performance is needed
and you are confident that your simulation is running as expected.

Examples :ref:`lst_Use_of_standard_module <lst_Use_of_standard_module>` and
:ref:`lst_Use_of_optimized_module <lst_Use_of_optimized_module>` demonstrate the
differences between standard and optimized modules, by executing two invalid
commands. A standard module checks all input values and raises exceptions when
invalid inputs are detected. An interactive Python session would catch these
exceptions and print proper error messages. In constrast, an optimized module
returns erroneous results and or simply crashes when such inputs are given.

.. _lst_Use_of_standard_module:

**Example**: *Use of standard simuPOP modules*

::

   >>> import simuPOP as sim
   >>> pop = sim.Population(10, loci=2)
   >>> pop.locusPos(10)
   Traceback (most recent call last):
     File "/var/folders/ys/gnzk0qbx5wbdgm531v82xxljv5yqy8/T/tmp6boewtoh", line 1, in <module>
       #begin_file log/standard.py
   IndexError: genoStru.h: 557 absolute locus index (10) out of range of 0 ~ 1
   >>> pop.individual(20).setAllele(1, 0)
   Traceback (most recent call last):
     File "/var/folders/ys/gnzk0qbx5wbdgm531v82xxljv5yqy8/T/tmp6boewtoh", line 1, in <module>
       #begin_file log/standard.py
   IndexError: population.h: 566 individual index (20)  out of range of 0 ~ 9

   now exiting runScriptInteractively...

`Download standard.py <standard.py>`_

.. index:: single: r

Example :ref:`lst_Use_of_optimized_module <lst_Use_of_optimized_module>` also
demonstrates how to use the :func:`setOptions` function in the :mod:`simuOpt`
module to control the choice of one of the six simuPOP modules. By specifying
one of the values ``short, long`` or ``binary``\ for option ``alleleType``, and
setting\ ``optimized`` to ``True`` or ``False``, the right flavor of module will
be chosen when simuPOP is loaded. In addition, option ``quiet`` can be used
suppress the banner message when simuPOP is loaded. An alternative method is to
set environmental variable ``SIMUALLELETYPE`` to ``short``, ``long`` or
``binary`` to use the standard short, long or binary module, and variable
``SIMUOPTIMIZED`` to use the optimized modules. Command line options
``--optimized`` can also be used.

.. _lst_Use_of_optimized_module:

**Example**: *Use of optimized simuPOP modules*

::

   % python
   >>> from simuOpt import setOptions
   >>> setOptions(optimized=True, alleleType='long', quiet=True)
   >>> import simuPOP as sim
   >>> pop = sim.Population(10, loci=[2])
   >>> pop.locusPos(10)
   1.2731974748756028e-313
   >>> pop.individual(20).setAllele(1, 0)
   Segmentation fault


Execution in multiple threads
-----------------------------

simuPOP is capable of executing in multiple threads but it by default only makes
use of one thread. If you have a multi-core CPU, it is often beneficial to set
the number of threads to 2 or more to take advantage of this feature. The
recommended number of threads is usually the number of cores of your CPU but you
might want to set it to a lower number to leave room for the execution of other
applications. The number of threads used in simuPOP can be controlled in the
following ways:

* If an environmental variable ``OMP_NUM_THREADS`` is set to a positive number,
  simuPOP will be started with specified number of threads.

* Before simuPOP is imported, you can set the number of threads using function
  :func:`simuOpt.setOptions`\ (``numThreads=x``) where ``x`` can be a positive
  number (number of threads) or ``0``, which is intepreted as the number of cores
  available for your computer.

The number of threads a simuPOP session is used will be displayed in the banner
message when simuPOP is imported, and can be retrieved through
:func:`moduleInfo`\ ``['threads']``.

Although simuPOP can usually benefit from the use of multiple cores, certain
features of your script might prevent the execution of simuPOP in multiple
threads. For example, if your script uses a sex mode of ``GLOBAL_SEX_SEQUENCE``
to set the sex of offspring according to the global sequence of sexes (e.g.
male, male, female), simuPOP will only use on thread to generate offspring
because it is not feasible to assign individual sex from a single source of list
across multiple threads.


Graphical user interface
------------------------

A complete graphical user interface (GUI) for users to interactively construct
evolutionary processes is still in the planning stage. However, some simuPOP
classes and functions can make use of a GUI to improve user interaction. For
example, a parameter input dialog can be constructed automatically from a
parameter specification list, and be used to accept user input if class
``simuOpt.Params`` is used to handle parameters. Other examples include a
progress bar :class:`simuPOP.utils.ProgressBar` and a dialog used by function
:func:`simuPOP.utils.viewVars` to display a large number of variables. The most
notable feature of the use of GUI in simuPOP is that **all functionalities can
be achieved without a GUI**. For examples, ``simuOpt.getParam`` will use a
terminal to accept user input interactively and
:class:`simuPOP.utils.ProgressBar` will turn to a text-based progress bar in the
non-GUI mode.

The use of GUI can be controlled either globally or Individually. First, a
global GUI parameter could be set by environmental variable ``SIMUGUI``,
function :func:`simuOpt.setOptions`\ (``gui``) or a command line option
``--gui`` of a simuPOP scripts. Allowed values include

* ``True``: This is the system default value. A GUI is used whenever possible.
  All GUI-capable functions support ``wxPython`` so a ``wxPython`` dialog will be
  used if ``wxPython`` is available. Otherwise, ``tkInter`` based dialogs or text-
  mode will be used.

* ``False``: no GUI will be used. All functions will use text-based
  implementation. Note that ``--gui=False`` is commonly used to run scripts in
  batch mode.

* ``wxPython``: Force the use of ``wxPython`` GUI toolkit.

* ``Tkinter``: Force the use of ``Tkinter`` GUI toolkit.

Individual classes and functions that could make use a GUI usually have their
own ``gui`` parameters, which can be set to override global GUI settings. For
example, you could force the use of a text-based progress bar by using
``ProgressBar(gui=False)``.


