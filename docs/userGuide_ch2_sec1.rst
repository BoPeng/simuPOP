What is simuPOP?
================

simuPOP is a **general-purpose individual-based forward-time population genetics
simulation environment** based on Python, a dynamic object-oriented programming
language that has been widely used in biological studies. More specifically,

* simuPOP is a **population genetics simulator** that simulates the evolution of
  populations. It uses a discrete generation model although overlapping
  generations could be simulated using nonrandom mating schemes.

* simuPOP explicitly models **populations with individuals** who have their own
  genotype, sex, and auxiliary information such as age. The evolution of a
  population is modeled by populating an offspring population from parents in the
  parental population.

* Unlike coalescent-based programs, simuPOP evolves populations **forward in
  time**, subject to arbitrary number of genetic and environmental forces such as
  mutation, recombination, migration and Population/subpopulation size changes.

* simuPOP is a **general-purpose** simulator that is designed to simulate
  arbitrary evolutionary processes. In contrast to competing applications that use
  command-line options or configuration files to direct the execution of a limited
  number of predefined evolutionary scenarios, users of simuPOP's scripting
  interface could make use of many of its unique features, such as customized
  chromosome types, arbitrary nonrandom mating schemes, virtual subpopulations,
  information fields and Python operators, to construct and study almost
  arbitrarily complex evolutionary scenarios. In addition, because simuPOP
  provides a large number of functions to manipulate populations, it can be used
  as an data manipulatation and analysis tool.

simuPOP is provided as a number of Python modules, which consist of a large
number of Python objects and functions, including Population, mating schemes,
operators (objects that manipulate populations) and simulators to coordinate the
evolutionary processes. It is the users' responsibility to write a Python script
to glue these pieces together and form a simulation. At a more user-friendly
level, an increasing number of functions and scripts contributed by simuPOP
users is available in the online simuPOP cookbook
(``http://simupop.sourceforge.net/cookbook``). They provide useful functions for
different applications (e.g. load and manipulate HapMap samples, import and
export files from another application) and allow users who are unfamiliar with
simuPOP to perform a large number of simulations ranging from basic population
genetics models to generating datasets under complex evolutionary scenarios.


