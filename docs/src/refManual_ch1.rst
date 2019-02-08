.. _front:

************
Front Matter
************


.. topic:: Abstract

   simuPOP is a general-purpose individual-based forward-time population genetics
   simulation environment. Unlike coalescent-based programs, simuPOP evolves
   populations forward in time, subject to arbitrary number of genetic and
   environmental forces such as mutation, recombination, migration and
   population/subpopulation size changes. In contrast to competing applications
   that use command-line options or configuration files to direct the execution of
   a limited number of predefined evolutionary scenarios, users of simuPOP's
   scripting interface could make use of many of its unique features, such as
   customized chromosome types, arbitrary nonrandom mating schemes, virtual
   subpopulations, information fields and Python operators, to construct and study
   almost arbitrarily complex evolutionary scenarios.

   simuPOP is provided as a number of Python modules, which consist of a large
   number of Python objects and functions, including population, mating schemes,
   operators (objects that manipulate populations) and simulators to coordinate the
   evolutionary processes. It is the users' responsibility to write a Python script
   to glue these pieces together and form a simulation. At a more user-friendly
   level, an increasing number of functions and scripts contributed by simuPOP
   users is available in the online simuPOP cookbook. They provide useful functions
   for different applications (e.g. load and manipulate HapMap samples, import and
   export files from another application) and allow users who are unfamiliar with
   simuPOP to perform a large number of simulations ranging from basic population
   genetics models to generating datasets under complex evolutionary scenarios.

   This document provides complete references to all classes and functions of
   simuPOP and its utility modules. Please refer to the *simuPOP user's guide* for
   a detailed introduction to simuPOP concepts, and a number of examples on how to
   use simuPOP to perform various simulations. All resources, including a pdf
   version of this guide and a mailing list can be found at the simuPOP homepage
   ``http://simupop.sourceforge.net``.

   **How to cite simuPOP:**

      Bo Peng and Marek Kimmel (2005) simuPOP: a forward-time population genetics
      simulation environment. *bioinformatics*, **21** (18): 3686-3687.

      Bo Peng and Christopher Amos (2008) Forward-time simulations of nonrandom mating
      populations using simuPOP. *bioinformatics*, **24** (11): 1408-1409.



