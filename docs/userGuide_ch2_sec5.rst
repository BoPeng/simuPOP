How to read this user's guide
=============================

This user's guide describes all simuPOP features using a lot of examples. The
first few chapters describes all classes in the simuPOP core. Chapter
:ref:`cha_simuPOP_Operators <cha_simuPOP_Operators>` describes almost all
simuPOP operators, divided largely by genetic models. Features listed in these
two chapters are generally implemented at the C++ level and are provided through
the :mod:`simuPOP` module. Chapter :ref:`cha_Utility_Modules
<cha_Utility_Modules>` describes features that are provided by various simuPOP
utility modules. These modules provide extensions to the simuPOP core that
improves the usability and userfriendliness of simuPOP. The next chapter
(Chapter :ref:`cha_A_real_example <cha_A_real_example>`) demonstrates how to
write a script to solve a real-world simulation problem. Because some sections
describe advanced features that are only used in the construction of highly
complex simulations, or implementation details that concern only advanced users,
new simuPOP users can safely skip these sections. **Sections that describe
advanced topics are marked by one or two asterisks (\*) after the section**
**titles**.

simuPOP is a comprehensive forward-time population genetics simulation
environment with many unique features. If you are new to simuPOP, you can go
through this guide quickly and understand what simuPOP is and what features it
provides. Then, you can read Chapter :ref:`cha_A_real_example
<cha_A_real_example>` and learn how to apply simuPOP in real-world problems.
After you play with simuPOP for a while and start to write simple scripts, you
can study relevant sections in details. The *simuPOP reference manual* will
become more and more useful when the complexity of your scripts grows.

Before we dive into the details of simuPOP, it is helpful to know a few name
conventions that simuPOP tries to follow. Generally speaking,

* All class names use the CapWords convention (e.g. :class:`Population`\ (),
  :class:`InitSex`\ ()) .

* All standalone functions (e.g. :func:`loadPopulation`\ () and ``initSex``),
  member functions (e.g. :meth:`Population.mergeSubPops`\ ()) and parameter names
  use the mixedCases style.

* Constants are written in all capital characters with underscores separating
  words (e.g. ``CHROMOSOME_X``, ``UNIFORM_DISTRIBUTION``). Their names instead of
  their actual values should be used because those values can change without
  notice.

* simuPOP uses the abbreviated form of the following words in function and
  parameter names:

     ``pop`` (population), ``pops`` (populations), ``pos`` (position),  ``info``
     (information), ``migr`` (migration), ``subPop`` (subpopulation and virtual
     subpopulation), ``subPops`` (subpopulations and virtual subpopulations), ``rep``
     (replicates), ``gen`` (generation), ``ops`` (operators), ``expr`` (expression),
     ``stmts`` (statements).


* simuPOP uses both singular and plural forms of parameters, according to the
  following rules:

* If a parameter only accept a single input, singular names such as ``field``,
    ``locus``, ``value``, and ``name`` are used.

* If a parameter accepts a list of values, plural names such as ``fields``,
    ``loci``, ``values`` and ``names`` are used. **Such parameters usually accept
    single inputs.** For example, ``loci=1`` can be used as a shortcut for
    ``loci=[1]`` and ``infoFields='x'`` can be used as a shortcut for
    ``infoFields=['x']``.

  The same rules also hold for function names. For example,
  :meth:`Population.addInfoFields`\ () accept a list of information fields but
  ``pop.addInfoFields('field')`` is also acceptable.


