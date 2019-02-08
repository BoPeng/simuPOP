Features
========

simuPOP offers a long list of features, many of which are unique among all
forward-time population genetics simulation programs. The most distinguishing
features include:

#. simuPOP provides three types of modules that use 1, 8 or 32/64 bits to store
   an allele. The binary module (1 bit) is suitable for simulating a large number
   of SNP markers, and the long module (32 or 64 bits depending on platform) is
   suitable for simulating some population genetics models such as the infinite
   allele mutation model.

#. [NEW in simuPOP 1.0.7] simuPOP provides modules to store a large number of
   rare variants in a compressed manner (the mutant module), and to store origin of
   each allele so that it is easy to track allelic lineage during evolution.

#. The core of simuPOP is implemented in C++ which is heavily optimized for
   large-scale simulations. simuPOP can be executed in multiple threads with
   boosted performance on modern multi-core CPUs.

#. In addition to autosomes and sex chromosomes, simuPOP supports arbitrary
   types of chromosomes through customizable genotype transmitters. Random maternal
   transmission of mitochondrial DNAs is supported as a special case of this
   feature.

#. An arbitrary number of float numbers, called **information fields**, can be
   attached to individuals of a population. For example, information field
   ``father_idx`` and ``mother_idx`` can be used to track an individual's parents,
   and ``pack_year`` can be used to simulate an environmental factor associated
   with smoking.

#. simuPOP does not impose a limit on the number of homologous sets of
   chromosomes, the size of the genome or populations. The size of your simulation
   is only limited by the physical memory of your computer.

#. During an evolutionary process, a population can hold more than one most-
   recent generation. Pedigrees can be sampled from such multi-generation
   populations.

#. An operator can be native (implemented in C++) or hybrid (Python-assisted). A
   hybrid operator calls a user-provided Python function to implement arbitrary
   genetic effects. For example, a hybrid mutator passes to-be-mutated alleles to a
   function and mutates these alleles according to the returned values.

#. simuPOP provides more than 60 operators that cover all important aspects of
   genetic studies. These include mutation (*e.g. k*-allele, stepwise, generalized
   stepwise and context-sensitive mutation models), migration (arbitrary, can
   create new subpopulation), recombination and gene conversion (uniform or
   nonuniform), selection (single-locus, additive, multiplicative or hybrid multi-
   locus models, support selection of both parents and offspring), penetrance
   (single, multi-locus or hybrid), ascertainment (casecontrol, affected sibpairs,
   random, nuclear and large Pedigrees), statistics calculation (including but not
   limited to allele, genotype, haplotype, heterozygote number and frequency;
   linkage disequilibrium measures, Hardy-Weinberg test), pedigree tracing,
   visualization (using R or other Python modules) and load/save in simuPOP's
   native format and many external formats such as Linkage.

#. Mating schemes can work on virtual subpopulations of a subpopulation. For
   example, positive assortative mating can be implemented by mating individuals
   with similar properties such as ancestry and overlapping generations could be
   simulated by copying individuals acorss generations. The number of offspring per
   mating event can be fixed or can follow a statistical distribution.

A number of forward-time simulation programs are available. If we exclude early
forward-time simulation applications developed primarily for teaching purposes,
notable forward-time simulation programs include *easyPOP*, *FPG*, *Nemo* and
*quantiNemo*, *genoSIM* and *genomeSIMLA*, *FreGene*, *GenomePop*, *ForwSim*,
and *ForSim*. These programs are designed with specific applications and
specific evolutionary scenarios in mind, and excel in what they are designed
for. For some applications, these programs may be easier to use than simuPOP.
For example, using a special look-ahead algorithm, *ForwSim* is among the
fastest programs to simulate a standard Wright-Fisher process, and should be
used if such a simulation is needed. However, these programs are not flexible
enough to be applied to problems outside of their designed application area. For
example, none of these programs can be used to study the evolution of a disease
predisposing mutant, a process that is of great importance in statistical
genetics and genetic epidemiology. Compared to such programs, simuPOP has the
following advantages:

* The scripting interface gives simuPOP the flexibility to create arbitrarily
  complex evolutionary scenarios. For example, it is easy to use simuPOP to
  explicitly introduce a disease predisposing mutant to an evolving population,
  trace the allele frequency of them, and restart the simulation if they got lost
  due to genetic drift.

* The Python interface allows users to define customized genetic effects in
  Python. In contrast, other programs either do not allow customized effects or
  force users to modify code at a lower (e.g. C++) level.

* simuPOP is the only application that embodies the concept of virtual
  subpopulation that allows evolutions at a finer scale. This is required for
  realistic simulations of complex evolutionary scenarios.

* simuPOP allows users to examine an evolutionary process very closely because
  all simuPOP objects are Python objects that can be assessed using their member
  functions. For example, users can keep track of genotype at particular loci
  during evolution. In contrast, other programs work more or less like a black box
  where only limited types of statistics can be outputted.


