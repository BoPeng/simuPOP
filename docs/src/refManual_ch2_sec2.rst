Virtual splitters
=================


class BaseVspSplitter
---------------------

.. class:: BaseVspSplitter

   This class is the base class of all virtual subpopulation (VSP)
   splitters, which provide ways to define groups of individuals in a
   subpopulation who share certain properties. A splitter defines a
   fixed number of named VSPs. They do not have to add up to the whole
   subpopulation, nor do they have to be distinct. After a splitter is
   assigned to a population, many functions and operators can be
   applied to individuals within specified VSPs.

   Each VSP has a name. A default name is determined by each splitter
   but you can also assign a name to each VSP. The name of a VSP can
   be retrieved by function ``BaseVspSplitter.name()`` or
   ``Population.subPopName()``.

   Only one VSP splitter can be assigned to a population, which
   defined VSPs for all its subpopulations. If different splitters are
   needed for different subpopulations, a :class:`CombinedSplitter`
   can be used.


   .. method:: BaseVspSplitter(names=[])


      This is a virtual class that cannot be instantiated.


   .. method:: BaseVspSplitter.clone()

      All VSP splitter defines a ``clone()`` function to create an
      identical copy of itself.

   .. method:: BaseVspSplitter.name(vsp)

      Return the name of VSP *vsp* (an index between ``0`` and
      ``numVirtualSubPop()``).

   .. method:: BaseVspSplitter.numVirtualSubPop()

      Return the number of VSPs defined by this splitter.

   .. method:: BaseVspSplitter.vspByName(name)

      Return the index of a virtual subpopulation from its name. If
      multiple virtual subpopulations share the same name, the first
      vsp is returned.


class SexSplitter
-----------------

.. class:: SexSplitter

   This splitter defines two VSPs by individual sex. The first VSP
   consists of all male individuals and the second VSP consists of all
   females in a subpopulation.


   .. method:: SexSplitter(names=[])


      Create a sex splitter that defines male and female VSPs. These
      VSPs are named ``Male`` and ``Female`` unless a new set of names
      are specified by parameter *names*.


   .. method:: SexSplitter.name(vsp)

      Return ``"Male"`` if *vsp=0* and ``"Female"`` otherwise, unless
      a new set of names are specified.

   .. method:: SexSplitter.numVirtualSubPop()

      Return ``2``.


class AffectionSplitter
-----------------------

.. class:: AffectionSplitter

   This class defines two VSPs according individual affection status.
   The first VSP consists of unaffected invidiauls and the second VSP
   consists of affected ones.


   .. method:: AffectionSplitter(names=[])


      Create a splitter that defined two VSPs by affection
      status.These VSPs are named ``Unaffected`` and ``Affected``
      unless a new set of names are specified by parameter *names*.


   .. method:: AffectionSplitter.name(vsp)

      Return ``"Unaffected"`` if *vsp=0* and ``"Affected"`` if
      *vsp=1*, unless a new set of names are specified.

   .. method:: AffectionSplitter.numVirtualSubPop()

      Return 2.


class InfoSplitter
------------------

.. class:: InfoSplitter

   This splitter defines VSPs according to the value of an information
   field of each indivdiual. A VSP is defined either by a value or a
   range of values.


   .. method:: InfoSplitter(field, values=[], cutoff=[], ranges=[], names=[])


      Create an infomration splitter using information field *field*.
      If parameter *values* is specified, each item in this list
      defines a VSP in which all individuals have this value at
      information field *field*. If a set of cutoff values are defined
      in parameter *cutoff*, individuals are grouped by intervals
      defined by these cutoff values. For example, ``cutoff=[1,2]``
      defines three VSPs with ``v < 1``, ``1 <= v < 2`` and ``v >=2``
      where ``v`` is the value of an individual at information field
      *field*. If parameter ``ranges`` is specified, each range
      defines a VSP. For example, ``ranges=[[1, 3], [2, 5]]`` defines
      two VSPs with ``1 <= v < 3`` and ``2 <= 3 < 5``. Of course, only
      one of the parameters *values*, *cutoff* and *ranges* should be
      defined, and values in *cutoff* should be distinct, and in an
      increasing order. A default set of names are given to each VSP
      unless a new set of names is given by parameter *names*.


   .. method:: InfoSplitter.name(vsp)

      Return the name of a VSP *vsp*, which is ``field = value`` if
      VSPs are defined by values in parameter *values*, or ``field <
      value`` (the first VSP), ``v1 <= field < v2`` and ``field >= v``
      (the last VSP) if VSPs are defined by cutoff values. A user-
      specified name, if specified, will be returned instead.

   .. method:: InfoSplitter.numVirtualSubPop()

      Return the number of VSPs defined by this splitter, which is the
      length parameter *values* or the length of *cutoff* plus one,
      depending on which parameter is specified.


class ProportionSplitter
------------------------

.. class:: ProportionSplitter

   This splitter divides subpopulations into several VSPs by
   proportion.


   .. method:: ProportionSplitter(proportions=[], names=[])


      Create a splitter that divides subpopulations by *proportions*,
      which should be a list of float numbers (between ``0`` and
      ``1``) that add up to ``1``. A default set of names are given to
      each VSP unless a new set of names is given by parameter
      *names*.


   .. method:: ProportionSplitter.name(vsp)

      Return the name of VSP *vsp*, which is ``"Prop p"`` where
      ``p=propotions[vsp]``. A user specified name will be returned if
      specified.

   .. method:: ProportionSplitter.numVirtualSubPop()

      Return the number of VSPs defined by this splitter, which is the
      length of parameter *proportions*.


class RangeSplitter
-------------------

.. class:: RangeSplitter

   This class defines a splitter that groups individuals in certain
   ranges into VSPs.


   .. method:: RangeSplitter(ranges, names=[])


      Create a splitter according to a number of individual ranges
      defined in *ranges*. For example, ``RangeSplitter(ranges=[[0,
      20], [40, 50]])`` defines two VSPs. The first VSP consists of
      individuals ``0``, ``1``, ..., ``19``, and the sceond VSP
      consists of individuals ``40``, ``41``, ..., ``49``. Note that a
      nested list has to be used even if only one range is defined. A
      default set of names are given to each VSP unless a new set of
      names is given by parameter *names*.


   .. method:: RangeSplitter.name(vsp)

      Return the name of VSP *vsp*, which is ``"Range [a, b)"`` where
      ``[a, b)`` is range ``ranges[vsp]``. A user specified name will
      be returned if specified.

   .. method:: RangeSplitter.numVirtualSubPop()

      Return the number of VSPs, which is the number of ranges defined
      in parameter *ranges*.


class GenotypeSplitter
----------------------

.. class:: GenotypeSplitter

   This class defines a VSP splitter that defines VSPs according to
   individual genotype at specified loci.


   .. method:: GenotypeSplitter(loci, alleles, phase=False, names=[])


      Create a splitter that defines VSPs by individual genotype at
      *loci* (can be indexes or names of one or more loci). Each list
      in a list *allele* defines a VSP, which is a list of allowed
      alleles at these *loci*. If only one VSP is defined, the outer
      list of the nested list can be ignored. If phase if true, the
      order of alleles in each list is significant. If more than one
      set of alleles are given, Individuals having either of them is
      qualified.

      For example, in a haploid population, ``loci=1, alleles=[0, 1]``
      defines a VSP with individuals having allele ``0`` or ``1`` at
      locus ``1``, ``alleles=[[0, 1], [2]]`` defines two VSPs with
      indivdiuals in the second VSP having allele ``2`` at locus
      ``1``. If multiple loci are involved, alleles at each locus need
      to be defined. For example, VSP defined by ``loci=[0, 1],
      alleles=[0, 1, 1, 1]`` consists of individuals having alleles
      ``[0, 1]`` or ``[1, 1]`` at loci ``[0, 1]``.

      In a haploid population, ``loci=1, alleles=[0, 1]`` defines a
      VSP with individuals having genotype ``[0, 1]`` or ``[1, 0]`` at
      locus ``1``. ``alleles[[0, 1], [2, 2]]`` defines two VSPs with
      indivdiuals in the second VSP having genotype ``[2, 2]`` at
      locus ``1``. If *phase* is set to ``True``, the first VSP will
      only has individuals with genotype ``[0, 1]``. In the multiple
      loci case, alleles should be arranged by haplotypes, for
      example, ``loci=[0, 1], alleles=[0, 0, 1, 1], phase=True``
      defines a VSP with individuals having genotype ``-0-0-, -1-1-``
      at loci ``0`` and ``1``. If ``phase=False`` (default), genotypes
      ``-1-1-, -0-0-``, ``-0-1-`` and ``-1-0-`` are all allowed.

      A default set of names are given to each VSP unless a new set of
      names is given by parameter *names*.


   .. method:: GenotypeSplitter.name(vsp)

      Return name of VSP *vsp*, which is ``"Genotype
      loc1,loc2:genotype"`` as defined by parameters *loci* and
      *alleles*. A user provided name will be returned if specified.

   .. method:: GenotypeSplitter.numVirtualSubPop()

      number of virtual subpops of subpopulation sp


class CombinedSplitter
----------------------

.. class:: CombinedSplitter

   This splitter takes several splitters and stacks their VSPs
   together. For example, if the first splitter defines ``3`` VSPs and
   the second splitter defines ``2``, the two VSPs from the second
   splitter become the fourth (index ``3``) and the fifth (index
   ``4``) VSPs of the combined splitter. In addition, a new set of
   VSPs could be defined as the union of one or more of the original
   VSPs. This splitter is usually used to define different types of
   VSPs to a population.


   .. method:: CombinedSplitter(splitters=[], vspMap=[], names=[])


      Create a combined splitter using a list of *splitters*. For
      example, ``CombinedSplitter([SexSplitter(),
      AffectionSplitter()])`` defines a combined splitter with four
      VSPs, defined by male (vsp ``0``), female (vsp ``1``),
      unaffected (vsp ``2``) and affected individuals (vsp ``3``).
      Optionally, a new set of VSPs could be defined by parameter
      *vspMap*. Each item in this parameter is a list of VSPs that
      will be combined to a single VSP. For example, ``vspMap=[(0, 2),
      (1, 3)]`` in the previous example will define two VSPs defined
      by male or unaffected, and female or affected individuals. VSP
      names are usually determined by splitters, but can also be
      specified using parameter *names*.


   .. method:: CombinedSplitter.name(vsp)

      Return the name of a VSP *vsp*, which is the name a VSP defined
      by one of the combined splitters unless a new set of names is
      specified. If a *vspMap* was used, names from different VSPs
      will be joined by ``"or"``.

   .. method:: CombinedSplitter.numVirtualSubPop()

      Return the number of VSPs defined by this splitter, which is the
      sum of the number of VSPs of all combined splitters.


class ProductSplitter
---------------------

.. class:: ProductSplitter

   This splitter takes several splitters and take their intersections
   as new VSPs. For example, if the first splitter defines ``3`` VSPs
   and the second splitter defines ``2``, ``6`` VSPs will be defined
   by splitting 3 VSPs defined by the first splitter each to two VSPs.
   This splitter is usually used to define finer VSPs from existing
   VSPs.


   .. method:: ProductSplitter(splitters=[], names=[])


      Create a product splitter using a list of *splitters*. For
      example, ``ProductSplitter([SexSplitter(),
      AffectionSplitter()])`` defines four VSPs by male unaffected,
      male affected, female unaffected, and female affected
      individuals. VSP names are usually determined by splitters, but
      can also be specified using parameter *names*.


   .. method:: ProductSplitter.name(vsp)

      Return the name of a VSP *vsp*, which is the names of indivdual
      VSPs separated by a comma, unless a new set of names is
      specified for each VSP.

   .. method:: ProductSplitter.numVirtualSubPop()

      Return the number of VSPs defined by this splitter, which is the
      sum of the number of VSPs of all combined splitters.


