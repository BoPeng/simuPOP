%feature("docstring") simuPOP::affectionSplitter "

Details:

    This class defines two VSPs according individual affection status.
    The first VSP consists of unaffected invidiauls and the second VSP
    consists of affected ones.

"; 

%feature("docstring") simuPOP::affectionSplitter::affectionSplitter "

Description:

    Create a splitter that defined two VSPs by affection status.

Usage:

    affectionSplitter()

"; 

%feature("docstring") simuPOP::affectionSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::affectionSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::affectionSplitter::numVirtualSubPop "

Description:

    Return 2.

Usage:

    x.numVirtualSubPop()

"; 

%ignore simuPOP::affectionSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::affectionSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::affectionSplitter::name "

Description:

    Return \"Unaffected\" if vsp=0 and \"Affected\" if vsp=1.

Usage:

    x.name(vsp)

"; 

%feature("docstring") simuPOP::alphaParentsChooser "

Applicability: all ploidy

Details:

    This parent chooser chooses two parents randomly, a male and a
    female, from their respective sex groups randomly. If selection is
    turned on, parents are chosen from their sex groups with
    probabilities that are proportional to their fitness values. This
    parents chooser also allows polygamous mating by reusing a parent
    multiple times when returning parents, and allows specification of
    a few alpha individuals who will be the only mating individuals in
    their sex group.

"; 

%feature("docstring") simuPOP::alphaParentsChooser::alphaParentsChooser "

Usage:

    alphaParentsChooser(alphaSex=Male, alphaNum=0,
      alphaField=string)

Details:

    Note: If selection is enabled, it works regularly on on-alpha sex,
    but works twice on alpha sex. That is to say, alphaNum alpha
    indiviudals are chosen selectively, and selected again during
    mating.

Arguments:

    replacement:    choose with (True, default) or without (False)
                    replacement. When choosing without replacement,
                    parents will be paired and can only mate once.
    replenish:      if set to true, one or both sex groups will be
                    replenished if they are exhausted.
    polySex:        Male (polygyny) or Female (polyandry) parent that
                    will have polyNum sex partners.
    polyNum:        Number of sex partners.
    alphaSex:       the sex of the alpha individual, i.e. alpha male
                    or alpha female who be the only mating individuals
                    in their sex group.
    alphaNum:       Number of alpha individuals. If infoField is not
                    given, alphaNum random individuals with alphaSex
                    will be chosen. If selection is enabled,
                    individuals with higher fitness values have higher
                    probability to be selected. There is by default no
                    alpha individual (alphaNum = 0).
    alphaField:     if an information field is given, individuals with
                    non-zero values at this information field are
                    alpha individuals. Note that these individuals
                    must have alphaSex.

"; 

%feature("docstring") simuPOP::alphaParentsChooser::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::alphaParentsChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::alphaParentsChooser::chooseParents(RawIndIterator basePtr);

%ignore simuPOP::alphaParentsChooser::numMale();

%ignore simuPOP::alphaParentsChooser::numFemale();

%feature("docstring") simuPOP::baseOperator "

Details:

    Operators are objects that act on populations. They can be applied
    to populations directly using their function forms, but they are
    usually managed and applied by a simulator. In the latter case,
    operators are passed to the evolve function of a simulator, and
    are applied repeatedly during the evolution of the simulator.
    The baseOperator class is the base class for all operators. It
    defines a common user interface that specifies at which
    generations, at which stage of a life cycle, to which populations
    and subpopulation an operator will be applied. These are achieved
    by a common set of parameters such as begin, end, step, at, stage
    for all operators. Note that a specific operator does not have to
    honor all these parameters. For example, a recombinator can only
    be applied during mating so it ignores the stage parameter.
    An operator can be applied to all or part of the generations
    during the evolution of a simulator. At the beginning of an
    evolution, a simulator is usually at the beginning of generation
    0. If it evolves 10 generations, it evolves generations 0, 1,
    ,,,., and 9 (10 generations) and stops at the begging of
    generation 10. A negative generation number a has generation
    number 10 + a, with -1 referring to the last evolved generation 9.
    Note that the starting generation number of a simulator can be
    changed by its setGen() member function.
    Output from an operator is usually directed to the standard output
    (sys.stdout). This can be configured using a output specification
    string, which can be '' for no output, '>' standard terminal
    output (default), or a filename prefixed by one or more '>'
    characters. In the case of '>filename' (or equivalently
    'filename'), the output from an operator is written to this file.
    However, if two operators write to the same file filename, or if
    an operator write to this file more than once, only the last write
    operation will succeed. In the case of '>>filename', file filename
    will be opened at the beginning of the evolution and closed at the
    end. Outputs from multiple operators are appended. >>>filename
    works similar to >>filename but filename, if it already exists at
    the beginning of an evolutionary process, will not be cleared.

"; 

%feature("docstring") simuPOP::baseOperator::baseOperator "

Usage:

    baseOperator(output, outputExpr, stage, begin, end, step, at,
      rep, subPop, infoFields)

Details:

    The following parameters can be specified by all operators.
    However, an operator can ignore some parameters and the exact
    meaning of a parameter can vary.

Arguments:

    output:         A string that specifies how output from an
                    operator is written, which can be '' (no output),
                    '>' (standard output), or 'filename' prefixed by
                    one or more '>'.
    outputExpr:     An expression that determines the output parameter
                    dynamically. This expression will be evaluated
                    against a population's local namespace each time
                    when an output filename is required. For example,
                    \"'>>out%s_%s.xml' % (gen, rep)\" will output to
                    >>out10_1.xml for replicate 1 at generation 10.
    stage:          Stage(s) of a life cycle at which an operator will
                    be applied. It can be PreMating, DuringMating,
                    PostMating and any of their combined stages
                    PrePostMating, PreDuringMatingDuringPostMating and
                    PreDuringPostMating. Note that all operators have
                    their default stage parameter and some of them
                    ignores this parameter because they can only be
                    applied at certain stage(s) of a life cycle.
    begin:          The starting generation at which an operator will
                    be applied. Default to 0. A negative number is
                    interpreted as a generation counted from the end
                    of an evolution (-1 being the last evolved
                    generation).
    end:            The last generation at which an operator will be
                    applied. Default to -1, namely the last
                    generation.
    step:           The number of generations between applicable
                    generations. Default to 1.
    at:             A list of applicable generations. Parameters
                    begin, end, and step will be ignored if this
                    parameter is specified.
    rep:            A list of applicable replicates. An empty list
                    (default) is interpreted as all replicates in a
                    simulator. Negative indexes such as -1 (last
                    replicate) is acceptable. rep=idx can be used as a
                    shortcut for rep=[idx].
    subPop:         A list of applicable (virtual) subpopulations,
                    such as subPop=[sp1, sp2, (sp2, vsp1)]. An empty
                    list (default) is interpreted as all
                    subpopulations. subPop=[sp1] can be simplied as
                    subPop=sp1. Negative indexes are not supported.
                    Suport for this parameter vary from operator to
                    operator. Some operators do not support virtual
                    subpopulations and some operators do not support
                    this parameter at all. Please refer to the
                    reference manual of individual operators for their
                    support for this parameter.
    infoFields:     A list of information fields that will be used by
                    an operator. You usually do not need to specify
                    this parameter because operators that use
                    information fields usually have default values for
                    this parameter.

"; 

%feature("docstring") simuPOP::baseOperator::~baseOperator "

Description:

    destroy an operator

Usage:

    x.~baseOperator()

"; 

%feature("docstring") simuPOP::baseOperator::clone "

Usage:

    x.clone()

Details:

    Return a cloned copy of an operator. This function is available to
    all operators.

"; 

%ignore simuPOP::baseOperator::isActive(UINT rep, UINT numRep, long gen, long end, bool repOnly=false);

%ignore simuPOP::baseOperator::setApplicableStage(int stage);

%ignore simuPOP::baseOperator::canApplyPreMating();

%ignore simuPOP::baseOperator::canApplyDuringMating();

%ignore simuPOP::baseOperator::canApplyPostMating();

%ignore simuPOP::baseOperator::isCompatible(const population &pop);

%ignore simuPOP::baseOperator::haploidOnly();

%ignore simuPOP::baseOperator::diploidOnly();

%ignore simuPOP::baseOperator::setHaploidOnly();

%ignore simuPOP::baseOperator::setDiploidOnly();

%ignore simuPOP::baseOperator::infoSize();

%ignore simuPOP::baseOperator::infoField(UINT idx);

%ignore simuPOP::baseOperator::formOffGenotype();

%ignore simuPOP::baseOperator::setFormOffGenotype(bool flag=true);

%feature("docstring") simuPOP::baseOperator::apply "

Usage:

    x.apply(pop)

Details:

    Apply an operator to population pop directly, without checking its
    applicability.

"; 

%ignore simuPOP::baseOperator::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%ignore simuPOP::baseOperator::getOstream(PyObject *dict=NULL, bool readable=false);

%ignore simuPOP::baseOperator::closeOstream();

%ignore simuPOP::baseOperator::atRepr();

%feature("docstring") simuPOP::baseOperator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::baseOperator::noOutput();

%feature("docstring") simuPOP::baseOperator::initialize "

Usage:

    x.initialize(pop)

"; 

%ignore simuPOP::baseOperator::applicableSubPops() const;

%feature("docstring") simuPOP::BernulliTrials "

Details:

    this class encapsulate behavior of a sequence of Bernulli trial.
    the main idea is that when doing a sequence of Bernulli trials of
    the same probability, we can use much quicker algorithms instead
    of doing n Bernulli trials For example, when N=10000, p=0.001. The
    usual way to do N Bin(p) trials is to do N randUnif(0,1)<p
    comparison. using the new method, we can use geometric
    distrubution to find the next true event. Also, for the cases of
    p=0.5, random bits are generated. This class maintain a two
    dimensional table: a vector of probabilities cross expected number
    of trials p1 p2 p3 p4 p5 trial 1 trial 2 ... trial N We expect
    that N is big (usually populaiton size) and p_i are small using
    fast BernulliTrial method for fix p, we can fill up this table
    very quickly column by column This class will provide easy access
    to row (each trial) or column (called each prob) of this table. if
    this table is accessed row by row (each trial), a internal index
    is used. if index exceeds N, trials will be generated all again.
    if trial will be called, e.g., N+2 times all the time, this
    treatment might not be very efficient.

"; 

%ignore simuPOP::BernulliTrials::BernulliTrials(RNG &rng);

%feature("docstring") simuPOP::BernulliTrials::BernulliTrials "

Usage:

    BernulliTrials(rng, prob, trials)

"; 

%feature("docstring") simuPOP::BernulliTrials::~BernulliTrials "

Usage:

    x.~BernulliTrials()

"; 

%ignore simuPOP::BernulliTrials::trialSize() const;

%feature("docstring") simuPOP::BernulliTrials::probSize "

Usage:

    x.probSize()

"; 

%ignore simuPOP::BernulliTrials::setParameter(const vectorf &prob, ULONG trials);

%feature("docstring") simuPOP::BernulliTrials::doTrial "

Description:

    generate the trial table, reset m_cur

Usage:

    x.doTrial()

"; 

%ignore simuPOP::BernulliTrials::curTrial();

%feature("docstring") simuPOP::BernulliTrials::trial "

Description:

    if necessary, do trail again.

Usage:

    x.trial()

"; 

%feature("docstring") simuPOP::BernulliTrials::trialSucc "

Usage:

    x.trialSucc(idx)

"; 

%feature("docstring") simuPOP::BernulliTrials::probFirstSucc "

Usage:

    x.probFirstSucc()

"; 

%feature("docstring") simuPOP::BernulliTrials::probNextSucc "

Usage:

    x.probNextSucc(pos)

"; 

%feature("docstring") simuPOP::BernulliTrials::trialFirstSucc "

Usage:

    x.trialFirstSucc(idx)

"; 

%feature("docstring") simuPOP::BernulliTrials::trialNextSucc "

Usage:

    x.trialNextSucc(idx, pos)

"; 

%feature("docstring") simuPOP::BernulliTrials::setTrialSucc "

Usage:

    x.setTrialSucc(idx, succ)

"; 

%feature("docstring") simuPOP::BernulliTrials::trialSuccRate "

Description:

    return the succ rate for one index, used for verification pruposes

Usage:

    x.trialSuccRate(index)

"; 

%feature("docstring") simuPOP::BernulliTrials::probSuccRate "

Description:

    return the succ rate for current trial, used for verification
    pruposes

Usage:

    x.probSuccRate()

"; 

%ignore simuPOP::BernulliTrials::probabilities();

%feature("docstring") simuPOP::cloneGenoTransmitter "

Details:

    This during mating operator copies parental genotype directly to
    offspring. This operator works for all mating schemes when one or
    two parents are involved. If both parents are passed, maternal
    genotype are copied.

"; 

%feature("docstring") simuPOP::cloneGenoTransmitter::cloneGenoTransmitter "

Usage:

    cloneGenoTransmitter(begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Details:

    Create a cloneGenoTransmitter.

"; 

%feature("docstring") simuPOP::cloneGenoTransmitter::clone "

Usage:

    x.clone()

Details:

    Return a cloned copy of an operator. This function is available to
    all operators.

"; 

%feature("docstring") simuPOP::cloneGenoTransmitter::transmitGenotype "

Usage:

    x.transmitGenotype(parent, parPloidy, offspring, ploidy)

Details:

    Transmit the parPloidy set of homologous chromosomes from parent
    to the ploidy set of homologous chromosomes of offspring.
    Customized chromosomes are not copied.

"; 

%feature("docstring") simuPOP::cloneGenoTransmitter::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::cloneGenoTransmitter::initialize "

Usage:

    x.initialize(pop)

"; 

%ignore simuPOP::cloneGenoTransmitter::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::CombinedAlleleIterator "

Details:

    this class implements a C++ iterator class that iterate through
    infomation fields in a (sub)population using 1. an IndIterator
    that will skip invisible individuals, or 2. a gapped iterator that
    will run faster. Note that 1, 2 should yield identical result, and
    2 should be used when there is no virtual subpopulation.q

"; 

%feature("docstring") simuPOP::CombinedAlleleIterator::CombinedAlleleIterator "

Usage:

    CombinedAlleleIterator()

"; 

%feature("docstring") simuPOP::CombinedAlleleIterator::ptr "

Usage:

    x.ptr()

"; 

%feature("docstring") simuPOP::combinedSplitter "

Details:

    This splitter takes several splitters and stacks their VSPs
    together. For example, if the first splitter defines 3 VSPs and
    the second splitter defines 2, the two VSPs from the second
    splitter becomes the fourth (index 3) and the fifth (index 4) VSPs
    of the combined splitter. This splitter is usually used to define
    different types of VSPs to a population.

"; 

%feature("docstring") simuPOP::combinedSplitter::combinedSplitter "

Usage:

    combinedSplitter(splitters=[])

Details:

    Create a combined splitter using a list of splitters. For example,
    combinedSplitter([sexSplitter(), affectionSplitter()]) defines a
    combined splitter with four VSPs.

"; 

%feature("docstring") simuPOP::combinedSplitter::~combinedSplitter "

Usage:

    x.~combinedSplitter()

"; 

%feature("docstring") simuPOP::combinedSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::combinedSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::combinedSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs defined by this splitter, which is the
    sum of the number of VSPs of all combined splitters.

"; 

%ignore simuPOP::combinedSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::combinedSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::combinedSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of a VSP vsp, which is the name a VSP defined by
    one of the combined splitters.

"; 

%feature("docstring") simuPOP::dumper "

Description:

    dump the content of a population.

"; 

%feature("docstring") simuPOP::dumper::dumper "

Description:

    dump a population

Usage:

    dumper(genotype=True, structure=True, ancGen=0, width=1,
      max=100, chrom=[], loci=[], subPop=[], indRange=[], output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=[], infoFields=[])

Arguments:

    genotype:       Whether or not display genotype
    structure:      Whether or not display genotypic structure
    width:          number of characters to display an allele. Default
                    to 1.
    ancGen:         how many ancestral generations to display
    chrom:          chromosome(s) to display
    loci:           loci to display
    subPop:         only display subpopulation(s)
    indRange:       range(s) of individuals to display
    max:            the maximum number of individuals to display.
                    Default to 100. This is to avoid careless dump of
                    huge populations.
    output:         output file. Default to the standard output.
    outputExpr:     and other parameters: refer to
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::dumper::clone "

Usage:

    x.clone()

Details:

    Return a cloned copy of an operator. This function is available to
    all operators.

"; 

%feature("docstring") simuPOP::dumper::apply "

Usage:

    x.apply(pop)

Details:

    Apply an operator to population pop directly, without checking its
    applicability.

"; 

%feature("docstring") simuPOP::dumper::~dumper "

Usage:

    x.~dumper()

"; 

%feature("docstring") simuPOP::dumper::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::Exception "

Description:

    exception handler. Exceptions will be passed to Python.

"; 

%feature("docstring") simuPOP::Exception::Exception "

Usage:

    Exception(msg)

Arguments:

    msg:            error message

"; 

%feature("docstring") simuPOP::Exception::message "

Description:

    return error message

Usage:

    x.message()

"; 

%feature("docstring") simuPOP::Exception::~Exception "

Usage:

    x.~Exception()

"; 

%ignore simuPOP::Expression;

%feature("docstring") simuPOP::Expression::Expression "

Usage:

    Expression(expr=\"\", stmts=\"\", locals=None)

"; 

%feature("docstring") simuPOP::Expression::~Expression "

Usage:

    x.~Expression()

"; 

%ignore simuPOP::Expression::Expression(const Expression &rhs);

%ignore simuPOP::Expression::setLocalDict(PyObject *dict);

%ignore simuPOP::Expression::empty();

%ignore simuPOP::Expression::setExpr(const string &expr="");

%ignore simuPOP::Expression::setStmts(const string &stmts="");

%ignore simuPOP::Expression::evaluate();

%ignore simuPOP::Expression::valueAsBool();

%ignore simuPOP::Expression::valueAsInt();

%ignore simuPOP::Expression::valueAsDouble();

%ignore simuPOP::Expression::valueAsString();

%ignore simuPOP::Expression::valueAsArray();

%ignore simuPOP::Expression::valueAsStrDict();

%ignore simuPOP::Expression::valueAsIntDict();

%ignore simuPOP::GenoStructure;

%ignore simuPOP::GenoStructure::GenoStructure();

%ignore simuPOP::GenoStructure::GenoStructure(UINT ploidy, const vectoru &loci, const vectoru &chromTypes, bool haplodiploid, const vectorf &lociPos, const vectorstr &chromNames, const vectorstr &alleleNames, const vectorstr &lociNames, const vectorstr &infoFields);

%feature("docstring") simuPOP::GenoStructure::~GenoStructure "

Description:

    destructor, do nothing.

Usage:

    x.~GenoStructure()

"; 

%ignore simuPOP::GenoStructure::locusPos(UINT locus) const;

%ignore simuPOP::GenoStructure::chromIndex(UINT ch) const;

%ignore simuPOP::GenoStructure::setChromTypes(const vectoru &chromTypes);

%feature("docstring") simuPOP::GenoStruTrait "

Details:

    All individuals in a population share the same genotypic
    properties such as number of chromosomes, number and position of
    loci, names of markers, chromosomes, and information fields. These
    properties are stored in this GenoStruTrait class and are
    accessible from individual, population, and simulator classes.
    Currently, a genotypic structure consists of
    *  Ploidy, namely the number of homologous sets of chromosomes, of
    a population. Haplodiploid population is also supported.
    *  Number of chromosomes and number of loci on each chromosome.
    *  Positions of loci, which determine the relative distance
    between loci on the same chromosome. No unit is assumed so these
    positions can be ordinal (1, 2, 3, ..., the default), in physical
    distance (bp, kb or mb), or in map distance (e.g. centiMorgan)
    depending on applications.
    *  Names of alleles. Although alleles at different loci usually
    have different names, simuPOP uses the same names for alleles
    across loci for simplicity.
    *  Names of loci and chromosomes.
    *  Names of information fields attached to each individual. In
    addition to basic property access functions, this class also
    provides some utility functions such as locusByName, which looks
    up a locus by its name.

"; 

%feature("docstring") simuPOP::GenoStruTrait::GenoStruTrait "

Usage:

    GenoStruTrait()

Details:

    A GenoStruTrait object is created with the creation of a
    population so it cannot be initialized directly.

"; 

%ignore simuPOP::GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru &loci, const vectoru &chromTypes, bool haplodiploid, const vectorf &lociPos, const vectorstr &chromNames, const vectorstr &alleleNames, const vectorstr &lociNames, const vectorstr &infoFields);

%ignore simuPOP::GenoStruTrait::setGenoStructure(GenoStructure &rhs);

%ignore simuPOP::GenoStruTrait::setGenoStruIdx(size_t idx);

%feature("docstring") simuPOP::GenoStruTrait::lociDist "

Usage:

    x.lociDist(loc1, loc2)

Details:

    Return the distance between loci loc1 and loc2 on the same
    chromosome. A negative value will be returned if loc1 is after
    loc2.

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociLeft "Obsolete or undocumented function."

%feature("docstring") simuPOP::GenoStruTrait::distLeft "Obsolete or undocumented function."

%feature("docstring") simuPOP::GenoStruTrait::lociCovered "Obsolete or undocumented function."

%ignore simuPOP::GenoStruTrait::gsAddChromFromStru(size_t idx) const;

%ignore simuPOP::GenoStruTrait::gsAddLociFromStru(size_t idx) const;

%ignore simuPOP::GenoStruTrait::gsRemoveLoci(const vectoru &loci, vectoru &kept);

%ignore simuPOP::GenoStruTrait::gsAddChrom(const vectorf &lociPos, const vectorstr &lociNames, const string &chromName, UINT chromType) const;

%ignore simuPOP::GenoStruTrait::gsAddLoci(const vectoru &chrom, const vectorf &pos, const vectorstr &names, vectoru &newIndex) const;

%ignore simuPOP::GenoStruTrait::genoStru() const;

%ignore simuPOP::GenoStruTrait::genoStruIdx() const;

%feature("docstring") simuPOP::GenoStruTrait::ploidy "

Usage:

    x.ploidy()

Details:

    return the number of homologous sets of chromosomes, specified by
    the ploidy parameter of the population function. Return 2 for a
    haplodiploid population because two sets of chromosomes are stored
    for both males and females in such a population.

"; 

%feature("docstring") simuPOP::GenoStruTrait::ploidyName "

Usage:

    x.ploidyName()

Details:

    return the ploidy name of this population, can be one of haploid,
    diploid, haplodiploid, triploid, tetraploid or #-ploid where # is
    the ploidy number.

"; 

%feature("docstring") simuPOP::GenoStruTrait::numLoci "

Usage:

    x.numLoci(chrom)

Details:

    return the number of loci on chromosome chrom, equivalent to
    numLoci()[chrom].

"; 

%feature("docstring") simuPOP::GenoStruTrait::numLoci "

Usage:

    x.numLoci()

Details:

    return the number of loci on all chromosomes.

"; 

%ignore simuPOP::GenoStruTrait::chromX() const;

%ignore simuPOP::GenoStruTrait::chromY() const;

%ignore simuPOP::GenoStruTrait::customizedChroms() const;

%feature("docstring") simuPOP::GenoStruTrait::sexChrom "Obsolete or undocumented function."

%feature("docstring") simuPOP::GenoStruTrait::isHaplodiploid "Obsolete or undocumented function."

%feature("docstring") simuPOP::GenoStruTrait::totNumLoci "

Usage:

    x.totNumLoci()

Details:

    return the total number of loci on all chromosomes.

"; 

%feature("docstring") simuPOP::GenoStruTrait::genoSize "Obsolete or undocumented function."

%feature("docstring") simuPOP::GenoStruTrait::locusPos "

Usage:

    x.locusPos(loc)

Details:

    return the position of locus loc specified by the lociPos
    parameter of the population function. An IndexError will be raised
    if the absolute index loc is greater than or equal to the total
    number of loci.

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociPos "

Usage:

    x.lociPos()

Details:

    return the positions of all loci, specified by the lociPos
    prameter of the population function. The default positions are 1,
    2, 3, 4, ... on each chromosome.

"; 

%feature("docstring") simuPOP::GenoStruTrait::numChrom "

Usage:

    x.numChrom()

Details:

    return the number of chromosomes.

"; 

%ignore simuPOP::GenoStruTrait::chromIndex() const;

%feature("docstring") simuPOP::GenoStruTrait::chromBegin "

Usage:

    x.chromBegin(chrom)

Details:

    return the index of the first locus on chromosome chrom.

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromEnd "

Usage:

    x.chromEnd(chrom)

Details:

    return the index of the last locus on chromosome chrom plus 1.

"; 

%feature("docstring") simuPOP::GenoStruTrait::absLocusIndex "

Usage:

    x.absLocusIndex(chrom, locus)

Details:

    return the absolute index of locus locus on chromosome chrom. An
    IndexError will be raised if chrom or locus is out of range. c.f.
    chromLocusPair.

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromLocusPair "

Usage:

    x.chromLocusPair(locus)

Details:

    return the chromosome and relative index of a locus using its
    absolute index locus. c.f. absLocusIndex.

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromName "

Usage:

    x.chromName(chrom)

Details:

    return the name of a chromosome chrom. Default to chrom# where #
    is the 1-based index of the chromosome.

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromNames "

Usage:

    x.chromNames()

Details:

    return a list of the names of all chromosomes.

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromType "

Usage:

    x.chromType(chrom)

Details:

    return the type of a chromosome chrom (Customized, Autosome,
    ChromosomeX, or ChromosomeY).

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromTypes "

Usage:

    x.chromTypes()

Details:

    return the type of all chromosomes (Customized, Autosome,
    ChromosomeX or ChromosomeY).

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromByName "

Usage:

    x.chromByName(name)

Details:

    return the index of a chromosome by its name.

"; 

%feature("docstring") simuPOP::GenoStruTrait::alleleName "

Usage:

    x.alleleName(allele)

Details:

    return the name of allele allele specified by the alleleNames
    parameter of the population function. If the name of an allele is
    not specified, its index ('0', '1', '2', etc) is returned. An
    IndexError will be raised if allele is larger than the maximum
    allowed allele state of this module ( MaxAllele()).

"; 

%feature("docstring") simuPOP::GenoStruTrait::alleleNames "

Usage:

    x.alleleNames()

Details:

    return a list of allele names given by the alleleNames parameter
    of the population function. This list does not have to cover all
    possible allele states of a population so alleleNames()[allele]
    might fail (use alleleNames(allele) instead).

"; 

%feature("docstring") simuPOP::GenoStruTrait::locusName "

Usage:

    x.locusName(loc)

Details:

    return the name of locus loc specified by the lociNames parameter
    of the population function. Default to locX-Y where X and Y are
    1-based chromosome and locus indexes (loc1-1, loc1-2, ... etc)

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociNames "

Usage:

    x.lociNames()

Details:

    return the names of all loci specified by the lociNames parameter
    of the population function.

"; 

%feature("docstring") simuPOP::GenoStruTrait::locusByName "

Usage:

    x.locusByName(name)

Details:

    return the index of a locus with name name. Raise a ValueError if
    no locus is found.

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociByNames "

Usage:

    x.lociByNames(names)

Details:

    return the indexes of loci with names names. Raise a ValueError if
    any of the loci cannot be found.

"; 

%feature("docstring") simuPOP::GenoStruTrait::hasInfoField "Obsolete or undocumented function."

%feature("docstring") simuPOP::GenoStruTrait::infoSize "Obsolete or undocumented function."

%feature("docstring") simuPOP::GenoStruTrait::infoFields "

Usage:

    x.infoFields()

Details:

    return a list of the names of all information fields of the
    population.

"; 

%feature("docstring") simuPOP::GenoStruTrait::infoField "

Usage:

    x.infoField(idx)

Details:

    return the name of information field idx.

"; 

%feature("docstring") simuPOP::GenoStruTrait::infoIdx "

Usage:

    x.infoIdx(name)

Details:

    return the index of information field name. Raise an IndexError if
    name is not one of the information fields.

"; 

%ignore simuPOP::GenoStruTrait::struAddInfoFields(const vectorstr &fields);

%ignore simuPOP::GenoStruTrait::struSetInfoFields(const vectorstr &fields);

%ignore simuPOP::GenoStruTrait::swap(GenoStruTrait &rhs);

%feature("docstring") simuPOP::genotypeSplitter "

Details:

    This class defines a VSP splitter that defines VSPs according to
    individual genotype at specified loci.

"; 

%feature("docstring") simuPOP::genotypeSplitter::genotypeSplitter "

Usage:

    genotypeSplitter(loci (or locus), alleles, phase=False)

Details:

    Create a splitter that defined VSPs by individual genotype at loci
    loci (or locus if only one locus is used). Each list in a list
    allele defines a VSP, which is a list of allowed alleles at these
    loci. If only one VSP is defined, the outer list of the nested
    list can be ignored. If phase if true, the order of alleles in
    each list is significant. If more than one set of alleles are
    given, individuals having either of them is qualified.
    For example, in a haploid population, locus=1, alleles=[0, 1]
    defines a VSP with individuals having allele 0 or 1 at locus 1,
    alleles=[[0, 1], [2]] defines two VSPs with indivdiuals in the
    second VSP having allele 2 at locus 1. If multiple loci are
    involved, alleles at each locus need to be defined. For example,
    VSP defined by loci=[0, 1], alleles=[0, 1, 1, 1] consists of
    individuals having alleles [0, 1] or [1, 1] at loci [0, 1].
    In a haploid population, locus=1, alleles=[0, 1] defines a VSP
    with individuals having genotype [0, 1] or [1, 0] at locus 1.
    alleles[[0, 1], [2, 2]] defines two VSPs with indivdiuals in the
    second VSP having genotype [2, 2] at locus 1. If phase is set to
    True, the first VSP will only has individuals with genotype [0,
    1]. In the multiple loci case, alleles should be arranged by
    haplotypes, for example, loci=[0, 1], alleles=[0, 0, 1, 1],
    phase=True defines a VSP with individuals having genotype -0-0-,
    -1-1- at loci 0 and 1. If phase=False (default), genotypes -1-1-,
    -0-0-, -0-1- and -1-0- are all allowed.

"; 

%feature("docstring") simuPOP::genotypeSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::genotypeSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::genotypeSplitter::numVirtualSubPop "

Description:

    number of virtual subpops of subpopulation sp

Usage:

    x.numVirtualSubPop()

"; 

%ignore simuPOP::genotypeSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::genotypeSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::genotypeSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return name of VSP vsp, which is \"Genotype loc1,loc2:genotype\" as
    defined by parameters loci and alleles.

"; 

%feature("docstring") simuPOP::gsmMutator "

Function form:

    GsmMutate

Description:

    generalized stepwise mutation model

Details:

    The Generalized Stepwise Mutation model (GSM) is an extension to
    the stepwise mutation model. This model assumes that alleles are
    represented by integer values and that a mutation either increases
    or decreases the allele value by a random value. In other words,
    in this model the change in the allelic state is drawn from a
    random distribution. A geometric generalized stepwise model uses a
    geometric distribution with parameter $ p $, which has mean $
    \\frac{p}{1-p} $ and variance $ \\frac{p}{\\left(1-p\\right)^{2}} $.
    gsmMutator implements both models. If you specify a Python
    function without a parameter, this mutator will use its return
    value each time a mutation occur; otherwise, a parameter $ p $
    should be provided and the mutator will act as a geometric
    generalized stepwise model.

"; 

%feature("docstring") simuPOP::gsmMutator::gsmMutator "

Description:

    create a gsmMutator

Usage:

    gsmMutator(rate=[], loci=[], maxAllele=0, incProb=0.5, p=0,
      func=None, output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
      end=-1, step=1, at=[], rep=[], subPop=[], infoFields=[])

Details:

    The GSM model is developed for allozymes. It provides better
    description for these kinds of evolutionary processes.  Please see
    class mutator for the descriptions of other parameters.

Arguments:

    incProb:        probability to increase allele state. Default to
                    0.5.
    func:           a function that returns the number of steps. This
                    function does not accept any parameter.

"; 

%feature("docstring") simuPOP::gsmMutator::~gsmMutator "

Usage:

    x.~gsmMutator()

"; 

%feature("docstring") simuPOP::gsmMutator::clone "

Description:

    deep copy of a gsmMutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::gsmMutator::mutate "

Description:

    mutate according to the GSM model

Usage:

    x.mutate(allele)

"; 

%feature("docstring") simuPOP::gsmMutator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the gsmMutator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::haplodiploidGenoTransmitter "

Applicability: haplodiploid only

Details:

    haplodiploid offspring generator mimics sex-determination in honey
    bees. Given a female (queen) parent and a male parent, the female
    is considered as diploid with two set of chromosomes, and the male
    is condiered as haploid. Actually, the first set of male
    chromosomes are used. During mating, female produce eggs, subject
    to potential recombination and gene conversion, while male sperm
    is identical to the parental chromosome. Female offspring has two
    sets of chromosomes, one from mother and one from father. Male
    offspring has one set of chromosomes from his mother.

"; 

%feature("docstring") simuPOP::haplodiploidGenoTransmitter::haplodiploidGenoTransmitter "

Usage:

    haplodiploidGenoTransmitter(begin=0, end=-1, step=1, at=[],
      rep=[], subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::haplodiploidGenoTransmitter::clone "

Usage:

    x.clone()

Details:

    Return a cloned copy of an operator. This function is available to
    all operators.

"; 

%feature("docstring") simuPOP::haplodiploidGenoTransmitter::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::haplodiploidGenoTransmitter::initialize "

Usage:

    x.initialize(pop)

"; 

%ignore simuPOP::haplodiploidGenoTransmitter::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::heteroMating "

Applicability: diploid only

Details:

    a heterogeneous mating scheme that applies a list of mating
    schemes to different (virtual) subpopulations.

"; 

%feature("docstring") simuPOP::heteroMating::heteroMating "

Description:

    create a heterogeneous Python mating scheme

Usage:

    heteroMating(matingSchemes, newSubPopSize=[],
      newSubPopSizeExpr=\"\", newSubPopSizeFunc=None,
      shuffleOffspring=True, subPop=[], weight=0)

Details:

    Parameter subpop, virtualSubPOp and weight of this mating scheme
    is ignored.

Arguments:

    matingSchemes:  A list of mating schemes. If parameter subPop of
                    an mating scheme is specified, it will be applied
                    to specific subpopulation. If virtualSubPop if
                    specified, it will be applied to specifc virtual
                    subpopulations.

"; 

%feature("docstring") simuPOP::heteroMating::~heteroMating "

Description:

    destructor

Usage:

    x.~heteroMating()

"; 

%ignore simuPOP::heteroMating::heteroMating(const heteroMating &rhs);

%feature("docstring") simuPOP::heteroMating::clone "

Description:

    deep copy of a Python mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::heteroMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the Python mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::heteroMating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%feature("docstring") simuPOP::ifElse "

Description:

    conditional operator

Details:

    This operator accepts
    *  an expression that will be evaluated when this operator is
    applied.
    *  an operator that will be applied if the expression is True
    (default to null).
    *  an operator that will be applied if the expression is False
    (default to null). When this operator is applied to a population,
    it will evaluate the expression and depending on its value, apply
    the supplied operator. Note that the begin, end, step, and at
    parameters of ifOp and elseOp will be ignored. For example, you
    can mimic the at parameter of an operator by ifElse('rep in
    [2,5,9]' operator). The real use of this machanism is to monitor
    the population statistics and act accordingly.

"; 

%feature("docstring") simuPOP::ifElse::ifElse "

Description:

    create a conditional operator

Usage:

    ifElse(cond, ifOp=None, elseOp=None, output=\">\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Arguments:

    cond:           expression that will be treated as a boolean
                    variable
    ifOp:           an operator that will be applied when cond is True
    elseOp:         an operator that will be applied when cond is
                    False

"; 

%feature("docstring") simuPOP::ifElse::~ifElse "

Description:

    destructor

Usage:

    x.~ifElse()

"; 

%ignore simuPOP::ifElse::ifElse(const ifElse &rhs);

%feature("docstring") simuPOP::ifElse::clone "Obsolete or undocumented function."

%ignore simuPOP::ifElse::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::ifElse::apply "

Description:

    apply the ifElse operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::ifElse::__repr__ "

Description:

    used by Python print function to print out the general information
    of the ifElse operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::indCompare;

%feature("docstring") simuPOP::indCompare::indCompare "

Usage:

    indCompare(idx)

"; 

%feature("docstring") simuPOP::IndexError "

Description:

    exception, thrown if index out of range

"; 

%feature("docstring") simuPOP::IndexError::IndexError "

Usage:

    IndexError(msg)

"; 

%feature("docstring") simuPOP::individual "

Details:

    A population consists of individuals with the same genotypic
    structure. An individual object cannot be created independently,
    but refences to inidividuals can be retrieved using member
    functions of a population object. In addition to structural
    information shared by all individuals in a population (provided by
    class genoStruTrait), the individual class provides member
    functions to get and set genotype, sex, affection status and
    information fields of an individual.  Genotypes of an individual
    are stored sequentially and can be accessed locus by locus, or in
    batch. The alleles are arranged by position, chromosome and
    ploidy. That is to say, the first allele on the first chromosome
    of the first homologous set is followed by alleles at other loci
    on the same chromsome, then markers on the second and later
    chromosomes, followed by alleles on the second homologous set of
    the chromosomes for a diploid individual. A consequence of this
    memory layout is that alleles at the same locus of a non-haploid
    individual are separated by individual::totNumLoci() loci. It is
    worth noting that access to invalid chromosomes, such as the Y
    chromosomes of female individuals, are not restricted.

"; 

%feature("docstring") simuPOP::individual::individual "

Usage:

    individual()

Details:

    An individual object cannot be created directly. It has to be
    accessed from a population object using functions such as
    population::individual(idx).

"; 

%ignore simuPOP::individual::individual(const individual &ind);

%feature("docstring") simuPOP::individual::~individual "

Description:

    destructor. Do nothing.

Usage:

    x.~individual()

"; 

%ignore simuPOP::individual::setGenoPtr(GenoIterator pos);

%ignore simuPOP::individual::setInfoPtr(InfoIterator pos);

%ignore simuPOP::individual::copyFrom(const individual &rhs);

%ignore simuPOP::individual::genoPtr() const;

%ignore simuPOP::individual::infoPtr() const;

%feature("docstring") simuPOP::individual::allele "

Usage:

    x.allele(idx)

Details:

    return the current allele at a locus, using its absolute index
    idx.

"; 

%feature("docstring") simuPOP::individual::allele "

Usage:

    x.allele(idx, p)

Details:

    return the current allele at locus idx on the p-th set of
    homologous chromosomes.

"; 

%feature("docstring") simuPOP::individual::allele "

Usage:

    x.allele(idx, p, chrom)

Details:

    return the current allele at locus idx on chromosome chrom of the
    p-th set of homologous chromosomes.

"; 

%feature("docstring") simuPOP::individual::alleleChar "Obsolete or undocumented function."

%feature("docstring") simuPOP::individual::alleleChar "Obsolete or undocumented function."

%feature("docstring") simuPOP::individual::alleleChar "Obsolete or undocumented function."

%feature("docstring") simuPOP::individual::setAllele "

Usage:

    x.setAllele(allele, idx)

Details:

    set allele allele to a locus, using its absolute index idx.

"; 

%feature("docstring") simuPOP::individual::setAllele "

Usage:

    x.setAllele(allele, idx, p)

Details:

    set allele allele to locus idx on the p-th homologous set of
    chromosomes.

"; 

%feature("docstring") simuPOP::individual::setAllele "

Usage:

    x.setAllele(allele, idx, p, chrom)

Details:

    set allele allele to locus idx on chromosome chrom of the p-th
    homologous set of chromosomes.

"; 

%feature("docstring") simuPOP::individual::genotype "

Usage:

    x.genotype()

Details:

    return an editable array (a carray of length
    totNumLoci()*ploidy()) that represents all alleles of an
    individual.

"; 

%feature("docstring") simuPOP::individual::genotype "

Usage:

    x.genotype(p)

Details:

    return an editable array (a carray of length totNumLoci()) that
    represents all alleles on the p-th homologous set of chromosomes.

"; 

%feature("docstring") simuPOP::individual::genotype "

Usage:

    x.genotype(p, chrom)

Details:

    return an editable array (a carrary of legnth numLoci(chrom)) that
    represents all alleles on chromosome chrom of the p-th homologous
    set of chromosomes.

"; 

%feature("docstring") simuPOP::individual::setGenotype "

Usage:

    x.setGenotype(geno)

Details:

    Fill the genotype of an individual using a list of alleles geno.
    geno will be reused if its length is less than
    totNumLoci()*ploidy().

"; 

%feature("docstring") simuPOP::individual::setGenotype "

Usage:

    x.setGenotype(geno, p)

Details:

    Fill the genotype of the p-th homologous set of chromosomes using
    a list of alleles geno. geno will be reused if its length is less
    than totNumLoci().

"; 

%feature("docstring") simuPOP::individual::setGenotype "

Usage:

    x.setGenotype(geno, p, chrom)

Details:

    Fill the genotype of chromosome chrom on the p-th homologous set
    of chromosomes using a list of alleles geno. geno will be reused
    if its length is less than mumLoci(chrom).

"; 

%feature("docstring") simuPOP::individual::sex "

Usage:

    x.sex()

Details:

    return the sex of an individual, 1 for male and 2 for female.

"; 

%feature("docstring") simuPOP::individual::sexChar "

Usage:

    x.sexChar()

Details:

    return the sex of an individual, M for male or F for female.

"; 

%feature("docstring") simuPOP::individual::setSex "

Usage:

    x.setSex(sex)

Details:

    set individual sex to Male or Female.

"; 

%feature("docstring") simuPOP::individual::affected "

Usage:

    x.affected()

Details:

    Return True if this individual is affected.

"; 

%feature("docstring") simuPOP::individual::affectedChar "

Usage:

    x.affectedChar()

Details:

    Return A if this individual is affected, or U otherwise.

"; 

%feature("docstring") simuPOP::individual::setAffected "

Usage:

    x.setAffected(affected)

Details:

    set affection status to affected (True or False).

"; 

%ignore simuPOP::individual::iteratable() const;

%ignore simuPOP::individual::setIteratable(bool iteratable);

%ignore simuPOP::individual::visible() const;

%ignore simuPOP::individual::setVisible(bool visible);

%feature("docstring") simuPOP::individual::info "

Usage:

    x.info(idx)

Details:

    Return the value of an information field idx (an index).

"; 

%feature("docstring") simuPOP::individual::intInfo "

Usage:

    x.intInfo(idx)

Details:

    Return the value of an information field idx (an index) as an
    integer number.

"; 

%feature("docstring") simuPOP::individual::info "

Usage:

    x.info(name)

Details:

    Return the value of an information field name.

"; 

%feature("docstring") simuPOP::individual::intInfo "

Usage:

    x.intInfo(name)

Details:

    Return the value of an information field name as an integer
    number.

"; 

%feature("docstring") simuPOP::individual::setInfo "

Usage:

    x.setInfo(value, idx)

Details:

    set the value of an information field idx (an index) to value.

"; 

%feature("docstring") simuPOP::individual::setInfo "

Usage:

    x.setInfo(value, name)

Details:

    set the value of an information field name to value.

"; 

%ignore simuPOP::individual::genoBegin() const;

%ignore simuPOP::individual::genoEnd() const;

%ignore simuPOP::individual::genoBegin(UINT p) const;

%ignore simuPOP::individual::genoEnd(UINT p) const;

%ignore simuPOP::individual::genoBegin(UINT p, UINT chrom) const;

%ignore simuPOP::individual::infoBegin() const;

%ignore simuPOP::individual::infoEnd() const;

%feature("docstring") simuPOP::individual::__cmp__ "

Description:

    a python function used to compare the individual objects

Usage:

    x.__cmp__(rhs)

"; 

%feature("docstring") simuPOP::individual::__repr__ "

Description:

    used by Python print function to print out the general information
    of the individual

Usage:

    x.__repr__()

"; 

%ignore simuPOP::individual::swap(individual &ind, bool swapContent=true);

%ignore simuPOP::individual::display(ostream &out, int width=1, const vectori &chrom=vectori(), const vectori &loci=vectori());

%feature("docstring") simuPOP::IndividualIterator "

Details:

    this class implements a C++ iterator class that iterate through
    individuals in a (sub)population. If allInds are true, the
    visiblility of individuals will not be checked. Note that
    individualIterator *will* iterate through only visible
    individuals, and allInds is only provided when we know in advance
    that all individuals are visible. This is a way to obtain better
    performance in simple cases.

"; 

%feature("docstring") simuPOP::IndividualIterator::IndividualIterator "

Usage:

    IndividualIterator()

"; 

%feature("docstring") simuPOP::IndividualIterator::valid "

Usage:

    x.valid()

"; 

%feature("docstring") simuPOP::IndividualIterator::rawIter "

Usage:

    x.rawIter()

"; 

%feature("docstring") simuPOP::infoEval "

Function form:

    infoEval

Details:

    Unlike operator pyEval and pyExec that work at the population
    level, in its local namespace, infoEval works at the individual
    level, working with individual information fields. is statement
    can change the value of existing information fields. Optionally,
    variables in population's local namespace can be used in the
    statement, but this should be used with caution.

"; 

%feature("docstring") simuPOP::infoEval::infoEval "

Description:

    evaluate Python statements with variables being an individual's
    information fields

Usage:

    infoEval(expr=\"\", stmts=\"\", subPops=[], usePopVars=False,
      exposePop=False, name=\"\", output=\">\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Details:

    The expression and statements will be executed for each
    individual, in a Python namespace (dictionary) where individual
    information fields are made available as variables. Population
    dictionary can be made avaialbe with option usePopVars. Changes to
    these variables will change the corresponding information fields
    of individuals. Please note that, 1. If population variables are
    used, and there are name conflicts between information fields and
    variables, population variables will be overridden by information
    fields, without any warning. 2. Information fields are float
    numbers. An exceptions will raise if an information field can not
    be converted to a float number. 3. This operator can be used in
    all stages. When it is used during-mating, it will act on each
    offspring.

Arguments:

    expr:           the expression to be evaluated. The result will be
                    sent to output.
    stmts:          the statement that will be executed before the
                    expression
    subPop:         a shortcut to subPops=[subPop]
    subPops:        subpopulations this operator will apply to.
                    Default to all.
    usePopVars:     if True, import variables from expose the current
                    population as a variable named pop
    exposePop:      if True, expose the current population as a
                    variable named pop
    name:           used to let pure Python operator to identify
                    themselves
    output:         default to >. I.e., output to standard output.
                    Note that because the expression will be executed
                    for each individual, the output can be large.

"; 

%feature("docstring") simuPOP::infoEval::~infoEval "

Usage:

    x.~infoEval()

"; 

%feature("docstring") simuPOP::infoEval::clone "

Description:

    deep copy of a infoEval operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::infoEval::apply "

Description:

    apply the infoEval operator

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::infoEval::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::infoEval::__repr__ "

Description:

    used by Python print function to print out the general information
    of the infoEval operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::infoEval::name "

Description:

    return the name of an expression

Usage:

    x.name()

Details:

    The name of a infoEval operator is given by an optional parameter
    name. It can be used to identify this infoEval operator in debug
    output, or in the dryrun mode of simulator::evolve.

"; 

%feature("docstring") simuPOP::infoExec "

Function form:

    infoExec

Description:

    execute a Python statement for each individual, using information
    fields

Details:

    This operator takes a list of statements and executes them. No
    value will be returned or outputted.

"; 

%feature("docstring") simuPOP::infoExec::infoExec "

Description:

    fields, optionally with variable in population's local namespace

Usage:

    infoExec(stmts=\"\", subPops=[], usePopVars=False,
      exposePop=False, name=\"\", output=\">\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Details:

    Please refer to class infoEval for parameter descriptions.

"; 

%feature("docstring") simuPOP::infoExec::~infoExec "

Usage:

    x.~infoExec()

"; 

%feature("docstring") simuPOP::infoExec::clone "

Description:

    deep copy of a infoExec operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::infoExec::__repr__ "

Description:

    used by Python print function to print out the general information
    of the infoExec operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::infoParentsChooser "

Applicability: all ploidy

Details:

    This parents chooser choose an individual randomly, but choose
    his/her spouse from a given set of information fields, which
    stores indexes of individuals in the same generation. A field will
    be ignored if its value is negative, or if sex is compatible.
    Depending on what indexes are stored in these information fields,
    this parent chooser can be used to implement consanguineous mating
    where close relatives are located for each individual, or certain
    non-random mating schemes where each individual can only mate with
    a small number of pre-determinable individuals. This parent
    chooser (currently) uses randomParentChooser to choose one parent
    and randomly choose another one from the information fields.
    Because of potentially non-even distribution of valid information
    fields, the overall process may not be as random as expected,
    especially when selection is applied. Note: if there is no valid
    individual, this parents chooser works like a double
    parentChooser.

"; 

%feature("docstring") simuPOP::infoParentsChooser::infoParentsChooser "

Usage:

    infoParentsChooser(infoFields=[], func=None, param=None,
      replacement=True)

Arguments:

    infoFields:     information fields that store index of matable
                    individuals.

"; 

%feature("docstring") simuPOP::infoParentsChooser::~infoParentsChooser "

Usage:

    x.~infoParentsChooser()

"; 

%feature("docstring") simuPOP::infoParentsChooser::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::infoParentsChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::infoParentsChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::InformationIterator "

Details:

    this class implements a C++ iterator class that iterate through
    infomation fields in a (sub)population using 1. an IndIterator
    that will skip invisible individuals, or 2. a gapped iterator that
    will run faster. Note that 1, 2 should yield identical result, and
    2 should be used when there is no virtual subpopulation.q

"; 

%feature("docstring") simuPOP::InformationIterator::InformationIterator "

Usage:

    InformationIterator()

"; 

%feature("docstring") simuPOP::infoSplitter "

Details:

    This splitter defines VSPs according to the value of an
    information field of each indivdiual. A VSP is defined either by a
    value or a range of values.

"; 

%feature("docstring") simuPOP::infoSplitter::infoSplitter "

Usage:

    infoSplitter(field, values=[], cutoff=[])

Details:

    Create an infomration splitter using information field field. If
    parameter values is specified, each item in this list defines a
    VSP in which all individuals have this value at information field
    field. If a set of cutoff values are defined in parameter cutoff,
    individuals are grouped by intervals defined by these cutoff
    values. For example, cutoff=[1,2] defines three VSPs with v < 1, 1
    <= v < 2 and v >=2 where v is the value of an individual at
    information field field. Of course, only one of the parameters
    values and cutoff should be defined, values in cutoff should be
    distinct, and in an increasing order.

"; 

%feature("docstring") simuPOP::infoSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::infoSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::infoSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs defined by this splitter, which is the
    length parameter values or the length of cutoff plus one,
    depending on which parameter is specified.

"; 

%ignore simuPOP::infoSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::infoSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::infoSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of a VSP vsp, which is field = value if VSPs are
    defined by values in parameter values, or field < value (the first
    VSP), v1 <= field < v2 and field >= v (the last VSP) if VSPs are
    defined by cutoff values.

"; 

%feature("docstring") simuPOP::inheritTagger "

Description:

    inherite tag from parents

Details:

    This during-mating operator will copy the tag (information field)
    from his/her parents. Depending on mode parameter, this tagger
    will obtain tag, value of the first specified information fields,
    from his/her father or mother (two tag fields), or both (first tag
    field from father, and second tag field from mother).
    An example may be tagging one or a few parents and examining, at
    the last generation, how many offspring they have.

"; 

%feature("docstring") simuPOP::inheritTagger::inheritTagger "

Description:

    create an inheritTagger that inherits a tag from one or both
    parents

Usage:

    inheritTagger(mode=TAG_Paternal, begin=0, end=-1, step=1, at=[],
      rep=[], subPop=[], output=\"\", outputExpr=\"\",
      infoFields=[\"paternal_tag\", \"maternal_tag\"])

Arguments:

    mode:           can be one of TAG_Paternal, TAG_Maternal, and
                    TAG_Both

"; 

%feature("docstring") simuPOP::inheritTagger::~inheritTagger "

Usage:

    x.~inheritTagger()

"; 

%feature("docstring") simuPOP::inheritTagger::__repr__ "

Description:

    used by Python print function to print out the general information
    of the inheritTagger

Usage:

    x.__repr__()

"; 

%ignore simuPOP::inheritTagger::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::inheritTagger::clone "

Description:

    deep copy of a inheritTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initByFreq "

Function form:

    InitByFreq

Details:

    This operator assigns alleles at all or part of loci with given
    allele frequencies. Alternatively, an individual can be
    initialized and be copied to all individuals in the same (virtual)
    subpopulations.

"; 

%feature("docstring") simuPOP::initByFreq::initByFreq "

Usage:

    initByFreq(alleleFreq=[], loci=[], ploidy=[],
      identicalInds=False, initSex=True, maleFreq=0.5, sex=[],
      stage=PreMating, begin=0, end=1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Details:

    This function creates an initializer that initialize individual
    genotypes randomly. alleleFreq specified the allele frequencies of
    allele 0, 1, ... respectively. These frequencies should add up to
    1. If loci, ploidy and/or subPop are specified, only specified
    loci, ploidy, and individuals in these (virtual) subpopulations
    will be initialized. If identicalInds is True, the first
    individual in each (virtual) subpopulation will be initialized
    randomly, and be copied to all other individuals in this (virtual)
    subpopulation. If a list of frequencies are given, they will be
    used for each (virtual) subpopulation. If initSex is True
    (default), initSex(maleFreq, sex) will be applied. This operator
    initializes all chromosomes, including unused genotype locations
    and customized chromosomes.

"; 

%feature("docstring") simuPOP::initByFreq::~initByFreq "

Usage:

    x.~initByFreq()

"; 

%feature("docstring") simuPOP::initByFreq::clone "

Description:

    deep copy of the operator initByFreq

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initByFreq::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator initByFreq

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::initByFreq::apply "

Description:

    apply this operator to population pop

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::initByValue "

Function form:

    InitByValue

Details:

    This operator initialize individuals by given values.

"; 

%feature("docstring") simuPOP::initByValue::initByValue "

Usage:

    initByValue(value=[], loci=[], ploidy=[], proportions=[],
      initSex=True, maleFreq=0.5, sex=[], stage=PreMating, begin=0,
      end=1, step=1, at=[], rep=[], subPop=[], infoFields=[])

Details:

    This function creates an initializer that initialize individual
    genotypes with given genotype value. If loci, ploidy and/or subPop
    are specified, only specified loci, ploidy, and individuals in
    these (virtual) subpopulations will be initialized. value can be
    used to initialize given loci, all loci, and all homologous copies
    of these loci. If proportions (a list of positive numbers that add
    up to 1) is given, value should be a list of values that will be
    assigned randomly according to their respective proportion. If a
    list of values are given without proportions, they will be used
    for each (virtual) subpopulations. If initSex is True (default),
    initSex(maleFreq, sex) will be applied. This operator initializes
    all chromosomes, including unused genotype locations and
    customized chromosomes.

"; 

%feature("docstring") simuPOP::initByValue::~initByValue "

Usage:

    x.~initByValue()

"; 

%feature("docstring") simuPOP::initByValue::clone "

Description:

    deep copy of the operator initByValue

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initByValue::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator initByValue

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::initByValue::apply "

Description:

    apply this operator to population pop

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::initSex "

Function form:

    InitSex

Details:

    This operator initialize sex of individuals, either randomly or
    use a list of sexes. For convenience, the function of this
    operator is included in other initializers such as initByFreq and
    initByValue so that you do not have to intiailize sexes separately
    from genotype.

"; 

%feature("docstring") simuPOP::initSex::initSex "

Usage:

    initSex(maleFreq=0.5, sex=[], stage=PreMating, begin=0, end=-1,
      step=1, at=[], rep=[], subPop=[], infoFields=[])

Details:

    Create an operator that initialize individual sex to Male or
    Female. By default, it assign sex to individuals randomly, with
    equal probability of having a male or a female. This probabability
    can be adjusted through parameter maleFreq. Alternatively, a fixed
    sequence of sexes can be assigned. For example, if sex=[Male,
    Female], individuals will be assigned Male and Female
    successively. Parameter maleFreq is ignored if sex is given. If a
    list of (virtual) subpopulation is specified in parameter subPop,
    only individuals in these subpopulations will be initialized.

"; 

%feature("docstring") simuPOP::initSex::~initSex "

Description:

    destructor

Usage:

    x.~initSex()

"; 

%feature("docstring") simuPOP::initSex::clone "

Description:

    deep copy of an initSex

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initSex::__repr__ "

Description:

    used by Python print function to print out the general information
    of the initSex

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::initSex::apply "

Description:

    apply this operator to population pop

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::IOError "

Description:

    exception, thrown if file io failure

"; 

%feature("docstring") simuPOP::IOError::IOError "

Usage:

    IOError(msg)

"; 

%ignore simuPOP::isAffected;

%feature("docstring") simuPOP::isAffected::isAffected "

Usage:

    isAffected()

"; 

%feature("docstring") simuPOP::kamMutator "

Function form:

    KamMutate

Description:

    K-Allele Model mutator.

Details:

    This mutator mutate an allele to another allelic state with equal
    probability. The specified mutation rate is actually the
    'probability to mutate'. So the mutation rate to any other allelic
    state is actually $ \\frac{rate}{K-1} $, where $ K $ is specified
    by parameter maxAllele.

"; 

%feature("docstring") simuPOP::kamMutator::kamMutator "

Description:

    create a K-Allele Model mutator

Usage:

    kamMutator(rate=[], loci=[], maxAllele=0, output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=[], subPop=[], infoFields=[])

Details:

    Please see class mutator for the descriptions of other parameters.

Arguments:

    rate:           mutation rate. It is the 'probability to mutate'.
                    The actual mutation rate to any of the other K-1
                    allelic states are rate/(K-1).
    maxAllele:      maximum allele that can be mutated to. For binary
                    libraries, allelic states will be [0, maxAllele].
                    Otherwise, they are [1, maxAllele].

"; 

%feature("docstring") simuPOP::kamMutator::~kamMutator "

Usage:

    x.~kamMutator()

"; 

%feature("docstring") simuPOP::kamMutator::mutate "

Description:

    mutate to a state other than current state with equal probability

Usage:

    x.mutate(allele)

"; 

%feature("docstring") simuPOP::kamMutator::clone "

Description:

    deep copy of a kamMutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::kamMutator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the kamMutator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::maPenetrance "

Function form:

    MaPenetrance

Description:

    multiple allele penetrance operator

Details:

    This is called 'multiple-allele' penetrance. It separates alleles
    into two groups: wildtype and diseased alleles. Wildtype alleles
    are specified by parameter wildtype and any other alleles are
    considered as diseased alleles. maPenetrance accepts an array of
    penetrance for AA, Aa, aa in the single-locus case, and a longer
    table for the multi-locus case. Penetrance is then set for any
    given genotype.

"; 

%feature("docstring") simuPOP::maPenetrance::maPenetrance "

Description:

    create a multiple allele penetrance operator (penetrance according
    to diseased or wildtype alleles)

Usage:

    maPenetrance(loci, penet, wildtype, ancestralGen=-1,
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Arguments:

    locus:          the locus index. The genotype of this locus will
                    be used to determine penetrance.
    loci:           the locus indexes. The genotypes of these loci
                    will be examed.
    penet:          an array of penetrance values of AA, Aa, aa. A is
                    the wild type group. In the case of multiple loci,
                    penetrance should be in the order of AABB, AABb,
                    AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb.
    wildtype:       an array of alleles in the wildtype group. Any
                    other alleles will be considered as in the
                    diseased allele group.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::maPenetrance::~maPenetrance "

Usage:

    x.~maPenetrance()

"; 

%feature("docstring") simuPOP::maPenetrance::clone "

Description:

    deep copy of a multi-allele penetrance operator

Usage:

    x.clone()

"; 

%ignore simuPOP::maPenetrance::penet(individual *ind);

%feature("docstring") simuPOP::maPenetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the multi-allele penetrance operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mapPenetrance "

Function form:

    MapPenetrance

Description:

    penetrance according to the genotype at one locus

Details:

    Assign penetrance using a table with keys 'X-Y' where X and Y are
    allele numbers.

"; 

%feature("docstring") simuPOP::mapPenetrance::mapPenetrance "

Description:

    create a map penetrance operator

Usage:

    mapPenetrance(loci, penet, phase=False, ancestralGen=-1,
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Arguments:

    locus:          the locus index. Shortcut to loci=[locus]
    loci:           the locus indexes. The genotypes of these loci
                    will be used to determine penetrance.
    penet:          a dictionary of penetrance. The genotype must be
                    in the form of 'a-b' for a single locus.
    phase:          if True, a/b and b/a will have different
                    penetrance values. Default to False.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::mapPenetrance::~mapPenetrance "

Usage:

    x.~mapPenetrance()

"; 

%feature("docstring") simuPOP::mapPenetrance::clone "

Description:

    deep copy of a map penetrance operator

Usage:

    x.clone()

"; 

%ignore simuPOP::mapPenetrance::penet(individual *ind);

%feature("docstring") simuPOP::mapPenetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the map penetrance operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mapQuanTrait "

Function form:

    MapQuanTrait

Description:

    quantitative trait according to genotype at one locus

Details:

    Assign quantitative trait using a table with keys 'X-Y' where X
    and Y are allele numbers. If parameter sigma is not zero, the
    return value is the sum of the trait plus $
    N\\left(0,\\sigma^{2}\\right) $. This random part is usually
    considered as the environmental factor of the trait.

"; 

%feature("docstring") simuPOP::mapQuanTrait::mapQuanTrait "

Description:

    create a map quantitative trait operator

Usage:

    mapQuanTrait(loci, qtrait, sigma=0, phase=False,
      ancestralGen=-1, stage=PostMating, begin=0, end=-1, step=1,
      at=[], rep=[], subPop=[], infoFields=[\"qtrait\"])

Arguments:

    locus:          the locus index. The quantitative trait is
                    determined by genotype at this locus.
    loci:           an array of locus indexes. The quantitative trait
                    is determined by genotypes at these loci.
    qtrait:         a dictionary of quantitative traits. The genotype
                    must be in the form of 'a-b'. This is the mean of
                    the quantitative trait. The actual trait value
                    will be $ N\\left(mean,\\sigma^{2}\\right) $. For
                    multiple loci, the form is 'a-b|c-d|e-f' etc.
    sigma:          standard deviation of the environmental factor $
                    N\\left(0,\\sigma^{2}\\right) $.
    phase:          if True, a/b and b/a will have different
                    quantitative trait values. Default to False.
    output:         and other parameters please refer to help
                    (baseOperator.__init__)

"; 

%feature("docstring") simuPOP::mapQuanTrait::~mapQuanTrait "

Usage:

    x.~mapQuanTrait()

"; 

%feature("docstring") simuPOP::mapQuanTrait::clone "

Description:

    deep copy of a map quantitative trait operator

Usage:

    x.clone()

"; 

%ignore simuPOP::mapQuanTrait::qtrait(individual *ind);

%feature("docstring") simuPOP::mapQuanTrait::__repr__ "

Description:

    used by Python print function to print out the general information
    of the map quantitative trait operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mapSelector "

Function form:

    MapSelector

Applicability: all ploidy

Description:

    selection according to the genotype at one or more loci

Details:

    This map selector implements selection according to genotype at
    one or more loci. A user provided dictionary (map) of genotypes
    will be used in this selector to set each individual's fitness
    value.

"; 

%feature("docstring") simuPOP::mapSelector::mapSelector "

Description:

    create a map selector

Usage:

    mapSelector(loci, fitness, phase=False, subPops=[],
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[\"fitness\"])

Arguments:

    locus:          the locus index. A shortcut to loci=[locus]
    loci:           the locus indexes. The genotypes at these loci
                    will be used to determine the fitness value.
    fitness:        a dictionary of fitness values. The genotype must
                    be in the form of 'a-b' for a single locus, and
                    'a-b|c-d|e-f' for multi-loci. In the haploid case,
                    the genotype should be specified in the form of
                    'a' for single locus, and 'a|b|c' for multi-locus
                    models.
    phase:          if True, genotypes a-b and b-a will have different
                    fitness values. Default to False.
    output:         and other parameters please refer to help
                    (baseOperator.__init__)

"; 

%feature("docstring") simuPOP::mapSelector::~mapSelector "

Usage:

    x.~mapSelector()

"; 

%feature("docstring") simuPOP::mapSelector::clone "

Description:

    deep copy of a map selector

Usage:

    x.clone()

"; 

%ignore simuPOP::mapSelector::indFitness(individual *ind, ULONG gen);

%feature("docstring") simuPOP::mapSelector::__repr__ "

Description:

    used by Python print function to print out the general information
    of the map selector

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::maQuanTrait "

Function form:

    MaQuanTrait

Description:

    multiple allele quantitative trait (quantitative trait according
    to disease or wildtype alleles)

Details:

    This is called 'multiple-allele' quantitative trait. It separates
    alleles into two groups: wildtype and diseased alleles. Wildtype
    alleles are specified by parameter wildtype and any other alleles
    are considered as diseased alleles. maQuanTrait accepts an array
    of fitness. Quantitative trait is then set for any given genotype.
    A standard normal distribution $ N\\left(0,\\sigma^{2}\\right) $ will
    be added to the returned trait value.

"; 

%feature("docstring") simuPOP::maQuanTrait::maQuanTrait "

Description:

    create a multiple allele quantitative trait operator

Usage:

    maQuanTrait(loci, qtrait, wildtype, sigma=[], ancestralGen=-1,
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[\"qtrait\"])

Details:

    Please refer to quanTrait for other parameter descriptions.

Arguments:

    qtrait:         an array of quantitative traits of AA, Aa, aa. A
                    is the wildtype group
    sigma:          an array of standard deviations for each of the
                    trait genotype (AA, Aa, aa)
    wildtype:       an array of alleles in the wildtype group. Any
                    other alleles will be considered as diseased
                    alleles. Default to [0].
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::maQuanTrait::~maQuanTrait "

Description:

    destructor

Usage:

    x.~maQuanTrait()

"; 

%feature("docstring") simuPOP::maQuanTrait::clone "

Description:

    deep copy of a multiple allele quantitative trait

Usage:

    x.clone()

"; 

%ignore simuPOP::maQuanTrait::qtrait(individual *ind);

%feature("docstring") simuPOP::maQuanTrait::__repr__ "

Description:

    used by Python print function to print out the general information
    of the multiple allele quantitative trait operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::maSelector "

Function form:

    MaSelect

Description:

    multiple allele selector (selection according to wildtype or
    diseased alleles)

Details:

    This is called 'multiple-allele' selector. It separates alleles
    into two groups: wildtype and diseased alleles. Wildtype alleles
    are specified by parameter wildtype and any other alleles are
    considered as diseased alleles. This selector accepts an array of
    fitness values:
    *  For single-locus, fitness is the fitness for genotypes AA, Aa,
    aa, while A stands for wildtype alleles.
    *  For a two-locus model, fitness is the fitness for genotypes
    AABB, AABb, AAbb, AaBB, AbBb, Aabb, aaBB, aaBb and aaBb.
    *  For a model with more than two loci, use a table of length $
    3^{n} $ in a order similar to the two-locus model.

"; 

%feature("docstring") simuPOP::maSelector::maSelector "

Description:

    create a multiple allele selector

Usage:

    maSelector(loci, fitness, wildtype, subPops=[], stage=PreMating,
      begin=0, end=-1, step=1, at=[], rep=[], subPop=[],
      infoFields=[\"fitness\"])

Details:

    Please refer to baseOperator for other parameter descriptions.

Arguments:

    fitness:        for the single locus case, fitness is an array of
                    fitness of AA, Aa, aa. A is the wildtype group. In
                    the case of multiple loci, fitness should be in
                    the order of AABB, AABb, AAbb, AaBB, AaBb, Aabb,
                    aaBB, aaBb, aabb.
    wildtype:       an array of alleles in the wildtype group. Any
                    other alleles are considered to be diseased
                    alleles. Default to [0].
    output:         and other parameters please refer to help
                    (baseOperator.__init__)

Note:

    * maSelector only works for diploid populations.
    * wildtype alleles at all loci are the same.

"; 

%feature("docstring") simuPOP::maSelector::~maSelector "

Usage:

    x.~maSelector()

"; 

%feature("docstring") simuPOP::maSelector::clone "

Description:

    deep copy of a maSelector

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::maSelector::indFitness "

Description:

    calculate/return the fitness value, currently assuming diploid

Usage:

    x.indFitness(ind, gen)

"; 

%feature("docstring") simuPOP::maSelector::__repr__ "

Description:

    used by Python print function to print out the general information
    of the maSelector

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mating "

Description:

    the base class of all mating schemes - a required parameter of
    simulator

Details:

    Mating schemes specify how to generate offspring from the current
    population. It must be provided when a simulator is created.
    Mating can perform the following tasks:
    *  change population/subpopulation sizes;
    *  randomly select parent(s) to generate offspring to populate the
    offspring generation;
    *  apply during-mating operators;
    *  apply selection if applicable.

"; 

%ignore simuPOP::mating::isCompatible(const population &pop) const;

%feature("docstring") simuPOP::mating::mating "

Description:

    create a mating scheme (do not use this base mating scheme, use
    one of its derived classes instead)

Usage:

    mating(newSubPopSize=[], newSubPopSizeExpr=\"\",
      newSubPopSizeFunc=None, subPop=[], weight=0)

Details:

    By default, a mating scheme keeps a constant population size,
    generates one offspring per mating event. These can be changed
    using certain parameters. newSubPopSize, newSubPopSizeExpr and
    newSubPopSizeFunc can be used to specify subpopulation sizes of
    the offspring generation.

Arguments:

    newSubPopSize:  an array of subpopulations sizes, should have the
                    same number of subpopulations as the current
                    population
    newSubPopSizeExpr:an expression that will be evaluated as an array
                    of new subpopulation sizes
    newSubPopSizeFunc:a function that takes parameters gen (generation
                    number) and oldsize (an array of current
                    population size) and return an array of
                    subpopulation sizes of the next generation. This
                    is usually easier to use than its expression
                    version of this parameter.
    subPop:         if this parameter is given, the mating scheme will
                    be applied only to the given (virtual)
                    subpopulation. This is only used in heteroMating
                    where mating schemes are passed to.
    weight:         When subPop is virtual, this is used to detemine
                    the number of offspring for this mating scheme.
                    Weight can be
                    * 0 (default) the weight will be proportional to
                    the current (virtual) subpopulation size. If other
                    virutal subpopulation has non-zero weight, this
                    virtual subpopulation will produce no offspring
                    (weight 0).
                    * any negative number -n: the size will be n*m
                    where m is the size of the (virtual) subpopulation
                    of the parental generation.
                    * any positive number n: the size will be
                    determined by weights from all (virtual)
                    subpopulations.

"; 

%ignore simuPOP::mating::mating(const mating &rhs);

%feature("docstring") simuPOP::mating::~mating "

Description:

    destructor

Usage:

    x.~mating()

"; 

%ignore simuPOP::mating::subPop() const;

%ignore simuPOP::mating::virtualSubPop() const;

%ignore simuPOP::mating::weight() const;

%feature("docstring") simuPOP::mating::clone "

Description:

    deep copy of a mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::mating::preparePopulation(population &pop);

%ignore simuPOP::mating::submitScratch(population &pop, population &scratch);

%ignore simuPOP::mating::mateSubPop(population &pop, SubPopID subPop, RawIndIterator offBegin, RawIndIterator offEnd, vector< baseOperator * > &ops);

%ignore simuPOP::mating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%ignore simuPOP::mating::prepareScratchPop(population &pop, population &scratch);

%feature("docstring") simuPOP::mendelianGenoTransmitter "

Applicability: diploid only

Details:

    Mendelian offspring generator accepts two parents and pass their
    genotype to a number of offspring following Mendelian's law.
    Basically, one of the paternal chromosomes is chosen randomly to
    form the paternal copy of the offspring, and one of the maternal
    chromosome is chosen randomly to form the maternal copy of the
    offspring. The number of offspring produced is controled by
    parameters numOffspring, numOffspringFunc, maxNumOffspring and
    mode. Recombination will not happen unless a during-mating
    operator recombinator is used.

"; 

%feature("docstring") simuPOP::mendelianGenoTransmitter::mendelianGenoTransmitter "

Usage:

    mendelianGenoTransmitter(begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::mendelianGenoTransmitter::clone "

Usage:

    x.clone()

Details:

    Return a cloned copy of an operator. This function is available to
    all operators.

"; 

%feature("docstring") simuPOP::mendelianGenoTransmitter::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::mendelianGenoTransmitter::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::mendelianGenoTransmitter::initialize "

Usage:

    x.initialize(pop)

"; 

%feature("docstring") simuPOP::mendelianGenoTransmitter::transmitGenotype "

Usage:

    x.transmitGenotype(parent, offspring, ploidy)

Details:

    Transmit genotype from parent to offspring, and fill the ploidy
    homologous set of chromosomes. This function does not set
    genotypes of customized chromosomes and handles sex chromosomes
    properly, according to offspring sex and ploidy.

"; 

%feature("docstring") simuPOP::mergeSubPops "

Function form:

    MergeSubPops

Description:

    merge subpopulations

Details:

    This operator merges subpopulations subPops to a single
    subpopulation. If subPops is ignored, all subpopulations will be
    merged.

"; 

%feature("docstring") simuPOP::mergeSubPops::mergeSubPops "

Description:

    merge subpopulations

Usage:

    mergeSubPops(subPops=[], stage=PreMating, begin=0, end=-1,
      step=1, at=[], rep=[], subPop=[], infoFields=[])

Arguments:

    subPops:        subpopulations to be merged. Default to all.

"; 

%feature("docstring") simuPOP::mergeSubPops::~mergeSubPops "

Description:

    destructor

Usage:

    x.~mergeSubPops()

"; 

%feature("docstring") simuPOP::mergeSubPops::clone "

Description:

    deep copy of a mergeSubPops operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mergeSubPops::apply "

Description:

    apply a mergeSubPops operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::mergeSubPops::__repr__ "

Description:

    used by Python print function to print out the general information
    of the mergeSubPops operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::migrator "

Description:

    migrate individuals from (virtual) subpopulations to other
    subpopulations

Details:

    Migrator is the only way to mix genotypes of several
    subpopulations because mating is strictly within subpopulations in
    simuPOP. Migrators are quite flexible in simuPOP in the sense that
    *  migration can happen from and to a subset of subpopulations.
    *  migration can be done by probability, proportion or by counts.
    In the case of probability, if the migration rate from
    subpopulation a to b is r, then everyone in subpopulation a will
    have this probability to migrate to b. In the case of proportion,
    exactly r*size_of_subPop_a individuals (chosen by random) will
    migrate to subpopulation b. In the last case, a given number of
    individuals will migrate.
    *  new subpopulation can be generated through migration. You
    simply need to migrate to a subpopulation with a new subpopulation
    number.

"; 

%feature("docstring") simuPOP::migrator::migrator "

Description:

    create a migrator

Usage:

    migrator(rate, mode=MigrByProbability, fromSubPop=[],
      toSubPop=[], stage=PreMating, begin=0, end=-1, step=1, at=[],
      rep=[], subPop=[], infoFields=[\"migrate_to\"])

Arguments:

    rate:           migration rate, can be a proportion or counted
                    number. Determined by parameter mode. rate should
                    be an m by n matrix. If a number is given, the
                    migration rate will be a m by n matrix of value r
    mode:           one of MigrByProbability (default),
                    MigrByProportion or MigrByCounts
    fromSubPop:     an array of 'from' subpopulations (a number) or
                    virtual subpopulations (a pair of numbers).
                    Default to all subpopulations. For example, if you
                    define a virtual subpopulation by sex, you can use
                    fromSubpop=[(0,0), 1] to choose migrants from the
                    first virtual subpopulation of subpopulation 0,
                    and from subpopulation 1. If a single number sp is
                    given, it is intepretted as [sp]. Note that
                    fromSubPop=(0, 1) (two subpopulation) is different
                    from fromSubPop=[(0,1)] (a virtual subpopulation).
    toSubPop:       an array of 'to' subpopulations. Default to all
                    subpopulations. If a single subpopulation is
                    specified, [] can be ignored.
    stage:          default to PreMating

Note:

    * The overall population size will not be changed. (Mating schemes
    can do that). If you would like to keep the subpopulation sizes
    after migration, you can use the newSubPopSize or
    newSubPopSizeExpr parameter of a mating scheme.
    * rate is a matrix with dimensions determined by fromSubPop and
    toSubPop. By default, rate is a matrix with element r(i,j), where
    r(i, j) is the migration rate, probability or count from
    subpopulation i to j. If fromSubPop and/or toSubPop are given,
    migration will only happen between these subpopulations. An
    extreme case is 'point migration', rate=[[r]], fromSubPop=a,
    toSubPop=b which migrate from subpopulation a to b with given rate
    r.

"; 

%feature("docstring") simuPOP::migrator::~migrator "

Description:

    destructor

Usage:

    x.~migrator()

"; 

%feature("docstring") simuPOP::migrator::clone "

Description:

    deep copy of a migrator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::migrator::rate "

Description:

    return migration rate

Usage:

    x.rate()

"; 

%feature("docstring") simuPOP::migrator::setRates "

Description:

    set migration rate

Usage:

    x.setRates(rate, mode)

Details:

    Format should be 0-0 0-1 0-2, 1-0 1-1 1-2, 2-0, 2-1, 2-2. For mode
    MigrByProbability or MigrByProportion, 0-0,1-1,2-2 will be set
    automatically regardless of input.

"; 

%feature("docstring") simuPOP::migrator::apply "

Description:

    apply the migrator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::migrator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the migrator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mitochondrialGenoTransmitter "

Details:

    This geno transmitter transmits some customized chromosomes as
    human mitochondrial chromosomes. It randomly inherit the first
    homologous copy of several customized chromosomes of the female
    parent.

"; 

%feature("docstring") simuPOP::mitochondrialGenoTransmitter::mitochondrialGenoTransmitter "

Usage:

    mitochondrialGenoTransmitter(chroms=[], begin=0, end=-1, step=1,
      at=[], rep=[], subPop=[], infoFields=[])

Details:

    chroms: if not given, all customized chromosomes.

"; 

%feature("docstring") simuPOP::mitochondrialGenoTransmitter::clone "

Usage:

    x.clone()

Details:

    Return a cloned copy of an operator. This function is available to
    all operators.

"; 

%feature("docstring") simuPOP::mitochondrialGenoTransmitter::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mitochondrialGenoTransmitter::initialize "

Usage:

    x.initialize(pop)

"; 

%ignore simuPOP::mitochondrialGenoTransmitter::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::mlPenetrance "

Function form:

    MlPenetrance

Description:

    penetrance according to the genotype according to a multiple loci
    multiplicative model

Details:

    This is the 'multiple-locus' penetrnace calculator. It accepts a
    list of penetrances and combine them according to the mode
    parameter, which takes one of the following values:
    *  PEN_Multiplicative: the penetrance is calculated as $ f=\\prod
    f_{i} $.
    *  PEN_Additive: the penetrance is calculated as $
    f=\\min\\left(1,\\sum f_{i}\\right) $. $ f $ will be set to 1 when $
    f<0 $. In this case, $ s_{i} $ are added, not $ f_{i} $ directly.
    *  PEN_Heterogeneity: the penetrance is calculated as $
    f=1-\\prod\\left(1-f_{i}\\right) $. Please refer to Neil Risch (1990)
    for detailed information about these models.

"; 

%feature("docstring") simuPOP::mlPenetrance::mlPenetrance "

Description:

    create a multiple locus penetrance operator

Usage:

    mlPenetrance(peneOps, mode=PEN_Multiplicative, ancestralGen=-1,
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Arguments:

    peneOps:        a list of penetrance operators
    mode:           can be one of PEN_Multiplicative, PEN_Additive,
                    and PEN_Heterogeneity

"; 

%feature("docstring") simuPOP::mlPenetrance::~mlPenetrance "

Usage:

    x.~mlPenetrance()

"; 

%feature("docstring") simuPOP::mlPenetrance::clone "

Description:

    deep copy of a multi-loci penetrance operator

Usage:

    x.clone()

"; 

%ignore simuPOP::mlPenetrance::penet(individual *ind);

%feature("docstring") simuPOP::mlPenetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the multiple-loci penetrance operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mlQuanTrait "

Function form:

    MlQuanTrait

Description:

    quantitative trait according to genotypes from a multiple loci
    multiplicative model

Details:

    Operator mlQuanTrait is a 'multiple-locus' quantitative trait
    calculator. It accepts a list of quantitative traits and combine
    them according to the mode parameter, which takes one of the
    following values
    *  QT_Multiplicative: the mean of the quantitative trait is
    calculated as $ f=\\prod f_{i} $.
    *  QT_Additive: the mean of the quantitative trait is calculated
    as $ f=\\sum f_{i} $. Note that all $ \\sigma_{i} $ (for $ f_{i} $)
    and $ \\sigma $ (for $ f $) will be considered. I.e, the trait
    value should be
    $ f=\\sum_{i}\\left(f_{i}+N\\left(0,\\sigma_{i}^{2}\\right)\\right)+\\sig
    ma^{2} $ for QT_Additive case. If this is not desired, you can set
    some of the $ \\sigma $ to zero.

"; 

%feature("docstring") simuPOP::mlQuanTrait::mlQuanTrait "

Description:

    create a multiple locus quantitative trait operator

Usage:

    mlQuanTrait(qtraits, mode=QT_Multiplicative, sigma=0,
      ancestralGen=-1, stage=PostMating, begin=0, end=-1, step=1,
      at=[], rep=[], subPop=[], infoFields=[\"qtrait\"])

Details:

    Please refer to quanTrait for other parameter descriptions.

Arguments:

    qtraits:        a list of quantitative traits
    mode:           can be one of QT_Multiplicative and QT_Additive

"; 

%feature("docstring") simuPOP::mlQuanTrait::~mlQuanTrait "

Usage:

    x.~mlQuanTrait()

"; 

%feature("docstring") simuPOP::mlQuanTrait::clone "

Description:

    deep copy of a multiple loci quantitative trait operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mlQuanTrait::qtrait "

Description:

    currently assuming diploid

Usage:

    x.qtrait(ind)

"; 

%feature("docstring") simuPOP::mlQuanTrait::__repr__ "

Description:

    used by Python print function to print out the general information
    of the multiple loci quantitative trait operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mlSelector "

Function form:

    MlSelect

Description:

    selection according to genotypes at multiple loci in a
    multiplicative model

Details:

    This selector is a 'multiple-locus model' selector. The selector
    takes a vector of selectors (can not be another mlSelector) and
    evaluate the fitness of an individual as the product or sum of
    individual fitness values. The mode is determined by parameter
    mode, which takes one of the following values
    *  SEL_Multiplicative: the fitness is calculated as $
    f=\\prod_{i}f_{i} $, where $ f_{i} $ is the single-locus fitness
    value.
    *  SEL_Additive: the fitness is calculated as $
    f=\\max\\left(0,1-\\sum_{i}(1-f_{i})\\right) $. $ f $ will be set to 0
    when $ f<0 $.

"; 

%feature("docstring") simuPOP::mlSelector::mlSelector "

Description:

    create a multiple-locus selector

Usage:

    mlSelector(selectors, mode=SEL_Multiplicative, subPops=[],
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[\"fitness\"])

Details:

    Please refer to mapSelector for other parameter descriptions.

Arguments:

    selectors:      a list of selectors

"; 

%feature("docstring") simuPOP::mlSelector::~mlSelector "

Usage:

    x.~mlSelector()

"; 

%feature("docstring") simuPOP::mlSelector::clone "

Description:

    deep copy of a mlSelector

Usage:

    x.clone()

"; 

%ignore simuPOP::mlSelector::indFitness(individual *ind, ULONG gen);

%feature("docstring") simuPOP::mlSelector::__repr__ "

Description:

    used by Python print function to print out the general information
    of the mlSelector

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mutator "

Description:

    Base class of all mutators.

Details:

    The base class of all functional mutators. It is not supposed to
    be called directly.
    Every mutator can specify rate (equal rate or different rates for
    different loci) and a vector of applicable loci (default to all
    but should have the same length as rate if rate has length greater
    than one).
    Maximum allele can be specified as well but more parameters, if
    needed, should be implemented by individual mutator classes.
    There are numbers of possible allelic states. Most theoretical
    studies assume an infinite number of allelic states to avoid any
    homoplasy. If it facilitates any analysis, this is however
    extremely unrealistic.

"; 

%feature("docstring") simuPOP::mutator::mutator "

Description:

    create a mutator, do not call this constructor directly

Usage:

    mutator(rate=[], loci=[], maxAllele=0, output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=[], subPop=[], infoFields=[])

Details:

    All mutators have the following common parameters. However, the
    actual meaning of these parameters may vary according to different
    models. The only differences between the following mutators are
    the way they actually mutate an allele, and corresponding input
    parameters. The number of mutation events at each locus is
    recorded and can be accessed from the mutationCount or
    mutationCounts functions.

Arguments:

    rate:           can be a number (uniform rate) or an array of
                    mutation rates (the same length as loci)
    loci:           a vector of locus indexes. Will be ignored only
                    when single rate is specified. Default to all
                    loci.
    maxAllele:      maximum allowed allele. Interpreted by each sub
                    mutator class. Default to pop.maxAllele().

"; 

%feature("docstring") simuPOP::mutator::~mutator "

Description:

    destructor

Usage:

    x.~mutator()

"; 

%feature("docstring") simuPOP::mutator::clone "

Description:

    deep copy of a mutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mutator::rate "

Description:

    return the mutation rate

Usage:

    x.rate()

"; 

%feature("docstring") simuPOP::mutator::setRate "

Description:

    set an array of mutation rates

Usage:

    x.setRate(rate, loci=[])

"; 

%feature("docstring") simuPOP::mutator::maxAllele "

Description:

    return maximum allowable allele number

Usage:

    x.maxAllele()

"; 

%feature("docstring") simuPOP::mutator::setMaxAllele "

Description:

    set maximum allowable allele

Usage:

    x.setMaxAllele(maxAllele)

"; 

%feature("docstring") simuPOP::mutator::mutationCount "

Description:

    return mutation count at locus

Usage:

    x.mutationCount(locus)

"; 

%feature("docstring") simuPOP::mutator::mutationCounts "

Description:

    return mutation counts

Usage:

    x.mutationCounts()

"; 

%feature("docstring") simuPOP::mutator::mutate "

Description:

    describe how to mutate a single allele

Usage:

    x.mutate(allele)

"; 

%feature("docstring") simuPOP::mutator::apply "

Description:

    apply a mutator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::noneOp "

Description:

    none operator

Details:

    This operator does nothing.

"; 

%feature("docstring") simuPOP::noneOp::noneOp "

Description:

    create a none operator

Usage:

    noneOp(output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
      end=0, step=1, at=[], rep=[], subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::noneOp::~noneOp "

Description:

    destructor

Usage:

    x.~noneOp()

"; 

%feature("docstring") simuPOP::noneOp::clone "Obsolete or undocumented function."

%ignore simuPOP::noneOp::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::noneOp::apply "

Description:

    apply the noneOp operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::noneOp::__repr__ "

Description:

    used by Python print function to print out the general information
    of the noneOp operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::offspringGenerator "

Applicability: all ploidy

Details:

    Offspring generators generate offspring from given parents.
    Generators differ from each other by how and how many offspring is
    generated at each mating event. Parameters mode, numOffspring,
    numOffspringParam and numOffspringFunc are used to specify how
    many offspring will be produced at each mating event. mode can be
    one of
    *  MATE_NumOffspring: a fixed number of offspring will be produced
    at all mating events .
    *  MATE_PyNumOffspring: A python function, specified by parameter
    numOffspringFunc, is called at each mating event to determine the
    number of offspring to produce.
    *  MATE_GeometricDistribution: a Geometric distribution with
    parameter numOffspring is used to determine the number of
    offspring of each family.
    *  MATE_PoissonDistribution: a Poisson distribution with parameter
    numOffspring is used to determine the number of offspring of each
    family.
    *  MATE_BinomialDistribution: a Binomial distribution with
    parameter numOffspring is used to determine the number of
    offspring of each family.
    *  MATE_UniformDistribution: a Uniform [a, b] distribution with
    parameter numOffspring (a) and numOffspringParam (b) is used to
    determine the number of offspring of each family. This is the base
    class of all offspring generators, and should not be used
    directly.

"; 

%feature("docstring") simuPOP::offspringGenerator::offspringGenerator "

Usage:

    offspringGenerator(ops, numParents=0, numOffspring=1.,
      numOffspringFunc=None, numOffspringParam=1,
      mode=MATE_NumOffspring, sexParam=0.5, sexMode=MATE_RandomSex)

Arguments:

    numOffspring:   Depending on , this paramter can be the number of
                    offspring to produce, or a paremter of a random
                    distribution.
    numOffspringFunc:a Python function that returns the number of
                    offspring at each mating event. The setting of
                    this parameter implies =MATE_PyNumOffspring.
    numOffspringParam:used when numOffspring is generated from a
                    binomial or random distribution.
    mode:           can be one of MATE_NumOffspring,
                    MATE_PyNumOffspring, MATE_GeometricDistribution,
                    MATE_PoissonDistribution,
                    MATE_BinomialDistribution, or
                    MATE_UniformDistribution.
    sexParam:       parameter that controls the sex distribution among
                    offspring. Its exact meaning is determined by
                    parameter sexMode. Default to 0.5.
    sexMode:        can be one of
                    * MATE_RandomSex Set sex to Male or Female with
                    probability 0.5. Parameter sexParam is ignored in
                    this case. This is the default mode.
                    * MATE_ProbOfMale Set an offspring to Male with
                    probability sexParam (default to 0.5)
                    * MATE_NumOfMale Set sexParam offspring to Male
                    * MATE_NumOfFemale Set sexParam offspring to
                    Female. If there are sex chromosomes, sex is
                    determined by sex chromosomes when sexMode id
                    MATE_RandomSex. Otherwise, some offspring will be
                    rejected so that offspring sex match what is
                    specified in other modes.
    transmitter:    is an during mating operator, that will be used if
                    no during mating operator is used to produce
                    offspring.

"; 

%feature("docstring") simuPOP::offspringGenerator::~offspringGenerator "

Usage:

    x.~offspringGenerator()

"; 

%ignore simuPOP::offspringGenerator::offspringGenerator(const offspringGenerator &rhs);

%feature("docstring") simuPOP::offspringGenerator::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::offspringGenerator::initialize(const population &pop, vector< baseOperator * > const &ops);

%ignore simuPOP::offspringGenerator::generateOffspring(population &pop, individual *dad, individual *mom, RawIndIterator &offBegin, RawIndIterator &offEnd, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::offspringGenerator::finalize "

Usage:

    x.finalize(pop)

"; 

%feature("docstring") simuPOP::offspringGenerator::initialized "

Usage:

    x.initialized()

"; 

%ignore simuPOP::offspringGenerator::numOffspring(int gen);

%ignore simuPOP::offspringGenerator::getSex(int count);

%ignore simuPOP::offspringGenerator::numParents() const;

%ignore simuPOP::OstreamManager;

%ignore simuPOP::OstreamManager::OstreamManager();

%feature("docstring") simuPOP::OstreamManager::~OstreamManager "

Usage:

    x.~OstreamManager()

"; 

%ignore simuPOP::OstreamManager::getOstream(const string &name, bool readable, bool realAppend, bool useString);

%ignore simuPOP::OstreamManager::hasOstream(const string &filename);

%ignore simuPOP::OstreamManager::listAll();

%ignore simuPOP::OstreamManager::closeAll();

%feature("docstring") simuPOP::OutOfMemory "

Description:

    exception, thrown if out of memory

"; 

%feature("docstring") simuPOP::OutOfMemory::OutOfMemory "

Usage:

    OutOfMemory(msg)

"; 

%feature("docstring") simuPOP::outputer "

Description:

    Base class of all operators that out information. different
    format.

Details:

    Bo Peng

"; 

%feature("docstring") simuPOP::outputer::outputer "

Description:

    constructor.

Usage:

    outputer(output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
      end=-1, step=1, at=[], rep=[], subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::outputer::~outputer "

Description:

    destructor

Usage:

    x.~outputer()

"; 

%feature("docstring") simuPOP::outputer::clone "

Usage:

    x.clone()

Details:

    Return a cloned copy of an operator. This function is available to
    all operators.

"; 

%feature("docstring") simuPOP::parentChooser "

Details:

    Parent choosers repeatedly choose parent(s) from a parental
    population, and pass them to offspring generators. A parent
    chooser can select one or two parents, which should match what is
    required by the offspring generator used. This is the base class
    of all parent choosers, and should not be used directly.

"; 

%feature("docstring") simuPOP::parentChooser::parentChooser "

Usage:

    parentChooser(numParents)

"; 

%feature("docstring") simuPOP::parentChooser::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::parentChooser::initialize(population &pop, SubPopID subPop);

%feature("docstring") simuPOP::parentChooser::finalize "

Usage:

    x.finalize(pop, subPop)

"; 

%ignore simuPOP::parentChooser::initialized() const;

%ignore simuPOP::parentChooser::numParents();

%ignore simuPOP::parentChooser::chooseParent(RawIndIterator basePtr);

%ignore simuPOP::parentChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::parentChooser::~parentChooser "

Usage:

    x.~parentChooser()

"; 

%feature("docstring") simuPOP::parentsTagger "

Description:

    tagging according to parents' indexes

Details:

    This during-mating operator set tag(), currently a pair of
    numbers, of each individual with indexes of his/her parents in the
    parental population. This information will be used by pedigree-
    related operators like affectedSibpairSample to track the pedigree
    information. Because parental population will be discarded or
    stored after mating, these index will not be affected by post-
    mating operators. This tagger record parental index to one or both
    *  one or two information fields. Default to father_idx and
    mother_idx. If only one parent is passed in a mating scheme (such
    as selfing), only the first information field is used. If two
    parents are passed, the first information field records paternal
    index, and the second records maternal index.
    *  a file. Indexes will be written to this file. This tagger will
    also act as a post-mating operator to add a new-line to this file.

"; 

%feature("docstring") simuPOP::parentsTagger::parentsTagger "

Description:

    create a parentsTagger

Usage:

    parentsTagger(begin=0, end=-1, step=1, at=[], rep=[], subPop=[],
      output=\"\", outputExpr=\"\", infoFields=[\"father_idx\",
      \"mother_idx\"])

"; 

%feature("docstring") simuPOP::parentsTagger::~parentsTagger "

Usage:

    x.~parentsTagger()

"; 

%feature("docstring") simuPOP::parentsTagger::clone "

Description:

    deep copy of a parentsTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::parentsTagger::__repr__ "

Description:

    used by Python print function to print out the general information
    of the parentsTagger

Usage:

    x.__repr__()

"; 

%ignore simuPOP::parentsTagger::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::parentsTagger::apply "

Usage:

    x.apply(pop)

Details:

    at the end of a generation, write population structure information
    to a file with a newline.

"; 

%feature("docstring") simuPOP::parentTagger "

Description:

    tagging according to parental indexes

Details:

    This during-mating operator set tag() each individual with indexes
    of his/her parent in the parental population. Because only one
    parent is recorded, this is recommended to be used for mating
    schemes that requires only one parent (such as selfMating). This
    tagger record indexes to information field parent_idx, and/or a
    given file. The usage is similar to parentsTagger.

"; 

%feature("docstring") simuPOP::parentTagger::parentTagger "

Description:

    create a parentTagger

Usage:

    parentTagger(begin=0, end=-1, step=1, at=[], rep=[], subPop=[],
      output=\"\", outputExpr=\"\", infoFields=[\"parent_idx\"])

"; 

%feature("docstring") simuPOP::parentTagger::~parentTagger "

Usage:

    x.~parentTagger()

"; 

%feature("docstring") simuPOP::parentTagger::clone "

Description:

    deep copy of a parentTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::parentTagger::__repr__ "

Description:

    used by Python print function to print out the general information
    of the parentTagger

Usage:

    x.__repr__()

"; 

%ignore simuPOP::parentTagger::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::parentTagger::apply "

Usage:

    x.apply(pop)

Details:

    at the end of a generation, write population structure information
    to a file with a newline.

"; 

%feature("docstring") simuPOP::pause "

Description:

    pause a simulator

Details:

    This operator pauses the evolution of a simulator at given
    generations or at a key stroke, using stopOnKeyStroke=True option.
    Users can use 'q' to stop an evolution. When a simulator is
    stopped, press any other key to resume the simulation or escape to
    a Python shell to examine the status of the simulation by pressing
    's'.
    There are two ways to use this operator, the first one is to pause
    the simulation at specified generations, using the usual operator
    parameters such as at. Another way is to pause a simulation with
    any key stroke, using the stopOnKeyStroke parameter. This feature
    is useful for a presentation or an interactive simulation. When
    's' is pressed, this operator expose the current population to the
    main Python dictionary as variable pop and enter an interactive
    Python session. The way current population is exposed can be
    controlled by parameter exposePop and popName. This feature is
    useful when you want to examine the properties of a population
    during evolution.

"; 

%feature("docstring") simuPOP::pause::pause "

Description:

    stop a simulation. Press 'q' to exit or any other key to continue.

Usage:

    pause(prompt=True, stopOnKeyStroke=False, exposePop=True,
      popName=\"pop\", output=\">\", outputExpr=\"\", stage=PostMating,
      begin=0, end=-1, step=1, at=[], rep=-1, subPop=[],
      infoFields=[])

Arguments:

    prompt:         if True (default), print prompt message.
    stopOnKeyStroke:if True, stop only when a key was pressed.
    exposePop:      whether or not expose pop to user namespace, only
                    useful when user choose 's' at pause. Default to
                    True.
    popName:        by which name the population is exposed. Default
                    to pop.

"; 

%feature("docstring") simuPOP::pause::~pause "

Description:

    destructor

Usage:

    x.~pause()

"; 

%feature("docstring") simuPOP::pause::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::pause::apply "

Description:

    apply the pause operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pause::__repr__ "

Description:

    used by Python print function to print out the general information
    of the pause operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pedigree "

"; 

%feature("docstring") simuPOP::pedigree::pedigree "

Usage:

    pedigree(pop, loci=[], infoFields=[], ancGen=-1)

Details:

    Create a pedigree object from a population, using a subset of
    loci, information fields and ancestral generations.

"; 

%feature("docstring") simuPOP::pedigree::locateRelatives "

Usage:

    x.locateRelatives(relType, relFields, gen=-1, relSex=AnySex,
      parentFields=[])

Details:

    This function locates relatives (of type relType, and sex relSex)
    of each individual and store their indexes in specified
    information fields relFields. The indexes of parents in the
    parental generation should be available in information fields
    parentFields (default to ['father_idx', 'mother_idx'] which are
    the information fields used by operator parentsTagger. This
    function currently only work for diploid populations.

Arguments:

    relType:        Relative type, which can be
                    * REL_Self set indexes of individual themselves.
                    * REL_Spouse locate spouses of individuals in the
                    current generation. A spouse is defined as two
                    individuals having an offspring with shared
                    parentFields. If more than one infoFields is
                    given, multiple spouses can be identified.
                    * REL_Offspring index of offspring in the
                    offspring generation. If only one parent is given,
                    only paternal or maternal relationship is
                    considered. For example,
                    parentFields=['father_idx'] will locate offspring
                    for all fathers.
                    * REL_FullSibling all siblings with the same
                    parents
                    * REL_Sibling all sibs with at least one shared
                    parent
    relFields:      information fields to hold relatives. The number
                    of these fields limits the number of relatives to
                    locate.
    gen:            Find relatives for individuals for how many
                    generations. Default to -1, meaning for all
                    generations. If a non-negative number is given, up
                    till generation gen will be processed.
    relSex:         Whether or not only locate relative or certain
                    sex. It can be AnySex (do not care, default),
                    MaleOnly, FemaleOnly, or OppositeSex (only locate
                    relatives of opposite sex.

"; 

%feature("docstring") simuPOP::pedigree::setIndexesOfRelatives "

Description:

    Trace a relative path in a population and record the result in the
    given information fields.

Usage:

    x.setIndexesOfRelatives(pathGen, pathFields, pathSex=[],
      resultFields=[])

Details:

    For example, setInfoWithRelatives(pathGen = [0, 1, 1, 0],
    pathFields = [['father_idx', 'mother_idx'], ['sib1', 'sib2'],
    ['off1', 'off2']], pathSex = [AnySex, MaleOnly, FemaleOnly],
    resultFields = ['cousin1', 'cousin2']) This function will 1.
    locate father_idx and mother_idx for each individual at generation
    0 (pathGen[0]) 2. find AnySex individuals referred by father_idx
    and mother_idx at generation 1 (pathGen[1]) 3. find informaton
    fields sib1 and sib2 from these parents 4. locate MaleOnly
    individuals referred by sib1 and sib2 from generation 1
    (pathGen[2]) 5. find information fields off1 and off2 from these
    individuals, and 6. locate FemaleOnly indiviudals referred by off1
    and from geneartion 0 (pathGen[3]) 7. Save index of these
    individuals to information fields cousin1 and cousin2 at
    genearation pathGen[0]. In short, this function locates father or
    mother's brother's daughters.

Arguments:

    pathGen:        A list of generations that form a relative path.
                    This array is one element longer than pathFields,
                    with gen_i, gen_i+1 indicating the current and
                    destinating generation of information fields
                    path_i.
    pathFields:     A list of list of information fields forming a
                    path to trace a certain type of relative.
    resultFields:   Where to store located relatives. Note that the
                    result will be saved in the starting generation
                    specified in pathGen[0], which is usually 0.
    pathSex:        (Optional) A list of sex choices, AnySex, Male,
                    Female or OppositeSex, that is used to choose
                    individuals at each step. Default to AnySex.

"; 

%feature("docstring") simuPOP::pedigreeTagger "

Details:

    Pedigree tagger is used to save a complete pedigree to a pedigree
    file during an evolution process. Because is destroyedof record
    individuals involved in an evolutioary process. This is a simple
    post-mating tagger that write given information fields to a file
    (or standard output).

"; 

%feature("docstring") simuPOP::pedigreeTagger::pedigreeTagger "

Usage:

    pedigreeTagger(begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], stage=PostMating, output=\">\", outputExpr=\"\",
      pedigreeFields=[])

"; 

%ignore simuPOP::pedigreeTagger::apply(population &pop);

%feature("docstring") simuPOP::penetrance "

Description:

    Base class of all penetrance operators.

Details:

    Penetrance is the probability that one will have the disease when
    he has certain genotype(s). An individual will be randomly marked
    as affected/unaffected according to his/her penetrance value. For
    example, an individual will have probability 0.8 to be affected if
    the penetrance is 0.8.
    Penetrance can be applied at any stage (default to DuringMating).
    When a penetrance operator is applied, it calculates the
    penetrance value of each offspring and assigns affected status
    accordingly. Penetrance can also be used PreMating or PostMating.
    In these cases, the affected status will be set to all individuals
    according to their penetrance values.
    Penetrance values are usually not saved. If you would like to know
    the penetrance value, you need to
    *  use addInfoField('penetrance') to the population to analyze.
    (Or use infoFields parameter of the population constructor), and
    *  use e.g., mlPenetrance(...., infoFields=['penetrance']) to add
    the penetrance field to the penetrance operator you use. You may
    choose a name other than 'penetrance' as long as the field names
    for the operator and population match. Penetrance functions can be
    applied to the current, all, or certain number of ancestral
    generations. This is controlled by the ancestralGen parameter,
    which is default to -1 (all available ancestral generations). You
    can set it to 0 if you only need affection status for the current
    generation, or specify a number n for the number of ancestral
    generations (n + 1 total generations) to process. Note that the
    ancestralGen parameter is ignored if the penetrance operator is
    used as a during mating operator.

"; 

%feature("docstring") simuPOP::penetrance::penetrance "

Description:

    create a penetrance operator

Usage:

    penetrance(ancestralGen=-1, stage=DuringMating, begin=0, end=-1,
      step=1, at=[], rep=[], subPop=[], infoFields=[])

Arguments:

    ancestralGen:   if this parameter is set to be 0, apply penetrance
                    to the current generation; if -1, apply to all
                    generations; otherwise, apply to the specified
                    numbers of ancestral generations.
    stage:          specify the stage this operator will be applied.
                    Default to DuringMating.
    infoFields:     If one field is specified, it will be used to
                    store penetrance values.

"; 

%feature("docstring") simuPOP::penetrance::~penetrance "

Description:

    destructor

Usage:

    x.~penetrance()

"; 

%feature("docstring") simuPOP::penetrance::clone "

Description:

    deep copy of a penetrance operator

Usage:

    x.clone()

"; 

%ignore simuPOP::penetrance::penet(individual *);

%feature("docstring") simuPOP::penetrance::apply "

Description:

    set penetrance to all individuals and record penetrance if
    requested

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::penetrance::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::penetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the penetrance operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pointMutator "

Function form:

    PointMutate

Description:

    point mutator

Details:

    Mutate specified individuals at specified loci to a spcified
    allele. I.e., this is a non-random mutator used to introduce
    diseases etc. pointMutator, as its name suggest, does point
    mutation. This mutator will turn alleles at loci on the first
    chromosome copy to toAllele for individual inds. You can specify
    atPloidy to mutate other, or all ploidy copies.

"; 

%feature("docstring") simuPOP::pointMutator::pointMutator "

Description:

    create a pointMutator

Usage:

    pointMutator(loci, toAllele, atPloidy=[], inds=[], output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=[], subPop=[], infoFields=[])

Details:

    Please see class mutator for the descriptions of other parameters.

Arguments:

    inds:           individuals who will mutate
    toAllele:       allele that will be mutate to

"; 

%feature("docstring") simuPOP::pointMutator::~pointMutator "

Description:

    destructor

Usage:

    x.~pointMutator()

"; 

%feature("docstring") simuPOP::pointMutator::clone "

Description:

    deep copy of a pointMutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pointMutator::apply "

Description:

    apply a pointMutator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pointMutator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the pointMutator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pointMutator::mutationCount "

Description:

    return mutation count at locus

Usage:

    x.mutationCount(locus)

"; 

%feature("docstring") simuPOP::pointMutator::mutationCounts "

Description:

    return mutation counts

Usage:

    x.mutationCounts()

"; 

%feature("docstring") simuPOP::polyParentsChooser "

Applicability: all ploidy

Details:

    This parent chooser chooses two parents randomly, a male and a
    female, from their respective sex groups randomly. If selection is
    turned on, parents are chosen from their sex groups with
    probabilities that are proportional to their fitness values. Note
    that selection is not allowed in the case of monopoly because this
    poses a particular order on individuals in the offspring
    generation. This parents chooser also allows polygamous mating by
    reusing a parent multiple times when returning parents, and allows
    specification of a few alpha individuals who will be the only
    mating individuals in their sex group.

"; 

%feature("docstring") simuPOP::polyParentsChooser::polyParentsChooser "

Usage:

    polyParentsChooser(polySex=Male, polyNum=1)

Details:

    Note: If selection is enabled, it works regularly on on-alpha sex,
    but works twice on alpha sex. That is to say, alphaNum alpha
    indiviudals are chosen selectively, and selected again during
    mating.

Arguments:

    polySex:        Male (polygyny) or Female (polyandry) parent that
                    will have polyNum sex partners.
    polyNum:        Number of sex partners.

"; 

%feature("docstring") simuPOP::polyParentsChooser::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::polyParentsChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::polyParentsChooser::chooseParents(RawIndIterator basePtr);

%ignore simuPOP::polyParentsChooser::numMale();

%ignore simuPOP::polyParentsChooser::numFemale();

%feature("docstring") simuPOP::population "

Details:

    A simuPOP population consists of individuals of the same genotypic
    structure, organized by generations, subpopulations and virtual
    subpopulations. It also contains a Python dictionary that is used
    to store arbitrary population variables. In addition to genotypic
    structured related functions provided by the genoStruTrait class,
    the population class provides a large number of member functions
    that can be used to
    *  Create, copy and compare populations.
    *  Manipulate subpopulations. A population can be divided into
    several subpopulations. Because individuals only mate with
    individuals within the same subpopulation, exchange of genetic
    information across subpopulations can only be done through
    migration. A number of functions are provided to access
    subpopulation structure information, and to merge and split
    subpopulations.
    *  Define and access virtual subpopulations. A virtual
    subpopulation splitter can be assigned to a population, which
    defines groups of individuals called virtual subpopulations (VSP)
    within each subpopulation.
    *  Access individuals individually, or through iterators that
    iterate through individuals in (virtual) subpopulations.
    *  Access genotype and information fields of individuals at the
    population level. From a population point of view, all genotypes
    are arranged sequentially individual by individual. Please refer
    to class individual for an introduction to genotype arragement of
    each individual.
    *  Store and access ancestral generations. A population can save
    arbitrary number of ancestral generations. It is possible to
    directly access an ancestor, or make an ancestral generation the
    current generation for more efficient access.
    *  Insert or remove loci, resize (shrink or expand) a population,
    sample from a population, or merge with other populations.
    *  Manipulate population variables and evaluate expressions in
    this local namespace.
    *  Save and load a population.

"; 

%feature("docstring") simuPOP::population::population "

Usage:

    population(size=[], ploidy=2, loci=[], chromTypes=[],
      lociPos=[], ancGen=0, chromNames=[], alleleNames=[],
      lociNames=[], subPopNames=[], infoFields=[])

Details:

    The following parameters are used to create a population object:

Arguments:

    size:           A list of subpopulation sizes. The length of this
                    list determines the number of subpopulations of
                    this population. If there is no subpopulation,
                    size=[popSize] can be written as size=popSize.
    ploidy:         Number of homologous sets of chromosomes. Default
                    to 2 (diploid). For efficiency considerations, all
                    chromosomes have the same number of homologous
                    sets, even if some customized chromosomes or some
                    individuals (e.g. males in a haplodiploid
                    population) have different numbers of homologous
                    sets. The first case is handled by setting
                    chromTypes of each chromosome. Only the
                    haplodiploid populations are handled for the
                    second case, for which ploidy=Haplodiploid should
                    be used.
    loci:           A list of numbers of loci on each chromosome. The
                    length of this parameter determines the number of
                    chromosomes. Default to [1], meaning one
                    chromosome with a single locus.
    chromTypes:     A list that specifies the type of each chromosome,
                    which can be Autosome, ChromosomeX, ChromosomeY,
                    or Customized. All chromosomes are assumed to be
                    autosomes if this parameter is ignored. Sex
                    chromosome can only be specified in a diploid
                    population where the sex of an individual is
                    determined by the existence of these chromosomes
                    using the XX (Female) and XY (Male) convention.
                    Both sex chromosomes have to be available and be
                    specified only once. Because chromosomes X and Y
                    are treated as two chromosomes, recombination on
                    the pseudo-autosomal regions of the sex chromsomes
                    is not supported. Customized chromosomes are
                    special chromosomes whose inheritance patterns are
                    undefined. They rely on user-defined functions and
                    operators to be passed from parents to offspring.
                    Multiple customized chromosomes have to be
                    arranged consecutively.
    lociPos:        Positions of all loci on all chromosome, as a list
                    of float numbers. Default to 1, 2, ... etc on each
                    chromosome. Positions on the same chromosome
                    should be ordered. A nested list that specifies
                    positions of loci on each chromosome is also
                    acceptable.
    ancGen:         Number of the most recent ancestral generations to
                    keep during evolution. Default to 0, which means
                    only the current generation will be kept. If it is
                    set to -1, all ancestral generations will be kept
                    in this population (and exhaust your computer RAM
                    quickly).
    chromNames:     A list of chromosome names. Default to chrom1,
                    chrom2, ... etc.
    alleleNames:    A list of allele names for all markers. For
                    example, alleleNames=('A','C','T','G') names
                    allele 0 -- 3'A', 'C', 'T', and 'G' respectively.
                    Note that simuPOP does not yet support locus-
                    specific allele names.
    lociNames:      A list or a matrix (separated by chromosomes) of
                    names for each locus. Default to \"locX-Y\" where X
                    and Y are 1-based chromosome and locus indexes,
                    respectively.
    subPopNames:    A list of subpopulation names. All subpopulations
                    will have name '' if this parameter is not
                    specified.
    infoFields:     Names of information fields (named float number)
                    that will be attached to each individual.

"; 

%ignore simuPOP::population::population(const population &rhs);

%feature("docstring") simuPOP::population::clone "

Usage:

    x.clone()

Details:

    Create a cloned copy of a population. Note that Python statement
    pop1 = pop only creates a reference to an existing population pop.

"; 

%feature("docstring") simuPOP::population::swap "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::~population "

Description:

    destroy a population

Usage:

    x.~population()

"; 

%feature("docstring") simuPOP::population::__repr__ "

Description:

    used by Python print function to print out the general information
    of the population

Usage:

    x.__repr__()

"; 

%ignore simuPOP::population::validate(const string &msg) const;

%ignore simuPOP::population::fitSubPopStru(const vectorlu &newSubPopSizes, const vectorstr &newSubPopNames);

%ignore simuPOP::population::hasActivatedVirtualSubPop() const;

%ignore simuPOP::population::hasActivatedVirtualSubPop(SubPopID subPop) const;

%ignore simuPOP::population::hasVirtualSubPop() const;

%ignore simuPOP::population::virtualSplitter() const;

%feature("docstring") simuPOP::population::setVirtualSplitter "

Usage:

    x.setVirtualSplitter(splitter)

Details:

    Set a VSP splitter to the population, which defines the same VSPs
    for all subpopulations. If different VSPs are needed for different
    subpopulations, a combinedSplitter can be used to make these VSPs
    available to all subpopulations.

"; 

%feature("docstring") simuPOP::population::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of virtual subpopulations (VSP) defined by a VSP
    splitter. Return 0 if no VSP is defined.

"; 

%feature("docstring") simuPOP::population::activateVirtualSubPop "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::deactivateVirtualSubPop "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::__cmp__ "

Description:

    a python function used to compare the population objects

Usage:

    x.__cmp__(rhs)

"; 

%feature("docstring") simuPOP::population::setSubPopStru "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::numSubPop "

Usage:

    x.numSubPop()

Details:

    Return the number of subpopulations in a population. Return 1 if
    there is no subpopulation structure.

"; 

%feature("docstring") simuPOP::population::subPopSize "

Usage:

    x.subPopSize(subPop)

Details:

    Return the size of a subpopulation (subPopSize(sp)) or a virtual
    subpopulation (subPopSize([sp, vsp])).

"; 

%feature("docstring") simuPOP::population::subPopByName "

Usage:

    x.subPopByName(name)

Details:

    Return the index of the first subpopulation with name name. An
    IndexError will be raised if subpopulations are not named, or if
    no subpopulation with name name is found. Virtual subpopulation
    name is not supported.

"; 

%feature("docstring") simuPOP::population::subPopName "

Usage:

    x.subPopName(subPop)

Details:

    Return the name of a subpopulation subPop, and 'unnamed' if no
    name is assigned to subPop. If subPop is a virtual subpopulation
    (specified by a (sp, vsp) pair), a combined name such as subPop1 -
    Male is returned.

"; 

%feature("docstring") simuPOP::population::subPopNames "

Usage:

    x.subPopNames()

Details:

    Return the names of all subpopulations (excluding virtual
    subpopulations). 'unnamed' will be returned for unnamed
    subpopulations.

"; 

%feature("docstring") simuPOP::population::setSubPopName "

Usage:

    x.setSubPopName(name, subPop)

Details:

    Assign a name name to subpopulation subPop. does not have to be
    unique.

"; 

%feature("docstring") simuPOP::population::subPopSizes "

Usage:

    x.subPopSizes()

Details:

    Return the sizes of all subpopulations in a list. Virtual
    subpopulations are not considered.

"; 

%feature("docstring") simuPOP::population::popSize "

Usage:

    x.popSize()

Details:

    Return the total number of individuals in all subpopulations.

"; 

%feature("docstring") simuPOP::population::absIndIndex "

Usage:

    x.absIndIndex(idx, subPop)

Details:

    return the absolute index of an individual idx in subpopulation
    subPop.

"; 

%feature("docstring") simuPOP::population::subPopIndPair "

Usage:

    x.subPopIndPair(idx)

Details:

    Return the subpopulation ID and relative index of an individual,
    given its absolute index idx.

"; 

%feature("docstring") simuPOP::population::subPopBegin "

Usage:

    x.subPopBegin(subPop)

Details:

    Return the index of the first individual in subpopulation subPop.
    An IndexError will be raised if subPop is out of range.

"; 

%feature("docstring") simuPOP::population::subPopEnd "

Usage:

    x.subPopEnd(subPop)

Details:

    Return the index of the last individual in subpopulation subPop
    plus 1, so that range(subPopBegin(subPop), subPopEnd(subPop) can
    iterate through the index of all individuals in subpopulation
    subPop.

"; 

%feature("docstring") simuPOP::population::ind "

Usage:

    x.ind(idx)

Details:

    Return a refernce to individual ind in the population.

"; 

%feature("docstring") simuPOP::population::ind "

Usage:

    x.ind(idx, subPop)

Details:

    Return a refernce to individual ind in subpopulation subPop.

"; 

%ignore simuPOP::population::ind(ULONG idx, UINT subPop=0) const;

%feature("docstring") simuPOP::population::ancestor "

Usage:

    x.ancestor(idx, gen)

Details:

    Return a reference to individual idx in ancestral generation gen.
    The correct individual will be returned even if the current
    generation is not the present one (see also useAncestralGen).

"; 

%feature("docstring") simuPOP::population::ancestor "

Description:

    refrence to an individual ind in an ancestral generation

Usage:

    x.ancestor(ind, gen)

"; 

%feature("docstring") simuPOP::population::ancestor "

Usage:

    x.ancestor(ind, subPop, gen)

Details:

    Return a reference to individual idx of subpopulation subPop in
    ancestral generation gen.

"; 

%feature("docstring") simuPOP::population::ancestor "

Description:

    refrence to an individual ind in a specified subpopulaton or an
    ancestral generation

Usage:

    x.ancestor(ind, subPop, gen)

"; 

%feature("docstring") simuPOP::population::individuals "

Usage:

    x.individuals()

Details:

    Return a Python iterator that can be used to iterate through all
    individuals in a population.

"; 

%feature("docstring") simuPOP::population::individuals "

Usage:

    x.individuals(subPop)

Details:

    Return an iterator that can be used to iterate through all
    individuals in a subpopulation (subPop=spID) or a virtual
    subpopulation (subPop=[spID, vspID]).

"; 

%ignore simuPOP::population::indOrdered() const;

%ignore simuPOP::population::setIndOrdered(bool s) const;

%ignore simuPOP::population::indBegin(IterationType type=VisibleInds);

%ignore simuPOP::population::indEnd(IterationType type=VisibleInds);

%ignore simuPOP::population::indBegin(UINT subPop, IterationType type=VisibleInds);

%ignore simuPOP::population::indEnd(UINT subPop, IterationType type=VisibleInds);

%ignore simuPOP::population::indBegin(IterationType type=VisibleInds) const;

%ignore simuPOP::population::indEnd(IterationType type=VisibleInds) const;

%ignore simuPOP::population::indEnd(UINT subPop, IterationType type) const;

%ignore simuPOP::population::rawIndBegin();

%ignore simuPOP::population::rawIndEnd();

%ignore simuPOP::population::rawIndBegin(UINT subPop);

%ignore simuPOP::population::rawIndEnd(UINT subPop);

%ignore simuPOP::population::rawIndBegin() const;

%ignore simuPOP::population::rawIndEnd() const;

%ignore simuPOP::population::rawIndEnd(UINT subPop) const;

%ignore simuPOP::population::alleleBegin(UINT locus);

%ignore simuPOP::population::alleleEnd(UINT locus);

%ignore simuPOP::population::alleleBegin(UINT locus, UINT subPop);

%ignore simuPOP::population::genoBegin(bool order);

%ignore simuPOP::population::genoEnd(bool order);

%ignore simuPOP::population::genoBegin(UINT subPop, bool order);

%ignore simuPOP::population::genoEnd(UINT subPop, bool order);

%ignore simuPOP::population::indGenoBegin(ULONG ind) const;

%ignore simuPOP::population::indGenoEnd(ULONG ind) const;

%feature("docstring") simuPOP::population::arrGenotype "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::arrGenotype "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::genotype "

Usage:

    x.genotype()

Details:

    Return an editable array of the genotype of all individuals in
    this population.

"; 

%feature("docstring") simuPOP::population::genotype "

Usage:

    x.genotype(subPop)

Details:

    Return an editable array of the genotype of all individuals in
    subpopulation subPop. Virtual subpopulation is unsupported.

"; 

%feature("docstring") simuPOP::population::setGenotype "

Usage:

    x.setGenotype(geno)

Details:

    Fill the genotype of all individuals of a population using a list
    of alleles geno. geno will be reused if its length is less than
    popSize()*totNumLoci()*ploidy().

"; 

%feature("docstring") simuPOP::population::setGenotype "

Usage:

    x.setGenotype(geno, subPop)

Details:

    Fill the genotype of all individuals of in (virtual) subpopulation
    subPop using a list of alleles geno. geno will be reused if its
    length is less than subPopSize(subPop)*totNumLoci()*ploidy().

"; 

%feature("docstring") simuPOP::population::setSubPopByIndInfo "

Usage:

    x.setSubPopByIndInfo(field)

Details:

    Rearrange individuals to their new subpopulations according to
    their integer values at information field field (value returned by
    individual::indInfo(field)). Individuals with negative values at
    this field will be removed. Existing subpopulation names are
    unchanged but new subpopulations will not assign a name
    ('unnamed').

"; 

%feature("docstring") simuPOP::population::splitSubPop "

Usage:

    x.splitSubPop(subPop, sizes)

Details:

    Split subpopulation subPop into subpopulations of given sizes,
    which should add up to the size of subpopulation subPop.
    Alternatively, sizes can be a list of proportions (add up to 1)
    from which the sizes of new subpopulations are determined. If
    subPop is not the last subpopulation, indexes of subpopulations
    after subPop are shifted. If subPop is named, the same name will
    be given to all split subpopulations.

"; 

%feature("docstring") simuPOP::population::removeSubPops "

Usage:

    x.removeSubPops(subPops)

Details:

    Remove subpopulations subPop and all their individuals. Indexes of
    subpopulations after removed subpopulations will be shifted.

"; 

%feature("docstring") simuPOP::population::removeIndividuals "

Usage:

    x.removeIndividuals(inds)

Details:

    remove individuals inds (absolute indexes) from the current
    population. A subpopulation will be kept even if all individuals
    from it are removed. This function only affects the current
    generation.

"; 

%feature("docstring") simuPOP::population::mergeSubPops "

Usage:

    x.mergeSubPops(subPops=[])

Details:

    Merge subpopulations subPops. If subPops is empty (default), all
    subpopulations will be merged. subPops do not have to be adjacent
    to each other. They will all be merged to the subpopulation with
    the smallest subpopulation ID. Indexes of the rest of the
    subpopulation may be changed.

"; 

%feature("docstring") simuPOP::population::addIndFromPop "

Usage:

    x.addIndFromPop(pop)

Details:

    Add all individuals, including ancestors, in pop to the current
    population. Two populations should have the same genotypic
    structures and number of ancestral generations. Subpopulations in
    population pop are kept.

"; 

%feature("docstring") simuPOP::population::addChromFromPop "

Usage:

    x.addChromFromPop(pop)

Details:

    Add chromosomes in population pop to the current population.
    Population pop should have the same number of individuals as the
    current population in the current and all ancestral generations.
    This function merges genotypes on the new chromosomes from
    population pop individual by individual.

"; 

%feature("docstring") simuPOP::population::addLociFromPop "

Usage:

    x.addLociFromPop(pop)

Details:

    Add loci from population pop, chromosome by chromosome. Added loci
    will be inserted according to their position. Their position and
    names should not overlap with any locus in the current population.
    Population pop should have the same number of individuals as the
    current population in the current and all ancestral generations.

"; 

%feature("docstring") simuPOP::population::addChrom "

Usage:

    x.addChrom(lociPos, lociNames=[], chromName=\"\",
      chromType=Autosome)

Details:

    Add chromosome chromName with given type chromType to a
    population, with loci lociNames inserted at position lociPos.
    lociPos should be ordered. lociNames and chromName should not
    exist in the current population. If they are not specified,
    simuPOP will try to assign default names, and raise a ValueError
    if the default names have been used.

"; 

%feature("docstring") simuPOP::population::addLoci "

Usage:

    x.addLoci(chrom, pos, names=[])

Details:

    Insert loci names at positions pos on chromosome chrom. These
    parameters should be lists of the same length, although names may
    be ignored, in which case random names will be given. Alleles at
    inserted loci are initialized with zero alleles. Note that loci
    have to be added to existing chromosomes. If loci on a new
    chromosome need to be added, function addChrom should be used.
    This function returns indexes of the inserted loci.

"; 

%feature("docstring") simuPOP::population::resize "

Usage:

    x.resize(newSubPopSizes, propagate=False)

Details:

    Resize population by giving new subpopulation sizes
    newSubPopSizes. Individuals at the end of some subpopulations will
    be removed if the new subpopulation size is smaller than the old
    one. New individuals will be appended to a subpopulation if the
    new size is larger. Their genotypes will be set to zero (default),
    or be copied from existing individuals if propagate is set to
    True. More specifically, if a subpopulation with 3 individuals is
    expanded to 7, the added individuals will copy genotypes from
    individual 1, 2, 3, and 1 respectively. Note that this function
    only resizes the current generation.

"; 

%feature("docstring") simuPOP::population::extract "

Usage:

    x.extract(field=None, loci=None, infoFields=None, ancGen=-1)

Details:

    Extract subsets of individuals, loci and/or information fields
    from the current population and create a new one. If information
    field field is not None, individuals with negative values at this
    information field will be removed, and others are put into
    subpopulations specified by this field. If loci is not None, only
    genotypes at loci are extracted. If infoFields is not None, only
    these information fields will be extracted. If ancGen is not -1
    (default, meaing all ancestral generations), only ancGen ancestral
    generations will be kept. As an advanced feature, field can be
    information field of a pedigree object ped. This allows extraction
    of individuals according to pedigrees identified in a pedigree
    object. This pedigree should have the same number of individuals
    in all generations.

"; 

%feature("docstring") simuPOP::population::removeLoci "

Usage:

    x.removeLoci(loci=[], keep=[])

Details:

    Remove loci (absolute indexes) and genotypes at these loci from
    the current population. Alternatively, a parameter keep can be
    used to specify loci that will not be removed.

"; 

%feature("docstring") simuPOP::population::push "

Usage:

    x.push(pop)

Details:

    Push population pop into the current population. Both populations
    should have the same genotypic structure. The current population
    is discarded if ancestralDepth (maximum number of ancestral
    generations to hold) is zero so no ancestral generation can be
    kept. Otherise, the current population will become the parental
    generation of pop, advancing the greatness level of all existing
    ancestral generations by one. If ancestralDepth is positive and
    there are already ancestralDepth ancestral generations (see also:
    ancestralGens()), the greatest ancestral generation will be
    discarded. In any case, population pop becomes invalid as all its
    individuals are absorbed by the current population.

"; 

%feature("docstring") simuPOP::population::ancestralGens "

Usage:

    x.ancestralGens()

Details:

    Return the actual number of ancestral generations stored in a
    population, which does not necessarily equal to the number set by
    setAncestralDepth().

"; 

%feature("docstring") simuPOP::population::setIndInfo "

Usage:

    x.setIndInfo(values, idx)

Details:

    Set information field idx (an index) of the current population to
    values. values will be reused if its length is smaller than
    popSize().

"; 

%feature("docstring") simuPOP::population::setIndInfo "

Usage:

    x.setIndInfo(values, name)

Details:

    Set information field name of the current population to values.
    values will be reused if its length is smaller than popSize().

"; 

%feature("docstring") simuPOP::population::setIndInfo "

Usage:

    x.setIndInfo(values, idx, subPop)

Details:

    Set information field idx (an index) of a subpopulation
    (subPop=sp) or a virtual subpopulation (subPop=[sp, vsp]) to
    values. values will be reused if its length is smaller than
    subPopSize(subPop).

"; 

%feature("docstring") simuPOP::population::setIndInfo "

Usage:

    x.setIndInfo(values, name, subPop)

Details:

    Set information field name of a subpopulation (subPop=sp) or a
    virtual subpopulation (subPop=[sp, vsp]) to values. values will be
    reused if its length is smaller than subPopSize(subPop).

"; 

%ignore simuPOP::population::infoBegin(UINT idx);

%ignore simuPOP::population::infoEnd(UINT idx);

%feature("docstring") simuPOP::population::indInfo "

Usage:

    x.indInfo(idx)

Details:

    Return the information field idx (an index) of all individuals as
    a list.

"; 

%feature("docstring") simuPOP::population::indInfo "

Usage:

    x.indInfo(name)

Details:

    Return the information field name of all individuals as a list.

"; 

%feature("docstring") simuPOP::population::indInfo "

Usage:

    x.indInfo(idx, subPop)

Details:

    Return the information field idx (an index) of all individuals in
    (virtual) subpopulation subPop as a list.

"; 

%feature("docstring") simuPOP::population::indInfo "

Usage:

    x.indInfo(name, subPop)

Details:

    Return the information field name of all individuals in (virtual)
    subpopulation subPop as a list.

"; 

%feature("docstring") simuPOP::population::addInfoField "

Usage:

    x.addInfoField(field, init=0)

Details:

    Add an information field field to a population and initialize its
    values to init.

"; 

%feature("docstring") simuPOP::population::addInfoFields "

Usage:

    x.addInfoFields(fields, init=0)

Details:

    Add a list of information fields fields to a population and
    initialize their values to init. If an information field alreay
    exists, it will be re-initialized.

"; 

%feature("docstring") simuPOP::population::setInfoFields "

Usage:

    x.setInfoFields(fields, init=0)

Details:

    Set information fields fields to a population and initialize them
    with value init. All existing information fields will be removed.

"; 

%feature("docstring") simuPOP::population::updateInfoFieldsFrom "

Usage:

    x.updateInfoFieldsFrom(fields, pop, fromFields=[], ancGen=-1)

Details:

    Update information fields fields from fromFields of another
    population (or pedigree) pop. Two populations should have the same
    number of individuals. If fromFields is not specified, it is
    assumed to be the same as fields. If ancGen is not -1, only the
    most recent ancGen generations are updated.

"; 

%feature("docstring") simuPOP::population::setAncestralDepth "

Usage:

    x.setAncestralDepth(depth)

Details:

    set the intended ancestral depth of a population to depth, which
    can be 0 (does not store any ancestral generation), -1 (store all
    ancestral generations), and a positive number (store depth
    ancestral generations.

"; 

%feature("docstring") simuPOP::population::useAncestralGen "

Usage:

    x.useAncestralGen(idx)

Details:

    Making ancestral generation idx (0 for current generation, 1 for
    parental generation, 2 for grand-parental generation, etc) the
    current generation. This is an efficient way to access population
    properties of an ancestral generation. useAncestralGen(0) should
    always be called afterward to restore the correct order of
    ancestral generations.

"; 

%ignore simuPOP::population::equalTo(const population &rhs);

%ignore simuPOP::population::sortIndividuals(bool infoOnly=false) const;

%feature("docstring") simuPOP::population::save "

Usage:

    x.save(filename)

Details:

    Save population to a file filename, which can be loaded by a
    global function LoadPopulation(filename).

"; 

%ignore simuPOP::population::load(const string &filename);

%ignore simuPOP::population::selectionOn() const;

%ignore simuPOP::population::selectionOn(UINT sp) const;

%feature("docstring") simuPOP::population::turnOffSelection "Obsolete or undocumented function."

%ignore simuPOP::population::turnOnSelection(UINT sp);

%ignore simuPOP::population::turnOnSelection();

%ignore simuPOP::population::rep();

%ignore simuPOP::population::setRep(int rep, bool setVar=true);

%ignore simuPOP::population::gen();

%ignore simuPOP::population::setGen(ULONG gen, bool setVar=true);

%feature("docstring") simuPOP::population::vars "

Usage:

    x.vars()

Details:

    return variables of a population as a Python dictionary.

"; 

%feature("docstring") simuPOP::population::vars "

Usage:

    x.vars(subPop)

Details:

    return a dictionary vars()[\"subPop\"][subPop]. subPop can be a
    number (subPop=spID), or a pair of numbers (subPop=(spID, vspID)).
    A ValueError will be raised if key 'subPop' does not exist in
    vars(), or if key subPop does not exist in vars()[\"subPop\"].

"; 

%ignore simuPOP::population::dict(int subPop=-1);

%ignore simuPOP::population::setDict(PyObject *dict);

%ignore simuPOP::population::hasVar(const string &name);

%ignore simuPOP::population::removeVar(const string &name);

%ignore simuPOP::population::setBoolVar(const string &name, const bool val);

%ignore simuPOP::population::setIntVar(const string &name, const int val);

%ignore simuPOP::population::setDoubleVar(const string &name, const double val);

%ignore simuPOP::population::setStringVar(const string &name, const string &val);

%ignore simuPOP::population::setIntVectorVar(const string &name, const vectori &val);

%ignore simuPOP::population::setDoubleVectorVar(const string &name, const vectorf &val);

%ignore simuPOP::population::setStrDictVar(const string &name, const strDict &val);

%ignore simuPOP::population::setIntDictVar(const string &name, const intDict &val);

%ignore simuPOP::population::setVar(const string &name, PyObject *val);

%ignore simuPOP::population::getVar(const string &name, bool nameError=true);

%ignore simuPOP::population::getVarAsBool(const string &name, bool nameError=true);

%ignore simuPOP::population::getVarAsInt(const string &name, bool nameError=true);

%ignore simuPOP::population::getVarAsDouble(const string &name, bool nameError=true);

%ignore simuPOP::population::getVarAsString(const string &name, bool nameError=true);

%ignore simuPOP::population::getVarAsStrDict(const string &name, bool nameError=true);

%ignore simuPOP::population::getVarAsIntDict(const string &name, bool nameError=true);

%ignore simuPOP::population::varsAsString() const;

%ignore simuPOP::population::varsFromString(const string &vars);

%feature("docstring") simuPOP::population::evaluate "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::execute "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::scramble "Obsolete or undocumented function."

%feature("docstring") simuPOP::proportionSplitter "

Details:

    This splitter divides subpopulations into several VSPs by
    proportion.

"; 

%feature("docstring") simuPOP::proportionSplitter::proportionSplitter "

Usage:

    proportionSplitter(proportions=[])

Details:

    Create a splitter that divides subpopulations by proportions,
    which should be a list of float numbers (between 0 and 1) that add
    up to 1.

"; 

%feature("docstring") simuPOP::proportionSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::proportionSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::proportionSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs defined by this splitter, which is the
    length of parameter proportions.

"; 

%ignore simuPOP::proportionSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::proportionSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::proportionSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of VSP vsp, which is \"Prop p\" where
    p=propotions[vsp].

"; 

%feature("docstring") simuPOP::pyEval "

Function form:

    PyEval

Description:

    evaluate an expression

Details:

    Python expressions/statements will be executed when pyEval is
    applied to a population by using parameters expr/stmts. Statements
    can also been executed when pyEval is created and destroyed or
    before expr is executed. The corresponding parameters are
    preStmts, postStmts and stmts. For example, operator varPlotter
    uses this feature to initialize R plots and save plots to a file
    when finished.

"; 

%feature("docstring") simuPOP::pyEval::pyEval "

Description:

    evaluate expressions/statments in the local namespace of a
    replicate

Usage:

    pyEval(expr=\"\", stmts=\"\", preStmts=\"\", postStmts=\"\",
      exposePop=False, name=\"\", output=\">\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Arguments:

    expr:           the expression to be evaluated. The result will be
                    sent to output.
    stmts:          the statement that will be executed before the
                    expression
    preStmts:       the statement that will be executed when the
                    operator is constructed
    postStmts:      the statement that will be executed when the
                    operator is destroyed
    exposePop:      if True, expose the current population as a
                    variable named pop
    name:           used to let pure Python operator to identify
                    themselves
    output:         default to >. I.e., output to standard output.

"; 

%feature("docstring") simuPOP::pyEval::~pyEval "

Usage:

    x.~pyEval()

"; 

%feature("docstring") simuPOP::pyEval::clone "

Description:

    deep copy of a pyEval operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyEval::apply "

Description:

    apply the pyEval operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyEval::__repr__ "

Description:

    used by Python print function to print out the general information
    of the pyEval operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyEval::name "

Description:

    return the name of an expression

Usage:

    x.name()

Details:

    The name of a pyEval operator is given by an optional parameter
    name. It can be used to identify this pyEval operator in debug
    output, or in the dryrun mode of simulator::evolve.

"; 

%feature("docstring") simuPOP::pyExec "

Function form:

    PyExec

Description:

    execute a Python statement

Details:

    This operator takes a list of statements and executes them. No
    value will be returned or outputted.

"; 

%feature("docstring") simuPOP::pyExec::pyExec "

Description:

    evaluate statments in the local replicate namespace, no return
    value

Usage:

    pyExec(stmts=\"\", preStmts=\"\", postStmts=\"\", exposePop=False,
      name=\"\", output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
      end=-1, step=1, at=[], rep=[], subPop=[], infoFields=[])

Details:

    Please refer to class pyEval for parameter descriptions.

"; 

%feature("docstring") simuPOP::pyExec::~pyExec "

Usage:

    x.~pyExec()

"; 

%feature("docstring") simuPOP::pyExec::clone "

Description:

    deep copy of a pyExec operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyExec::__repr__ "

Description:

    used by Python print function to print out the general information
    of the pyExec operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyIndIterator "

Details:

    this class implements a Python itertor class that can be used to
    iterate through individuals in a (sub)population. If allInds are
    true, visiblility of individuals will not be checked. Note that
    individualIterator *will* iterate through only visible
    individuals, and allInds is only provided when we know in advance
    that all individuals are visible. This is a way to obtain better
    performance in simple cases. An instance of this class is returned
    by population::individuals() and population::individuals(subPop)

"; 

%feature("docstring") simuPOP::pyIndIterator::pyIndIterator "

Usage:

    pyIndIterator(begin, end, allInds, allVisibles)

"; 

%feature("docstring") simuPOP::pyIndIterator::~pyIndIterator "

Usage:

    x.~pyIndIterator()

"; 

%feature("docstring") simuPOP::pyIndIterator::__iter__ "

Usage:

    x.__iter__()

"; 

%feature("docstring") simuPOP::pyIndIterator::next "

Usage:

    x.next()

"; 

%feature("docstring") simuPOP::pyIndOperator "

Description:

    individual operator

Details:

    This operator is similar to a pyOperator but works at the
    individual level. It expects a function that accepts an
    individual, optional genotype at certain loci, and an optional
    parameter. When it is applied, it passes each individual to this
    function. When infoFields is given, this function should return an
    array to fill these infoFields. Otherwise, True or False is
    expected. More specifically, func can be
    *  func(ind) when neither loci nor param is given.
    *  func(ind, genotype) when loci is given.
    *  func(ind, param) when param is given.
    *  func(ind, genotype, param) when both loci and param are given.

"; 

%feature("docstring") simuPOP::pyIndOperator::pyIndOperator "

Description:

    a Pre- or PostMating Python operator that apply a function to each
    individual

Usage:

    pyIndOperator(func, loci=[], param=None, stage=PostMating,
      formOffGenotype=False, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Arguments:

    func:           a Python function that accepts an individual and
                    optional genotype and parameters.
    param:          any Python object that will be passed to func
                    after pop parameter. Multiple parameters can be
                    passed as a tuple.
    infoFields:     if given, func is expected to return an array of
                    the same length and fill these infoFields of an
                    individual.

"; 

%feature("docstring") simuPOP::pyIndOperator::~pyIndOperator "

Description:

    destructor

Usage:

    x.~pyIndOperator()

"; 

%ignore simuPOP::pyIndOperator::pyIndOperator(const pyIndOperator &rhs);

%feature("docstring") simuPOP::pyIndOperator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::pyIndOperator::apply "

Description:

    apply the pyIndOperator operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyIndOperator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the pyIndOperator operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyMating "

Applicability: all ploidy

Description:

    a Python mating scheme

Details:

    This hybrid mating scheme does not have to involve a python
    function. It requires a parent chooser, and an offspring
    generator. The parent chooser chooses parent(s) and pass them to
    the offspring generator to produce offspring.

"; 

%feature("docstring") simuPOP::pyMating::pyMating "

Description:

    create a Python mating scheme

Usage:

    pyMating(chooser, generator, newSubPopSize=[],
      newSubPopSizeExpr=\"\", newSubPopSizeFunc=None, subPop=[],
      weight=0)

Arguments:

    chooser:        a parent chooser that chooses parent(s) from the
                    parental generation.
    generator:      an offspring generator that produce offspring of
                    given parents.

"; 

%feature("docstring") simuPOP::pyMating::~pyMating "

Description:

    destructor

Usage:

    x.~pyMating()

"; 

%ignore simuPOP::pyMating::pyMating(const pyMating &rhs);

%feature("docstring") simuPOP::pyMating::clone "

Description:

    deep copy of a Python mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the Python mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::pyMating::mateSubPop(population &pop, SubPopID subPop, RawIndIterator offBegin, RawIndIterator offEnd, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::pyMigrator "

Description:

    a more flexible Python migrator

Details:

    This migrator can be used in two ways
    *  define a function that accepts a generation number and returns
    a migration rate matrix. This can be used in various migration
    rate cases.
    *  define a function that accepts individuals etc, and returns the
    new subpopulation ID. More specifically, func can be
    *  func(ind) when neither loci nor param is given.
    *  func(ind, genotype) when loci is given.
    *  func(ind, param) when param is given.
    *  func(ind, genotype, param) when both loci and param are given.

"; 

%feature("docstring") simuPOP::pyMigrator::pyMigrator "

Description:

    create a hybrid migrator

Usage:

    pyMigrator(rateFunc=None, indFunc=None, mode=MigrByProbability,
      fromSubPop=[], toSubPop=[], loci=[], param=None,
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[\"migrate_to\"])

Arguments:

    rateFunc:       a Python function that accepts a generation
                    number, current subpopulation sizes, and returns a
                    migration rate matrix. The migrator then migrate
                    like a usual migrator.
    indFunc:        a Python function that accepts an individual,
                    optional genotypes and parameters, then returns a
                    subpopulation ID. This method can be used to
                    separate a population according to individual
                    genotype.
    stage:          default to PreMating

"; 

%feature("docstring") simuPOP::pyMigrator::~pyMigrator "

Description:

    destructor

Usage:

    x.~pyMigrator()

"; 

%ignore simuPOP::pyMigrator::pyMigrator(const pyMigrator &rhs);

%feature("docstring") simuPOP::pyMigrator::clone "

Description:

    deep copy of a pyMigrator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyMigrator::apply "

Description:

    apply a pyMigrator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyMigrator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the pyMigrator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyMutator "

Function form:

    PyMutate

Description:

    A hybrid mutator.

Details:

    Parameters such as mutation rate of this operator are set just
    like others and you are supposed to provide a Python function to
    return a new allele state given an old state. pyMutator will
    choose an allele as usual and call your function to mutate it to
    another allele.

"; 

%feature("docstring") simuPOP::pyMutator::pyMutator "

Description:

    create a pyMutator

Usage:

    pyMutator(rate=[], loci=[], maxAllele=0, func=None, output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=[], subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::pyMutator::~pyMutator "

Description:

    destructor

Usage:

    x.~pyMutator()

"; 

%ignore simuPOP::pyMutator::pyMutator(const pyMutator &rhs);

%feature("docstring") simuPOP::pyMutator::clone "

Description:

    deep copy of a pyMutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyMutator::mutate "

Description:

    mutate according to the mixed model

Usage:

    x.mutate(allele)

"; 

%feature("docstring") simuPOP::pyMutator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the pyMutator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyOperator "

Description:

    A python operator that directly operate a population.

Details:

    This operator accepts a function that can take the form of
    *  func(pop) when stage=PreMating or PostMating, without setting
    param;
    *  func(pop, param) when stage=PreMating or PostMating, with
    param;
    *  func(pop, off, dad, mom) when stage=DuringMating and
    passOffspringOnly=False, without setting param;
    *  func(off) when stage=DuringMating and passOffspringOnly=True,
    and without setting param;
    *  func(pop, off, dad, mom, param) when stage=DuringMating and
    passOffspringOnly=False, with param;
    *  func(off, param) when stage=DuringMating and
    passOffspringOnly=True, with param. For Pre- and PostMating
    usages, a population and an optional parameter is passed to the
    given function. For DuringMating usages, population, offspring,
    its parents and an optional parameter are passed to the given
    function. Arbitrary operations can be applied to the population
    and offspring (if stage=DuringMating).

"; 

%feature("docstring") simuPOP::pyOperator::pyOperator "

Description:

    Python operator, using a function that accepts a population
    object.

Usage:

    pyOperator(func, param=None, stage=PostMating,
      formOffGenotype=False, passOffspringOnly=False, begin=0, end=-1,
      step=1, at=[], rep=[], subPop=[], infoFields=[])

Arguments:

    func:           a Python function. Its form is determined by other
                    parameters.
    param:          any Python object that will be passed to func
                    after pop parameter. Multiple parameters can be
                    passed as a tuple.
    formOffGenotype:This option tells the mating scheme this operator
                    will set the genotype of offspring (valid only for
                    stage=DuringMating). By default
                    (formOffGenotype=False), a mating scheme will set
                    the genotype of offspring before it is passed to
                    the given Python function. Otherwise, a 'blank'
                    offspring will be passed.
    passOffspringOnly:if True, pyOperator will expect a function of form
                    func(off [,param]), instead of func(pop, off, dad,
                    mom [, param]) which is used when
                    passOffspringOnly is False. Because many during-
                    mating pyOperator only need access to offspring,
                    this will improve efficiency. Default to False.

Note:

    * Output to output or outputExpr is not supported. That is to say,
    you have to open/close/append to files explicitly in the Python
    function. Because files specified by output or outputExpr are
    controlled (opened/closed) by simulators, they should not be
    manipulated in a pyOperator operator.
    * This operator can be applied Pre-, During- or Post- Mating and
    is applied PostMating by default. For example, if you would like
    to examine the fitness values set by a selector, a PreMating
    Python operator should be used.

"; 

%feature("docstring") simuPOP::pyOperator::~pyOperator "

Description:

    destructor

Usage:

    x.~pyOperator()

"; 

%ignore simuPOP::pyOperator::pyOperator(const pyOperator &rhs);

%feature("docstring") simuPOP::pyOperator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::pyOperator::apply "

Description:

    apply the pyOperator operator to one population

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::pyOperator::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::pyOperator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the pyOperator operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyOutput "

Description:

    Output a given string.

Details:

    A common usage is to output a new line for the last replicate.

"; 

%feature("docstring") simuPOP::pyOutput::pyOutput "

Description:

    Create a pyOutput operator that outputs a given string.

Usage:

    pyOutput(str=\"\", output=\">\", outputExpr=\"\", stage=PostMating,
      begin=0, end=-1, step=1, at=[], rep=[], subPop=[],
      infoFields=[])

Arguments:

    str:            string to be outputted

"; 

%feature("docstring") simuPOP::pyOutput::apply "

Description:

    simply output some info

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyOutput::~pyOutput "

Usage:

    x.~pyOutput()

"; 

%feature("docstring") simuPOP::pyOutput::clone "

Usage:

    x.clone()

Details:

    Return a cloned copy of an operator. This function is available to
    all operators.

"; 

%feature("docstring") simuPOP::pyOutput::setString "

Description:

    set output string.

Usage:

    x.setString(str)

"; 

%feature("docstring") simuPOP::pyOutput::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyParentsChooser "

Applicability: all ploidy

Details:

    This parents chooser accept a Python generator function that
    yields repeatedly an index (relative to each subpopulation) of a
    parent, or indexes of two parents as a Python list of tuple. The
    generator function is responsible for handling sex or selection if
    needed.

"; 

%feature("docstring") simuPOP::pyParentsChooser::pyParentsChooser "

Usage:

    pyParentsChooser(parentsGenerator)

Arguments:

    parentsGenerator:A Python generator function

"; 

%ignore simuPOP::pyParentsChooser::pyParentsChooser(const pyParentsChooser &rhs);

%feature("docstring") simuPOP::pyParentsChooser::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::pyParentsChooser::initialize(population &pop, SubPopID sp);

%feature("docstring") simuPOP::pyParentsChooser::finalize "

Usage:

    x.finalize(pop, sp)

"; 

%feature("docstring") simuPOP::pyParentsChooser::~pyParentsChooser "

Description:

    destructor

Usage:

    x.~pyParentsChooser()

"; 

%ignore simuPOP::pyParentsChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::pyPenetrance "

Function form:

    PyPenetrance

Description:

    assign penetrance values by calling a user provided function

Details:

    For each individual, the penetrance is determined by a user-
    defined penetrance function func. This function takes genetypes at
    specified loci, and optionally values of specified information
    fields. The return value is considered as the penetrance for this
    individual. More specifically, func can be
    *  func(geno) if infoFields has length 0 or 1.
    *  func(geno, fields) when infoFields has more than 1 fields. Both
    parameters should be an list.

"; 

%feature("docstring") simuPOP::pyPenetrance::pyPenetrance "

Description:

    provide locus and penetrance for 11, 12, 13 (in the form of
    dictionary)

Usage:

    pyPenetrance(loci, func, ancestralGen=-1, stage=DuringMating,
      begin=0, end=-1, step=1, at=[], rep=[], subPop=[],
      infoFields=[])

Arguments:

    loci:           the genotypes at these loci will be passed to the
                    provided Python function in the form of loc1_1,
                    loc1_2, loc2_1, loc2_2, ... if the individuals are
                    diploid.
    func:           a user-defined Python function that accepts an
                    array of genotypes at specified loci and return a
                    penetrance value. The return value should be
                    between 0 and 1.
    infoFields:     if specified, the first field should be the
                    information field to save calculated penetrance
                    value. The values of the rest of the information
                    fields (if available) will also be passed to the
                    user defined penetrance function.
    output:         and other parameters please refer to help
                    (baseOperator.__init__)

"; 

%feature("docstring") simuPOP::pyPenetrance::~pyPenetrance "

Description:

    destructor

Usage:

    x.~pyPenetrance()

"; 

%ignore simuPOP::pyPenetrance::pyPenetrance(const pyPenetrance &rhs);

%feature("docstring") simuPOP::pyPenetrance::clone "

Description:

    deep copy of a Python penetrance operator

Usage:

    x.clone()

"; 

%ignore simuPOP::pyPenetrance::penet(individual *ind);

%feature("docstring") simuPOP::pyPenetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the Python penetrance operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::pyPopIterator;

%feature("docstring") simuPOP::pyPopIterator::pyPopIterator "

Usage:

    pyPopIterator(begin, end)

"; 

%feature("docstring") simuPOP::pyPopIterator::~pyPopIterator "

Usage:

    x.~pyPopIterator()

"; 

%feature("docstring") simuPOP::pyPopIterator::__iter__ "

Usage:

    x.__iter__()

"; 

%feature("docstring") simuPOP::pyPopIterator::next "

Usage:

    x.next()

"; 

%feature("docstring") simuPOP::pyQuanTrait "

Function form:

    PyQuanTrait

Description:

    quantitative trait using a user provided function

Details:

    For each individual, a user provided function is used to calculate
    quantitative trait.

"; 

%feature("docstring") simuPOP::pyQuanTrait::pyQuanTrait "

Description:

    create a Python quantitative trait operator

Usage:

    pyQuanTrait(loci, func, ancestralGen=-1, stage=PostMating,
      begin=0, end=-1, step=1, at=[], rep=[], subPop=[],
      infoFields=[\"qtrait\"])

Details:

    Please refer to quanTrait for other parameter descriptions.

Arguments:

    loci:           The genotypes at these loci will be passed to
                    func.
    func:           a Python function that accepts genotypes at
                    specified loci and returns the quantitative trait
                    value.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::pyQuanTrait::~pyQuanTrait "

Usage:

    x.~pyQuanTrait()

"; 

%ignore simuPOP::pyQuanTrait::pyQuanTrait(const pyQuanTrait &rhs);

%feature("docstring") simuPOP::pyQuanTrait::clone "

Description:

    deep copy of a Python quantitative trait operator

Usage:

    x.clone()

"; 

%ignore simuPOP::pyQuanTrait::qtrait(individual *ind);

%feature("docstring") simuPOP::pyQuanTrait::__repr__ "

Description:

    used by Python print function to print out the general information
    of the Python quantitative trait operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pySelector "

Function form:

    PySelect

Description:

    selection using user provided function

Details:

    This selector assigns fitness values by calling a user provided
    function. It accepts a list of loci and a Python function func.
    For each individual, this operator will pass the genotypes at
    these loci, generation number, and optionally values at some
    information fields to this function. The return value is treated
    as the fitness value. The genotypes are arranged in the order of
    0-0,0-1,1-0,1-1 etc. where X-Y represents locus X - ploidy Y. More
    specifically, func can be
    *  func(geno, gen) if infoFields has length 0 or 1.
    *  func(geno, gen, fields) when infoFields has more than 1 fields.
    Values of fields 1, 2, ... will be passed. Both geno and fields
    should be a list.

"; 

%feature("docstring") simuPOP::pySelector::pySelector "

Description:

    create a Python hybrid selector

Usage:

    pySelector(loci, func, subPops=[], stage=PreMating, begin=0,
      end=-1, step=1, at=[], rep=[], subPop=[],
      infoFields=[\"fitness\"])

Arguments:

    loci:           susceptibility loci. The genotype at these loci
                    will be passed to func.
    func:           a Python function that accepts genotypes at
                    specified loci, generation number, and optionally
                    information fields. It returns the fitness value.
    output:         and other parameters please refer to help
                    (baseOperator.__init__)
    infoFields:     if specified, the first field should be the
                    information field to save calculated fitness value
                    (should be 'fitness' in most cases). The values of
                    the rest of the information fields (if available)
                    will also be passed to the user defined penetrance
                    function.

"; 

%feature("docstring") simuPOP::pySelector::~pySelector "

Description:

    destructor

Usage:

    x.~pySelector()

"; 

%ignore simuPOP::pySelector::pySelector(const pySelector &rhs);

%feature("docstring") simuPOP::pySelector::clone "

Description:

    deep copy of a pySelector

Usage:

    x.clone()

"; 

%ignore simuPOP::pySelector::indFitness(individual *ind, ULONG gen);

%feature("docstring") simuPOP::pySelector::__repr__ "

Description:

    used by Python print function to print out the general information
    of the pySelector

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyTagger "

Description:

    Python tagger.

Details:

    This tagger takes some information fields from both parents, pass
    to a Python function and set the individual field with the return
    value. This operator can be used to trace the inheritance of trait
    values.

"; 

%feature("docstring") simuPOP::pyTagger::pyTagger "

Description:

    creates a pyTagger that works on specified information fields

Usage:

    pyTagger(func=None, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], output=\"\", outputExpr=\"\", infoFields=[])

Arguments:

    infoFields:     information fields. The user should gurantee the
                    existence of these fields.
    func:           a Pyton function that returns a list to assign the
                    information fields. e.g., if fields=['A', 'B'],
                    the function will pass values of fields 'A' and
                    'B' of father, followed by mother if there is one,
                    to this function. The return value is assigned to
                    fields 'A' and 'B' of the offspring. The return
                    value has to be a list even if only one field is
                    given.

"; 

%feature("docstring") simuPOP::pyTagger::~pyTagger "

Usage:

    x.~pyTagger()

"; 

%ignore simuPOP::pyTagger::pyTagger(const pyTagger &rhs);

%feature("docstring") simuPOP::pyTagger::clone "

Description:

    deep copy of a pyTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyTagger::__repr__ "

Description:

    used by Python print function to print out the general information
    of the pyTagger

Usage:

    x.__repr__()

"; 

%ignore simuPOP::pyTagger::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::quanTrait "

Description:

    base class of quantitative trait

Details:

    Quantitative trait is the measure of certain phenotype for given
    genotype. Quantitative trait is similar to penetrance in that the
    consequence of penetrance is binary: affected or unaffected; while
    it is continuous for quantitative trait.
    In simuPOP, different operators or functions were implemented to
    calculate quantitative traits for each individual and store the
    values in the information fields specified by the user (default to
    qtrait). The quantitative trait operators also accept the
    ancestralGen parameter to control the number of generations for
    which the qtrait information field will be set.

"; 

%feature("docstring") simuPOP::quanTrait::quanTrait "

Description:

    create a quantitative trait operator

Usage:

    quanTrait(ancestralGen=-1, stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=[], subPop=[], infoFields=[\"qtrait\"])

"; 

%feature("docstring") simuPOP::quanTrait::~quanTrait "

Description:

    destructor

Usage:

    x.~quanTrait()

"; 

%feature("docstring") simuPOP::quanTrait::clone "

Description:

    deep copy of a quantitative trait operator

Usage:

    x.clone()

"; 

%ignore simuPOP::quanTrait::qtrait(individual *);

%feature("docstring") simuPOP::quanTrait::apply "

Description:

    set qtrait to all individual

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::quanTrait::__repr__ "

Description:

    used by Python print function to print out the general information
    of the quantitative trait operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::randomParentChooser "

Applicability: all ploidy

Details:

    This parent chooser chooses a parent randomly from the parental
    generation. If selection is turned on, parents are chosen with
    probabilities that are proportional to their fitness values. Sex
    is not considered. Parameter replacement determines if a parent
    can be chosen multiple times. Note that selection is not allowed
    when replacement=false because this poses a particular order on
    individuals in the offspring generation.

"; 

%feature("docstring") simuPOP::randomParentChooser::randomParentChooser "

Usage:

    randomParentChooser(replacement=True)

Arguments:

    replacement:    if replacement is false, a parent can not be
                    chosen more than once.

"; 

%feature("docstring") simuPOP::randomParentChooser::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::randomParentChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::randomParentChooser::chooseParent(RawIndIterator basePtr);

%feature("docstring") simuPOP::randomParentsChooser "

Applicability: all ploidy

Details:

    This parent chooser chooses two parents randomly, a male and a
    female, from their respective sex groups randomly. If selection is
    turned on, parents are chosen from their sex groups with
    probabilities that are proportional to their fitness values. If
    replacement = False, each parent can only be used once.

"; 

%feature("docstring") simuPOP::randomParentsChooser::randomParentsChooser "

Usage:

    randomParentsChooser(replacement=True)

Details:

    Note: If selection is enabled, it works regularly on on-alpha sex,
    but works twice on alpha sex. That is to say, alphaNum alpha
    indiviudals are chosen selectively, and selected again during
    mating.

"; 

%feature("docstring") simuPOP::randomParentsChooser::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::randomParentsChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::randomParentsChooser::chooseParents(RawIndIterator basePtr);

%ignore simuPOP::randomParentsChooser::numMale();

%ignore simuPOP::randomParentsChooser::numFemale();

%feature("docstring") simuPOP::rangeSplitter "

Details:

    This class defines a splitter that groups individuals in certain
    ranges into VSPs.

"; 

%feature("docstring") simuPOP::rangeSplitter::rangeSplitter "

Usage:

    rangeSplitter(ranges)

Details:

    Create a splitter according to a number of individual ranges
    defined in ranges. For example, rangeSplitter(ranges=[[0, 20],
    [40, 50]]) defines two VSPs. The first VSP consists of individuals
    0, 1, ..., 19, and the sceond VSP consists of individuals 40, 41,
    ..., 49. Note that a nested list has to be used even if only one
    range is defined.

"; 

%feature("docstring") simuPOP::rangeSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::rangeSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::rangeSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs, which is the number of ranges defined
    in parameter ranges.

"; 

%ignore simuPOP::rangeSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::rangeSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::rangeSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of VSP vsp, which is \"Range [a, b]\" where [a, b]
    is range ranges[vsp].

"; 

%feature("docstring") simuPOP::recombinator "

Details:

    In simuPOP, only one recombinator is provided. Recombination
    events between loci a/b and b/c are independent, otherwise there
    will be some linkage between loci. Users need to specify physical
    recombination rate between adjacent loci. In addition, for the
    recombinator
    *  it only works for diploid (and for females in haplodiploid)
    populations.
    *  the recombination rate must be comprised between 0.0 and 0.5. A
    recombination rate of 0.0 means that the loci are completely
    linked, and thus behave together as a single linked locus. A
    recombination rate of 0.5 is equivalent to free of recombination.
    All other values between 0.0 and 0.5 will represent various
    linkage intensities between adjacent pairs of loci. The
    recombination rate is equivalent to 1-linkage and represents the
    probability that the allele at the next locus is randomly drawn.
    *  it works for selfing. I.e., when only one parent is provided,
    it will be recombined twice, producing both maternal and paternal
    chromosomes of the offspring.
    *  conversion is allowed. Note that conversion will nullify many
    recombination events, depending on the parameters chosen.

"; 

%feature("docstring") simuPOP::recombinator::recombinator "

Description:

    recombine chromosomes from parents

Usage:

    recombinator(intensity=-1, rate=[], loci=[], convProb=0,
      convMode=CONVERT_NumMarkers, convParam=1., begin=0, end=-1,
      step=1, at=[], rep=[], subPop=[], infoFields=[])

Arguments:

    intensity:      intensity of recombination. The actual
                    recombination rate between two loci is determined
                    by intensity*locus distance (between them).
    rate:           recombination rate regardless of locus distance
                    after all afterLoci. It can also be an array of
                    recombination rates. Should have the same length
                    as afterLoci or totNumOfLoci(). The recombination
                    rates are independent of locus distance.
    afterLoci:      an array of locus indexes. Recombination will
                    occur after these loci. If rate is also specified,
                    they should have the same length. Default to all
                    loci (but meaningless for those loci located at
                    the end of a chromosome). If this parameter is
                    given, it should be ordered, and can not include
                    loci at the end of a chromosome.
    maleIntensity:  recombination intensity for male individuals. If
                    given, parameter intensity will be considered as
                    female intensity.
    maleRate:       recombination rate for male individuals. If given,
                    parameter rate will be considered as female
                    recombination rate.
    maleAfterLoci:  if given, males will recombine at different
                    locations.
    convProb:       The probability of conversion event among all
                    recombination events. When a recombination event
                    happens, it may become a recombination event if
                    the Holliday junction is resolved/repaired
                    successfully, or a conversion event if the
                    junction is not resolved/repaired. The default
                    convProb is 0, meaning no conversion event at all.
                    Note that the ratio of conversion to recombination
                    events varies greatly from study to study, ranging
                    from 0.1 to 15 (Chen et al, Nature Review
                    Genetics, 2007). This translate to 0.1/0.9~0.1 to
                    15/16~0.94 of this parameter. When convProb is 1,
                    all recombination events will be conversion
                    events.
    convMode:       conversion mode, determines how track length is
                    determined.
                    * CONVERT_NumMarkers Converts a fixed number of
                    markers.
                    * CONVERT_GeometricDistribution An geometric
                    distribution is used to determine how many markers
                    will be converted.
                    * CONVERT_TractLength Converts a fixed length of
                    tract.
                    * CONVERT_ExponentialDistribution An exponential
                    distribution with parameter convLen will be used
                    to determine track length.
    convParam:      Parameter for the conversion process. The exact
                    meaning of this parameter is determined by
                    convMode. Note that
                    * conversion tract length is usually short, and is
                    estimated to be between 337 and 456 bp, with
                    overall range between maybe 50 - 2500 bp.
                    * simuPOP does not impose a unit for marker
                    distance so your choice of convParam needs to be
                    consistent with your unit. In the HapMap dataset,
                    cM is usually assumed and marker distances are
                    around 10kb (0.001cM ~- 1kb). Gene conversion can
                    largely be ignored. This is important when you use
                    distance based conversion mode such as
                    CONVERT_TrackLength or
                    CONVERT_ExponentialDistribution.
                    * After a track length is determined, if a second
                    recombination event happens within this region,
                    the track length will be shortened. Note that
                    conversion is identical to double recombination
                    under this context.
    haplodiploid:   If set to true, the first copy of paternal
                    chromosomes is copied directly as the paternal
                    chromosomes of the offspring. This is because
                    haplodiploid male has only one set of chromosome.

Note:

    There is no recombination between sex chromosomes of male
    individuals if sexChrom()=True. This may change later if the
    exchanges of genes between pseudoautosomal regions of XY need to
    be modeled.

"; 

%feature("docstring") simuPOP::recombinator::~recombinator "

Usage:

    x.~recombinator()

"; 

%feature("docstring") simuPOP::recombinator::clone "

Description:

    deep copy of a recombinator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::recombinator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the recombinator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::recombinator::recCounts "

Description:

    return recombination counts (only valid in standard modules)

Usage:

    x.recCounts()

"; 

%feature("docstring") simuPOP::recombinator::convCount "

Description:

    return the count of conversion of a certain size (only valid in
    standard modules)

Usage:

    x.convCount(size)

"; 

%feature("docstring") simuPOP::recombinator::convCounts "

Description:

    return the count of conversions of all sizes (only valid in
    standard modules)

Usage:

    x.convCounts()

"; 

%feature("docstring") simuPOP::recombinator::initialize "

Usage:

    x.initialize(pop)

"; 

%feature("docstring") simuPOP::recombinator::transmitGenotype "

Usage:

    x.transmitGenotype(parent, offspring, ploidy)

"; 

%ignore simuPOP::recombinator::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad, individual *mom);

%feature("docstring") simuPOP::repList "

Details:

    A class to specify replicate list. The reason why I cannot simple
    use vectori() is that users have got used to use a single number
    to specify a single replicate.

"; 

%feature("docstring") simuPOP::repList::repList "

Usage:

    repList(reps=[])

"; 

%feature("docstring") simuPOP::repList::match "

Usage:

    x.match(rep, numRep)

"; 

%feature("docstring") simuPOP::resizeSubPops "

Function form:

    ResizeSubPops

Description:

    resize subpopulations

Details:

    This operator resize subpopulations subPops to a another size. If
    subPops is ignored, all subpopulations will be resized. If the new
    size is smaller than the original one, the remaining individuals
    are discarded. If the new size if greater, individuals will be
    copied again if propagate is true, and be empty otherwise.

"; 

%feature("docstring") simuPOP::resizeSubPops::resizeSubPops "

Description:

    resize subpopulations

Usage:

    resizeSubPops(newSizes=[], subPops=[], propagate=True,
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Arguments:

    newSizes:       of the specified (or all) subpopulations.
    subPops:        subpopulations to be resized. Default to all.
    propagate:      if true (default) and the new size if greater than
                    the original size, individuals will be copied
                    over.

"; 

%feature("docstring") simuPOP::resizeSubPops::~resizeSubPops "

Description:

    destructor

Usage:

    x.~resizeSubPops()

"; 

%feature("docstring") simuPOP::resizeSubPops::clone "

Description:

    deep copy of a resizeSubPops operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::resizeSubPops::apply "

Description:

    apply a resizeSubPops operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::resizeSubPops::__repr__ "

Description:

    used by Python print function to print out the general information
    of the resizeSubPops operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::RNG "

Description:

    random number generator

Details:

    This random number generator class wraps around a number of random
    number generators from GNU Scientific Library. You can obtain and
    change system random number generator through the rng() function.
    Or create a separate random number generator and use it in your
    script.

"; 

%feature("docstring") simuPOP::RNG::RNG "

Description:

    RNG used by simuPOP.

Usage:

    RNG(rng=None, seed=0)

"; 

%feature("docstring") simuPOP::RNG::~RNG "

Usage:

    x.~RNG()

"; 

%feature("docstring") simuPOP::RNG::setRNG "

Description:

    choose an random number generator, or set seed to the current RNG

Usage:

    x.setRNG(rng=None, seed=0)

Arguments:

    rng:            name of the RNG. If rng is not given,
                    environmental variable GSL_RNG_TYPE will be used
                    if it is available. Otherwise, RNGmt19937 will be
                    used.
    seed:           random seed. If not given, /dev/urandom,
                    /dev/random, system time will be used, depending
                    on availability, in that order. Note that windows
                    system does not have /dev so system time is used.

"; 

%feature("docstring") simuPOP::RNG::name "

Description:

    return RNG name

Usage:

    x.name()

"; 

%feature("docstring") simuPOP::RNG::seed "

Description:

    return the seed of this RNG

Usage:

    x.seed()

"; 

%feature("docstring") simuPOP::RNG::maxSeed "

Description:

    return the maximum allowed seed value

Usage:

    x.maxSeed()

"; 

%feature("docstring") simuPOP::RNG::setSeed "

Description:

    if seed is 0, method described in setRNG is used.

Usage:

    x.setSeed(seed)

"; 

%ignore simuPOP::RNG::generateRandomSeed();

%feature("docstring") simuPOP::RNG::max "

Description:

    Maximum value of this RNG.

Usage:

    x.max()

"; 

%feature("docstring") simuPOP::RNG::__repr__ "

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::RNG::randGet "

Description:

    return a random number in the range of [0, 2, ... max()-1]

Usage:

    x.randGet()

"; 

%feature("docstring") simuPOP::RNG::randBit "

Usage:

    x.randBit()

"; 

%feature("docstring") simuPOP::RNG::randInt "

Description:

    return a random number in the range of [0, 1, 2, ... n-1]

Usage:

    x.randInt(n)

"; 

%ignore simuPOP::RNG::randIntArray(ULONG n, ULONG size, ULONG *vec);

%feature("docstring") simuPOP::RNG::randGeometric "

Description:

    Geometric distribution.

Usage:

    x.randGeometric(p)

"; 

%feature("docstring") simuPOP::RNG::randUniform01 "

Description:

    Uniform distribution [0,1).

Usage:

    x.randUniform01()

"; 

%feature("docstring") simuPOP::RNG::randNormal "

Description:

    Normal distribution.

Usage:

    x.randNormal(m, v)

"; 

%feature("docstring") simuPOP::RNG::randExponential "

Usage:

    x.randExponential(v)

"; 

%ignore simuPOP::RNG::randUniform01Array(ULONG size, double *vec);

%feature("docstring") simuPOP::RNG::randBinomial "

Description:

    Binomial distribution B(n, p).

Usage:

    x.randBinomial(n, p)

"; 

%feature("docstring") simuPOP::RNG::randMultinomial "

Description:

    Multinomial distribution.

Usage:

    x.randMultinomial(N, p, n)

"; 

%feature("docstring") simuPOP::RNG::randMultinomialVal "

Usage:

    x.randMultinomialVal(N, p)

"; 

%feature("docstring") simuPOP::RNG::randPoisson "

Description:

    Poisson distribution.

Usage:

    x.randPoisson(p)

"; 

%feature("docstring") simuPOP::RNG::pvalChiSq "

Description:

    right hand side (single side) p-value for ChiSq value

Usage:

    x.pvalChiSq(chisq, df)

"; 

%feature("docstring") simuPOP::savePopulation "

Description:

    save population to a file

"; 

%feature("docstring") simuPOP::savePopulation::savePopulation "

Description:

    save population

Usage:

    savePopulation(output=\"\", outputExpr=\"\", format=\"\",
      compress=True, stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=[], subPop=[], infoFields=[])

Arguments:

    output:         output filename.
    outputExpr:     An expression that will be evalulated dynamically
                    to determine file name. Parameter output will be
                    ignored if this parameter is given.
    format:         obsolete parameter
    compress:       obsolete parameter

"; 

%feature("docstring") simuPOP::savePopulation::~savePopulation "

Usage:

    x.~savePopulation()

"; 

%feature("docstring") simuPOP::savePopulation::clone "

Usage:

    x.clone()

Details:

    Return a cloned copy of an operator. This function is available to
    all operators.

"; 

%feature("docstring") simuPOP::savePopulation::apply "

Usage:

    x.apply(pop)

Details:

    Apply an operator to population pop directly, without checking its
    applicability.

"; 

%feature("docstring") simuPOP::savePopulation::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::selector "

Description:

    A base selection operator for all selectors.

Details:

    Genetic selection is tricky to simulate since there are many
    different fitness values and many different ways to apply
    selection. simuPOP employs an 'ability-to-mate' approach. Namely,
    the probability that an individual will be chosen for mating is
    proportional to its fitness value. More specifically,
    *  PreMating selectors assign fitness values to each individual,
    and mark part or all subpopulations as under selection.
    *  during sexless mating (e.g. binomialSelection mating scheme),
    individuals are chosen at probabilities that are proportional to
    their fitness values. If there are $ N $ individuals with fitness
    values $ f_{i},i=1,...,N $, individual $ i $ will have probability
    $ \\frac{f_{i}}{\\sum_{j}f_{j}} $ to be chosen and passed to the
    next generation.
    *  during randomMating, males and females are separated. They are
    chosen from their respective groups in the same manner as
    binomialSelection and mate.
    All of the selection operators, when applied, will set an
    information field fitness (configurable) and then mark part or all
    subpopulations as under selection. (You can use different
    selectors to simulate various selection intensities for different
    subpopulations). Then, a 'selector-aware' mating scheme can select
    individuals according to their fitness information fields. This
    implies that
    *  only mating schemes can actually select individuals.
    *  a selector has to be a PreMating operator. This is not a
    problem when you use the operator form of the selector since its
    default stage is PreMating. However, if you use the function form
    of the selector in a pyOperator, make sure to set the stage of
    pyOperator to PreMating.

Note:

    You can not apply two selectors to the same subpopulation, because
    only one fitness value is allowed for each individual.

"; 

%feature("docstring") simuPOP::selector::selector "

Description:

    create a selector

Usage:

    selector(subPops=[], stage=PreMating, begin=0, end=-1, step=1,
      at=[], rep=[], subPop=[], infoFields=[\"fitness\"])

Arguments:

    subPop:         a shortcut to subPops=[subPop]
    subPops:        subpopulations that the selector will apply to.
                    Default to all.

"; 

%feature("docstring") simuPOP::selector::~selector "

Description:

    destructor

Usage:

    x.~selector()

"; 

%feature("docstring") simuPOP::selector::clone "

Description:

    deep copy of a selector

Usage:

    x.clone()

"; 

%ignore simuPOP::selector::indFitness(individual *, ULONG gen);

%feature("docstring") simuPOP::selector::apply "

Description:

    set fitness to all individuals. No selection will happen!

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::selector::__repr__ "

Description:

    used by Python print function to print out the general information
    of the selector

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::selfingGenoTransmitter "

Applicability: diploid only

Details:

    selfing offspring generator works similarly as a mendelian
    offspring generator but a single parent produces both the paternal
    and maternal copy of the offspring chromosomes. This offspring
    generator accepts a dipload parent. A random copy of the parental
    chromosomes is chosen randomly to form the parental copy of the
    offspring chromosome, and is chosen randomly again to form the
    maternal copy of the offspring chromosome.

"; 

%feature("docstring") simuPOP::selfingGenoTransmitter::selfingGenoTransmitter "

Usage:

    selfingGenoTransmitter(begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::selfingGenoTransmitter::clone "

Usage:

    x.clone()

Details:

    Return a cloned copy of an operator. This function is available to
    all operators.

"; 

%feature("docstring") simuPOP::selfingGenoTransmitter::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::selfingGenoTransmitter::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::sequentialParentChooser "

Applicability: all ploidy

Details:

    This parent chooser chooses a parent linearly, regardless of sex
    or fitness values (selection is not considered).

"; 

%feature("docstring") simuPOP::sequentialParentChooser::sequentialParentChooser "

Usage:

    sequentialParentChooser()

"; 

%feature("docstring") simuPOP::sequentialParentChooser::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::sequentialParentChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::sequentialParentChooser::chooseParent(RawIndIterator basePtr);

%feature("docstring") simuPOP::sequentialParentsChooser "

Applicability: all ploidy

Details:

    This parents chooser chooses two parents sequentially. The parents
    are chosen from their respective sex groups. Selection is not
    considered.

"; 

%feature("docstring") simuPOP::sequentialParentsChooser::sequentialParentsChooser "

Usage:

    sequentialParentsChooser()

"; 

%feature("docstring") simuPOP::sequentialParentsChooser::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::sequentialParentsChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::sequentialParentsChooser::chooseParents(RawIndIterator basePtr);

%ignore simuPOP::sequentialParentsChooser::numMale();

%ignore simuPOP::sequentialParentsChooser::numFemale();

%feature("docstring") simuPOP::setAncestralDepth "

Description:

    set ancestral depth

Details:

    This operator set the number of ancestral generations to keep in a
    population. It is usually called like setAncestral(at=[-2]) to
    start recording ancestral generations to a population at the end
    of the evolution. This is useful when constructing pedigree trees
    from a population.

"; 

%feature("docstring") simuPOP::setAncestralDepth::setAncestralDepth "

Description:

    create a setAncestralDepth operator

Usage:

    setAncestralDepth(depth, output=\">\", outputExpr=\"\",
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::setAncestralDepth::~setAncestralDepth "

Description:

    destructor

Usage:

    x.~setAncestralDepth()

"; 

%feature("docstring") simuPOP::setAncestralDepth::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::setAncestralDepth::apply "

Description:

    apply the setAncestralDepth operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::setAncestralDepth::__repr__ "

Description:

    used by Python print function to print out the general information
    of the setAncestralDepth operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::sexSplitter "

Details:

    This splitter defines two VSPs by individual sex. The first VSP
    consists of all male individuals and the second VSP consists of
    all females in a subpopulation.

"; 

%feature("docstring") simuPOP::sexSplitter::sexSplitter "

Description:

    Create a sex splitter that defines male and female VSPs.

Usage:

    sexSplitter()

"; 

%feature("docstring") simuPOP::sexSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::sexSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::sexSplitter::numVirtualSubPop "

Description:

    Return 2.

Usage:

    x.numVirtualSubPop()

"; 

%ignore simuPOP::sexSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::sexSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::sexSplitter::name "

Description:

    Return \"Male\" if vsp=0 and \"Female\" otherwise.

Usage:

    x.name(vsp)

"; 

%ignore simuPOP::SharedVariables;

%ignore simuPOP::SharedVariables::SharedVariables();

%ignore simuPOP::SharedVariables::SharedVariables(const SharedVariables &rhs);

%ignore simuPOP::SharedVariables::swap(SharedVariables &rhs);

%feature("docstring") simuPOP::SharedVariables::~SharedVariables "

Description:

    hence call this destructore.

Usage:

    x.~SharedVariables()

"; 

%ignore simuPOP::SharedVariables::clear();

%ignore simuPOP::SharedVariables::setDict(PyObject *dict);

%ignore simuPOP::SharedVariables::setVar(const string &name, const PyObject *val);

%ignore simuPOP::SharedVariables::getVar(const string &name, bool nameError=true);

%feature("docstring") simuPOP::SharedVariables::hasVar "

Usage:

    x.hasVar(name)

"; 

%feature("docstring") simuPOP::SharedVariables::removeVar "

Description:

    remove variable

Usage:

    x.removeVar(name)

"; 

%ignore simuPOP::SharedVariables::setBoolVar(const string &name, const bool val);

%ignore simuPOP::SharedVariables::setIntVar(const string &name, const int val);

%ignore simuPOP::SharedVariables::setDoubleVar(const string &name, const double val);

%ignore simuPOP::SharedVariables::setStringVar(const string &name, const string &val);

%ignore simuPOP::SharedVariables::setIntVectorVar(const string &name, const vectori &val);

%ignore simuPOP::SharedVariables::setDoubleVectorVar(const string &name, const vectorf &val);

%ignore simuPOP::SharedVariables::setStrDictVar(const string &name, const strDict &val);

%ignore simuPOP::SharedVariables::setIntDictVar(const string &name, const intDict &val);

%ignore simuPOP::SharedVariables::getVarAsBool(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsInt(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsDouble(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsString(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsStrDict(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsIntDict(const string &name, bool nameError=true);

%feature("docstring") simuPOP::SharedVariables::dict "

Usage:

    x.dict()

"; 

%ignore simuPOP::SharedVariables::asString() const;

%feature("docstring") simuPOP::SharedVariables::fromString "

Usage:

    x.fromString(vars)

"; 

%feature("docstring") simuPOP::simulator "

Details:

    A simuPOP simulator is responsible for evolving one or more
    replicates of a population forward in time, subject to various
    operators. Populations in a simulator are created as identical
    copies of a population and will become different after evolution.
    A mating scheme needs to be specified, which will be used to
    generate offspring generations during evolution. A number of
    functions are provided to access simulator properties, access
    populations and their variables, copy, save and load a simulator.
    The most important member function of a simulator is evolve, which
    evolves populations forward in time, subject to various operators.
    Because populations in a simulator have to keep the same genotypic
    structure, several functions are provided to change ancestral
    depth and information fields of all populations. These functions
    cannot be replaced by similar calls to all populations in a
    simulator because the genotypic structure of the simulator itself
    needs to be updated.

"; 

%feature("docstring") simuPOP::simulator::simulator "

Usage:

    simulator(pop, matingScheme, rep=1)

Details:

    Create a simulator with rep replicates of population pop.
    Population pop will be copied rep times (default to 1), while
    keeping the passed population intact. A mating scheme matingScheme
    will be used to evolve these populations.

"; 

%feature("docstring") simuPOP::simulator::~simulator "

Usage:

    x.~simulator()

"; 

%ignore simuPOP::simulator::simulator(const simulator &rhs);

%feature("docstring") simuPOP::simulator::clone "

Usage:

    x.clone()

Details:

    Clone a simulator, along with all its populations. Note that
    Python assign statement simu1 = simu only creates a symbolic link
    to an existing simulator.

"; 

%feature("docstring") simuPOP::simulator::gen "

Usage:

    x.gen()

Details:

    Return the current generation number, which is the initial
    generation number (0, or some value set by setGen(gen)) plus the
    total number of generations evolved.

"; 

%feature("docstring") simuPOP::simulator::setGen "

Usage:

    x.setGen(gen)

Details:

    Set the current generation number of a simulator to gen.

"; 

%feature("docstring") simuPOP::simulator::numRep "

Usage:

    x.numRep()

Details:

    Return the number of replicates.

"; 

%feature("docstring") simuPOP::simulator::pop "

Usage:

    x.pop(rep)

Details:

    Return a reference to the rep-th population of a simulator. The
    reference will become invalid once the simulator starts evolving
    or becomes invalid (removed). Modifying the returned object is
    discouraged because it will change the population within the
    simulator. If an independent copy of the population is needed, use
    simu.population(rep). clone().

"; 

%feature("docstring") simuPOP::simulator::extract "

Usage:

    x.extract(rep)

Details:

    Extract the rep-th population from a simulator. This will reduce
    the number of populations in this simulator by one.

"; 

%feature("docstring") simuPOP::simulator::populations "

Usage:

    x.populations()

Details:

    Return a Python iterator that can be used to iterate through all
    populations in a simulator.

"; 

%ignore simuPOP::simulator::setPopulation(population &pop, UINT rep);

%feature("docstring") simuPOP::simulator::evolve "

Usage:

    x.evolve(ops, preOps=[], postOps=[], gen=-1, dryrun=False)

Details:

    Evolve all populations gen generations, subject to operators
    opspreOps and postOps. Operators preOps are applied to all
    populations (subject to applicability restrictions of the
    operators, imposed by the rep parameter of these operators) before
    evolution. They are usually used to initialize populations.
    Operators postOps are applied to all populations after the
    evolution.
    Operators ops are applied during the life cycle of each
    generation. Depending on the stage of these operators, they can be
    applied before-, during-, and/or post-mating. These operators can
    be applied at all or some of the generations, depending the begin,
    end, step, and at parameters of these operators. Populations in a
    simulator are evolved one by one. At each generation, the
    applicability of these operators are determined. Pre-mating
    operators are applied to a population first. A mating scheme is
    then used to populate an offspring generation, using applicable
    during-mating operators. After an offspring generation is
    successfully generated and becomes the current generation,
    applicable post-mating operators are applied to it. Because the
    order at which operators are applied can be important, and the
    stage(s) at which operators are applied are not always clear, a
    parameter dryRun can be used. If set to True, this function will
    print out the order at which all operators are applied, without
    actually evolving the populations.
    Parameter gen can be set to a positive number, which is the number
    of generations to evolve. If gen is negative (default), the
    evolution will continue indefinitely, until all replicates are
    stopped by a special kind of operators called terminators. At the
    end of the evolution, the generations that each replicates have
    evolved are returned.

"; 

%ignore simuPOP::simulator::apply(const vectorop ops, bool dryrun=false);

%feature("docstring") simuPOP::simulator::addInfoField "

Usage:

    x.addInfoField(field, init=0)

Details:

    Add an information field field to all populations in a simulator,
    and update the genotypic structure of the simulator itself. The
    information field will be initialized by value init.

"; 

%feature("docstring") simuPOP::simulator::addInfoFields "

Usage:

    x.addInfoFields(fields, init=0)

Details:

    Add information fields fields to all populations in a simulator,
    and update the genotypic structure of the simulator itself. The
    information field will be initialized by value init.

"; 

%feature("docstring") simuPOP::simulator::setAncestralDepth "

Usage:

    x.setAncestralDepth(depth)

Details:

    Set ancestral depth of all populations in a simulator.

"; 

%feature("docstring") simuPOP::simulator::setMatingScheme "

Usage:

    x.setMatingScheme(matingScheme)

Details:

    Set a new mating scheme matingScheme to a simulator.

"; 

%feature("docstring") simuPOP::simulator::vars "

Usage:

    x.vars(rep)

Details:

    Return the local namespace of the rep-th population, equivalent to
    x.population(rep). vars().

"; 

%feature("docstring") simuPOP::simulator::vars "

Usage:

    x.vars(rep, subPop)

Details:

    Return a dictionary of subpopulation variables in the local
    namespace of the rep-th population, equivalent to
    x.population(rep).vars(subPop).

"; 

%feature("docstring") simuPOP::simulator::save "

Usage:

    x.save(filename)

Details:

    Save a simulator to file filename, which can be loaded by a global
    function LoadSimulator.

"; 

%ignore simuPOP::simulator::load(string filename);

%feature("docstring") simuPOP::simulator::__repr__ "

Description:

    of the simulator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::smmMutator "

Function form:

    SmmMutate

Description:

    The stepwise mutation model.

Details:

    The Stepwise Mutation Model (SMM) assumes that alleles are
    represented by integer values and that a mutation either increases
    or decreases the allele value by one. For variable number tandem
    repeats(VNTR) loci, the allele value is generally taken as the
    number of tandem repeats in the DNA sequence.

"; 

%feature("docstring") simuPOP::smmMutator::smmMutator "

Description:

    create a SMM mutator

Usage:

    smmMutator(rate=[], loci=[], maxAllele=0, incProb=0.5,
      output=\">\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=[], subPop=[], infoFields=[])

Details:

    The SMM is developed for allozymes. It provides better description
    for these kinds of evolutionary processes.  Please see class
    mutator for the descriptions of other parameters.

Arguments:

    incProb:        probability to increase allele state. Default to
                    0.5.

"; 

%feature("docstring") simuPOP::smmMutator::~smmMutator "

Usage:

    x.~smmMutator()

"; 

%ignore simuPOP::smmMutator::mutate(AlleleRef allele);

%feature("docstring") simuPOP::smmMutator::clone "

Description:

    deep copy of a smmMutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::smmMutator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the smmMutator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::splitSubPop "

Function form:

    SplitSubPop

Description:

    split a subpopulation

Details:

    

"; 

%feature("docstring") simuPOP::splitSubPop::splitSubPop "

Description:

    split a subpopulation

Usage:

    splitSubPop(which=0, sizes=[], proportions=[], randomize=True,
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[\"migrate_to\"])

Details:

    Split a subpopulation by sizes or proportions. Individuals are
    randomly (by default) assigned to the resulting subpopulations.
    Because mating schemes may introduce certain order to individuals,
    randomization ensures that split subpopulations have roughly even
    distribution of genotypes.

Arguments:

    which:          which subpopulation to split. If there is no
                    subpopulation structure, use 0 as the first (and
                    only) subpopulation.
    sizes:          new subpopulation sizes. The sizes should be added
                    up to the original subpopulation (subpopulation
                    which) size.
    proportions:    proportions of new subpopulations. Should be added
                    up to 1.
    randomize:      Whether or not randomize individuals before
                    population split. Default to True.

"; 

%feature("docstring") simuPOP::splitSubPop::~splitSubPop "

Description:

    destructor

Usage:

    x.~splitSubPop()

"; 

%feature("docstring") simuPOP::splitSubPop::clone "

Description:

    deep copy of a splitSubPop operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::splitSubPop::apply "

Description:

    apply a splitSubPop operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::splitSubPop::__repr__ "

Description:

    used by Python print function to print out the general information
    of the splitSubPop operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::stat "

Function form:

    Stat

Description:

    calculate statistics

Details:

    Operator stat calculates various basic statistics for the
    population and sets variables in the local namespace. Other
    operators or functions can refer to the results from the namespace
    after stat is applied. Stat is the function form of the operator.
    Note that these statistics are dependent to each other. For
    example, heterotype and allele frequencies of related loci will be
    automatically calculated if linkage diseqilibrium is requested.

"; 

%feature("docstring") simuPOP::stat::stat "

Description:

    create an stat operator

Usage:

    stat(popSize=False, numOfMale=False, numOfMale_param={},
      numOfAffected=False, numOfAffected_param={}, numOfAlleles=[],
      numOfAlleles_param={}, alleleFreq=[], alleleFreq_param={},
      heteroFreq=[], expHetero=[], expHetero_param={}, homoFreq=[],
      genoFreq=[], genoFreq_param={}, haploFreq=[], LD=[],
      LD_param={}, association=[], association_param={}, Fst=[],
      Fst_param={}, relGroups=[], relLoci=[], rel_param={},
      relBySubPop=False, relMethod=[], relMinScored=10,
      hasPhase=False, midValues=False, output=\"\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=[],
      subPop=[], infoFields=[])

Arguments:

    popSize:        whether or not calculate population and virtual
                    subpopulation sizes. This parameter will set the
                    following variables:
                    * numSubPop the number of subpopulations.
                    * subPopSize an array of subpopulation sizes.
                    * virtualSubPopSize (optional) an array of virtual
                    subpopulation sizes. If a subpopulation does not
                    have any virtual subpopulation, the subpopulation
                    size is returned.
                    * popSize, subPop[sp]['popSize'] the
                    population/subpopulation size.
    numOfMale:      whether or not count the numbers or proportions of
                    males and females. This parameter can set the
                    following variables by user's specification:
                    * numOfMale, subPop[sp]['numOfMale'] the number of
                    males in the population/subpopulation.
                    * numOfFemale, subPop[sp]['numOfFemale'] the
                    number of females in the population/subpopulation.
                    * propOfMale, subPop[sp]['propOfMale'] the
                    proportion of males in the
                    population/subpopulation.
                    * propOfFemale, subPop[sp]['propOfFemale'] the
                    proportion of females in the
                    population/subpopulation.
    numOfMale_param:a dictionary of parameters of numOfMale
                    statistics. Can be one or more items choosen from
                    the following options: numOfMale, propOfMale,
                    numOfFemale, and propOfFemale.
    numOfAffected:  whether or not count the numbers or proportions of
                    affected and unaffected individuals. This
                    parameter can set the following variables by
                    user's specification:
                    * numOfAffected, subPop[sp]['numOfAffected'] the
                    number of affected individuals in the
                    population/subpopulation.
                    * numOfUnaffected, subPop[sp]['numOfUnAffected']
                    the number of unaffected individuals in the
                    population/subpopulation.
                    * propOfAffected, subPop[sp]['propOfAffected'] the
                    proportion of affected individuals in the
                    population/subpopulation.
                    * propOfUnaffected, subPop[sp]['propOfUnAffected']
                    the proportion of unaffected individuals in the
                    population/subpopulation.
    numOfAffected_param:a dictionary of parameters of numOfAffected
                    statistics. Can be one or more items choosen from
                    the following options: numOfAffected,
                    propOfAffected, numOfUnaffected, propOfUnaffected.
    numOfAlleles:   an array of loci at which the numbers of distinct
                    alleles will be counted (numOfAlleles=[loc1, loc2,
                    ...] where loc1 etc. are absolute locus indexes).
                    This is done through the calculation of allele
                    frequencies. Therefore, allele frequencies will
                    also be calculated if this statistics is
                    requested. This parameter will set the following
                    variables (carray objects of the numbers of
                    alleles for all loci). Unrequested loci will have
                    0 distinct alleles.
                    * numOfAlleles, subPop[sp]['numOfAlleles'] the
                    number of distinct alleles at each locus.
                    (Calculated only at requested loci.)
    numOfAlleles_param:a dictionary of parameters of numOfAlleles
                    statistics. Can be one or more items choosen from
                    the following options: numOfAffected,
                    propOfAffected, numOfUnaffected, propOfUnaffected.
    alleleFreq:     an array of loci at which all allele frequencies
                    will be calculated (alleleFreq=[loc1, loc2, ...]
                    where loc1 etc. are loci where allele frequencies
                    will be calculated). This parameter will set the
                    following variables (carray objects); for example,
                    alleleNum[1][2] will be the number of allele 2 at
                    locus 1:
                    * alleleNum[a], subPop[sp]['alleleNum'][a]
                    * alleleFreq[a], subPop[sp]['alleleFreq'][a].
    alleleFreq_param:a dictionary of parameters of alleleFreq
                    statistics. Can be one or more items choosen from
                    the following options: numOfAlleles, alleleNum,
                    and alleleFreq.
    genoFreq:       an array of loci at which all genotype frequencies
                    will be calculated (genoFreq=[loc1, loc2, ...].
                    You may use parameter genoFreq_param to control if
                    a/b and b/a are the same genotype. This parameter
                    will set the following dictionary variables. Note
                    that unlike list used for alleleFreq etc., the
                    indexes a, b of genoFreq[loc][a][b] are dictionary
                    keys, so you will get a KeyError when you used a
                    wrong key. You can get around this problem by
                    using expressions like
                    genoNum[loc].setDefault(a,{}).
                    * genoNum[loc][allele1][allele2] and
                    subPop[sp]['genoNum'][loc][allele1][allele2], the
                    number of genotype allele1-allele2 at locus loc.
                    * genoFreq[loc][allele1][allele2] and
                    subPop[sp]['genoFreq'][loc][allele1][allele2], the
                    frequency of genotype allele1-allele2 at locus
                    loc.
                    * genoFreq_param a dictionary of parameters of
                    phase = 0 or 1.
    heteroFreq:     an array of loci at which observed
                    heterozygosities will be calculated
                    (heteroFreq=[loc1, loc2, ...]). For each locus,
                    the number and frequency of allele specific and
                    overall heterozygotes will be calcuated and stored
                    in four population variables. For example,
                    heteroNum[loc][1] stores number of heterozygotes
                    at locus loc, with respect to allele 1, which is
                    the number of all genotype 1x or x1 where does not
                    equal to 1. All other genotypes such as 02 are
                    considered as homozygotes when heteroFreq[loc][1]
                    is calculated. The overall number of heterozygotes
                    (HeteroNum[loc]) is the number of genotype xy if x
                    does not equal to y.
                    * HeteroNum[loc], subPop[sp]['HeteroNum'][loc],
                    the overall heterozygote count.
                    * HeteroFreq[loc], subPop[sp]['HeteroFreq'][loc],
                    the overall heterozygote frequency.
                    * heteroNum[loc][allele],
                    subPop[sp]['heteroNum'][loc][allele], allele-
                    specific heterozygote counts.
                    * heteroFreq[loc][allele],
                    subPop[sp]['heteroFreq'][loc][allele], allele-
                    specific heterozygote frequency.
    homoFreq:       an array of loci to calculate observed
                    homozygosities and expected homozygosities
                    (homoFreq=[loc1, loc2, ...]). This parameter will
                    calculate the numbers and frequencies of
                    homozygotes xx and set the following variables:
                    * homoNum[loc], subPop[sp]['homoNum'][loc].
                    * homoFreq[loc], subPop[sp]['homoFreq'][loc].
    expHetero:      an array of loci at which the expected
                    heterozygosities will be calculated
                    (expHetero=[loc1, loc2, ...]). The expected
                    heterozygosity is calculated by
                    $ h_{exp}=1-p_{i}^{2}, $ where $ p_i $ is the
                    allele frequency of allele $ i $. The following
                    variables will be set:
                    * expHetero[loc], subPop[sp]['expHetero'][loc].
    expHetero_param:a dictionary of parameters of expHetero
                    statistics. Can be one or more items choosen from
                    the following options: subpop and midValues.
    haploFreq:      a matrix of haplotypes (allele sequences on
                    different loci) to count. For example, haploFreq =
                    [ [ 0,1,2 ], [1,2] ] will count all haplotypes on
                    loci 0, 1 and 2; and all haplotypes on loci 1, 2.
                    If only one haplotype is specified, the outer []
                    can be omitted. I.e., haploFreq=[0,1] is
                    acceptable. The following dictionary variables
                    will be set with keys 0-1-2 etc. For example,
                    haploNum['1-2']['5-6'] is the number of allele
                    pair 5, 6 (on loci 1 and 2 respectively) in the
                    population.
                    * haploNum[haplo] and
                    subPop[sp]['haploNum'][haplo], the number of
                    allele sequencies on loci haplo.
                    * haploFreq[haplo],
                    subPop[sp]['haploFreq'][haplo], the frequency of
                    allele sequencies on loci haplo.
    LD:             calculate linkage disequilibria $ LD $, $ LD' $
                    and $ r^{2} $, given LD=[ [loc1, loc2], [ loc1,
                    loc2, allele1, allele2], ... ]. For each item
                    [loc1, loc2, allele1, allele2], $ D $, $ D' $ and
                    $ r^{2} $ will be calculated based on allele1 at
                    loc1 and allele2 at loc2. If only two loci are
                    given, the LD values are averaged over all allele
                    pairs. For example, for allele $ A $ at locus 1
                    and allele $ B $ at locus 2,
                    $ D = P_{AB}-P_{A}P_{B} $
                    $ D' = D/D_{max} $
                    $ D_{max} = \\min\\left(P_{A}\\left(1-P_{B}\\right),\\l
                    eft(1-P_{A}\\right)P_{B}\\right) \\textrm{if }D>0 \\\\ \\
                    min\\left(P_{A}P_{B},\\left(1-P_{A}\\right)\\left(1-P_
                    {B}\\right)\\right) \\textrm{if }D<0 $
                    $ r^{2} = \\frac{D^{2}}{P_{A}\\left(1-P_{A}\\right)P_
                    {B}\\left(1-P_{B}\\right)} $ If only one item is
                    specified, the outer [] can be ignored. I.e.,
                    LD=[loc1, loc2] is acceptable. This parameter will
                    set the following variables. Please note that the
                    difference between the data structures used for ld
                    and LD.
                    * ld['loc1-loc2']['allele1-allele2'],
                    subPop[sp]['ld']['loc1-loc2']['allele1-allele2'].
                    * ld_prime['loc1-loc2']['allele1-allele2'], subPop
                    [sp]['ld_prime']['loc1-loc2']['allele1-allele2'].
                    * r2['loc1-loc2']['allele1-allele2'],
                    subPop[sp]['r2']['loc1-loc2']['allele1-allele2'].
                    * LD[loc1][loc2], subPop[sp]['LD'][loc1][loc2].
                    * LD_prime[loc1][loc2],
                    subPop[sp]['LD_prime'][loc1][loc2].
                    * R2[loc1][loc2], subPop[sp]['R2'][loc1][loc2].
    LD_param:       a dictionary of parameters of LD statistics. Can
                    have key stat which is a list of statistics to
                    calculate. Default to all. If any statistics is
                    specified, only those specified will be
                    calculated. For example, you may use
                    LD_param={LD_prime} to calculate D' only, where
                    LD_prime is a shortcut for 'stat':['LD_prime'].
                    Other parameters that you may use are:
                    * subPop whether or not calculate statistics for
                    subpopulations.
                    * midValues whether or not keep intermediate
                    results.
    association:    association measures
    association_param:a dictionary of parameters of association
                    statistics. Can be one or more items choosen from
                    the following options: ChiSq, ChiSq_P, UC_U, and
                    CramerV.
    Fst:            calculate $ F_{st} $, $ F_{is} $, $ F_{it} $. For
                    example, Fst = [0,1,2] will calculate $ F_{st} $,
                    $ F_{is} $, $ F_{it} $ based on alleles at loci 0,
                    1, 2. The locus-specific values will be used to
                    calculate AvgFst, which is an average value over
                    all alleles (Weir & Cockerham, 1984). Terms and
                    values that match Weir & Cockerham are:
                    * $ F $ ( $ F_{IT} $) the correlation of genes
                    within individuals (inbreeding);
                    * $ \\theta $ ( $ F_{ST} $) the correlation of
                    genes of difference individuals in the same
                    population (will evaluate for each subpopulation
                    and the whole population)
                    * $ f $ ( $ F_{IS} $) the correlation of genes
                    within individuals within populations. This
                    parameter will set the following variables:
                    * Fst[loc], Fis[loc], Fit[loc]
                    * AvgFst, AvgFis, AvgFit.
    Fst_param:      a dictionary of parameters of Fst statistics. Can
                    be one or more items choosen from the following
                    options: Fst, Fis, Fit, AvgFst, AvgFis, and
                    AvgFit.
    relMethod:      method used to calculate relatedness. Can be
                    either REL_Queller or REL_Lynch. The relatedness
                    values between two individuals, or two groups of
                    individuals are calculated according to Queller &
                    Goodnight (1989) (method=REL_Queller) and Lynch et
                    al. (1999) (method=REL_Lynch). The results are
                    pairwise relatedness values, in the form of a
                    matrix. Original group or subpopulation numbers
                    are discarded. There is no subpopulation level
                    relatedness value.
    relGroups:      calculate pairwise relatedness between groups. Can
                    be in the form of either [[1,2,3],[5,6,7],[8,9]]
                    or [2,3,4]. The first one specifies groups of
                    individuals, while the second specifies
                    subpopulations. By default, relatedness between
                    subpopulations is calculated.
    relLoci:        loci on which relatedness values are calculated
    rel_param:      a dictionary of parameters of relatedness
                    statistics. Can be one or more items choosen from
                    the following options: Fst, Fis, Fit, AvgFst,
                    AvgFis, and AvgFit.
    hasPhase:       if a/b and b/a are the same genotype. Default to
                    False.
    midValues:      whether or not post intermediate results. Default
                    to False. For example, Fst will need to calculate
                    allele frequencise. If midValues is set to True,
                    allele frequencies will be posted as well. This
                    will be helpful in debugging and sometimes in
                    deriving statistics.

"; 

%feature("docstring") simuPOP::stat::~stat "

Usage:

    x.~stat()

"; 

%feature("docstring") simuPOP::stat::clone "

Description:

    deep copy of a stat operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::stat::apply "

Description:

    apply the stat operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::stat::__repr__ "

Description:

    used by Python print function to print out the general information
    of the stat operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::statAlleleFreq;

%feature("docstring") simuPOP::statAlleleFreq::statAlleleFreq "

Usage:

    statAlleleFreq(atLoci=[], param={})

"; 

%feature("docstring") simuPOP::statAlleleFreq::~statAlleleFreq "

Description:

    destructor, nested vectors have to be cleared manually

Usage:

    x.~statAlleleFreq()

"; 

%feature("docstring") simuPOP::statAlleleFreq::addLocus "

Usage:

    x.addLocus(locus, post, subPop, numOfAlleles)

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleNumAll "

Usage:

    x.alleleNumAll()

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleNumVec "

Usage:

    x.alleleNumVec(loc)

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleNum "

Usage:

    x.alleleNum(allele, loc)

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleFreqAll "

Usage:

    x.alleleFreqAll()

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleFreqVec "

Usage:

    x.alleleFreqVec(loc)

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleFreq "

Usage:

    x.alleleFreq(allele, loc)

"; 

%feature("docstring") simuPOP::statAlleleFreq::numOfAlleles "

Usage:

    x.numOfAlleles()

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleles "

Usage:

    x.alleles(loc)

"; 

%feature("docstring") simuPOP::statAlleleFreq::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statAssociation;

%feature("docstring") simuPOP::statAssociation::statAssociation "

Usage:

    statAssociation(alleleFreq, haploFreq, Association=[], param={})

"; 

%feature("docstring") simuPOP::statAssociation::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statExpHetero;

%feature("docstring") simuPOP::statExpHetero::statExpHetero "

Usage:

    statExpHetero(alleleFreq, expHetero=[], param={})

"; 

%feature("docstring") simuPOP::statExpHetero::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statFst;

%feature("docstring") simuPOP::statFst::statFst "

Usage:

    statFst(alleleFreq, heteroFreq, Fst=[], param={})

"; 

%feature("docstring") simuPOP::statFst::Fst "

Usage:

    x.Fst()

"; 

%feature("docstring") simuPOP::statFst::Fis "

Usage:

    x.Fis()

"; 

%feature("docstring") simuPOP::statFst::Fit "

Usage:

    x.Fit()

"; 

%feature("docstring") simuPOP::statFst::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statGenoFreq;

%feature("docstring") simuPOP::statGenoFreq::statGenoFreq "

Usage:

    statGenoFreq(genoFreq=[], param={})

"; 

%feature("docstring") simuPOP::statGenoFreq::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statHaploFreq;

%feature("docstring") simuPOP::statHaploFreq::statHaploFreq "

Usage:

    statHaploFreq(haploFreq=[])

"; 

%feature("docstring") simuPOP::statHaploFreq::~statHaploFreq "

Usage:

    x.~statHaploFreq()

"; 

%feature("docstring") simuPOP::statHaploFreq::addHaplotype "

Usage:

    x.addHaplotype(haplo, post=False)

"; 

%feature("docstring") simuPOP::statHaploFreq::numOfHaplotypes "

Usage:

    x.numOfHaplotypes(haplo)

"; 

%feature("docstring") simuPOP::statHaploFreq::haploNum "

Usage:

    x.haploNum(haplo)

"; 

%feature("docstring") simuPOP::statHaploFreq::haploFreq "

Usage:

    x.haploFreq(haplo)

"; 

%feature("docstring") simuPOP::statHaploFreq::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statHeteroFreq;

%feature("docstring") simuPOP::statHeteroFreq::statHeteroFreq "

Usage:

    statHeteroFreq(heteroFreq=[], homoFreq=[])

"; 

%feature("docstring") simuPOP::statHeteroFreq::addLocus "

Usage:

    x.addLocus(locus, post=False)

"; 

%feature("docstring") simuPOP::statHeteroFreq::heteroNum "

Usage:

    x.heteroNum(allele, loc)

"; 

%feature("docstring") simuPOP::statHeteroFreq::heteroFreq "

Usage:

    x.heteroFreq(allele, loc)

"; 

%feature("docstring") simuPOP::statHeteroFreq::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statLD;

%feature("docstring") simuPOP::statLD::statLD "

Usage:

    statLD(alleleFreq, haploFreq, LD=[], LD_param={})

"; 

%feature("docstring") simuPOP::statLD::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statNumOfAffected;

%feature("docstring") simuPOP::statNumOfAffected::statNumOfAffected "

Usage:

    statNumOfAffected(numOfAffected=False, param={})

"; 

%feature("docstring") simuPOP::statNumOfAffected::~statNumOfAffected "

Usage:

    x.~statNumOfAffected()

"; 

%feature("docstring") simuPOP::statNumOfAffected::activate "

Usage:

    x.activate(yes=True)

"; 

%feature("docstring") simuPOP::statNumOfAffected::numOfAffected "

Usage:

    x.numOfAffected()

"; 

%feature("docstring") simuPOP::statNumOfAffected::numOfUnaffected "

Usage:

    x.numOfUnaffected()

"; 

%feature("docstring") simuPOP::statNumOfAffected::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statNumOfAlleles;

%feature("docstring") simuPOP::statNumOfAlleles::statNumOfAlleles "

Usage:

    statNumOfAlleles(calc, atLoci=[], param={})

"; 

%feature("docstring") simuPOP::statNumOfAlleles::~statNumOfAlleles "

Usage:

    x.~statNumOfAlleles()

"; 

%feature("docstring") simuPOP::statNumOfAlleles::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statNumOfMale;

%feature("docstring") simuPOP::statNumOfMale::statNumOfMale "

Usage:

    statNumOfMale(numOfMale=False, param={})

"; 

%feature("docstring") simuPOP::statNumOfMale::activate "

Usage:

    x.activate(yes=True)

"; 

%feature("docstring") simuPOP::statNumOfMale::numOfMale "

Usage:

    x.numOfMale()

"; 

%feature("docstring") simuPOP::statNumOfMale::numOfFemale "

Usage:

    x.numOfFemale()

"; 

%feature("docstring") simuPOP::statNumOfMale::apply "

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::stator "

Description:

    base class of all the statistics calculator

Details:

    Operator stator calculates various basic statistics for the
    population and set variables in the local namespace. Other
    operators or functions can refer to the results from the namespace
    after stat is applied.

"; 

%feature("docstring") simuPOP::stator::stator "

Description:

    create a stator

Usage:

    stator(output=\"\", outputExpr=\"\", stage=PostMating, begin=0,
      end=-1, step=1, at=[], rep=[], subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::stator::~stator "

Description:

    destructor

Usage:

    x.~stator()

"; 

%feature("docstring") simuPOP::stator::clone "

Description:

    deep copy of a stator

Usage:

    x.clone()

"; 

%ignore simuPOP::statPopSize;

%feature("docstring") simuPOP::statPopSize::statPopSize "

Usage:

    statPopSize(popSize=False)

"; 

%feature("docstring") simuPOP::statPopSize::activate "

Usage:

    x.activate()

"; 

%feature("docstring") simuPOP::statPopSize::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statRelatedness;

%feature("docstring") simuPOP::statRelatedness::statRelatedness "

Description:

    calculate relatedness measures between elements in groups

Usage:

    statRelatedness(alleleFreq, groups=[], useSubPop=False, loci=[],
      method=[], minScored=10, param={})

Arguments:

    groups:         can be [ [1,2,3],[4,5,6],[7,8,9]] as three groups
                    of individuals; or [ 1 3 4] as three
                    subpopulations. To specify between individual
                    relatedness, use [[1],[2],[3]] (the first form).
                    If this parameter is ignored, this operator
                    calculate relatedness between all subpopulations.
    method:         can be REL_Queller, REL_Lynch, REL_IR, REL_D2 or
                    REL_Rel. Please refer to the manual for details.

"; 

%feature("docstring") simuPOP::statRelatedness::relQueller "

Usage:

    x.relQueller(ind1, ind2)

"; 

%feature("docstring") simuPOP::statRelatedness::relLynch "

Usage:

    x.relLynch(ind1, ind2)

"; 

%feature("docstring") simuPOP::statRelatedness::relIR "

Usage:

    x.relIR(ind1, locus)

"; 

%feature("docstring") simuPOP::statRelatedness::relD2 "

Usage:

    x.relD2(ind1, locus)

"; 

%feature("docstring") simuPOP::statRelatedness::relRel "

Usage:

    x.relRel(ind1, ind2, locus)

"; 

%feature("docstring") simuPOP::statRelatedness::groupRelatedness "

Description:

    for group i and locus j otherwise

Usage:

    x.groupRelatedness(pop, i, j, method)

"; 

%feature("docstring") simuPOP::statRelatedness::apply "

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::StopEvolution "

Description:

    all replicates.

"; 

%feature("docstring") simuPOP::StopEvolution::StopEvolution "

Usage:

    StopEvolution(msg)

"; 

%feature("docstring") simuPOP::StopIteration "

Description:

    exception, thrown if out of memory

"; 

%feature("docstring") simuPOP::StopIteration::StopIteration "

Usage:

    StopIteration(msg)

"; 

%ignore simuPOP::StreamElem;

%ignore simuPOP::StreamElem::StreamElem(const string &name, bool readable, bool realAppend, bool useString);

%ignore simuPOP::StreamElem::StreamElem(const StreamElem &rhs);

%feature("docstring") simuPOP::StreamElem::~StreamElem "

Description:

    destructor. Close stream and delete m_stream pointer.

Usage:

    x.~StreamElem()

"; 

%ignore simuPOP::StreamElem::makeReadable();

%ignore simuPOP::StreamElem::makeAppend(bool append);

%ignore simuPOP::StreamElem::stream();

%ignore simuPOP::StreamElem::type();

%ignore simuPOP::StreamElem::info();

%ignore simuPOP::StreamElem::append();

%ignore simuPOP::StreamProvider;

%ignore simuPOP::StreamProvider::StreamProvider(const string &output, const string &outputExpr);

%feature("docstring") simuPOP::StreamProvider::~StreamProvider "

Usage:

    x.~StreamProvider()

"; 

%ignore simuPOP::StreamProvider::setOutput(const string &output, const string &outputExpr);

%ignore simuPOP::StreamProvider::noOutput();

%ignore simuPOP::StreamProvider::getOstream(PyObject *dict=NULL, bool readable=false);

%ignore simuPOP::StreamProvider::closeOstream();

%feature("docstring") simuPOP::subPopList "

Details:

    A class to specify (virtual) subpopulation list. Using a dedicated
    class allows users to specify a single subpopulation, or a list of
    (virutal) subpoulations easily.

"; 

%feature("docstring") simuPOP::subPopList::subPopList "

Usage:

    subPopList(subPops=[])

"; 

%feature("docstring") simuPOP::subPopList::empty "

Usage:

    x.empty()

"; 

%feature("docstring") simuPOP::subPopList::size "

Usage:

    x.size()

"; 

%feature("docstring") simuPOP::subPopList::push_back "

Usage:

    x.push_back(subPop)

"; 

%feature("docstring") simuPOP::subPopList::begin "

Usage:

    x.begin()

"; 

%feature("docstring") simuPOP::subPopList::end "

Usage:

    x.end()

"; 

%feature("docstring") simuPOP::SystemError "

Description:

    exception, thrown if system error occurs

"; 

%feature("docstring") simuPOP::SystemError::SystemError "

Usage:

    SystemError(msg)

"; 

%feature("docstring") simuPOP::tagger "

Description:

    base class of tagging individuals

Details:

    This is a during-mating operator that tags individuals with
    various information. Potential usages are:
    *  recording the parental information to track pedigree;
    *  tagging an individual/allele and monitoring its spread in the
    population etc.

"; 

%feature("docstring") simuPOP::tagger::tagger "

Description:

    create a tagger, default to be always active but no output

Usage:

    tagger(output=\"\", outputExpr=\"\", stage=DuringMating, begin=0,
      end=-1, step=1, at=[], rep=[], subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::tagger::~tagger "

Description:

    destructor

Usage:

    x.~tagger()

"; 

%feature("docstring") simuPOP::tagger::clone "

Description:

    deep copy of a \\ tagger

Usage:

    x.clone()

"; 

%ignore simuPOP::tagger::apply(population &pop);

%feature("docstring") simuPOP::terminateIf "

Details:

    This operator evaluates an expression in a population's local
    namespace and terminate the evolution of this population, or the
    whole simulator, if the return value of this expression is True.
    Termination caused by an operator will stop the execution of all
    operators after it. Because a life-cycle is considered to be
    complete if mating is complete, the evolved generations (return
    value from simulator::evolve) of a terminated replicate is
    determined by when the last evolution cycle is terminated.

"; 

%feature("docstring") simuPOP::terminateIf::terminateIf "

Usage:

    terminateIf(condition=\"\", stopAll=False, message=\"\", output=\"\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=[], subPop=[], infoFields=[])

Details:

    Create a terminator with an expression condition, which will be
    evalulated in a population's local namespace when the operator is
    applied to this population. If the return value of condition is
    True, the evolution of the population will be terminated. If
    stopAll is set to True, the evolution of all replicates of the
    simulator will be terminated. If this operator is allowed to write
    to an output or outputExpr (both default to \"\"), the generation
    number, preceeded with an optional message will be written to it.

"; 

%feature("docstring") simuPOP::terminateIf::clone "

Description:

    deep copy of a terminateIf terminator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::terminateIf::__repr__ "

Description:

    used by Python print function to print out the general information
    of the terminateIf terminator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::terminateIf::apply "

Usage:

    x.apply(pop)

Details:

    Apply an operator to population pop directly, without checking its
    applicability.

"; 

%feature("docstring") simuPOP::terminateIf::~terminateIf "

Usage:

    x.~terminateIf()

"; 

%feature("docstring") simuPOP::ticToc "

Function form:

    TicToc

Description:

    timer operator

Details:

    This operator, when called, output the difference between current
    and the last called clock time. This can be used to estimate
    execution time of each generation. Similar information can also be
    obtained from turnOnDebug(DBG_PROFILE), but this operator has the
    advantage of measuring the duration between several generations by
    setting step parameter.

"; 

%feature("docstring") simuPOP::ticToc::ticToc "

Description:

    create a timer

Usage:

    ticToc(output=\">\", outputExpr=\"\", stage=PreMating, begin=0,
      end=-1, step=1, at=[], rep=[], subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::ticToc::~ticToc "

Description:

    destructor

Usage:

    x.~ticToc()

"; 

%feature("docstring") simuPOP::ticToc::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::ticToc::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::ticToc::__repr__ "

Description:

    used by Python print function to print out the general information
    of the ticToc operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::turnOffDebug "

Function form:

    TurnOffDebug

Description:

    set debug off

Details:

    Turn off debug.

"; 

%feature("docstring") simuPOP::turnOffDebug::turnOffDebug "

Description:

    create a turnOffDebug operator

Usage:

    turnOffDebug(code, stage=PreMating, begin=0, end=-1, step=1,
      at=[], rep=[], subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::turnOffDebug::~turnOffDebug "

Description:

    destructor

Usage:

    x.~turnOffDebug()

"; 

%feature("docstring") simuPOP::turnOffDebug::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::turnOffDebug::apply "

Description:

    apply the turnOffDebug operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::turnOffDebug::__repr__ "

Description:

    used by Python print function to print out the general information
    of the turnOffDebug operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::turnOnDebug "

Function form:

    TurnOnDebug

Description:

    set debug on

Details:

    Turn on debug. There are several ways to turn on debug information
    for non-optimized modules, namely
    *  set environment variable SIMUDEBUG.
    *  use simuOpt.setOptions(debug) function.
    *  use TurnOnDebug or TurnOnDebugByName function.
    *  use this turnOnDebug operator The advantage of using this
    operator is that you can turn on debug at given generations.

"; 

%feature("docstring") simuPOP::turnOnDebug::turnOnDebug "

Description:

    create a turnOnDebug operator

Usage:

    turnOnDebug(code, stage=PreMating, begin=0, end=-1, step=1,
      at=[], rep=[], subPop=[], infoFields=[])

"; 

%feature("docstring") simuPOP::turnOnDebug::~turnOnDebug "

Description:

    destructor

Usage:

    x.~turnOnDebug()

"; 

%feature("docstring") simuPOP::turnOnDebug::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::turnOnDebug::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::turnOnDebug::__repr__ "

Description:

    used by Python print function to print out the general information
    of the turnOnDebug operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::TypeError "

Description:

    exception, thrown if type mismatch

"; 

%feature("docstring") simuPOP::TypeError::TypeError "

Usage:

    TypeError(msg)

"; 

%feature("docstring") simuPOP::ValueError "

Description:

    exception, thrown if value of range etc

"; 

%feature("docstring") simuPOP::ValueError::ValueError "

Usage:

    ValueError(msg)

"; 

%feature("docstring") simuPOP::vspID "

Details:

    A class to specify virtual subpopulation, which is composed of a
    subPopulation ID and a virtual subpopulation ID.

"; 

%feature("docstring") simuPOP::vspID::vspID "

Usage:

    vspID(subPop)

"; 

%feature("docstring") simuPOP::vspID::subPop "

Usage:

    x.subPop()

"; 

%feature("docstring") simuPOP::vspID::virtualSubPop "

Usage:

    x.virtualSubPop()

"; 

%feature("docstring") simuPOP::vspID::isVirtual "

Usage:

    x.isVirtual()

"; 

%feature("docstring") simuPOP::vspSplitter "

Details:

    This class is the base class of all virtual subpopulation (VSP)
    splitters, which provide ways to define groups of individuals in a
    subpopulation who share certain properties. A splitter defines a
    fixed number of named VSPs. They do not have to add up to the
    whole subpopulation, nor do they have to be distinct. After a
    splitter is assigned to a population, many functions and operators
    can be applied to individuals within specified VSPs. Only one VSP
    splitter can be assigned to a population, which defined VSPs for
    all its subpopulations. It different splitters are needed for
    different subpopulations, a combinedSplitter should be.

"; 

%feature("docstring") simuPOP::vspSplitter::vspSplitter "

Usage:

    vspSplitter()

Details:

    This is a virtual class that cannot be instantiated.

"; 

%feature("docstring") simuPOP::vspSplitter::clone "

Usage:

    x.clone()

Details:

    All VSP splitter defines a clone() function to create an identical
    copy of itself.

"; 

%feature("docstring") simuPOP::vspSplitter::~vspSplitter "

Usage:

    x.~vspSplitter()

"; 

%ignore simuPOP::vspSplitter::activatedSubPop() const;

%ignore simuPOP::vspSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::vspSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs defined by this splitter.

"; 

%ignore simuPOP::vspSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::vspSplitter::deactivate(population &pop, SubPopID subPop);

%feature("docstring") simuPOP::vspSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of VSP vsp (an index between 0 and
    numVirtualSubPop()).

"; 

%feature("docstring") simuPOP::Weightedsampler "

"; 

%feature("docstring") simuPOP::Weightedsampler::Weightedsampler "

Usage:

    Weightedsampler(rng, weight=[], fast=True)

"; 

%feature("docstring") simuPOP::Weightedsampler::~Weightedsampler "

Usage:

    x.~Weightedsampler()

"; 

%feature("docstring") simuPOP::Weightedsampler::set "

Usage:

    x.set(weight)

"; 

%feature("docstring") simuPOP::Weightedsampler::biSearch "

Usage:

    x.biSearch(a)

"; 

%feature("docstring") simuPOP::Weightedsampler::get "

Usage:

    x.get()

"; 

%feature("docstring") simuPOP::Weightedsampler::q "

Usage:

    x.q()

"; 

%feature("docstring") simuPOP::Weightedsampler::a "

Usage:

    x.a()

"; 

%ignore simuPOP::countAlleles(population &pop, int subpop, const vectori &loci, const vectori &alleles, vectorlu &numAllele);

%ignore simuPOP::getExpectedAlleles(population &pop, vectorf &expFreq, const vectori &loci, const vectori &alleles, vectoru &expAlleles);

%feature("docstring") simuPOP::FreqTrajectoryStoch "

Usage:

    FreqTrajectoryStoch(curGen=0, freq=0, N=0, NtFunc=None,
      fitness=[], fitnessFunc=None, minMutAge=0, maxMutAge=100000,
      ploidy=2, restartIfFail=False, maxAttempts=1000,
      allowFixation=False)

"; 

%ignore simuPOP::MarginalFitness(unsigned nLoci, const vectorf &fitness, const vectorf &freq);

%feature("docstring") simuPOP::FreqTrajectoryMultiStoch "

Usage:

    FreqTrajectoryMultiStoch(curGen=0, freq=[], N=0, NtFunc=None,
      fitness=[], fitnessFunc=None, minMutAge=0, maxMutAge=100000,
      ploidy=2, restartIfFail=False, maxAttempts=1000)

"; 

%feature("docstring") simuPOP::ForwardFreqTrajectory "

Usage:

    ForwardFreqTrajectory(curGen=0, endGen=0, curFreq=[], freq=[],
      N=[], NtFunc=None, fitness=[], fitnessFunc=None, migrRate=0,
      ploidy=2, maxAttempts=1000)

"; 

%feature("docstring") simuPOP::FreqTrajectorySelSim "

Usage:

    FreqTrajectorySelSim(sel, Ne, freq, dom_h, selection)

"; 

%feature("docstring") simuPOP::FreqTrajectoryForward "

Usage:

    FreqTrajectoryForward(lowbound, highbound, disAge, grate, N0,
      seleCo)

"; 

%feature("docstring") simuPOP::ApplyDuringMatingOperator "Obsolete or undocumented function."

%feature("docstring") simuPOP::LoadPopulation "

Usage:

    LoadPopulation(file)

Details:

    load a population from a file.

"; 

%feature("docstring") simuPOP::LoadSimulator "

Description:

    load a simulator from a file with the specified mating scheme. The
    file format is by default determined by file extension
    (format=\"auto\"). Otherwise, format can be one of txt, bin, or xml.

Usage:

    LoadSimulator(file, matingScheme)

"; 

%ignore simuPOP::haploKey(const vectori &seq);

%feature("docstring") simuPOP::TurnOnDebug "

Description:

    set debug codes. Default to turn on all debug codes. Only
    available in non-optimized modules.

Usage:

    TurnOnDebug(code=DBG_ALL)

"; 

%feature("docstring") simuPOP::TurnOnDebug "

Usage:

    TurnOnDebug(code)

"; 

%feature("docstring") simuPOP::TurnOffDebug "

Description:

    turn off debug information. Default to turn off all debug codes.
    Only available in non-optimized modules.

Usage:

    TurnOffDebug(code=DBG_ALL)

"; 

%ignore simuPOP::debug(DBG_CODE code);

%feature("docstring") simuPOP::ListDebugCode "

Description:

    list all debug codes

Usage:

    ListDebugCode()

"; 

%ignore simuPOP::dbgString(DBG_CODE code);

%ignore simuPOP::simuPOP_kbhit();

%ignore simuPOP::simuPOP_getch();

%ignore simuPOP::PyObj_As_Bool(PyObject *obj, bool &val);

%ignore simuPOP::PyObj_As_Int(PyObject *obj, int &val);

%ignore simuPOP::PyObj_As_Double(PyObject *obj, double &val);

%ignore simuPOP::PyObj_As_String(PyObject *obj, string &val);

%ignore simuPOP::PyObj_As_Array(PyObject *obj, vectorf &val);

%ignore simuPOP::PyObj_As_IntArray(PyObject *obj, vectori &val);

%ignore simuPOP::PyObj_As_Matrix(PyObject *obj, matrix &val);

%ignore simuPOP::PyObj_As_StrDict(PyObject *obj, strDict &val);

%ignore simuPOP::PyObj_As_IntDict(PyObject *obj, intDict &val);

%ignore simuPOP::PyObj_Is_IntNumArray(PyObject *obj);

%ignore simuPOP::PyObj_Is_DoubleNumArray(PyObject *obj);

%ignore simuPOP::PyObj_Is_AlleleNumArray(PyObject *obj);

%ignore simuPOP::Double_Vec_As_NumArray(vectorf::iterator begin, vectorf::iterator end);

%ignore simuPOP::Allele_Vec_As_NumArray(GenoIterator begin, GenoIterator end);

%ignore simuPOP::NumArray_Size(PyObject *obj);

%ignore simuPOP::NumArray_Data(PyObject *obj);

%ignore simuPOP::mainVars();

%ignore simuPOP::moduleVars();

%ignore simuPOP::pyPopObj(void *p);

%ignore simuPOP::pyIndObj(void *p);

%ignore simuPOP::ostreamManager();

%feature("docstring") simuPOP::rng "

Description:

    return the currently used random number generator

Usage:

    rng()

"; 

%feature("docstring") simuPOP::SetRNG "

Description:

    set random number generator. If seed=0 (default), a random seed
    will be given. If rng=\"\", seed will be set to the current random
    number generator.

Usage:

    SetRNG(rng=\"\", seed=0)

"; 

%feature("docstring") simuPOP::AvailableRNGs "

Description:

    list the names of all available random number generators

Usage:

    AvailableRNGs()

"; 

%feature("docstring") simuPOP::simuRev "

Description:

    return the revision number of this simuPOP module. Can be used to
    test if a feature is available.

Usage:

    simuRev()

"; 

%feature("docstring") simuPOP::simuVer "

Description:

    return the version of this simuPOP module

Usage:

    simuVer()

"; 

%feature("docstring") simuPOP::ModuleCompiler "

Description:

    return the compiler used to compile this simuPOP module

Usage:

    ModuleCompiler()

"; 

%feature("docstring") simuPOP::ModuleDate "

Description:

    return the date when this simuPOP module is compiled

Usage:

    ModuleDate()

"; 

%feature("docstring") simuPOP::ModulePyVersion "

Description:

    return the Python version this simuPOP module is compiled for

Usage:

    ModulePyVersion()

"; 

%feature("docstring") simuPOP::ModulePlatForm "

Description:

    return the platform on which this simuPOP module is compiled

Usage:

    ModulePlatForm()

"; 

%ignore simuPOP::initialize();

%feature("docstring") simuPOP::Optimized "

Description:

    return True if this simuPOP module is optimized

Usage:

    Optimized()

"; 

%feature("docstring") simuPOP::Limits "

Description:

    print out system limits

Usage:

    Limits()

"; 

%feature("docstring") simuPOP::AlleleType "

Description:

    return the allele type of the current module. Can be binary,
    short, or long.

Usage:

    AlleleType()

"; 

%feature("docstring") simuPOP::MaxAllele "

Usage:

    MaxAllele()

Details:

    return the maximum allowed allele state of the current simuPOP
    module, which is 1 for binary modules, 255 for short modules and
    65535 for long modules.

"; 

%ignore simuPOP::cnull();

%feature("docstring") simuPOP::setLogOutput "

Description:

    set the standard output (default to standard Python output)

Usage:

    setLogOutput(filename=\"\")

"; 

%ignore simuPOP::isGzipped(const string &filename);

%ignore simuPOP::fileExtension(const string &filename);

%ignore std::pow3(unsigned n);

%feature("docstring") simuPOP::population::dvars "

Usage:

    x.dvars()

Details:

    Return a wrapper of Python dictionary returned by vars() so that
    dictionary keys can be accessed as attributes. For example
    pop.dvars().alleleFreq is equivalent to pop.vars()[\"alleleFreq\"].

"; 

%feature("docstring") simuPOP::simulator::dvars "

Usage:

    x.dvars(rep)

Details:

    Return a wrapper of Python dictionary returned by vars(rep) so
    that dictionary keys can be accessed as attributes. For example
    simu.dvars(1).alleleFreq is equivalent to
    simu.vars(1)[\"alleleFreq\"].

"; 

%feature("docstring") simuPOP::population::dvars "

Usage:

    x.dvars(subPop)

Details:

    Return a wrapper of Python dictionary returned by vars(subPop) so
    that dictionary keys can be accessed as attributes.

"; 

%feature("docstring") simuPOP::simulator::dvars "

Usage:

    x.dvars(rep, subPop)

Details:

    Return a wrapper of Python dictionary returned by vars(rep,
    subPop) so that dictionary keys can be accessed as attributes.

"; 
