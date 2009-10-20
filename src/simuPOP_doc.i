%feature("docstring") simuPOP::affectionSplitter "

Details:

    This class defines two VSPs according individual affection status.
    The first VSP consists of unaffected invidiauls and the second VSP
    consists of affected ones.

"; 

%feature("docstring") simuPOP::affectionSplitter::affectionSplitter "

Usage:

    affectionSplitter(names=[])

Details:

    Create a splitter that defined two VSPs by affection status.These
    VSPs are named Unaffected and Affected unless a new set of names
    are specified by parameter names.

"; 

%feature("docstring") simuPOP::affectionSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::affectionSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::affectionSplitter::numVirtualSubPop "

Description:

    Return 2.

Usage:

    x.numVirtualSubPop()

"; 

%ignore simuPOP::affectionSplitter::contains(const population &pop, ULONG ind, vspID vsp);

%ignore simuPOP::affectionSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::affectionSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::affectionSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return \"Unaffected\" if vsp=0 and \"Affected\" if vsp=1, unless a new
    set of names are specified.

"; 

%feature("docstring") simuPOP::alphaParentsChooser "

Details:

    This parent chooser mimicks some animal populations where only
    certain individuals (usually males) can mate. Alpha individuals
    can be chosen either randomly (with natural selection) or
    according to an information field. After the alpha individuals are
    selected, the parent chooser works identical to a random mating
    scheme, except that one of the parents are chosen from these alpha
    individuals.

"; 

%feature("docstring") simuPOP::alphaParentsChooser::alphaParentsChooser "

Usage:

    alphaParentsChooser(alphaSex=Male, alphaNum=0, alphaField=\"\",
      selectionField=\"fitness\")

Details:

    Create a parent chooser that chooses father (if alphaSex is Male)
    or mother (if alphaSex is Female) from a selected group of alpha
    individuals. If alphaNum is given, alpha individuals are chosen
    randomly or according to individual fitness if natural selection
    is enabled. If alphaField is given, individuals with non-zero
    values at this information field are considered as alpha
    individuals. After alpha individuals are selected, alphaSex parent
    will be chosen from the alpha individuals randomly or according to
    individual fitness. The other parents are chosen randomly.

"; 

%feature("docstring") simuPOP::alphaParentsChooser::clone "

Description:

    Deep copy of an alpha parents chooser.

Usage:

    x.clone()

"; 

%ignore simuPOP::alphaParentsChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::alphaParentsChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::baseOperator "

Details:

    Operators are objects that act on populations. They can be applied
    to populations directly using their function forms, but they are
    usually managed and applied by a simulator. In the latter case,
    operators are passed to the evolve function of a simulator, and
    are applied repeatedly during the evolution of the simulator.  The
    baseOperator class is the base class for all operators. It defines
    a common user interface that specifies at which generations, at
    which stage of a life cycle, to which populations and
    subpopulations an operator is applied. These are achieved by a
    common set of parameters such as begin, end, step, at, stage for
    all operators. Note that a specific operator does not have to
    honor all these parameters. For example, a recombinator can only
    be applied during mating so it ignores the stage parameter.  An
    operator can be applied to all or part of the generations during
    the evolution of a simulator. At the beginning of an evolution, a
    simulator is usually at the beginning of generation 0. If it
    evolves 10 generations, it evolves generations 0, 1, ,,,., and 9
    (10 generations) and stops at the begging of generation 10. A
    negative generation number a has generation number 10 + a, with -1
    referring to the last evolved generation 9. Note that the starting
    generation number of a simulator can be changed by its setGen()
    member function.  Output from an operator is usually directed to
    the standard output (sys.stdout). This can be configured using a
    output specification string, which can be '' for no output, '>'
    standard terminal output (default), a filename prefixed by one or
    more '>' characters or a Python expression indicated by a leading
    exclamation mark ('!expr'). In the case of '>filename' (or
    equivalently 'filename'), the output from an operator is written
    to this file. However, if two operators write to the same file
    filename, or if an operator writes to this file more than once,
    only the last write operation will succeed. In the case of
    '>>filename', file filename will be opened at the beginning of the
    evolution and closed at the end. Outputs from multiple operators
    are appended. >>>filename works similar to >>filename but
    filename, if it already exists at the beginning of an evolutionary
    process, will not be cleared. If the output specification is
    prefixed by an exclamation mark, the string after the mark is
    considered as a Python expression. When an operator is applied to
    a population, this expression will be evaluated within the
    population's local namespace to obtain a population specific
    output specification. As an advanced feature, a Python function
    can be assigned to this parameter. Output strings will be sent to
    this function for processing.

"; 

%feature("docstring") simuPOP::baseOperator::baseOperator "

Usage:

    baseOperator(output, begin, end, step, at, reps, subPops,
      infoFields)

Details:

    The following parameters can be specified by all operators.
    However, an operator can ignore some parameters and the exact
    meaning of a parameter can vary.

Arguments:

    output:         A string that specifies how output from an
                    operator is written, which can be '' (no output),
                    '>' (standard output), 'filename' prefixed by one
                    or more '>', or an Python expression prefixed by
                    an exclamation mark ('!expr'). Alternatively, a
                    Python function can be given to handle outputs.
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
                    parameter is specified. A single generation number
                    is also acceptable.
    reps:           A list of applicable replicates. A common default
                    value AllAvail is interpreted as all replicates in
                    a simulator. Negative indexes such as -1 (last
                    replicate) is acceptable. rep=idx can be used as a
                    shortcut for rep=[idx].
    subPops:        A list of applicable (virtual) subpopulations,
                    such as subPops=[sp1, sp2, (sp2, vsp1)].
                    subPops=[sp1] can be simplied as subPops=sp1.
                    Negative indexes are not supported. A common
                    default value (AllAvail) of this parameter
                    reprents all subpopulations of the population
                    being aplied. Suport for this parameter vary from
                    operator to operator and some operators do not
                    support virtual subpopulations at all. Please
                    refer to the reference manual of individual
                    operators for their support for this parameter.
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

%ignore simuPOP::baseOperator::isActive(UINT rep, long gen, long end, const vector< bool > &activeRep, bool repOnly=false);

%ignore simuPOP::baseOperator::isActive(UINT rep, long gen);

%ignore simuPOP::baseOperator::isCompatible(const population &pop);

%ignore simuPOP::baseOperator::haploidOnly();

%ignore simuPOP::baseOperator::diploidOnly();

%ignore simuPOP::baseOperator::setHaploidOnly();

%ignore simuPOP::baseOperator::setDiploidOnly();

%ignore simuPOP::baseOperator::infoSize();

%ignore simuPOP::baseOperator::infoField(UINT idx);

%ignore simuPOP::baseOperator::infoFields();

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

%feature("docstring") simuPOP::baseOperator::description "Obsolete or undocumented function."

%ignore simuPOP::baseOperator::noOutput();

%feature("docstring") simuPOP::baseOperator::initializeIfNeeded "

Usage:

    x.initializeIfNeeded(pop)

"; 

%feature("docstring") simuPOP::baseOperator::initialize "Obsolete or undocumented function."

%ignore simuPOP::baseOperator::applicableSubPops() const;

%feature("docstring") simuPOP::basePenetrance "

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
    accordingly. Penetrance can also be used pre- or post mating. In
    these cases, the affected status will be set to all individuals
    according to their penetrance values.
    Penetrance values are usually not saved. If you would like to know
    the penetrance value, you need to
    *   use addInfoField('penetrance') to the population to analyze.
    (Or use infoFields parameter of the population constructor), and
    *   use e.g., mlPenetrance(...., infoFields=['penetrance']) to add
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

%feature("docstring") simuPOP::basePenetrance::basePenetrance "

Description:

    create a penetrance operator

Usage:

    basePenetrance(ancestralGen=-1, begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=[])

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

%feature("docstring") simuPOP::basePenetrance::~basePenetrance "

Description:

    destructor

Usage:

    x.~basePenetrance()

"; 

%feature("docstring") simuPOP::basePenetrance::clone "

Description:

    deep copy of a penetrance operator

Usage:

    x.clone()

"; 

%ignore simuPOP::basePenetrance::penet(individual *);

%feature("docstring") simuPOP::basePenetrance::apply "

Description:

    set penetrance to all individuals and record penetrance if
    requested

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::basePenetrance::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::basePenetrance::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::BernulliTrials "

Details:

    this class encapsulate behavior of a sequence of Bernulli trial.
    the main idea is that when doing a sequence of Bernulli trials of
    the same probability, we can use much quicker algorithms instead
    of doing n Bernulli trials  For example, when N=10000, p=0.001.
    The usual way to do N Bin(p) trials is to do N randUnif(0,1)<p
    comparison.  using the new method, we can use geometric
    distrubution to find the next true event.  Also, for the cases of
    p=0.5, random bits are generated.  This class maintain a two
    dimensional table: a vector of probabilities cross expected number
    of trials  p1 p2 p3 p4 p5 trial 1 trial 2 ... trial N  We expect
    that N is big (usually populaiton size) and p_i are small  using
    fast BernulliTrial method for fix p, we can fill up this table
    very quickly column by column  This class will provide easy access
    to row (each trial) or column (called each prob) of this table.
    if this table is accessed row by row (each trial), a internal
    index is used.  if index exceeds N, trials will be generated all
    again. if trial will be called, e.g., N+2 times all the time, this
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
    genotype are copied. This genotype transmitter does not copy
    genotype on customized chromosomes.

"; 

%feature("docstring") simuPOP::cloneGenoTransmitter::cloneGenoTransmitter "

Usage:

    cloneGenoTransmitter(output=\"\", begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=AllAvail)

Details:

    Create a clone genotype transmitter (a during-mating operator)
    that copies genotypes from parents to offspring. If two parents
    are specified, genotypes are copied maternally. After genotype
    transmission, offspring sex is copied from parental sex even if
    sex has been determined by an offspring generator. All or
    specified information fields (parameter infoFields, default to
    AllAvail) will also be copied from parent to offspring. Parameters
    subPops is ignored.

"; 

%feature("docstring") simuPOP::cloneGenoTransmitter::clone "

Description:

    Deep copy of a clone genotype transmitter.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::cloneGenoTransmitter::description "Obsolete or undocumented function."

%ignore simuPOP::cloneGenoTransmitter::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::CombinedAlleleIterator "

Details:

    This class implements a C++ iterator class that iterate through
    all alleles in a (virtual) (sub)population using 1. an IndIterator
    that will skip invisible individuals and invalid alleles, or 2. a
    gapped iterator that will run faster, in the case that a): no
    virtual subpopulation b): not sex chromosomes c): not haplodiploid

"; 

%feature("docstring") simuPOP::CombinedAlleleIterator::CombinedAlleleIterator "

Usage:

    CombinedAlleleIterator()

"; 

%feature("docstring") simuPOP::CombinedAlleleIterator::valid "

Usage:

    x.valid()

"; 

%feature("docstring") simuPOP::CombinedAlleleIterator::ptr "

Usage:

    x.ptr()

"; 

%feature("docstring") simuPOP::CombinedAlleleIterator::advance "

Usage:

    x.advance(it, p, valid)

"; 

%feature("docstring") simuPOP::combinedSplitter "

Details:

    This splitter takes several splitters and stacks their VSPs
    together. For example, if the first splitter defines 3 VSPs and
    the second splitter defines 2, the two VSPs from the second
    splitter become the fourth (index 3) and the fifth (index 4) VSPs
    of the combined splitter. In addition, a new set of VSPs could be
    defined as the union of one or more of the original VSPs. This
    splitter is usually used to define different types of VSPs to a
    population.

"; 

%feature("docstring") simuPOP::combinedSplitter::combinedSplitter "

Usage:

    combinedSplitter(splitters=[], vspMap=[], names=[])

Details:

    Create a combined splitter using a list of splitters. For example,
    combinedSplitter([sexSplitter(), affectionSplitter()]) defines a
    combined splitter with four VSPs, defined by male (vsp 0), female
    (vsp 1), unaffected (vsp 2) and affected individuals (vsp 3).
    Optionally, a new set of VSPs could be defined by parameter
    vspMap. Each item in this parameter is a list of VSPs that will be
    combined to a single VSP. For example, vspMap=[(0, 2), (1, 3)] in
    the previous example will define two VSPs defined by male or
    unaffected, and female or affected individuals. VSP names are
    usually determined by splitters, but can also be specified using
    parameter names.

"; 

%ignore simuPOP::combinedSplitter::combinedSplitter(const combinedSplitter &rhs);

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

%ignore simuPOP::combinedSplitter::contains(const population &pop, ULONG ind, vspID vsp);

%ignore simuPOP::combinedSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::combinedSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::combinedSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of a VSP vsp, which is the name a VSP defined by
    one of the combined splitters unless a new set of names is
    specified.

"; 

%feature("docstring") simuPOP::contextMutator "

Function form:

    ContextMutate

Details:

    This context-dependent mutator accepts a list of mutators and use
    one of them to mutate an allele depending on the context of the
    mutated allele.

"; 

%feature("docstring") simuPOP::contextMutator::contextMutator "

Usage:

    contextMutator(rates=[], loci=AllAvail, mutators=[],
      contexts=[], mapIn=[], mapOut=[], output=\">\", begin=0, end=-1,
      step=1, at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a mutator that choose one of the specified mutators to
    mutate an allele when a mutation event happens. The mutators are
    choosen according to the context of the mutated allele, which is
    specified as a list of alleles to the left and right of an allele
    ( contexts). For example, contexts=[(0,0), (0,1), (1,1)] indicates
    which mutators should be used to mutate allele X in the context of
    0X0, 0X1, and 1X1. A context can include more than one alleles at
    both left and right sides of a mutated allele but all contexts
    should have the same (even) number of alleles. If an allele does
    not have full context (e.g. when a locus is the first locus on a
    chromosome), unavailable alleles will be marked as -1. There
    should be a mutator for each context but an additional mutator can
    be specified as the default mutator for unmatched contexts. If
    parameters mapIn is specified, both mutated allele and its context
    alleles will be mapped. Most parameters, including loci, mapIn,
    mapOut, rep, and subPops of mutators specified in parameter
    mutators are ignored. This mutator by default applies to all loci
    unless parameter loci is specified. Please refer to classes
    mutator and baseOperator for descriptions of other parameters.

"; 

%feature("docstring") simuPOP::contextMutator::clone "

Description:

    deep copy of a context-dependentMutator

Usage:

    x.clone()

"; 

%ignore simuPOP::contextMutator::initialize(population &pop);

%ignore simuPOP::contextMutator::mutate(AlleleRef allele, UINT locus);

%feature("docstring") simuPOP::contextMutator::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::controlledOffspringGenerator "

Details:

    This offspring generator populates an offspring population and
    controls allele frequencies at specified loci. At each generation,
    expected allele frequencies at these loci are passed from a user
    defined allele frequency trajectory function. The offspring
    population is populated in two steps. At the first step, only
    families with disease alleles are accepted until until the
    expected number of disease alleles are met. At the second step,
    only families with wide type alleles are accepted to populate the
    rest of the offspring generation. This method is described in
    detail in \"Peng et al, (2007) Forward-time simulations of
    populations with complex human diseases, PLoS Genetics\".

"; 

%feature("docstring") simuPOP::controlledOffspringGenerator::controlledOffspringGenerator "

Usage:

    controlledOffspringGenerator(loci, alleles, freqFunc, ops=[],
      numOffspring=1, sexMode=RandomSex)

Details:

    Create an offspring generator that selects offspring so that
    allele frequency at specified loci in the offspring generation
    reaches specified allele frequency. At the beginning of each
    generation, expected allele frequency of alleles at loci is
    returned from a user-defined trajectory function freqFunc. If
    there is no subpopulation, this function should return a list of
    frequencies for each locus. If there are multiple subpopulations,
    freqFunc can return a list of allele frequencies for all
    subpopulations or combined frequencies that ignore population
    structure. In the former case, allele frequencies should be
    arranged by loc0_sp0, loc1_sp0, ... loc0_sp1, loc1_sp1, ..., and
    so on. In the latter case, overall expected number of alleles are
    scattered to each subpopulation in proportion to existing number
    of alleles in each subpopulation, using a multinomial
    distribution.  After the expected alleles are calculated, this
    offspring generator accept and reject families according to their
    genotype at loci until allele frequecies reach their expected
    values. The rest of the offspring generation is then filled with
    families without only wild type alleles at these loci.  This
    offspring generator is derived from class offspringGenerator.
    Please refer to class offspringGenerator for a detailed
    description of parameters ops, numOffspring and sexMode.

"; 

%ignore simuPOP::controlledOffspringGenerator::controlledOffspringGenerator(const controlledOffspringGenerator &rhs);

%ignore simuPOP::controlledOffspringGenerator::initialize(const population &pop, SubPopID subPop, vector< baseOperator * > const &ops);

%ignore simuPOP::controlledOffspringGenerator::generateOffspring(population &pop, individual *dad, individual *mom, RawIndIterator &offBegin, RawIndIterator &offEnd, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::controlledOffspringGenerator::clone "

Description:

    Deep copy of a controlled random mating scheme.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::dumper "

Details:

    This operator dumps the content of a population in a human
    readable format. Because this output format is not structured and
    can not be imported back to simuPOP, this operator is usually used
    to dump a small population to a terminal for demonstration and
    debugging purposes.

"; 

%feature("docstring") simuPOP::dumper::dumper "

Usage:

    dumper(genotype=True, structure=True, ancGen=0, width=1,
      max=100, loci=[], output=\">\", begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a operator that dumps the genotype structure (if structure
    is True) and genotype (if genotype is True) to an output ( default
    to standard terminal output). Because a population can be large,
    this operator will only output the first 100 (parameter max)
    individuals of the present generation (parameter ancGen). All loci
    will be outputed unless parameter loci are used to specify a
    subset of loci. If a list of (virtual) subpopulations are
    specified, this operator will only output individuals in these
    outputs. Please refer to class baseOperator for a detailed
    explanation for common parameters such as output and stage.

"; 

%feature("docstring") simuPOP::dumper::clone "

Description:

    Deep copy of a dumper operator.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::dumper::apply "

Description:

    Apply a dumper operator to population pop.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::dumper::~dumper "

Description:

    destructor.

Usage:

    x.~dumper()

"; 

%feature("docstring") simuPOP::dumper::description "Obsolete or undocumented function."

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

%ignore simuPOP::Expression::setExpr(const string &expr=string());

%ignore simuPOP::Expression::setStmts(const string &stmts=string());

%ignore simuPOP::Expression::evaluate();

%ignore simuPOP::Expression::valueAsBool();

%ignore simuPOP::Expression::valueAsInt();

%ignore simuPOP::Expression::valueAsDouble();

%ignore simuPOP::Expression::valueAsString();

%ignore simuPOP::Expression::valueAsArray();

%feature("docstring") simuPOP::floatList "

"; 

%feature("docstring") simuPOP::floatList::floatList "

Usage:

    floatList(values=[])

"; 

%ignore simuPOP::floatList::elems() const;

%feature("docstring") simuPOP::floatListFunc "

"; 

%feature("docstring") simuPOP::floatListFunc::floatListFunc "

Usage:

    floatListFunc(values=[])

"; 

%ignore simuPOP::floatListFunc::empty() const;

%ignore simuPOP::floatListFunc::size() const;

%ignore simuPOP::floatListFunc::func() const;

%ignore simuPOP::GenoStructure;

%ignore simuPOP::GenoStructure::GenoStructure();

%ignore simuPOP::GenoStructure::GenoStructure(UINT ploidy, const vectoru &loci, const vectoru &chromTypes, bool haplodiploid, const vectorf &lociPos, const vectorstr &chromNames, const matrixstr &alleleNames, const vectorstr &lociNames, const vectorstr &infoFields);

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
    *   Ploidy, namely the number of homologous sets of chromosomes,
    of a population. Haplodiploid population is also supported.
    *   Number of chromosomes and number of loci on each chromosome.
    *   Positions of loci, which determine the relative distance
    between loci on the same chromosome. No unit is assumed so these
    positions can be ordinal (1, 2, 3, ..., the default), in physical
    distance (bp, kb or mb), or in map distance (e.g. centiMorgan)
    depending on applications.
    *   Names of alleles. Although alleles at different loci usually
    have different names, simuPOP uses the same names for alleles
    across loci for simplicity.
    *   Names of loci and chromosomes.
    *   Names of information fields attached to each individual. In
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

%ignore simuPOP::GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru &loci, const vectoru &chromTypes, bool haplodiploid, const vectorf &lociPos, const vectorstr &chromNames, const matrixstr &alleleNames, const vectorstr &lociNames, const vectorstr &infoFields);

%ignore simuPOP::GenoStruTrait::setGenoStructure(const GenoStructure &rhs);

%ignore simuPOP::GenoStruTrait::setGenoStruIdx(size_t idx);

%feature("docstring") simuPOP::GenoStruTrait::lociDist "

Usage:

    x.lociDist(locus1, locus2)

Details:

    Return the distance between loci locus1 and locus2 on the same
    chromosome. A negative value will be returned if locus1 is after
    locus2.

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociLeft "Obsolete or undocumented function."

%feature("docstring") simuPOP::GenoStruTrait::distLeft "Obsolete or undocumented function."

%feature("docstring") simuPOP::GenoStruTrait::lociCovered "Obsolete or undocumented function."

%ignore simuPOP::GenoStruTrait::gsAddChromFromStru(size_t idx) const;

%ignore simuPOP::GenoStruTrait::gsAddLociFromStru(size_t idx, vectoru &index1, vectoru &index2) const;

%ignore simuPOP::GenoStruTrait::gsRemoveLoci(const vectoru &loci, vectoru &kept);

%ignore simuPOP::GenoStruTrait::gsAddChrom(const vectorf &lociPos, const vectorstr &lociNames, const string &chromName, const matrixstr &alleleNames, UINT chromType) const;

%ignore simuPOP::GenoStruTrait::gsSetAlleleNames(const uintList &loci, const matrixstr &alleleNames);

%ignore simuPOP::GenoStruTrait::gsAddLoci(const vectoru &chrom, const vectorf &pos, const vectorstr &lociNames, const matrixstr &alleleNames, vectoru &newIndex) const;

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

    x.locusPos(locus)

Details:

    return the position of locus locus specified by the lociPos
    parameter of the population function. An IndexError will be raised
    if the absolute index locus is greater than or equal to the total
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

    x.alleleName(allele, locus=0)

Details:

    return the name of allele allele at lcous specified by the
    alleleNames parameter of the population function. locus could be
    ignored if alleles at all loci share the same names. If the name
    of an allele is unspecified, its index ('0', '1', '2', etc) is
    returned. An IndexError will be raised if allele is larger than
    the maximum allowed allele state of this module (ModuleMaxAllele).

"; 

%ignore simuPOP::GenoStruTrait::allAlleleNames() const;

%feature("docstring") simuPOP::GenoStruTrait::alleleNames "

Usage:

    x.alleleNames(locus=0)

Details:

    return a list of allele names at given by the alleleNames
    parameter of the population function. locus could be ignored if
    alleles at all loci share the same names. This list does not have
    to cover all possible allele states of a population so
    alleleNames()[allele] might fail (use alleleNames(allele)
    instead).

"; 

%feature("docstring") simuPOP::GenoStruTrait::locusName "

Usage:

    x.locusName(locus)

Details:

    return the name of locus locus specified by the lociNames
    parameter of the population function. An empty string will be
    returned if no name has been given to locus locus.

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociNames "

Usage:

    x.lociNames()

Details:

    return the names of all loci specified by the lociNames parameter
    of the population function. An empty list will be returned if
    lociNames was not specified.

"; 

%feature("docstring") simuPOP::GenoStruTrait::locusByName "

Usage:

    x.locusByName(name)

Details:

    return the index of a locus with name name. Raise a ValueError if
    no locus is found. Note that empty strings are used for loci
    without name but you cannot lookup such loci using this function.

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

%ignore simuPOP::GenoStruTrait::gsAddInfoFields(const vectorstr &fields);

%ignore simuPOP::GenoStruTrait::gsSetInfoFields(const vectorstr &fields);

%ignore simuPOP::GenoStruTrait::swap(GenoStruTrait &rhs);

%ignore simuPOP::GenoStruTrait::incGenoStruRef() const;

%ignore simuPOP::GenoStruTrait::decGenoStruRef() const;

%feature("docstring") simuPOP::genoTransmitter "

Details:

    This during mating operator is the base class of all genotype
    transmitters. It is made available to users because it provides a
    few member functions that can be used by derived transmitters, and
    by customized Python during mating operators.

"; 

%feature("docstring") simuPOP::genoTransmitter::genoTransmitter "

Usage:

    genoTransmitter(output=\"\", begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a base genotype transmitter.

"; 

%feature("docstring") simuPOP::genoTransmitter::clone "

Description:

    Deep copy of a base genotype transmitter.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::genoTransmitter::clearChromosome "

Usage:

    x.clearChromosome(ind, ploidy, chrom)

Details:

    Clear (set alleles to zero) chromosome chrom on the ploidy-th
    homologous set of chromosomes of individual ind. It is equivalent
    to ind.setGenotype([0], ploidy, chrom).

"; 

%feature("docstring") simuPOP::genoTransmitter::copyChromosome "

Usage:

    x.copyChromosome(parent, parPloidy, offspring, ploidy, chrom)

Details:

    Transmit chromosome chrom on the parPloidy set of homologous
    chromosomes from parent to the ploidy set of homologous
    chromosomes of offspring. It is equivalent to
    offspring.setGenotype(parent.genotype(parPloidy, chrom), polidy,
    chrom).

"; 

%feature("docstring") simuPOP::genoTransmitter::copyChromosomes "

Usage:

    x.copyChromosomes(parent, parPloidy, offspring, ploidy)

Details:

    Transmit the parPloidy set of homologous chromosomes from parent
    to the ploidy set of homologous chromosomes of offspring.
    Customized chromosomes are not copied. It is equivalent to
    offspring.setGenotype(parent.genotype(parPloidy), ploidy).

"; 

%feature("docstring") simuPOP::genoTransmitter::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::genoTransmitter::initialize "

Usage:

    x.initialize(pop)

Details:

    Initialize a base genotype operator for a population. This
    function should be called before any other functions are used to
    transmit genotype.

"; 

%ignore simuPOP::genoTransmitter::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::genotypeSplitter "

Details:

    This class defines a VSP splitter that defines VSPs according to
    individual genotype at specified loci.

"; 

%feature("docstring") simuPOP::genotypeSplitter::genotypeSplitter "

Usage:

    genotypeSplitter(loci, alleles, phase=False, names=[])

Details:

    Create a splitter that defines VSPs by individual genotype at loci
    (a locus index or a list of loci indexes). Each list in a list
    allele defines a VSP, which is a list of allowed alleles at these
    loci. If only one VSP is defined, the outer list of the nested
    list can be ignored. If phase if true, the order of alleles in
    each list is significant. If more than one set of alleles are
    given, individuals having either of them is qualified.  For
    example, in a haploid population, loci=1, alleles=[0, 1] defines a
    VSP with individuals having allele 0 or 1 at locus 1, alleles=[[0,
    1], [2]] defines two VSPs with indivdiuals in the second VSP
    having allele 2 at locus 1. If multiple loci are involved, alleles
    at each locus need to be defined. For example, VSP defined by
    loci=[0, 1], alleles=[0, 1, 1, 1] consists of individuals having
    alleles [0, 1] or [1, 1] at loci [0, 1].  In a haploid population,
    loci=1, alleles=[0, 1] defines a VSP with individuals having
    genotype [0, 1] or [1, 0] at locus 1. alleles[[0, 1], [2, 2]]
    defines two VSPs with indivdiuals in the second VSP having
    genotype [2, 2] at locus 1. If phase is set to True, the first VSP
    will only has individuals with genotype [0, 1]. In the multiple
    loci case, alleles should be arranged by haplotypes, for example,
    loci=[0, 1], alleles=[0, 0, 1, 1], phase=True defines a VSP with
    individuals having genotype -0-0-, -1-1- at loci 0 and 1. If
    phase=False (default), genotypes -1-1-, -0-0-, -0-1- and -1-0- are
    all allowed.  A default set of names are given to each VSP unless
    a new set of names is given by parameter names.

"; 

%feature("docstring") simuPOP::genotypeSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::genotypeSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::genotypeSplitter::numVirtualSubPop "

Description:

    number of virtual subpops of subpopulation sp

Usage:

    x.numVirtualSubPop()

"; 

%ignore simuPOP::genotypeSplitter::contains(const population &pop, ULONG ind, vspID vsp);

%ignore simuPOP::genotypeSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::genotypeSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::genotypeSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return name of VSP vsp, which is \"Genotype loc1,loc2:genotype\" as
    defined by parameters loci and alleles. A user provided name will
    be returned if specified.

"; 

%feature("docstring") simuPOP::haplodiploidGenoTransmitter "

Details:

    A genotype transmitter (during-mating operator) for haplodiploid
    populations. The female parent is considered as diploid and the
    male parent is considered as haploid (only the first homologous
    copy is valid). If the offspring is Female, she will get a random
    copy of two homologous chromosomes of her mother, and get the only
    paternal copy from her father. If the offspring is Male, he will
    only get a set of chromosomes from his mother.

"; 

%feature("docstring") simuPOP::haplodiploidGenoTransmitter::haplodiploidGenoTransmitter "

Usage:

    haplodiploidGenoTransmitter(output=\"\", begin=0, end=-1, step=1,
      at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a haplodiploid genotype transmitter (during-mating
    operator) that transmit parental genotypes from parents to
    offspring in a haplodiploid population. Parameters subPops and
    infoFields are ignored.

"; 

%feature("docstring") simuPOP::haplodiploidGenoTransmitter::clone "

Description:

    Deep copy of a haplodiploid transmitter.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::haplodiploidGenoTransmitter::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::haplodiploidGenoTransmitter::initialize "Obsolete or undocumented function."

%ignore simuPOP::haplodiploidGenoTransmitter::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::heteroMating "

Details:

    A heterogeneous mating scheme that applies a list of mating
    schemes to different (virtual) subpopulations.

"; 

%feature("docstring") simuPOP::heteroMating::heteroMating "

Usage:

    heteroMating(matingSchemes, subPopSize=[],
      shuffleOffspring=True)

Details:

    Create a heterogeneous mating scheme that will apply a list of
    homogeneous mating schemes matingSchemes to different (virtual)
    subpopulations. The size of the offspring generation is determined
    by parameter subPopSize, which can be a list of subpopulation
    sizes or a Python function that returns a list of subpopulation
    sizes at each generation. Please refer to homoMating for a
    detailed explanation of this parameter.  Each mating scheme
    defined in matingSchemes can be applied to one or more (virtual)
    subpopulation. If parameter subPop is not specified, a mating
    scheme will be applied to all subpopulations. If a (virtual)
    subpopulation is specified, a mating scheme will be applied to a
    specific (virtual) subpopulation. A special case is when subPop is
    given as (-1, vsp). In this case, the mating scheme will be
    applied to virtual subpopulation vsp in all subpopulations.  If
    multiple mating schemes are applied to the same subpopulation, a
    weight (parameter weight) can be given to each mating scheme to
    determine how many offspring it will produce. The default for all
    mating schemes are 0. In this case, the number of offspring each
    mating scheme produces is proportional to the size of its parental
    (virtual) subpopulation. If all weights are negative, the numbers
    of offspring are determined by the multiplication of the absolute
    values of the weights and their respective parental (virtual)
    subpopulation sizes. If all weights are positive, the number of
    offspring produced by each mating scheme is proportional to these
    weights. Mating schemes with zero weight in this case will produce
    no offspring. If both negative and positive weights are present,
    negative weights are processed before positive ones.  If multiple
    mating schemes are applied to the same subpopulation, offspring
    produced by these mating schemes are shuffled randomly. If this is
    not desired, you can turn off offspring shuffling by setting
    parameter shuffleOffspring to False.

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

    deep copy of a heterogeneous mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::heteroMating::mate(population &pop, population &scratch, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::homoMating "

Details:

    A homogeneous mating scheme that uses a parent chooser to choose
    parents from a prental generation, and an offspring generator to
    generate offspring from chosen parents. It can be either used
    directly, or within a heterogeneous mating scheme. In the latter
    case, it can be applied to a (virtual) subpopulation.

"; 

%feature("docstring") simuPOP::homoMating::homoMating "

Usage:

    homoMating(chooser, generator, subPopSize=[], subPop=[],
      weight=0)

Details:

    Create a homogeneous mating scheme using a parent chooser chooser
    and an offspring generator generator.  If this mating scheme is
    used directly in a simulator, it will be responsible for creating
    an offspring population according to parameter subPopSize. This
    parameter can be a list of subpopulation sizes (or a number if
    there is only one subpopulation) or a Python function. The
    function should take two parameters, a generation number and a
    population object which is the parental population just before
    mating. The return value of this function should be a list of
    subpopulation sizes for the offspring generation. A single number
    can be returned if there is only one subpopulation. If a function
    is passed to this parameter, it will be called at each generation
    to determine the size of the offspring generation. The passed
    parental population is usually used to determine offspring
    population size from parental population size but nothing stops
    you from modifying this parental population to prepare it for
    mating.  If this mating shcme is used within a heterogeneous
    mating scheme. Parameters subPop and weight are used to determine
    which (virtual) subpopulation this mating scheme will be applied
    to, and how many offspring this mating scheme will produce. Please
    refer to mating scheme heteroMating for the use of these two
    parameters.

"; 

%feature("docstring") simuPOP::homoMating::~homoMating "

Description:

    destructor

Usage:

    x.~homoMating()

"; 

%ignore simuPOP::homoMating::homoMating(const homoMating &rhs);

%feature("docstring") simuPOP::homoMating::clone "

Description:

    Deep copy of a homogeneous mating scheme.

Usage:

    x.clone()

"; 

%ignore simuPOP::homoMating::subPop() const;

%ignore simuPOP::homoMating::virtualSubPop() const;

%ignore simuPOP::homoMating::weight() const;

%ignore simuPOP::homoMating::mateSubPop(population &pop, SubPopID subPop, RawIndIterator offBegin, RawIndIterator offEnd, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::idTagger "

Function form:

    TagID

Details:

    An idTagger gives a unique ID for each individual it is applies
    to. These ID can be used to uniquely identify an individual in a
    multi-generational population and be used to reliably reconstruct
    a pedigree.  To ensure uniqueness across populations, a single
    source of ID is used for this operator. Individual IDs are
    assigned consecutively starting from 0. If you would like to reset
    the sequence or start from a different number, you can call the
    reset(startID) function of any idTagger.  An idTagger is usually
    used during-mating to assign ID to each offspring. However, if it
    is applied directly to a population, it will assign unique IDs to
    all individuals in this population. This property is usually used
    in the preOps parameter of function simulator.evolve to assign
    initial ID to a population.

"; 

%feature("docstring") simuPOP::idTagger::idTagger "

Usage:

    idTagger(begin=0, end=-1, step=1, at=[], reps=AllAvail,
      subPops=AllAvail, output=\"\", infoFields=\"ind_id\")

Details:

    Create an idTagger that assign an unique ID for each individual it
    is applied to. The IDs are created sequentially and are stored in
    an information field specified in parameter infoFields (default to
    ind_id). This operator is considered a during-mating operator but
    it can be used to set ID for all individuals of a population when
    it is directly applied to the population. Because the information
    field is supposed to record a unique ID for the whole population,
    and because the IDs are increasingly assigned, this operator will
    raise a RuntimeError if parental IDs are the same, or are larger
    than the ID to be assigned to an offspring.

"; 

%feature("docstring") simuPOP::idTagger::~idTagger "

Usage:

    x.~idTagger()

"; 

%feature("docstring") simuPOP::idTagger::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::idTagger::reset "

Usage:

    x.reset(startID=0)

Details:

    Reset the global individual ID number so that idTaggers will start
    from id (default to 0) again.

"; 

%feature("docstring") simuPOP::idTagger::apply "

Usage:

    x.apply(pop)

Details:

    Set an unique ID to all individuals with zero ID.

"; 

%ignore simuPOP::idTagger::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::idTagger::clone "

Description:

    deep copy of an idTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::ifElse "

Details:

    This operator accepts an expression that will be evaluated when
    this operator is applied. A list of if-operators will be applied
    when the expression returns True. Otherwise a list of else-
    operators will be applied.

"; 

%feature("docstring") simuPOP::ifElse::ifElse "

Usage:

    ifElse(cond, ifOps=[], elseOps=[], output=\">\", begin=0, end=-1,
      step=1, at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a conditional operator that will apply operators ifOps if
    condition cond is met and elseOps otherwise. The replicate and
    generation applicability parameters (begin, end, step, at and rep)
    of the ifOps and elseOps are ignored because their applicability
    is determined by the ifElse operator.

"; 

%feature("docstring") simuPOP::ifElse::~ifElse "

Description:

    destructor

Usage:

    x.~ifElse()

"; 

%feature("docstring") simuPOP::ifElse::clone "Obsolete or undocumented function."

%ignore simuPOP::ifElse::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::ifElse::apply "

Description:

    apply the ifElse operator to population pop.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::ifElse::description "Obsolete or undocumented function."

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
    information fields of an individual.    Genotypes of an individual
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

    x.allele(idx, ploidy=-1, chrom=-1)

Details:

    return the current allele at a locus, using its absolute index
    idx. If a ploidy ploidy and/or a chromosome indexes are given, idx
    is relative to the beginning of specified homologous copy of
    chromosomes (if chrom=-1) or the beginning of the specified
    homologous copy of specified chromosome (if chrom >= 0).

"; 

%feature("docstring") simuPOP::individual::alleleChar "Obsolete or undocumented function."

%feature("docstring") simuPOP::individual::setAllele "

Usage:

    x.setAllele(allele, idx, ploidy=-1, chrom=-1)

Details:

    set allele allele to a locus, using its absolute index idx. If a
    ploidy ploidy and/or a chromosome indexes are given, idx is
    relative to the beginning of specified homologous copy of
    chromosomes (if chrom=-1) or the beginning of the specified
    homologous copy of specified chromosome (if chrom >= 0).

"; 

%feature("docstring") simuPOP::individual::genotype "

Usage:

    x.genotype(ploidy=-1, chrom=-1)

Details:

    return an editable array (a carray of length totNumLoci()) that
    represents all alleles on the p-th homologous set of chromosomes.
    If ploidy or chrom is given, only alleles on the chrom-th
    chromosome (or all chromosomes if chrom = -1) of ploidy-th
    homologous copy of chromosomes will be returned.

"; 

%feature("docstring") simuPOP::individual::setGenotype "

Usage:

    x.setGenotype(geno, ploidy=-1, chrom=-1)

Details:

    Fill the genotype of an individual using a list of alleles geno.
    If parameters ploidy and/or chrom are specified, alleles will be
    copied to only all or specified chromosome on selected homologous
    copy of chromosomes. geno will be reused if its length is less
    than number of alleles to be filled.

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

    x.info(field)

Details:

    Return the value of an information field filed (by index or name).
    ind.info(name) is equivalent to ind.name although the function
    form allows the use of indexes of information fieldes.

"; 

%ignore simuPOP::individual::intInfo(const uintString &field) const;

%feature("docstring") simuPOP::individual::__getattr__ "

Usage:

    x.__getattr__(field)

Details:

    read info as attribute

"; 

%feature("docstring") simuPOP::individual::__setattr__ "

Usage:

    x.__setattr__(field, value)

Details:

    write info as attribute

"; 

%feature("docstring") simuPOP::individual::setInfo "

Usage:

    x.setInfo(value, field)

Details:

    set the value of an information field field (by index or name) to
    value. ind.setInfo(value, field) is equivalent to ind.field =
    value although the function form allows the use of indexes of
    information fieldes.

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

%ignore simuPOP::individual::swap(individual &ind, bool swapContent=true);

%ignore simuPOP::individual::display(ostream &out, int width=1, const vectoru &loci=vectoru());

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

    InfoEval

Details:

    Unlike operator pyEval and pyExec that work at the population
    level, in a population's local namespace, operator infoEval works
    at the individual level, working with individual information
    fields. When this operator is applied to a population, information
    fields of eligible individuals are put into either a temporary
    dictionary or in the local namespace of the population. A Python
    expression is then evaluated for each individual. The result is
    written to an output.

Note:

    Unlike operator ``infoExec``, individual information fields are
    not updated after this operator is applied to a population.This
    operator tends to generate a large amount of output so use it is
    with caution.

"; 

%feature("docstring") simuPOP::infoEval::infoEval "

Usage:

    infoEval(expr=\"\", stmts=\"\", usePopVars=False, exposeInd=\"\",
      output=\">\", begin=0, end=-1, step=1, at=[], reps=AllAvail,
      subPops=AllAvail, infoFields=[])

Details:

    Create an operator that evaluate a Python expression expr using
    individual information fields as variables. For each eligible
    individual (individuals in (virtual) subpopulations specified by
    parameter subPops, default to all individuals), its information
    fields are copied either to a temporary namespace (default) or the
    population's local namespace (if usePopVars is True). If exposeInd
    is not empty, the individual itself will be exposed in this
    namespace as a variable with name specified by exposeInd. In the
    usePopVars=True case, any population variable whose name matches
    an information field or exposeInd will be silently overridden.  A
    Python expression (expr) is evaluated for each individual. The
    results are converted to strings and are written to an output
    specified by parameter output. Optionally, a statement (or several
    statements separated by newline) can be executed before expr is
    evaluated.  This operator is by default applied post-mating. If it
    stage is set to DuringMating, it will be applied to all offspring,
    regardless of subPops settings.

Note:

    Although expr is evaluated in individual or population level local
    namespaces, it can also access a global namespace which is the
    module namespace of your script. However, using module level
    variables and functions in this operator is discouraged.

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

%feature("docstring") simuPOP::infoEval::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::infoExec "

Function form:

    InfoExec

Details:

    Operator infoExec is similar to infoEval in that it works at the
    individual level, using individual information fields as
    variables. The difference is that instead of evaluating an
    expression and outputing its result, this operator execute one or
    more statements and update individual information fields from the
    namespace after the specified statements are execuated.

"; 

%feature("docstring") simuPOP::infoExec::infoExec "

Usage:

    infoExec(stmts=\"\", usePopVars=False, exposeInd=\"\", output=\"\",
      begin=0, end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=[])

Details:

    Create an operator that executes Python statements stmts using
    individual information fields as variables. For each eligible
    individual (individuals in (virtual) subpopulations specified by
    parameter subPops, default to all individuals), its information
    fields are copied either to a temporary namespace (default) or the
    population's local namespace (if usePopVars is True). If exposeInd
    is not empty, the individual itself will be exposed in this
    namespace as a variable with name specified by exposeInd. In the
    usePopVars=True case, any population variable whose name matches
    an information field or exposeInd will be silently overridden.
    One or more python statements (stmts) are executed for each
    individual. Information fields of these individuals are then
    updated from the corresponding variables. For example, a=1 will
    set information field a of all individuals to 1, a=b will set
    information field a of all individuals to information field b or a
    population variable b if b is not an information field but a
    population variable (needs usePopVars=True), and a=ind.sex() will
    set information field a of all individuals to its sex (needs
    exposeInd='ind'.  This operator is by default applied post-mating.
    If it stage is set to DuringMating, it will be applied to all
    offspring, regardless of subPops settings.

Note:

    Although stmts are executed in individual or population level
    local namespaces, they also have access to a global namespace
    which is the module namespace of your script. However, using
    module level variables and functions in stmts is discouraged.

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

%feature("docstring") simuPOP::infoExec::apply "

Description:

    apply the infoExec operator

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::infoExec::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::infoExec::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::infoParentsChooser "

Details:

    This parent chooser chooses an individual randomly, and then
    his/her spouse his/her spouse from a given set of information
    fields, which stores indexes of individuals in the same
    generation. An information field will be ignored if its value is
    negative, or if sex is incompatible.  Depending on what indexes
    are stored in these information fields, this parent chooser can be
    used to implement different types of mating schemes where
    selection of spouse is limited. For example, a consanguineous
    mating scheme can be implemeneted using this mating scheme if
    certain type of relatives are located for each individual, and are
    used for mating.  This parent chooser uses randomParentChooser to
    choose one parent and randomly choose another one from the
    information fields. Natural selection is supported during the
    selection of the first parent. Because of potentially uneven
    distribution of valid information fields, the overall process may
    not be as random as expected.

"; 

%feature("docstring") simuPOP::infoParentsChooser::infoParentsChooser "

Usage:

    infoParentsChooser(infoFields=[], func=None, param=None,
      selectionField=\"fitness\")

Details:

    Create a information parent chooser a parent randomly (with
    replacement, and with selection if natural selection is enabled),
    and then his/her spouse from indexes stored in infoFields. If a
    Python function func is specified, it will be called before
    parents are chosen. This function accepts the parental population
    and an optional parameter param and is usually used to locate
    qualified spouse for each parent. The return value of this
    function is ignored.

"; 

%feature("docstring") simuPOP::infoParentsChooser::clone "

Description:

    Deep copy of a infomation parent chooser.

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

%feature("docstring") simuPOP::InformationIterator::valid "

Usage:

    x.valid()

"; 

%feature("docstring") simuPOP::infoSplitter "

Details:

    This splitter defines VSPs according to the value of an
    information field of each indivdiual. A VSP is defined either by a
    value or a range of values.

"; 

%feature("docstring") simuPOP::infoSplitter::infoSplitter "

Usage:

    infoSplitter(field, values=[], cutoff=[], ranges=[], names=[])

Details:

    Create an infomration splitter using information field field. If
    parameter values is specified, each item in this list defines a
    VSP in which all individuals have this value at information field
    field. If a set of cutoff values are defined in parameter cutoff,
    individuals are grouped by intervals defined by these cutoff
    values. For example, cutoff=[1,2] defines three VSPs with v < 1, 1
    <= v < 2 and v >=2 where v is the value of an individual at
    information field field. If parameter ranges is specified, each
    range defines a VSP. For example, ranges=[[1, 3], [2, 5]] defines
    two VSPs with 1 <= v < 3 and 2 <= 3 < 5. Of course, only one of
    the parameters values, cutoff and ranges should be defined, and
    values in cutoff should be distinct, and in an increasing order. A
    default set of names are given to each VSP unless a new set of
    names is given by parameter names.

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

%ignore simuPOP::infoSplitter::contains(const population &pop, ULONG ind, vspID vsp);

%ignore simuPOP::infoSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::infoSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::infoSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of a VSP vsp, which is field = value if VSPs are
    defined by values in parameter values, or field < value (the first
    VSP), v1 <= field < v2 and field >= v (the last VSP) if VSPs are
    defined by cutoff values. A user-specified name, if specified,
    will be returned instead.

"; 

%feature("docstring") simuPOP::inheritTagger "

Details:

    An inheritance tagger passes values of parental information
    field(s) to the corresponding fields of offspring. If there are
    two parental values from parents of a sexual mating event, a
    parameter mode is used to specify how to assign offspring
    information fields.

"; 

%feature("docstring") simuPOP::inheritTagger::inheritTagger "

Usage:

    inheritTagger(mode=Paternal, begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, output=\"\", infoFields=[])

Details:

    Creates an inheritance tagger that passes values of parental
    information fields (parameter infoFields) to the corresponding
    fields of offspring. If there is only one parent, values at the
    specified information fields are copied directly. If there are two
    parents, parameter mode specifies how to pass them to an
    offspring. More specifically,
    *   mode=Maternal Passing the value from mother.
    *   mode=Paternal Passing the value from father.
    *   mode=Mean Passing the average of two values.
    *   mode=Maximum Passing the maximum value of two values.
    *   mode=Minumum Passing the minimum value of two values.
    *   mode=Summation Passing the summation of two values.
    *   mode=Multiplication Passing the multiplication of two values.
    An RuntimeError will be raised if any of the parents does not
    exist. This operator does not support parameter subPops and does
    not output any information.

"; 

%feature("docstring") simuPOP::inheritTagger::~inheritTagger "

Usage:

    x.~inheritTagger()

"; 

%feature("docstring") simuPOP::inheritTagger::description "Obsolete or undocumented function."

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

    initByFreq(alleleFreq=[], loci=AllAvail, ploidy=AllAvail,
      identicalInds=False, begin=0, end=1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    This function creates an initializer that initializes individual
    genotypes randomly, using allele frequencies specified in
    parameter alleleFreq. Elements in alleleFreq specifies the allele
    frequencies of allele 0, 1, ... respectively. These frequencies
    should add up to 1. If loci, ploidy and/or subPop are specified,
    only specified loci, ploidy, and individuals in these (virtual)
    subpopulations will be initialized. If identicalInds is True, the
    first individual in each (virtual) subpopulation will be
    initialized randomly, and be copied to all other individuals in
    this (virtual) subpopulation. If a list of frequencies are given,
    they will be used for each (virtual) subpopulation. This operator
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

%feature("docstring") simuPOP::initByFreq::description "Obsolete or undocumented function."

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

    initByValue(value=[], loci=AllAvail, ploidy=AllAvail,
      proportions=[], freq=[], begin=0, end=1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    This function creates an initializer that initializes individual
    genotypes with given genotype value. If loci, ploidy and/or subPop
    are specified, only specified loci, ploidy, and individuals in
    these (virtual) subpopulations will be initialized. value can be
    used to initialize given loci, all loci, and all homologous copies
    of these loci. If freq (a list of positive numbers that add up to
    1) is given, value should be a list of values that will be
    assigned randomly according to their respective proportion.
    Althernatively, you can use parameter proportions which assign
    values randomly, but with exact proportion. If a list of values
    are given without frequencies or proportions, they will be used
    for each (virtual) subpopulations. This operator initializes all
    chromosomes, including unused genotype locations and customized
    chromosomes.

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

%feature("docstring") simuPOP::initByValue::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::initByValue::apply "

Description:

    apply this operator to population pop

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::initInfo "

Function form:

    InitInfo

Details:

    This operator initializes given information fields with a sequence
    of values, or a user-provided function such as random.random.

"; 

%feature("docstring") simuPOP::initInfo::initInfo "

Usage:

    initInfo(values, begin=0, end=-1, step=1, at=[], reps=AllAvail,
      subPops=AllAvail, infoFields=[])

Details:

    Create an operator that initialize individual information fields
    infoFields using a sequence of values or a user-defined function.
    If a list of values are given, it will be used sequentially for
    all individuals. The values will be reused if its length is less
    than the number of individuals. The values will be assigned
    repeatedly regardless of subpopulation boundaries. If a Python
    function is given, it will be called, without any argument,
    whenever a value is needed. If a list of (virtual) subpopulation
    is specified in parameter subPop, only individuals in these
    subpopulations will be initialized.

"; 

%feature("docstring") simuPOP::initInfo::~initInfo "

Description:

    destructor

Usage:

    x.~initInfo()

"; 

%feature("docstring") simuPOP::initInfo::clone "

Description:

    deep copy of an initInfo operator.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initInfo::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::initInfo::apply "

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

    initSex(maleFreq=0.5, maleProp=-1, sex=[], begin=0, end=-1,
      step=1, at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create an operator that initialize individual sex to Male or
    Female. By default, it assign sex to individuals randomly, with
    equal probability of having a male or a female. This probabability
    can be adjusted through parameter maleFreq or be made to exact
    proportions by specifying parameter maleProp. Alternatively, a
    fixed sequence of sexes can be assigned. For example, if
    sex=[Male, Female], individuals will be assigned Male and Female
    successively. Parameter maleFreq or maleProp are ignored if sex is
    given. If a list of (virtual) subpopulation is specified in
    parameter subPop, only individuals in these subpopulations will be
    initialized. Note that the sex sequence, if used, is assigned
    repeatedly regardless of subpopulation boundaries.

"; 

%feature("docstring") simuPOP::initSex::~initSex "

Description:

    destructor

Usage:

    x.~initSex()

"; 

%feature("docstring") simuPOP::initSex::clone "

Description:

    deep copy of an initSex operator.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initSex::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::initSex::apply "

Description:

    apply this operator to population pop

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::intList "

Details:

    A class to specify replicate list. The reason why I cannot simple
    use vectori() is that users have got used to use a single number
    to specify a single replicate.

"; 

%feature("docstring") simuPOP::intList::intList "

Usage:

    intList(obj=None)

"; 

%ignore simuPOP::intList::intList(const vectori &reps);

%ignore simuPOP::intList::elems() const;

%ignore simuPOP::intList::allAvail();

%ignore simuPOP::intList::match(UINT rep, const vector< bool > *activeRep=NULL);

%feature("docstring") simuPOP::kamMutator "

Function form:

    KamMutate

Details:

    This mutator implements a k-allele mutation model that assumes k
    allelic states (alleles 0, 1, 2, ..., k-1) at each locus. When a
    mutation event happens, it mutates an allele to any other states
    with equal probability.

"; 

%feature("docstring") simuPOP::kamMutator::kamMutator "

Usage:

    kamMutator(k, rates=[], loci=AllAvail, mapIn=[], mapOut=[],
      output=\">\", begin=0, end=-1, step=1, at=[], reps=AllAvail,
      subPops=AllAvail, infoFields=[])

Details:

    Create a k-allele mutator that mutates alleles to one of the other
    k-1 alleles with equal probability. This mutator by default
    applies to all loci unless parameter loci is specified. A single
    mutation rate will be used for all loci if a single value of
    parameter rates is given. Otherwise, a list of mutation rates can
    be specified for each locus in parameter loci. Please refer to
    classes mutator and baseOperator for descriptions of other
    parameters.

"; 

%feature("docstring") simuPOP::kamMutator::~kamMutator "

Usage:

    x.~kamMutator()

"; 

%ignore simuPOP::kamMutator::mutate(AlleleRef allele, UINT locus);

%feature("docstring") simuPOP::kamMutator::clone "

Description:

    deep copy of a kamMutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::kamMutator::description "Obsolete or undocumented function."

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

    maPenetrance(loci, penetrance, wildtype=AllAvail, ancGen=-1,
      begin=0, end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=[])

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

%feature("docstring") simuPOP::maPenetrance::description "Obsolete or undocumented function."

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

    mapPenetrance(loci, penetrance, phase=False, ancGen=-1, begin=0,
      end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=[])

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

%feature("docstring") simuPOP::mapPenetrance::description "Obsolete or undocumented function."

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

    mapQuanTrait(loci, qtrait, sigma=0, phase=False, ancGen=-1,
      begin=0, end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=AllAvail)

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

%feature("docstring") simuPOP::mapQuanTrait::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::mapSelector "

Applicability: all ploidy

Details:

    This selector assigns individual fitness values using a user-
    specified dictionary.

"; 

%feature("docstring") simuPOP::mapSelector::mapSelector "

Usage:

    mapSelector(loci, fitness, begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=AllAvail)

Details:

    Create a selector that assigns individual fitness values using a
    dictionary fitness with genotype at loci as keys, and fitness as
    values. For each individual (parents if this operator is applied
    before mating, and offspring if this operator is applied during
    mating), genotypes at loci are collected one by one (e.g. p0_loc0,
    p1_loc0, p0_loc1, p1_loc1... for a diploid individual) and are
    looked up in the dictionary. If a genotype cannot be found, it
    will be looked up again without phase information (e.g. (1,0) will
    match key (0,1)). If the genotype still can not be found, a
    ValueError will be returned.

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

%feature("docstring") simuPOP::mapSelector::description "Obsolete or undocumented function."

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

    maQuanTrait(loci, qtrait, wildtype, sigma=[], ancGen=-1,
      begin=0, end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=AllAvail)

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

%feature("docstring") simuPOP::maQuanTrait::description "Obsolete or undocumented function."

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
    considered as diseased alleles.  This selector accepts an array of
    fitness values:
    *   For single-locus, fitness is the fitness for genotypes AA, Aa,
    aa, while A stands for wildtype alleles.
    *   For a two-locus model, fitness is the fitness for genotypes
    AABB, AABb, AAbb, AaBB, AbBb, Aabb, aaBB, aaBb and aaBb.
    *   For a model with more than two loci, use a table of length $
    3^{n} $ in a order similar to the two-locus model.

"; 

%feature("docstring") simuPOP::maSelector::maSelector "

Description:

    create a multiple allele selector

Usage:

    maSelector(loci, fitness, wildtype=0, begin=0, end=-1, step=1,
      at=[], reps=AllAvail, subPops=AllAvail, infoFields=AllAvail)

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

%feature("docstring") simuPOP::maSelector::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::mating "Obsolete or undocumented function."

%ignore simuPOP::mating::isCompatible(const population &pop) const;

%feature("docstring") simuPOP::mating::mating "

Usage:

    mating(subPopSize=[])

Details:

    Create a mating scheme. subPopSize can be used to determine
    subpopulatio sizes of an offspring generation.

"; 

%feature("docstring") simuPOP::mating::~mating "

Description:

    destructor

Usage:

    x.~mating()

"; 

%feature("docstring") simuPOP::mating::clone "

Description:

    deep copy of a mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::mating::submitScratch(population &pop, population &scratch);

%ignore simuPOP::mating::mateSubPop(population &pop, SubPopID subPop, RawIndIterator offBegin, RawIndIterator offEnd, vector< baseOperator * > &ops);

%ignore simuPOP::mating::mate(population &pop, population &scratch, vector< baseOperator * > &ops);

%ignore simuPOP::mating::prepareScratchPop(population &pop, population &scratch);

%feature("docstring") simuPOP::matrixMutator "

Function form:

    MatrixMutate

Details:

    A matrix mutator mutates alleles 0, 1, ..., n-1 using a n by n
    matrix, which specifies the probability at which each allele
    mutates to another. Conceptually speaking, this mutator goes
    through all mutable allele and mutate it to another state
    according to probabilities in the corresponding row of the rate
    matrix. Only one mutation rate matrix can be specified which will
    be used for all specified loci. #

"; 

%feature("docstring") simuPOP::matrixMutator::matrixMutator "

Usage:

    matrixMutator(rate, loci=AllAvail, mapIn=[], mapOut=[],
      output=\">\", begin=0, end=-1, step=1, at=[], reps=AllAvail,
      subPops=AllAvail, infoFields=[])

Details:

    Create a mutator that mutates alleles 0, 1, ..., n-1 using a n by
    n matrix rate. Item (i,j) of this matrix specifies the probability
    at which allele i mutates to allele j. Diagnal items (i, i) are
    ignored because they are automatically determined by other
    probabilities. Only one mutation rate matrix can be specified
    which will be used for all loci in the applied population, or loci
    specified by parameter loci. Please refer to classes mutator and
    baseOperator for detailed explanation of other parameters.

"; 

%feature("docstring") simuPOP::matrixMutator::~matrixMutator "

Description:

    destructor.

Usage:

    x.~matrixMutator()

"; 

%ignore simuPOP::matrixMutator::mutate(AlleleRef allele, UINT locus);

%feature("docstring") simuPOP::matrixMutator::clone "

Description:

    deep copy of a matrixMutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::matrixMutator::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::mendelianGenoTransmitter "

Details:

    This Mendelian offspring generator accepts two parents and pass
    their genotypes to an offspring following Mendel's laws. Sex
    chromosomes are handled according to the sex of the offspring,
    which is usually determined in advance by an offspring generator.
    Customized chromosomes are not handled.

"; 

%feature("docstring") simuPOP::mendelianGenoTransmitter::mendelianGenoTransmitter "

Usage:

    mendelianGenoTransmitter(output=\"\", begin=0, end=-1, step=1,
      at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a Mendelian genotype transmitter (a during-mating operator)
    that transmits genotypes from parents to offspring following
    Mendel's laws. Autosomes and sex chromosomes are handled but
    customized chromosomes are ignored. Parameters subPops and
    infoFields are ignored.

"; 

%feature("docstring") simuPOP::mendelianGenoTransmitter::clone "

Description:

    Deep copy of a Mendelian genotype transmitter.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mendelianGenoTransmitter::description "Obsolete or undocumented function."

%ignore simuPOP::mendelianGenoTransmitter::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::mendelianGenoTransmitter::initialize "

Usage:

    x.initialize(pop)

Details:

    Initialize a base genotype operator for a population. This
    function should be called before function transmitGenotype is used
    to transmit genotype.

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

Details:

    This operator merges subpopulations subPops to a single
    subpopulation. If subPops is ignored, all subpopulations will be
    merged. Virtual subpopulations are not allowed in subPops.

"; 

%feature("docstring") simuPOP::mergeSubPops::mergeSubPops "

Usage:

    mergeSubPops(subPops=AllAvail, name=\"\", begin=0, end=-1, step=1,
      at=[], reps=AllAvail, infoFields=[])

Details:

    Create an operator that merges subpopulations subPops to a single
    subpopulation. If subPops is not given, all subpopulations will be
    merged. The merged subpopulation will take the name of the first
    subpopulation being merged unless a new name is given.  This
    operator is by default applied pre-mating (parameter stage).
    Please refer to operator baseOperator for a detailed explanation
    for all parameters.

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

%feature("docstring") simuPOP::mergeSubPops::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::migrator "

Details:

    This operator migrates individuals from (virtual) subpopulations
    to other subpopulations, according to either pre-specified
    destination subpopulation stored in an information field, or
    randomly according to a migration matrix.  In the former case,
    values in a specified information field (default to migrate_to)
    are considered as destination subpopulation for each individual.
    If subPops is given, only individuals in specified (virtual)
    subpopulations will be migrated where others will stay in their
    original subpopulation. Negative values are not allowed in this
    information field because they do not represent a valid
    destination subpopulation ID.  In the latter case, a migration
    matrix is used to randomly assign destination subpoulations to
    each individual. The elements in this matrix can be probabilities
    to migrate, proportions of individuals to migrate, or exact number
    of individuals to migrate.  By default, the migration matrix
    should have m by m elements if there are m subpopulations. Element
    (i, j) in this matrix represents migration probability, rate or
    count from subpopulation i to j. If subPops (length m) and/or
    toSubPops (length n) are given, the matrix should have m by n
    elements, corresponding to specified source and destination
    subpopulations. Subpopulations in subPops can be virtual
    subpopulations, which makes it possible to migrate, for example,
    males and females at different rates from a subpopulation. If a
    subpopulation in toSubPops does not exist, it will be created. In
    case that all individuals from a subpopulation are migrated, the
    empty subpopulation will be kept.  If migration is applied by
    probability, the row of the migration matrix corresponding to a
    source subpopulation is intepreted as probabilities to migrate to
    each destination subpopulation. Each individual's detination
    subpopulation is assigned randomly according to these
    probabilities. Note that the probability of staying at the present
    subpopulation is automatically calculated so the corresponding
    matrix elements are ignored.  If migration is applied by
    proportion, the row of the migration matrix corresponding to a
    source subpopulation is intepreted as proportions to migrate to
    each destination subpopulation. The number of migrants to each
    destination subpopulation is determined before random
    indidividuals are chosen to migrate.  If migration is applied by
    counts, the row of the migration matrix corresponding to a source
    subpopulation is intepreted as number of individuals to migrate to
    each detination subpopulation. The migrants are chosen randomly.
    This operator goes through all source (virtual) subpopulations and
    assign detination subpopulation of each individual to an
    information field. An RuntimeError will be raised if an individual
    is assigned to migrate more than once. This might happen if you
    are migrating from two overlapping virtual subpopulations.

"; 

%feature("docstring") simuPOP::migrator::migrator "

Usage:

    migrator(rate=[], mode=ByProbability, toSubPops=AllAvail,
      begin=0, end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=\"migrate_to\")

Details:

    Create a migrator that moves individuals from source (virtual)
    subpopulations subPops (default to migrate from all
    subpopulations) to destination subpopulations toSubPops (default
    to all subpopulations), according to existing values in an
    information field infoFields[0], or randomly according to a
    migration matrix rate. In the latter case, the size of the matrix
    should match the number of source and destination subpopulations.
    Depending on the value of parameter mode, elements in the
    migration matrix (rate) are interpreted as either the
    probabilities to migrate from source to destination subpopulations
    (mode = ByProbability), proportions of individuals in the source
    (virtual) subpopulations to the destination subpopulations (mode =
    ByProportion), numbers of migrants in the source (virtual)
    subpopulations (mode = ByCounts), or ignored completely (mode =
    ByIndInfo). In the last case, parameter subPops is respected (only
    individuals in specified (virtual) subpopulations will migrate)
    but toSubPops is ignored.  This operator is by default applied
    pre-mating (parameter stage). Please refer to operator
    baseOperator for a detailed explanation for all parameters.

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

%ignore simuPOP::migrator::setRates(int mode, const subPopList &fromSubPops, const vectoru &toSubPops);

%feature("docstring") simuPOP::migrator::apply "

Description:

    apply the migrator to populaiton pop.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::migrator::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::mitochondrialGenoTransmitter "

Details:

    This geno transmitter assumes that the first homologous copy of
    several (or all) Customized chromosomes are copies of
    mitochondrial chromosomes. It transmits these chromosomes randomly
    from the female parent to offspring. If this transmitter is
    applied to populations with more than one homologous copies of
    chromosomes, it transmits the first homologous copy of chromosomes
    and clears alleles (set to zero) on other homologous copies.

"; 

%feature("docstring") simuPOP::mitochondrialGenoTransmitter::mitochondrialGenoTransmitter "

Usage:

    mitochondrialGenoTransmitter(output=\"\", chroms=[], begin=0,
      end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=[])

Details:

    Createa a mitochondrial genotype transmitter that treats all
    Customized chromosomes, or a list of chromosomes specified by
    chroms, as human mitochondrial chromosomes. These chromosomes
    should have the same length and the same number of loci. This
    operator transmits these chromosomes randomly from the female
    parent to offspring of both sexes.

"; 

%feature("docstring") simuPOP::mitochondrialGenoTransmitter::clone "

Description:

    Deep copy of a mitochondrial genotype transmitter.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mitochondrialGenoTransmitter::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::mitochondrialGenoTransmitter::initialize "Obsolete or undocumented function."

%ignore simuPOP::mitochondrialGenoTransmitter::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::mixedMutator "

Function form:

    MixedMutate

Details:

    This mixed mutator accepts a list of mutators and use one of them
    to mutate an allele when an mutation event happens.

"; 

%feature("docstring") simuPOP::mixedMutator::mixedMutator "

Usage:

    mixedMutator(rates=[], loci=AllAvail, mutators=[], prob=[],
      mapIn=[], mapOut=[], context=0, output=\">\", begin=0, end=-1,
      step=1, at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a mutator that randomly chooses one of the specified
    mutators to mutate an allele when a mutation event happens. The
    mutators are choosen according to a list of probabilities ( prob)
    that should add up to 1. The passed and returned alleles might be
    changed if parameters mapIn and mapOut are used. Most parameters,
    including loci, mapIn, mapOut, rep, and subPops of mutators
    specified in parameter mutators are ignored. This mutator by
    default applies to all loci unless parameter loci is specified.
    Please refer to classes mutator and baseOperator for descriptions
    of other parameters.

"; 

%feature("docstring") simuPOP::mixedMutator::clone "

Description:

    deep copy of a mixedMutator

Usage:

    x.clone()

"; 

%ignore simuPOP::mixedMutator::initialize(population &pop);

%ignore simuPOP::mixedMutator::mutate(AlleleRef allele, UINT locus);

%feature("docstring") simuPOP::mixedMutator::description "Obsolete or undocumented function."

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
    *   PEN_Multiplicative: the penetrance is calculated as $ f=\\prod
    f_{i} $.
    *   PEN_Additive: the penetrance is calculated as $
    f=\\min\\left(1,\\sum f_{i}\\right) $. $ f $ will be set to 1 when $
    f<0 $. In this case, $ s_{i} $ are added, not $ f_{i} $ directly.
    *   PEN_Heterogeneity: the penetrance is calculated as $
    f=1-\\prod\\left(1-f_{i}\\right) $. Please refer to Neil Risch (1990)
    for detailed information about these models.

"; 

%feature("docstring") simuPOP::mlPenetrance::mlPenetrance "

Description:

    create a multiple locus penetrance operator

Usage:

    mlPenetrance(peneOps, mode=Multiplicative, ancGen=-1, begin=0,
      end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=[])

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

%feature("docstring") simuPOP::mlPenetrance::description "Obsolete or undocumented function."

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
    *   Multiplicative: the mean of the quantitative trait is
    calculated as $ f=\\prod f_{i} $.
    *   Additive: the mean of the quantitative trait is calculated as
    $ f=\\sum f_{i} $. Note that all $ \\sigma_{i} $ (for $ f_{i} $) and
    $ \\sigma $ (for $ f $) will be considered. I.e, the trait value
    should be
    $ f=\\sum_{i}\\left(f_{i}+N\\left(0,\\sigma_{i}^{2}\\right)\\right)+\\sig
    ma^{2} $ for Additive case. If this is not desired, you can set
    some of the $ \\sigma $ to zero.

"; 

%feature("docstring") simuPOP::mlQuanTrait::mlQuanTrait "

Description:

    create a multiple locus quantitative trait operator

Usage:

    mlQuanTrait(qtraits, mode=Multiplicative, sigma=0, ancGen=-1,
      begin=0, end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=AllAvail)

Details:

    Please refer to quanTrait for other parameter descriptions.

Arguments:

    qtraits:        a list of quantitative traits
    mode:           can be one of Multiplicative and Additive

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

%feature("docstring") simuPOP::mlQuanTrait::description "Obsolete or undocumented function."

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
    *   Multiplicative: the fitness is calculated as $
    f=\\prod_{i}f_{i} $, where $ f_{i} $ is the single-locus fitness
    value.
    *   Additive: the fitness is calculated as $
    f=\\max\\left(0,1-\\sum_{i}(1-f_{i})\\right) $. $ f $ will be set to 0
    when $ f<0 $.

"; 

%feature("docstring") simuPOP::mlSelector::mlSelector "

Description:

    create a multiple-locus selector

Usage:

    mlSelector(selectors, mode=Multiplicative, begin=0, end=-1,
      step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=AllAvail)

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

%feature("docstring") simuPOP::mlSelector::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::mutator "

Details:

    Class mutator is the base class of all mutators. It handles all
    the work of picking an allele at specified loci from certain
    (virtual) subpopulation with certain probability, and calling a
    derived mutator to mutate the allele. Alleles can be changed
    before and after mutation if existing allele numbers do not match
    those of a mutation model.

"; 

%feature("docstring") simuPOP::mutator::mutator "

Usage:

    mutator(rates=[], loci=AllAvail, mapIn=[], mapOut=[], context=0,
      output=\">\", begin=0, end=-1, step=1, at=[], reps=AllAvail,
      subPops=AllAvail, infoFields=[])

Details:

    A mutator mutates alleles from one state to another with given
    probability. This base mutator does not perform any mutation but
    it defines common behaviors of all mutators.  By default, a
    mutator mutates all alleles in all populations of a simulator at
    all generations. A number of parameters can be used to restrict
    mutations to certain generations (parameters begin, end, step and
    at), replicate populations (parameter rep), (virtual)
    subpopulations (parameter subPops) and loci (parameter loci).
    Please refer to class baseOperator for a detailed explanation of
    these parameters.  Parameter rate or its equivalence specifies the
    probability that a a mutation event happens. The exact form and
    meaning of rate is mutator-specific. If a single rate is
    specified, it will be applied to all loci. If a list of mutation
    rates are given, they will be applied to each locus specified in
    parameter loci. Note that not all mutators allow specification of
    multiple mutation rate, especially when the mutation rate itself
    is a list or matrix.  Alleles at a locus are non-negative numbers
    0, 1, ... up to the maximum allowed allele for the loaded module
    (1 for binary, 255 for short and 65535 for long modules). Whereas
    some general mutation models treat alleles as numbers, other
    models assume specific interpretation of alleles. For example, an
    acgtMutator assumes alleles 0, 1, 2 and 3 as nucleotides A, C, G
    and T. Using a mutator that is incompatible with your simulation
    will certainly yield erroneous results.  If your simulation
    assumes different alleles with a mutation model, you can map an
    allele to the allele used in the model and map the mutated allele
    back. This is achieved using a mapIn list with its i-th item being
    the corresponding allele of real allele i, and a mapOut list with
    its i-th item being the real allele of allele i assumed in the
    model. For example mapIn=[0, 0, 1] and mapOut=[1, 2] would allow
    the use of a snpMutator to mutate between alleles 1 and 2, instead
    of 0 and 1. Parameters mapIn and mapOut also accept a user-defined
    Python function that returns a corresponding allele for a given
    allele. This allows easier mapping between a large number of
    alleles and advanced models such as random emission of alleles.
    Some mutation models are context dependent. Namely, how an allele
    mutates will depend on its adjecent alleles. Whereas most simuPOP
    mutators are context independent, some of them accept a parameter
    context which is the number of alleles to the left and right of
    the mutated allele. For example context=1 will make the alleles to
    the immediate left and right to a mutated allele available to a
    mutator. These alleles will be mapped in if parameter mapIn is
    defined. How exactly a mutator makes use of these information is
    mutator dependent.

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

%ignore simuPOP::mutator::setRate(const vectorf &rates, const uintList &loci);

%ignore simuPOP::mutator::mutRate(UINT loc);

%ignore simuPOP::mutator::mutate(AlleleRef allele, UINT locus);

%feature("docstring") simuPOP::mutator::fillContext "

Description:

    a rarely used feature, performance should be a secondary
    consideration.

Usage:

    x.fillContext(pop, ptr, locus)

"; 

%ignore simuPOP::mutator::setContext(int context);

%ignore simuPOP::mutator::context();

%feature("docstring") simuPOP::mutator::apply "

Description:

    Apply a mutator.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::noneOp "

Details:

    This operator does nothing when it is applied to a population. It
    is usually used as a placeholder when an operator is needed
    syntactically.

"; 

%feature("docstring") simuPOP::noneOp::noneOp "

Usage:

    noneOp(output=\">\", begin=0, end=0, step=1, at=[], reps=AllAvail,
      subPops=AllAvail, infoFields=[])

Details:

    Create a noneOp.

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

%feature("docstring") simuPOP::noneOp::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::offspringGenerator "

Details:

    An offspring generator generates offspring from parents chosen by
    a parent chooser. It is responsible for creating a certain number
    of offspring, determinning their sex, and transmitting genotypes
    from parents to offspring.

"; 

%feature("docstring") simuPOP::offspringGenerator::offspringGenerator "

Usage:

    offspringGenerator(ops, numOffspring=1, sexMode=RandomSex)

Details:

    Create a basic offspring generator. This offspring generator uses
    ops genotype transmitters to transmit genotypes from parents to
    offspring. It expects numParents from an upstream parents chooser
    and raises an RuntimeError if incorrect number of parents are
    passed. If both one and two parents can be handled, 0 should be
    specified for this parameter.  A number of during-mating operators
    (parameter ops) can be used to, among other possible duties such
    as setting information fields of offspring, transmit genotype from
    parents to offspring. Additional during-mating operators passed
    from the simulator.evolve() function will be applied afterwards.
    This general offspring generator does not have any default during-
    mating operator but all stock mating schemes use an offspring
    generator with a default operator. For example, a
    mendelianOffspringGenerator is used by randomMating to trasmit
    genotypes. Note that applicability parameters begin, step, end, at
    and reps could be used in these operators but negative population
    and generation indexes are unsupported.  Parameter numOffspring is
    used to control the number of offspring per mating event, or in
    another word the number of offspring in each family. It can be a
    number, a function, or a mode parameter followed by some optional
    arguments. If a number is given, given number of offspring will be
    generated at each mating event. If a Python function is given, it
    will be called each time when a mating event happens. Current
    generation number will be passed to this function, and its return
    value will be considered the number of offspring. In the last
    case, a tuple (or a list) in one of the following forms:
    (GeometricDistribution, p), (PoissonDistribution, p),
    (BinomialDistribution, p, N), or (UniformDistribution, a, b) can
    be given. The number of offspring will be determined randomly
    following these statistical distributions. Please refer to the
    simuPOP user's guide for a detailed description of these
    distributions and their parameters.  Parameter sexMode is used to
    control the sex of each offspring. Its default value is usually
    RandomSex which assign Male or Female to each individual randomly,
    with equal probabilities. If NoSex is given, all individuals will
    be Male. sexMode can also be one of (ProbOfMale, p), (NumOfMale,
    n), and (NumOfFemale, n). The first case specifies the probability
    of male for each offspring. The next two cases specifies the
    number of male or female individuals in each family, respectively.
    If n is greater than or equal to the number of offspring in this
    family, all offspring in this family will be Male or Female.

"; 

%feature("docstring") simuPOP::offspringGenerator::~offspringGenerator "

Usage:

    x.~offspringGenerator()

"; 

%feature("docstring") simuPOP::offspringGenerator::clone "

Description:

    Make a deep copy of this offspring generator.

Usage:

    x.clone()

"; 

%ignore simuPOP::offspringGenerator::initialize(const population &pop, SubPopID subPop, vector< baseOperator * > const &ops);

%ignore simuPOP::offspringGenerator::generateOffspring(population &pop, individual *dad, individual *mom, RawIndIterator &offBegin, RawIndIterator &offEnd, vector< baseOperator * > &ops);

%ignore simuPOP::offspringGenerator::finalize(const population &pop);

%ignore simuPOP::offspringGenerator::initialized();

%ignore simuPOP::offspringGenerator::numOffspring(int gen);

%ignore simuPOP::offspringGenerator::getSex(int count);

%feature("docstring") simuPOP::opList "

"; 

%feature("docstring") simuPOP::opList::opList "

Usage:

    opList(ops=[])

"; 

%ignore simuPOP::opList::opList(const opList &rhs);

%feature("docstring") simuPOP::opList::~opList "

Usage:

    x.~opList()

"; 

%ignore simuPOP::opList::begin();

%ignore simuPOP::opList::end();

%ignore simuPOP::opList::size() const;

%ignore simuPOP::opList::empty() const;

%ignore simuPOP::opList::elems() const;

%ignore simuPOP::OstreamManager;

%ignore simuPOP::OstreamManager::OstreamManager();

%feature("docstring") simuPOP::OstreamManager::~OstreamManager "

Usage:

    x.~OstreamManager()

"; 

%ignore simuPOP::OstreamManager::getOstream(const string &name, bool readable, bool realAppend, bool useString);

%ignore simuPOP::OstreamManager::hasOstream(const string &filename);

%ignore simuPOP::OstreamManager::listAll();

%ignore simuPOP::OstreamManager::closeOstream(const string &filename);

%ignore simuPOP::OstreamManager::closeAll();

%feature("docstring") simuPOP::parentChooser "

Details:

    A parent chooser repeatedly chooses parent(s) from a parental
    population and pass them to an offspring generator. A parent
    chooser can select one or two parents, which should be matched by
    the offspring generator. This class is the base class of all
    parent choosers, and should not be used directly.

"; 

%feature("docstring") simuPOP::parentChooser::parentChooser "

Usage:

    parentChooser(selectionField=\"\")

"; 

%feature("docstring") simuPOP::parentChooser::clone "

Description:

    Deep copy of a parent chooser.

Usage:

    x.clone()

"; 

%ignore simuPOP::parentChooser::initialize(population &pop, SubPopID subPop);

%ignore simuPOP::parentChooser::finalize(population &pop, SubPopID subPop);

%ignore simuPOP::parentChooser::initialized() const;

%ignore simuPOP::parentChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::parentChooser::~parentChooser "

Description:

    destructor

Usage:

    x.~parentChooser()

"; 

%feature("docstring") simuPOP::parentsTagger "

Details:

    This tagging operator records the indexes of parents (relative to
    the parental generation) of each offspring in specified
    information fields ( default to father_idx and mother_idx). Only
    one information field should be specified if an asexsual mating
    scheme is used so there is one parent for each offspring.
    Information recorded by this operator is intended to be used to
    look up parents of each individual in multi-generational
    population.

"; 

%feature("docstring") simuPOP::parentsTagger::parentsTagger "

Usage:

    parentsTagger(begin=0, end=-1, step=1, at=[], reps=AllAvail,
      subPops=AllAvail, output=\"\", infoFields=[\"father_idx\",
      \"mother_idx\"])

Details:

    Create a parents tagger that records the indexes of parents of
    each offspring when it is applied to an offspring during-mating.
    If two information fields are specified (parameter infoFields,
    with default value ['father_idx', 'mother_idx']), they are used to
    record the indexes of each individual's father and mother. Value
    -1 will be assigned if any of the parent is missing. If only one
    information field is given, it will be used to record the index of
    the first valid parent (father if both parents are valid). This
    operator ignores parameters stage, output, and subPops.

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

%feature("docstring") simuPOP::parentsTagger::description "Obsolete or undocumented function."

%ignore simuPOP::parentsTagger::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::pause "

Details:

    This operator pauses the evolution of a simulator at given
    generations or at a key stroke. When a simulator is stopped, you
    can go to a Python shell to examine the status of an evolutionary
    process, resume or stop the evolution.

"; 

%feature("docstring") simuPOP::pause::pause "

Usage:

    pause(stopOnKeyStroke=False, prompt=True, output=\">\", begin=0,
      end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=[])

Details:

    Create an operator that pause the evolution of a population when
    it is applied to this population. If stopOnKeyStroke is False
    (default), it will always pause a population when it is applied,
    if this parameter is set to True, the operator will pause a
    population if *any* key has been pressed. If a specific character
    is set, the operator will stop when this key has been pressed.
    This allows, for example, the use of several pause operators to
    pause different populations.  After a population has been paused,
    a message will be displayed (unless prompt is set to False) and
    tells you how to proceed. You can press 's' to stop the evolution
    of this population, 'S' to stop the evolution of all populations,
    or 'p' to enter a Python shell. The current population will be
    available in this Python shell as \"pop_X_Y\" when X is generation
    number and Y is replicate number. The evolution will continue
    after you exit this interactive Python shell.

Note:

    Ctrl-C will be intercepted even if a specific character is
    specified in parameter stopOnKeyStroke.

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

%feature("docstring") simuPOP::pause::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::pedigree "

Details:

    The pedigree class is derived from the population class. Unlike a
    population class that emphasizes on individual properties, the
    pedigree class emphasizes on relationship between individuals.  A
    pedigree class can be created from a population, or loaded from a
    disk file, which is usually saved by an operator during a previous
    evolutionary process. Depending on how a pedigree is saved, sex
    and affection status information may be missing.

"; 

%feature("docstring") simuPOP::pedigree::pedigree "

Usage:

    pedigree(pop, loci=[], infoFields=[], ancGen=-1,
      fatherField=\"father_idx\", motherField=\"mother_idx\")

Details:

    Create a pedigree object from a population, using a subset of loci
    (parameter loci, default to no loci), information fields
    (parameter infoFields, default to no information field except for
    parentFields), and ancestral generations (parameter ancGen,
    default to all ancestral generations). By default, information
    field father_idx and mother_idx are used to locate parents. If
    individuals in a pedigree has only one parent, the information
    field that stores parental indexes should be specified in
    parameter fatherField or motherField. The other field should be
    set to an empty string.

"; 

%ignore simuPOP::pedigree::pedigree(const pedigree &rhs);

%feature("docstring") simuPOP::pedigree::clone "

Usage:

    x.clone()

Details:

    Create a cloned copy of a pedigree.

"; 

%feature("docstring") simuPOP::pedigree::numParents "

Usage:

    x.numParents()

Details:

    Return the number of parents each individual has. This function
    returns the number of information fields used to store parental
    indexes, even if one of the fields are unused.

"; 

%feature("docstring") simuPOP::pedigree::father "

Usage:

    x.father(idx, subPop)

Details:

    Return the index of the father of individual idx in subpopulation
    subPop in the parental generation. Return -1 if this individual
    has no father (fatherField is empty or the value of information
    field is negative).

"; 

%feature("docstring") simuPOP::pedigree::mother "

Usage:

    x.mother(idx, subPop)

Details:

    Return the index of the mother of individual idx in subpopulation
    subPop in the parental generation. Return -1 if this individual
    has no mother (motherField is empty or the value of information
    field is negative).

"; 

%feature("docstring") simuPOP::pedigree::locateRelatives "

Usage:

    x.locateRelatives(relType=[], relFields=[], ancGen=-1)

Details:

    This function locates relatives (of type relType) of each
    individual and store their indexes in specified information fields
    relFields. The length of relFields determines how many relatives
    an individual can have.  Parameter relType specifies what type of
    relative to locate. It can be Self, Spouse (having at least one
    common offspring), Offspring, FullSibling (having common father
    and mother), Sibling (having at least one common parent) or
    SpouseAndOffspring (One spouse and their common offspring).
    Optionally, you can specify the sex of relatives you would like to
    locate, in the form of relType=(type, sexChoice). sexChoice can be
    AnySex (default), MaleOnly, FemaleOnly, SameSex or OppositeSex.
    sexChoice for SpouseAndOffspring only refer to sex of offspring.
    This function will by default go through all ancestral generations
    and locate relatives for all individuals. This can be changed by
    setting parameter ancGen to the greatest ancestral generation you
    would like to process.

"; 

%feature("docstring") simuPOP::pedigree::traceRelatives "

Usage:

    x.traceRelatives(pathGen, pathFields, pathSex=[],
      resultFields=[])

Details:

    Trace a relative path in a population and record the result in the
    given information fields resultFields. This function is used to
    locate more distant relatives based on the relatives located by
    function locateRelatives. For example, after siblings and
    offspring of all individuals are located, you can locate mother's
    sibling's offspring using a relative path, and save their indexes
    in each individuals information fields resultFields.  A relative
    path consits of three pieces of information specified by three
    parameters. Parameter pathGen specifies starting, intermediate and
    ending generations. pathFields specifies which information fields
    to look for at each step, and pathSex specifies sex choices at
    each generation, which should be a list of AnySex, MaleOnly,
    FemaleOnly, SameSex and OppsiteSex. The default value for this
    paramter is AnySex at all steps. The length of pathGen should be
    one more than pathFields, and pathSex if pathSex is given.  For
    example, if pathGen=[0, 1, 1, 0], pathFields = [['father_idx',
    'mother_idx'], ['sib1', 'sib2'], ['off1', 'off2']], and pathSex =
    [AnySex, MaleOnly, FemaleOnly], this function will locate
    father_idx and mother_idx for each individual at generation 0,
    find all individuals referred by father_idx and mother_idx at
    generation 1, find informaton fields sib1 and sib2 from these
    parents and locate male individuals referred by these two
    information fields. Finally, the information fields off1 and off2
    from these siblings are located and are used to locate their
    female offspring at the present geneartion. The results are father
    or mother's brother's daughters. Their indexes will be saved in
    each individuals information fields resultFields. Note that this
    function will locate and set relatives for individuals only at the
    starting generation specified at pathGen[0].

"; 

%feature("docstring") simuPOP::pedigreeMating "

Details:

    A pedigree mating scheme that evolves a population following a
    pedigree object.

"; 

%feature("docstring") simuPOP::pedigreeMating::pedigreeMating "

Usage:

    pedigreeMating(ped, generator, setSex=False, setAffection=False,
      copyFields=[])

Details:

    Creates a mating scheme that evolve a population following a
    pedigree object ped. Considering this pedigree as a population
    with N ancestral generations, the starting population is the
    greatest ancestral generation of ped. The mating scheme creates an
    offspring generation that match the size of generation N-1 and
    chooses parents according to the parents of individuals at this
    generation. Depending on the gen parameter of the simulator, the
    process continues generation by generation for N generations if
    gen >= N), or gen generations if gen < N. During the evolution, an
    offspring generator generator is used to produce one offspring at
    a time, regardless of the numOffspring setting of this offspring
    generator. If individuals in pedigree ped has only one parent, the
    offspring generator should be compatible.  By default, the
    pedigree mating scheme does not set offspring sex and affection
    status using sex and affection status of corresponding individuals
    in the pedigree. However, if such information is valid in the
    pedigree object ped, you can set parameters setSex and/or
    setAffection to True to set sex and/of affection status to
    offspring during the evolutionary process. Similarly, you can
    specify some information fields in copyFields to copy some
    information fields from pedigree to the evolving population. Note
    that these information will be copied also to the starting
    population (from the greatest ancestral generation in ped).

"; 

%feature("docstring") simuPOP::pedigreeMating::~pedigreeMating "

Description:

    destructor

Usage:

    x.~pedigreeMating()

"; 

%ignore simuPOP::pedigreeMating::pedigreeMating(const pedigreeMating &rhs);

%feature("docstring") simuPOP::pedigreeMating::clone "

Description:

    deep copy of a Python mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::pedigreeMating::prepareScratchPop(population &pop, population &scratch);

%ignore simuPOP::pedigreeMating::mate(population &pop, population &scratch, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::pedigreeTagger "

Details:

    This tagging operator records the ID of parents of each offspring
    in specified information fields (default to father_id and
    mother_id). Only one information field should be specified if an
    asexsual mating scheme is used so there is one parent for each
    offspring. Information recorded by this operator is intended to be
    used to record full pedigree information of an evolutionary
    process.

"; 

%feature("docstring") simuPOP::pedigreeTagger::pedigreeTagger "

Usage:

    pedigreeTagger(idField=\"ind_id\", output=\"\", begin=0, end=-1,
      step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=[\"father_id\", \"mother_id\"])

Details:

    Create a pedigree tagger that records the ID of parents of each
    offspring when it is applied to an offspring during-mating. If two
    information fields are specified (parameter infoFields, with
    default value ['father_id', 'mother_id']), they are used to record
    the ID of each individual's father and mother stored in the
    idField (default to ind_id) field of the parents. Value -1 will be
    assigned if any of the parent is missing. If only one information
    field is given, it will be used to record the ID of the first
    valid parent (father if both pedigree are valid).  This operator
    by default does not send any output but will output the ID of
    offspring, father, and mother if a valid output stream is
    specified. The output will be in the format of off_id father_id
    mother_id. father_id or mother_id will be ignored if only one
    parent is involved. This operator ignores parameter stage, and
    subPops.

"; 

%feature("docstring") simuPOP::pedigreeTagger::~pedigreeTagger "

Usage:

    x.~pedigreeTagger()

"; 

%feature("docstring") simuPOP::pedigreeTagger::clone "

Description:

    deep copy of a pedigreeTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pedigreeTagger::description "Obsolete or undocumented function."

%ignore simuPOP::pedigreeTagger::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::pointMutator "

Function form:

    PointMutate

Details:

    A point mutator is different from all other mutators because
    mutations in this mutator do not happen randomly. Instead, it
    happens to specific loci and mutate an allele to a specific state,
    regardless of its original state. This mutator is usually used to
    introduce a mutant to a population.

"; 

%feature("docstring") simuPOP::pointMutator::pointMutator "

Usage:

    pointMutator(loci, allele, ploidy=0, inds=[], output=\">\",
      begin=0, end=-1, step=1, at=[], reps=AllAvail, subPops=0,
      infoFields=[])

Details:

    Create a point mutator that mutates alleles at specified loci to a
    given allele of individuals inds. If there are multiple alleles at
    a locus (e.g. individuals in a diploid population), only the first
    allele is mutated unless indexes of alleles are listed in
    parameter ploidy. This operator is by default applied to
    individuals in the first subpopulation but you can apply it to a
    different or more than one (virtual) subpopulations using
    parameter *subPops* (``AllAvail`` is also accepted). Please refer
    to class baseOperator for detailed descriptions of other
    parameters.

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

%feature("docstring") simuPOP::pointMutator::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::polyParentsChooser "

Details:

    This parent chooser is similar to random parents chooser but
    instead of selecting a new pair of parents each time, one of the
    parents in this parent chooser will mate with several spouses
    before he/she is replaced. This mimicks multi-spouse mating
    schemes such as polygyny or polyandry in some populations. Natural
    selection is supported for both sexes.

"; 

%feature("docstring") simuPOP::polyParentsChooser::polyParentsChooser "

Usage:

    polyParentsChooser(polySex=Male, polyNum=1,
      selectionField=\"fitness\")

Details:

    Create a multi-spouse parents chooser where each father (if
    polySex is Male) or mother (if polySex is Female) has polyNum
    spouses. The parents are chosen with replacement. If natural
    selection is enabled, the probability that an individual is chosen
    is proportional to his/her fitness value among all individuals
    with the same sex. Selection will be disabled if specified
    information field selectionField (default to \"fitness\") does not
    exist.

"; 

%feature("docstring") simuPOP::polyParentsChooser::clone "

Description:

    Deep copy of a parent chooser.

Usage:

    x.clone()

"; 

%ignore simuPOP::polyParentsChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::polyParentsChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::population "

Details:

    A simuPOP population consists of individuals of the same genotypic
    structure, organized by generations, subpopulations and virtual
    subpopulations. It also contains a Python dictionary that is used
    to store arbitrary population variables.  In addition to genotypic
    structured related functions provided by the GenoStruTrait class,
    the population class provides a large number of member functions
    that can be used to
    *   Create, copy and compare populations.
    *   Manipulate subpopulations. A population can be divided into
    several subpopulations. Because individuals only mate with
    individuals within the same subpopulation, exchange of genetic
    information across subpopulations can only be done through
    migration. A number of functions are provided to access
    subpopulation structure information, and to merge and split
    subpopulations.
    *   Define and access virtual subpopulations. A virtual
    subpopulation splitter can be assigned to a population, which
    defines groups of individuals called virtual subpopulations (VSP)
    within each subpopulation.
    *   Access individuals individually, or through iterators that
    iterate through individuals in (virtual) subpopulations.
    *   Access genotype and information fields of individuals at the
    population level. From a population point of view, all genotypes
    are arranged sequentially individual by individual. Please refer
    to class individual for an introduction to genotype arragement of
    each individual.
    *   Store and access ancestral generations. A population can save
    arbitrary number of ancestral generations. It is possible to
    directly access an ancestor, or make an ancestral generation the
    current generation for more efficient access.
    *   Insert or remove loci, resize (shrink or expand) a population,
    sample from a population, or merge with other populations.
    *   Manipulate population variables and evaluate expressions in
    this local namespace.
    *   Save and load a population.

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
                    chromosome. lociPos should be arranged chromosome
                    by chromosome. If lociPos are not in order within
                    a chromosome, they will be re-arranged along with
                    corresponding lociNames (if specified).
    ancGen:         Number of the most recent ancestral generations to
                    keep during evolution. Default to 0, which means
                    only the current generation will be kept. If it is
                    set to -1, all ancestral generations will be kept
                    in this population (and exhaust your computer RAM
                    quickly).
    chromNames:     A list of chromosome names. Default to '' for all
                    chromosomes.
    alleleNames:    A list or a nested list of allele names. If a list
                    of alleles is given, it will be used for all loci
                    in this population. For example,
                    alleleNames=('A','C','T','G') gives names A, C, T,
                    and G to alleles 0, 1, 2, and 3 respectively. If a
                    nested list of names is given, it should specify
                    alleles names for all loci.
    lociNames:      A list of names for each locus. It can be empty or
                    a list of unique names for each locus. If loci are
                    not specified in order, loci names will be
                    rearranged according to their position on the
                    chromosome.
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

%ignore simuPOP::population::validate(const string &msg) const;

%ignore simuPOP::population::fitSubPopStru(const vectoru &newSubPopSizes, const vectorstr &newSubPopNames);

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

%feature("docstring") simuPOP::population::fitGenoStru "Obsolete or undocumented function."

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

    x.subPopSize(subPop=[])

Details:

    Return the size of a subpopulation (subPopSize(sp)) or a virtual
    subpopulation (subPopSize([sp, vsp])). If no subpop is given, it
    is the same as popSize().

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

    Return the \"spName - vspName\" (virtual named subpopulation), \"\"
    (unnamed non-virtual subpopulation), \"spName\" (named
    subpopulation) or \"vspName\" (unnamed virtual subpopulation),
    depending on whether subPopulation is named or if subPop is
    virtual.

"; 

%feature("docstring") simuPOP::population::subPopNames "

Usage:

    x.subPopNames()

Details:

    Return the names of all subpopulations (excluding virtual
    subpopulations). An empty string will be returned for unnamed
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

    x.ind(idx, subPop=[])

Details:

    Return a refernce to individual idx in the population (if
    subPop=[], default) or a subpopulation (if subPop=sp). Virtual
    subpopulation is not supported.

"; 

%feature("docstring") simuPOP::population::indByID "

Usage:

    x.indByID(id, ancGen=-1, idField=\"ind_id\")

Details:

    Return a reference to individual with id stored in information
    field idField (default to ind_id). This function by default search
    the present and all ancestral generations (ancGen=-1), but you can
    specify a specific generation if you know which generation to
    search (ancGen=0 for present generation, ancGen=1 for parental
    generation, and so on). If no individual with id is found,an
    IndexError will be raised.

"; 

%ignore simuPOP::population::ind(ULONG idx, vspID subPop=vspID()) const;

%feature("docstring") simuPOP::population::ancestor "

Usage:

    x.ancestor(idx, gen, subPop=[])

Details:

    Return a reference to individual idx in ancestral generation gen.
    The correct individual will be returned even if the current
    generation is not the present one (see also useAncestralGen). If a
    valid subPop is specified, index is relative to that subPop.
    Virtual subpopulation is not supported.

"; 

%ignore simuPOP::population::ancestor(ULONG idx, UINT gen, vspID subPop=vspID()) const;

%feature("docstring") simuPOP::population::individuals "

Usage:

    x.individuals(subPop=[])

Details:

    Return an iterator that can be used to iterate through all
    individuals in a population (if subPop=[], default), or a
    (virtual) subpopulation (subPop=spID or (spID, vspID)).

"; 

%ignore simuPOP::population::indOrdered() const;

%ignore simuPOP::population::setIndOrdered(bool s) const;

%ignore simuPOP::population::indIterator(IterationType type=VisibleInds);

%ignore simuPOP::population::indIterator(UINT subPop, IterationType type=VisibleInds);

%ignore simuPOP::population::indIterator(IterationType type=VisibleInds) const;

%ignore simuPOP::population::rawIndBegin();

%ignore simuPOP::population::rawIndEnd();

%ignore simuPOP::population::rawIndBegin(UINT subPop);

%ignore simuPOP::population::rawIndEnd(UINT subPop);

%ignore simuPOP::population::rawIndBegin() const;

%ignore simuPOP::population::alleleIterator(UINT locus);

%ignore simuPOP::population::alleleIterator(UINT locus, UINT subPop);

%ignore simuPOP::population::genoBegin(bool order);

%ignore simuPOP::population::genoEnd(bool order);

%ignore simuPOP::population::genoBegin(UINT subPop, bool order);

%ignore simuPOP::population::genoEnd(UINT subPop, bool order);

%ignore simuPOP::population::indGenoBegin(ULONG ind) const;

%ignore simuPOP::population::indGenoEnd(ULONG ind) const;

%feature("docstring") simuPOP::population::genotype "

Usage:

    x.genotype(subPop=[])

Details:

    Return an editable array of the genotype of all individuals in a
    population (if subPop=[], default), or individuals in a
    subpopulation subPop. Virtual subpopulation is unsupported.

"; 

%feature("docstring") simuPOP::population::setGenotype "

Usage:

    x.setGenotype(geno, subPop=[])

Details:

    Fill the genotype of all individuals in a population (if
    subPop=[]) or in a (virtual) subpopulation subPop (if subPop=sp or
    (sp, vsp)) using a list of alleles geno. geno will be reused if
    its length is less than subPopSize(subPop)*totNumLoci()*ploidy().

"; 

%feature("docstring") simuPOP::population::setSubPopByIndInfo "

Usage:

    x.setSubPopByIndInfo(field)

Details:

    Rearrange individuals to their new subpopulations according to
    their integer values at information field field (value returned by
    individual::indInfo(field)). Individuals with negative values at
    this field will be removed. Existing subpopulation names are kept.
    New subpopulations will have empty names.

"; 

%feature("docstring") simuPOP::population::splitSubPop "

Usage:

    x.splitSubPop(subPop, sizes, names=[])

Details:

    Split subpopulation subPop into subpopulations of given sizes,
    which should add up to the size of subpopulation subPop. If subPop
    is not the last subpopulation, indexes of subpopulations after
    subPop are shifted. If subPop is named, the same name will be
    given to all new subpopulations unless a new set of names are
    specified for these subpopulations. This function returns the IDs
    of split subpopulations.

"; 

%feature("docstring") simuPOP::population::removeSubPops "

Usage:

    x.removeSubPops(subPops)

Details:

    Remove subpopulation(s) subPop and all their individuals. Indexes
    of subpopulations after removed subpopulations will be shifted.

"; 

%feature("docstring") simuPOP::population::removeIndividuals "

Usage:

    x.removeIndividuals(inds)

Details:

    remove individual(s) inds (absolute indexes) from the current
    population. A subpopulation will be kept even if all individuals
    from it are removed. This function only affects the current
    generation.

"; 

%feature("docstring") simuPOP::population::mergeSubPops "

Usage:

    x.mergeSubPops(subPops=[], name=\"\")

Details:

    Merge subpopulations subPops. If subPops is empty (default), all
    subpopulations will be merged. subPops do not have to be adjacent
    to each other. They will all be merged to the subpopulation with
    the smallest subpopulation ID. Indexes of the rest of the
    subpopulation may be changed. A new name can be assigned to the
    merged subpopulation through parameter name (an empty name will be
    ignored). This function returns the ID of the merged
    subpopulation.

"; 

%feature("docstring") simuPOP::population::addIndFrom "

Usage:

    x.addIndFrom(pop)

Details:

    Add all individuals, including ancestors, in pop to the current
    population. Two populations should have the same genotypic
    structures and number of ancestral generations. Subpopulations in
    population pop are kept.

"; 

%feature("docstring") simuPOP::population::addChromFrom "

Usage:

    x.addChromFrom(pop)

Details:

    Add chromosomes in population pop to the current population.
    Population pop should have the same number of individuals as the
    current population in the current and all ancestral generations.
    This function merges genotypes on the new chromosomes from
    population pop individual by individual.

"; 

%feature("docstring") simuPOP::population::addLociFrom "

Usage:

    x.addLociFrom(pop)

Details:

    Add loci from population pop, chromosome by chromosome. Added loci
    will be inserted according to their position. Their position and
    names should not overlap with any locus in the current population.
    Population pop should have the same number of individuals as the
    current population in the current and all ancestral generations.

"; 

%feature("docstring") simuPOP::population::addChrom "

Usage:

    x.addChrom(lociPos, lociNames=[], chromName=\"\", alleleNames=[],
      chromType=Autosome)

Details:

    Add chromosome chromName with given type chromType to a
    population, with loci lociNames inserted at position lociPos.
    lociPos should be ordered. lociNames and chromName should not
    exist in the current population. Allele names could be specified
    for all loci (a list of names) or differently for each locus (a
    nested list of names), using parameter alleleNames. Empty loci
    names will be used if lociNames is not specified.

"; 

%feature("docstring") simuPOP::population::addLoci "

Usage:

    x.addLoci(chrom, pos, lociNames=[], alleleNames=[])

Details:

    Insert loci lociNames at positions pos on chromosome chrom. These
    parameters should be lists of the same length, although names may
    be ignored, in which case empty strings will be assumed. Single-
    value input is allowed for parameter chrom and pos if only one
    locus is added. Alleles at inserted loci are initialized with zero
    alleles. Note that loci have to be added to existing chromosomes.
    If loci on a new chromosome need to be added, function addChrom
    should be used. Optionally, allele names could be specified either
    for all loci (a single list) or each loci (a nested list). This
    function returns indexes of the inserted loci.

"; 

%feature("docstring") simuPOP::population::resize "

Usage:

    x.resize(sizes, propagate=False)

Details:

    Resize population by giving new subpopulation sizes sizes.
    Individuals at the end of some subpopulations will be removed if
    the new subpopulation size is smaller than the old one. New
    individuals will be appended to a subpopulation if the new size is
    larger. Their genotypes will be set to zero (default), or be
    copied from existing individuals if propagate is set to True. More
    specifically, if a subpopulation with 3 individuals is expanded to
    7, the added individuals will copy genotypes from individual 1, 2,
    3, and 1 respectively. Note that this function only resizes the
    current generation.

"; 

%feature("docstring") simuPOP::population::extract "

Usage:

    x.extract(field=\"\", loci=AllAvail, infoFields=AllAvail,
      ancGen=-1, ped=None, pedFields=[])

Details:

    Extract subsets of individuals, loci and/or information fields
    from the current population and create a new population. By
    default, all genotypes and information fields for all individuals
    in all ancestral generations are extracted. If an valid (non-
    empty) information field (field) is given, individuals with
    negative values at this field will be removed and others are put
    into subpopulations specified by this field. The extracted
    population will keep the original subpopulation names if two
    populations have the same number of subpopulations. If a list of
    loci is specified, only genotypes at specified loci are extracted.
    If a list of infoFields is specified, only these information
    fields are extracted. If ancGen is not -1 (default, meaing all
    ancestral generations), only ancGen ancestral generations will be
    extracted. As an advanced feature, field can be information field
    of a pedigree object ped. This allows extraction of individuals
    according to pedigrees identified in a pedigree object. Additional
    information fields from pedFields can be copied to the extracted
    population. This pedigree should have the same number of
    individuals in all generations.

"; 

%feature("docstring") simuPOP::population::removeLoci "

Usage:

    x.removeLoci(loci=[], keep=[])

Details:

    Remove loci (absolute indexes) and genotypes at these loci from
    the current population. Alternatively, a parameter keep can be
    used to specify loci that will not be removed.

"; 

%feature("docstring") simuPOP::population::recodeAlleles "

Usage:

    x.recodeAlleles(alleles, loci=AllAvail, alleleNames=[])

Details:

    Recode alleles at loci (default to all loci in a population) to
    other values according to parameter alleles. This parameter can a
    list of new allele numbers for alleles 0, 1, 2, ... (allele x will
    be recoded to newAlleles[x]) or a Python function. In the latter
    case, each allele and the index of the locus it resides are passed
    to this function. The return value will become the new allele. A
    new list of allele names could be specified for these loci.
    Different sets of names could be specified for each locus if a
    nested list of names are given. This function recode alleles for
    all subpopulations in all ancestral generations.

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
    generation of pop, a the greatness level of all existing ancestral
    generations by one. If ancestralDepth is positive and there are
    already ancestralDepth ancestral generations (see also:
    ancestralGens()), the greatest ancestral generation will be
    discarded. In any case, population pop becomes invalid as all its
    individuals are absorbed by the current population.

"; 

%ignore simuPOP::population::curAncestralGen() const;

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

    x.setIndInfo(values, field, subPop=[])

Details:

    Set information field field (specified by index or name) of all
    individuals (if subPop=[], default), or individuals in a (virtual)
    subpopulation (subPop=sp or (sp, vsp)) to values. values will be
    reused if its length is smaller than the size of the population or
    (virtual) subpopulation.

"; 

%ignore simuPOP::population::infoBegin(UINT idx);

%ignore simuPOP::population::infoEnd(UINT idx);

%feature("docstring") simuPOP::population::indInfo "

Usage:

    x.indInfo(field, subPop=[])

Details:

    Return the values (as a list) of information field field (by index
    or name) of all individuals (if subPop=[], default), or
    individuals in a (virtual) subpopulation (if subPop=sp or (sp,
    vsp)).

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
    ancestral generations. If there exists more than depth ancestral
    generations (if depth > 0), extra ancestral generations are
    removed.

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

%ignore simuPOP::population::gen() const;

%ignore simuPOP::population::setGen(ULONG gen, bool setVar=true);

%feature("docstring") simuPOP::population::vars "

Usage:

    x.vars(subPop=[])

Details:

    return variables of a population as a Python dictionary. If a
    valid subpopulation subPop is specified, a dictionary
    vars()[\"subPop\"][subPop] is returned. A ValueError will be raised
    if key subPop does not exist in vars(), or if key subPop does not
    exist in vars()[\"subPop\"].

"; 

%ignore simuPOP::population::dict(int subPop=-1);

%ignore simuPOP::population::getVars();

%ignore simuPOP::population::setDict(PyObject *dict);

%ignore simuPOP::population::varsAsString() const;

%ignore simuPOP::population::varsFromString(const string &vars);

%feature("docstring") simuPOP::population::evaluate "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::execute "Obsolete or undocumented function."

%feature("docstring") simuPOP::productSplitter "

Details:

    This splitter takes several splitters and take their intersections
    as new VSPs. For example, if the first splitter defines 3 VSPs and
    the second splitter defines 2, 6 VSPs will be defined by splitting
    3 VSPs defined by the first splitter each to two VSPs. This
    splitter is usually used to define finer VSPs from existing VSPs.

"; 

%feature("docstring") simuPOP::productSplitter::productSplitter "

Usage:

    productSplitter(splitters=[], names=[])

Details:

    Create a product splitter using a list of splitters. For example,
    productSplitter([sexSplitter(), affectionSplitter()]) defines four
    VSPs by male unaffected, male affected, female unaffected, and
    female affected individuals. VSP names are usually determined by
    splitters, but can also be specified using parameter names.

"; 

%ignore simuPOP::productSplitter::productSplitter(const productSplitter &rhs);

%feature("docstring") simuPOP::productSplitter::~productSplitter "

Usage:

    x.~productSplitter()

"; 

%feature("docstring") simuPOP::productSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::productSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::productSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs defined by this splitter, which is the
    sum of the number of VSPs of all combined splitters.

"; 

%ignore simuPOP::productSplitter::contains(const population &pop, ULONG ind, vspID vsp);

%ignore simuPOP::productSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::productSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::productSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of a VSP vsp, which is the names of indivdual VSPs
    separated by a comma, unless a new set of names is specified for
    each VSP.

"; 

%feature("docstring") simuPOP::proportionSplitter "

Details:

    This splitter divides subpopulations into several VSPs by
    proportion.

"; 

%feature("docstring") simuPOP::proportionSplitter::proportionSplitter "

Usage:

    proportionSplitter(proportions=[], names=[])

Details:

    Create a splitter that divides subpopulations by proportions,
    which should be a list of float numbers (between 0 and 1) that add
    up to 1. A default set of names are given to each VSP unless a new
    set of names is given by parameter names.

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

%ignore simuPOP::proportionSplitter::contains(const population &pop, ULONG ind, vspID vsp);

%ignore simuPOP::proportionSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::proportionSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::proportionSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of VSP vsp, which is \"Prop p\" where
    p=propotions[vsp]. A user specified name will be returned if
    specified.

"; 

%feature("docstring") simuPOP::pyEval "

Function form:

    PyEval

Details:

    A pyEval operator evaluates a Python expression in a population's
    local namespace when it is applied to this population. The result
    is written to an output specified by parameter output.

"; 

%feature("docstring") simuPOP::pyEval::pyEval "

Usage:

    pyEval(expr=\"\", stmts=\"\", exposePop=\"\", output=\">\", begin=0,
      end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=[])

Details:

    Crete a pyEval operator that evaluates a Python expression expr in
    a population's local namespace when it is applied to this
    population. If Python statements stmts is given (a single or
    multi-line string), the statement will be executed before expr. If
    exposePop is set to an non-empty string, the current population
    will be exposed in its own local namespace as a variable with this
    name. This allows the execution of expressions such as
    'pop.individual(0).allele(0)'. The result of expr will be sent to
    an output stream specified by parameter output. The exposed
    population variable will be removed after expr is evaluated.
    Please refer to class baseOperator for other parameters.

Note:

    Although the statements and expressions are evaluated in a
    population's local namespace, they have access to a **global**
    namespace which is the module global namespace. It is therefore
    possible to refer to any module variable in these expressions.
    Such mixed use of local and global variables is, however, strongly
    discouraged.

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

%feature("docstring") simuPOP::pyEval::evaluate "

Usage:

    x.evaluate(pop)

Details:

    Evaluate the expression and optional statements in the local
    namespace of population pop and return its result as a string.

"; 

%feature("docstring") simuPOP::pyEval::apply "

Description:

    Apply the pyEval operator to population pop.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyEval::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::pyExec "

Function form:

    PyExec

Details:

    This operator executes given Python statements in a population's
    local namespace when it is applied to this population.

"; 

%feature("docstring") simuPOP::pyExec::pyExec "

Usage:

    pyExec(stmts=\"\", exposePop=\"\", output=\">\", begin=0, end=-1,
      step=1, at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a pyExec operator that executes statements stmts in a
    population's local namespace when it is applied to this
    population. If exposePop is given, current population will be
    exposed in its local namespace as a variable named by exposePop.
    Although multiple statements can be executed, it is recommended
    that you use this operator to execute short statements and use
    pyOperator for more complex once. Note that exposed population
    variable will be removed after the statements are executed.

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

%feature("docstring") simuPOP::pyExec::description "Obsolete or undocumented function."

%ignore simuPOP::pyFunc;

%feature("docstring") simuPOP::pyFunc::pyFunc "

Usage:

    pyFunc(func)

"; 

%feature("docstring") simuPOP::pyFunc::isValid "

Usage:

    x.isValid()

"; 

%feature("docstring") simuPOP::pyFunc::func "

Usage:

    x.func()

"; 

%feature("docstring") simuPOP::pyIndIterator "

Details:

    this class implements a Python itertor class that can be used to
    iterate through individuals in a (sub)population. If allInds are
    true, visiblility of individuals will not be checked. Note that
    individualIterator *will* iterate through only visible
    individuals, and allInds is only provided when we know in advance
    that all individuals are visible. This is a way to obtain better
    performance in simple cases.  An instance of this class is
    returned by population::individuals() and
    population::individuals(subPop)

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

%feature("docstring") simuPOP::pyMutator "

Function form:

    PyMutate

Details:

    This hybrid mutator accepts a Python function that determines how
    to mutate an allele when an mutation event happens.

"; 

%feature("docstring") simuPOP::pyMutator::pyMutator "

Usage:

    pyMutator(rates=[], loci=AllAvail, func=None, context=0,
      mapIn=[], mapOut=[], output=\">\", begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a hybrid mutator that uses a user-provided function to
    mutate an allele when a mutation event happens. This function
    (parameter func) accepts the allele to be mutated and return a
    mutated allele. If context is specified, the context alleles to
    the left and to the right of the mutated alleles will be passed to
    this function as the second parameter. Invalid context alleles
    (e.g. left allele to the first locus of a chromosome) will be
    marked by -1. The passed, returned and context alleles might be
    changed if parameters mapIn and mapOut are used although allele
    mappings, if needed, are usually handled in func as well. This
    mutator by default applies to all loci unless parameter loci is
    specified. A single mutation rate will be used for all loci if a
    single value of parameter rates is given. Otherwise, a list of
    mutation rates can be specified for each locus in parameter loci.
    Please refer to classes mutator and baseOperator for descriptions
    of other parameters.

"; 

%feature("docstring") simuPOP::pyMutator::clone "

Description:

    deep copy of a pyMutator

Usage:

    x.clone()

"; 

%ignore simuPOP::pyMutator::mutate(AlleleRef allele, UINT locus);

%feature("docstring") simuPOP::pyMutator::description "Obsolete or undocumented function."

%ignore simuPOP::pyObject;

%feature("docstring") simuPOP::pyObject::pyObject "

Usage:

    pyObject(obj)

"; 

%feature("docstring") simuPOP::pyObject::~pyObject "

Usage:

    x.~pyObject()

"; 

%feature("docstring") simuPOP::pyObject::object "

Usage:

    x.object()

"; 

%feature("docstring") simuPOP::pyObject::isValid "

Usage:

    x.isValid()

"; 

%feature("docstring") simuPOP::pyOperator "

Details:

    An operator that calls a user-defined function when it is applied
    to a population (pre- or post-mating) or offsprings (during-
    mating). The function can have have parameters pop when the
    operator is applied pre- or post-mating, pop, off, dad, mom when
    the operator is applied during-mating. An optional parameter can
    be passed if parameter param is given. In the during-mating case,
    parameters pop, dad and mom can be ignored if offspringOnly is set
    to True.

"; 

%feature("docstring") simuPOP::pyOperator::pyOperator "

Usage:

    pyOperator(func, param=None, offspringOnly=False, begin=0,
      end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=[])

Details:

    Create a pure-Python operator that calls a user-defined function
    when it is applied. Depending on parameters stage, param, and
    offspringOnly, the function should have one of the following
    forms:
    *   func(pop) if used pre- or post-mating without param.
    *   func(pop, param) if used pre- or post-mating with param.
    *   func(pop, off, dad, mom) if used during mating with param.
    *   func(pop, off, dad, mom, param) if used during mating with
    param.
    *   func(off) if used during mating, with offspringOnly=True and
    without param.
    *   func(off, param) if used during mating with offspringOnly=True
    and with param. where pop is the population to which the operator
    is applied, off is the offspring of dad and mom, and param is the
    parameter param specified when the operator is created. When this
    operator is applied during mating, it can be used in the ops
    parameter of a mating scheme, or used in the ops parameter of
    simulator.evolve and be applied after an offspring has been
    created. Please refer to the simuPOP user's guide for a detailed
    explanation.  This operator does not support parameters output,
    subPops and infoFields. If certain output is needed, it should be
    handled in the user defined function func. Because the status of
    files used by other operators through parameter output is
    undetermined during evolution, they should not be open or closed
    in this Python operator.

"; 

%feature("docstring") simuPOP::pyOperator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::pyOperator::apply "

Usage:

    x.apply(pop)

Details:

    Apply the pyOperator operator to population pop. Calling this
    function is equivalent to call func with parameter pop and
    optional parameter param.

"; 

%ignore simuPOP::pyOperator::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::pyOperator::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::pyOutput "

Details:

    This operator outputs a given string when it is applied to a
    population.

"; 

%feature("docstring") simuPOP::pyOutput::pyOutput "

Usage:

    pyOutput(msg=\"\", output=\">\", begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Creates a pyOutput operator that outputs a string msg to output
    (default to standard terminal output) when it is applied to a
    population. Please refer to class baseOperator for a detailed
    description of common operator parameters such as stage, begin and
    output.

"; 

%feature("docstring") simuPOP::pyOutput::apply "

Description:

    simply output some info

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyOutput::~pyOutput "

Description:

    destructor

Usage:

    x.~pyOutput()

"; 

%feature("docstring") simuPOP::pyOutput::clone "

Description:

    Deep copy of a pyOutput operator.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyOutput::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::pyParentsChooser "

Details:

    This parents chooser accept a Python generator function that
    repeatedly yields one or two parents, which can be references to
    individual objects or indexes relative to each subpopulation. The
    parent chooser calls the generator function with parental
    population and a subpopulation index for each subpopulation and
    retrieves parents repeatedly using the iterator interface of the
    generator function.  This parent chooser does not support virtual
    subpopulation directly. A ValueError will be raised if this parent
    chooser is applied to a virtual subpopulation. However, because
    virtual subpopulations are defined in the passed parental
    population, it is easy to return parents from a particular virtual
    subpopulation using virtual subpopulation related functions.

"; 

%feature("docstring") simuPOP::pyParentsChooser::pyParentsChooser "

Usage:

    pyParentsChooser(parentsGenerator)

Details:

    Create a Python parent chooser using a Python generator function
    parentsGenerator. This function should accept a population object
    (the parental population) and a subpopulation number and return
    the reference or index (relative to subpopulation) of a parent or
    a pair of parents repeatedly using the iterator interface of the
    generator function.

"; 

%ignore simuPOP::pyParentsChooser::pyParentsChooser(const pyParentsChooser &rhs);

%feature("docstring") simuPOP::pyParentsChooser::clone "

Description:

    Deep copy of a python parent chooser.

Usage:

    x.clone()

"; 

%ignore simuPOP::pyParentsChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::pyParentsChooser::finalize(population &pop, SubPopID sp);

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
    individual.  More specifically, func can be
    *   func(geno) if infoFields has length 0 or 1.
    *   func(geno, fields) when infoFields has more than 1 fields.
    Both parameters should be an list.

"; 

%feature("docstring") simuPOP::pyPenetrance::pyPenetrance "

Description:

    provide locus and penetrance for 11, 12, 13 (in the form of
    dictionary)

Usage:

    pyPenetrance(loci, func, ancGen=-1, begin=0, end=-1, step=1,
      at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

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

%feature("docstring") simuPOP::pyPenetrance::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::pyPopIterator "

Details:

    This class implements a Python itertor class that can be used to
    iterate through populations in a population.

"; 

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

    pyQuanTrait(loci, func, ancGen=-1, begin=0, end=-1, step=1,
      at=[], reps=AllAvail, subPops=AllAvail, infoFields=AllAvail)

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

%ignore simuPOP::pyQuanTrait::pyQuanTrait(const pyQuanTrait &rhs);

%feature("docstring") simuPOP::pyQuanTrait::clone "

Description:

    deep copy of a Python quantitative trait operator

Usage:

    x.clone()

"; 

%ignore simuPOP::pyQuanTrait::qtrait(individual *ind);

%feature("docstring") simuPOP::pyQuanTrait::description "Obsolete or undocumented function."

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
    *   func(geno, gen) if infoFields has length 0 or 1.
    *   func(geno, gen, fields) when infoFields has more than 1
    fields. Values of fields 1, 2, ... will be passed. Both geno and
    fields should be a list.

"; 

%feature("docstring") simuPOP::pySelector::pySelector "

Description:

    create a Python hybrid selector

Usage:

    pySelector(loci, func, begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=AllAvail)

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

%ignore simuPOP::pySelector::pySelector(const pySelector &rhs);

%feature("docstring") simuPOP::pySelector::clone "

Description:

    deep copy of a pySelector

Usage:

    x.clone()

"; 

%ignore simuPOP::pySelector::indFitness(individual *ind, ULONG gen);

%feature("docstring") simuPOP::pySelector::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::pyTagger "

Details:

    A Python tagger takes some information fields from both parents,
    pass them to a user provided Python function and set the offspring
    individual fields with the return values.

"; 

%feature("docstring") simuPOP::pyTagger::pyTagger "

Usage:

    pyTagger(func=None, begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, output=\"\", infoFields=[])

Details:

    Create a hybrid tagger that passes parental information fields
    (parameter infoFields) to an user provided function func and use
    its return values to assign corresponding information fields of
    offspring. If more than one parent are available, maternal values
    are passed after paternal values. For example, if infoFields=['A',
    'B'], the user-defined function should expect an array of size 4,
    with paternal values at fields 'A', 'B', followed by maternal
    values at these fields. The return value of this function should
    be a list, although a single value will be accepted if only one
    information field is specified. This operator ignores parameters
    stage, output and subPops.

"; 

%feature("docstring") simuPOP::pyTagger::clone "

Description:

    deep copy of a pyTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyTagger::description "Obsolete or undocumented function."

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

    quanTrait(ancGen=-1, begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=AllAvail)

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

%feature("docstring") simuPOP::quanTrait::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::randomParentChooser "

Details:

    This parent chooser chooses a parent randomly from a (virtual)
    parental subpopulation. Parents are chosen with or without
    replacement. If parents are chosen with replacement, a parent can
    be selected multiple times. If natural selection is enabled, the
    probability that an individual is chosen is proportional to
    his/her fitness value stored in an information field
    selectionField (default to \"fitness\"). If parents are chosen
    without replacement, a parent can be chosen only once. An
    RuntimeError will be raised if all parents are exhausted.
    Selection is disabled in the without-replacement case.

"; 

%feature("docstring") simuPOP::randomParentChooser::randomParentChooser "

Usage:

    randomParentChooser(replacement=True, selectionField=\"fitness\")

Details:

    Create a random parent chooser that choose parents with or without
    replacement (parameter replacement, default to True). If selection
    is enabled and information field selectionField exists in the
    passed population, the probability that a parent is chosen is
    proportional to his/her fitness value stored in selectionField.

"; 

%feature("docstring") simuPOP::randomParentChooser::clone "

Description:

    Deep copy of a random parent chooser.

Usage:

    x.clone()

"; 

%ignore simuPOP::randomParentChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::randomParentChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::randomParentsChooser "

Details:

    This parent chooser chooses two parents, a male and a female,
    randomly from a (virtual) parental subpopulation. Parents are
    chosen with or without replacement from their respective sex
    group. If parents are chosen with replacement, a parent can be
    selected multiple times. If natural selection is enabled, the
    probability that an individual is chosen is proportional to
    his/her fitness value among all individuals with the same sex.
    Selection will be disabled if specified information field
    selectionField (default to \"fitness\") does not exist.If parents
    are chosen without replacement, a parent can be chosen only once.
    An RuntimeError will be raised if all males or females are
    exhausted. Selection is disabled in the without-replacement case.

"; 

%feature("docstring") simuPOP::randomParentsChooser::randomParentsChooser "

Usage:

    randomParentsChooser(replacement=True, selectionField=\"fitness\")

Details:

    Create a random parents chooser that choose two parents with or
    without replacement (parameter replacement, default to True). If
    selection is enabled and information field selectionField exists
    in the passed population, the probability that a parent is chosen
    is proportional to his/her fitness value stored in selectionField.

"; 

%feature("docstring") simuPOP::randomParentsChooser::clone "

Description:

    Deep copy of a random parents chooser.

Usage:

    x.clone()

"; 

%ignore simuPOP::randomParentsChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::randomParentsChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::rangeSplitter "

Details:

    This class defines a splitter that groups individuals in certain
    ranges into VSPs.

"; 

%feature("docstring") simuPOP::rangeSplitter::rangeSplitter "

Usage:

    rangeSplitter(ranges, names=[])

Details:

    Create a splitter according to a number of individual ranges
    defined in ranges. For example, rangeSplitter(ranges=[[0, 20],
    [40, 50]]) defines two VSPs. The first VSP consists of individuals
    0, 1, ..., 19, and the sceond VSP consists of individuals 40, 41,
    ..., 49. Note that a nested list has to be used even if only one
    range is defined. A default set of names are given to each VSP
    unless a new set of names is given by parameter names.

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

%ignore simuPOP::rangeSplitter::contains(const population &pop, ULONG ind, vspID vsp);

%ignore simuPOP::rangeSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::rangeSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::rangeSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of VSP vsp, which is \"Range [a, b]\" where [a, b]
    is range ranges[vsp]. A user specified name will be returned if
    specified.

"; 

%feature("docstring") simuPOP::recombinator "

Details:

    A genotype transmitter (during-mating operator) that transmits
    parental chromosomes to offspring, subject to recombination and
    gene conversion. This can be used to replace
    mendelianGenoTransmitter and selfingGenoTransmitter. It does not
    work in haplodiploid populations, although a customized genotype
    transmitter that makes uses this operator could be defined. Please
    refer to the simuPOP user's guide or online cookbook for details.
    Recombination could be applied to all adjacent markers or after
    specified loci. Recombination rate between two adjacent markers
    could be specified directly, or calculated using physical distance
    between them. In the latter case, a recombination intensity is
    multiplied by physical distance between markers.  Gene conversion
    is interpreted as double-recombination events. That is to say, if
    a recombination event happens, it has a certain probability (can
    be 1) to become a conversion event, namely triggering another
    recombination event down the chromosome. The length of the
    converted chromosome can be controlled in a number of ways.

Note:

    simuPOP does not assume any unit to loci positions so
    recombination intensity could be explained differntly (e.g. cM/Mb,
    Morgan/Mb) depending on your intepretation of loci positions. For
    example, if basepair is used for loci position, intensity=10^-8
    indicates 10^-8 per basepair, which is equivalent to 10^-2 per Mb
    or 1 cM/Mb. If Mb is used for physical positions, the same
    recombination intensity could be achieved by intensity=0.01.

"; 

%feature("docstring") simuPOP::recombinator::recombinator "

Usage:

    recombinator(rates=[], intensity=-1, loci=AllAvail,
      convMode=NoConversion, output=\"\", begin=0, end=-1, step=1,
      at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a recombinator (a mendelian genotype transmitter with
    recombination and gene conversion) that passes genotypes from
    parents (or a parent in case of self-fertilization) to offspring.
    Recombination happens by default between all adjacent markers but
    can be limited to a given set of loci. Each locus in this list
    specifies a recombination point between the locus and the locus
    immediately before it. Loci that are the first locus on each
    chromosome are ignored.  If a single recombination rate (parameter
    rates) is specified, it will used for all loci (all loci or loci
    specified by parameter loci), regardless of physical distances
    between adjacent loci.  If a list of recombination rates are
    specified in rates, a parameter loci with the same length should
    also be specified. Different recombination rates can then be used
    after these loci (between specified loci and their immediate
    neighbor to the right).  A recombination intensity (intensity) can
    be used to specify recombination rates that are proportional to
    physical distances between adjacent markers. If the physical
    distance between two markers is d, the recombination rate between
    them will be intensity * d. No unit is assume for loci position
    and recombination intensity.  Gene conversion is controlled using
    parameter convMode, which can be
    *   NoConversion: no gene conversion (default).
    *   (NumMarkers, prob, n): With probability prob, convert a fixed
    number (n) of markers if a recombination event happens.
    *   (GeometricDistribution, prob, p): With probability prob,
    convert a random number of markers if a recombination event
    happens. The number of markes converted follows a geometric
    distribution with probability p.
    *   (TractLength, prob, n): With probability prob, convert a
    region of fixed tract length (n) if a recombination event happens.
    The actual number of markers converted depends on loci positions
    of surrounding loci. The starting position of this tract is the
    middle of two adjacent markers. For example, if four loci are
    located at 0, 1, 2, 3 respectively, a conversion event happens
    between 0 and 1, with a tract length 2 will start at 0.5 and end
    at 2.5, covering the second and third loci.
    *   (ExponentialDistribution, prob, p): With probability prob,
    convert a region of random tract length if a recombination event
    happens. The distribution of tract length follows a exponential
    distribution with probability p. The actual number of markers
    converted depends on loci positions of surrounding loci. simuPOP
    uses this probabilistic model of gene conversion because when a
    recombination event happens, it may become a recombination event
    if the if the Holliday junction is resolved/repaired successfully,
    or a conversion event if the junction is not resolved/repaired.
    The probability, however, is more commonly denoted by the ratio of
    conversion to recombination events in the literature. This ratio
    varies greatly from study to study, ranging from 0.1 to 15 (Chen
    et al, Nature Review Genetics, 2007). This translate to
    0.1/0.9~0.1 to 15/16~0.94 of the gene conversion probability.  A
    recombinator usually does not send any output. However, if an
    information field is given (parameter infoFields), this operator
    will treat this information field as an unique ID of parents and
    offspring and output all recombination events in the format of
    offspring_id parent_id starting_ploidy loc1 loc2 ... where
    starting_ploidy indicates which homologous copy genotype
    replication starts from (0 or 1), loc1, loc2 etc are loci after
    which recombination events happens. If there are multiple
    chromosomes on the genome, you will see a lot of (fake)
    recombination events because of independent segregation of
    chromosomes. Such a record will be generated for each set of
    homologous chromosomes so an diploid offspring will have two lines
    of output. Note that individual IDs need to be set (using a
    idTagger operator) before this recombinator is applied.

Note:

    conversion tract length is usually short, and is estimated to be
    between 337 and 456 bp, with overall range between maybe 50 - 2500
    bp. This is usually not enough to convert, for example, two
    adjacent markers from the HapMap dataset.There is no recombination
    between sex chromosomes (Chromosomes X and Y), although
    recombination is possible between pesudoautosomal regions on these
    chromosomes. If such a feature is required, you will have to
    simulate the pesudoautosomal regions as separate chromosomes.

"; 

%feature("docstring") simuPOP::recombinator::clone "

Description:

    deep copy of a recombinator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::recombinator::~recombinator "

Usage:

    x.~recombinator()

"; 

%feature("docstring") simuPOP::recombinator::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::recombinator::initialize "

Usage:

    x.initialize(pop)

Details:

    Initialize a recombinator for the genotypic structure of
    population pop. This function should be called before a
    recombinator is explicitly applied to a population.

"; 

%feature("docstring") simuPOP::recombinator::transmitGenotype "

Usage:

    x.transmitGenotype(parent, offspring, ploidy)

Details:

    This function transmits genotypes from a parent to the ploidy-th
    homologous set of chromosomes of an offspring. It can be used, for
    example, by a customized genotype transmitter to use sex-specific
    recombination rates to transmit parental genotypes to offspring.

"; 

%ignore simuPOP::recombinator::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad, individual *mom);

%feature("docstring") simuPOP::resizeSubPops "

Function form:

    ResizeSubPops

Details:

    This operator resizes subpopulations to specified sizes.
    Individuals are added or removed depending on the new
    subpopulation sizes.

"; 

%feature("docstring") simuPOP::resizeSubPops::resizeSubPops "

Usage:

    resizeSubPops(subPops=AllAvail, sizes=[], proportions=[],
      propagate=True, begin=0, end=-1, step=1, at=[], reps=AllAvail,
      infoFields=[])

Details:

    Resize given subpopulations subPops to new sizes size, or sizes
    proportional to original sizes (parameter proportions). All
    subpopulations will be resized if subPops is not specified. If the
    new size of a subpopulation is smaller than its original size,
    extra individuals will be removed. If the new size is larger, new
    individuals with empty genotype will be inserted, unless parameter
    propagate is set to True (default). In this case, existing
    individuals will be copied sequentially, and repeatedly if needed.
    This operator is by default applied pre-mating (parameter stage).
    Please refer to operator baseOperator for a detailed explanation
    for all parameters.

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

%feature("docstring") simuPOP::resizeSubPops::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::RNG "

Details:

    This random number generator class wraps around a number of random
    number generators from GNU Scientific Library. You can obtain and
    change the RNG used by the current simuPOP module through the
    GetRNG() function, or create a separate random number generator
    and use it in your script.

"; 

%feature("docstring") simuPOP::RNG::RNG "

Usage:

    RNG(name=None, seed=0)

Details:

    Create a RNG object using specified name and seed. If rng is not
    given, environmental variable GSL_RNG_TYPE will be used if it is
    available. Otherwise, generator mt19937 will be used. If seed is
    not given, /dev/urandom, /dev/random, or other system random
    number source will be used to guarantee that random seeds are used
    even if more than one simuPOP sessions are started simultaneously.

"; 

%feature("docstring") simuPOP::RNG::~RNG "

Usage:

    x.~RNG()

"; 

%feature("docstring") simuPOP::RNG::setRNG "

Usage:

    x.setRNG(name=None, seed=0)

Details:

    Use another underlying RNG for the current RNG object. The
    handling of parameters rng and seed is the same as RNG::RNG(name,
    seed).

"; 

%feature("docstring") simuPOP::RNG::setSeed "

Usage:

    x.setSeed(seed=0)

Details:

    Set random seed for this random number generator. If seed is 0,
    method described in setRNG is used.

"; 

%feature("docstring") simuPOP::RNG::name "

Usage:

    x.name()

Details:

    Return the name of the current random number generator.

"; 

%feature("docstring") simuPOP::RNG::seed "

Usage:

    x.seed()

Details:

    Return the seed used to initialize the RNG. This can be used to
    repeat a previous session.

"; 

%ignore simuPOP::RNG::generateRandomSeed();

%feature("docstring") simuPOP::RNG::randUniform "

Usage:

    x.randUniform()

Details:

    Generate a random number following a rng_uniform [0, 1)
    distribution.

"; 

%feature("docstring") simuPOP::RNG::randBit "Obsolete or undocumented function."

%feature("docstring") simuPOP::RNG::randInt "

Usage:

    x.randInt(n)

Details:

    return a random number in the range of [0, 1, 2, ... n-1]

"; 

%feature("docstring") simuPOP::RNG::randNormal "

Usage:

    x.randNormal(mu, sigma)

Details:

    Generate a random number following a normal distribution with mean
    mu and standard deviation sigma.

"; 

%feature("docstring") simuPOP::RNG::randExponential "

Usage:

    x.randExponential(mu)

Details:

    Generate a random number following a exponential distribution with
    parameter mu.

"; 

%feature("docstring") simuPOP::RNG::randGamma "

Usage:

    x.randGamma(a, b)

Details:

    Generate a random number following a gamma distribution with
    parameters a and b.

"; 

%feature("docstring") simuPOP::RNG::randChisq "

Usage:

    x.randChisq(nu)

Details:

    Generate a random number following a Chi-squared distribution with
    nu degrees of freedom.

"; 

%feature("docstring") simuPOP::RNG::randGeometric "

Usage:

    x.randGeometric(p)

Details:

    Generate a random number following a geometric distribution with
    parameter p.

"; 

%feature("docstring") simuPOP::RNG::randBinomial "

Usage:

    x.randBinomial(n, p)

Details:

    Generate a random number following a binomial distribution with
    parameters n and p.

"; 

%feature("docstring") simuPOP::RNG::randPoisson "

Usage:

    x.randPoisson(mu)

Details:

    Generate a random number following a Poisson distribution with
    parameter mu.

"; 

%feature("docstring") simuPOP::RNG::randMultinomial "

Usage:

    x.randMultinomial(N, p)

Details:

    Generate a random number following a multinomial distribution with
    parameters N and p (a list of probabilities).

"; 

%ignore simuPOP::RNG::randomShuffle(T begin, T end) const;

%ignore simuPOP::RNG_func;

%feature("docstring") simuPOP::RNG_func::RNG_func "

Usage:

    RNG_func(rng)

"; 

%feature("docstring") simuPOP::RuntimeError "

Description:

    exception, thrown if a runtime error occurs

"; 

%feature("docstring") simuPOP::RuntimeError::RuntimeError "

Usage:

    RuntimeError(msg)

"; 

%feature("docstring") simuPOP::savePopulation "

Details:

    An operator that save populations to specified files.

"; 

%feature("docstring") simuPOP::savePopulation::savePopulation "

Usage:

    savePopulation(output=\"\", begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create an operator that saves a population to output when it is
    applied to the population. This operator supports all output
    specifications ('', 'filename', 'filename' prefixed by one or more
    '>' characters, and '!expr') but output from different operators
    will always replace existing files (effectively ignore '>'
    specification). Parameter subPops is ignored. Please refer to
    class baseOperator for a detailed description about common
    operator parameters such as stage and begin.

"; 

%feature("docstring") simuPOP::savePopulation::~savePopulation "

Description:

    destructor.

Usage:

    x.~savePopulation()

"; 

%feature("docstring") simuPOP::savePopulation::clone "

Description:

    Deep copy of a savePopulation operator.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::savePopulation::apply "

Description:

    Apply operator to population pop.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::savePopulation::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::selector "

Details:

    Genetic selection is tricky to simulate since there are many
    different fitness values and many different ways to apply
    selection. simuPOP employs an 'ability-to-mate' approach. Namely,
    the probability that an individual will be chosen for mating is
    proportional to its fitness value. More specifically,
    *   PreMating selectors assign fitness values to each individual,
    and mark part or all subpopulations as under selection.
    *   during sexless mating (e.g. binomialSelection mating scheme),
    individuals are chosen at probabilities that are proportional to
    their fitness values. If there are $ N $ individuals with fitness
    values $ f_{i},i=1,...,N $, individual $ i $ will have probability
    $ \\frac{f_{i}}{\\sum_{j}f_{j}} $ to be chosen and passed to the
    next generation.
    *   during randomMating, males and females are separated. They are
    chosen from their respective groups in the same manner as
    binomialSelection and mate.
    All of the selection operators, when applied, will set an
    information field fitness (configurable) and then mark part or all
    subpopulations as under selection. (You can use different
    selectors to simulate various selection intensities for different
    subpopulations). Then, a 'selector-aware' mating scheme can select
    individuals according to their fitness information fields. This
    implies that
    *   only mating schemes can actually select individuals.
    *   a selector has to be a PreMating operator. This is not a
    problem when you use the operator form of the selector since its
    default stage is PreMating. However, if you use the function form
    of the selector in a pyOperator, make sure to set the stage of
    pyOperator to PreMating.

Note:

    You can not apply two selectors to the same subpopulation, because
    only one fitness value is allowed for each individual.

"; 

%feature("docstring") simuPOP::selector::selector "

Usage:

    selector(begin=0, end=-1, step=1, at=[], reps=AllAvail,
      subPops=AllAvail, infoFields=AllAvail)

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

%feature("docstring") simuPOP::selector::applyDuringMating "

Description:

    apply the operator during mating.

Usage:

    x.applyDuringMating(pop, offspring, dad=None, mom=None)

"; 

%feature("docstring") simuPOP::selector::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::selfingGenoTransmitter "

Details:

    A genotype transmitter (during-mating operator) that transmits
    parental genotype of a parent through self-fertilization. That is
    to say, the offspring genotype is formed according to Mendel's
    laws, only that a parent serves as both maternal and paternal
    parents.

"; 

%feature("docstring") simuPOP::selfingGenoTransmitter::selfingGenoTransmitter "

Usage:

    selfingGenoTransmitter(output=\"\", begin=0, end=-1, step=1,
      at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a self-fertilization genotype transmitter that transmits
    genotypes of a parent to an offspring through self-fertilization.
    Cutsomized chromosomes are not handled. Parameters subPops and
    infoFields are ignored.

"; 

%feature("docstring") simuPOP::selfingGenoTransmitter::clone "

Description:

    Deep copy of a selfing genotype transmitter.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::selfingGenoTransmitter::description "Obsolete or undocumented function."

%ignore simuPOP::selfingGenoTransmitter::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::sequentialParentChooser "

Details:

    This parent chooser chooses a parent from a parental (virtual)
    subpopulation sequentially. Sex and selection is not considered.
    If the last parent is reached, this parent chooser will restart
    from the beginning of the (virtual) subpopulation.

"; 

%feature("docstring") simuPOP::sequentialParentChooser::sequentialParentChooser "

Usage:

    sequentialParentChooser()

Details:

    Create a parent chooser that chooses a parent from a parental
    (virtual) subpopulation sequentially.

"; 

%feature("docstring") simuPOP::sequentialParentChooser::clone "

Description:

    Deep copy of a sequential parent chooser.

Usage:

    x.clone()

"; 

%ignore simuPOP::sequentialParentChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::sequentialParentChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::sequentialParentsChooser "

Details:

    This parent chooser chooses two parents (a father and a mother)
    sequentially from their respective sex groups. Selection is not
    considered. If all fathers (or mothers) are exhausted, this parent
    chooser will choose fathers (or mothers) from the beginning of the
    (virtual) subpopulation again.

"; 

%feature("docstring") simuPOP::sequentialParentsChooser::sequentialParentsChooser "

Usage:

    sequentialParentsChooser()

Details:

    Create a parent chooser that chooses two parents sequentially from
    a parental (virtual) subpopulation.

"; 

%feature("docstring") simuPOP::sequentialParentsChooser::clone "

Description:

    Deep copy of a sequential parents chooser.

Usage:

    x.clone()

"; 

%ignore simuPOP::sequentialParentsChooser::initialize(population &pop, SubPopID sp);

%ignore simuPOP::sequentialParentsChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::setAncestralDepth "

Details:

    This operator sets the number of ancestral generations to keep
    during the evolution of a population. This is usually used to
    start storing ancestral generations at the end of an evolutionary
    process. A typical usage is setAncestralDepth(1, at=-1) which will
    cause the parental generation of the present population to be
    stored at the last generation of an evolutionary process.

"; 

%feature("docstring") simuPOP::setAncestralDepth::setAncestralDepth "

Usage:

    setAncestralDepth(depth, output=\">\", begin=0, end=-1, step=1,
      at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a setAncestralDepth operator that sets the ancestral depth
    of an population. It basically calls the
    population.setAncestralDepth member function of a population.

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

    apply the setAncestralDepth operator to population pop.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::setAncestralDepth::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::sexSplitter "

Details:

    This splitter defines two VSPs by individual sex. The first VSP
    consists of all male individuals and the second VSP consists of
    all females in a subpopulation.

"; 

%feature("docstring") simuPOP::sexSplitter::sexSplitter "

Usage:

    sexSplitter(names=[])

Details:

    Create a sex splitter that defines male and female VSPs. These
    VSPs are named Male and Female unless a new set of names are
    specified by parameter names.

"; 

%feature("docstring") simuPOP::sexSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::sexSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const;

%feature("docstring") simuPOP::sexSplitter::numVirtualSubPop "

Description:

    Return 2.

Usage:

    x.numVirtualSubPop()

"; 

%ignore simuPOP::sexSplitter::contains(const population &pop, ULONG ind, vspID vsp);

%ignore simuPOP::sexSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::sexSplitter::deactivate(population &pop, SubPopID sp);

%feature("docstring") simuPOP::sexSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return \"Male\" if vsp=0 and \"Female\" otherwise, unless a new set of
    names are specified.

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

%ignore simuPOP::SharedVariables::setIntDefDictVar(const string &name, const intDict &val);

%ignore simuPOP::SharedVariables::setTupleDefDictVar(const string &name, const tupleDict &val);

%ignore simuPOP::SharedVariables::getVarAsBool(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsInt(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsDouble(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsString(const string &name, bool nameError=true);

%feature("docstring") simuPOP::SharedVariables::dict "

Usage:

    x.dict()

"; 

%ignore simuPOP::SharedVariables::asString() const;

%feature("docstring") simuPOP::SharedVariables::fromString "

Usage:

    x.fromString(vars)

"; 

%ignore simuPOP::simpleStmt;

%feature("docstring") simuPOP::simpleStmt::simpleStmt "

Usage:

    simpleStmt(stmt)

"; 

%feature("docstring") simuPOP::simpleStmt::var "

Usage:

    x.var()

"; 

%feature("docstring") simuPOP::simpleStmt::operation "

Usage:

    x.operation()

"; 

%feature("docstring") simuPOP::simpleStmt::value "

Usage:

    x.value()

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
    For convenience, member functions are provided to set virtual
    splitter, add information field and set ancestral depth to all
    populations in a simulator.

"; 

%feature("docstring") simuPOP::simulator::simulator "

Usage:

    simulator(pop, matingScheme, rep=1)

Details:

    Create a simulator with rep replicates of population pop.
    Population pop will be copied rep times (default to 1), while
    keeping the passed population intact. A mating scheme matingScheme
    will be used to evolve these populations.

Note:

    Population pop is copied to a simulator so the input population
    will be kept untouched.

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

%feature("docstring") simuPOP::simulator::add "

Usage:

    x.add(pop)

Details:

    Add a population pop to the end of an existing simulator. This
    creates an cloned copy of pop in the simulator so the evolution of
    the simulator will not change pop.

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

%feature("docstring") simuPOP::simulator::evolve "

Usage:

    x.evolve(initOps=[], preOps=[], duringOps=[], postOps=[],
      finalOps=[], gen=-1, dryrun=False)

Details:

    Evolve all populations gen generations, subject to several lists
    of operators which are applied at different stages of an
    evolutionary process. Operators initOps are applied to all
    populations (subject to applicability restrictions of the
    operators, imposed by the rep parameter of these operators) before
    evolution. They are used to initialize populations before
    evolution. Operators finalOps are applied to all populations after
    the evolution.  Operators preOps, duringOps and postOps are
    applied during the life cycle of each generation. These operators
    can be applied at all or some of the generations, to all or some
    of the evolving populations, depending the begin, end, step, at
    and reps parameters of these operators. These operators are
    applied in the order at which they are specified. Populations in a
    simulator are evolved one by one. At each generation, operators
    preOps are applied to the parental generations. A mating scheme is
    then used to populate an offspring generation. For each offspring,
    his or her sex is determined before during-mating operators of the
    mating scheme are used to transmit parental genotypes. During-
    mating operators specified in parameters duringOps are applied
    afterwards. An offspring will be discarded if any of the during-
    mating operator fails (return False). After an offspring
    generation is successfully generated and becomes the current
    generation, operators postOps are applied to the offspring
    generation. If any of the preOps and postOps fails (return False),
    the evolution of a population will be stopped. The generation
    number of a population is increased by one if an offspring
    generation has been successfully populated even if a post-during
    operator fails.  Parameter gen can be set to a positive number,
    which is the number of generations to evolve. Because a simulator
    always starts at the beginning of a generation g (e.g. 0), a
    simulator will stop at the beginning (instead of the end) of
    generation g + gen (e.g. gen). If gen is negative (default), the
    evolution will continue indefinitely, until all replicates are
    stopped by operators that return False at some point (these
    operators are called terminators). At the end of the evolution,
    the generations that each replicates have evolved are returned.
    Note that finalOps are applied to all applicable population,
    including those that have stopped before others.  The last
    parameter dryrun, if set to True, will print a description of this
    evolutionary process. It can help you understand what exactly will
    happen at each generation during an evolutionary process.

"; 

%ignore simuPOP::simulator::apply(const opList &ops, bool dryrun=false);

%feature("docstring") simuPOP::simulator::setMatingScheme "

Usage:

    x.setMatingScheme(matingScheme)

Details:

    Set a new mating scheme matingScheme to a simulator.

"; 

%feature("docstring") simuPOP::simulator::vars "

Usage:

    x.vars(rep, subPop=[])

Details:

    Return the local namespace of the rep-th population, equivalent to
    x.population(rep).vars(subPop).

"; 

%feature("docstring") simuPOP::simulator::__cmp__ "

Description:

    Note that mating schemes are not tested.

Usage:

    x.__cmp__(rhs)

"; 

%feature("docstring") simuPOP::simulator::save "

Usage:

    x.save(filename)

Details:

    Save a simulator to file filename, which can be loaded by a global
    function LoadSimulator.

"; 

%ignore simuPOP::simulator::load(string filename);

%feature("docstring") simuPOP::smmMutator "

Function form:

    SmmMutate

Details:

    A stepwise mutation model treats alleles at a locus as the number
    of tandem repeats of microsatellite or minisatellite markers. When
    a mutation event happens, the number of repeats (allele) either
    increase or decrease. A standard stepwise mutation model increases
    of decreases an allele by 1 with equal probability. More complex
    models (generalized stepwise mutation model) are also allowed.
    Note that an allele cannot be mutated beyond boundaries (0 and
    maximum allowed allele).

"; 

%feature("docstring") simuPOP::smmMutator::smmMutator "

Usage:

    smmMutator(rates=[], loci=AllAvail, incProb=0.5, maxAllele=0,
      mutStep=[], mapIn=[], mapOut=[], output=\">\", begin=0, end=-1,
      step=1, at=[], reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a stepwise mutation mutator that mutates an allele by
    increasing or decreasing it. This mutator by default applies to
    all loci unless parameter loci is specified. A single mutation
    rate will be used for all loci if a single value of parameter
    rates is given. Otherwise, a list of mutation rates can be
    specified for each locus in parameter loci.  When a mutation event
    happens, this operator increases or decreases an allele by mutStep
    steps. Acceptable input of parameter mutStep include
    *   A number: This is the default mode with default value 1.
    *   (GeometricDistribution, p): The number of steps follows a a
    geometric distribution with parameter p.
    *   A Python function: This user defined function accepts the
    allele being mutated and return the steps to mutate. The mutation
    process is usually neutral in the sense that mutating up and down
    is equally likely. You can adjust parameter incProb to change this
    behavior.  If you need to use other generalized stepwise mutation
    models, you can implement them using a pyMutator. If performance
    becomes a concern, I may add them to this operator if provided
    with a reliable reference.

"; 

%feature("docstring") simuPOP::smmMutator::~smmMutator "

Usage:

    x.~smmMutator()

"; 

%ignore simuPOP::smmMutator::mutate(AlleleRef allele, UINT locus);

%feature("docstring") simuPOP::smmMutator::clone "

Description:

    deep copy of a smmMutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::smmMutator::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::splitSubPops "

Function form:

    SplitSubPops

Details:

    Split a given list of subpopulations according to either sizes of
    the resulting subpopulations, proportion of individuals, or an
    information field. The resulting subpopulations will have the same
    name as the original subpopulation.

"; 

%feature("docstring") simuPOP::splitSubPops::splitSubPops "

Usage:

    splitSubPops(subPops=AllAvail, sizes=[], proportions=[],
      names=[], randomize=True, begin=0, end=-1, step=1, at=[],
      reps=AllAvail, infoFields=[])

Details:

    Split a list of subpopulations subPops into finer subpopulations.
    A single subpopulation is acceptable but virtual subpopulations
    are not allowed. All subpopulations will be split if subPops is
    not specified.  The subpopulations can be split in three ways:
    *   If parameter sizes is given, each subpopulation will be split
    into subpopulations with given size. The sizes should add up to
    the size of all orignal subpopulations.
    *   If parameter proportions is given, each subpopulation will be
    split into subpopulations with corresponding proportion of
    individuals. proportions should add up to 1.
    *   If an information field is given (parameter infoFields),
    individuals having the same value at this information field will
    be grouped into a subpopulation. The number of resulting
    subpopulations is determined by the number of distinct values at
    this information field. If parameter randomize is True (default),
    individuals will be randomized before a subpopulation is split.
    This is designed to remove artificial order of individuals
    introduced by, for example, some non- random mating schemes. Note
    that, however, the original individual order is not guaranteed
    even if this parameter is set to False.  Unless the last
    subpopulation is split, the indexes of existing subpopulations
    will be changed. If a subpopulation has a name, this name will
    become the name for all subpopulations separated from this
    subpopulation. Optionally, you can assign names to the new
    subpopulations using a list of names specified in parameter names.
    Because the same set of names will be used for all subpopulations,
    this parameter is not recommended when multiple subpopulations are
    split.  This operator is by default applied pre-mating (parameter
    stage). Please refer to operator baseOperator for a detailed
    explanation for all parameters.

Note:

    Unlike operator migrator, this operator does not require an
    information field such as migrate_to.

"; 

%feature("docstring") simuPOP::splitSubPops::~splitSubPops "

Description:

    destructor

Usage:

    x.~splitSubPops()

"; 

%feature("docstring") simuPOP::splitSubPops::clone "

Description:

    deep copy of a splitSubPops operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::splitSubPops::apply "

Description:

    apply a splitSubPops operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::splitSubPops::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::stat "

Function form:

    Stat

Details:

    Operator stat calculates various statistics of the population
    being applied and sets variables in its local namespace. Other
    operators or functions can retrieve results from or evalulate
    expressions in this local namespace after stat is applied.

"; 

%feature("docstring") simuPOP::stat::stat "

Usage:

    stat(popSize=False, numOfMale=False, numOfAffected=False,
      alleleFreq=[], heteroFreq=[], homoFreq=[], genoFreq=[],
      haploFreq=[], sumOfInfo=[], meanOfInfo=[], varOfInfo=[],
      maxOfInfo=[], minOfInfo=[], LD=[], association=[],
      neutrality=[], structure=[], HWE=[], vars=AllAvail, suffix=\"\",
      output=\"\", begin=0, end=-1, step=1, at=[], reps=AllAvail,
      subPops=AllAvail, infoFields=[])

Details:

    Create a stat operator that calculates specified statistics of a
    population when it is applied to this population. This operator is
    by default applied after mating (parameter stage) and can be
    applied to specified replicates (parameter rep) at specified
    generations (parameter begin, end, step, and at). This operator
    does not produce any output (ignore parameter output) after
    statistics are calculated. Instead, it stores results in the local
    namespace of the population being applied. Other operators can
    retrieve these variables or evalulate expression directly in this
    local namespace. Please refer to operator baseOperator for a
    detailed explanation of these common operator parameters.  stat
    supports parameter subPops. It usually calculate the same set of
    statistics for all subpopulations (subPops=subPopList()). If a
    list of (virtual) subpopulations are specified, statistics for
    only specified subpopulations will be calculated. However,
    different statistics treat this parameter differently and it is
    very important to check its reference before you use subPops for
    any statistics.  Calculated statistics are saved as variables in a
    population's local namespace. These variables can be numbers,
    lists or dictionaries and can be retrieved using functions
    population.vars() or population.dvars(). A special default
    dictionary (defdict) is used for dictionaries whose keys are
    determined dynamically. Accessing elements of such a dictionary
    with an invalid key will yield value 0 instead of a KeyError. If
    the same variables are calculated for one or more (virtual)
    subpopulation, the variables are stored in
    vars()['subPop'][sp]['var'] where sp is a subpopulation ID (sp) or
    a tuple of virtual subpopulation ID ((sp, vsp)).
    population.vars(sp) and population.dvars(sp) provide shortcuts to
    these variables.  Operator stat outputs a number of most useful
    variables for each type of statistic. For example, alleleFreq
    calculates both allele counts and allele frequencies and it by
    default sets variable alleleFreq (dvars().alleleFreq) for all or
    specified subpopulations. If this does not fit your need, you can
    use parameter vars to output additional parameters, or limit the
    output of existing parameters. More specifically, for this
    particular statistic, the available variables are 'alleleFreq',
    'alleleNum', 'alleleFreq_sp' ('alleleFreq' in each subpopulation),
    and 'alleleNum_sp' ('alleleNum' in each subpopulation). You can
    set vars=['alleleNum_sp'] to output only subpopulation specific
    allele count. An optional suffix (parameter suffix) can be used to
    append a suffix to default parameter names. This parameter can be
    used, for example, to calculate and store the same statistics for
    different subpopulations (e.g. pairwise Fst).  Operator stat
    supports the following statistics:  popSize: If popSize=True,
    number of individuals in all or specified subpopulations
    (parameter subPops) will be set to the following variables:
    *   popSize (default): Number of individuals in all or specified
    subpopulations. Because subPops does not have to cover all
    individuals, it may not be the actual population size.
    *   popSize_sp: Size of (virtual) subpopulation sp.
    *   subPopSize (default): A list of (virtual) subpopulation sizes.
    This variable is easier to use than accessing popSize from each
    (virtual) subpopulation.numOfMale: If numOfMale=True, number of
    male individuals in all or specified (virtual) subpopulations will
    be set to the following variables:
    *   numOfMale (default): Total number of male individuals in all
    or specified (virtual) subpopulations.
    *   numOfMale (default): Total number of female individuals in all
    or specified (virtual) subpopulations.
    *   propOfMale: Proportion of male individuals.
    *   propOfFemale: Proportion of female individuals.
    *   numOfMale_sp: Number of male individuals in each (virtual)
    subpopulation.
    *   numOfFemale_sp: Number of female individuals in each (virtual)
    subpopulation.
    *   propOfMale_sp: Proportion of male individuals in each
    (virtual) subpopulation.
    *   propOfFemale_sp: Proportion of female individuals in each
    (virtual) subpopulation.numOfAffected: If numOfAffected=True,
    number of affected individuals in all or specified (virtual)
    subpopulations will be set to the following variables:
    *   numOfAffected (default): Total number of affected individuals
    in all or specified (virtual) subpopulations.
    *   numOfAffected (default): Total number of unaffected
    individuals in all or specified (virtual) subpopulations.
    *   propOfAffected: Proportion of affected individuals.
    *   propOfUnaffected: Proportion of unaffected individuals.
    *   numOfAffected_sp: Number of affected individuals in each
    (virtual) subpopulation.
    *   numOfUnaffected_sp: Number of unaffected individuals in each
    (virtual) subpopulation.
    *   propOfAffected_sp: Proportion of affected individuals in each
    (virtual) subpopulation.
    *   propOfUnaffected_sp: Proportion of unaffected individuals in
    each (virtual) subpopulation.alleleFreq: This parameter accepts a
    list of loci (by indexes), at which allele frequencies will be
    calculated. This statistic outputs the following variables, all of
    which are dictionary (with loci indexes as keys) of default
    dictionaries (with alleles as keys). For example,
    alleleFreq[loc][a] returns 0 if allele a does not exist.
    *   alleleFreq (default): alleleFreq[loc][a] is the frequency of
    allele a at locus for all or specified (virtual) subpopulations.
    *   alleleNum (default): alleleNum[loc][a] is the number of allele
    a at locus for all or specified (virtual) subpopulations.
    *   alleleFreq_sp: Allele frequency in each (virtual)
    subpopulation.
    *   alleleNum_sp: Allele count in each (virtual)
    subpopulation.heteroFreq and homoFreq: These parameters accept a
    list of loci (by indexes), at which the number and frequency of
    homozygotes and/or heterozygotes will be calculated. These
    statistics are only available for diploid populations. The
    following variables will be outputted:
    *   heteroFreq (default for parameter heteroFreq): A dictionary of
    proportion of heterozygotes in all or specified (virtual)
    subpopulations, with loci indexes as dictionary keys.
    *   homoFreq (default for parameter homoFreq): A dictionary of
    proportion of homozygotes in all or specified (virtual)
    subpopulations.
    *   heteroNum: A dictionary of number of heterozygotes in all or
    specified (virtual) subpopulations.
    *   homoNum: A dictionary of number of homozygotes in all or
    specified (virtual) subpopulations.
    *   heteroFreq_sp: A dictionary of proportion of heterozygotes in
    each (virtual) subpopulation.
    *   homoFreq_sp: A dictionary of proportion of homozygotes in each
    (virtual) subpopulation.
    *   heteroNum_sp: A dictionary of number of heterozygotes in each
    (virtual) subpopulation.
    *   homoNum_sp: A dictionary of number of homozygotes in each
    (virtual) subpopulation.genoFreq: This parameter accept a list of
    loci (by index) at which number and frequency of all genotypes are
    outputed as a dictionary (indexed by loci indexes) of default
    dictionaries (indexed by tuples of possible indexes). This
    statistic is available for all population types with genotype
    defined as ordered alleles at a locus. The length of genotype
    equals the number of homologous copies of chromosomes (ploidy) of
    a population. Genotypes for males or females on sex chromosomes or
    in haplodiploid populations will have different length. Because
    genotypes are ordered, (1, 0) and (0, 1) (two possible genotypes
    in a diploid population) are considered as different genotypes.
    This statistic outputs the following variables:
    *   genoFreq (default): A dictionary (by loci indexes) of default
    dictionaries (by genotype) of genotype frequencies. For example,
    genoFreq[1][(1, 0)] is the frequency of genotype (1, 0) at locus
    1.
    *   genoNum (default): A dictionary of default dictionaries of
    genotype counts of all or specified (virtual) subpopulations.
    *   genoFreq_sp: genotype frequency in each specified (virtual)
    subpopulation.
    *   genoFreq_sp: genotype count in each specified (virtual)
    subpopulation.haploFreq: This parameter accepts one or more lists
    of loci (by index) at which number and frequency of haplotypes are
    outputted as default dictionaries. [(1,2)] can be abbreviated to
    (1,2). For example, using parameter haploFreq=(1,2,4), all
    haplotypes at loci 1, 2 and 4 are counted. This statistic saves
    results to dictionary (with loci index as keys) of default
    dictionaries (with haplotypes as keys) such as
    haploFreq[(1,2,4)][(1,1,0)] (frequency of haplotype (1,1,0) at
    loci (1,2,3)). This statistic works for all population types.
    Number of haplotypes for each individual equals to his/her ploidy
    number. Haplodiploid populations are supported in the sense that
    the second homologous copy of the haplotype is not counted for
    male individuals. This statistic outputs the following variables:
    *   haploFreq (default): A dictionary (with tuples of loci indexes
    as keys) of default dictionaries of haplotype frequencies. For
    example, haploFreq[(0, 1)][(1,1)] records the frequency of
    haplotype (1,1) at loci (0, 1) in all or specified (virtual)
    subpopulations.
    *   haploNum (default): A dictionary of default dictionaries of
    haplotype counts in all or specified (virtual) subpopulations.
    *   haploFreq_sp: Halptype frequencies in each (virtual)
    subpopulation.
    *   haploNum_sp: Halptype count in each (virtual)
    subpopulation.sumOfinfo, meanOfInfo, varOfInfo, maxOfInfo and
    minOfInfo: Each of these five parameters accepts a list of
    information fields. For each information field, the sum, mean,
    variance, maximum or minimal (depending on the specified
    parameter(s)) of this information field at iddividuals in all or
    specified (virtual) subpopulations will be calculated. The results
    will be put into the following population variables:
    *   sumOfInfo (default for sumOfInfo): A dictionary of the sum of
    specified information fields of individuals in all or specified
    (virtual) subpopulations. This dictionary is indexed by names of
    information fields.
    *   meanOfInfo (default for meanOfInfo): A dictionary of the mean
    of information fields of all individuals.
    *   varOfInfo (default for varOfInfo): A dictionary of the sample
    variance of information fields of all individuals.
    *   maxOfInfo (default for maxOfInfo): A dictionary of the maximum
    value of information fields of all individuals.
    *   minOfInfo (default for minOfInfo): A dictionary of the minimal
    value of information fields of all individuals.
    *   sumOfInfo_sp: A dictionary of the sum of information fields of
    individuals in each subpopulation.
    *   meanOfInfo_sp: A dictionary of the mean of information fields
    of individuals in each subpopulation.
    *   varOfInfo_sp: A dictionary of the sample variance of
    information fields of individuals in each subpopulation.
    *   maxOfInfo_sp: A dictionary of the maximum value of information
    fields of individuals in each subpopulation.
    *   minOfInfo_sp: A dictionary of the minimal value of information
    fields of individuals in each subpopulation.LD: Parameter LD
    accepts one or a list of loci pairs (e.g. LD=[[0,1], [2,3]]) with
    optional primary alleles at both loci (e.g. LD=[0,1,0,0]). For
    each pair of loci, this operator calculates linkage disequilibrium
    and optional association statistics between two loci. When primary
    alleles are specified, signed linkage disequilibrium values are
    calculated with non-primary alleles are combined. Otherwise,
    absolute values of diallelic measures are combined to yield
    positive measure of LD. Association measures are calculated from a
    m by n contigency of haplotype counts (m=n=2 if primary alleles
    are specified). Please refer to the simuPOP user's guide for
    detailed information. This statistic sets the following variables:
    *   LD (default) Basic LD measure for haplotypes in all or
    specified (virtual) subpopulations. Signed if primary alleles are
    specified.
    *   LD_prime (default) Lewontin's D' measure for haplotypes in all
    or specified (virtual) subpopulations. Signed if primary alleles
    are specified.
    *   R2 (default) Correlation LD measure for haplotypes in all or
    specified (virtual) subpopulations.
    *   LD_ChiSq ChiSq statistics for a contigency table with
    frequencies of haplotypes in all or specified (virtual)
    subpopulations.
    *   LD_ChiSq_p Single side p-value for the ChiSq statistic.
    Degrees of freedom is determined by number of alleles at both loci
    and the specification of primary alleles.
    *   CramerV Normalized ChiSq statistics.
    *   LD_sp Basic LD measure for haplotypes in each (virtual)
    subpopulation.
    *   LD_prime_sp Lewontin's D' measure for haplotypes in each
    (virtual) subpopulation.
    *   R2_sp R2 measure for haplotypes in each (virtual)
    subpopulation.
    *   LD_ChiSq_sp ChiSq statistics for each (virtual) subpopulation.
    *   LD_ChiSq_p_sp p value for the ChiSq statistics for each
    (virtual) subpopulation.
    *   CramerV_sp Cramer V statistics for each (virtual)
    subpopulation.association: Parameter association accepts a list of
    loci. At each locus, one or more statistical tests will be
    performed to test association between this locus and individual
    affection status. Currently, simuPOP provides the following tests:
    *   An allele-based Chi-square test using alleles counts. This
    test can be applied to loci with more than two alleles, and to
    haploid populations.
    *   A genotype-based Chi-square test using genotype counts. This
    test can be applied to loci with more than two alleles (more than
    3 genotypes) in diploid populations. aA and Aa are considered to
    be the same genotype.
    *   A genotype-based Cochran-Armitage trend test. This test can
    only be applied to diallelic loci in diploid populations. A
    codominant model is assumed. This statistic sets the following
    variables:
    *   Allele_ChiSq A dictionary of allele-based Chi-Square
    statistics for each locus, using cases and controls in all or
    specified (virtual) subpopulations.
    *   Allele_ChiSq_p (default) A dictionary of p-values of the
    corresponding Chi-square statistics.
    *   Geno_ChiSq A dictionary of genotype-based Chi-Square
    statistics for each locus, using cases and controls in all or
    specified (virtual) subpopulations.
    *   Geno_ChiSq_p A dictionary of p-values of the corresponding
    genotype-based Chi-square test.
    *   Armitage_p A dictionary of p-values of the Cochran-Armitage
    tests, using cases and controls in all or specified (virtual)
    subpopulations.
    *   Allele_ChiSq_sp A dictionary of allele-based Chi-Square
    statistics for each locus, using cases and controls from each
    subpopulation.
    *   Allele_ChiSq_p_sp A dictionary of p-values of allele-based
    Chi-square tests, using cases and controls from each (virtual)
    subpopulation.
    *   Geno_ChiSq_sp A dictionary of genotype-based Chi-Square tests
    for each locus, using cases and controls from each subpopulation.
    *   Geno_ChiSq_p_sp A dictionary of p-values of genotype-based
    Chi-Square tests, using cases and controls from each
    subpopulation.
    *   Armitage_p_sp A dictionary of p-values of the Cochran-
    Armitage tests, using cases and controls from each
    subpopulation.neutrality: This parameter performs neutrality tests
    (detection of natural selection) on specified loci. It currently
    only outputs Pi, which is the average number of pairwise
    difference between loci. This statistic outputs the following
    variables:
    *   Pi Mean pairwise difference between all sequences from all or
    specified (virtual) subpopulations.
    *   Pi_sp Mean paiewise difference between all sequences in each
    (virtual) subpopulation.structure: Parameter structure accepts a
    list of loci at which statistics that measure population structure
    are calculated. This parameter currently supports the following
    statistics:
    *   Weir and Cockerham's Fst (1984). This is the most widely used
    estimator of Wright's fixation index and can be used to measure
    population differentiation. However, this method is designed to
    estimate Fst from samples of larger populations and might not be
    appropriate for the calculation of Fst of large populations.
    *   Nei's Gst (1973). The Gst estimator is another estimator for
    Wright's fixation index but it is extended for multi-allele (more
    than two alleles) and multi-loci cases. This statistics should be
    used if you would like to obtain a true Fst value of a large
    population.
    *   F_st (default) The WC84 Fst statistic estimated for all
    specified loci.
    *   F_is The WC84 Fis statistic estimated for all specified loci.
    *   F_it The WC84 Fit statistic estimated for all specified loci.
    *   f_st A dictionary of locus level WC84 Fst values.
    *   f_is A dictionary of locus level WC84 Fis values.
    *   f_it A dictionary of locus level WC84 Fit values.
    *   G_st Nei's Gst statistic estimated for all specified loci.
    *   g_st A dictionary of Nei's Gst statistic estimated for each
    locus.HWE: Parameter HWE accepts a list of loci at which exact
    two-side tests for Hardy-Weinberg equilibrium will be performed.
    This statistic is only available for diallelic loci in diploid
    populations. It outputs the following variables:
    *   HWE (default) A dictionary of p-values of HWE tests using
    genotypes in all or specified (virtual) subpopulations.
    *   HWE_sp A dictionary of p-values of HWS tests using genotypes
    in each (virtual) subpopulation.

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

%feature("docstring") simuPOP::stat::description "Obsolete or undocumented function."

%ignore simuPOP::statAlleleFreq;

%feature("docstring") simuPOP::statAlleleFreq::statAlleleFreq "

Usage:

    statAlleleFreq(loci, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statAlleleFreq::~statAlleleFreq "

Description:

    destructor, nested vectors have to be cleared manually

Usage:

    x.~statAlleleFreq()

"; 

%feature("docstring") simuPOP::statAlleleFreq::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statAssociation;

%feature("docstring") simuPOP::statAssociation::statAssociation "

Usage:

    statAssociation(loci, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statAssociation::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statGenoFreq;

%feature("docstring") simuPOP::statGenoFreq::statGenoFreq "

Usage:

    statGenoFreq(genoFreq, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statGenoFreq::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statHaploFreq;

%feature("docstring") simuPOP::statHaploFreq::statHaploFreq "

Usage:

    statHaploFreq(haploFreq, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statHaploFreq::~statHaploFreq "

Usage:

    x.~statHaploFreq()

"; 

%feature("docstring") simuPOP::statHaploFreq::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statHeteroFreq;

%feature("docstring") simuPOP::statHeteroFreq::statHeteroFreq "

Usage:

    statHeteroFreq(heteroFreq, homoFreq, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statHeteroFreq::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statHWE;

%feature("docstring") simuPOP::statHWE::statHWE "

Usage:

    statHWE(loci, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statHWE::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statInfo;

%feature("docstring") simuPOP::statInfo::statInfo "

Usage:

    statInfo(sumOfInfo, meanOfInfo, varOfInfo, maxOfInfo, minOfInfo,
      subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statInfo::~statInfo "

Usage:

    x.~statInfo()

"; 

%feature("docstring") simuPOP::statInfo::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statLD;

%feature("docstring") simuPOP::statLD::statLD "

Usage:

    statLD(LD, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statLD::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statNeutrality;

%feature("docstring") simuPOP::statNeutrality::statNeutrality "

Usage:

    statNeutrality(loci, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statNeutrality::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statNumOfAffected;

%feature("docstring") simuPOP::statNumOfAffected::statNumOfAffected "

Usage:

    statNumOfAffected(numOfAffected, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statNumOfAffected::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statNumOfMale;

%feature("docstring") simuPOP::statNumOfMale::statNumOfMale "

Usage:

    statNumOfMale(numOfMale, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statNumOfMale::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statPopSize;

%feature("docstring") simuPOP::statPopSize::statPopSize "

Usage:

    statPopSize(popSize, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statPopSize::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statStructure;

%feature("docstring") simuPOP::statStructure::statStructure "

Usage:

    statStructure(Fst, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statStructure::apply "

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

%ignore simuPOP::StreamProvider::StreamProvider(const string &output, const pyFunc &func);

%feature("docstring") simuPOP::StreamProvider::~StreamProvider "

Usage:

    x.~StreamProvider()

"; 

%ignore simuPOP::StreamProvider::noOutput();

%ignore simuPOP::StreamProvider::getOstream(PyObject *dict=NULL, bool readable=false);

%ignore simuPOP::StreamProvider::closeOstream();

%feature("docstring") simuPOP::stringFunc "

"; 

%feature("docstring") simuPOP::stringFunc::stringFunc "

Usage:

    stringFunc(value)

"; 

%ignore simuPOP::stringFunc::value() const;

%ignore simuPOP::stringFunc::func() const;

%ignore simuPOP::stringFunc::empty() const;

%feature("docstring") simuPOP::stringList "

"; 

%feature("docstring") simuPOP::stringList::stringList "

Usage:

    stringList(str=None)

"; 

%ignore simuPOP::stringList::stringList(const string &str);

%ignore simuPOP::stringList::obtainFrom(const stringList &items, const char *allowedItems[], const char *defaultItems[]);

%ignore simuPOP::stringList::allAvail() const;

%ignore simuPOP::stringList::empty() const;

%ignore simuPOP::stringList::contains(const string &str) const;

%ignore simuPOP::stringList::push_back(const string &str);

%ignore simuPOP::stringList::elems() const;

%feature("docstring") simuPOP::stringMatrix "

"; 

%feature("docstring") simuPOP::stringMatrix::stringMatrix "

Usage:

    stringMatrix(str=None)

"; 

%ignore simuPOP::stringMatrix::empty() const;

%ignore simuPOP::stringMatrix::elems() const;

%feature("docstring") simuPOP::subPopList "

Details:

    A class to specify (virtual) subpopulation list. Using a dedicated
    class allows users to specify a single subpopulation, or a list of
    (virutal) subpoulations easily.

"; 

%feature("docstring") simuPOP::subPopList::subPopList "

Usage:

    subPopList(obj=None)

"; 

%ignore simuPOP::subPopList::subPopList(const vectorvsp &subPops);

%ignore simuPOP::subPopList::allAvail();

%ignore simuPOP::subPopList::empty() const;

%ignore simuPOP::subPopList::size() const;

%feature("docstring") simuPOP::subPopList::__len__ "

Usage:

    x.__len__()

"; 

%ignore simuPOP::subPopList::push_back(const vspID subPop);

%ignore simuPOP::subPopList::contains(const vspID subPop) const;

%ignore simuPOP::subPopList::begin() const;

%ignore simuPOP::subPopList::end() const;

%ignore simuPOP::subPopList::useSubPopsFrom(const population &pop);

%feature("docstring") simuPOP::summaryTagger "

Details:

    A summary tagger summarize values of one or more parental
    information field to another information field of an offspring. If
    mating is sexual, two sets of parental values will be involved.

"; 

%feature("docstring") simuPOP::summaryTagger::summaryTagger "

Usage:

    summaryTagger(mode=Mean, begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, output=\"\", infoFields=[])

Details:

    Creates a summary tagger that summarize values of one or more
    parental information field (infoFields[:-1]) to an offspring
    information field (infoFields[-1]). A parameter mode specifies how
    to pass summarize parental values. More specifically,
    *   mode=Mean Passing the average of values.
    *   mode=Maximum Passing the maximum value of values.
    *   mode=Minumum Passing the minimum value of values.
    *   mode=Summation Passing the sum of values.
    *   mode=Multiplication Passing the multiplication of values. This
    operator does not support parameter subPops and does not output
    any information.

"; 

%feature("docstring") simuPOP::summaryTagger::~summaryTagger "

Usage:

    x.~summaryTagger()

"; 

%feature("docstring") simuPOP::summaryTagger::description "Obsolete or undocumented function."

%ignore simuPOP::summaryTagger::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::summaryTagger::clone "

Description:

    deep copy of a summaryTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::SystemError "

Description:

    exception, thrown if system error occurs

"; 

%feature("docstring") simuPOP::SystemError::SystemError "

Usage:

    SystemError(msg)

"; 

%feature("docstring") simuPOP::terminateIf "

Details:

    This operator evaluates an expression in a population's local
    namespace and terminate the evolution of this population, or the
    whole simulator, if the return value of this expression is True.
    Termination caused by an operator will stop the execution of all
    operators after it. The generation at which the population is
    terminated will be counted in the evolved generations (return
    value from simulator::evolve) if termination happens after mating.

"; 

%feature("docstring") simuPOP::terminateIf::terminateIf "

Usage:

    terminateIf(condition=\"\", stopAll=False, message=\"\", output=\"\",
      begin=0, end=-1, step=1, at=[], reps=AllAvail, subPops=AllAvail,
      infoFields=[])

Details:

    Create a terminator with an expression condition, which will be
    evalulated in a population's local namespace when the operator is
    applied to this population. If the return value of condition is
    True, the evolution of the population will be terminated. If
    stopAll is set to True, the evolution of all replicates of the
    simulator will be terminated. If this operator is allowed to write
    to an output (default to \"\"), the generation number, proceeded
    with an optional message.

"; 

%feature("docstring") simuPOP::terminateIf::clone "

Description:

    deep copy of a terminateIf terminator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::terminateIf::description "Obsolete or undocumented function."

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

Details:

    This operator, when called, output the difference between current
    and the last called clock time. This can be used to estimate
    execution time of each generation. Similar information can also be
    obtained from turnOnDebug(\"DBG_PROFILE\"), but this operator has
    the advantage of measuring the duration between several
    generations by setting step parameter.

"; 

%feature("docstring") simuPOP::ticToc::ticToc "

Usage:

    ticToc(output=\">\", begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    Create a ticToc operator that outputs the elapsed since the last
    time it was applied, and the overall time since it was created.

"; 

%feature("docstring") simuPOP::ticToc::~ticToc "

Description:

    destructor

Usage:

    x.~ticToc()

"; 

%feature("docstring") simuPOP::ticToc::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::ticToc::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::ticToc::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::turnOffDebug "

Function form:

    TurnOffDebug

Details:

    Turn off certain debug information. Please refer to operator
    turnOnDebug for detailed usages.

"; 

%feature("docstring") simuPOP::turnOffDebug::turnOffDebug "

Usage:

    turnOffDebug(code=\"DBG_ALL\", begin=0, end=-1, step=1, at=[],
      reps=AllAvail, subPops=AllAvail, infoFields=[])

Details:

    create a turnOffDebug operator that turns off debug information
    code when it is applied to a population.

"; 

%feature("docstring") simuPOP::turnOffDebug::~turnOffDebug "

Description:

    destructor

Usage:

    x.~turnOffDebug()

"; 

%feature("docstring") simuPOP::turnOffDebug::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::turnOffDebug::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::turnOffDebug::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::turnOnDebug "

Function form:

    TurnOnDebug

Details:

    Turn on debug. There are several ways to turn on debug information
    for non-optimized modules, namely
    *   set environment variable SIMUDEBUG.
    *   use simuOpt.setOptions(debug) function.
    *   use function TurnOnDebug
    *   use the turnOnDebug operator The advantage of using an
    operator is that you can turn on debug at given generations.

"; 

%feature("docstring") simuPOP::turnOnDebug::turnOnDebug "

Usage:

    turnOnDebug(code, begin=0, end=-1, step=1, at=[], reps=AllAvail,
      subPops=AllAvail, infoFields=[])

Details:

    create a turnOnDebug operator that turns on debug information code
    when it is applied to a population.

"; 

%feature("docstring") simuPOP::turnOnDebug::~turnOnDebug "

Description:

    destructor

Usage:

    x.~turnOnDebug()

"; 

%feature("docstring") simuPOP::turnOnDebug::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::turnOnDebug::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::turnOnDebug::description "Obsolete or undocumented function."

%feature("docstring") simuPOP::uintList "

"; 

%feature("docstring") simuPOP::uintList::uintList "

Usage:

    uintList(obj=None)

"; 

%ignore simuPOP::uintList::uintList(const vectoru &values);

%ignore simuPOP::uintList::allAvail() const;

%ignore simuPOP::uintList::elems() const;

%feature("docstring") simuPOP::uintListFunc "

"; 

%feature("docstring") simuPOP::uintListFunc::uintListFunc "

Usage:

    uintListFunc(values=[])

"; 

%ignore simuPOP::uintListFunc::func() const;

%ignore simuPOP::uintListFunc::empty() const;

%feature("docstring") simuPOP::uintString "

"; 

%feature("docstring") simuPOP::uintString::uintString "

Usage:

    uintString(value)

"; 

%ignore simuPOP::uintString::empty() const;

%ignore simuPOP::uintString::name() const;

%ignore simuPOP::uintString::value() const;

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

%feature("docstring") simuPOP::vspID::valid "

Usage:

    x.valid()

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
    can be applied to individuals within specified VSPs.  Each VSP has
    a name. A default name is determined by each splitter but you can
    also assign a name to each VSP. The name of a VSP can be retrieved
    by function population.subPopName.  Only one VSP splitter can be
    assigned to a population, which defined VSPs for all its
    subpopulations. If different splitters are needed for different
    subpopulations, a combinedSplitter should be used.

"; 

%feature("docstring") simuPOP::vspSplitter::vspSplitter "

Usage:

    vspSplitter(names=[])

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

%ignore simuPOP::vspSplitter::contains(const population &pop, ULONG ind, vspID vsp);

%ignore simuPOP::vspSplitter::activate(population &pop, SubPopID subPop, SubPopID virtualSubPop, IterationType type);

%ignore simuPOP::vspSplitter::deactivate(population &pop, SubPopID subPop);

%feature("docstring") simuPOP::vspSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of VSP vsp (an index between 0 and
    numVirtualSubPop()).

"; 

%feature("docstring") simuPOP::weightedSampler "

"; 

%feature("docstring") simuPOP::weightedSampler::weightedSampler "

Usage:

    weightedSampler(rng)

"; 

%feature("docstring") simuPOP::weightedSampler::~weightedSampler "

Usage:

    x.~weightedSampler()

"; 

%feature("docstring") simuPOP::weightedSampler::set "

Description:

    set parameters

Usage:

    x.set(weight)

"; 

%feature("docstring") simuPOP::weightedSampler::set "

Description:

    set parameters for the second case.

Usage:

    x.set(weight, N)

"; 

%feature("docstring") simuPOP::weightedSampler::get "

Usage:

    x.get()

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

%feature("docstring") simuPOP::DebugCodes "

Usage:

    DebugCodes()

Details:

    Return names of all debug codes

"; 

%feature("docstring") simuPOP::TurnOnDebug "

Usage:

    TurnOnDebug(code=\"\")

Details:

    Set debug code code. More than one code could be specified using a
    comma separated string. Name of available codes are available from
    DebugCodes.

"; 

%feature("docstring") simuPOP::TurnOffDebug "

Usage:

    TurnOffDebug(code=\"DBG_ALL\")

Details:

    Turn off debug code code. More than one code could be specified
    using a comma separated string. Default to turn off all debug
    codes.

"; 

%ignore simuPOP::debug(DBG_CODE code);

%ignore simuPOP::simuPOP_kbhit();

%ignore simuPOP::simuPOP_getch();

%ignore simuPOP::PyObj_As_Bool(PyObject *obj, bool &val);

%ignore simuPOP::PyObj_As_Int(PyObject *obj, long int &val);

%ignore simuPOP::PyObj_As_Double(PyObject *obj, double &val);

%ignore simuPOP::PyObj_As_String(PyObject *obj, string &val);

%ignore simuPOP::PyObj_As_Array(PyObject *obj, vectorf &val);

%ignore simuPOP::PyObj_As_IntArray(PyObject *obj, vectori &val);

%ignore simuPOP::PyObj_Is_IntNumArray(PyObject *obj);

%ignore simuPOP::PyObj_Is_DoubleNumArray(PyObject *obj);

%ignore simuPOP::PyObj_Is_AlleleNumArray(PyObject *obj);

%ignore simuPOP::Double_Vec_As_NumArray(vectorf::iterator begin, vectorf::iterator end);

%ignore simuPOP::Int_Vec_As_NumArray(vectori::iterator begin, vectori::iterator end);

%ignore simuPOP::Allele_Vec_As_NumArray(GenoIterator begin, GenoIterator end);

%ignore simuPOP::NumArray_Size(PyObject *obj);

%ignore simuPOP::NumArray_Data(PyObject *obj);

%ignore simuPOP::mainVars();

%ignore simuPOP::moduleVars();

%ignore simuPOP::pyPopObj(void *p);

%ignore simuPOP::pyIndObj(void *p);

%ignore simuPOP::pyIndPointer(PyObject *p);

%ignore simuPOP::ostreamManager();

%feature("docstring") simuPOP::CloseOutput "

Usage:

    CloseOutput(output=\"\")

Details:

    Output files specified by '>' are closed immediately after they
    are written. Those specified by '>>' and '>>>' are closed by a
    simulator after simulator.evolve(). However, these files will be
    kept open if the operators are applied directly to a population
    using the operators' function form. In this case, function
    closeOutput can be used to close a specific file output, and close
    all unclosed files if output is unspecified. An exception will be
    raised if output does not exist or it has already been closed.

"; 

%ignore simuPOP::chisqTest(const vector< vectoru > &table, double &chisq, double &chisq_p);

%ignore simuPOP::armitageTrendTest(const vector< vectoru > &table, const vectorf &weight);

%ignore simuPOP::hweTest(const vectoru &cnt);

%ignore simuPOP::propToCount(const vectorf &prop, ULONG N, vectoru &count);

%feature("docstring") simuPOP::GetRNG "

Description:

    return the currently used random number generator

Usage:

    GetRNG()

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

%feature("docstring") simuPOP::ModuleInfo "

Usage:

    ModuleInfo()

Details:

    Return a dictionary with information regarding the currently
    loaded simuPOP module. This dictionary has the following keys:
    *   revision: revision number.
    *   version: simuPOP version string.
    *   optimized: Is this module optimized (True or False).
    *   alleleType: Allele type of the module (short, long or binary).
    *   maxAllele: the maximum allowed allele state, which is 1 for
    binary modules, 255 for short modules and 65535 for long modules.
    *   compiler: the compiler that compiles this module.
    *   date: date on which this module is compiled.
    *   python: version of python.
    *   platform: platform of the module.
    *   maxNumSubPop: maximum number of subpopulations.
    *   maxIndex: maximum index size (limits population size * total
    number of marker).
    *   debug: A list of effective debugging codes.

"; 

%ignore simuPOP::initialize();

%ignore simuPOP::cnull();

%ignore std::pow3(unsigned n);

%feature("docstring") simuPOP::population::dvars "

Usage:

    x.dvars(subPop=[])

Details:

    Return a wrapper of Python dictionary returned by vars(subPop) so
    that dictionary keys can be accessed as attributes.

"; 

%feature("docstring") simuPOP::simulator::dvars "

Usage:

    x.dvars(rep, subPop=[])

Details:

    Return a wrapper of Python dictionary returned by vars(rep,
    subPop) so that dictionary keys can be accessed as attributes.

"; 

