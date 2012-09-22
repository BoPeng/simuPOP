%feature("docstring") simuPOP::AffectionSplitter "

Details:

    This class defines two VSPs according individual affection status.
    The first VSP consists of unaffected invidiauls and the second VSP
    consists of affected ones.

"; 

%feature("docstring") simuPOP::AffectionSplitter::AffectionSplitter "

Usage:

    AffectionSplitter(names=[])

Details:

    Create a splitter that defined two VSPs by affection status.These
    VSPs are named Unaffected and Affected unless a new set of names
    are specified by parameter names.

"; 

%feature("docstring") simuPOP::AffectionSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::AffectionSplitter::size(const Population &pop, size_t subPop, size_t virtualSubPop) const;

%feature("docstring") simuPOP::AffectionSplitter::numVirtualSubPop "

Description:

    Return 2.

Usage:

    x.numVirtualSubPop()

"; 

%ignore simuPOP::AffectionSplitter::contains(const Population &pop, size_t ind, vspID vsp) const;

%ignore simuPOP::AffectionSplitter::activate(const Population &pop, size_t subPop, size_t virtualSubPop);

%feature("docstring") simuPOP::AffectionSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return \"Unaffected\" if vsp=0 and \"Affected\" if vsp=1, unless a new
    set of names are specified.

"; 

%feature("docstring") simuPOP::BaseMutator "

Details:

    Class mutator is the base class of all mutators. It handles all
    the work of picking an allele at specified loci from certain
    (virtual) subpopulation with certain probability, and calling a
    derived mutator to mutate the allele. Alleles can be changed
    before and after mutation if existing allele numbers do not match
    those of a mutation model.

"; 

%feature("docstring") simuPOP::BaseMutator::BaseMutator "

Usage:

    BaseMutator(rates=[], loci=ALL_AVAIL, mapIn=[], mapOut=[],
      context=0, output=\">\", begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=\"ind_id\",
      lineageMode=FROM_INFO)

Details:

    A mutator mutates alleles from one state to another with given
    probability. This base mutator does not perform any mutation but
    it defines common behaviors of all mutators.  By default, a
    mutator mutates all alleles in all populations of a simulator at
    all generations. A number of parameters can be used to restrict
    mutations to certain generations (parameters begin, end, step and
    at), replicate populations (parameter rep), (virtual)
    subpopulations (parameter subPops) and loci (parameter loci).
    Parameter loci can be a list of loci indexes, names or ALL_AVAIL.
    Please refer to class BaseOperator for a detailed explanation of
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
    AcgtMutator assumes alleles 0, 1, 2 and 3 as nucleotides A, C, G
    and T. Using a mutator that is incompatible with your simulation
    will certainly yield erroneous results.  If your simulation
    assumes different alleles with a mutation model, you can map an
    allele to the allele used in the model and map the mutated allele
    back. This is achieved using a mapIn list with its i-th item being
    the corresponding allele of real allele i, and a mapOut list with
    its i-th item being the real allele of allele i assumed in the
    model. For example mapIn=[0, 0, 1] and mapOut=[1, 2] would allow
    the use of a SNPMutator to mutate between alleles 1 and 2, instead
    of 0 and 1. Parameters mapIn and mapOut also accept a user-defined
    Python function that returns a corresponding allele for a given
    allele. This allows easier mapping between a large number of
    alleles and advanced models such as random emission of alleles.
    If a valid information field is specified for parameter infoFields
    (default to ind_id) for modules with lineage allele type, the
    lineage of the mutated alleles will be the ID (stored in the first
    field of infoFields) of individuals that harbor the mutated
    alleles if lineageMode is set to FROM_INFO (default). If
    lineageMode is set to FROM_INFO_SIGNED, the IDs will be assigned a
    sign depending on the ploidy the mutation happens (1 for ploidy 0,
    -1 for ploidy 1, etc). The lineage information will be transmitted
    along with the alleles so this feature allows you to track the
    source of mutants during evolution.A  Some mutation models are
    context dependent. Namely, how an allele mutates will depend on
    its adjecent alleles. Whereas most simuPOP mutators are context
    independent, some of them accept a parameter context which is the
    number of alleles to the left and right of the mutated allele. For
    example context=1 will make the alleles to the immediate left and
    right to a mutated allele available to a mutator. These alleles
    will be mapped in if parameter mapIn is defined. How exactly a
    mutator makes use of these information is mutator dependent.

"; 

%feature("docstring") simuPOP::BaseMutator::~BaseMutator "

Description:

    destructor

Usage:

    x.~BaseMutator()

"; 

%feature("docstring") simuPOP::BaseMutator::clone "Obsolete or undocumented function."

%ignore simuPOP::BaseMutator::setRate(const vectorf &rates, const lociList &loci);

%ignore simuPOP::BaseMutator::mutRate(size_t loc) const;

%ignore simuPOP::BaseMutator::mutate(Allele, size_t) const;

%ignore simuPOP::BaseMutator::fillContext(const Population &pop, IndAlleleIterator ptr, size_t locus) const;

%ignore simuPOP::BaseMutator::setContext(size_t context);

%ignore simuPOP::BaseMutator::context() const;

%feature("docstring") simuPOP::BaseMutator::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::BaseOperator "

Details:

    Operators are objects that act on populations. They can be applied
    to populations directly using their function forms, but they are
    usually managed and applied by a simulator. In the latter case,
    operators are passed to the evolve function of a simulator, and
    are applied repeatedly during the evolution of the simulator.  The
    BaseOperator class is the base class for all operators. It defines
    a common user interface that specifies at which generations, at
    which stage of a life cycle, to which populations and
    subpopulations an operator is applied. These are achieved by a
    common set of parameters such as begin, end, step, at, stage for
    all operators. Note that a specific operator does not have to
    honor all these parameters. For example, a Recombinator can only
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

%feature("docstring") simuPOP::BaseOperator::BaseOperator "

Usage:

    BaseOperator(output, begin, end, step, at, reps, subPops,
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
                    value ALL_AVAIL is interpreted as all replicates
                    in a simulator. Negative indexes such as -1 (last
                    replicate) is acceptable. rep=idx can be used as a
                    shortcut for rep=[idx].
    subPops:        A list of applicable (virtual) subpopulations,
                    such as subPops=[sp1, sp2, (sp2, vsp1)].
                    subPops=[sp1] can be simplied as subPops=sp1.
                    Negative indexes are not supported. A common
                    default value (ALL_AVAIL) of this parameter
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

%feature("docstring") simuPOP::BaseOperator::~BaseOperator "

Description:

    destroy an operator

Usage:

    x.~BaseOperator()

"; 

%feature("docstring") simuPOP::BaseOperator::clone "

Usage:

    x.clone()

Details:

    Return a cloned copy of an operator. This function is available to
    all operators.

"; 

%ignore simuPOP::BaseOperator::isActive(size_t rep, ssize_t gen, ssize_t end, const vector< bool > &activeRep, bool repOnly=false) const;

%ignore simuPOP::BaseOperator::isActive(ssize_t rep, ssize_t gen) const;

%ignore simuPOP::BaseOperator::infoSize() const;

%ignore simuPOP::BaseOperator::infoField(size_t idx) const;

%ignore simuPOP::BaseOperator::infoFields() const;

%feature("docstring") simuPOP::BaseOperator::apply "

Usage:

    x.apply(pop)

Details:

    Apply an operator to population pop directly, without checking its
    applicability.

"; 

%ignore simuPOP::BaseOperator::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%ignore simuPOP::BaseOperator::getOstream(PyObject *dict=NULL, bool readable=false) const;

%ignore simuPOP::BaseOperator::closeOstream() const;

%ignore simuPOP::BaseOperator::applicability(bool subPops=true, bool gen=true) const;

%feature("docstring") simuPOP::BaseOperator::describe "Obsolete or undocumented function."

%ignore simuPOP::BaseOperator::noOutput() const;

%ignore simuPOP::BaseOperator::applicableSubPops(const Population &pop) const;

%ignore simuPOP::BaseOperator::applicableToAllOffspring() const;

%ignore simuPOP::BaseOperator::applicableToOffspring(const Population &pop, RawIndIterator offspring) const;

%ignore simuPOP::BaseOperator::parallelizable() const;

%ignore simuPOP::BaseOperator::initialize(const Individual &ind) const;

%ignore simuPOP::BaseOperator::initializeIfNeeded(const Individual &ind) const;

%feature("docstring") simuPOP::BasePenetrance "

Details:

    A penetrance model models the probability that an individual has a
    certain disease provided that he or she has certain genetic
    (genotype) and environmental (information field) riske factors. A
    penetrance operator calculates this probability according to
    provided information and set his or her affection status randomly.
    For example, an individual will have probability 0.8 to be
    affected if the penetrance is 0.8. This class is the base class to
    all penetrance operators and defines a common interface for all
    penetrance operators.  A penetrance operator can be applied at any
    stage of an evolutionary cycle. If it is applied before or after
    mating, it will set affection status of all parents and offspring,
    respectively. If it is applied during mating, it will set the
    affection status of each offspring. You can also apply a
    penetrance operator to an individual using its applyToIndividual
    member function.  By default, a penetrance operator assigns
    affection status of individuals but does not save the actual
    penetrance value. However, if an information field is specified,
    penetrance values will be saved to this field for future analysis.
    When a penetrance operator is applied to a population, it is only
    applied to the current generation. You can, however, use parameter
    ancGens to set affection status for all ancestral generations
    (ALL_AVAIL), or individuals in specified generations if a list of
    ancestral generations is specified. Note that this parameter is
    ignored if the operator is applied during mating.

"; 

%feature("docstring") simuPOP::BasePenetrance::BasePenetrance "

Usage:

    BasePenetrance(ancGens=UNSPECIFIED, begin=0, end=-1, step=1,
      at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a base penetrance operator. This operator assign individual
    affection status in the present generation (default). If ALL_AVAIL
    or a list of ancestral generations are spcified in parameter
    ancGens, individuals in specified ancestral generations will be
    processed. A penetrance operator can be applied to specified
    (virtual) subpopulations (parameter subPops) and replicates
    (parameter reps). If an informatio field is given, penetrance
    value will be stored in this information field of each individual.

"; 

%feature("docstring") simuPOP::BasePenetrance::~BasePenetrance "

Description:

    destructor

Usage:

    x.~BasePenetrance()

"; 

%feature("docstring") simuPOP::BasePenetrance::clone "Obsolete or undocumented function."

%ignore simuPOP::BasePenetrance::penet(Population *, RawIndIterator) const;

%feature("docstring") simuPOP::BasePenetrance::apply "

Description:

    set penetrance to all individuals and record penetrance if
    requested

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::BasePenetrance::applyToIndividual "

Usage:

    x.applyToIndividual(ind, pop=None)

Details:

    Apply the penetrance operator to a single individual ind and set
    his or her affection status. A population reference can be passed
    if the penetrance model depends on population properties such as
    generation number. This function returns the affection status.

"; 

%ignore simuPOP::BasePenetrance::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::BasePenetrance::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::BaseQuanTrait "

Details:

    A quantitative trait in simuPOP is simply an information field. A
    quantitative trait model simply assigns values to one or more
    information fields (called trait fields) of each individual
    according to its genetic (genotype) and environmental (information
    field) factors. It can be applied at any stage of an evolutionary
    cycle. If a quantitative trait operator is applied before or after
    mating, it will set the trait fields of all parents and offspring.
    If it is applied during mating, it will set the trait fields of
    each offspring.  When a quantitative trait operator is applied to
    a population, it is only applied to the current generation. You
    can, however, use parameter ancGen=-1 to set the trait field of
    all ancestral generations, or a generation index to apply to only
    ancestral generation younger than ancGen. Note that this parameter
    is ignored if the operator is applied during mating.

"; 

%feature("docstring") simuPOP::BaseQuanTrait::BaseQuanTrait "

Usage:

    BaseQuanTrait(ancGens=UNSPECIFIED, begin=0, end=-1, step=1,
      at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a base quantitative trait operator. This operator assigns
    one or more quantitative traits to trait fields in the present
    generation (default). If ALL_AVAIL or a list of ancestral
    generations are specified, this operator will be applied to
    individuals in these generations as well. A quantitative trait
    operator can be applied to specified (virtual) subpopulations
    (parameter subPops) and replicates (parameter reps).

"; 

%feature("docstring") simuPOP::BaseQuanTrait::~BaseQuanTrait "

Description:

    destructor

Usage:

    x.~BaseQuanTrait()

"; 

%feature("docstring") simuPOP::BaseQuanTrait::clone "Obsolete or undocumented function."

%ignore simuPOP::BaseQuanTrait::qtrait(Individual *, size_t, vectorf &) const;

%feature("docstring") simuPOP::BaseQuanTrait::apply "

Description:

    set qtrait to all individual

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::BaseQuanTrait::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::BaseQuanTrait::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::BaseSelector "

Details:

    This class is the base class to all selectors, namely operators
    that perform natural selection. It defines a common interface for
    all selectors.  A selector can be applied before mating or during
    mating. If a selector is applied to one or more (virtual)
    subpopulations of a parental population before mating, it sets
    individual fitness values to all involved parents to an
    information field (default to fitness). When a mating scheme that
    supports natural selection is applied to the parental population,
    it will select parents with probabilities that are proportional to
    individual fitness stored in an information field (default to
    fitness). Individual fitness is considered relative fitness and
    can be any non-negative number. This simple process has some
    implications that can lead to advanced usages of natural selection
    in simuPOP:
    *   It is up to the mating scheme how to handle individual
    fitness. Some mating schemes do not support natural selection at
    all.
    *   A mating scheme performs natural selection according to
    fitness values stored in an information field. It does not care
    how these values are set. For example, fitness values can be
    inherited from a parent using a tagging operator, or set directly
    using a Python operator.
    *   A mating scheme can treat any information field as fitness
    field. If an specified information field does not exist, or if all
    individuals have the same fitness values (e.g. 0), the mating
    scheme selects parents randomly.
    *   Multiple selectors can be applied to the same parental
    generation. individual fitness is determined by the last fitness
    value it is assigned.
    *   A selection operator can be applied to virtual subpopulations
    and set fitness values only to part of the individuals.
    *   individuals with zero fitness in a subpopulation with anyone
    having a positive fitness value will not be selected to produce
    offspring. This can sometimes lead to unexpected behaviors. For
    example, if you only assign fitness value to part of the
    individuals in a subpopulation, the rest of them will be
    effectively discarded. If you migrate individuals with valid
    fitness values to a subpopulation with all individuals having zero
    fitness, the migrants will be the only mating parents.
    *   It is possible to assign multiple fitness values to different
    information fields so that different homogeneous mating schemes
    can react to different fitness schemes when they are used in a
    heterogeneous mating scheme.
    *   You can apply a selector to the offspring generation using the
    postOps parameter of Simulator.evolve, these fitness values will
    be used when the offspring generation becomes parental generation
    in the next generation. Alternatively, a selector can be used as a
    during mating operator. In this case, it caculates fitness value
    for each offspring which will be treated as absolute fitness,
    namely the probability for each offspring to survive. This process
    uses the fact that an individual will be discarded when any of the
    during mating operators returns False. It is important to remember
    that:
    *   individual fitness needs to be between 0 and 1 in this case.
    *   Fitness values are not stored so the population does not need
    an information field fitness.
    *   This method applies natural selection to offspring instead of
    parents. These two implementation can be identical or different
    depending on the mating scheme used.
    *   Seleting offspring is less efficient than the selecting
    parents, especially when fitness values are low.
    *   Parameter subPops are applied to the offspring population and
    is used to judge if an operator should be applied. It thus does
    not make sense to apply a selector to a virtual subpopulation with
    affected individuals.

"; 

%feature("docstring") simuPOP::BaseSelector::BaseSelector "

Usage:

    BaseSelector(output=\"\", begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)

Details:

    Create a base selector object. This operator should not be created
    directly.

"; 

%feature("docstring") simuPOP::BaseSelector::~BaseSelector "

Description:

    destructor

Usage:

    x.~BaseSelector()

"; 

%feature("docstring") simuPOP::BaseSelector::clone "Obsolete or undocumented function."

%ignore simuPOP::BaseSelector::indFitness(Population &, RawIndIterator) const;

%feature("docstring") simuPOP::BaseSelector::apply "Obsolete or undocumented function."

%ignore simuPOP::BaseSelector::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::BaseSelector::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::BaseVspSplitter "

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
    by function BaseVspSplitter.name() or Population.subPopName().
    Only one VSP splitter can be assigned to a population, which
    defined VSPs for all its subpopulations. If different splitters
    are needed for different subpopulations, a CombinedSplitter can be
    used.

"; 

%feature("docstring") simuPOP::BaseVspSplitter::BaseVspSplitter "

Usage:

    BaseVspSplitter(names=[])

Details:

    This is a virtual class that cannot be instantiated.

"; 

%feature("docstring") simuPOP::BaseVspSplitter::clone "

Usage:

    x.clone()

Details:

    All VSP splitter defines a clone() function to create an identical
    copy of itself.

"; 

%feature("docstring") simuPOP::BaseVspSplitter::~BaseVspSplitter "

Usage:

    x.~BaseVspSplitter()

"; 

%ignore simuPOP::BaseVspSplitter::activatedSubPop() const;

%ignore simuPOP::BaseVspSplitter::size(const Population &pop, size_t subPop, size_t virtualSubPop) const;

%feature("docstring") simuPOP::BaseVspSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs defined by this splitter.

"; 

%ignore simuPOP::BaseVspSplitter::contains(const Population &pop, size_t ind, vspID vsp) const;

%ignore simuPOP::BaseVspSplitter::activate(const Population &pop, size_t subPop, size_t virtualSubPop);

%ignore simuPOP::BaseVspSplitter::deactivate(size_t subPop);

%feature("docstring") simuPOP::BaseVspSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of VSP vsp (an index between 0 and
    numVirtualSubPop()).

"; 

%feature("docstring") simuPOP::BaseVspSplitter::vspByName "

Usage:

    x.vspByName(name)

Details:

    Return the index of a virtual subpopulation from its name. If
    multiple virtual subpopulations share the same name, the first vsp
    is returned.

"; 

%feature("docstring") simuPOP::Bernullitrials "

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
    of trials p1 p2 p3 p4 p5  trial 1 trial 2 ... trial N  We expect
    that N is big (usually populaiton size) and p_i are small  using
    fast bernulliTrial method for fix p, we can fill up this table
    very quickly column by column  This class will provide easy access
    to row (each trial) or column (called each prob) of this table.
    if this table is accessed row by row (each trial), a internal
    index is used.  if index exceeds N, trials will be generated all
    again. if trial will be called, e.g., N+2 times all the time, this
    treatment might not be very efficient.

"; 

%ignore simuPOP::Bernullitrials::Bernullitrials(RNG &);

%feature("docstring") simuPOP::Bernullitrials::Bernullitrials "

Usage:

    Bernullitrials(, prob, trials=0)

"; 

%feature("docstring") simuPOP::Bernullitrials::~Bernullitrials "

Usage:

    x.~Bernullitrials()

"; 

%ignore simuPOP::Bernullitrials::trialSize() const;

%feature("docstring") simuPOP::Bernullitrials::probSize "

Usage:

    x.probSize()

"; 

%ignore simuPOP::Bernullitrials::setParameter(const vectorf &prob, size_t trials=0);

%feature("docstring") simuPOP::Bernullitrials::doTrial "

Description:

    generate the trial table, reset m_cur

Usage:

    x.doTrial()

"; 

%ignore simuPOP::Bernullitrials::curTrial();

%feature("docstring") simuPOP::Bernullitrials::trial "

Description:

    if necessary, do trail again.

Usage:

    x.trial()

"; 

%feature("docstring") simuPOP::Bernullitrials::trialSucc "

Usage:

    x.trialSucc(idx)

"; 

%feature("docstring") simuPOP::Bernullitrials::probFirstSucc "

Usage:

    x.probFirstSucc()

"; 

%feature("docstring") simuPOP::Bernullitrials::probNextSucc "

Usage:

    x.probNextSucc(pos)

"; 

%feature("docstring") simuPOP::Bernullitrials::trialFirstSucc "

Usage:

    x.trialFirstSucc(idx)

"; 

%feature("docstring") simuPOP::Bernullitrials::trialNextSucc "

Usage:

    x.trialNextSucc(idx, pos)

"; 

%feature("docstring") simuPOP::Bernullitrials::setTrialSucc "

Usage:

    x.setTrialSucc(idx, succ)

"; 

%feature("docstring") simuPOP::Bernullitrials::trialSuccRate "

Description:

    return the succ rate for one index, used for verification pruposes

Usage:

    x.trialSuccRate(index)

"; 

%feature("docstring") simuPOP::Bernullitrials::probSuccRate "

Description:

    return the succ rate for current trial, used for verification
    pruposes

Usage:

    x.probSuccRate()

"; 

%ignore simuPOP::Bernullitrials::probabilities();

%ignore simuPOP::BinomialNumOffModel;

%feature("docstring") simuPOP::BinomialNumOffModel::BinomialNumOffModel "

Usage:

    BinomialNumOffModel(N, mu)

"; 

%feature("docstring") simuPOP::BinomialNumOffModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::BinomialNumOffModel::getNumOff "

Usage:

    x.getNumOff(ssize_t)

"; 

%feature("docstring") simuPOP::CloneGenoTransmitter "

Details:

    This during mating operator copies parental genotype directly to
    offspring. This operator works for all mating schemes when one or
    two parents are involved. If both parents are passed, maternal
    genotype are copied. In addition to genotypes on all non-
    customized or specified chromosomes, sex and information fields
    are by default also coped copied from parent to offspring.

"; 

%feature("docstring") simuPOP::CloneGenoTransmitter::CloneGenoTransmitter "

Usage:

    CloneGenoTransmitter(output=\"\", chroms=ALL_AVAIL, begin=0,
      end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL,
      infoFields=ALL_AVAIL)

Details:

    Create a clone genotype transmitter (a during-mating operator)
    that copies genotypes from parents to offspring. If two parents
    are specified, genotypes are copied maternally. After genotype
    transmission, offspring sex and affection status is copied from
    the parent even if sex has been determined by an offspring
    generator. All or specified information fields (parameter
    infoFields, default to ALL_AVAIL) will also be copied from parent
    to offspring. Parameters subPops is ignored. This operator by
    default copies genotypes on all autosome and sex chromosomes
    (excluding customized chromosomes), unless a parameter chroms is
    used to specify which chromosomes to copy. This operator also
    copies allelic lineage when it is executed in a module with
    lineage allele type.

"; 

%feature("docstring") simuPOP::CloneGenoTransmitter::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::CloneGenoTransmitter::describe "Obsolete or undocumented function."

%ignore simuPOP::CloneGenoTransmitter::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%ignore simuPOP::CloneGenoTransmitter::parallelizable() const;

%ignore simuPOP::CombinedAlleleIterator;

%feature("docstring") simuPOP::CombinedAlleleIterator::CombinedAlleleIterator "

Usage:

    CombinedAlleleIterator()

"; 

%feature("docstring") simuPOP::CombinedAlleleIterator::valid "

Usage:

    x.valid()

"; 

%feature("docstring") simuPOP::CombinedAlleleIterator::assignIfDiffer "

Usage:

    x.assignIfDiffer(a)

"; 

%feature("docstring") simuPOP::CombinedAlleleIterator::value "

Usage:

    x.value()

"; 

%feature("docstring") simuPOP::CombinedAlleleIterator::ptr "

Usage:

    x.ptr()

"; 

%feature("docstring") simuPOP::CombinedAlleleIterator::advance "

Usage:

    x.advance(T, T::pointer, it, p, valid)

"; 

%feature("docstring") simuPOP::CombinedParentsChooser "

Details:

    This parent chooser accepts two parent choosers. It takes one
    parent from each parent chooser and return them as father and
    mother. Because two parent choosers do not have to choose parents
    from the same virtual subpopulation, this parent chooser allows
    you to choose parents from different subpopulations.

"; 

%feature("docstring") simuPOP::CombinedParentsChooser::CombinedParentsChooser "

Usage:

    CombinedParentsChooser(fatherChooser, motherChooser)

Details:

    Create a Python parent chooser using two parent choosers
    fatherChooser and motherChooser. It takes one parent from each
    parent chooser and return them as father and mother. If two valid
    parents are returned, the first valid parent (father) will be used
    for fatherChooser, the second valid parent (mother) will be used
    for motherChooser. Although these two parent choosers are supposed
    to return a father and a mother respectively, the sex of returned
    parents are not checked so it is possible to return parents with
    the same sex using this parents chooser.

"; 

%ignore simuPOP::CombinedParentsChooser::CombinedParentsChooser(const CombinedParentsChooser &rhs);

%feature("docstring") simuPOP::CombinedParentsChooser::~CombinedParentsChooser "

Usage:

    x.~CombinedParentsChooser()

"; 

%feature("docstring") simuPOP::CombinedParentsChooser::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::CombinedParentsChooser::describe "Obsolete or undocumented function."

%ignore simuPOP::CombinedParentsChooser::parallelizable() const;

%ignore simuPOP::CombinedParentsChooser::initialize(Population &pop, size_t sp);

%ignore simuPOP::CombinedParentsChooser::finalize(Population &pop, size_t sp);

%ignore simuPOP::CombinedParentsChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::CombinedSplitter "

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

%feature("docstring") simuPOP::CombinedSplitter::CombinedSplitter "

Usage:

    CombinedSplitter(splitters=[], vspMap=[], names=[])

Details:

    Create a combined splitter using a list of splitters. For example,
    CombinedSplitter([SexSplitter(), AffectionSplitter()]) defines a
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

%ignore simuPOP::CombinedSplitter::CombinedSplitter(const CombinedSplitter &rhs);

%feature("docstring") simuPOP::CombinedSplitter::~CombinedSplitter "

Usage:

    x.~CombinedSplitter()

"; 

%feature("docstring") simuPOP::CombinedSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::CombinedSplitter::size(const Population &pop, size_t subPop, size_t virtualSubPop) const;

%feature("docstring") simuPOP::CombinedSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs defined by this splitter, which is the
    sum of the number of VSPs of all combined splitters.

"; 

%ignore simuPOP::CombinedSplitter::contains(const Population &pop, size_t ind, vspID vsp) const;

%ignore simuPOP::CombinedSplitter::activate(const Population &pop, size_t subPop, size_t virtualSubPop);

%feature("docstring") simuPOP::CombinedSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of a VSP vsp, which is the name a VSP defined by
    one of the combined splitters unless a new set of names is
    specified. If a vspMap was used, names from different VSPs will be
    joined by \"or\".

"; 

%ignore simuPOP::ConstNumOffModel;

%feature("docstring") simuPOP::ConstNumOffModel::ConstNumOffModel "

Usage:

    ConstNumOffModel(numOff)

"; 

%feature("docstring") simuPOP::ConstNumOffModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::ConstNumOffModel::getNumOff "

Usage:

    x.getNumOff(ssize_t)

"; 

%feature("docstring") simuPOP::ConstNumOffModel::parallelizable "

Usage:

    x.parallelizable()

"; 

%feature("docstring") simuPOP::ContextMutator "

Details:

    This context-dependent mutator accepts a list of mutators and use
    one of them to mutate an allele depending on the context of the
    mutated allele.

"; 

%feature("docstring") simuPOP::ContextMutator::ContextMutator "

Usage:

    ContextMutator(rates=[], loci=ALL_AVAIL, mutators=[],
      contexts=[], mapIn=[], mapOut=[], output=\">\", begin=0, end=-1,
      step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL,
      infoFields=\"ind_id\", lineageMode=FROM_INFO)

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
    mutator and BaseOperator for descriptions of other parameters.

"; 

%feature("docstring") simuPOP::ContextMutator::clone "Obsolete or undocumented function."

%ignore simuPOP::ContextMutator::mutate(Allele allele, size_t locus) const;

%feature("docstring") simuPOP::ContextMutator::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::ControlledOffspringGenerator "

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
    detail in \"Peng et al, (2007) PLoS Genetics\".

"; 

%feature("docstring") simuPOP::ControlledOffspringGenerator::ControlledOffspringGenerator "

Usage:

    ControlledOffspringGenerator(loci, alleles, freqFunc, ops=[],
      numOffspring=1, sexMode=RANDOM_SEX)

Details:

    Create an offspring generator that selects offspring so that
    allele frequency at specified loci in the offspring generation
    reaches specified allele frequency. At the beginning of each
    generation, expected allele frequency of alleles at loci is
    returned from a user-defined trajectory function freqFunc.
    Parameter loci can be a list of loci indexes, names, or ALL_AVAIL.
    If there is no subpopulation, this function should return a list
    of frequencies for each locus. If there are multiple
    subpopulations, freqFunc can return a list of allele frequencies
    for all subpopulations or combined frequencies that ignore
    population structure. In the former case, allele frequencies
    should be arranged by loc0_sp0, loc1_sp0, ... loc0_sp1, loc1_sp1,
    ..., and so on. In the latter case, overall expected number of
    alleles are scattered to each subpopulation in proportion to
    existing number of alleles in each subpopulation, using a
    multinomial distribution.  After the expected alleles are
    calculated, this offspring generator accept and reject families
    according to their genotype at loci until allele frequecies reach
    their expected values. The rest of the offspring generation is
    then filled with families without only wild type alleles at these
    loci.  This offspring generator is derived from class
    OffspringGenerator. Please refer to class OffspringGenerator for a
    detailed description of parameters ops, numOffspring and sexMode.

"; 

%ignore simuPOP::ControlledOffspringGenerator::ControlledOffspringGenerator(const ControlledOffspringGenerator &rhs);

%ignore simuPOP::ControlledOffspringGenerator::initialize(const Population &pop, size_t subPop);

%ignore simuPOP::ControlledOffspringGenerator::generateOffspring(Population &pop, Population &offPop, Individual *dad, Individual *mom, RawIndIterator &offBegin, RawIndIterator &offEnd);

%feature("docstring") simuPOP::ControlledOffspringGenerator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::ControlledOffspringGenerator::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::DiscardIf "

Details:

    This operator discards individuals according to either an
    expression that evaluates according to individual information
    field, or a Python function that accepts individual and its
    information fields.

"; 

%feature("docstring") simuPOP::DiscardIf::DiscardIf "

Usage:

    DiscardIf(cond, exposeInd=\"\", output=\"\", begin=0, end=-1,
      step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create an operator that discard individuals according to an
    expression or the return value of a Python function (parameter
    cond). This operator can be applied to a population before or
    after mating, or to offspring during mating. If an expression is
    passed to cond, it will be evalulated with each individual's
    information fields (see operator InfoEval for details). If
    exposeInd is non-empty, individuals will be available for
    evaluation in the expression as an variable with name spacied by
    exposeInd. If the expression is evaluated to be True, individuals
    (if applied before or after mating) or offspring (if applied
    during mating) will be removed or discard. If a function is passed
    to cond, it should accept paramters ind and pop or names of
    information fields when it is applied to a population (pre or post
    mating), or parameters off, dad, mom, pop (parental population),
    or names of information fields if the operator is applied during
    mating. Individuals will be discarded if this function returns
    True. A constant expression (e.g. True) is also acceptable).
    Because this operator supports parameter subPops, only individuals
    belonging to specified (virtual) subpopulations will be screened.

"; 

%feature("docstring") simuPOP::DiscardIf::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::DiscardIf::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::DiscardIf::apply "Obsolete or undocumented function."

%ignore simuPOP::DiscardIf::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad, Individual *mom) const;

%feature("docstring") simuPOP::DiscardIf::~DiscardIf "

Usage:

    x.~DiscardIf()

"; 

%feature("docstring") simuPOP::Dumper "

Details:

    This operator dumps the content of a population in a human
    readable format. Because this output format is not structured and
    can not be imported back to simuPOP, this operator is usually used
    to dump a small population to a terminal for demonstration and
    debugging purposes.

"; 

%feature("docstring") simuPOP::Dumper::Dumper "

Usage:

    Dumper(genotype=True, structure=True, ancGens=UNSPECIFIED,
      width=1, max=100, loci=[], output=\">\", begin=0, end=-1, step=1,
      at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a operator that dumps the genotype structure (if structure
    is True) and genotype (if genotype is True) to an output ( default
    to standard terminal output). Because a population can be large,
    this operator will only output the first 100 (parameter max)
    individuals of the present generation (parameter ancGens). All
    loci will be outputed unless parameter loci are used to specify a
    subset of loci. If a list of (virtual) subpopulations are
    specified, this operator will only output individuals in these
    outputs. Please refer to class BaseOperator for a detailed
    explanation for common parameters such as output and stage.

"; 

%feature("docstring") simuPOP::Dumper::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::Dumper::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::Dumper::~Dumper "

Description:

    destructor.

Usage:

    x.~Dumper()

"; 

%feature("docstring") simuPOP::Dumper::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::Exception "

Description:

    exception handler. Exceptions will be passed to Python.

"; 

%feature("docstring") simuPOP::Exception::Exception "

Description:

    constructor

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

%ignore simuPOP::Expression::expr() const;

%ignore simuPOP::Expression::stmts() const;

%ignore simuPOP::Expression::setLocalDict(PyObject *dict) const;

%ignore simuPOP::Expression::empty();

%ignore simuPOP::Expression::setExpr(const string &expr=string());

%ignore simuPOP::Expression::setStmts(const string &stmts=string());

%ignore simuPOP::Expression::evaluate() const;

%ignore simuPOP::Expression::valueAsBool() const;

%ignore simuPOP::Expression::valueAsInt() const;

%ignore simuPOP::Expression::valueAsDouble() const;

%ignore simuPOP::Expression::valueAsString() const;

%ignore simuPOP::Expression::valueAsArray() const;

%feature("docstring") simuPOP::FiniteSitesMutator "

Details:

    This is an infite site mutation model in mutational space. The
    alleles in the population is assumed to be locations of mutants. A
    mutation rate is given that mutate alleles in 'regions'. If number
    of mutants for an individual exceed the number of loci, 10 loci
    will be added to everyone in the population.

"; 

%feature("docstring") simuPOP::FiniteSitesMutator::FiniteSitesMutator "

Usage:

    FiniteSitesMutator(rate, ranges, model=1, output=\"\", begin=0,
      end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL,
      infoFields=\"ind_id\", lineageMode=FROM_INFO)

Details:

    This operator accepts a list of ranges which is the 'real range'
    of each chromosome. Mutation happens with muation rate rate and
    mutants will be recorded to the population (instead of alleles).
    By default, this mutator assumes an finite-allele model where all
    mutations are allowed and if a mutant (allele 1) is mutated, it
    will be mutated to allele 0 (back mutation). Alternatively (model
    = 2), an infinite-sites mutation model can be used where mutations
    can happen only at a new locus. Mutations happen at a locus with
    existing mutants will be moved to a random locus without existing
    mutant. A warning message will be printed if there is no vacant
    locus available. If a valid output is given, mutants will be
    outputted in the format of \"gen mutant ind type\" where type is 0
    for forward (0->1), 1 for backward (1->0), 2 for relocated
    mutations, and 3 for ignored mutation because no vacent locus is
    available. The second mode has the advantage that all mutants in
    the simulated population can be traced to a single mutation event.
    If the regions are reasonably wide and mutation rates are low,
    these two mutation models should yield similar results.

"; 

%feature("docstring") simuPOP::FiniteSitesMutator::~FiniteSitesMutator "

Description:

    destructor.

Usage:

    x.~FiniteSitesMutator()

"; 

%feature("docstring") simuPOP::FiniteSitesMutator::apply "

Usage:

    x.apply(pop)

Details:

    Apply an operator to population pop directly, without checking its
    applicability.

"; 

%feature("docstring") simuPOP::FiniteSitesMutator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::FiniteSitesMutator::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::floatList "

"; 

%feature("docstring") simuPOP::floatList::floatList "

Usage:

    floatList(obj=None)

"; 

%ignore simuPOP::floatList::floatList(double val);

%ignore simuPOP::floatList::elems() const;

%feature("docstring") simuPOP::floatListFunc "

"; 

%feature("docstring") simuPOP::floatListFunc::floatListFunc "

Usage:

    floatListFunc(func)

"; 

%ignore simuPOP::floatListFunc::floatListFunc(double val);

%ignore simuPOP::floatListFunc::empty() const;

%ignore simuPOP::floatListFunc::size() const;

%ignore simuPOP::floatListFunc::func() const;

%feature("docstring") simuPOP::floatMatrix "

"; 

%feature("docstring") simuPOP::floatMatrix::floatMatrix "

Usage:

    floatMatrix(obj=None)

"; 

%ignore simuPOP::floatMatrix::empty() const;

%ignore simuPOP::floatMatrix::elems() const;

%ignore simuPOP::FuncNumOffModel;

%feature("docstring") simuPOP::FuncNumOffModel::FuncNumOffModel "

Usage:

    FuncNumOffModel(func)

"; 

%feature("docstring") simuPOP::FuncNumOffModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::FuncNumOffModel::getNumOff "

Usage:

    x.getNumOff(gen)

"; 

%feature("docstring") simuPOP::FuncNumOffModel::reset "

Usage:

    x.reset()

"; 

%ignore simuPOP::FuncSexModel;

%feature("docstring") simuPOP::FuncSexModel::FuncSexModel "

Usage:

    FuncSexModel(func)

"; 

%feature("docstring") simuPOP::FuncSexModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::FuncSexModel::getSex "

Usage:

    x.getSex(count)

"; 

%feature("docstring") simuPOP::FuncSexModel::reset "

Usage:

    x.reset()

"; 

%feature("docstring") simuPOP::FuncSexModel::parallelizable "

Usage:

    x.parallelizable()

"; 

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
    accessible from both Individual and Population classes. Currently,
    a genotypic structure consists of
    *   Ploidy, namely the number of homologous sets of chromosomes,
    of a population. Haplodiploid population is also supported.
    *   Number of chromosomes and number of loci on each chromosome.
    *   Positions of loci, which determine the relative distance
    between loci on the same chromosome. No unit is assumed so these
    positions can be ordinal (1, 2, 3, ..., the default), in physical
    distance (bp, kb or mb), or in map distance (e.g. centiMorgan)
    depending on applications.
    *   Names of alleles, which can either be shared by all loci or be
    specified for each locus.
    *   Names of loci and chromosomes.
    *   Names of information fields attached to each individual. In
    addition to basic property access functions, this class provides
    some utility functions such as locusByName, which looks up a locus
    by its name.

"; 

%feature("docstring") simuPOP::GenoStruTrait::GenoStruTrait "

Usage:

    GenoStruTrait()

Details:

    A GenoStruTrait object is created with the construction of a
    Population object and cannot be initialized directly.

"; 

%ignore simuPOP::GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru &loci, const vectoru &chromTypes, bool haplodiploid, const vectorf &lociPos, const vectorstr &chromNames, const matrixstr &alleleNames, const vectorstr &lociNames, const vectorstr &infoFields);

%ignore simuPOP::GenoStruTrait::setGenoStructure(const GenoStructure &rhs);

%ignore simuPOP::GenoStruTrait::setGenoStruIdx(size_t idx);

%ignore simuPOP::GenoStruTrait::swapGenoStru(GenoStruTrait &rhs);

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

%ignore simuPOP::GenoStruTrait::gsRemoveLoci(const vectoru &kept);

%ignore simuPOP::GenoStruTrait::gsAddChrom(const vectorf &lociPos, const vectorstr &lociNames, const string &chromName, const matrixstr &alleleNames, size_t chromType) const;

%ignore simuPOP::GenoStruTrait::gsSetAlleleNames(const lociList &loci, const matrixstr &alleleNames);

%ignore simuPOP::GenoStruTrait::gsAddLoci(const vectoru &chrom, const vectorf &pos, const vectorstr &lociNames, const matrixstr &alleleNames, vectoru &newIndex) const;

%ignore simuPOP::GenoStruTrait::genoStru() const;

%ignore simuPOP::GenoStruTrait::genoStruIdx() const;

%feature("docstring") simuPOP::GenoStruTrait::ploidy "

Usage:

    x.ploidy()

Details:

    return the number of homologous sets of chromosomes, specified by
    the ploidy parameter of the Population function. Return 2 for a
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

    return the number of loci on chromosome chrom.

"; 

%feature("docstring") simuPOP::GenoStruTrait::numLoci "

Usage:

    x.numLoci()

Details:

    return a list of the number of loci on all chromosomes.

"; 

%ignore simuPOP::GenoStruTrait::chromX() const;

%ignore simuPOP::GenoStruTrait::chromY() const;

%ignore simuPOP::GenoStruTrait::customizedChroms() const;

%ignore simuPOP::GenoStruTrait::mitochondrial() const;

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
    parameter of the Population function.

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociPos "

Usage:

    x.lociPos()

Details:

    return the positions of all loci, specified by the lociPos
    prameter of the Population function. The default positions are 1,
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

    return the absolute index of locus locus on chromosome chrom. c.f.
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

    return the name of a chromosome chrom.

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

    return the type of a chromosome chrom (CUSTOMIZED, AUTOSOME,
    CHROMOSOME_X, CHROMOSOME_Y or MITOCHONDRIAL.

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromTypes "

Usage:

    x.chromTypes()

Details:

    return the type of all chromosomes (CUSTOMIZED, AUTOSOME,
    CHROMOSOME_X, CHROMOSOME_Y, or MITOCHONDRIAL).

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
    alleleNames parameter of the Population function. locus could be
    ignored if alleles at all loci share the same names. If the name
    of an allele is unspecified, its value ('0', '1', '2', etc) is
    returned.

"; 

%ignore simuPOP::GenoStruTrait::allAlleleNames() const;

%feature("docstring") simuPOP::GenoStruTrait::alleleNames "

Usage:

    x.alleleNames(locus=0)

Details:

    return a list of allele names at given by the alleleNames
    parameter of the Population function. locus could be ignored if
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
    parameter of the Population function. An empty string will be
    returned if no name has been given to locus locus.

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociNames "

Usage:

    x.lociNames()

Details:

    return the names of all loci specified by the lociNames parameter
    of the Population function. An empty list will be returned if
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

%feature("docstring") simuPOP::GenoTransmitter "

Details:

    This during mating operator is the base class of all genotype
    transmitters. It is made available to users because it provides a
    few member functions that can be used by derived transmitters, and
    by customized Python during mating operators.

"; 

%feature("docstring") simuPOP::GenoTransmitter::GenoTransmitter "

Usage:

    GenoTransmitter(output=\"\", begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a base genotype transmitter.

"; 

%feature("docstring") simuPOP::GenoTransmitter::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::GenoTransmitter::clearChromosome "

Usage:

    x.clearChromosome(ind, ploidy, chrom)

Details:

    Clear (set alleles to zero) chromosome chrom on the ploidy-th
    homologous set of chromosomes of individual ind. It is equivalent
    to ind.setGenotype([0], ploidy, chrom), except that it also clears
    allele lineage if it is executed in a module with lineage allele
    type.

"; 

%feature("docstring") simuPOP::GenoTransmitter::copyChromosome "

Usage:

    x.copyChromosome(parent, parPloidy, offspring, ploidy, chrom)

Details:

    Transmit chromosome chrom on the parPloidy set of homologous
    chromosomes from parent to the ploidy set of homologous
    chromosomes of offspring. It is equivalent to
    offspring.setGenotype(parent.genotype(parPloidy, chrom), polidy,
    chrom), except that it also copies allelic lineage when it is
    executed in a module with lineage allele type.

"; 

%feature("docstring") simuPOP::GenoTransmitter::copyChromosomes "

Usage:

    x.copyChromosomes(parent, parPloidy, offspring, ploidy)

Details:

    Transmit the parPloidy set of homologous chromosomes from parent
    to the ploidy set of homologous chromosomes of offspring.
    Customized chromosomes are not copied. It is equivalent to
    offspring.setGenotype(parent.genotype(parPloidy), ploidy), except
    that it also copies allelic lineage when it is executed in a
    module with lineage allele type.

"; 

%feature("docstring") simuPOP::GenoTransmitter::describe "Obsolete or undocumented function."

%ignore simuPOP::GenoTransmitter::initialize(const Individual &ind) const;

%ignore simuPOP::GenoTransmitter::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%ignore simuPOP::GenoTransmitter::initializeIfNeeded(const Individual &ind) const;

%feature("docstring") simuPOP::GenotypeSplitter "

Details:

    This class defines a VSP splitter that defines VSPs according to
    individual genotype at specified loci.

"; 

%feature("docstring") simuPOP::GenotypeSplitter::GenotypeSplitter "

Usage:

    GenotypeSplitter(loci, alleles, phase=False, names=[])

Details:

    Create a splitter that defines VSPs by individual genotype at loci
    (can be indexes or names of one or more loci). Each list in a list
    allele defines a VSP, which is a list of allowed alleles at these
    loci. If only one VSP is defined, the outer list of the nested
    list can be ignored. If phase if true, the order of alleles in
    each list is significant. If more than one set of alleles are
    given, Individuals having either of them is qualified.  For
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

%feature("docstring") simuPOP::GenotypeSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::GenotypeSplitter::size(const Population &pop, size_t subPop, size_t virtualSubPop) const;

%feature("docstring") simuPOP::GenotypeSplitter::numVirtualSubPop "

Description:

    number of virtual subpops of subpopulation sp

Usage:

    x.numVirtualSubPop()

"; 

%ignore simuPOP::GenotypeSplitter::contains(const Population &pop, size_t ind, vspID vsp) const;

%ignore simuPOP::GenotypeSplitter::activate(const Population &pop, size_t subPop, size_t virtualSubPop);

%feature("docstring") simuPOP::GenotypeSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return name of VSP vsp, which is \"Genotype loc1,loc2:genotype\" as
    defined by parameters loci and alleles. A user provided name will
    be returned if specified.

"; 

%ignore simuPOP::GeometricNumOffModel;

%feature("docstring") simuPOP::GeometricNumOffModel::GeometricNumOffModel "

Usage:

    GeometricNumOffModel(p)

"; 

%feature("docstring") simuPOP::GeometricNumOffModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::GeometricNumOffModel::getNumOff "

Usage:

    x.getNumOff(ssize_t)

"; 

%ignore simuPOP::GlobalSeqSexModel;

%feature("docstring") simuPOP::GlobalSeqSexModel::GlobalSeqSexModel "

Usage:

    GlobalSeqSexModel(sex)

"; 

%feature("docstring") simuPOP::GlobalSeqSexModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::GlobalSeqSexModel::getSex "

Usage:

    x.getSex(UINT)

"; 

%feature("docstring") simuPOP::GlobalSeqSexModel::reset "

Usage:

    x.reset()

"; 

%feature("docstring") simuPOP::GlobalSeqSexModel::parallelizable "

Usage:

    x.parallelizable()

"; 

%feature("docstring") simuPOP::HaplodiploidGenoTransmitter "

Details:

    A genotype transmitter (during-mating operator) for haplodiploid
    populations. The female parent is considered as diploid and the
    male parent is considered as haploid (only the first homologous
    copy is valid). If the offspring is FEMALE, she will get a random
    copy of two homologous chromosomes of her mother, and get the only
    paternal copy from her father. If the offspring is MALE, he will
    only get a set of chromosomes from his mother.

"; 

%feature("docstring") simuPOP::HaplodiploidGenoTransmitter::HaplodiploidGenoTransmitter "

Usage:

    HaplodiploidGenoTransmitter(output=\"\", begin=0, end=-1, step=1,
      at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a haplodiploid genotype transmitter (during-mating
    operator) that transmit parental genotypes from parents to
    offspring in a haplodiploid population. Parameters subPops and
    infoFields are ignored. This operator also copies allelic lineage
    when it is executed in a module with lineage allele type.

"; 

%feature("docstring") simuPOP::HaplodiploidGenoTransmitter::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::HaplodiploidGenoTransmitter::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::HaplodiploidGenoTransmitter::initialize "Obsolete or undocumented function."

%ignore simuPOP::HaplodiploidGenoTransmitter::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::HeteroMating "

Details:

    A heterogeneous mating scheme that applies a list of homogeneous
    mating schemes to different (virtual) subpopulations.

"; 

%feature("docstring") simuPOP::HeteroMating::HeteroMating "

Usage:

    HeteroMating(matingSchemes, subPopSize=[],
      shuffleOffspring=True)

Details:

    Create a heterogeneous mating scheme that will apply a list of
    homogeneous mating schemes matingSchemes to different (virtual)
    subpopulations. The size of the offspring generation is determined
    by parameter subPopSize, which can be a list of subpopulation
    sizes or a Python function that returns a list of subpopulation
    sizes at each generation. Please refer to class MatingScheme for a
    detailed explanation of this parameter.  Each mating scheme
    defined in matingSchemes can be applied to one or more (virtual)
    subpopulation. If parameter subPops is not specified, a mating
    scheme will be applied to all subpopulations. If a list of
    (virtual) subpopulation is specified, the mating scheme will be
    applied to specific (virtual) subpopulations.  If multiple mating
    schemes are applied to the same subpopulation, a weight (parameter
    weight) can be given to each mating scheme to determine how many
    offspring it will produce. The default for all mating schemes are
    0. In this case, the number of offspring each mating scheme
    produces is proportional to the size of its parental (virtual)
    subpopulation. If all weights are negative, the numbers of
    offspring are determined by the multiplication of the absolute
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

%feature("docstring") simuPOP::HeteroMating::~HeteroMating "

Description:

    destructor

Usage:

    x.~HeteroMating()

"; 

%ignore simuPOP::HeteroMating::HeteroMating(const HeteroMating &rhs);

%feature("docstring") simuPOP::HeteroMating::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::HeteroMating::describe "Obsolete or undocumented function."

%ignore simuPOP::HeteroMating::mate(Population &pop, Population &scratch);

%feature("docstring") simuPOP::HomoMating "

Details:

    A homogeneous mating scheme that uses a parent chooser to choose
    parents from a prental generation, and an offspring generator to
    generate offspring from chosen parents. It can be either used
    directly, or within a heterogeneous mating scheme. In the latter
    case, it can be applied to a (virtual) subpopulation.

"; 

%feature("docstring") simuPOP::HomoMating::HomoMating "

Usage:

    HomoMating(chooser, generator, subPopSize=[], subPops=ALL_AVAIL,
      weight=0)

Details:

    Create a homogeneous mating scheme using a parent chooser chooser
    and an offspring generator generator.  If this mating scheme is
    used directly in a simulator, it will be responsible for creating
    an offspring population according to parameter subPopSize. This
    parameter can be a list of subpopulation sizes (or a number if
    there is only one subpopulation) or a Python function which will
    be called at each generation to determine the subpopulation sizes
    of the offspring generation. Please refer to class MatingScheme
    for details about this parameter.  If this mating shcme is used
    within a heterogeneous mating scheme. Parameters subPops and
    weight are used to determine which (virtual) subpopulations this
    mating scheme will be applied to, and how many offspring this
    mating scheme will produce. Please refer to mating scheme
    HeteroMating for the use of these two parameters.

"; 

%feature("docstring") simuPOP::HomoMating::~HomoMating "

Description:

    destructor

Usage:

    x.~HomoMating()

"; 

%ignore simuPOP::HomoMating::HomoMating(const HomoMating &rhs);

%feature("docstring") simuPOP::HomoMating::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::HomoMating::describe "Obsolete or undocumented function."

%ignore simuPOP::HomoMating::subPops() const;

%ignore simuPOP::HomoMating::weight() const;

%ignore simuPOP::HomoMating::mateSubPop(Population &pop, Population &offPop, size_t subPop, RawIndIterator offBegin, RawIndIterator offEnd);

%feature("docstring") simuPOP::IdTagger "

Details:

    An IdTagger gives a unique ID for each individual it is applies
    to. These ID can be used to uniquely identify an individual in a
    multi-generational population and be used to reliably reconstruct
    a Pedigree.  To ensure uniqueness across populations, a single
    source of ID is used for this operator. individual IDs are
    assigned consecutively starting from 1. Value 1 instead of 0 is
    used because most software applications use 0 as missing values
    for parentship. If you would like to reset the sequence or start
    from a different number, you can call the reset(startID) function
    of any IdTagger.  An IdTagger is usually used during-mating to
    assign ID to each offspring. However, if it is applied directly to
    a population, it will assign unique IDs to all individuals in this
    population. This property is usually used in the preOps parameter
    of function Simulator.evolve to assign initial ID to a population.

"; 

%feature("docstring") simuPOP::IdTagger::IdTagger "

Usage:

    IdTagger(begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, output=\"\", infoFields=\"ind_id\")

Details:

    Create an IdTagger that assign an unique ID for each individual it
    is applied to. The IDs are created sequentially and are stored in
    an information field specified in parameter infoFields (default to
    ind_id). This operator is considered a during-mating operator but
    it can be used to set ID for all individuals of a population when
    it is directly applied to the population.

"; 

%feature("docstring") simuPOP::IdTagger::~IdTagger "

Usage:

    x.~IdTagger()

"; 

%feature("docstring") simuPOP::IdTagger::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::IdTagger::reset "

Usage:

    x.reset(startID=1)

Details:

    Reset the global individual ID number so that IdTaggers will start
    from id (default to 1) again.

"; 

%feature("docstring") simuPOP::IdTagger::apply "Obsolete or undocumented function."

%ignore simuPOP::IdTagger::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::IdTagger::clone "Obsolete or undocumented function."

%ignore simuPOP::IdTagger::parallelizable() const;

%feature("docstring") simuPOP::IfElse "

Details:

    This operator uses a condition, which can be a fixed condition, an
    expression or a user-defined function, to determine which
    operators to be applied when this operator is applied. A list of
    if-operators will be applied when the condition is True. Otherwise
    a list of else-operators will be applied.

"; 

%feature("docstring") simuPOP::IfElse::IfElse "

Usage:

    IfElse(cond, ifOps=[], elseOps=[], output=\">\", begin=0, end=-1,
      step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a conditional operator that will apply operators ifOps if
    condition cond is met and elseOps otherwise. If a Python
    expression (a string) is given to parameter cond, the expression
    will be evalulated in each population's local namespace when this
    operator is applied. When a Python function is specified, it
    accepts parameter pop when it is applied to a population, and one
    or more parameters pop, off, dad or mom when it is applied during
    mating. The return value of this function should be True or False.
    Otherwise, parameter cond will be treated as a fixed condition
    (converted to True or False) upon which one set of operators is
    always applied. The applicability of ifOps and elseOps are
    controlled by parameters begin, end, step, at and rep of both the
    IfElse operator and individual operators but ifOps and elseOps
    opeartors does not support negative indexes for replicate and
    generation numbers.

"; 

%feature("docstring") simuPOP::IfElse::~IfElse "

Description:

    destructor

Usage:

    x.~IfElse()

"; 

%feature("docstring") simuPOP::IfElse::clone "Obsolete or undocumented function."

%ignore simuPOP::IfElse::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::IfElse::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::IfElse::describe "Obsolete or undocumented function."

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

%feature("docstring") simuPOP::Individual "

Details:

    A Population consists of individuals with the same genotypic
    structure. An Individual object cannot be created independently,
    but refences to inidividuals can be retrieved using member
    functions of a Population object. In addition to structural
    information shared by all individuals in a population (provided by
    class GenoStruTrait), the Individual class provides member
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
    individual are separated by Individual::totNumLoci() loci. It is
    worth noting that access to invalid chromosomes, such as the Y
    chromosomes of female individuals, is not restricted.

"; 

%feature("docstring") simuPOP::Individual::Individual "

Usage:

    Individual()

Details:

    An Individual object cannot be created directly. It has to be
    accessed from a Population object using functions such as
    Population::Individual(idx).

"; 

%ignore simuPOP::Individual::Individual(const Individual &ind);

%feature("docstring") simuPOP::Individual::~Individual "

Description:

    destructor. Do nothing.

Usage:

    x.~Individual()

"; 

%ignore simuPOP::Individual::setGenoPtr(GenoIterator pos);

%ignore simuPOP::Individual::setInfoPtr(InfoIterator pos);

%ignore simuPOP::Individual::copyFrom(const Individual &rhs);

%ignore simuPOP::Individual::genoPtr() const;

%ignore simuPOP::Individual::infoPtr() const;

%feature("docstring") simuPOP::Individual::allele "

Usage:

    x.allele(idx, ploidy=-1, chrom=-1)

Details:

    return the current allele at a locus, using its absolute index
    idx. If a ploidy ploidy and/or a chromosome indexes is given, idx
    is relative to the beginning of specified homologous copy of
    chromosomes (if chrom=-1) or the beginning of the specified
    homologous copy of specified chromosome (if chrom >= 0).

"; 

%feature("docstring") simuPOP::Individual::alleleChar "

Usage:

    x.alleleChar(idx, ploidy=-1, chrom=-1)

Details:

    return the name of allele(idx, ploidy, chrom). If idx is invalid
    (e.g. second homologus copy of chromosome Y), '_' is returned.

"; 

%feature("docstring") simuPOP::Individual::setAllele "

Usage:

    x.setAllele(allele, idx, ploidy=-1, chrom=-1)

Details:

    set allele allele to a locus, using its absolute index idx. If a
    ploidy ploidy and/or a chromosome indexes are given, idx is
    relative to the beginning of specified homologous copy of
    chromosomes (if chrom=-1) or the beginning of the specified
    homologous copy of specified chromosome (if chrom >= 0).

"; 

%feature("docstring") simuPOP::Individual::alleleLineage "

Usage:

    x.alleleLineage(idx, ploidy=-1, chrom=-1)

Details:

    return the lineage of the allele at a locus, using its absolute
    index idx. If a ploidy ploidy and/or a chromosome indexes is
    given, idx is relative to the beginning of specified homologous
    copy of chromosomes (if chrom=-1) or the beginning of the
    specified homologous copy of specified chromosome (if chrom >= 0).
    This function returns 0 for modules without lineage information.

"; 

%feature("docstring") simuPOP::Individual::setAlleleLineage "

Usage:

    x.setAlleleLineage(lineage, idx, ploidy=-1, chrom=-1)

Details:

    set lineage lineage to an allele, using its absolute index idx. If
    a ploidy ploidy and/or a chromosome indexes are given, idx is
    relative to the beginning of specified homologous copy of
    chromosomes (if chrom=-1) or the beginning of the specified
    homologous copy of specified chromosome (if chrom >= 0). This
    function does nothing for modules without lineage information.

"; 

%feature("docstring") simuPOP::Individual::genotype "

Usage:

    x.genotype(ploidy=ALL_AVAIL, chroms=ALL_AVAIL)

Details:

    return an editable array (a carray object) that represents all
    alleles of an individual. If ploidy or chroms is given, only
    alleles on the specified chromosomes and homologous copy of
    chromosomes will be returned. If multiple chromosomes are
    specified, there should not be gaps between chromosomes. This
    function ignores type of chromosomes so it will return unused
    alleles for sex and mitochondrial chromosomes.

"; 

%feature("docstring") simuPOP::Individual::mutants "

Usage:

    x.mutants(ploidy=ALL_AVAIL, chroms=ALL_AVAIL)

Details:

    return an itertor that iterate through all mutants (non-zero
    alleles) of an individual. Each mutant is presented as a tuple of
    (index, value) where index is the index of mutant ranging from
    zero to totNumLoci() * ploidy() - 1, so you will have to adjust
    indexes to check multiple alleles at a locus. If ploidy or chroms
    is given, only alleles on the specified chromosomes and homologous
    copy of chromosomes will be iterated. If multiple chromosomes are
    specified, there should not be gaps between chromosomes. This
    function ignores type of chromosomes so it will return unused
    alleles for sex and mitochondrial chromosomes.

"; 

%feature("docstring") simuPOP::Individual::lineage "

Usage:

    x.lineage(ploidy=ALL_AVAIL, chroms=ALL_AVAIL)

Details:

    return an editable array (a carray_lineage object) that represents
    the lineages of all alleles of an individual. If ploidy or chroms
    is given, only lineages on the specified chromosomes and
    homologous copy of chromosomes will be returned. If multiple
    chromosomes are specified, there should not be gaps between
    chromosomes. This function ignores type of chromosomes so it will
    return lineage of unused alleles for sex and mitochondrial
    chromosomes. A None object will be returned for modules without
    lineage information.

"; 

%ignore simuPOP::Individual::genoAtLoci(const lociList &loci);

%feature("docstring") simuPOP::Individual::setGenotype "

Usage:

    x.setGenotype(geno, ploidy=ALL_AVAIL, chroms=ALL_AVAIL)

Details:

    Fill the genotype of an individual using a list of alleles geno.
    If parameters ploidy and/or chroms are specified, alleles will be
    copied to only all or specified chromosomes on selected homologous
    copies of chromosomes. geno will be reused if its length is less
    than number of alleles to be filled. This function ignores type of
    chromosomes so it will set genotype for unused alleles for sex and
    mitochondrial chromosomes.

"; 

%feature("docstring") simuPOP::Individual::setLineage "

Usage:

    x.setLineage(lineage, ploidy=ALL_AVAIL, chroms=ALL_AVAIL)

Details:

    Fill the lineage of an individual using a list of IDs lineage. If
    parameters ploidy and/or chroms are specified, lineages will be
    copied to only all or specified chromosomes on selected homologous
    copies of chromosomes. lineage will be reused if its length is
    less than number of allelic lineage to be filled. This function
    ignores type of chromosomes so it will set lineage to unused
    alleles for sex and mitochondrial chromosomes. It does nothing for
    modules without lineage information.

"; 

%feature("docstring") simuPOP::Individual::sex "

Usage:

    x.sex()

Details:

    return the sex of an individual, 1 for male and 2 for female.

"; 

%feature("docstring") simuPOP::Individual::setSex "

Usage:

    x.setSex(sex)

Details:

    set individual sex to MALE or FEMALE.

"; 

%feature("docstring") simuPOP::Individual::affected "

Usage:

    x.affected()

Details:

    Return True if this individual is affected.

"; 

%feature("docstring") simuPOP::Individual::setAffected "

Usage:

    x.setAffected(affected)

Details:

    set affection status to affected (True or False).

"; 

%ignore simuPOP::Individual::visible() const;

%ignore simuPOP::Individual::setVisible(bool visible) const;

%ignore simuPOP::Individual::marked() const;

%ignore simuPOP::Individual::setMarked(bool mark=true) const;

%feature("docstring") simuPOP::Individual::info "

Usage:

    x.info(field)

Details:

    Return the value of an information field filed (by index or name).
    ind.info(name) is equivalent to ind.name although the function
    form allows the use of indexes of information fieldes.

"; 

%ignore simuPOP::Individual::intInfo(const uintString &field) const;

%feature("docstring") simuPOP::Individual::setInfo "

Usage:

    x.setInfo(value, field)

Details:

    set the value of an information field field (by index or name) to
    value. ind.setInfo(value, field) is equivalent to ind.field =
    value although the function form allows the use of indexes of
    information fieldes.

"; 

%ignore simuPOP::Individual::genoBegin() const;

%ignore simuPOP::Individual::genoEnd() const;

%ignore simuPOP::Individual::genoBegin(size_t p) const;

%ignore simuPOP::Individual::genoEnd(size_t p) const;

%ignore simuPOP::Individual::genoBegin(size_t p, size_t chrom) const;

%ignore simuPOP::Individual::infoBegin() const;

%ignore simuPOP::Individual::infoEnd() const;

%feature("docstring") simuPOP::Individual::cmp "

Description:

    a python function used to compare the individual objects

Usage:

    x.__cmp__(rhs)

"; 

%ignore simuPOP::Individual::swap(Individual &ind, bool swapContent=true);

%ignore simuPOP::Individual::display(ostream &out, int width=1, const vectoru &loci=vectoru());

%feature("docstring") simuPOP::IndividualIterator "

Details:

    this class implements a C++ iterator class that iterate through
    individuals in a (sub)population. If allInds are true, the
    visiblility of individuals will not be checked. Note that
    IndividualIteratorwill iterate through only visible individuals,
    and allInds is only provided when we know in advance that all
    individuals are visible. This is a way to obtain better
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

%feature("docstring") simuPOP::InfoEval "

Details:

    Unlike operator PyEval and PyExec that work at the population
    level, in a population's local namespace, operator InfoEval works
    at the individual level, working with individual information
    fields. When this operator is applied to a population, information
    fields of eligible individuals are put into the local namespace of
    the population. A Python expression is then evaluated for each
    individual. The result is written to an output.

"; 

%feature("docstring") simuPOP::InfoEval::InfoEval "

Usage:

    InfoEval(expr=\"\", stmts=\"\", usePopVars=False, exposeInd=\"\",
      output=\">\", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=[])

Details:

    Create an operator that evaluate a Python expression expr using
    individual information fields and population variables as
    variables. If exposeInd is not empty, the individual itself will
    be exposed in the population's local namespace as a variable with
    name specified by exposeInd.  A Python expression (expr) is
    evaluated for each individual. The results are converted to
    strings and are written to an output specified by parameter
    output. Optionally, a statement (or several statements separated
    by newline) can be executed before expr is evaluated. The
    evaluation of this statement may change the value of information
    fields.  This operator is by default applied post-mating. If it
    stage is set to DuringMating, it will be applied to all offspring,
    regardless of subPops settings.  Parameter usePopVars is obsolete
    because population variables are always usable in such
    expressions.

"; 

%feature("docstring") simuPOP::InfoEval::~InfoEval "

Usage:

    x.~InfoEval()

"; 

%feature("docstring") simuPOP::InfoEval::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::InfoEval::apply "Obsolete or undocumented function."

%ignore simuPOP::InfoEval::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::InfoEval::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::InfoExec "

Details:

    Operator InfoExec is similar to InfoEval in that it works at the
    individual level, using individual information fields as
    variables. This is usually used to change the value of information
    fields. For example, \"b=a*2\" will set the value of information
    field b to a*a for all individuals.

"; 

%feature("docstring") simuPOP::InfoExec::InfoExec "

Usage:

    InfoExec(stmts=\"\", usePopVars=False, exposeInd=\"\", output=\"\",
      begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=[])

Details:

    Create an operator that executes Python statements stmts using
    individual information fields and population variables as
    variables. If exposeInd is not empty, the individual itself will
    be exposed in the population's local namespace as a variable with
    name specified by exposeInd.  One or more python statements
    (stmts) are executed for each individual. Information fields of
    these individuals are then updated from the corresponding
    variables. For example, a=1 will set information field a of all
    individuals to 1, a=b will set information field a of all
    individuals to information field b or a population variable b if b
    is not an information field but a population variable, and
    a=ind.sex() will set information field a of all individuals to its
    sex (needs exposeInd='ind'.  This operator is by default applied
    post-mating. If it stage is set to DuringMating, it will be
    applied to all offspring, regardless of subPops settings.
    Parameter usePopVars is obsolete because population variables will
    always be usable.

"; 

%feature("docstring") simuPOP::InfoExec::~InfoExec "

Usage:

    x.~InfoExec()

"; 

%feature("docstring") simuPOP::InfoExec::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::InfoExec::apply "Obsolete or undocumented function."

%ignore simuPOP::InfoExec::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::InfoExec::describe "Obsolete or undocumented function."

%ignore simuPOP::InformationIterator;

%feature("docstring") simuPOP::InformationIterator::InformationIterator "

Usage:

    InformationIterator()

"; 

%feature("docstring") simuPOP::InformationIterator::valid "

Usage:

    x.valid()

"; 

%feature("docstring") simuPOP::InfoSplitter "

Details:

    This splitter defines VSPs according to the value of an
    information field of each indivdiual. A VSP is defined either by a
    value or a range of values.

"; 

%feature("docstring") simuPOP::InfoSplitter::InfoSplitter "

Usage:

    InfoSplitter(field, values=[], cutoff=[], ranges=[], names=[])

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

%feature("docstring") simuPOP::InfoSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::InfoSplitter::size(const Population &pop, size_t subPop, size_t virtualSubPop) const;

%feature("docstring") simuPOP::InfoSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs defined by this splitter, which is the
    length parameter values or the length of cutoff plus one,
    depending on which parameter is specified.

"; 

%ignore simuPOP::InfoSplitter::contains(const Population &pop, size_t ind, vspID vsp) const;

%ignore simuPOP::InfoSplitter::activate(const Population &pop, size_t subPop, size_t virtualSubPop);

%feature("docstring") simuPOP::InfoSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of a VSP vsp, which is field = value if VSPs are
    defined by values in parameter values, or field < value (the first
    VSP), v1 <= field < v2 and field >= v (the last VSP) if VSPs are
    defined by cutoff values. A user-specified name, if specified,
    will be returned instead.

"; 

%feature("docstring") simuPOP::InheritTagger "

Details:

    An inheritance tagger passes values of parental information
    field(s) to the corresponding fields of offspring. If there are
    two parental values from parents of a sexual mating event, a
    parameter mode is used to specify how to assign offspring
    information fields.

"; 

%feature("docstring") simuPOP::InheritTagger::InheritTagger "

Usage:

    InheritTagger(mode=PATERNAL, begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, output=\"\", infoFields=[])

Details:

    Creates an inheritance tagger that passes values of parental
    information fields (parameter infoFields) to the corresponding
    fields of offspring. If there is only one parent, values at the
    specified information fields are copied directly. If there are two
    parents, parameter mode specifies how to pass them to an
    offspring. More specifically,
    *   mode=MATERNAL Passing the value from mother.
    *   mode=PATERNAL Passing the value from father.
    *   mode=MEAN Passing the average of two values.
    *   mode=MAXIMUM Passing the maximum value of two values.
    *   mode=MINIMUM Passing the minimum value of two values.
    *   mode=SUMMATION Passing the summation of two values.
    *   mode=MULTIPLICATION Passing the multiplication of two values.
    An RuntimeError will be raised if any of the parents does not
    exist. This operator does not support parameter subPops and does
    not output any information.

"; 

%feature("docstring") simuPOP::InheritTagger::~InheritTagger "

Usage:

    x.~InheritTagger()

"; 

%feature("docstring") simuPOP::InheritTagger::describe "Obsolete or undocumented function."

%ignore simuPOP::InheritTagger::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::InheritTagger::clone "Obsolete or undocumented function."

%ignore simuPOP::InheritTagger::parallelizable() const;

%feature("docstring") simuPOP::InitGenotype "

Details:

    This operator assigns alleles at all or part of loci with given
    allele frequencies, proportions or values. This operator
    initializes all chromosomes, including unused genotype locations
    and customized chromosomes.

"; 

%feature("docstring") simuPOP::InitGenotype::InitGenotype "

Usage:

    InitGenotype(freq=[], genotype=[], prop=[], haplotypes=[],
      loci=ALL_AVAIL, ploidy=ALL_AVAIL, begin=0, end=1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    This function creates an initializer that initializes individual
    genotypes with random alleles or haplotypes with specified
    frequencies (parameter freq) or proportions (parameter prop). If
    parameter haplotypes is not specified, freq specifies the allele
    frequencies of alleles 0, 1, ... respectively. Alternatively, you
    can use parameter prop to specified the exact proportions of
    alleles 0, 1, ..., although alleles with small proportions might
    not be assigned at all. Values of parameter prob or prop should
    add up to 1. If parameter haplotypes is specified, it should
    contain a list of haplotypes and parameter prob or prop specifies
    frequencies or proportions of each haplotype. If loci, ploidy
    and/or subPop are specified, only specified loci, ploidy, and
    individuals in these (virtual) subpopulations will be initialized.
    Parameter loci can be a list of loci indexes, names or ALL_AVAIL.
    If the length of a haplotype is not enough to fill all loci, the
    haplotype will be reused. If a list (or a single) haplotypes are
    specified without freq or prop, they are used with equal
    probability.  In the last case, if a sequence of genotype is
    specified, it will be used repeatedly to initialize all alleles
    sequentially. This works similar to function
    Population.setGenotype() except that you can limit the
    initialization to certain loci and ploidy.

"; 

%feature("docstring") simuPOP::InitGenotype::~InitGenotype "

Usage:

    x.~InitGenotype()

"; 

%feature("docstring") simuPOP::InitGenotype::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::InitGenotype::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::InitGenotype::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::InitInfo "

Details:

    This operator initializes given information fields with a sequence
    of values, or a user-provided function such as random.random.

"; 

%feature("docstring") simuPOP::InitInfo::InitInfo "

Usage:

    InitInfo(values, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=[])

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

%feature("docstring") simuPOP::InitInfo::~InitInfo "

Description:

    destructor

Usage:

    x.~InitInfo()

"; 

%feature("docstring") simuPOP::InitInfo::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::InitInfo::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::InitInfo::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::InitLineage "

Details:

    This operator assigns lineages at all or part of loci with given
    values. This operator initializes all chromosomes, including
    unused lineage locations and customized chromosomes.

"; 

%feature("docstring") simuPOP::InitLineage::InitLineage "

Usage:

    InitLineage(lineage=[], mode=PER_LOCI, loci=ALL_AVAIL,
      ploidy=ALL_AVAIL, begin=0, end=1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=\"ind_id\")

Details:

    This function creates an initializer that initializes lineages
    with either a specified set of values or from the field infoFields
    (default to ind_id), whose value will be saved as the lineage of
    modified alleles. If a list of values is specified in parameter
    lineage, each value in this list is applied to one or more alleles
    so that each locus (PER_LOCI), alleles on each chromosome
    (PER_CHROMOSOME), on chromosomes of each ploidy (PER_PLOIDY), or
    for each individual (PER_INDIVIDUAL) have the same lineage. A
    single value is allowed and values in lineage will be re-used if
    not enough values are provided. If a valid field is specified
    (default to ind_id), the value of this field will be used for all
    alleles of each individual if mode is set to FROM_INFO, or be
    adjusted to produce positive values for alleles on the frist
    ploidy, and negative values for the second ploidy (and so on) if
    mode equals to FROM_INFO_SIGNED. If loci, ploidy and/or subPops
    are specified, only specified loci, ploidy, and individuals in
    these (virtual) subpopulations will be initialized.

"; 

%feature("docstring") simuPOP::InitLineage::~InitLineage "

Usage:

    x.~InitLineage()

"; 

%feature("docstring") simuPOP::InitLineage::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::InitLineage::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::InitLineage::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::InitSex "

Details:

    This operator initializes sex of individuals, either randomly or
    use a list of sexes.

"; 

%feature("docstring") simuPOP::InitSex::InitSex "

Usage:

    InitSex(maleFreq=0.5, maleProp=-1, sex=[], begin=0, end=-1,
      step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create an operator that initializes individual sex to MALE or
    FEMALE. By default, it assigns sex to individuals randomly, with
    equal probability of having a male or a female. This probabability
    can be adjusted through parameter maleFreq or be made to exact
    proportions by specifying parameter maleProp. Alternatively, a
    fixed sequence of sexes can be assigned. For example, if
    sex=[MALE, FEMALE], individuals will be assigned MALE and FEMALE
    successively. Parameter maleFreq or maleProp are ignored if sex is
    given. If a list of (virtual) subpopulation is specified in
    parameter subPop, only individuals in these subpopulations will be
    initialized. Note that the sex sequence, if used, is assigned
    repeatedly regardless of (virtual) subpopulation boundaries so
    that you can assign sex to all individuals in a population.

"; 

%feature("docstring") simuPOP::InitSex::~InitSex "

Description:

    destructor

Usage:

    x.~InitSex()

"; 

%feature("docstring") simuPOP::InitSex::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::InitSex::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::InitSex::apply "Obsolete or undocumented function."

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

%ignore simuPOP::intList::match(ssize_t rep, const vector< bool > *activeRep=NULL) const;

%feature("docstring") simuPOP::intMatrix "

"; 

%feature("docstring") simuPOP::intMatrix::intMatrix "

Usage:

    intMatrix(obj=None)

"; 

%ignore simuPOP::intMatrix::empty() const;

%ignore simuPOP::intMatrix::elems() const;

%feature("docstring") simuPOP::KAlleleMutator "

Details:

    This mutator implements a k-allele mutation model that assumes k
    allelic states (alleles 0, 1, 2, ..., k-1) at each locus. When a
    mutation event happens, it mutates an allele to any other states
    with equal probability.

"; 

%feature("docstring") simuPOP::KAlleleMutator::KAlleleMutator "

Usage:

    KAlleleMutator(k, rates=[], loci=ALL_AVAIL, mapIn=[], mapOut=[],
      output=\">\", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=\"ind_id\", lineageMode=FROM_INFO)

Details:

    Create a k-allele mutator that mutates alleles to one of the other
    k-1 alleles with equal probability. This mutator by default
    applies to all loci unless parameter loci is specified. A single
    mutation rate will be used for all loci if a single value of
    parameter rates is given. Otherwise, a list of mutation rates can
    be specified for each locus in parameter loci. If the mutated
    allele is larger than or equal to k, it will not be mutated. A
    warning message will be displayed if debugging code DBG_WARNING is
    turned on. Please refer to classes mutator and BaseOperator for
    descriptions of other parameters.

"; 

%feature("docstring") simuPOP::KAlleleMutator::~KAlleleMutator "

Usage:

    x.~KAlleleMutator()

"; 

%ignore simuPOP::KAlleleMutator::mutate(Allele allele, size_t locus) const;

%feature("docstring") simuPOP::KAlleleMutator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::KAlleleMutator::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::lociList "

"; 

%feature("docstring") simuPOP::lociList::lociList "

Usage:

    lociList(obj=Py_True)

"; 

%ignore simuPOP::lociList::lociList(const vectoru &values);

%feature("docstring") simuPOP::lociList::empty "

Usage:

    x.empty()

"; 

%ignore simuPOP::lociList::size() const;

%ignore simuPOP::lociList::name(size_t i) const;

%ignore simuPOP::lociList::allAvail() const;

%ignore simuPOP::lociList::unspecified() const;

%feature("docstring") simuPOP::lociList::dynamic "

Usage:

    x.dynamic()

"; 

%ignore simuPOP::lociList::elems(const GenoStruTrait *trait) const;

%feature("docstring") simuPOP::MaPenetrance "

Details:

    This operator is called a 'multi-allele' penetrance operator
    because it groups multiple alleles into two groups: wildtype and
    non-wildtype alleles. Alleles in each allele group are assumed to
    have the same effect on individual penetrance. If we denote all
    wildtype alleles as A, and all non-wildtype alleles a, this
    operator assign Individual penetrance according to genotype AA,
    Aa, aa in the diploid case, and A and a in the haploid case.

"; 

%feature("docstring") simuPOP::MaPenetrance::MaPenetrance "

Usage:

    MaPenetrance(loci, penetrance, wildtype=0, ancGens=UNSPECIFIED,
      begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=[])

Details:

    Creates a multi-allele penetrance operator that groups multiple
    alleles into a wildtype group (with alleles wildtype, default to
    [0]), and a non-wildtype group. A list of penetrance values is
    specified through parameter penetrance, for genotypes at one or
    more loci. Parameter loci can be a list of loci indexes, names or
    ALL_AVAIL. If we denote wildtype alleles using capital letters A,
    B ... and non-wildtype alleles using small letters a, b ..., the
    penetrance values should be for
    *   genotypes A and a for the haploid single-locus case,
    *   genotypes AB, Ab, aB and bb for haploid two=locus cases,
    *   genotypes AA, Aa and aa for diploid single-locus cases,
    *   genotypes AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, and
    aabb for diploid two-locus cases,
    *   and in general 2**n for diploid and 3**n for haploid cases if
    there are n loci. This operator does not support haplodiploid
    populations and sex chromosomes.

"; 

%feature("docstring") simuPOP::MaPenetrance::~MaPenetrance "

Usage:

    x.~MaPenetrance()

"; 

%feature("docstring") simuPOP::MaPenetrance::clone "Obsolete or undocumented function."

%ignore simuPOP::MaPenetrance::penet(Population *pop, RawIndIterator ind) const;

%feature("docstring") simuPOP::MaPenetrance::describe "Obsolete or undocumented function."

%ignore simuPOP::MaPenetrance::parallelizable() const;

%feature("docstring") simuPOP::MapPenetrance "

Details:

    This penetrance operator assigns individual affection status using
    a user-specified penetrance dictionary.

"; 

%feature("docstring") simuPOP::MapPenetrance::MapPenetrance "

Usage:

    MapPenetrance(loci, penetrance, ancGens=UNSPECIFIED, begin=0,
      end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL,
      infoFields=[])

Details:

    Create a penetrance operator that get penetrance value from a
    dictionary penetrance with genotype at loci as keys, and
    penetrance as values. For each individual, genotypes at loci are
    collected one by one (e.g. p0_loc0, p1_loc0, p0_loc1, p1_loc1...
    for a diploid individual) and are looked up in the dictionary.
    Parameter loci can be a list of loci indexes, names or ALL_AVAIL.
    If a genotype cannot be found, it will be looked up again without
    phase information (e.g. (1,0) will match key (0,1)). If the
    genotype still can not be found, a ValueError will be raised. This
    operator supports sex chromosomes and haplodiploid populations. In
    these cases, only valid genotypes should be used to generator the
    dictionary keys.

"; 

%feature("docstring") simuPOP::MapPenetrance::~MapPenetrance "

Usage:

    x.~MapPenetrance()

"; 

%feature("docstring") simuPOP::MapPenetrance::clone "Obsolete or undocumented function."

%ignore simuPOP::MapPenetrance::penet(Population *pop, RawIndIterator ind) const;

%feature("docstring") simuPOP::MapPenetrance::describe "Obsolete or undocumented function."

%ignore simuPOP::MapPenetrance::parallelizable() const;

%feature("docstring") simuPOP::MapSelector "

Details:

    This selector assigns individual fitness values using a user-
    specified dictionary. This operator can be applied to populations
    with arbitrary number of homologous chromosomes.

"; 

%feature("docstring") simuPOP::MapSelector::MapSelector "

Usage:

    MapSelector(loci, fitness, begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)

Details:

    Create a selector that assigns individual fitness values using a
    dictionary fitness with genotype at loci as keys, and fitness as
    values. Parameter loci can be a list of indexes, loci names or
    ALL_AVAIL. For each individual (parents if this operator is
    applied before mating, and offspring if this operator is applied
    during mating), genotypes at loci are collected one by one (e.g.
    p0_loc0, p1_loc0, p0_loc1, p1_loc1... for a diploid individual,
    with number of alleles varying for sex and mitochondrial DNAs) and
    are looked up in the dictionary. If a genotype cannot be found, it
    will be looked up again without phase information (e.g. (1,0) will
    match key (0,1)). If the genotype still can not be found, a
    ValueError will be raised. This operator supports sex chromosomes
    and haplodiploid populations. In these cases, only valid genotypes
    should be used to generator the dictionary keys.

"; 

%feature("docstring") simuPOP::MapSelector::~MapSelector "

Usage:

    x.~MapSelector()

"; 

%feature("docstring") simuPOP::MapSelector::clone "Obsolete or undocumented function."

%ignore simuPOP::MapSelector::indFitness(Population &pop, RawIndIterator ind) const;

%feature("docstring") simuPOP::MapSelector::describe "Obsolete or undocumented function."

%ignore simuPOP::MapSelector::parallelizable() const;

%feature("docstring") simuPOP::MaSelector "

Details:

    This operator is called a 'multi-allele' selector because it
    groups multiple alleles into two groups: wildtype and non-wildtype
    alleles. Alleles in each allele group are assumed to have the same
    effect on individual fitness. If we denote all wildtype alleles as
    A, and all non-wildtype alleles a, this operator assign individual
    fitness according to genotype AA, Aa, aa in the diploid case, and
    A and a in the haploid case.

"; 

%feature("docstring") simuPOP::MaSelector::MaSelector "

Usage:

    MaSelector(loci, fitness, wildtype=0, begin=0, end=-1, step=1,
      at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)

Details:

    Creates a multi-allele selector that groups multiple alleles into
    a wildtype group (with alleles wildtype, default to [0]), and a
    non-wildtype group. A list of fitness values is specified through
    parameter fitness, for genotypes at one or more loci. Parameter
    loci can be a list of indexes, loci names or ALL_AVAIL. If we
    denote wildtype alleles using capital letters A, B ... and non-
    wildtype alleles using small letters a, b ..., the fitness values
    should be for
    *   genotypes A and a for the haploid single-locus case,
    *   genotypes AB, Ab, aB and bb for haploid two=locus cases,
    *   genotypes AA, Aa and aa for diploid single-locus cases,
    *   genotypes AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, and
    aabb for diploid two-locus cases,
    *   and in general 2**n for diploid and 3**n for haploid cases if
    there are n loci. This operator does not support haplodiploid
    populations, sex and mitochondrial chromosomes.

"; 

%feature("docstring") simuPOP::MaSelector::~MaSelector "

Usage:

    x.~MaSelector()

"; 

%feature("docstring") simuPOP::MaSelector::clone "Obsolete or undocumented function."

%ignore simuPOP::MaSelector::indFitness(Population &pop, RawIndIterator ind) const;

%feature("docstring") simuPOP::MaSelector::describe "Obsolete or undocumented function."

%ignore simuPOP::MaSelector::parallelizable() const;

%feature("docstring") simuPOP::MatingScheme "

Details:

    This mating scheme is the base class of all mating schemes. It
    evolves a population generation by generation but does not
    actually transmit genotype.

"; 

%feature("docstring") simuPOP::MatingScheme::MatingScheme "

Usage:

    MatingScheme(subPopSize=[])

Details:

    Create a base mating scheme that evolves a population without
    transmitting genotypes. At each generation, this mating scheme
    creates an offspring generation according to parameter subPopSize,
    which can be a list of subpopulation sizes (or a number if there
    is only one subpopulation) or a Python function which will be
    called at each generation, just before mating, to determine the
    subpopulation sizes of the offspring generation. The function
    should be defined with one or both parameters of gen and pop where
    gen is the current generation number and pop is the parental
    population just before mating. The return value of this function
    should be a list of subpopulation sizes for the offspring
    generation. A single number can be returned if there is only one
    subpopulation. The passed parental population is usually used to
    determine offspring population size from parental population size
    but you can also modify this population to prepare for mating. A
    common practice is to split and merge parental populations in this
    function so that you demographic related information and actions
    could be implemented in the same function.

"; 

%feature("docstring") simuPOP::MatingScheme::~MatingScheme "

Description:

    destructor

Usage:

    x.~MatingScheme()

"; 

%feature("docstring") simuPOP::MatingScheme::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::MatingScheme::describe "Obsolete or undocumented function."

%ignore simuPOP::MatingScheme::submitScratch(Population &pop, Population &scratch);

%ignore simuPOP::MatingScheme::mateSubPop(Population &, Population &, size_t, RawIndIterator, RawIndIterator);

%ignore simuPOP::MatingScheme::mate(Population &pop, Population &scratch);

%ignore simuPOP::MatingScheme::prepareScratchPop(Population &pop, Population &scratch);

%ignore simuPOP::MatingScheme::subPopSizeSpecified();

%feature("docstring") simuPOP::MatrixMutator "

Details:

    A matrix mutator mutates alleles 0, 1, ..., n-1 using a n by n
    matrix, which specifies the probability at which each allele
    mutates to another. Conceptually speaking, this mutator goes
    through all mutable allele and mutate it to another state
    according to probabilities in the corresponding row of the rate
    matrix. Only one mutation rate matrix can be specified which will
    be used for all specified loci.

"; 

%feature("docstring") simuPOP::MatrixMutator::MatrixMutator "

Usage:

    MatrixMutator(rate, loci=ALL_AVAIL, mapIn=[], mapOut=[],
      output=\">\", begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=\"ind_id\", lineageMode=FROM_INFO)

Details:

    Create a mutator that mutates alleles 0, 1, ..., n-1 using a n by
    n matrix rate. Item (i,j) of this matrix specifies the probability
    at which allele i mutates to allele j. Diagnal items (i, i) are
    ignored because they are automatically determined by other
    probabilities. Only one mutation rate matrix can be specified
    which will be used for all loci in the applied population, or loci
    specified by parameter loci. If alleles other than 0, 1, ..., n-1
    exist in the population, they will not be mutated although a
    warning message will be given if debugging code DBG_WARNING is
    turned on. Please refer to classes mutator and BaseOperator for
    detailed explanation of other parameters.

"; 

%feature("docstring") simuPOP::MatrixMutator::~MatrixMutator "

Description:

    destructor.

Usage:

    x.~MatrixMutator()

"; 

%ignore simuPOP::MatrixMutator::mutate(Allele allele, size_t locus) const;

%feature("docstring") simuPOP::MatrixMutator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::MatrixMutator::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::MendelianGenoTransmitter "

Details:

    This Mendelian offspring generator accepts two parents and pass
    their genotypes to an offspring following Mendel's laws. Sex
    chromosomes are handled according to the sex of the offspring,
    which is usually determined in advance by an offspring generator.
    Customized chromosomes are not handled.

"; 

%feature("docstring") simuPOP::MendelianGenoTransmitter::MendelianGenoTransmitter "

Usage:

    MendelianGenoTransmitter(output=\"\", begin=0, end=-1, step=1,
      at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a Mendelian genotype transmitter (a during-mating operator)
    that transmits genotypes from parents to offspring following
    Mendel's laws. Autosomes and sex chromosomes are handled but
    customized chromosomes are ignored. Parameters subPops and
    infoFields are ignored. This operator also copies allelic lineage
    when it is executed in a module with lineage allele type.

"; 

%feature("docstring") simuPOP::MendelianGenoTransmitter::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::MendelianGenoTransmitter::describe "Obsolete or undocumented function."

%ignore simuPOP::MendelianGenoTransmitter::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::MendelianGenoTransmitter::initialize "Obsolete or undocumented function."

%feature("docstring") simuPOP::MendelianGenoTransmitter::transmitGenotype "

Usage:

    x.transmitGenotype(parent, offspring, ploidy)

Details:

    Transmit genotype from parent to offspring, and fill the ploidy
    homologous set of chromosomes. This function does not set
    genotypes of customized chromosomes and handles sex chromosomes
    properly, according to offspring sex and ploidy.

"; 

%ignore simuPOP::MendelianGenoTransmitter::parallelizable() const;

%feature("docstring") simuPOP::MergeSubPops "

Details:

    This operator merges subpopulations subPops to a single
    subpopulation. If subPops is ignored, all subpopulations will be
    merged. Virtual subpopulations are not allowed in subPops.

"; 

%feature("docstring") simuPOP::MergeSubPops::MergeSubPops "

Usage:

    MergeSubPops(subPops=ALL_AVAIL, name=\"\", begin=0, end=-1,
      step=1, at=[], reps=ALL_AVAIL, infoFields=[])

Details:

    Create an operator that merges subpopulations subPops to a single
    subpopulation. If subPops is not given, all subpopulations will be
    merged. The merged subpopulation will take the name of the first
    subpopulation being merged unless a new name is given.  This
    operator is by default applied pre-mating (parameter stage).
    Please refer to operator BaseOperator for a detailed explanation
    for all parameters.

"; 

%feature("docstring") simuPOP::MergeSubPops::~MergeSubPops "

Description:

    destructor

Usage:

    x.~MergeSubPops()

"; 

%feature("docstring") simuPOP::MergeSubPops::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::MergeSubPops::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::MergeSubPops::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::Migrator "

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
    information field. Unexpected results may happen if individuals
    migrate from overlapping virtual subpopulations.

"; 

%feature("docstring") simuPOP::Migrator::Migrator "

Usage:

    Migrator(rate=[], mode=BY_PROBABILITY, toSubPops=ALL_AVAIL,
      begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=\"migrate_to\")

Details:

    Create a Migrator that moves individuals from source (virtual)
    subpopulations subPops (default to migrate from all
    subpopulations) to destination subpopulations toSubPops (default
    to all subpopulations), according to existing values in an
    information field infoFields[0], or randomly according to a
    migration matrix rate. In the latter case, the size of the matrix
    should match the number of source and destination subpopulations.
    Depending on the value of parameter mode, elements in the
    migration matrix (rate) are interpreted as either the
    probabilities to migrate from source to destination subpopulations
    (mode = BY_PROBABILITY), proportions of individuals in the source
    (virtual) subpopulations to the destination subpopulations (mode =
    BY_PROPORTION), numbers of migrants in the source (virtual)
    subpopulations (mode = BY_COUNTS), or ignored completely (mode =
    BY_IND_INFO). In the last case, parameter subPops is respected
    (only individuals in specified (virtual) subpopulations will
    migrate) but toSubPops is ignored.  This operator is by default
    applied pre-mating (parameter stage). Please refer to operator
    BaseOperator for a detailed explanation for all parameters.

"; 

%feature("docstring") simuPOP::Migrator::~Migrator "

Description:

    destructor

Usage:

    x.~Migrator()

"; 

%feature("docstring") simuPOP::Migrator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::Migrator::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::Migrator::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::MitochondrialGenoTransmitter "

Details:

    This geno transmitter transmits the first homologous copy of a
    Mitochondrial chromosome. If no mitochondrial chromosome is
    present, it assumes that the first homologous copy of several (or
    all) Customized chromosomes are copies of mitochondrial
    chromosomes. This operator transmits the mitochondrial chromosome
    from the female parent to offspring for sexsual reproduction, and
    any parent to offspring for asexual reproduction. If there are
    multiple chromosomes, the organelles are selected randomly. If
    this transmitter is applied to populations with more than one
    homologous copies of chromosomes, it transmits the first
    homologous copy of chromosomes and clears alleles (set to zero) on
    other homologous copies.

"; 

%feature("docstring") simuPOP::MitochondrialGenoTransmitter::MitochondrialGenoTransmitter "

Usage:

    MitochondrialGenoTransmitter(output=\"\", chroms=ALL_AVAIL,
      begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=[])

Details:

    Createa a mitochondrial genotype transmitter that treats the
    Mitochondiral chromosome, or Customized chromosomes if no
    Mitochondrial chromosome is specified, or a list of chromosomes
    specified by chroms, as human mitochondrial chromosomes. These
    chromosomes should have the same length and the same number of
    loci. This operator transmits these chromosomes randomly from the
    female parent to offspring of both sexes. It also copies allelic
    lineage when it is executed in a module with lineage allele type.

"; 

%feature("docstring") simuPOP::MitochondrialGenoTransmitter::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::MitochondrialGenoTransmitter::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::MitochondrialGenoTransmitter::initialize "Obsolete or undocumented function."

%ignore simuPOP::MitochondrialGenoTransmitter::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%ignore simuPOP::MitochondrialGenoTransmitter::parallelizable() const;

%feature("docstring") simuPOP::MixedMutator "

Details:

    This mixed mutator accepts a list of mutators and use one of them
    to mutate an allele when an mutation event happens.

"; 

%feature("docstring") simuPOP::MixedMutator::MixedMutator "

Usage:

    MixedMutator(rates=[], loci=ALL_AVAIL, mutators=[], prob=[],
      mapIn=[], mapOut=[], context=0, output=\">\", begin=0, end=-1,
      step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL,
      infoFields=\"ind_id\", lineageMode=FROM_INFO)

Details:

    Create a mutator that randomly chooses one of the specified
    mutators to mutate an allele when a mutation event happens. The
    mutators are choosen according to a list of probabilities ( prob)
    that should add up to 1. The passed and returned alleles might be
    changed if parameters mapIn and mapOut are used. Most parameters,
    including loci, mapIn, mapOut, rep, and subPops of mutators
    specified in parameter mutators are ignored. This mutator by
    default applies to all loci unless parameter loci is specified.
    Please refer to classes mutator and BaseOperator for descriptions
    of other parameters.

"; 

%feature("docstring") simuPOP::MixedMutator::clone "Obsolete or undocumented function."

%ignore simuPOP::MixedMutator::mutate(Allele allele, size_t locus) const;

%feature("docstring") simuPOP::MixedMutator::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::MlPenetrance "

Details:

    This penetrance operator is created by a list of penetrance
    operators. When it is applied to an individual, it applies these
    penetrance operators to the individual, obtain a list of
    penetrance values, and compute a combined penetrance value from
    them and assign affection status accordingly. ADDITIVE,
    multiplicative, and a heterogeneour multi-locus model are
    supported. Please refer to Neil Rish (1989) \"Linkage Strategies
    for  Genetically Complex Traits\" for some analysis of these
    models.

"; 

%feature("docstring") simuPOP::MlPenetrance::MlPenetrance "

Usage:

    MlPenetrance(ops, mode=MULTIPLICATIVE, ancGens=UNSPECIFIED,
      begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a multiple-locus penetrance operator from a list penetrance
    operator ops. When this operator is applied to an individual
    (parents when used before mating and offspring when used during
    mating), it applies these operators to the individual and obtain a
    list of (usually single-locus) penetrance values. These penetrance
    values are combined to a single penetrance value using
    *   Prod(f_i), namely the product of individual penetrance if mode
    = MULTIPLICATIVE,
    *   sum(f_i) if mode = ADDITIVE, and
    *   1-Prod(1 - f_i) if mode = HETEROGENEITY 0 or 1 will be
    returned if the combined penetrance value is less than zero or
    greater than 1.  Applicability parameters (begin, end, step, at,
    reps, subPops) could be used in both MlSelector and selectors in
    parameter ops, but parameters in MlSelector will be interpreted
    first.

"; 

%feature("docstring") simuPOP::MlPenetrance::~MlPenetrance "

Usage:

    x.~MlPenetrance()

"; 

%feature("docstring") simuPOP::MlPenetrance::clone "Obsolete or undocumented function."

%ignore simuPOP::MlPenetrance::penet(Population *pop, RawIndIterator ind) const;

%feature("docstring") simuPOP::MlPenetrance::describe "Obsolete or undocumented function."

%ignore simuPOP::MlPenetrance::parallelizable() const;

%feature("docstring") simuPOP::MlSelector "

Details:

    This selector is created by a list of selectors. When it is
    applied to an individual, it applies these selectors to the
    individual, obtain a list of fitness values, and compute a
    combined fitness value from them. ADDITIVE, multiplicative, and a
    heterogeneour multi-locus model are supported.

"; 

%feature("docstring") simuPOP::MlSelector::MlSelector "

Usage:

    MlSelector(ops, mode=MULTIPLICATIVE, begin=0, end=-1, step=1,
      at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)

Details:

    Create a multiple-locus selector from a list selection operator
    selectors. When this operator is applied to an individual (parents
    when used before mating and offspring when used during mating), it
    applies these operators to the individual and obtain a list of
    (usually single-locus) fitness values. These fitness values are
    combined to a single fitness value using
    *   Prod(f_i), namely the product of individual fitness if mode =
    MULTIPLICATIVE,
    *   1-sum(1 - f_i) if mode = ADDITIVE,
    *   1-Prod(1 - f_i) if mode = HETEROGENEITY, and
    *   exp(- sum(1 - f_i)) if mode = EXPONENTIAL, zero will be
    returned if the combined fitness value is less than zero.
    Applicability parameters (begin, end, step, at, reps, subPops)
    could be used in both MlSelector and selectors in parameter ops,
    but parameters in MlSelector will be interpreted first.

"; 

%feature("docstring") simuPOP::MlSelector::~MlSelector "

Usage:

    x.~MlSelector()

"; 

%feature("docstring") simuPOP::MlSelector::clone "Obsolete or undocumented function."

%ignore simuPOP::MlSelector::indFitness(Population &pop, RawIndIterator ind) const;

%feature("docstring") simuPOP::MlSelector::describe "Obsolete or undocumented function."

%ignore simuPOP::MlSelector::parallelizable() const;

%feature("docstring") simuPOP::NoneOp "

Details:

    This operator does nothing when it is applied to a population. It
    is usually used as a placeholder when an operator is needed
    syntactically.

"; 

%feature("docstring") simuPOP::NoneOp::NoneOp "

Usage:

    NoneOp(output=\">\", begin=0, end=0, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a NoneOp.

"; 

%feature("docstring") simuPOP::NoneOp::~NoneOp "

Description:

    destructor

Usage:

    x.~NoneOp()

"; 

%feature("docstring") simuPOP::NoneOp::clone "Obsolete or undocumented function."

%ignore simuPOP::NoneOp::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::NoneOp::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::NoneOp::describe "Obsolete or undocumented function."

%ignore simuPOP::NoSexModel;

%feature("docstring") simuPOP::NoSexModel::NoSexModel "

Usage:

    NoSexModel()

"; 

%feature("docstring") simuPOP::NoSexModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::NoSexModel::getSex "

Usage:

    x.getSex(UINT)

"; 

%feature("docstring") simuPOP::NoSexModel::parallelizable "

Usage:

    x.parallelizable()

"; 

%ignore simuPOP::NumOfFemalesSexModel;

%feature("docstring") simuPOP::NumOfFemalesSexModel::NumOfFemalesSexModel "

Usage:

    NumOfFemalesSexModel(numOfFemales)

"; 

%feature("docstring") simuPOP::NumOfFemalesSexModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::NumOfFemalesSexModel::getSex "

Usage:

    x.getSex(count)

"; 

%feature("docstring") simuPOP::NumOfFemalesSexModel::parallelizable "

Usage:

    x.parallelizable()

"; 

%ignore simuPOP::NumOffModel;

%feature("docstring") simuPOP::NumOffModel::NumOffModel "

Usage:

    NumOffModel()

"; 

%feature("docstring") simuPOP::NumOffModel::~NumOffModel "

Usage:

    x.~NumOffModel()

"; 

%feature("docstring") simuPOP::NumOffModel::getNumOff "

Usage:

    x.getNumOff(gen)

"; 

%feature("docstring") simuPOP::NumOffModel::reset "

Usage:

    x.reset()

"; 

%feature("docstring") simuPOP::NumOffModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::NumOffModel::parallelizable "

Usage:

    x.parallelizable()

"; 

%ignore simuPOP::NumOfMalesSexModel;

%feature("docstring") simuPOP::NumOfMalesSexModel::NumOfMalesSexModel "

Usage:

    NumOfMalesSexModel(numOfMales)

"; 

%feature("docstring") simuPOP::NumOfMalesSexModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::NumOfMalesSexModel::getSex "

Usage:

    x.getSex(count)

"; 

%feature("docstring") simuPOP::NumOfMalesSexModel::parallelizable "

Usage:

    x.parallelizable()

"; 

%feature("docstring") simuPOP::OffspringGenerator "

Details:

    An offspring generator generates offspring from parents chosen by
    a parent chooser. It is responsible for creating a certain number
    of offspring, determinning their sex, and transmitting genotypes
    from parents to offspring.

"; 

%feature("docstring") simuPOP::OffspringGenerator::OffspringGenerator "

Usage:

    OffspringGenerator(ops, numOffspring=1, sexMode=RANDOM_SEX)

Details:

    Create a basic offspring generator. This offspring generator uses
    ops genotype transmitters to transmit genotypes from parents to
    offspring.  A number of during-mating operators (parameter ops)
    can be used to, among other possible duties such as setting
    information fields of offspring, transmit genotype from parents to
    offspring. This general offspring generator does not have any
    default during-mating operator but all stock mating schemes use an
    offspring generator with a default operator. For example, a
    mendelianOffspringGenerator is used by RandomMating to trasmit
    genotypes. Note that applicability parameters begin, step, end, at
    and reps could be used in these operators but negative population
    and generation indexes are unsupported.  Parameter numOffspring is
    used to control the number of offspring per mating event, or in
    another word the number of offspring in each family. It can be a
    number, a Python function or generator, or a mode parameter
    followed by some optional arguments. If a number is given, given
    number of offspring will be generated at each mating event. If a
    Python function is given, it will be called each time when a
    mating event happens. When a generator function is specified, it
    will be called for each subpopulation to provide number of
    offspring for all mating events during the populating of this
    subpopulation. Current generation number will be passed to this
    function or generator function if parameter \"gen\" is used in this
    function. In the last case, a tuple (or a list) in one of the
    following forms can be given:
    *   (GEOMETRIC_DISTRIBUTION, p)
    *   (POISSON_DISTRIBUTION, p), p > 0
    *   (BINOMIAL_DISTRIBUTION, p, N), 0 < p <=1, N > 0
    *   (UNIFORM_DISTRIBUTION, a, b), 0 <= a <= b. In this case, the
    number of offspring will be determined randomly following the
    specified statistical distributions. Because families with zero
    offspring are silently ignored, the distribution of the observed
    number of offspring per mating event (excluding zero) follows
    zero-truncated versions of these distributions.  Parameter
    numOffspring specifies the number of offspring per mating event
    but the actual surviving offspring can be less than specified.
    More spefically, if any during-mating operator returns False, an
    offspring will be discarded so the actually number of offspring of
    a mating event will be reduced. This is essentially how during-
    mating selector works.  Parameter sexMode is used to control the
    sex of each offspring. Its default value is usually RANDOM_SEX
    which assign MALE or FEMALE to each individual randomly, with
    equal probabilities. If NO_SEX is given, offspring sex will not be
    changed. sexMode can also be one of
    *   (PROB_OF_MALES, p) where p is the probability of male for each
    offspring,
    *   (NUM_OF_MALES, n) where n is the number of males in a mating
    event. If n is greater than or equal to the number of offspring in
    this family, all offspring in this family will be MALE.
    *   (NUM_OF_FEMALES, n) where n is the number of females in a
    mating event,
    *   (SEQUENCE_OF_SEX, s1, s2 ...) where s1, s2 etc are MALE or
    FEMALE. The sequence will be used for each mating event. It will
    be reused if the number of offspring in a mating event is greater
    than the length of sequence.
    *   (GLOBAL_SEQUENCE_OF_SEX, s1, s2, ...) where s1, s2 etc are
    MALE or FEMALE. The sequence will be used across mating events. It
    will be reused if the number of offspring in a subpopulation is
    greater than the length of sequence. Finally, parameter sexMode
    accepts a function or a generator function. A function will be
    called whenever an offspring is produced. A generator will be
    created at each subpopulation and will be used to produce sex for
    all offspring in this subpopulation. No parameter is accepted.

"; 

%ignore simuPOP::OffspringGenerator::OffspringGenerator(const OffspringGenerator &rhs);

%feature("docstring") simuPOP::OffspringGenerator::~OffspringGenerator "

Usage:

    x.~OffspringGenerator()

"; 

%feature("docstring") simuPOP::OffspringGenerator::clone "Obsolete or undocumented function."

%ignore simuPOP::OffspringGenerator::initialize(const Population &pop, size_t subPop);

%ignore simuPOP::OffspringGenerator::generateOffspring(Population &pop, Population &offPop, Individual *dad, Individual *mom, RawIndIterator &offBegin, RawIndIterator &offEnd);

%ignore simuPOP::OffspringGenerator::finalize(const Population &);

%feature("docstring") simuPOP::OffspringGenerator::describe "Obsolete or undocumented function."

%ignore simuPOP::OffspringGenerator::initialized();

%ignore simuPOP::OffspringGenerator::numOffspring(ssize_t gen);

%ignore simuPOP::OffspringGenerator::getSex(UINT count);

%ignore simuPOP::OffspringGenerator::parallelizable() const;

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

%ignore simuPOP::opList::begin() const;

%ignore simuPOP::opList::end() const;

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

%feature("docstring") simuPOP::ParentChooser "

Details:

    A parent chooser repeatedly chooses parent(s) from a parental
    population and pass them to an offspring generator. A parent
    chooser can select one or two parents, which should be matched by
    the offspring generator. This class is the base class of all
    parent choosers, and should not be used directly.

"; 

%feature("docstring") simuPOP::ParentChooser::ParentChooser "

Usage:

    ParentChooser(selectionField=\"\")

"; 

%feature("docstring") simuPOP::ParentChooser::clone "Obsolete or undocumented function."

%ignore simuPOP::ParentChooser::initialize(Population &, size_t);

%ignore simuPOP::ParentChooser::finalize(Population &, size_t);

%feature("docstring") simuPOP::ParentChooser::describe "Obsolete or undocumented function."

%ignore simuPOP::ParentChooser::parallelizable() const;

%ignore simuPOP::ParentChooser::initialized() const;

%ignore simuPOP::ParentChooser::chooseParents(RawIndIterator);

%feature("docstring") simuPOP::ParentChooser::~ParentChooser "

Description:

    destructor

Usage:

    x.~ParentChooser()

"; 

%feature("docstring") simuPOP::ParentsTagger "

Details:

    This tagging operator records the indexes of parents (relative to
    the parental generation) of each offspring in specified
    information fields ( default to father_idx and mother_idx). Only
    one information field should be specified if an asexsual mating
    scheme is used so there is one parent for each offspring.
    Information recorded by this operator is intended to be used to
    look up parents of each individual in multi-generational
    Population.

"; 

%feature("docstring") simuPOP::ParentsTagger::ParentsTagger "

Usage:

    ParentsTagger(begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, output=\"\", infoFields=[\"father_idx\",
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

%feature("docstring") simuPOP::ParentsTagger::~ParentsTagger "

Usage:

    x.~ParentsTagger()

"; 

%feature("docstring") simuPOP::ParentsTagger::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::ParentsTagger::describe "Obsolete or undocumented function."

%ignore simuPOP::ParentsTagger::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%ignore simuPOP::ParentsTagger::parallelizable() const;

%feature("docstring") simuPOP::Pause "

Details:

    This operator pauses the evolution of a simulator at given
    generations or at a key stroke. When a simulator is stopped, you
    can go to a Python shell to examine the status of an evolutionary
    process, resume or stop the evolution.

"; 

%feature("docstring") simuPOP::Pause::Pause "

Usage:

    Pause(stopOnKeyStroke=False, prompt=True, output=\">\", begin=0,
      end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL,
      infoFields=[])

Details:

    Create an operator that pause the evolution of a population when
    it is applied to this population. If stopOnKeyStroke is False
    (default), it will always pause a population when it is applied,
    if this parameter is set to True, the operator will pause a
    population if any key has been pressed. If a specific character is
    set, the operator will stop when this key has been pressed. This
    allows, for example, the use of several pause operators to pause
    different populations.  After a population has been paused, a
    message will be displayed (unless prompt is set to False) and
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

%feature("docstring") simuPOP::Pause::~Pause "

Description:

    destructor

Usage:

    x.~Pause()

"; 

%feature("docstring") simuPOP::Pause::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pause::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pause::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pedigree "

Details:

    The pedigree class is derived from the population class. Unlike a
    population class that emphasizes on individual properties, the
    pedigree class emphasizes on relationship between individuals. An
    unique ID for all individuals is needed to create a pedigree
    object from a population object. Compared to the Population class,
    a Pedigree object is optimized for access individuals by their
    IDs, regardless of population structure and ancestral generations.
    Note that the implementation of some algorithms rely on the fact
    that parental IDs are smaller than their offspring because
    individual IDs are assigned sequentially during evolution.
    Pedigrees with manually assigned IDs should try to obey such a
    rule.

"; 

%feature("docstring") simuPOP::Pedigree::Pedigree "

Usage:

    Pedigree(pop, loci=[], infoFields=[], ancGens=ALL_AVAIL,
      idField=\"ind_id\", fatherField=\"father_id\",
      motherField=\"mother_id\", stealPop=False)

Details:

    Create a pedigree object from a population, using a subset of loci
    (parameter loci, can be a list of loci indexes, names, or
    ALL_AVAIL, default to no locus), information fields (parameter
    infoFields, default to no information field besides idField,
    fatherField and motherField), and ancestral generations (parameter
    ancGens, default to all ancestral generations). By default,
    information field father_id (parameter fatherField) and mother_id
    (parameter motherField) are used to locate parents identified by
    ind_id (parameter idField), which should store an unique ID for
    all individuals. Multiple individuls with the same ID are allowed
    and will be considered as the same individual, but a warning will
    be given if they actually differ in genotype or information
    fields. Operators IdTagger and PedigreeTagger are usually used to
    assign such IDs, although function sampling.indexToID could be
    used to assign unique IDs and construct parental IDs from index
    based relationship recorded by operator ParentsTagger. A pedigree
    object could be constructed with one or no parent but certain
    functions such as relative tracking will not be available for such
    pedigrees. In case that your are no longer using your population
    object, you could steal the content from the population by setting
    stealPop to True.

"; 

%ignore simuPOP::Pedigree::Pedigree(const Pedigree &rhs);

%ignore simuPOP::Pedigree::idIdx() const;

%ignore simuPOP::Pedigree::fatherOf(size_t id) const;

%ignore simuPOP::Pedigree::motherOf(size_t id) const;

%feature("docstring") simuPOP::Pedigree::clone "

Usage:

    x.clone()

Details:

    Create a cloned copy of a Pedigree.

"; 

%feature("docstring") simuPOP::Pedigree::save "

Usage:

    x.save(filename, infoFields=[], loci=[])

Details:

    Save a pedigree to file filename. This function goes through all
    individuals of a pedigree and outputs in each line the ID of
    individual, IDs of his or her parents, sex ('M' or 'F'), affection
    status ('A' or 'U'), values of specified information fields
    infoFields and genotypes at specified loci (parameter loci, which
    can be a list of loci indexes, names, or ALL_AVAIL). Allele
    numbers, instead of their names are outputed. Two columns are used
    for each locus if the population is diploid. This file can be
    loaded using function loadPedigree although additional information
    such as names of information fields need to be specified. This
    format differs from a .ped file used in some genetic analysis
    software in that there is no family ID and IDs of all individuals
    have to be unique. Note that parental IDs will be set to zero if
    the parent is not in the pedigree object. Therefore, the parents
    of individuals in the top-most ancestral generation will always be
    zero.

"; 

%feature("docstring") simuPOP::Pedigree::indByID "

Usage:

    x.indByID(id)

Details:

    Return a reference to individual with id. An IndexError will be
    raised if no individual with id is found. An float id is
    acceptable as long as it rounds closely to an integer.

"; 

%ignore simuPOP::Pedigree::indByID(size_t id) const;

%feature("docstring") simuPOP::Pedigree::numParents "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pedigree::locateRelatives "

Usage:

    x.locateRelatives(relType, resultFields=[], sex=ANY_SEX,
      affectionStatus=ANY_AFFECTION_STATUS, ancGens=ALL_AVAIL)

Details:

    This function locates relatives (of type relType) of each
    individual and store their IDs in information fields relFields.
    The length of relFields determines how many relatives an
    individual can have.  Parameter relType specifies what type of
    relative to locate, which can be
    *   SPOUSE locate spouses with whom an individual has at least one
    common offspring.
    *   OUTBRED_SPOUSE locate non-slibling spouses, namely spouses
    with no shared parent.
    *   OFFSPRING all offspring of each individual.
    *   COMMON_OFFSPRING common offspring between each individual and
    its spouse (located by SPOUSE or OUTBRED_SPOUSE). relFields should
    consist of an information field for spouse and m-1 fields for
    offspring where m is the number of fields.
    *   FULLSIBLING siblings with common father and mother,
    *   SIBLING siblings with at least one common parent. Optionally,
    you can specify the sex and affection status of relatives you
    would like to locate, using parameters sex and affectionStatus.
    sex can be ANY_SEX (default), MALE_ONLY, FEMALE_ONLY, SAME_SEX or
    OPPOSITE_SEX, and affectionStatus can be AFFECTED, UNAFFECTED or
    ANY_AFFECTION_STATUS (default). Only relatives with specified
    properties will be located.  This function will by default go
    through all ancestral generations and locate relatives for all
    individuals. This can be changed by setting parameter ancGens to
    certain ancestral generations you would like to process.

"; 

%feature("docstring") simuPOP::Pedigree::traceRelatives "

Usage:

    x.traceRelatives(fieldPath, sex=[], affectionStatus=[],
      resultFields=[], ancGens=ALL_AVAIL)

Details:

    Trace a relative path in a population and record the result in the
    given information fields resultFields. This function is used to
    locate more distant relatives based on the relatives located by
    function locateRelatives. For example, after siblings and
    offspring of all individuals are located, you can locate mother's
    sibling's offspring using a relative path, and save their indexes
    in each individuals information fields resultFields.  A relative
    path consits of a fieldPath that specifies which information
    fields to look for at each step, a sex specifies sex choices at
    each generation, and a affectionStatus that specifies affection
    status at each generation. fieldPath should be a list of
    information fields, sex and affectionStatus are optional. If
    specified, they should be a list of ANY_SEX, MALE_ONLY,
    FEMALE_ONLY, SAME_SEX and OppsiteSex for parameter sex, and a list
    of UNAFFECTED, AFFECTED and ANY_AFFECTION_STATUS for parameter
    affectionStatus.  For example, if fieldPath = [['father_id',
    'mother_id'], ['sib1', 'sib2'], ['off1', 'off2']], and sex =
    [ANY_SEX, MALE_ONLY, FEMALE_ONLY], this function will locate
    father_id and mother_id for each individual, find all individuals
    referred by father_id and mother_id, find informaton fields sib1
    and sib2 from these parents and locate male individuals referred
    by these two information fields. Finally, the information fields
    off1 and off2 from these siblings are located and are used to
    locate their female offspring. The results are father or mother's
    brother's daughters. Their indexes will be saved in each
    individuals information fields resultFields. If a list of
    ancestral generations is given in parameter ancGens is given, only
    individuals in these ancestral generations will be processed.

"; 

%feature("docstring") simuPOP::Pedigree::individualsWithRelatives "

Usage:

    x.individualsWithRelatives(infoFields, sex=[],
      affectionStatus=[], subPops=ALL_AVAIL, ancGens=ALL_AVAIL)

Details:

    Return a list of IDs of individuals who have non-negative values
    at information fields infoFields. Additional requirements could be
    specified by parameters sex and affectionStatus. sex can be
    ANY_SEX (default), MALE_ONLY, FEMALE_ONLY, SAME_SEX or
    OPPOSITE_SEX, and affectionStatus can be AFFECTED, UNAFFECTED or
    ANY_AFFECTION_STATUS (default). This function by default check all
    individuals in all ancestral generations, but you could limit the
    search using parameter subPops (a list of (virtual)
    subpopulations) and ancestral generations ancGens. Relatives fall
    out of specified subpopulations and ancestral generaions will be
    considered invalid.

"; 

%feature("docstring") simuPOP::Pedigree::identifyFamilies "

Usage:

    x.identifyFamilies(pedField=\"\", subPops=ALL_AVAIL,
      ancGens=ALL_AVAIL)

Details:

    This function goes through all individuals in a pedigree and group
    related individuals into families. If an information field
    pedField is given, indexes of families will be assigned to this
    field of each family member. The return value is a list of family
    sizes corresponding to families 0, 1, 2, ... etc. If a list of
    (virtual) subpopulations (parameter subPops) or ancestral
    generations are specified (parameter ancGens), the search will be
    limited to individuals in these subpopulations and generations.

"; 

%feature("docstring") simuPOP::Pedigree::identifyAncestors "

Usage:

    x.identifyAncestors(IDs=ALL_AVAIL, subPops=ALL_AVAIL,
      ancGens=ALL_AVAIL)

Details:

    If a list of individuals (IDs) is given, this function traces
    backward in time and find all ancestors of these individuals. If
    IDs is ALL_AVAIL, ancestors of all individuals in the present
    generation will be located. If a list of (virtual) subpopulations
    (subPops) or ancestral geneartions (ancGens) is given, the search
    will be limited to individuals in these subpopulations and
    generations. This could be used to, for example, find all fathers
    of IDs. This function returns a list of IDs, which includes valid
    specified IDs. Invalid IDs will be silently ignored. Note that
    parameters subPops and ancGens will limit starting IDs if IDs is
    set to ALL_AVAIL, but specified IDs will not be trimmed according
    to these parameters.

"; 

%feature("docstring") simuPOP::Pedigree::identifyOffspring "

Usage:

    x.identifyOffspring(IDs=[], subPops=ALL_AVAIL,
      ancGens=ALL_AVAIL)

Details:

    This function traces forward in time and find all offspring of
    individuals specified in parameter IDs. If a list of (virtual)
    subpopulations (subPops) or ancestral geneartions (ancGens) is
    given, the search will be limited to individuals in these
    subpopulations and generations. This could be used to, for
    example, find all male offspring of IDs. This function returns a
    list of IDs, which includes valid starting IDs. Invalid IDs are
    silently ignored. Note that parameters subPops and ancGens will
    limit search result but will not be used to trim specified IDs.

"; 

%feature("docstring") simuPOP::Pedigree::removeIndividuals "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pedigree::removeSubPops "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pedigree::push "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pedigree::addChrom "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pedigree::addChromFrom "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pedigree::addIndFrom "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pedigree::mergeSubPops "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pedigree::resize "Obsolete or undocumented function."

%feature("docstring") simuPOP::Pedigree::setSubPopByIndInfo "Obsolete or undocumented function."

%feature("docstring") simuPOP::PedigreeMating "

Details:

    This mating scheme evolves a population following an existing
    pedigree structure. If the Pedigree object has N ancestral
    generations and a present generation, it can be used to evolve a
    population for N generations, starting from the topmost ancestral
    generation. At the k-th generation, this mating scheme produces an
    offspring generation according to subpopulation structure of the
    N-k-1 ancestral generation in the pedigree object (e.g. producing
    the offspring population of generation 0 according to the N-1
    ancestral generation of the pedigree object ). For each offspring,
    this mating scheme copies individual ID and sex from the
    corresponing individual in the pedigree object. It then locates
    the parents of each offspring using their IDs in the pedigree
    object. A list of during mating operators are then used to
    transmit parental genotype to the offspring. The population being
    evolved must have an information field 'ind_id'.

"; 

%feature("docstring") simuPOP::PedigreeMating::PedigreeMating "

Usage:

    PedigreeMating(ped, ops, idField=\"ind_id\")

Details:

    Creates a pedigree mating scheme that evolves a population
    according to Pedigree object ped. The evolved population should
    contain individuals with ID (at information field idField, default
    to 'ind_id') that match those individual in the topmost ancestral
    generation who have offspring. After parents of each individuals
    are determined from their IDs, a list of during-mating operators
    ops are applied to transmit genotypes. The return value of these
    operators are not checked.

"; 

%feature("docstring") simuPOP::PedigreeMating::~PedigreeMating "

Usage:

    x.~PedigreeMating()

"; 

%feature("docstring") simuPOP::PedigreeMating::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::PedigreeMating::describe "Obsolete or undocumented function."

%ignore simuPOP::PedigreeMating::mate(Population &pop, Population &scratch);

%feature("docstring") simuPOP::PedigreeMating::parallelizable "

Usage:

    x.parallelizable()

"; 

%feature("docstring") simuPOP::PedigreeTagger "

Details:

    This tagging operator records the ID of parents of each offspring
    in specified information fields (default to father_id and
    mother_id). Only one information field should be specified if an
    asexsual mating scheme is used so there is one parent for each
    offspring. Information recorded by this operator is intended to be
    used to record full pedigree information of an evolutionary
    process.

"; 

%feature("docstring") simuPOP::PedigreeTagger::PedigreeTagger "

Usage:

    PedigreeTagger(idField=\"ind_id\", output=\"\", outputFields=[],
      outputLoci=[], begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=[\"father_id\", \"mother_id\"])

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
    by default does not send any output. If a valid output stream is
    given (should be in the form of '>>filename' so that output will
    be concatenated), this operator will output the ID of offspring,
    IDs of his or her parent(s), sex and affection status of
    offspring, and values at specified information fields
    (outputFields) and loci (outputLoci) in the format of off_id
    father_id mother_id M/F A/U fields genotype. father_id or
    mother_id will be ignored if only one parent is involved. This
    file format can be loaded using function loadPedigree.  Because
    only offspring will be outputed, individuals in the top-most
    ancestral generation will not be outputed. This is usually not a
    problem because individuals who have offspring in the next
    generation will be constructed by function loadPedigree, although
    their information fields and genotype will be missing. If you
    would like to create a file with complete pedigree information,
    you can apply this operator before evolution in the initOps
    parameter of functions Population.evolve or Simulator.evolve. This
    will output all individuals in the initial population (the top-
    most ancestral population after evolution) in the same format.
    Note that sex, affection status and genotype can be changed by
    other operators so this operator should usually be applied after
    all other operators are applied.

"; 

%feature("docstring") simuPOP::PedigreeTagger::~PedigreeTagger "

Usage:

    x.~PedigreeTagger()

"; 

%feature("docstring") simuPOP::PedigreeTagger::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::PedigreeTagger::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::PedigreeTagger::apply "Obsolete or undocumented function."

%ignore simuPOP::PedigreeTagger::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%ignore simuPOP::PedigreeTagger::parallelizable() const;

%feature("docstring") simuPOP::PointMutator "

Details:

    A point mutator is different from all other mutators because
    mutations in this mutator do not happen randomly. Instead, it
    happens to specific loci and mutate an allele to a specific state,
    regardless of its original state. This mutator is usually used to
    introduce a mutant to a population.

"; 

%feature("docstring") simuPOP::PointMutator::PointMutator "

Usage:

    PointMutator(loci, allele, ploidy=[], 0, inds=[], output=\">\",
      begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=0,
      infoFields=\"ind_id\", lineageMode=FROM_INFO)

Details:

    Create a point mutator that mutates alleles at specified loci to a
    given allele of individuals inds. If there are multiple alleles at
    a locus (e.g. individuals in a diploid population), only the first
    allele is mutated unless indexes of alleles are listed in
    parameter ploidy. This operator is by default applied to
    individuals in the first subpopulation but you can apply it to a
    different or more than one (virtual) subpopulations using
    parameter subPops (AllAvail is also accepted). Please refer to
    class BaseOperator for detailed descriptions of other parameters.

"; 

%feature("docstring") simuPOP::PointMutator::~PointMutator "

Description:

    destructor

Usage:

    x.~PointMutator()

"; 

%feature("docstring") simuPOP::PointMutator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::PointMutator::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::PointMutator::describe "Obsolete or undocumented function."

%ignore simuPOP::PoissonNumOffModel;

%feature("docstring") simuPOP::PoissonNumOffModel::PoissonNumOffModel "

Usage:

    PoissonNumOffModel(mu)

"; 

%feature("docstring") simuPOP::PoissonNumOffModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::PoissonNumOffModel::getNumOff "

Usage:

    x.getNumOff(ssize_t)

"; 

%feature("docstring") simuPOP::PolyParentsChooser "

Details:

    This parent chooser is similar to random parents chooser but
    instead of selecting a new pair of parents each time, one of the
    parents in this parent chooser will mate with several spouses
    before he/she is replaced. This mimicks multi-spouse mating
    schemes such as polygyny or polyandry in some populations. Natural
    selection is supported for both sexes.

"; 

%feature("docstring") simuPOP::PolyParentsChooser::PolyParentsChooser "

Usage:

    PolyParentsChooser(polySex=MALE, polyNum=1,
      selectionField=\"fitness\")

Details:

    Create a multi-spouse parents chooser where each father (if
    polySex is MALE) or mother (if polySex is FEMALE) has polyNum
    spouses. The parents are chosen with replacement. If individual
    fitness values are assigned (stored to information field
    selectionField, default to \"fitness\"), the probability that an
    individual is chosen is proportional to his/her fitness value
    among all individuals with the same sex.

"; 

%feature("docstring") simuPOP::PolyParentsChooser::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::PolyParentsChooser::describe "Obsolete or undocumented function."

%ignore simuPOP::PolyParentsChooser::parallelizable() const;

%ignore simuPOP::PolyParentsChooser::initialize(Population &pop, size_t sp);

%ignore simuPOP::PolyParentsChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::Population "

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
    to class Individual for an introduction to genotype arragement of
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

%feature("docstring") simuPOP::Population::Population "

Usage:

    Population(size=[], ploidy=2, loci=[], chromTypes=[],
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
                    second case, for which ploidy=HAPLODIPLOID should
                    be used.
    loci:           A list of numbers of loci on each chromosome. The
                    length of this parameter determines the number of
                    chromosomes. If there is only one chromosome,
                    numLoci instead of [numLoci] can be used.
    chromTypes:     A list that specifies the type of each chromosome,
                    which can be AUTOSOME, CHROMOSOME_X, CHROMOSOME_Y,
                    or CUSTOMIZED. All chromosomes are assumed to be
                    autosomes if this parameter is ignored. Sex
                    chromosome can only be specified in a diploid
                    population where the sex of an individual is
                    determined by the existence of these chromosomes
                    using the XX (FEMALE) and XY (MALE) convention.
                    Both sex chromosomes have to be available and be
                    specified only once. Because chromosomes X and Y
                    are treated as two chromosomes, recombination on
                    the pseudo-autosomal regions of the sex chromsomes
                    is not supported. CUSTOMIZED chromosomes are
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

%ignore simuPOP::Population::Population(const Population &rhs);

%feature("docstring") simuPOP::Population::clone "

Usage:

    x.clone()

Details:

    Create a cloned copy of a population. Note that Python statement
    pop1 = pop only creates a reference to an existing population pop.

"; 

%feature("docstring") simuPOP::Population::swap "

Usage:

    x.swap(rhs)

Details:

    Swap the content of two population objects, which can be handy in
    some particular circumstances. For example, you could swap out a
    population in a simulator.

"; 

%feature("docstring") simuPOP::Population::~Population "

Description:

    destroy a population

Usage:

    x.~Population()

"; 

%ignore simuPOP::Population::validate(const string &msg) const;

%ignore simuPOP::Population::fitSubPopStru(const vectoru &newSubPopSizes, const vectorstr &newSubPopNames);

%ignore simuPOP::Population::hasActivatedVirtualSubPop() const;

%ignore simuPOP::Population::hasActivatedVirtualSubPop(size_t subPop) const;

%ignore simuPOP::Population::hasVirtualSubPop() const;

%feature("docstring") simuPOP::Population::virtualSplitter "

Usage:

    x.virtualSplitter()

Details:

    Return the virtual splitter associated with the population, None
    will be returned if there is no splitter.

"; 

%feature("docstring") simuPOP::Population::setVirtualSplitter "

Usage:

    x.setVirtualSplitter(splitter)

Details:

    Set a VSP splitter to the population, which defines the same VSPs
    for all subpopulations. If different VSPs are needed for different
    subpopulations, a CombinedSplitter can be used to make these VSPs
    available to all subpopulations.

"; 

%feature("docstring") simuPOP::Population::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of virtual subpopulations (VSP) defined by a VSP
    splitter. Return 0 if no VSP is defined.

"; 

%feature("docstring") simuPOP::Population::activateVirtualSubPop "Obsolete or undocumented function."

%feature("docstring") simuPOP::Population::deactivateVirtualSubPop "Obsolete or undocumented function."

%feature("docstring") simuPOP::Population::cmp "

Description:

    a python function used to compare the population objects

Usage:

    x.__cmp__(rhs)

"; 

%feature("docstring") simuPOP::Population::fitGenoStru "Obsolete or undocumented function."

%feature("docstring") simuPOP::Population::setSubPopStru "Obsolete or undocumented function."

%feature("docstring") simuPOP::Population::numSubPop "

Usage:

    x.numSubPop()

Details:

    Return the number of subpopulations in a population. Return 1 if
    there is no subpopulation structure.

"; 

%feature("docstring") simuPOP::Population::subPopSize "

Usage:

    x.subPopSize(subPop=[], ancGen=-1)

Details:

    Return the size of a subpopulation (subPopSize(sp)) or a virtual
    subpopulation (subPopSize([sp, vsp])) in the current generation
    (default) or a specified ancestral generation ancGen. If no subpop
    is given, it is the same as popSize(ancGen). Population and
    virtual subpopulation names can be used.
    <group>2-subpopsize</grouplociList()>

"; 

%feature("docstring") simuPOP::Population::subPopByName "

Usage:

    x.subPopByName(name)

Details:

    Return the index of the first subpopulation with name name. An
    IndexError will be raised if subpopulations are not named, or if
    no subpopulation with name name is found. Virtual subpopulation
    name is not supported.

"; 

%feature("docstring") simuPOP::Population::subPopName "

Usage:

    x.subPopName(subPop)

Details:

    Return the \"spName - vspName\" (virtual named subpopulation), \"\"
    (unnamed non-virtual subpopulation), \"spName\" (named
    subpopulation) or \"vspName\" (unnamed virtual subpopulation),
    depending on whether subpopulation is named or if subPop is
    virtual.

"; 

%feature("docstring") simuPOP::Population::subPopNames "

Usage:

    x.subPopNames()

Details:

    Return the names of all subpopulations (excluding virtual
    subpopulations). An empty string will be returned for unnamed
    subpopulations.

"; 

%feature("docstring") simuPOP::Population::setSubPopName "

Usage:

    x.setSubPopName(name, subPop)

Details:

    Assign a name name to subpopulation subPop. Note that
    subpopulation names do not have to be unique.

"; 

%feature("docstring") simuPOP::Population::subPopSizes "

Usage:

    x.subPopSizes(ancGen=-1)

Details:

    Return the sizes of all subpopulations at the current generation
    (default) or specified ancestral generation ancGen. Virtual
    subpopulations are not considered.

"; 

%feature("docstring") simuPOP::Population::popSize "

Usage:

    x.popSize(ancGen=-1)

Details:

    Return the total number of individuals in all subpopulations of
    the current generation (default) or the an ancestral generation
    ancGen.

"; 

%feature("docstring") simuPOP::Population::absIndIndex "

Usage:

    x.absIndIndex(idx, subPop)

Details:

    return the absolute index of an individual idx in subpopulation
    subPop.

"; 

%feature("docstring") simuPOP::Population::subPopIndPair "

Usage:

    x.subPopIndPair(idx)

Details:

    Return the subpopulation ID and relative index of an individual,
    given its absolute index idx.

"; 

%feature("docstring") simuPOP::Population::subPopBegin "

Usage:

    x.subPopBegin(subPop)

Details:

    Return the index of the first individual in subpopulation subPop.

"; 

%feature("docstring") simuPOP::Population::subPopEnd "

Usage:

    x.subPopEnd(subPop)

Details:

    Return the index of the last individual in subpopulation subPop
    plus 1, so that range(subPopBegin(subPop), subPopEnd(subPop) can
    iterate through the index of all individuals in subpopulation
    subPop.

"; 

%ignore simuPOP::Population::individual(size_t idx, vspID subPop=vspID());

%feature("docstring") simuPOP::Population::individual "

Usage:

    x.individual(idx, subPop=[])

Details:

    Return a refernce to individual idx in the population (if
    subPop=[], default) or a subpopulation (if subPop=sp). Virtual
    subpopulation is not supported. Note that a float idx is
    acceptable as long as it rounds closely to an integer.

"; 

%feature("docstring") simuPOP::Population::indByID "

Usage:

    x.indByID(id, ancGens=ALL_AVAIL, idField=\"ind_id\")

Details:

    Return a reference to individual with id stored in information
    field idField (default to ind_id). This function by default search
    the present and all ancestral generations (ancGen=ALL_AVAIL), but
    you can limit the search in specific generations if you know which
    generations to search (ancGens=[0,1] for present and parental
    generations) or UNSPECIFIED to search only the current generation.
    If no individual with id is found, an IndexError will be raised. A
    float id is acceptable as long as it rounds closely to an integer.
    Note that this function uses a dynamic searching algorithm which
    tends to be slow. If you need to look for multiple individuals
    from a static population, you might want to convert a population
    object to a pedigree object and use function Pedigree.indByID.

"; 

%ignore simuPOP::Population::individual(double idx, vspID subPop=vspID()) const;

%feature("docstring") simuPOP::Population::ancestor "

Usage:

    x.ancestor(idx, gen, subPop=[])

Details:

    Return a reference to individual idx in ancestral generation gen.
    The correct individual will be returned even if the current
    generation is not the present one (see also useAncestralGen). If a
    valid subPop is specified, index is relative to that subPop.
    Virtual subpopulation is not supported. Note that a float idx is
    acceptable as long as it rounds closely to an integer.

"; 

%ignore simuPOP::Population::ancestor(double idx, ssize_t gen, vspID subPop=vspID()) const;

%feature("docstring") simuPOP::Population::individuals "

Usage:

    x.individuals(subPop=[])

Details:

    Return an iterator that can be used to iterate through all
    individuals in a population (if subPop=[], default), or a
    (virtual) subpopulation (subPop=spID or (spID, vspID)). If you
    would like to iterate through multiple subpopulations in multiple
    ancestral generations, please use function
    Population.allIndividuals().

"; 

%ignore simuPOP::Population::indOrdered() const;

%ignore simuPOP::Population::setIndOrdered(bool s) const;

%ignore simuPOP::Population::indIterator();

%ignore simuPOP::Population::indIterator(size_t subPop);

%ignore simuPOP::Population::indIterator() const;

%ignore simuPOP::Population::rawIndBegin();

%ignore simuPOP::Population::rawIndEnd();

%ignore simuPOP::Population::rawIndBegin(size_t subPop);

%ignore simuPOP::Population::rawIndEnd(size_t subPop);

%ignore simuPOP::Population::rawIndBegin() const;

%ignore simuPOP::Population::alleleIterator(size_t locus);

%ignore simuPOP::Population::alleleIterator(size_t locus, size_t subPop);

%ignore simuPOP::Population::alleleIterator(size_t locus) const;

%ignore simuPOP::Population::genoBegin(bool order);

%ignore simuPOP::Population::genoEnd(bool order);

%ignore simuPOP::Population::genoBegin(size_t subPop, bool order);

%ignore simuPOP::Population::genoEnd(size_t subPop, bool order);

%ignore simuPOP::Population::indGenoBegin(size_t ind) const;

%ignore simuPOP::Population::indGenoEnd(size_t ind) const;

%feature("docstring") simuPOP::Population::genotype "

Usage:

    x.genotype(subPop=[])

Details:

    Return an editable array of the genotype of all individuals in a
    population (if subPop=[], default), or individuals in a
    subpopulation subPop. Virtual subpopulation is unsupported.

"; 

%feature("docstring") simuPOP::Population::mutants "

Usage:

    x.mutants(subPop=[])

Details:

    Return an iterator that iterate through mutants of all individuals
    in a population (if subPop=[], default), or individuals in a
    subpopulation subPop. Virtual subpopulation is unsupported. Each
    mutant is presented as a tuple of (index, value) where index is
    the index of mutant (from 0 to totNumLoci()*ploidy()) so you will
    have to adjust its value to check multiple alleles at a locus.
    This function ignores type of chromosomes so non-zero alleles in
    unused alleles of sex and mitochondrial chromosomes are also
    iterated.

"; 

%feature("docstring") simuPOP::Population::lineage "

Usage:

    x.lineage(subPop=[])

Details:

    Return an editable array of the lineage of alleles for all
    individuals in a population (if subPop=[], default), or
    individuals in a subpopulation subPop. Virtual subpopulation is
    unsupported. This function returns None for modules without
    lineage information.

"; 

%feature("docstring") simuPOP::Population::setGenotype "

Usage:

    x.setGenotype(geno, subPop=[])

Details:

    Fill the genotype of all individuals in a population (if
    subPop=[]) or in a (virtual) subpopulation subPop (if subPop=sp or
    (sp, vsp)) using a list of alleles geno. geno will be reused if
    its length is less than subPopSize(subPop)*totNumLoci()*ploidy().

"; 

%feature("docstring") simuPOP::Population::setLineage "

Usage:

    x.setLineage(geno, subPop=[])

Details:

    Fill the lineage of all individuals in a population (if subPop=[])
    or in a (virtual) subpopulation subPop (if subPop=sp or (sp, vsp))
    using a list of IDs lineage. lineage will be reused if its length
    is less than subPopSize(subPop)*totNumLoci()*ploidy(). This
    function returns directly for modules without lineage information.

"; 

%feature("docstring") simuPOP::Population::sortIndividuals "

Usage:

    x.sortIndividuals(infoFields)

Details:

    Sort individuals according to values at specified information
    fields (infoFields). Individuals will be sorted at an increasing
    order.

"; 

%feature("docstring") simuPOP::Population::setSubPopByIndInfo "

Usage:

    x.setSubPopByIndInfo(field)

Details:

    Rearrange individuals to their new subpopulations according to
    their integer values at information field field (value returned by
    Individual::info(field)). individuals with negative values at this
    field will be removed. Existing subpopulation names are kept. New
    subpopulations will have empty names.

"; 

%feature("docstring") simuPOP::Population::splitSubPop "

Usage:

    x.splitSubPop(subPop, sizes, names=[])

Details:

    Split subpopulation subPop into subpopulations of given sizes,
    which should add up to the size of subpopulation subPop or 1, in
    which case sizes are treated as proportions. If subPop is not the
    last subpopulation, indexes of subpopulations after subPop are
    shifted. If subPop is named, the same name will be given to all
    new subpopulations unless a new set of names are specified for
    these subpopulations. This function returns the IDs of split
    subpopulations.

"; 

%feature("docstring") simuPOP::Population::removeSubPops "

Usage:

    x.removeSubPops(subPops)

Details:

    Remove (virtual) subpopulation(s) subPops and all their
    individuals. This function can be used to remove complete
    subpopulations (with shifted subpopulation indexes) or individuals
    belonging to virtual subpopulations of a subpopulation. In the
    latter case, the subpopulations are kept even if all individuals
    have been removed. This function only handles the present
    generation.

"; 

%ignore simuPOP::Population::removeMarkedIndividuals();

%feature("docstring") simuPOP::Population::removeIndividuals "

Usage:

    x.removeIndividuals(indexes=[], IDs=[], idField=\"ind_id\",
      filter=None)

Details:

    remove individual(s) by absolute indexes (parameter index) or
    their IDs (parameter IDs), or using a filter function (paramter
    filter). If indexes are used, only individuals at the current
    generation will be removed. If IDs are used, all individuals with
    one of the IDs at information field idField (default to \"ind_id\")
    will be removed. Although \"ind_id\" usually stores unique IDs of
    individuals, this function is frequently used to remove groups of
    individuals with the same value at an information field. An
    IndexError will be raised if an index is out of bound, but no
    error will be given if an invalid ID is specified. In the last
    case, a user-defined function should be provided. This function
    should accept parameter \"ind\" or one or more of the information
    fields. All individuals, including ancestors if there are multiple
    ancestral generations, will be passed to this function.
    Individuals that returns True will be removed. This function does
    not affect subpopulation structure in the sense that a
    subpopulation will be kept even if all individuals from it are
    removed.

"; 

%feature("docstring") simuPOP::Population::mergeSubPops "

Usage:

    x.mergeSubPops(subPops=ALL_AVAIL, name=\"\")

Details:

    Merge subpopulations subPops. If subPops is ALL_AVAIL (default),
    all subpopulations will be merged. subPops do not have to be
    adjacent to each other. They will all be merged to the
    subpopulation with the smallest subpopulation ID. Indexes of the
    rest of the subpopulation may be changed. A new name can be
    assigned to the merged subpopulation through parameter name (an
    empty name will be ignored). This function returns the ID of the
    merged subpopulation.

"; 

%feature("docstring") simuPOP::Population::addIndFrom "

Usage:

    x.addIndFrom(pop)

Details:

    Add all individuals, including ancestors, in pop to the current
    population. Two populations should have the same genotypic
    structures and number of ancestral generations. Subpopulations in
    population pop are kept.

"; 

%feature("docstring") simuPOP::Population::addChromFrom "

Usage:

    x.addChromFrom(pop)

Details:

    Add chromosomes in population pop to the current population.
    population pop should have the same number of individuals as the
    current population in the current and all ancestral generations.
    This function merges genotypes on the new chromosomes from
    population pop individual by individual.

"; 

%feature("docstring") simuPOP::Population::addLociFrom "

Usage:

    x.addLociFrom(pop)

Details:

    Add loci from population pop, chromosome by chromosome. Added loci
    will be inserted according to their position. Their position and
    names should not overlap with any locus in the current population.
    population pop should have the same number of individuals as the
    current population in the current and all ancestral generations.
    Allele lineages are also copied from pop in modules with lineage
    information.

"; 

%feature("docstring") simuPOP::Population::addChrom "

Usage:

    x.addChrom(lociPos, lociNames=[], chromName=\"\", alleleNames=[],
      chromType=AUTOSOME)

Details:

    Add chromosome chromName with given type chromType to a
    population, with loci lociNames inserted at position lociPos.
    lociPos should be ordered. lociNames and chromName should not
    exist in the current population. Allele names could be specified
    for all loci (a list of names) or differently for each locus (a
    nested list of names), using parameter alleleNames. Empty loci
    names will be used if lociNames is not specified. The newly added
    alleles will have zero lineage in modules wiht lineage
    information.

"; 

%feature("docstring") simuPOP::Population::addLoci "

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
    function returns indexes of the inserted loci. Newly inserted
    alleles will have zero lineage in modules with lineage
    information.

"; 

%feature("docstring") simuPOP::Population::resize "

Usage:

    x.resize(sizes, propagate=False)

Details:

    Resize population by giving new subpopulation sizes sizes.
    individuals at the end of some subpopulations will be removed if
    the new subpopulation size is smaller than the old one. New
    individuals will be appended to a subpopulation if the new size is
    larger. Their genotypes will be set to zero (default), or be
    copied from existing individuals if propagate is set to True. More
    specifically, if a subpopulation with 3 individuals is expanded to
    7, the added individuals will copy genotypes from individual 1, 2,
    3, and 1 respectively. Note that this function only resizes the
    current generation.

"; 

%feature("docstring") simuPOP::Population::extractSubPops "

Usage:

    x.extractSubPops(subPops=ALL_AVAIL, rearrange=False)

Details:

    Extract a list of (virtual) subpopulations from a population and
    create a new population. If rearrange is False (default),
    structure and names of extracted subpopulations are kept although
    extracted subpopulations can have fewer individuals if they are
    created from extracted virtual subpopulations. (e.g. it is
    possible to extract all male individuals from a subpopulation
    using a SexSplitter()). If rearrange is True, each (virtual)
    subpopulation in subPops becomes a new subpopulation in the
    extracted population in the order at which they are specified.
    Because each virtual subpopulation becomes a subpopulation, this
    function could be used, for example, to separate male and female
    individuals to two subpopulations ( subPops=[(0,0), (0,1)]). If
    overlapping (virtual) subpopulations are specified, individuals
    will be copied multiple times. This function only extract
    individuals from the present generation.

"; 

%ignore simuPOP::Population::extractMarkedIndividuals() const;

%feature("docstring") simuPOP::Population::extractIndividuals "

Usage:

    x.extractIndividuals(indexes=[], IDs=[], idField=\"ind_id\",
      filter=None)

Details:

    Extract individuals with given absolute indexes (parameter
    indexes), IDs (parameter IDs, stored in information field idField,
    default to ind_id), or a filter function (parameter filter). If a
    list of absolute indexes are specified, the present generation
    will be extracted and form a one-generational population. If a
    list of IDs are specified, this function will look through all
    ancestral generations and extract individuals with given ID.
    Individuals with shared IDs are allowed. In the last case, a user-
    defined Python function should be provided. This function should
    accept parameter \"ind\" or one or more of the information fields.
    All individuals, including ancestors if there are multiple
    ancestral generations, will be passed to this function.
    Individuals that returns True will be extracted. Extracted
    individuals will be in their original ancestral generations and
    subpopulations, even if some subpopulations or generations are
    empty. An IndexError will be raised if an index is out of bound
    but no error will be given if an invalid ID is encountered.

"; 

%ignore simuPOP::Population::extract(const lociList &extractedLoci, const stringList &infoFieldList, const subPopList &subPops=subPopList(), const uintList &ancGens=uintList()) const;

%feature("docstring") simuPOP::Population::removeLoci "

Usage:

    x.removeLoci(loci=UNSPECIFIED, keep=UNSPECIFIED)

Details:

    Remove loci (absolute indexes or names) and genotypes at these
    loci from the current population. Alternatively, a parameter keep
    can be used to specify loci that will not be removed.

"; 

%feature("docstring") simuPOP::Population::recodeAlleles "

Usage:

    x.recodeAlleles(alleles, loci=ALL_AVAIL, alleleNames=[])

Details:

    Recode alleles at loci (can be a list of loci indexes or names, or
    all loci in a population (ALL_AVAIL)) to other values according to
    parameter alleles. This parameter can a list of new allele numbers
    for alleles 0, 1, 2, ... (allele x will be recoded to
    newAlleles[x], x outside of the range of newAlleles will not be
    recoded, although a warning will be given if DBG_WARNING is
    defined) or a Python function, which should accept one or both
    parameters allele (existing allele) and locus (index of locus).
    The return value will become the new allele. This function is
    intended to recode some alleles without listing all alleles in a
    list. It will be called once for each existing allele so it is not
    possible to recode an allele to different alleles. A new list of
    allele names could be specified for these loci. Different sets of
    names could be specified for each locus if a nested list of names
    are given. This function recode alleles for all subpopulations in
    all ancestral generations.

"; 

%feature("docstring") simuPOP::Population::push "

Usage:

    x.push(pop)

Details:

    Push population pop into the current population. Both populations
    should have the same genotypic structure. The current population
    is discarded if ancestralDepth (maximum number of ancestral
    generations to hold) is zero so no ancestral generation can be
    kept. Otherise, the current population will become the parental
    generation of pop. If ancGen of a population is positive and there
    are already ancGen ancestral generations (c.f. ancestralGens()),
    the greatest ancestral generation will be discarded. In any case,
    Populationpop becomes invalid as all its individuals are absorbed
    by the current population.

"; 

%feature("docstring") simuPOP::Population::curAncestralGen "Obsolete or undocumented function."

%feature("docstring") simuPOP::Population::ancestralGens "

Usage:

    x.ancestralGens()

Details:

    Return the actual number of ancestral generations stored in a
    population, which does not necessarily equal to the number set by
    setAncestralDepth().

"; 

%ignore simuPOP::Population::clearInfo();

%ignore simuPOP::Population::markIndividuals(vspID subPop, bool mark) const;

%feature("docstring") simuPOP::Population::setIndInfo "

Usage:

    x.setIndInfo(values, field, subPop=[])

Details:

    Set information field field (specified by index or name) of all
    individuals (if subPop=[], default), or individuals in a (virtual)
    subpopulation (subPop=sp or (sp, vsp)) to values. values will be
    reused if its length is smaller than the size of the population or
    (virtual) subpopulation.

"; 

%ignore simuPOP::Population::infoBegin(size_t idx);

%ignore simuPOP::Population::infoEnd(size_t idx);

%feature("docstring") simuPOP::Population::indInfo "

Usage:

    x.indInfo(field, subPop=[])

Details:

    Return the values (as a list) of information field field (by index
    or name) of all individuals (if subPop=[], default), or
    individuals in a (virtual) subpopulation (if subPop=sp or (sp,
    vsp)).

"; 

%feature("docstring") simuPOP::Population::addInfoFields "

Usage:

    x.addInfoFields(fields, init=0)

Details:

    Add a list of information fields fields to a population and
    initialize their values to init. If an information field alreay
    exists, it will be re-initialized.

"; 

%feature("docstring") simuPOP::Population::setInfoFields "

Usage:

    x.setInfoFields(fields, init=0)

Details:

    Set information fields fields to a population and initialize them
    with value init. All existing information fields will be removed.

"; 

%feature("docstring") simuPOP::Population::removeInfoFields "

Usage:

    x.removeInfoFields(fields)

Details:

    Remove information fields fields from a population.

"; 

%feature("docstring") simuPOP::Population::updateInfoFieldsFrom "

Usage:

    x.updateInfoFieldsFrom(fields, pop, fromFields=[],
      ancGens=ALL_AVAIL)

Details:

    Update information fields fields from fromFields of another
    population (or Pedigree) pop. Two populations should have the same
    number of individuals. If fromFields is not specified, it is
    assumed to be the same as fields. If ancGens is not ALL_AVAIL,
    only the specified ancestral generations are updated.

"; 

%feature("docstring") simuPOP::Population::setAncestralDepth "

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

%ignore simuPOP::Population::keepAncestralGens(const uintList &ancGens);

%feature("docstring") simuPOP::Population::useAncestralGen "

Usage:

    x.useAncestralGen(idx)

Details:

    Making ancestral generation idx (0 for current generation, 1 for
    parental generation, 2 for grand-parental generation, etc) the
    current generation. This is an efficient way to access Population
    properties of an ancestral generation. useAncestralGen(0) should
    always be called afterward to restore the correct order of
    ancestral generations.

"; 

%ignore simuPOP::Population::syncIndPointers(bool infoOnly=false) const;

%feature("docstring") simuPOP::Population::save "

Usage:

    x.save(filename)

Details:

    Save population to a file filename, which can be loaded by a
    global function loadPopulation(filename).

"; 

%ignore simuPOP::Population::load(const string &filename);

%feature("docstring") simuPOP::Population::vars "

Usage:

    x.vars(subPop=[])

Details:

    return variables of a population as a Python dictionary. If a
    valid subpopulation subPop is specified, a dictionary
    vars()[\"subPop\"][subPop] is returned. A ValueError will be raised
    if key subPop does not exist in vars(), or if key subPop does not
    exist in vars()[\"subPop\"].

"; 

%ignore simuPOP::Population::dict(vspID subPop=vspID());

%ignore simuPOP::Population::getVars() const;

%ignore simuPOP::Population::setDict(PyObject *dict);

%ignore simuPOP::Population::varsAsString() const;

%ignore simuPOP::Population::varsFromString(const string &vars);

%feature("docstring") simuPOP::Population::evaluate "Obsolete or undocumented function."

%feature("docstring") simuPOP::Population::execute "Obsolete or undocumented function."

%ignore simuPOP::Population::rep();

%ignore simuPOP::Population::setRep(size_t rep);

%ignore simuPOP::Population::gen() const;

%ignore simuPOP::Population::setGen(size_t gen);

%ignore simuPOP::ProbOfMalesSexModel;

%feature("docstring") simuPOP::ProbOfMalesSexModel::ProbOfMalesSexModel "

Usage:

    ProbOfMalesSexModel(probOfMales)

"; 

%feature("docstring") simuPOP::ProbOfMalesSexModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::ProbOfMalesSexModel::getSex "

Usage:

    x.getSex(UINT)

"; 

%feature("docstring") simuPOP::ProbOfMalesSexModel::parallelizable "

Usage:

    x.parallelizable()

"; 

%feature("docstring") simuPOP::ProductSplitter "

Details:

    This splitter takes several splitters and take their intersections
    as new VSPs. For example, if the first splitter defines 3 VSPs and
    the second splitter defines 2, 6 VSPs will be defined by splitting
    3 VSPs defined by the first splitter each to two VSPs. This
    splitter is usually used to define finer VSPs from existing VSPs.

"; 

%feature("docstring") simuPOP::ProductSplitter::ProductSplitter "

Usage:

    ProductSplitter(splitters=[], names=[])

Details:

    Create a product splitter using a list of splitters. For example,
    ProductSplitter([SexSplitter(), AffectionSplitter()]) defines four
    VSPs by male unaffected, male affected, female unaffected, and
    female affected individuals. VSP names are usually determined by
    splitters, but can also be specified using parameter names.

"; 

%ignore simuPOP::ProductSplitter::ProductSplitter(const ProductSplitter &rhs);

%feature("docstring") simuPOP::ProductSplitter::~ProductSplitter "

Usage:

    x.~ProductSplitter()

"; 

%feature("docstring") simuPOP::ProductSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::ProductSplitter::size(const Population &pop, size_t subPop, size_t virtualSubPop) const;

%feature("docstring") simuPOP::ProductSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs defined by this splitter, which is the
    sum of the number of VSPs of all combined splitters.

"; 

%ignore simuPOP::ProductSplitter::contains(const Population &pop, size_t ind, vspID vsp) const;

%ignore simuPOP::ProductSplitter::activate(const Population &pop, size_t subPop, size_t virtualSubPop);

%feature("docstring") simuPOP::ProductSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of a VSP vsp, which is the names of indivdual VSPs
    separated by a comma, unless a new set of names is specified for
    each VSP.

"; 

%feature("docstring") simuPOP::ProportionSplitter "

Details:

    This splitter divides subpopulations into several VSPs by
    proportion.

"; 

%feature("docstring") simuPOP::ProportionSplitter::ProportionSplitter "

Usage:

    ProportionSplitter(proportions=[], names=[])

Details:

    Create a splitter that divides subpopulations by proportions,
    which should be a list of float numbers (between 0 and 1) that add
    up to 1. A default set of names are given to each VSP unless a new
    set of names is given by parameter names.

"; 

%feature("docstring") simuPOP::ProportionSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::ProportionSplitter::size(const Population &pop, size_t subPop, size_t virtualSubPop) const;

%feature("docstring") simuPOP::ProportionSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs defined by this splitter, which is the
    length of parameter proportions.

"; 

%ignore simuPOP::ProportionSplitter::contains(const Population &pop, size_t ind, vspID vsp) const;

%ignore simuPOP::ProportionSplitter::activate(const Population &pop, size_t subPop, size_t virtualSubPop);

%feature("docstring") simuPOP::ProportionSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of VSP vsp, which is \"Prop p\" where
    p=propotions[vsp]. A user specified name will be returned if
    specified.

"; 

%feature("docstring") simuPOP::PyEval "

Details:

    A PyEval operator evaluates a Python expression in a population's
    local namespace when it is applied to this population. The result
    is written to an output specified by parameter output.

"; 

%feature("docstring") simuPOP::PyEval::PyEval "

Usage:

    PyEval(expr=\"\", stmts=\"\", exposePop=\"\", output=\">\", begin=0,
      end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=Py_False,
      infoFields=[])

Details:

    Create a PyEval operator that evaluates a Python expression expr
    in a population's local namespaces when it is applied to this
    population. This namespace can either be the population's local
    namespace (pop.vars()), or namespaces subPop[sp] for (virtual)
    subpop (pop.vars(subpop)) in specified subPops. If Python
    statements stmts is given (a single or multi-line string), the
    statement will be executed before expr. If exposePop is set to an
    non-empty string, the current population will be exposed in its
    own local namespace as a variable with this name. This allows the
    execution of expressions such as 'pop.individual(0).allele(0)'.
    The result of expr will be sent to an output stream specified by
    parameter output. The exposed population variable will be removed
    after expr is evaluated. Please refer to class BaseOperator for
    other parameters.

Note:

    Although the statements and expressions are evaluated in a
    population's local namespace, they have access to a global
    namespace which is the module global namespace. It is therefore
    possible to refer to any module variable in these expressions.
    Such mixed use of local and global variables is, however, strongly
    discouraged.

"; 

%feature("docstring") simuPOP::PyEval::~PyEval "

Usage:

    x.~PyEval()

"; 

%feature("docstring") simuPOP::PyEval::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyEval::evaluate "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyEval::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyEval::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyExec "

Details:

    This operator executes given Python statements in a population's
    local namespace when it is applied to this population.

"; 

%feature("docstring") simuPOP::PyExec::PyExec "

Usage:

    PyExec(stmts=\"\", exposePop=\"\", output=\">\", begin=0, end=-1,
      step=1, at=[], reps=ALL_AVAIL, subPops=Py_False, infoFields=[])

Details:

    Create a PyExec operator that executes statements stmts in a
    population's local namespace when it is applied to this
    population. This namespace can either be the population's local
    namespace (pop.vars()), or namespaces subPop[sp] for each
    (virtual) subpop (pop.vars(subpop)) in specified subPops. If
    exposePop is given, current population will be exposed in its
    local namespace as a variable named by exposePop. Although
    multiple statements can be executed, it is recommended that you
    use this operator to execute short statements and use PyOperator
    for more complex once. Note that exposed population variables will
    be removed after the statements are executed.

"; 

%feature("docstring") simuPOP::PyExec::~PyExec "

Usage:

    x.~PyExec()

"; 

%feature("docstring") simuPOP::PyExec::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyExec::describe "Obsolete or undocumented function."

%ignore simuPOP::pyFunc;

%feature("docstring") simuPOP::pyFunc::pyFunc "

Usage:

    pyFunc(func)

"; 

%feature("docstring") simuPOP::pyFunc::numArgs "

Description:

    return number of arguments this function accepts. This function
    does not count tuple parameters.

Usage:

    x.numArgs()

"; 

%feature("docstring") simuPOP::pyFunc::name "

Usage:

    x.name()

"; 

%feature("docstring") simuPOP::pyFunc::arg "

Usage:

    x.arg(arg)

"; 

%feature("docstring") simuPOP::pyFunc::isValid "

Usage:

    x.isValid()

"; 

%feature("docstring") simuPOP::pyFunc::func "

Usage:

    x.func()

"; 

%ignore simuPOP::pyGenerator;

%feature("docstring") simuPOP::pyGenerator::pyGenerator "

Usage:

    pyGenerator(gen=None)

"; 

%feature("docstring") simuPOP::pyGenerator::isValid "

Usage:

    x.isValid()

"; 

%feature("docstring") simuPOP::pyGenerator::set "

Usage:

    x.set(gen)

"; 

%feature("docstring") simuPOP::pyGenerator::~pyGenerator "

Usage:

    x.~pyGenerator()

"; 

%feature("docstring") simuPOP::pyGenerator::next "

Usage:

    x.next()

"; 

%feature("docstring") simuPOP::pyIndIterator "

Details:

    this class implements a Python itertor class that can be used to
    iterate through individuals in a (sub)population. If allInds are
    true, visiblility of individuals will not be checked. Otherwise, a
    functor will be used to check if indiviudals belong to a specified
    virtual subpopulation.  An instance of this class is returned by
    population::Individuals() and Population::Individuals(subPop)

"; 

%feature("docstring") simuPOP::pyIndIterator::pyIndIterator "

Usage:

    pyIndIterator(begin, end, allInds, func)

"; 

%feature("docstring") simuPOP::pyIndIterator::~pyIndIterator "

Usage:

    x.~pyIndIterator()

"; 

%feature("docstring") simuPOP::pyIndIterator::iter "

Usage:

    x.__iter__()

"; 

%feature("docstring") simuPOP::pyIndIterator::next "

Usage:

    x.next()

"; 

%feature("docstring") simuPOP::pyIndIterator::next "

Usage:

    x.__next__()

"; 

%feature("docstring") simuPOP::PyMlPenetrance "

Details:

    This penetrance operator is a multi-locus Python penetrance
    operator that assigns penetrance values by combining locus and
    genotype specific penetrance values. It differs from a
    PyPenetrance in that the python function is responsible for
    penetrance values values for each gentoype type at each locus,
    which can potentially be random, and locus or gentoype-specific.

"; 

%feature("docstring") simuPOP::PyMlPenetrance::PyMlPenetrance "

Usage:

    PyMlPenetrance(func, mode=MULTIPLICATIVE, loci=ALL_AVAIL,
      ancGens=UNSPECIFIED, output=\"\", begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a penetrance operator that assigns individual affection
    status according to penetrance values combined from locus-specific
    penetrance values that are determined by a Python call-back
    function. The callback function accepts parameter loc, alleles
    (both optional) and returns location- or genotype-specific
    penetrance values that can be constant or random. The penetrance
    values for each genotype will be cached so the same penetrance
    values will be assigned to genotypes with previously assigned
    values. Note that a function that does not examine the genotype
    naturally assumes a dominant model where genotypes with one or two
    mutants have the same penetrance value. Because genotypes at a
    locus are passed separately and in no particular order, this
    function is also responsible for assigning consistent fitness
    values for genotypes at the same locus (a class is usually used).
    This operator currently ignores chromosome types so unused alleles
    will be passed for loci on sex or mitochondrial chromosomes. This
    operator also ignores the phase of genotype so genotypes (a,b) and
    (b,a) are assumed to have the same fitness effect.   Individual
    penetrance will be combined in ADDITIVE, MULTIPLICATIVE, or
    HETEROGENEITY mode from penetrance values of loci with at least
    one non-zero allele (See MlPenetrance for details).

"; 

%feature("docstring") simuPOP::PyMlPenetrance::~PyMlPenetrance "

Usage:

    x.~PyMlPenetrance()

"; 

%feature("docstring") simuPOP::PyMlPenetrance::clone "Obsolete or undocumented function."

%ignore simuPOP::PyMlPenetrance::penet(Population *pop, RawIndIterator ind) const;

%feature("docstring") simuPOP::PyMlPenetrance::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyMlSelector "

Details:

    This selector is a multi-locus Python selector that assigns
    fitness to individuals by combining locus and genotype specific
    fitness values. It differs from a PySelector in that the python
    function is responsible for assigning fitness values for each
    gentoype type at each locus, which can potentially be random, and
    locus or gentoype-specific.

"; 

%feature("docstring") simuPOP::PyMlSelector::PyMlSelector "

Usage:

    PyMlSelector(func, mode=EXPONENTIAL, loci=ALL_AVAIL, output=\"\",
      begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=ALL_AVAIL)

Details:

    Create a selector that assigns individual fitness values by
    combining locus-specific fitness values that are determined by a
    Python call-back function. The callback function accepts parameter
    loc, alleles (both optional) and returns location- or genotype-
    specific fitness values that can be constant or random. The
    fitness values for each genotype will be cached so the same
    fitness values will be assigned to genotypes with previously
    assigned values. Note that a function that does not examine the
    genotype naturally assumes a dominant model where genotypes with
    one or two mutants have the same fitness effect. Because genotypes
    at a locus are passed separately and in no particular order, this
    function is also responsible for assigning consistent fitness
    values for genotypes at the same locus (a class is usually used).
    This operator currently ignores chromosome types so unused alleles
    will be passed for loci on sex or mitochondrial chromosomes. It
    also ignores phase of genotype so it will use the same fitness
    value for genotype (a,b) and (b,a).   Individual fitness will be
    combined in ADDITIVE, MULTIPLICATIVE, HETEROGENEITY, or
    EXPONENTIAL mode from fitness values of loci with at least one
    non-zero allele (See MlSelector for details). If an output is
    given, location, genotype, fitness and generation at which the new
    genotype is assgined the value will be written to the output, in
    the format of 'loc a1 a2 fitness gen' for loci on autosomes of
    diploid populations.

"; 

%feature("docstring") simuPOP::PyMlSelector::~PyMlSelector "

Usage:

    x.~PyMlSelector()

"; 

%feature("docstring") simuPOP::PyMlSelector::clone "Obsolete or undocumented function."

%ignore simuPOP::PyMlSelector::indFitness(Population &pop, RawIndIterator ind) const;

%feature("docstring") simuPOP::PyMlSelector::describe "Obsolete or undocumented function."

%ignore simuPOP::PyMlSelector::apply(Population &pop) const;

%feature("docstring") simuPOP::pyMutantIterator "

Details:

    this class implements a Python itertor class that can be used to
    iterate through individuals in a (sub)population. If allInds are
    true, visiblility of individuals will not be checked. Otherwise, a
    functor will be used to check if indiviudals belong to a specified
    virtual subpopulation.  An instance of this class is returned by
    population::Individuals() and Population::Individuals(subPop)

"; 

%feature("docstring") simuPOP::pyMutantIterator::pyMutantIterator "

Usage:

    pyMutantIterator(base, begin, end, step)

"; 

%feature("docstring") simuPOP::pyMutantIterator::~pyMutantIterator "

Usage:

    x.~pyMutantIterator()

"; 

%feature("docstring") simuPOP::pyMutantIterator::iter "

Usage:

    x.__iter__()

"; 

%feature("docstring") simuPOP::pyMutantIterator::next "

Usage:

    x.next()

"; 

%feature("docstring") simuPOP::pyMutantIterator::next "

Usage:

    x.__next__()

"; 

%feature("docstring") simuPOP::PyMutator "

Details:

    This hybrid mutator accepts a Python function that determines how
    to mutate an allele when an mutation event happens.

"; 

%feature("docstring") simuPOP::PyMutator::PyMutator "

Usage:

    PyMutator(rates=[], loci=ALL_AVAIL, func=None, context=0,
      mapIn=[], mapOut=[], output=\">\", begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=\"ind_id\",
      lineageMode=FROM_INFO)

Details:

    Create a hybrid mutator that uses a user-provided function to
    mutate an allele when a mutation event happens. This function
    (parameter func) accepts the allele to be mutated as parameter
    allele and optional array of alleles as parameter context, which
    are context alleles the left and right of the mutated allele.
    Invalid context alleles (e.g. left allele to the first locus of a
    chromosome) will be marked by -1. The return value of this
    function will be used to mutate the passed allele. The passed,
    returned and context alleles might be altered if parameter mapIn
    and mapOut are used. This mutator by default applies to all loci
    unless parameter loci is specified. A single mutation rate will be
    used for all loci if a single value of parameter rates is given.
    Otherwise, a list of mutation rates can be specified for each
    locus in parameter loci. Please refer to classes mutator and
    BaseOperator for descriptions of other parameters.

"; 

%feature("docstring") simuPOP::PyMutator::clone "Obsolete or undocumented function."

%ignore simuPOP::PyMutator::mutate(Allele allele, size_t locus) const;

%feature("docstring") simuPOP::PyMutator::describe "Obsolete or undocumented function."

%ignore simuPOP::pyObject;

%feature("docstring") simuPOP::pyObject::pyObject "

Usage:

    pyObject(obj, defToNone=False)

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

%feature("docstring") simuPOP::PyOperator "

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

%feature("docstring") simuPOP::PyOperator::PyOperator "

Usage:

    PyOperator(func, param=None, begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a pure-Python operator that calls a user-defined function
    when it is applied. If this operator is applied before or after
    mating, your function should have form func(pop) or func(pop,
    param) where pop is the population to which the operator is
    applied, param is the value specified in parameter param. param
    will be ignored if your function only accepts one parameter.
    Althernatively, the function should have form func(ind) with
    optional parameters pop and param. In this case, the function will
    be called for all individuals, or individuals in subpopulations
    subPops. Individuals for which the function returns False will be
    removed from the population. This operator can therefore perform
    similar functions as operator DiscardIf.  If this operator is
    applied during mating, your function should accept parameters pop,
    off (or ind), dad, mom and param where pop is the parental
    population, and off or ind, dad, and mom are offspring and their
    parents for each mating event, and param is an optional parameter.
    If subPops are provided, only offspring in specified (virtual)
    subpopulations are acceptable.  This operator does not support
    parameters output, and infoFields. If certain output is needed, it
    should be handled in the user defined function func. Because the
    status of files used by other operators through parameter output
    is undetermined during evolution, they should not be open or
    closed in this Python operator.

"; 

%feature("docstring") simuPOP::PyOperator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyOperator::apply "Obsolete or undocumented function."

%ignore simuPOP::PyOperator::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::PyOperator::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyOutput "

Details:

    This operator outputs a given string when it is applied to a
    population.

"; 

%feature("docstring") simuPOP::PyOutput::PyOutput "

Usage:

    PyOutput(msg=\"\", output=\">\", begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Creates a PyOutput operator that outputs a string msg to output
    (default to standard terminal output) when it is applied to a
    population. Please refer to class BaseOperator for a detailed
    description of common operator parameters such as stage, begin and
    output.

"; 

%feature("docstring") simuPOP::PyOutput::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyOutput::~PyOutput "

Description:

    destructor

Usage:

    x.~PyOutput()

"; 

%feature("docstring") simuPOP::PyOutput::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyOutput::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyParentsChooser "

Details:

    This parent chooser accepts a Python generator function that
    repeatedly yields one or two parents, which can be references to
    individual objects or indexes relative to each subpopulation. The
    parent chooser calls the generator function with parental
    population and a subpopulation index for each subpopulation and
    retrieves parents repeatedly using the iterator interface of the
    generator function.  This parent chooser does not support virtual
    subpopulation directly. However, because virtual subpopulations
    are defined in the passed parental population, it is easy to
    return parents from a particular virtual subpopulation using
    virtual subpopulation related functions.

"; 

%feature("docstring") simuPOP::PyParentsChooser::PyParentsChooser "

Usage:

    PyParentsChooser(generator)

Details:

    Create a Python parent chooser using a Python generator function
    parentsGenerator. This function should accept one or both of
    parameters pop (the parental population) and subPop (index of
    subpopulation) and return the reference or index (relative to
    subpopulation) of a parent or a pair of parents repeatedly using
    the iterator interface of the generator function.

"; 

%ignore simuPOP::PyParentsChooser::PyParentsChooser(const PyParentsChooser &rhs);

%feature("docstring") simuPOP::PyParentsChooser::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyParentsChooser::describe "Obsolete or undocumented function."

%ignore simuPOP::PyParentsChooser::parallelizable() const;

%ignore simuPOP::PyParentsChooser::initialize(Population &pop, size_t sp);

%ignore simuPOP::PyParentsChooser::finalize(Population &pop, size_t sp);

%feature("docstring") simuPOP::PyParentsChooser::~PyParentsChooser "

Description:

    destructor

Usage:

    x.~PyParentsChooser()

"; 

%ignore simuPOP::PyParentsChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::PyPenetrance "

Details:

    This penetrance operator assigns penetrance values by calling a
    user provided function. It accepts a list of loci (parameter
    loci), and a Python function func which should be defined with one
    or more of parameters geno, gen, ind, pop, or names of information
    fields. When this operator is applied to a population, it passes
    genotypes at specified loci, generation number, a reference to an
    individual, a reference to the current population (usually used to
    retrieve population variables) and values at specified information
    fields to respective parameters of this function. The returned
    penetrance values will be used to determine the affection status
    of each individual.

"; 

%feature("docstring") simuPOP::PyPenetrance::PyPenetrance "

Usage:

    PyPenetrance(func, loci=[], ancGens=UNSPECIFIED, begin=0,
      end=-1, step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL,
      infoFields=[])

Details:

    Create a Python hybrid penetrance operator that passes genotype at
    specified loci, values at specified information fields (if
    requested), and a generation number to a user-defined function
    func. Parameter loci can be a list of loci indexes, names, or
    ALL_AVAIL. The return value will be treated as Individual
    penetrance.

"; 

%feature("docstring") simuPOP::PyPenetrance::clone "Obsolete or undocumented function."

%ignore simuPOP::PyPenetrance::penet(Population *pop, RawIndIterator ind) const;

%feature("docstring") simuPOP::PyPenetrance::describe "Obsolete or undocumented function."

%ignore simuPOP::PyPenetrance::parallelizable() const;

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

%feature("docstring") simuPOP::pyPopIterator::iter "

Usage:

    x.__iter__()

"; 

%feature("docstring") simuPOP::pyPopIterator::next "

Usage:

    x.next()

"; 

%feature("docstring") simuPOP::pyPopIterator::next "

Usage:

    x.__next__()

"; 

%feature("docstring") simuPOP::PyQuanTrait "

Details:

    This quantitative trait operator assigns a trait field by calling
    a user provided function. It accepts a list of loci (parameter
    loci), and a Python function func which should be defined with one
    or more of parameters geno, gen, ind, or names of information
    fields. When this operator is applied to a population, it passes
    genotypes at specified loci, generation number, a reference to an
    individual, and values at specified information fields to
    respective parameters of this function. The return values will be
    assigned to specified trait fields.

"; 

%feature("docstring") simuPOP::PyQuanTrait::PyQuanTrait "

Usage:

    PyQuanTrait(func, loci=[], ancGens=UNSPECIFIED, begin=0, end=-1,
      step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a Python hybrid quantitative trait operator that passes
    genotype at specified loci, optional values at specified
    information fields (if requested), and an optional generation
    number to a user-defined function func. Parameter loci can be a
    list of loci indexes, names, or ALL_AVAIL. The return value will
    be assigned to specified trait fields (infoField). If only one
    trait field is specified, a number or a sequence of one element is
    acceptable. Otherwise, a sequence of values will be accepted and
    be assigned to each trait field.

"; 

%feature("docstring") simuPOP::PyQuanTrait::clone "Obsolete or undocumented function."

%ignore simuPOP::PyQuanTrait::qtrait(Individual *ind, size_t gen, vectorf &traits) const;

%feature("docstring") simuPOP::PyQuanTrait::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::PySelector "

Details:

    This selector assigns fitness values by calling a user provided
    function. It accepts a list of loci (parameter loci) and a Python
    function func which should be defined with one or more of
    parameters geno, gen, ind, pop or names of information fields.
    Parameter loci can be a list of loci indexes, names or ALL_AVAIL.
    When this operator is applied to a population, it passes genotypes
    at specified loci, generation number, a reference to an
    individual, a reference to the current population (usually used to
    retrieve population variable), and values at specified information
    fields to respective parameters of this function. The returned
    value will be used to determine the fitness of each individual.

"; 

%feature("docstring") simuPOP::PySelector::PySelector "

Usage:

    PySelector(func, loci=[], begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=ALL_AVAIL)

Details:

    Create a Python hybrid selector that passes genotype at specified
    loci, values at specified information fields (if requested) and a
    generation number to a user-defined function func. The return
    value will be treated as individual fitness.

"; 

%feature("docstring") simuPOP::PySelector::clone "Obsolete or undocumented function."

%ignore simuPOP::PySelector::indFitness(Population &pop, RawIndIterator ind) const;

%feature("docstring") simuPOP::PySelector::describe "Obsolete or undocumented function."

%ignore simuPOP::PySelector::parallelizable() const;

%feature("docstring") simuPOP::PyTagger "

Details:

    A Python tagger takes some information fields from both parents,
    pass them to a user provided Python function and set the offspring
    individual fields with the return values.

"; 

%feature("docstring") simuPOP::PyTagger::PyTagger "

Usage:

    PyTagger(func=None, begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, output=\"\", infoFields=[])

Details:

    Create a hybrid tagger that provides an user provided function
    func with values of specified information fields (determined by
    parameter names of this function) of parents and assign
    corresponding information fields of offspring with its return
    value. If more than one parent are available, maternal values are
    passed after paternal values. For example, if a function func(A,
    B) is passed, this operator will send two tuples with parental
    values of information fields 'A' and 'B' to this function and
    assign its return values to fields 'A' and 'B' of each offspring.
    The return value of this function should be a list, although a
    single value will be accepted if only one information field is
    specified. This operator ignores parameters stage, output and
    subPops.

"; 

%feature("docstring") simuPOP::PyTagger::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::PyTagger::describe "Obsolete or undocumented function."

%ignore simuPOP::PyTagger::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%ignore simuPOP::PyTagger::parallelizable() const;

%feature("docstring") simuPOP::RandomParentChooser "

Details:

    This parent chooser chooses a parent randomly from a (virtual)
    parental subpopulation. Parents are chosen with or without
    replacement. If parents are chosen with replacement, a parent can
    be selected multiple times. If individual fitness values are
    assigned to individuals ( stored in an information field
    selectionField (default to \"fitness\"), individuals will be chosen
    at a probability proportional to his or her fitness value. If
    parents are chosen without replacement, a parent can be chosen
    only once. An RuntimeError will be raised if all parents are
    exhausted. Natural selection is disabled in the without-
    replacement case.

"; 

%feature("docstring") simuPOP::RandomParentChooser::RandomParentChooser "

Usage:

    RandomParentChooser(replacement=True, selectionField=\"fitness\",
      sexChoice=ANY_SEX)

Details:

    Create a random parent chooser that choose parents with or without
    replacement (parameter replacement, default to True). If selection
    is enabled and information field selectionField exists in the
    passed population, the probability that a parent is chosen is
    proportional to his/her fitness value stored in selectionField.
    This parent chooser by default chooses parent from all individuals
    (ANY_SEX), but it can be made to select only male (MALE_ONLY) or
    female (FEMALE_ONLY) individuals by setting parameter sexChoice.

"; 

%feature("docstring") simuPOP::RandomParentChooser::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::RandomParentChooser::describe "Obsolete or undocumented function."

%ignore simuPOP::RandomParentChooser::parallelizable() const;

%ignore simuPOP::RandomParentChooser::initialize(Population &pop, size_t sp);

%ignore simuPOP::RandomParentChooser::chooseParents(RawIndIterator basePtr);

%feature("docstring") simuPOP::RandomParentsChooser "

Details:

    This parent chooser chooses two parents, a male and a female,
    randomly from a (virtual) parental subpopulation. Parents are
    chosen with or without replacement from their respective sex
    group. If parents are chosen with replacement, a parent can be
    selected multiple times. If individual fitness values are assigned
    (stored in information field selectionField, default to \"fitness\",
    the probability that an individual is chosen is proportional to
    his/her fitness value among all individuals with the same sex. If
    parents are chosen without replacement, a parent can be chosen
    only once. An RuntimeError will be raised if all males or females
    are exhausted. Natural selection is disabled in the without-
    replacement case.

"; 

%feature("docstring") simuPOP::RandomParentsChooser::RandomParentsChooser "

Usage:

    RandomParentsChooser(replacement=True, selectionField=\"fitness\")

Details:

    Create a random parents chooser that choose two parents with or
    without replacement (parameter replacement, default to True). If
    selection is enabled and information field selectionField exists
    in the passed population, the probability that a parent is chosen
    is proportional to his/her fitness value stored in selectionField.

"; 

%feature("docstring") simuPOP::RandomParentsChooser::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::RandomParentsChooser::describe "Obsolete or undocumented function."

%ignore simuPOP::RandomParentsChooser::parallelizable() const;

%ignore simuPOP::RandomParentsChooser::initialize(Population &pop, size_t sp);

%ignore simuPOP::RandomParentsChooser::chooseParents(RawIndIterator basePtr);

%ignore simuPOP::RandomSexModel;

%feature("docstring") simuPOP::RandomSexModel::RandomSexModel "

Usage:

    RandomSexModel()

"; 

%feature("docstring") simuPOP::RandomSexModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::RandomSexModel::getSex "

Usage:

    x.getSex(UINT)

"; 

%feature("docstring") simuPOP::RandomSexModel::parallelizable "

Usage:

    x.parallelizable()

"; 

%feature("docstring") simuPOP::RangeSplitter "

Details:

    This class defines a splitter that groups individuals in certain
    ranges into VSPs.

"; 

%feature("docstring") simuPOP::RangeSplitter::RangeSplitter "

Usage:

    RangeSplitter(ranges, names=[])

Details:

    Create a splitter according to a number of individual ranges
    defined in ranges. For example, RangeSplitter(ranges=[[0, 20],
    [40, 50]]) defines two VSPs. The first VSP consists of individuals
    0, 1, ..., 19, and the sceond VSP consists of individuals 40, 41,
    ..., 49. Note that a nested list has to be used even if only one
    range is defined. A default set of names are given to each VSP
    unless a new set of names is given by parameter names.

"; 

%feature("docstring") simuPOP::RangeSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::RangeSplitter::size(const Population &pop, size_t subPop, size_t virtualSubPop) const;

%feature("docstring") simuPOP::RangeSplitter::numVirtualSubPop "

Usage:

    x.numVirtualSubPop()

Details:

    Return the number of VSPs, which is the number of ranges defined
    in parameter ranges.

"; 

%ignore simuPOP::RangeSplitter::contains(const Population &pop, size_t ind, vspID vsp) const;

%ignore simuPOP::RangeSplitter::activate(const Population &pop, size_t subPop, size_t virtualSubPop);

%feature("docstring") simuPOP::RangeSplitter::name "

Usage:

    x.name(vsp)

Details:

    Return the name of VSP vsp, which is \"Range [a, b)\" where [a, b)
    is range ranges[vsp]. A user specified name will be returned if
    specified.

"; 

%feature("docstring") simuPOP::Recombinator "

Details:

    A genotype transmitter (during-mating operator) that transmits
    parental chromosomes to offspring, subject to recombination and
    gene conversion. This can be used to replace
    MendelianGenoTransmitter and SelfingGenoTransmitter. It does not
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

%feature("docstring") simuPOP::Recombinator::Recombinator "

Usage:

    Recombinator(rates=[], intensity=-1, loci=ALL_AVAIL,
      convMode=NO_CONVERSION, output=\"\", begin=0, end=-1, step=1,
      at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a Recombinator (a mendelian genotype transmitter with
    recombination and gene conversion) that passes genotypes from
    parents (or a parent in case of self-fertilization) to offspring.
    Recombination happens by default between all adjacent markers but
    can be limited to a given set of loci, which can be a list of loci
    indexes, names or ALL_AVAIL. Each locus in this list specifies a
    recombination point between the locus and the locus immediately
    after it. Loci that are the last locus on each chromosome are
    ignored.  If a single recombination rate (parameter rates) is
    specified, it will used for all loci (all loci or loci specified
    by parameter loci), regardless of physical distances between
    adjacent loci.  If a list of recombination rates are specified in
    rates, a parameter loci with the same length should also be
    specified. Different recombination rates can then be used after
    these loci (between specified loci and their immediate neighbor to
    the right).  A recombination intensity (intensity) can be used to
    specify recombination rates that are proportional to physical
    distances between adjacent markers. If the physical distance
    between two markers is d, the recombination rate between them will
    be intensity * d. No unit is assume for loci position and
    recombination intensity.  Gene conversion is controlled using
    parameter convMode, which can be
    *   NoConversion: no gene conversion (default).
    *   (NUM_MARKERS, prob, n): With probability prob, convert a fixed
    number (n) of markers if a recombination event happens.
    *   (GEOMETRIC_DISTRIBUTION, prob, p): With probability prob,
    convert a random number of markers if a recombination event
    happens. The number of markes converted follows a geometric
    distribution with probability p.
    *   (TRACT_LENGTH, prob, n): With probability prob, convert a
    region of fixed tract length (n) if a recombination event happens.
    The actual number of markers converted depends on loci positions
    of surrounding loci. The starting position of this tract is the
    middle of two adjacent markers. For example, if four loci are
    located at 0, 1, 2, 3 respectively, a conversion event happens
    between 0 and 1, with a tract length 2 will start at 0.5 and end
    at 2.5, covering the second and third loci.
    *   (EXPONENTIAL_DISTRIBUTION, prob, p): With probability prob,
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
    Recombinator usually does not send any output. However, if an
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
    IdTagger operator) before this Recombinator is applied.  In
    addition to genotypes, this operator also copies alleleic lineage
    if it is executed in a module with lineage allele type.

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

%feature("docstring") simuPOP::Recombinator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::Recombinator::~Recombinator "

Usage:

    x.~Recombinator()

"; 

%feature("docstring") simuPOP::Recombinator::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::Recombinator::initialize "Obsolete or undocumented function."

%feature("docstring") simuPOP::Recombinator::transmitGenotype "

Usage:

    x.transmitGenotype(parent, offspring, ploidy)

Details:

    This function transmits genotypes from a parent to the ploidy-th
    homologous set of chromosomes of an offspring. It can be used, for
    example, by a customized genotype transmitter to use sex-specific
    recombination rates to transmit parental genotypes to offspring.

"; 

%ignore simuPOP::Recombinator::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad, Individual *mom) const;

%ignore simuPOP::Recombinator::parallelizable() const;

%feature("docstring") simuPOP::ResizeSubPops "

Details:

    This operator resizes subpopulations to specified sizes.
    individuals are added or removed depending on the new
    subpopulation sizes.

"; 

%feature("docstring") simuPOP::ResizeSubPops::ResizeSubPops "

Usage:

    ResizeSubPops(subPops=ALL_AVAIL, sizes=[], proportions=[],
      propagate=True, begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
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
    Please refer to operator BaseOperator for a detailed explanation
    for all parameters.

"; 

%feature("docstring") simuPOP::ResizeSubPops::~ResizeSubPops "

Description:

    destructor

Usage:

    x.~ResizeSubPops()

"; 

%feature("docstring") simuPOP::ResizeSubPops::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::ResizeSubPops::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::ResizeSubPops::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::RevertFixedSites "

Details:

    This operator checks all or specified loci of a population and
    revert all mutants at loci to wildtype alleles if they are fixed
    in the population. If a list of (virtual) subpopulations are
    specified, alleles are reverted if they are fixed in each
    subpopulation, regardless if the alleles are fixed in other
    subpopulations.

"; 

%feature("docstring") simuPOP::RevertFixedSites::RevertFixedSites "

Usage:

    RevertFixedSites(loci=ALL_AVAIL, output=\"\", begin=0, end=-1,
      step=1, at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL,
      infoFields=\"ind_id\")

Details:

    Create an operator to set all alleles to zero at specified
    (parameter loci) or all loci if they are fixed (having no zero-
    allele) at these loci. If parameter subPops are specified, only
    individuals in these subpopulations are considered.

"; 

%feature("docstring") simuPOP::RevertFixedSites::~RevertFixedSites "

Description:

    destructor

Usage:

    x.~RevertFixedSites()

"; 

%feature("docstring") simuPOP::RevertFixedSites::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::RevertFixedSites::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::RevertFixedSites::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::RNG "

Details:

    This random number generator class wraps around a number of random
    number generators from GNU Scientific Library. You can obtain and
    change the RNG used by the current simuPOP module through the
    getRNG() function, or create a separate random number generator
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
    Names of supported random number generators are available from
    moduleInfo()['availableRNGs'].

"; 

%ignore simuPOP::RNG::RNG(const RNG &);

%feature("docstring") simuPOP::RNG::~RNG "

Usage:

    x.~RNG()

"; 

%feature("docstring") simuPOP::RNG::set "

Usage:

    x.set(name=None, seed=0)

Details:

    Replace the existing random number generator using RNGname with
    seed seed. If seed is 0, a random seed will be used. If name is
    empty, use the existing RNG but reset the seed.

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

    Generate a random number following a gamma distribution with a
    shape parameters a and scale parameter b.

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

%feature("docstring") simuPOP::RNG::randTruncatedPoisson "

Usage:

    x.randTruncatedPoisson(mu)

Details:

    Generate a positive random number following a zero-truncated
    Poisson distribution with parameter mu.

"; 

%feature("docstring") simuPOP::RNG::randTruncatedBinomial "

Usage:

    x.randTruncatedBinomial(n, p)

Details:

    Generate a positive random number following a zero-truncated
    binomial distribution with parameters n and p.

"; 

%feature("docstring") simuPOP::RNG::randMultinomial "

Usage:

    x.randMultinomial(N, p)

Details:

    Generate a random number following a multinomial distribution with
    parameters N and p (a list of probabilities).

"; 

%ignore simuPOP::RNG::randomShuffle(T begin, T end) const;

%ignore simuPOP::RNGfunc;

%feature("docstring") simuPOP::RNGfunc::RNGfunc "

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

%feature("docstring") simuPOP::SavePopulation "

Details:

    An operator that save populations to specified files.

"; 

%feature("docstring") simuPOP::SavePopulation::SavePopulation "

Usage:

    SavePopulation(output=\"\", begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create an operator that saves a population to output when it is
    applied to the population. This operator supports all output
    specifications ('', 'filename', 'filename' prefixed by one or more
    '>' characters, and '!expr') but output from different operators
    will always replace existing files (effectively ignore '>'
    specification). Parameter subPops is ignored. Please refer to
    class BaseOperator for a detailed description about common
    operator parameters such as stage and begin.

"; 

%feature("docstring") simuPOP::SavePopulation::~SavePopulation "

Description:

    destructor.

Usage:

    x.~SavePopulation()

"; 

%feature("docstring") simuPOP::SavePopulation::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::SavePopulation::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::SavePopulation::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::SelfingGenoTransmitter "

Details:

    A genotype transmitter (during-mating operator) that transmits
    parental genotype of a parent through self-fertilization. That is
    to say, the offspring genotype is formed according to Mendel's
    laws, only that a parent serves as both maternal and paternal
    parents.

"; 

%feature("docstring") simuPOP::SelfingGenoTransmitter::SelfingGenoTransmitter "

Usage:

    SelfingGenoTransmitter(output=\"\", begin=0, end=-1, step=1,
      at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a self-fertilization genotype transmitter that transmits
    genotypes of a parent to an offspring through self-fertilization.
    Cutsomized chromosomes are not handled. Parameters subPops and
    infoFields are ignored. This operator also copies allelic lineage
    when it is executed in a module with lineage allele type.

"; 

%feature("docstring") simuPOP::SelfingGenoTransmitter::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::SelfingGenoTransmitter::describe "Obsolete or undocumented function."

%ignore simuPOP::SelfingGenoTransmitter::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%ignore simuPOP::SeqSexModel;

%feature("docstring") simuPOP::SeqSexModel::SeqSexModel "

Usage:

    SeqSexModel(sex)

"; 

%feature("docstring") simuPOP::SeqSexModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::SeqSexModel::getSex "

Usage:

    x.getSex(count)

"; 

%feature("docstring") simuPOP::SeqSexModel::parallelizable "

Usage:

    x.parallelizable()

"; 

%feature("docstring") simuPOP::SequentialParentChooser "

Details:

    This parent chooser chooses a parent from a parental (virtual)
    subpopulation sequentially. Natural selection is not considered.
    If the last parent is reached, this parent chooser will restart
    from the beginning of the (virtual) subpopulation.

"; 

%feature("docstring") simuPOP::SequentialParentChooser::SequentialParentChooser "

Usage:

    SequentialParentChooser(sexChoice=ANY_SEX)

Details:

    Create a parent chooser that chooses a parent from a parental
    (virtual) subpopulation sequentially. Parameter choice can be
    ANY_SEX (default), MALE_ONLY and FEMALE_ONLY. In the latter two
    cases, only male or female individuals are selected. A
    RuntimeError will be raised if there is no male or female
    individual from the population.

"; 

%feature("docstring") simuPOP::SequentialParentChooser::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::SequentialParentChooser::describe "Obsolete or undocumented function."

%ignore simuPOP::SequentialParentChooser::initialize(Population &pop, size_t sp);

%ignore simuPOP::SequentialParentChooser::chooseParents(RawIndIterator basePtr);

%ignore simuPOP::SexModel;

%feature("docstring") simuPOP::SexModel::SexModel "

Usage:

    SexModel()

"; 

%feature("docstring") simuPOP::SexModel::~SexModel "

Usage:

    x.~SexModel()

"; 

%feature("docstring") simuPOP::SexModel::getSex "

Usage:

    x.getSex(count)

"; 

%feature("docstring") simuPOP::SexModel::reset "

Usage:

    x.reset()

"; 

%feature("docstring") simuPOP::SexModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::SexModel::parallelizable "

Usage:

    x.parallelizable()

"; 

%feature("docstring") simuPOP::SexSplitter "

Details:

    This splitter defines two VSPs by individual sex. The first VSP
    consists of all male individuals and the second VSP consists of
    all females in a subpopulation.

"; 

%feature("docstring") simuPOP::SexSplitter::SexSplitter "

Usage:

    SexSplitter(names=[])

Details:

    Create a sex splitter that defines male and female VSPs. These
    VSPs are named Male and Female unless a new set of names are
    specified by parameter names.

"; 

%feature("docstring") simuPOP::SexSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::SexSplitter::size(const Population &pop, size_t subPop, size_t virtualSubPop) const;

%feature("docstring") simuPOP::SexSplitter::numVirtualSubPop "

Description:

    Return 2.

Usage:

    x.numVirtualSubPop()

"; 

%ignore simuPOP::SexSplitter::contains(const Population &pop, size_t ind, vspID vsp) const;

%ignore simuPOP::SexSplitter::activate(const Population &pop, size_t subPop, size_t virtualSubPop);

%feature("docstring") simuPOP::SexSplitter::name "

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

    destructor I can not clear dict here since a resize of g_vars will
    copy this object and hence call this destructore.

Usage:

    x.~SharedVariables()

"; 

%ignore simuPOP::SharedVariables::clear();

%ignore simuPOP::SharedVariables::setDict(PyObject *dict);

%ignore simuPOP::SharedVariables::setVar(const string &name, const PyObject *val);

%ignore simuPOP::SharedVariables::getVar(const string &name, bool nameError=true) const;

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

%ignore simuPOP::SharedVariables::setVar(const string &name, const bool val);

%ignore simuPOP::SharedVariables::getVarAsBool(const string &name, bool nameError=true) const;

%ignore simuPOP::SharedVariables::getVarAsInt(const string &name, bool nameError=true) const;

%ignore simuPOP::SharedVariables::getVarAsDouble(const string &name, bool nameError=true) const;

%ignore simuPOP::SharedVariables::getVarAsString(const string &name, bool nameError=true) const;

%ignore simuPOP::SharedVariables::getVarAsIntDict(const string &name, uintDict &res, bool nameError=true) const;

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

    simpleStmt(stmt, indVar)

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

%feature("docstring") simuPOP::Simulator "

Details:

    A simuPOP simulator is responsible for evolving one or more
    populations forward in time, subject to various operators.
    Populations in a simulator are created from one or more replicates
    of specified populations. A number of functions are provided to
    access and manipulate populations, and most importantly, to evolve
    them.

"; 

%feature("docstring") simuPOP::Simulator::Simulator "

Usage:

    Simulator(pops, rep=1, stealPops=True)

Details:

    Create a simulator with rep (default to 1) replicates of
    populations pops, which is a list of populations although a single
    population object is also acceptable. Contents of passed
    populations are by default moved to the simulator to avoid
    duplication of potentially large population objects, leaving empty
    populations behind. This behavior can be changed by setting
    stealPops to False, in which case populations are copied to the
    simulator.

"; 

%feature("docstring") simuPOP::Simulator::~Simulator "

Usage:

    x.~Simulator()

"; 

%ignore simuPOP::Simulator::Simulator(const Simulator &rhs);

%feature("docstring") simuPOP::Simulator::clone "

Usage:

    x.clone()

Details:

    Clone a simulator, along with all its populations. Note that
    Python assign statement simu1 = simu only creates a symbolic link
    to an existing simulator.

"; 

%feature("docstring") simuPOP::Simulator::numRep "

Usage:

    x.numRep()

Details:

    Return the number of replicates.

"; 

%feature("docstring") simuPOP::Simulator::population "

Usage:

    x.population(rep)

Details:

    Return a reference to the rep-th population of a simulator. The
    reference will become invalid once the simulator starts evolving
    or becomes invalid (removed). If an independent copy of the
    population is needed, you can use population.clone() to create a
    cloned copy or simulator.extract() to remove the population from
    the simulator.

"; 

%feature("docstring") simuPOP::Simulator::add "

Usage:

    x.add(pop, stealPop=True)

Details:

    Add a population pop to the end of an existing simulator. This
    function by default moves pop to the simulator, leaving an empty
    population for passed population object. If steal is set to False,
    the population will be copied to the simulator, and thus
    unchanged.

"; 

%feature("docstring") simuPOP::Simulator::extract "

Usage:

    x.extract(rep)

Details:

    Extract the rep-th population from a simulator. This will reduce
    the number of populations in this simulator by one.

"; 

%feature("docstring") simuPOP::Simulator::populations "

Usage:

    x.populations()

Details:

    Return a Python iterator that can be used to iterate through all
    populations in a simulator.

"; 

%ignore simuPOP::Simulator::describe(bool format=true) const;

%feature("docstring") simuPOP::Simulator::evolve "

Usage:

    x.evolve(initOps=[], preOps=[], matingScheme=MatingScheme,
      postOps=[], finalOps=[], gen=-1, dryrun=False)

Details:

    Evolve all populations gen generations, subject to several lists
    of operators which are applied at different stages of an
    evolutionary process. Operators initOps are applied to all
    populations (subject to applicability restrictions of the
    operators, imposed by the rep parameter of these operators) before
    evolution. They are used to initialize populations before
    evolution. Operators finalOps are applied to all populations after
    the evolution.  Operators preOps, and postOps are applied during
    the life cycle of each generation. These operators can be applied
    at all or some of the generations, to all or some of the evolving
    populations, depending the begin, end, step, at and reps
    parameters of these operators. These operators are applied in the
    order at which they are specified. populations in a simulator are
    evolved one by one. At each generation, operators preOps are
    applied to the parental generations. A mating scheme is then used
    to populate an offspring generation. For each offspring, his or
    her sex is determined before during-mating operators of the mating
    scheme are used to transmit parental genotypes. After an offspring
    generation is successfully generated and becomes the current
    generation, operators postOps are applied to the offspring
    generation. If any of the preOps and postOps fails (return False),
    the evolution of a population will be stopped. The generation
    number of a population, which is the variable \"gen\" in each
    populations local namespace, is increased by one if an offspring
    generation has been successfully populated even if a post-mating
    operator fails. Another variable \"rep\" will also be set to
    indicate the index of each population in the simulator. Note that
    populations in a simulator does not have to have the same
    generation number. You could reset a population's generation
    number by changing this variable.  Parameter gen can be set to a
    positive number, which is the number of generations to evolve. If
    a simulator starts at the beginning of a generation g (for example
    0), a simulator will stop at the beginning (instead of the end) of
    generation g + gen (for example gen). If gen is negative
    (default), the evolution will continue indefinitely, until all
    replicates are stopped by operators that return False at some
    point (these operators are called terminators). At the end of the
    evolution, the generations that each replicates have evolved are
    returned. Note that finalOps are applied to all applicable
    population, including those that have stopped before others.  If
    parameter dryrun is set to True, this function will print a
    description of the evolutionary process generated by function
    describeEvolProcess() and exits.

"; 

%ignore simuPOP::Simulator::apply(const opList &ops);

%feature("docstring") simuPOP::Simulator::vars "

Usage:

    x.vars(rep, subPop=[])

Details:

    Return the local namespace of the rep-th population, equivalent to
    x.Population(rep).vars(subPop).

"; 

%feature("docstring") simuPOP::Simulator::cmp "

Description:

    a Pyton function used to compare the simulator objects Note that
    mating schemes are not tested.

Usage:

    x.__cmp__(rhs)

"; 

%feature("docstring") simuPOP::SplitSubPops "

Details:

    Split a given list of subpopulations according to either sizes of
    the resulting subpopulations, proportion of individuals, or an
    information field. The resulting subpopulations will have the same
    name as the original subpopulation.

"; 

%feature("docstring") simuPOP::SplitSubPops::SplitSubPops "

Usage:

    SplitSubPops(subPops=ALL_AVAIL, sizes=[], proportions=[],
      names=[], randomize=True, begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, infoFields=[])

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
    stage). Please refer to operator BaseOperator for a detailed
    explanation for all parameters.

Note:

    Unlike operator Migrator, this operator does not require an
    information field such as migrate_to.

"; 

%feature("docstring") simuPOP::SplitSubPops::~SplitSubPops "

Description:

    destructor

Usage:

    x.~SplitSubPops()

"; 

%feature("docstring") simuPOP::SplitSubPops::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::SplitSubPops::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::SplitSubPops::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::Stat "

Details:

    Operator Stat calculates various statistics of the population
    being applied and sets variables in its local namespace. Other
    operators or functions can retrieve results from or evalulate
    expressions in this local namespace after Stat is applied.

"; 

%feature("docstring") simuPOP::Stat::Stat "

Usage:

    Stat(popSize=False, numOfMales=False, numOfAffected=False,
      numOfSegSites=[], numOfMutants=[], alleleFreq=[], heteroFreq=[],
      homoFreq=[], genoFreq=[], haploFreq=[], haploHeteroFreq=[],
      haploHomoFreq=[], sumOfInfo=[], meanOfInfo=[], varOfInfo=[],
      maxOfInfo=[], minOfInfo=[], LD=[], association=[],
      neutrality=[], structure=[], HWE=[], effectiveSize=[],
      vars=ALL_AVAIL, suffix=\"\", output=\"\", begin=0, end=-1, step=1,
      at=[], reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a Stat operator that calculates specified statistics of a
    population when it is applied to this population. This operator is
    by default applied after mating (parameter stage) and can be
    applied to specified replicates (parameter rep) at specified
    generations (parameter begin, end, step, and at). This operator
    does not produce any output (ignore parameter output) after
    statistics are calculated. Instead, it stores results in the local
    namespace of the population being applied. Other operators can
    retrieve these variables or evalulate expression directly in this
    local namespace. Please refer to operator BaseOperator for a
    detailed explanation of these common operator parameters.   Stat
    supports parameter subPops. It usually calculate the same set of
    statistics for all subpopulations (subPops=subPopList()). If a
    list of (virtual) subpopulations are specified, statistics for
    only specified subpopulations will be calculated. However,
    different statistics treat this parameter differently and it is
    very important to check its reference before you use subPops for
    any statistics.  Calculated statistics are saved as variables in a
    population's local namespace. These variables can be numbers,
    lists or dictionaries and can be retrieved using functions
    Population.vars() or Population.dvars(). A special default
    dictionary (defdict) is used for dictionaries whose keys are
    determined dynamically. Accessing elements of such a dictionary
    with an invalid key will yield value 0 instead of a KeyError. If
    the same variables are calculated for one or more (virtual)
    subpopulation, the variables are stored in
    vars()['subPop'][sp]['var'] where sp is a subpopulation ID (sp) or
    a tuple of virtual subpopulation ID ((sp, vsp)).
    Population.vars(sp) and Population.dvars(sp) provide shortcuts to
    these variables.  Operator Stat outputs a number of most useful
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
    different subpopulations (e.g. pairwise Fst).  Operator Stat
    supports the following statistics:  popSize: If popSize=True,
    number of individuals in all or specified subpopulations
    (parameter subPops) will be set to the following variables:
    *   popSize (default): Number of individuals in all or specified
    subpopulations. Because subPops does not have to cover all
    individuals, it may not be the actual population size.
    *   popSize_sp: Size of (virtual) subpopulation sp.
    *   subPopSize (default): A list of (virtual) subpopulation sizes.
    This variable is easier to use than accessing popSize from each
    (virtual) subpopulation.numOfMales: If numOfMales=True, number of
    male individuals in all or specified (virtual) subpopulations will
    be set to the following variables:
    *   numOfMales (default): Total number of male individuals in all
    or specified (virtual) subpopulations.
    *   numOfFemales (default): Total number of female individuals in
    all or specified (virtual) subpopulations.
    *   propOfMales: Proportion of male individuals.
    *   propOfFemales: Proportion of female individuals.
    *   numOfMales_sp: Number of male individuals in each (virtual)
    subpopulation.
    *   numOfFemales_sp: Number of female individuals in each
    (virtual) subpopulation.
    *   propOfMales_sp: Proportion of male individuals in each
    (virtual) subpopulation.
    *   propOfFemales_sp: Proportion of female individuals in each
    (virtual) subpopulation.numOfAffected: If numOfAffected=True,
    number of affected individuals in all or specified (virtual)
    subpopulations will be set to the following variables:
    *   numOfAffected (default): Total number of affected individuals
    in all or specified (virtual) subpopulations.
    *   numOfUnaffected (default): Total number of unaffected
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
    each (virtual) subpopulation.numOfSegSites: Parameter
    numOfSegSites accepts a list of loci (loci indexes, names, or
    ALL_AVAIL) and count the number of loci with at least two
    different alleles (segregating sites) or loci with only one non-
    zero allele (no zero allele, not segragating) for individuals in
    all or specified (virtual) subpopulations. This parameter sets
    variables
    *   numOfSegSites (default): Number of segregating sites in all or
    specified (virtual) subpopulations.
    *   numOfSegSites_sp: Number of segregating sites in each
    (virtual) subpopulation.
    *   numOfFixedSites: Number of sites with one non-zero allele in
    all or specified (virtual) subpopulations.
    *   numOfFixedSites_sp: Number of sites with one non-zero allele
    in in each (virtual) subpopulations.numOfMutants: Parameter
    numOfMutants accepts a list of loci (loci indexes, names, or
    ALL_AVAIL) and count the number of mutants (non-zero alleles) for
    individuals in all or specified (virtual) subpopulations. It sets
    variables
    *   numOfMutants (default): Number of mutants in all or specified
    (virtual) subpopulations.
    *   numOfMutants_sp: Number of mutants in each (virtual)
    subpopulations.alleleFreq: This parameter accepts a list of loci
    (loci indexes, names, or ALL_AVAIL), at which allele frequencies
    will be calculated. This statistic outputs the following
    variables, all of which are dictionary (with loci indexes as keys)
    of default dictionaries (with alleles as keys). For example,
    alleleFreq[loc][a] returns 0 if allele a does not exist.
    *   alleleFreq (default): alleleFreq[loc][a] is the frequency of
    allele a at locus for all or specified (virtual) subpopulations.
    *   alleleNum (default): alleleNum[loc][a] is the number of allele
    a at locus for all or specified (virtual) subpopulations.
    *   alleleFreq_sp: Allele frequency in each (virtual)
    subpopulation.
    *   alleleNum_sp: Allele count in each (virtual)
    subpopulation.heteroFreq and homoFreq: These parameters accept a
    list of loci (by indexes or names), at which the number and
    frequency of homozygotes and/or heterozygotes will be calculated.
    These statistics are only available for diploid populations. The
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
    loci (by indexes or names) at which number and frequency of all
    genotypes are outputed as a dictionary (indexed by loci indexes)
    of default dictionaries (indexed by tuples of possible indexes).
    This statistic is available for all population types with genotype
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
    subpopulation.haploHeteroFreq and haploHomoFreq: These parameters
    accept a list of haplotypes (list of loci), at which the number
    and frequency of haplotype homozygotes and/or heterozygotes will
    be calculated. Note that these statistics are observed count of
    haplotype heterozygote. The following variables will be outputted:
    *   haploHeteroFreq (default for parameter haploHeteroFreq): A
    dictionary of proportion of haplotype heterozygotes in all or
    specified (virtual) subpopulations, with haplotype indexes as
    dictionary keys.
    *   haploHomoFreq (default for parameter haploHomoFreq): A
    dictionary of proportion of homozygotes in all or specified
    (virtual) subpopulations.
    *   haploHeteroNum: A dictionary of number of heterozygotes in all
    or specified (virtual) subpopulations.
    *   haploHomoNum: A dictionary of number of homozygotes in all or
    specified (virtual) subpopulations.
    *   haploHeteroFreq_sp: A dictionary of proportion of
    heterozygotes in each (virtual) subpopulation.
    *   haploHomoFreq_sp: A dictionary of proportion of homozygotes in
    each (virtual) subpopulation.
    *   haploHeteroNum_sp: A dictionary of number of heterozygotes in
    each (virtual) subpopulation.
    *   haploHomoNum_sp: A dictionary of number of homozygotes in each
    (virtual) subpopulation.sumOfinfo, meanOfInfo, varOfInfo,
    maxOfInfo and minOfInfo: Each of these five parameters accepts a
    list of information fields. For each information field, the sum,
    mean, variance, maximum or minimal (depending on the specified
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
    loci, which can be a list of indexes, names, or ALL_AVAIL. At each
    locus, one or more statistical tests will be performed to test
    association between this locus and individual affection status.
    Currently, simuPOP provides the following tests:
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
    (detection of natural selection) on specified loci, which can be a
    list of loci indexes, names or ALL_AVAIL. It currently only
    outputs Pi, which is the average number of pairwise difference
    between loci. This statistic outputs the following variables:
    *   Pi Mean pairwise difference between all sequences from all or
    specified (virtual) subpopulations.
    *   Pi_sp Mean paiewise difference between all sequences in each
    (virtual) subpopulation.structure: Parameter structure accepts a
    list of loci at which statistics that measure population structure
    are calculated. structure accepts a list of loci indexes, names or
    ALL_AVAIL. This parameter currently supports the following
    statistics:
    *   Weir and Cockerham's Fst (1984). This is the most widely used
    estimator of Wright's fixation index and can be used to measure
    Population differentiation. However, this method is designed to
    estimate Fst from samples of larger populations and might not be
    appropriate for the calculation of Fst of large populations.
    *   Nei's Gst (1973). The Gst estimator is another estimator for
    Wright's fixation index but it is extended for multi-allele (more
    than two alleles) and multi-loci cases. This statistics should be
    used if you would like to obtain a true Fst value of a large
    Population. Nei's Gst uses only allele frequency information so it
    is available for all population type (haploid, diploid etc). Weir
    and Cockerham's Fst uses heterozygosity frequency so it is best
    for autosome of diploid populations. For non-diploid population,
    sex, and mitochondrial DNAs, simuPOP uses expected heterozygosity
    (1 - sum p_i^2) when heterozygosity is needed. These statistics
    output the following variables:
    *   F_st (default) The WC84 Fst statistic estimated for all *
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
    populations. HWE can be a list of loci indexes, names or
    ALL_AVAIL. This statistic outputs the following variables:
    *   HWE (default) A dictionary of p-values of HWE tests using
    genotypes in all or specified (virtual) subpopulations.
    *   HWE_sp A dictionary of p-values of HWS tests using genotypes
    in each (virtual) subpopulation.effectiveSize: Parameter
    effectiveSize accepts a list of loci at which the effective
    population size for the whole or specified (virtual)
    subpopulations is calculated. effectiveSize can be a list of loci
    indexes, names or ALL_AVAIL. This statistic outputs the following
    variables:
    *   Ne_waples89 (default) Effective population size, 2.5% and
    97.5% confidence interval as a list of size three, estimated using
    a temporal method as described in Waples 1989, Genetics. Because
    this is a temporal method, Ne_waples89 is set to census population
    size when it is first calculated, and is set to the temporal
    effecitive size between the present and last call afterward. The
    number of generations between samples for each estimate is
    therefore controlled by applicability parameters such as step and
    at. Because this implementation uses census population size as s_0
    and s_t, it is recommended that you apply this statistics to a
    random sample that retains population variables from the
    population from which the sample is drawn, or to a virtual
    subpopulation (e.g. using a rangeSplitter), to estimate Ne from
    samples. simuPOP assumes the samples are from a population with
    unknown size and only returns Ne under plan under study sampling
    plan 2.
    *   Ne_tempoFS (default) Effective population size, 2.5% and 97.5%
    confidence interval as a list of size three, estimated using a
    temporal method as described in Jorde & Ryman (2007), and as
    implemented by software tempoFS
    (http://www.zoologi.su.se/~ryman/). This variable is set to census
    population size when the operator is first applied, and to the
    temporal effective size between the present and the last
    generation when it is applied afterward. It does not make use of
    more than two samples, and assumes unknown census population size
    (sampling plan 2). A value of -1 implies an infinite estimate of
    effective size.

"; 

%feature("docstring") simuPOP::Stat::~Stat "

Usage:

    x.~Stat()

"; 

%feature("docstring") simuPOP::Stat::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::Stat::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::Stat::apply "Obsolete or undocumented function."

%ignore simuPOP::statAlleleFreq;

%feature("docstring") simuPOP::statAlleleFreq::statAlleleFreq "

Usage:

    statAlleleFreq(loci, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statAlleleFreq::describe "

Usage:

    x.describe(format=True)

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

%feature("docstring") simuPOP::statAssociation::describe "

Usage:

    x.describe(format=True)

"; 

%feature("docstring") simuPOP::statAssociation::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statEffectiveSize;

%feature("docstring") simuPOP::statEffectiveSize::statEffectiveSize "

Usage:

    statEffectiveSize(loci, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statEffectiveSize::describe "

Usage:

    x.describe(format=True)

"; 

%feature("docstring") simuPOP::statEffectiveSize::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statGenoFreq;

%feature("docstring") simuPOP::statGenoFreq::statGenoFreq "

Usage:

    statGenoFreq(genoFreq, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statGenoFreq::describe "

Usage:

    x.describe(format=True)

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

%feature("docstring") simuPOP::statHaploFreq::describe "

Usage:

    x.describe(format=True)

"; 

%feature("docstring") simuPOP::statHaploFreq::~statHaploFreq "

Usage:

    x.~statHaploFreq()

"; 

%feature("docstring") simuPOP::statHaploFreq::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statHaploHomoFreq;

%feature("docstring") simuPOP::statHaploHomoFreq::statHaploHomoFreq "

Usage:

    statHaploHomoFreq(haploHeteroFreq, haploHomoFreq, subPops, vars,
      suffix)

"; 

%feature("docstring") simuPOP::statHaploHomoFreq::~statHaploHomoFreq "

Usage:

    x.~statHaploHomoFreq()

"; 

%feature("docstring") simuPOP::statHaploHomoFreq::describe "

Usage:

    x.describe(format=True)

"; 

%feature("docstring") simuPOP::statHaploHomoFreq::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statHeteroFreq;

%feature("docstring") simuPOP::statHeteroFreq::statHeteroFreq "

Usage:

    statHeteroFreq(heteroFreq, homoFreq, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statHeteroFreq::describe "

Usage:

    x.describe(format=True)

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

%feature("docstring") simuPOP::statHWE::describe "

Usage:

    x.describe(format=True)

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

%feature("docstring") simuPOP::statInfo::describe "

Usage:

    x.describe(format=True)

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

%feature("docstring") simuPOP::statLD::describe "

Usage:

    x.describe(format=True)

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

%feature("docstring") simuPOP::statNeutrality::describe "

Usage:

    x.describe(format=True)

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

%feature("docstring") simuPOP::statNumOfAffected::describe "

Usage:

    x.describe(format=True)

"; 

%feature("docstring") simuPOP::statNumOfAffected::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statNumOfMales;

%feature("docstring") simuPOP::statNumOfMales::statNumOfMales "

Usage:

    statNumOfMales(numOfMales, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statNumOfMales::apply "

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::statNumOfMales::describe "

Usage:

    x.describe(format=True)

"; 

%ignore simuPOP::statNumOfMutants;

%feature("docstring") simuPOP::statNumOfMutants::statNumOfMutants "

Usage:

    statNumOfMutants(loci, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statNumOfMutants::describe "

Usage:

    x.describe(format=True)

"; 

%feature("docstring") simuPOP::statNumOfMutants::apply "

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statNumOfSegSites;

%feature("docstring") simuPOP::statNumOfSegSites::statNumOfSegSites "

Usage:

    statNumOfSegSites(loci, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statNumOfSegSites::describe "

Usage:

    x.describe(format=True)

"; 

%feature("docstring") simuPOP::statNumOfSegSites::apply "

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

%feature("docstring") simuPOP::statPopSize::describe "

Usage:

    x.describe(format=True)

"; 

%ignore simuPOP::statStructure;

%feature("docstring") simuPOP::statStructure::statStructure "

Usage:

    statStructure(Fst, subPops, vars, suffix)

"; 

%feature("docstring") simuPOP::statStructure::describe "

Usage:

    x.describe(format=True)

"; 

%feature("docstring") simuPOP::statStructure::apply "

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::StepwiseMutator "

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

%feature("docstring") simuPOP::StepwiseMutator::StepwiseMutator "

Usage:

    StepwiseMutator(rates=[], loci=ALL_AVAIL, incProb=0.5,
      maxAllele=0, mutStep=[], mapIn=[], mapOut=[], output=\">\",
      begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=\"ind_id\", lineageMode=FROM_INFO)

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
    *   (GEOMETRIC_DISTRIBUTION, p): The number of steps follows a a
    geometric distribution with parameter p.
    *   A Python function: This user defined function accepts the
    allele being mutated and return the steps to mutate. The mutation
    process is usually neutral in the sense that mutating up and down
    is equally likely. You can adjust parameter incProb to change this
    behavior.  If you need to use other generalized stepwise mutation
    models, you can implement them using a PyMutator. If performance
    becomes a concern, I may add them to this operator if provided
    with a reliable reference.

"; 

%feature("docstring") simuPOP::StepwiseMutator::~StepwiseMutator "

Usage:

    x.~StepwiseMutator()

"; 

%ignore simuPOP::StepwiseMutator::mutate(Allele allele, size_t locus) const;

%feature("docstring") simuPOP::StepwiseMutator::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::StepwiseMutator::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::StopEvolution "

Description:

    exception, throw if an operator would like to stop all replicates.

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

%ignore simuPOP::stringList::pushback(const string &str);

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

%ignore simuPOP::subPopList::allAvail() const;

%ignore simuPOP::subPopList::empty() const;

%ignore simuPOP::subPopList::size() const;

%feature("docstring") simuPOP::subPopList::len "

Usage:

    x.__len__()

"; 

%ignore simuPOP::subPopList::pushback(const vspID subPop);

%ignore simuPOP::subPopList::contains(const vspID subPop) const;

%ignore simuPOP::subPopList::overlap(const size_t subPop) const;

%ignore simuPOP::subPopList::begin() const;

%ignore simuPOP::subPopList::end() const;

%feature("docstring") simuPOP::subPopList::expandFrom "

Description:

    expand ALL_AVAIL and [(ALL_AVAIL, vsp), ...] according to pop

Usage:

    x.expandFrom(pop)

"; 

%feature("docstring") simuPOP::SummaryTagger "

Details:

    A summary tagger summarize values of one or more parental
    information field to another information field of an offspring. If
    mating is sexual, two sets of parental values will be involved.

"; 

%feature("docstring") simuPOP::SummaryTagger::SummaryTagger "

Usage:

    SummaryTagger(mode=MEAN, begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, output=\"\", infoFields=[])

Details:

    Creates a summary tagger that summarize values of one or more
    parental information field (infoFields[:-1]) to an offspring
    information field (infoFields[-1]). A parameter mode specifies how
    to pass summarize parental values. More specifically,
    *   mode=MEAN Passing the average of values.
    *   mode=MAXIMUM Passing the maximum value of values.
    *   mode=Minumum Passing the minimum value of values.
    *   mode=SUMMATION Passing the sum of values.
    *   mode=MULTIPLICATION Passing the multiplication of values. This
    operator does not support parameter subPops and does not output
    any information.

"; 

%feature("docstring") simuPOP::SummaryTagger::~SummaryTagger "

Usage:

    x.~SummaryTagger()

"; 

%feature("docstring") simuPOP::SummaryTagger::describe "Obsolete or undocumented function."

%ignore simuPOP::SummaryTagger::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::SummaryTagger::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::SystemError "

Description:

    exception, thrown if system error occurs

"; 

%feature("docstring") simuPOP::SystemError::SystemError "

Usage:

    SystemError(msg)

"; 

%feature("docstring") simuPOP::TerminateIf "

Details:

    This operator evaluates an expression in a population's local
    namespace and terminate the evolution of this population, or the
    whole simulator, if the return value of this expression is True.
    Termination caused by an operator will stop the execution of all
    operators after it. The generation at which the population is
    terminated will be counted in the evolved generations (return
    value from Simulator::evolve) if termination happens after mating.

"; 

%feature("docstring") simuPOP::TerminateIf::TerminateIf "

Usage:

    TerminateIf(condition=\"\", stopAll=False, message=\"\", output=\"\",
      begin=0, end=-1, step=1, at=[], reps=ALL_AVAIL,
      subPops=ALL_AVAIL, infoFields=[])

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

%feature("docstring") simuPOP::TerminateIf::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::TerminateIf::describe "Obsolete or undocumented function."

%feature("docstring") simuPOP::TerminateIf::apply "Obsolete or undocumented function."

%feature("docstring") simuPOP::TerminateIf::~TerminateIf "

Usage:

    x.~TerminateIf()

"; 

%feature("docstring") simuPOP::TicToc "

Details:

    This operator, when called, output the difference between current
    and the last called clock time. This can be used to estimate
    execution time of each generation. Similar information can also be
    obtained from turnOnDebug(\"DBG_PROFILE\"), but this operator has
    the advantage of measuring the duration between several
    generations by setting step parameter. As an advanced feature that
    mainly used for performance testing, this operator accepts a
    parameter stopAfter (seconds), and will stop the evolution of a
    population if the overall time exceeds stopAfter. Note that
    elapsed time is only checked when this operator is applied to a
    population so it might not be able to stop the evolution process
    right after stopAfter seconds. This operator can also be applied
    during mating. Note that to avoid excessive time checking, this
    operator does not always check system time accurately.

"; 

%feature("docstring") simuPOP::TicToc::TicToc "

Usage:

    TicToc(output=\">\", stopAfter=0, begin=0, end=-1, step=1, at=[],
      reps=ALL_AVAIL, subPops=ALL_AVAIL, infoFields=[])

Details:

    Create a TicToc operator that outputs the elapsed since the last
    time it was applied, and the overall time since the first time
    this operator is applied.

"; 

%feature("docstring") simuPOP::TicToc::~TicToc "

Description:

    destructor

Usage:

    x.~TicToc()

"; 

%feature("docstring") simuPOP::TicToc::clone "Obsolete or undocumented function."

%feature("docstring") simuPOP::TicToc::apply "Obsolete or undocumented function."

%ignore simuPOP::TicToc::applyDuringMating(Population &pop, Population &offPop, RawIndIterator offspring, Individual *dad=NULL, Individual *mom=NULL) const;

%feature("docstring") simuPOP::TicToc::describe "Obsolete or undocumented function."

%ignore simuPOP::TicToc::parallelizable() const;

%feature("docstring") simuPOP::uintList "

"; 

%feature("docstring") simuPOP::uintList::uintList "

Usage:

    uintList(obj=Py_True)

"; 

%ignore simuPOP::uintList::uintList(const vectoru &values);

%ignore simuPOP::uintList::allAvail() const;

%feature("docstring") simuPOP::uintList::unspecified "

Usage:

    x.unspecified()

"; 

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

%ignore simuPOP::UniformNumOffModel;

%feature("docstring") simuPOP::UniformNumOffModel::UniformNumOffModel "

Usage:

    UniformNumOffModel(low, high)

"; 

%feature("docstring") simuPOP::UniformNumOffModel::clone "

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::UniformNumOffModel::getNumOff "

Usage:

    x.getNumOff(ssize_t)

"; 

%feature("docstring") simuPOP::ValueError "

Description:

    exception, thrown if value of range etc

"; 

%feature("docstring") simuPOP::ValueError::ValueError "

Usage:

    ValueError(msg)

"; 

%ignore simuPOP::vspFunctor;

%feature("docstring") simuPOP::vspFunctor::vspFunctor "

Usage:

    vspFunctor()

"; 

%feature("docstring") simuPOP::vspID "

Details:

    A class to specify virtual subpopulation, which is composed of a
    subpopulation ID and a virtual subpopulation ID.

"; 

%feature("docstring") simuPOP::vspID::vspID "

Usage:

    vspID(id)

Details:

    Create a subpopulation id. Accept id as well as names.

"; 

%ignore simuPOP::vspID::vspID(const vectori &subPop, bool allAvailSP=false, bool allAvailVSP=false, const string &spName=string(), const string &vspName=string());

%feature("docstring") simuPOP::vspID::~vspID "

Usage:

    x.~vspID()

"; 

%ignore simuPOP::vspID::subPop() const;

%ignore simuPOP::vspID::virtualSubPop() const;

%ignore simuPOP::vspID::valid() const;

%ignore simuPOP::vspID::isVirtual() const;

%ignore simuPOP::vspID::allAvailSP() const;

%ignore simuPOP::vspID::allAvailVSP() const;

%feature("docstring") simuPOP::vspID::resolve "

Usage:

    x.resolve(pop)

"; 

%ignore simuPOP::vspID::spName() const;

%ignore simuPOP::vspID::vspName() const;

%feature("docstring") simuPOP::WeightedSampler "

Details:

    A random number generator that returns 0, 1, ..., k-1 with
    probabilites that are proportional to their weights. For example,
    a weighted sampler with weights 4, 3, 2 and 1 will return numbers
    0, 1, 2 and 3 with probabilities 0.4, 0.3, 0.2 and 0.1,
    respectively. If an additional parameter N is specified, the
    weighted sampler will return exact proportions of numbers if N
    numbers are returned. The version without additional parameter is
    similar to the sample(prob, replace=FALSE) function of the R
    statistical package.

"; 

%feature("docstring") simuPOP::WeightedSampler::WeightedSampler "

Usage:

    WeightedSampler(weights=[], N=0)

Details:

    Creates a weighted sampler that returns 0, 1, ... k-1 when a list
    of k weights are specified (weights). weights do not have to add
    up to 1. If a non-zero N is specified, exact proportions of
    numbers will be returned in N returned numbers.

"; 

%feature("docstring") simuPOP::WeightedSampler::~WeightedSampler "

Description:

    destructor

Usage:

    x.~WeightedSampler()

"; 

%ignore simuPOP::WeightedSampler::set(IT first, IT last, size_t N=0);

%feature("docstring") simuPOP::WeightedSampler::draw "

Usage:

    x.draw()

Details:

    Returns a random number between 0 and k-1 with probabilities that
    are proportional to specified weights.

"; 

%feature("docstring") simuPOP::WeightedSampler::drawSamples "

Usage:

    x.drawSamples(n=1)

Details:

    Returns a list of n random numbers

"; 

%feature("docstring") simuPOP::applyDuringMatingOperator "Obsolete or undocumented function."

%feature("docstring") simuPOP::loadPedigree "

Usage:

    loadPedigree(file, idField=\"ind_id\", fatherField=\"father_id\",
      motherField=\"mother_id\", ploidy=2, loci=[], chromTypes=[],
      lociPos=[], chromNames=[], alleleNames=[], lociNames=[],
      subPopNames=[], infoFields=[])

Details:

    Load a pedigree from a file saved by operator PedigreeTagger or
    function Pedigree.save. This file contains the ID of each
    offspring and their parent(s) and optionally sex ('M' or 'F'),
    affection status ('A' or 'U'), values of information fields and
    genotype at some loci. IDs of each individual and their parents
    are loaded to information fields idField, fatherField and
    motherField. Only numeric IDs are allowed, and individual IDs must
    be unique across all generations.  Because this file does not
    contain generation information, generations to which offspring
    belong are determined by the parent-offspring relationships.
    Individuals without parents are assumed to be in the top-most
    ancestral generation. This is the case for individuals in the top-
    most ancestral generation if the file is saved by function
    Pedigree.save(), and for individuals who only appear as another
    individual's parent, if the file is saved by operator
    PedigreeTagger. The order at which offsprng is specified is not
    important because this function essentially creates a top-most
    ancestral generation using IDs without parents, and creates the
    next generation using offspring of these parents, and so on until
    all generations are recreated. That is to say, if you have a
    mixture of pedigrees with different generations, they will be
    lined up from the top most ancestral generation.  If individual
    sex is not specified, sex of of parents are determined by their
    parental roles (father or mother) but the sex of individuals in
    the last generation can not be determined so they will all be
    males. If additional information fields are given, their names
    have to be specified using parameter infoFields. The rest of the
    columns are assued to be alleles, arranged ploidy consecutive
    columns for each locus. If paraemter loci is not specified, the
    number of loci is calculated by number of columns divided by
    ploidy (default to 2). All loci are assumed to be on one
    chromosome unless parameter loci is used to specified number of
    loci on each chromosome. Additional parameters such as ploidy,
    chromTypes, lociPos, chromNames, alleleNames, lociNames could be
    used to specified the genotype structured of the loaded pedigree.
    Please refer to class Population for details about these
    parameters.

"; 

%feature("docstring") simuPOP::loadPopulation "

Usage:

    loadPopulation(file)

Details:

    load a population from a file saved by Population::save().

"; 

%feature("docstring") simuPOP::describeEvolProcess "

Usage:

    describeEvolProcess(initOps=[], preOps=[],
      matingScheme=MatingScheme, postOps=[], finalOps=[], gen=-1,
      numRep=1)

Details:

    This function takes the same parameters as Simulator.evolve and
    output a description of how an evolutionary process will be
    executed. It is recommended that you call this function if you
    have any doubt how your simulation will proceed.

"; 

%feature("docstring") simuPOP::turnOnDebug "

Usage:

    turnOnDebug(code=\"\")

Details:

    Set debug code code. More than one code could be specified using a
    comma separated string. Name of available codes are available from
    moduleInfo()['debug'].keys().

"; 

%feature("docstring") simuPOP::turnOffDebug "

Usage:

    turnOffDebug(code=\"DBG_ALL\")

Details:

    Turn off debug code code. More than one code could be specified
    using a comma separated string. Default to turn off all debug
    codes.

"; 

%ignore simuPOP::debug(DBG_CODE code);

%ignore simuPOP::repeatedWarning(const string &message);

%ignore simuPOP::initClock();

%feature("docstring") simuPOP::elapsedTime "

Usage:

    elapsedTime(name)

"; 

%feature("docstring") simuPOP::setOptions "

Usage:

    setOptions(numThreads=-1, name=None, seed=0)

Details:

    First argument is to set number of thread in openMP. The number of
    threads can be be positive, integer (number of threads) or 0,
    which implies all available cores, or a number set by
    environmental variable OMP_NUM_THREADS. Second and third argument
    is to set the type or seed of existing random number generator
    using RNGname with seed. If using openMP, it sets the type or seed
    of random number generator of each thread.

"; 

%ignore simuPOP::numThreads();

%ignore simuPOP::fetchAndIncrement(ATOMICLONG *val);

%ignore simuPOP::parallelSort(T1 start, T1 end, T2 cmp);

%ignore simuPOP::simuPOPkbhit();

%ignore simuPOP::simuPOPgetch();

%ignore simuPOP::PyObjAsBool(PyObject *obj, bool &val);

%ignore simuPOP::PyObjAsInt(PyObject *obj, long &val);

%ignore simuPOP::PyObjAsDouble(PyObject *obj, double &val);

%ignore simuPOP::PyObjAsString(PyObject *obj, string &val);

%ignore simuPOP::PyObjAsArray(PyObject *obj, vectorf &val);

%ignore simuPOP::PyObjAsIntArray(PyObject *obj, vectori &val);

%ignore simuPOP::AlleleVecAsNumArray(GenoIterator begin, GenoIterator end);

%ignore simuPOP::LineageVecAsNumArray(LineageIterator begin, LineageIterator end);

%ignore simuPOP::PyObjAsString(PyObject *str);

%ignore simuPOP::mainVars();

%ignore simuPOP::moduleVars();

%ignore simuPOP::pyPopObj(void *p);

%ignore simuPOP::pyIndObj(void *p);

%ignore simuPOP::pyIndPointer(PyObject *p);

%ignore simuPOP::pyPopPointer(PyObject *p);

%ignore simuPOP::shorten(const string &val, size_t length=40);

%ignore simuPOP::ostreamManager();

%feature("docstring") simuPOP::closeOutput "

Usage:

    closeOutput(output=\"\")

Details:

    Output files specified by '>' are closed immediately after they
    are written. Those specified by '>>' and '>>>' are closed by a
    simulator after Simulator.evolve(). However, these files will be
    kept open if the operators are applied directly to a population
    using the operators' function form. In this case, function
    closeOutput can be used to close a specific file output, and close
    all unclosed files if output is unspecified. An exception will be
    raised if output does not exist or it has already been closed.

"; 

%feature("docstring") simuPOP::getRNG "

Description:

    return the currently used random number generator

Usage:

    getRNG()

"; 

%ignore simuPOP::chisqTest(const vector< vectoru > &table, double &chisq, double &chisq_p);

%ignore simuPOP::armitageTrendTest(const vector< vectoru > &table, const vectorf &weight);

%ignore simuPOP::hweTest(const vectoru &cnt);

%ignore simuPOP::propToCount(IT first, IT last, size_t N, vectoru &count);

%ignore simuPOP::formatDescription(const string &text);

%feature("docstring") simuPOP::moduleInfo "

Usage:

    moduleInfo()

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
    *   wordsize: size of word, can be either 32 or 64.
    *   alleleBits: the number of bits used to store an allele
    *   maxNumSubPop: maximum number of subpopulations.
    *   maxIndex: maximum index size (limits population size * total
    number of marker).
    *   debug: A dictionary with debugging codes as keys and the
    status of each debugging code (True or False) as their values.

"; 

%ignore simuPOP::initialize();

%ignore simuPOP::cnull();

%ignore std::powthree(unsigned n);

%feature("docstring") simuPOP::Population::dvars "

Usage:

    x.dvars(subPop=[])

Details:

    Return a wrapper of Python dictionary returned by vars(subPop) so
    that dictionary keys can be accessed as attributes.

"; 

%feature("docstring") simuPOP::Simulator::dvars "

Usage:

    x.dvars(rep, subPop=[])

Details:

    Return a wrapper of Python dictionary returned by vars(rep,
    subPop) so that dictionary keys can be accessed as attributes.

"; 

