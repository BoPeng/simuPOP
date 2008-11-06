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

%ignore simuPOP::affectionSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const ;

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

%feature("docstring") simuPOP::affectionTagger "

Description:

    Tagging affection status.

Details:

    This is a simple post-mating tagger that write affection status to
    a file. By default, 1 for unaffected, 2 for affected.

"; 

%feature("docstring") simuPOP::affectionTagger::affectionTagger "

Usage:

    affectionTagger(code=[], begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, stage=PostMating, output=\">\", outputExpr=\"\",
      infoFields=[])

Arguments:

    code:           code for Male and Female, default to 1 and 2,
                    respectively. This is used by Linkage format.

"; 

%feature("docstring") simuPOP::affectionTagger::apply "

Description:

    add a newline

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::alphaMating "

Applicability: diploid only

Description:

    Only a number of alpha individuals can mate with individuals of
    opposite sex.

Details:

    This mating scheme is composed of an random parents chooser with
    alpha individuals, and a Mendelian offspring generator. That is to
    say, a certain number of alpha individual (male or female) are
    determined by alphaNum or an information field. Then, only these
    alpha individuals are able to mate with random individuals of
    opposite sex.

"; 

%feature("docstring") simuPOP::alphaMating::alphaMating "

Usage:

    alphaMating(alphaSex=Male, alphaNum=0, alphaField=string,
      numOffspring=1., numOffspringFunc=None, maxNumOffspring=0,
      mode=MATE_NumOffspring, sexParam=0.5, sexMode=MATE_RandomSex,
      newSubPopSize=[], newSubPopSizeFunc=None, newSubPopSizeExpr=\"\",
      subPop=[], weight=0)

Details:

    Please refer to class mating for descriptions of other parameters.
    Note: If selection is enabled, it works regularly on on-alpha sex,
    but works twice on alpha sex. That is to say, alphaNum alpha
    indiviudals are chosen selectively, and selected again during
    mating.

Arguments:

    alphaSex:       the sex of the alpha individual, i.e. alpha male
                    or alpha female who be the only mating individuals
                    in their sex group.
    alphaNum:       Number of alpha individuals. If infoField is not
                    given, alphaNum random individuals with alphaSex
                    will be chosen. If selection is enabled,
                    individuals with higher+ fitness values have
                    higher probability to be selected. There is by
                    default no alpha individual (alphaNum = 0).
    alphaField:     if an information field is given, individuals with
                    non-zero values at this information field are
                    alpha individuals. Note that these individuals
                    must have alphaSex.

"; 

%feature("docstring") simuPOP::alphaMating::clone "

Description:

    deep copy of a random mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::alphaMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the random mating scheme

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::baseOperator "

Description:

    base class of all classes that manipulate populations

Details:

    Operators are objects that act on populations. They can be applied
    to populations directly using their function forms, but they are
    usually managed and applied by a simulator.
    There are three kinds of operators:
    *  built-in: written in C++, the fastest. They do not interact
    with Python shell except that some of them set variables that are
    accessible from Python.
    *  hybrid: written in C++ but calls a Python function during
    execution. Less efficient. For example, a hybrid mutator
    pyMutator will go through a population and mutate alleles with
    given mutation rate. How exactly the allele will be mutated is
    determined by a user-provided Python function. More specifically,
    this operator will pass the current allele to a user-provided
    Python function and take its return value as the mutant allele.
    *  pure Python: written in Python. The same speed as Python. For
    example, a varPlotter can plot Python variables that are set by
    other operators. Usually, an individual or a population object is
    passed to a user-provided Python function. Because arbitrary
    operations can be performed on the passed object, this operator is
    very flexible. Operators can be applied at different stages of the
    life cycle of a generation. It is possible for an operator to
    apply multiple times in a life cycle. For example, a
    savePopulation operator might be applied before and after mating
    to trace parental information. More specifically, operators can be
    applied at pre-, during-, post-mating, or a combination of these
    stages. Applicable stages are usually set by default but you can
    change it by setting stage=(PreMating|PostMating|DuringMating|PreP
    ostMating|PreDuringMating|DuringPostMating) parameter. Some
    operators ignore stage parameter because they only work at one
    stage.
    Operators do not have to be applied at all generations. You can
    specify starting and/or ending generations (parameters start,
    end), gaps between applicable generations (parameter step), or
    specific generations (parameter at). For example, you might want
    to start applying migrations after certain burn-in generations, or
    calculate certain statistics only sparsely. Generation numbers can
    be counted from the last generation, using negative generation
    numbers.
    Most operators are applied to every replicate of a simulator
    during evolution. Operators can have outputs, which can be
    standard (terminal) or a file. Output can vary with replicates
    and/or generations, and outputs from different operators can be
    accumulated to the same file to form table-like outputs.
    Filenames can have the following format:
    *  'filename' this file will be overwritten each time. If two
    operators output to the same file, only the last one will succeed;
    *  '>filename' the same as 'filename';
    *  '>>filename' the file will be created at the beginning of
    evolution ( simulator::evolve) and closed at the end. Outputs from
    several operators are appended;
    *  '>>>filename' the same as '>>filename' except that the file
    will not be cleared at the beginning of evolution if it is not
    empty;
    *  '>' standard output (terminal);
    *  '' suppress output. The output filename does not have to be
    fixed. If parameter outputExpr is used (parameter output will be
    ignored), it will be evaluated when a filename is needed. This is
    useful when you need to write different files for different
    replicates/generations.

"; 

%feature("docstring") simuPOP::baseOperator::baseOperator "

Description:

    common interface for all operators (this base operator does
    nothing by itself.)

Usage:

    baseOperator(output, outputExpr, stage, begin, end, step, at,
      rep, infoFields)

Arguments:

    begin:          the starting generation. Default to 0. A negative
                    number is allowed.
    end:            stop applying after this generation. A negative
                    numbers is allowed.
    step:           the number of generations between active
                    generations. Default to 1.
    at:             an array of active generations. If given, stage,
                    begin, end, and step will be ignored.
    rep:            applicable replicates. It can be a valid replicate
                    number, REP_ALL (all replicates, default), or
                    REP_LAST (only the last replicate). REP_LAST is
                    useful in adding newlines to a table output.
    output:         a string of the output filename. Different
                    operators will have different default output (most
                    commonly '>' or '').
    outputExpr:     an expression that determines the output filename
                    dynamically. This expression will be evaluated
                    against a population's local namespace each time
                    when an output filename is required. For example,
                    \"'>>out%s_%s.xml' % (gen, rep)\"  will output to
                    >>>out1_1.xml for replicate 1 at generation 1.

Note:

    * Negative generation numbers are allowed for parameters begin,
    end and at. They are interpreted as endGen + gen + 1. For example,
    begin = -2 in simu.evolve(..., end=20) starts at generation 19.
    * REP_ALL, REP_LAST are special constant that can only be used in
    the constructor of an operator. That is to say, explicit test of
    rep() == REP_LAST will not work.

Example:

Testsrc_operator.log Common features of all operators 

"; 

%feature("docstring") simuPOP::baseOperator::~baseOperator "

Description:

    destroy an operator

Usage:

    x.~baseOperator()

"; 

%feature("docstring") simuPOP::baseOperator::clone "

Description:

    deep copy of an operator

Usage:

    x.clone()

"; 

%ignore simuPOP::baseOperator::isActive(UINT rep, UINT numRep, long gen, long end, bool repOnly=false);

%ignore simuPOP::baseOperator::applicableReplicate();

%ignore simuPOP::baseOperator::setApplicableReplicate(int rep);

%ignore simuPOP::baseOperator::setActiveGenerations(int begin=0, int end=-1, int step=1, vectorl at=vectorl());

%ignore simuPOP::baseOperator::setApplicableStage(int stage);

%ignore simuPOP::baseOperator::canApplyPreMating();

%ignore simuPOP::baseOperator::canApplyDuringMating();

%ignore simuPOP::baseOperator::canApplyPostMating();

%ignore simuPOP::baseOperator::canApplyPreOrPostMating();

%ignore simuPOP::baseOperator::isCompatible(const population &pop);

%feature("docstring") simuPOP::baseOperator::haploidOnly "

Description:

    determine if the operator can be applied only for haploid
    population

Usage:

    x.haploidOnly()

"; 

%feature("docstring") simuPOP::baseOperator::diploidOnly "

Description:

    determine if the operator can be applied only for diploid
    population

Usage:

    x.diploidOnly()

"; 

%ignore simuPOP::baseOperator::setHaploidOnly();

%ignore simuPOP::baseOperator::setDiploidOnly();

%feature("docstring") simuPOP::baseOperator::infoSize "

Description:

    get the length of information fields for this operator

Usage:

    x.infoSize()

"; 

%feature("docstring") simuPOP::baseOperator::infoField "

Description:

    get the information field specified by user (or by default)

Usage:

    x.infoField(idx)

"; 

%ignore simuPOP::baseOperator::formOffGenotype();

%ignore simuPOP::baseOperator::setFormOffGenotype(bool flag=true);

%ignore simuPOP::baseOperator::applyWithScratch(population &pop, population &scratch, int stage);

%feature("docstring") simuPOP::baseOperator::apply "

Description:

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::baseOperator::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%ignore simuPOP::baseOperator::setOutput(string output="", string outputExpr="");

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

%feature("docstring") simuPOP::baseRandomMating "

Applicability: diploid only

Details:

    This base class defines a general random mating scheme that makes
    full use of a general random parents chooser, and a Mendelian
    offspring generator. A general random parents chooser allows
    selection without replacement, polygemous parents selection (a
    parent with more than one partners), and the definition of several
    alpha individuals. Direct use of this mating scheme is not
    recommended.  randomMating, monogemousMating, polygemousMating,
    alphaMating are all special cases of this mating scheme. They
    should be used whenever possible.

"; 

%feature("docstring") simuPOP::baseRandomMating::baseRandomMating "

Usage:

    baseRandomMating(replacement=True, replenish=False,
      polySex=Male, polyNum=1, alphaSex=Male, alphaNum=0,
      alphaField=string, numOffspring=1., numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, sexParam=0.5,
      sexMode=MATE_RandomSex, newSubPopSize=[], newSubPopSizeExpr=\"\",
      newSubPopSizeFunc=None, contWhenUniSex=True, subPop=[],
      weight=0)

Arguments:

    replacement:    If set to True, a parent can be chosen to mate
                    again. Default to False.
    replenish:      In case that replacement=True, whether or not
                    replenish a sex group when it is exhausted.
    polySex:        sex of polygamous mating. Male for polygyny,
                    Female for polyandry.
    polyNum:        Number of sex partners.
    alphaSex:       the sex of the alpha individual, i.e. alpha male
                    or alpha female who be the only mating individuals
                    in their sex group.
    alphaNum:       Number of alpha individuals. If infoField is not
                    given, alphaNum random individuals with alphaSex
                    will be chosen. If selection is enabled,
                    individuals with higher+ fitness values have
                    higher probability to be selected. There is by
                    default no alpha individual (alphaNum = 0).
    alphaField:     if an information field is given, individuals with
                    non-zero values at this information field are
                    alpha individuals. Note that these individuals
                    must have alphaSex.

"; 

%feature("docstring") simuPOP::baseRandomMating::~baseRandomMating "

Description:

    destructor

Usage:

    x.~baseRandomMating()

"; 

%feature("docstring") simuPOP::baseRandomMating::clone "

Description:

    deep copy of a random mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::baseRandomMating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::baseRandomMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the random mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::baseRandomMating::mateSubPop(population &pop, SubPopID subPop, RawIndIterator offBegin, RawIndIterator offEnd, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::BernulliTrials "

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

%feature("docstring") simuPOP::BernulliTrials::trialSize "

Description:

    return size of trial

Usage:

    x.trialSize()

"; 

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

%feature("docstring") simuPOP::binomialSelection "

Applicability: all ploidy

Description:

    a mating scheme that uses binomial selection, regardless of sex

Details:

    No sex information is involved (binomial random selection).
    Offspring is chosen from parental generation by random or
    according to the fitness values. In this mating scheme,
    *  numOffspring protocol is honored;
    *  population size changes are allowed;
    *  selection is possible;
    *  haploid population is allowed.

"; 

%feature("docstring") simuPOP::binomialSelection::binomialSelection "

Description:

    create a binomial selection mating scheme

Usage:

    binomialSelection(numOffspring=1., numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, sexParam=0.5,
      sexMode=MATE_RandomSex, newSubPopSize=[], newSubPopSizeExpr=\"\",
      newSubPopSizeFunc=None, subPop=[], weight=0)

Details:

    Please refer to class mating for parameter descriptions.

"; 

%feature("docstring") simuPOP::binomialSelection::~binomialSelection "

Description:

    destructor

Usage:

    x.~binomialSelection()

"; 

%feature("docstring") simuPOP::binomialSelection::clone "

Description:

    deep copy of a binomial selection mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::binomialSelection::__repr__ "

Description:

    used by Python print function to print out the general information
    of the binomial selection mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::binomialSelection::mateSubPop(population &pop, SubPopID subPop, RawIndIterator offBegin, RawIndIterator offEnd, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::cloneMating "

Applicability: all ploidy

Description:

    a clone mating that copy everyone from parental to offspring
    generation.

Details:

    Note that
    *  selection is not considered (fitness is ignored)
    *  sequentialParentMating is used. If offspring (virtual)
    subpopulation size is smaller than parental subpopulation size,
    not all parents will be cloned. If offspring (virtual)
    subpopulation size is larger, some parents will be cloned more
    than once.
    *  numOffspring interface is respected.
    *  during mating operators are applied.

"; 

%feature("docstring") simuPOP::cloneMating::cloneMating "

Description:

    create a binomial selection mating scheme

Usage:

    cloneMating(numOffspring=1., numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, sexParam=0.5,
      sexMode=MATE_RandomSex, newSubPopSize=[], newSubPopSizeExpr=\"\",
      newSubPopSizeFunc=None, subPop=[], weight=0)

Details:

    Please refer to class mating for parameter descriptions.

"; 

%feature("docstring") simuPOP::cloneMating::~cloneMating "

Description:

    destructor

Usage:

    x.~cloneMating()

"; 

%feature("docstring") simuPOP::cloneMating::clone "

Description:

    deep copy of a binomial selection mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::cloneMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the binomial selection mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::cloneMating::mateSubPop(population &pop, SubPopID subPop, RawIndIterator offBegin, RawIndIterator offEnd, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::cloneOffspringGenerator "

Applicability: all ploidy

Details:

    clone offspring generator copies parental geneotype to a number of
    offspring. Only one parent is accepted. The number of offspring
    produced is controled by parameters numOffspring,
    numOffspringFunc, maxNumOffspring and mode. Parameters sexParam
    and sexMode is ignored.

"; 

%feature("docstring") simuPOP::cloneOffspringGenerator::cloneOffspringGenerator "

Usage:

    cloneOffspringGenerator(numOffspring=1, numOffspringFunc=None,
      maxNumOffspring=1, mode=MATE_NumOffspring, sexParam=0.5,
      sexMode=MATE_RandomSex)

Arguments:

    sexParam:       ignored because sex is copied from the parent.
    sexMode:        ignored because sex is copied from the parent.

"; 

%feature("docstring") simuPOP::cloneOffspringGenerator::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::cloneOffspringGenerator::generateOffspring(population &pop, individual *dad, individual *mom, RawIndIterator &offBegin, RawIndIterator &offEnd, vector< baseOperator * > &ops);

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

%ignore simuPOP::combinedSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const ;

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

%feature("docstring") simuPOP::consanguineousMating "

Applicability: diploid only

Description:

    a mating scheme of consanguineous mating

Details:

    In this mating scheme, a parent is choosen randomly and mate with
    a relative that has been located and written to a number of
    information fields.

"; 

%feature("docstring") simuPOP::consanguineousMating::consanguineousMating "

Description:

    create a consanguineous mating scheme

Usage:

    consanguineousMating(relativeFields=[], func=None, param=None,
      replacement=False, replenish=True, numOffspring=1.,
      numOffspringFunc=None, maxNumOffspring=0,
      mode=MATE_NumOffspring, sexParam=0.5, sexMode=MATE_RandomSex,
      newSubPopSize=[], newSubPopSizeFunc=None, newSubPopSizeExpr=\"\",
      contWhenUniSex=True, subPop=[], weight=0)

Details:

    This mating scheme randomly choose a parent and then choose
    his/her spouse from indexes stored in infoFields.  Please refer to
    infoParentsChooser and  mendelianOffspringGenerator for other
    parameters.

Arguments:

    relativeFields: The information fields that stores indexes to
                    other individuals in a population. If more than
                    one valid (positive value) indexes exist, a random
                    index will be chosen. (c.f.  infoParentsChooser )
                    If there is no individual having any valid index,
                    the second parent will be chosen randomly from the
                    whole population.
    func:           A python function that can be used to prepare the
                    indexes of these information fields. For example,
                    functions population::locateRelatives and/or
                    population::setIndexesOfRelatives can be used to
                    locate certain types of relatives of each
                    individual.
    param:          An optional parameter that can be passed to func.

"; 

%feature("docstring") simuPOP::consanguineousMating::~consanguineousMating "

Description:

    destructor

Usage:

    x.~consanguineousMating()

"; 

%ignore simuPOP::consanguineousMating::consanguineousMating(const consanguineousMating &rhs);

%feature("docstring") simuPOP::consanguineousMating::clone "

Description:

    deep copy of a consanguineous mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::consanguineousMating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::consanguineousMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the consanguineous mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::consanguineousMating::preparePopulation(population &pop);

%ignore simuPOP::consanguineousMating::mateSubPop(population &pop, SubPopID subPop, RawIndIterator offBegin, RawIndIterator offEnd, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::continueIf "

Description:

    terminate according to a condition failure

Details:

    The same as  terminateIf but continue if the condition is True.

"; 

%feature("docstring") simuPOP::continueIf::continueIf "

Description:

    create a  continueIf terminator

Usage:

    continueIf(condition=\"\", message=\"\", var=\"terminate\", output=\"\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::continueIf::clone "

Description:

    deep copy of a  continueIf terminator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::continueIf::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  continueIf terminator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::continueIf::apply "

Description:

    apply this operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::continueIf::~continueIf "

Usage:

    x.~continueIf()

"; 

%feature("docstring") simuPOP::controlledMating "

Applicability: diploid only

Description:

    a controlled mating scheme

Details:

    This is an experimental mating scheme that uses a frequency range
    to control the allele frequency of the offspring generation at
    given loci. When allele frequencies at the offspring generation
    does not fall into the given range, the offspring generation is
    regenerated. Any mating scheme can be used with this mating scheme
    by passing through parameter matingScheme.

"; 

%feature("docstring") simuPOP::controlledMating::controlledMating "

Description:

    control allele frequencies at a locus

Usage:

    controlledMating(matingScheme, loci, alleles, freqFunc,
      range=0.01)

Arguments:

    matingScheme:   a mating scheme
    loci:           loci at which allele frequency is controlled. Note
                    that controlling the allele frequencies at several
                    loci may take a long time.
    alleles:        alleles to control at each locus. Should have the
                    same length as loci.
    freqFunc:       frequency boundaries. If the length of the
                    returned value equals the size of loci, the range
                    for loci will be [value0, value0+range], [value1,
                    value1+range] etc. If the length of the returned
                    value is 2 times the size of loci, it will be
                    interpreted as [low1, high1, low2, high2, ...].

"; 

%ignore simuPOP::controlledMating::controlledMating(const controlledMating &rhs);

%feature("docstring") simuPOP::controlledMating::~controlledMating "

Description:

    destructor

Usage:

    x.~controlledMating()

"; 

%ignore simuPOP::controlledMating::submitScratch(population &pop, population &scratch);

%feature("docstring") simuPOP::controlledMating::clone "

Description:

    deep copy of a controlled mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::controlledMating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::controlledMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the controlled mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::controlledMating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%feature("docstring") simuPOP::controlledRandomMating "

Applicability: diploid only

Description:

    a controlled random mating scheme

Details:

    This is the controlled random mating scheme described in  Peng
    2007 (PLoS Genetics) . Basically, a freqFunc is passed to this
    mating scheme and set the allele frequencies of given alleles at
    given loci at the offspring generation.
    The offspring generation is conceptually populated in two steps.
    At the first step, only families with disease alleles are accepted
    until the expected number of disease alleles are met. At the
    second step, only families with wide type alleles are accepted to
    populate the rest of the offspring generation.

"; 

%feature("docstring") simuPOP::controlledRandomMating::controlledRandomMating "

Description:

    create a controlled random mating scheme

Usage:

    controlledRandomMating(loci, alleles, freqFunc, acceptScheme=0,
      numOffspring=1., sexParam=0.5, sexMode=MATE_RandomSex,
      numOffspringFunc=None, maxNumOffspring=0,
      mode=MATE_NumOffspring, newSubPopSize=[],
      newSubPopSizeFunc=None, newSubPopSizeExpr=\"\",
      contWhenUniSex=True, subPop=[], weight=0)

Details:

    Please refer to class mating for descriptions of other parameters.

Arguments:

    loci:           loci at which allele frequencies are monitored
                    (controlled)
    alleles:        alleles at given loci. It should have the same
                    length as loci
    freqFunc:       a Python function that accepts a generation number
                    and returns expected allele frequencies at given
                    loci
    acceptScheme:   internal use only

"; 

%ignore simuPOP::controlledRandomMating::controlledRandomMating(const controlledRandomMating &rhs);

%feature("docstring") simuPOP::controlledRandomMating::~controlledRandomMating "

Description:

    destructor

Usage:

    x.~controlledRandomMating()

"; 

%feature("docstring") simuPOP::controlledRandomMating::clone "

Description:

    deep copy of a controlled random mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::controlledRandomMating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::controlledRandomMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the controlled random mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::controlledRandomMating::submitScratch(population &pop, population &scratch);

%ignore simuPOP::controlledRandomMating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%feature("docstring") simuPOP::dumper "

Description:

    dump the content of a population.

"; 

%feature("docstring") simuPOP::dumper::dumper "

Description:

    dump a population

Usage:

    dumper(alleleOnly=False, infoOnly=False, ancestralPops=False,
      dispWidth=1, max=100, chrom=[], loci=[], subPop=[], indRange=[],
      output=\">\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, infoFields=[])

Arguments:

    alleleOnly:     only display allele
    infoOnly:       only display genotypic information
    dispWidth:      number of characters to display an allele. Default
                    to 1.
    ancestralPops:  whether or not display ancestral populations.
                    Default to False.
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

Description:

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::dumper::alleleOnly "

Description:

    only show alleles (not structure, gene information?

Usage:

    x.alleleOnly()

"; 

%feature("docstring") simuPOP::dumper::setAlleleOnly "

Usage:

    x.setAlleleOnly(alleleOnly)

"; 

%feature("docstring") simuPOP::dumper::infoOnly "

Description:

    only show info

Usage:

    x.infoOnly()

"; 

%feature("docstring") simuPOP::dumper::setInfoOnly "

Usage:

    x.setInfoOnly(infoOnly)

"; 

%feature("docstring") simuPOP::dumper::apply "

Description:

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

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

%feature("docstring") simuPOP::Expression::setLocalDict "

Description:

    set local dictionary

Usage:

    x.setLocalDict(dict)

"; 

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

%ignore simuPOP::GenoStructure::locusPos(UINT locus) const ;

%ignore simuPOP::GenoStructure::chromIndex(UINT ch) const ;

%ignore simuPOP::GenoStructure::setChromTypes(const vectoru &chromTypes);

%feature("docstring") simuPOP::GenoStruTrait "

Details:

    All individuals in a population share the same genotypic
    properties such as number of chromosomes, number and position of
    loci, names of markers, chromosomes, and information fields. These
    properties are stored in this  GenoStruTrait class and are
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
    have different names,  simuPOP uses the same names for alleles
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

    A  GenoStruTrait object is created with the creation of a
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

%ignore simuPOP::GenoStruTrait::gsAddChromFromStru(size_t idx) const ;

%ignore simuPOP::GenoStruTrait::gsAddLociFromStru(size_t idx) const ;

%ignore simuPOP::GenoStruTrait::gsRemoveLoci(const vectoru &loci, vectoru &kept);

%ignore simuPOP::GenoStruTrait::gsAddChrom(const vectorf &lociPos, const vectorstr &lociNames, const string &chromName, UINT chromType) const ;

%ignore simuPOP::GenoStruTrait::gsAddLoci(const vectoru &chrom, const vectorf &pos, const vectorstr &names, vectoru &newIndex) const ;

%ignore simuPOP::GenoStruTrait::genoStru() const ;

%ignore simuPOP::GenoStruTrait::genoStruIdx() const ;

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

%ignore simuPOP::GenoStruTrait::chromX() const ;

%ignore simuPOP::GenoStruTrait::chromY() const ;

%ignore simuPOP::GenoStruTrait::mitochondrial() const ;

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
    parameter of the population function. An  IndexError will be
    raised if the absolute index loc is greater than or equal to the
    total number of loci.

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

%ignore simuPOP::GenoStruTrait::chromIndex() const ;

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

    return the type of a chromosome chrom (1 for Autosome, 2 for
    ChromosomeX, 3 for ChromosomeY, and 4 for Mitochondrial).

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromTypes "

Usage:

    x.chromTypes()

Details:

    return the type of all chromosomes (1 for Autosome, 2 for
    ChromosomeX, 3 for ChromosomeY, and 4 for Mitochondrial).

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
    not specified, its index ('0', '1', '2', etc) is returned.

"; 

%feature("docstring") simuPOP::GenoStruTrait::alleleNames "

Usage:

    x.alleleNames()

Details:

    return a list of allele names given by the alleleNames parameter
    of the population function. This list does not have to cover all
    possible allele states of a population so  alleleNames()[allele]
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

    return the index of a locus with name name. Raise a  ValueError if
    no locus is found.

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociByNames "

Usage:

    x.lociByNames(names)

Details:

    return the indexes of loci with names names. Raise a  ValueError
    if any of the loci cannot be found.

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

    return the index of information field name. Raise an  IndexError
    if name is not one of the information fields.

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
    1]. In the multiple loci, alleles should be arranged by
    haplotypes, for example, loci=[0, 1], alleles=[0, 0, 1, 1] defines
    a VSP with individuals having genotype -0-0-, -1-1- or -1-1-,
    -0-0- at loci 0 and 1. If haplotypes -0-1- should be allowed, it
    should be added to explicitly, such as using alleles=[0, 0, 1, 1,
    0, 1, 0, 1].

"; 

%feature("docstring") simuPOP::genotypeSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::genotypeSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const ;

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
    geometric distribution with parameter  $ p $, which has mean  $
    \\frac{p}{1-p} $ and variance  $ \\frac{p}{\\left(1-p\\right)^{2}} $.
    gsmMutator implements both models. If you specify a Python
    function without a parameter, this mutator will use its return
    value each time a mutation occur; otherwise, a parameter  $ p $
    should be provided and the mutator will act as a geometric
    generalized stepwise model.

"; 

%feature("docstring") simuPOP::gsmMutator::gsmMutator "

Description:

    create a  gsmMutator

Usage:

    gsmMutator(rate=[], loci=[], maxAllele=0, incProb=0.5, p=0,
      func=None, output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
      end=-1, step=1, at=[], rep=REP_ALL, infoFields=[])

Details:

    The GSM model is developed for allozymes. It provides better
    description for these kinds of evolutionary processes.  Please see
    class mutator for the descriptions of other parameters.

Arguments:

    incProb:        probability to increase allele state. Default to
                    0.5.
    func:           a function that returns the number of steps. This
                    function does not accept any parameter.

Example:

Testsrc_gsmMutator.log Operator <tt> gsmMutator</tt>

"; 

%feature("docstring") simuPOP::gsmMutator::~gsmMutator "

Usage:

    x.~gsmMutator()

"; 

%feature("docstring") simuPOP::gsmMutator::clone "

Description:

    deep copy of a  gsmMutator

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
    of the  gsmMutator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::haplodiploidMating "

Applicability: haplodiploid only

Description:

    haplodiploid mating scheme of many hymemopterans

Details:

    This mating scheme is composed of an alphaParentChooser and a
    haplodiploidOffspringGenerator. The alphaParentChooser chooses a
    single Female randomly or from a given information field. This
    female will mate with random males from the colony. The offspring
    will have one of the two copies of chromosomes from the female
    parent, and the first copy of chromosomes from the male parent.
    Note that if a recombinator is used, it should disable
    recombination of male parent.

"; 

%feature("docstring") simuPOP::haplodiploidMating::haplodiploidMating "

Usage:

    haplodiploidMating(alphaSex=Female, alphaNum=1,
      alphaField=string, numOffspring=1., numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, sexParam=0.5,
      sexMode=MATE_RandomSex, newSubPopSize=[],
      newSubPopSizeFunc=None, newSubPopSizeExpr=\"\", subPop=[],
      weight=0)

Details:

    Please refer to class mating for descriptions of other parameters.

Arguments:

    alphaSex:       sex of the alpha individual. Default to Female.
    alphaNum:       Number of alpha individual. Default to one.
    alphaField:     information field that identifies the queen of the
                    colony. By default, a random female will be
                    chosen.

"; 

%feature("docstring") simuPOP::haplodiploidMating::~haplodiploidMating "

Description:

    destructor

Usage:

    x.~haplodiploidMating()

"; 

%feature("docstring") simuPOP::haplodiploidMating::clone "

Description:

    deep copy of a random mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::haplodiploidMating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::haplodiploidMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the random mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::haplodiploidMating::mateSubPop(population &pop, SubPopID subPop, RawIndIterator offBegin, RawIndIterator offEnd, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::haplodiploidOffspringGenerator "

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

%feature("docstring") simuPOP::haplodiploidOffspringGenerator::haplodiploidOffspringGenerator "

Usage:

    haplodiploidOffspringGenerator(numOffspring=1,
      numOffspringFunc=None, maxNumOffspring=1,
      mode=MATE_NumOffspring, sexParam=0.5, sexMode=MATE_RandomSex)

"; 

%feature("docstring") simuPOP::haplodiploidOffspringGenerator::copyParentalGenotype "

Usage:

    x.copyParentalGenotype(parent, it, ploidy, count)

"; 

%feature("docstring") simuPOP::haplodiploidOffspringGenerator::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::haplodiploidOffspringGenerator::generateOffspring(population &pop, individual *dad, individual *mom, RawIndIterator &offBegin, RawIndIterator &offEnd, vector< baseOperator * > &ops);

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
    can mimic the at parameter of an operator by  ifElse('rep in
    [2,5,9]' operator). The real use of this machanism is to monitor
    the population statistics and act accordingly.

"; 

%feature("docstring") simuPOP::ifElse::ifElse "

Description:

    create a conditional operator

Usage:

    ifElse(cond, ifOp=None, elseOp=None, output=\">\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

Arguments:

    cond:           expression that will be treated as a boolean
                    variable
    ifOp:           an operator that will be applied when cond is True
    elseOp:         an operator that will be applied when cond is
                    False

Example:

Testsrc_ifElse.log Operator <tt> ifElse</tt>

"; 

%feature("docstring") simuPOP::ifElse::~ifElse "

Description:

    destructor

Usage:

    x.~ifElse()

"; 

%ignore simuPOP::ifElse::ifElse(const ifElse &rhs);

%feature("docstring") simuPOP::ifElse::clone "

Description:

    deep copy of an  ifElse operator

Usage:

    x.clone()

"; 

%ignore simuPOP::ifElse::applyWithScratch(population &pop, population &scratch, int stage);

%ignore simuPOP::ifElse::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::ifElse::apply "

Description:

    apply the  ifElse operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::ifElse::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  ifElse operator

Usage:

    x.__repr__()

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
    class genoStruTrait), an individual class provides member
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
    individual are separated by  individual::totNumLoci() loci. It is
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

%ignore simuPOP::individual::genoPtr() const ;

%ignore simuPOP::individual::infoPtr() const ;

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

    return an editable array (a carray of length  totNumLoci()) that
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
    than  totNumLoci().

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

%feature("docstring") simuPOP::individual::unaffected "Obsolete or undocumented function."

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

%ignore simuPOP::individual::iteratable() const ;

%ignore simuPOP::individual::setIteratable(bool iteratable);

%ignore simuPOP::individual::visible() const ;

%ignore simuPOP::individual::setVisible(bool visible);

%feature("docstring") simuPOP::individual::subPopID "Obsolete or undocumented function."

%feature("docstring") simuPOP::individual::setSubPopID "Obsolete or undocumented function."

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

%ignore simuPOP::individual::genoBegin() const ;

%ignore simuPOP::individual::genoEnd() const ;

%ignore simuPOP::individual::genoBegin(UINT p) const ;

%ignore simuPOP::individual::genoEnd(UINT p) const ;

%ignore simuPOP::individual::genoBegin(UINT p, UINT chrom) const ;

%ignore simuPOP::individual::infoBegin() const ;

%ignore simuPOP::individual::infoEnd() const ;

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

    Unlike operator  pyEval and  pyExec that work at the population
    level, in its local namespace,  infoEval works at the individual
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
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

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

    deep copy of a  infoEval operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::infoEval::apply "

Description:

    apply the  infoEval operator

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::infoEval::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::infoEval::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  infoEval operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::infoEval::name "

Description:

    return the name of an expression

Usage:

    x.name()

Details:

    The name of a  infoEval operator is given by an optional parameter
    name. It can be used to identify this  infoEval operator in debug
    output, or in the dryrun mode of  simulator::evolve.

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
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

Details:

    Please refer to class  infoEval for parameter descriptions.

"; 

%feature("docstring") simuPOP::infoExec::~infoExec "

Usage:

    x.~infoExec()

"; 

%feature("docstring") simuPOP::infoExec::clone "

Description:

    deep copy of a  infoExec operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::infoExec::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  infoExec operator

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
    chooser (currently) uses  randomParentChooser to choose one parent
    and randomly choose another one from the information fields.
    Because of potentially non-even distribution of valid information
    fields, the overall process may not be as random as expected,
    especially when selection is applied. Note: if there is no valid
    individual, this parents chooser works like a double
    parentChooser.

"; 

%feature("docstring") simuPOP::infoParentsChooser::infoParentsChooser "

Usage:

    infoParentsChooser(infoFields=[], replacement=True,
      replenish=False)

Arguments:

    infoFields:     information fields that store index of matable
                    individuals.
    replacement:    if replacement is false, a parent can not be
                    chosen more than once.
    replenish:      if all parent has been chosen, choose from the
                    whole parental population again.

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

%ignore simuPOP::infoSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const ;

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

%feature("docstring") simuPOP::infoTagger "

Description:

    Tagging information fields.

Details:

    This is a simple post-mating tagger that write given information
    fields to a file (or standard output).

"; 

%feature("docstring") simuPOP::infoTagger::infoTagger "

Usage:

    infoTagger(begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      stage=PostMating, output=\">\", outputExpr=\"\", infoFields=[])

"; 

%feature("docstring") simuPOP::infoTagger::apply "

Description:

    add a newline

Usage:

    x.apply(pop)

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

    create an  inheritTagger that inherits a tag from one or both
    parents

Usage:

    inheritTagger(mode=TAG_Paternal, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, output=\"\", outputExpr=\"\",
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
    of the  inheritTagger

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::inheritTagger::applyDuringMating "

Description:

    apply the  inheritTagger

Usage:

    x.applyDuringMating(pop, offspring, dad=None, mom=None)

"; 

%feature("docstring") simuPOP::inheritTagger::clone "

Description:

    deep copy of a  inheritTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initByFreq "

Function form:

    InitByFreq

Description:

    initialize genotypes by given allele frequencies, and sex by male
    frequency

Details:

    This operator assigns alleles at loci with given allele
    frequencies. By default, all individuals will be assigned with
    random alleles. If identicalInds=True, an individual is assigned
    with random alleles and is then copied to all others. If subPop or
    indRange is given, multiple arrays of alleleFreq can be given to
    given different frequencies for different subpopulation or
    individual ranges.

"; 

%feature("docstring") simuPOP::initByFreq::initByFreq "

Description:

    randomly assign alleles according to given allele frequencies

Usage:

    initByFreq(alleleFreq=[], identicalInds=False, subPop=[],
      indRange=[], loci=[], atPloidy=-1, maleFreq=0.5, sex=[],
      stage=PreMating, begin=0, end=1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

Arguments:

    alleleFreq:     an array of allele frequencies. The sum of all
                    frequencies must be 1; or for a matrix of allele
                    frequencies, each row corresponses to a
                    subpopulation or range.
    identicalInds:  whether or not make individual genotypes identical
                    in all subpopulations. If True, this operator will
                    randomly generate genotype for an individual and
                    spread it to the whole subpopulation in the given
                    range.
    sex:            an array of sex [Male, Female, Male...] for
                    individuals. The length of sex will not be
                    checked. If it is shorter than the number of
                    individuals, sex will be reused from the
                    beginning.
    stage:          default to PreMating.

Example:

Testsrc_initByFreq.log Operator <tt> initByFreq</tt>

"; 

%feature("docstring") simuPOP::initByFreq::~initByFreq "

Usage:

    x.~initByFreq()

"; 

%feature("docstring") simuPOP::initByFreq::clone "

Description:

    deep copy of the operator  initByFreq

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initByFreq::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator  initByFreq

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

Description:

    initialize genotype by value and then copy to all individuals

Details:

    Operator  initByValue gets one copy of chromosomes or the whole
    genotype (or of those corresponds to loci) of an individual and
    copy them to all or a subset of individuals. This operator assigns
    given alleles to specified individuals. Every individual will have
    the same genotype. The parameter combinations should be
    *  value - subPop/indRange: individual in subPop or in range(s)
    will be assigned genotype value;
    *  subPop/indRange: subPop or indRange should have the same length
    as value. Each item of value will be assigned to each subPop or
    indRange.

"; 

%feature("docstring") simuPOP::initByValue::initByValue "

Description:

    initialize a population by given alleles

Usage:

    initByValue(value=[], loci=[], atPloidy=-1, subPop=[],
      indRange=[], proportions=[], maleFreq=0.5, sex=[],
      stage=PreMating, begin=0, end=1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

Arguments:

    value:          an array of genotypes of one individual, having
                    the same length as the length of loci() or
                    loci()*ploidy() or pop.genoSize() (whole genotype)
                    or totNumLoci() (one copy of chromosomes). This
                    parameter can also be an array of arrays of
                    genotypes of one individual. If value is an array
                    of values, it should have the length one, number
                    of subpopulations, or the length of ranges of
                    proportions.
    proportions:    an array of percentages for each item in value. If
                    given, assign given genotypes randomly.
    maleFreq:       male frequency
    sex:            an array of sex [Male, Female, Male...] for
                    individuals. The length of sex will not be
                    checked. If length of sex is shorter than the
                    number of individuals, sex will be reused from the
                    beginning.
    stages:         default to PreMating.

Example:

Testsrc_initByValue.log Operator <tt> initByValue</tt>

"; 

%feature("docstring") simuPOP::initByValue::~initByValue "

Usage:

    x.~initByValue()

"; 

%feature("docstring") simuPOP::initByValue::clone "

Description:

    deep copy of the operator  initByValue

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initByValue::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator  initByValue

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::initByValue::apply "

Description:

    apply this operator to population pop

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::initializer "

Description:

    initialize alleles at the start of a generation

Details:

    Initializers are used to initialize populations before evolution.
    They are set to be PreMating operators by default.  simuPOP
    provides three initializers. One assigns alleles by random, one
    assigns a fixed set of genotypes, and the last one calls a user-
    defined function.

"; 

%feature("docstring") simuPOP::initializer::initializer "

Description:

    create an initializer. Default to be always active.

Usage:

    initializer(subPop=[], indRange=[], loci=[], atPloidy=-1,
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

Arguments:

    subPop:         an array specifies applicable subpopulations
    indRange:       a [begin, end] pair of the range of absolute
                    indexes of individuals, for example, ([1,2]); or
                    an array of [begin, end] pairs, such as
                    ([[1,4],[5,6]]). This is how you can initialize
                    individuals differently within subpopulations.
                    Note that ranges are in the form of [a,b). I.e.,
                    range [4,6] will intialize individual 4, 5, but
                    not 6. As a shortcut for [4,5], you can use [4] to
                    specify one individual.
    loci:           a vector of locus indexes at which initialization
                    will be done. If empty, apply to all loci.
    locus:          a shortcut to loci
    atPloidy:       initialize which copy of chromosomes. Default to
                    all.

"; 

%feature("docstring") simuPOP::initializer::~initializer "

Description:

    destructor

Usage:

    x.~initializer()

"; 

%feature("docstring") simuPOP::initializer::clone "

Description:

    deep copy of an initializer

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initializer::__repr__ "

Description:

    used by Python print function to print out the general information
    of the initializer

Usage:

    x.__repr__()

"; 

%ignore simuPOP::initializer::setRanges(population &pop);

%feature("docstring") simuPOP::initSex "

Function form:

    InitSex

Details:

    An operator to initialize individual sex. For convenience, this
    operator is included by other initializers such as  initByFreq,
    initByValue, or  pyInit.

"; 

%feature("docstring") simuPOP::initSex::initSex "

Description:

    initialize individual sex.

Usage:

    initSex(maleFreq=0.5, sex=[], subPop=[], indRange=[], loci=[],
      atPloidy=-1, stage=PreMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, infoFields=[])

Arguments:

    maleFreq:       male frequency. Default to 0.5. Sex will be
                    initialized with this parameter.
    sex:            a list of sexes (Male or Female) and will be
                    applied to individuals in in turn. If specified,
                    parameter maleFreq is ignored.

"; 

%feature("docstring") simuPOP::initSex::~initSex "

Description:

    destructor

Usage:

    x.~initSex()

"; 

%feature("docstring") simuPOP::initSex::clone "

Description:

    deep copy of an  initSex

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initSex::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  initSex

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
    state is actually  $ \\frac{rate}{K-1} $, where  $ K $ is specified
    by parameter maxAllele.

"; 

%feature("docstring") simuPOP::kamMutator::kamMutator "

Description:

    create a K-Allele Model mutator

Usage:

    kamMutator(rate=[], loci=[], maxAllele=0, output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, infoFields=[])

Details:

    Please see class mutator for the descriptions of other parameters.

Arguments:

    rate:           mutation rate. It is the 'probability to mutate'.
                    The actual mutation rate to any of the other K-1
                    allelic states are rate/(K-1).
    maxAllele:      maximum allele that can be mutated to. For binary
                    libraries, allelic states will be [0, maxAllele].
                    Otherwise, they are [1, maxAllele].

Example:

Testsrc_kamMutator.log Operator <tt> kamMutator</tt>

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

    deep copy of a  kamMutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::kamMutator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  kamMutator

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
    considered as diseased alleles.  maPenetrance accepts an array of
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
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
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

Example:

Testsrc_maPenetrance.log Operator <tt> maPenetrance</tt>

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

%feature("docstring") simuPOP::maPenetrance::penet "

Description:

    currently assuming diploid

Usage:

    x.penet(ind)

"; 

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
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
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

Example:

Testsrc_mapPenetrance.log Operator <tt> mapPenetrance</tt>

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
    return value is the sum of the trait plus  $
    N\\left(0,\\sigma^{2}\\right) $. This random part is usually
    considered as the environmental factor of the trait.

"; 

%feature("docstring") simuPOP::mapQuanTrait::mapQuanTrait "

Description:

    create a map quantitative trait operator

Usage:

    mapQuanTrait(loci, qtrait, sigma=0, phase=False,
      ancestralGen=-1, stage=PostMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, infoFields=[\"qtrait\"])

Arguments:

    locus:          the locus index. The quantitative trait is
                    determined by genotype at this locus.
    loci:           an array of locus indexes. The quantitative trait
                    is determined by genotypes at these loci.
    qtrait:         a dictionary of quantitative traits. The genotype
                    must be in the form of 'a-b'. This is the mean of
                    the quantitative trait. The actual trait value
                    will be  $ N\\left(mean,\\sigma^{2}\\right) $. For
                    multiple loci, the form is 'a-b|c-d|e-f' etc.
    sigma:          standard deviation of the environmental factor  $
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

%feature("docstring") simuPOP::mapQuanTrait::qtrait "

Description:

    currently assuming diploid

Usage:

    x.qtrait(ind)

"; 

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
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[\"fitness\"])

Arguments:

    locus:          the locus index. A shortcut to  loci=[locus]
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

Example:

Testsrc_mapSelector.log Operator <tt> mapSelector</tt>

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

%feature("docstring") simuPOP::mapSelector::indFitness "

Description:

    calculate/return the fitness value, currently assuming diploid

Usage:

    x.indFitness(ind, gen)

"; 

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
    are considered as diseased alleles.  maQuanTrait accepts an array
    of fitness. Quantitative trait is then set for any given genotype.
    A standard normal distribution  $ N\\left(0,\\sigma^{2}\\right) $
    will be added to the returned trait value.

"; 

%feature("docstring") simuPOP::maQuanTrait::maQuanTrait "

Description:

    create a multiple allele quantitative trait operator

Usage:

    maQuanTrait(loci, qtrait, wildtype, sigma=[], ancestralGen=-1,
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[\"qtrait\"])

Details:

    Please refer to  quanTrait for other parameter descriptions.

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

%feature("docstring") simuPOP::maQuanTrait::qtrait "

Description:

    currently assuming diploid

Usage:

    x.qtrait(ind)

"; 

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
    *  For a model with more than two loci, use a table of length  $
    3^{n} $ in a order similar to the two-locus model.

"; 

%feature("docstring") simuPOP::maSelector::maSelector "

Description:

    create a multiple allele selector

Usage:

    maSelector(loci, fitness, wildtype, subPops=[], stage=PreMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[\"fitness\"])

Details:

    Please refer to  baseOperator for other parameter descriptions.

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

    *  maSelector only works for diploid populations.
    * wildtype alleles at all loci are the same.

Example:

Testsrc_maSelector.log Operator <tt> maSelector</tt>

"; 

%feature("docstring") simuPOP::maSelector::~maSelector "

Usage:

    x.~maSelector()

"; 

%feature("docstring") simuPOP::maSelector::clone "

Description:

    deep copy of a  maSelector

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
    of the  maSelector

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

%ignore simuPOP::mating::isCompatible(const population &pop) const ;

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
                    subpopulation. This is only used in  heteroMating
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

Example:

Testsrc_mating.log Demographic models and control of number of offspring per mating event 

"; 

%ignore simuPOP::mating::mating(const mating &rhs);

%feature("docstring") simuPOP::mating::~mating "

Description:

    destructor

Usage:

    x.~mating()

"; 

%ignore simuPOP::mating::subPop() const ;

%ignore simuPOP::mating::virtualSubPop() const ;

%ignore simuPOP::mating::weight() const ;

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

%feature("docstring") simuPOP::mating::submitScratch "

Description:

    a common submit procedure is defined.

Usage:

    x.submitScratch(pop, scratch)

"; 

%ignore simuPOP::mating::mateSubPop(population &pop, SubPopID subPop, RawIndIterator offBegin, RawIndIterator offEnd, vector< baseOperator * > &ops);

%ignore simuPOP::mating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%ignore simuPOP::mating::prepareScratchPop(population &pop, population &scratch);

%feature("docstring") simuPOP::mendelianOffspringGenerator "

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

%feature("docstring") simuPOP::mendelianOffspringGenerator::mendelianOffspringGenerator "

Usage:

    mendelianOffspringGenerator(numOffspring=1,
      numOffspringFunc=None, maxNumOffspring=1,
      mode=MATE_NumOffspring, sexParam=0.5, sexMode=MATE_RandomSex)

"; 

%feature("docstring") simuPOP::mendelianOffspringGenerator::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::mendelianOffspringGenerator::initialize(const population &pop, vector< baseOperator * > const &ops);

%feature("docstring") simuPOP::mendelianOffspringGenerator::formOffspringGenotype "

Description:

    does not set sex if count == -1.

Usage:

    x.formOffspringGenotype(parent, it, ploidy, count)

Arguments:

    count:          index of offspring, used to set offspring sex

"; 

%ignore simuPOP::mendelianOffspringGenerator::generateOffspring(population &pop, individual *dad, individual *mom, RawIndIterator &offBegin, RawIndIterator &offEnd, vector< baseOperator * > &ops);

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

    mergeSubPops(subPops=[], removeEmptySubPops=False,
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

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

    deep copy of a  mergeSubPops operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mergeSubPops::apply "

Description:

    apply a  mergeSubPops operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::mergeSubPops::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  mergeSubPops operator

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
    simuPOP. Migrators are quite flexible in  simuPOP in the sense
    that
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
      rep=REP_ALL, infoFields=[])

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
    *  PEN_Multiplicative: the penetrance is calculated as  $ f=\\prod
    f_{i} $.
    *  PEN_Additive: the penetrance is calculated as  $
    f=\\min\\left(1,\\sum f_{i}\\right) $.  $ f $ will be set to 1 when  $
    f<0 $. In this case,  $ s_{i} $ are added, not  $ f_{i} $
    directly.
    *  PEN_Heterogeneity: the penetrance is calculated as  $
    f=1-\\prod\\left(1-f_{i}\\right) $. Please refer to Neil Risch (1990)
    for detailed information about these models.

"; 

%feature("docstring") simuPOP::mlPenetrance::mlPenetrance "

Description:

    create a multiple locus penetrance operator

Usage:

    mlPenetrance(peneOps, mode=PEN_Multiplicative, ancestralGen=-1,
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

Arguments:

    peneOps:        a list of penetrance operators
    mode:           can be one of PEN_Multiplicative, PEN_Additive,
                    and PEN_Heterogeneity

Example:

Testsrc_mlPenetrance.log Operator <tt> mlPenetrance</tt>

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

%feature("docstring") simuPOP::mlPenetrance::penet "

Description:

    currently assuming diploid

Usage:

    x.penet(ind)

"; 

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

    Operator  mlQuanTrait is a 'multiple-locus' quantitative trait
    calculator. It accepts a list of quantitative traits and combine
    them according to the mode parameter, which takes one of the
    following values
    *  QT_Multiplicative: the mean of the quantitative trait is
    calculated as  $ f=\\prod f_{i} $.
    *  QT_Additive: the mean of the quantitative trait is calculated
    as  $ f=\\sum f_{i} $. Note that all  $ \\sigma_{i} $ (for  $ f_{i}
    $) and  $ \\sigma $ (for  $ f $) will be considered. I.e, the trait
    value should be
    $ f=\\sum_{i}\\left(f_{i}+N\\left(0,\\sigma_{i}^{2}\\right)\\right)+\\sig
    ma^{2} $ for QT_Additive case. If this is not desired, you can set
    some of the  $ \\sigma $ to zero.

"; 

%feature("docstring") simuPOP::mlQuanTrait::mlQuanTrait "

Description:

    create a multiple locus quantitative trait operator

Usage:

    mlQuanTrait(qtraits, mode=QT_Multiplicative, sigma=0,
      ancestralGen=-1, stage=PostMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, infoFields=[\"qtrait\"])

Details:

    Please refer to  quanTrait for other parameter descriptions.

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
    takes a vector of selectors (can not be another  mlSelector) and
    evaluate the fitness of an individual as the product or sum of
    individual fitness values. The mode is determined by parameter
    mode, which takes one of the following values
    *  SEL_Multiplicative: the fitness is calculated as  $
    f=\\prod_{i}f_{i} $, where  $ f_{i} $ is the single-locus fitness
    value.
    *  SEL_Additive: the fitness is calculated as  $
    f=\\max\\left(0,1-\\sum_{i}(1-f_{i})\\right) $.  $ f $ will be set to
    0 when  $ f<0 $.

"; 

%feature("docstring") simuPOP::mlSelector::mlSelector "

Description:

    create a multiple-locus selector

Usage:

    mlSelector(selectors, mode=SEL_Multiplicative, subPops=[],
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[\"fitness\"])

Details:

    Please refer to  mapSelector for other parameter descriptions.

Arguments:

    selectors:      a list of selectors

Example:

Testsrc_mlSelector.log Operator <tt> mlSelector</tt>

"; 

%feature("docstring") simuPOP::mlSelector::~mlSelector "

Usage:

    x.~mlSelector()

"; 

%feature("docstring") simuPOP::mlSelector::clone "

Description:

    deep copy of a  mlSelector

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mlSelector::indFitness "

Description:

    calculate/return the fitness value, currently assuming diploid

Usage:

    x.indFitness(ind, gen)

"; 

%feature("docstring") simuPOP::mlSelector::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  mlSelector

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::monogamousMating "

Applicability: diploid only

Description:

    a mating scheme of monogamy

Details:

    This mating scheme is identical to random mating except that
    parents are chosen without replacement. Under this mating scheme,
    offspring share the same mother must share the same father. In
    case that all parental pairs are exhausted, parameter
    replenish=True allows for the replenishment of one or both sex
    groups.

"; 

%feature("docstring") simuPOP::monogamousMating::monogamousMating "

Usage:

    monogamousMating(replenish=False, numOffspring=1.,
      numOffspringFunc=None, maxNumOffspring=0,
      mode=MATE_NumOffspring, sexParam=0.5, sexMode=MATE_RandomSex,
      newSubPopSize=[], newSubPopSizeFunc=None, newSubPopSizeExpr=\"\",
      contWhenUniSex=True, subPop=[], weight=0)

Details:

    replenish This parameter allows replenishment of one or both
    parental sex groups in case that they are are exhausted. Default
    to False. Please refer to class mating for descriptions of other
    parameters.

"; 

%feature("docstring") simuPOP::monogamousMating::clone "

Description:

    deep copy of a random mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::monogamousMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the random mating scheme

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
      rep=REP_ALL, infoFields=[])

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

%feature("docstring") simuPOP::noMating "

Applicability: all ploidy

Description:

    a mating scheme that does nothing

Details:

    In this scheme, there is
    *  no mating. Parent generation will be considered as offspring
    generation.
    *  no subpopulation change. During-mating operators will be
    applied, but the return values are not checked. I.e.,
    subpopulation size parameters will be ignored although some
    during-mating operators might be applied. Note that because the
    offspring population is the same as parental population, this
    mating scheme can not be used with other mating schemes in a
    heterogeneous mating scheme.  cloneMating is recommended for that
    purpose.

"; 

%feature("docstring") simuPOP::noMating::noMating "

Description:

    creat a scheme with no mating

Usage:

    noMating(numOffspring=1.0, numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, sexParam=0.5,
      sexMode=MATE_RandomSex, newSubPopSize=[], newSubPopSizeExpr=\"\",
      newSubPopSizeFunc=None, subPop=[], weight=0)

Note:

    All parameters are ignored!

"; 

%feature("docstring") simuPOP::noMating::~noMating "

Description:

    destructor

Usage:

    x.~noMating()

"; 

%feature("docstring") simuPOP::noMating::clone "

Description:

    deep copy of a scheme with no mating

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::noMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the scheme with no mating

Usage:

    x.__repr__()

"; 

%ignore simuPOP::noMating::submitScratch(population &pop, population &scratch);

%ignore simuPOP::noMating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

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
      end=0, step=1, at=[], rep=REP_ALL, infoFields=[])

Example:

Testsrc_noneOp.log Operator <tt> noneOp</tt>

"; 

%feature("docstring") simuPOP::noneOp::~noneOp "

Description:

    destructor

Usage:

    x.~noneOp()

"; 

%feature("docstring") simuPOP::noneOp::clone "

Description:

    deep copy of a  noneOp operator

Usage:

    x.clone()

"; 

%ignore simuPOP::noneOp::applyWithScratch(population &pop, population &scratch, int stage);

%ignore simuPOP::noneOp::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::noneOp::apply "

Description:

    apply the  noneOp operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::noneOp::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  noneOp operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::offspringGenerator "

Applicability: all ploidy

Details:

    Offspring generators generate offspring from given parents.
    Generators differ from each other by how and how many offspring is
    generated at each mating event. Parameters mode, numOffspring,
    maxNumOffspring and numOffspringFunc are used to specify how many
    offspring will be produced at each mating event. mode can be one
    of
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
    *  MATE_UniformDistribution: a Uniform  [a, b]  distribution with
    parameter numOffspring (a) and maxNumOffspring (b) is used to
    determine the number of offspring of each family. This is the base
    class of all offspring generators, and should not be used
    directly.

"; 

%feature("docstring") simuPOP::offspringGenerator::offspringGenerator "

Usage:

    offspringGenerator(numOffspring, numOffspringFunc,
      maxNumOffspring, mode, sexParam, sexMode)

Arguments:

    numOffspring:   Depending on , this paramter can be the number of
                    offspring to produce, or a paremter of a random
                    distribution.
    numOffspringFunc:a Python function that returns the number of
                    offspring at each mating event. The setting of
                    this parameter implies =MATE_PyNumOffspring.
    maxNumOffspring:used when numOffspring is generated from a
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
                    Female.

Note:

    : Parameter sexMode and sexParam are ignored if sex chromosome is
    defined. Offspring sex is defined by genotype in this case.

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

%ignore simuPOP::offspringGenerator::generateOffspring(population &pop, individual *dad, individual *mom, RawIndIterator &offBegin, RawIndIterator &offEnd, vector< baseOperator * > &ops);

%ignore simuPOP::offspringGenerator::fixedFamilySize() const ;

%ignore simuPOP::offspringGenerator::initialize(const population &pop, vector< baseOperator * > const &ops);

%feature("docstring") simuPOP::offspringGenerator::finalize "

Usage:

    x.finalize(pop)

"; 

%ignore simuPOP::offspringGenerator::numOffspring(int gen);

%feature("docstring") simuPOP::offspringGenerator::getSex "

Usage:

    x.getSex(count)

Arguments:

    count:          the index of offspring

"; 

%ignore simuPOP::offspringGenerator::initialized() const ;

%ignore simuPOP::offspringGenerator::setNumParents(int numParents);

%ignore simuPOP::offspringGenerator::numParents() const ;

%ignore simuPOP::OstreamManager;

%ignore simuPOP::OstreamManager::OstreamManager();

%feature("docstring") simuPOP::OstreamManager::~OstreamManager "

Usage:

    x.~OstreamManager()

"; 

%feature("docstring") simuPOP::OstreamManager::getOstream "

Description:

    if the stream does not exist, create one and return.

Usage:

    x.getOstream(name, readable, realAppend, useString)

"; 

%feature("docstring") simuPOP::OstreamManager::hasOstream "

Description:

    this is mostly for debug purposes

Usage:

    x.hasOstream(filename)

"; 

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
      end=-1, step=1, at=[], rep=REP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::outputer::~outputer "

Description:

    destructor

Usage:

    x.~outputer()

"; 

%feature("docstring") simuPOP::outputer::clone "

Description:

    deep copy of an operator

Usage:

    x.clone()

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

%ignore simuPOP::parentChooser::initialized() const ;

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

    create a  parentsTagger

Usage:

    parentsTagger(begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      output=\"\", outputExpr=\"\", infoFields=[\"father_idx\",
      \"mother_idx\"])

"; 

%feature("docstring") simuPOP::parentsTagger::~parentsTagger "

Usage:

    x.~parentsTagger()

"; 

%feature("docstring") simuPOP::parentsTagger::clone "

Description:

    deep copy of a  parentsTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::parentsTagger::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  parentsTagger

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::parentsTagger::applyDuringMating "

Description:

    apply the  parentsTagger

Usage:

    x.applyDuringMating(pop, offspring, dad=None, mom=None)

"; 

%feature("docstring") simuPOP::parentsTagger::apply "

Description:

    with a newline.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::parentTagger "

Description:

    tagging according to parental indexes

Details:

    This during-mating operator set tag() each individual with indexes
    of his/her parent in the parental population. Because only one
    parent is recorded, this is recommended to be used for mating
    schemes that requires only one parent (such as  selfMating). This
    tagger record indexes to information field parent_idx, and/or a
    given file. The usage is similar to  parentsTagger.

"; 

%feature("docstring") simuPOP::parentTagger::parentTagger "

Description:

    create a  parentTagger

Usage:

    parentTagger(begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      output=\"\", outputExpr=\"\", infoFields=[\"parent_idx\"])

"; 

%feature("docstring") simuPOP::parentTagger::~parentTagger "

Usage:

    x.~parentTagger()

"; 

%feature("docstring") simuPOP::parentTagger::clone "

Description:

    deep copy of a  parentTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::parentTagger::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  parentTagger

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::parentTagger::applyDuringMating "

Description:

    apply the  parentTagger

Usage:

    x.applyDuringMating(pop, offspring, dad=None, mom=None)

"; 

%feature("docstring") simuPOP::parentTagger::apply "

Description:

    with a newline.

Usage:

    x.apply(pop)

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
      begin=0, end=-1, step=1, at=[], rep=REP_LAST, infoFields=[])

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

%feature("docstring") simuPOP::pause::clone "

Description:

    deep copy of a pause operator

Usage:

    x.clone()

"; 

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
    the information fields used by operator  parentsTagger. This
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

    For example,  setInfoWithRelatives(pathGen = [0, 1, 1, 0],
    pathFields = [['father_idx', 'mother_idx'], ['sib1', 'sib2'],
    ['off1', 'off2']], pathSex = [AnySex, MaleOnly, FemaleOnly],
    resultFields = ['cousin1', 'cousin2'])  This function will 1.
    locate father_idx and mother_idx for each individual at generation
    0 (pathGen[0]) 2. find AnySex individuals referred by father_idx
    and mother_idx at generation 1 (pathGen[1]) 3. find informaton
    fields sib1 and sib2 from these parents 4. locate MaleOnly
    individuals referred by sib1 and sib2 from generation 1
    (pathGen[2]) 5. find information fields off1 and off2 from these
    individuals, and 6. locate FemaleOnly indiviudals referred by off1
    and  from geneartion 0 (pathGen[3]) 7. Save index of these
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
    *  use e.g.,  mlPenetrance(...., infoFields=['penetrance']) to add
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
      step=1, at=[], rep=REP_ALL, infoFields=[])

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

%feature("docstring") simuPOP::penetrance::penet "

Description:

    calculate/return penetrance etc.

Usage:

    x.penet()

"; 

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
    diseases etc.  pointMutator, as its name suggest, does point
    mutation. This mutator will turn alleles at loci on the first
    chromosome copy to toAllele for individual inds. You can specify
    atPloidy to mutate other, or all ploidy copies.

"; 

%feature("docstring") simuPOP::pointMutator::pointMutator "

Description:

    create a  pointMutator

Usage:

    pointMutator(loci, toAllele, atPloidy=[], inds=[], output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, infoFields=[])

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

    deep copy of a  pointMutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pointMutator::apply "

Description:

    apply a  pointMutator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pointMutator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  pointMutator

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

%feature("docstring") simuPOP::polygamousMating "

Applicability: diploid only

Description:

    a mating scheme of polygymy or polyandry

Details:

    This mating scheme is composed of a random parents chooser that
    allows for polygamous mating, and a mendelian offspring generator.
    In this mating scheme, a male (or female) parent will have more
    than one sex partner (numPartner). Parents returned from this
    parents chooser will yield the same male (or female) parents, each
    with varying partners.

"; 

%feature("docstring") simuPOP::polygamousMating::polygamousMating "

Usage:

    polygamousMating(polySex=Male, polyNum=1, replacement=False,
      replenish=False, numOffspring=1., numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, sexParam=0.5,
      sexMode=MATE_RandomSex, newSubPopSize=[],
      newSubPopSizeFunc=None, newSubPopSizeExpr=\"\",
      contWhenUniSex=True, subPop=[], weight=0)

Arguments:

    polySex:        sex of polygamous mating. Male for polygyny,
                    Female for polyandry.
    polyNum:        Number of sex partners.
    replacement:    If set to True, a parent can be chosen to mate
                    again. Default to False.
    replenish:      In case that replacement=True, whether or not
                    replenish a sex group when it is exhausted. Please
                    refer to class mating for descriptions of other
                    parameters.

"; 

%feature("docstring") simuPOP::polygamousMating::clone "

Description:

    deep copy of a random mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::polygamousMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the random mating scheme

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::population "

Details:

    A  simuPOP population consists of individuals of the same
    genotypic structure, organized by generations, subpopulations and
    virtual subpopulations. It also contains a Python dictionary that
    is used to store arbitrary population variables. In addition to
    genotypic structured related functions provided by the
    genoStruTrait class, the population class provides a large number
    of member functions that can be used to
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
      lociPos=[], ancestralGens=0, chromNames=[], alleleNames=[],
      lociNames=[], infoFields=[])

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
                    sets, even if some chromosomes (e.g.
                    mitochondrial) or some individuals (e.g. males in
                    a haplodiploid population) have different numbers
                    of homologous sets. The first case is handled by
                    setting chromTypes of each chromosome. Only the
                    haplodiploid populations are handled for the
                    second case, for which ploidy=Haplodiploid should
                    be used.
    loci:           A list of numbers of loci on each chromosome. The
                    length of this parameter determines the number of
                    chromosomes. Default to [1], meaning one
                    chromosome with a single locus.
    chromTypes:     A list that specifies the type of each chromosome,
                    which can be Autosome, ChromosomeX, ChromosomeY,
                    or Mitochondrial. All chromosomes are assumed to
                    be autosomes if this parameter is ignored. Sex
                    chromosome can only be specified in a diploid
                    population where the sex of an individual is
                    determined by the existence of these chromosomes
                    using the XX (Female) and XY (Male) convention.
                    Both sex chromosomes have to be available and be
                    specified only once. Because chromosomes X and Y
                    are treated as two chromosomes, recombination on
                    the pseudo-autosomal regions of the sex chromsomes
                    is not supported. A Mitochondrial chromosome only
                    exists in females and is inherited maternally.
    lociPos:        Positions of all loci on all chromosome, as a list
                    of float numbers. Default to 1, 2, ... etc on each
                    chromosome. Positions on the same chromosome
                    should be ordered. A nested list that specifies
                    positions of loci on each chromosome is also
                    acceptable.
    ancestralGens:  Number of the most recent ancestral generations to
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
                    Note that  simuPOP does not yet support locus-
                    specific allele names.
    lociNames:      A list or a matrix (separated by chromosomes) of
                    names for each locus. Default to \"locX-Y\" where X
                    and Y are 1-based chromosome and locus indexes,
                    respectively.
    infoFields:     Names of information fields (named float number)
                    that will be attached to each individual.

"; 

%ignore simuPOP::population::population(const population &rhs);

%feature("docstring") simuPOP::population::clone "

Usage:

    x.clone(keepAncestralPops=-1)

Details:

    Copy a population, with the option to keep all (default), no, or a
    given number of ancestral generations (keepAncestralPops = -1, 0,
    or a positive number, respectively). Note that Python statement
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

%ignore simuPOP::population::validate(const string &msg) const ;

%ignore simuPOP::population::fitSubPopStru(const vectorlu &newSubPopSizes);

%ignore simuPOP::population::hasActivatedVirtualSubPop() const ;

%ignore simuPOP::population::hasActivatedVirtualSubPop(SubPopID subPop) const ;

%ignore simuPOP::population::hasVirtualSubPop() const ;

%ignore simuPOP::population::virtualSplitter() const ;

%feature("docstring") simuPOP::population::setVirtualSplitter "

Usage:

    x.setVirtualSplitter(splitter)

Details:

    Set a VSP splitter to the population, which defines the same VSPs
    for all subpopulations. If different VSPs are needed for different
    subpopulations, a  combinedSplitter can be used to make these VSPs
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

%feature("docstring") simuPOP::population::virtualSubPopName "

Usage:

    x.virtualSubPopName(subPop)

Details:

    Return the name of a virtual subpopulation subPop (specified by a
    (sp, vsp) pair). Because VSP names are the same across all
    subpopulations, a single VSP index is also acceptable.

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
    An  IndexError will be raised if subPop is out of range.

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

    x.ind(idx, subPop=0)

Details:

    Return a refernce to individual ind in subpopulation subPop.

"; 

%ignore simuPOP::population::ind(ULONG idx, UINT subPop=0) const ;

%feature("docstring") simuPOP::population::ancestor "

Usage:

    x.ancestor(idx, gen)

Details:

    Return a reference to individual idx in ancestral generation gen.
    The correct individual will be returned even if the current
    generation is not the present one (see useAncestralGen).

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
    subpopulation (subPop=[spID,  vspID]).

"; 

%ignore simuPOP::population::indOrdered();

%ignore simuPOP::population::setIndOrdered(bool s);

%ignore simuPOP::population::indBegin(IterationType type=VisibleInds);

%ignore simuPOP::population::indEnd(IterationType type=VisibleInds);

%ignore simuPOP::population::indBegin(UINT subPop, IterationType type=VisibleInds);

%ignore simuPOP::population::indEnd(UINT subPop, IterationType type=VisibleInds);

%ignore simuPOP::population::indBegin(IterationType type=VisibleInds) const ;

%ignore simuPOP::population::indEnd(IterationType type=VisibleInds) const ;

%ignore simuPOP::population::indEnd(UINT subPop, IterationType type) const ;

%ignore simuPOP::population::rawIndBegin();

%ignore simuPOP::population::rawIndEnd();

%ignore simuPOP::population::rawIndBegin(UINT subPop);

%ignore simuPOP::population::rawIndEnd(UINT subPop);

%ignore simuPOP::population::rawIndBegin() const ;

%ignore simuPOP::population::rawIndEnd() const ;

%ignore simuPOP::population::rawIndEnd(UINT subPop) const ;

%ignore simuPOP::population::alleleBegin(UINT locus);

%ignore simuPOP::population::alleleEnd(UINT locus);

%ignore simuPOP::population::alleleBegin(UINT locus, UINT subPop);

%ignore simuPOP::population::genoBegin(bool order);

%ignore simuPOP::population::genoEnd(bool order);

%ignore simuPOP::population::genoBegin(UINT subPop, bool order);

%ignore simuPOP::population::genoEnd(UINT subPop, bool order);

%ignore simuPOP::population::indGenoBegin(ULONG ind) const ;

%ignore simuPOP::population::indGenoEnd(ULONG ind) const ;

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
    subpopulation subPop.

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

    Fill the genotype of all individuals of in subpopulation subPop
    using a list of alleles geno. geno will be reused if its length is
    less than subPopSize(subPop)*totNumLoci()*ploidy().

"; 

%feature("docstring") simuPOP::population::setIndSubPopID "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::setIndSubPopIDWithID "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::setSubPopByIndID "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::splitSubPop "

Usage:

    x.splitSubPop(subPop, sizes, keepOrder=True)

Details:

    Split subpopulation subPop into subpopulations of given sizes,
    which should add up to the size of subpopulation subPop.
    Alternatively, sizes can be a list of proportions (add up to 1)
    from which the sizes of new subpopulations are determined. By
    default, subpopulation indexes will be adjusted so that
    individuals can keep their original order. That is to say, if
    subpopulation 1 of a population having four subpopulations is
    split into three subpopulation, the new subpopulation ID would be
    0, 1.1->1, 1.2->2, 1.3->3, 2->4, 3->5. If keepOrder is set to
    False, the subpopulation IDs of existing subpopulations will not
    be changed so the new subpopulation IDs of the previous example
    would be 0, 1.1->1, 2, 3, 1.2->4, 1.3->5.

"; 

%feature("docstring") simuPOP::population::removeEmptySubPops "

Usage:

    x.removeEmptySubPops()

Details:

    remove empty subpopulations by adjusting subpopulation IDs.

"; 

%feature("docstring") simuPOP::population::removeSubPops "

Usage:

    x.removeSubPops(subPops)

Details:

    Remove all individuals from subpopulations subPops. The removed
    subpopulations will have size zero, and can be removed by function
    removeEmptySubPops.

"; 

%feature("docstring") simuPOP::population::removeIndividuals "

Usage:

    x.removeIndividuals(inds)

Details:

    remove individuals inds (absolute indexes) from the current
    population. A subpopulation will be kept even if all individuals
    from it are removed.

"; 

%feature("docstring") simuPOP::population::mergeSubPops "

Usage:

    x.mergeSubPops(subPops=[])

Details:

    Merge subpopulations subPops. If subPops is empty (default), all
    subpopulations will be merged. Subpopulations subPops do not have
    to be adjacent to each other. The ID of the first subpopulation in
    parameter subPops will become the ID of the new large
    subpopulation. Other subpopulations will keep their IDs although
    their sizes become zero. Function removeEmptySubPops can be used
    to remove these empty subpopulation.

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
    simuPOP will try to assign default names, and raise a  ValueError
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

%ignore simuPOP::population::reorderSubPops(const vectoru &order=vectoru(), const vectoru &rank=vectoru(), bool removeEmptySubPops=false);

%ignore simuPOP::population::newPopByIndID(int keepAncestralPops=-1, const vectori &id=vectori(), bool removeEmptySubPops=false);

%feature("docstring") simuPOP::population::removeLoci "

Usage:

    x.removeLoci(loci)

Details:

    Remove loci (absolute indexes) and genotypes at these loci from
    the current population.

"; 

%feature("docstring") simuPOP::population::pushAndDiscard "Obsolete or undocumented function."

%feature("docstring") simuPOP::population::ancestralGens "

Usage:

    x.ancestralGens()

Details:

    Return the actual number of ancestral generations stored in a
    population, which does not necessarily equal to the number set by
    setAncestralDepth().

"; 

%ignore simuPOP::population::curAncestralGen() const ;

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
    values will be reused if its length is smaller than  popSize().

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

    Add information fields fields to a population and initialize their
    values to init. If an information field alreay exists, it will be
    re-initialized.

"; 

%feature("docstring") simuPOP::population::setInfoFields "

Usage:

    x.setInfoFields(fields, init=0)

Details:

    Set information fields fields to a population and initialize them
    with value init. All existing information fields will be removed.

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
    always be called to restore the correct order of ancestral
    generations.

"; 

%ignore simuPOP::population::equalTo(const population &rhs);

%ignore simuPOP::population::sortIndividuals(bool infoOnly=false);

%feature("docstring") simuPOP::population::save "

Usage:

    x.save(filename)

Details:

    Save population to a file filename. The population can be restored
    from this file, using a global function LoadPopulation(filename).

"; 

%ignore simuPOP::population::load(const string &filename);

%ignore simuPOP::population::selectionOn() const ;

%ignore simuPOP::population::selectionOn(UINT sp) const ;

%feature("docstring") simuPOP::population::turnOffSelection "Obsolete or undocumented function."

%ignore simuPOP::population::turnOnSelection(UINT sp);

%ignore simuPOP::population::turnOnSelection();

%ignore simuPOP::population::rep();

%ignore simuPOP::population::setRep(int rep, bool setVar=true);

%ignore simuPOP::population::gen();

%ignore simuPOP::population::setGen(ULONG gen, bool setVar=true);

%feature("docstring") simuPOP::population::vars "

Usage:

    x.vars(subPop=-1)

Details:

    return variables of a population. If subPop is given, return a
    dictionary for specified subpopulation.

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

%ignore simuPOP::population::varsAsString() const ;

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

%ignore simuPOP::proportionSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const ;

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

    Python expressions/statements will be executed when  pyEval is
    applied to a population by using parameters expr/stmts. Statements
    can also been executed when  pyEval is created and destroyed or
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
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

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

    deep copy of a  pyEval operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyEval::apply "

Description:

    apply the  pyEval operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyEval::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  pyEval operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyEval::name "

Description:

    return the name of an expression

Usage:

    x.name()

Details:

    The name of a  pyEval operator is given by an optional parameter
    name. It can be used to identify this  pyEval operator in debug
    output, or in the dryrun mode of  simulator::evolve.

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
      end=-1, step=1, at=[], rep=REP_ALL, infoFields=[])

Details:

    Please refer to class  pyEval for parameter descriptions.

"; 

%feature("docstring") simuPOP::pyExec::~pyExec "

Usage:

    x.~pyExec()

"; 

%feature("docstring") simuPOP::pyExec::clone "

Description:

    deep copy of a  pyExec operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyExec::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  pyExec operator

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
    by  population::individuals() and population::individuals(subPop)

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

    This operator is similar to a  pyOperator but works at the
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
      formOffGenotype=False, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, infoFields=[])

Arguments:

    func:           a Python function that accepts an individual and
                    optional genotype and parameters.
    param:          any Python object that will be passed to func
                    after pop parameter. Multiple parameters can be
                    passed as a tuple.
    infoFields:     if given, func is expected to return an array of
                    the same length and fill these infoFields of an
                    individual.

Example:

Testsrc_pyIndOperator.log Applying a <tt> pyIndOperator</tt>

"; 

%feature("docstring") simuPOP::pyIndOperator::~pyIndOperator "

Description:

    destructor

Usage:

    x.~pyIndOperator()

"; 

%ignore simuPOP::pyIndOperator::pyIndOperator(const pyIndOperator &rhs);

%feature("docstring") simuPOP::pyIndOperator::clone "

Description:

    deep copy of a  pyIndOperator operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyIndOperator::apply "

Description:

    apply the  pyIndOperator operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyIndOperator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  pyIndOperator operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyInit "

Function form:

    PyInit

Description:

    A python operator that uses a user-defined function to initialize
    individuals.

Details:

    This is a hybrid initializer. Users of this operator must supply a
    Python function with parameters allele, ploidy and subpopulation
    indexes (index, ploidy, subPop), and return an allele value. This
    operator will loop through all individuals in each subpopulation
    and call this function to initialize populations. The arrange of
    parameters allows different initialization scheme for each
    subpopulation.

"; 

%feature("docstring") simuPOP::pyInit::pyInit "

Description:

    initialize populations using given user function

Usage:

    pyInit(func, subPop=[], loci=[], atPloidy=-1, indRange=[],
      maleFreq=0.5, sex=[], stage=PreMating, begin=0, end=1, step=1,
      at=[], rep=REP_ALL, infoFields=[])

Arguments:

    func:           a Python function with parameter (index, ploidy,
                    subPop), where
                    * index is the allele index ranging from 0 to
                    totNumLoci-1;
                    * ploidy is the index of the copy of chromosomes;
                    * subPop is the subpopulation index. The return
                    value of this function should be an integer.
    loci:           a vector of locus indexes. If empty, apply to all
                    loci.
    locus:          a shortcut to loci.
    atPloidy:       initialize which copy of chromosomes. Default to
                    all.
    stage:          default to PreMating.

Example:

Testsrc_pyInit.log Operator <tt> pyInit</tt>

"; 

%feature("docstring") simuPOP::pyInit::~pyInit "

Usage:

    x.~pyInit()

"; 

%ignore simuPOP::pyInit::pyInit(const pyInit &rhs);

%feature("docstring") simuPOP::pyInit::clone "

Description:

    deep copy of the operator  pyInit

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyInit::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator  pyInit

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyInit::apply "

Description:

    apply this operator to population pop

Usage:

    x.apply(pop)

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
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

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

    deep copy of a  pyMigrator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyMigrator::apply "

Description:

    apply a  pyMigrator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyMigrator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  pyMigrator

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
    return a new allele state given an old state.  pyMutator will
    choose an allele as usual and call your function to mutate it to
    another allele.

"; 

%feature("docstring") simuPOP::pyMutator::pyMutator "

Description:

    create a  pyMutator

Usage:

    pyMutator(rate=[], loci=[], maxAllele=0, func=None, output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, infoFields=[])

Example:

Testsrc_pyMutator.log Operator <tt> pyMutator</tt>

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

    deep copy of a  pyMutator

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
    of the  pyMutator

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

Example:

Testsrc_pyOperator.log Operator <tt> pyOperator</tt>

"; 

%feature("docstring") simuPOP::pyOperator::pyOperator "

Description:

    Python operator, using a function that accepts a population
    object.

Usage:

    pyOperator(func, param=None, stage=PostMating,
      formOffGenotype=False, passOffspringOnly=False, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, infoFields=[])

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
    passOffspringOnly:if True,  pyOperator will expect a function of
                    form func(off [,param]), instead of func(pop, off,
                    dad, mom [, param]) which is used when
                    passOffspringOnly is False. Because many during-
                    mating  pyOperator only need access to offspring,
                    this will improve efficiency. Default to False.

Note:

    * Output to output or outputExpr is not supported. That is to say,
    you have to open/close/append to files explicitly in the Python
    function. Because files specified by output or outputExpr are
    controlled (opened/closed) by simulators, they should not be
    manipulated in a  pyOperator operator.
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

%feature("docstring") simuPOP::pyOperator::clone "

Description:

    deep copy of a  pyOperator operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyOperator::apply "

Description:

    apply the  pyOperator operator to one population

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::pyOperator::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::pyOperator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  pyOperator operator

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

    Create a  pyOutput operator that outputs a given string.

Usage:

    pyOutput(str=\"\", output=\">\", outputExpr=\"\", stage=PostMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, infoFields=[])

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

Description:

    deep copy of an operator

Usage:

    x.clone()

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
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, infoFields=[])

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

Example:

Testsrc_pyPenetrance.log Operator <tt> pyPenetrance</tt>

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

%feature("docstring") simuPOP::pyPenetrance::penet "

Description:

    currently assuming diploid

Usage:

    x.penet(ind)

"; 

%feature("docstring") simuPOP::pyPenetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the Python penetrance operator

Usage:

    x.__repr__()

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
      begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[\"qtrait\"])

Details:

    Please refer to  quanTrait for other parameter descriptions.

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

%feature("docstring") simuPOP::pyQuanTrait::qtrait "

Description:

    currently assuming diploid

Usage:

    x.qtrait(ind)

"; 

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
      end=-1, step=1, at=[], rep=REP_ALL, infoFields=[\"fitness\"])

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

Example:

Testsrc_pySelector.log Operator <tt> pySelector</tt>

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

    deep copy of a  pySelector

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pySelector::indFitness "

Description:

    calculate/return the fitness value, currently assuming diploid

Usage:

    x.indFitness(ind, gen)

"; 

%feature("docstring") simuPOP::pySelector::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  pySelector

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

    creates a  pyTagger that works on specified information fields

Usage:

    pyTagger(func=None, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      output=\"\", outputExpr=\"\", infoFields=[])

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

    deep copy of a  pyTagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyTagger::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  pyTagger

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyTagger::applyDuringMating "

Description:

    apply the  pyTagger

Usage:

    x.applyDuringMating(pop, offspring, dad=None, mom=None)

"; 

%feature("docstring") simuPOP::quanTrait "

Description:

    base class of quantitative trait

Details:

    Quantitative trait is the measure of certain phenotype for given
    genotype. Quantitative trait is similar to penetrance in that the
    consequence of penetrance is binary: affected or unaffected; while
    it is continuous for quantitative trait.
    In  simuPOP, different operators or functions were implemented to
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
      step=1, at=[], rep=REP_ALL, infoFields=[\"qtrait\"])

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

%feature("docstring") simuPOP::quanTrait::qtrait "

Description:

    calculate/return quantitative trait etc.

Usage:

    x.qtrait()

"; 

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

%feature("docstring") simuPOP::randomMating "

Applicability: diploid only

Description:

    a mating scheme of basic sexually random mating

Details:

    In this scheme, sex information is considered for each individual,
    and ploidy is always 2. Within each subpopulation, males and
    females are randomly chosen. Then randomly get one copy of
    chromosomes from father and mother. If only one sex exists in a
    subpopulation, a parameter (contWhenUniSex) can be set to
    determine the behavior. Default to continuing without warning.

"; 

%feature("docstring") simuPOP::randomMating::randomMating "

Usage:

    randomMating(numOffspring=1., numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, sexParam=0.5,
      sexMode=MATE_RandomSex, newSubPopSize=[],
      newSubPopSizeFunc=None, newSubPopSizeExpr=\"\",
      contWhenUniSex=True, subPop=[], weight=0)

Details:

    Please refer to class mating for descriptions of other parameters.

Arguments:

    contWhenUniSex: continue when there is only one sex in the
                    population. Default to True.

"; 

%feature("docstring") simuPOP::randomMating::clone "

Description:

    deep copy of a random mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::randomMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the random mating scheme

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
    can be chosen multiple times. In case that replacement=false,
    paremeter replenish=true allows restart of the process if all
    parents are exhausted. Note that selection is not allowed when
    replacement=false because this poses a particular order on
    individuals in the offspring generation.

"; 

%feature("docstring") simuPOP::randomParentChooser::randomParentChooser "

Usage:

    randomParentChooser(replacement=True, replenish=False)

Arguments:

    replacement:    if replacement is false, a parent can not be
                    chosen more than once.
    replenish:      if all parent has been chosen, choose from the
                    whole parental population again.

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
    parameter replacement is false, a chosen pair of parents can no
    longer be selected. This feature can be used to simulate monopoly.
    If replenish is true, a sex group can be replenished when it is
    exhausted. Note that selection is not allowed in the case of
    monopoly because this poses a particular order on individuals in
    the offspring generation. This parents chooser also allows
    polygamous mating by reusing a parent multiple times when
    returning parents, and allows specification of a few alpha
    individuals who will be the only mating individuals in their sex
    group.

"; 

%feature("docstring") simuPOP::randomParentsChooser::randomParentsChooser "

Usage:

    randomParentsChooser(replacement=True, replenish=False,
      polySex=Male, polyNum=1, alphaSex=Male, alphaNum=0,
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
    defined in ranges. For example,  rangeSplitter(ranges=[[0, 20],
    [40, 50]]) defines two VSPs. The first VSP consists of individuals
    0, 1, ..., 19, and the sceond VSP consists of individuals 40, 41,
    ..., 49. Note that a nested list has to be used even if only one
    range is defined.

"; 

%feature("docstring") simuPOP::rangeSplitter::clone "Obsolete or undocumented function."

%ignore simuPOP::rangeSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const ;

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

Description:

    recombination and conversion

Details:

    In  simuPOP, only one recombinator is provided. Recombination
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

    recombinator(intensity=-1, rate=[], afterLoci=[],
      maleIntensity=-1, maleRate=[], maleAfterLoci=[], convProb=0,
      convMode=CONVERT_NumMarkers, convParam=1., begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, infoFields=[])

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
                    *  simuPOP does not impose a unit for marker
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

Example:

Testsrc_recombinator.log Operator <tt>recombinator</tt>

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

%feature("docstring") simuPOP::recombinator::recCount "

Description:

    return recombination count at a locus (only valid in standard
    modules)

Usage:

    x.recCount(locus)

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

%ignore simuPOP::recombinator::applyDuringMating(population &pop, RawIndIterator offspring, individual *dad=NULL, individual *mom=NULL);

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
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

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

    deep copy of a  resizeSubPops operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::resizeSubPops::apply "

Description:

    apply a  resizeSubPops operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::resizeSubPops::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  resizeSubPops operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::RNG "

Description:

    random number generator

Details:

    This random number generator class wraps around a number of random
    number generators from GNU Scientific Library. You can obtain and
    change system random number generator through the  rng() function.
    Or create a separate random number generator and use it in your
    script.

"; 

%feature("docstring") simuPOP::RNG::RNG "

Description:

    RNG used by  simuPOP.

Usage:

    RNG(rng=None, seed=0)

"; 

%feature("docstring") simuPOP::RNG::~RNG "

Usage:

    x.~RNG()

"; 

%feature("docstring") simuPOP::RNG::setRNG "

Description:

    choose an random number generator, or set seed to the current  RNG

Usage:

    x.setRNG(rng=None, seed=0)

Arguments:

    rng:            name of the  RNG. If rng is not given,
                    environmental variable GSL_RNG_TYPE will be used
                    if it is available. Otherwise,  RNGmt19937 will be
                    used.
    seed:           random seed. If not given, /dev/urandom,
                    /dev/random, system time will be used, depending
                    on availability, in that order. Note that windows
                    system does not have /dev so system time is used.

"; 

%feature("docstring") simuPOP::RNG::name "

Description:

    return  RNG name

Usage:

    x.name()

"; 

%feature("docstring") simuPOP::RNG::seed "

Description:

    return the seed of this  RNG

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

    Maximum value of this  RNG.

Usage:

    x.max()

"; 

%feature("docstring") simuPOP::RNG::__repr__ "

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::RNG::randGet "

Description:

    return a random number in the range of [0, 2, ...  max()-1]

Usage:

    x.randGet()

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
      rep=REP_ALL, infoFields=[])

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

Description:

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::savePopulation::apply "

Description:

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

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
    selection.  simuPOP employs an 'ability-to-mate' approach. Namely,
    the probability that an individual will be chosen for mating is
    proportional to its fitness value. More specifically,
    *  PreMating selectors assign fitness values to each individual,
    and mark part or all subpopulations as under selection.
    *  during sexless mating (e.g.  binomialSelection mating scheme),
    individuals are chosen at probabilities that are proportional to
    their fitness values. If there are  $ N $ individuals with fitness
    values  $ f_{i},i=1,...,N $, individual  $ i $ will have
    probability  $ \\frac{f_{i}}{\\sum_{j}f_{j}} $ to be chosen and
    passed to the next generation.
    *  during  randomMating, males and females are separated. They are
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
    of the selector in a  pyOperator, make sure to set the stage of
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
      at=[], rep=REP_ALL, infoFields=[\"fitness\"])

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

%feature("docstring") simuPOP::selfingOffspringGenerator "

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

%feature("docstring") simuPOP::selfingOffspringGenerator::selfingOffspringGenerator "

Usage:

    selfingOffspringGenerator(numOffspring=1, numOffspringFunc=None,
      maxNumOffspring=1, mode=MATE_NumOffspring, sexParam=0.5,
      sexMode=MATE_RandomSex)

"; 

%feature("docstring") simuPOP::selfingOffspringGenerator::clone "

Usage:

    x.clone()

"; 

%ignore simuPOP::selfingOffspringGenerator::generateOffspring(population &pop, individual *parent, individual *, RawIndIterator &offBegin, RawIndIterator &offEnd, vector< baseOperator * > &ops);

%feature("docstring") simuPOP::selfMating "

Applicability: diploid only

Description:

    a mating scheme of selfing

Details:

    In this mating scheme, a parent is choosen randomly, acts both as
    father and mother in the usual random mating. The parent is chosen
    randomly, regardless of sex. If selection is turned on, the
    probability that an individual is chosen is proportional to
    his/her fitness.

"; 

%feature("docstring") simuPOP::selfMating::selfMating "

Description:

    create a self mating scheme

Usage:

    selfMating(numOffspring=1., numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, sexParam=0.5,
      sexMode=MATE_RandomSex, newSubPopSize=[],
      newSubPopSizeFunc=None, newSubPopSizeExpr=\"\",
      contWhenUniSex=True, subPop=[], weight=0)

Details:

    Please refer to class mating for descriptions of other parameters.

Arguments:

    contWhenUniSex: continue when there is only one sex in the
                    population. Default to True.

"; 

%feature("docstring") simuPOP::selfMating::~selfMating "

Description:

    destructor

Usage:

    x.~selfMating()

"; 

%feature("docstring") simuPOP::selfMating::clone "

Description:

    deep copy of a self mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::selfMating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::selfMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the self mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::selfMating::mateSubPop(population &pop, SubPopID subPop, RawIndIterator offBegin, RawIndIterator offEnd, vector< baseOperator * > &ops);

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

    create a  setAncestralDepth operator

Usage:

    setAncestralDepth(depth, output=\">\", outputExpr=\"\",
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

"; 

%feature("docstring") simuPOP::setAncestralDepth::~setAncestralDepth "

Description:

    destructor

Usage:

    x.~setAncestralDepth()

"; 

%feature("docstring") simuPOP::setAncestralDepth::clone "

Description:

    deep copy of a  setAncestralDepth operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::setAncestralDepth::apply "

Description:

    apply the  setAncestralDepth operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::setAncestralDepth::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  setAncestralDepth operator

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

%ignore simuPOP::sexSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const ;

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

%feature("docstring") simuPOP::sexTagger "

Description:

    Tagging sex status.

Details:

    This is a simple post-mating tagger that write sex status to a
    file. By default, 1 for Male, 2 for Female.

"; 

%feature("docstring") simuPOP::sexTagger::sexTagger "

Usage:

    sexTagger(code=[], begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      stage=PostMating, output=\">\", outputExpr=\"\", infoFields=[])

Arguments:

    code:           code for Male and Female, default to 1 and 2,
                    respectively. This is used by Linkage format.

"; 

%feature("docstring") simuPOP::sexTagger::apply "

Description:

    add a newline

Usage:

    x.apply(pop)

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

%ignore simuPOP::SharedVariables::asString() const ;

%feature("docstring") simuPOP::SharedVariables::fromString "

Usage:

    x.fromString(vars)

"; 

%feature("docstring") simuPOP::simulator "

Description:

    simulator manages several replicates of a population, evolve them
    using given mating scheme and operators

Details:

    Simulators combine three important components of  simuPOP:
    population, mating scheme and operator together. A simulator is
    created with an instance of population, a replicate number rep and
    a mating scheme. It makes rep number of replicates of this
    population and control the evolutionary process of them.
    The most important function of a simulator is  evolve(). It
    accepts an array of operators as its parameters, among which,
    preOps and postOps will be applied to the populations at the
    beginning and the end of evolution, respectively, whereas ops will
    be applied at every generation.
    A simulator separates operators into pre-, during-, and post-
    mating operators. During evolution, a simulator first apply all
    pre-mating operators and then call the mate() function of the
    given mating scheme, which will call during-mating operators
    during the birth of each offspring. After mating is completed,
    post-mating operators are applied to the offspring in the order at
    which they appear in the operator list.
    Simulators can evolve a given number of generations (the end
    parameter of evolve), or evolve indefinitely until a certain type
    of operators called terminator terminates it. In this case, one or
    more terminators will check the status of evolution and determine
    if the simulation should be stopped. An obvious example of such a
    terminator is a fixation-checker.
    A simulator can be saved to a file in the format of 'txt', 'bin',
    or 'xml'. This allows you to stop a simulator and resume it at
    another time or on another machine.

"; 

%feature("docstring") simuPOP::simulator::simulator "

Description:

    create a simulator

Usage:

    simulator(pop, matingScheme, stopIfOneRepStops=False,
      applyOpToStoppedReps=False, rep=1)

Arguments:

    population:     a population created by population() function.
                    This population will be copied rep times to the
                    simulator. Its content will not be changed.
    matingScheme:   a mating scheme
    rep:            number of replicates. Default to 1.
    applyOpToStoppedReps:If set, the simulator will continue to apply
                    operators to all stopped replicates until all
                    replicates are marked 'stopped'.
    stopIfOneRepStops:If set, the simulator will stop evolution if one
                    replicate stops.

"; 

%feature("docstring") simuPOP::simulator::~simulator "

Description:

    destroy a simulator along with all its populations

Usage:

    x.~simulator()

Note:

    pop = simulator::population() returns temporary reference to an
    internal population. After a simulator evolves another genertion
    or after the simulator is destroyed, this referenced population
    should not be used.

"; 

%ignore simuPOP::simulator::simulator(const simulator &rhs);

%feature("docstring") simuPOP::simulator::clone "

Description:

    deep copy of a simulator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::simulator::addInfoField "

Description:

    add an information field to all replicates

Usage:

    x.addInfoField(field, init=0)

Details:

    Add an information field to all replicate, and to the simulator
    itself. This is important because all populations must have the
    same genotypic information as the simulator. Adding an information
    field to one or more of the replicates will compromise the
    integrity of the simulator.

Arguments:

    field:          information field to be added

"; 

%feature("docstring") simuPOP::simulator::addInfoFields "

Description:

    add information fields to all replicates

Usage:

    x.addInfoFields(fields, init=0)

Details:

    Add given information fields to all replicate, and to the
    simulator itself.

"; 

%feature("docstring") simuPOP::simulator::setAncestralDepth "

Description:

    set ancestral depth of all replicates

Usage:

    x.setAncestralDepth(depth)

"; 

%feature("docstring") simuPOP::simulator::pop "

Description:

    Return a reference to the rep replicate of this simulator.

Usage:

    x.pop(rep)

Arguments:

    rep:            the index number of replicate which will be
                    accessed

Note:

    The returned reference is temporary in the sense that the refered
    population will be invalid after another round of evolution. If
    you would like to get a persistent population, please use
    getPopulation(rep).

"; 

%feature("docstring") simuPOP::simulator::getPopulation "

Description:

    return a copy of population rep

Usage:

    x.getPopulation(rep, destructive=False)

Details:

    By default return a cloned copy of population rep of the
    simulator. If destructive==True, the population is extracted from
    the simulator, leaving a defunct simulator.

Arguments:

    rep:            the index number of the replicate which will be
                    obtained
    destructive:    if true, destroy the copy of population within
                    this simulator. Default to false.
                    getPopulation(rep, true) is a more efficient way
                    to get hold of a population when the simulator
                    will no longer be used.

"; 

%feature("docstring") simuPOP::simulator::setMatingScheme "

Description:

    set a new mating scheme

Usage:

    x.setMatingScheme(matingScheme)

"; 

%ignore simuPOP::simulator::setPopulation(population &pop, UINT rep);

%ignore simuPOP::simulator::curRep() const ;

%feature("docstring") simuPOP::simulator::numRep "

Description:

    return the number of replicates

Usage:

    x.numRep()

"; 

%feature("docstring") simuPOP::simulator::gen "

Description:

    return the current generation number

Usage:

    x.gen()

"; 

%feature("docstring") simuPOP::simulator::setGen "

Description:

    set the current generation. Usually used to reset a simulator.

Usage:

    x.setGen(gen)

Arguments:

    gen:            new generation index number

"; 

%feature("docstring") simuPOP::simulator::step "

Description:

    evolve steps generation

Usage:

    x.step(ops=[], preOps=[], postOps=[], steps=1, dryrun=False)

"; 

%feature("docstring") simuPOP::simulator::evolve "

Description:

    evolve all replicates of the population, subject to operators

Usage:

    x.evolve(ops, preOps=[], postOps=[], end=-1, gen=-1,
      dryrun=False)

Details:

    Evolve to the end generation unless end=-1. An operator
    (terminator) may stop the evolution earlier.
    ops will be applied to each replicate of the population in the
    order of:
    *  all pre-mating opertors
    *  during-mating operators called by the mating scheme at the
    birth of each offspring
    *  all post-mating operators If any pre- or post-mating operator
    fails to apply, that replicate will be stopped. The behavior of
    the simulator will be determined by flags applyOpToStoppedReps and
    stopIfOneRepStopss.

Arguments:

    ops:            operators that will be applied at each generation,
                    if they are active at that generation. (Determined
                    by the begin, end, step and at parameters of the
                    operator.)
    preOps:         operators that will be applied before evolution.
                    evolve() function will not check if they are
                    active.
    postOps:        operators that will be applied after evolution.
                    evolve() function will not check if they are
                    active.
    gen:            generations to evolve. Default to -1. In this
                    case, there is no ending generation and a
                    simulator will only be ended by a terminator. Note
                    that simu.gen() refers to the begining of a
                    generation, and starts at 0.
    dryrun:         dryrun mode. Default to False.

Note:

    When gen = -1, you can not specify negative generation parameters
    to operators. How would an operator know which genertion is the -1
    genertion if no ending genertion is given?

"; 

%ignore simuPOP::simulator::apply(const vectorop ops, bool dryrun=false);

%ignore simuPOP::simulator::setStopIfOneRepStops(bool on=true);

%ignore simuPOP::simulator::stopIfOneRepStops();

%ignore simuPOP::simulator::setApplyOpToStoppedReps(bool on=true);

%ignore simuPOP::simulator::applyOpToStoppedReps();

%feature("docstring") simuPOP::simulator::vars "

Description:

    Return the local namespace of population rep, equivalent to
    x.population(rep).vars(subPop).

Usage:

    x.vars(rep, subPop=-1)

"; 

%feature("docstring") simuPOP::simulator::saveSimulator "

Description:

    save simulator in 'txt', 'bin' or 'xml' format

Usage:

    x.saveSimulator(filename, format=\"\", compress=True)

Arguments:

    filename:       filename to save the simulator. Default to simu.
    format:         obsolete parameter
    compress:       obsolete parameter

"; 

%ignore simuPOP::simulator::loadSimulator(string filename, string format="");

%feature("docstring") simuPOP::simulator::__repr__ "

Description:

    used by Python print function to print out the general information
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
      step=1, at=[], rep=REP_ALL, infoFields=[])

Details:

    The SMM is developed for allozymes. It provides better description
    for these kinds of evolutionary processes.  Please see class
    mutator for the descriptions of other parameters.

Arguments:

    incProb:        probability to increase allele state. Default to
                    0.5.

Example:

Testsrc_smmMutator.log Operator <tt> smmMutator</tt>

"; 

%feature("docstring") simuPOP::smmMutator::~smmMutator "

Usage:

    x.~smmMutator()

"; 

%ignore simuPOP::smmMutator::mutate(AlleleRef allele);

%feature("docstring") simuPOP::smmMutator::clone "

Description:

    deep copy of a  smmMutator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::smmMutator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  smmMutator

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

    splitSubPop(which=0, sizes=[], proportions=[], keepOrder=True,
      randomize=True, stage=PreMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, infoFields=[])

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

Example:

Testsrc_splitSubPop.log Operator <tt> splitSubPop</tt>

"; 

%feature("docstring") simuPOP::splitSubPop::~splitSubPop "

Description:

    destructor

Usage:

    x.~splitSubPop()

"; 

%feature("docstring") simuPOP::splitSubPop::clone "

Description:

    deep copy of a  splitSubPop operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::splitSubPop::apply "

Description:

    apply a  splitSubPop operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::splitSubPop::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  splitSubPop operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::spread "

Function form:

    Spread

Description:

    copy the genotype of an individual to all individuals

Details:

    Function Spread(ind, subPop) spreads the genotypes of ind to all
    individuals in an array of subpopulations. The default value of
    subPop is the subpopulation where ind resides.

"; 

%feature("docstring") simuPOP::spread::spread "

Description:

    copy genotypes of ind to all individuals in subPop

Usage:

    spread(ind, subPop=[], stage=PreMating, begin=0, end=1, step=1,
      at=[], rep=REP_ALL, infoFields=[])

Example:

Testsrc_spread.log Operator <tt>spread</tt>

"; 

%feature("docstring") simuPOP::spread::~spread "

Usage:

    x.~spread()

"; 

%feature("docstring") simuPOP::spread::clone "

Description:

    deep copy of the operator spread

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::spread::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator spread

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::spread::apply "

Description:

    apply this operator to population pop

Usage:

    x.apply(pop)

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
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

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
                    the number of all genotype 1x or x1 where  does
                    not equal to 1. All other genotypes such as 02 are
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
                    $ h_{exp}=1-p_{i}^{2}, $ where  $ p_i $ is the
                    allele frequency of allele  $ i $. The following
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
    LD:             calculate linkage disequilibria  $ LD $,  $ LD' $
                    and  $ r^{2} $, given LD=[ [loc1, loc2], [ loc1,
                    loc2, allele1, allele2], ... ]. For each item
                    [loc1, loc2, allele1, allele2],  $ D $,  $ D' $
                    and  $ r^{2} $ will be calculated based on allele1
                    at loc1 and allele2 at loc2. If only two loci are
                    given, the LD values are averaged over all allele
                    pairs. For example, for allele  $ A $ at locus 1
                    and allele  $ B $ at locus 2,
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
    Fst:            calculate  $ F_{st} $,  $ F_{is} $,  $ F_{it} $.
                    For example, Fst = [0,1,2] will calculate  $
                    F_{st} $,  $ F_{is} $,  $ F_{it} $ based on
                    alleles at loci 0, 1, 2. The locus-specific values
                    will be used to calculate AvgFst, which is an
                    average value over all alleles (Weir & Cockerham,
                    1984). Terms and values that match Weir &
                    Cockerham are:
                    *  $ F $ ( $ F_{IT} $) the correlation of genes
                    within individuals (inbreeding);
                    *  $ \\theta $ ( $ F_{ST} $) the correlation of
                    genes of difference individuals in the same
                    population (will evaluate for each subpopulation
                    and the whole population)
                    *  $ f $ ( $ F_{IS} $) the correlation of genes
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
      end=-1, step=1, at=[], rep=REP_ALL, infoFields=[])

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

%feature("docstring") simuPOP::StopIteration "

Description:

    exception, thrown if out of memory

"; 

%feature("docstring") simuPOP::StopIteration::StopIteration "

Usage:

    StopIteration(msg)

"; 

%ignore simuPOP::StreamElem;

%feature("docstring") simuPOP::StreamElem::StreamElem "

Usage:

    StreamElem(name, readable, realAppend, useString)

Arguments:

    name:           filename
    readable:       iostream or just ostream
    realAppend:     whether or not keep old content when open an
                    existing file
    useString:      use a stringstream rather than a file.

"; 

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

%feature("docstring") simuPOP::StreamElem::info "

Description:

    mostly for debug purposes

Usage:

    x.info()

"; 

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

    tagger(output=\"\", outputExpr=\"\", begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, infoFields=[])

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

%feature("docstring") simuPOP::tagger::apply "

Description:

    add a newline

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::terminateIf "

Description:

    terminate according to a condition

Details:

    This operator terminates the evolution under certain conditions.
    For example,  terminateIf(condition='alleleFreq[0][1]<0.05',
    begin=100) terminates the evolution if the allele frequency of
    allele 1 at locus 0 is less than 0.05. Of course, to make this
    opertor work, you will need to use a stat operator before it so
    that variable alleleFreq exists in the local namespace.
    When the value of condition is True, a shared variable
    var=\"terminate\" will be set to the current generation.

"; 

%feature("docstring") simuPOP::terminateIf::terminateIf "

Description:

    create a  terminateIf terminator

Usage:

    terminateIf(condition=\"\", message=\"\", var=\"terminate\",
      output=\"\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::terminateIf::clone "

Description:

    deep copy of a  terminateIf terminator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::terminateIf::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  terminateIf terminator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::terminateIf::apply "

Description:

    apply the  terminateIf terminator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::terminateIf::~terminateIf "

Usage:

    x.~terminateIf()

"; 

%feature("docstring") simuPOP::terminator "

Description:

    Base class of all terminators.

Details:

    Teminators are used to see if an evolution is running as expected,
    and terminate the evolution if a certain condition fails.

"; 

%feature("docstring") simuPOP::terminator::terminator "

Description:

    create a terminator

Usage:

    terminator(message=\"\", output=\">\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      infoFields=[])

Arguments:

    message:        a message that will be displayed when the
                    evolution is terminated.

"; 

%feature("docstring") simuPOP::terminator::~terminator "

Description:

    destructor

Usage:

    x.~terminator()

"; 

%feature("docstring") simuPOP::terminator::clone "

Description:

    deep copy of a terminator

Usage:

    x.clone()

"; 

%ignore simuPOP::terminator::message();

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
      end=-1, step=1, at=[], rep=REP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::ticToc::~ticToc "

Description:

    destructor

Usage:

    x.~ticToc()

"; 

%feature("docstring") simuPOP::ticToc::clone "

Description:

    deep copy of a  ticToc operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::ticToc::apply "

Description:

    apply the  ticToc operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::ticToc::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  ticToc operator

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

    create a  turnOffDebug operator

Usage:

    turnOffDebug(code, stage=PreMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::turnOffDebug::~turnOffDebug "

Description:

    destructor

Usage:

    x.~turnOffDebug()

"; 

%feature("docstring") simuPOP::turnOffDebug::clone "

Description:

    deep copy of a  turnOffDebug operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::turnOffDebug::apply "

Description:

    apply the  turnOffDebug operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::turnOffDebug::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  turnOffDebug operator

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
    *  use this  turnOnDebug operator The advantage of using this
    operator is that you can turn on debug at given generations.

"; 

%feature("docstring") simuPOP::turnOnDebug::turnOnDebug "

Description:

    create a  turnOnDebug operator

Usage:

    turnOnDebug(code, stage=PreMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::turnOnDebug::~turnOnDebug "

Description:

    destructor

Usage:

    x.~turnOnDebug()

"; 

%feature("docstring") simuPOP::turnOnDebug::clone "

Description:

    deep copy of a  turnOnDebug operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::turnOnDebug::apply "

Description:

    apply the  turnOnDebug operator to one population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::turnOnDebug::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  turnOnDebug operator

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
    different subpopulations, a  combinedSplitter should be.

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

    All VSP splitter defines a  clone() function to create an
    identical copy of itself.

"; 

%feature("docstring") simuPOP::vspSplitter::~vspSplitter "

Usage:

    x.~vspSplitter()

"; 

%ignore simuPOP::vspSplitter::activatedSubPop() const ;

%ignore simuPOP::vspSplitter::size(const population &pop, SubPopID subPop, SubPopID virtualSubPop) const ;

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

    LoadSimulator(file, mate, format=\"auto\")

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

%feature("docstring") simuPOP::mainVars "

Description:

    get main dictionary (user namespace)

Usage:

    mainVars()

"; 

%feature("docstring") simuPOP::moduleVars "

Description:

    get module dictionary (it is different than mainDict!

Usage:

    moduleVars()

"; 

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

%feature("docstring") simuPOP::setRNG "

Description:

    for backward compatibilit, will remove later

Usage:

    setRNG(rng=\"\", seed=0)

"; 

%feature("docstring") simuPOP::ListAllRNG "

Description:

    list the names of all available random number generators

Usage:

    ListAllRNG()

"; 

%ignore simuPOP::listAllRNG();

%feature("docstring") simuPOP::simuRev "

Description:

    return the revision number of this  simuPOP module. Can be used to
    test if a feature is available.

Usage:

    simuRev()

"; 

%feature("docstring") simuPOP::simuVer "

Description:

    return the version of this  simuPOP module

Usage:

    simuVer()

"; 

%feature("docstring") simuPOP::ModuleCompiler "

Description:

    return the compiler used to compile this  simuPOP module

Usage:

    ModuleCompiler()

"; 

%feature("docstring") simuPOP::ModuleDate "

Description:

    return the date when this  simuPOP module is compiled

Usage:

    ModuleDate()

"; 

%feature("docstring") simuPOP::ModulePyVersion "

Description:

    return the Python version this  simuPOP module is compiled for

Usage:

    ModulePyVersion()

"; 

%feature("docstring") simuPOP::ModulePlatForm "

Description:

    return the platform on which this  simuPOP module is compiled

Usage:

    ModulePlatForm()

"; 

%ignore simuPOP::initialize();

%feature("docstring") simuPOP::Optimized "

Description:

    return True if this  simuPOP module is optimized

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

    return the maximum allowed allele state of the current  simuPOP
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

    x.dvars(subPop=-1)

Details:

    Return a wrapper of Python dictionary returned by vars(subPop) so
    that dictionary keys can be accessed as attributes. For example
    pop.dvars().alleleFreq is equivalent to pop.vars()[\"alleleFreq\"].

"; 

%feature("docstring") simuPOP::simulator::dvars "

Usage:

    x.dvars(rep, subPop=-1)

Details:

    Return a wrapper of Python dictionary returned by vars(rep,
    subPop) so that dictionary keys can be accessed as attributes. For
    example simu.dvars(1).alleleFreq is equivalent to
    simu.vars(1)[\"alleleFreq\"].

"; 

