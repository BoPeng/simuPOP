%feature("docstring") simuPOP::affectedSibpairSample "

Function form:

    AffectedSibpairSample

Function form:

    AffectedSibpairSample

Description:

    draw an affected sibling pair  sample

Details:

    Special preparation for the  population is needed in order to use
    this operator. Obviously, to obtain affected sibling pairs, we
    need to know the parents and the affectedness status of each
    individual. Furthermore, to get parental genotype, the  population
    should have ancestralDepth at least 1. The most important problem,
    however, comes from the  mating scheme we are using.
    randomMating() is usually used for diploid populations. The
    realrandom mating requires that a  mating will generate only one
    offspring. Since parents are chosen with replacement, a parent can
    have multiple offspring with different parents. On the other hand,
    it is very unlikely that two offspring will have the same parents.
    The probability of having a sibling for an offspring is  $
    \\frac{1}{N^{2}} $ (if do not consider selection). Therefore, we
    will have to allow multiple offspring per  mating at the cost of
    small effective  population size.
    All these requirements come at a cost: multiple ancestral
    populations, determining affectedness status and tagging will slow
    down evolution; multiple offspring will reduce effective
    population size. Fortunately,  simuPOP is flexible enough to let
    all these happen only at the last several generations.

Details:

    Special preparation for the  population is needed in order to use
    this operator. Obviously, to obtain affected sibling pairs, we
    need to know the parents and the affectedness status of each
    individual. Furthermore, to get parental genotype, the  population
    should have ancestralDepth at least 1. The most important problem,
    however, comes from the  mating scheme we are using.
    randomMating() is usually used for diploid populations. The
    realrandom mating requires that a  mating will generate only one
    offspring. Since parents are chosen with replacement, a parent can
    have multiple offspring with different parents. On the other hand,
    it is very unlikely that two offspring will have the same parents.
    The probability of having a sibling for an offspring is  $
    \\frac{1}{N^{2}} $ (if do not consider selection). Therefore, we
    will have to allow multiple offspring per  mating at the cost of
    small effective  population size.
    All these requirements come at a cost: multiple ancestral
    populations, determining affectedness status and tagging will slow
    down evolution; multiple offspring will reduce effective
    population size. Fortunately,  simuPOP is flexible enough to let
    all these happen only at the last several generations.

"; 

%feature("docstring") simuPOP::affectedSibpairSample::affectedSibpairSample "

Description:

    draw an affected sibling pair  sample

Usage:

    affectedSibpairSample(size=[], chooseUnaffected=False,
      countOnly=False, name=\"sample\", nameExpr=\"\", times=1, saveAs=\"\",
      saveAsExpr=\"\", format=\"auto\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[\"father_idx\", \"mother_idx\"])

Details:

    Please refer to class  sample for other parameter descriptions.

Details:

    Please refer to class  sample for other parameter descriptions.

Arguments:

    size:           the number of affected sibling pairs to be
                    sampled. Can be a number or an array. If a number
                    is given, it is the total number of sibpairs,
                    ignoring  population structure. Otherwise, given
                    number of sibpairs are sampled from
                    subpopulations. If size is unspecified, this
                    operator will return all affected sibpairs.
    chooseUnaffected:instead of affected sibpairs, choose unaffected
                    families.
    countOnly:      set variables about number of affected sibpairs,
                    do not actually draw the  sample

"; 

%feature("docstring") simuPOP::affectedSibpairSample::~affectedSibpairSample "

Description:

    destructor

Usage:

    x.~affectedSibpairSample()

"; 

%feature("docstring") simuPOP::affectedSibpairSample::clone "

Description:

    deep copy of a  affectedSibpairSample operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::affectedSibpairSample::prepareSample "

Description:

    preparation before drawing a  sample

Usage:

    x.prepareSample(pop)

"; 

%feature("docstring") simuPOP::affectedSibpairSample::drawsample "

Description:

    draw a  sample

Usage:

    x.drawsample(pop)

"; 

%feature("docstring") simuPOP::affectedSibpairSample::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  affectedSibpairSample operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::baseOperator "

Description:

    base class of all classes that manipulate populations

Details:

    Operators are objects that act on populations. They can be applied
    to populations directly using their function forms, but they are
    usually managed and applied by a  simulator.
    Operators can be applied at different stages of the life cycle of
    a generation. More specifically, they can be applied at pre-,
    during-, post-mating, or a combination of these stages. Applicable
    stages are usually set by default but you can change it by setting
    stage=(PreMating|PostMating|DuringMating|PrePostMating) parameter.
    Some operators ignore stage parameter because they only work at
    one stage.
    Operators do not have to be applied at all generations. You can
    specify starting/ending generation, gaps between applicable
    generations, or even specific generations. For example, you might
    want to start applying migrations after certain burn-in
    generations, or calculate certain statistics only sparsely.
    Operators can have outputs. Output can be standard (terminal) or a
    file, which can vary with replicates and/or generations. Outputs
    from different operators can be accumulated to the same file to
    form table-like outputs.
    Operators are applied to every replicate of a  simulator by
    default. However, you can apply operators to one or a group of
    replicates using parameter rep or grp.
    Filenames can have the following format:
    * 'filename' this file will be overwritten each time. If two
    operators output to the same file, only the last one will succeed;
    * '>filename' the same as 'filename';
    * '>>filename' the file will be created at the beginning of
    evolution ( simulator::evolve) and closed at the end. Output from
    several operators is allowed;
    * '>>>filename' the same as '>>filename' except that the file will
    not be cleared at the beginning of evolution if it is not empty;
    * '>' standard output (terminal);
    * '' supress output.

Details:

    Operators are objects that act on populations. They can be applied
    to populations directly using their function forms, but they are
    usually managed and applied by a  simulator.
    Operators can be applied at different stages of the life cycle of
    a generation. More specifically, they can be applied at pre-,
    during-, post-mating, or a combination of these stages. Applicable
    stages are usually set by default but you can change it by setting
    stage=(PreMating|PostMating|DuringMating|PrePostMating) parameter.
    Some operators ignore stage parameter because they only work at
    one stage.
    Operators do not have to be applied at all generations. You can
    specify starting/ending generation, gaps between applicable
    generations, or even specific generations. For example, you might
    want to start applying migrations after certain burn-in
    generations, or calculate certain statistics only sparsely.
    Operators can have outputs. Output can be standard (terminal) or a
    file, which can vary with replicates and/or generations. Outputs
    from different operators can be accumulated to the same file to
    form table-like outputs.
    Operators are applied to every replicate of a  simulator by
    default. However, you can apply operators to one or a group of
    replicates using parameter rep or grp.
    Filenames can have the following format:
    * 'filename' this file will be overwritten each time. If two
    operators output to the same file, only the last one will succeed;
    * '>filename' the same as 'filename';
    * '>>filename' the file will be created at the beginning of
    evolution ( simulator::evolve) and closed at the end. Output from
    several operators is allowed;
    * '>>>filename' the same as '>>filename' except that the file will
    not be cleared at the beginning of evolution if it is not empty;
    * '>' standard output (terminal);
    * '' supress output.

"; 

%feature("docstring") simuPOP::baseOperator::baseOperator "

Description:

    common interface for all operators (this base operator does
    nothing by itself.)

Usage:

    baseOperator(output, outputExpr, stage, begin, end, step, at,
      rep, grp, infoFields)

Arguments:

    begin:          the starting generation. Default to 0. Negative
                    numbers are allowed.
    end:            stop applying after this generation. Negative
                    numbers are allowed.
    step:           the number of generations between active
                    generations. Default to 1.
    at:             an array of active generations. If given, stage,
                    begin, end, and step will be ignored.
    rep:            applicable replicates. It can be a valid replicate
                    number, REP_ALL (all replicates, default), or
                    REP_LAST (only the last replicate). REP_LAST is
                    useful in adding newlines to a table output.
    grp:            applicable group. Default to GRP_ALL. A group
                    number for each replicate is set by
                    simulator.__init__ or  simulator::setGroup().
    output:         a string of the output filename. Different
                    operators will have different default output (most
                    commonly '>' or '').
    outputExpr:     an expression that determines the output filename
                    dynamically. This expression will be evaluated
                    against a population's local namespace each time
                    when an output filename is required. For example,
                    \"'>>out%s_%s.xml' % (gen, rep)\"  will output to
                    >>>out1_1.xml  for replicate 1 at generation 1.

Note:

    Negative generation numbers are allowed for begin, end and at.
    They are intepretted as endGen + gen + 1. For example, begin = -2
    in simu.evolve(..., end=20) starts at generation 19.

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

%ignore simuPOP::baseOperator::isActive(UINT rep, UINT numRep, long gen, long end, int grp, bool repOnly=false);

%feature("docstring") simuPOP::baseOperator::applicableGroup "

Description:

    return applicable group

Usage:

    x.applicableGroup()

"; 

%feature("docstring") simuPOP::baseOperator::setApplicableGroup "

Description:

    set applicable group

Usage:

    x.setApplicableGroup(grp=GRP_ALL)

Details:

    Default to GRP_ALL (applicable to all groups). Otherwise, the
    operator is applicable to only one group of replicates. Groups can
    be set in  simulator::setGroup().

Details:

    Default to GRP_ALL (applicable to all groups). Otherwise, the
    operator is applicable to only one group of replicates. Groups can
    be set in  simulator::setGroup().

"; 

%feature("docstring") simuPOP::baseOperator::applicableReplicate "

Description:

    return applicable replicate

Usage:

    x.applicableReplicate()

"; 

%feature("docstring") simuPOP::baseOperator::setApplicableReplicate "

Description:

    set applicable replicate

Usage:

    x.setApplicableReplicate(rep)

"; 

%feature("docstring") simuPOP::baseOperator::setActiveGenerations "

Description:

    set applicable generation parameters: begin, end, step and at

Usage:

    x.setActiveGenerations(begin=0, end=-1, step=1, at=[])

"; 

%feature("docstring") simuPOP::baseOperator::setApplicableStage "

Description:

    set applicable stage. Another way to set stage parameter.

Usage:

    x.setApplicableStage(stage)

"; 

%feature("docstring") simuPOP::baseOperator::canApplyPreMating "

Description:

    set if this operator can be applied pre-mating

Usage:

    x.canApplyPreMating()

"; 

%feature("docstring") simuPOP::baseOperator::canApplyDuringMating "

Description:

    set if this operator can be applied during-mating

Usage:

    x.canApplyDuringMating()

"; 

%feature("docstring") simuPOP::baseOperator::canApplyPostMating "

Description:

    set if this operator can be applied post-mating

Usage:

    x.canApplyPostMating()

"; 

%feature("docstring") simuPOP::baseOperator::canApplyPreOrPostMating "

Description:

    set of this operator can be applied pre- or post-mating

Usage:

    x.canApplyPreOrPostMating()

"; 

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

%feature("docstring") simuPOP::baseOperator::MPIReady "

Description:

    determine if this operator can be used in a MPI module

Usage:

    x.MPIReady()

"; 

%ignore simuPOP::baseOperator::setHaploidOnly();

%ignore simuPOP::baseOperator::setDiploidOnly();

%ignore simuPOP::baseOperator::setMPIReady();

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

    apply to one  population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::baseOperator::applyDuringMating(population &pop, population::IndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::baseOperator::setOutput "

Description:

    set ouput stream, if was not set during construction

Usage:

    x.setOutput(output=\"\", outputExpr=\"\")

"; 

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

%feature("docstring") simuPOP::BernulliTrials "

Description:

    this class encapsulate behavior of a sequence of Bernulli trial.
    the main idea is that when doing a sequence of Bernulli trials of
    the same probability, we can use much quicker algorithms instead
    of doing n Bernulli trials

Details:

    For example, when N=10000, p=0.001. The usual way to do N Bin(p)
    trials is to do N randUnif(0,1)<p comparison.using the new method,
    we can use geometric distrubution to find the next true
    event.Also, for the cases of p=0.5, random bits are generated.This
    class maintain a two dimensional table: a vector of probabilities
    cross expected number of trialsp1 p2 p3 p4 p5 trial 1 trial 2 ...
    trial NWe expect that N is big (usually populaiton size) and p_i
    are smallusing fast BernulliTrial method for fix p, we can fill up
    this table very quickly column by columnThis class will provide
    easy access to row (each trial) or column (called each prob) of
    this table.if this table is accessed row by row (each trial), a
    internal index is used.if index exceeds N, trials will be
    generated all again. if trial will be called, e.g., N+2 times all
    the time, this treatment might not be very efficient.

Details:

    For example, when N=10000, p=0.001. The usual way to do N Bin(p)
    trials is to do N randUnif(0,1)<p comparison.using the new method,
    we can use geometric distrubution to find the next true
    event.Also, for the cases of p=0.5, random bits are generated.This
    class maintain a two dimensional table: a vector of probabilities
    cross expected number of trialsp1 p2 p3 p4 p5 trial 1 trial 2 ...
    trial NWe expect that N is big (usually populaiton size) and p_i
    are smallusing fast BernulliTrial method for fix p, we can fill up
    this table very quickly column by columnThis class will provide
    easy access to row (each trial) or column (called each prob) of
    this table.if this table is accessed row by row (each trial), a
    internal index is used.if index exceeds N, trials will be
    generated all again. if trial will be called, e.g., N+2 times all
    the time, this treatment might not be very efficient.

"; 

%ignore simuPOP::BernulliTrials::BernulliTrials(RNG &rng);

%feature("docstring") simuPOP::BernulliTrials::BernulliTrials "

Description:

    simuPOP::BernulliTrials::BernulliTrials

Usage:

    BernulliTrials(rng, prob, trials)

"; 

%feature("docstring") simuPOP::BernulliTrials::~BernulliTrials "

Description:

    simuPOP::BernulliTrials::~BernulliTrials

Usage:

    x.~BernulliTrials()

"; 

%ignore simuPOP::BernulliTrials::trialSize() const;

%feature("docstring") simuPOP::BernulliTrials::probSize "

Description:

    simuPOP::BernulliTrials::probSize

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

Description:

    simuPOP::BernulliTrials::trialSucc

Usage:

    x.trialSucc(idx)

"; 

%feature("docstring") simuPOP::BernulliTrials::probFirstSucc "

Description:

    simuPOP::BernulliTrials::probFirstSucc

Usage:

    x.probFirstSucc()

"; 

%feature("docstring") simuPOP::BernulliTrials::probNextSucc "

Description:

    simuPOP::BernulliTrials::probNextSucc

Usage:

    x.probNextSucc(pos)

"; 

%feature("docstring") simuPOP::BernulliTrials::trialFirstSucc "

Description:

    simuPOP::BernulliTrials::trialFirstSucc

Usage:

    x.trialFirstSucc(idx)

"; 

%feature("docstring") simuPOP::BernulliTrials::trialNextSucc "

Description:

    simuPOP::BernulliTrials::trialNextSucc

Usage:

    x.trialNextSucc(idx, pos)

"; 

%feature("docstring") simuPOP::BernulliTrials::setTrialSucc "

Description:

    simuPOP::BernulliTrials::setTrialSucc

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

Description:

    a  mating scheme that uses binomial selection, regardless of sex

Details:

    No sex information is involved (binomial random selection).
    Offspring is chosen from parental generation by random or
    according to the fitness values. In this  mating scheme,
    * numOffspring protocol is honored;
    *  population size changes are allowed;
    * selection is possible;
    * haploid populaton is allowed.

Details:

    No sex information is involved (binomial random selection).
    Offspring is chosen from parental generation by random or
    according to the fitness values. In this  mating scheme,
    * numOffspring protocol is honored;
    *  population size changes are allowed;
    * selection is possible;
    * haploid populaton is allowed.

"; 

%feature("docstring") simuPOP::binomialSelection::binomialSelection "

Description:

    create a binomial selection  mating scheme

Usage:

    binomialSelection(numOffspring=1., numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
      newSubPopSizeExpr=\"\", newSubPopSizeFunc=None)

Details:

    Please refer to class  mating for parameter descriptions.

Details:

    Please refer to class  mating for parameter descriptions.

"; 

%feature("docstring") simuPOP::binomialSelection::~binomialSelection "

Description:

    destructor

Usage:

    x.~binomialSelection()

"; 

%feature("docstring") simuPOP::binomialSelection::clone "

Description:

    deep copy of a binomial selection  mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::binomialSelection::__repr__ "

Description:

    used by Python print function to print out the general information
    of the binomial selection  mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::binomialSelection::submitScratch(population &pop, population &scratch);

%ignore simuPOP::binomialSelection::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%feature("docstring") simuPOP::caseControlSample "

Function form:

    CaseControlSample

Function form:

    CaseControlSample

Description:

    draw a case-control  sample from a  population

Details:

    This operator will randomly choose cases affected individuals and
    controls unaffected individuals as a  sample. The affected status
    is usually set by  penetrance functions/operators. The  sample
    populations will have two subpopulations: cases and controls.
    You may specify the number of cases and the number of controls
    from each subpopulation using the array form of the parameters.
    The  sample population will still have only two subpoulations
    (cases/controls) though.
    A special case of this sampling scheme occurs when one of or both
    cases and controls are omitted (zeros). In this case, all cases
    and/or controls are chosen. If both parameters are omitted, the
    sample is effectively the same  population with affected and
    unaffected individuals separated into two subpopulations.

Details:

    This operator will randomly choose cases affected individuals and
    controls unaffected individuals as a  sample. The affected status
    is usually set by  penetrance functions/operators. The  sample
    populations will have two subpopulations: cases and controls.
    You may specify the number of cases and the number of controls
    from each subpopulation using the array form of the parameters.
    The  sample population will still have only two subpoulations
    (cases/controls) though.
    A special case of this sampling scheme occurs when one of or both
    cases and controls are omitted (zeros). In this case, all cases
    and/or controls are chosen. If both parameters are omitted, the
    sample is effectively the same  population with affected and
    unaffected individuals separated into two subpopulations.

"; 

%feature("docstring") simuPOP::caseControlSample::caseControlSample "

Description:

    draw cases and controls as a  sample

Usage:

    caseControlSample(cases=[], controls=[], spSample=False,
      name=\"sample\", nameExpr=\"\", times=1, saveAs=\"\", saveAsExpr=\"\",
      format=\"auto\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Details:

    Please refer to class  sample for other parameter descriptions.

Details:

    Please refer to class  sample for other parameter descriptions.

Arguments:

    cases:          the number of cases, or an array of the numbers of
                    cases from each subpopulation
    controls:       the number of controls, or an array of the numbers
                    of controls from each subpopulation

"; 

%feature("docstring") simuPOP::caseControlSample::~caseControlSample "

Description:

    destructor

Usage:

    x.~caseControlSample()

"; 

%feature("docstring") simuPOP::caseControlSample::clone "

Description:

    deep copy of a  caseControlSample operator

Usage:

    x.clone()

"; 

%ignore simuPOP::caseControlSample::prepareSample(population &pop);

%ignore simuPOP::caseControlSample::drawsample(population &pop);

%feature("docstring") simuPOP::caseControlSample::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  caseControlSample operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::continueIf "

Description:

    terminate according to a condition failure

Details:

    The same as  terminateIf but continue if the condition is True.

Details:

    The same as  terminateIf but continue if the condition is True.

"; 

%feature("docstring") simuPOP::continueIf::continueIf "

Description:

    create a  continueIf terminator

Usage:

    continueIf(condition=\"\", message=\"\", var=\"terminate\", output=\"\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

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

    apply the  continueIf terminator???

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::continueIf::~continueIf "

Description:

    simuPOP::continueIf::~continueIf

Usage:

    x.~continueIf()

"; 

%feature("docstring") simuPOP::controlledBinomialSelection "

Description:

    a controlled binomial random selection  mating scheme

Details:

    This is the controlled binomial random selection  mating scheme
    described in  Peng 2007 (PLoS Genetics) . Basically, a freqFunc is
    passed to this  mating scheme and set the allele frequencies of
    given alleles at given loci at the offspring generation.
    The offspring generation is conceptually populated in two steps.
    At the first step, only families with disease alleles are accepted
    until the expected number of disease alleles are met. At the
    second step, only families with wide type alleles are accepted to
    populate the rest of the offspring generation.

Details:

    This is the controlled binomial random selection  mating scheme
    described in  Peng 2007 (PLoS Genetics) . Basically, a freqFunc is
    passed to this  mating scheme and set the allele frequencies of
    given alleles at given loci at the offspring generation.
    The offspring generation is conceptually populated in two steps.
    At the first step, only families with disease alleles are accepted
    until the expected number of disease alleles are met. At the
    second step, only families with wide type alleles are accepted to
    populate the rest of the offspring generation.

"; 

%feature("docstring") simuPOP::controlledBinomialSelection::controlledBinomialSelection "

Description:

    create a controlled binomial random selection  mating scheme

Usage:

    controlledBinomialSelection(loci, alleles, freqFunc,
      numOffspring=1., numOffspringFunc=None, maxNumOffspring=0,
      mode=MATE_NumOffspring, newSubPopSize=[], newSubPopSizeExpr=\"\",
      newSubPopSizeFunc=None)

Details:

    Please refer to class  mating for descriptions of parameters.

Details:

    Please refer to class  mating for descriptions of parameters.

"; 

%ignore simuPOP::controlledBinomialSelection::controlledBinomialSelection(const controlledBinomialSelection &rhs);

%feature("docstring") simuPOP::controlledBinomialSelection::~controlledBinomialSelection "

Description:

    destructor

Usage:

    x.~controlledBinomialSelection()

"; 

%feature("docstring") simuPOP::controlledBinomialSelection::clone "

Description:

    deep copy of a controlled binomial random selection  mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::controlledBinomialSelection::__repr__ "

Description:

    used by Python print function to print out the general information
    of the controlled binomial random selection  mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::controlledBinomialSelection::submitScratch(population &pop, population &scratch);

%ignore simuPOP::controlledBinomialSelection::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%feature("docstring") simuPOP::controlledMating "

Description:

    a controlled  mating scheme

Details:

    This is an experimental  mating scheme that uses a frequency range
    to control the allele frequency of the offspring generation at
    given loci. When allele frequencies at the offspring generation
    does not fall into the given range, the offspring generation is
    regenerated. Any  mating scheme can be used with this  mating
    scheme by passing through parameter matingScheme.

Details:

    This is an experimental  mating scheme that uses a frequency range
    to control the allele frequency of the offspring generation at
    given loci. When allele frequencies at the offspring generation
    does not fall into the given range, the offspring generation is
    regenerated. Any  mating scheme can be used with this  mating
    scheme by passing through parameter matingScheme.

"; 

%feature("docstring") simuPOP::controlledMating::controlledMating "

Description:

    control allele frequencies at a locus

Usage:

    controlledMating(matingScheme, loci, alleles, freqFunc,
      range=0.01)

Arguments:

    matingScheme:   a  mating scheme
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

%feature("docstring") simuPOP::controlledMating::clone "

Description:

    deep copy of a controlled  mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::controlledMating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::controlledMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the controlled  mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::controlledMating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%feature("docstring") simuPOP::controlledRandomMating "

Description:

    a controlled random  mating scheme

Details:

    This is the controlled random  mating scheme described in  Peng
    2007 (PLoS Genetics) . Basically, a freqFunc is passed to this
    mating scheme and set the allele frequencies of given alleles at
    given loci at the offspring generation.
    The offspring generation is conceptually populated in two steps.
    At the first step, only families with disease alleles are accepted
    until the expected number of disease alleles are met. At the
    second step, only families with wide type alleles are accepted to
    populate the rest of the offspring generation.

Details:

    This is the controlled random  mating scheme described in  Peng
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

    create a controlled random  mating scheme

Usage:

    controlledRandomMating(loci, alleles, freqFunc, acceptScheme=0,
      numOffspring=1., numOffspringFunc=None, maxNumOffspring=0,
      mode=MATE_NumOffspring, newSubPopSize=[],
      newSubPopSizeFunc=None, newSubPopSizeExpr=\"\",
      contWhenUniSex=True)

Details:

    Please refer to class  mating for descriptions of other
    parameters.

Details:

    Please refer to class  mating for descriptions of other
    parameters.

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

    deep copy of a controlled random  mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::controlledRandomMating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::controlledRandomMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the controlled random  mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::controlledRandomMating::submitScratch(population &pop, population &scratch);

%ignore simuPOP::controlledRandomMating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%feature("docstring") simuPOP::dumper "

Description:

    dump the content of a  population.

"; 

%feature("docstring") simuPOP::dumper::dumper "

Description:

    dump  population

Usage:

    dumper(alleleOnly=False, infoOnly=False, ancestralPops=False,
      dispWidth=1, max=100, chrom=[], loci=[], subPop=[], indRange=[],
      output=\">\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    alleleOnly:     only display allele
    infoOnly:       only display info
    dispWidth:      width of allele display, default to 1
    ancestralPops:  whether or not display ancestral populations,
                    default to False
    chrom:          chromsoome(s) to display
    loci:           loci to display
    subPop:         only display subPop(s)
    indRange:       range(s) of individuals to display
    max:            max number of individuals to display, default to
                    100. This is to avoid careless dump of huge
                    populations.
    output:         output file, default to standard output.
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

Description:

    simuPOP::dumper::setAlleleOnly

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

Description:

    simuPOP::dumper::setInfoOnly

Usage:

    x.setInfoOnly(infoOnly)

"; 

%feature("docstring") simuPOP::dumper::apply "

Description:

    apply to one  population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::dumper::~dumper "

Description:

    simuPOP::dumper::~dumper

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

Description:

    simuPOP::Exception::~Exception

Usage:

    x.~Exception()

"; 

%ignore simuPOP::Expression;

%feature("docstring") simuPOP::Expression::Expression "

Description:

    simuPOP::Expression::Expression

Usage:

    Expression(expr=\"\", stmts=\"\", locals=None)

"; 

%feature("docstring") simuPOP::Expression::~Expression "

Description:

    simuPOP::Expression::~Expression

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

%ignore simuPOP::GappedIterator;

%feature("docstring") simuPOP::GappedIterator::GappedIterator "

Description:

    simuPOP::GappedIterator::GappedIterator

Usage:

    GappedIterator()

"; 

%ignore simuPOP::GappedIterator::GappedIterator(pointer p, difference_type s=1);

%feature("docstring") simuPOP::GappedIterator::~GappedIterator "

Description:

    simuPOP::GappedIterator::~GappedIterator

Usage:

    x.~GappedIterator()

"; 

%feature("docstring") simuPOP::GappedIterator::ptr "

Description:

    get pointer so that  ptr()+1 will be the next locus.

Usage:

    x.ptr()

"; 

%ignore simuPOP::GappedIterator::step();

%ignore simuPOP::GenoStructure;

%ignore simuPOP::GenoStructure::GenoStructure();

%ignore simuPOP::GenoStructure::GenoStructure(UINT ploidy, const vectoru &loci, bool sexChrom, const vectorf &lociPos, const vectorstr &alleleNames, const vectorstr &lociNames, UINT maxAllele, const vectorstr &infoFields, const vectori &chromMap);

%feature("docstring") simuPOP::GenoStructure::~GenoStructure "

Description:

    destructor, do nothing.

Usage:

    x.~GenoStructure()

"; 

%ignore simuPOP::GenoStructure::locusPos(UINT locus) const;

%ignore simuPOP::GenoStructure::chromIndex(UINT ch) const ;

%feature("docstring") simuPOP::GenoStruTrait "

Description:

    genotypic structure related functions, can be accessed from both
    individuals and populations

Details:

    Genotypic structure refers to the number of chromosomes,
    positions, the number of loci on each chromosome, and allele and
    locus names etc. All individuals in a  population share the same
    genotypic structure. Because class  GenoStruTrait is inherited by
    class  population, class  individual, and class  simulator,
    functions provided in this class can be accessed at the
    individual,  population and  simulator levels.

Details:

    Genotypic structure refers to the number of chromosomes,
    positions, the number of loci on each chromosome, and allele and
    locus names etc. All individuals in a  population share the same
    genotypic structure. Because class  GenoStruTrait is inherited by
    class  population, class  individual, and class  simulator,
    functions provided in this class can be accessed at the
    individual,  population and  simulator levels.

"; 

%feature("docstring") simuPOP::GenoStruTrait::GenoStruTrait "

Description:

    This object can not be created directly. It is created by a
    population.

Usage:

    GenoStruTrait()

Example:

>>> # create a population, most parameters have default values
>>> pop = population(size=5, ploidy=2, loci=[5,10],
...     lociPos=[range(0,5),range(0,20,2)],
...     alleleNames=['A','C','T','G'],
...     subPop=[2,3], maxAllele=3)
>>> print pop.popSize()
5
>>> print pop.ploidy()
2
>>> print pop.ploidyName()
diploid
>>> print pop.numChrom()
2
>>> print pop.locusPos(2)
2.0
>>> print pop.alleleName(1)
C
>>> # get the fourth individual of the population
>>> ind = pop.individual(3)
>>> # access genotypic structure info
>>> print ind.ploidy()
2
>>> print ind.numChrom()
2
>>> print ind.numLoci(0)
5
>>> print ind.genoSize()
30
>>> # and from simulator level
>>> simu = simulator(pop, randomMating(), rep=3)
>>> print simu.numChrom()
2
>>>


"; 

%ignore simuPOP::GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru &loci, bool sexChrom, const vectorf &lociPos, const vectorstr &alleleNames, const vectorstr &lociNames, UINT maxAllele, const vectorstr &infoFields, const vectori &chromMap);

%ignore simuPOP::GenoStruTrait::setGenoStructure(GenoStructure &rhs);

%ignore simuPOP::GenoStruTrait::setGenoStruIdx(size_t idx);

%ignore simuPOP::GenoStruTrait::mergeGenoStru(size_t idx, bool byChromosome) const ;

%ignore simuPOP::GenoStruTrait::removeLociFromGenoStru(const vectoru &remove=vectoru(), const vectoru &keep=vectoru());

%ignore simuPOP::GenoStruTrait::insertBeforeLociToGenoStru(const vectoru &idx, const vectorf &pos, const vectorstr &names) const ;

%ignore simuPOP::GenoStruTrait::insertAfterLociToGenoStru(const vectoru &idx, const vectorf &pos, const vectorstr &names) const ;

%ignore simuPOP::GenoStruTrait::genoStru() const;

%ignore simuPOP::GenoStruTrait::genoStruIdx() const;

%feature("docstring") simuPOP::GenoStruTrait::ploidy "

Description:

    return ploidy, the number of homologous sets of chromosomes

Usage:

    x.ploidy()

"; 

%feature("docstring") simuPOP::GenoStruTrait::ploidyName "

Description:

    return ploidy name, haploid, diploid, or triploid etc.

Usage:

    x.ploidyName()

"; 

%feature("docstring") simuPOP::GenoStruTrait::numLoci "

Description:

    return the number of loci on chromosome chrom, equivalent to
    numLoci()[chrom]

Usage:

    x.numLoci(chrom)

"; 

%feature("docstring") simuPOP::GenoStruTrait::numLoci "

Description:

    return the number of loci on all chromosomes

Usage:

    x.numLoci()

"; 

%feature("docstring") simuPOP::GenoStruTrait::sexChrom "

Description:

    determine whether or not the last chromosome is sex chromosome

Usage:

    x.sexChrom()

"; 

%feature("docstring") simuPOP::GenoStruTrait::totNumLoci "

Description:

    return the total number of loci on all chromosomes

Usage:

    x.totNumLoci()

"; 

%feature("docstring") simuPOP::GenoStruTrait::genoSize "

Description:

    return the total number of loci times ploidy

Usage:

    x.genoSize()

"; 

%feature("docstring") simuPOP::GenoStruTrait::locusPos "

Description:

    return the position of a locus

Usage:

    x.locusPos(locus)

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociPos "

Description:

    return loci positions

Usage:

    x.lociPos()

"; 

%feature("docstring") simuPOP::GenoStruTrait::arrLociPos "

Description:

    return an (editable) array of loci positions of all loci

Usage:

    x.arrLociPos()

Note:

    Modifying loci position directly using this function is strongly
    discouraged.

"; 

%feature("docstring") simuPOP::GenoStruTrait::arrLociPos "

Description:

    return an (editable) array of loci positions on a given chromosome

Usage:

    x.arrLociPos(chrom)

Note:

    Modifying loci position directly using this function is strongly
    discouraged.

"; 

%feature("docstring") simuPOP::GenoStruTrait::numChrom "

Description:

    return the number of chromosomes

Usage:

    x.numChrom()

"; 

%ignore simuPOP::GenoStruTrait::chromIndex() const;

%feature("docstring") simuPOP::GenoStruTrait::chromBegin "

Description:

    return the index of the first locus on a chromosome

Usage:

    x.chromBegin(chrom)

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromEnd "

Description:

    return the index of the last locus on a chromosome plus 1

Usage:

    x.chromEnd(chrom)

Note:

    From the description of this function, the returned value may not
    be a valid index. (This is consistant with Python ranges.)

"; 

%feature("docstring") simuPOP::GenoStruTrait::absLocusIndex "

Description:

    return the absolute index of a locus on a chromosome

Usage:

    x.absLocusIndex(chrom, locus)

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromLocusPair "

Description:

    return a (chrom, locus) pair of an absolute locus index

Usage:

    x.chromLocusPair(locus)

"; 

%feature("docstring") simuPOP::GenoStruTrait::alleleName "

Description:

    return the name of an allele (if previously specified). Default to
    allele index.

Usage:

    x.alleleName(allele)

"; 

%feature("docstring") simuPOP::GenoStruTrait::alleleNames "

Description:

    return an array of allele names

Usage:

    x.alleleNames()

"; 

%feature("docstring") simuPOP::GenoStruTrait::locusName "

Description:

    return the name of a locus

Usage:

    x.locusName(loc)

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociNames "

Description:

    return names of all loci

Usage:

    x.lociNames()

"; 

%feature("docstring") simuPOP::GenoStruTrait::locusByName "

Description:

    return the index of a locus by its locus name

Usage:

    x.locusByName(name)

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociByNames "

Description:

    return an array of locus indices by locus names

Usage:

    x.lociByNames(names)

"; 

%feature("docstring") simuPOP::GenoStruTrait::maxAllele "

Description:

    return the maximum allele value for all loci. Default to maximum
    allowed allele state.

Usage:

    x.maxAllele()

Details:

    Maximum allele value has to be 1 for binary modules. maxAllele is
    the maximum possible allele value, which allows maxAllele+1
    alleles 0, 1, ..., maxAllele.

Details:

    Maximum allele value has to be 1 for binary modules. maxAllele is
    the maximum possible allele value, which allows maxAllele+1
    alleles 0, 1, ..., maxAllele.

"; 

%ignore simuPOP::GenoStruTrait::setMaxAllele(UINT maxAllele);

%feature("docstring") simuPOP::GenoStruTrait::hasInfoField "

Description:

    determine if an information field exists

Usage:

    x.hasInfoField(name)

"; 

%feature("docstring") simuPOP::GenoStruTrait::infoSize "

Description:

    obtain the number of information fields

Usage:

    x.infoSize()

"; 

%feature("docstring") simuPOP::GenoStruTrait::infoFields "

Description:

    return an array of all information fields

Usage:

    x.infoFields()

"; 

%feature("docstring") simuPOP::GenoStruTrait::infoField "

Description:

    obtain the name of information field idx

Usage:

    x.infoField(idx)

"; 

%feature("docstring") simuPOP::GenoStruTrait::infoIdx "

Description:

    return the index of the field name, return -1 if not found

Usage:

    x.infoIdx(name)

"; 

%ignore simuPOP::GenoStruTrait::struAddInfoFields(const vectorstr &fields);

%ignore simuPOP::GenoStruTrait::struSetInfoFields(const vectorstr &fields);

%ignore simuPOP::GenoStruTrait::swap(GenoStruTrait &rhs);

%ignore simuPOP::GenoStruTrait::chromMap() const;

%feature("docstring") simuPOP::gsmMutator "

Function form:

    GsmMutate

Function form:

    GsmMutate

Description:

    generalized stepwise mutation model

Details:

    Generalized Stepwise Mutation model (GSM) is an extension to
    stepwise mutation model. This model assumes that alleles are
    represented by integer values and that a mutation either increases
    or decreases the allele value by a random value. In other words,
    in this model the change in the allelic state is drawn from a
    random distribution. A geometric generalized stepwise model uses a
    geometric distribution with parameter  $ p $, which has mean  $
    \\frac{p}{1-p} $ and variance  $ \\frac{p}{\\left(1-p\\right)^{2}} $.
    gsmMutator implements both models. If you specify a Python
    function without a parameter, this  mutator will use its return
    value each time a mutation occur; otherwise, a parameter  $ p $
    should be provided and the  mutator will act as a geometric
    generalized stepwise model.

Details:

    Generalized Stepwise Mutation model (GSM) is an extension to
    stepwise mutation model. This model assumes that alleles are
    represented by integer values and that a mutation either increases
    or decreases the allele value by a random value. In other words,
    in this model the change in the allelic state is drawn from a
    random distribution. A geometric generalized stepwise model uses a
    geometric distribution with parameter  $ p $, which has mean  $
    \\frac{p}{1-p} $ and variance  $ \\frac{p}{\\left(1-p\\right)^{2}} $.
    gsmMutator implements both models. If you specify a Python
    function without a parameter, this  mutator will use its return
    value each time a mutation occur; otherwise, a parameter  $ p $
    should be provided and the  mutator will act as a geometric
    generalized stepwise model.

Example:

>>> simu = simulator(population(size=3, loci=[3,5]), noMating())
>>> simu.step([
...     initByFreq( [.2,.3,.5]),
...     gsmMutator(rate=1, p=.8, incProb=.8),
...     dumper(alleleOnly=True, stage=PrePostMating)])
individual info: 
sub population 0:
   0: MU   0  2  1   0  2  1  2  1 |   2  2  1   2  0  2  1  2 
   1: FU   2  2  2   0  1  1  0  2 |   2  2  2   0  1  2  1  1 
   2: FU   2  1  1   2  2  2  1  1 |   2  2  0   0  0  2  2  2 
End of individual info.


No ancenstral population recorded.
individual info: 
sub population 0:
   0: MU   2  3  3   2  3  2  3  2 |   3  3  0   3  1  3  2  1 
   1: FU   1  3  1   1  2  0  2  3 |   4  0  1   1  2  3  3  2 
   2: FU   3  0  2   5  3  0  3  4 |   3  5  1   0  2  3  3  4 
End of individual info.


No ancenstral population recorded.
True
>>> 
>>> import random
>>> def rndInt():
...   return random.randrange(3,6)
... 
>>> simu.step([
...     initByFreq( [.2,.3,.5]),
...     gsmMutator(rate=1, func=rndInt, incProb=.8),
...     dumper(alleleOnly=True, stage=PrePostMating)])
individual info: 
sub population 0:
   0: FU   1  1  2   2  1  0  2  0 |   2  1  1   1  2  2  2  2 
   1: FU   2  1  2   1  2  2  2  2 |   1  2  1   1  1  2  1  2 
   2: FU   2  2  2   1  2  2  2  1 |   2  2  2   2  2  2  2  0 
End of individual info.


No ancenstral population recorded.
individual info: 
sub population 0:
   0: FU   0  0  7   7  0  5  7  5 |   0  4  6   6  0  0  6  7 
   1: FU   0  4  0   4  0  7  0  7 |   4  6  4   5  6  6  4  7 
   2: FU   7  5  5   6  7  0  5  6 |   7  0  5   5  7  5  0  5 
End of individual info.


No ancenstral population recorded.
True
>>> 
>>>


"; 

%feature("docstring") simuPOP::gsmMutator::gsmMutator "

Description:

    create a  gsmMutator

Usage:

    gsmMutator(rate=[], loci=[], maxAllele=0, incProb=0.5, p=0,
      func=None, output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
      end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Details:

    The generalized stepwise mutation model (GMM) is developed for
    allozymes. It provides better description for these kinds of
    evolutionary processes. Please see  mutator for the description of
    other parameters.

Details:

    The generalized stepwise mutation model (GMM) is developed for
    allozymes. It provides better description for these kinds of
    evolutionary processes. Please see  mutator for the description of
    other parameters.

Arguments:

    incProb:        probability to increase allele state. Default to
                    0.5.
    func:           return number of steps. No parameter.???

"; 

%feature("docstring") simuPOP::gsmMutator::~gsmMutator "

Description:

    simuPOP::gsmMutator::~gsmMutator

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

%feature("docstring") simuPOP::ifElse "

Description:

    conditional operator

Details:

    This operator accepts
    * an expression that will be evaluated when this operator is
    applied;
    * an operator that will be applied if the expression is True
    (default to null);
    * an operator that will be applied if the expression is False
    (default to null). When this operator is applied to a  population,
    it will evaluate the expression and depending on its value, apply
    the supplied operator. Note that the begin, at, step, and at
    parameters of ifOp and elseOp will be ignored. For example, you
    can mimic the at parameter of an operator by  ifElse('rep in
    [2,5,9]' operator). The real use of this machanism is to monitor
    the  population statistics and act accordingly.

Details:

    This operator accepts
    * an expression that will be evaluated when this operator is
    applied;
    * an operator that will be applied if the expression is True
    (default to null);
    * an operator that will be applied if the expression is False
    (default to null). When this operator is applied to a  population,
    it will evaluate the expression and depending on its value, apply
    the supplied operator. Note that the begin, at, step, and at
    parameters of ifOp and elseOp will be ignored. For example, you
    can mimic the at parameter of an operator by  ifElse('rep in
    [2,5,9]' operator). The real use of this machanism is to monitor
    the  population statistics and act accordingly.

"; 

%feature("docstring") simuPOP::ifElse::ifElse "

Description:

    simuPOP::ifElse::ifElse

Usage:

    ifElse(cond, ifOp=None, elseOp=None, output=\">\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    cond:           expression that will be treated as a bool variable
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

%feature("docstring") simuPOP::ifElse::clone "

Description:

    deep copy of an  ifElse operator

Usage:

    x.clone()

"; 

%ignore simuPOP::ifElse::applyWithScratch(population &pop, population &scratch, int stage);

%ignore simuPOP::ifElse::applyDuringMating(population &pop, population::IndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::ifElse::apply "

Description:

    apply the  ifElse operator to one  population

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

Description:

    simuPOP::IndexError::IndexError

Usage:

    IndexError(msg)

"; 

%feature("docstring") simuPOP::individual "

Description:

    individuals with genotype, affection status, sex etc.

Details:

    Individuals are the building blocks of populations, each having
    the following  individual information:
    * shared genotypic structure information
    * genotype
    * sex, affection status, subpopulation ID
    * optional information fields Individual genotypes are arranged by
    locus, chromosome, ploidy, in that order, and can be accessed from
    a single index. For example, for a diploid  individual with two
    loci on the first chromosome, one locus on the second, its
    genotype is arranged as  1-1-1 1-1-2 1-2-1 2-1-1 2-1-2 2-2-1
    where x-y-z represents ploidy x chromosome y and locus z. An
    allele 2-1-2 can be accessed by allele(4) (by absolute index),
    allele(1, 1) (by index and ploidy) or allele(1, 1, 0) (by index,
    ploidy and chromosome).

Details:

    Individuals are the building blocks of populations, each having
    the following  individual information:
    * shared genotypic structure information
    * genotype
    * sex, affection status, subpopulation ID
    * optional information fields Individual genotypes are arranged by
    locus, chromosome, ploidy, in that order, and can be accessed from
    a single index. For example, for a diploid  individual with two
    loci on the first chromosome, one locus on the second, its
    genotype is arranged as  1-1-1 1-1-2 1-2-1 2-1-1 2-1-2 2-2-1
    where x-y-z represents ploidy x chromosome y and locus z. An
    allele 2-1-2 can be accessed by allele(4) (by absolute index),
    allele(1, 1) (by index and ploidy) or allele(1, 1, 0) (by index,
    ploidy and chromosome).

"; 

%feature("docstring") simuPOP::individual::individual "

Description:

    Individuals are created by populations automatically. Do not call
    this function directly.

Usage:

    individual()

Example:

>>> pop = population(500, loci=[2, 5, 10])
>>> # get an individual
>>> ind = pop.individual(9)
>>> # oops, wrong index
>>> ind = pop.individual(3)
>>> # you can access genotypic structure info
>>> print ind.ploidy()
2
>>> print ind.numChrom()
3
>>> # ...
>>> # as well as genotype
>>> print ind.allele(1) 
0
>>> ind.setAllele(1,5)
>>> print ind.allele(1)
0
>>> # you can also use an overloaded function
>>> # with a second parameter being the ploidy index
>>> print ind.allele(1,1) # second locus at the second copy of chromosome
0
>>> # other information
>>> print ind.affected()
False
>>> print ind.affectedChar()
U
>>> ind.setAffected(1)
>>> print ind.affectedChar()
A
>>> print ind.sexChar()
M
>>>


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

%feature("docstring") simuPOP::individual::arrGenotype "

Description:

    return an editable array (a carray of length
    totNumLoci()*ploidy()) of genotypes of an  individual

Usage:

    x.arrGenotype()

Details:

    This function returns the whole genotype. Although this function
    is not as easy to use as other functions that access alleles, it
    is the fastest one since you can read/write genotype directly.

Details:

    This function returns the whole genotype. Although this function
    is not as easy to use as other functions that access alleles, it
    is the fastest one since you can read/write genotype directly.

"; 

%feature("docstring") simuPOP::individual::arrGenotype "

Description:

    return a carray with the genotype of the p-th copy of the
    chromosomes

Usage:

    x.arrGenotype(p)

"; 

%feature("docstring") simuPOP::individual::arrGenotype "

Description:

    return a carray with the genotype of the ch-th chromosome of the
    p-th copy

Usage:

    x.arrGenotype(p, ch)

"; 

%feature("docstring") simuPOP::individual::arrInfo "

Description:

    return a carray of all information fields (of size infosSize()) of
    this  individual

Usage:

    x.arrInfo()

"; 

%feature("docstring") simuPOP::individual::allele "

Description:

    return the allele at locus index

Usage:

    x.allele(index)

Arguments:

    index:          absolute index from the beginning of the genotype,
                    ranging from 0 to  totNumLoci()*ploidy()

"; 

%feature("docstring") simuPOP::individual::allele "

Description:

    return the allele at locus index of the p-th copy of the
    chromosomes

Usage:

    x.allele(index, p)

Arguments:

    index:          index from the begining of the p-th set of the
                    chromosomes, ranging from 0 to  totNumLoci()
    p:              index of the ploidy

"; 

%feature("docstring") simuPOP::individual::allele "

Description:

    return the allele at locus index of the ch-th chromosome of the
    p-th chromosome set

Usage:

    x.allele(index, p, ch)

Arguments:

    index:          index from the begining of chromosome ch of ploidy
                    p, ranging from 0 to  numLoci(ch)
    p:              index of the polidy
    ch:             index of the chromosome in the p-th chromosome set

"; 

%feature("docstring") simuPOP::individual::alleleChar "

Description:

    return the name of allele(index)

Usage:

    x.alleleChar(index)

"; 

%feature("docstring") simuPOP::individual::alleleChar "

Description:

    return the name of allele(index, p)

Usage:

    x.alleleChar(index, p)

"; 

%feature("docstring") simuPOP::individual::alleleChar "

Description:

    return the name of allele(idx, p, ch)

Usage:

    x.alleleChar(index, p, ch)

"; 

%feature("docstring") simuPOP::individual::setAllele "

Description:

    set the allele at locus index

Usage:

    x.setAllele(allele, index)

Arguments:

    allele:         allele to be set
    index:          index from the begining of genotype, ranging from
                    0 to  totNumLoci()*ploidy()

"; 

%feature("docstring") simuPOP::individual::setAllele "

Description:

    set the allele at locus index of the p-th copy of the chromosomes

Usage:

    x.setAllele(allele, index, p)

Arguments:

    allele:         allele to be set
    index:          index from the begining of the poloidy p, ranging
                    from 0 to  totNumLoci(p)
    p:              index of the poloidy

"; 

%feature("docstring") simuPOP::individual::setAllele "

Description:

    set the allele at locus index of the ch-th chromosome in the p-th
    chromosome set

Usage:

    x.setAllele(allele, index, p, ch)

Arguments:

    allele:         allele to be set
    index:          index from the begining of the chromosome, ranging
                    from 0 to numLoci(ch)
    p:              index of the ploidy
    ch:             index of the chromosome in ploidy p

"; 

%feature("docstring") simuPOP::individual::sex "

Description:

    return the sex of an  individual, 1 for males and 2 for females.

Usage:

    x.sex()

"; 

%feature("docstring") simuPOP::individual::sexChar "

Description:

    return the sex of an  individual, M or F

Usage:

    x.sexChar()

"; 

%feature("docstring") simuPOP::individual::setSex "

Description:

    set sex. sex can be Male of Female.

Usage:

    x.setSex(sex)

"; 

%feature("docstring") simuPOP::individual::affected "

Description:

    whether or not an  individual is affected

Usage:

    x.affected()

"; 

%feature("docstring") simuPOP::individual::unaffected "

Description:

    equals to not  affected()

Usage:

    x.unaffected()

"; 

%feature("docstring") simuPOP::individual::affectedChar "

Description:

    return A or U for affection status

Usage:

    x.affectedChar()

"; 

%feature("docstring") simuPOP::individual::setAffected "

Description:

    set affection status

Usage:

    x.setAffected(affected)

"; 

%feature("docstring") simuPOP::individual::subPopID "

Description:

    return the ID of the subpopulation to which this  individual
    blongs

Usage:

    x.subPopID()

Note:

    subPopID is not set by default. It only corresponds to the
    subpopulation in which this  individual resides after
    pop::setIndSubPopID is called.

"; 

%feature("docstring") simuPOP::individual::setSubPopID "

Description:

    set new subpopulation ID, pop.rearrangeByIndID will move this
    individual to that  population

Usage:

    x.setSubPopID(id)

"; 

%feature("docstring") simuPOP::individual::info "

Description:

    get information field idx

Usage:

    x.info(idx)

Arguments:

    idx:            index of the information field

"; 

%feature("docstring") simuPOP::individual::info "

Description:

    get information field name

Usage:

    x.info(name)

Details:

    Equivalent to info(infoIdx(name)).

Details:

    Equivalent to info(infoIdx(name)).

Arguments:

    name:           name of the information field

"; 

%feature("docstring") simuPOP::individual::setInfo "

Description:

    set information field by idx

Usage:

    x.setInfo(value, idx)

"; 

%feature("docstring") simuPOP::individual::setInfo "

Description:

    set information field by name

Usage:

    x.setInfo(value, name)

"; 

%ignore simuPOP::individual::genoBegin() const;

%ignore simuPOP::individual::genoEnd() const;

%ignore simuPOP::individual::genoBegin(UINT p) const ;

%ignore simuPOP::individual::genoEnd(UINT p) const ;

%ignore simuPOP::individual::genoBegin(UINT p, UINT chrom) const ;

%ignore simuPOP::individual::infoBegin() const;

%ignore simuPOP::individual::infoEnd() const;

%feature("docstring") simuPOP::individual::__cmp__ "

Description:

    a python function used to compare the  individual objects

Usage:

    x.__cmp__(rhs)

"; 

%feature("docstring") simuPOP::individual::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  individual

Usage:

    x.__repr__()

"; 

%ignore simuPOP::individual::swap(individual &ind, bool swapContent=true);

%ignore simuPOP::individual::shallowCopied() const;

%ignore simuPOP::individual::setShallowCopied(bool shallowCopied);

%ignore simuPOP::individual::display(ostream &out, int width=1, const vectori &chrom=vectori(), const vectori &loci=vectori());

%feature("docstring") simuPOP::individualIterator "

Description:

    this class implements a Python itertor class that can be used to
    iterate through individuals in a  population.

Details:

    an instance of this class is returned by
    population::individuals() and population::individuals(subPop)

Details:

    an instance of this class is returned by
    population::individuals() and population::individuals(subPop)

"; 

%feature("docstring") simuPOP::individualIterator::individualIterator "

Description:

    simuPOP::individualIterator::individualIterator

Usage:

    individualIterator(pop, s, e)

"; 

%feature("docstring") simuPOP::individualIterator::~individualIterator "

Description:

    simuPOP::individualIterator::~individualIterator

Usage:

    x.~individualIterator()

"; 

%feature("docstring") simuPOP::individualIterator::__iter__ "

Description:

    simuPOP::individualIterator::__iter__

Usage:

    x.__iter__()

"; 

%feature("docstring") simuPOP::individualIterator::next "

Description:

    simuPOP::individualIterator::next

Usage:

    x.next()

"; 

%feature("docstring") simuPOP::inheritTagger "

Description:

    inherite tag from parents.

Details:

    This during-mating operator will copy the tag information from
    his/her parents. Depending on mode parameter, this  tagger will
    obtain tag from his/her father (two tag fields), mother (two tag
    fields) or both (first tag field from both father and mother).An
    example may be tagging one or a few parents and see, at the last
    generation, how many offspring they have.

Details:

    This during-mating operator will copy the tag information from
    his/her parents. Depending on mode parameter, this  tagger will
    obtain tag from his/her father (two tag fields), mother (two tag
    fields) or both (first tag field from both father and mother).An
    example may be tagging one or a few parents and see, at the last
    generation, how many offspring they have.

"; 

%feature("docstring") simuPOP::inheritTagger::inheritTagger "

Description:

    create an  inheritTagger, default to be always active

Usage:

    inheritTagger(mode=TAG_Paternal, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[\"paternal_tag\",
      \"maternal_tag\"])

Arguments:

    mode:           can be one of TAG_Paternal, TAG_Maternal, and
                    TAG_Both

"; 

%feature("docstring") simuPOP::inheritTagger::~inheritTagger "

Description:

    simuPOP::inheritTagger::~inheritTagger

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

Function form:

    InitByFreq

Description:

    initialize genotypes by given allele frequencies, and sex by male
    frequency

Details:

    This operator accepts alleleFreq or alleleFreqs.??? The first one
    ignores subpopulation structure while the second one gives
    different initial allele frequencies to different subpop or
    ranges. Allele frequencies can differ in subpop. Sex is also
    assigned randomly.

Details:

    This operator accepts alleleFreq or alleleFreqs.??? The first one
    ignores subpopulation structure while the second one gives
    different initial allele frequencies to different subpop or
    ranges. Allele frequencies can differ in subpop. Sex is also
    assigned randomly.

"; 

%feature("docstring") simuPOP::initByFreq::initByFreq "

Description:

    randomly assign alleles according to given allele frequencies

Usage:

    initByFreq(alleleFreq=[], identicalInds=False, subPop=[],
      indRange=[], loci=[], atPloidy=-1, maleFreq=0.5, sex=[],
      stage=PreMating, begin=0, end=1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    alleleFreq:     an array of allele frequencies. The sum of all the
                    frequencies must be 1; or for a matrix of allele
                    frequencies, each row corresponses to a
                    subpopulation.
    subPop:         an array specifies applicable subpopulations
    indRange:       a [begin, end] pair of the range of absolute
                    indices of individuals, for example, ([1,2]); or
                    an array of [begin, end] pairs, such as
                    ([[1,4],[5,6]]). This is how you can initialize
                    individuals differently within subpopulations.
                    Note that ranges are in the form of [a,b). I.e.,
                    range [4,6] will intialize  individual 4, 5, but
                    not 6. As a shortcut for [4,5], you can use [4] to
                    specify one  individual.
    identicalInds:  whether or not make  individual genotypies
                    identical in all subpopulation. If True, this
                    operator will randomly generate genotype for an
                    individual and  spread it to the whole
                    subpopulation in the given range.
    loci:           a vector of locus indices at which initialization
                    will be done. If empty, apply to all loci.
    locus:          a shortcut to loci
    atPloidy:       initialize which copy of chromosomes. Default to
                    all.
    maleFreq:       male frequency. Default to 0.5. Sex will be
                    initialized with this parameter.
    sex:            an array of sex [Male, Female, Male...] for
                    individuals. The length of sex will not be
                    checked. If it is shorter than the number of
                    individuals, sex will be reused from the
                    beginning.
    stage:          default to PreMating

Example:

>>> simu = simulator( 
...     population(subPop=[2,3], loci=[5,7], maxAllele=1),
...     randomMating(), rep=1)
>>> simu.step([
...     initByFreq(alleleFreq=[ [.2,.8],[.8,.2]]),
...     dumper(alleleOnly=True)
...   ])
individual info: 
sub population 0:
   0: MU 11111 1011101 | 10111 1111011 
   1: FU 10111 1111011 | 11111 1011101 
sub population 1:
   2: FU 00000 0000010 | 10010 0100010 
   3: MU 00000 0000010 | 01100 0000000 
   4: FU 00000 0000100 | 01000 0100000 
End of individual info.


No ancenstral population recorded.
True
>>>


"; 

%feature("docstring") simuPOP::initByFreq::~initByFreq "

Description:

    simuPOP::initByFreq::~initByFreq

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

    apply operator  initByFreq???

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::initByValue "

Function form:

    InitByValue

Function form:

    InitByValue

Description:

    initialize genotype by value and then copy to all individuals

Details:

    initByValue operator gets one copy of chromosomes or the whole
    genotype (or of those corresponds to loci) of an  individual and
    copy them to all or a subset of individuals. This operator assign
    given alleles to specified individuals. Every  individual will
    have the same genotype. The parameter combinations should be
    * value - subPop/indRange:  individual in subPop or in range(s)
    will be assigned genotype 'value';
    * subPop/indRange: subPop or indRange should have the same length
    as values. Each item of values will be assigned to each subPop or
    indRange.

Details:

    initByValue operator gets one copy of chromosomes or the whole
    genotype (or of those corresponds to loci) of an  individual and
    copy them to all or a subset of individuals. This operator assign
    given alleles to specified individuals. Every  individual will
    have the same genotype. The parameter combinations should be
    * value - subPop/indRange:  individual in subPop or in range(s)
    will be assigned genotype 'value';
    * subPop/indRange: subPop or indRange should have the same length
    as values. Each item of values will be assigned to each subPop or
    indRange.

"; 

%feature("docstring") simuPOP::initByValue::initByValue "

Description:

    initialize populations by given alleles

Usage:

    initByValue(value=[], loci=[], atPloidy=-1, subPop=[],
      indRange=[], proportions=[], maleFreq=0.5, sex=[],
      stage=PreMating, begin=0, end=1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    value:          an array of genotypes of one  individual, having
                    the same length as the length of loci() or
                    loci()*ploidy() or pop.genoSize() (whole genotype)
                    or totNumLoci() (one copy of chromosome). This
                    parameter can also be an array of arrays of
                    genotypes of one  individual. Should have length
                    one or equal to subpop or ranges or proportion. If
                    value is an array of values, it should have the
                    same length as subpop, indRange or proportions.
    loci:           a vector of loci indices. If empty, apply to all
                    loci.
    locus:          a shortcut to loci
    atPloidy:       initialize which copy of chromosomes. Default to
                    all.
    subPop:         an array of applicable subpopulations. If values
                    are given, should have equal length to values.
    indRange:       a [begin, end] pair of range of individuals; or an
                    array of [begin, end] pairs.
    proportions:    an array of percentages for each item in values.
                    If given, assign given genotypes randomly.
    maleFreq:       male frequency
    sex:            an array of sex [Male, Female, Male...] for
                    individuals. The length of sex will not be
                    checked. If length of sex is shorter than the
                    number of individuals, sex will be reused from the
                    beginning.
    stages:         default to PreMating

Example:

>>> simu = simulator( population(subPop=[2,3], loci=[5,7]),
...     randomMating(), rep=1)
>>> simu.step([
...     initByValue([1]*5 + [2]*7 + [3]*5 +[4]*7),
...     dumper(alleleOnly=True)])
individual info: 
sub population 0:
   0: MU   3  3  3  3  3   2  2  2  2  2  2  2 |   1  1  1  1  1   4  4  4  4  4  4  4 
   1: FU   3  3  3  3  3   2  2  2  2  2  2  2 |   1  1  1  1  1   2  2  2  2  2  2  2 
sub population 1:
   2: FU   3  3  3  3  3   2  2  2  2  2  2  2 |   1  1  1  1  1   4  4  4  4  4  4  4 
   3: MU   3  3  3  3  3   2  2  2  2  2  2  2 |   1  1  1  1  1   2  2  2  2  2  2  2 
   4: MU   1  1  1  1  1   4  4  4  4  4  4  4 |   3  3  3  3  3   2  2  2  2  2  2  2 
End of individual info.


No ancenstral population recorded.
True
>>>


"; 

%feature("docstring") simuPOP::initByValue::~initByValue "

Description:

    simuPOP::initByValue::~initByValue

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

    apply operator  initByValue???

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

Details:

    Initializers are used to initialize populations before evolution.
    They are set to be PreMating operators by default.  simuPOP
    provides three initializers. One assigns alleles by random, one
    assigns a fixed set of genotypes, and the last one calls a user-
    defined function.

"; 

%feature("docstring") simuPOP::initializer::initializer "

Description:

    create an  initializer. default to be always active

Usage:

    initializer(subPop=[], indRange=[], loci=[], atPloidy=-1,
      maleFreq=0.5, sex=[], stage=PreMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::initializer::~initializer "

Description:

    destructor

Usage:

    x.~initializer()

"; 

%feature("docstring") simuPOP::initializer::clone "

Description:

    deep copy of an  initializer

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initializer::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  initializer

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::initializer::setRanges "

Description:

    set the range of a  populationpop???

Usage:

    x.setRanges(pop)

"; 

%feature("docstring") simuPOP::initializer::initSexIter "

Description:

    ???

Usage:

    x.initSexIter()

"; 

%feature("docstring") simuPOP::initializer::nextSex "

Description:

    ???

Usage:

    x.nextSex()

"; 

%feature("docstring") simuPOP::IOError "

Description:

    exception, thrown if file io failure

"; 

%feature("docstring") simuPOP::IOError::IOError "

Description:

    simuPOP::IOError::IOError

Usage:

    IOError(msg)

"; 

%ignore simuPOP::isAffected;

%feature("docstring") simuPOP::isAffected::isAffected "

Description:

    simuPOP::isAffected::isAffected

Usage:

    isAffected()

"; 

%feature("docstring") simuPOP::kamMutator "

Function form:

    KamMutate

Function form:

    KamMutate

Description:

    K-Allele Model  mutator.

Details:

    This  mutator mutate an allele to another allelic state with equal
    probability. The specified mutation rate is actually the
    'probability to mutate'. So the mutation rate to any other allelic
    state is actually  $ (rate/(K-1)) $, where  $ K $ is specified by
    parameter maxAllele. You can also specify states for this
    mutator. If the state parameter is given, all alleles must be one
    of the states, and mutation will happen among them. states is
    defaulted to 1-maxAllele.???

Details:

    This  mutator mutate an allele to another allelic state with equal
    probability. The specified mutation rate is actually the
    'probability to mutate'. So the mutation rate to any other allelic
    state is actually  $ (rate/(K-1)) $, where  $ K $ is specified by
    parameter maxAllele. You can also specify states for this
    mutator. If the state parameter is given, all alleles must be one
    of the states, and mutation will happen among them. states is
    defaulted to 1-maxAllele.???

"; 

%feature("docstring") simuPOP::kamMutator::kamMutator "

Description:

    create a K-Allele Model  mutator

Usage:

    kamMutator(rate=[], loci=[], maxAllele=0, output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    rate:           mutation rate. It is the 'probability to mutate'.
                    The actual mutation rate to any of the other K-1
                    allelic states are rates/(K-1).
    loci:           a vector of loci indices. Will be ignored only
                    when single rate is specified. Default to all
                    loci.
    maxAllele:      maximum allele that can be mutated to. For binary
                    libraries allelic states will be [0, maxAllele].
                    Otherwise, they are [1, maxAllele].

Example:

>>> simu = simulator(population(size=5, loci=[3,5]), noMating())
>>> simu.step([
...     kamMutator( rate=[.2,.6,.5], loci=[0,2,6], maxAllele=9),
...     dumper(alleleOnly=True)])
individual info: 
sub population 0:
   0: MU   0  0  0   0  0  0  0  0 |   0  0  6   0  0  0  1  0 
   1: MU   0  0  1   0  0  0  0  0 |   0  0  8   0  0  0  0  0 
   2: MU   0  0  8   0  0  0  0  0 |   0  0  1   0  0  0  0  0 
   3: MU   6  0  2   0  0  0  7  0 |   0  0  0   0  0  0  5  0 
   4: MU   0  0  0   0  0  0  5  0 |   0  0  8   0  0  0  1  0 
End of individual info.


No ancenstral population recorded.
True
>>>


"; 

%feature("docstring") simuPOP::kamMutator::~kamMutator "

Description:

    simuPOP::kamMutator::~kamMutator

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

%feature("docstring") simuPOP::largePedigreeSample "

Description:

    draw a large pedigree  sample

"; 

%feature("docstring") simuPOP::largePedigreeSample::largePedigreeSample "

Description:

    draw a large pedigree  sample

Usage:

    largePedigreeSample(size=[], minTotalSize=0, maxOffspring=5,
      minPedSize=5, minAffected=0, countOnly=False, name=\"sample\",
      nameExpr=\"\", times=1, saveAs=\"\", saveAsExpr=\"\", format=\"auto\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[\"father_idx\", \"mother_idx\"])

Details:

    Please refer to class  sample for other parameter descriptions.

Details:

    Please refer to class  sample for other parameter descriptions.

Arguments:

    minTotalSize:   the minimum number of individuals in the  sample
    maxOffspring:   the maximum number of offspring a parent may have
    minPedSize:     minimal pedigree size, default to 5
    minAffected:    minimal number of affected individuals in each
                    pedigree, default to 0
    countOnly:      set variables about number of affected sibpairs,
                    do not actually draw the  sample.

"; 

%feature("docstring") simuPOP::largePedigreeSample::~largePedigreeSample "

Description:

    destructor

Usage:

    x.~largePedigreeSample()

"; 

%feature("docstring") simuPOP::largePedigreeSample::clone "

Description:

    deep copy of a  largePedigreeSample operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::largePedigreeSample::prepareSample "

Description:

    preparation before drawing a  sample

Usage:

    x.prepareSample(pop)

"; 

%feature("docstring") simuPOP::largePedigreeSample::drawsample "

Description:

    draw a a large pedigree  sample

Usage:

    x.drawsample(pop)

"; 

%feature("docstring") simuPOP::largePedigreeSample::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  largePedigreeSample operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::maPenetrance "

Function form:

    MaPenetrance

Function form:

    MaPenetrance

Description:

    multiple allele  penetrance operator

Details:

    This is called 'multiple-alleles'???  penetrance. It separates
    alleles into two groups: wildtype and disease alleles. Wildtype
    alleles are specified by parameter wildtype and any other alleles
    are considered as diseased alleles.  maPenetrance accepts an array
    of fitness for AA, Aa, aa in the single-locus case, and a longer
    table for multi-locus case. Penetrance is then set for any given
    genotype.

Details:

    This is called 'multiple-alleles'???  penetrance. It separates
    alleles into two groups: wildtype and disease alleles. Wildtype
    alleles are specified by parameter wildtype and any other alleles
    are considered as diseased alleles.  maPenetrance accepts an array
    of fitness for AA, Aa, aa in the single-locus case, and a longer
    table for multi-locus case. Penetrance is then set for any given
    genotype.

"; 

%feature("docstring") simuPOP::maPenetrance::maPenetrance "

Description:

    create a multiple allele  penetrance operator ( penetrance
    according to diseased or wildtype alleles)

Usage:

    maPenetrance(loci, penet, wildtype, ancestralGen=-1,
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    locus:          the locus index. The genotype of this locus will
                    be examed.???
    loci:           the loci indices. The genotypes of these loci will
                    be examed.
    penetrance:     an array of  penetrance values of AA, Aa, aa. A is
                    the wild type group. In the case of multiple loci,
                    fitness should be in the order of AABB, AABb,
                    AAbb, AaBB, AaBb, Aabb, aaBB, aaBb, aabb.
    wildtype:       an array of alleles in the wildtype group. Any
                    other alleles will be considered as in the disease
                    allele group.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)???

"; 

%feature("docstring") simuPOP::maPenetrance::~maPenetrance "

Description:

    simuPOP::maPenetrance::~maPenetrance

Usage:

    x.~maPenetrance()

"; 

%feature("docstring") simuPOP::maPenetrance::clone "

Description:

    deep copy of a multi-allele  penetrance operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::maPenetrance::penet "

Description:

    currently assuming diploid???

Usage:

    x.penet(ind)

"; 

%feature("docstring") simuPOP::maPenetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the multi-allele  penetrance operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mapPenetrance "

Function form:

    MapPenetrance

Function form:

    MapPenetrance

Description:

    penetrance according to the genotype at one locus

Details:

    Assign  penetrance using a table with keys 'X-Y' where X and Y are
    allele numbers.

Details:

    Assign  penetrance using a table with keys 'X-Y' where X and Y are
    allele numbers.

"; 

%feature("docstring") simuPOP::mapPenetrance::mapPenetrance "

Description:

    create a map  penetrance operator

Usage:

    mapPenetrance(loci, penet, phase=False, ancestralGen=-1,
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    locus:          the locus index. The genotype of this locus will
                    be examed.???
    loci:           the loci indices. The genotypes of these loci will
                    be examed.
    penetrance:     a dictionary of  penetrance. The genotype must be
                    in the form of 'a-b' for a single locus.
    phase:          if True, a/b and b/a will have different
                    penetrance values. Default to False.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)???

"; 

%feature("docstring") simuPOP::mapPenetrance::~mapPenetrance "

Description:

    simuPOP::mapPenetrance::~mapPenetrance

Usage:

    x.~mapPenetrance()

"; 

%feature("docstring") simuPOP::mapPenetrance::clone "

Description:

    deep copy of a map  penetrance operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mapPenetrance::penet "

Description:

    currently assuming diploid???

Usage:

    x.penet(ind)

"; 

%feature("docstring") simuPOP::mapPenetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the map  penetrance operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mapQuanTrait "

Function form:

    MapQuanTrait

Function form:

    MapQuanTrait

Description:

    quantitative trait according to genotype at one locus

Details:

    Assign quantitative trait using a table with keys 'X-Y' where X
    and Y are allele numbers. If parameter sigma is not zero, the
    returned value is the sum of the trait plus  $
    N\\left(0,\\sigma^{2}\\right) $. This random part is usually
    considered as the environmental factor of the trait.

Details:

    Assign quantitative trait using a table with keys 'X-Y' where X
    and Y are allele numbers. If parameter sigma is not zero, the
    returned value is the sum of the trait plus  $
    N\\left(0,\\sigma^{2}\\right) $. This random part is usually
    considered as the environmental factor of the trait.

"; 

%feature("docstring") simuPOP::mapQuanTrait::mapQuanTrait "

Description:

    create a map quantitative trait operator

Usage:

    mapQuanTrait(loci, qtrait, sigma=0, phase=False,
      ancestralGen=-1, stage=PostMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[\"qtrait\"])

Arguments:

    locus:          the locus index. The quantitative trait is
                    determined by genotype at this locus.
    loci:           an array of locus indices. The quantitative trait
                    is determined by genotype at these loci.
    qtrait:         a dictionary of quantitative traits. The genotype
                    must be in the form of 'a-b'. This is the mean of
                    the quantitative trait. The actual trait value
                    will be N(mean, sigma^2). For multiple loci, the
                    form is 'a-b|c-d|e-f' etc.
    sigma:          standard deviation of the environmental factor
                    N(0, sigma^2).
    phase:          if True, a/b and b/a will have different
                    quantitative trait values. Default to False.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::mapQuanTrait::~mapQuanTrait "

Description:

    simuPOP::mapQuanTrait::~mapQuanTrait

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

Function form:

    MapSelector

Description:

    selection according to the genotype at one locus

Details:

    This map  selector implements selection at one locus. A user
    provided dictionary (map) of genotypes will be used in this
    selector to set each individual's fitness value.

Details:

    This map  selector implements selection at one locus. A user
    provided dictionary (map) of genotypes will be used in this
    selector to set each individual's fitness value.

"; 

%feature("docstring") simuPOP::mapSelector::mapSelector "

Description:

    create a map  selector

Usage:

    mapSelector(loci, fitness, phase=False, subPops=[],
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[\"fitness\"])

Arguments:

    locus:          the locus index. The genotype of this locus will
                    be examed.???
    loci:           the locus indices. The genotypes of these loci
                    will be examed.
    fitness:        a dictionary of fitness values. The genotype must
                    be in the form of 'a-b' for a single locus, and
                    'a-b|c-d|e-f' for multi-loci.
    phase:          if True, genotypes a-b and b-a will have different
                    fitness values. Default to false.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)???

"; 

%feature("docstring") simuPOP::mapSelector::~mapSelector "

Description:

    simuPOP::mapSelector::~mapSelector

Usage:

    x.~mapSelector()

"; 

%feature("docstring") simuPOP::mapSelector::clone "

Description:

    deep copy of a map  selector

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
    of the map  selector

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::maQuanTrait "

Function form:

    MaQuanTrait

Function form:

    MaQuanTrait

Description:

    multiple allele quantitative trait (quantitative trait according
    to disease or wildtype alleles)

Details:

    This is called 'multiple-allele' quantitative trait. It separates
    alleles into two groups: wildtype and disease susceptibility
    alleles. Wildtype alleles are specified by parameter wildtype and
    any other alleles are considered as disease alleles.  maQuanTrait
    accepts an array of fitness. Quantitative trait is then set for
    any given genotype. A standard normal distribution  $
    N\\left(0,\\sigma^{2}\\right) $ will be added to the returned trait
    value.

Details:

    This is called 'multiple-allele' quantitative trait. It separates
    alleles into two groups: wildtype and disease susceptibility
    alleles. Wildtype alleles are specified by parameter wildtype and
    any other alleles are considered as disease alleles.  maQuanTrait
    accepts an array of fitness. Quantitative trait is then set for
    any given genotype. A standard normal distribution  $
    N\\left(0,\\sigma^{2}\\right) $ will be added to the returned trait
    value.

"; 

%feature("docstring") simuPOP::maQuanTrait::maQuanTrait "

Description:

    create a multiple allele quantitative trait operator

Usage:

    maQuanTrait(loci, qtrait, wildtype, sigma=[], ancestralGen=-1,
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[\"qtrait\"])

Details:

    Please refer to  quanTrait for other parameter descriptions.

Details:

    Please refer to  quanTrait for other parameter descriptions.

Arguments:

    qtrait:         an array of quantitative traits of AA, Aa, aa. A
                    is the wild type group
    sigma:          an array of standard deviations for each of the
                    trait genotype (AA, Aa, aa)
    wildtype:       an array of alleles in the wildtype group. Any
                    other alleles will be considered as disease
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

Function form:

    MaSelect

Description:

    multiple allele  selector (selection according to wildtype or
    diseased alleles)

Details:

    This is called 'multiple-allele'  selector. It separate alleles
    into two groups: wildtype and disease alleles. Wildtype alleles
    are specified by parameter wildtype and any other alleles are
    considered as diseased alleles.

Details:

    This is called 'multiple-allele'  selector. It separate alleles
    into two groups: wildtype and disease alleles. Wildtype alleles
    are specified by parameter wildtype and any other alleles are
    considered as diseased alleles.

"; 

%feature("docstring") simuPOP::maSelector::maSelector "

Description:

    create a multiple allele  selector

Usage:

    maSelector(loci, fitness, wildtype, subPops=[], stage=PreMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[\"fitness\"])

Details:

    Please refer to  mapSelector for other parameter descriptions.

Details:

    Please refer to  mapSelector for other parameter descriptions.

Arguments:

    fitness:        for the single locus case, fitness is an array of
                    fitness of AA, Aa, aa. A is the wildtype group. In
                    the case of multiple loci, fitness should be in
                    the order of AABB, AABb, AAbb, AaBB, AaBb, Aabb,
                    aaBB, aaBb, aabb.
    wildtype:       an array of alleles in the wildtype group. Any
                    other alleles are considered to be diseased
                    alleles. Default to [0].
    output:         and other parameters please refer to
                    help(baseOperator.__init__)???

Note:

    *  maSelector only works for diploid populations now.
    * wildtype at all loci are the same.

"; 

%feature("docstring") simuPOP::maSelector::~maSelector "

Description:

    simuPOP::maSelector::~maSelector

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

    the base class of all  mating schemes - a required parameter of
    simulator

Details:

    Mating schemes specify how to generate offspring from the current
    population. It must be provided when a  simulator is created.
    Mating can perform the following tasks:
    * change population/subpopulation sizes;
    * randomly select parent(s) to generate offspring to populate the
    offspring generation;
    * apply during-mating operators;
    * apply selection if applicable.

Details:

    Mating schemes specify how to generate offspring from the current
    population. It must be provided when a  simulator is created.
    Mating can perform the following tasks:
    * change population/subpopulation sizes;
    * randomly select parent(s) to generate offspring to populate the
    offspring generation;
    * apply during-mating operators;
    * apply selection if applicable.

"; 

%ignore simuPOP::mating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::mating::mating "

Description:

    create a  mating scheme (do not use this base  mating scheme, use
    one of its derived classes instead)

Usage:

    mating(numOffspring=1.0, numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
      newSubPopSizeExpr=\"\", newSubPopSizeFunc=None)

Details:

    By default, a  mating scheme keeps a constant  population size,
    generates one offspring per  mating event. These can be changed
    using certain parameters. newSubPopSize, newSubPopSizeExpr and
    newSubPopSizeFunc can be used to specify subpopulation sizes of
    the offspring generation. mode, numOffspring, maxNumOffspring can
    be used to specify how many offspring will be produced at each
    mating event. This mode parameter can be one of
    * MATE_NumOffspring: a fixed number of offspring at all  mating
    events at this generation. If numOffspring is given, all
    generations use this fixed number. If numOffspringFunc is given,
    the number of offspring at each generation is determined by the
    value returned from this function.
    * MATE_NumOffspringEachFamily: each family can have its own number
    of offspring. Usually, numOffspringFunc is used to determine the
    number of offspring of each family. If numOffspring is used,
    MATE_NumOffspringEachFamily is equivalent to MATE_NumOffspring.
    * MATE_GeometricDistribution: a Geometric distribution with
    parameter numOffspring is used to determine the number of
    offspring of each family.
    * MATE_PoissonDistribution: a Poisson distribution with parameter
    numOffspring is used to determine the number of offspring of each
    family.
    * MATE_BinomialDistribution: a Binomial distribution with
    parameter numOffspring is used to determine the number of
    offspring of each family.
    * MATE_UniformDistribution: a Uniform distribution  [a, b]  with
    parameter numOffspring (a) and maxNumOffspring (b) is used to
    determine the number of offspring of each family.

Details:

    By default, a  mating scheme keeps a constant  population size,
    generates one offspring per  mating event. These can be changed
    using certain parameters. newSubPopSize, newSubPopSizeExpr and
    newSubPopSizeFunc can be used to specify subpopulation sizes of
    the offspring generation. mode, numOffspring, maxNumOffspring can
    be used to specify how many offspring will be produced at each
    mating event. This mode parameter can be one of
    * MATE_NumOffspring: a fixed number of offspring at all  mating
    events at this generation. If numOffspring is given, all
    generations use this fixed number. If numOffspringFunc is given,
    the number of offspring at each generation is determined by the
    value returned from this function.
    * MATE_NumOffspringEachFamily: each family can have its own number
    of offspring. Usually, numOffspringFunc is used to determine the
    number of offspring of each family. If numOffspring is used,
    MATE_NumOffspringEachFamily is equivalent to MATE_NumOffspring.
    * MATE_GeometricDistribution: a Geometric distribution with
    parameter numOffspring is used to determine the number of
    offspring of each family.
    * MATE_PoissonDistribution: a Poisson distribution with parameter
    numOffspring is used to determine the number of offspring of each
    family.
    * MATE_BinomialDistribution: a Binomial distribution with
    parameter numOffspring is used to determine the number of
    offspring of each family.
    * MATE_UniformDistribution: a Uniform distribution  [a, b]  with
    parameter numOffspring (a) and maxNumOffspring (b) is used to
    determine the number of offspring of each family.

Arguments:

    numOffspring:   the number of offspring or p for a random
                    distribution. Default to 1. This parameter
                    determines the number of offspring that a  mating
                    event will produce. Therefore, it determines the
                    family size.
    numOffspringFunc:a Python function that returns the number of
                    offspring or p
    maxNumOffspring:used when numOffspring is generated from a
                    binomial distribution
    mode:           can be one of MATE_NumOffspring,
                    MATE_NumOffspringEachFamily,
                    MATE_GeometricDistribution,
                    MATE_PoissonDistribution,
                    MATE_BinomialDistribution,
                    MATE_UniformDistribution.
    newSubPopSize:  an array of subpopulaitons sizes
    newSubPopSizeExpr:an expression that will return the new
                    subpopulation size
    newSubPopSizeFunc:a function that accepts an int
                    parameter(generation), an array of current
                    population size and return an array of
                    subpopulation sizes. This is usually easier to use
                    than its expression version of this parameter.

Example:

>>> # arbitrary demographic model
>>> def lin_inc(gen, oldsize=[]):
...     return [10+gen]*5
... 
>>> simu = simulator(
...     population(subPop=[5]*5, loci=[1]),
...     randomMating(newSubPopSizeFunc=lin_inc)
... )
>>> simu.evolve(
...     ops = [
...         stat(popSize=True),
...         pyEval(r'\"%d %d\\\\n\"%(gen, subPop[0][\"popSize\"])'),
...     ],
...     end=5
... )
0 10
1 11
2 12
3 13
4 14
5 15
True
>>> 
>>> #
>>> # control the number of offspring per mating event
>>> # famSizes is only defined when DBG_MATING is defined
>>> TurnOnDebug(DBG_MATING)
>>> simu = simulator(population(50, loci=[1]),
...     randomMating(numOffspring=2, 
...         maxNumOffspring=5,
...         mode=MATE_UniformDistribution))
>>> simu.step(ops=[])
True
>>> print simu.population(0).dvars().famSizes
[5, 3, 2, 5, 2, 2, 5, 3, 3, 3, 4, 3, 2, 5, 3]
>>> TurnOffDebug(DBG_MATING)
Debug code DBG_MATING is turned off. cf. ListDebugCode(), TurnOnDebug().
>>>


"; 

%ignore simuPOP::mating::mating(const mating &rhs);

%feature("docstring") simuPOP::mating::~mating "

Description:

    destructor

Usage:

    x.~mating()

"; 

%feature("docstring") simuPOP::mating::clone "

Description:

    deep copy of a  mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::mating::submitScratch(population &pop, population &scratch);

%ignore simuPOP::mating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%ignore simuPOP::mating::fixedFamilySize();

%ignore simuPOP::mating::numOffspring(int gen);

%ignore simuPOP::mating::resetNumOffspring();

%ignore simuPOP::mating::prepareScratchPop(population &pop, population &scratch);

%feature("docstring") simuPOP::mergeSubPops "

Function form:

    MergeSubPops

Function form:

    MergeSubPops

Description:

    merge subpopulations

Details:

    This operator merges subpopulations subPops (the only
    parameter???) to a single subpopulation. If subPops is ignored,
    all subpopulations will be merged.

Details:

    This operator merges subpopulations subPops (the only
    parameter???) to a single subpopulation. If subPops is ignored,
    all subpopulations will be merged.

"; 

%feature("docstring") simuPOP::mergeSubPops::mergeSubPops "

Description:

    merge subpopulations

Usage:

    mergeSubPops(subPops=[], removeEmptySubPops=False,
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    subPops:        subpopulatinos to be merged. Default to all.

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

    migrate individuals from a (sub) population to another (sub)
    population

Details:

    Migrator is the only way to mix genotypes of several
    subpopulations because  mating is strictly within subpopulations
    in  simuPOP. Migrators are quite flexible in  simuPOP in the sense
    that
    * Migration can happen from and to a subset of subpopulations.
    * Migration can be done by probability, proportion or by counts.
    In the case of probability, if the migration rate from
    subpopulation a to b is r, then everyone in subpopulation a will
    have this probability to migrate to b. In the case of proportion,
    exactly r*size_of_subPop_a individuals (chosen by random) will
    migrate to subpopulation b. In the last case, a given number of
    individuals will migrate.
    * New subpopulation can be generated through migration. You simply
    need to migrate to a new subpopulation number.

Details:

    Migrator is the only way to mix genotypes of several
    subpopulations because  mating is strictly within subpopulations
    in  simuPOP. Migrators are quite flexible in  simuPOP in the sense
    that
    * Migration can happen from and to a subset of subpopulations.
    * Migration can be done by probability, proportion or by counts.
    In the case of probability, if the migration rate from
    subpopulation a to b is r, then everyone in subpopulation a will
    have this probability to migrate to b. In the case of proportion,
    exactly r*size_of_subPop_a individuals (chosen by random) will
    migrate to subpopulation b. In the last case, a given number of
    individuals will migrate.
    * New subpopulation can be generated through migration. You simply
    need to migrate to a new subpopulation number.

"; 

%feature("docstring") simuPOP::migrator::migrator "

Description:

    create a  migrator

Usage:

    migrator(rate, mode=MigrByProbability, fromSubPop=[],
      toSubPop=[], stage=PreMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    rate:           migration rate, can be a proportion or counted
                    number. Determined by parameter mode. rate should
                    be an m by n matrix. If a number is given, the
                    migration rate will be r*ones(m,n)???.
    mode:           one of MigrByProbability (default),
                    MigrByProportion or MigrByCounts
    fromSubPop:     an array of 'from' subpopulations. Default to all.
                    If a single subpop is specified, [] can be
                    ignored. I.e., [a] is equvalent to a.
    toSubPop:       an array of 'to' subpopulations. Default to all
                    subpopulations. If a single subpop is specified,
                    [] can be ignored.
    stage:          default to PreMating

Note:

    * The overall  population size will not be changed. (Mating
    schemes can do that). If you would like to keep the subpopulation
    size after migration, you can use the newSubPopSize or
    newSubPopSizeExpr parameter of a  mating scheme.
    * rate is a matrix with dimensions determined by fromSubPop and
    toSubPop. By default, rate is a matrix with element r(i,j), where
    r(i, j) is the migration rate, probability or count from
    subpopulation i to j. If fromSubPop and/or toSubPop are given,
    migration will only happen between these subpopulations. An
    extreme case is 'point migration', rate=[[r]], fromSubPop=a,
    toSubPop=b which migrate from subpopulation a to b with given rate
    r.???

"; 

%feature("docstring") simuPOP::migrator::~migrator "

Description:

    destructor

Usage:

    x.~migrator()

"; 

%feature("docstring") simuPOP::migrator::clone "

Description:

    deep copy of a  migrator

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

Details:

    Format should be 0-0 0-1 0-2, 1-0 1-1 1-2, 2-0, 2-1, 2-2. For mode
    MigrByProbability or MigrByProportion, 0-0,1-1,2-2 will be set
    automatically regardless of input.

"; 

%feature("docstring") simuPOP::migrator::apply "

Description:

    apply the  migrator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::migrator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  migrator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mlPenetrance "

Function form:

    MlPenetrance

Function form:

    MlPenetrance

Description:

    penetrance according to the genotype according to a multiple loci
    multiplicative model

Details:

    mlPentrance is the 'multiple-loci'??? penetrnace calculator. It
    accepts a list of penetrances and combine them according to the
    mode parameter, which takes one of the following values:
    * PEN_Multiplicative: the  penetrance is calculated as  $ f=\\prod
    f_{i} $.
    * PEN_Additive: the  penetrance is calculated as  $
    f=\\min\\left(1,\\sum f_{i}\\right) $.  $ f $ will be set to 1 when  $
    f<0 $. In this case,  $ s_{i} $??? are added, not  $ f_{i} $
    directly.
    * PEN_Heterogeneity: the  penetrance is calculated as  $
    f=1-\\prod\\left(1-f_{i}\\right) $. Please refer to Neil Risch (1990)
    for detailed information about these models.

Details:

    mlPentrance is the 'multiple-loci'??? penetrnace calculator. It
    accepts a list of penetrances and combine them according to the
    mode parameter, which takes one of the following values:
    * PEN_Multiplicative: the  penetrance is calculated as  $ f=\\prod
    f_{i} $.
    * PEN_Additive: the  penetrance is calculated as  $
    f=\\min\\left(1,\\sum f_{i}\\right) $.  $ f $ will be set to 1 when  $
    f<0 $. In this case,  $ s_{i} $??? are added, not  $ f_{i} $
    directly.
    * PEN_Heterogeneity: the  penetrance is calculated as  $
    f=1-\\prod\\left(1-f_{i}\\right) $. Please refer to Neil Risch (1990)
    for detailed information about these models.

"; 

%feature("docstring") simuPOP::mlPenetrance::mlPenetrance "

Description:

    create a multiple loci  penetrance operator using a multiplicative
    model

Usage:

    mlPenetrance(peneOps, mode=PEN_Multiplicative, ancestralGen=-1,
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    peneOps:        a list of selectors???
    mode:           can be one of PEN_Multiplicative, PEN_Additive,
                    and PEN_Heterogeneity

"; 

%feature("docstring") simuPOP::mlPenetrance::~mlPenetrance "

Description:

    simuPOP::mlPenetrance::~mlPenetrance

Usage:

    x.~mlPenetrance()

"; 

%feature("docstring") simuPOP::mlPenetrance::clone "

Description:

    deep copy of a multi-loci  penetrance operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mlPenetrance::penet "

Description:

    currently assuming diploid???

Usage:

    x.penet(ind)

"; 

%feature("docstring") simuPOP::mlPenetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the multiple-loci  penetrance operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mlQuanTrait "

Function form:

    MlQuanTrait

Function form:

    MlQuanTrait

Description:

    quantitative trait according to genotypes from a multiple loci
    multiplicative model

Details:

    mlQuanTrait is a 'multiple-loci' quantitative trait calculator. It
    accepts a list of quantitative traits and combine them according
    to the mode parameter, which takes one of the following values
    * QT_Multiplicative: the mean of the quantitative trait is
    calculated as  $ f=\\prod f_{i} $.
    * QT_Additive: the mean of the quantitative trait is calculated as
    $ f=\\sum f_{i} $. Note that all  $ \\sigma_{i} $ (for  $ f_{i} $)
    and  $ \\sigma $ (for  $ f $) will all be considered. I.e, the
    trait value should be
    $ f=\\sum_{i}\\left(f_{i}+N\\left(0,\\sigma_{i}^{2}\\right)\\right)+\\sig
    ma^{2} $ for QT_Additive case. If this is not desired, you can set
    some of the  $ \\sigma $ to zero.

Details:

    mlQuanTrait is a 'multiple-loci' quantitative trait calculator. It
    accepts a list of quantitative traits and combine them according
    to the mode parameter, which takes one of the following values
    * QT_Multiplicative: the mean of the quantitative trait is
    calculated as  $ f=\\prod f_{i} $.
    * QT_Additive: the mean of the quantitative trait is calculated as
    $ f=\\sum f_{i} $. Note that all  $ \\sigma_{i} $ (for  $ f_{i} $)
    and  $ \\sigma $ (for  $ f $) will all be considered. I.e, the
    trait value should be
    $ f=\\sum_{i}\\left(f_{i}+N\\left(0,\\sigma_{i}^{2}\\right)\\right)+\\sig
    ma^{2} $ for QT_Additive case. If this is not desired, you can set
    some of the  $ \\sigma $ to zero.

"; 

%feature("docstring") simuPOP::mlQuanTrait::mlQuanTrait "

Description:

    multiple loci quantitative trait using a multiplicative model

Usage:

    mlQuanTrait(qtraits, mode=QT_Multiplicative, sigma=0,
      ancestralGen=-1, stage=PostMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[\"qtrait\"])

Details:

    Please refer to  quanTrait for other parameter descriptions.

Details:

    Please refer to  quanTrait for other parameter descriptions.

Arguments:

    qtraits:        a list of quantitative traits
    mode:           can be one of QT_Multiplicative and QT_Additive

"; 

%feature("docstring") simuPOP::mlQuanTrait::~mlQuanTrait "

Description:

    simuPOP::mlQuanTrait::~mlQuanTrait

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

Function form:

    MlSelect

Description:

    selection according to genotypes at multiple loci in a
    multiplicative model

Details:

    This  selector is a 'multiple-loci model'  selector. The  selector
    takes a vector of selectors (can not be another  mlSelector) and
    evaluate the fitness of an  individual as the the product or sum
    of  individual fitness values. The mode is determined by parameter
    mode, which takes the value
    * SEL_Multiplicative: the fitness is calculated as  $
    f=\\prod_{i}f_{i} $.
    * SEL_Additive: the fitness is calculated as  $ f=\\max\\left(0,1-\\s
    um_{i}(1-f_{i})\\right)=\\max\\left(0,1-\\sum_{i}s_{i}\\right) $.  $ f
    $ will be set to 0 when  $ f<0 $. In this case,  $ s_{i} $ are
    added, not  $ f_{i} $ directly.

Details:

    This  selector is a 'multiple-loci model'  selector. The  selector
    takes a vector of selectors (can not be another  mlSelector) and
    evaluate the fitness of an  individual as the the product or sum
    of  individual fitness values. The mode is determined by parameter
    mode, which takes the value
    * SEL_Multiplicative: the fitness is calculated as  $
    f=\\prod_{i}f_{i} $.
    * SEL_Additive: the fitness is calculated as  $ f=\\max\\left(0,1-\\s
    um_{i}(1-f_{i})\\right)=\\max\\left(0,1-\\sum_{i}s_{i}\\right) $.  $ f
    $ will be set to 0 when  $ f<0 $. In this case,  $ s_{i} $ are
    added, not  $ f_{i} $ directly.

"; 

%feature("docstring") simuPOP::mlSelector::mlSelector "

Description:

    create a multi-loci  selector

Usage:

    mlSelector(selectors, mode=SEL_Multiplicative, subPops=[],
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[\"fitness\"])

Details:

    Please refer to  mapSelector for other parameter descriptions.

Details:

    Please refer to  mapSelector for other parameter descriptions.

Arguments:

    selectors:      a list of selectors

"; 

%feature("docstring") simuPOP::mlSelector::~mlSelector "

Description:

    simuPOP::mlSelector::~mlSelector

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

%feature("docstring") simuPOP::mutator "

Description:

    mutator class

Details:

    The base class of all functional mutators. It is not supposed to
    be called directly.Every  mutator can specify rate (equal rate or
    different rates for different loci) and a vector of applicable
    loci (default to all but should have the same length as rate if
    rate has length greater than one).Maximum allele can be specified
    as well but more parameter, if needed, should be implemented by
    individual mutator classes.There are number of possible allelic
    states. Most theoretical studies assume an infinite number of
    allelic states to avoid any homoplasy. If it facilitates any
    analysis, this is however extremely unrealistic.

Details:

    The base class of all functional mutators. It is not supposed to
    be called directly.Every  mutator can specify rate (equal rate or
    different rates for different loci) and a vector of applicable
    loci (default to all but should have the same length as rate if
    rate has length greater than one).Maximum allele can be specified
    as well but more parameter, if needed, should be implemented by
    individual mutator classes.There are number of possible allelic
    states. Most theoretical studies assume an infinite number of
    allelic states to avoid any homoplasy. If it facilitates any
    analysis, this is however extremely unrealistic.

"; 

%feature("docstring") simuPOP::mutator::mutator "

Description:

    create a  mutator

Usage:

    mutator(rate=[], loci=[], maxAllele=0, output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Details:

    All mutators have the following common parameters. However, the
    actual meaning of these parameters may vary according to different
    model. The only differences between the following mutators are
    they way they actually mutate an allele, and corresponding input
    parameters. Mutators record the number of mutation events at each
    locus.

Details:

    All mutators have the following common parameters. However, the
    actual meaning of these parameters may vary according to different
    model. The only differences between the following mutators are
    they way they actually mutate an allele, and corresponding input
    parameters. Mutators record the number of mutation events at each
    locus.

Arguments:

    rate:           can be a number (uniform rate) or an array of
                    mutation rates (the same length as loci)
    loci:           a vector of loci indices. Will be ignored only
                    when single rate is specified. Default to all
                    loci.
    maxAllele:      maximum allowable allele. Interpreted by each sub
                    mutaor class. Default to  pop.maxAllele().

"; 

%feature("docstring") simuPOP::mutator::~mutator "

Description:

    destructor

Usage:

    x.~mutator()

"; 

%feature("docstring") simuPOP::mutator::clone "

Description:

    deep copy of a  mutator

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

    apply a  mutator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::noMating "

Description:

    a  mating scheme that does nothing

Details:

    In this scheme, there is
    * no  mating. Parent generation will be considered as offspring
    generation.
    * no subpopulation change. During-mating operators will be
    applied, but the return values are not checked. I.e., subpopsizes
    will be ignored although some during-mating operators may be
    applied.

Details:

    In this scheme, there is
    * no  mating. Parent generation will be considered as offspring
    generation.
    * no subpopulation change. During-mating operators will be
    applied, but the return values are not checked. I.e., subpopsizes
    will be ignored although some during-mating operators may be
    applied.

"; 

%feature("docstring") simuPOP::noMating::noMating "

Description:

    creat a scheme with no  mating

Usage:

    noMating(numOffspring=1.0, numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
      newSubPopSizeExpr=\"\", newSubPopSizeFunc=None)

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

    deep copy of a scheme with no  mating

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::noMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the scheme with no  mating

Usage:

    x.__repr__()

"; 

%ignore simuPOP::noMating::submitScratch(population &pop, population &scratch);

%ignore simuPOP::noMating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%feature("docstring") simuPOP::noneOp "

Description:

    none operator

"; 

%feature("docstring") simuPOP::noneOp::noneOp "

Usage:

    noneOp(output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
      end=0, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Details:

    This operator does nothing.

Details:

    This operator does nothing.

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

%ignore simuPOP::noneOp::applyDuringMating(population &pop, population::IndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::noneOp::apply "

Description:

    apply the  noneOp operator to one  population

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

%feature("docstring") simuPOP::nuclearFamilySample "

Description:

    draw a nuclear family  sample

"; 

%feature("docstring") simuPOP::nuclearFamilySample::nuclearFamilySample "

Description:

    draw a nuclear family  sample

Usage:

    nuclearFamilySample(size=[], minTotalSize=0, maxOffspring=5,
      minPedSize=5, minAffected=0, countOnly=False, name=\"sample\",
      nameExpr=\"\", times=1, saveAs=\"\", saveAsExpr=\"\", format=\"auto\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[\"father_idx\", \"mother_idx\"])

Details:

    Please refer to class  sample for parameter descriptions.

Details:

    Please refer to class  sample for parameter descriptions.

"; 

%feature("docstring") simuPOP::nuclearFamilySample::~nuclearFamilySample "

Description:

    destructor

Usage:

    x.~nuclearFamilySample()

"; 

%feature("docstring") simuPOP::nuclearFamilySample::clone "

Description:

    deep copy of a  nuclearFamilySample operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::nuclearFamilySample::prepareSample "

Description:

    preparation before drawing a  sample

Usage:

    x.prepareSample(pop)

"; 

%feature("docstring") simuPOP::nuclearFamilySample::drawsample "

Description:

    draw a nuclear family  sample

Usage:

    x.drawsample(pop)

"; 

%feature("docstring") simuPOP::nuclearFamilySample::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  nuclearFamilySample operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::offspringGenerator;

%feature("docstring") simuPOP::offspringGenerator::offspringGenerator "

Description:

    create an offspring generator, save information from pop and ops
    to speed up the calls to generateOffspring

Usage:

    offspringGenerator(pop, ops)

"; 

%feature("docstring") simuPOP::offspringGenerator::generateOffspring "

Description:

    generate numOff offspring, or until reach offEnd

Usage:

    x.generateOffspring(pop, dad, mom, numOff, offBegin)

Details:

    This is because offBegin+numOff may go beyond the subpopulation
    boundary.

Details:

    This is because offBegin+numOff may go beyond the subpopulation
    boundary.

"; 

%feature("docstring") simuPOP::offspringGenerator::copyOffspring "

Description:

    copy numOff offspring, or until reach offEnd

Usage:

    x.copyOffspring(pop, par, numOff, offBegin)

Details:

    This is because offBegin+numOff may go beyond the subpopulation
    boundary.

Details:

    This is because offBegin+numOff may go beyond the subpopulation
    boundary.

"; 

%ignore simuPOP::OstreamManager;

%ignore simuPOP::OstreamManager::OstreamManager();

%feature("docstring") simuPOP::OstreamManager::~OstreamManager "

Description:

    simuPOP::OstreamManager::~OstreamManager

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

Description:

    simuPOP::OutOfMemory::OutOfMemory

Usage:

    OutOfMemory(msg)

"; 

%feature("docstring") simuPOP::outputer "

Description:

    outputer is a (special) subclass of Operator that will output
    files with different format.

Details:

    Bo Peng

Details:

    Bo Peng

"; 

%feature("docstring") simuPOP::outputer::outputer "

Description:

    constructor. default to be always active.

Usage:

    outputer(output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
      end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

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

%feature("docstring") simuPOP::outputHelper "

Description:

    simuPOP::outputHelper

"; 

%feature("docstring") simuPOP::outputHelper::outputHelper "

Description:

    simuPOP::outputHelper::outputHelper

Usage:

    outputHelper(str=\"\\n\", output=\">\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::outputHelper::apply "

Description:

    simply output some info

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::outputHelper::~outputHelper "

Description:

    simuPOP::outputHelper::~outputHelper

Usage:

    x.~outputHelper()

"; 

%feature("docstring") simuPOP::outputHelper::clone "

Description:

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::outputHelper::setString "

Description:

    set output string.

Usage:

    x.setString(str)

"; 

%feature("docstring") simuPOP::outputHelper::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::parentsTagger "

Description:

    tagging according to parents' indices

Details:

    This during-mating operator set \\c tag(), currently a pair of
    numbers, of each  individual with indices of his/her parents in
    the parental  population. This information will be used by
    pedigree-related operators like  affectedSibpairSample to track
    the pedigree information. Since parental  population will be
    discarded or stored after  mating, and tagging information will be
    passed with individuals, mating/population change etc. will not
    interfere with this simple tagging system.

Details:

    This during-mating operator set \\c tag(), currently a pair of
    numbers, of each  individual with indices of his/her parents in
    the parental  population. This information will be used by
    pedigree-related operators like  affectedSibpairSample to track
    the pedigree information. Since parental  population will be
    discarded or stored after  mating, and tagging information will be
    passed with individuals, mating/population change etc. will not
    interfere with this simple tagging system.

"; 

%feature("docstring") simuPOP::parentsTagger::parentsTagger "

Description:

    create a  parentsTagger, default to be always active

Usage:

    parentsTagger(begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[\"father_idx\", \"mother_idx\"])

"; 

%feature("docstring") simuPOP::parentsTagger::~parentsTagger "

Description:

    simuPOP::parentsTagger::~parentsTagger

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

%feature("docstring") simuPOP::pause "

Description:

    pause a  simulator

Details:

    This operator pauses the evolution of a  simulator at given
    generations or at a key stroke, using stopOnKeyStroke=True option.
    Users can use 'q' to stop an evolution. When a  simulator is
    stopped, press any other key to resume the simulation or escape to
    a Python shell to examine the status of the simulation by press
    's'.
    There are two ways to use this operator, the first one is to
    pause the simulation at specified generations, using the usual
    operator parameters such as at. Another way is to  pause a
    simulation with any key stroke, using the stopOnKeyStroke
    parameter. This feature is useful for a presentation or an
    interactive simulation. When 's' is pressed, this operator expose
    the current  population to the main Python dictionary as variable
    pop and enter an interactive Python session. The way current
    population is exposed can be controlled by parameter exposePop and
    popName. This feature is useful when you want to examine the
    properties of a  population during evolution.

Details:

    This operator pauses the evolution of a  simulator at given
    generations or at a key stroke, using stopOnKeyStroke=True option.
    Users can use 'q' to stop an evolution. When a  simulator is
    stopped, press any other key to resume the simulation or escape to
    a Python shell to examine the status of the simulation by press
    's'.
    There are two ways to use this operator, the first one is to
    pause the simulation at specified generations, using the usual
    operator parameters such as at. Another way is to  pause a
    simulation with any key stroke, using the stopOnKeyStroke
    parameter. This feature is useful for a presentation or an
    interactive simulation. When 's' is pressed, this operator expose
    the current  population to the main Python dictionary as variable
    pop and enter an interactive Python session. The way current
    population is exposed can be controlled by parameter exposePop and
    popName. This feature is useful when you want to examine the
    properties of a  population during evolution.

"; 

%feature("docstring") simuPOP::pause::pause "

Description:

    stop a simulation. Press 'q' to exit or any other key to continue.

Usage:

    pause(prompt=True, stopOnKeyStroke=False, exposePop=True,
      popName=\"pop\", output=\">\", outputExpr=\"\", stage=PostMating,
      begin=0, end=-1, step=1, at=[], rep=REP_LAST, grp=GRP_ALL,
      infoFields=[])

Arguments:

    prompt:         if True (default), print prompt message
    stopOnKeyStroke:if True, stop only when a key was pressed
    exposePop:      whether or not expose pop to user namespace, only
                    useful when user choose 's' at  pause. Default to
                    True.
    popName:        by which name the  population is exposed. Default
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

    deep copy of a  pause operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pause::apply "

Description:

    apply the  pause operator to one  population

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pause::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  pause operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::penetrance "

Description:

    basic class of a  penetrance operator

Details:

    Penetrance is the probability that one will have the disease when
    he has certain genotype(s). Calculation and the parameter set of
    penetrance are similar to those of fitness. An  individual will be
    randomly marked as affected/unaffected according to his
    penetrance value.??? For example, an  individual will have
    probability 0.8 to be affected if the  penetrance is 0.8.
    Penetrance can be applied at any stage (default to DuringMating).
    It will be calculated during  mating, and then the affected status
    will be set for each offspring. Penetrance can also be used as
    PreMating, PostMating or even PrePostMating??? operator. In these
    cases, the affected status will be set to all individuals
    according to their  penetrance values. It is also possible to
    store  penetrance in a given information field specified by
    infoFields parameter (e.g. infoFields=['penetrance']). This is
    useful to check the  penetrance values at a later time.
    Affected status will be used for statistical purpose, and most
    importantly, ascertainment. They will be calculated along with
    fitness although they might not be used at every generation. You
    can use two operators: one for fitness/selection, active at every
    generation; one for affected status, active only at
    ascertainments, to avoid unnecessary calculation of the affected
    status.
    Pentrance values are used to set the affectedness of individuals,
    and are usually not saved. If you would like to know the
    penetrance value, you need to
    * use addInfoField('penetrance') to the  population to analyze.
    (Or use infoFields parameter of the  population constructor), and
    * use e.g.,  mlPenetrance(...., infoFields=['penetrance']) to add
    the  penetrance field to the  penetrance operator you use. You may
    choose a name other than 'penetrance' as long as the field names
    for the operator and  population match. Penetrance functions can
    be applied to the current, all, or certain number of ancestral
    generations. This is controlled by the ancestralGen parameter,
    which is default to -1 (all available ancestral generations). You
    can set it to 0 if you only need affection??? status for the
    current generation, or specify a number n for the number of
    ancestral generations (n + 1 total generations) to process. Note
    that ancestralGen parameter is ignored if the  penetrance operator
    is used as a during  mating operator.

Details:

    Penetrance is the probability that one will have the disease when
    he has certain genotype(s). Calculation and the parameter set of
    penetrance are similar to those of fitness. An  individual will be
    randomly marked as affected/unaffected according to his
    penetrance value.??? For example, an  individual will have
    probability 0.8 to be affected if the  penetrance is 0.8.
    Penetrance can be applied at any stage (default to DuringMating).
    It will be calculated during  mating, and then the affected status
    will be set for each offspring. Penetrance can also be used as
    PreMating, PostMating or even PrePostMating??? operator. In these
    cases, the affected status will be set to all individuals
    according to their  penetrance values. It is also possible to
    store  penetrance in a given information field specified by
    infoFields parameter (e.g. infoFields=['penetrance']). This is
    useful to check the  penetrance values at a later time.
    Affected status will be used for statistical purpose, and most
    importantly, ascertainment. They will be calculated along with
    fitness although they might not be used at every generation. You
    can use two operators: one for fitness/selection, active at every
    generation; one for affected status, active only at
    ascertainments, to avoid unnecessary calculation of the affected
    status.
    Pentrance values are used to set the affectedness of individuals,
    and are usually not saved. If you would like to know the
    penetrance value, you need to
    * use addInfoField('penetrance') to the  population to analyze.
    (Or use infoFields parameter of the  population constructor), and
    * use e.g.,  mlPenetrance(...., infoFields=['penetrance']) to add
    the  penetrance field to the  penetrance operator you use. You may
    choose a name other than 'penetrance' as long as the field names
    for the operator and  population match. Penetrance functions can
    be applied to the current, all, or certain number of ancestral
    generations. This is controlled by the ancestralGen parameter,
    which is default to -1 (all available ancestral generations). You
    can set it to 0 if you only need affection??? status for the
    current generation, or specify a number n for the number of
    ancestral generations (n + 1 total generations) to process. Note
    that ancestralGen parameter is ignored if the  penetrance operator
    is used as a during  mating operator.

"; 

%feature("docstring") simuPOP::penetrance::penetrance "

Description:

    create a  penetrance operator

Usage:

    penetrance(ancestralGen=-1, stage=DuringMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Details:

    default to be always active.

Details:

    default to be always active.

Arguments:

    ancestralGen:   if this parameter is set to be 0, then apply
                    penetrance to the current generation; if -1, apply
                    to all generations; otherwise, apply to the
                    specified number of ancestral generations
    stage:          specify the stage this operator will be applied,
                    default to DuringMating.
    infoFields:     If one field is specified, it will be used to
                    store  penetrance values.???

"; 

%feature("docstring") simuPOP::penetrance::~penetrance "

Description:

    destructor

Usage:

    x.~penetrance()

"; 

%feature("docstring") simuPOP::penetrance::clone "

Description:

    deep copy of a  penetrance operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::penetrance::penet "

Description:

    calculate/return  penetrance etc.

Usage:

    x.penet()

"; 

%feature("docstring") simuPOP::penetrance::apply "

Description:

    set  penetrance to all individuals and record  penetrance if
    requested

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::penetrance::applyDuringMating "

Description:

    set  penetrance to all individuals

Usage:

    x.applyDuringMating(pop, offspring, dad=None, mom=None)

"; 

%feature("docstring") simuPOP::penetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  penetrance operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pointMutator "

Function form:

    PointMutate

Function form:

    PointMutate

Description:

    point  mutator

Details:

    Mutate specified individuals at a specified loci to a spcified
    allele. I.e., this is a non-random  mutator used to introduce
    diseases etc.  pointMutator, as its name suggest, does point
    mutation. This  mutator will turn alleles at loci on the first
    chromosome copy to toAllele for  individualinds. You can specify
    atPloidy to mutate other, or all ploidy copies.

Details:

    Mutate specified individuals at a specified loci to a spcified
    allele. I.e., this is a non-random  mutator used to introduce
    diseases etc.  pointMutator, as its name suggest, does point
    mutation. This  mutator will turn alleles at loci on the first
    chromosome copy to toAllele for  individualinds. You can specify
    atPloidy to mutate other, or all ploidy copies.

"; 

%feature("docstring") simuPOP::pointMutator::pointMutator "

Description:

    create a  pointMutator

Usage:

    pointMutator(loci, toAllele, atPloidy=[], inds=[], output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Details:

    Please see  mutator for the description of other parameters.

Details:

    Please see  mutator for the description of other parameters.

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

%feature("docstring") simuPOP::population "

Description:

    a collection of individuals with the same genotypic structure

Details:

    A  simuPOP population consists of individuals of the same
    genotypic structure, which refers to the number of chromosomes,
    number and position of loci on each chromosome etc. The most
    important components of a  population are:
    * subpopulation. A  population is divided into subpopulations
    (unstructured  population has a single subpopulation, which is the
    whole  population itself). Subpopulation structure limits the
    usually random exchange of genotypes between individuals by
    disallowing  mating between individuals from different
    subpopulations. In the presence of subpopualtion structure,
    exchange of genetic information across subpopulations can only be
    done through migration. Note that in  simuPOP there is no sub-
    subpopulation or family in subpopulations.
    * variables. Every  population has its own variable space, or
    local namespace in  simuPOP term. This namespace is a Python
    dictionary that is attached to each  population and can be exposed
    to the users through  vars() or dvars() function. Many functions
    and operators work and store their results in this namespace. For
    example, function Stat sets variables such as alleleFreq[loc], and
    you can access it via pop.dvars().alleleFreq[loc][allele].
    * ancestral generations. A  population can save arbitrary number
    of ancestral generations. During evolution, the latest several (or
    all) ancestral generations are saved. Functions to switch between
    ancestral generations are provided so that one can examine and
    modify ancestral generations. Other important concepts like
    information fields are explained in class  individual.

Details:

    A  simuPOP population consists of individuals of the same
    genotypic structure, which refers to the number of chromosomes,
    number and position of loci on each chromosome etc. The most
    important components of a  population are:
    * subpopulation. A  population is divided into subpopulations
    (unstructured  population has a single subpopulation, which is the
    whole  population itself). Subpopulation structure limits the
    usually random exchange of genotypes between individuals by
    disallowing  mating between individuals from different
    subpopulations. In the presence of subpopualtion structure,
    exchange of genetic information across subpopulations can only be
    done through migration. Note that in  simuPOP there is no sub-
    subpopulation or family in subpopulations.
    * variables. Every  population has its own variable space, or
    local namespace in  simuPOP term. This namespace is a Python
    dictionary that is attached to each  population and can be exposed
    to the users through  vars() or dvars() function. Many functions
    and operators work and store their results in this namespace. For
    example, function Stat sets variables such as alleleFreq[loc], and
    you can access it via pop.dvars().alleleFreq[loc][allele].
    * ancestral generations. A  population can save arbitrary number
    of ancestral generations. During evolution, the latest several (or
    all) ancestral generations are saved. Functions to switch between
    ancestral generations are provided so that one can examine and
    modify ancestral generations. Other important concepts like
    information fields are explained in class  individual.

Note:

    Although a large number of member functions are provided, most of
    the operations are performed by operators. These functions will
    only be useful when you need to manipulate a  population
    explicitly.

"; 

%feature("docstring") simuPOP::population::population "

Description:

    Create a  population object with given size and genotypic
    structure.

Usage:

    population(size=0, ploidy=2, loci=[], sexChrom=False,
      lociPos=[], subPop=[], ancestralDepth=0, alleleNames=[],
      lociNames=[], maxAllele=MaxAllele, infoFields=[], chromMap=[])

Details:

    This is techniquely the __init__ function of the  population
    object.

Details:

    This is techniquely the __init__ function of the  population
    object.

Arguments:

    size:           population size. Can be ignored if subPop is
                    specified. In that case, size is the sum of
                    subPop. Default to 0.
    ploidy:         number of sets of homologous copies of
                    chromosomes. Default to 2 (diploid).
    loci:           an array of numbers of loci on each chromosome.
                    The length of parameter loci determines the number
                    of chromosomes. Default to [1], meaning one
                    chromosome with a single locus.
                    The last chromosome can be sex chromosome. In this
                    case, the maximum number of loci on X and Y should
                    be provided. I.e., if there are 3 loci on Y
                    chromosme and 5 on X chromosome, use 5.
    sexChrom:       true or false. Diploid  population only. If true,
                    the last homologous chromosome will be treated as
                    sex chromosome. (XY for male and XX for female.)
                    If X and Y have different number of loci, number
                    of loci of the longer one of the last (sex)
                    chromosome should be specified in loci.
    lociPos:        a 1-d or 2-d array specifying positions of loci on
                    each chromosome. You can use a nested array to
                    specify loci position for each chromosome. For
                    example, you can use lociPos=[1,2,3] when loci=[3]
                    or lociPos=[[1,2],[1.5,3,5]] for loci=[2,3].
                    simuPOP does not assume a unit for these
                    positions, although they are usually intepreted as
                    centiMorgans. The default values are 1, 2, etc. on
                    each chromosome.
    subPop:         an array of subpopulation sizes. Default value is
                    [size] which means a single subpopulation of the
                    whole  population. If both size and subPop are
                    provided, subPop should add up to size.
    ancestralDepth: number of most recent ancestral generations to
                    keep during evolution. Default to 0, which means
                    only the current generation will be available. You
                    can set it to a positive number m to keep the
                    latest m generations in the  population, or -1 to
                    keep all ancestral populations. Note that keeping
                    track of all ancestral generations may quickly
                    exhaust your computer RAM. If you really need to
                    do that, use  savePopulation operator to save each
                    generation to a file is a much better choice.
    alleleNames:    an array of allele names. For example, for a locus
                    with alleles A,C,T,G, you can specify alleleNames
                    as ('A','C','T','G').
    lociNames:      an array or a matrix (separated by chromosomes) of
                    names for each locus. Default to \"locX-Y\" where X
                    is chromosome index and Y is locus number, both
                    starting from 1.
    maxAllele:      maximum allele number. Default to the maximum
                    allowed allele state of the current library. This
                    will set a cap for all loci. For  individual
                    locus, you can specify maxAllele in mutation
                    models, which can be smaller than global maxAllele
                    but not larger. Note that this number is the
                    number of allele states minus 1 since allele
                    number starts from 0.
    infoFields:     name of information fields that will be attached
                    to each  individual. For example, if you need to
                    record the parents of each  individual using
                    operator \\c parentTagger(), you will need two
                    fields father_idx and mother_idx.

Example:

>>> # use of population function
>>> # a Wright-Fisher population
>>> WF = population(size=100, ploidy=1, loci=[1])
>>> 
>>> # a diploid population of size 10
>>> # there are two chromosomes with 5 and 7 loci respectively
>>> pop = population(size=10, ploidy=2, loci=[5, 7], subPop=[2, 8])
>>> 
>>> # a population with SNP markers (with names A,C,T,G)
>>> # range() are python functions
>>> pop = population(size=5, ploidy=2, loci=[5,10],
...     lociPos=[range(0,5),range(0,20,2)],
...     alleleNames=['A','C','T','G'],
...     subPop=[2,3], maxAllele=3)
>>> 
>>> #
>>> # population structure functions
>>> print pop.popSize()
5
>>> print pop.numSubPop()
2
>>> print pop.subPopSize(0)
2
>>> print pop.subPopSizes()
(2, 3)
>>> print pop.subPopBegin(1)
2
>>> print pop.subPopEnd(1)
5
>>> print pop.subPopIndPair(3)
(1, 1)
>>> print pop.absIndIndex(1,1)
3
>>> 
>>> #
>>> # functions of setting population structure
>>> pop.setIndSubPopID([1,2,2,3,1])
>>> pop.setSubPopByIndID()
>>> pop.removeLoci(keep=range(2,7))
>>> Dump(pop)
Ploidy:         	2
Number of chrom:	2
Number of loci: 	3 2 
Maximum allele state:	3
Loci positions: 
		2 3 4 
		0 2 
Loci names: 
		loc1-3 loc1-4 loc1-5 
		loc2-1 loc2-2 
population size:	5
Number of subPop:	4
Subpop sizes:   	0 2 2 1 
Number of ancestral populations:	0
individual info: 
sub population 1:
   0: MU AAA AA | AAA AA 
   1: MU AAA AA | AAA AA 
sub population 2:
   2: MU AAA AA | AAA AA 
   3: MU AAA AA | AAA AA 
sub population 3:
   4: MU AAA AA | AAA AA 
End of individual info.


No ancenstral population recorded.
>>> 
>>> #
>>> # save and load population
>>> # save it in various formats, default format is \"txt\"
>>> pop = population(1000, loci=[2, 5, 10])
>>> pop.savePopulation(\"pop.txt\")
>>> pop.savePopulation(\"pop.txt\", compress=False)
>>> pop.savePopulation(\"pop.xml\", format=\"xml\")
>>> pop.savePopulation(\"pop.bin\", format=\"bin\")
>>> 
>>> # load it in another population
>>> pop1 = LoadPopulation(\"pop.xml\", format=\"xml\")
>>>


"; 

%ignore simuPOP::population::population(const population &rhs);

%feature("docstring") simuPOP::population::clone "

Description:

    deep copy of a  population. (In python, pop1 = pop will only
    create a reference to pop.)

Usage:

    x.clone(keepAncestralPops=-1)

Details:

    This function by default copies all ancestral generations, but you
    can copy only one (current, keepAncestralPops=0), or specified
    number of ancestral generations.

Details:

    This function by default copies all ancestral generations, but you
    can copy only one (current, keepAncestralPops=0), or specified
    number of ancestral generations.

"; 

%feature("docstring") simuPOP::population::swap "

Description:

    swap the content of two populations

Usage:

    x.swap(rhs)

"; 

%feature("docstring") simuPOP::population::~population "

Description:

    destroy a  population

Usage:

    x.~population()

"; 

%feature("docstring") simuPOP::population::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  population

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::population::__cmp__ "

Description:

    a python function used to compare the  population objects

Usage:

    x.__cmp__(rhs)

"; 

%feature("docstring") simuPOP::population::setSubPopStru "

Description:

    set population/subpopulation structure given subpopulation sizes

Usage:

    x.setSubPopStru(newSubPopSizes, allowPopSizeChange=False)

Arguments:

    subPopSize:     an array of subpopulation sizes. The  population
                    may or may not change according to parameter
                    allowPopSizeChange if the sum of subPopSize does
                    not match popSize.
    allowPopSizeChange:if this parameter is true,  population will be
                    resized.

"; 

%feature("docstring") simuPOP::population::numSubPop "

Description:

    number of subpopulations in a  population

Usage:

    x.numSubPop()

"; 

%feature("docstring") simuPOP::population::subPopSize "

Description:

    return size of a subpopulation subPop

Usage:

    x.subPopSize(subPop)

Arguments:

    subPop:         index of subpopulation (start from 0)

"; 

%feature("docstring") simuPOP::population::subPopSizes "

Description:

    return an array of all subpopulation sizes

Usage:

    x.subPopSizes()

"; 

%feature("docstring") simuPOP::population::popSize "

Description:

    total  population size

Usage:

    x.popSize()

"; 

%feature("docstring") simuPOP::population::absIndIndex "

Description:

    return the absolute index of an  individual in a subpopulation

Usage:

    x.absIndIndex(ind, subPop)

Arguments:

    index:          index of an  individual in a subpopulation subPop
    subPop:         subpopulation index (start from 0)

"; 

%feature("docstring") simuPOP::population::subPopIndPair "

Description:

    return the (sp, idx) pair from an absolute index of an  individual

Usage:

    x.subPopIndPair(ind)

"; 

%feature("docstring") simuPOP::population::subPopBegin "

Description:

    index of the first  individual of a subpopulation subPop

Usage:

    x.subPopBegin(subPop)

"; 

%feature("docstring") simuPOP::population::subPopEnd "

Description:

    return the value of the index of the last  individual of a
    subpopulation subPop plus 1

Usage:

    x.subPopEnd(subPop)

"; 

%feature("docstring") simuPOP::population::ind "

Description:

    refernce to  individualind in subpopulation subPop

Usage:

    x.ind(ind, subPop=0)

Details:

    This function is named  individual in the Python interface.

Details:

    This function is named  individual in the Python interface.

Arguments:

    ind:            individual index within subPop
    subPop:         subpopulation index

"; 

%feature("docstring") simuPOP::population::individuals "

Description:

    return an iterator that can be used to iterate through all
    individuals

Usage:

    x.individuals()

Details:

    Typical usage is
    for ind in pop.individuals():

Details:

    Typical usage is
    for ind in pop.individuals():

"; 

%feature("docstring") simuPOP::population::individuals "

Description:

    return an iterator that can be used to iterate through all
    individuals in subpopulation subPop

Usage:

    x.individuals(subPop)

"; 

%ignore simuPOP::population::ind(ULONG ind, UINT subPop=0) const ;

%ignore simuPOP::population::shallowCopied();

%ignore simuPOP::population::setShallowCopied(bool s);

%ignore simuPOP::population::infoOrdered();

%ignore simuPOP::population::setInfoOrdered(bool s);

%ignore simuPOP::population::indBegin();

%ignore simuPOP::population::indEnd();

%ignore simuPOP::population::indBegin(UINT subPop);

%ignore simuPOP::population::indEnd(UINT subPop);

%ignore simuPOP::population::alleleBegin(UINT locus, bool order);

%ignore simuPOP::population::alleleEnd(UINT locus, bool order);

%ignore simuPOP::population::alleleBegin(UINT locus, UINT subPop, bool order);

%ignore simuPOP::population::genoBegin(bool order);

%ignore simuPOP::population::genoEnd(bool order);

%ignore simuPOP::population::genoBegin(UINT subPop, bool order);

%ignore simuPOP::population::genoEnd(UINT subPop, bool order);

%ignore simuPOP::population::indGenoBegin(ULONG ind) const ;

%ignore simuPOP::population::indGenoEnd(ULONG ind) const ;

%feature("docstring") simuPOP::population::arrGenotype "

Description:

    get the whole genotypes

Usage:

    x.arrGenotype(order)

Details:

    Return an editable array of all genotypes of the  population. You
    need to know how these genotypes are organized to safely
    read/write genotype directly.

Details:

    Return an editable array of all genotypes of the  population. You
    need to know how these genotypes are organized to safely
    read/write genotype directly.

Arguments:

    order:          if order is true, individuals will be ordered such
                    that pop.individual(x). arrGenotype() ==
                    pop.arrGenotype()[x*pop. genoSize():(x+1)*pop.
                    genoSize()].

"; 

%feature("docstring") simuPOP::population::arrGenotype "

Description:

    get the whole genotypes of individuals in a subpopulation

Usage:

    x.arrGenotype(subPop, order)

Details:

    Return an editable array of all genotype in a subpopulation.

Details:

    Return an editable array of all genotype in a subpopulation.

Arguments:

    subPop:         index of subpopulation (start from 0)
    order:          if order is true, individuals will be ordered.

"; 

%feature("docstring") simuPOP::population::setIndSubPopID "

Description:

    set subpopulation ID with given ID

Usage:

    x.setIndSubPopID(id)

Details:

    Set subpopulation ID of each  individual with given ID.
    Individuals can be rearranged afterwards using setSubPopByIndID.

Details:

    Set subpopulation ID of each  individual with given ID.
    Individuals can be rearranged afterwards using setSubPopByIndID.

Arguments:

    id:             an array of the same length of  population size,
                    resprenting subpopulation ID of each  individual.

"; 

%feature("docstring") simuPOP::population::setIndSubPopIDWithID "

Description:

    set subpopulation ID of each  individual with their current
    subpopulation ID

Usage:

    x.setIndSubPopIDWithID()

"; 

%feature("docstring") simuPOP::population::setSubPopByIndID "

Description:

    adjust subpopulation according to  individual subpopulation ID.

Usage:

    x.setSubPopByIndID(id=[])

Details:

    Rearrange individuals to their new subpopulations according to
    their subpopulation ID (or the new given ID). Order within each
    subpopulation is not respected.

Details:

    Rearrange individuals to their new subpopulations according to
    their subpopulation ID (or the new given ID). Order within each
    subpopulation is not respected.

Arguments:

    id:             new subpopulation ID, if given, current
                    individual subpopulation ID will be ignored.

Note:

    Individual with negative info will be removed!

"; 

%feature("docstring") simuPOP::population::splitSubPop "

Description:

    split a subpopulation into subpopulations of given sizes

Usage:

    x.splitSubPop(which, sizes, subPopID=[])

Details:

    The sum of given sizes should be equal to the size of the split
    subpopulation. Subpopulation IDs can be specified. The
    subpopulation IDs of non-split subpopulations will be kept. For
    example, if subpopulation 1 of 0 1 2 3 is split into three parts,
    the new subpop id will be 0 (1 4 5) 2 3.

Details:

    The sum of given sizes should be equal to the size of the split
    subpopulation. Subpopulation IDs can be specified. The
    subpopulation IDs of non-split subpopulations will be kept. For
    example, if subpopulation 1 of 0 1 2 3 is split into three parts,
    the new subpop id will be 0 (1 4 5) 2 3.

Note:

    subpop with negative ID will be removed. So, you can shrink one
    subpop by splitting and setting one of the new subpop with
    negative ID.

"; 

%feature("docstring") simuPOP::population::splitSubPopByProportion "

Description:

    split a subpopulation into subpopulations of given proportions

Usage:

    x.splitSubPopByProportion(which, proportions, subPopID=[])

Details:

    The sum of given proportions should add up to one. Subpopulation
    IDs can be specified.

Details:

    The sum of given proportions should add up to one. Subpopulation
    IDs can be specified.

Note:

    subpop with negative ID will be removed. So, you can shrink one
    subpop by splitting and setting one of the new subpop with
    negative ID.

"; 

%feature("docstring") simuPOP::population::removeEmptySubPops "

Description:

    remove empty subpopulations by adjusting subpopulation IDs

Usage:

    x.removeEmptySubPops()

"; 

%feature("docstring") simuPOP::population::removeSubPops "

Description:

    remove subpopulations and adjust subpopulation IDs so that there
    will be no 'empty' subpopulation left

Usage:

    x.removeSubPops(subPops=[], shiftSubPopID=True,
      removeEmptySubPops=False)

Details:

    Remove specified subpopulations (and all individuals within). If
    shiftSubPopID is false, subPopID will be kept intactly.

Details:

    Remove specified subpopulations (and all individuals within). If
    shiftSubPopID is false, subPopID will be kept intactly.

"; 

%feature("docstring") simuPOP::population::removeIndividuals "

Description:

    remove individuals. If a valid subPop is given, remove individuals
    from this subpopulation. Indexes in inds will be treated as
    relative indexes.

Usage:

    x.removeIndividuals(inds=[], subPop=-1,
      removeEmptySubPops=False)

"; 

%feature("docstring") simuPOP::population::mergeSubPops "

Description:

    merge given subpopulations

Usage:

    x.mergeSubPops(subPops=[], removeEmptySubPops=False)

Details:

    Merge subpopulations, the first subpopulation ID (the first one in
    array subPops) will be used as the ID of the new subpopulation.
    That is to say, all merged subpopulations will take the ID of the
    first one. The subpopulation ID of the empty subpopulations will
    be kept (so that other subpopulations are unaffected, unless they
    are removed by removeEmptySubPops = True).

Details:

    Merge subpopulations, the first subpopulation ID (the first one in
    array subPops) will be used as the ID of the new subpopulation.
    That is to say, all merged subpopulations will take the ID of the
    first one. The subpopulation ID of the empty subpopulations will
    be kept (so that other subpopulations are unaffected, unless they
    are removed by removeEmptySubPops = True).

"; 

%feature("docstring") simuPOP::population::mergePopulation "

Description:

    merge populations by individuals

Usage:

    x.mergePopulation(pop, newSubPopSizes=[], keepAncestralPops=-1)

Details:

    Merge individuals from pop to the current  population. Two
    populations should have the same genotypic structures. By default,
    subpopulations of the merged populations are kept. I.e., if you
    merge two populations with one subpopulation, the resulting
    population will have two subpopulations. All ancestral generations
    are also merged by default.

Details:

    Merge individuals from pop to the current  population. Two
    populations should have the same genotypic structures. By default,
    subpopulations of the merged populations are kept. I.e., if you
    merge two populations with one subpopulation, the resulting
    population will have two subpopulations. All ancestral generations
    are also merged by default.

Arguments:

    newSubPopSizes: subpopulation sizes can be specified. The overall
                    size should be the combined size of the two
                    populations. Because this parameter will be used
                    for all ancestral generations, it may fail if
                    ancestral generations have different sizes. To
                    avoid this problem, you can run mergePopulation
                    without this parameter, and then adjust
                    subpopulation sizes generation by generation.
    keepAncestralPops:ancestral populations to merge, default to all
                    (-1)

Note:

    Population variables are not copied to pop.

"; 

%feature("docstring") simuPOP::population::mergePopulationByLoci "

Description:

    merge populations by loci

Usage:

    x.mergePopulationByLoci(pop, newNumLoci=[], newLociPos=[],
      byChromosome=False)

Details:

    Two populations should have the same number of individuals. This
    also holds for any ancestral generations. By default, chromosomes
    of pop are appended to the current  population. You can change
    this arrangement in two ways
    * specify new chromosome structure using parameter newLoci and
    newLociPos. Loci from new and old populations are still in their
    original order, but chromosome number and positions can be changed
    in this way.
    * specify byChromosome=true so that chromosomes will be merged one
    by one. In this case, loci position of two popualtions are
    important because loci will be arranged in the order of loci
    position; and identical loci position of two loci in two
    populations will lead to error.

Details:

    Two populations should have the same number of individuals. This
    also holds for any ancestral generations. By default, chromosomes
    of pop are appended to the current  population. You can change
    this arrangement in two ways
    * specify new chromosome structure using parameter newLoci and
    newLociPos. Loci from new and old populations are still in their
    original order, but chromosome number and positions can be changed
    in this way.
    * specify byChromosome=true so that chromosomes will be merged one
    by one. In this case, loci position of two popualtions are
    important because loci will be arranged in the order of loci
    position; and identical loci position of two loci in two
    populations will lead to error.

Arguments:

    newNumLoci:     the new number of loci for the combined genotypic
                    structure.
    newLociPos:     the new loci position if number of loci on each
                    chromosomes are changed with newNumLoci. New loci
                    positions should be in order on the new
                    chromosomes.
    byChromosome:   merge chromosome by chromosome, loci are ordered
                    by loci position Default to False.

Note:

    * Information fields are not merged.
    * All ancestral generations are merged because all individuals in
    a  population have to have the same genotypic structure.

"; 

%feature("docstring") simuPOP::population::insertBeforeLoci "

Description:

    insert loci at given locations

Usage:

    x.insertBeforeLoci(idx, pos, names=[])

Details:

    Insert loci at some given locations. Alleles at inserted loci are
    initialized with zero allele.

Details:

    Insert loci at some given locations. Alleles at inserted loci are
    initialized with zero allele.

Arguments:

    idx:            an array of locus index. The loci will be inserted
                    before each index. If you need to append to the
                    last locus, use insertAfterLoci instead. If your
                    index is the first locus of a chromosome, the
                    inserted locus will become the first locus of that
                    chromosome. If you need to insert multiple loci
                    before a locus, repeat that locus number.
    pos:            an array of locus positions. The positions of the
                    appended loci have to be between adjacent markers.
    names:          an array of locus names. If this parameter is not
                    given, some unique names such as \"insX_Y\" will be
                    given.

"; 

%feature("docstring") simuPOP::population::insertBeforeLocus "

Description:

    insert an locus at given location.

Usage:

    x.insertBeforeLocus(idx, pos, name=string)

Details:

    insertBeforeLocus(idx, pos, name) is a shortcut to
    insertBeforeLoci([idx], [pos], [name])

Details:

    insertBeforeLocus(idx, pos, name) is a shortcut to
    insertBeforeLoci([idx], [pos], [name])

"; 

%feature("docstring") simuPOP::population::insertAfterLoci "

Description:

    append loci at given locations

Usage:

    x.insertAfterLoci(idx, pos, names=[])

Details:

    Append loci at some given locations. Alleles at inserted loci are
    initialized with zero allele.

Details:

    Append loci at some given locations. Alleles at inserted loci are
    initialized with zero allele.

Arguments:

    idx:            an array of locus index. The loci will be added
                    after each index. If you need to append to the
                    first locus of a chromosome, use insertBeforeLoci
                    instead. If your index is the last locus of a
                    chromosome, the appended locus will become the
                    last locus of that chromosome. If you need to
                    append multiple loci after a locus, repeat that
                    locus number.
    pos:            an array of locus positions. The positions of the
                    appended loci have to be between adjacent markers.
    names:          an array of locus names. If this parameter is not
                    given, some unique names such as \"insX_Y\" will be
                    given.

"; 

%feature("docstring") simuPOP::population::insertAfterLocus "

Description:

    append an locus at a given location

Usage:

    x.insertAfterLocus(idx, pos, name=string)

Details:

    insertAfterLocus(idx, pos, name) is a shortcut to
    insertAfterLoci([idx], [pos], [name]).

Details:

    insertAfterLocus(idx, pos, name) is a shortcut to
    insertAfterLoci([idx], [pos], [name]).

"; 

%feature("docstring") simuPOP::population::resize "

Description:

    resize  population

Usage:

    x.resize(newSubPopSizes, propagate=False)

Details:

    Resize  population by giving new subpopulation sizes.

Details:

    Resize  population by giving new subpopulation sizes.

Arguments:

    newSubPopSizes: an array of new subpopulation sizes. If there is
                    only one subpopulation, use [newPopSize].
    propagate:      if propagate is true, copy individuals to new
                    comers. I.e., 1, 2, 3 ==> 1, 2, 3, 1, 2, 3, 1

Note:

    This function only resizes the current generation.

"; 

%feature("docstring") simuPOP::population::reorderSubPops "

Description:

    reorder subpopulations by order or by rank

Usage:

    x.reorderSubPops(order=[], rank=[], removeEmptySubPops=False)

Arguments:

    order:          new order of the subpopulations. For examples, 3 2
                    0 1 means subpop3, subpop2, subpop0, subpop1 will
                    be the new layout.
    rank:           you may also specify a new rank for each
                    subpopulation. For example, 3,2,0,1 means the
                    original subpopulations will have new IDs 3,2,0,1,
                    respectively. To achive order 3,2,0,1, the rank
                    should be 1 0 2 3.

"; 

%feature("docstring") simuPOP::population::newPopByIndID "

Usage:

    x.newPopByIndID(keepAncestralPops=-1, id=[],
      removeEmptySubPops=False)

Details:

    Form a new  population according to  individual subpopulation ID.
    Individuals with negative subpopulation ID will be removed.

Details:

    Form a new  population according to  individual subpopulation ID.
    Individuals with negative subpopulation ID will be removed.

"; 

%feature("docstring") simuPOP::population::removeLoci "

Description:

    remove some loci from the current  population. Only one of the two
    parameters can be specified.

Usage:

    x.removeLoci(remove=[], keep=[])

"; 

%feature("docstring") simuPOP::population::newPopWithPartialLoci "

Description:

    obtain a new  population with selected loci

Usage:

    x.newPopWithPartialLoci(remove=[], keep=[])

Details:

    Copy current  population to a new one with selected loci keep or
    remove specified loci remove (no change on the current
    population), equivalent to
    y=x.clone
    y.removeLoci(remove, keep)

Details:

    Copy current  population to a new one with selected loci keep or
    remove specified loci remove (no change on the current
    population), equivalent to
    y=x.clone
    y.removeLoci(remove, keep)

"; 

%feature("docstring") simuPOP::population::pushAndDiscard "

Description:

    absorb rhs population as the current generation of a  population

Usage:

    x.pushAndDiscard(rhs, force=False)

Details:

    This function is used by a  simulator to push offspring generation
    rhs to the current  population, while the current  population is
    pushed back as an ancestral  population (if ancestralDepath() !=
    0). Because rhs population is swapped in, rhs will be empty after
    this operation.

Details:

    This function is used by a  simulator to push offspring generation
    rhs to the current  population, while the current  population is
    pushed back as an ancestral  population (if ancestralDepath() !=
    0). Because rhs population is swapped in, rhs will be empty after
    this operation.

"; 

%feature("docstring") simuPOP::population::ancestralDepth "

Description:

    ancestral depth of the current  population

Usage:

    x.ancestralDepth()

Note:

    The returned value is the number of ancestral generations exist in
    the  population, not necessarily equals to the number set by
    setAncestralDepth().

"; 

%feature("docstring") simuPOP::population::ancestralGen "

Description:

    currently used ancestral  population (0 for the latest generation)

Usage:

    x.ancestralGen()

Details:

    Current ancestral  population activated by  useAncestralPop().
    There can be several ancestral generations in a  population. 0
    (current), 1 (parental) etc. When useAncestralPop(gen) is used,
    current generation is set to one of the parental generations,
    which is the information returned by this function.
    useAncestralPop(0) should always be used to set a  population to
    its usual ancestral order after operations to the ancestral
    generation are done.

Details:

    Current ancestral  population activated by  useAncestralPop().
    There can be several ancestral generations in a  population. 0
    (current), 1 (parental) etc. When useAncestralPop(gen) is used,
    current generation is set to one of the parental generations,
    which is the information returned by this function.
    useAncestralPop(0) should always be used to set a  population to
    its usual ancestral order after operations to the ancestral
    generation are done.

"; 

%feature("docstring") simuPOP::population::setIndInfo "

Description:

    set  individual information for the given information field
    (index),

Usage:

    x.setIndInfo(values, idx, order=True)

Arguments:

    values:         an array that has the same length as  population
                    size.
    idx:            index to the information field.
    order:          if true, info will be in the order of individuals

"; 

%feature("docstring") simuPOP::population::setIndInfo "

Description:

    set  individual information for the given information field (name)

Usage:

    x.setIndInfo(values, name, order=True)

Details:

    x.setIndInfo(values, name) is equivalent to the idx version
    x.setIndInfo(values, x.infoIdx(name)).

Details:

    x.setIndInfo(values, name) is equivalent to the idx version
    x.setIndInfo(values, x.infoIdx(name)).

"; 

%ignore simuPOP::population::infoBegin(UINT idx, bool order);

%ignore simuPOP::population::infoEnd(UINT idx, bool order);

%ignore simuPOP::population::infoBegin(UINT index, UINT subPop, bool order);

%feature("docstring") simuPOP::population::indInfo "

Description:

    get information field idx of all individuals

Usage:

    x.indInfo(idx, order)

Arguments:

    idx:            index of the information field
    order:          if true, sort returned vector in  individual order

"; 

%feature("docstring") simuPOP::population::indInfo "

Description:

    get information field name of all individuals

Usage:

    x.indInfo(name, order)

Arguments:

    name:           name of the information field
    order:          if true, sort returned vector in  individual order

"; 

%feature("docstring") simuPOP::population::indInfo "

Description:

    get information field idx of all individuals in a subpopulation
    subPop

Usage:

    x.indInfo(idx, subPop, order)

Arguments:

    idx:            index of the information field
    subPop:         subpopulation index
    order:          if true, sort returned vector in  individual order

"; 

%feature("docstring") simuPOP::population::indInfo "

Description:

    get information field name of all individuals in a subpopulation
    subPop

Usage:

    x.indInfo(name, subPop, order)

Arguments:

    name:           name of the information field
    subPop:         subpopulation index
    order:          if true, sort returned vector in  individual order

"; 

%feature("docstring") simuPOP::population::arrIndInfo "

Description:

    get an editable array (Python list) of all information fields

Usage:

    x.arrIndInfo(order)

Details:

    The length of the array is  infoSize()*popSize().

Details:

    The length of the array is  infoSize()*popSize().

Arguments:

    order:          whether or not the list has the same order as
                    individuals

"; 

%feature("docstring") simuPOP::population::arrIndInfo "

Description:

    get an editable array (Python list) of all information fields

Usage:

    x.arrIndInfo(subPop, order)

Arguments:

    order:          whether or not the list has the same order as
                    individuals

"; 

%feature("docstring") simuPOP::population::addInfoField "

Description:

    add an information field to a  population

Usage:

    x.addInfoField(field, init=0)

Arguments:

    field:          new information field. If it already exists, it
                    will be re-initialized.
    init:           initial value for the new field.

"; 

%feature("docstring") simuPOP::population::addInfoFields "

Description:

    add one or more information fields to a  population

Usage:

    x.addInfoFields(fields, init=0)

Arguments:

    fields:         an array of new information fields. If one or more
                    of the fields alreay exist, they will be re-
                    initialized.
    init:           initial value for the new fields.

"; 

%feature("docstring") simuPOP::population::setInfoFields "

Description:

    set information fields for an existing  population. The existing
    fields will be removed.

Usage:

    x.setInfoFields(fields, init=0)

Arguments:

    fields:         an array of fields
    init:           initial value for the new fields.

"; 

%feature("docstring") simuPOP::population::setAncestralDepth "

Description:

    set ancestral depth.

Usage:

    x.setAncestralDepth(depth)

Arguments:

    depth:          0 for none, -1 for unlimited, a positive number
                    sets the number of ancestral generations to save.

"; 

%feature("docstring") simuPOP::population::useAncestralPop "

Description:

    use an ancestral generation. 0 for the latest generation.

Usage:

    x.useAncestralPop(idx)

Arguments:

    idx:            Index of the ancestral generation. 0 for current,
                    1 for parental, etc. idx can not exceed ancestral
                    depth (see  setAncestralDepth).

"; 

%ignore simuPOP::population::equalTo(const population &rhs);

%ignore simuPOP::population::adjustGenoPosition(bool order);

%ignore simuPOP::population::adjustInfoPosition(bool order);

%feature("docstring") simuPOP::population::savePopulation "

Description:

    save  population to a file

Usage:

    x.savePopulation(filename, format=\"auto\", compress=True)

Arguments:

    filename:       save to filename
    format:         format to save. Can be one of the following:
                    'txt', 'bin', or 'xml', or 'auto' which is
                    determined by the extension of filename.

"; 

%ignore simuPOP::population::loadPopulation(const string &filename, const string &format="auto");

%ignore simuPOP::population::selectionOn();

%ignore simuPOP::population::selectionOn(UINT sp);

%feature("docstring") simuPOP::population::turnOffSelection "

Description:

    turn off selection for all subpopulations

Usage:

    x.turnOffSelection()

Details:

    This is only used when you would like to apply two selectors.
    Maybe using two different information fields.

Details:

    This is only used when you would like to apply two selectors.
    Maybe using two different information fields.

"; 

%ignore simuPOP::population::turnOnSelection(UINT sp);

%ignore simuPOP::population::turnOnSelection();

%feature("docstring") simuPOP::population::rep "

Description:

    current replicate in a  simulator which is not meaningful for a
    stand-alone  population

Usage:

    x.rep()

"; 

%ignore simuPOP::population::setRep(int rep, bool setVar=true);

%feature("docstring") simuPOP::population::grp "

Description:

    current group ID in a  simulator which is not meaningful for a
    stand-alone  population.

Usage:

    x.grp()

"; 

%ignore simuPOP::population::setGrp(int grp, bool setVar=true);

%feature("docstring") simuPOP::population::gen "

Description:

    current generation during evolution

Usage:

    x.gen()

"; 

%ignore simuPOP::population::setGen(ULONG gen, bool setVar=true);

%feature("docstring") simuPOP::population::vars "

Description:

    return variables of a  population. If subPop is given, return a
    dictionary for specified subpopulation.

Usage:

    x.vars(subPop=-1)

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

%feature("docstring") simuPOP::population::evaluate "

Description:

    evaluate a python statment/expression in the population's local
    namespace

Usage:

    x.evaluate(expr=\"\", stmts=\"\")

Details:

    This function evaluates a Python statment( stmts )/expression(
    expr ) and return its result as a string. Optionally run
    statement( stmts ) first.

Details:

    This function evaluates a Python statment( stmts )/expression(
    expr ) and return its result as a string. Optionally run
    statement( stmts ) first.

"; 

%feature("docstring") simuPOP::population::execute "

Description:

    execute a statement (can be a multi-line string) in the
    population's local namespace

Usage:

    x.execute(stmts=\"\")

"; 

%ignore simuPOP::population::rearrangeLoci(const vectoru &newNumLoci, const vectorf &newLociPos);

%feature("docstring") simuPOP::pyEval "

Function form:

    PyEval

Function form:

    PyEval

Description:

    evaluate an expression

Details:

    Python expressions/statements will be executed when  pyEval is
    applied to a  population by using parameters expr/stmts.
    Statements can also been executed when  pyEval is created and
    destroyed or before expr is executed. The corresponding parameters
    are preStmts, postStmts and stmts. For example, operator
    varPlotter uses this feature to initialize R plots and save plots
    to a file when finished.

Details:

    Python expressions/statements will be executed when  pyEval is
    applied to a  population by using parameters expr/stmts.
    Statements can also been executed when  pyEval is created and
    destroyed or before expr is executed. The corresponding parameters
    are preStmts, postStmts and stmts. For example, operator
    varPlotter uses this feature to initialize R plots and save plots
    to a file when finished.

"; 

%feature("docstring") simuPOP::pyEval::pyEval "

Description:

    evaluate expressions/statments in the local namespace of a
    replicate

Usage:

    pyEval(expr=\"\", stmts=\"\", preStmts=\"\", postStmts=\"\",
      exposePop=False, name=\"\", output=\">\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    expr:           the expression to be evaluated. Its result will be
                    sent to output.
    stmts:          the statement that will be executed before the
                    expression
    preStmts:       the statement that will be executed when the
                    operator is constructed
    postStmts:      the statement that will be executed when the
                    operator is destroyed
    exposePop:      if true, expose current  population as variable
                    pop
    name:           used to let pure Python operator to identify
                    themselves
    output:         default to >. I.e., output to standard output.

"; 

%feature("docstring") simuPOP::pyEval::~pyEval "

Description:

    simuPOP::pyEval::~pyEval

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

Details:

    The name of a  pyEval operator is given by an optional parameter
    name. It can be used to identify this  pyEval operator in debug
    output, or in the dryrun mode of  simulator::evolve.

"; 

%feature("docstring") simuPOP::pyExec "

Function form:

    PyExec

Function form:

    PyExec

Description:

    execute a Python statement

Details:

    This operator takes a list of statements and execute them. No
    value will be returned or outputted.

Details:

    This operator takes a list of statements and execute them. No
    value will be returned or outputted.

"; 

%feature("docstring") simuPOP::pyExec::pyExec "

Description:

    evaluate statments in the local replicate namespace, no return
    value

Usage:

    pyExec(stmts=\"\", preStmts=\"\", postStmts=\"\", exposePop=False,
      name=\"\", output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
      end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    stmts:          the statements (a single or multi-line string)
                    that will be executed when this operator is
                    applied.
    preStmts:       the statement that will be executed when the
                    operator is constructed
    postStmts:      the statement that will be executed when the
                    operator is destroyed
    exposePop:      if true, expose current  population as variable
                    pop
    default:        to >. I.e., output to standard output.

"; 

%feature("docstring") simuPOP::pyExec::~pyExec "

Description:

    simuPOP::pyExec::~pyExec

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

%feature("docstring") simuPOP::pyIndOperator "

Description:

    individual operator

Details:

    This operator is similar to a  pyOperator but works at the
    individual level. It expects a function that accepts an
    individual, optional genotype at certain loci, and an optional
    parameter. When it is applied, it passes each  individual to this
    function. When infoFields is given, this function should return an
    array to fill these infoFields. Otherwise, True/False is
    expected.More specifically, func can be
    * func(ind) when neither loci nor param is given.
    * func(ind, genotype) when loci is given
    * func(ind, param) when param is given
    * func(ind, genotype, param) when both loci and param are given.

Details:

    This operator is similar to a  pyOperator but works at the
    individual level. It expects a function that accepts an
    individual, optional genotype at certain loci, and an optional
    parameter. When it is applied, it passes each  individual to this
    function. When infoFields is given, this function should return an
    array to fill these infoFields. Otherwise, True/False is
    expected.More specifically, func can be
    * func(ind) when neither loci nor param is given.
    * func(ind, genotype) when loci is given
    * func(ind, param) when param is given
    * func(ind, genotype, param) when both loci and param are given.

"; 

%feature("docstring") simuPOP::pyIndOperator::pyIndOperator "

Description:

    a Pre- or PostMating Python operator that apply a function to each
    individual

Usage:

    pyIndOperator(func, loci=[], param=None, stage=PostMating,
      formOffGenotype=False, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    func:           a Python function that accepts an  individual and
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

%feature("docstring") simuPOP::pyIndOperator::clone "

Description:

    deep copy of a  pyIndOperator operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyIndOperator::apply "

Description:

    apply the  pyIndOperator operator to one  population

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

Function form:

    PyInit

Description:

    a hybrid  initializer???

Details:

    pyInit is a hybrid  initializer. User should define a function
    with parameters allele, ploidy and subpopulation indices, and
    return an allele value. Users of this operator must supply a
    Python function with parameter (index, ploidy, subpop). This
    operator will loop through all  individual in each subpopulation
    and call this function to initialize populations. The arrange of
    parameters allows different initialization scheme for each subpop.

Details:

    pyInit is a hybrid  initializer. User should define a function
    with parameters allele, ploidy and subpopulation indices, and
    return an allele value. Users of this operator must supply a
    Python function with parameter (index, ploidy, subpop). This
    operator will loop through all  individual in each subpopulation
    and call this function to initialize populations. The arrange of
    parameters allows different initialization scheme for each subpop.

"; 

%feature("docstring") simuPOP::pyInit::pyInit "

Description:

    initialize populations using given user function

Usage:

    pyInit(func, subPop=[], loci=[], atPloidy=-1, indRange=[],
      maleFreq=0.5, sex=[], stage=PreMating, begin=0, end=1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    func:           a Python function with parameter (index, ploidy,
                    subpop), where
                    * index is the allele index ranging from 0 to
                    totNumLoci(-1),
                    * ploidy is the index of the copy of chromosomes)
                    * subpop is the subpopulation index. The return
                    value of this function should be an integer.
    loci:           a vector of loci indices. If empty, apply to all
                    loci.
    locus:          a shortcut to loci
    atPloidy:       initialize which copy of chromosomes. Default to
                    all.
    stage:          default to PreMating

Example:

>>> def initAllele(ind, p, sp):
...   return sp + ind + p
... 
>>> simu = simulator( 
...     population(subPop=[2,3], loci=[5,7]),
...     randomMating(), rep=1)
>>> simu.step([
...     pyInit(func=initAllele),
...     dumper(alleleOnly=True, dispWidth=2)])
individual info: 
sub population 0:
   0: MU   1  2  3  4  5   6  7  8  9 10 11 12 |   0  1  2  3  4   5  6  7  8  9 10 11 
   1: FU   1  2  3  4  5   5  6  7  8  9 10 11 |   1  2  3  4  5   6  7  8  9 10 11 12 
sub population 1:
   2: FU   1  2  3  4  5   7  8  9 10 11 12 13 |   1  2  3  4  5   6  7  8  9 10 11 12 
   3: FU   1  2  3  4  5   6  7  8  9 10 11 12 |   1  2  3  4  5   7  8  9 10 11 12 13 
   4: FU   2  3  4  5  6   6  7  8  9 10 11 12 |   1  2  3  4  5   7  8  9 10 11 12 13 
End of individual info.


No ancenstral population recorded.
True
>>>


"; 

%feature("docstring") simuPOP::pyInit::~pyInit "

Description:

    simuPOP::pyInit::~pyInit

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

    apply operator  pyInit???

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyMating "

Description:

    a Python  mating scheme

Details:

    Hybird  mating scheme. This  mating scheme takes a Python function
    that accepts both the parental and offspring populations and this
    function is responsible for setting genotype, sex of the offspring
    generation. During-mating operators, if needed, have to be applied
    from this function as well. Note that the subpopulaton size
    parameters are honored and the passed offspring generation has the
    desired (sub) population sizes. Parameters that control the number
    of offspring of each family are ignored.
    This is likely an extremely slow  mating scheme and should be used
    for experimental uses only. When a  mating scheme is tested, it is
    recommended to implement it at the C++ level.

Details:

    Hybird  mating scheme. This  mating scheme takes a Python function
    that accepts both the parental and offspring populations and this
    function is responsible for setting genotype, sex of the offspring
    generation. During-mating operators, if needed, have to be applied
    from this function as well. Note that the subpopulaton size
    parameters are honored and the passed offspring generation has the
    desired (sub) population sizes. Parameters that control the number
    of offspring of each family are ignored.
    This is likely an extremely slow  mating scheme and should be used
    for experimental uses only. When a  mating scheme is tested, it is
    recommended to implement it at the C++ level.

"; 

%feature("docstring") simuPOP::pyMating::pyMating "

Description:

    create a Python  mating scheme

Usage:

    pyMating(func=None, newSubPopSize=[], newSubPopSizeExpr=\"\",
      newSubPopSizeFunc=None)

Arguments:

    func:           a Python function that accepts two parameters: the
                    parental and the offspring populations. The
                    offspring  population is empty, and this function
                    is responsible for setting genotype, sex etc. of
                    individuals in the offspring generation.

"; 

%feature("docstring") simuPOP::pyMating::~pyMating "

Description:

    destructor

Usage:

    x.~pyMating()

"; 

%feature("docstring") simuPOP::pyMating::clone "

Description:

    deep copy of a Python  mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::pyMating::pyMating(const pyMating &rhs);

%feature("docstring") simuPOP::pyMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the Python  mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::pyMating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%feature("docstring") simuPOP::pyMigrator "

Description:

    a more flexible Python  migrator

Details:

    This  migrator can be used in two ways
    * define a function that accepts a generation number and returns a
    migration rate matrix. This can be used in the varying migration
    rate cases.
    * define a function that accepts individuals etc, and returns the
    new subpopulation ID. More specifically, func can be
    * func(ind) when neither loci nor param is given.
    * func(ind, genotype) when loci is given.
    * func(ind, param) when param is given.
    * func(ind, genotype, param) when both loci and param are given.

Details:

    This  migrator can be used in two ways
    * define a function that accepts a generation number and returns a
    migration rate matrix. This can be used in the varying migration
    rate cases.
    * define a function that accepts individuals etc, and returns the
    new subpopulation ID. More specifically, func can be
    * func(ind) when neither loci nor param is given.
    * func(ind, genotype) when loci is given.
    * func(ind, param) when param is given.
    * func(ind, genotype, param) when both loci and param are given.

"; 

%feature("docstring") simuPOP::pyMigrator::pyMigrator "

Description:

    create a hybrid  migrator

Usage:

    pyMigrator(rateFunc=None, mode=MigrByProbability, fromSubPop=[],
      toSubPop=[], indFunc=None, loci=[], param=None, stage=PreMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[])

Arguments:

    rateFunc:       a Python function that accepts a generation
                    number, current subpopulation sizes, and returns a
                    migration rate matrix. The  migrator then migrate
                    like a usual  migrator.
    indFunc:        a Python function that accepts an  individual,
                    optional genotype and parameter, then returns a
                    subpopulation id. This method can be used to
                    separate a  population according to  individual
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

Function form:

    PyMutate

Description:

    hybrid  mutator

Details:

    Hybrid  mutator. Mutation rate etc. are set just like others and
    you are supposed to provide a Python function to return a new
    allele state given an old state.  pyMutator will choose an allele
    as usual and call your function to mutate it to another allele.

Details:

    Hybrid  mutator. Mutation rate etc. are set just like others and
    you are supposed to provide a Python function to return a new
    allele state given an old state.  pyMutator will choose an allele
    as usual and call your function to mutate it to another allele.

"; 

%feature("docstring") simuPOP::pyMutator::pyMutator "

Description:

    create a  pyMutator

Usage:

    pyMutator(rate=[], loci=[], maxAllele=0, func=None, output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Example:

>>> def mut(x):
...   return 8
... 
>>> simu = simulator(population(size=3, loci=[3,5]), noMating())
>>> simu.step([
...   pyMutator(rate=.5, loci=[3,4,5], func=mut),
...   dumper(alleleOnly=True)])
individual info: 
sub population 0:
   0: MU   0  0  0   0  8  8  0  0 |   0  0  0   8  0  8  0  0 
   1: MU   0  0  0   0  0  8  0  0 |   0  0  0   0  0  8  0  0 
   2: MU   0  0  0   0  8  8  0  0 |   0  0  0   8  0  0  0  0 
End of individual info.


No ancenstral population recorded.
True
>>>


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

    the one and only Python operator???

Details:

    This operator accepts a function that can take the form of
    * func(pop) when stage=PreMating or PostMating, without setting
    param;
    * func(pop, param) when stage=PreMating or PostMating, with param;
    * func(pop, off, dad, mom) when stage=DuringMating and
    passOffspringOnly=False, without setting param;
    * func(off) when stage=DuringMating and passOffspringOnly=True,
    and without setting param;
    * func(pop, off, dad, mom, param) when stage=DuringMating and
    passOffspringOnly=False, with param;
    * func(off, param) when stage=DuringMating and
    passOffspringOnly=True, with param. For Pre- and PostMating
    usages, a  population and an optional parameter is passed to the
    given function. For DuringMating usages,  population, offspring,
    its parents and an optional parameter are passed to the given
    function. Arbitrary operations can be applied to the  population
    and offspring (if stage=DuringMating).

Details:

    This operator accepts a function that can take the form of
    * func(pop) when stage=PreMating or PostMating, without setting
    param;
    * func(pop, param) when stage=PreMating or PostMating, with param;
    * func(pop, off, dad, mom) when stage=DuringMating and
    passOffspringOnly=False, without setting param;
    * func(off) when stage=DuringMating and passOffspringOnly=True,
    and without setting param;
    * func(pop, off, dad, mom, param) when stage=DuringMating and
    passOffspringOnly=False, with param;
    * func(off, param) when stage=DuringMating and
    passOffspringOnly=True, with param. For Pre- and PostMating
    usages, a  population and an optional parameter is passed to the
    given function. For DuringMating usages,  population, offspring,
    its parents and an optional parameter are passed to the given
    function. Arbitrary operations can be applied to the  population
    and offspring (if stage=DuringMating).

Example:

>>> def dynaMutator(pop, param):
...   ''' this mutator mutate common loci with low mutation rate
...   and rare loci with high mutation rate, as an attempt to
...   bring allele frequency of these loci at an equal level.'''
...   # unpack parameter
...   (cutoff, mu1, mu2) = param;
...   Stat(pop, alleleFreq=range( pop.totNumLoci() ) )
...   for i in range( pop.totNumLoci() ):
...     # 1-freq of wild type = total disease allele frequency
...     if 1-pop.dvars().alleleFreq[i][1] < cutoff:
...       KamMutate(pop, maxAllele=2, rate=mu1, loci=[i])
...     else:
...       KamMutate(pop, maxAllele=2, rate=mu2, loci=[i])
...   return True
... 
>>> pop = population(size=10000, ploidy=2, loci=[2, 3])
>>> 
>>> simu = simulator(pop, randomMating())
>>> 
>>> simu.evolve(
...   preOps = [ 
...     initByFreq( [.6, .4], loci=[0,2,4]),
...     initByFreq( [.8, .2], loci=[1,3]) ],
...   ops = [ 
...     pyOperator( func=dynaMutator, param=(.5, .1, 0) ),
...     stat(alleleFreq=range(5)),
...     pyEval(r'\"%f\\\\t%f\\\\n\"%(alleleFreq[0][1],alleleFreq[1][1])', step=10)
...     ],
...   end = 30
... )        
0.387100	0.199450
0.392450	0.198800
0.409650	0.208250
0.406100	0.206300
True
>>>


"; 

%feature("docstring") simuPOP::pyOperator::pyOperator "

Description:

    Python operator, using a function that accepts a  population
    object.

Usage:

    pyOperator(func, param=None, stage=PostMating,
      formOffGenotype=False, passOffspringOnly=False, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    func:           a Python function. Its form is determined by other
                    parameters.
    param:          any Python object that will be passed to func
                    after pop parameter. Multiple parameters can be
                    passed as a tuple.
    formOffGenotype:This option tells the  mating scheme this operator
                    will set the genotype of offspring (valid only for
                    stage=DuringMating). By default
                    (formOffGenotype=False), a  mating scheme will set
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
    function.
    * This operator can be applied Pre-, During- or Post- mating and
    is applied PostMating by default. For example, if you would like
    to examine the fitness values set by a  selector, a PreMating
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

    apply the  pyOperator operator to one  population

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::pyOperator::applyDuringMating(population &pop, population::IndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::pyOperator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  pyOperator operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyPenetrance "

Function form:

    PyPenetrance

Function form:

    PyPenetrance

Description:

    assign  penetrance values by calling a user provided function

Details:

    For each  individual, users provide a function to calculate
    penetrance. This method is very flexible but will be slower than
    previous operators since a function will be called for each
    individual.

Details:

    For each  individual, users provide a function to calculate
    penetrance. This method is very flexible but will be slower than
    previous operators since a function will be called for each
    individual.

"; 

%feature("docstring") simuPOP::pyPenetrance::pyPenetrance "

Description:

    provide locus and  penetrance for 11, 12, 13 (in the form of
    dictionary)

Usage:

    pyPenetrance(loci, func, ancestralGen=-1, stage=DuringMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[])

Arguments:

    loci:           disease susceptibility loci. The genotypes at
                    these loci will be passed to the provided Python
                    function in the form of loc1_1, loc1_2, loc2_1,
                    loc2_2, ... if the individuals are diploid.
    func:           a user-defined Python function that accepts an
                    array of genotypes at susceptibility loci and
                    return a  penetrance value. The returned value
                    should be between 0 and 1.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)???

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

    deep copy of a Python  penetrance operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyPenetrance::penet "

Description:

    currently assuming diploid???

Usage:

    x.penet(ind)

"; 

%feature("docstring") simuPOP::pyPenetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the Python  penetrance operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyQuanTrait "

Function form:

    PyQuanTrait

Function form:

    PyQuanTrait

Description:

    quantitative trait using a user provided function

Details:

    For each  individual, a user provided function is used to
    calculate quantitative trait.

Details:

    For each  individual, a user provided function is used to
    calculate quantitative trait.

"; 

%feature("docstring") simuPOP::pyQuanTrait::pyQuanTrait "

Description:

    create a Python quantitative trait operator

Usage:

    pyQuanTrait(loci, func, ancestralGen=-1, stage=PostMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[\"qtrait\"])

Details:

    Please refer to  quanTrait for other parameter descriptions.

Details:

    Please refer to  quanTrait for other parameter descriptions.

Arguments:

    loci:           susceptibility loci. The genotypes at these loci
                    will be passed to func.
    func:           a Python function that accepts genotypes at
                    susceptibility loci and returns the quantitative
                    trait value.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::pyQuanTrait::~pyQuanTrait "

Description:

    simuPOP::pyQuanTrait::~pyQuanTrait

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

%feature("docstring") simuPOP::pySample "

Function form:

    PySample

Function form:

    PySample

Description:

    Python sampler.

Details:

    A Python sampler that generate a  sample with given individuals.

Details:

    A Python sampler that generate a  sample with given individuals.

"; 

%feature("docstring") simuPOP::pySample::pySample "

Description:

    create a Python sampler

Usage:

    pySample(keep, keepAncestralPops=-1, name=\"sample\", nameExpr=\"\",
      times=1, saveAs=\"\", saveAsExpr=\"\", format=\"auto\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Details:

    Please refer to class  sample for other parameter descriptions.

Details:

    Please refer to class  sample for other parameter descriptions.

Arguments:

    keep:           subpopulation IDs of all individuals
    keepAncestralPop:the number of ancestral populations that will be
                    kept. If -1, keep all ancestral populations
                    (default). If 0, no ancestral  population will be
                    kept.

"; 

%feature("docstring") simuPOP::pySample::~pySample "

Description:

    destructor

Usage:

    x.~pySample()

"; 

%ignore simuPOP::pySample::pySample(const pySample &rhs);

%feature("docstring") simuPOP::pySample::clone "

Description:

    deep copy of a Python sampler

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pySample::drawsample "

Description:

    draw a Python  sample

Usage:

    x.drawsample(pop)

"; 

%feature("docstring") simuPOP::pySample::__repr__ "

Description:

    used by Python print function to print out the general information
    of the Python sampler

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pySelector "

Function form:

    PySelect

Function form:

    PySelect

Description:

    selection using user provided function

Details:

    pySelector assigns fitness values by calling a user provided
    function. It accepts a list of susceptibility loci and a Python
    function. For each  individual, this operator will pass the
    genotypes at these loci (in the order of 0-0,0-1,1-0,1-1 etc.
    where X-Y represents locus X - ploidy Y, in the case of diploid
    population), generation number,??? and expect a returned fitness
    value. This, at least in theory, can accommodate all selection
    scenarios.

Details:

    pySelector assigns fitness values by calling a user provided
    function. It accepts a list of susceptibility loci and a Python
    function. For each  individual, this operator will pass the
    genotypes at these loci (in the order of 0-0,0-1,1-0,1-1 etc.
    where X-Y represents locus X - ploidy Y, in the case of diploid
    population), generation number,??? and expect a returned fitness
    value. This, at least in theory, can accommodate all selection
    scenarios.

"; 

%feature("docstring") simuPOP::pySelector::pySelector "

Description:

    create a Python hybrid  selector

Usage:

    pySelector(loci, func, subPops=[], stage=PreMating, begin=0,
      end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[\"fitness\"])

Arguments:

    loci:           susceptibility loci. The genotype at these loci
                    will be passed to func.
    func:           a Python function that accepts genotypes at
                    susceptibility loci generation number, and return
                    fitness value.???
    output:         and other parameters please refer to
                    help(baseOperator.__init__)???

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

%feature("docstring") simuPOP::pySubset "

Function form:

    PySubset

Function form:

    PySubset

Description:

    shrink  population

Details:

    This operator shrinks a  population according to a given array or
    the subPopID() value of each indvidual. Subpopulations are kept
    intact.

Details:

    This operator shrinks a  population according to a given array or
    the subPopID() value of each indvidual. Subpopulations are kept
    intact.

"; 

%feature("docstring") simuPOP::pySubset::pySubset "

Description:

    create a  pySubset operator

Usage:

    pySubset(keep=[], stage=PostMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    keep:           an array of subpopulation IDs for each
                    individual.

"; 

%feature("docstring") simuPOP::pySubset::~pySubset "

Description:

    destructor

Usage:

    x.~pySubset()

"; 

%feature("docstring") simuPOP::pySubset::clone "

Description:

    deep copy of a  pySubset operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pySubset::apply "

Description:

    apply the  pySubset operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pySubset::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  pySubset operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyTagger "

Details:

    This  tagger takes some information fields from both parents, pass
    to a Python function and set the  individual field with the
    returned value.
    This operator can be used to trace the inheritance of trait
    values.

Details:

    This  tagger takes some information fields from both parents, pass
    to a Python function and set the  individual field with the
    returned value.
    This operator can be used to trace the inheritance of trait
    values.

"; 

%feature("docstring") simuPOP::pyTagger::pyTagger "

Description:

    simuPOP::pyTagger::pyTagger

Usage:

    pyTagger(func=None, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    infoFields:     information fields. The user should gurantee the
                    existence of these fields.
    func:           a Pyton function that returns a list to assign the
                    information fields. e.g., if fields=['A', 'B'],
                    the function will pass values of fields 'A' and
                    'B' of father, followed by mother if there is one,
                    to this function. The returned value is assigned
                    to fields 'A' and 'B' of the offspring. The
                    returned value has to be a list even if only one
                    field is given.

"; 

%feature("docstring") simuPOP::pyTagger::~pyTagger "

Description:

    simuPOP::pyTagger::~pyTagger

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

    basic class of quantitative trait

Details:

    Quantitative trait is the measure of certain phenotype for given
    genotype. Quantitative trait is similar to  penetrance in that the
    consequence of  penetrance is binary: affected or unaffected;
    while it is continuous for quantitative trait.
    In  simuPOP, different operators/functions were implemented to
    calculate quantitative traits for each  individual and store the
    values in the information fields specified by user (default to
    qtrait). The quantitative trait operators also accept the
    ancestralGen parameter to control the number of generations for
    which the qtrait information field will be set.

Details:

    Quantitative trait is the measure of certain phenotype for given
    genotype. Quantitative trait is similar to  penetrance in that the
    consequence of  penetrance is binary: affected or unaffected;
    while it is continuous for quantitative trait.
    In  simuPOP, different operators/functions were implemented to
    calculate quantitative traits for each  individual and store the
    values in the information fields specified by user (default to
    qtrait). The quantitative trait operators also accept the
    ancestralGen parameter to control the number of generations for
    which the qtrait information field will be set.

"; 

%feature("docstring") simuPOP::quanTrait::quanTrait "

Description:

    create a quantitative trait operator, default to be always active

Usage:

    quanTrait(ancestralGen=-1, stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[\"qtrait\"])

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

    set qtrait to all  individual

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

Description:

    a  mating scheme of basic sexually random  mating

Details:

    In this scheme, sex information is considered for each
    individual, and ploidy is always 2. Within each subpopulation,
    males and females are randomly chosen. Then randomly get one copy
    of chromosomes from father and mother. When only one sex exists in
    a subpopulation, a parameter (contWhenUniSex) can be set to
    determine the behavior. Default to continuing without warning.

Details:

    In this scheme, sex information is considered for each
    individual, and ploidy is always 2. Within each subpopulation,
    males and females are randomly chosen. Then randomly get one copy
    of chromosomes from father and mother. When only one sex exists in
    a subpopulation, a parameter (contWhenUniSex) can be set to
    determine the behavior. Default to continuing without warning.

"; 

%feature("docstring") simuPOP::randomMating::randomMating "

Description:

    create a random  mating scheme

Usage:

    randomMating(numOffspring=1., numOffspringFunc=None,
      maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
      newSubPopSizeFunc=None, newSubPopSizeExpr=\"\",
      contWhenUniSex=True)

Arguments:

    numOffspring:   number of offspring or p in some modes
    numOffspringFunc:a python function that determines the number of
                    offspring or p
    maxNumOffspring:used when numOffspring is generated from a
                    binomial distribution
    mode:           can be one of MATE_NumOffspring,
                    MATE_NumOffspringEachFamily,
                    MATE_GeometricDistribution,
                    MATE_PoissonDistribution,
                    MATE_BinomialDistribution
    newSubPopSize:  an array of subpopulation sizes, should have the
                    same number of subpopulations as the current
                    population
    newSubPopSizeExpr:an expression that will be evaluated as an array
                    of subpopulation sizes
    newSubPopSizeFunc:an function that have parameter gen and oldSize
                    (current subpopulation size)
    contWhenUniSex: continue when there is only one sex in the
                    population, default to true
                    Please refer to class  mating for descriptions of
                    other parameters.

"; 

%feature("docstring") simuPOP::randomMating::~randomMating "

Description:

    destructor

Usage:

    x.~randomMating()

"; 

%feature("docstring") simuPOP::randomMating::clone "

Description:

    deep copy of a random  mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::randomMating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::randomMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the random  mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::randomMating::submitScratch(population &pop, population &scratch);

%ignore simuPOP::randomMating::mate(population &pop, population &scratch, vector< baseOperator * > &ops, bool submit);

%feature("docstring") simuPOP::randomSample "

Function form:

    RandomSample

Function form:

    RandomSample

Description:

    randomly draw a  sample from a  population

Details:

    This operator will randomly choose size individuals (or  size[i]
    individuals from subpopulation i) and return a new  population.
    The function form of this operator returns the samples directly.
    The operator keeps samples in an array name in the local
    namespace. You may access them through dvars() or vars()
    functions.
    The original subpopulation structure/boundary is kept in the
    samples.

Details:

    This operator will randomly choose size individuals (or  size[i]
    individuals from subpopulation i) and return a new  population.
    The function form of this operator returns the samples directly.
    The operator keeps samples in an array name in the local
    namespace. You may access them through dvars() or vars()
    functions.
    The original subpopulation structure/boundary is kept in the
    samples.

"; 

%feature("docstring") simuPOP::randomSample::randomSample "

Description:

    draw a random  sample, regardless of the affected status

Usage:

    randomSample(size=[], name=\"sample\", nameExpr=\"\", times=1,
      saveAs=\"\", saveAsExpr=\"\", format=\"auto\", stage=PostMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[])

Details:

    Please refer to class  sample for other parameter descriptions.

Details:

    Please refer to class  sample for other parameter descriptions.

Arguments:

    size:           size of the  sample. It can be either a number
                    which represents the overall  sample size,
                    regardless of the  population structure; or an
                    array which represents the number of samples drawn
                    from each subpopulation.

Note:

    Ancestral populations will not be copied to the samples.

"; 

%feature("docstring") simuPOP::randomSample::~randomSample "

Description:

    destructor

Usage:

    x.~randomSample()

"; 

%feature("docstring") simuPOP::randomSample::clone "

Description:

    deep copy of a  randomSample operator

Usage:

    x.clone()

"; 

%ignore simuPOP::randomSample::prepareSample(population &pop);

%ignore simuPOP::randomSample::drawsample(population &pop);

%feature("docstring") simuPOP::randomSample::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  randomSample operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::recombinator "

Description:

    recombination

Details:

    In  simuPOP, only one  recombinator is provided. Recombination
    events between loci a/b and b/c are independent, otherwise there
    will be some linkage between loci, users need to specify physical
    recombination rate between adjacent loci. In addition, for the
    recombinator
    * it only works for diploid (and for females in haplodiploid)
    populations.
    * the recombination rate must be comprised between 0.0 and 0.5. A
    recombination rate of 0.0 means that the loci are completely
    linked, and thus behave together as a single linked locus. A
    recombination rate of 0.5 is equivalent to free recombination. All
    other values between 0.0 and 0.5 will represent various linkage
    intensities between adjacent pairs of loci. The recombination rate
    is equivalent to 1-linkage and represents the probability that the
    allele at the next locus is randomly drawn.

Details:

    In  simuPOP, only one  recombinator is provided. Recombination
    events between loci a/b and b/c are independent, otherwise there
    will be some linkage between loci, users need to specify physical
    recombination rate between adjacent loci. In addition, for the
    recombinator
    * it only works for diploid (and for females in haplodiploid)
    populations.
    * the recombination rate must be comprised between 0.0 and 0.5. A
    recombination rate of 0.0 means that the loci are completely
    linked, and thus behave together as a single linked locus. A
    recombination rate of 0.5 is equivalent to free recombination. All
    other values between 0.0 and 0.5 will represent various linkage
    intensities between adjacent pairs of loci. The recombination rate
    is equivalent to 1-linkage and represents the probability that the
    allele at the next locus is randomly drawn.

"; 

%feature("docstring") simuPOP::recombinator::recombinator "

Description:

    recombine chromosomes from parents

Usage:

    recombinator(intensity=-1, rate=[], afterLoci=[],
      maleIntensity=-1, maleRate=[], maleAfterLoci=[], begin=0,
      end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    intensity:      intensity of recombination. The actually
                    recombination rate between two loci is determined
                    by intensity*locus distance between them.
    rate:           recombination rate regardless of locus distance
                    after all afterLoci. It can also be an array of
                    recombination rates. Should have the same length
                    as afterLoci or totNumOfLoci(). If totNumLoci, the
                    last item can be ignored.??? The recombination
                    rates are independent of locus distance.
    afterLoci:      an array of locus indices. Recombination will
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
                    locations. This is rarely used.???

Note:

    There is no recombination between sex chromosomes of male
    individuals if sexChrom()=True.??? This may change later if the
    exchanges of genes between pseudoautosomal regions of XY need to
    be modeled.

Example:

>>> simu = simulator(population(4, loci=[4,5,6], 
...     infoFields=['father_idx', 'mother_idx']),
...     randomMating())
>>> simu.step([
...   parentsTagger(),
...   ],
...   preOps = [initByFreq([.2,.2,.4,.2]), dumper(alleleOnly=True) ],
...   postOps = [ dumper(alleleOnly=True)]
... )
individual info: 
sub population 0:
   0: FU   3  0  1  3   1  2  0  2  2   3  0  2  3  0  1 |   2  0  0  3   2  2  2  3  2   2  3  2  3  2  0 
   1: MU   2  3  2  0   2  2  1  2  3   0  2  2  2  1  1 |   1  3  2  0   2  3  2  3  3   1  1  2  1  3  0 
   2: MU   2  1  3  0   2  0  1  3  3   2  1  3  2  0  2 |   2  2  1  2   3  2  3  3  3   3  3  2  2  2  2 
   3: MU   2  2  0  0   2  2  1  3  0   3  1  0  3  2  0 |   2  3  1  3   2  3  2  2  2   0  1  2  0  1  2 
End of individual info.


No ancenstral population recorded.
individual info: 
sub population 0:
   0: FU   1  3  2  0   2  2  1  2  3   1  1  2  1  3  0 |   3  0  1  3   1  2  0  2  2   3  0  2  3  0  1 
   1: FU   2  3  1  3   2  2  1  3  0   0  1  2  0  1  2 |   2  0  0  3   2  2  2  3  2   2  3  2  3  2  0 
   2: MU   2  2  1  2   2  0  1  3  3   3  3  2  2  2  2 |   3  0  1  3   1  2  0  2  2   2  3  2  3  2  0 
   3: MU   2  2  0  0   2  3  2  2  2   3  1  0  3  2  0 |   3  0  1  3   2  2  2  3  2   2  3  2  3  2  0 
End of individual info.


No ancenstral population recorded.
True
>>> simu.step([
...   parentsTagger(),
...   recombinator(rate=[1,1,1], afterLoci=[2,6,10])
...   ],
...   postOps = [ dumper(alleleOnly=True)]
... )
individual info: 
sub population 0:
   0: MU   1  3  2  3   2  2  1  2  2   1  1  2  3  0  1 |   3  0  1  0   2  3  2  3  2   2  3  0  3  2  0 
   1: FU   2  3  1  3   2  2  1  3  2   0  1  2  3  2  0 |   2  2  0  3   2  2  2  2  2   2  3  0  3  2  0 
   2: MU   3  0  1  0   2  2  1  2  2   1  1  2  3  0  1 |   2  2  1  3   1  2  0  3  3   3  3  2  3  2  0 
   3: FU   2  0  0  3   2  2  1  3  2   0  1  2  3  2  0 |   3  0  1  2   1  2  0  3  3   2  3  2  2  2  2 
End of individual info.


No ancenstral population recorded.
True
>>>


"; 

%feature("docstring") simuPOP::recombinator::~recombinator "

Description:

    simuPOP::recombinator::~recombinator

Usage:

    x.~recombinator()

"; 

%feature("docstring") simuPOP::recombinator::clone "

Description:

    deep copy of a  recombinator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::recombinator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  recombinator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::recombinator::recCount "

Description:

    return recombination count

Usage:

    x.recCount(locus)

"; 

%feature("docstring") simuPOP::recombinator::recCounts "

Description:

    return recombination counts

Usage:

    x.recCounts()

"; 

%feature("docstring") simuPOP::recombinator::applyDuringMating "

Description:

    apply the  recombinator during  mating???

Usage:

    x.applyDuringMating(pop, offspring, dad=None, mom=None)

"; 

%feature("docstring") simuPOP::RNG "

Description:

    simuPOP::RNG

"; 

%feature("docstring") simuPOP::RNG::RNG "

Description:

    simuPOP::RNG::RNG

Usage:

    RNG(rng=None, seed=0)

"; 

%feature("docstring") simuPOP::RNG::~RNG "

Description:

    simuPOP::RNG::~RNG

Usage:

    x.~RNG()

"; 

%ignore simuPOP::RNG::setRNG(const char *rng=NULL, unsigned long seed=0);

%feature("docstring") simuPOP::RNG::name "

Description:

    return  RNG name

Usage:

    x.name()

"; 

%feature("docstring") simuPOP::RNG::seed "

Description:

    simuPOP::RNG::seed

Usage:

    x.seed()

"; 

%feature("docstring") simuPOP::RNG::maxSeed "

Description:

    simuPOP::RNG::maxSeed

Usage:

    x.maxSeed()

"; 

%feature("docstring") simuPOP::RNG::setSeed "

Description:

    simuPOP::RNG::setSeed

Usage:

    x.setSeed(seed)

"; 

%feature("docstring") simuPOP::RNG::generateRandomSeed "

Description:

    simuPOP::RNG::generateRandomSeed

Usage:

    x.generateRandomSeed()

"; 

%feature("docstring") simuPOP::RNG::max "

Description:

    simuPOP::RNG::max

Usage:

    x.max()

"; 

%feature("docstring") simuPOP::RNG::__repr__ "

Description:

    simuPOP::RNG::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::RNG::randGet "

Description:

    return min, 1, 2, ... max

Usage:

    x.randGet()

"; 

%feature("docstring") simuPOP::RNG::randInt "

Description:

    return [0, 1, 2, ... n-1]

Usage:

    x.randInt(n)

"; 

%feature("docstring") simuPOP::RNG::randIntArray "

Description:

    return [0, 1, 2, ... n-1]

Usage:

    x.randIntArray(n, size, vec)

"; 

%feature("docstring") simuPOP::RNG::randGeometric "

Description:

    return geometric with parameter p (k>=1)

Usage:

    x.randGeometric(p)

"; 

%feature("docstring") simuPOP::RNG::randUniform01 "

Description:

    return double [0,1)

Usage:

    x.randUniform01()

"; 

%feature("docstring") simuPOP::RNG::randNormal "

Description:

    return double -inf, inf, v is standard deviation

Usage:

    x.randNormal(m, v)

"; 

%feature("docstring") simuPOP::RNG::randUniform01Array "

Description:

    return double [0,1)

Usage:

    x.randUniform01Array(size, vec)

"; 

%feature("docstring") simuPOP::RNG::randBinomial "

Description:

    binomial distribution

Usage:

    x.randBinomial(n, p)

"; 

%feature("docstring") simuPOP::RNG::randMultinomial "

Description:

    simuPOP::RNG::randMultinomial

Usage:

    x.randMultinomial(N, p, n)

"; 

%feature("docstring") simuPOP::RNG::randMultinomialVal "

Description:

    simuPOP::RNG::randMultinomialVal

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

%feature("docstring") simuPOP::sample "

Description:

    basic class of other  sample operator

Details:

    Ascertainment/sampling refers to ways to select individuals from a
    population. In  simuPOP, ascerntainment operators form separate
    populations in a population's namespace. All the ascertainment
    operators work like this except for  pySubset which shrink the
    population itself.
    Individuals in sampled populations may or may not keep their
    original order but their indices in the whole  population are
    stored in a information field oldindex. That is to say, you can
    use ind.info('oldindex') to check the original position of an
    individual.
    Two forms of  sample size specification are supported: with or
    without subpopulation structure. For example, the size parameter
    of  randomSample can be a number or an array (which has the length
    of the number of subpopulations). If a number is given, a  sample
    will be drawn from the whole  population, regardless of the
    population structure. If an array is given, individuals will be
    drawn from each subpopulation sp according to size[sp].
    An important special case of  sample size specification occurs
    when size=[] (default). In this case, usually all qualified
    individuals will be returned.
    The function forms of these operators are a little different from
    others. They do return a value: an array of samples.

Details:

    Ascertainment/sampling refers to ways to select individuals from a
    population. In  simuPOP, ascerntainment operators form separate
    populations in a population's namespace. All the ascertainment
    operators work like this except for  pySubset which shrink the
    population itself.
    Individuals in sampled populations may or may not keep their
    original order but their indices in the whole  population are
    stored in a information field oldindex. That is to say, you can
    use ind.info('oldindex') to check the original position of an
    individual.
    Two forms of  sample size specification are supported: with or
    without subpopulation structure. For example, the size parameter
    of  randomSample can be a number or an array (which has the length
    of the number of subpopulations). If a number is given, a  sample
    will be drawn from the whole  population, regardless of the
    population structure. If an array is given, individuals will be
    drawn from each subpopulation sp according to size[sp].
    An important special case of  sample size specification occurs
    when size=[] (default). In this case, usually all qualified
    individuals will be returned.
    The function forms of these operators are a little different from
    others. They do return a value: an array of samples.

"; 

%feature("docstring") simuPOP::sample::sample "

Description:

    draw a  sample

Usage:

    sample(name=\"sample\", nameExpr=\"\", times=1, saveAs=\"\",
      saveAsExpr=\"\", format=\"auto\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Details:

    Please refer to baseOperator::__init__ for other parameters.

Details:

    Please refer to baseOperator::__init__ for other parameters.

Arguments:

    name:           name of the  sample in local namespace. This
                    variable is an array of populations of size times.
                    Default to  sample. If name='' is set, samples
                    will not be saved in local namespace.
    nameExpr:       expression version of parameter name. If both name
                    and nameExpr are empty, do not store pop. This
                    expression will be evaluated dynamically in
                    population's local namespace.
    times:          how many times to  sample from the  population.
                    This is usually 1, but we may want to take several
                    random samples.
    saveAs:         filename to save the samples
    saveAsExpr:     expression version of parameter saveAs. It will be
                    evaluated dynamically in population's local
                    namespace.
    format:         format to save the samples

"; 

%feature("docstring") simuPOP::sample::~sample "

Description:

    destructor

Usage:

    x.~sample()

"; 

%feature("docstring") simuPOP::sample::clone "

Description:

    deep copy of a  sample operator

Usage:

    x.clone()

"; 

%ignore simuPOP::sample::prepareSample(population &);

%ignore simuPOP::sample::drawsample(population &pop);

%feature("docstring") simuPOP::sample::samples "

Description:

    return the samples

Usage:

    x.samples(pop)

"; 

%feature("docstring") simuPOP::sample::apply "

Description:

    apply the  sample operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::sample::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  sample operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::sample::saveIndIndex "

Description:

    save the index of each  individual to a field (usually oldindex)

Usage:

    x.saveIndIndex(pop, indexField=\"oldindex\")

"; 

%feature("docstring") simuPOP::sample::resetParentalIndex "

Description:

    reset father_idx and mother_idx

Usage:

    x.resetParentalIndex(pop, fatherField=\"father_idx\",
      motherField=\"mother_idx\", indexField=\"oldindex\")

"; 

%feature("docstring") simuPOP::sample::findOffspringAndSpouse "

Description:

    find offspring and spouse

Usage:

    x.findOffspringAndSpouse(pop, ancestralDepth, maxOffspring,
      fatherField, motherField, spouseField, offspringField)

"; 

%feature("docstring") simuPOP::sample::resetSubPopID "

Description:

    set all subpopulation IDs to -1 (remove)

Usage:

    x.resetSubPopID(pop)

"; 

%feature("docstring") simuPOP::savePopulation "

Description:

    save  population to a file

"; 

%feature("docstring") simuPOP::savePopulation::savePopulation "

Description:

    simuPOP::savePopulation::savePopulation

Usage:

    savePopulation(output=\"\", outputExpr=\"\", format=\"bin\",
      compress=True, stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::savePopulation::~savePopulation "

Description:

    simuPOP::savePopulation::~savePopulation

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

    apply to one  population. It does not check if the operator is
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

    genetic selection

Details:

    Genetic selection is tricky to simulate since there are many
    different fitness values and many different ways to apply
    selection.  simuPOP employs an 'ability-to-mate' approach. Namely,
    the probability that an  individual will be chosen for  mating is
    proportional to its fitness value. More specifically,
    * PreMating selectors assign fitness values to each  individual.
    * During sexless  mating (e.g.  binomialSelection???), individuals
    are chosen at probabilities that are proportional to their fitness
    values. If there are  $ N $ individuals with fitness values  $
    f_{i},i=1,...,N $,  individual $ i $ will have probability  $
    \\frac{f_{i}}{\\sum_{j}f_{j}} $ to be chosen and passed to the next
    generation.
    * During  randomMating, males and females are separated. They are
    chosen from their respective groups in the same manner and mate.
    It is not very clear that our method agrees with the traditional
    'average number of offspring' definition of fitness. (Note that
    this concept is very difficult to simulate since we do not know
    who will determine the number of offspring if two parents are
    involved.)All of the selection operators, when applied, will set a
    variable fitness and an indicator so that 'selector-aware'  mating
    scheme can select individuals according to these values. Hence,
    two consequences are stated below:
    *  selector alone can not do selection! Only  mating schemes can
    actually select individuals.
    *  selector has to be PreMating operator. This is not a problem
    when you use the operator form of the selectors since their
    default stage is PreMating. However, if you use the function form
    of these selectors in a  pyOperator, make sure to set the stage of
    pyOperator to PreMating.

Details:

    Genetic selection is tricky to simulate since there are many
    different fitness values and many different ways to apply
    selection.  simuPOP employs an 'ability-to-mate' approach. Namely,
    the probability that an  individual will be chosen for  mating is
    proportional to its fitness value. More specifically,
    * PreMating selectors assign fitness values to each  individual.
    * During sexless  mating (e.g.  binomialSelection???), individuals
    are chosen at probabilities that are proportional to their fitness
    values. If there are  $ N $ individuals with fitness values  $
    f_{i},i=1,...,N $,  individual $ i $ will have probability  $
    \\frac{f_{i}}{\\sum_{j}f_{j}} $ to be chosen and passed to the next
    generation.
    * During  randomMating, males and females are separated. They are
    chosen from their respective groups in the same manner and mate.
    It is not very clear that our method agrees with the traditional
    'average number of offspring' definition of fitness. (Note that
    this concept is very difficult to simulate since we do not know
    who will determine the number of offspring if two parents are
    involved.)All of the selection operators, when applied, will set a
    variable fitness and an indicator so that 'selector-aware'  mating
    scheme can select individuals according to these values. Hence,
    two consequences are stated below:
    *  selector alone can not do selection! Only  mating schemes can
    actually select individuals.
    *  selector has to be PreMating operator. This is not a problem
    when you use the operator form of the selectors since their
    default stage is PreMating. However, if you use the function form
    of these selectors in a  pyOperator, make sure to set the stage of
    pyOperator to PreMating.

"; 

%feature("docstring") simuPOP::selector::selector "

Description:

    create a  selector

Usage:

    selector(subPops=[], stage=PreMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[\"fitness\"])

Arguments:

    subPop:         a shortcut to subPops=[subPop]
    subPops:        subpopulations that the  selector will apply to.
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

    deep copy of a  selector

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::selector::indFitness "

Description:

    calculate/return the fitness value???

Usage:

    x.indFitness(, gen)

"; 

%feature("docstring") simuPOP::selector::apply "

Description:

    set fitness to all individuals???

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::selector::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  selector

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::setAncestralDepth "

Description:

    set ancestral depth

Details:

    This operator set the number of ancestral generations to keep in a
    population. It is usually called like setAncestral(at=[-2]) to
    start recording ancestral generations to a  population at the end
    of the evolution. This is useful when constructing pedigree trees
    from a  population.

Details:

    This operator set the number of ancestral generations to keep in a
    population. It is usually called like setAncestral(at=[-2]) to
    start recording ancestral generations to a  population at the end
    of the evolution. This is useful when constructing pedigree trees
    from a  population.

"; 

%feature("docstring") simuPOP::setAncestralDepth::setAncestralDepth "

Description:

    create a  setAncestralDepth operator

Usage:

    setAncestralDepth(depth, output=\">\", outputExpr=\"\",
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

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

    apply the  setAncestralDepth operator to one  population

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

%ignore simuPOP::SharedVariables::getVar(const string &name, bool nameError=true);

%feature("docstring") simuPOP::SharedVariables::hasVar "

Description:

    simuPOP::SharedVariables::hasVar

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

Description:

    simuPOP::SharedVariables::dict

Usage:

    x.dict()

"; 

%ignore simuPOP::SharedVariables::asString() const;

%feature("docstring") simuPOP::SharedVariables::fromString "

Description:

    simuPOP::SharedVariables::fromString

Usage:

    x.fromString(vars)

"; 

%feature("docstring") simuPOP::simulator "

Description:

    simulator manages several replicates of a  population, evolve them
    using given  mating scheme and operators

Details:

    Simulators combine three important components of  simuPOP:
    population,  mating scheme and operators together. A  simulator is
    created with an instance of  population, a replicate number rep
    and a  mating scheme. It makes rep number of replicates of this
    population and control the evolution process of them.
    The most important function of a  simulator is  evolve(). It
    accepts an array of operators as its parameters, among which,
    preOps and postOps will be applied to the populations at the
    beginning and the end of evolution, respectively, whereas ops will
    be applied at every generation.
    Simulators separate operators into pre-, during-, and post-mating
    operators. During evolution, a  simulator first apply all pre-
    mating operators and then call the mate() function of the given
    mating scheme, which will call during-mating operators during the
    birth of each offspring. After  mating is completed, post-mating
    operators are applied to the offspring in the order at which they
    appear in the operator list.
    Operators can be applied to a specific replicate, a group of
    replicates, or specific generations, determined by the rep, grp,
    begin, end, step, and at parameters.
    Simulators can evolve a given number of generations (the end
    parameter of evolve), or evolve indefinitely until a certain type
    of operators called terminators terminates it. In this case, one
    or more terminators will check the status of evolution and
    determine if the simulation should be stopped. An obvious example
    of such a  terminator is a fixation-checker.
    Finally, a  simulator can be saved to a file in the format of
    'txt', 'bin', or 'xml'. So we can stop a simulation and resume it
    at another time or on another machine. It is also a good idea to
    save a snapshot of a simulation every several hundred generations.

Details:

    Simulators combine three important components of  simuPOP:
    population,  mating scheme and operators together. A  simulator is
    created with an instance of  population, a replicate number rep
    and a  mating scheme. It makes rep number of replicates of this
    population and control the evolution process of them.
    The most important function of a  simulator is  evolve(). It
    accepts an array of operators as its parameters, among which,
    preOps and postOps will be applied to the populations at the
    beginning and the end of evolution, respectively, whereas ops will
    be applied at every generation.
    Simulators separate operators into pre-, during-, and post-mating
    operators. During evolution, a  simulator first apply all pre-
    mating operators and then call the mate() function of the given
    mating scheme, which will call during-mating operators during the
    birth of each offspring. After  mating is completed, post-mating
    operators are applied to the offspring in the order at which they
    appear in the operator list.
    Operators can be applied to a specific replicate, a group of
    replicates, or specific generations, determined by the rep, grp,
    begin, end, step, and at parameters.
    Simulators can evolve a given number of generations (the end
    parameter of evolve), or evolve indefinitely until a certain type
    of operators called terminators terminates it. In this case, one
    or more terminators will check the status of evolution and
    determine if the simulation should be stopped. An obvious example
    of such a  terminator is a fixation-checker.
    Finally, a  simulator can be saved to a file in the format of
    'txt', 'bin', or 'xml'. So we can stop a simulation and resume it
    at another time or on another machine. It is also a good idea to
    save a snapshot of a simulation every several hundred generations.

"; 

%feature("docstring") simuPOP::simulator::simulator "

Description:

    create a  simulator

Usage:

    simulator(pop, matingScheme, stopIfOneRepStops=False,
      applyOpToStoppedReps=False, rep=1, grp=[])

Arguments:

    population:     a  population created by population() function.
                    This  population will be copied rep times to the
                    simulator. Its content will not be changed.
    matingScheme:   a  mating scheme
    rep:            number of replicates. Default to 1.
    grp:            group number for each replicate. Operators can be
                    applied to a group of replicates using its grp
                    parameter.
    applyOpToStoppedReps:If set, the  simulator will continue to apply
                    operators to all stopped replicates until all
                    replicates are marked 'stopped'.
    stopIfOneRepStops:If set, the  simulator will stop evolution if one
                    replicate stops.

"; 

%feature("docstring") simuPOP::simulator::~simulator "

Description:

    destroy a  simulator along with all its populations

Usage:

    x.~simulator()

Note:

    pop = simulator::population() returns temporary reference to an
    internal  population. After a  simulator evolves another genertion
    or after the  simulator is destroyed, this referenced  population
    should not be used.

"; 

%feature("docstring") simuPOP::simulator::clone "

Description:

    deep copy of a  simulator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::simulator::addInfoField "

Description:

    add an information field to all replicates

Usage:

    x.addInfoField(field, init=0)

Details:

    Add an information field to all replicate, and to the  simulator
    itself. This is important because all populations must have the
    same genotypic information as the  simulator. Adding an
    information field to one or more of the replicates will compromise
    the integrity of the  simulator.

Details:

    Add an information field to all replicate, and to the  simulator
    itself. This is important because all populations must have the
    same genotypic information as the  simulator. Adding an
    information field to one or more of the replicates will compromise
    the integrity of the  simulator.

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

    the rep replicate of this  simulator

Usage:

    x.pop(rep)

Details:

    This function is named  population in the Python interface.

Details:

    This function is named  population in the Python interface.

Arguments:

    rep:            the index number of replicate which will be
                    accessed

Note:

    The returned reference is temporary in the sense that the refered
    population will be invalid after another round of evolution.
    Therefore, the use of this function should be limited to
    immediateafterretrival. If you would like to get a persistent
    population, please use getPopulation(rep).

"; 

%feature("docstring") simuPOP::simulator::getPopulation "

Description:

    return a copy of  populationrep

Usage:

    x.getPopulation(rep, destructive=False)

Details:

    return a temporary reference of one of the populations.
    'Reference' means that the changes to the referred  population
    will reflect to the one in  simulator. 'Temporary' means that the
    referred  population might be invalid after evolution.

Details:

    return a temporary reference of one of the populations.
    'Reference' means that the changes to the referred  population
    will reflect to the one in  simulator. 'Temporary' means that the
    referred  population might be invalid after evolution.

Arguments:

    rep:            the index number of the replicate which will be
                    obtained
    destructive:    if true, destroy the copy of  population within
                    this  simulator. Default to false.
                    getPopulation(rep, true) is a more efficient way
                    to get hold of a  population when the  simulator
                    will no longer be used.

"; 

%feature("docstring") simuPOP::simulator::setMatingScheme "

Description:

    set  mating scheme

Usage:

    x.setMatingScheme(matingScheme)

"; 

%ignore simuPOP::simulator::setPopulation(population &pop, UINT rep);

%ignore simuPOP::simulator::curRep() const;

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

%ignore simuPOP::simulator::grp();

%feature("docstring") simuPOP::simulator::group "

Description:

    return group indices

Usage:

    x.group()

"; 

%feature("docstring") simuPOP::simulator::setGroup "

Description:

    set groups for replicates

Usage:

    x.setGroup(grp)

"; 

%feature("docstring") simuPOP::simulator::setGen "

Description:

    set the current generation. Usually used to reset a  simulator.

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

    evolve all replicates of the  population, subject to operators

Usage:

    x.evolve(ops, preOps=[], postOps=[], end=-1, dryrun=False)

Details:

    Evolve to the end generation unless an operator ( terminator)
    stops it earlier.
    ops will be applied in the order of:
    * all pre-mating opertors
    * during-mating operators called by the  mating scheme at the
    birth of each offspring
    * all post-mating operators If any pre- or post-mating operator
    fails to apply, that replicate will be stopped. The behavior of
    the  simulator will be determined by flags applyOpToStoppedReps
    and stopIfOneRepStopss. This is exactly how terminators work.

Details:

    Evolve to the end generation unless an operator ( terminator)
    stops it earlier.
    ops will be applied in the order of:
    * all pre-mating opertors
    * during-mating operators called by the  mating scheme at the
    birth of each offspring
    * all post-mating operators If any pre- or post-mating operator
    fails to apply, that replicate will be stopped. The behavior of
    the  simulator will be determined by flags applyOpToStoppedReps
    and stopIfOneRepStopss. This is exactly how terminators work.

Arguments:

    ops:            operators that will be applied at each generation,
                    if they are active at that generation. (Determined
                    by the begin, end, step and at parameters of the
                    operator.)
    preOps:         operators that will be applied before evolution.
                    evolve() function will not check if they are
                    active.
    postOps:        operators that will be applied after evolution
    end:            ending generation. Default to -1. In this case,
                    there is no ending generation and a  simulator
                    will only be ended by a  terminator. Otherwise, it
                    should be a number greater than current generation
                    number.
    dry:            run mode. Default to False.

Note:

    When end = -1, you can not specify negative generation parameters
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

    get  simulator namespace. If rep > 0 is given, return the
    namespace of replicate rep

Usage:

    x.vars(rep, subPop=-1)

"; 

%feature("docstring") simuPOP::simulator::saveSimulator "

Description:

    save  simulator in 'txt', 'bin' or 'xml' format

Usage:

    x.saveSimulator(filename, format=\"auto\", compress=True)

Details:

    The default format is 'txt' but the output is not supposed to be
    read. 'bin' has smaller size and should be used for large
    populations. 'xml' format is most verbose and should be used when
    you would like to convert  simuPOP populations to other formats.

Details:

    The default format is 'txt' but the output is not supposed to be
    read. 'bin' has smaller size and should be used for large
    populations. 'xml' format is most verbose and should be used when
    you would like to convert  simuPOP populations to other formats.

Arguments:

    filename:       filename to save the  simulator. Default to simu.
    format:         format to save. Default to auto. I.e., determine
                    the format by file extensions.
    compress:       whether or not compress the file in 'gzip' format

"; 

%ignore simuPOP::simulator::loadSimulator(string filename, string format="auto");

%feature("docstring") simuPOP::simulator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  simulator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::smmMutator "

Function form:

    SmmMutate

Function form:

    SmmMutate

Description:

    stepwise mutation model

Details:

    Stepwise Mutation Model (SMM) assumes that alleles are represented
    by integer values and that a mutation either increases or
    decreases the allele value by one. For variable number tandem
    repeats loci (VNTR), the allele value is generally taken as the
    number of tandem repeats in the DNA sequence.

Details:

    Stepwise Mutation Model (SMM) assumes that alleles are represented
    by integer values and that a mutation either increases or
    decreases the allele value by one. For variable number tandem
    repeats loci (VNTR), the allele value is generally taken as the
    number of tandem repeats in the DNA sequence.

"; 

%feature("docstring") simuPOP::smmMutator::smmMutator "

Description:

    create a SMM  mutator

Usage:

    smmMutator(rate=[], loci=[], maxAllele=0, incProb=0.5,
      output=\">\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Details:

    The stepwise mutation model (SMM) is developed for allozymes. It
    provides better description for these kinds of evolutionary
    processes. Please see  mutator for the description of other
    parameters.

Details:

    The stepwise mutation model (SMM) is developed for allozymes. It
    provides better description for these kinds of evolutionary
    processes. Please see  mutator for the description of other
    parameters.

Arguments:

    incProb:        probability to increase allele state. Default to
                    0.5.

Example:

>>> simu = simulator(population(size=3, loci=[3,5]), noMating())
>>> simu.step([
...     initByFreq( [.2,.3,.5]),
...     smmMutator(rate=1,  incProb=.8),
...     dumper(alleleOnly=True, stage=PrePostMating)])
individual info: 
sub population 0:
   0: FU   0  2  0   0  2  2  2  0 |   2  2  0   2  2  0  0  2 
   1: FU   1  0  1   2  0  2  0  2 |   2  1  0   2  2  2  0  2 
   2: FU   0  1  2   2  0  2  2  2 |   0  2  1   0  0  0  2  2 
End of individual info.


No ancenstral population recorded.
individual info: 
sub population 0:
   0: FU   1  3  1   1  3  3  3  1 |   3  3  1   3  1  1  1  3 
   1: FU   2  1  2   3  1  1  1  3 |   1  2  1   3  1  3  0  3 
   2: FU   0  0  3   3  0  3  3  1 |   1  3  2   1  1  1  3  3 
End of individual info.


No ancenstral population recorded.
True
>>>


"; 

%feature("docstring") simuPOP::smmMutator::~smmMutator "

Description:

    simuPOP::smmMutator::~smmMutator

Usage:

    x.~smmMutator()

"; 

%feature("docstring") simuPOP::smmMutator::mutate "

Description:

    mutate according to the SMM model ???

Usage:

    x.mutate(allele)

"; 

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

Function form:

    SplitSubPop

Description:

    split a subpopulation

Details:

    

Details:

    

"; 

%feature("docstring") simuPOP::splitSubPop::splitSubPop "

Description:

    split a subpopulation or the whole  population as subpopulation 0

Usage:

    splitSubPop(which=0, sizes=[], proportions=[], subPopID=[],
      randomize=True, stage=PreMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    which:          which subpopulation to split. If there is no
                    subpopulation structure, use 0 as the first (and
                    only) subpopulation.
    sizes:          new subpopulation sizes. The sizes should be added
                    up to the original subpopulation (subpopulation
                    which) size.
    proportions:    proportions of new subpopulations. Should be added
                    up to 1. Optionally, you can use one of subPopID
                    or proportions to split. ???
    subPopID:       new subpopulation IDs. Otherwise, the operator
                    will automatically set new subpopulation IDs to
                    new subpopulations. If given, should have the same
                    length as subPop or proportions.??? Since subpop
                    with negative id will be removed. You can remove
                    part of a subpop by setting a new negative ID.???

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

Function form:

    Spread

Description:

    initialize genotype by value and then copy to all individuals

Details:

    Spread(ind, subPop) spreads the genotype of ind to all individuals
    in an array of subpopulations. The default value of subPop is the
    subpopulation where ind resides.

Details:

    Spread(ind, subPop) spreads the genotype of ind to all individuals
    in an array of subpopulations. The default value of subPop is the
    subpopulation where ind resides.

"; 

%feature("docstring") simuPOP::spread::spread "

Description:

    copy genotypes of ind to all individuals in subPop

Usage:

    spread(ind, subPop=[], stage=PreMating, begin=0, end=1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::spread::~spread "

Description:

    simuPOP::spread::~spread

Usage:

    x.~spread()

"; 

%feature("docstring") simuPOP::spread::clone "

Description:

    deep copy of the operator  spread

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::spread::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator  spread

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::spread::apply "

Description:

    apply operator  spread???

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::stat "

Function form:

    Stat

Function form:

    Stat

Description:

    calculate statistics

Details:

    Operator  stat calculatse various basic statistics for the
    population and sets variables in the local namespace. Other
    operators/functions can refer to the results from the namespace
    after  stat is applied. Stat is the function form of the operator.
    Note that these statistics are dependent to each other. For
    example, heterotype and allele frequencies of related loci will be
    automatically calculated if linkage diseqilibrium is requested.

Details:

    Operator  stat calculatse various basic statistics for the
    population and sets variables in the local namespace. Other
    operators/functions can refer to the results from the namespace
    after  stat is applied. Stat is the function form of the operator.
    Note that these statistics are dependent to each other. For
    example, heterotype and allele frequencies of related loci will be
    automatically calculated if linkage diseqilibrium is requested.

"; 

%feature("docstring") simuPOP::stat::stat "

Description:

    create an  stat operator

Usage:

    stat(popSize=False, numOfMale=False, numOfMale_param={},
      numOfAffected=False, numOfAffected_param={}, numOfAlleles=[],
      numOfAlleles_param={}, alleleFreq=[], alleleFreq_param={},
      heteroFreq=[], expHetero=[], expHetero_param={}, homoFreq=[],
      genoFreq=[], haploFreq=[], LD=[], LD_param={}, association=[],
      association_param={}, Fst=[], Fst_param={}, relGroups=[],
      relLoci=[], rel_param={}, relBySubPop=False, relMethod=[],
      relMinScored=10, hasPhase=False, midValues=False, output=\"\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Details:

    If only one item is specified, the outer [] can be ignored. I.e.,
    LD=[loc1, loc2] is acceptable. This parameter will set the
    following variables. Please note that the difference between the
    data structures used for ld and LD. The names are potentially very
    confusing but I have no better idea.
    * ld['loc1-loc2']['allele1-allele2'], subPop[sp]['ld']['loc1-
    loc2']['allele1-allele2']
    * ld_prime['loc1-loc2']['allele1-allele2'],
    subPop[sp]['ld_prime']['loc1-loc2']['allele1-allele2']
    * r2['loc1-loc2']['allele1-allele2'], subPop[sp]['r2']['loc1-
    loc2']['allele1-allele2']
    * LD[loc1][loc2], subPop[sp]['LD'][loc1][loc2]
    * LD_prime[loc1][loc2], subPop[sp]['LD_prime'][loc1][loc2]
    * R2[loc1][loc2], subPop[sp]['R2'][loc1][loc2]

Details:

    If only one item is specified, the outer [] can be ignored. I.e.,
    LD=[loc1, loc2] is acceptable. This parameter will set the
    following variables. Please note that the difference between the
    data structures used for ld and LD. The names are potentially very
    confusing but I have no better idea.
    * ld['loc1-loc2']['allele1-allele2'], subPop[sp]['ld']['loc1-
    loc2']['allele1-allele2']
    * ld_prime['loc1-loc2']['allele1-allele2'],
    subPop[sp]['ld_prime']['loc1-loc2']['allele1-allele2']
    * r2['loc1-loc2']['allele1-allele2'], subPop[sp]['r2']['loc1-
    loc2']['allele1-allele2']
    * LD[loc1][loc2], subPop[sp]['LD'][loc1][loc2]
    * LD_prime[loc1][loc2], subPop[sp]['LD_prime'][loc1][loc2]
    * R2[loc1][loc2], subPop[sp]['R2'][loc1][loc2]

Arguments:

    popSize:        whether or not calculate  population sizes. This
                    parameter will set the following variables:
                    * numSubPop the number of subpopulations
                    * subPopSize an array of subpopulation sizes. Not
                    available for subpopulations.
                    * popSize, subPop[sp]['popSize']
                    population/subpopulation size.
    numOfMale:      whether or not count the numbers/proportions of
                    males and females. This parameter can set the
                    following variables by user's specification:
                    * numOfMale, subPop[sp]['numOfMale'] the number of
                    males in the population/subpopulation
                    * numOfFemale, subPop[sp]['numOfFemale'] the
                    number of females in the population/subpopulation.
                    * propOfMale, subPop[sp]['propOfMale'] the
                    proportion of males in the
                    population/subpopulation
                    * propOfFemale, subPop[sp]['propOfFemale'] the
                    proportion of females in the
                    population/subpopulation
    numOfMale_param:a dictionary of parameters of numOfMale
                    statistics. Can be one or more items choosen from
                    the following options: numOfMale, propOfMale,
                    numOfFemale, and propOfFemale.
    numOfAffected:  whether or not count the numbers/proportions of
                    affected and unaffected individuals. This
                    parameter can set the following variables by
                    user's specification:
                    * numOfAffected, subPop[sp]['numOfAffected'] the
                    number of affected individuals in the
                    population/subpopulation
                    * numOfUnaffected, subPop[sp]['numOfUnAffected']
                    the number of unaffected individuals in the
                    population/subpopulation
                    * propOfAffected, subPop[sp]['propOfAffected'] the
                    proportion of affected individuals in the
                    population/subpopulation
                    * propOfUnaffected, subPop[sp]['propOfUnAffected']
                    the proportion of unaffected individuals in the
                    population/subpopulation
    numOfAffected_param:a dictionary of parameters of numOfAffected
                    statistics. Can be one or more items choosen from
                    the following options: numOfAffected,
                    propOfAffected, numOfUnaffected, propOfUnaffected.
    numOfAlleles:   an array of loci at which the numbers of distinct
                    alleles will be counted (numOfAlleles=[loc1, loc2,
                    ...] where loc1 etc. are absolute locus indices).
                    This is done through the calculation of allele
                    frequencies. Therefore, allele frequencies will
                    also be calculated if this statistics is
                    requested. This parameter will set the following
                    variables (carray objects of the numbers of
                    alleles for allloci. Unrequested loci will have 0
                    distinct alleles.):
                    * numOfAlleles, subPop[sp]['numOfAlleles'], number
                    of distinct alleles at each locus. (Calculated
                    only at requested loci.)
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
                    will be calculated (genoFreq=[loc1, loc2, ...]
                    where loc1 etc. are loci where genotype
                    frequencies will be calculated). All the genotypes
                    in the  population will be counted. You may use
                    hasPhase to set if a/b and b/a are the same
                    genotype. This parameter will set the following
                    dictionary variables. Note that unlike list used
                    for alleleFreq etc., the indices a, b of
                    genoFreq[a][b] are dictionary keys, so you will
                    get a KeyError when you used a wrong key. Usually,
                    genoNum.setDefault(a,{}) is preferred.
                    * genoNum[a][geno] and
                    subPop[sp]['genoNum'][a][geno], the number of
                    genotype geno at allele a. geno has the form x-y.
                    * genoFreq[a][geno] and
                    subPop[sp]['genoFreq'][a][geno], the frequency of
                    genotype geno at allele a.
    heteroFreq:     an array of loci to calaulate observed
                    heterozygosities and expected heterozygosities
                    (heteroFreq=[loc1, loc2, ...]). This parameter
                    will set the following variables (arrays of
                    observed heterozygosities). Note that
                    heteroNum[loc][1] is the number of heterozygote
                    1x,  $ x \\ne 1 $. Numbers and frequencies
                    (proportions) of heterozygotes are calculated for
                    each allele. HeteroNum[loc] and HeterFreq[loc] are
                    the overall heterozygosity number and frequency.
                    I.e., the number/frequency of genotype xy,  $ x
                    \\ne y $. From this number, we can easily derive
                    the number of homozygosity.
                    * HeteroNum[loc], subPop[sp]['HeteroNum'][loc],
                    the overall heterozygote number
                    * HeteroFreq[loc], subPop[sp]['HeteroFreq'][loc],
                    the overall heterozygote frequency
                    * heteroNum[loc][allele],
                    subPop[sp]['heteroNum'][loc][allele]
                    * heteroFreq[loc][allele],
                    subPop[sp]['heteroFreq'][loc][allele]
    homoFreq:       an array of loci to calaulate observed
                    homozygosities and expected homozygosities
                    (homoFreq=[loc1, loc2, ...]). This parameter will
                    calculate the numbers and frequencies of
                    homozygotes xx and set the following variables:
                    * homoNum[loc], subPop[sp]['homoNum'][loc],
                    * homoFreq[loc], subPop[sp]['homoFreq'][loc].
    expHetero:      an array of loci at which the expected
                    heterozygosities will be calculated
                    (expHetero=[loc1, loc2, ...]). The expected
                    heterozygosity is calculated by
                    $ h_{exp}=1-p_{i}^{2}. $ The following variables
                    will be set:
                    * expHetero[loc], subPop[sp]['expHetero'][loc].
    expHetero_param:a dictionary of parameters of expHetero
                    statistics. Can be one or more items choosen from
                    the following options: subpop and midValues.
    haploFreq:      a matrix of haplotypes (allele sequences on
                    different loci) to count. For example, haploFreq =
                    [ [ 0,1,2 ], [1,2] ] will count all haplotypes on
                    loci 0,1 and 2; and all haplotypes on loci 1, 2.
                    If only one haplotype is specified, the outer []
                    can be omitted. I.e., haploFreq=[0,1] is
                    acceptable. The following dictionary variables
                    will be set with keys 0-1-2 etc. For example,
                    haploNum['1-2']['5-6'] is the number of allele
                    pair 5,6 (on loci 1 and 2 respectively) in the
                    population.
                    * haploNum[haplo] and
                    subPop[sp]['haploNum'][haplo], the number of
                    allele sequencies on loci haplo.
                    * haploFreq[haplo],
                    subPop[sp]['haploFreq'][haplo], the frequency of
                    allele sequencies on loci haplo.
    LD:             calculate linkage disequilibria  $ LD $,  $ LD' $
                    and  $ r^{2} $, given LD=[ [loc1, loc2], [ loc1,
                    loc2, allele1, allele2], ... ] For each item
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
                    {B}\\left(1-P_{B}\\right)} $
    LD_param:       a dictionary of parameters of LD statistics. Can
                    have key  stat which is a list of statistics to
                    calculate. Default to all. If any statistics is
                    specified, only those specified will be
                    calculated. For example, you may use
                    LD_param={LD_prime} to calculate D' only, where
                    LD_prime is a shortcut for 'stat':['LD_prime'].
                    Other parameters that you may use are:
                    * subPop, whether or not calculate statistics for
                    subpopulations
                    * midValues, whether or not keep intermediate
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
                    Cockerham:
                    *  $ F $ ( $ F_{IT} $) the correlation of genes
                    within individuals (inbreeding);
                    *  $ \\theta $ ( $ F_{ST} $) the correlation of
                    genes of difference individuals in the same
                    population (will evaluate for each subpopulation
                    and the whole  population)
                    *  $ f $ ( $ F_{IS} $) the correlation of genes
                    within individuals within populations. Population
                    refers to subpopulations in  simuPOP term.??? This
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
                    are discarded. relatedness[grp1][grp2] is the
                    relatedness value between grp1 and grp2. There is
                    no subpopulation level relatedness values.
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

Description:

    simuPOP::stat::~stat

Usage:

    x.~stat()

"; 

%feature("docstring") simuPOP::stat::clone "

Description:

    deep copy of a  stat operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::stat::apply "

Description:

    apply the  stat operator

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::stat::__repr__ "

Description:

    used by Python print function to print out the general information
    of the  stat operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::statAlleleFreq;

%feature("docstring") simuPOP::statAlleleFreq::statAlleleFreq "

Description:

    simuPOP::statAlleleFreq::statAlleleFreq

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

Description:

    simuPOP::statAlleleFreq::addLocus

Usage:

    x.addLocus(locus, post, subPop, numOfAlleles)

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleNumAll "

Description:

    simuPOP::statAlleleFreq::alleleNumAll

Usage:

    x.alleleNumAll()

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleNumVec "

Description:

    simuPOP::statAlleleFreq::alleleNumVec

Usage:

    x.alleleNumVec(loc)

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleNum "

Description:

    simuPOP::statAlleleFreq::alleleNum

Usage:

    x.alleleNum(allele, loc)

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleFreqAll "

Description:

    simuPOP::statAlleleFreq::alleleFreqAll

Usage:

    x.alleleFreqAll()

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleFreqVec "

Description:

    simuPOP::statAlleleFreq::alleleFreqVec

Usage:

    x.alleleFreqVec(loc)

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleleFreq "

Description:

    simuPOP::statAlleleFreq::alleleFreq

Usage:

    x.alleleFreq(allele, loc)

"; 

%feature("docstring") simuPOP::statAlleleFreq::numOfAlleles "

Description:

    simuPOP::statAlleleFreq::numOfAlleles

Usage:

    x.numOfAlleles()

"; 

%feature("docstring") simuPOP::statAlleleFreq::alleles "

Description:

    simuPOP::statAlleleFreq::alleles

Usage:

    x.alleles(loc)

"; 

%feature("docstring") simuPOP::statAlleleFreq::apply "

Description:

    simuPOP::statAlleleFreq::apply

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statAssociation;

%feature("docstring") simuPOP::statAssociation::statAssociation "

Description:

    simuPOP::statAssociation::statAssociation

Usage:

    statAssociation(alleleFreq, haploFreq, Association=[], param={})

"; 

%feature("docstring") simuPOP::statAssociation::apply "

Description:

    simuPOP::statAssociation::apply

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statExpHetero;

%feature("docstring") simuPOP::statExpHetero::statExpHetero "

Description:

    simuPOP::statExpHetero::statExpHetero

Usage:

    statExpHetero(alleleFreq, expHetero=[], param={})

"; 

%feature("docstring") simuPOP::statExpHetero::apply "

Description:

    simuPOP::statExpHetero::apply

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statFst;

%feature("docstring") simuPOP::statFst::statFst "

Description:

    simuPOP::statFst::statFst

Usage:

    statFst(alleleFreq, heteroFreq, Fst=[], param={})

"; 

%feature("docstring") simuPOP::statFst::Fst "

Description:

    simuPOP::statFst::Fst

Usage:

    x.Fst()

"; 

%feature("docstring") simuPOP::statFst::Fis "

Description:

    simuPOP::statFst::Fis

Usage:

    x.Fis()

"; 

%feature("docstring") simuPOP::statFst::Fit "

Description:

    simuPOP::statFst::Fit

Usage:

    x.Fit()

"; 

%feature("docstring") simuPOP::statFst::apply "

Description:

    simuPOP::statFst::apply

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statGenoFreq;

%feature("docstring") simuPOP::statGenoFreq::statGenoFreq "

Description:

    simuPOP::statGenoFreq::statGenoFreq

Usage:

    statGenoFreq(genoFreq=[], phase=True)

"; 

%feature("docstring") simuPOP::statGenoFreq::apply "

Description:

    simuPOP::statGenoFreq::apply

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statHaploFreq;

%feature("docstring") simuPOP::statHaploFreq::statHaploFreq "

Description:

    simuPOP::statHaploFreq::statHaploFreq

Usage:

    statHaploFreq(haploFreq=[])

"; 

%feature("docstring") simuPOP::statHaploFreq::~statHaploFreq "

Description:

    simuPOP::statHaploFreq::~statHaploFreq

Usage:

    x.~statHaploFreq()

"; 

%feature("docstring") simuPOP::statHaploFreq::addHaplotype "

Description:

    simuPOP::statHaploFreq::addHaplotype

Usage:

    x.addHaplotype(haplo, post=False)

"; 

%feature("docstring") simuPOP::statHaploFreq::numOfHaplotypes "

Description:

    simuPOP::statHaploFreq::numOfHaplotypes

Usage:

    x.numOfHaplotypes(haplo)

"; 

%feature("docstring") simuPOP::statHaploFreq::haploNum "

Description:

    simuPOP::statHaploFreq::haploNum

Usage:

    x.haploNum(haplo)

"; 

%feature("docstring") simuPOP::statHaploFreq::haploFreq "

Description:

    simuPOP::statHaploFreq::haploFreq

Usage:

    x.haploFreq(haplo)

"; 

%feature("docstring") simuPOP::statHaploFreq::apply "

Description:

    simuPOP::statHaploFreq::apply

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statHeteroFreq;

%feature("docstring") simuPOP::statHeteroFreq::statHeteroFreq "

Description:

    simuPOP::statHeteroFreq::statHeteroFreq

Usage:

    statHeteroFreq(heteroFreq=[], homoFreq=[])

"; 

%feature("docstring") simuPOP::statHeteroFreq::addLocus "

Description:

    simuPOP::statHeteroFreq::addLocus

Usage:

    x.addLocus(locus, post=False)

"; 

%feature("docstring") simuPOP::statHeteroFreq::heteroNum "

Description:

    simuPOP::statHeteroFreq::heteroNum

Usage:

    x.heteroNum(allele, loc)

"; 

%feature("docstring") simuPOP::statHeteroFreq::heteroFreq "

Description:

    simuPOP::statHeteroFreq::heteroFreq

Usage:

    x.heteroFreq(allele, loc)

"; 

%feature("docstring") simuPOP::statHeteroFreq::apply "

Description:

    simuPOP::statHeteroFreq::apply

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statLD;

%feature("docstring") simuPOP::statLD::statLD "

Description:

    simuPOP::statLD::statLD

Usage:

    statLD(alleleFreq, haploFreq, LD=[], LD_param={})

"; 

%feature("docstring") simuPOP::statLD::apply "

Description:

    simuPOP::statLD::apply

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statNumOfAffected;

%feature("docstring") simuPOP::statNumOfAffected::statNumOfAffected "

Description:

    simuPOP::statNumOfAffected::statNumOfAffected

Usage:

    statNumOfAffected(numOfAffected=False, param={})

"; 

%feature("docstring") simuPOP::statNumOfAffected::~statNumOfAffected "

Description:

    simuPOP::statNumOfAffected::~statNumOfAffected

Usage:

    x.~statNumOfAffected()

"; 

%feature("docstring") simuPOP::statNumOfAffected::activate "

Description:

    simuPOP::statNumOfAffected::activate

Usage:

    x.activate(yes=True)

"; 

%feature("docstring") simuPOP::statNumOfAffected::numOfAffected "

Description:

    simuPOP::statNumOfAffected::numOfAffected

Usage:

    x.numOfAffected()

"; 

%feature("docstring") simuPOP::statNumOfAffected::numOfUnaffected "

Description:

    simuPOP::statNumOfAffected::numOfUnaffected

Usage:

    x.numOfUnaffected()

"; 

%feature("docstring") simuPOP::statNumOfAffected::apply "

Description:

    simuPOP::statNumOfAffected::apply

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statNumOfAlleles;

%feature("docstring") simuPOP::statNumOfAlleles::statNumOfAlleles "

Description:

    simuPOP::statNumOfAlleles::statNumOfAlleles

Usage:

    statNumOfAlleles(calc, atLoci=[], param={})

"; 

%feature("docstring") simuPOP::statNumOfAlleles::~statNumOfAlleles "

Description:

    simuPOP::statNumOfAlleles::~statNumOfAlleles

Usage:

    x.~statNumOfAlleles()

"; 

%feature("docstring") simuPOP::statNumOfAlleles::apply "

Description:

    simuPOP::statNumOfAlleles::apply

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statNumOfMale;

%feature("docstring") simuPOP::statNumOfMale::statNumOfMale "

Description:

    simuPOP::statNumOfMale::statNumOfMale

Usage:

    statNumOfMale(numOfMale=False, param={})

"; 

%feature("docstring") simuPOP::statNumOfMale::activate "

Description:

    simuPOP::statNumOfMale::activate

Usage:

    x.activate(yes=True)

"; 

%feature("docstring") simuPOP::statNumOfMale::numOfMale "

Description:

    simuPOP::statNumOfMale::numOfMale

Usage:

    x.numOfMale()

"; 

%feature("docstring") simuPOP::statNumOfMale::numOfFemale "

Description:

    simuPOP::statNumOfMale::numOfFemale

Usage:

    x.numOfFemale()

"; 

%feature("docstring") simuPOP::statNumOfMale::apply "

Description:

    simuPOP::statNumOfMale::apply

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::stator "

Description:

    basic class of all the statistics

Details:

    Operator  stator calculate various basic statistics for the
    population and set variables in the local namespace. Other
    operators/functions can refer to the results from the namespace
    after  stat is applied. Stat is the function form of the operator.
    ????

Details:

    Operator  stator calculate various basic statistics for the
    population and set variables in the local namespace. Other
    operators/functions can refer to the results from the namespace
    after  stat is applied. Stat is the function form of the operator.
    ????

"; 

%feature("docstring") simuPOP::stator::stator "

Description:

    create a  stator

Usage:

    stator(output=\"\", outputExpr=\"\", stage=PostMating, begin=0,
      end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::stator::~stator "

Description:

    destructor

Usage:

    x.~stator()

"; 

%feature("docstring") simuPOP::stator::clone "

Description:

    deep copy of a  stator

Usage:

    x.clone()

"; 

%ignore simuPOP::statPopSize;

%feature("docstring") simuPOP::statPopSize::statPopSize "

Description:

    simuPOP::statPopSize::statPopSize

Usage:

    statPopSize(popSize=False)

"; 

%feature("docstring") simuPOP::statPopSize::activate "

Description:

    simuPOP::statPopSize::activate

Usage:

    x.activate()

"; 

%feature("docstring") simuPOP::statPopSize::apply "

Description:

    simuPOP::statPopSize::apply

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
                    subpopulations. To specify between  individual
                    relatedness, use [[1],[2],[3]] (the first form).
                    If this parameter is ignored, this operator
                    calculate relatedness between all subpopulations.
    method:         can be REL_Queller, REL_Lynch, REL_IR, REL_D2 or
                    REL_Rel. Please refer to the manual for details.

"; 

%feature("docstring") simuPOP::statRelatedness::relQueller "

Description:

    simuPOP::statRelatedness::relQueller

Usage:

    x.relQueller(ind1, ind2)

"; 

%feature("docstring") simuPOP::statRelatedness::relLynch "

Description:

    simuPOP::statRelatedness::relLynch

Usage:

    x.relLynch(ind1, ind2)

"; 

%feature("docstring") simuPOP::statRelatedness::relIR "

Description:

    simuPOP::statRelatedness::relIR

Usage:

    x.relIR(ind1, locus)

"; 

%feature("docstring") simuPOP::statRelatedness::relD2 "

Description:

    simuPOP::statRelatedness::relD2

Usage:

    x.relD2(ind1, locus)

"; 

%feature("docstring") simuPOP::statRelatedness::relRel "

Description:

    simuPOP::statRelatedness::relRel

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

Description:

    simuPOP::statRelatedness::apply

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::StopIteration "

Description:

    exception, thrown if out of memory

"; 

%feature("docstring") simuPOP::StopIteration::StopIteration "

Description:

    simuPOP::StopIteration::StopIteration

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

Description:

    simuPOP::StreamProvider::~StreamProvider

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

Description:

    simuPOP::SystemError::SystemError

Usage:

    SystemError(msg)

"; 

%feature("docstring") simuPOP::tagger "

Description:

    basic class of tagging individuals

Details:

    tagger is a during  mating operator that tag individuals with
    various information. Potential usages are:
    * recording parental information to track pedigree;
    * tagging an individual/allele and monitor its  spread in the
    population etc.

Details:

    tagger is a during  mating operator that tag individuals with
    various information. Potential usages are:
    * recording parental information to track pedigree;
    * tagging an individual/allele and monitor its  spread in the
    population etc.

"; 

%feature("docstring") simuPOP::tagger::tagger "

Description:

    create a  tagger, default to be always active but no output

Usage:

    tagger(begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[])

"; 

%feature("docstring") simuPOP::tagger::~tagger "

Description:

    destructor

Usage:

    x.~tagger()

"; 

%feature("docstring") simuPOP::tagger::clone "

Description:

    deep copy of a \\  tagger

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::terminateIf "

Description:

    terminate according to a condition

Details:

    This operator terminates the evolution under certain conditions.
    For example,  terminateIf(condition='alleleFreq[0][1]<0.05',
    begin=100) terminates the evolution if the allele frequency of
    allele 1 at locus 0 is less than 0.05. Of course, to make this
    opertor work, you will need to use a  stat operator before it so
    that variable alleleFreq exists in the local namespace.
    When the condition is true, a shared variable var=\"terminate\" will
    be set to the current generation.

Details:

    This operator terminates the evolution under certain conditions.
    For example,  terminateIf(condition='alleleFreq[0][1]<0.05',
    begin=100) terminates the evolution if the allele frequency of
    allele 1 at locus 0 is less than 0.05. Of course, to make this
    opertor work, you will need to use a  stat operator before it so
    that variable alleleFreq exists in the local namespace.
    When the condition is true, a shared variable var=\"terminate\" will
    be set to the current generation.

"; 

%feature("docstring") simuPOP::terminateIf::terminateIf "

Description:

    create a  terminateIf terminator

Usage:

    terminateIf(condition=\"\", message=\"\", var=\"terminate\",
      output=\"\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

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

Description:

    simuPOP::terminateIf::~terminateIf

Usage:

    x.~terminateIf()

"; 

%feature("docstring") simuPOP::terminator "

Description:

    terminate the evolution

Details:

    These operators are used to see if an evolution is running as
    expected, and terminate the evolution if a certain condition
    fails.

Details:

    These operators are used to see if an evolution is running as
    expected, and terminate the evolution if a certain condition
    fails.

"; 

%feature("docstring") simuPOP::terminator::terminator "

Description:

    create a  terminator, default to be always active

Usage:

    terminator(message=\"\", output=\">\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::terminator::~terminator "

Description:

    destructor

Usage:

    x.~terminator()

"; 

%feature("docstring") simuPOP::terminator::clone "

Description:

    deep copy of a  terminator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::terminator::message "

Description:

    return the message to print when terminated???

Usage:

    x.message()

"; 

%feature("docstring") simuPOP::ticToc "

Function form:

    TicToc

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
      end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

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

    apply the  ticToc operator to one  population

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

Function form:

    TurnOffDebug

Description:

    set debug off

Details:

    Turn off debug.

Details:

    Turn off debug.

"; 

%feature("docstring") simuPOP::turnOffDebug::turnOffDebug "

Description:

    simuPOP::turnOffDebug::turnOffDebug

Usage:

    turnOffDebug(code, stage=PreMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

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

    apply the  turnOffDebug operator to one  population

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

Function form:

    TurnOnDebug

Description:

    set debug on

Details:

    Turn on debug. There are several ways to turn on debug information
    for non-optimized modules, namely
    * set environment variable SIMUDEBUG
    * use simuOpt.setOptions(debug) function, or
    * use TurnOnDebug or TurnOnDebugByName function
    * use this  turnOnDebug operator The advantage of using this
    operator is that you can turn on debug at given generations.

Details:

    Turn on debug. There are several ways to turn on debug information
    for non-optimized modules, namely
    * set environment variable SIMUDEBUG
    * use simuOpt.setOptions(debug) function, or
    * use TurnOnDebug or TurnOnDebugByName function
    * use this  turnOnDebug operator The advantage of using this
    operator is that you can turn on debug at given generations.

"; 

%feature("docstring") simuPOP::turnOnDebug::turnOnDebug "

Description:

    simuPOP::turnOnDebug::turnOnDebug

Usage:

    turnOnDebug(code, stage=PreMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

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

    apply the  turnOnDebug operator to one  population

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

Description:

    simuPOP::TypeError::TypeError

Usage:

    TypeError(msg)

"; 

%feature("docstring") simuPOP::ValueError "

Description:

    exception, thrown if value of range etc

"; 

%feature("docstring") simuPOP::ValueError::ValueError "

Description:

    simuPOP::ValueError::ValueError

Usage:

    ValueError(msg)

"; 

%feature("docstring") simuPOP::Weightedsampler "

Description:

    simuPOP::Weightedsampler

"; 

%feature("docstring") simuPOP::Weightedsampler::Weightedsampler "

Description:

    simuPOP::Weightedsampler::Weightedsampler

Usage:

    Weightedsampler(rng, weight=[], fast=True)

"; 

%feature("docstring") simuPOP::Weightedsampler::~Weightedsampler "

Description:

    simuPOP::Weightedsampler::~Weightedsampler

Usage:

    x.~Weightedsampler()

"; 

%feature("docstring") simuPOP::Weightedsampler::set "

Description:

    simuPOP::Weightedsampler::set

Usage:

    x.set(weight)

"; 

%feature("docstring") simuPOP::Weightedsampler::biSearch "

Description:

    simuPOP::Weightedsampler::biSearch

Usage:

    x.biSearch(a)

"; 

%feature("docstring") simuPOP::Weightedsampler::get "

Description:

    simuPOP::Weightedsampler::get

Usage:

    x.get()

"; 

%feature("docstring") simuPOP::Weightedsampler::q "

Description:

    simuPOP::Weightedsampler::q

Usage:

    x.q()

"; 

%feature("docstring") simuPOP::Weightedsampler::a "

Description:

    simuPOP::Weightedsampler::a

Usage:

    x.a()

"; 

%ignore simuPOP::countAlleles(population &pop, int subpop, const vectori &loci, const vectori &alleles, vectorlu &numAllele);

%ignore simuPOP::getExpectedAlleles(population &pop, vectorf &expFreq, const vectori &loci, const vectori &alleles, vectoru &expAlleles);

%feature("docstring") simuPOP::FreqTrajectoryStoch "

Description:

    simuPOP::FreqTrajectoryStoch

Usage:

    FreqTrajectoryStoch(curGen=0, freq=0, N=0, NtFunc=None,
      fitness=[], fitnessFunc=None, minMutAge=0, maxMutAge=100000,
      ploidy=2, restartIfFail=False, maxAttempts=1000,
      allowFixation=False)

"; 

%ignore simuPOP::MarginalFitness(unsigned nLoci, const vectorf &fitness, const vectorf &freq);

%feature("docstring") simuPOP::FreqTrajectoryMultiStoch "

Description:

    simuPOP::FreqTrajectoryMultiStoch

Usage:

    FreqTrajectoryMultiStoch(curGen=0, freq=[], N=0, NtFunc=None,
      fitness=[], fitnessFunc=None, minMutAge=0, maxMutAge=100000,
      ploidy=2, restartIfFail=False, maxAttempts=1000)

"; 

%feature("docstring") simuPOP::FreqTrajectorySelSim "

Description:

    simuPOP::FreqTrajectorySelSim

Usage:

    FreqTrajectorySelSim(sel, Ne, freq, dom_h, selection)

"; 

%feature("docstring") simuPOP::FreqTrajectoryForward "

Description:

    simuPOP::FreqTrajectoryForward

Usage:

    FreqTrajectoryForward(lowbound, highbound, disAge, grate, N0,
      seleCo)

"; 

%feature("docstring") simuPOP::LoadPopulation "

Description:

    load a  population from a file. The file format is by default
    determined by file extension (format=\"auto\"). Otherwise, format
    can be one of txt, bin, or xml.

Usage:

    LoadPopulation(file, format=\"auto\")

"; 

%feature("docstring") simuPOP::testGetinfoFromInd "

Description:

    get info through ind.info()

Usage:

    testGetinfoFromInd(pop)

"; 

%feature("docstring") simuPOP::testGetinfoFromPop "

Description:

    get info through GappedInfoIterator

Usage:

    testGetinfoFromPop(pop, order)

"; 

%feature("docstring") simuPOP::LoadSimulator "

Description:

    load a  simulator from a file with the specified  mating scheme.
    The file format is by default determined by file extension
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

Description:

    simuPOP::TurnOnDebug

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

    show all debug codes (print to cout)

Usage:

    ListDebugCode()

"; 

%ignore simuPOP::dbgString(DBG_CODE code);

%ignore simuPOP::simuPOP_kbhit();

%ignore simuPOP::simuPOP_getch();

%feature("docstring") simuPOP::testGappedIterator "

Description:

    simuPOP::testGappedIterator

Usage:

    testGappedIterator()

"; 

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

%ignore simuPOP::Int_Vec_As_NumArray(vectori::iterator begin, vectori::iterator end);

%ignore simuPOP::Double_Vec_As_NumArray(vectorf::iterator begin, vectorf::iterator end);

%ignore simuPOP::Allele_Vec_As_NumArray(GenoIterator begin, GenoIterator end);

%ignore simuPOP::Info_Vec_As_NumArray(InfoIterator begin, InfoIterator end);

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

%feature("docstring") simuPOP::setRNG "

Description:

    for backward compatibilit, will remove later

Usage:

    setRNG(rng=\"\", seed=0)

"; 

%feature("docstring") simuPOP::ListAllRNG "

Description:

    list all available random number generators

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

%feature("docstring") simuPOP::compileCompiler "

Description:

    return the compiler used to compile this  simuPOP module

Usage:

    compileCompiler()

"; 

%feature("docstring") simuPOP::compileDate "

Description:

    return the date when this  simuPOP module is compiled

Usage:

    compileDate()

"; 

%feature("docstring") simuPOP::compilePyVersion "

Description:

    return the Python version this  simuPOP module is compiled for

Usage:

    compilePyVersion()

"; 

%feature("docstring") simuPOP::compilePlatForm "

Description:

    return the platform on which this  simuPOP module is compiled

Usage:

    compilePlatForm()

"; 

%ignore simuPOP::initialize();

%feature("docstring") simuPOP::optimized "

Description:

    return True if this  simuPOP module is optimized

Usage:

    optimized()

"; 

%feature("docstring") simuPOP::mpi "

Description:

    simuPOP::mpi

Usage:

    mpi()

"; 

%feature("docstring") simuPOP::limits "

Description:

    print out system limits

Usage:

    limits()

"; 

%feature("docstring") simuPOP::mpiRank "

Description:

    simuPOP::mpiRank

Usage:

    mpiRank()

"; 

%ignore simuPOP::mpiSize();

%ignore simuPOP::mpiBarrier();

%feature("docstring") simuPOP::alleleType "

Description:

    return the allele type of the current module. Can be binary,
    short, or long.

Usage:

    alleleType()

"; 

%feature("docstring") simuPOP::maxAllele "

Description:

    return  $ 1 $,  $ 2^8-1 $,  $ 2^{16}-1 $ for binary, short, or
    long allele modules, respectively

Usage:

    maxAllele()

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

