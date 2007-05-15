%feature("docstring") simuPOP::affectedSibpairSample "

Description:

    thrink population accroding to some outside value

"; 

%feature("docstring") simuPOP::affectedSibpairSample::affectedSibpairSample "

Description:

    draw cases and controls

Usage:

    affectedSibpairSample(size=[], chooseUnaffected=False,
      countOnly=False, name=\"sample\", nameExpr=\"\", times=1, saveAs=\"\",
      saveAsExpr=\"\", format=\"auto\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[\"father_idx\", \"mother_idx\"])

Arguments:

    size:           number of affected sibpairs to be sampled. Can be
                    a number or an array. If a number is given, it is
                    the total number of sibpairs, ignoring population
                    structure. Otherwise, given number of sibpairs are
                    sampled from subpopulations. If size is
                    unspecified, this operator will return all
                    affected sibpairs.
    countOnly:      set variables about number of affected sibpairs,
                    do not actually draw the sample
    name:           variable name of the sampled population (will be
                    put in pop local namespace)
    nameExpr:       expression version of name. If both name and
                    nameExpr is empty, do not store pop.
    times:          how many times to run the sample process? This is
                    usually one, but we may want to take several
                    random samples.
    saveAs:         filename to save the population.
    saveAsExpr:     expression for save filename
    format:         to save sample(s)
    stage:          and other parameters please see
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::affectedSibpairSample::~affectedSibpairSample "

Description:

    destructor

Usage:

    x.~affectedSibpairSample()

"; 

%feature("docstring") simuPOP::affectedSibpairSample::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::affectedSibpairSample::prepareSample "

Description:

    simuPOP::affectedSibpairSample::prepareSample

Usage:

    x.prepareSample(pop)

"; 

%feature("docstring") simuPOP::affectedSibpairSample::drawsample "

Details:

    collect all families

Usage:

    x.drawsample(pop)

"; 

%feature("docstring") simuPOP::affectedSibpairSample::__repr__ "

Description:

    simuPOP::affectedSibpairSample::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::BernulliTrials "

Description:

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

%ignore simuPOP::BernulliTrials::trialSize() const ;

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

    get a trial corresponding to m_prob.

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

Details:

    binomial random selectionNo sex. Choose one individual from last
    generation.1. numOffspring protocol is honored 2. population size
    changes are allowed 3. selection is possible.So this works just
    like a sexless random mating. If ploidy is one, this is
    chromosomal mating.

"; 

%feature("docstring") simuPOP::binomialSelection::binomialSelection "

Description:

    constructor

Usage:

    binomialSelection(numOffspring=1., *numOffspringFunc=NULL,
      maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
      newSubPopSizeExpr=\"\", *newSubPopSizeFunc=NULL)

"; 

%feature("docstring") simuPOP::binomialSelection::~binomialSelection "

Description:

    destructor

Usage:

    x.~binomialSelection()

"; 

%feature("docstring") simuPOP::binomialSelection::clone "

Description:

    clone() const. The same as  mating::clone() const.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::binomialSelection::__repr__ "

Description:

    return name of the mating type

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::binomialSelection::submitScratch "

Description:

    simuPOP::binomialSelection::submitScratch

Usage:

    x.submitScratch(pop, scratch)

"; 

%feature("docstring") simuPOP::binomialSelection::mate "

Description:

    do the mating.

Usage:

    x.mate(pop, scratch, ops, submit)

Arguments:

    pop:            population
    scratch:        scratch population, will be used in this mating
                    scheme.
    ops:            during mating operators

"; 

%feature("docstring") simuPOP::caseControlSample "

Description:

    thrink population accroding to some outside value

"; 

%feature("docstring") simuPOP::caseControlSample::caseControlSample "

Description:

    draw cases and controls

Usage:

    caseControlSample(cases=[], controls=[], spSample=False,
      name=\"sample\", nameExpr=\"\", times=1, saveAs=\"\", saveAsExpr=\"\",
      format=\"auto\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    cases:          number of cases, or an array of number of cases
                    from each subpopulation.
    controls:       number of controls, or an array of number of
                    controls from each subpopulation.
    name:           variable name of the sampled population (will be
                    put in pop local namespace)
    nameExpr:       expression version of name. If both name and
                    nameExpr is empty, do not store pop.
    times:          how many times to run the sample process? This is
                    usually one, but we may want to take several
                    random samples.
    saveAs:         filename to save the population.
    saveAsExpr:     expression for save filename
    format:         to save sample(s)
    stage:          and other parameters please see
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::caseControlSample::~caseControlSample "

Description:

    destructor

Usage:

    x.~caseControlSample()

"; 

%feature("docstring") simuPOP::caseControlSample::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::caseControlSample::prepareSample "

Description:

    simuPOP::caseControlSample::prepareSample

Usage:

    x.prepareSample(pop)

"; 

%feature("docstring") simuPOP::caseControlSample::drawsample "

Description:

    simuPOP::caseControlSample::drawsample

Usage:

    x.drawsample(pop)

"; 

%feature("docstring") simuPOP::caseControlSample::__repr__ "

Description:

    simuPOP::caseControlSample::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::continueIf "

Description:

    terminate according to a condition which can be, e.g.
    any(alleleNum0) == 0 all(alleleNum1) > 0.5 alleleNum0{2} == 0 etc.
    When the condition is true, a shared variable var=\"terminate\" will
    be set to current generation.

"; 

%feature("docstring") simuPOP::continueIf::continueIf "

Description:

    simuPOP::continueIf::continueIf

Usage:

    continueIf(condition=\"\", message=\"\", var=\"terminate\", output=\"\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::continueIf::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::continueIf::__repr__ "

Description:

    simuPOP::continueIf::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::continueIf::apply "

Description:

    check all alleles in vector allele if they are fixed.

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

Details:

    binomial random selectionNo sex. Choose one individual from last
    generation.1. numOffspring protocol is honored 2. population size
    changes are allowed 3. selection is possible.So this works just
    like a sexless random mating. If ploidy is one, this is
    chromosomal mating.

"; 

%feature("docstring") simuPOP::controlledBinomialSelection::controlledBinomialSelection "

Description:

    constructor

Usage:

    controlledBinomialSelection(loci, alleles, *freqFunc,
      numOffspring=1., *numOffspringFunc=NULL, maxNumOffspring=0,
      mode=MATE_NumOffspring, newSubPopSize=[], newSubPopSizeExpr=\"\",
      *newSubPopSizeFunc=NULL)

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

    clone() const. The same as  mating::clone() const.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::controlledBinomialSelection::__repr__ "

Description:

    return name of the mating type

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::controlledBinomialSelection::submitScratch "

Description:

    simuPOP::controlledBinomialSelection::submitScratch

Usage:

    x.submitScratch(pop, scratch)

"; 

%feature("docstring") simuPOP::controlledBinomialSelection::mate "

Description:

    do the mating.

Usage:

    x.mate(pop, scratch, ops, submit)

Arguments:

    pop:            population
    scratch:        scratch population, will be used in this mating
                    scheme.
    ops:            during mating operators

"; 

%feature("docstring") simuPOP::controlledMating "

Details:

    controlled mating

"; 

%feature("docstring") simuPOP::controlledMating::controlledMating "

Description:

    Controlled mating, control allele frequency at a locus.

Usage:

    controlledMating(matingScheme, loci, alleles, *freqFunc,
      range=0.01)

Arguments:

    mating:         a mating scheme.
    loci:           loci at which allele frequency is controlled. Note
                    that controlling the allele frequencies at several
                    loci may take a long time.
    alleles:        alleles to control at each loci. Should have the
                    same length as loci
    freqFunc:       frequency boundaries. If the length of the return
                    value equals the size of loci, the range for loci
                    will be [value0, value0+range], [value1,
                    value1+range] etc. If the length of the return
                    value is 2 times size of loci, it will be
                    interpreted as [low1, high1, low2, high2 ...]

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

    clone() const. Generate a copy of itself and return pointer this
    is to make sure the object is persistent and will not be freed by
    python.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::controlledMating::isCompatible "

Description:

    check if the mating type is compatible with population structure

Details:

    possible things to check:
    * need certain types of individual (age, sex etc)
    * need resizeable population...

Usage:

    x.isCompatible(pop)

"; 

%feature("docstring") simuPOP::controlledMating::__repr__ "

Description:

    return name of the mating type

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::controlledMating::mate "

Description:

    mate: This is not supposed to be called for base mating class.

Usage:

    x.mate(pop, scratch, ops, submit)

Arguments:

    pop:            population
    scratch:        scratch population
    ops:            during mating operators

"; 

%feature("docstring") simuPOP::controlledRandomMating "

Details:

    basic sexual random mating.Within each subpopulation, choose male
    and female randomly randmly get one copy of chromosome from
    father/mother.require: sexed individual; ploidy == 2apply during
    mating operators and put into the next generation.if
    ignoreParentsSex is set, parents will be chosen regardless of
    sex.Otherwise, male and female will be collected and be chosen
    randomly.If there is no male or female in a subpopulation, if
    m_UseSameSexIfUniSex is true, an warning will be generated and
    same sex mating (?) will be used otherwise,  randomMating will
    return false.if there is no during mating operator to copy
    alleles, a direct copy will be used.

"; 

%feature("docstring") simuPOP::controlledRandomMating::controlledRandomMating "

Description:

    create a random mating scheme

Usage:

    controlledRandomMating(loci, alleles, *freqFunc, acceptScheme=0,
      numOffspring=1., *numOffspringFunc=NULL, maxNumOffspring=0,
      mode=MATE_NumOffspring, newSubPopSize=[],
      *newSubPopSizeFunc=NULL, newSubPopSizeExpr=\"\",
      contWhenUniSex=True)

Arguments:

    numOffspring:   
    number:         of offspring or p in some modes
    numOffspringFunc:
    a:              python function that determine number of offspring
                    or p depending on mode
    maxNumOffspring:used when mode=MATE_BinomialDistribution
    mode:           one of MATE_NumOffspring ,
                    MATE_NumOffspringEachFamily,
                    MATE_GeometricDistribution,
                    MATE_PoissonDistribution,
                    MATE_BinomialDistribution
    newSubPopSize:  an array of subpopulation sizes, should have the
                    same number of subpopulations as current
                    population
    newSubPopSizeExpr:an expression that will be evaluated as an array
                    of subpop sizes
    newSubPopSizeFunc:an function that have parameter gen and oldSize
                    (current subpop size).
    contWhenUniSex: continue when there is only one sex in the
                    population, default to true

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

    clone() const. Generate a copy of itself and return pointer this
    is to make sure the object is persistent and will not be freed by
    python.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::controlledRandomMating::isCompatible "

Description:

    check if the mating type is compatible with population structure

Details:

    possible things to check:
    * need certain types of individual (age, sex etc)
    * need resizeable population...

Usage:

    x.isCompatible(pop)

"; 

%feature("docstring") simuPOP::controlledRandomMating::__repr__ "

Description:

    return name of the mating type

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::controlledRandomMating::submitScratch "

Description:

    simuPOP::controlledRandomMating::submitScratch

Usage:

    x.submitScratch(pop, scratch)

"; 

%feature("docstring") simuPOP::controlledRandomMating::mate "

Description:

    do the mating. parameters see  mating::mate .

Details:

    Within each subpopulation, choose male and female randomly randmly
    get one copy of chromosome from father/mother.require: sexed
    individual; ploidy == 2apply during mating operators and put into
    the next generation.Otherwise, male and female will be collected
    and be chosen randomly.
    * If there is no male or female in a subpopulation,
    * if m_contWhenUniSex is true, an warning will be generated and
    same sex mating (?) will be used
    * otherwise,  controlledRandomMating will return false.

Usage:

    x.mate(pop, scratch, ops, submit)

"; 

%feature("docstring") simuPOP::dumper "

Description:

    dump the content of a population.

"; 

%feature("docstring") simuPOP::dumper::dumper "

Description:

    dump population

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

    this function is very important

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

Details:

    dump population structuredump all genotypic info

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

    simuPOP::dumper::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::Exception "

Description:

    exception handler. Exceptions will be passed to Python.

"; 

%feature("docstring") simuPOP::Exception::Exception "

Description:

    constructor msg error message

Usage:

    Exception(msg)

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

    Expression(expr=\"\", stmts=\"\", *locals=NULL)

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

%feature("docstring") simuPOP::Expression::evaluate "

Description:

    python expression

Usage:

    x.evaluate()

"; 

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

%feature("docstring") simuPOP::GenoStruTrait "

Description:

    genoStruTrait

Details:

    A trait class that maintain a static array of geno structure, and
    provide interfaces around a  GenoStructure Index.

"; 

%feature("docstring") simuPOP::GenoStruTrait::GenoStruTrait "

Description:

    constructor, but m_genoStruIdx will be set later.

Usage:

    GenoStruTrait()

"; 

%ignore simuPOP::GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru &loci, bool sexChrom, const vectorf &lociPos, const vectorstr &alleleNames, const vectorstr &lociNames, UINT maxAllele, const vectorstr &infoFields, const vectori &chromMap);

%feature("docstring") simuPOP::GenoStruTrait::setGenoStructure "

Description:

    set an existing geno structure, simply use it This is NOT
    efficient! (but has to be used when, for example, loading a
    structure from file

Usage:

    x.setGenoStructure(rhs)

"; 

%ignore simuPOP::GenoStruTrait::setGenoStruIdx(size_t idx);

%ignore simuPOP::GenoStruTrait::mergeGenoStru(size_t idx) const ;

%ignore simuPOP::GenoStruTrait::removeLociFromGenoStru(const vectoru &remove=vectoru(), const vectoru &keep=vectoru());

%ignore simuPOP::GenoStruTrait::insertBeforeLociToGenoStru(const vectoru &idx, const vectorf &pos, const vectorstr &names) const ;

%ignore simuPOP::GenoStruTrait::insertAfterLociToGenoStru(const vectoru &idx, const vectorf &pos, const vectorstr &names) const ;

%ignore simuPOP::GenoStruTrait::genoStru() const ;

%ignore simuPOP::GenoStruTrait::genoStruIdx() const ;

%feature("docstring") simuPOP::GenoStruTrait::ploidy "

Description:

    return ploidy

Usage:

    x.ploidy()

"; 

%feature("docstring") simuPOP::GenoStruTrait::ploidyName "

Description:

    return ploidy

Usage:

    x.ploidyName()

"; 

%feature("docstring") simuPOP::GenoStruTrait::numLoci "

Description:

    number of loci on chromosome chrom

Usage:

    x.numLoci(chrom)

"; 

%feature("docstring") simuPOP::GenoStruTrait::numLoci "

Description:

    simuPOP::GenoStruTrait::numLoci

Usage:

    x.numLoci()

"; 

%feature("docstring") simuPOP::GenoStruTrait::sexChrom "

Description:

    whether or not the last chromosome is sex chromosome

Usage:

    x.sexChrom()

"; 

%feature("docstring") simuPOP::GenoStruTrait::totNumLoci "

Description:

    return totNumLoci (STATIC)

Usage:

    x.totNumLoci()

"; 

%feature("docstring") simuPOP::GenoStruTrait::genoSize "

Description:

    return totNumLoci * ploidy

Usage:

    x.genoSize()

"; 

%feature("docstring") simuPOP::GenoStruTrait::locusPos "

Description:

    locus distance.

Usage:

    x.locusPos(locus)

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociPos "

Description:

    simuPOP::GenoStruTrait::lociPos

Usage:

    x.lociPos()

"; 

%feature("docstring") simuPOP::GenoStruTrait::arrLociPos "

Description:

    expose loci distance

Usage:

    x.arrLociPos()

"; 

%feature("docstring") simuPOP::GenoStruTrait::arrLociPos "

Description:

    expose loci distance of a chromosome

Usage:

    x.arrLociPos(chrom)

"; 

%feature("docstring") simuPOP::GenoStruTrait::numChrom "

Description:

    number of chromosome

Usage:

    x.numChrom()

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromIndex "

Description:

    chromosome index

Usage:

    x.chromIndex()

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromBegin "

Description:

    chromosome index of chromosome chrom

Usage:

    x.chromBegin(chrom)

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromEnd "

Description:

    chromosome index of chromosome chrom

Usage:

    x.chromEnd(chrom)

"; 

%feature("docstring") simuPOP::GenoStruTrait::absLocusIndex "

Description:

    convert from relative locus (on chromsome) to absolute locus (no
    chromosome structure)

Usage:

    x.absLocusIndex(chrom, locus)

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromLocusPair "

Description:

    return chrom, locus pair from an absolute locus position.

Usage:

    x.chromLocusPair(locus)

"; 

%feature("docstring") simuPOP::GenoStruTrait::alleleName "

Description:

    return allele name

Usage:

    x.alleleName(allele)

"; 

%feature("docstring") simuPOP::GenoStruTrait::alleleNames "

Description:

    allele names

Usage:

    x.alleleNames()

"; 

%feature("docstring") simuPOP::GenoStruTrait::locusName "

Description:

    return locus name

Usage:

    x.locusName(loc)

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociNames "

Description:

    simuPOP::GenoStruTrait::lociNames

Usage:

    x.lociNames()

"; 

%feature("docstring") simuPOP::GenoStruTrait::locusByName "

Description:

    return the index of a locus by locus name

Usage:

    x.locusByName(name)

"; 

%feature("docstring") simuPOP::GenoStruTrait::lociByNames "

Description:

    return an array of locus index by loci names

Usage:

    x.lociByNames(names)

"; 

%feature("docstring") simuPOP::GenoStruTrait::maxAllele "

Description:

    simuPOP::GenoStruTrait::maxAllele

Usage:

    x.maxAllele()

"; 

%feature("docstring") simuPOP::GenoStruTrait::setMaxAllele "

Description:

    simuPOP::GenoStruTrait::setMaxAllele

Usage:

    x.setMaxAllele(maxAllele)

"; 

%feature("docstring") simuPOP::GenoStruTrait::hasInfoField "

Description:

    get info length

Usage:

    x.hasInfoField(name)

"; 

%feature("docstring") simuPOP::GenoStruTrait::infoSize "

Description:

    simuPOP::GenoStruTrait::infoSize

Usage:

    x.infoSize()

"; 

%feature("docstring") simuPOP::GenoStruTrait::infoFields "

Description:

    simuPOP::GenoStruTrait::infoFields

Usage:

    x.infoFields()

"; 

%feature("docstring") simuPOP::GenoStruTrait::infoField "

Description:

    simuPOP::GenoStruTrait::infoField

Usage:

    x.infoField(idx)

"; 

%feature("docstring") simuPOP::GenoStruTrait::infoIdx "

Description:

    return the index of field name, return -1 if not found.

Usage:

    x.infoIdx(name)

"; 

%ignore simuPOP::GenoStruTrait::struAddInfoField(const string &field);

%ignore simuPOP::GenoStruTrait::struSetInfoFields(const vectorstr &fields);

%feature("docstring") simuPOP::GenoStruTrait::swap "

Description:

    simuPOP::GenoStruTrait::swap

Usage:

    x.swap(rhs)

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromMap "

Description:

    simuPOP::GenoStruTrait::chromMap

Usage:

    x.chromMap()

"; 

%feature("docstring") simuPOP::gsmMutator "

Description:

    stepwise mutation model.

Details:

    Generalized Stepwise mutation model (GSM) assumes that alleles are
    represented by integer values and that a mutation either increases
    or decreases the allele value by a random value.

"; 

%feature("docstring") simuPOP::gsmMutator::gsmMutator "

Details:

    The generalized stepwise mutation model (GMM) is developed for
    allozymes. It provides better description for these kinds of
    evolutionary processes.

Usage:

    gsmMutator(rate=[], atLoci=[], maxAllele=0, incProb=0.5, p=0,
      *func=NULL, output=\">\", outputExpr=\"\", stage=PostMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[])

Arguments:

    rate::          mutation rate
    incProb:        probability to increase allele state. Default to
                    0.5
    atLoci:         and other parameters: refer to help(mutator),
                    help(baseOperator.__init__)
    func:           return number of steps. no parameter

"; 

%feature("docstring") simuPOP::gsmMutator::~gsmMutator "

Description:

    simuPOP::gsmMutator::~gsmMutator

Usage:

    x.~gsmMutator()

"; 

%feature("docstring") simuPOP::gsmMutator::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::gsmMutator::mutate "

Description:

    how to mutate a single allele. this is usually the only function
    that need to be defined by the subclasses.

Usage:

    x.mutate(allele)

"; 

%feature("docstring") simuPOP::gsmMutator::__repr__ "

Description:

    simuPOP::gsmMutator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::ifElse "

Description:

    simuPOP::ifElse

"; 

%feature("docstring") simuPOP::ifElse::ifElse "

Description:

    conditional operator

Usage:

    ifElse(cond, *ifOp=NULL, *elseOp=NULL, output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    cond:           expression, will be treated as bool variable.
    ifOp:           if operator, be called when expr is true
    elseOp:         else operator, be called when expr is false

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::ifElse::applyWithScratch "

Description:

    simply output some info providing interface to apply operator
    before during or after mating.

Usage:

    x.applyWithScratch(pop, scratch, stage)

"; 

%feature("docstring") simuPOP::ifElse::applyDuringMating "

Description:

    give pop, offspring, pop and mom.

Usage:

    x.applyDuringMating(pop, offspring, *dad=NULL, *mom=NULL)

"; 

%feature("docstring") simuPOP::ifElse::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::ifElse::__repr__ "

Description:

    simuPOP::ifElse::__repr__

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

    Basic individual class.

Details:

    class individual with
    * genotypic information
    * shared genotypic structure info (through a  GenoStructure
    pointer)
    * flags about sex, affected status
    * an internal info field other individuals will be derived from
    this class, adding age info etc.Note  that
    * individual DOES NOT manage memory. It will use a pointer passed
    from class population. This causes A LOT of trouble and I have not
    evaluated how much benefic I get.
    * operator = uses shallow copy. This is required by sort algorithm
    since otherwise individuals are non-copiable. However, in
    population memory management, it is showtimes required that
    genotypic information within one subPop should go together. This
    is done by a shollow_copied flag for each individual and for all
    individuals. population might have to re-arrange individuals to
    solve this problem.
    * output of individual can be adjusted by setOutputDelimeter.
    Usage info: (for population classes developers)
    * for individuals are created, you are responsible to set its
    genotypic pointer and genotypic information. This is done by
    *  setSubPopID()  and  subPopID()  can be used for any temporary
    purpose.

"; 

%feature("docstring") simuPOP::individual::individual "

Description:

    default constructor,

Usage:

    individual()

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

%feature("docstring") simuPOP::individual::copyFrom "

Description:

    Deep copy! Important!

Usage:

    x.copyFrom(rhs)

"; 

%ignore simuPOP::individual::genoPtr() const ;

%ignore simuPOP::individual::infoPtr() const ;

%feature("docstring") simuPOP::individual::arrGenotype "

Description:

    return genotype as python Numeric.array object This is the whole
    genotype (all)

Usage:

    x.arrGenotype()

"; 

%feature("docstring") simuPOP::individual::arrGenotype "

Description:

    return genotype as python Numeric.array object This is the p'th
    copy of chromosomes

Usage:

    x.arrGenotype(p)

"; 

%feature("docstring") simuPOP::individual::arrGenotype "

Description:

    return genotype as python Numeric.array object This is the ch
    chromosome of the pth copy of chromosome

Usage:

    x.arrGenotype(p, ch)

"; 

%feature("docstring") simuPOP::individual::arrInfo "

Description:

    simuPOP::individual::arrInfo

Usage:

    x.arrInfo()

"; 

%feature("docstring") simuPOP::individual::allele "

Description:

    get allele from an index

Usage:

    x.allele(index)

Arguments:

    index:          index from the beginning of genotypic info

"; 

%feature("docstring") simuPOP::individual::allele "

Description:

    get allele from an index, on the pth set of chromosome

Usage:

    x.allele(index, p)

Arguments:

    index:          index from the begining of the p'th set of
                    chromosomes.
    p:              on p'th set of chromosomes, default to 0

"; 

%feature("docstring") simuPOP::individual::allele "

Description:

    simuPOP::individual::allele

Usage:

    x.allele(index, p, ch)

"; 

%feature("docstring") simuPOP::individual::alleleChar "

Description:

    simuPOP::individual::alleleChar

Usage:

    x.alleleChar(index)

"; 

%feature("docstring") simuPOP::individual::alleleChar "

Description:

    get allele from an index, on the pth set of chromosome

Usage:

    x.alleleChar(index, p)

Arguments:

    index:          index from the begining of the p'th set of
                    chromosomes.
    p:              on p'th set of chromosomes, p=0 by default

"; 

%feature("docstring") simuPOP::individual::setAllele "

Description:

    set allele from an index.

Usage:

    x.setAllele(allele, index)

Arguments:

    index:          index from the begining of genotype

"; 

%feature("docstring") simuPOP::individual::setAllele "

Description:

    simuPOP::individual::setAllele

Usage:

    x.setAllele(allele, index, p, ch)

"; 

%feature("docstring") simuPOP::individual::sex "

Description:

    sex?

Usage:

    x.sex()

"; 

%feature("docstring") simuPOP::individual::sexChar "

Description:

    return M or F for sex, for display purpose

Usage:

    x.sexChar()

"; 

%feature("docstring") simuPOP::individual::setSex "

Description:

    set sex

Usage:

    x.setSex(sex)

"; 

%feature("docstring") simuPOP::individual::affected "

Description:

    affected?

Usage:

    x.affected()

"; 

%feature("docstring") simuPOP::individual::unaffected "

Description:

    unaffected?

Usage:

    x.unaffected()

"; 

%feature("docstring") simuPOP::individual::affectedChar "

Description:

    return A or U for affected/Unaffected, for display purpose

Usage:

    x.affectedChar()

"; 

%feature("docstring") simuPOP::individual::setAffected "

Description:

    set affected status

Usage:

    x.setAffected(affected)

"; 

%feature("docstring") simuPOP::individual::subPopID "

Description:

    get subpop id

Usage:

    x.subPopID()

"; 

%feature("docstring") simuPOP::individual::setSubPopID "

Description:

    set subpop if

Usage:

    x.setSubPopID(id)

"; 

%feature("docstring") simuPOP::individual::info "

Description:

    get info

Usage:

    x.info(idx)

"; 

%feature("docstring") simuPOP::individual::setInfo "

Description:

    set info

Usage:

    x.setInfo(value, idx)

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

    simuPOP::individual::__cmp__

Usage:

    x.__cmp__(rhs)

"; 

%feature("docstring") simuPOP::individual::__repr__ "

Description:

    simuPOP::individual::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::individual::swap "

Description:

    swap individuals

Details:

    The default behavior is swapping all info, but not the position of
    genotypic info. If swapContent is false, pointer to genotypic info
    is swapped instead. This will lead to better performance for
    swapping but may affected performance of allele counting.

Usage:

    x.swap(ind, swapContent=True)

Arguments:

    ind:            individual to be swapped in
    swapContent:    swapContent or only the pointers. The guideline is
                    that if we swap individuals across subpopulation,
                    we should swap content. Otherwise, swap pointers.
                    (There is no order right now within subpopulation
                    so the later case is rare, at best.

"; 

%ignore simuPOP::individual::shallowCopied() const ;

%ignore simuPOP::individual::setShallowCopied(bool shallowCopied);

%ignore simuPOP::individual::display(ostream &out, int width=1, const vectori &chrom=vectori(), const vectori &loci=vectori());

%feature("docstring") simuPOP::individualIterator "

Description:

    this class implements a Python itertor class that can be used to
    iterate through individuals in a population. an instance of this
    class is returned by  population::individuals() and
    population::individuals(subPop)

"; 

%feature("docstring") simuPOP::individualIterator::individualIterator "

Description:

    simuPOP::individualIterator::individualIterator

Usage:

    individualIterator(*pop, s, e)

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

    inherite tag from parents. If both parents have tags, use fathers.

"; 

%feature("docstring") simuPOP::inheritTagger::inheritTagger "

Description:

    constructor. default to be always active.

Usage:

    inheritTagger(mode=TAG_Paternal, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[\"paternal_tag\",
      \"maternal_tag\"])

"; 

%feature("docstring") simuPOP::inheritTagger::~inheritTagger "

Description:

    simuPOP::inheritTagger::~inheritTagger

Usage:

    x.~inheritTagger()

"; 

%feature("docstring") simuPOP::inheritTagger::__repr__ "

Description:

    simuPOP::inheritTagger::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::inheritTagger::applyDuringMating "

Description:

    give pop, offspring, pop and mom.

Usage:

    x.applyDuringMating(pop, offspring, *dad=NULL, *mom=NULL)

"; 

%feature("docstring") simuPOP::inheritTagger::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initByFreq "

Description:

    initialize genotype by allele frequency and sex by male frequency

"; 

%feature("docstring") simuPOP::initByFreq::initByFreq "

Description:

    randomly assign alleles according to allele frequency

Details:

    This operator randomly assign alleles according to given allele
    frequency. Allele frequencies can differ by subpop. Sex is also
    assigned randomly.

Usage:

    initByFreq(alleleFreq=[], identicalInds=False, subPop=[],
      indRange=intMatrix, atLoci=[], atPloidy=-1, maleFreq=0.5,
      sex=[], stage=PreMating, begin=0, end=1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    alleleFreq:     an array of allele frequencies. Must add up to 1;
                    or a matrix of allele frequencies, each row
                    corresponse a subpopulation.
    subPop:         an array of applicable subpopulations. default to
                    all
    indRange:       a [begin, end] pair of range of individuals; or an
                    array of [begin, end] pairs.
    identicalInds:  whether or not make individual genotype identical
                    in all subpopulation. If true, this operator will
                    randomly generate genotype for an individual and
                    spread it to the whole subpopulation.
    atLoci:         a vector of loci indices. If empty, apply to all
                    loci
    atPloidy:       initialize which copy of chromosomes. Default to
                    all.
    maleFreq:       male frequency. Default to 0.5.
    sex:            an arry of sex [Male, Female, Male]... for
                    individuals. The length of sex will not be
                    checked. If length of sex is shorter than number
                    of individuals, sex will be reused from the
                    beginning.
    stages:         is set to PreMating. Other parameters please see
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::initByFreq::~initByFreq "

Description:

    simuPOP::initByFreq::~initByFreq

Usage:

    x.~initByFreq()

"; 

%feature("docstring") simuPOP::initByFreq::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initByFreq::__repr__ "

Description:

    simuPOP::initByFreq::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::initByFreq::apply "

Details:

    initialize m_ranges

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::initByValue "

Description:

    initialize genotype by value and then copy to all individuals

"; 

%feature("docstring") simuPOP::initByValue::initByValue "

Description:

    simuPOP::initByValue::initByValue

Usage:

    initByValue(value=intMatrix, atLoci=[], atPloidy=-1, subPop=[],
      indRange=intMatrix, proportions=[], maleFreq=0.5, sex=[],
      stage=PreMating, begin=0, end=1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::initByValue::~initByValue "

Description:

    simuPOP::initByValue::~initByValue

Usage:

    x.~initByValue()

"; 

%feature("docstring") simuPOP::initByValue::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initByValue::__repr__ "

Description:

    simuPOP::initByValue::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::initByValue::apply "

Details:

    fixme: check length of src?atLoci is in effect

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::initializer "

Description:

    initialize alleles at the start of generation.

Details:

    Bo Peng

"; 

%feature("docstring") simuPOP::initializer::initializer "

Description:

    constructor. default to be always active.

Usage:

    initializer(subPop=[], indRange=intMatrix, atLoci=[],
      atPloidy=-1, maleFreq=0.5, sex=[], stage=PreMating, begin=0,
      end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::initializer::~initializer "

Description:

    destructor

Usage:

    x.~initializer()

"; 

%feature("docstring") simuPOP::initializer::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::initializer::__repr__ "

Description:

    simuPOP::initializer::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::initializer::setRanges "

Description:

    simuPOP::initializer::setRanges

Usage:

    x.setRanges(pop)

"; 

%feature("docstring") simuPOP::initializer::initSexIter "

Description:

    simuPOP::initializer::initSexIter

Usage:

    x.initSexIter()

"; 

%feature("docstring") simuPOP::initializer::nextSex "

Description:

    simuPOP::initializer::nextSex

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

Description:

    K-Allele Model mutator.

Details:

    Under this model, there are K (here refers as maxAllele) possible
    allele states, and any allele has a constant probability
    (rate/(K-1)) of mutation towards any of the K-1 allelic states.

Note:

    the theoretical mutation rate is rates/(K-1) towards any of the
    K-1 allelic states. So rates is actually the probability to
    mutate!

"; 

%feature("docstring") simuPOP::kamMutator::kamMutator "

Description:

    K-Allele Model mutator.

Usage:

    kamMutator(rate=[], atLoci=[], maxAllele=0, output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    rate:           mutation rate. It is 'probability to mutate'. The
                    actual mutation rate to any of the other K-1
                    allelic states are rates/(K-1)!
    atLoci:         and other parameters: refer to help(mutator),
                    help(baseOperator.__init__)
    maxAllele:      maxAllele that can be mutated to. For binary
                    libraries allelic states will be [0, maxAllele].
                    For others, they are [1, maxAllele]

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::kamMutator::__repr__ "

Description:

    simuPOP::kamMutator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::largePedigreeSample "

Description:

    simuPOP::largePedigreeSample

"; 

%feature("docstring") simuPOP::largePedigreeSample::largePedigreeSample "

Description:

    draw cases and controls

Usage:

    largePedigreeSample(size=[], minTotalSize=0, maxOffspring=5,
      minPedSize=5, minAffected=0, countOnly=False, name=\"sample\",
      nameExpr=\"\", times=1, saveAs=\"\", saveAsExpr=\"\", format=\"auto\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[\"father_idx\", \"mother_idx\"])

Arguments:

    size:           number of affected sibpairs to be sampled. Can be
                    a number or an array. If a number is given, it is
                    the total number of sibpairs, ignoring population
                    structure. Otherwise, given number of sibpairs are
                    sampled from subpopulations. If size is
                    unspecified, this operator will return all
                    affected sibpairs.
    countOnly:      set variables about number of affected sibpairs,
                    do not actually draw the sample
    name:           variable name of the sampled population (will be
                    put in pop local namespace)
    nameExpr:       expression version of name. If both name and
                    nameExpr is empty, do not store pop.
    times:          how many times to run the sample process? This is
                    usually one, but we may want to take several
                    random samples.
    saveAs:         filename to save the population.
    saveAsExpr:     expression for save filename
    format:         to save sample(s)
    stage:          and other parameters please see
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::largePedigreeSample::~largePedigreeSample "

Description:

    destructor

Usage:

    x.~largePedigreeSample()

"; 

%feature("docstring") simuPOP::largePedigreeSample::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::largePedigreeSample::prepareSample "

Description:

    simuPOP::largePedigreeSample::prepareSample

Usage:

    x.prepareSample(pop)

"; 

%feature("docstring") simuPOP::largePedigreeSample::drawsample "

Details:

    collect all families

Usage:

    x.drawsample(pop)

"; 

%feature("docstring") simuPOP::largePedigreeSample::__repr__ "

Description:

    simuPOP::largePedigreeSample::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::maPenetrance "

Description:

    penetrance according to genotype at one locus

Details:

    multiple allele selector. This selector group alleles to disease
    and wild type and return penetrance to AA,Aa,aa. (A is wildtype).

"; 

%feature("docstring") simuPOP::maPenetrance::maPenetrance "

Description:

    create a multiple allele selector (penetrance according to
    diseased or wildtype alleles)

Usage:

    maPenetrance(loci, penet, wildtype, ancestralGen=-1,
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    locus:          the locus index. The genotype of this locus will
                    be axamed.
    loci:           the loci index.
    penetrance:     an array of penetrance of AA,Aa,aa. A is the wild
                    type group. In the case of multiple loci, fitness
                    should be in the order of BB Bb bb AA 1 2 3 Aa 4 5
                    6 aa 7 8 9
    wildtype:       an array of alleles in the wildtype group.
                    Anything else is disease allele., default = [0]
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::maPenetrance::~maPenetrance "

Description:

    simuPOP::maPenetrance::~maPenetrance

Usage:

    x.~maPenetrance()

"; 

%feature("docstring") simuPOP::maPenetrance::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::maPenetrance::penet "

Description:

    currently assuming diploid

Details:

    get genotype of ind

Usage:

    x.penet(*ind)

"; 

%feature("docstring") simuPOP::maPenetrance::__repr__ "

Description:

    simuPOP::maPenetrance::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mapPenetrance "

Description:

    penetrance according to genotype at one locus

Details:

    map selector. Assign penetrance value according to a given
    dictionary.

"; 

%feature("docstring") simuPOP::mapPenetrance::mapPenetrance "

Description:

    create a map penetrance function (penetrance according to genotype
    at one locus

Usage:

    mapPenetrance(loci, penet, phase=False, ancestralGen=-1,
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    locus:          the locus index. The genotype of this locus will
                    be axamed.
    loci:           the loci index. The genotype of this locus will be
                    axamed.
    penetrance:     a dictionary of penetrance. The genotype must be
                    in the form of 'a-b' for single locus.
    phase:          if true, a/b and b/a will have different
                    penetrance value. Default to false.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::mapPenetrance::~mapPenetrance "

Description:

    simuPOP::mapPenetrance::~mapPenetrance

Usage:

    x.~mapPenetrance()

"; 

%feature("docstring") simuPOP::mapPenetrance::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mapPenetrance::penet "

Description:

    currently assuming diploid

Details:

    get genotype of ind

Usage:

    x.penet(*ind)

"; 

%feature("docstring") simuPOP::mapPenetrance::__repr__ "

Description:

    simuPOP::mapPenetrance::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mapQuanTrait "

Description:

    quantitative trait according to genotype at one locus

Details:

    map selector. Assign qtrait value according to a given dictionary.

"; 

%feature("docstring") simuPOP::mapQuanTrait::mapQuanTrait "

Description:

    create a map selector (quantitative trait according to genotype at
    one locus

Usage:

    mapQuanTrait(loci, qtrait, sigma=0, phase=False,
      ancestralGen=-1, stage=PostMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[\"qtrait\"])

Arguments:

    locus:          the locus index. The genotype of this locus will
                    be axamed.
    loci:           the loci.
    qtrait:         a dictionary of qtrait. The genotype must be in
                    the form of 'a-b'. This is the mean of
                    quantitative trait. The actual trait value will be
                    N(mean, sigma^2) For multiple loci, the form is
                    'a-b|c-d|e-f' etc.
    sigma:          standard deviation of the environmental facotr
                    N(0,sigma^2).
    phase:          if true, a/b and b/a will have different qtrait
                    value. Default to false.
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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mapQuanTrait::qtrait "

Description:

    currently assuming diploid

Details:

    get genotype of ind

Usage:

    x.qtrait(*ind)

"; 

%feature("docstring") simuPOP::mapQuanTrait::__repr__ "

Description:

    simuPOP::mapQuanTrait::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mapSelector "

Description:

    selection according to genotype at one locus

Details:

    map selector. Assign fitness value according to a given
    dictionary.

"; 

%feature("docstring") simuPOP::mapSelector::mapSelector "

Description:

    create a map selector (selection according to genotype at one
    locus

Usage:

    mapSelector(loci, fitness, phase=False, stage=PreMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[\"fitness\"])

Arguments:

    locus:          the locus index. The genotype of this locus will
                    be axamed.
    loci:           the loci index. The genotype of this locus will be
                    axamed.
    fitness:        a dictionary of fitness. The genotype must be in
                    the form of 'a-b' for single locus, and
                    'a-b|c-d|e-f' for multi-locus..
    phase:          if true, a/b and b/a will have different fitness
                    value. Default to false.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::mapSelector::~mapSelector "

Description:

    simuPOP::mapSelector::~mapSelector

Usage:

    x.~mapSelector()

"; 

%feature("docstring") simuPOP::mapSelector::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mapSelector::indFitness "

Description:

    currently assuming diploid

Details:

    get genotype of ind

Usage:

    x.indFitness(*ind, gen)

"; 

%feature("docstring") simuPOP::mapSelector::__repr__ "

Description:

    simuPOP::mapSelector::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::maQuanTrait "

Description:

    quantitative trait according to genotype at one locus

Details:

    multiple allele selector. This selector group alleles to disease
    and wild type and return qtrait to AA,Aa,aa. (A is wildtype).

"; 

%feature("docstring") simuPOP::maQuanTrait::maQuanTrait "

Description:

    create a multiple allele selector (quantitative trait according to
    diseased or wildtype alleles)

Usage:

    maQuanTrait(loci, qtrait, wildtype, sigma=[], ancestralGen=-1,
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[\"qtrait\"])

Arguments:

    locus:          the locus index. The genotype of this locus will
                    be axamed.
    qtrait:         an array of qtrait of AA,Aa,aa. A is the wild type
                    group.
    sigma:          an array of standard deviation for each of the
                    trait genotype (AA, Aa, aa)
    wildtype:       an array of alleles in the wildtype group.
                    Anything else is disease allele. default = [0]
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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::maQuanTrait::qtrait "

Description:

    currently assuming diploid

Details:

    get genotype of ind

Usage:

    x.qtrait(*ind)

"; 

%feature("docstring") simuPOP::maQuanTrait::__repr__ "

Description:

    simuPOP::maQuanTrait::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::maSelector "

Description:

    selection according to genotype at one locus

Details:

    multiple allele selector. This selector group alleles to disease
    and wild type and return fitness to AA,Aa,aa. (A is wildtype).

"; 

%feature("docstring") simuPOP::maSelector::maSelector "

Description:

    create a multiple allele selector (selection according to diseased
    or wildtype alleles) Note that  maSelector only work for diploid
    population now.

Usage:

    maSelector(loci, fitness, wildtype, stage=PreMating, begin=0,
      end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[\"fitness\"])

Arguments:

    locus:          the locus index. The genotype of this locus will
                    be axamed.
    loci:           the loci index.
    fitness:        For the single locus case, fitness is an array of
                    fitness of AA,Aa,aa. A is the wild type group. In
                    the case of multiple loci, fitness should be in
                    the order of BB Bb bb AA 1 2 3 Aa 4 5 6 aa 7 8 9
                    The length for such table is 3^(#loci).
    wildtype:       an array of alleles in the wildtype group.
                    Anything else is disease allele. default = [0]
                    NOTE that wildtype at all loci are the same.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::maSelector::~maSelector "

Description:

    simuPOP::maSelector::~maSelector

Usage:

    x.~maSelector()

"; 

%feature("docstring") simuPOP::maSelector::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::maSelector::indFitness "

Description:

    currently assuming diploid

Details:

    get genotype of ind

Usage:

    x.indFitness(*ind, gen)

"; 

%feature("docstring") simuPOP::maSelector::__repr__ "

Description:

    simuPOP::maSelector::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mating "

Details:

    The mating classes describe various mating scheme --- a required
    parameter of simulator.

"; 

%feature("docstring") simuPOP::mating::isCompatible "

Description:

    check if the mating type is compatible with population structure

Details:

    possible things to check:
    * need certain types of individual (age, sex etc)
    * need resizeable population...

Usage:

    x.isCompatible(pop)

"; 

%feature("docstring") simuPOP::mating::mating "

Description:

    constructor

Usage:

    mating(numOffspring=1.0, *numOffspringFunc=NULL,
      maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
      newSubPopSizeExpr=\"\", *newSubPopSizeFunc=NULL)

"; 

%feature("docstring") simuPOP::mating::mating "

Description:

    simuPOP::mating::mating

Usage:

    mating(rhs)

"; 

%feature("docstring") simuPOP::mating::~mating "

Description:

    destructor

Usage:

    x.~mating()

"; 

%feature("docstring") simuPOP::mating::clone "

Description:

    clone() const. Generate a copy of itself and return a pointer

Details:

    This function is important because Python automatically release an
    object after it is used.For example: will fail since  mating() is
    released after the first line being executed.With the help of
    clone() const, the C++ implementation can avoid this problem by

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mating::__repr__ "

Description:

    return name of the mating type used primarily in logging.

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mating::submitScratch "

Description:

    simuPOP::mating::submitScratch

Usage:

    x.submitScratch(pop, scratch)

"; 

%feature("docstring") simuPOP::mating::mate "

Description:

    mate: This is not supposed to be called for base mating class.

Usage:

    x.mate(pop, scratch, ops, submit)

Arguments:

    pop:            population
    scratch:        scratch population
    ops:            during mating operators

"; 

%feature("docstring") simuPOP::mating::fixedFamilySize "

Description:

    simuPOP::mating::fixedFamilySize

Usage:

    x.fixedFamilySize()

"; 

%feature("docstring") simuPOP::mating::numOffspring "

Description:

    simuPOP::mating::numOffspring

Usage:

    x.numOffspring(gen)

"; 

%feature("docstring") simuPOP::mating::resetNumOffspring "

Description:

    simuPOP::mating::resetNumOffspring

Usage:

    x.resetNumOffspring()

"; 

%feature("docstring") simuPOP::mating::prepareScratchPop "

Description:

    dealing with pop/subPop size change, copy of structure etc.

Usage:

    x.prepareScratchPop(pop, scratch)

"; 

%feature("docstring") simuPOP::mergeSubPops "

Description:

    merge subpopulations

"; 

%feature("docstring") simuPOP::mergeSubPops::mergeSubPops "

Description:

    merge subpopulations

Usage:

    mergeSubPops(subPops=[], removeEmptySubPops=False,
      stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    subPops:        subpops to be merged, default to all subpops.

"; 

%feature("docstring") simuPOP::mergeSubPops::~mergeSubPops "

Description:

    destructor

Usage:

    x.~mergeSubPops()

"; 

%feature("docstring") simuPOP::mergeSubPops::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mergeSubPops::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::mergeSubPops::__repr__ "

Description:

    simuPOP::mergeSubPops::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::migrator "

Details:

    Bo Peng

"; 

%feature("docstring") simuPOP::migrator::migrator "

Description:

    create a migrator

Usage:

    migrator(rate, mode=MigrByProbability, fromSubPop=[],
      toSubPop=[], stage=PreMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    rate:           migration rate, proportion or count. Determined by
                    parameter mode. rate should be a m by n matrix. If
                    a number is given, the migration rate will be
                    r*ones(m,n).
    mode:           one of MigrByProbability (default),
                    MigrByProportion or MigrByCounts
    fromSubPop:     an array of 'from' subpops, default to all
                    subpopulations. If a single subpop is specified,
                    [] can be ignored. I.e., [a] is equvalent to a.
    toSubPop:       an array of 'to' subpops, default to all
                    subpopulations. If a single subpop is specified,
                    [] can be ignored.
    stage:          is default to PreMating. For details about other
                    parameters, please refer to
                    help(baseOperator.__init__) rate is a matrix with
                    dimensions determined by fromSubPop and toSubPop.
                    By default, rate is a matrix with element (i,j)
                    being the migration rate, probability or count
                    from subpop i to subpop j. If fromSubPop and/or
                    toSubPop are given, migration only happen between
                    these subpopulations. An extreme case is 'point
                    migration'rate=[[r]], fromSubPop=a,
                    toSubPop=bwhich migrate from subpop a to b with
                    given rate r.

"; 

%feature("docstring") simuPOP::migrator::~migrator "

Description:

    destructor

Usage:

    x.~migrator()

"; 

%feature("docstring") simuPOP::migrator::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::migrator::rate "

Description:

    return rate

Usage:

    x.rate()

"; 

%feature("docstring") simuPOP::migrator::setRates "

Description:

    set migration rate

Details:

    format 0-0 0-1 0-2, 1-0 1-1 1-2, 2-0, 2-1, 2-2. for mode 1 or 2,
    00,11,22 will be set automatically. regardless of input.

Usage:

    x.setRates(rate, mode)

"; 

%feature("docstring") simuPOP::migrator::apply "

Details:

    2nd, or 3rd methodcreate a vector and assign indices, then random
    shuffle and assign infofor all subPop.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::migrator::__repr__ "

Description:

    simuPOP::migrator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mlPenetrance "

Description:

    penetrance according to genotype at multiple loci multiplicative
    model

Details:

    multiple loci selector. This selector takes several selectors and
    multiply their penetrance values... e.g. mlmpenetrance(
    [mappenetrance(...), mapenetrance(...) ])

"; 

%feature("docstring") simuPOP::mlPenetrance::mlPenetrance "

Description:

    multiple loci selector using a multiplicative model.

Usage:

    mlPenetrance(peneOps, mode=PEN_Multiplicative, ancestralGen=-1,
      stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    selectors:      a list of selectors.
    mode:           one of PEN_Multiplicative, PEN_Additive,
                    PEN_Heterogeneity

"; 

%feature("docstring") simuPOP::mlPenetrance::~mlPenetrance "

Description:

    simuPOP::mlPenetrance::~mlPenetrance

Usage:

    x.~mlPenetrance()

"; 

%feature("docstring") simuPOP::mlPenetrance::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mlPenetrance::penet "

Description:

    currently assuming diploid

Usage:

    x.penet(*ind)

"; 

%feature("docstring") simuPOP::mlPenetrance::__repr__ "

Description:

    simuPOP::mlPenetrance::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mlQuanTrait "

Description:

    quantitative trait according to genotype at multiple loci
    multiplicative model

Details:

    multiple loci selector. This selector takes several selectors and
    multiply their qtrait values... e.g. mlmquanTrait(
    [mapquanTrait(...), maquanTrait(...) ])

"; 

%feature("docstring") simuPOP::mlQuanTrait::mlQuanTrait "

Description:

    multiple loci selector using a multiplicative model.

Usage:

    mlQuanTrait(qtraits, mode=QT_Multiplicative, sigma=0,
      ancestralGen=-1, stage=PostMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[\"qtrait\"])

Arguments:

    qtraits:        a list of qtraits.

"; 

%feature("docstring") simuPOP::mlQuanTrait::~mlQuanTrait "

Description:

    simuPOP::mlQuanTrait::~mlQuanTrait

Usage:

    x.~mlQuanTrait()

"; 

%feature("docstring") simuPOP::mlQuanTrait::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mlQuanTrait::qtrait "

Description:

    currently assuming diploid

Usage:

    x.qtrait(*ind)

"; 

%feature("docstring") simuPOP::mlQuanTrait::__repr__ "

Description:

    simuPOP::mlQuanTrait::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mlSelector "

Description:

    selection according to genotype at multiple loci multiplicative
    model

Details:

    multiple loci selector. This selector takes several selectors and
    multiply their fitness values... e.g. mlmselector(
    [mapselector(...), maselector(...) ])

"; 

%feature("docstring") simuPOP::mlSelector::mlSelector "

Description:

    multiple loci selector using a multiplicative model.

Usage:

    mlSelector(selectors, mode=SEL_Multiplicative, stage=PreMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[\"fitness\"])

Arguments:

    selectors:      a list of selectors.

"; 

%feature("docstring") simuPOP::mlSelector::~mlSelector "

Description:

    simuPOP::mlSelector::~mlSelector

Usage:

    x.~mlSelector()

"; 

%feature("docstring") simuPOP::mlSelector::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mlSelector::indFitness "

Description:

    currently assuming diploid

Details:

    fixme

Usage:

    x.indFitness(*ind, gen)

"; 

%feature("docstring") simuPOP::mlSelector::__repr__ "

Description:

    simuPOP::mlSelector::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::mutator "

Description:

    mutator class.

Details:

    Do not use this class directly. It just provide interface for real
    mutators.Every mutator can specify rate (equal rate) or rates
    (different rate for different loci) and a vector of applicable
    loci (default to all but should have the same length with rates if
    rates have length greater than one).max allele can be specified as
    well but more parameter, if needed, should be implemented by
    individual mutator classes.Number of possible allelic states: Most
    theoretical studies assume an infinite number of allelic states to
    avoid any homoplasy. If it facilitates analysis, this is however
    extremely unrealistic.Bo Peng

"; 

%feature("docstring") simuPOP::mutator::mutator "

Description:

    create a mutator All mutators have the following common
    parameters. However, the actual meaning of these parameters may
    vary according to different model. Check the manual for details!!!
    (help(kamMutator) for example.)

Usage:

    mutator(rate=[], atLoci=[], maxAllele=0, output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    rate:           single rate for all applicable loci (atLoci). Will
                    be ignored if rates is specified; or it can be an
                    array of rates, the same length as atLoci.
    atLoci:         a vector of loci index. Can be ignored only when
                    single rate is specified. Default to all loci.
    maxAllele:      max allowable allele. Interpreted by each sub
                    mutaor class. Default to pop.maxAllele().

"; 

%feature("docstring") simuPOP::mutator::~mutator "

Description:

    destructor

Usage:

    x.~mutator()

"; 

%feature("docstring") simuPOP::mutator::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mutator::rate "

Description:

    return mutation rate

Usage:

    x.rate()

"; 

%feature("docstring") simuPOP::mutator::setRate "

Description:

    set an array of rates

Usage:

    x.setRate(rate, atLoci=[])

"; 

%feature("docstring") simuPOP::mutator::maxAllele "

Description:

    return max allowable allele number

Usage:

    x.maxAllele()

"; 

%feature("docstring") simuPOP::mutator::setMaxAllele "

Description:

    simuPOP::mutator::setMaxAllele

Usage:

    x.setMaxAllele(maxAllele)

"; 

%feature("docstring") simuPOP::mutator::mutationCount "

Description:

    return mutation count

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

    how to mutate a single allele. this is usually the only function
    that need to be defined by the subclasses.

Usage:

    x.mutate(allele)

"; 

%feature("docstring") simuPOP::mutator::apply "

Description:

    apply!

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::noMating "

Details:

    No mating. No subpopulation change. During mating operator will be
    applied, but the return values are not checked.

"; 

%feature("docstring") simuPOP::noMating::noMating "

Description:

    constructor, no new subPopsize parameter

Usage:

    noMating()

"; 

%feature("docstring") simuPOP::noMating::~noMating "

Description:

    destructor

Usage:

    x.~noMating()

"; 

%feature("docstring") simuPOP::noMating::clone "

Description:

    clone() const. The same as  mating::clone() const.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::noMating::__repr__ "

Description:

    return name of the mating type

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::noMating::submitScratch "

Description:

    simuPOP::noMating::submitScratch

Usage:

    x.submitScratch(pop, scratch)

"; 

%feature("docstring") simuPOP::noMating::mate "

Description:

    do the mating. --- no mating :-)

Details:

    All individuals will be passed to during mating operators but no
    one will die (ignore during mating failing signal).

Usage:

    x.mate(pop, scratch, ops, submit)

"; 

%feature("docstring") simuPOP::noneOp "

Description:

    simuPOP::noneOp

"; 

%feature("docstring") simuPOP::noneOp::noneOp "

Description:

    do nothing

Usage:

    noneOp(output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
      end=0, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::noneOp::~noneOp "

Description:

    destructor

Usage:

    x.~noneOp()

"; 

%feature("docstring") simuPOP::noneOp::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::noneOp::applyWithScratch "

Description:

    simply output some info providing interface to apply operator
    before during or after mating.

Usage:

    x.applyWithScratch(pop, scratch, stage)

"; 

%feature("docstring") simuPOP::noneOp::applyDuringMating "

Description:

    give pop, offspring, pop and mom.

Usage:

    x.applyDuringMating(pop, offspring, *dad=NULL, *mom=NULL)

"; 

%feature("docstring") simuPOP::noneOp::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::noneOp::__repr__ "

Description:

    simuPOP::noneOp::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::nuclearFamilySample "

Description:

    simuPOP::nuclearFamilySample

"; 

%feature("docstring") simuPOP::nuclearFamilySample::nuclearFamilySample "

Description:

    draw nuclear family

Usage:

    nuclearFamilySample(size=[], minTotalSize=0, maxOffspring=5,
      minPedSize=5, minAffected=0, countOnly=False, name=\"sample\",
      nameExpr=\"\", times=1, saveAs=\"\", saveAsExpr=\"\", format=\"auto\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[\"father_idx\", \"mother_idx\"])

Arguments:

    name:           variable name of the sampled population (will be
                    put in pop local namespace)
    nameExpr:       expression version of name. If both name and
                    nameExpr is empty, do not store pop.
    times:          how many times to run the sample process? This is
                    usually one, but we may want to take several
                    random samples.
    saveAs:         filename to save the population.
    saveAsExpr:     expression for save filename
    format:         to save sample(s)
    stage:          and other parameters please see
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::nuclearFamilySample::~nuclearFamilySample "

Description:

    destructor

Usage:

    x.~nuclearFamilySample()

"; 

%feature("docstring") simuPOP::nuclearFamilySample::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::nuclearFamilySample::prepareSample "

Description:

    simuPOP::nuclearFamilySample::prepareSample

Usage:

    x.prepareSample(pop)

"; 

%feature("docstring") simuPOP::nuclearFamilySample::drawsample "

Details:

    collect all families

Usage:

    x.drawsample(pop)

"; 

%feature("docstring") simuPOP::nuclearFamilySample::__repr__ "

Description:

    simuPOP::nuclearFamilySample::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::NullStreamBuf "

Description:

    create a null stream buf that discard everything

"; 

%feature("docstring") simuPOP::NullStreamBuf::NullStreamBuf "

Description:

    simuPOP::NullStreamBuf::NullStreamBuf

Usage:

    NullStreamBuf()

"; 

%feature("docstring") simuPOP::offspringGenerator "

Description:

    the default method to generate offspring from parents This part is
    separated from the mating schemes, because mating schemes usually
    only differ by the way parents are choosing. input: parents,
    output: offsprings

"; 

%feature("docstring") simuPOP::offspringGenerator::offspringGenerator "

Description:

    constructor, save information from pop and ops to speed up the
    calls to generateOffspring

Usage:

    offspringGenerator(pop, ops)

"; 

%feature("docstring") simuPOP::offspringGenerator::generateOffspring "

Description:

    generate numOff offspring, or until reach offEnd this is because
    offBegin+numOff may go beyond subpopulation boundary. return the
    ending iterator

Details:

    apply all during mating operators

Usage:

    x.generateOffspring(pop, *dad, *mom, numOff, offBegin)

"; 

%feature("docstring") simuPOP::offspringGenerator::copyOffspring "

Description:

    generate numOff offspring, or until reach offEnd this is because
    offBegin+numOff may go beyond subpopulation boundary. return the
    ending iterator

Details:

    use deep copy!!!!!!!

Usage:

    x.copyOffspring(pop, *par, numOff, offBegin)

"; 

%feature("docstring") simuPOP::Operator "

Description:

    base class of all classes that manipulate populations.

Details:

    Operators are object that act on populations. They can be applied
    to populations directly using  apply() member function, but most
    of the time they are managed and applied by a simulator.Operators
    can be applied at different stage(s) of a life cycle. More
    specifically, at pre-, duing- or post mating stage(s). Note that
    it is possible for an operator to apply multiple times in a life
    cycle. For example, an save-to-file operator might be applied
    before and after mating to trace parental information.Operators do
    not have to be applied at all generations. You can specify
    starting genertion, ending generation, gaps between applicable
    generations, or even specific generations to apply. For example,
    you might want to start applying migrations after certain heat-up
    generation; or you want to calculate every 10
    generations.Operators can have outputs. Output can be standard
    output (terminal) or a file, which can be constant, or change with
    generation or replicate. Different operators can append to the
    same file to form table-like outputs.filename can have the
    following format:1. 'filename' this file will be closed after each
    use. I.e., if several operators output to the same file, only the
    last one will succeed.2. '>filename' the same as 'filaname'3.
    '>>filename' The file will be created at the beginning of
    evolution ( simulator::evolve) and close at the end. Several
    operators can output to this file to form a table.4. '>>>filename'
    The same as '>>filename' except that the file will not be cleared
    at the beginning of evolution if it is not empty.5. '>' out put to
    standard output.6. '' supress output.Most operators are applied to
    every replicate of a simulator during evolution. However, you can
    apply opertors to one or a group of replicates only. For example,
    you can initialize different replicates with different initial
    values and then start evolution. c.f.  simulator::setGroup .Please
    refer to help(baseOperator) and help(baseOperator.__init__) for
    detailed information about member functions and parameters.Bo Peng

"; 

%feature("docstring") simuPOP::Operator::Operator "

Description:

    create an operator (this function is not supposed to be called
    directly)

Usage:

    Operator(output, outputExpr, stage, begin, end, step, at, rep,
      grp, infoFields)

Arguments:

    output:         a string of output filename. Different operators
                    will have different default output (most commonly
                    '>' or '')
    outputExpr:     an expression that determines output filename
                    dynamically.
    begin:          start generation. default to 1. negative number is
                    interpreted as endGeneration + begin
    end:            stop applying after this generation. negative
                    number is allowed
    step:           number of generations between active generations.
                    default to 1
    at:             an array of active generations. If given, stage,
                    begin, end, step will be ignored.
    rep:            applicable replicate. It can be replicate number 0
                    ~ (number of replicate -1), REP_ALL (all
                    replicates) or REP_LAST (only to the last
                    replicate). Usually default to REP_ALL.
    grp:            applicable group, default to GRP_ALL. A group
                    number for each replicate is set by
                    simulator.__init__ or  simulator::setGroup(). grp,
                    if not GRP_ALL, will be compared to the group
                    number of this replicate before applying. DEVONLY{
                    DO NOT SET DEFAULT PARAMETE. This will force all
                    derived classes to pay attention to parameter
                    number. }

"; 

%feature("docstring") simuPOP::Operator::~Operator "

Description:

    destroy an operator

Usage:

    x.~Operator()

"; 

%feature("docstring") simuPOP::Operator::clone "

Description:

    this function is very important

Details:

    use of parameter start, end, every, at, group, rep

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::Operator::isActive "

Description:

    judge if this operator is active

Details:

    Judge if this operator is active under the conditions like current
    replicate, group number of current replicate, current generation.
    ending generation etc.

Usage:

    x.isActive(rep, numRep, gen, end, grp, repOnly=False)

Note:

    will be called by simulator before applying.

"; 

%feature("docstring") simuPOP::Operator::applicableGroup "

Description:

    return applicable group

Usage:

    x.applicableGroup()

"; 

%feature("docstring") simuPOP::Operator::setApplicableGroup "

Description:

    set applicable group.

Details:

    GRP_ALL is the default value (applicable to all groups. )
    Otherwise, the operator is applicable to ONE group of replicates.
    groups can be set in  simulator::setGroup()

Usage:

    x.setApplicableGroup(grp=GRP_ALL)

"; 

%feature("docstring") simuPOP::Operator::applicableReplicate "

Description:

    return applicable replicate

Usage:

    x.applicableReplicate()

"; 

%feature("docstring") simuPOP::Operator::setApplicableReplicate "

Description:

    set applicable replicate.

Usage:

    x.setApplicableReplicate(rep)

"; 

%feature("docstring") simuPOP::Operator::setActiveGenerations "

Description:

    set applicable generation parrameters: stage, begin, end, step and
    at

Details:

    set certain m_flags to speed up using this machanism. pre, during,
    post-mating methods

Usage:

    x.setActiveGenerations(begin=0, end=-1, step=1, at=[])

"; 

%feature("docstring") simuPOP::Operator::setApplicableStage "

Description:

    set m_stage settings. This is usually not usable since the m_stage
    setting are set by default for each  Operator.

Usage:

    x.setApplicableStage(stage)

"; 

%feature("docstring") simuPOP::Operator::canApplyPreMating "

Description:

    Can this operator be applied pre-mating?

Usage:

    x.canApplyPreMating()

"; 

%feature("docstring") simuPOP::Operator::canApplyDuringMating "

Description:

    Can this operator be applied uring mating?

Usage:

    x.canApplyDuringMating()

"; 

%feature("docstring") simuPOP::Operator::canApplyPostMating "

Description:

    Can this operator be applied post-mating?

Usage:

    x.canApplyPostMating()

"; 

%feature("docstring") simuPOP::Operator::canApplyPreOrPostMating "

Description:

    can be applied pre or post mating.

Usage:

    x.canApplyPreOrPostMating()

"; 

%feature("docstring") simuPOP::Operator::isCompatible "

Description:

    simuPOP::Operator::isCompatible

Usage:

    x.isCompatible(pop)

"; 

%feature("docstring") simuPOP::Operator::haploidOnly "

Description:

    simuPOP::Operator::haploidOnly

Usage:

    x.haploidOnly()

"; 

%feature("docstring") simuPOP::Operator::diploidOnly "

Description:

    simuPOP::Operator::diploidOnly

Usage:

    x.diploidOnly()

"; 

%feature("docstring") simuPOP::Operator::MPIReady "

Description:

    simuPOP::Operator::MPIReady

Usage:

    x.MPIReady()

"; 

%feature("docstring") simuPOP::Operator::setHaploidOnly "

Description:

    simuPOP::Operator::setHaploidOnly

Usage:

    x.setHaploidOnly()

"; 

%feature("docstring") simuPOP::Operator::setDiploidOnly "

Description:

    simuPOP::Operator::setDiploidOnly

Usage:

    x.setDiploidOnly()

"; 

%feature("docstring") simuPOP::Operator::setMPIReady "

Description:

    simuPOP::Operator::setMPIReady

Usage:

    x.setMPIReady()

"; 

%feature("docstring") simuPOP::Operator::infoSize "

Description:

    get the number of information fields for this operator

Usage:

    x.infoSize()

"; 

%feature("docstring") simuPOP::Operator::infoField "

Description:

    get the information field specified by user (or by default)

Usage:

    x.infoField(idx)

"; 

%ignore simuPOP::Operator::formOffGenotype();

%ignore simuPOP::Operator::setFormOffGenotype(bool flag=true);

%feature("docstring") simuPOP::Operator::applyWithScratch "

Description:

    providing interface to apply operator before during or after
    mating.

Usage:

    x.applyWithScratch(pop, scratch, stage)

"; 

%feature("docstring") simuPOP::Operator::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::Operator::applyDuringMating "

Description:

    give pop, offspring, pop and mom.

Details:

    separator, persistant files, $gen etc substitution.

Usage:

    x.applyDuringMating(pop, offspring, *dad=NULL, *mom=NULL)

"; 

%feature("docstring") simuPOP::Operator::setOutput "

Description:

    ostream, if not set during construction.

Usage:

    x.setOutput(output=\"\", outputExpr=\"\")

"; 

%feature("docstring") simuPOP::Operator::getOstream "

Description:

    get output stream. This function is not exposed to user.

Usage:

    x.getOstream(*dict=NULL, readable=False)

"; 

%feature("docstring") simuPOP::Operator::closeOstream "

Description:

    close ostream and delete ostream pointer.. if it is a ofstream.

Usage:

    x.closeOstream()

"; 

%ignore simuPOP::Operator::atRepr();

%feature("docstring") simuPOP::Operator::__repr__ "

Description:

    simuPOP::Operator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::Operator::noOutput "

Description:

    if output=\">\". used internally

Usage:

    x.noOutput()

"; 

%ignore simuPOP::OstreamManager;

%feature("docstring") simuPOP::OstreamManager::OstreamManager "

Description:

    OStream Manager///////////////////////////////////////////////////
    ////////.

Usage:

    OstreamManager()

"; 

%feature("docstring") simuPOP::OstreamManager::~OstreamManager "

Description:

    simuPOP::OstreamManager::~OstreamManager

Usage:

    x.~OstreamManager()

"; 

%ignore simuPOP::OstreamManager::getOstream(const string &name, bool readable, bool realAppend, bool useString);

%ignore simuPOP::OstreamManager::hasOstream(const string &filename);

%ignore simuPOP::OstreamManager::listAll();

%feature("docstring") simuPOP::OstreamManager::closeAll "

Description:

    close all files and clean the map

Usage:

    x.closeAll()

"; 

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

    outputer is a (special) subclass of  Operator that will output
    files with different format.

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

    this function is very important

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

    outputHelper(str=\"\\\\n\", output=\">\", outputExpr=\"\",
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

    this function is very important

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

    simuPOP::outputHelper::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::parentsTagger "

Description:

    inherite tag from parents. If both parents have tags, use fathers.

"; 

%feature("docstring") simuPOP::parentsTagger::parentsTagger "

Description:

    constructor. default to be always active. string can be any string
    (m_Delimiter will be ignored for this class.) r will be replicate
    number g will be generation number.

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::parentsTagger::__repr__ "

Description:

    simuPOP::parentsTagger::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::parentsTagger::applyDuringMating "

Description:

    give pop, offspring, pop and mom.

Usage:

    x.applyDuringMating(pop, offspring, *dad=NULL, *mom=NULL)

"; 

%feature("docstring") simuPOP::pause "

Description:

    simuPOP::pause

"; 

%feature("docstring") simuPOP::pause::pause "

Description:

    stop simulation. press q to exit and any other key to continue

Usage:

    pause(prompt=True, stopOnKeyStroke=False, exposePop=True,
      popName=\"pop\", output=\">\", outputExpr=\"\", stage=PostMating,
      begin=0, end=-1, step=1, at=[], rep=REP_LAST, grp=GRP_ALL,
      infoFields=[])

Arguments:

    prompt:         if true (default), print prompt message.
    stopOnKeyStroke:if true, goon if no key was pressed
    exposePop:      whether or not expose pop to user namespace, only
                    useful when user choose 's' at pause. Default to
                    true.
    popName:        by which name the population is exposed? default
                    to 'pop'

"; 

%feature("docstring") simuPOP::pause::~pause "

Description:

    destructor

Usage:

    x.~pause()

"; 

%feature("docstring") simuPOP::pause::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pause::apply "

Description:

    simply output some info

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pause::__repr__ "

Description:

    simuPOP::pause::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::penetrance "

Description:

    penetrance

Details:

    Please refer to the user's guide for details.

"; 

%feature("docstring") simuPOP::penetrance::penetrance "

Description:

    If one field is specified, it will be used to store penetrance
    values. default to post mating.

Usage:

    penetrance(ancestralGen=-1, stage=DuringMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::penetrance::~penetrance "

Description:

    destructor

Usage:

    x.~penetrance()

"; 

%feature("docstring") simuPOP::penetrance::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::penetrance::penet "

Description:

    calculate/return penetrance etc

Usage:

    x.penet(*)

"; 

%feature("docstring") simuPOP::penetrance::apply "

Description:

    set pentrance to all individuals and record penetrance if
    requested.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::penetrance::applyDuringMating "

Description:

    set penetrance to all individual

Usage:

    x.applyDuringMating(pop, offspring, *dad=NULL, *mom=NULL)

"; 

%feature("docstring") simuPOP::penetrance::__repr__ "

Description:

    simuPOP::penetrance::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pointMutator "

Description:

    point mutator

Details:

    mutate specified individuals at specified loci to spcified allele.
    I.e., this is a non-random mutator used to introduce disease etc.

"; 

%feature("docstring") simuPOP::pointMutator::pointMutator "

Description:

    mutate once

Usage:

    pointMutator(atLoci, toAllele, atPloidy=[], inds=[], output=\">\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    atLoci:         a vector of loci index.
    inds:           mutate 'inds' individuals
    toAllele:       mutate to 'toAllele'

"; 

%feature("docstring") simuPOP::pointMutator::~pointMutator "

Description:

    destructor

Usage:

    x.~pointMutator()

"; 

%feature("docstring") simuPOP::pointMutator::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pointMutator::apply "

Description:

    apply!

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pointMutator::__repr__ "

Description:

    simuPOP::pointMutator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pointMutator::mutationCount "

Description:

    return mutation count

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

    simuPOP populations consists of individuals of the same genotypic
    structure, which refers to number of chromosome, number and
    position of loci on each chromosome etc. The most important
    components of a population are:
    * subpopulation. A population is divided into subpopulations
    (unstructured population has a single subpopulation, which is the
    whole population itself). Subpopulation structure limits the
    usually random exchange of genotypes between individuals
    disallowing mating between individuals from different
    subpopulation. In the presence of subpopualtion structure,
    exchange of genetic information across subpopulations can only be
    done through migration. Note that  simuPOP uses one-level
    population structure. That is to say, there is no sub-
    subpopulation or families in subpopulations.
    * variables. Every population has its own variable space, or
    localnamespaces  in  simuPOP term. This namespace is a Python
    dictionary that is attached to each population and can be exposed
    to the users through  vars()  or dvars()  function. Many functions
    and operators work in these namespaces and store their results in
    them. For example, function Stat  set variables like
    alleleFreq[loc]  and you can access them via
    pop.dvars().alleleFreq[loc][allele] .
    * ancestral generations. A population can save arbitrary number of
    ancestral generations. During evolution, the latest severl (or
    all) ancestral generations are saved. Functions to make a certain
    ancestray generation current  is provided so that one can examine
    and modify ancestral generations. Other important concepts like
    {information fields} are explained in class individual.
    It is worth noting that although a large number of member
    functions are provided, most of the operations are performed by
    operators . These functions will only be useful when you need to
    manipulate a population explicitly.

"; 

%feature("docstring") simuPOP::population::population "

Description:

    Create a population object with given size and genotypic
    structure.

Details:

    FIXME: Details of constructure is missing.This is techniquely the
    __init__  function of the population object.

Usage:

    population(size=0, ploidy=2, loci=[], sexChrom=False,
      lociPos=[], subPop=[], ancestralDepth=0, alleleNames=[],
      lociNames=[], maxAllele=MaxAllele, infoFields=[], chromMap=[])

Arguments:

    size:           population size. Can be ignored if subPop  is
                    specified. In that case, size  is the sum of
                    subPop . Default to 0 .
    ploidy:         number of sets of chromosomes. Default to 2
                    (diploid).
    loci:           an array of numbers of loci on each chromosome.
                    The length of parameter loci  determines the
                    number of chromosomes. Default to [1], meaning one
                    chromosome with a single locus.
                    The last chromosome can be sex chromosome. In this
                    case, the maximum number of loci on X and Y should
                    be provided. I.e., if there are 3 loci on Y
                    chromosme and 5 on X chromosome, use 5 .
    sexChrom:       true  or false . Diploid population only. If true
                    , the last homologous chromosome will be treated
                    as sex chromosome. (XY for male and XX for
                    female.) If X and Y have different number of loci,
                    number of loci of the longer one of the last (sex)
                    chromosome should be specified in loci .
    lociPos:        a 1-d or 2-d array specifying positions of loci on
                    each chromosome. You can use a nested array to
                    specify loci position for each chromosome. For
                    example, you can use lociPos=[1,2,3]  when
                    loci=[3]  or lociPos=[[1,2],[1.5,3,5]]  for
                    loci=[2,3] .  simuPOP does not assume a unit for
                    these locations, although they are usually
                    intepreted as centiMorgans. The default values are
                    1 , 2 , etc. on each chromosome.
    subPop:         an array of subpopulation sizes. Default value is
                    [size]  which means a single subpopulation of the
                    whole population. If both size  and subPop  are
                    are provided, subPop  should add up to size .
    ancestralDepth: number of most recent ancestral generations to
                    keep during evolution. Default to 0 , which means
                    only the current generation will be available. You
                    can set it to a positive number m to keep the
                    latest m generations in the population, or -1  to
                    keep all ancestral populations. Note that keeping
                    track of all ancestral populations may quickly
                    exhaust your computer RAM. If you really need to
                    do that, use  savePopulation  operator to save
                    each generation to a file is a much better choice.
    alleleNames:    an array of allele names. The first element should
                    be given to invalid/unknown allele. For example,
                    for a locus with alleles A,C,T,G, you can specify
                    alleleNames  as ('_','A','C','T','G') . Note that
                    simuPOP uses 1,2,3,4  internally and these names
                    will only be used for display purpose.
    lociNames:      an array or a matrix (separated by chromosomes) of
                    names for each loci. Default to \"locX-X\"  where
                    X-X  is a chromosome-loci index starting from 1.
    maxAllele:      maximum allele number. Default to the max allowed
                    allele states of current library (standard or long
                    allele version) maximum allele state for the whole
                    population. This will set a cap for all loci. For
                    individual locus, you can specify maxAllele in
                    mutation models, which can be smaller than global
                    maxAllele but not larger. Note that this number is
                    the number of allele states minus 1 since allele
                    number starts from 0.
    infoFields:     name of information fields that will be attached
                    to each individual. For example, if you need to
                    record the parents of each individual you will
                    need two, if you need to record the age of
                    individual, you need an additional one. Other
                    possibilities include offspring IDs etc. Note that
                    you have to plan this ahead of time since, for
                    example, tagger will need to know what info unit
                    to use. Default to none .

Examples:

>>> # a Wright-Fisher population
>>> WF = population(size=100, ploidy=1, loci=[1])
>>> 
>>> # a diploid population of size 10
>>> # there are two chromosomes with 5 and 7 loci respectively
>>> pop = population(size=10, ploidy=2, loci=[5, 7], subPop=[2, 8])
>>> 
>>> # a population with SNP markers (with names A,C,T,G
>>> #  range() are python functions
>>> pop = population(size=5, ploidy=2, loci=[5,10],
...     lociPos=[range(0,5),range(0,20,2)],
...     alleleNames=['A','C','T','G'],
...     subPop=[2,3], maxAllele=3)
>>>  


"; 

%ignore simuPOP::population::population(const population &rhs);

%feature("docstring") simuPOP::population::clone "

Description:

    deep copy of a population. (In python, pop1 = pop  will only
    create a reference to pop .)

Usage:

    x.clone(keepAncestralPops=-1)

"; 

%ignore simuPOP::population::swap(population &rhs);

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

%feature("docstring") simuPOP::population::__cmp__ "

Description:

    a python function used to compare the population class

Usage:

    x.__cmp__(rhs)

"; 

%feature("docstring") simuPOP::population::setSubPopStru "

Description:

    set population/subpopulation structure given subpopulation sizes

Usage:

    x.setSubPopStru(newSubPopSizes, allowPopSizeChange=False)

Arguments:

    subPopSize:     an array of subpopulation sizes. The population
                    may or may not change according to parameter
                    allowPopSizeChange  if the sum of subPopSize  does
                    not match popSize .
    allowPopSizeChange:if this parameter is true , popSize  can be
                    changed to the sum of subPopSize .

"; 

%feature("docstring") simuPOP::population::numSubPop "

Description:

    number of subpopulations in a population

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

    obtain total population size

Usage:

    x.popSize()

"; 

%feature("docstring") simuPOP::population::absIndIndex "

Description:

    return the absolute index of an individual in a subpopulation

Usage:

    x.absIndIndex(ind, subPop)

Arguments:

    index:          index of an individual in a subpopulation subPop
    subPop:         subpopulation index (start from 0 )

"; 

%feature("docstring") simuPOP::population::subPopIndPair "

Description:

    return the (sp, idx)  pair from an absolute index of an individual

Usage:

    x.subPopIndPair(ind)

"; 

%feature("docstring") simuPOP::population::subPopBegin "

Description:

    index of the first individual of a subpopulation

Usage:

    x.subPopBegin(subPop)

Arguments:

    subPop:         subpopulation index

"; 

%feature("docstring") simuPOP::population::subPopEnd "

Description:

    return the value of the index of the last individual of a
    subpopulation plus 1

Usage:

    x.subPopEnd(subPop)

Arguments:

    subPop:         subpopulation index

Note:

    as with all ...End functions, the returning index is out of the
    range so the actual range is [xxxBegin, xxxEnd). This agrees with
    all STL conventions and Python range.

"; 

%feature("docstring") simuPOP::population::ind "

Description:

    refernce to individual ind  in subpopulation subPop

Details:

    This function is named individual  in the Python interface.

Usage:

    x.ind(ind, subPop=0)

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

"; 

%feature("docstring") simuPOP::population::individuals "

Description:

    return an iterator that can be used to iterate through all
    individuals in subpopulation subPop

Usage:

    x.individuals(subPop)

"; 

%feature("docstring") simuPOP::population::ind "

Description:

    refernce to individual ind  in subpopulation subPop

Details:

    Return individual  from subpopulation . This function is named
    individual  in the Python interface.

Usage:

    x.ind(ind, subPop=0)

"; 

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

Details:

    Return an editable array of all genotypes of the population. You
    need to know how these genotypes are organized to safely
    read/write genotype directly. Individuals will be in order before
    exposing their genotypes.

Usage:

    x.arrGenotype(order)

Arguments:

    order:          if order is true , respect order; otherwise, do
                    not repect population structure.

"; 

%feature("docstring") simuPOP::population::arrGenotype "

Description:

    get the whole genotypes

Details:

    Return an editable array of all genotype of a subpopulation.
    Individuals will be in order before exposing their genotypes.

Usage:

    x.arrGenotype(subPop, order)

Arguments:

    subPop:         index of subpopulation (start from 0)
    order:          if order is true , keep order; otherwise, respect
                    subpop  structure.

"; 

%feature("docstring") simuPOP::population::exposeAffectedness "

Description:

    expose the (non-editable) effectedness fields of all individuals

Usage:

    x.exposeAffectedness(name=\"affected\")

"; 

%feature("docstring") simuPOP::population::setIndSubPopID "

Description:

    individual::setSubPopID,  individual::subPopID

Details:

    This function set information field for all individuals. Info can
    be used to sort individuals and set subpopulations. Therefore, the
    following code will do a migration:setIndSubPopID([
    an_array_of_info_value ]) setSubPopByIndID()

Usage:

    x.setIndSubPopID(id)

Arguments:

    info:           an array of info values, should have length of pop
                    size

"; 

%feature("docstring") simuPOP::population::setIndSubPopIDWithID "

Description:

    set information field of all individuals using their respective
    subpopulation ID

Usage:

    x.setIndSubPopIDWithID()

"; 

%feature("docstring") simuPOP::population::setSubPopByIndID "

Description:

    adjust subpopulation according to individual information values

Details:

    Assume individuals has subpopulation index in their information
    values and put them into corresponding subpopulations.

Usage:

    x.setSubPopByIndID(id=[])

Arguments:

    info:           optional info that will be used to set sub pop

"; 

%feature("docstring") simuPOP::population::splitSubPop "

Description:

    split a subpopulation into subpopulations of given sizes

Details:

    The sum of given sizes should be equal to the size of the split
    subpopulation. Subpopulation IDs can be specified. The
    subpopulation IDs of non-split subpopulations will be kept. For
    example, if subpopulation 1 of 0 1 2 3 is split into three parts,
    the new subpop id will be 0 (1 4 5) 2 3.

Usage:

    x.splitSubPop(which, sizes, subPopID=[])

"; 

%feature("docstring") simuPOP::population::splitSubPopByProportion "

Description:

    split a subpopulation into subpopulations of given proportions

Details:

    The sum of given proportions should add up to one. Subpopulation
    IDs can be specified.

Usage:

    x.splitSubPopByProportion(which, proportions, subPopID=[])

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
    will be no 'empty'  subpopulation left

Details:

    Remove specified subpopulations (and all individuals within). If
    shiftSubPopID  is false, subPopID  will be kept intactly.

Usage:

    x.removeSubPops(subPops=[], shiftSubPopID=True,
      removeEmptySubPops=False)

"; 

%feature("docstring") simuPOP::population::removeIndividuals "

Description:

    remove individuals

Details:

    If a valid subPop  is given, remove individuals from this
    subpopulation.

Usage:

    x.removeIndividuals(inds=[], subPop=-1,
      removeEmptySubPops=False)

"; 

%feature("docstring") simuPOP::population::mergeSubPops "

Description:

    merge populations

Details:

    Merge subpopulations, the first subpopulation ID (the first one in
    array subPops ) will be used as the ID of the new population.
    I.e., all subpopulations will take the ID of the first one.

Usage:

    x.mergeSubPops(subPops=[], removeEmptySubPops=False)

"; 

%feature("docstring") simuPOP::population::mergePopulation "

Description:

    merge populations by individuals

Details:

    Merge individuals from pop  to the current population. Two
    populations should have the same genotypic structures. By default,
    subpopulations of the merged populations are kept. I.e., if you
    merge two populations with one subpopulation, the resulting
    population will have two subpopulations. All ancestral generations
    are also merged.

Usage:

    x.mergePopulation(pop, newSubPopSizes=[], keepAncestralPops=-1)

Arguments:

    newSubPopSizes: subpopulation sizes can be specified. The overall
                    size should be the combined size of the two
                    populations. Because this parameter will be used
                    for all ancestral generations, it may fail if
                    ancestral generations have different sizes. To
                    avoid this problem, you may run mergePopulation
                    without this parameter, and then adjust
                    subpopulation sizes generation by generation.
    keepAncestralPops:ancestral populations to merge, default to all (-1
                    )

"; 

%feature("docstring") simuPOP::population::mergePopulationByLoci "

Description:

    merge populations by loci

Details:

    Two populations should have the same number of individuals. This
    also holds for any ancestral generations. By default, chromosomes
    of pop  are added to the current population. You can also specify
    new chromosome structure using parameter newLoci .

Usage:

    x.mergePopulationByLoci(pop, newNumLoci=[], newLociPos=[])

Arguments:

    newLoci:        the new number of loci for the combined genotypic
                    structure.

"; 

%feature("docstring") simuPOP::population::insertBeforeLoci "

Description:

    insert loci at given locations

Details:

    Insert loci at some given locations. In an inserted location,
    alleles will be zero.

Usage:

    x.insertBeforeLoci(idx, pos, names=[])

Arguments:

    idx:            an array of locus index. The loci will be inserted
                    before  each index. If you need to append to the
                    last locus, please use insertAfterLoci . If your
                    index is the first locus of a chromosome, the
                    inserted locus will become the first of that
                    chromosome. If you need to insert multiple loci
                    before a locus, please repeat that locus number.
    pos:            an array of locus positions. You need to make sure
                    that the position will make the inserted locus
                    between adjacent markers.
    names:          an array of locus names. If this parameter is not
                    given, some unique names such as \"insX_X\", will be
                    given.

"; 

%feature("docstring") simuPOP::population::insertBeforeLocus "

Description:

    insert an locus at given location.

Details:

    insertBeforeLocus(idx, pos, name)  is a shortcut to
    insertBeforeLoci([idx], [pos], [name])

Usage:

    x.insertBeforeLocus(idx, pos, name=string)

"; 

%feature("docstring") simuPOP::population::insertAfterLoci "

Description:

    append loci at given locations

Details:

    Append loci at some given locations. In an appended location,
    alleles will be zero.

Usage:

    x.insertAfterLoci(idx, pos, names=[])

Arguments:

    idx:            an array of locus index. The loci will be added
                    after  each index. If you need to append to the
                    first locus, please use insertBeforeLoci . If your
                    index is the last locus of a chromosome, the
                    appended locus will become the last of that
                    chromosome. If you need to append multiple loci
                    after a locus, please repeat that locus number.
    pos:            an array of locus positions. You need to make sure
                    that the position will make the appended locus
                    between adjacent markers.
    names:          an array of locus names. If this parameter is not
                    given, some unique names such as \"insX_X\", will be
                    given.

"; 

%feature("docstring") simuPOP::population::insertAfterLocus "

Description:

    append an locus at a given location

Details:

    insertAfterLocus(idx, pos, name)  is a shortcut to
    insertAfterLoci([idx], [pos], [name]) .

Usage:

    x.insertAfterLocus(idx, pos, name=string)

"; 

%feature("docstring") simuPOP::population::resize "

Description:

    resize population

Details:

    Resize population by giving new subpopulation sizes.

Usage:

    x.resize(newSubPopSizes, propagate=False)

Arguments:

    newSubPopSizes: an array of new subpopulation sizes. If there is
                    only one subpopulation, use [newPopSize] .
    propagate:      if propagate  is true , copy individuals to new
                    comers. I.e., 1, 2, 3 ==> 1, 2, 3, 1, 2, 3, 1

"; 

%feature("docstring") simuPOP::population::reorderSubPops "

Description:

    reorder subpopulations by order or by rank

Usage:

    x.reorderSubPops(order=[], rank=[], removeEmptySubPops=False)

Arguments:

    order:          new order of the subpopulations. For examples, 3 2
                    0 1 means subpop3 , subpop2 , subpop0 , subpop1
                    will be the new layout.
    rank:           you may also specify a new rank for each
                    subpopulation. For example, 3,2,0,1 means the
                    original subpopulations will have new IDs 3,2,0,1,
                    respectively. To achive order 3,2,0,1, the rank
                    should be 1 0 2 3.

"; 

%feature("docstring") simuPOP::population::newPopByIndID "

Details:

    Form a new population according to the parameter information.
    Information can be given directly as
    * keepAncestralPops=-1 : keep all
    * keepAncestralPops=0 : only current
    * keepAncestralPops=1 : keep one ...

Usage:

    x.newPopByIndID(keepAncestralPops=-1, id=[],
      removeEmptySubPops=False)

"; 

%feature("docstring") simuPOP::population::removeLoci "

Description:

    remove some loci from the current population. Loci that will be
    removed or kept can be specified.

Usage:

    x.removeLoci(remove=[], keep=[])

"; 

%feature("docstring") simuPOP::population::newPopWithPartialLoci "

Description:

    obtain a new population with selected loci

Details:

    Copy current population to a new one with selected loci and remove
    specified loci. (No change on the current population.)

Usage:

    x.newPopWithPartialLoci(remove=[], keep=[])

"; 

%feature("docstring") simuPOP::population::pushAndDiscard "

Description:

    (almost) internal use only

Usage:

    x.pushAndDiscard(rhs, force=False)

"; 

%feature("docstring") simuPOP::population::ancestralDepth "

Description:

    ancestral depth of the current population

Usage:

    x.ancestralDepth()

"; 

%feature("docstring") simuPOP::population::ancestralGen "

Description:

    return the current ancestral generation index

Usage:

    x.ancestralGen()

"; 

%feature("docstring") simuPOP::population::setIndInfo "

Description:

    set individual information for the given information field
    (index), for the given subpopulation (default to all population)

Usage:

    x.setIndInfo(values, idx)

"; 

%feature("docstring") simuPOP::population::setIndInfo "

Description:

    set individual information for the given information field (name),
    for the given subpopulation (default to all population)

Usage:

    x.setIndInfo(values, name)

"; 

%ignore simuPOP::population::infoBegin(UINT idx, bool order);

%ignore simuPOP::population::infoEnd(UINT idx, bool order);

%ignore simuPOP::population::infoBegin(UINT index, UINT subPop, bool order);

%feature("docstring") simuPOP::population::indInfo "

Description:

    Get information idx  of all individuals.

Usage:

    x.indInfo(idx, order)

Arguments:

    idx:            index in all information fields
    order:          if true, sort returned vector in individual order

"; 

%feature("docstring") simuPOP::population::indInfo "

Description:

    Get information name  of all individuals.

Usage:

    x.indInfo(name, order)

Arguments:

    name:           name of the information field
    order:          if true, sort returned vector in individual order

"; 

%feature("docstring") simuPOP::population::indInfo "

Description:

    Get information idx  of all individuals in a subpopulation subPop
    .

Usage:

    x.indInfo(idx, subPop, order)

Arguments:

    idx:            index in all information fields
    subPop:         subpopulation index
    order:          if true, sort returned vector in individual order

"; 

%feature("docstring") simuPOP::population::indInfo "

Description:

    Get information name  of all individuals in a subpopulation subPop
    .

Usage:

    x.indInfo(name, subPop, order)

Arguments:

    name:           name of the information field
    subPop:         subpopulation index
    order:          if true, sort returned vector in individual order

"; 

%feature("docstring") simuPOP::population::arrIndInfo "

Description:

    Get an editable array (Python list) of all information fields.

Usage:

    x.arrIndInfo(order)

Arguments:

    order:          whether or not the list has the same order as
                    individuals

"; 

%feature("docstring") simuPOP::population::addInfoField "

Details:

    add an information field to a population. For example,
    \"penetrance\", for the analysis of an existing population. If the
    field already exists in the geno structure, matches local info
    fields with the structure.

Usage:

    x.addInfoField(field, init=0)

"; 

%feature("docstring") simuPOP::population::addInfoFields "

Description:

    add one or more information fields to a population

Usage:

    x.addInfoFields(fields, init=0)

"; 

%feature("docstring") simuPOP::population::setInfoFields "

Description:

    set fields

Usage:

    x.setInfoFields(fields, init=0)

"; 

%feature("docstring") simuPOP::population::setAncestralDepth "

Description:

    set ancestral depth, can be -1

Usage:

    x.setAncestralDepth(depth)

"; 

%feature("docstring") simuPOP::population::ancestralPop "

Description:

    current ancestral population selected (0 for the latest
    generation)

Usage:

    x.ancestralPop()

"; 

%feature("docstring") simuPOP::population::useAncestralPop "

Description:

    use an ancestral population. 0  for the latest generation.

Usage:

    x.useAncestralPop(idx)

"; 

%feature("docstring") simuPOP::population::equalTo "

Description:

    compare two populations

Usage:

    x.equalTo(rhs)

"; 

%ignore simuPOP::population::adjustGenoPosition(bool order);

%ignore simuPOP::population::adjustInfoPosition(bool order);

%feature("docstring") simuPOP::population::savePopulation "

Description:

    global function loadPopulation

Usage:

    x.savePopulation(filename, format=\"auto\", compress=True)

Arguments:

    filename:       save to filename
    format:         format to save. Can be one of the following:
                    'text', 'bin', or 'xml'. The default format is
                    'text' but the output is not supposed to be read.
                    'bin' has smaller size than the other two and
                    should be used for large populations. 'xml' is the
                    most readable format and should be used when you
                    would like to convert  simuPOP populations to
                    other formats.

"; 

%ignore simuPOP::population::loadPopulation(const string &filename, const string &format="auto");

%feature("docstring") simuPOP::population::rep "

Description:

    current replicate in a simulator

Usage:

    x.rep()

"; 

%ignore simuPOP::population::setRep(int rep, bool setVar=true);

%feature("docstring") simuPOP::population::grp "

Description:

    current group ID in a simulator

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

    return variables of a population. If subPop  is given, return a
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

%ignore simuPOP::population::varsAsString() const ;

%ignore simuPOP::population::varsFromString(const string &vars);

%feature("docstring") simuPOP::population::evaluate "

Description:

    evaluate a python statment/expression

Details:

    This function evaluates a python statment/expression and return
    its result as a string. Optionally run statement first.

Usage:

    x.evaluate(expr=\"\", stmts=\"\")

"; 

%feature("docstring") simuPOP::population::execute "

Description:

    evaluate a statement (can be multi-line string)

Usage:

    x.execute(stmts=\"\")

"; 

%feature("docstring") simuPOP::population::rearrangeLoci "

Description:

    rearrange loci on chromosomes, e.g. combine two chromosomes into
    one

Details:

    This is used by mergeByLoci

Usage:

    x.rearrangeLoci(newNumLoci, newLociPos)

"; 

%feature("docstring") simuPOP::pyEval "

Description:

    evaluate an expression.

"; 

%feature("docstring") simuPOP::pyEval::pyEval "

Description:

    evaluate expr/statments in local replicate namespace

Usage:

    pyEval(expr=\"\", stmts=\"\", preStmts=\"\", postStmts=\"\",
      exposePop=False, name=\"\", output=\">\", outputExpr=\"\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    expr:           expression to be evaluated. Its result will be
                    sent to output.
    stmts:          statement that will be executed before the
                    expression.
    preStmts:       statement that will be executed when the operaot
                    is constructed.
    postStmts:      statement that will be executed when the operator
                    is destroyed.
    exposePop:      if true, expose current pop as variable \"pop\"
    name:           used to let pure python operator to identify
                    themselves
    output:         default to \">\" . i.e., output to standard output.
                    for usage of other parameters, see
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::pyEval::~pyEval "

Description:

    simuPOP::pyEval::~pyEval

Usage:

    x.~pyEval()

"; 

%feature("docstring") simuPOP::pyEval::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyEval::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyEval::__repr__ "

Description:

    simuPOP::pyEval::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyEval::name "

Description:

    simuPOP::pyEval::name

Usage:

    x.name()

"; 

%feature("docstring") simuPOP::pyExec "

Description:

    evaluate an expression.

"; 

%feature("docstring") simuPOP::pyExec::pyExec "

Description:

    evaluate statments in local replicate namespace, no return value

Usage:

    pyExec(stmts=\"\", preStmts=\"\", postStmts=\"\", exposePop=False,
      name=\"\", output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
      end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    stmts:          statements (a single or multi-line string) that
                    will be evaluated before the expression.
    preStmts:       statement that will be executed when the operaot
                    is constructed.
    postStmts:      statement that will be executed when the operator
                    is destroyed.
    exposePop:      if true, expose current pop as variable \"pop\"
    output:         default to \">\" . i.e., output to standard output.
                    for usage of other parameters, see
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::pyExec::~pyExec "

Description:

    simuPOP::pyExec::~pyExec

Usage:

    x.~pyExec()

"; 

%feature("docstring") simuPOP::pyExec::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyExec::__repr__ "

Description:

    simuPOP::pyExec::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyIndOperator "

Description:

    simuPOP::pyIndOperator

"; 

%feature("docstring") simuPOP::pyIndOperator::pyIndOperator "

Description:

    individual operator that apply a function to each individual

Usage:

    pyIndOperator(*func, *param=NULL, stage=PostMating,
      formOffGenotype=False, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    func:           a python function that accept a individual and
                    perform whatever operation it wants to.
    para:           any python object that will be passed to func
                    after pop parameter. Multiple parameter can be
                    passed as a tuple.
    formOffGenotype:if stage=DuringMating, set this parameter to false
                    will disallow random mating to set genotype. Note:
                    (FIXME) output to output or outputExpr is not yet
                    supported. Ideally, this func will take two
                    parameters with pop and then a filehandle to
                    output, however, differentiating output, append
                    etc is too troublesome right now.

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyIndOperator::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyIndOperator::__repr__ "

Description:

    simuPOP::pyIndOperator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyInit "

Description:

    simuPOP::pyInit

"; 

%feature("docstring") simuPOP::pyInit::pyInit "

Description:

    initialize populations using given user function.

Details:

    User of this operator must supply a Python function with parameter
    (index, ploidy, subpop). This operator will loop through all
    individual in each subpop and call this function to initialize
    populations.The arrange of parameters allows different
    initialization scheme for each subpop.

Usage:

    pyInit(*func, subPop=[], atLoci=[], atPloidy=-1,
      indRange=intMatrix, maleFreq=0.5, sex=[], stage=PreMating,
      begin=0, end=1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[])

Arguments:

    func:           a python function with parameter (index, ploidy,
                    subpop) index is the allele index (0 ~
                    totNumLoci()-1), ploidy (index to copy of
                    chromosomes), subpop (subpop number). The return
                    value of this function should be a integer.
    atLoci:         a vector of loci indices. If empty, apply to all
                    loci
    atPloidy:       initialize which copy of chromosomes. Default to
                    all.
    stage:          is et to PreMating. Other parameters please refer
                    to help(baseOperator.__init__)

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyInit::__repr__ "

Description:

    simuPOP::pyInit::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyInit::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyMating "

Details:

    python mating. Parental and offspring generation, along with
    during mating operators, are passed to a python function. All
    mating are done there, and the resulting population be
    returned.This process will be slow and should be used mainly for
    prototyping or demonstration purposes.

"; 

%feature("docstring") simuPOP::pyMating::pyMating "

Description:

    constructor, no new subPopsize parameter

Usage:

    pyMating(*func=NULL, newSubPopSize=[], newSubPopSizeExpr=\"\",
      *newSubPopSizeFunc=NULL)

"; 

%feature("docstring") simuPOP::pyMating::~pyMating "

Description:

    destructor

Usage:

    x.~pyMating()

"; 

%feature("docstring") simuPOP::pyMating::clone "

Description:

    clone() const. The same as  mating::clone() const.

Usage:

    x.clone()

"; 

%ignore simuPOP::pyMating::pyMating(const pyMating &rhs);

%feature("docstring") simuPOP::pyMating::__repr__ "

Description:

    return name of the mating type

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyMating::mate "

Description:

    do the mating. --- no mating :-)

Details:

    All individuals will be passed to during mating operators but no
    one will die (ignore during mating failing signal).

Usage:

    x.mate(pop, scratch, ops, submit)

"; 

%feature("docstring") simuPOP::pyMigrator "

Description:

    migrate using given info vector

Details:

    You can use directmigrator to accomplish any migration: that is to
    say you directly specify subpopulation numbers for each individual
    and this operator will do the rest.

"; 

%feature("docstring") simuPOP::pyMigrator::pyMigrator "

Description:

    create a directmigrator

Details:

    This operator accept a one-dimensional Numeric Python int array.
    (created by Numeric.array ). The contend of the array will be
    considered as subpopulation id.

Usage:

    pyMigrator(*subPopID=NULL, stage=PreMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    subPopID:       a 1-d array (list, typle, carray). Must has length
                    greater or equal to population size.
    stage:          is default to PreMating, please refer to
                    help(baseOperator.__init__) for details about
                    other parameters.

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyMigrator::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyMigrator::__repr__ "

Description:

    simuPOP::pyMigrator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyMutator "

Description:

    mixed mutation model . has not been implemented.

"; 

%feature("docstring") simuPOP::pyMutator::pyMutator "

Description:

    simuPOP::pyMutator::pyMutator

Usage:

    pyMutator(rate=[], atLoci=[], maxAllele=0, *func=NULL,
      output=\">\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::pyMutator::~pyMutator "

Description:

    simuPOP::pyMutator::~pyMutator

Usage:

    x.~pyMutator()

"; 

%ignore simuPOP::pyMutator::pyMutator(const pyMutator &rhs);

%feature("docstring") simuPOP::pyMutator::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyMutator::mutate "

Description:

    how to mutate a single allele. this is usually the only function
    that need to be defined by the subclasses.

Usage:

    x.mutate(allele)

"; 

%feature("docstring") simuPOP::pyMutator::__repr__ "

Description:

    simuPOP::pyMutator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyOperator "

Description:

    simuPOP::pyOperator

"; 

%feature("docstring") simuPOP::pyOperator::pyOperator "

Description:

    python operator, using a function that accept a population object

Usage:

    pyOperator(*func, *param=NULL, stage=PostMating,
      formOffGenotype=False, passOffspringOnly=False, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    func:           a python function that accept a population and
                    perform whatever operation it wants to.
    para:           any python object that will be passed to func
                    after pop parameter. Multiple parameter can be
                    passed as a tuple.
    formOffGenotype:if stage=DuringMating, set this parameter to false
                    will disallow random mating to set genotype.
    passOffspringOnly:Default to false. If true,  pyOperator will expect
                    a function of form func(off, param), instead of
                    func(pop, off, dad, mon, param) when
                    passOffspringOnly is false. Since many
                    duringMating  pyOperator only need access to
                    offspring, this will imporve efficiency. Note:
                    (FIXME) output to output or outputExpr is not yet
                    supported. Ideally, this func will take two
                    parameters with pop and then a filehandle to
                    output, however, differentiating output, append
                    etc is too troublesome right now.

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyOperator::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyOperator::applyDuringMating "

Description:

    give pop, offspring, pop and mom.

Usage:

    x.applyDuringMating(pop, offspring, *dad=NULL, *mom=NULL)

"; 

%feature("docstring") simuPOP::pyOperator::__repr__ "

Description:

    simuPOP::pyOperator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyPenetrance "

Description:

    penetrance using user supplied function

Details:

    Assign penetrance value by calling a user supplied function

"; 

%feature("docstring") simuPOP::pyPenetrance::pyPenetrance "

Description:

    provide locus and penetrance for 11, 12, 13 (in the form of
    dictionary)

Usage:

    pyPenetrance(loci, *func, ancestralGen=-1, stage=DuringMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[])

Arguments:

    loci:           susceptibility loci. The genotype at these loci
                    will be passed to func.
    func:           a Python function that accept genotypes at
                    susceptibility loci and return penetrance value.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyPenetrance::penet "

Description:

    currently assuming diploid

Usage:

    x.penet(*ind)

"; 

%feature("docstring") simuPOP::pyPenetrance::__repr__ "

Description:

    simuPOP::pyPenetrance::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyQuanTrait "

Description:

    quantitative trait using user supplied function

Details:

    Assign qtrait value by calling a user supplied function

"; 

%feature("docstring") simuPOP::pyQuanTrait::pyQuanTrait "

Description:

    provide locus and qtrait for 11, 12, 13 (in the form of
    dictionary)

Usage:

    pyQuanTrait(loci, *func, ancestralGen=-1, stage=PostMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[\"qtrait\"])

Arguments:

    loci:           susceptibility loci. The genotype at these loci
                    will be passed to func.
    func:           a Python function that accept genotypes at
                    susceptibility loci and return qtrait value.
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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyQuanTrait::qtrait "

Description:

    currently assuming diploid

Usage:

    x.qtrait(*ind)

"; 

%feature("docstring") simuPOP::pyQuanTrait::__repr__ "

Description:

    simuPOP::pyQuanTrait::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pySample "

Description:

    thrink population accroding to some outside value

"; 

%feature("docstring") simuPOP::pySample::pySample "

Description:

    create a python sampler

Usage:

    pySample(*keep, keepAncestralPops=-1, name=\"sample\",
      nameExpr=\"\", times=1, saveAs=\"\", saveAsExpr=\"\", format=\"auto\",
      stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
      grp=GRP_ALL, infoFields=[])

Arguments:

    keep:           a carray of the length of population. its values
                    will be assigned to info.
    keepAncestralPop::-1 (all), 0 (no), 1(one ancestral pop) and so on.
                    and other parameters please see
                    help(baseOperator.__init__)

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pySample::drawsample "

Description:

    simuPOP::pySample::drawsample

Usage:

    x.drawsample(pop)

"; 

%feature("docstring") simuPOP::pySample::__repr__ "

Description:

    simuPOP::pySample::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pySelector "

Description:

    selection using user supplied function

Details:

    Assign fitness value by calling a user supplied function

"; 

%feature("docstring") simuPOP::pySelector::pySelector "

Description:

    provide locus and fitness for 11, 12, 13 (in the form of
    dictionary)

Usage:

    pySelector(loci, *func, stage=PreMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[\"fitness\"])

Arguments:

    loci:           susceptibility loci. The genotype at these loci
                    will be passed to func.
    func:           a Python function that accept genotypes at
                    susceptibility loci generation number, and return
                    fitness value.
    output:         and other parameters please refer to
                    help(baseOperator.__init__)

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pySelector::indFitness "

Description:

    currently assuming diploid

Usage:

    x.indFitness(*ind, gen)

"; 

%feature("docstring") simuPOP::pySelector::__repr__ "

Description:

    simuPOP::pySelector::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pySubset "

Description:

    thrink population accroding to some outside value

"; 

%feature("docstring") simuPOP::pySubset::pySubset "

Description:

    create a directmigrator

Usage:

    pySubset(keep=[], stage=PostMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    keep:           a carray of the length of population. its values
                    will be assigned to info.  and other parameters
                    please see help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::pySubset::~pySubset "

Description:

    destructor

Usage:

    x.~pySubset()

"; 

%feature("docstring") simuPOP::pySubset::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pySubset::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pySubset::__repr__ "

Description:

    simuPOP::pySubset::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyTagger "

Details:

    This tagger take some information fields from both parents, pass
    to a python function and set individual field with the return
    value.This operator can be used to trace the inheritance of trait
    values.

"; 

%feature("docstring") simuPOP::pyTagger::pyTagger "

Description:

    simuPOP::pyTagger::pyTagger

Usage:

    pyTagger(*func=NULL, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    infoFields:     information fields. The user should gurantee the
                    existence of these fields.
    func:           a pyton function that return a list to assign the
                    information fields. e.g. if fields=['A', 'B'], the
                    function will pass values of fields 'A' and 'B' of
                    father, followed by mother if there is one, to
                    this function. The returned value is assigned to
                    fields 'A' and 'B' of the offspring. The returned
                    value has to be a list even if only one field is
                    given.

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyTagger::__repr__ "

Description:

    simuPOP::pyTagger::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyTagger::applyDuringMating "

Description:

    give pop, offspring, pop and mom.

Usage:

    x.applyDuringMating(pop, offspring, *dad=NULL, *mom=NULL)

"; 

%feature("docstring") simuPOP::PythonCoutBuf "

Description:

    create a stream buf that print to python sys.stdout cout will be
    redirected to this buf to really output to python console.

"; 

%feature("docstring") simuPOP::PythonCoutBuf::PythonCoutBuf "

Description:

    simuPOP::PythonCoutBuf::PythonCoutBuf

Usage:

    PythonCoutBuf()

"; 

%feature("docstring") simuPOP::quanTrait "

Description:

    quantitative trait

Details:

    Genetic quantitative trait is tricky to simulate. In  simuPOP, I
    employee an ability (fitness) to mate approach. Namely, the
    probability that an individual will be chosen for mating is
    proportional to its fitness value. More specifically,
    * PreMating selectors assign fitness values to each individual.
    * Sexless mating (e.g.  binomialSelection) : individuals are
    chosen at probabilities that are proportional to their fitness
    values. More specifically, if there are N individuals with fitness
    values  $f_i, i=1,...,N $, individual  $i$ will have probability
    $ \\\\frac{f_i}{\\\\sum_{j=1}^N f_j} $ to be chosen to be passed to the
    next generation.
    * Random mating with sex (e.g. randommating): males and females
    are separated and each are chosen as described above.Please refer
    to the user's guide for details.

"; 

%feature("docstring") simuPOP::quanTrait::quanTrait "

Description:

    constructor. default to be always active.

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::quanTrait::qtrait "

Description:

    calculate/return quantitative trait etc

Usage:

    x.qtrait(*)

"; 

%feature("docstring") simuPOP::quanTrait::apply "

Description:

    set qtrait to all individual

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::quanTrait::__repr__ "

Description:

    simuPOP::quanTrait::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::randomMating "

Details:

    basic sexual random mating.Within each subpopulation, choose male
    and female randomly randmly get one copy of chromosome from
    father/mother.require: sexed individual; ploidy == 2apply during
    mating operators and put into the next generation.if
    ignoreParentsSex is set, parents will be chosen regardless of
    sex.Otherwise, male and female will be collected and be chosen
    randomly.If there is no male or female in a subpopulation, if
    m_UseSameSexIfUniSex is true, an warning will be generated and
    same sex mating (?) will be used otherwise,  randomMating will
    return false.if there is no during mating operator to copy
    alleles, a direct copy will be used.

"; 

%feature("docstring") simuPOP::randomMating::randomMating "

Description:

    create a random mating scheme

Usage:

    randomMating(numOffspring=1., *numOffspringFunc=NULL,
      maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
      *newSubPopSizeFunc=NULL, newSubPopSizeExpr=\"\",
      contWhenUniSex=True)

Arguments:

    numOffspring:   
    number:         of offspring or p in some modes
    numOffspringFunc:
    a:              python function that determine number of offspring
                    or p depending on mode
    maxNumOffspring:used when mode=MATE_BinomialDistribution
    mode:           one of MATE_NumOffspring ,
                    MATE_NumOffspringEachFamily,
                    MATE_GeometricDistribution,
                    MATE_PoissonDistribution,
                    MATE_BinomialDistribution
    newSubPopSize:  an array of subpopulation sizes, should have the
                    same number of subpopulations as current
                    population
    newSubPopSizeExpr:an expression that will be evaluated as an array
                    of subpop sizes
    newSubPopSizeFunc:an function that have parameter gen and oldSize
                    (current subpop size).
    contWhenUniSex: continue when there is only one sex in the
                    population, default to true

"; 

%feature("docstring") simuPOP::randomMating::~randomMating "

Description:

    destructor

Usage:

    x.~randomMating()

"; 

%feature("docstring") simuPOP::randomMating::clone "

Description:

    clone() const. Generate a copy of itself and return pointer this
    is to make sure the object is persistent and will not be freed by
    python.

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::randomMating::isCompatible "

Description:

    check if the mating type is compatible with population structure

Details:

    possible things to check:
    * need certain types of individual (age, sex etc)
    * need resizeable population...

Usage:

    x.isCompatible(pop)

"; 

%feature("docstring") simuPOP::randomMating::__repr__ "

Description:

    return name of the mating type

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::randomMating::submitScratch "

Description:

    simuPOP::randomMating::submitScratch

Usage:

    x.submitScratch(pop, scratch)

"; 

%feature("docstring") simuPOP::randomMating::mate "

Description:

    do the mating. parameters see  mating::mate .

Details:

    Within each subpopulation, choose male and female randomly randmly
    get one copy of chromosome from father/mother.require: sexed
    individual; ploidy == 2apply during mating operators and put into
    the next generation.Otherwise, male and female will be collected
    and be chosen randomly.
    * If there is no male or female in a subpopulation,
    * if m_contWhenUniSex is true, an warning will be generated and
    same sex mating (?) will be used
    * otherwise,  randomMating will return false.

Usage:

    x.mate(pop, scratch, ops, submit)

"; 

%feature("docstring") simuPOP::randomSample "

Description:

    thrink population accroding to some outside value

"; 

%feature("docstring") simuPOP::randomSample::randomSample "

Description:

    draw random sample, regardless of affected status

Usage:

    randomSample(size=[], name=\"sample\", nameExpr=\"\", times=1,
      saveAs=\"\", saveAsExpr=\"\", format=\"auto\", stage=PostMating,
      begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
      infoFields=[])

Arguments:

    size:           size of sample. It can be either a number,
                    representing the overall sample size, regardless
                    of population strucutre; or an array, representing
                    number of samples drawn from each subpopulation.
    stage:          and other parameters please see
                    help(baseOperator.__init__)
    name:           variable name of the sampled population (will be
                    put in pop local namespace)
    nameExpr:       expression version of name. If both name and
                    nameExpr is empty, do not store pop.
    times:          how many times to run the sample process? This is
                    usually one, but we may want to take several
                    random samples.
    saveAs:         filename to save the population.
    saveAsExpr:     expression for save filename
    format:         to save sample(s)
    stage:          and other parameters please see
                    help(baseOperator.__init__)

Note:

    ancestral populations will not be copied to the samples

"; 

%feature("docstring") simuPOP::randomSample::~randomSample "

Description:

    destructor

Usage:

    x.~randomSample()

"; 

%feature("docstring") simuPOP::randomSample::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::randomSample::prepareSample "

Description:

    value checking

Usage:

    x.prepareSample(pop)

"; 

%feature("docstring") simuPOP::randomSample::drawsample "

Description:

    simuPOP::randomSample::drawsample

Usage:

    x.drawsample(pop)

"; 

%feature("docstring") simuPOP::randomSample::__repr__ "

Description:

    simuPOP::randomSample::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::recombinator "

Description:

    Recombination.

Details:

    * only works for diploids (and for females in haplodiploids)
    population.
    * Free recombination between loci. Loci behave completely
    independently.
    * otherwise there will be some linkage between loci, user need to
    specify physical recombination rate between adjacent loci (ie
    between locus n and n+1)
    * The recombination rate must be comprised between 0.0 and 0.5.
    * A recombination rate of 0.0 means that the loci are completely
    linked and thus behave together as a single linked locus.
    * A recombination rate of 0.5 is equivalent to free recombination.
    * All values in between will represent various linkage intensities
    between adjacent pairs of loci. The recombination rate is
    equivalent to 1-linkage and represents the probability that the
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

    intensity:      recombination rate per unit of loci distance.
                    I.e., the really recombination rate between two
                    loci is determined by intensity*loci distance
                    between them.
    rate:           recombination rate regardless of loci distance; it
                    can also be an array of recombination rates. Must
                    be the same length as afterLoci or totNumOfLoci().
                    If totNumLociThe last item can be ignored.
    afterLoci:      an array of loci index. If rates is also
                    specified, they should have the same length.
                    Default to all loci (but meaningless for those
                    loci locate at the end of chromosome.) If given,
                    afterLoci should be ordered, and can not include
                    loci at the end of a chromosome.  recombination
                    rate for male individuals. If given, parameter
                    rate will be considered as female rate.
                    recombination intensity for male individuals. If
                    given, parameter intensity will be considered as
                    female intensity.  if given, male will recombine
                    at different locations. This is rarely used.

Note:

    rate gives recombination rate PER unit, rates gives plain rates.
    there is no recombination between sex chromosomes of male
    individuals. (sexChrom=True).

"; 

%feature("docstring") simuPOP::recombinator::~recombinator "

Description:

    simuPOP::recombinator::~recombinator

Usage:

    x.~recombinator()

"; 

%feature("docstring") simuPOP::recombinator::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::recombinator::__repr__ "

Description:

    simuPOP::recombinator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::recombinator::prepareRecRates "

Description:

    this function takes intensity, rate, afterLoci, ... inputs and
    return a bernulli trailer and a recBeforeLoci vector.

Details:

    get loci distance * rate and then recombinant pointsget loci
    distance * rate and then recombinant pointsinitialize
    recombination counter, This will count recombination events after
    each locus.

Usage:

    x.prepareRecRates(pop, intensity, rate, afterLoci, sexChrom,
      recBeforeLoci, vecP)

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

%feature("docstring") simuPOP::recombinator::recombine "

Description:

    simuPOP::recombinator::recombine

Usage:

    x.recombine(*parent, offspring, offPloidy, bt, recBeforeLoci,
      setSex=False)

"; 

%feature("docstring") simuPOP::recombinator::applyDuringMating "

Description:

    give pop, offspring, pop and mom.

Usage:

    x.applyDuringMating(pop, offspring, *dad=NULL, *mom=NULL)

"; 

%feature("docstring") simuPOP::RNG "

Description:

    Random number generator///////////////////////////////////////////
    ////////////////.

"; 

%feature("docstring") simuPOP::RNG::RNG "

Description:

    Random number generator///////////////////////////////////////////
    ////////////////.

Usage:

    RNG(*rng=NULL, seed=0)

"; 

%feature("docstring") simuPOP::RNG::~RNG "

Description:

    simuPOP::RNG::~RNG

Usage:

    x.~RNG()

"; 

%feature("docstring") simuPOP::RNG::setRNG "

Description:

    choose an random number generator. This can be done by setting
    GSL_RNG_TYPE as well.

Usage:

    x.setRNG(*rng=NULL, seed=0)

"; 

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

    x.randIntArray(n, size, *vec)

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

    x.randUniform01Array(size, *vec)

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

    sample from population and save samples sample operator will
    generate a new subpopulation in pop namespace.

"; 

%feature("docstring") simuPOP::sample::sample "

Description:

    create a sample

Usage:

    sample(name=\"sample\", nameExpr=\"\", times=1, saveAs=\"\",
      saveAsExpr=\"\", format=\"auto\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    name:           variable name of the sampled population (will be
                    put in pop local namespace)
    nameExpr:       expression version of name. If both name and
                    nameExpr is empty, do not store pop.
    times:          how many times to run the sample process? This is
                    usually one, but we may want to take several
                    random samples.
    saveAs:         filename to save the population.
    saveAsExpr:     expression for save filename
    format:         to save sample(s)
    stage:          and other parameters please see
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::sample::~sample "

Description:

    destructor

Usage:

    x.~sample()

"; 

%feature("docstring") simuPOP::sample::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::sample::prepareSample "

Description:

    simuPOP::sample::prepareSample

Usage:

    x.prepareSample()

"; 

%feature("docstring") simuPOP::sample::drawsample "

Description:

    simuPOP::sample::drawsample

Usage:

    x.drawsample(pop)

"; 

%feature("docstring") simuPOP::sample::samples "

Description:

    return the samples

Usage:

    x.samples(pop)

"; 

%feature("docstring") simuPOP::sample::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::sample::__repr__ "

Description:

    simuPOP::sample::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::sample::saveIndIndex "

Description:

    service functions that can be used by everyone save the idx of
    each individual to a filed usually 'oldindex'

Usage:

    x.saveIndIndex(pop, indexField=\"oldindex\")

"; 

%feature("docstring") simuPOP::sample::resetParentalIndex "

Description:

    reset father_idx and mother_idx using the samed

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

    set all sub pop id to -1 (remove)

Usage:

    x.resetSubPopID(pop)

"; 

%feature("docstring") simuPOP::savePopulation "

Description:

    save population to a file

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::savePopulation::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::savePopulation::__repr__ "

Description:

    simuPOP::savePopulation::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::selector "

Description:

    selection

Details:

    Genetic selection is tricky to simulate. In  simuPOP, I employee
    an ability (fitness) to mate approach. Namely, the probability
    that an individual will be chosen for mating is proportional to
    its fitness value. More specifically,
    * PreMating selectors assign fitness values to each individual.
    * Sexless mating (e.g.  binomialSelection) : individuals are
    chosen at probabilities that are proportional to their fitness
    values. More specifically, if there are N individuals with fitness
    values  $f_i, i=1,...,N $, individual  $i$ will have probability
    $ \\\\frac{f_i}{\\\\sum_{j=1}^N f_j} $ to be chosen to be passed to the
    next generation.
    * Random mating with sex (e.g. randommating): males and females
    are separated and each are chosen as described above.Please refer
    to the user's guide for details.

"; 

%feature("docstring") simuPOP::selector::selector "

Description:

    constructor. default to be always active.

Usage:

    selector(stage=PreMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[\"fitness\"])

"; 

%feature("docstring") simuPOP::selector::~selector "

Description:

    destructor

Usage:

    x.~selector()

"; 

%feature("docstring") simuPOP::selector::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::selector::indFitness "

Description:

    calculate/return w11 etc

Usage:

    x.indFitness(*, gen)

"; 

%feature("docstring") simuPOP::selector::apply "

Description:

    set fitness to all individual

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::selector::__repr__ "

Description:

    simuPOP::selector::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::setAncestralDepth "

Description:

    simuPOP::setAncestralDepth

"; 

%feature("docstring") simuPOP::setAncestralDepth::setAncestralDepth "

Description:

    timer if called, output time passed since last calling time.

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::setAncestralDepth::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::setAncestralDepth::__repr__ "

Description:

    simuPOP::setAncestralDepth::__repr__

Usage:

    x.__repr__()

"; 

%ignore simuPOP::SharedVariables;

%ignore simuPOP::SharedVariables::SharedVariables();

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

%ignore simuPOP::SharedVariables::removeVar(const string &name);

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

%ignore simuPOP::SharedVariables::asString() const ;

%ignore simuPOP::SharedVariables::fromString(const string &vars);

%feature("docstring") simuPOP::simulator "

Description:

    simulator manage several replicates of a population, evolve them
    using given mating scheme and operators.

Details:

    simulators combine three important components of  simuPOP:
    population, mating scheme and operators together. A simulator is
    created with an instance of population, a replicate number and a
    mating scheme. It makes 'rep' replicates of this population and
    control the evolution process of these populations.The most
    important functions of a simulator is of course  evolve(). It
    accepts arrays of operators as its parameters, among which,
    'preOps' and 'postOps' will be applied to the populations at the
    begining and end of evolution, whereas 'operators'will be applied
    at every generation.simulators separates operators into pre-,
    during- and post- mating operators. During evolution, simulator
    first apply all pre-mating operators and then call the mate()
    function of the given mating scheme, which will call during-mating
    operators during the birth of each offsrping. After the mating is
    finished, post-mating operators are applied in the order they
    apprear in the operator list.Since operators can apply to specific
    replicate or replicates group, and might not be active at all
    time, the isActive(m_curRep, m_numRep, m_gen, end,  grp())
    function of each operator is called before it is applied to the
    populations.simulators can evolve a given number of generations
    (the 'end' parameter of evolve), or evolve indefinitely using a
    certain type of operators called terminators. In this case, one or
    more terminators will check the status of evolution and determine
    if the simulation should be stopped. An obvious example of such a
    terminator is a fixation-checker.Finally, a simulator can be saved
    to a file in the format of 'text', 'bin', or 'xml'. This enables
    us to stop a simulation and ressume it at another time or on
    another machine. It is also a good idea to save a snapshot of a
    simulation every several generations.DEVONLY{This is a template
    class that manage the whole simulation process. It glues three
    major components of  simuPOP together: population, mating,
    Operation. More specifically, a simuPOP<population> object manages
    several copies of population (or one of its subclasses) plus a
    'scratch population'; it has a mating object that knows how to
    generate next gen; it controls the evolution process by applying
    pre- during- and post- mating Operators during evolution. }

"; 

%feature("docstring") simuPOP::simulator::simulator "

Description:

    create a simulator

Details:

    DEVONLY{ m_curRep, gen are reference to glocal shared variables. }

Usage:

    simulator(pop, matingScheme, rep=1, grp=[])

Arguments:

    population:     a population created by population() function.
                    This population will be copied to the simulator so
                    its content will not be changed.
    mate:           a mating scheme
    rep:            number of replicates. default to 1
    grp:            grp number for each replicate. For example, we can
                    seprate all replicates into two groups and give
                    them differnt initial values before evolution.

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
    should *not* be used.

"; 

%feature("docstring") simuPOP::simulator::clone "

Description:

    simuPOP::simulator::clone

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::simulator::pop "

Description:

    the 'rep' replicate of this simulator

Usage:

    x.pop(rep)

Arguments:

    rep:            number of replicate.

Note:

    The returned reference is temporary in the sense that the refered
    population will be invalid after another round of evolution.
    Therefore, the use of this function should be limited to
    'immediate after retrival'. If you would like to get a persistent
    population, use getPopulation(rep)

"; 

%feature("docstring") simuPOP::simulator::getPopulation "

Details:

    this function returns a copy of population rep

Usage:

    x.getPopulation(rep)

Arguments:

    rep:            number of replicate.

"; 

%feature("docstring") simuPOP::simulator::setMatingScheme "

Description:

    simuPOP::simulator::setMatingScheme

Usage:

    x.setMatingScheme(matingScheme)

"; 

%feature("docstring") simuPOP::simulator::setPopulation "

Description:

    simuPOP::simulator::setPopulation

Usage:

    x.setPopulation(pop, rep)

"; 

%ignore simuPOP::simulator::curRep() const ;

%feature("docstring") simuPOP::simulator::numRep "

Description:

    simuPOP::simulator::numRep

Usage:

    x.numRep()

"; 

%feature("docstring") simuPOP::simulator::gen "

Description:

    simuPOP::simulator::gen

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

    set generation number

Usage:

    x.setGen(gen)

Arguments:

    gen:            new generation number

Note:

    this will also set shared variable gen

"; 

%feature("docstring") simuPOP::simulator::step "

Description:

    evolve one step

Usage:

    x.step(ops=[], preOps=[], postOps=[], steps=1, dryrun=False)

"; 

%feature("docstring") simuPOP::simulator::evolve "

Description:

    evolve till 'end' generation subject to given operators

Details:

    'operators' will be applied in the order of:all pre-mating
    opertorsduring-mating operators called by the mating scheme at the
    birth of each offspringall post-mating operatorsIf any pre or
    post-mating opertor fails to apply, that replicate will be
    stopped. This is exactly how terminators work.

Usage:

    x.evolve(ops, preOps=[], postOps=[], end=-1, dryrun=False)

Arguments:

    ops:            operators that will be applied at all generations.
                    Of course they might not be active at all
                    generations.
    preOps:         operators that will be applied before evolution.
                    evolve() function will *not* check if they are
                    active.
    postOps:        operators that will be applied after the last
                    generation, or after a replicate is terminated.
    end:            ending generation. Should be -1 (no ending
                    generation) or a number greater than current
                    generation number. When end=-1, simulator can only
                    be stopped by terminators.

Note:

    When end=-1, you can not specify negative generation parameters to
    operators. How would an operator know which genertion is the -1
    genertion if no ending genertion is given?

"; 

%feature("docstring") simuPOP::simulator::apply "

Description:

    apply some ops, geneartion of population does not change No mating
    is allowed.

Details:

    pre-mating oeprators are applied before post-mating operators. no
    during-mating operators are allowed.

Usage:

    x.apply(ops, dryrun=False)

Arguments:

    ops:            operators that will be applied at all generations.
                    Of course they might not be active at all
                    generations.

"; 

%feature("docstring") simuPOP::simulator::setStopIfOneRepStop "

Description:

    stop if one replicate stops or not

Usage:

    x.setStopIfOneRepStop(on=True)

Arguments:

    on:             turn on or off stopIfOneRepStop if set, the
                    simulator will stop evolution if one replicate
                    stops. This is sometimes useful.

"; 

%feature("docstring") simuPOP::simulator::stopIfOneRepStop "

Description:

    simuPOP::simulator::stopIfOneRepStop

Usage:

    x.stopIfOneRepStop()

"; 

%feature("docstring") simuPOP::simulator::setApplyOpToStoppedReps "

Description:

    apply ops even if rep stops

Usage:

    x.setApplyOpToStoppedReps(on=True)

Arguments:

    on:             turn on or off applyOpToStoppedReps flag if set,
                    the simulator will continue to apply operators to
                    all stopped repicates until all replicates are
                    marked stopped. This is sometimes useful.

"; 

%feature("docstring") simuPOP::simulator::applyOpToStoppedReps "

Description:

    simuPOP::simulator::applyOpToStoppedReps

Usage:

    x.applyOpToStoppedReps()

"; 

%feature("docstring") simuPOP::simulator::vars "

Description:

    get simulator namespace, if rep > 0 is given, return replicate rep
    namespace

Usage:

    x.vars(rep, subPop=-1)

"; 

%feature("docstring") simuPOP::simulator::saveSimulator "

Description:

    save simulator in 'text','bin' or 'xml' format

Usage:

    x.saveSimulator(filename, format=\"auto\", compress=True)

Arguments:

    filename:       save to filename
    format:         format to save. Can be one of 'text', 'bin', 'xml'
                    The default format is 'text' but the output is not
                    suppored to be read. 'bin' has smaller size and
                    should be used for large populations. 'xml' format
                    is most readable and should be used when you would
                    like to convert  simuPOP populations to other
                    formats.

"; 

%ignore simuPOP::simulator::loadSimulator(string filename, string format="auto");

%feature("docstring") simuPOP::simulator::__repr__ "

Description:

    simuPOP::simulator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::smmMutator "

Description:

    stepwise mutation model.

Details:

    Stepwise mutation model (SMM) assumes that alleles are represented
    by integer values and that a mutation either increases or
    decreases the allele value by one. For variable number tandem
    repeats loci (VNTR), the allele value is generally taken as the
    number of tandem repeats in the DNA sequence.

"; 

%feature("docstring") simuPOP::smmMutator::smmMutator "

Details:

    The stepwise mutation model (SMM) is developed for allozymes. It
    provides better description for these kinds of evolutionary
    processes.

Usage:

    smmMutator(rate=[], atLoci=[], maxAllele=0, incProb=0.5,
      output=\">\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    rate::          mutation rate
    incProb:        probability to increase allele state. Default to 1
    atLoci:         and other parameters: refer to help(mutator),
                    help(baseOperator.__init__)

"; 

%feature("docstring") simuPOP::smmMutator::~smmMutator "

Description:

    simuPOP::smmMutator::~smmMutator

Usage:

    x.~smmMutator()

"; 

%feature("docstring") simuPOP::smmMutator::mutate "

Description:

    how to mutate a single allele. this is usually the only function
    that need to be defined by the subclasses.

Usage:

    x.mutate(allele)

"; 

%feature("docstring") simuPOP::smmMutator::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::smmMutator::__repr__ "

Description:

    simuPOP::smmMutator::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::splitSubPop "

Description:

    split subpopulation

"; 

%feature("docstring") simuPOP::splitSubPop::splitSubPop "

Description:

    split a subpopulation (or whole population as subpop 0)

Usage:

    splitSubPop(which=0, sizes=[], proportions=[], subPopID=[],
      randomize=True, stage=PreMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    which:          which subpop to split (if there is no subpop
                    structure, 0 is the only subpop)
    subPop:         new subpop size, should add up to the size of
                    subpop to be splitted
    proportions:    proportions of new subpop. (use one of subPop or
                    proportions). Should add up to one.
    subPopID:       optional. new subpop IDs. If given, should have
                    the same length as subPop or proportions. SInce
                    subpop with negative id will be removed. You can
                    remove part of a subpop by setting a new negative
                    id.

"; 

%feature("docstring") simuPOP::splitSubPop::~splitSubPop "

Description:

    destructor

Usage:

    x.~splitSubPop()

"; 

%feature("docstring") simuPOP::splitSubPop::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::splitSubPop::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::splitSubPop::__repr__ "

Description:

    simuPOP::splitSubPop::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::spread "

Description:

    initialize genotype by value and then copy to all individuals

"; 

%feature("docstring") simuPOP::spread::spread "

Description:

    simuPOP::spread::spread

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::spread::__repr__ "

Description:

    simuPOP::spread::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::spread::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::stat "

Description:

    simuPOP::stat

"; 

%feature("docstring") simuPOP::stat::stat "

Description:

    create an stat

Usage:

    stat(popSize=False, numOfMale=False, numOfAffected=False,
      numOfAlleles=[], numOfAlleles_param=strDict, alleleFreq=[],
      alleleFreq_param=strDict, heteroFreq=[], expHetero=[],
      expHetero_param=strDict, homoFreq=[], genoFreq=[],
      haploFreq=intMatrix, LD=intMatrix, LD_param=strDict,
      association=intMatrix, association_param=strDict, Fst=[],
      Fst_param=strDict, relGroups=intMatrix, relLoci=[],
      rel_param=strDict, relBySubPop=False, relMethod=[],
      relMinScored=10, hasPhase=False, midValues=False, output=\"\",
      outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
      rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    popSize:        whether or not calculate population sizes. will
                    set numSubPop, subPopSize, popSize,
                    subPop[sp]['popSize']
    numOfMale:      whether or not count number of male and female,
                    will set numOfMale and numOfFemale for all
                    population/subpopulations
    numOfAffected:  whether or not count number of affected
                    individuals. Will set numOfAffected, and
                    numOfUnaffected.
    numOfAlleles:   an array of loci at which number of alleles will
                    be counted (0 is excluded). Note that number of
                    alleles will be automatically set if alleleFreq is
                    counted.
    alleleFreq:     an array of loci at which all alleles will be
                    counted.
    genoFreq:       an array of loci at which all genotype will be
                    counted each item is the locus index followed by
                    allele pairs.
    heteroFreq:     an array of loci at which the observed proportion
                    of individuausl heterozygous will be applyd for
                    each allele. Expected heterozygosity will also be
                    calculuate and put in heteroFreq[locus][0] (since
                    allele 0 is not used.)
    homoFreq:       an array of loci at which homozygosity number and
                    frequcies will be calculated
    expHetero:      an array of loci at which expected heterozygosity
                    will be calculated.
    haploFreq:      a matrix of haplotypes (allele sequence) to count
                    format: haploFreq = [ [ 0,1,2 ], [1,2] ]All
                    haplotypes on loci 012, 12 will be counted. If
                    only one haplotype is specified, the outer [] can
                    be ommited. I.e., haploFreq=[0,1] is acceptable.
    LD:             apply LD, LD' and r2, given LD=[ [locus1 locus2],
                    [ locus1 locus2 allele1 allele2], ,...] If two
                    numbers are given, D, D' and r2 overall possible
                    allele pairs will be calculated and saved as AvgLD
                    etc. If four numbers are given, D, D' and r2 using
                    specified alleles are provided. If only one item
                    is specified, the outer [] can be ignored. I.e.,
                    LD=[locus1 locus2] is acceptable.
    LD_param::      a dictionary of parameters to LD statistics. Can
                    have key stat: a list of statistics to calculate.
                    default to all. if any statistics is specified,
                    only those specified will be calculated. i.e.:
                    LD_param={'stat':['LD']} LD: True/False, shortcut
                    for 'stat':['LD'] LD_prime: True/False, shortcut
                    for 'stat':['LD_prime'] ... subPop: True/False:
                    whether or not calculate statistics for
                    subpopulations midValues: True/False: whether or
                    not keep intermediate results
    Fst:            calculate Fst. Fis and Fit will be given as a side
                    product. format Fst = [ 0, 1 ] Fst calculated at
                    locus 0 and 1. Since any allele can be used to
                    calculate Fst, Fst[0] will be the average over all
                    alleles as suggested by Weir and Cockerham.
    relGroups:      calculated pairwise relatedness between groups.
                    relGroups can be in the form of either
                    [[1,2,3],[4,5],[7,8]] (gourps of individuals) or
                    [1,3,4] (use subpopulations).
    relMethod:      method used to calculate relatedness. Can be
                    either REL_Queller or REL_Lynch.
    relLoci:        loci on which relatedness calues are calculated.
    hasPhase:       if a/b and b/a are the same genotype. default is
                    false.
    midValues:      whether or not post intermediate results. Default
                    to false. For example, Fst will need to calculate
                    allele frequency, if midValues is set to true,
                    allele frequencies will be posted as well. This
                    will help debuging and sometimes derived
                    statistics.
    others:         there is NO output for this operator. other
                    parameters please see help(baseOperator.__init__)

Note:

    previous provide 'search for all allele/genotype' option but this
    has proven to be troublesome. In this version, everything should
    be explicitly specified.

"; 

%feature("docstring") simuPOP::stat::~stat "

Description:

    simuPOP::stat::~stat

Usage:

    x.~stat()

"; 

%feature("docstring") simuPOP::stat::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::stat::apply "

Description:

    count various statistics. use m_alleles etc to save (potentially)
    time to resize all these variables.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::stat::__repr__ "

Description:

    simuPOP::stat::__repr__

Usage:

    x.__repr__()

"; 

%ignore simuPOP::statAlleleFreq;

%feature("docstring") simuPOP::statAlleleFreq::statAlleleFreq "

Description:

    simuPOP::statAlleleFreq::statAlleleFreq

Usage:

    statAlleleFreq(atLoci=[], param=strDict)

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

    statAssociation(alleleFreq, haploFreq, Association=intMatrix,
      param=strDict)

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

    statExpHetero(alleleFreq, expHetero=[], param=strDict)

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

    statFst(alleleFreq, heteroFreq, Fst=[], param=strDict)

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

Details:

    go through a single allele for all individual, all diploidneed to
    replace previous values

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::statHaploFreq;

%feature("docstring") simuPOP::statHaploFreq::statHaploFreq "

Description:

    simuPOP::statHaploFreq::statHaploFreq

Usage:

    statHaploFreq(haploFreq=intMatrix)

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

    statLD(alleleFreq, haploFreq, LD=intMatrix, LD_param=strDict)

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

    statNumOfAffected(numOfAffected=False)

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

    statNumOfAlleles(calc, atLoci=[], param=strDict)

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

    statNumOfMale(numOfMale=False)

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

    NOTE: the default output for stator is \"\", i.e., no output i.e.,
    stator will write to shared variables and unless specified by
    output=\">\" etc, no output will be generated. this class will also
    list ALL statistics and its names?

"; 

%feature("docstring") simuPOP::stator::stator "

Description:

    constructor. default to be always active. default to have NO
    output (shared variables will be set.) phase: if we treat Aa!=aA,
    default is false, i.e, Aa=aA.

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

    this function is very important

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

    statRelatedness(alleleFreq, groups=intMatrix, useSubPop=False,
      loci=[], method=[], minScored=10, param=strDict)

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

%feature("docstring") simuPOP::StreamElem::StreamElem "

Description:

    Stream element, can be of different
    types///////////////////////////////////////////////////////////.

Usage:

    StreamElem(name, readable, realAppend, useString)

"; 

%feature("docstring") simuPOP::StreamElem::StreamElem "

Description:

    copy constructor, we need to clear rhs.m_stream to avoid closing
    file too early this is techniquely advanced (and dangerous)

Usage:

    StreamElem(rhs)

"; 

%feature("docstring") simuPOP::StreamElem::~StreamElem "

Description:

    destructor

Usage:

    x.~StreamElem()

"; 

%feature("docstring") simuPOP::StreamElem::makeReadable "

Description:

    the file was write-only, re-open it as read-write

Usage:

    x.makeReadable()

"; 

%feature("docstring") simuPOP::StreamElem::makeAppend "

Description:

    change the append status

Usage:

    x.makeAppend(append)

"; 

%ignore simuPOP::StreamElem::stream();

%ignore simuPOP::StreamElem::type();

%ignore simuPOP::StreamElem::info();

%ignore simuPOP::StreamElem::append();

%ignore simuPOP::StreamProvider;

%feature("docstring") simuPOP::StreamProvider::StreamProvider "

Description:

    Stream provider///////////////////////////////////////////////////
    ////////.

Usage:

    StreamProvider(output, outputExpr)

"; 

%feature("docstring") simuPOP::StreamProvider::~StreamProvider "

Description:

    simuPOP::StreamProvider::~StreamProvider

Usage:

    x.~StreamProvider()

"; 

%ignore simuPOP::StreamProvider::setOutput(const string &output, const string &outputExpr);

%ignore simuPOP::StreamProvider::noOutput();

%ignore simuPOP::StreamProvider::getOstream(PyObject *dict=NULL, bool readable=false);

%feature("docstring") simuPOP::StreamProvider::closeOstream "

Description:

    close ostream and delete ostream pointer.. if it is a ofstream.

Usage:

    x.closeOstream()

"; 

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

    tagger is a during mating operator that tag individual with
    various information. Potential usages are 1. record parenting
    information to track pedigree. 2. tag a individual/allele and
    monitor its spread in the population etc. 3...

Details:

    Bo Peng

"; 

%feature("docstring") simuPOP::tagger::tagger "

Description:

    constructor. default to be always active but no output.

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::terminateIf "

Description:

    terminate according to a condition which can be, e.g.
    any(alleleNum0) == 0 all(alleleNum1) > 0.5 alleleNum0{2} == 0 etc.
    When the condition is true, a shared variable var=\"terminate\" will
    be set to current generation.

"; 

%feature("docstring") simuPOP::terminateIf::terminateIf "

Description:

    simuPOP::terminateIf::terminateIf

Usage:

    terminateIf(condition=\"\", message=\"\", var=\"terminate\",
      output=\"\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::terminateIf::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::terminateIf::__repr__ "

Description:

    simuPOP::terminateIf::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::terminateIf::apply "

Description:

    check all alleles in vector allele if they are fixed.

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

    simuPOP::terminator

"; 

%feature("docstring") simuPOP::terminator::terminator "

Description:

    constructor. default to be always active.

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::terminator::message "

Description:

    simuPOP::terminator::message

Usage:

    x.message()

"; 

%feature("docstring") simuPOP::ticToc "

Description:

    simuPOP::ticToc

"; 

%feature("docstring") simuPOP::ticToc::ticToc "

Description:

    timer if called, output time passed since last calling time.

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

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::ticToc::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::ticToc::__repr__ "

Description:

    simuPOP::ticToc::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::turnOffDebugOp "

Description:

    simuPOP::turnOffDebugOp

"; 

%feature("docstring") simuPOP::turnOffDebugOp::turnOffDebugOp "

Description:

    turn on debug

Usage:

    turnOffDebugOp(code, stage=PreMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::turnOffDebugOp::~turnOffDebugOp "

Description:

    destructor

Usage:

    x.~turnOffDebugOp()

"; 

%feature("docstring") simuPOP::turnOffDebugOp::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::turnOffDebugOp::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::turnOffDebugOp::__repr__ "

Description:

    simuPOP::turnOffDebugOp::__repr__

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::turnOnDebugOp "

Description:

    simuPOP::turnOnDebugOp

"; 

%feature("docstring") simuPOP::turnOnDebugOp::turnOnDebugOp "

Description:

    turn on debug

Usage:

    turnOnDebugOp(code, stage=PreMating, begin=0, end=-1, step=1,
      at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

"; 

%feature("docstring") simuPOP::turnOnDebugOp::~turnOnDebugOp "

Description:

    destructor

Usage:

    x.~turnOnDebugOp()

"; 

%feature("docstring") simuPOP::turnOnDebugOp::clone "

Description:

    this function is very important

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::turnOnDebugOp::apply "

Description:

    apply to one population, does not check if the oeprator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::turnOnDebugOp::__repr__ "

Description:

    simuPOP::turnOnDebugOp::__repr__

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

    FIXME: consider adopting R's implementation. They may be quicker.

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

%ignore simuPOP::countAlleles(population &pop, int subpop, const vectori &loci, const vectori &alleles, vectorlu &alleleNum);

%ignore simuPOP::getExpectedAlleles(population &pop, vectorf &expFreq, const vectori &loci, const vectori &alleles, vectoru &expAlleles);

%feature("docstring") simuPOP::FreqTrajectoryStoch "

Description:

    simuPOP::FreqTrajectoryStoch

Usage:

    FreqTrajectoryStoch(curGen, freq, N, *NtFunc, fitness,
      *fitnessFunc, minMutAge, maxMutAge, ploidy, restartIfFail,
      maxAttempts, allowFixation)

"; 

%ignore simuPOP::fitOfGeno(unsigned loc, const vectori &allgeno, const vectorf &fitness, const vectorf::const_iterator &freq);

%ignore simuPOP::interFitness(unsigned nLoci, const vectorf &fitness, const vectorf::const_iterator &freq, vectorf &sAll);

%ignore simuPOP::MarginalFitness(unsigned nLoci, const vectorf &fitness, const vectorf &freq);

%feature("docstring") simuPOP::FreqTrajectoryMultiStoch "

Description:

    simuPOP::FreqTrajectoryMultiStoch

Usage:

    FreqTrajectoryMultiStoch(curGen, freq, N, *NtFunc, fitness,
      *fitnessFunc, minMutAge, maxMutAge, ploidy, restartIfFail,
      maxAttempts)

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

    simuPOP::LoadPopulation

Usage:

    LoadPopulation(file, format)

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

    simuPOP::LoadSimulator

Usage:

    LoadSimulator(file, mate, format)

"; 

%ignore simuPOP::haploKey(const vectori &seq);

%feature("docstring") simuPOP::TurnOnDebug "

Description:

    set debug code, default to turn all code on

Usage:

    TurnOnDebug(code)

"; 

%feature("docstring") simuPOP::TurnOnDebugWithName "

Description:

    set debug code, using name

Usage:

    TurnOnDebugWithName(code)

"; 

%feature("docstring") simuPOP::TurnOffDebug "

Description:

    turn off debug, default to turn all code off

Usage:

    TurnOffDebug(code)

"; 

%ignore simuPOP::debug(DBG_CODE code);

%feature("docstring") simuPOP::ListDebugCode "

Description:

    show all dbg codes (print to cout)

Usage:

    ListDebugCode()

"; 

%ignore simuPOP::dbgString(DBG_CODE code);

%ignore simuPOP::simuPOP_kbhit(void);

%ignore simuPOP::simuPOP_getch(void);

%ignore simuPOP::PyObj_As_Bool(PyObject *obj, bool &val);

%ignore simuPOP::PyObj_As_Int(PyObject *obj, int &val);

%ignore simuPOP::PyObj_As_Double(PyObject *obj, double &val);

%ignore simuPOP::PyObj_As_String(PyObject *obj, string &val);

%ignore simuPOP::PyObj_As_StrDict(PyObject *obj, strDict &val);

%ignore simuPOP::PyObj_As_Array(PyObject *obj, vectorf &val);

%ignore simuPOP::PyObj_As_IntArray(PyObject *obj, vectori &val);

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

%ignore simuPOP::save_none(string &str);

%ignore simuPOP::load_none(const string &str, size_t &offset);

%ignore simuPOP::save_int(string &str, PyObject *args);

%ignore simuPOP::load_int(const string &str, size_t &offset);

%ignore simuPOP::save_long(string &str, PyObject *args);

%ignore simuPOP::load_long(const string &str, size_t &offset);

%ignore simuPOP::save_float(string &str, PyObject *args);

%ignore simuPOP::load_float(const string &str, size_t &offset);

%ignore simuPOP::save_string(string &str, PyObject *args);

%ignore simuPOP::load_string(const string &str, size_t &offset);

%ignore simuPOP::saveObj(string &str, PyObject *args);

%ignore simuPOP::loadObj(const string &vars, size_t &offset);

%ignore simuPOP::save_dict(string &str, PyObject *args);

%ignore simuPOP::load_dict(const string &vars, size_t &offset);

%ignore simuPOP::save_list(string &str, PyObject *args);

%ignore simuPOP::load_list(const string &vars, size_t &offset);

%ignore simuPOP::save_tuple(string &str, PyObject *args);

%ignore simuPOP::load_tuple(const string &vars, size_t &offset);

%ignore simuPOP::mainVars();

%ignore simuPOP::moduleVars();

%ignore simuPOP::pyPopObj(void *p);

%ignore simuPOP::pyIndObj(void *p);

%ignore simuPOP::ostreamManager();

%feature("docstring") simuPOP::rng "

Description:

    currently, return a global  RNG.

Usage:

    rng()

"; 

%feature("docstring") simuPOP::SetRNG "

Description:

    set random number generator

Usage:

    SetRNG(r, seed)

"; 

%feature("docstring") simuPOP::setRNG "

Description:

    for backward compatibilit, will remove later

Usage:

    setRNG(r, seed)

"; 

%feature("docstring") simuPOP::ListAllRNG "

Description:

    list all available  RNG.

Usage:

    ListAllRNG()

"; 

%feature("docstring") simuPOP::listAllRNG "

Description:

    for backward compatibility, will remove later

Usage:

    listAllRNG()

"; 

%ignore simuPOP::gsl_error_handler(const char *reason, const char *, int, int gsl_errno);

%feature("docstring") simuPOP::g_cnull "

Description:

    null stream

Usage:

    g_cnull(g_nullStreamBuf)

"; 

%ignore simuPOP::cnull();

%feature("docstring") simuPOP::setLogOutput "

Description:

    set standard output to (default standard Python output)

Usage:

    setLogOutput(filename)

"; 

%feature("docstring") simuPOP::simuRev "

Description:

    Global debug and initialization related functions/////////////////
    ////////////////////////////////////////// return svn revision.

Usage:

    simuRev()

"; 

%feature("docstring") simuPOP::simuVer "

Description:

    return version infomation of simuPOP

Usage:

    simuVer()

"; 

%feature("docstring") simuPOP::optimized "

Description:

    simuPOP::optimized

Usage:

    optimized()

"; 

%feature("docstring") simuPOP::mpi "

Description:

    simuPOP::mpi

Usage:

    mpi()

"; 

%feature("docstring") simuPOP::mpiRank "

Description:

    simuPOP::mpiRank

Usage:

    mpiRank()

"; 

%feature("docstring") simuPOP::mpiSize "

Description:

    simuPOP::mpiSize

Usage:

    mpiSize()

"; 

%ignore simuPOP::mpiBarrier();

%feature("docstring") simuPOP::alleleType "

Description:

    simuPOP::alleleType

Usage:

    alleleType()

"; 

%feature("docstring") simuPOP::compileCompiler "

Description:

    simuPOP::compileCompiler

Usage:

    compileCompiler()

"; 

%feature("docstring") simuPOP::compileDate "

Description:

    simuPOP::compileDate

Usage:

    compileDate()

"; 

%feature("docstring") simuPOP::compilePyVersion "

Description:

    simuPOP::compilePyVersion

Usage:

    compilePyVersion()

"; 

%feature("docstring") simuPOP::compilePlatForm "

Description:

    simuPOP::compilePlatForm

Usage:

    compilePlatForm()

"; 

%ignore simuPOP::isGzipped(const string &filename);

%ignore simuPOP::fileExtension(const string &filename);

%ignore simuPOP::initialize();

%feature("docstring") simuPOP::testGappedIterator "

Description:

    simuPOP::testGappedIterator

Usage:

    testGappedIterator()

"; 

%ignore std::pow3(unsigned n);

%feature("docstring") simuPOP::newcarrayobject "

Description:

    simuPOP::newcarrayobject

Usage:

    newcarrayobject(*buf, type, size)

"; 

%ignore simuPOP::newcarrayiterobject(GenoIterator begin, GenoIterator end);

%ignore simuPOP::is_carrayobject(PyObject *);

%ignore simuPOP::carray_length(PyObject *a);

%ignore simuPOP::carray_itemsize(PyObject *a);

%ignore simuPOP::carray_type(PyObject *a);

%ignore simuPOP::carray_data(PyObject *a);

%ignore simuPOP::initcarray(void);

%ignore simuPOP::pow3(unsigned n);

