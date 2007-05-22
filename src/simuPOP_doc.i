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

Usage:

    x.drawsample(pop)

Details:

    collect all families

"; 

%feature("docstring") simuPOP::affectedSibpairSample::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

Description:

    a mating scheme that uses binomial selection, regardless of sex

Details:

    No sex information is involved (binomial random selection).
    Offspring is chosen from parental generation by random or
    according to the fitness values. In this mating scheme,
    * numOffspring  protocol is honored;
    * population size changes are allowed;
    * selection is possible;
    * haploid populaton is allowed.

"; 

%feature("docstring") simuPOP::binomialSelection::binomialSelection "

Description:

    create a binomial selection mating scheme

Usage:

    binomialSelection(numOffspring=1., *numOffspringFunc=NULL,
      maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
      newSubPopSizeExpr=\"\", *newSubPopSizeFunc=NULL)

Details:

    Please refer to class mating  for parameter descriptions.

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

%ignore simuPOP::binomialSelection::submitScratch(population &pop, population &scratch);

%ignore simuPOP::binomialSelection::mate(population &pop, population &scratch, vector< Operator * > &ops, bool submit);

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

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::continueIf::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

Description:

    a controlled binomial random selection mating scheme

Details:

    This is the controlled binomial random selection mating scheme
    described in  Peng 2007 (PLoS Genetics)  . Basically, a freqFunc
    is passed to this mating scheme and set the allele frequencies of
    given alleles at given loci at the offspring generation.
    The offspring generation is conceptually populated in two steps.
    At the first step, only families with disease alleles are accepted
    until the expected number of disease alleles are met. At the
    second step, only families with wide type alleles are accepted to
    populate the rest of the offspring generation.

"; 

%feature("docstring") simuPOP::controlledBinomialSelection::controlledBinomialSelection "

Description:

    create a controlled binomial random selection mating scheme

Usage:

    controlledBinomialSelection(loci, alleles, *freqFunc,
      numOffspring=1., *numOffspringFunc=NULL, maxNumOffspring=0,
      mode=MATE_NumOffspring, newSubPopSize=[], newSubPopSizeExpr=\"\",
      *newSubPopSizeFunc=NULL)

Details:

    Please refer to class mating  for descriptions of parameters.

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

    deep copy of a controlled binomial random selection mating scheme

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::controlledBinomialSelection::__repr__ "

Description:

    used by Python print function to print out the general information
    of the controlled binomial random selection mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::controlledBinomialSelection::submitScratch(population &pop, population &scratch);

%ignore simuPOP::controlledBinomialSelection::mate(population &pop, population &scratch, vector< Operator * > &ops, bool submit);

%feature("docstring") simuPOP::controlledMating "

Description:

    a controlled mating scheme

Details:

    This is an experimental mating scheme that uses a frequency range
    to control the allele frequency of the offspring generation at
    given loci. When allele frequencies at the offspring generation
    does not fall into the given range, the offspring generation is
    regenerated. Any mating scheme can be used with this mating scheme
    by passing through parameter matingScheme .

"; 

%feature("docstring") simuPOP::controlledMating::controlledMating "

Description:

    control allele frequencies at a locus

Usage:

    controlledMating(matingScheme, loci, alleles, *freqFunc,
      range=0.01)

Arguments:

    matingScheme:   a mating scheme
    loci:           loci at which allele frequency is controlled. Note
                    that controlling the allele frequencies at several
                    loci may take a long time.
    alleles:        alleles to control at each locus. Should have the
                    same length as loci .
    freqFunc:       frequency boundaries. If the length of the
                    returned value equals the size of loci , the range
                    for loci will be [value0, value0+range] , [value1,
                    value1+range]  etc. If the length of the returned
                    value is 2 times the size of loci , it will be
                    interpreted as [low1, high1, low2, high2, ...] .

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

%ignore simuPOP::controlledMating::mate(population &pop, population &scratch, vector< Operator * > &ops, bool submit);

%feature("docstring") simuPOP::controlledRandomMating "

Description:

    a controlled random mating scheme

Details:

    This is the controlled random mating scheme described in  Peng
    2007 (PLoS Genetics)  . Basically, a freqFunc  is passed to this
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

    controlledRandomMating(loci, alleles, *freqFunc, acceptScheme=0,
      numOffspring=1., *numOffspringFunc=NULL, maxNumOffspring=0,
      mode=MATE_NumOffspring, newSubPopSize=[],
      *newSubPopSizeFunc=NULL, newSubPopSizeExpr=\"\",
      contWhenUniSex=True)

Arguments:

    loci:           loci at which allele frequencies are monitored
                    (controlled)
    alleles:        alleles at given loci. It should have the same
                    length as loci
    freqFunc:       a Python function that accepts a generation number
                    and returns expected allele frequencies at given
                    loci
    acceptScheme:   internal use only

Details:

    Please refer to class mating  for descriptions of other
    parameters.

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

%ignore simuPOP::controlledRandomMating::mate(population &pop, population &scratch, vector< Operator * > &ops, bool submit);

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

Usage:

    x.apply(pop)

Details:

    dump population structuredump all genotypic info

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

    genotypic structure related functions, can be accessed from both
    individuals and populations

Details:

    Genotypic structure refers to the number of chromosomes,
    positions, the number of loci on each chromosome, and allele and
    locus names etc. All individuals in a population share the same
    genotypic structure and functions provided in this class can be
    accessed from individual, population and simulator levels.

"; 

%feature("docstring") simuPOP::GenoStruTrait::GenoStruTrait "

Description:

    Creat a  GenoStruTrait  class, but m_genoStruIdx  will be set
    later.

Usage:

    GenoStruTrait()

"; 

%ignore simuPOP::GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru &loci, bool sexChrom, const vectorf &lociPos, const vectorstr &alleleNames, const vectorstr &lociNames, UINT maxAllele, const vectorstr &infoFields, const vectori &chromMap);

%ignore simuPOP::GenoStruTrait::setGenoStructure(GenoStructure &rhs);

%ignore simuPOP::GenoStruTrait::setGenoStruIdx(size_t idx);

%ignore simuPOP::GenoStruTrait::mergeGenoStru(size_t idx) const ;

%ignore simuPOP::GenoStruTrait::removeLociFromGenoStru(const vectoru &remove=vectoru(), const vectoru &keep=vectoru());

%ignore simuPOP::GenoStruTrait::insertBeforeLociToGenoStru(const vectoru &idx, const vectorf &pos, const vectorstr &names) const ;

%ignore simuPOP::GenoStruTrait::insertAfterLociToGenoStru(const vectoru &idx, const vectorf &pos, const vectorstr &names) const ;

%ignore simuPOP::GenoStruTrait::genoStru() const ;

%ignore simuPOP::GenoStruTrait::genoStruIdx() const ;

%feature("docstring") simuPOP::GenoStruTrait::ploidy "

Description:

    return ploidy, the number of homologous sets of chromosomes

Usage:

    x.ploidy()

"; 

%feature("docstring") simuPOP::GenoStruTrait::ploidyName "

Description:

    return ploidy name, haploid , diploid , or triploid  etc.

Usage:

    x.ploidyName()

"; 

%feature("docstring") simuPOP::GenoStruTrait::numLoci "

Description:

    return the number of loci on chromosome chrom , equals to
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

    return an array of loci positions on a given chromosome

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

%feature("docstring") simuPOP::GenoStruTrait::chromIndex "

Description:

    return an array of chromosome indices

Usage:

    x.chromIndex()

"; 

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

"; 

%feature("docstring") simuPOP::GenoStruTrait::absLocusIndex "

Description:

    return the absolute index of a locus on a chromosome

Usage:

    x.absLocusIndex(chrom, locus)

"; 

%feature("docstring") simuPOP::GenoStruTrait::chromLocusPair "

Description:

    return (chrom, locus)  pair of an absolute locus index

Usage:

    x.chromLocusPair(locus)

"; 

%feature("docstring") simuPOP::GenoStruTrait::alleleName "

Description:

    return the name of an allele (if previously specified)

Usage:

    x.alleleName(allele)

"; 

%feature("docstring") simuPOP::GenoStruTrait::alleleNames "

Description:

    return an array of allelic names

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

    return the maximum allele value for all loci

Usage:

    x.maxAllele()

"; 

%feature("docstring") simuPOP::GenoStruTrait::setMaxAllele "

Description:

    set the maximum allele value for all loci

Usage:

    x.setMaxAllele(maxAllele)

Details:

    Maximum allele value has to be 1  for binary modules. maxAllele
    is the maximum possible allele value, which allows maxAllele+1
    alleles 0, 1, ..., maxAllele .

"; 

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

    return the index of the field name , return -1  if not found

Usage:

    x.infoIdx(name)

"; 

%ignore simuPOP::GenoStruTrait::struAddInfoFields(const vectorstr &fields);

%ignore simuPOP::GenoStruTrait::struSetInfoFields(const vectorstr &fields);

%ignore simuPOP::GenoStruTrait::swap(GenoStruTrait &rhs);

%feature("docstring") simuPOP::GenoStruTrait::chromMap "

Description:

    return the distribution of chromosomes across multiple nodes (MPI
    version of  simuPOP only)

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

Details:

    The generalized stepwise mutation model (GMM) is developed for
    allozymes. It provides better description for these kinds of
    evolutionary processes.

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

    used by Python print function to print out the general information
    of the operator

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

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::ifElse::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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
    the following individual information:
    * shared genotypic structure information
    * genotype
    * sex, affection status, subpopulation ID
    * optional information fields Individual genotypes are arranged by
    locus, chromosome, ploidy, in that order, and can be accessed from
    a single index. For example, for a diploid individual with two
    loci on the first chromosome, one locus on the second, its
    genotype is arranged as  1-1-1 1-1-2 1-2-1 2-1-1 2-1-2 2-2-1
    where x-y-z  represents ploidy x  chromosome y  and locus z . An
    allele 2-1-2  can be accessed by allele(4)  (by absolute index),
    allele(2, 1)  (by index and ploidy) or allele(1, 1, 0)  (by index,
    ploidy and chromosome).

"; 

%feature("docstring") simuPOP::individual::individual "

Description:

    Individuals are created by populations automatically. Do not call
    this function directly.

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

%ignore simuPOP::individual::copyFrom(const individual &rhs);

%ignore simuPOP::individual::genoPtr() const ;

%ignore simuPOP::individual::infoPtr() const ;

%feature("docstring") simuPOP::individual::arrGenotype "

Description:

    return an editable array (a Python list of length
    totNumLoci()*ploidy() ) of genotypes of an individual

Usage:

    x.arrGenotype()

Details:

    This function returns the whole genotype. Although this function
    is not as easy to use as other functions that access alleles, it
    is the fastest one since you can read/write genotype directly.

"; 

%feature("docstring") simuPOP::individual::arrGenotype "

Description:

    return only the p-th  copy of the chromosomes

Usage:

    x.arrGenotype(p)

"; 

%feature("docstring") simuPOP::individual::arrGenotype "

Description:

    return only the ch-th  chromosome of the p-th  copy

Usage:

    x.arrGenotype(p, ch)

"; 

%feature("docstring") simuPOP::individual::arrInfo "

Description:

    return an editable array of all information fields (a Python list
    of length infosSize() )

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
                    ranging from 0  to  totNumLoci()*ploidy()

"; 

%feature("docstring") simuPOP::individual::allele "

Description:

    return the allele at locus index  of the p-th  copy of the
    chromosomes

Usage:

    x.allele(index, p)

Arguments:

    index:          index from the begining of the p-th  set of the
                    chromosomes, ranging from 0  to  totNumLoci()
    p:              index of the ploidy

"; 

%feature("docstring") simuPOP::individual::allele "

Description:

    return the allele at locus index  of the ch-th  chromosome of the
    p-th  chromosome set

Usage:

    x.allele(index, p, ch)

Arguments:

    index:          index from the begining of chromosome ch  of
                    ploidy p , ranging from 0  to  numLoci(ch)
    p:              index of the polidy
    ch:             index of the chromosome in the p-th  chromosome
                    set

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
                    0  to  totNumLoci()*ploidy()

"; 

%feature("docstring") simuPOP::individual::setAllele "

Description:

    set the allele at locus index  of the p-th  copy of the
    chromosomes

Usage:

    x.setAllele(allele, index, p)

Arguments:

    allele:         allele to be set
    index:          index from the begining of the poloidy p , ranging
                    from 0  to  totNumLoci(p)
    p:              index of the poloidy

"; 

%feature("docstring") simuPOP::individual::setAllele "

Description:

    set the allele at locus index  of the ch-th  chromosome in the
    p-th  chromosome set

Usage:

    x.setAllele(allele, index, p, ch)

Arguments:

    allele:         allele to be set
    index:          index from the begining of the chromosome, ranging
                    from 0  to numLoci(ch)
    p:              index of the ploidy
    ch:             index of the chromosome in ploidy p

"; 

%feature("docstring") simuPOP::individual::sex "

Description:

    return the sex of an individual, 1  for males and 2  for females.
    However, this is not guranteed so please use  sexChar() .

Usage:

    x.sex()

"; 

%feature("docstring") simuPOP::individual::sexChar "

Description:

    return the sex of an individual, M  or F

Usage:

    x.sexChar()

"; 

%feature("docstring") simuPOP::individual::setSex "

Description:

    set the sex. You should use setSex(Male)  or setSex(Female)
    instead of 1  and 2 .

Usage:

    x.setSex(sex)

"; 

%feature("docstring") simuPOP::individual::affected "

Description:

    whether or not an individual is affected

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

    return A  or U  for affection status

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

    return the ID of the subpopulation to which this individual blongs

Usage:

    x.subPopID()

Note:

    subPopID  is not set by default. It only corresponds to the
    subpopulation in which this individual resides after
    pop::setIndSubPopID  is called.

"; 

%feature("docstring") simuPOP::individual::setSubPopID "

Description:

    set new subpopulation ID, pop.rearrangeByIndID  will move this
    individual to that population

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

Arguments:

    name:           name of the information field

Details:

    Equivalent to info(infoIdx(name)) .

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

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::inheritTagger::applyDuringMating(population &pop, population::IndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::inheritTagger::clone "

Description:

    deep copy of an operator

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

Details:

    This operator randomly assign alleles according to given allele
    frequency. Allele frequencies can differ by subpop. Sex is also
    assigned randomly.

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

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::initByFreq::apply "

Usage:

    x.apply(pop)

Details:

    initialize m_ranges

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

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::initByValue::apply "

Usage:

    x.apply(pop)

Details:

    fixme: check length of src?atLoci is in effect

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

    used by Python print function to print out the general information
    of the operator

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

    used by Python print function to print out the general information
    of the operator

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

Usage:

    x.drawsample(pop)

Details:

    collect all families

"; 

%feature("docstring") simuPOP::largePedigreeSample::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

Usage:

    x.penet(*ind)

Details:

    get genotype of ind

"; 

%feature("docstring") simuPOP::maPenetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

Usage:

    x.penet(*ind)

Details:

    get genotype of ind

"; 

%feature("docstring") simuPOP::mapPenetrance::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mapQuanTrait::qtrait "

Description:

    currently assuming diploid

Usage:

    x.qtrait(*ind)

Details:

    get genotype of ind

"; 

%feature("docstring") simuPOP::mapQuanTrait::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mapSelector::indFitness "

Description:

    currently assuming diploid

Usage:

    x.indFitness(*ind, gen)

Details:

    get genotype of ind

"; 

%feature("docstring") simuPOP::mapSelector::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::maQuanTrait::qtrait "

Description:

    currently assuming diploid

Usage:

    x.qtrait(*ind)

Details:

    get genotype of ind

"; 

%feature("docstring") simuPOP::maQuanTrait::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::maSelector::indFitness "

Description:

    currently assuming diploid

Usage:

    x.indFitness(*ind, gen)

Details:

    get genotype of ind

"; 

%feature("docstring") simuPOP::maSelector::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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
    * change population/subpopulation sizes;
    * randomly select parent(s) to generate offspring to fill the next
    generation;
    * during-mating  operators are applied to all offspring;
    * apply selection if applicable.

"; 

%ignore simuPOP::mating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::mating::mating "

Description:

    create a mating scheme

Usage:

    mating(numOffspring=1.0, *numOffspringFunc=NULL,
      maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
      newSubPopSizeExpr=\"\", *newSubPopSizeFunc=NULL)

Arguments:

    numOffspring:   number of offspring or p  for a random
                    distribution. Default to 1. This parameter
                    determines the number of offspring that a mating
                    event will produce. Therefore, it determines the
                    family size.
    numOffspringFunc:a python function that returns the number of
                    offspring or p
    maxNumOffspring:used when numOffspring  is generated from a
                    binomial distribution
    mode:           can be one of MATE_NumOffspring,
                    MATE_NumOffspringEachFamily,
                    MATE_GeometricDistribution,
                    MATE_PoissonDistribution,
                    MATE_BinomialDistribution,
                    MATE_UniformDistribution .
    newSubPopSize:  an array of subpopulaitons sizes
    newSubPopSizeExpr:an expression that will return the new
                    subpopulation size
    newSubPopSizeFunc:a function that accepts an int
                    parameter(generation), an array of current
                    population size and return an array of
                    subpopulation sizes. This is usually easier to use
                    than its expression version of this parameter.

Details:

    By default, a mating scheme keeps a constant population size,
    generate one offspring per mating event. These can be changed
    using a variety of parameters. First, newSubPopSize ,
    newSubPopSizeExpr  and newSubPopSizeFunc  can be used to specify
    subpopulation sizes of the offspring generation. mode ,
    numOffspring , maxNumOffspring  can be used to specify how many
    offspring will be produced for each mating event. This mode
    parameter ...

Please refer to the reference manual for more details.

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

%ignore simuPOP::mating::submitScratch(population &pop, population &scratch);

%ignore simuPOP::mating::mate(population &pop, population &scratch, vector< Operator * > &ops, bool submit);

%ignore simuPOP::mating::fixedFamilySize();

%ignore simuPOP::mating::numOffspring(int gen);

%ignore simuPOP::mating::resetNumOffspring();

%ignore simuPOP::mating::prepareScratchPop(population &pop, population &scratch);

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mergeSubPops::apply "

Description:

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::mergeSubPops::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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
                    help(baseOperator.__init__)

Details:

    rate is a matrix with dimensions determined by fromSubPop and
    toSubPop. By default, rate is a matrix with element (i,j) being
    the migration rate, probability or count from subpop i to subpop
    j. If fromSubPop and/or toSubPop are given, migration only happen
    between these subpopulations. An extreme case is 'point
    migration'rate=[[r]], fromSubPop=a, toSubPop=bwhich migrate from
    subpop a to b with given rate r.

"; 

%feature("docstring") simuPOP::migrator::~migrator "

Description:

    destructor

Usage:

    x.~migrator()

"; 

%feature("docstring") simuPOP::migrator::clone "

Description:

    deep copy of an operator

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

Usage:

    x.setRates(rate, mode)

Details:

    format 0-0 0-1 0-2, 1-0 1-1 1-2, 2-0, 2-1, 2-2. for mode 1 or 2,
    00,11,22 will be set automatically. regardless of input.

"; 

%feature("docstring") simuPOP::migrator::apply "

Usage:

    x.apply(pop)

Details:

    2nd, or 3rd methodcreate a vector and assign indices, then random
    shuffle and assign infofor all subPop.

"; 

%feature("docstring") simuPOP::migrator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

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

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::mlSelector::indFitness "

Description:

    currently assuming diploid

Usage:

    x.indFitness(*ind, gen)

Details:

    fixme

"; 

%feature("docstring") simuPOP::mlSelector::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

Description:

    a mating scheme that does nothing

Details:

    In this scheme, there is
    * no mating. Parent generation will be considered as offspring
    generation.
    * no subpopulation change. During-mating  operators will be
    applied, but the return values are not checked. I.e., subpopsizes
    will be ignored although some during-mating operators may be
    applied.

"; 

%feature("docstring") simuPOP::noMating::noMating "

Description:

    creat a scheme with no mating

Usage:

    noMating()

Note:

    There is no new subPopsize  parameter.

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

%ignore simuPOP::noMating::mate(population &pop, population &scratch, vector< Operator * > &ops, bool submit);

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

    providing interface to apply operator before during or after
    mating

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

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::noneOp::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

Usage:

    x.drawsample(pop)

Details:

    collect all families

"; 

%feature("docstring") simuPOP::nuclearFamilySample::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

%ignore simuPOP::offspringGenerator;

%feature("docstring") simuPOP::offspringGenerator::offspringGenerator "

Description:

    create an offspring generator, save information from pop  and ops
    to speed up the calls to generateOffspring

Usage:

    offspringGenerator(pop, ops)

"; 

%feature("docstring") simuPOP::offspringGenerator::generateOffspring "

Description:

    generate numOff  offspring, or until reach offEnd

Usage:

    x.generateOffspring(pop, *dad, *mom, numOff, offBegin)

Details:

    This is because offBegin+numOff  may go beyond the subpopulation
    boundary.

"; 

%feature("docstring") simuPOP::offspringGenerator::copyOffspring "

Description:

    copy numOff  offspring, or until reach offEnd

Usage:

    x.copyOffspring(pop, *par, numOff, offBegin)

Details:

    This is because offBegin+numOff  may go beyond the subpopulation
    boundary.

"; 

%feature("docstring") simuPOP::Operator "

Description:

    base class of all classes that manipulate populations

Details:

    Operators are objects that act on populations. They can be applied
    to populations directly using their function forms, but they are
    usually managed and applied by a simulator.
    Operators can be applied at different stages of the life cycle of
    a generation. More specifically, they can be applied at pre- ,
    during- , post-mating , or a combination of these stages.
    Applicable stages are usually set by default but you can change it
    by setting stage=(PreMating|PostMating|DuringMating|PrePostMating)
    parameter. Some operators ignore stage  parameter because they
    only work at one stage.
    Operators do not have to be applied at all generations. You can
    specify starting/ending generation, gaps between applicable
    generations, or even specific generations. For example, you might
    want to start applying migrations after certain burn-in
    generations, or calculate certain statistics only sparsely.
    Operators can have outputs. Output can be standard (terminal) or a
    file, which can vary with replicates and/or generations. Outputs
    from different operators can be accumulated to the same file to
    form table-like outputs.
    Operators are applied to every replicate of a simulator by
    default. However, you can apply operators to one or a group of
    replicates using parameter rep  or grp .
    Filenames can have the following format:
    * 'filename'  this file will be overwritten each time. If two
    operators output to the same file, only the last one will succeed;
    * '>filename'  the same as 'filename' ;
    * '>>filename'  the file will be created at the beginning of
    evolution ( simulator::evolve ) and closed at the end. Output from
    several operators is allowed;
    * '>>>filename'  the same as '>>filename'  except that the file
    will not be cleared at the beginning of evolution if it i...

Please refer to the reference manual for more details.

"; 

%feature("docstring") simuPOP::Operator::Operator "

Description:

    common interface for all operators (this base operator does
    nothing by itself.)

Usage:

    Operator(output, outputExpr, stage, begin, end, step, at, rep,
      grp, infoFields)

Arguments:

    begin:          the starting generation. Default to 0 . Negative
                    numbers are allowed.
    end:            stop applying after this generation. Negative
                    numbers are allowed.
    step:           the number of generations between active
                    generations. Default to 1 .
    at:             an array of active generations. If given, stage ,
                    begin , end , and step  will be ignored.
    rep:            applicable replicates. It can be a valid replicate
                    number, REP_ALL  (all replicates, default), or
                    REP_LAST  (only the last replicate). REP_LAST  is
                    useful in adding newlines to a table output.
    grp:            applicable group. Default to GRP_ALL . A group
                    number for each replicate is set by
                    simulator.__init__  or  simulator::setGroup() .
    output:         a string of the output filename. Different
                    operators will have different default output
                    (most commonly '>'  or '' ).
    outputExpr:     an expression that determines the output filename
                    dynamically. This expression will be evaluated
                    against a population's local namespace each time
                    when an output filename is required. For example,
                    \"'>>out%s_%s.xml' % (gen, rep)\"   will output to
                    >>>out1_1.xml   for replicate 1  at generation 1 .

Note:

    Negative generation numbers are allowed for begin , end  and at .
    They are intepretted as endGen + gen + 1 . For example, begin = -2
    in simu.evolve(..., end=20)  starts at generation 19 .

"; 

%feature("docstring") simuPOP::Operator::~Operator "

Description:

    destroy an operator

Usage:

    x.~Operator()

"; 

%feature("docstring") simuPOP::Operator::clone "

Description:

    deep copy of an operator

Usage:

    x.clone()

"; 

%ignore simuPOP::Operator::isActive(UINT rep, UINT numRep, long gen, long end, int grp, bool repOnly=false);

%feature("docstring") simuPOP::Operator::applicableGroup "

Description:

    return applicable group

Usage:

    x.applicableGroup()

"; 

%feature("docstring") simuPOP::Operator::setApplicableGroup "

Description:

    set applicable group

Usage:

    x.setApplicableGroup(grp=GRP_ALL)

Details:

    Default to GRP_ALL  (applicable to all groups). Otherwise, the
    operator is applicable to only one  group of replicates. Groups
    can be set in  simulator::setGroup() .

"; 

%feature("docstring") simuPOP::Operator::applicableReplicate "

Description:

    return applicable replicate

Usage:

    x.applicableReplicate()

"; 

%feature("docstring") simuPOP::Operator::setApplicableReplicate "

Description:

    set applicable replicate

Usage:

    x.setApplicableReplicate(rep)

"; 

%feature("docstring") simuPOP::Operator::setActiveGenerations "

Description:

    set applicable generation parameters: begin , end , step  and at

Usage:

    x.setActiveGenerations(begin=0, end=-1, step=1, at=[])

"; 

%feature("docstring") simuPOP::Operator::setApplicableStage "

Description:

    set applicable stage. Another way to set stage  parameter.

Usage:

    x.setApplicableStage(stage)

"; 

%feature("docstring") simuPOP::Operator::canApplyPreMating "

Description:

    set if this operator can be applied pre-mating

Usage:

    x.canApplyPreMating()

"; 

%feature("docstring") simuPOP::Operator::canApplyDuringMating "

Description:

    set if this operator can be applied during-mating

Usage:

    x.canApplyDuringMating()

"; 

%feature("docstring") simuPOP::Operator::canApplyPostMating "

Description:

    set if this operator can be applied post-mating

Usage:

    x.canApplyPostMating()

"; 

%feature("docstring") simuPOP::Operator::canApplyPreOrPostMating "

Description:

    set of this operator can be applied pre-  or post-mating

Usage:

    x.canApplyPreOrPostMating()

"; 

%ignore simuPOP::Operator::isCompatible(const population &pop);

%feature("docstring") simuPOP::Operator::haploidOnly "

Description:

    determine if the operator can be applied only for haploid
    population

Usage:

    x.haploidOnly()

"; 

%feature("docstring") simuPOP::Operator::diploidOnly "

Description:

    determine if the operator can be applied only for diploid
    population

Usage:

    x.diploidOnly()

"; 

%feature("docstring") simuPOP::Operator::MPIReady "

Description:

    determine if this operator can be used in a MPI module

Usage:

    x.MPIReady()

"; 

%ignore simuPOP::Operator::setHaploidOnly();

%ignore simuPOP::Operator::setDiploidOnly();

%ignore simuPOP::Operator::setMPIReady();

%feature("docstring") simuPOP::Operator::infoSize "

Description:

    get the length of information fields for this operator

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

%ignore simuPOP::Operator::applyWithScratch(population &pop, population &scratch, int stage);

%feature("docstring") simuPOP::Operator::apply "

Description:

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::Operator::applyDuringMating(population &pop, population::IndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::Operator::setOutput "

Description:

    set ouput stream, if was not set during construction

Usage:

    x.setOutput(output=\"\", outputExpr=\"\")

"; 

%ignore simuPOP::Operator::getOstream(PyObject *dict=NULL, bool readable=false);

%ignore simuPOP::Operator::closeOstream();

%ignore simuPOP::Operator::atRepr();

%feature("docstring") simuPOP::Operator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::Operator::noOutput();

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::parentsTagger::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::parentsTagger::applyDuringMating(population &pop, population::IndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::pause "

Details:

    Pause the evolution of a simulator at given generations or at key
    stroke, using stopOnKeyStroke=True  option. When a simulator is
    stopped, user can resume simulation by pressing '/??' or escape to
    a python shell to examine the status of the simulation.

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

    prompt:         if true (default), print prompt message. findpause
    stopOnKeyStroke:if true, go on if no key was pressed
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

    used by Python print function to print out the general information
    of the pause

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

    used by Python print function to print out the general information
    of the operator

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

    used by Python print function to print out the general information
    of the operator

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
    structure, which refers to the number of chromosomes, number and
    position of loci on each chromosome etc. The most important
    components of a population are:
    * subpopulation. A population is divided into subpopulations
    (unstructured population has a single subpopulation, which is the
    whole population itself). Subpopulation structure limits the
    usually random exchange of genotypes between individuals
    disallowing mating between individuals from different
    subpopulations. In the presence of subpopualtion structure,
    exchange of genetic information across subpopulations can only be
    done through migration. Note that  simuPOP uses one-level
    population structure. I.e., there is no sub-subpopulation or
    families in subpopulations.
    * variables. Every population has its own variable space, or
    localnamespaces  in  simuPOP term. This namespace is a Python
    dictionary that is attached to each population and can be exposed
    to the users through  vars()  or dvars()  function. Many functions
    and operators work and store their results in these namespaces.
    For example, function Stat  set variables like alleleFreq[loc]
    and you can access them via pop.dvars().alleleFreq[loc][allele] .
    * ancestral generations. A population can save arbitrary number of
    ancestral generations. During evolution, the latest several (or
    all) ancestral generations are saved. Functions to make a certain
    ancestray generation current  is provided so that one can examine
    and modify ancestral generations. Other important concepts like
    informationfields  are explained in class individual.
    Note that although a large number of member functions are
    provided, most of the operations are performed by operators .
    These function...

Please refer to the reference manual for more details.

"; 

%feature("docstring") simuPOP::population::population "

Description:

    Create a population object with given size and genotypic
    structure.

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
      ...

Please refer to the reference manual for more details.

"; 

%ignore simuPOP::population::population(const population &rhs);

%feature("docstring") simuPOP::population::clone "

Description:

    deep copy of a population. (In python, pop1 = pop  will only
    create a reference to pop .)

Usage:

    x.clone(keepAncestralPops=-1)

"; 

%feature("docstring") simuPOP::population::swap "

Description:

    swap the content of two populations

Usage:

    x.swap(rhs)

"; 

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

    a python function used to compare the population objects

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
    allowPopSizeChange:if this parameter is true , population will be
                    resized.

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

    As with all ...End functions, the returning index is out of the
    range so that the actual range is [xxxBegin, xxxEnd). This agrees
    with all STL conventions and Python range.

"; 

%feature("docstring") simuPOP::population::ind "

Description:

    refernce to individual ind  in subpopulation subPop

Usage:

    x.ind(ind, subPop=0)

Arguments:

    ind:            individual index within subPop
    subPop:         subpopulation index

Details:

    This function is named individual  in the Python interface.

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

Arguments:

    order:          if order is true , respect order; otherwise, do
                    not repect population structure.

Details:

    Return an editable array of all genotypes of the population. You
    need to know how these genotypes are organized to safely
    read/write genotype directly. Individuals will be in order before
    exposing their genotypes.

"; 

%feature("docstring") simuPOP::population::arrGenotype "

Description:

    get the whole genotypes

Usage:

    x.arrGenotype(subPop, order)

Arguments:

    subPop:         index of subpopulation (start from 0)
    order:          if order is true , keep order; otherwise, respect
                    subpop  structure.

Details:

    Return an editable array of all genotype of a subpopulation.
    Individuals will be in order before exposing their genotypes.

"; 

%feature("docstring") simuPOP::population::setIndSubPopID "

Description:

    set subpopulation ID with given ID

Usage:

    x.setIndSubPopID(id)

Arguments:

    id:             an array of the same length of population size,
                    resprenting subpopulation ID of each individual.

Details:

    Set subpopulation ID of each individual with given ID. Individuals
    can be rearranged afterwards using setSubPopByIndID .

"; 

%feature("docstring") simuPOP::population::setIndSubPopIDWithID "

Description:

    set subpopulation ID of each individual with their current
    subpopulation ID

Usage:

    x.setIndSubPopIDWithID()

"; 

%feature("docstring") simuPOP::population::setSubPopByIndID "

Description:

    adjust subpopulation according to individual subpopulation ID.

Usage:

    x.setSubPopByIndID(id=[])

Arguments:

    id:             new subpopulation ID, if given, current individual
                    subpopulation ID will be ignored.

Details:

    Rearrange individuals to their new subpopulations according to
    their subpopulation ID (or the new given id ). Order within each
    subpopulation is not respected.

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

Note:

    subpop  with negative ID will be removed. So, you can shrink one
    subpop  by splitting and setting one of the new subpop  with
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

Note:

    subpop  with negative ID will be removed. So, you can shrink one
    subpop  by splitting and setting one of the new subpop  with
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
    will be no 'empty'  subpopulation left

Usage:

    x.removeSubPops(subPops=[], shiftSubPopID=True,
      removeEmptySubPops=False)

Details:

    Remove specified subpopulations (and all individuals within). If
    shiftSubPopID  is false, subPopID  will be kept intactly.

"; 

%feature("docstring") simuPOP::population::removeIndividuals "

Description:

    remove individuals

Usage:

    x.removeIndividuals(inds=[], subPop=-1,
      removeEmptySubPops=False)

Details:

    If a valid subPop  is given, remove individuals from this
    subpopulation.

"; 

%feature("docstring") simuPOP::population::mergeSubPops "

Description:

    merge given subpopulations

Usage:

    x.mergeSubPops(subPops=[], removeEmptySubPops=False)

Details:

    Merge subpopulations, the first subpopulation ID (the first one in
    array subPops ) will be used as the ID of the new subpopulation.
    That is to say, all subpopulations will take the ID of the first
    one.

"; 

%feature("docstring") simuPOP::population::mergePopulation "

Description:

    merge populations by individuals

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

Details:

    Merge individuals from pop  to the current population. Two
    populations should have the same genotypic structures. By default,
    subpopulations of the merged populations are kept. I.e., if you
    merge two populations with one subpopulation, the resulting
    population will have two subpopulations. All ancestral generations
    are also merged.

Note:

    Population variables are not copied to pop .

"; 

%feature("docstring") simuPOP::population::mergePopulationByLoci "

Description:

    merge populations by loci

Usage:

    x.mergePopulationByLoci(pop, newNumLoci=[], newLociPos=[])

Arguments:

    newLoci:        the new number of loci for the combined genotypic
                    structure.

Details:

    Two populations should have the same number of individuals. This
    also holds for any ancestral generations. By default, chromosomes
    of pop  are added to the current population. You can also specify
    new chromosome structure using parameter newLoci .

Note:

    * Information fields are not merged.
    * All ancestral generations will be merged because all individuals
    in a population have to have the same genotypic structure.

"; 

%feature("docstring") simuPOP::population::insertBeforeLoci "

Description:

    insert loci at given locations

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

Details:

    Insert loci at some given locations. In an inserted location,
    alleles will be zero.

"; 

%feature("docstring") simuPOP::population::insertBeforeLocus "

Description:

    insert an locus at given location.

Usage:

    x.insertBeforeLocus(idx, pos, name=string)

Details:

    insertBeforeLocus(idx, pos, name)  is a shortcut to
    insertBeforeLoci([idx], [pos], [name])

"; 

%feature("docstring") simuPOP::population::insertAfterLoci "

Description:

    append loci at given locations

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

Details:

    Append loci at some given locations. In an appended location,
    alleles will be zero.

"; 

%feature("docstring") simuPOP::population::insertAfterLocus "

Description:

    append an locus at a given location

Usage:

    x.insertAfterLocus(idx, pos, name=string)

Details:

    insertAfterLocus(idx, pos, name)  is a shortcut to
    insertAfterLoci([idx], [pos], [name]) .

"; 

%feature("docstring") simuPOP::population::resize "

Description:

    resize population

Usage:

    x.resize(newSubPopSizes, propagate=False)

Arguments:

    newSubPopSizes: an array of new subpopulation sizes. If there is
                    only one subpopulation, use [newPopSize] .
    propagate:      if propagate  is true , copy individuals to new
                    comers. I.e., 1, 2, 3 ==> 1, 2, 3, 1, 2, 3, 1

Details:

    Resize population by giving new subpopulation sizes.

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
                    0 1 means subpop3 , subpop2 , subpop0 , subpop1
                    will be the new layout.
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

    Form a new population according to the parameter information.
    Information can be given directly as
    * keepAncestralPops=-1 : keep all
    * keepAncestralPops=0 : only current
    * keepAncestralPops=1 : keep one ...

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

Usage:

    x.newPopWithPartialLoci(remove=[], keep=[])

Details:

    Copy current population to a new one with selected loci and remove
    specified loci. (No change on the current population.)

"; 

%feature("docstring") simuPOP::population::pushAndDiscard "

Description:

    Absorb rhs  population as the current generation of a population.

Usage:

    x.pushAndDiscard(rhs, force=False)

Details:

    This is mainly used by a simulator to push offspring generation
    rhs  to the current population, while the current population is
    pushed back as an ancestral population (if ancestralDepath() != 0
    ). Because rhs  population is swapped in, rhs  will be empty after
    this operation.

"; 

%feature("docstring") simuPOP::population::ancestralDepth "

Description:

    ancestral depth of the current population

Usage:

    x.ancestralDepth()

Note:

    The returned value is the number of ancestral generations exists
    in the population, not necessarily equal to the number set by
    setAncestralDepth().

"; 

%feature("docstring") simuPOP::population::ancestralGen "

Description:

    currently used ancestral population (0 for the latest generation)

Usage:

    x.ancestralGen()

Details:

    Current ancestral population activated by useAncestralPop. There
    can be several ancestral generations in a population. 0
    (current), 1  (parental) etc. When useAncestralPop(gen) is used,
    current generation is set to one of the parental generation, which
    is the information returned by this function. useAncestralPop(0)
    should always be used to set a population to its usual ancestral
    order.

"; 

%feature("docstring") simuPOP::population::setIndInfo "

Description:

    set individual information for the given information field
    (index),

Usage:

    x.setIndInfo(values, idx)

Arguments:

    values:         an array that has the same length as population
                    size.
    idx:            index to the information field.

"; 

%feature("docstring") simuPOP::population::setIndInfo "

Description:

    set individual information for the given information field (name)

Usage:

    x.setIndInfo(values, name)

Details:

    setIndInfo using field name, x.setIndInfo(values, name)  is
    equivalent to the idx  version x.setIndInfo(values,
    x.infoIdx(name)) .

"; 

%ignore simuPOP::population::infoBegin(UINT idx, bool order);

%ignore simuPOP::population::infoEnd(UINT idx, bool order);

%ignore simuPOP::population::infoBegin(UINT index, UINT subPop, bool order);

%feature("docstring") simuPOP::population::indInfo "

Description:

    get information idx  of all individuals.

Usage:

    x.indInfo(idx, order)

Arguments:

    idx:            index in all information fields
    order:          if true, sort returned vector in individual order

"; 

%feature("docstring") simuPOP::population::indInfo "

Description:

    get information name  of all individuals.

Usage:

    x.indInfo(name, order)

Arguments:

    name:           name of the information field
    order:          if true, sort returned vector in individual order

"; 

%feature("docstring") simuPOP::population::indInfo "

Description:

    get information idx  of all individuals in a subpopulation subPop
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

    get information name  of all individuals in a subpopulation subPop
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

    get an editable array (Python list) of all information fields

Usage:

    x.arrIndInfo(order)

Arguments:

    order:          whether or not the list has the same order as
                    individuals

"; 

%feature("docstring") simuPOP::population::addInfoField "

Description:

    add an information field to a population.

Usage:

    x.addInfoField(field, init=0)

Arguments:

    field:          new information field. If it already exists, it
                    will be re-initialized.
    init:           inital value for the new field.

"; 

%feature("docstring") simuPOP::population::addInfoFields "

Description:

    add one or more information fields to a population

Usage:

    x.addInfoFields(fields, init=0)

Arguments:

    fields:         new information fields. If one or more of the
                    fields alreay exist, they will simply be re-
                    initialized.
    init:           initial value for the new fields.

"; 

%feature("docstring") simuPOP::population::setInfoFields "

Description:

    set fields

Usage:

    x.setInfoFields(fields, init=0)

Arguments:

    fields:         an array of fields
    init:           initial value for the new fields.

"; 

%feature("docstring") simuPOP::population::setAncestralDepth "

Description:

    set ancestral depth, can be -1

Usage:

    x.setAncestralDepth(depth)

Arguments:

    depth:          0  for none, -1  for unlimited, a positive number
                    sets the number of ancestral generations to save.

"; 

%feature("docstring") simuPOP::population::useAncestralPop "

Description:

    use an ancestral generation. 0  for the latest generation.

Usage:

    x.useAncestralPop(idx)

Arguments:

    idx:            Index of the ancestral generation. 0  for current,
                    1  for parental, etc. idx can not exceed ancestral
                    depth (see  setAncestralDepth).

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

Details:

    Replication number is not meaningful for a stand-alone population.

"; 

%ignore simuPOP::population::setRep(int rep, bool setVar=true);

%feature("docstring") simuPOP::population::grp "

Description:

    current group ID in a simulator

Usage:

    x.grp()

Details:

    Group number is not meaningful for a stand-alone population.

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

Usage:

    x.evaluate(expr=\"\", stmts=\"\")

Details:

    This function evaluates a python statment/expression and return
    its result as a string. Optionally run statement first.

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

Usage:

    x.rearrangeLoci(newNumLoci, newLociPos)

Details:

    This is used by mergeByLoci .

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

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyEval::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    used by Python print function to print out the general information
    of the operator

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
                    will disallow random mating to set genotype.

Details:

    Note: (FIXME) output to output or outputExpr is not yet supported.
    Ideally, this func will take two parameters with pop and then a
    filehandle to output, however, differentiating output, append etc
    is too troublesome right now.

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

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyIndOperator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

Details:

    User of this operator must supply a Python function with parameter
    (index, ploidy, subpop). This operator will loop through all
    individual in each subpop and call this function to initialize
    populations.The arrange of parameters allows different
    initialization scheme for each subpop.

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

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::pyInit::apply "

Description:

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyMating "

Description:

    a Python mating scheme

Details:

    Hybird mating scheme. This mating scheme takes a Python function
    that accepts both the parental and offspring populations and this
    function is responsible for setting genotype, sex of the offspring
    generation. During-mating operators, if needed, have to be applied
    from this function as well. Note that the subpopulaton size
    parameters are honored and the passed offspring generation has the
    desired (sub)population sizes. Parameters that control the number
    of offspring of each family are ignored.
    This is likely an extremely slow mating scheme and should be used
    for experimental uses only. When a mating scheme is tested, it is
    recommended to implement it at the C++ level.

"; 

%feature("docstring") simuPOP::pyMating::pyMating "

Description:

    create a Python mating scheme

Usage:

    pyMating(*func=NULL, newSubPopSize=[], newSubPopSizeExpr=\"\",
      *newSubPopSizeFunc=NULL)

Arguments:

    func:           a Python function that accepts two parameters: the
                    parental and the offspring populations. The
                    offspring population is empty, and this function
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

    deep copy of a Python mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::pyMating::pyMating(const pyMating &rhs);

%feature("docstring") simuPOP::pyMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the Python mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::pyMating::mate(population &pop, population &scratch, vector< Operator * > &ops, bool submit);

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

Usage:

    pyMigrator(*subPopID=NULL, stage=PreMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    subPopID:       a 1-d array (list, typle, carray). Must has length
                    greater or equal to population size.
    stage:          is default to PreMating, please refer to
                    help(baseOperator.__init__) for details about
                    other parameters.

Details:

    This operator accept a one-dimensional Numeric Python int array.
    (created by Numeric.array ). The contend of the array will be
    considered as subpopulation id.

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyMigrator::apply "

Description:

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pyMigrator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    used by Python print function to print out the general information
    of the operator

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
                    offspring, this will imporve efficiency.

Details:

    Note: (FIXME) output to output or outputExpr is not yet supported.
    Ideally, this func will take two parameters with pop and then a
    filehandle to output, however, differentiating output, append etc
    is too troublesome right now.

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

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%ignore simuPOP::pyOperator::applyDuringMating(population &pop, population::IndIterator offspring, individual *dad=NULL, individual *mom=NULL);

%feature("docstring") simuPOP::pyOperator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

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

    used by Python print function to print out the general information
    of the operator

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

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

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

    used by Python print function to print out the general information
    of the operator

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

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::pySubset::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::pyTagger::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%ignore simuPOP::pyTagger::applyDuringMating(population &pop, population::IndIterator offspring, individual *dad=NULL, individual *mom=NULL);

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

    deep copy of an operator

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

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::randomMating "

Description:

    a mating scheme of basic sexually random mating

Details:

    In this scheme, sex information is considered for each individual,
    and ploidy is always 2. Within each subpopulation, males and
    females are randomly chosen. Then randomly get one copy of
    chromosomes from father and mother. When only one sex exists in a
    subpopulation, a parameter (contWhenUniSex ) can be set to
    determine the behavior. Default to continuing without warning.

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

    numOffspring:   number of offspring or p  in some modes
    numOffspringFunc:a python function that determines the number of
                    offspring or p
    maxNumOffspring:used when numOffspring  is generated from a
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
    newSubPopSizeFunc:an function that have parameter gen  and oldSize
                    (current subpopulation size)
    contWhenUniSex: continue when there is only one sex in the
                    population, default to true
                    Please refer to class mating  for descriptions of
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

    deep copy of a random mating scheme

Usage:

    x.clone()

"; 

%ignore simuPOP::randomMating::isCompatible(const population &pop) const ;

%feature("docstring") simuPOP::randomMating::__repr__ "

Description:

    used by Python print function to print out the general information
    of the random mating scheme

Usage:

    x.__repr__()

"; 

%ignore simuPOP::randomMating::submitScratch(population &pop, population &scratch);

%ignore simuPOP::randomMating::mate(population &pop, population &scratch, vector< Operator * > &ops, bool submit);

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

    used by Python print function to print out the general information
    of the operator

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

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::recombinator::prepareRecRates "

Description:

    this function takes intensity, rate, afterLoci, ... inputs and
    return a bernulli trailer and a recBeforeLoci vector.

Usage:

    x.prepareRecRates(pop, intensity, rate, afterLoci, sexChrom,
      recBeforeLoci, vecP)

Details:

    get loci distance * rate and then recombinant pointsget loci
    distance * rate and then recombinant pointsinitialize
    recombination counter, This will count recombination events after
    each locus.

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

%ignore simuPOP::recombinator::applyDuringMating(population &pop, population::IndIterator offspring, individual *dad=NULL, individual *mom=NULL);

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

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::sample::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

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

    used by Python print function to print out the general information
    of the operator

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

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::setAncestralDepth::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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
    using given mating scheme and operators

Details:

    Simulators combine three important components of  simuPOP:
    population , mating scheme and operators together. A simulator is
    created with an instance of population , a replicate number rep
    and a mating scheme. It makes rep  number of replicates of this
    population and control the evolution process of them.
    The most important function of a simulator is  evolve() . It
    accepts an array of operators as its parameters, among which,
    preOps  and postOps  will be applied to the populations at the
    beginning and the end of evolution, respectively, whereas ops
    will be applied at every generation.
    Simulators separate operators into pre- , during-  and post-mating
    operators. During evolution, a simulator first apply all pre-
    mating operators and then call the mate()  function of the given
    mating scheme, which will call during-mating operators during the
    birth of each offspring. After mating is completed, post-mating
    operators are applied to the offspring in the order at which they
    appear in the operator list.
    Operators can be applied to specific replicate, a group of
    replicates, or specific generations, determined by the rep , grp ,
    begin , end , step , and at  parameters.
    Simulators can evolve a given number of generations (the end
    parameter of evolve ), or evolve indefinitely until a certain type
    of operators called terminators terminates it. In this case, one
    or more terminators will check the status of evolution and
    determine if the simulation should be stopped. An obvious example
    of such a terminator is a fixation-checker.
    Finally, a simulator can be saved to a file in the format of 'txt'
    , 'bin' , or 'xml' . So we can stop a simulation and resume it at
    another time or on another machine. It is a...

Please refer to the reference manual for more details.

"; 

%feature("docstring") simuPOP::simulator::simulator "

Description:

    create a simulator

Usage:

    simulator(pop, matingScheme, stopIfOneRepStops=False,
      applyOpToStoppedReps=False, rep=1, grp=[])

Arguments:

    population:     a population created by population()  function.
                    This population will be copied rep  times to the
                    simulator. Its content will not be changed.
    matingScheme:   a mating scheme
    rep:            number of replicates. Default to 1 .
    grp:            group number for each replicate. Operators can be
                    applied to a group of replicates using its grp
                    parameter.
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

    pop = simulator::population()  returns temporary reference to an
    internal population. After a simulator evolves another genertion
    or after the simulator is destroyed, this referenced population
    should not  be used.

"; 

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

Arguments:

    field:          information field to be added

Details:

    Add an information field to all replicate, and to the simulator
    itself. This is important because all populations must have the
    same genotypic information as the simulator. Adding an information
    field to one or more of the replicates will compromise the
    integrity of the simulator.

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

    the rep  replicate of this simulator

Usage:

    x.pop(rep)

Arguments:

    rep:            the index number of replicate which will be
                    accessed

Details:

    This function is named population  in the Python interface.

Note:

    The returned reference is temporary in the sense that the refered
    population will be invalid after another round of evolution.
    Therefore, the use of this function should be limited to
    immediateafterretrival . If you would like to get a persistent
    population, please use getPopulation(rep) .

"; 

%feature("docstring") simuPOP::simulator::getPopulation "

Description:

    return a copy of population rep

Usage:

    x.getPopulation(rep)

Arguments:

    rep:            the index number of the replicate which will be
                    obtained

Details:

    return a temporary reference of one of the populations.
    'Reference'  means that the changes to the referred population
    will reflect to the one in simulator. 'Temporary'  means that the
    referred population might be invalid after evolution.

"; 

%feature("docstring") simuPOP::simulator::setMatingScheme "

Description:

    set mating scheme

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

    return current generation number

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

    set current generation. Usually used to reset a simulator.

Usage:

    x.setGen(gen)

Arguments:

    gen:            new generation index number

"; 

%feature("docstring") simuPOP::simulator::step "

Description:

    evolve steps  generation

Usage:

    x.step(ops=[], preOps=[], postOps=[], steps=1, dryrun=False)

"; 

%feature("docstring") simuPOP::simulator::evolve "

Description:

    evolve all replicates of the population, subject to operators

Usage:

    x.evolve(ops, preOps=[], postOps=[], end=-1, dryrun=False)

Arguments:

    ops:            operators that will be applied at each generation,
                    if they are active at that generation. (Determined
                    by the begin , end , step  and at  parameters of
                    the operator.)
    preOps:         operators that will be applied before evolution.
                    evolve()  function will not  check if they are
                    active.
    postOps:        operators that will be applied after evolution
    end:            ending generation. Default to -1 . In this case,
                    there is no ending generation and a simulator will
                    only be ended by a terminator. Otherwise, it
                    should be a number greater than current generation
                    number.

Details:

    Evolve to the end  generation unless an operator (terminator)
    stops it earlier.
    ops  will be applied in the order of:
    * all pre-mating opertors
    * during-mating operators called by the mating scheme at the birth
    of each offspring
    * all post-mating operators If any pre- or post-mating operator
    fails to apply, that replicate will be stopped. The behavior of
    the simulator will be determined by flags applyOpToStoppedReps
    and stopIfOneRepStopss . This is exactly how terminators work.

Note:

    When end = -1 , you can not specify negative generation parameters
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

    get simulator namespace. If rep > 0  is given, return the
    namespace of replicate rep

Usage:

    x.vars(rep, subPop=-1)

"; 

%feature("docstring") simuPOP::simulator::saveSimulator "

Description:

    save simulator in 'txt' , 'bin'  or 'xml'  format

Usage:

    x.saveSimulator(filename, format=\"auto\", compress=True)

Arguments:

    filename:       filename to save the simulator. Default to simu .
    format:         format to save. Default to auto . I.e., determine
                    the format by file extensions.
    compress:       whether or not compress the file in 'gzip'  format

Details:

    The default format is 'txt'  but the output is not supposed to be
    read. 'bin'  has smaller size and should be used for large
    populations. 'xml'  format is most verbose and should be used when
    you would like to convert  simuPOP  populations to other formats.

"; 

%ignore simuPOP::simulator::loadSimulator(string filename, string format="auto");

%feature("docstring") simuPOP::simulator::__repr__ "

Description:

    used by Python print function to print out the general information
    of the simulator

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

Usage:

    smmMutator(rate=[], atLoci=[], maxAllele=0, incProb=0.5,
      output=\">\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
      step=1, at=[], rep=REP_ALL, grp=GRP_ALL, infoFields=[])

Arguments:

    rate::          mutation rate
    incProb:        probability to increase allele state. Default to 1
    atLoci:         and other parameters: refer to help(mutator),
                    help(baseOperator.__init__)

Details:

    The stepwise mutation model (SMM) is developed for allozymes. It
    provides better description for these kinds of evolutionary
    processes.

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

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::splitSubPop::apply "

Description:

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::splitSubPop::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    used by Python print function to print out the general information
    of the operator

Usage:

    x.__repr__()

"; 

%feature("docstring") simuPOP::spread::apply "

Description:

    apply to one population. It does not check if the operator is
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
                    counted
    heteroFreq:     an array of loci at which the observed proportion
                    of individuausl heterozygous will be applyd for
                    each allele. Expected heterozygosity will also be
                    calculuate and put in heteroFreq[locus][0] (since
                    allele 0 is not...

Please refer to the reference manual for more details.

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

    used by Python print function to print out the general information
    of the operator

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

Usage:

    x.apply(pop)

Details:

    go through a single allele for all individual, all diploidneed to
    replace previous values

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

    deep copy of an operator

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

    deep copy of an operator

Usage:

    x.clone()

"; 

%feature("docstring") simuPOP::terminateIf::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    deep copy of an operator

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

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::ticToc::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::turnOffDebugOp::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

    apply to one population. It does not check if the operator is
    activated.

Usage:

    x.apply(pop)

"; 

%feature("docstring") simuPOP::turnOnDebugOp::__repr__ "

Description:

    used by Python print function to print out the general information
    of the operator

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

