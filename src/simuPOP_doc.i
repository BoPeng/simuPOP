

%feature("docstring") arraydescr "

";

%feature("docstring") arrayobject "

";

%feature("docstring") simuPOP::BasicPenetrance "

penetrance according to genotype at one locus

Details:
  basic selector. Assign penetrance value according to a given
  dictionary.

";

%feature("docstring")  simuPOP::BasicPenetrance::BasicPenetrance " 

create a basic penetrance function (penetrance according to genotype
at one locus

Usage:
  basicPenetrance(locus, &penetrance, hasPhase=False,
    exposePenetrance=False, stage=DuringMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  penetrance:  a dictionary of penetrance. The genotype must be in the
      form of'a-b'.
      
  hasPhase:  if true, a/b and b/a will have different penetrance
      value. Default to false.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::BasicPenetrance::~BasicPenetrance " 

Usage:
  x.~BasicPenetrance()

";

%feature("docstring")  simuPOP::BasicPenetrance::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::BasicPenetrance::penet " 

currently assuming diploid

Usage:
  x.penet(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::BasicPenetrance::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::BasicQuanTrait "

quantitative trait according to genotype at one locus

Details:
  basic selector. Assign qtrait value according to a given dictionary.

";

%feature("docstring")  simuPOP::BasicQuanTrait::BasicQuanTrait " 

create a basic selector (quantitative trait according to genotype at
one locus

Usage:
  basicQuanTrait(locus, &qtrait, sigma, hasPhase=False,
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  qtrait:  a dictionary of qtrait. The genotype must be in the form of'a
      -b'. This is the mean of quantitative trait. The actual
      value will be N(mean, sigma)
      
  sigma:  sigma of the normal distribution of the trait.
      
  hasPhase:  if true, a/b and b/a will have different qtrait value.
      Default to false.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::BasicQuanTrait::~BasicQuanTrait " 

Usage:
  x.~BasicQuanTrait()

";

%feature("docstring")  simuPOP::BasicQuanTrait::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::BasicQuanTrait::qtrait " 

currently assuming diploid

Usage:
  x.qtrait(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::BasicQuanTrait::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::BasicSelector "

selection according to genotype at one locus

Details:
  basic selector. Assign fitness value according to a given
  dictionary.

";

%feature("docstring")  simuPOP::BasicSelector::BasicSelector " 

create a basic selector (selection according to genotype at one
locus

Usage:
  basicSelector(locus, &fitness, hasPhase=False, stage=PreMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
    sep=\"\\\\t\")

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  fitness:  a dictionary of fitness. The genotype must be in the form
      of'a-b'.
      
  hasPhase:  if true, a/b and b/a will have different fitness value.
      Default to false.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::BasicSelector::~BasicSelector " 

Usage:
  x.~BasicSelector()

";

%feature("docstring")  simuPOP::BasicSelector::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::BasicSelector::fitness " 

currently assuming diploid

Usage:
  x.fitness(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::BasicSelector::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::BasicStat "

";

%feature("docstring")  simuPOP::BasicStat::BasicStat " 

create an basicStat

Usage:
  basicStat(popSize=False, numOfMale=False, numOfAffected=False,
    numOfAlleles=vectori, alleleFreq=vectori, heteroFreq=vectori,
    genoFreq=vectori, haploFreq=intMatrix, LD=intMatrix, Fst=vectori,
    hasPhase=False, output=\"\", outputExpr=\"\", headers=[],
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  popSize:  whether or not apply population sizes, will set numSubPop,
      subPopSize, popSize, subPop[sp]['p
      opSize']
      
  numOfMale:  whether or not count number of male and female, will set
      numOfMale and numOfFemale for all population/subpopulations
      
  numOfAffected:  whether or not count number of affected individuals.
      Will set numOfAffected, and numOfUnaffected.
      
  numOfAlleles:  an array of loci at which number of alleles will be
      counted (0 is excluded). Note that number of alleles will be
      automatically set if alleleFreq is counted.
      
  alleleFreq:  an array of loci at which all alleles will be counted.
      
  genoFreq:  an array of loci at which all genotype will be counted
      
  each item is the locus index followed by allele pairs.

  heteroFreq:  an array of loci at which the observed proportion of
      individuausl heterozygous will be applyd for each allele.
      Expected heterozygosity will also be calculuate and put in
      heteroFreq[locus][0] (since allele 0 is not used.)
      
  haploFreq:  a matrix of haplotypes (allele sequence on the same
      chromosome) to count
      
  format: haploFreq = [ [ 0,1,2,1,2,3 ], [1,2,1,2] ]
  the first half of each item is loci number. the second half is
  allele number
  [0,1,2,1,2,3] is counting haplotype 1,2,3 (on the same chromosome)
  at loci 0,1,2.

  LD:  apply LD, LD'and r2, given
      LD=[ [ locus1 locus2 allele1 allele2], [locus1 locus2 allele1
      allele2],...]
      
  Fst:  calculate Fst. Fis and Fit will be given as a side product.
      
  format Fst = [ 0, 1 ] Fst calculated at locus 0 and 1. Since any
  allele can be used to calculate Fst, Fst[0] will be the average over
  all alleles as suggested by Weir and Cockerham.

  hasPhase:  if a/b and b/a are the same genotype. default is false.
      
  others:  there is NO output for this operator. other parameters
      please see help(baseOperator.__init__)

Note:
  previous provide'search for all
  allele/genotype'option but this has proven to be
  troublesome. In this version, everything should be explicitly
  specified.

";

%feature("docstring")  simuPOP::BasicStat::~BasicStat " 

Usage:
  x.~BasicStat()

";

%feature("docstring")  simuPOP::BasicStat::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::BasicStat::apply " 

count various statistics. use m_alleles etc to save (potentially)
time to resize all these variables.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::BasicStat::__repr__ " 

Usage:
  x.__repr__()

";

%ignore simuPOP::BernulliTrials;

%ignore simuPOP::BernulliTrials::BernulliTrials(RNG &rng);

%ignore simuPOP::BernulliTrials::BernulliTrials(RNG &rng, const vectorf &prob, ULONG trials);

%ignore simuPOP::BernulliTrials::~BernulliTrials();

%ignore simuPOP::BernulliTrials::size();

%feature("docstring")  simuPOP::BernulliTrials::prob " 

Usage:
  x.prob()

";

%ignore simuPOP::BernulliTrials::setParameter(const vectorf &prob, ULONG trials);

%ignore simuPOP::BernulliTrials::doTrial();

%ignore simuPOP::BernulliTrials::trial();

%feature("docstring")  simuPOP::BernulliTrials::succ " 

return succeed trials for p[index] fail when m_cur is not 0. (i.e.,
has retrieve the table through trial()C
PPONLY

Usage:
  x.succ(index)

";

%ignore simuPOP::BernulliTrials::probabilities();

%ignore simuPOP::BernulliTrials::numTrials();

%feature("docstring") simuPOP::BinomialSelection "

Details:
  binomial random selection
  No sex. Choose one individual from last generation. Subpopulations
  are dealt separately.

";

%feature("docstring")  simuPOP::BinomialSelection::BinomialSelection " 

constructor

Usage:
  binomialSelection(numOffsprings=1, newSubPopSize=[],
    newSubPopSizeExpr=\"\", *newSubPopSizeFunc=NULL)

";

%feature("docstring")  simuPOP::BinomialSelection::~BinomialSelection " 

destructor

Usage:
  x.~BinomialSelection()

";

%feature("docstring")  simuPOP::BinomialSelection::clone " 

 clone() const. The same as Mating::clone() const.

Usage:
  x.clone()

See Also:
   Mating::clone() const

";

%feature("docstring")  simuPOP::BinomialSelection::__repr__ " 

return name of the mating type

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::BinomialSelection::mate " 

do the mating.

Usage:
  x.mate(&pop, &scratch, &ops)

Arguments:

  pop:  population
      
  scratch:  scratch population, will be used in this mating scheme.
      
  ops:  during mating operators

Details:

  determine if mate() will generate offspring genotype
  use deep copy!!!!!!!

Value:
  return false when mating fails.

";

%feature("docstring") simuPOP::CaseControlSample "

thrink population accroding to some outside value

";

%feature("docstring")  simuPOP::CaseControlSample::CaseControlSample " 

draw cases and controls

Usage:
  caseControlSample(cases=0, controls=0, &name=\"sample\",
    &nameExpr=\"\", times=1, &saveAs=\"\", &saveAsExpr=\"\",
    &format=\"bin\", stage=PostMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  cases:  number of cases
      
  controls:  number of controls
      
  stage:  and other parameters please see help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::CaseControlSample::~CaseControlSample " 

destructor

Usage:
  x.~CaseControlSample()

";

%feature("docstring")  simuPOP::CaseControlSample::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::CaseControlSample::drawSample " 

Usage:
  x.drawSample(&pop)

";

%feature("docstring")  simuPOP::CaseControlSample::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::Dumper "

dump the content of a population.

";

%feature("docstring")  simuPOP::Dumper::Dumper " 

dump population

Usage:
  dumper(alleleOnly=False, infoOnly=False, dispWidth=1, max=100,
    output=\">\", outputExpr=\"\", headers=[], stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
    sep=\"\\\\t\")

Arguments:

  alleleOnly:  only display allele
      
  infoOnly:  only display info
      
  dispWidth:  width of allele display, default to 1
      
  max:  max number of individuals to display, default to 100. This is
      to avoid careless dump of huge populations.
      
  output:  output file, default to standard output.
      
  outputExpr:  and other parameters: refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::Dumper::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::Dumper::alleleOnly " 

only show alleles (not structure, gene information?

Usage:
  x.alleleOnly()

";

%feature("docstring")  simuPOP::Dumper::setAlleleOnly " 

Usage:
  x.setAlleleOnly(alleleOnly)

";

%feature("docstring")  simuPOP::Dumper::infoOnly " 

only show info

Usage:
  x.infoOnly()

";

%feature("docstring")  simuPOP::Dumper::setInfoOnly " 

Usage:
  x.setInfoOnly(infoOnly)

";

%feature("docstring")  simuPOP::Dumper::apply " 

Usage:
  x.apply(&pop)

Details:
  dump population structure
  dump all genotypic info

";

%feature("docstring")  simuPOP::Dumper::~Dumper " 

Usage:
  x.~Dumper()

";

%feature("docstring")  simuPOP::Dumper::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::Exception "

exception handler. Exceptions will be passed to Python.

";

%feature("docstring")  simuPOP::Exception::Exception " 

constructor str error message

Usage:
  exception(msg)

";

%feature("docstring")  simuPOP::Exception::message " 

return error message

Usage:
  x.message()

";

%feature("docstring")  simuPOP::Exception::~Exception " 

Usage:
  x.~Exception()

";

%ignore simuPOP::Expression;

%ignore simuPOP::Expression::Expression(const string &expr="", const string &stmts="", PyObject *locals=NULL);

%ignore simuPOP::Expression::~Expression();

%feature("docstring")  simuPOP::Expression::Expression " 

Copy constructor, need to be defined because of ref count issue.

Usage:
  expression(&rhs)

";

%ignore simuPOP::Expression::setLocalDict(PyObject *dict);

%ignore simuPOP::Expression::empty();

%ignore simuPOP::Expression::setExpr(const string &expr="");

%ignore simuPOP::Expression::setStmts(const string &stmts="");

%feature("docstring")  simuPOP::Expression::evaluate " 

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

%ignore simuPOP::GappedAlleleIterator;

%ignore simuPOP::GappedAlleleIterator::GappedAlleleIterator();

%ignore simuPOP::GappedAlleleIterator::GappedAlleleIterator(pointer p, difference_type s=1);

%ignore simuPOP::GappedAlleleIterator::GappedAlleleIterator(const GappedAlleleIterator &rhs);

%ignore simuPOP::GappedAlleleIterator::~GappedAlleleIterator();

%feature("docstring")  simuPOP::GappedAlleleIterator::ptr " 

get pointer so that ptr()+1 will be the next
locus.

Usage:
  x.ptr()

";

%ignore simuPOP::GappedAlleleIterator::step();

%ignore simuPOP::GenoStructure;

%ignore simuPOP::GenoStructure::GenoStructure();

%ignore simuPOP::GenoStructure::GenoStructure(UINT ploidy, const vectoru &loci, const matrix &lociDist, const vectorstr &alleleNames, UINT maxAllele);

%ignore simuPOP::GenoStructure::GenoStructure(const GenoStructure &rhs);

%feature("docstring")  simuPOP::GenoStructure::~GenoStructure " 

destructor, do nothing.

Usage:
  x.~GenoStructure()

";

%feature("docstring") simuPOP::GenoStruTrait "

genoStruTrait

Details:
  A trait class that provide interfaces around a GenoStructurep
  ointer.

";

%feature("docstring")  simuPOP::GenoStruTrait::GenoStruTrait " 

constructor, but m_genoStru will be set later.

Usage:
  genoStruTrait()

";

%ignore simuPOP::GenoStruTrait::setGenoStructure(GenoStructure *stru);

%ignore simuPOP::GenoStruTrait::destroyGenoStructure();

%ignore simuPOP::GenoStruTrait::genoStru();

%feature("docstring")  simuPOP::GenoStruTrait::ploidy " 

return ploidy

Usage:
  x.ploidy()

";

%feature("docstring")  simuPOP::GenoStruTrait::ploidyName " 

return ploidy

Usage:
  x.ploidyName()

";

%feature("docstring")  simuPOP::GenoStruTrait::numLoci " 

number of loci on chromosome chrom

Usage:
  x.numLoci(chrom)

";

%feature("docstring")  simuPOP::GenoStruTrait::totNumLoci " 

return totNumLoci (STATIC)

Usage:
  x.totNumLoci()

";

%feature("docstring")  simuPOP::GenoStruTrait::genoSize " 

return totNumLoci * ploidy

Usage:
  x.genoSize()

";

%feature("docstring")  simuPOP::GenoStruTrait::locusDist " 

locus distance.

Usage:
  x.locusDist(locus)

";

%feature("docstring")  simuPOP::GenoStruTrait::arrLociDist " 

return loci distance as python Numeric.array object

Usage:
  x.arrLociDist()

";

%feature("docstring")  simuPOP::GenoStruTrait::numChrom " 

number of chromosome

Usage:
  x.numChrom()

";

%feature("docstring")  simuPOP::GenoStruTrait::chromBegin " 

chromosome index of chromosome chrom

Usage:
  x.chromBegin(chrom)

";

%feature("docstring")  simuPOP::GenoStruTrait::chromEnd " 

chromosome index of chromosome chrom

Usage:
  x.chromEnd(chrom)

";

%feature("docstring")  simuPOP::GenoStruTrait::absLocusIndex " 

convert from relative locus (on chromsome) to absolute locus (no
chromosome structure)

Usage:
  x.absLocusIndex(chrom, locus)

";

%feature("docstring")  simuPOP::GenoStruTrait::chromLocusPair " 

return chrom, locus pair from an absolute locus position.

Usage:
  x.chromLocusPair(locus)

";

%feature("docstring")  simuPOP::GenoStruTrait::alleleName " 

return allele name

Usage:
  x.alleleName(allele)

";

%feature("docstring")  simuPOP::GenoStruTrait::maxAllele " 

Usage:
  x.maxAllele()

";

%feature("docstring") simuPOP::GSMMutator "

stepwise mutation model.

Details:
  Generalized Stepwise mutation model (GSM) assumes that alleles are
  represented by integer values and that a mutation either increases
  or decreases the allele value by a random value.

See Also:
  Kimura&Ohta 1978

";

%feature("docstring")  simuPOP::GSMMutator::GSMMutator " 

Usage:
  gSMMutator(rate=0, rates=[], atLoci=vectori, maxAllele=0,
    incProb=0.5, p=0, *func=NULL, output=\">\", outputExpr=\"\",
    headers=[], stage=PostMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  rate:  or rates: mutation rate
      
  incProb:  probability to increase allele state. Default to 1
      
  atLoci:  and other parameters: refer to help(mutator),
      help(baseOperator.__init__)
      
  func:  return number of steps. no parameter

Details:
  The generalized stepwise mutation model (GMM) is developed for
  allozymes. It provides better description for these kinds of
  evolutionary processes.

";

%feature("docstring")  simuPOP::GSMMutator::~GSMMutator " 

Usage:
  x.~GSMMutator()

";

%feature("docstring")  simuPOP::GSMMutator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::GSMMutator::mutate " 

how to mutate a single allele. this is usually the only function
that need to be defined by the subclasses.

Usage:
  x.mutate(*allele)

";

%feature("docstring")  simuPOP::GSMMutator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::IfElse "

";

%feature("docstring")  simuPOP::IfElse::IfElse " 

conditional operator

Usage:
  ifElse(&cond, *ifOp=NULL, *elseOp=NULL, output=\">\",
    outputExpr=\"\", headers=[], stage=PostMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  cond:  expression, will be treated as bool variable.
      
  ifOp:  if operator, be called when expr is true
      
  elseOp:  else operator, be called when expr is false

";

%feature("docstring")  simuPOP::IfElse::~IfElse " 

destructor

Usage:
  x.~IfElse()

";

%ignore simuPOP::IfElse::IfElse(const IfElse< Pop > &rhs);

%feature("docstring")  simuPOP::IfElse::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::IfElse::applyWithScratch " 

simply output some info providing interface to apply operator before
during or after mating.

Usage:
  x.applyWithScratch(&pop, &scratch, stage)

";

%feature("docstring")  simuPOP::IfElse::applyDuringMating " 

give pop, offspring, pop and mom.

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring")  simuPOP::IfElse::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::IfElse::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::IndexError "

exception, thrown if index out of range

";

%feature("docstring")  simuPOP::IndexError::IndexError " 

Usage:
  indexError(msg)

";

%feature("docstring") simuPOP::Individual "

Basic individual class.

Details:
  class individual withgenotypic information
  shared genotypic structure info (through a GenoStructurep
  ointer)
  flags about sex, affected status
  an internal info field and a tag field of any type (template)
  other individuals will be derived from this class, adding age info
  etc.
   Notethat IndividualDOES NOT manage
  memory. It will use a pointer passed from class Population.
   This causes A LOT of trouble and I have not evaluated how much
  benefic I get.
   Individualis a template class taking a tag parameter. tag can be a
  integer, a pair or any object to track individual information. The
  prolem with our template based design + SWIG makes instantiation of
  mutiple tag types difficult. :-(
  operator = uses shallow copy. This is required by sort algorithm
  since otherwise individuals are non-copiable. However, in population
  memory management, it is showtimes required that genotypic
  information within one subPop should go together. This is done by a
  shollow_copied flag for each individual and for all individuals. Popu
  lationmight have to re-arrange individuals to solve this problem.
  output of individual can be adjusted by setOutputDelimeter.
  Usage info: (for population classes developers)for individuals are
  created, you are responsible to set its genotypic pointer and
  genotypic information. This is done by
   setInfo()and info()can be
  used for any temporarypurpose.
   Tagare usually std::pair(int,int) used by class Taggert
  o track informations like pedigree structure.

";

%feature("docstring")  simuPOP::Individual::Individual " 

default constructor, Tag field need a default constructor

Usage:
  individual()

";

%ignore simuPOP::Individual::Individual(const Individual< Tag > &ind);

%feature("docstring")  simuPOP::Individual::~Individual " 

destructor. Do nothing.

Usage:
  x.~Individual()

";

%ignore simuPOP::Individual::setGenoPtr(ALLELE *pos);

%feature("docstring")  simuPOP::Individual::copyFrom " 

Deep copy! Important!

Usage:
  x.copyFrom(&rhs)

Details:
  can not do this. Imaging that we copy one individual to scratch
  population which has its own genoStructure!

";

%ignore simuPOP::Individual::genoPtr();

%feature("docstring")  simuPOP::Individual::arrAlleles " 

return genotype as python Numeric.array object This is the whole
genotype (all)

Usage:
  x.arrAlleles()

";

%feature("docstring")  simuPOP::Individual::arrAlleles " 

return genotype as python Numeric.array object This is the p't
h copy of chromosomes

Usage:
  x.arrAlleles(p)

";

%feature("docstring")  simuPOP::Individual::arrAlleles " 

return genotype as python Numeric.array object This is the ch
chromosome of the pth copy of chromosome

Usage:
  x.arrAlleles(p, ch)

";

%feature("docstring")  simuPOP::Individual::allele " 

get allele from an index

Usage:
  x.allele(index)

Arguments:

  index:  index from the beginning of genotypic info

";

%feature("docstring")  simuPOP::Individual::allele " 

get allele from an index, on the pth set of chromosome

Usage:
  x.allele(index, p)

Arguments:

  index:  index from the begining of the p't
      h set of chromosomes.
      
  p:  on p'th set of chromosomes, default to 0

";

%feature("docstring")  simuPOP::Individual::allele " 

Usage:
  x.allele(index, p, ch)

";

%feature("docstring")  simuPOP::Individual::alleleChar " 

Usage:
  x.alleleChar(index)

";

%feature("docstring")  simuPOP::Individual::alleleChar " 

get allele from an index, on the pth set of chromosome

Usage:
  x.alleleChar(index, p)

Arguments:

  index:  index from the begining of the p't
      h set of chromosomes.
      
  p:  on p'th set of chromosomes, p=0 by default

";

%feature("docstring")  simuPOP::Individual::alleleChar " 

get allele from an index, on the pth set of chromosome

Usage:
  x.alleleChar(index, p, ch)

Arguments:

  index:  index from the begining of the p't
      h set of chromosomes.
      
  p:  on p'th set of chromosomes, p=0 by default

";

%feature("docstring")  simuPOP::Individual::setAllele " 

set allele from an index.

Usage:
  x.setAllele(allele, index)

Arguments:

  index:  index from the begining of genotype

";

%feature("docstring")  simuPOP::Individual::setAllele " 

set allele from an index.

Usage:
  x.setAllele(allele, index, p)

Arguments:

  allele:  allele to set
      
  index:  index from the begining of genotype
      
  p:  on p'th set of chromosome, p=0 by default

";

%feature("docstring")  simuPOP::Individual::tag " 

return tag

Usage:
  x.tag()

";

%feature("docstring")  simuPOP::Individual::setTag " 

set tag

Usage:
  x.setTag(tag)

";

%feature("docstring")  simuPOP::Individual::sex " 

sex?

Usage:
  x.sex()

";

%feature("docstring")  simuPOP::Individual::sexChar " 

return M or F for sex, for display purpose

Usage:
  x.sexChar()

";

%feature("docstring")  simuPOP::Individual::setSex " 

set sex

Usage:
  x.setSex(sex)

";

%feature("docstring")  simuPOP::Individual::affected " 

affected?

Usage:
  x.affected()

";

%feature("docstring")  simuPOP::Individual::affectedChar " 

return A or U for affected/Unaffected, for display purpose

Usage:
  x.affectedChar()

";

%feature("docstring")  simuPOP::Individual::setAffected " 

set affected status

Usage:
  x.setAffected(affected)

";

%feature("docstring")  simuPOP::Individual::info " 

get info

Usage:
  x.info()

";

%feature("docstring")  simuPOP::Individual::setInfo " 

set info

Usage:
  x.setInfo(info)

";

%ignore simuPOP::Individual::genoBegin();

%ignore simuPOP::Individual::genoEnd();

%ignore simuPOP::Individual::genoBegin(UINT p);

%ignore simuPOP::Individual::genoEnd(UINT p);

%ignore simuPOP::Individual::genoBegin(UINT p, UINT chrom);

%ignore simuPOP::Individual::genoEnd(UINT p, UINT chrom);

%feature("docstring")  simuPOP::Individual::swap " 

swap individuals

Usage:
  x.swap(&ind, swapContent=True)

Arguments:

  ind:  individual to be swapped in
      
  swapContent:  swapContent or only the pointers.
      
  The guideline is that if we swap individuals across subPopulation,
  we should swap content. Otherwise, swap pointers. (There is no order
  right now within subPopulation so the later case is rare, at best.

Details:
  The default behavior is swapping all info, but not the position of
  genotypic info. If swapContent is false, pointer to genotypic info
  is swapped instead. This will lead to better performance for
  swapping but may affected performance of allele counting.

";

%ignore simuPOP::Individual::shallowCopied();

%ignore simuPOP::Individual::setShallowCopied(bool shallowCopied);

%ignore simuPOP::Individual::shallowCopiedFlagOn();

%ignore simuPOP::Individual::clearShallowCopiedFlag();

%feature("docstring")  simuPOP::Individual::dispWidth " 

display width

Usage:
  x.dispWidth()

";

%feature("docstring")  simuPOP::Individual::setDispWidth " 

set seperator

Usage:
  x.setDispWidth(width)

";

%feature("docstring") simuPOP::IndividualWithAge "

individual with age info

";

%feature("docstring")  simuPOP::IndividualWithAge::age " 

get age

Usage:
  x.age()

";

%feature("docstring")  simuPOP::IndividualWithAge::setAge " 

set age

Usage:
  x.setAge(age)

";

%feature("docstring")  simuPOP::IndividualWithAge::IndividualWithAge " 

default constructor

Usage:
  individualWithAge()

";

%feature("docstring")  simuPOP::IndividualWithAge::IndividualWithAge " 

copy constructor. will be a shallow copied one

Usage:
  individualWithAge(&ind)

";

%feature("docstring")  simuPOP::IndividualWithAge::~IndividualWithAge " 

destructor

Usage:
  x.~IndividualWithAge()

";

%feature("docstring")  simuPOP::IndividualWithAge::swap " 

swap info. See Individual::swap

Usage:
  x.swap(&ind, swapContent=True)

";

%feature("docstring")  simuPOP::IndividualWithAge::serialize " 

Usage:
  x.serialize(&ar, version)

";

%feature("docstring") simuPOP::InheritTagger "

inherite tag from parents. If both parents have tags, use fathers.

";

%feature("docstring")  simuPOP::InheritTagger::InheritTagger " 

constructor. default to be always active. string can be any string
(m_Delimiter will be ignored for this class.) r will be replicate
number g will be generation number.

Usage:
  inheritTagger(begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::InheritTagger::~InheritTagger " 

Usage:
  x.~InheritTagger()

";

%feature("docstring")  simuPOP::InheritTagger::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::InheritTagger::applyDuringMating " 

give pop, offspring, pop and mom.

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring")  simuPOP::InheritTagger::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring") simuPOP::InitByFreq "

initialize genotype by allele frequency and sex by male frequency

";

%feature("docstring")  simuPOP::InitByFreq::InitByFreq " 

randomly assign alleles according to allele frequency

Usage:
  initByFreq(&alleleFreq=[], &alleleFreqs=[], identicalInds=False,
    &subPop=[], &indRange=vectori, indRanges=intMatrix, atLoci=[],
    maleFreq=0.5, stage=PreMating, begin=0, end=1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  alleleFreq:  an array of allele frequencies. Must add up to 1.
      
  alleleFreqs:  a matrix of allele frequencies, each row corresponse a
      subpopulation. alleleFreq will be ignored if this is specified.
      
  subPop:  an array of applicable subpopulations. default to all
      
  indRange:  NOT implemented yet.
      
  identicalInds:  whether or not make individual genotype identical in
      all subpopulation. If true, this operator will randomly generate
      genotype for an individual and spread it to the whole
      subpopulation.
      
  atLoci:  a vector of loci indices. If empty, apply to all loci
      
  maleFreq:  male frequency. Default to 0.5.
      
  stages:  is set to PreMating. Other parameters please see
      help(baseOperator.__init__)

Details:
  This operator randomly assign alleles according to given allele
  frequency. Allele frequencies can differ by subpop. Sex is also
  assigned randomly.

";

%feature("docstring")  simuPOP::InitByFreq::~InitByFreq " 

Usage:
  x.~InitByFreq()

";

%feature("docstring")  simuPOP::InitByFreq::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::InitByFreq::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::InitByFreq::apply " 

Usage:
  x.apply(&pop)

Details:
  initialize m_ranges

";

%feature("docstring") simuPOP::InitByValue "

initialize genotype by value and then copy to all individuals

";

%feature("docstring")  simuPOP::InitByValue::InitByValue " 

Usage:
  initByValue(value=vectori, values=intMatrix, atLoci=[], subPop=[],
    indRange=vectori, indRanges=intMatrix, &proportions=[],
    maleFreq=0.5, stage=PreMating, begin=0, end=1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::InitByValue::~InitByValue " 

Usage:
  x.~InitByValue()

";

%feature("docstring")  simuPOP::InitByValue::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::InitByValue::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::InitByValue::apply " 

Usage:
  x.apply(&pop)

Details:
  atLoci is in effect

";

%feature("docstring") simuPOP::Initializer "

initialize alleles at the start of generation.

Details:
  Bo Peng

";

%feature("docstring")  simuPOP::Initializer::Initializer " 

constructor. default to be always active.

Usage:
  initializer(&subPop=[], &indRange=vectori, &indRanges=intMatrix,
    &atLoci=[], maleFreq=0.5, stage=PreMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::Initializer::~Initializer " 

destructor

Usage:
  x.~Initializer()

";

%feature("docstring")  simuPOP::Initializer::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::Initializer::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::Initializer::setRanges " 

Usage:
  x.setRanges(&pop)

";

%feature("docstring") simuPOP::IOError "

exception, thrown if file io failure

";

%feature("docstring")  simuPOP::IOError::IOError " 

Usage:
  iOError(msg)

";

%ignore simuPOP::isAffected;

%feature("docstring")  simuPOP::isAffected::isAffected " 

Usage:
  isAffected()

";

%ignore simuPOP::isMale;

%feature("docstring")  simuPOP::isMale::isMale " 

Usage:
  isMale()

";

%feature("docstring") simuPOP::KAMMutator "

K-Allele Model mutator.

Details:
  Under this model, there are K (here refers as maxAllele) possible
  allele states, and any allele has a constant probability
  (rate/(K-1)) of mutation towards any of the K-1 allelic states.

Note:
  the theoretical mutation rate is rates/(K-1) towards any of the K-1
  allelic states. So rates is actually the probability to mutate!

See Also:
  Crow&Kimura 1970

";

%feature("docstring")  simuPOP::KAMMutator::KAMMutator " 

K-Allele Model mutator.

Usage:
  kAMMutator(rate=0, rates=[], atLoci=vectori, maxAllele=0,
    states=[], output=\">\", outputExpr=\"\", headers=[],
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  rate:  and rates are mutation rate. They are'p
      robability to mutate'. The actual
      mutation rate to any of the other K-1 allelic states are
      rates/(K-1)!
      
  states:  allelic states. Default to is 1 - maxAllele . Given states
      should be in order!
      
  atLoci:  and other parameters: refer to help(mutator),
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::KAMMutator::~KAMMutator " 

Usage:
  x.~KAMMutator()

";

%feature("docstring")  simuPOP::KAMMutator::mutate " 

mutate to a state other than current state with equal probability

Usage:
  x.mutate(*allele)

";

%feature("docstring")  simuPOP::KAMMutator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::KAMMutator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::MAPenetrance "

penetrance according to genotype at one locus

Details:
  multiple allele selector. This selector group alleles to disease and
  wild type and return penetrance to AA,Aa,aa. (A is wildtype).

";

%feature("docstring")  simuPOP::MAPenetrance::MAPenetrance " 

create a multiple allele selector (penetrance according to diseased
or wildtype alleles)

Usage:
  mAPenetrance(locus, &penetrance, &wildtype,
    exposePenetrance=False, stage=DuringMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  penetrance:  an array of penetrance of AA,Aa,aa. A is the wild type
      group.
      
  wildtype:  an array of alleles in the wildtype group. Anything else
      is disease allele.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::MAPenetrance::~MAPenetrance " 

Usage:
  x.~MAPenetrance()

";

%feature("docstring")  simuPOP::MAPenetrance::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::MAPenetrance::penet " 

currently assuming diploid

Usage:
  x.penet(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::MAPenetrance::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::MAQuanTrait "

quantitative trait according to genotype at one locus

Details:
  multiple allele selector. This selector group alleles to disease and
  wild type and return qtrait to AA,Aa,aa. (A is wildtype).

";

%feature("docstring")  simuPOP::MAQuanTrait::MAQuanTrait " 

create a multiple allele selector (quantitative trait according to
diseased or wildtype alleles)

Usage:
  mAQuanTrait(locus, &qtrait, &sigma, &wildtype, stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
    sep=\"\\\\t\")

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  qtrait:  an array of qtrait of AA,Aa,aa. A is the wild type group.
      
  sigma:  sigma for each of the trait genotype
      
  wildtype:  an array of alleles in the wildtype group. Anything else
      is disease allele.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::MAQuanTrait::~MAQuanTrait " 

destructor

Usage:
  x.~MAQuanTrait()

";

%feature("docstring")  simuPOP::MAQuanTrait::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::MAQuanTrait::qtrait " 

currently assuming diploid

Usage:
  x.qtrait(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::MAQuanTrait::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::MASelector "

selection according to genotype at one locus

Details:
  multiple allele selector. This selector group alleles to disease and
  wild type and return fitness to AA,Aa,aa. (A is wildtype).

";

%feature("docstring")  simuPOP::MASelector::MASelector " 

create a multiple allele selector (selection according to diseased
or wildtype alleles)

Usage:
  mASelector(locus, &fitness, &wildtype, stage=PreMating, begin=0,
    end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  fitness:  an array of fitness of AA,Aa,aa. A is the wild type group.
      
  wildtype:  an array of alleles in the wildtype group. Anything else
      is disease allele.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::MASelector::~MASelector " 

Usage:
  x.~MASelector()

";

%feature("docstring")  simuPOP::MASelector::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::MASelector::fitness " 

currently assuming diploid

Usage:
  x.fitness(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::MASelector::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::Mating "

Details:
  The Matingclasses describe various mating scheme --- a
  required parameter of Simulator.

";

%feature("docstring")  simuPOP::Mating::isCompatible " 

check if the mating type is compatible with population structure

Usage:
  x.isCompatible(&pop)

Details:
  possible things to check:need certain types of individual (age,
  sex etc)
  need resizeable population...

";

%feature("docstring")  simuPOP::Mating::Mating " 

constructor

Usage:
  mating(numOffsprings=1, newSubPopSize=[], newSubPopSizeExpr=\"\",
    *newSubPopSizeFunc=NULL)

";

%feature("docstring")  simuPOP::Mating::Mating " 

Usage:
  mating(&rhs)

";

%feature("docstring")  simuPOP::Mating::~Mating " 

destructor

Usage:
  x.~Mating()

";

%feature("docstring")  simuPOP::Mating::clone " 

 clone() const. Generate a copy of itself and return a pointer

Usage:
  x.clone()

Details:
  This function is important because Python automatically release an
  object after it is used.
  For example:will fail since Mating() is released after the first
  line being executed.
  With the help of clone() const, the C++
  implementation can avoid this problem by

";

%feature("docstring")  simuPOP::Mating::__repr__ " 

return name of the mating type used primarily in logging.

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::Mating::mate " 

mate: This is not supposed to be called for base Matingc
lass.

Usage:
  x.mate(&pop, &scratch, &ops)

Arguments:

  pop:   Population
      
  scratch:  scratch population
      
  ops:  during mating operators

Value:
  return false when mating fail.

";

%feature("docstring")  simuPOP::Mating::numOffsprings " 

Usage:
  x.numOffsprings()

";

%feature("docstring")  simuPOP::Mating::setNumOffsprings " 

Usage:
  x.setNumOffsprings(numOffsprings)

";

%feature("docstring")  simuPOP::Mating::formOffGenotype " 

whether or not to generate offspring genotype this is true when none
of the during-mating operator can do this.

Usage:
  x.formOffGenotype(&ops)

";

%feature("docstring")  simuPOP::Mating::prepareScratchPop " 

dealing with pop/subPop size change, copy of structure etc.

Usage:
  x.prepareScratchPop(&pop, &scratch)

Details:
  force new subPop size?

";

%feature("docstring") simuPOP::Migrator "

Details:
  Bo Peng

";

%feature("docstring")  simuPOP::Migrator::Migrator " 

create a migrator

Usage:
  migrator(&rates, mode=MigrByProbability, fromSubPop=[],
    toSubPop=[], stage=PreMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  rates:  migration rate, proportion or count. Determined by parameter
      mode.
      
  mode:  one of MigrByProbability (default), MigrByProportion or
      MigrByCounts
      
  fromSubPop:  an array of'f
      rom'subpops, default to all subpopulations
      
  toSubPop:  an array of't
      o'subpops, default to all subpopulations
      
  stage:  is default to PreMating. For details about other parameters,
      please refer to help(baseOperator.__init__)
      
  rates is a matrix with dimensions determined by fromSubPop and
  toSubPop. By default, rates is a matrix with element (i,j) being the
  migration rate, probability or count from subpop i to subpop j. If
  fromSubPop and/or toSubPop are given, migration only happen between
  these subpopulations. An extreme case is'p
  oint migration'
  rates=[[r]], fromSubPop=[a], toSubPop[b]
  which migrate from subpop a to b with given rate r.

";

%feature("docstring")  simuPOP::Migrator::~Migrator " 

destructor

Usage:
  x.~Migrator()

";

%feature("docstring")  simuPOP::Migrator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::Migrator::rates " 

return rates

Usage:
  x.rates()

";

%feature("docstring")  simuPOP::Migrator::setRates " 

set migration rates

Usage:
  x.setRates(&rates, mode)

Details:
  format 0-0 0-1 0-2, 1-0 1-1 1-2, 2-0, 2-1, 2-2. for mode 1 or 2,
  00,11,22 will be set automatically. regardless of input.
  set r[i][i]--- may need to extend rate (to add i->i
  )

";

%feature("docstring")  simuPOP::Migrator::apply " 

Usage:
  x.apply(&pop)

Details:
  2nd, or 3rd method
  create a vector and assign indices, then random shuffle and assign
  info
  for all subPop.

";

%feature("docstring")  simuPOP::Migrator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::MLPenetrance "

penetrance according to genotype at multiple loci multiplicative
model

Details:
  multiple loci selector. This selector takes several selectors and
  multiply their penetrance values... e.g. mlmPenetrance(
  [basicPenetrance(...), maPenetrance(...) ])

";

%feature("docstring")  simuPOP::MLPenetrance::MLPenetrance " 

multiple loci selector using a multiplicative model.

Usage:
  mLPenetrance(peneOps, mode=PEN_Multiplicative,
    exposePenetrance=False, stage=DuringMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  selectors:  a list of selectors.

";

%feature("docstring")  simuPOP::MLPenetrance::~MLPenetrance " 

Usage:
  x.~MLPenetrance()

";

%feature("docstring")  simuPOP::MLPenetrance::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::MLPenetrance::penet " 

currently assuming diploid

Usage:
  x.penet(*ind)

";

%feature("docstring")  simuPOP::MLPenetrance::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::MLQuanTrait "

quantitative trait according to genotype at multiple loci
multiplicative model

Details:
  multiple loci selector. This selector takes several selectors and
  multiply their qtrait values... e.g. mlmQuanTrait(
  [basicQuanTrait(...), maQuanTrait(...) ])

";

%feature("docstring")  simuPOP::MLQuanTrait::MLQuanTrait " 

multiple loci selector using a multiplicative model.

Usage:
  mLQuanTrait(qtraits, mode=QT_Multiplicative, stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
    sep=\"\\\\t\")

Arguments:

  qtraits:  a list of qtraits.

";

%feature("docstring")  simuPOP::MLQuanTrait::~MLQuanTrait " 

Usage:
  x.~MLQuanTrait()

";

%feature("docstring")  simuPOP::MLQuanTrait::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::MLQuanTrait::qtrait " 

currently assuming diploid

Usage:
  x.qtrait(*ind)

";

%feature("docstring")  simuPOP::MLQuanTrait::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::MLSelector "

selection according to genotype at multiple loci multiplicative
model

Details:
  multiple loci selector. This selector takes several selectors and
  multiply their fitness values... e.g. mlmSelector(
  [basicSelector(...), maSelector(...) ])

";

%feature("docstring")  simuPOP::MLSelector::MLSelector " 

multiple loci selector using a multiplicative model.

Usage:
  mLSelector(selectors, mode=SEL_Multiplicative, stage=PreMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
    sep=\"\\\\t\")

Arguments:

  selectors:  a list of selectors.

";

%feature("docstring")  simuPOP::MLSelector::~MLSelector " 

Usage:
  x.~MLSelector()

";

%feature("docstring")  simuPOP::MLSelector::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::MLSelector::fitness " 

currently assuming diploid

Usage:
  x.fitness(*ind)

";

%feature("docstring")  simuPOP::MLSelector::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::Mutator "

mutator class.

Details:
  Do not use this class directly. It just provide interface for real
  mutators.
  Every mutator can specify rate (equal rate) or rates (different rate
  for different loci) and a vector of applicable loci (default to all
  but should have the same length with rates if rates have length
  greater than one).
  max allele can be specified as well but more parameter, if needed,
  should be implemented by individual mutator classes.
  Number of possible allelic states: Most theoretical studies assume
  an infinite number of allelic states to avoid any homoplasy. If it
  facilitates analysis, this is however extremely unrealistic.
  Bo Peng

";

%feature("docstring")  simuPOP::Mutator::Mutator " 

create a mutator All mutators have the following common parameters.
However, the actual meaning of these parameters may vary according
to different model. Check the manual for details!!!
(help(kamMutator) for example.)

Usage:
  mutator(rate=-1, rates=[], atLoci=vectori, maxAllele=0,
    output=\">\", outputExpr=\"\", headers=[], stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
    sep=\"\\\\t\")

Arguments:

  rate:  single rate for all applicable loci (atLoci). Will be ignored
      if rates is specified.
      
  rates:  a vector of rates, the same length as atLoci.
      
  atLoci:  a vector of loci index. Can be ignored only when single
      rate is specified. Default to all loci.
      
  maxAllele:  max allowable allele. Interpreted by each sub mutaor
      class. Default to pop.maxAllele().

";

%feature("docstring")  simuPOP::Mutator::~Mutator " 

destructor

Usage:
  x.~Mutator()

";

%feature("docstring")  simuPOP::Mutator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::Mutator::rates " 

return mutation rate(s)

Usage:
  x.rates()

";

%feature("docstring")  simuPOP::Mutator::setRate " 

set an equal mutation rate

Usage:
  x.setRate(rate, atLoci=vectori)

";

%feature("docstring")  simuPOP::Mutator::setRates " 

set an array of rates

Usage:
  x.setRates(rates, atLoci=vectori)

";

%feature("docstring")  simuPOP::Mutator::maxAllele " 

return max allowable allele number

Usage:
  x.maxAllele()

";

%feature("docstring")  simuPOP::Mutator::setMaxAllele " 

Usage:
  x.setMaxAllele(maxAllele)

";

%feature("docstring")  simuPOP::Mutator::mutate " 

how to mutate a single allele. this is usually the only function
that need to be defined by the subclasses.

Usage:
  x.mutate(*allele)

";

%feature("docstring")  simuPOP::Mutator::apply " 

apply!

Usage:
  x.apply(&pop)

";

%feature("docstring") simuPOP::NoMating "

Details:
  No matin. parent generation will be the offspring generation during
  mating operators will be applied though.

";

%feature("docstring")  simuPOP::NoMating::NoMating " 

constructor, no new subPopsize parameter

Usage:
  noMating()

";

%feature("docstring")  simuPOP::NoMating::~NoMating " 

destructor

Usage:
  x.~NoMating()

";

%feature("docstring")  simuPOP::NoMating::clone " 

 clone() const. The same as Mating::clone() const.

Usage:
  x.clone()

See Also:
   Mating::clone() const

";

%feature("docstring")  simuPOP::NoMating::__repr__ " 

return name of the mating type

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::NoMating::mate " 

do the mating. --- no mating :-)

Usage:
  x.mate(&pop, &scratch, &ops)

Details:
  All individuals will be passed to during mating operators but no one
  will die (ignore during mating failing signal).

";

%feature("docstring") simuPOP::NoneOp "

";

%feature("docstring")  simuPOP::NoneOp::NoneOp " 

do nothing

Usage:
  noneOp(output=\">\", outputExpr=\"\", headers=[],
    stage=PostMating, begin=0, end=0, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::NoneOp::~NoneOp " 

destructor

Usage:
  x.~NoneOp()

";

%feature("docstring")  simuPOP::NoneOp::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::NoneOp::applyWithScratch " 

simply output some info providing interface to apply operator before
during or after mating.

Usage:
  x.applyWithScratch(&pop, &scratch, stage)

";

%feature("docstring")  simuPOP::NoneOp::applyDuringMating " 

give pop, offspring, pop and mom.

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring")  simuPOP::NoneOp::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::NoneOp::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::NullStreamBuf "

create a null stream buf that discard everything

";

%feature("docstring")  simuPOP::NullStreamBuf::NullStreamBuf " 

Usage:
  nullStreamBuf()

";

%feature("docstring")  simuPOP::NullStreamBuf::NullStreamBuf " 

Usage:
  nullStreamBuf()

";

%feature("docstring") simuPOP::Operator "

base class of all classes that manipulate populations.

Details:
  Operators are object that act on populations. They can be applied to
  populations directly using apply() member function, but most of the
  time they are managed and applied by a simulator.
  Operators can be applied at different stage(s) of a life cycle. More
  specifically, at pre-, duing- or post mating stage(s). Note that it
  is possible for an operator to apply multiple times in a life cycle.
  For example, an save-to-file operator might be applied before and
  after mating to trace parental information.
  Operators do not have to be applied at all generations. You can
  specify starting genertion, ending generation, gaps between
  applicable generations, or even specific generations to apply. For
  example, you might want to start applying migrations after certain
  heat-up generation; or you want to calculate every 10 generations.
  Operators can have outputs. Output can be standard output (terminal)
  or a file, which can be constant, or change with generation or
  replicate. Different operators can append to the same file to form
  table-like outputs. Headers of such table can be specified or
  automatically generated.
  filename can have the following format:
  1.'filename't
  his file will be closed after each use. I.e., if several operators
  output to the same file, only the last one will succeed.
  2.'>filename't
  he same as'filaname'
  3.'>>filename'T
  he file will be created at the beginning of evolution
  (simulator::evolve) and close at the end. Several operators can
  output to this file to form a table.
  4.'>>>f
  ilename'The same as'>>f
  ilename'except that the file will not be cleared at the
  beginning of evolution if it is not empty.
  5.'>'out put to
  standard output.
  6.''supress output.
  Most operators are applied to every replicate of a simulator during
  evolution. However, you can apply opertors to one or a group of
  replicates only. For example, you can initialize different
  replicates with different initial values and then start evolution.
  c.f. simulator::setGroup .
  Please refer to help(baseOperator) and help(baseOperator.__init__)
  for detailed information about member functions and parameters.
  Bo Peng

";

%feature("docstring")  simuPOP::Operator::Operator " 

create an operator (this function is not supposed to be called
directly)

Usage:
  operator(output, outputExpr, headers, stage, begin, end, step, at,
    rep, grp, sep)

Arguments:

  output:  a string of output filename. Different operators will have
      different default output (most commonly'>'o
      r'')
      
  outputExpr:  an expression that determines output filename
      dynamically.
      
  headers:  an array of headers for each item to output. Most
      operators provide default headers.
      
  begin:  start generation. default to 1. negative number is
      interpreted as endGeneration + begin
      
  end:  stop applying after this generation. negative number is
      allowed
      
  step:  number of generations between active generations. default to
      1
      
  at:  an array of active generations. If given, stage, begin, end,
      step will be ignored.
      
  rep:  applicable replicate. It can be replicate number 0 ~ (number
      of replicate -1), REP_ALL (all replicates) or REP_LAST (only to
      the last replicate). Usually default to REP_ALL.
      
  grp:  applicable group, default to GRP_ALL. A group number for each
      replicate is set by simulator.__init__ or simulator::setGroup().
      grp, if not GRP_ALL, will be compared to the group number of
      this replicate before applying.
      
  sep:  separator string. Usually default to''.
       Operators that have array like output should use this to
      separate output.

";

%feature("docstring")  simuPOP::Operator::~Operator " 

destroy an operator

Usage:
  x.~Operator()

";

%feature("docstring")  simuPOP::Operator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::Operator::isActive " 

judge if this operator is active

Usage:
  x.isActive(rep, numRep, gen, end, grp, repOnly=False)

Details:
  Judge if this operator is active under the conditions like current
  replicate, group number of current replicate, current generation.
  ending generation etc.

Note:
  will be called by simulator before applying.

";

%feature("docstring")  simuPOP::Operator::applicableGroup " 

return applicable group

Usage:
  x.applicableGroup()

";

%feature("docstring")  simuPOP::Operator::setApplicableGroup " 

set applicable group.

Usage:
  x.setApplicableGroup(grp=GRP_ALL)

Details:
  GRP_ALL is the default value (applicable to all groups. ) Otherwise,
  the operator is applicable to ONE group of replicates. groups can be
  set in Simulator::setGroup()

";

%feature("docstring")  simuPOP::Operator::applicableReplicate " 

return applicable replicate

Usage:
  x.applicableReplicate()

";

%feature("docstring")  simuPOP::Operator::setApplicableReplicate " 

set applicable replicate.

Usage:
  x.setApplicableReplicate(rep)

";

%feature("docstring")  simuPOP::Operator::setActiveGenerations " 

set applicable generation parrameters: stage, begin, end, step and
at

Usage:
  x.setActiveGenerations(begin=0, end=-1, step=1, at=[])

Details:
  set certain m_flags to speed up using this machanism.

";

%feature("docstring")  simuPOP::Operator::setApplicableStage " 

set m_stage settings. This is usually not usable since the m_stage
setting are set by default for each Operator.

Usage:
  x.setApplicableStage(stage)

";

%feature("docstring")  simuPOP::Operator::canApplyPreMating " 

Can this operator be applied pre-mating?

Usage:
  x.canApplyPreMating()

";

%feature("docstring")  simuPOP::Operator::canApplyDuringMating " 

Can this operator be applied uring mating?

Usage:
  x.canApplyDuringMating()

";

%feature("docstring")  simuPOP::Operator::canApplyPostMating " 

Can this operator be applied post-mating?

Usage:
  x.canApplyPostMating()

";

%feature("docstring")  simuPOP::Operator::canApplyPreOrPostMating " 

can be applied pre or post mating.

Usage:
  x.canApplyPreOrPostMating()

";

%ignore simuPOP::Operator::formOffGenotype();

%ignore simuPOP::Operator::setFormOffGenotype(bool flag=true);

%feature("docstring")  simuPOP::Operator::applyWithScratch " 

providing interface to apply operator before during or after mating.

Usage:
  x.applyWithScratch(&pop, &scratch, stage)

";

%feature("docstring")  simuPOP::Operator::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::Operator::applyDuringMating " 

give pop, offspring, pop and mom.

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring")  simuPOP::Operator::setOutput " 

ostream, if not set during construction.

Usage:
  x.setOutput(output=\"\", outputExpr=\"\")

";

%feature("docstring")  simuPOP::Operator::separator " 

what is the separator?

Usage:
  x.separator()

";

%feature("docstring")  simuPOP::Operator::setSeparator " 

set separator

Usage:
  x.setSeparator(sep)

";

%feature("docstring")  simuPOP::Operator::getOstream " 

get output stream. This function is not exposed to user.

Usage:
  x.getOstream(*dict=NULL, readable=False)

";

%feature("docstring")  simuPOP::Operator::closeOstream " 

close ostream and delete ostream pointer.. if it is a ofstream.

Usage:
  x.closeOstream()

";

%feature("docstring")  simuPOP::Operator::outputHeader " 

before evoluion, header function will be called to output header?
This should be controlled in Simulator::evolve.

Usage:
  x.outputHeader(&pop)

";

%feature("docstring")  simuPOP::Operator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::Operator::headers " 

return headers

Usage:
  x.headers()

";

%feature("docstring")  simuPOP::Operator::hasHeader " 

if headers are given

Usage:
  x.hasHeader()

";

%feature("docstring")  simuPOP::Operator::setHeaders " 

Usage:
  x.setHeaders(&headers)

";

%feature("docstring")  simuPOP::Operator::noOutput " 

if output=\">\".
 used internally

Usage:
  x.noOutput()

";

%ignore simuPOP::OstreamManager;

%feature("docstring")  simuPOP::OstreamManager::OstreamManager " 

OStream
Manager///////////////////////////////////////////////////////////.

Usage:
  ostreamManager()

";

%ignore simuPOP::OstreamManager::~OstreamManager();

%ignore simuPOP::OstreamManager::getOstream(const string &name, bool readable, bool realAppend, bool useString);

%ignore simuPOP::OstreamManager::hasOstream(const string &filename);

%ignore simuPOP::OstreamManager::listAll();

%feature("docstring")  simuPOP::OstreamManager::closeAll " 

close all files and clean the map

Usage:
  x.closeAll(closeAppend)

";

%feature("docstring") simuPOP::OutOfMemory "

exception, thrown if out of memory

";

%feature("docstring")  simuPOP::OutOfMemory::OutOfMemory " 

Usage:
  outOfMemory(msg)

";

%feature("docstring") simuPOP::Outputer "

 Outputeris a (special) subclass of Operatort
hat will output files with different format.

Details:
  Bo Peng

";

%feature("docstring")  simuPOP::Outputer::Outputer " 

constructor. default to be always active.

Usage:
  outputer(output=\">\", outputExpr=\"\", headers=[],
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::Outputer::~Outputer " 

destructor

Usage:
  x.~Outputer()

";

%feature("docstring")  simuPOP::Outputer::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring") simuPOP::OutputHelper "

";

%feature("docstring")  simuPOP::OutputHelper::OutputHelper " 

Usage:
  outputHelper(str=\"\\\\n\", output=\">\", outputExpr=\"\",
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::OutputHelper::apply " 

simply output some info

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::OutputHelper::outputHeader " 

Header.

Usage:
  x.outputHeader(&pop)

";

%feature("docstring")  simuPOP::OutputHelper::~OutputHelper " 

Usage:
  x.~OutputHelper()

";

%feature("docstring")  simuPOP::OutputHelper::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::OutputHelper::setString " 

set output string.

Usage:
  x.setString(str)

";

%feature("docstring")  simuPOP::OutputHelper::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::ParentsTagger "

inherite tag from parents. If both parents have tags, use fathers.

";

%feature("docstring")  simuPOP::ParentsTagger::ParentsTagger " 

constructor. default to be always active. string can be any string
(m_Delimiter will be ignored for this class.) r will be replicate
number g will be generation number.

Usage:
  parentsTagger(begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::ParentsTagger::~ParentsTagger " 

Usage:
  x.~ParentsTagger()

";

%feature("docstring")  simuPOP::ParentsTagger::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::ParentsTagger::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::ParentsTagger::applyDuringMating " 

give pop, offspring, pop and mom.

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring") simuPOP::Pause "

";

%feature("docstring")  simuPOP::Pause::Pause " 

stop simulation. press q to exit and any other key to continue

Usage:
  pause(prompt=True, stopOnKeyStroke=False, output=\">\",
    outputExpr=\"\", headers=[], stage=PostMating, begin=0, end=-1,
    step=1, at=[], rep=REP_LAST, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  prompt:  if true (default), print prompt message.
      
  stopOnKeyStroke:  if true, goon if no key was pressed

";

%feature("docstring")  simuPOP::Pause::~Pause " 

destructor

Usage:
  x.~Pause()

";

%feature("docstring")  simuPOP::Pause::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::Pause::apply " 

simply output some info

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::Pause::outputHeader " 

Header.

Usage:
  x.outputHeader(&pop)

";

%feature("docstring")  simuPOP::Pause::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::Penetrance "

penetrance

Details:
  Please refer to the user's guide for
  details.

";

%feature("docstring")  simuPOP::Penetrance::Penetrance " 

constructor. default to be always active. default to post mating

Usage:
  penetrance(exposePenetrance=False, stage=DuringMating, begin=0,
    end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::Penetrance::~Penetrance " 

destructor

Usage:
  x.~Penetrance()

";

%feature("docstring")  simuPOP::Penetrance::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::Penetrance::penet " 

calculate/return penetrance etc

Usage:
  x.penet(*ind)

";

%feature("docstring")  simuPOP::Penetrance::apply " 

set pentrance to all individuals and record penetrance if requested.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::Penetrance::applyDuringMating " 

set penetrance to all individual

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring")  simuPOP::Penetrance::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::Population "

a collection of individuals with subPopulation structure

Details:
  population objects are essential to simuPOP.
   They are composed of subpopulations each with certain number of
  individuals. Several important feature/limitation of population
  should be emphasized:
  1. Populations are composed of individuals of the same type. I.e.,
  all individuals have the same genotypic structure (number of loci,
  chromosomes etc).
  2. Populationobject has one layer of structure --
  subPopulations. A subPopulation of the whole population is implied
  if no structure is specified.
  3. Matingis within subpopulation. Migration can change
  subpopulation number and size but not the total population size. Mati
  ngcan change population and subpopulation size but no the number of
  subpopulations. These are rather reasonable limitations though).
  4. Although loci are separated by chromosomes, they are referred to
  mostly by their absolute indices, starting from 0. (So you do not
  have to specify chromosome number and locus index each time you need
  to refer to a locus.) Some utility functions are provided to convert
  between absolute and relative indices. Similar functions are
  provided to convert between absolute and relative individual
  indices.
  5. Populations have their own associated variables that can be
  retrieved through vars() function. pop.eval and pop.exec can
  manipuate variables inside.
  A population object provides the following functions. Please refer
  to member functions of population class for further details.
  1. access to genotypic structure, at both population level and
  individual level.
  2. access to individuals information and their genotypes.
  3. methods to manipulate individuals within/between subpopulations.
  4. methods to manupulated associated variables
  6. save/load population, in txt, bin, xml formats.

";

%feature("docstring")  simuPOP::Population::Population " 

create a population object with given size and genotypic structure

Usage:
  population(size=1, ploidy=2, &loci=[], &lociDist=[], &subPop=[],
    &alleleNames=[], maxAllele=MAXALLELE)

Arguments:

  size:  population size. Default to 1. It will be set to the sum of
      subPop if subPop is specified.
      
  ploidy:  number of sets of chromosomes. Default to 2 (diploid).
      
  loci:  an array of numbers of loci on each chromosome. If not
      specified, assume a single locus on one chromosome. Number of
      chromosomes is determined by the size of this array.
      
  lociDist:  a matrix of loci distance on each chromosome. An array of
      distance should be given for each chromosome in the form of (
      (loc1, loc2...), (loc1,...) ). The default values are 1, 2, etc,
      on each chromosome.
      
  subPop:  an array of subPopulation sizes. Default value is [size]
      which means a single subpopulation of the while population. If
      both size and subPop are given, sum of subPop should agree with
      size.
      
  alleleNames:  an array of allele names. The first element should be
      given to invalid/unknown allele. For example, for SNP loci with
      alleles A,C,T,G, you can specify alleleName as ('_','A','C','T','G')
      . Note that simuPOPuses 1,2,3,4 internally and
      these names will only be used during display.
      
  maxAllele:  maximum allele number. (Mutators can specify individual
      maxAllele for each locus.)

Value:
  no return value. Exceptionwill be thrown is
  wrong parameters are given.

See Also:
  simulator, baseOperator, mating schemes

Examples:
  >>> # a Wright-Fisher population
  ... WF = population(size=100, ploidy=1, loci=[1])
  >>> # a diploid population of size 10
  ... # there are two chromosomes with 5 and 7 loci respectively
  ... pop = population(size=10, ploidy=2, loci=[5, 7], subPop=[2, 8])
  >>> # a population with SNP markers (with names A,C,T,G
  ... #  range() are python functions
  ... pop = population(size=5, ploidy=2, loci=[5,10],
  ...     lociDist=[range(0,5),range(0,20,2)],
  ...     alleleNames=['_','A','C','T','G'],
  ...     subPop=[2,3], maxAllele=4)
  >>> InitByFreq(pop, [.25]*4)
  >>>  

";

%ignore simuPOP::Population::Population(const Population &rhs);

%feature("docstring")  simuPOP::Population::~Population " 

destroy a population

Usage:
  x.~Population()

";

%feature("docstring")  simuPOP::Population::__repr__ " 

Usage:
  x.__repr__()

";

%ignore simuPOP::Population::setSubPopulation(const vectorlu &subPopsize, bool allowPopSizeChange=true);

%ignore simuPOP::Population::initIndividuals();

%ignore simuPOP::Population::copySubPopStrucFrom(Population &rhs);

%feature("docstring")  simuPOP::Population::numSubPop " 

number of sub populations.

Usage:
  x.numSubPop()

Value:
  number of subPopulations (>=1)

";

%feature("docstring")  simuPOP::Population::subPopSize " 

get size of subPopulation subPop

Usage:
  x.subPopSize(subPop)

Arguments:

  subPop:  index of subPopulation (start from 0)

Value:
  size of subPopulation subPop

";

%feature("docstring")  simuPOP::Population::subPopSizes " 

get size of all subPopulations

Usage:
  x.subPopSizes()

Value:
  an array of size of subPopulations

";

%feature("docstring")  simuPOP::Population::popSize " 

get population size

Usage:
  x.popSize()

Value:
  total number of individuals in this population

";

%feature("docstring")  simuPOP::Population::absIndIndex " 

absolute index of individual at a subPopulation

Usage:
  x.absIndIndex(index, subPop)

Arguments:

  index:  index of individual at subPopulation subPop
      
  subPop:  subpopulation index

Value:
  absolute index of individual at subPopulation subPop

See Also:
   subPopIndPair

";

%feature("docstring")  simuPOP::Population::subPopIndPair " 

subPop and relative index of an individual

Usage:
  x.subPopIndPair(ind)

";

%feature("docstring")  simuPOP::Population::subPopBegin " 

beginning index for subPopulation subPop

Usage:
  x.subPopBegin(subPop)

Arguments:

  subPop:  subpopulation index

Details:

   absIndIndex(index, subPop)= subPopBeg
  in(subPop)+ index

Value:
  beginning index of this subpopulation

See Also:
   absIndIndex

";

%feature("docstring")  simuPOP::Population::subPopEnd " 

ending index for subPopulation subPop

Usage:
  x.subPopEnd(subPop)

Arguments:

  subPop:  subpopulation index

Value:
  ending index of this subpopulation (not in this subpop)

Note:
  as with all ...End functions, the returning index is out of the
  range so the actual range is [xxxBegin, xxxEnd). This agrees with
  all STL conventions.

See Also:
   absIndIndex

";

%feature("docstring")  simuPOP::Population::individual " 

reference to individual ind

Usage:
  x.individual(ind)

Arguments:

  ind:  absolute index of an individual

Value:
  reference to an individual

";

%feature("docstring")  simuPOP::Population::individual " 

refernce to individual ind in subPopulation subPop

Usage:
  x.individual(ind, subPop)

Arguments:

  ind:  individual index within subPop
      
  subPop:  subPopulation index

Value:
  reference to an individual

";

%ignore simuPOP::Population::indBegin();

%ignore simuPOP::Population::indEnd();

%ignore simuPOP::Population::indBegin(UINT subPop);

%ignore simuPOP::Population::indEnd(UINT subPop);

%ignore simuPOP::Population::alleleBegin(UINT locus);

%ignore simuPOP::Population::alleleEnd(UINT locus);

%ignore simuPOP::Population::alleleBegin(UINT locus, UINT subPop);

%ignore simuPOP::Population::alleleEnd(UINT locus, UINT subPop);

%ignore simuPOP::Population::begin();

%ignore simuPOP::Population::end();

%ignore simuPOP::Population::begin(UINT subPop);

%ignore simuPOP::Population::end(UINT subPop);

%ignore simuPOP::Population::genoBegin(ULONG ind);

%ignore simuPOP::Population::genoEnd(ULONG ind);

%ignore simuPOP::Population::genoBegin(ULONG ind, UINT subPop);

%ignore simuPOP::Population::genoEnd(ULONG ind, UINT subPop);

%feature("docstring")  simuPOP::Population::arrGenotype " 

get the whole genotype. individuals will be in order before exposing
their genotypes.

Usage:
  x.arrGenotype()

";

%feature("docstring")  simuPOP::Population::arrGenotype " 

get the whole genotype. individuals will be in order before exposing
their genotypes.

Usage:
  x.arrGenotype(subPop)

";

%feature("docstring")  simuPOP::Population::exposeInfo " 

Usage:
  x.exposeInfo()

Details:
  brief return individual info in pop namespace

";

%feature("docstring")  simuPOP::Population::exposeAffectedness " 

Usage:
  x.exposeAffectedness()

Details:
  brief return individual affected status in pop namespace

";

%feature("docstring")  simuPOP::Population::setInfo " 

Usage:
  x.setInfo(info)

Arguments:

  info:  an array of info values, should have length of pop size

Details:
  This function set info field of all individuals. Info can be used to
  sort individuals and set subpopulations. Therefore, the following
  code will do a migration:
  setInfo([ an_array_of_info_value ])
  rearrangeIndByInfo()
   setSubPopByIndInfo()
  or ignore setInfo by running:
  rearrangeIndByInfo([ an_array_of_info_value ])
   setSubPopByIndInfo()

See Also:
   Individual::setInfo, Individual::info,
   info

";

%feature("docstring")  simuPOP::Population::setIndInfoWithSubPopID " 

set individual info with their subPopulation id.

Usage:
  x.setIndInfoWithSubPopID()

Details:
  set individual info by subpop id.

";

%feature("docstring")  simuPOP::Population::setSubPopByIndInfo " 

adjust subPopulation according to individual info values

Usage:
  x.setSubPopByIndInfo()

Details:
  assume individual has subpop index in their info value and put them
  into corresponding subpopulations. This function assume that
  individuals are already sorted by their info value so it will not
  move individuals around.

  rebuild index

Note:
  Assume individuals are ordered by their info value!

See Also:
   setInfo, rearrangeByIndInfo

";

%feature("docstring")  simuPOP::Population::rearrangeByIndInfo " 

rearrange individual to other subPops

Usage:
  x.rearrangeByIndInfo(info=vectori)

Arguments:

  info:  if not empty, set info to all individuals c.f. setInf
      o(info)

Details:
  rearrange individuals by their info value or given info

";

%feature("docstring")  simuPOP::Population::shrinkByIndInfo " 

Usage:
  x.shrinkByIndInfo(info=vectori)

Details:
  resize population according to info

";

%feature("docstring")  simuPOP::Population::newPopByIndInfo " 

Usage:
  x.newPopByIndInfo(info=vectori)

Details:
  form a new population according to info

";

%ignore simuPOP::Population::swap(ULONG id1, ULONG id2, bool swapContent=true);

%feature("docstring")  simuPOP::Population::equalTo " 

compare two populations

Usage:
  x.equalTo(&rhs)

";

%ignore simuPOP::Population::isActive();

%ignore simuPOP::Population::activate();

%ignore simuPOP::Population::deactivate();

%ignore simuPOP::Population::adjustGenoPosition(bool deep=false);

%feature("docstring")  simuPOP::Population::savePopulation " 

save population to a file

Usage:
  x.savePopulation(&filename, &format=\"auto\")

Arguments:

  filename:  save to filename
      
  format:  format to save. Can be one of't
      ext','b
      in','x
      ml'
      
  The default format is'text'b
  ut the output is not suppored to be read.'b
  in'has smaller size and should be used for large populations.'x
  ml'format is most readable and should be used when you would
  like to convert simuPOPpopulations to other formats.

See Also:
  global function loadPopulation

";

%ignore simuPOP::Population::loadPopulation(const string &filename, const string &format="auto");

%feature("docstring")  simuPOP::Population::rep " 

Usage:
  x.rep()

";

%ignore simuPOP::Population::setRep(int rep, bool setVar=true);

%feature("docstring")  simuPOP::Population::grp " 

Usage:
  x.grp()

";

%ignore simuPOP::Population::setGrp(UINT grp, bool setVar=true);

%feature("docstring")  simuPOP::Population::vars " 

return variables of this population if subPop is given, return
dictionary for specified subpopulation.

Usage:
  x.vars(subPop=-1)

";

%ignore simuPOP::Population::dict(int subPop=-1);

%ignore simuPOP::Population::setDict(PyObject *dict);

%ignore simuPOP::Population::hasVar(const string &name);

%ignore simuPOP::Population::removeVar(const string &name);

%ignore simuPOP::Population::setBoolVar(const string &name, const bool val);

%ignore simuPOP::Population::setIntVar(const string &name, const int val);

%ignore simuPOP::Population::setDoubleVar(const string &name, const double val);

%ignore simuPOP::Population::setStringVar(const string &name, const string &val);

%ignore simuPOP::Population::setStrDictVar(const string &name, const strDict &val);

%ignore simuPOP::Population::setIntDictVar(const string &name, const intDict &val);

%ignore simuPOP::Population::setVar(const string &name, PyObject *val);

%ignore simuPOP::Population::setDoubleNumArrayVar(const string &name, int dim, double *buf, bool copyOver=true);

%ignore simuPOP::Population::setIntNumArrayVar(const string &name, int dim, int *buf, bool copyOver=true);

%ignore simuPOP::Population::getVar(const string &name, bool nameError=true);

%ignore simuPOP::Population::getVarAsBool(const string &name, bool nameError=true);

%ignore simuPOP::Population::getVarAsInt(const string &name, bool nameError=true);

%ignore simuPOP::Population::getVarAsDouble(const string &name, bool nameError=true);

%ignore simuPOP::Population::getVarAsString(const string &name, bool nameError=true);

%ignore simuPOP::Population::getVarAsStrDict(const string &name, bool nameError=true);

%ignore simuPOP::Population::getVarAsIntDict(const string &name, bool nameError=true);

%ignore simuPOP::Population::getVarAsDoubleNumArray(const string &name, double *&buf, bool nameError=true);

%ignore simuPOP::Population::getVarAsIntNumArray(const string &name, int *&buf, bool nameError=true);

%feature("docstring")  simuPOP::Population::evaluate " 

evaluate python statment/expressions

Usage:
  x.evaluate(&expr=\"\", &stmts=\"\")

Details:
  this function evaluate python expressions and return as string
  representing the result

";

%feature("docstring")  simuPOP::Population::execute " 

Usage:
  x.execute(&stmts=\"\")

";

%feature("docstring") simuPOP::PyEval "

evaluate an expression.

";

%feature("docstring")  simuPOP::PyEval::PyEval " 

evaluate expr/statments in local replicate namespace

Usage:
  pyEval(&expr=\"\", &stmts=\"\", &preStmts=\"\", &postStmts=\"\",
    exposePop=False, &name=\"\", output=\">\", outputExpr=\"\",
    headers=[], stage=PostMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  expr:  expression to be evaluated. Its result will be sent to
      output.
      
  stmts:  statement that will be executed before the expression.
      
  preStmts:  statement that will be executed when the operaot is
      constructed.
      
  postStmts:  statement that will be executed when the operator is
      destroyed.
      
  exposePop:  if true, expose current pop as variable\"p
      op\"
      
  name:  used to let pure python operator to identify themselves
      
  output:  default to\">\". i.e., output to standard output. for usage of other
      parameters, see help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::PyEval::~PyEval " 

Usage:
  x.~PyEval()

";

%feature("docstring")  simuPOP::PyEval::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::PyEval::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::PyEval::outputHeader " 

before evoluion, header function will be called to output header?
This should be controlled in Simulator::evolve.

Usage:
  x.outputHeader(&pop)

";

%feature("docstring")  simuPOP::PyEval::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::PyEval::name " 

Usage:
  x.name()

";

%feature("docstring") simuPOP::PyExec "

evaluate an expression.

";

%feature("docstring")  simuPOP::PyExec::PyExec " 

evaluate statments in local replicate namespace, no return value

Usage:
  pyExec(&stmts=\"\", &preStmts=\"\", &postStmts=\"\",
    exposePop=False, &name=\"\", output=\">\", outputExpr=\"\",
    headers=[], stage=PostMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  stmts:  statements (a single or multi-line string) that will be
      evaluated before the expression.
      
  preStmts:  statement that will be executed when the operaot is
      constructed.
      
  postStmts:  statement that will be executed when the operator is
      destroyed.
      
  exposePop:  if true, expose current pop as variable\"p
      op\"
      
  output:  default to\">\". i.e., output to standard output. for usage of other
      parameters, see help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::PyExec::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::PyExec::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::PyInit "

";

%feature("docstring")  simuPOP::PyInit::PyInit " 

initialize populations using given user function.

Usage:
  pyInit(*func, subPop=[], atLoci=[], indRange=vectori,
    indRanges=intMatrix, maleFreq=0.5, stage=PreMating, begin=0,
    end=1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  func:  a python function with parameter (index, ploidy, subpop)
      index is the allele index (0 ~ totNumLoci()-1), ploidy (index to
      copy of chromosomes), subpop (subpop number). The return value
      of this function should be a integer.
      
  stage:  is et to PreMating. Other parameters please refer to
      help(baseOperator.__init__)

Details:
  User of this operator must supply a Python function with parameter
  (index, ploidy, subpop). This operator will loop through all
  individual in each subpop and call this function to initialize
  populations.
  The arrange of parameters allows different initialization scheme for
  each subpop.

";

%feature("docstring")  simuPOP::PyInit::~PyInit " 

Usage:
  x.~PyInit()

";

%ignore simuPOP::PyInit::PyInit(const PyInit &rhs);

%feature("docstring")  simuPOP::PyInit::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::PyInit::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::PyInit::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring") simuPOP::PyMigrator "

migrate using given info vector

Details:
  You can use directMigrator to accomplish any migration: that is to
  say you directly specify subpopulation numbers for each individual
  and this operator will do the rest.

";

%feature("docstring")  simuPOP::PyMigrator::PyMigrator " 

create a directMigrator

Usage:
  pyMigrator(*subPopID=NULL, stage=PreMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  subPopID:  a 1-d Numeric Python int array. Must has length greater
      or equal to population size.
      
  stage:  is default to PreMating, please refer to
      help(baseOperator.__init__) for details about other parameters.

Details:
  This operator accept a one-dimensional Numeric Python int array.
  (created by Numeric.array ). The contend of the array will be
  considered as subpopulation id.

";

%feature("docstring")  simuPOP::PyMigrator::~PyMigrator " 

destructor

Usage:
  x.~PyMigrator()

";

%ignore simuPOP::PyMigrator::PyMigrator(const PyMigrator &rhs);

%feature("docstring")  simuPOP::PyMigrator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::PyMigrator::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::PyMigrator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::PyMutator "

mixed mutation model . has not been implemented.

";

%feature("docstring")  simuPOP::PyMutator::PyMutator " 

Usage:
  pyMutator(rate=0, rates=[], atLoci=vectori, maxAllele=0,
    *func=NULL, output=\">\", outputExpr=\"\", headers=[],
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::PyMutator::~PyMutator " 

Usage:
  x.~PyMutator()

";

%ignore simuPOP::PyMutator::PyMutator(const PyMutator &rhs);

%feature("docstring")  simuPOP::PyMutator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::PyMutator::mutate " 

how to mutate a single allele. this is usually the only function
that need to be defined by the subclasses.

Usage:
  x.mutate(*allele)

";

%feature("docstring")  simuPOP::PyMutator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::PyPenetrance "

penetrance using user supplied function

Details:
  Assign penetrance value by calling a user supplied function

";

%feature("docstring")  simuPOP::PyPenetrance::PyPenetrance " 

provide locus and penetrance for 11, 12, 13 (in the form of
dictionary)

Usage:
  pyPenetrance(loci, *func, exposePenetrance=False,
    stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  loci:  susceptibility loci. The genotype at these loci will be
      passed to func.
      
  func:  a Python function that accept genotypes at susceptibility
      loci and return penetrance value.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::PyPenetrance::~PyPenetrance " 

destructor

Usage:
  x.~PyPenetrance()

";

%ignore simuPOP::PyPenetrance::PyPenetrance(const PyPenetrance &rhs);

%feature("docstring")  simuPOP::PyPenetrance::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::PyPenetrance::penet " 

currently assuming diploid

Usage:
  x.penet(*ind)

";

%feature("docstring")  simuPOP::PyPenetrance::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::PyQuanTrait "

quantitative trait using user supplied function

Details:
  Assign qtrait value by calling a user supplied function

";

%feature("docstring")  simuPOP::PyQuanTrait::PyQuanTrait " 

provide locus and qtrait for 11, 12, 13 (in the form of dictionary)

Usage:
  pyQuanTrait(loci, *func, stage=PostMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  loci:  susceptibility loci. The genotype at these loci will be
      passed to func.
      
  func:  a Python function that accept genotypes at susceptibility
      loci and return qtrait value.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::PyQuanTrait::~PyQuanTrait " 

Usage:
  x.~PyQuanTrait()

";

%ignore simuPOP::PyQuanTrait::PyQuanTrait(const PyQuanTrait &rhs);

%feature("docstring")  simuPOP::PyQuanTrait::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::PyQuanTrait::qtrait " 

currently assuming diploid

Usage:
  x.qtrait(*ind)

";

%feature("docstring")  simuPOP::PyQuanTrait::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::PySample "

thrink population accroding to some outside value

";

%feature("docstring")  simuPOP::PySample::PySample " 

create a python sampler

Usage:
  pySample(*keep, &name=\"sample\", &nameExpr=\"\", times=1,
    &saveAs=\"\", &saveAsExpr=\"\", &format=\"bin\", stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
    sep=\"\\\\t\")

Arguments:

  keep:  a carray of the length of population. its values will be
      assigned to info.  and other parameters please see
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::PySample::~PySample " 

destructor

Usage:
  x.~PySample()

";

%ignore simuPOP::PySample::PySample(const PySample &rhs);

%feature("docstring")  simuPOP::PySample::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::PySample::drawSample " 

Usage:
  x.drawSample(&pop)

";

%feature("docstring")  simuPOP::PySample::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::PySelector "

selection using user supplied function

Details:
  Assign fitness value by calling a user supplied function

";

%feature("docstring")  simuPOP::PySelector::PySelector " 

provide locus and fitness for 11, 12, 13 (in the form of dictionary)

Usage:
  pySelector(loci, *func, stage=PreMating, begin=0, end=-1, step=1,
    at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  loci:  susceptibility loci. The genotype at these loci will be
      passed to func.
      
  func:  a Python function that accept genotypes at susceptibility
      loci and return fitness value.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::PySelector::~PySelector " 

destructor

Usage:
  x.~PySelector()

";

%ignore simuPOP::PySelector::PySelector(const PySelector &rhs);

%feature("docstring")  simuPOP::PySelector::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::PySelector::fitness " 

currently assuming diploid

Usage:
  x.fitness(*ind)

";

%feature("docstring")  simuPOP::PySelector::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::PySubset "

thrink population accroding to some outside value

";

%feature("docstring")  simuPOP::PySubset::PySubset " 

create a directMigrator

Usage:
  pySubset(*keep=NULL, stage=PostMating, begin=0, end=-1, step=1,
    at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  keep:  a carray of the length of population. its values will be
      assigned to info.  and other parameters please see
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::PySubset::~PySubset " 

destructor

Usage:
  x.~PySubset()

";

%ignore simuPOP::PySubset::PySubset(const PySubset &rhs);

%feature("docstring")  simuPOP::PySubset::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::PySubset::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::PySubset::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::PythonCoutBuf "

create a stream buf that print to python sys.stdout cout will be
redirected to this buf to really output to python console.

";

%feature("docstring")  simuPOP::PythonCoutBuf::PythonCoutBuf " 

Usage:
  pythonCoutBuf()

";

%feature("docstring")  simuPOP::PythonCoutBuf::PythonCoutBuf " 

Usage:
  pythonCoutBuf()

";

%feature("docstring") simuPOP::QuanTrait "

quantitative trait

Details:
  Genetic quantitative trait is tricky to simulate. In simuPOP,
   I employee an ability (fitness) to mate approach. Namely, the
  probability that an individual will be chosen for mating is
  proportional to its fitness value. More specifically,
  PreMating selectors assign fitness values to each individual.
  
  Sexless mating (e.g. binomialSelection) : individuals are chosen at
  probabilities that are proportional to their fitness values. More
  specifically, if there are N individuals with fitness values $f_i, i=
  1,...,N $, individual $i$w
  ill have probability $ \\\\frac{f_i}{\\\\sum_{j=1}^N f_j} $t
  o be chosen to be passed to the next generation.
  
  Random mating with sex (e.g. randomMating): males and females are
  separated and each are chosen as described above.
  
  Please refer to the user's guide for
  details.

";

%feature("docstring")  simuPOP::QuanTrait::QuanTrait " 

constructor. default to be always active.

Usage:
  quanTrait(stage=PostMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::QuanTrait::~QuanTrait " 

destructor

Usage:
  x.~QuanTrait()

";

%feature("docstring")  simuPOP::QuanTrait::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::QuanTrait::qtrait " 

calculate/return quantitative trait etc

Usage:
  x.qtrait(*ind)

";

%feature("docstring")  simuPOP::QuanTrait::apply " 

set qtrait to all individual

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::QuanTrait::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::RandomMating "

Details:
  basic sexual random mating.
  Within each subPopulation, choose male and female randomly randmly
  get one copy of chromosome from father/mother.
  require: sexed individual; ploidy == 2
  apply during mating operators and put into the next generation.
  if ignoreParentsSex is set, parents will be chosen regardless of
  sex.
  Otherwise, male and female will be collected and be chosen randomly.
  If there is no male or female in a subPopulation, if
  m_UseSameSexIfUniSex is true, an warning will be generated and same
  sex mating (?) will be used otherwise, RandomMatingw
  ill return false.
  if there is no during mating operator to copy alleles, a direct copy
  will be used.

";

%feature("docstring")  simuPOP::RandomMating::RandomMating " 

constructor

Usage:
  randomMating(numOffsprings=1, newSubPopSize=[],
    *newSubPopSizeFunc=NULL, newSubPopSizeExpr=\"\",
    contWhenUniSex=True)

";

%feature("docstring")  simuPOP::RandomMating::~RandomMating " 

destructor

Usage:
  x.~RandomMating()

";

%feature("docstring")  simuPOP::RandomMating::clone " 

 clone() const. Generate a copy of itself and return pointer this
is to make sure the object is persistent and will not be freed by
python.

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::RandomMating::isCompatible " 

check if the mating type is compatible with population structure

Usage:
  x.isCompatible(&pop)

Details:
  possible things to check:need certain types of individual (age,
  sex etc)
  need resizeable population...

";

%feature("docstring")  simuPOP::RandomMating::__repr__ " 

return name of the mating type

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::RandomMating::mate " 

do the mating. parameters see Mating::mate.

Usage:
  x.mate(&pop, &scratch, &ops)

Details:
  Within each subPopulation, choose male and female randomly randmly
  get one copy of chromosome from father/mother.
  require: sexed individual; ploidy == 2
  apply during mating operators and put into the next generation.
  Otherwise, male and female will be collected and be chosen randomly.
  If there is no male or female in a subPopulation,
  if m_contWhenUniSex is true, an warning will be generated and same
  sex mating (?) will be used
  otherwise, RandomMatingwill return false.
  
  determine if mate() will generate offspring genotype
  if selection is on
  random mating happens within each subPopulation
  now, all individuals of needToFind sex is collected
  if selection is on
  normalize and cusum fitness

";

%feature("docstring") simuPOP::RandomSample "

thrink population accroding to some outside value

";

%feature("docstring")  simuPOP::RandomSample::RandomSample " 

draw random sample, regardless of affected status

Usage:
  randomSample(size=0, &name=\"sample\", &nameExpr=\"\", times=1,
    &saveAs=\"\", &saveAsExpr=\"\", &format=\"bin\", stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
    sep=\"\\\\t\")

Arguments:

  size:  size of sample
      
  stage:  and other parameters please see help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::RandomSample::~RandomSample " 

destructor

Usage:
  x.~RandomSample()

";

%feature("docstring")  simuPOP::RandomSample::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::RandomSample::drawSample " 

Usage:
  x.drawSample(&pop)

";

%feature("docstring")  simuPOP::RandomSample::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::Recombinator "

initialize alleles at the start of generation.

Details:
  objects Recombinator
  only works for diploids (and for females in haplodiploids) Populati
  on.
  
  Free recombination between loci. Loci behave completely
  independently.
  
  otherwise there will be some linkage between loci, user need to
  specify physical recombination rate between adjacent loci (ie
  between locus n and n+1)
  
  The recombination rate must be comprised between 0.0 and 1.0.
  
  A recombination rate of 0.0 means that the loci are completely
  linked and thus behave together as a single linked locus.
  
  A recombination rate of 1.0 is equivalent to free recombination.
  
  All values in between will represent various linkage intensities
  between adjacent pairs of loci. The recombination rate is equivalent
  to 1-linkage and represents the probability that the allele at the
  next locus is randomly drawn.

";

%feature("docstring")  simuPOP::Recombinator::Recombinator " 

recombine chromosomes from parents

Usage:
  recombinator(rate=-1, rates=[], afterLoci=[], begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  rate:  recombination rate per unit. I.e., the really recombination
      rate between two loci is determined by rate and the loci
      distance between them.
      
  rates:  an array of recombination rates. Must be the same length as
      afterLoci or totNumOfLoci(). The last item can be ignored.
      
  afterLoci:  an array of loci index. If rates is also specified, they
      should have the same length. Default to all loci (but
      meaningless for those loci locate at the end of chromosome.)

";

%feature("docstring")  simuPOP::Recombinator::~Recombinator " 

Usage:
  x.~Recombinator()

";

%feature("docstring")  simuPOP::Recombinator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::Recombinator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::Recombinator::prepareRecRates " 

clculate rates

Usage:
  x.prepareRecRates(&pop)

Details:
  get loci distance * rate and then recombinant points
  get loci distance * rate and then recombinant points

";

%feature("docstring")  simuPOP::Recombinator::applyDuringMating " 

give pop, offspring, pop and mom.

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring") simuPOP::RNG "

";

%feature("docstring")  simuPOP::RNG::RNG " 

Random number generator/////////////////////////////////////////////

Usage:
  rNG(*rng=NULL, seed=0)

";

%feature("docstring")  simuPOP::RNG::~RNG " 

Usage:
  x.~RNG()

";

%feature("docstring")  simuPOP::RNG::setRNG " 

choose an random number generator. This can be done by setting
GSL_RNG_TYPE as well.

Usage:
  x.setRNG(*rng=NULL, seed=0)

";

%feature("docstring")  simuPOP::RNG::name " 

return RNGname

Usage:
  x.name()

";

%feature("docstring")  simuPOP::RNG::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::RNG::randInt " 

return [0, 1, 2, ... n-1]

Usage:
  x.randInt(n)

";

%feature("docstring")  simuPOP::RNG::randIntArray " 

return [0, 1, 2, ... n-1]

Usage:
  x.randIntArray(n, size, *vec)

";

%feature("docstring")  simuPOP::RNG::randGeometric " 

return geometric with parameter p

Usage:
  x.randGeometric(p)

";

%feature("docstring")  simuPOP::RNG::randUniform01 " 

return double [0,1)

Usage:
  x.randUniform01()

";

%feature("docstring")  simuPOP::RNG::randNormal " 

return double -inf, inf, v is standard deviation

Usage:
  x.randNormal(m, v)

";

%feature("docstring")  simuPOP::RNG::randUniform01Array " 

return double [0,1)

Usage:
  x.randUniform01Array(size, *vec)

";

%feature("docstring")  simuPOP::RNG::randIntByFreq " 

return [0, n-1] according to freq freq should be cumulative [p_0,
p_1, ... p_n-2, p_n-1 ] p_(n-1) should be 1. p_n-1 should exist
since we then can deal with the case when n=1, p0=1. the new
implementation uses binary search instead of linear search to
imporve efficiency (less than 4 comparison)

Usage:
  x.randIntByFreq(n, *freq)

";

%feature("docstring")  simuPOP::RNG::randIntByFreqArray " 

if fast, freq should not be incremental

Usage:
  x.randIntByFreqArray(n, *freq, beg, end, offset=0)

";

%feature("docstring")  simuPOP::RNG::randBinomial " 

binomial distribution

Usage:
  x.randBinomial(n, p)

";

%feature("docstring") simuPOP::Sample "

sample from population and save samples sample operator will
generate a new subpopulation in pop namespace.

";

%feature("docstring")  simuPOP::Sample::Sample " 

create a sample

Usage:
  sample(&name=\"sample\", &nameExpr=\"\", times=1, &saveAs=\"\",
    &saveAsExpr=\"\", &format=\"bin\", stage=PostMating, begin=0,
    end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  name:  variable name of the sampled population (will be put in pop
      local namespace)
      
  nameExpr:  expression version of name. If both name and nameExpr is
      empty, do not store pop.
      
  times:  how many times to run the sample process? This is usually
      one, but we may want to take several random samples.
      
  saveAs:  filename to save the population.
      
  saveAsExpr:  expression for save filename
      
  format:  to save sample(s)
      
  stage:  and other parameters please see help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::Sample::~Sample " 

destructor

Usage:
  x.~Sample()

";

%feature("docstring")  simuPOP::Sample::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::Sample::drawSample " 

Usage:
  x.drawSample(&pop)

";

%feature("docstring")  simuPOP::Sample::samples " 

return the samples

Usage:
  x.samples(&pop)

";

%feature("docstring")  simuPOP::Sample::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::Sample::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::SavePopulation "

save population to a file

";

%feature("docstring")  simuPOP::SavePopulation::SavePopulation " 

Usage:
  savePopulation(output=\"\", outputExpr=\"\", headers=[],
    format=\"bin\", stage=PostMating, begin=-1, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::SavePopulation::~SavePopulation " 

Usage:
  x.~SavePopulation()

";

%feature("docstring")  simuPOP::SavePopulation::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::SavePopulation::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::SavePopulation::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::Selector "

selection

Details:
  Genetic selection is tricky to simulate. In simuPOP,
   I employee an ability (fitness) to mate approach. Namely, the
  probability that an individual will be chosen for mating is
  proportional to its fitness value. More specifically,
  PreMating selectors assign fitness values to each individual.
  
  Sexless mating (e.g. binomialSelection) : individuals are chosen at
  probabilities that are proportional to their fitness values. More
  specifically, if there are N individuals with fitness values $f_i, i=
  1,...,N $, individual $i$w
  ill have probability $ \\\\frac{f_i}{\\\\sum_{j=1}^N f_j} $t
  o be chosen to be passed to the next generation.
  
  Random mating with sex (e.g. randomMating): males and females are
  separated and each are chosen as described above.
  
  Please refer to the user's guide for
  details.

";

%feature("docstring")  simuPOP::Selector::Selector " 

constructor. default to be always active.

Usage:
  selector(stage=PreMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::Selector::~Selector " 

destructor

Usage:
  x.~Selector()

";

%feature("docstring")  simuPOP::Selector::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::Selector::fitness " 

calculate/return w11 etc

Usage:
  x.fitness(*ind)

";

%feature("docstring")  simuPOP::Selector::apply " 

set fitness to all individual

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::Selector::__repr__ " 

Usage:
  x.__repr__()

";

%ignore simuPOP::SharedVariables;

%ignore simuPOP::SharedVariables::SharedVariables();

%ignore simuPOP::SharedVariables::SharedVariables(PyObject *dict, bool ownVars);

%ignore simuPOP::SharedVariables::~SharedVariables();

%feature("docstring")  simuPOP::SharedVariables::clear " 

Usage:
  x.clear()

";

%feature("docstring")  simuPOP::SharedVariables::setDict " 

Usage:
  x.setDict(*dict)

";

%feature("docstring")  simuPOP::SharedVariables::setVar " 

setvars C++ ==>Python

Usage:
  x.setVar(&name, *val)

Details:
  if the size if enough, get the item

";

%ignore simuPOP::SharedVariables::getVar(const string &name, bool nameError=true);

%feature("docstring")  simuPOP::SharedVariables::hasVar " 

Usage:
  x.hasVar(&name)

";

%feature("docstring")  simuPOP::SharedVariables::removeVar " 

remove variable

Usage:
  x.removeVar(&name)

";

%ignore simuPOP::SharedVariables::setBoolVar(const string &name, const bool val);

%ignore simuPOP::SharedVariables::setIntVar(const string &name, const int val);

%ignore simuPOP::SharedVariables::setDoubleVar(const string &name, const double val);

%ignore simuPOP::SharedVariables::setStringVar(const string &name, const string &val);

%ignore simuPOP::SharedVariables::setStrDictVar(const string &name, const strDict &val);

%ignore simuPOP::SharedVariables::setIntDictVar(const string &name, const intDict &val);

%ignore simuPOP::SharedVariables::setDoubleNumArrayVar(const string &name, int dim, double *buf, bool copyOver=true);

%ignore simuPOP::SharedVariables::setIntNumArrayVar(const string &name, int dim, int *buf, bool copyOver=true);

%ignore simuPOP::SharedVariables::getVarAsBool(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsInt(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsDouble(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::numArrayOwnMem(const string &name);

%ignore simuPOP::SharedVariables::getVarAsString(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsStrDict(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsIntDict(const string &name, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsDoubleNumArray(const string &name, double *&buf, bool nameError=true);

%ignore simuPOP::SharedVariables::getVarAsIntNumArray(const string &name, int *&buf, bool nameError=true);

%feature("docstring")  simuPOP::SharedVariables::dict " 

Usage:
  x.dict()

";

%feature("docstring") simuPOP::Simulator "

simulator manage several replicates of a population, evolve them
using given mating scheme and operators.

Details:
  Simulators combine three important components of simuPOP:
  population, mating scheme and operators together. A simulator is
  created with an instance of population, a replicate number and a
  mating scheme. It makes'rep'r
  eplicates of this population and control the evolution process of
  these populations.
  The most important functions of a simulator is of course evolve().
  It accepts arrays of operators as its parameters, among which,'p
  reOps'and'pos
  tOps'will be applied to the populations at the begining and
  end of evolution, whereas'operators'w
  ill be applied at every generation.
  Simulators separates operators into pre-, during- and post- mating
  operators. During evolution, simulator first apply all pre-mating
  operators and then call the mate() function of the given mating
  scheme, which will call during-mating operators during the birth of
  each offsrping. After the mating is finished, post-mating operators
  are applied in the order they apprear in the operator list.
  Since operators can apply to specific replicate or replicates group,
  and might not be active at all time, the isActive(m_curRep,
  m_numRep, m_gen, end, grp()) function of each
  operator is called before it is applied to the populations.
  Simulators can evolve a given number of generations (the'e
  nd'parameter of evolve), or evolve indefinitely using a certain
  type of operators called terminators. In this case, one or more
  terminators will check the status of evolution and determine if the
  simulation should be stopped. An obvious example of such a
  terminator is a fixation-checker.
  Finally, a simulator can be saved to a file in the format of't
  ext','bin',
   or'xml'. This enables us
  to stop a simulation and ressume it at another time or on another
  machine. It is also a good idea to save a snapshot of a simulation
  every several generations.

";

%feature("docstring")  simuPOP::Simulator::Simulator " 

create a simulator

Usage:
  simulator(&population, &mate, varName=\"simuVars\", rep=1,
    grp=vectori)

Arguments:

  population:  a population created by population() function. This
      population will be copied to the simulator so its content will
      not be changed.
      
  mate:  a mating scheme
      
  rep:  number of replicates. default to 1
      
  grp:  grp number for each replicate. For example, we can seprate all
      replicates into two groups and give them differnt initial values
      before evolution.

Details:

Value:
  none

See Also:
   population, mating

";

%feature("docstring")  simuPOP::Simulator::~Simulator " 

destroy a simulator along with all its populations

Usage:
  x.~Simulator()

Note:
  pop = simulator::population() returns temporary reference to an
  internal population. After a simulator evolves another genertion or
  after the simulator is destroyed, this referenced population should
  *not* be used.

";

%feature("docstring")  simuPOP::Simulator::setMatingScheme " 

set another mating scheme.

Usage:
  x.setMatingScheme(&mate)

";

%feature("docstring")  simuPOP::Simulator::population " 

the'rep'rep
licate of this simulator

Usage:
  x.population(rep=0)

Arguments:

  rep:  number of replicate. This is maybe the only place in si
      muPOPto use 1 based index. E.g. we have population 1-numRep(),
      not 0 ~ ( numRep()-1).

Value:
  reference to a population.

Note:
  The returned reference is temporary in the sense that the refered
  population will be invalid after another round of evolution.
  Therefore, the use of this function should be limited to'i
  mmediate after retrival'.

";

%feature("docstring")  simuPOP::Simulator::setPopulation " 

Usage:
  x.setPopulation(&pop, rep)

";

%feature("docstring")  simuPOP::Simulator::curRep " 

Usage:
  x.curRep()

";

%feature("docstring")  simuPOP::Simulator::numRep " 

Usage:
  x.numRep()

";

%feature("docstring")  simuPOP::Simulator::gen " 

Usage:
  x.gen()

";

%feature("docstring")  simuPOP::Simulator::grp " 

Usage:
  x.grp()

";

%feature("docstring")  simuPOP::Simulator::grps " 

return group indices

Usage:
  x.grps()

";

%feature("docstring")  simuPOP::Simulator::setGroup " 

set groups for replicates

Usage:
  x.setGroup(&grp)

";

%feature("docstring")  simuPOP::Simulator::setGen " 

set generation number

Usage:
  x.setGen(gen)

Arguments:

  gen:  new generation number

Note:
  this will also set shared variable gen

";

%feature("docstring")  simuPOP::Simulator::step " 

evolve one step

Usage:
  x.step(operators, preOps=[], postOps=[], steps=1)

See Also:
  simulator::evolve()

";

%feature("docstring")  simuPOP::Simulator::evolve " 

evolve till'end'g
eneration subject to given operators

Usage:
  x.evolve(ops, preOps=[], postOps=[], end=-1, dryrun=False,
    showHeader=False)

Arguments:

  ops:  operators that will be applied at all generations. Of course
      they might not be active at all generations.
      
  preOps:  operators that will be applied before evolution. evolve()
      function will *not* check if they are active.
      
  postOps:  operators that will be applied after the last generation,
      or after a replicate is terminated.
      
  end:  ending generation. Should be -1 (no ending generation) or a
      number greater than current generation number. When end=-1,
      simulator can only be stopped by terminators.
      
  showHeader:  whether or not headers of operators should be outputed.
      This is used when we stop and resume evolution and would like
      supress intermediate heading output.

Details:
  
  'operators'will be applied in the order
  of:
  all pre-mating opertors
  during-mating operators called by the mating scheme at the birth of
  each offspring
  all post-mating operators
  If any pre or post-mating opertor fails to apply, that replicate
  will be stopped. This is exactly how terminators work.

  if apply to stopped reps, do it.

Value:
  True if evolution finishs successfully.

Note:
  When end=-1, you can not specify negative generation parameters to
  operators. How would an operator know which genertion is the -1
  genertion if no ending genertion is given?

See Also:
  simulator::step()

";

%feature("docstring")  simuPOP::Simulator::apply " 

apply some ops, geneartion of Populationd
oes not change No mating is allowed.

Usage:
  x.apply(ops, dryrun=False, showHeader=False)

Arguments:

  ops:  operators that will be applied at all generations. Of course
      they might not be active at all generations.
      
  showHeader:  whether or not headers of operators should be outputed.
      This is used when we stop and resume evolution and would like
      supress intermediate heading output.

Details:
  
  pre-mating oeprators are applied before post-mating operators. no
  during-mating operators are allowed.
  output header if available.

Value:
  True if evolution finishs successfully.

";

%feature("docstring")  simuPOP::Simulator::stopIfOneRepStop " 

stop if one replicate stops or not

Usage:
  x.stopIfOneRepStop(on=True)

Arguments:

  on:  turn on or off stopIfOneRepStop
      
  if set, the simulator will stop evolution if one replicate stops.
  This is sometimes useful.

";

%feature("docstring")  simuPOP::Simulator::applyOpToStoppedReps " 

apply ops even if rep stops

Usage:
  x.applyOpToStoppedReps(on=True)

Arguments:

  on:  turn on or off applyOpToStoppedReps flag
      
  if set, the simulator will continue to apply operators to all
  stopped repicates until all replicates are marked stopped. This is
  sometimes useful.

";

%feature("docstring")  simuPOP::Simulator::initVars " 

Usage:
  x.initVars()

";

%feature("docstring")  simuPOP::Simulator::removeVars " 

Usage:
  x.removeVars()

";

%feature("docstring")  simuPOP::Simulator::vars " 

get simulator namespace, if rep>0
 is given, return replicate rep namespace

Usage:
  x.vars(rep=-1, subPop=-1)

";

%feature("docstring")  simuPOP::Simulator::saveSimulator " 

save simulator in'text','b
in'or'xml'f
ormat

Usage:
  x.saveSimulator(filename, format=\"auto\")

Arguments:

  filename:  save to filename
      
  format:  format to save. Can be one of't
      ext','b
      in','x
      ml'
      
  The default format is'text'b
  ut the output is not suppored to be read.'b
  in'has smaller size and should be used for large populations.'x
  ml'format is most readable and should be used when you would
  like to convert simuPOPpopulations to other formats.

See Also:
  global function loadSimulator

";

%ignore simuPOP::Simulator::loadSimulator(string filename, string format="data");

%feature("docstring")  simuPOP::Simulator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::SMMMutator "

stepwise mutation model.

Details:
  Stepwise mutation model (SMM) assumes that alleles are represented
  by integer values and that a mutation either increases or decreases
  the allele value by one. For variable number tandem repeats loci
  (VNTR), the allele value is generally taken as the number of tandem
  repeats in the DNA sequence.

See Also:
  Kimura&Ohta 1978

";

%feature("docstring")  simuPOP::SMMMutator::SMMMutator " 

Usage:
  sMMMutator(rate=0, rates=[], atLoci=vectori, maxAllele=0,
    incProb=0.5, output=\">\", outputExpr=\"\", headers=[],
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

Arguments:

  rate:  or rates: mutation rate
      
  incProb:  probability to increase allele state. Default to 1
      
  atLoci:  and other parameters: refer to help(mutator),
      help(baseOperator.__init__)

Details:
  The stepwise mutation model (SMM) is developed for allozymes. It
  provides better description for these kinds of evolutionary
  processes.

";

%feature("docstring")  simuPOP::SMMMutator::~SMMMutator " 

Usage:
  x.~SMMMutator()

";

%feature("docstring")  simuPOP::SMMMutator::mutate " 

how to mutate a single allele. this is usually the only function
that need to be defined by the subclasses.

Usage:
  x.mutate(*allele)

";

%feature("docstring")  simuPOP::SMMMutator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::SMMMutator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::Spread "

initialize genotype by value and then copy to all individuals

";

%feature("docstring")  simuPOP::Spread::Spread " 

Usage:
  spread(ind, subPop=[], stage=PreMating, begin=0, end=1, step=1,
    at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::Spread::~Spread " 

Usage:
  x.~Spread()

";

%feature("docstring")  simuPOP::Spread::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::Spread::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::Spread::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%ignore simuPOP::statAlleleFreq;

%feature("docstring")  simuPOP::statAlleleFreq::statAlleleFreq " 

Usage:
  statAlleleFreq(&atLoci=vectori)

";

%feature("docstring")  simuPOP::statAlleleFreq::addLocus " 

Usage:
  x.addLocus(locus)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleNum " 

Usage:
  x.alleleNum(allele, loc)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleFreq " 

Usage:
  x.alleleFreq(allele, loc)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleNum " 

Usage:
  x.alleleNum(allele, loc, subPop)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleFreq " 

Usage:
  x.alleleFreq(allele, loc, subPop)

";

%feature("docstring")  simuPOP::statAlleleFreq::numOfAlleles " 

Usage:
  x.numOfAlleles()

";

%feature("docstring")  simuPOP::statAlleleFreq::numOfAlleles " 

Usage:
  x.numOfAlleles(subPop)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleles " 

Usage:
  x.alleles(loc)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleles " 

Usage:
  x.alleles(loc, subPop)

";

%feature("docstring")  simuPOP::statAlleleFreq::apply " 

Usage:
  x.apply(&pop)

";

%ignore simuPOP::statFst;

%feature("docstring")  simuPOP::statFst::statFst " 

Usage:
  statFst(&alleleFreq, &heteroFreq, &Fst=vectori)

";

%feature("docstring")  simuPOP::statFst::Fst " 

Usage:
  x.Fst()

";

%feature("docstring")  simuPOP::statFst::Fis " 

Usage:
  x.Fis()

";

%feature("docstring")  simuPOP::statFst::Fit " 

Usage:
  x.Fit()

";

%feature("docstring")  simuPOP::statFst::Fst " 

Usage:
  x.Fst(loc)

";

%feature("docstring")  simuPOP::statFst::Fis " 

Usage:
  x.Fis(loc)

";

%feature("docstring")  simuPOP::statFst::Fit " 

Usage:
  x.Fit(loc)

";

%feature("docstring")  simuPOP::statFst::apply " 

Usage:
  x.apply(&pop)

";

%ignore simuPOP::statGenoFreq;

%feature("docstring")  simuPOP::statGenoFreq::statGenoFreq " 

Usage:
  statGenoFreq(&genoFreq=vectori, phase=False)

";

%feature("docstring")  simuPOP::statGenoFreq::apply " 

Usage:
  x.apply(&pop)

Details:
  go through a single allele for all individual, all diploid
  need to replace previous values

";

%ignore simuPOP::statHaploFreq;

%feature("docstring")  simuPOP::statHaploFreq::statHaploFreq " 

Usage:
  statHaploFreq(&haploFreq=intMatrix)

";

%feature("docstring")  simuPOP::statHaploFreq::addHaplotype " 

Usage:
  x.addHaplotype(&haplo)

";

%feature("docstring")  simuPOP::statHaploFreq::numOfHaplotypes " 

Usage:
  x.numOfHaplotypes(&haplo)

";

%feature("docstring")  simuPOP::statHaploFreq::numOfHaplotypes " 

Usage:
  x.numOfHaplotypes(&haplo, subPop)

";

%feature("docstring")  simuPOP::statHaploFreq::haploNum " 

Usage:
  x.haploNum(&haplo)

";

%feature("docstring")  simuPOP::statHaploFreq::haploFreq " 

Usage:
  x.haploFreq(&haplo)

";

%feature("docstring")  simuPOP::statHaploFreq::haploNum " 

Usage:
  x.haploNum(&haplo, subPop)

";

%feature("docstring")  simuPOP::statHaploFreq::haploFreq " 

Usage:
  x.haploFreq(&haplo, subPop)

";

%feature("docstring")  simuPOP::statHaploFreq::apply " 

Usage:
  x.apply(&pop)

";

%ignore simuPOP::statHeteroFreq;

%feature("docstring")  simuPOP::statHeteroFreq::statHeteroFreq " 

Usage:
  statHeteroFreq(&alleleFreq, &heteroFreq=vectori)

";

%feature("docstring")  simuPOP::statHeteroFreq::addLocus " 

Usage:
  x.addLocus(locus)

";

%feature("docstring")  simuPOP::statHeteroFreq::heteroNum " 

Usage:
  x.heteroNum(allele, loc)

";

%feature("docstring")  simuPOP::statHeteroFreq::heteroFreq " 

Usage:
  x.heteroFreq(allele, loc)

";

%feature("docstring")  simuPOP::statHeteroFreq::heteroNum " 

Usage:
  x.heteroNum(allele, loc, subPop)

";

%feature("docstring")  simuPOP::statHeteroFreq::heteroFreq " 

Usage:
  x.heteroFreq(allele, loc, subPop)

";

%feature("docstring")  simuPOP::statHeteroFreq::apply " 

Usage:
  x.apply(&pop)

";

%ignore simuPOP::statLD;

%feature("docstring")  simuPOP::statLD::statLD " 

Usage:
  statLD(&alleleFreq, &haploFreq, &LD=intMatrix)

";

%feature("docstring")  simuPOP::statLD::apply " 

Usage:
  x.apply(&pop)

";

%ignore simuPOP::statNumOfAffected;

%feature("docstring")  simuPOP::statNumOfAffected::statNumOfAffected " 

Usage:
  statNumOfAffected(numOfAffected=False)

";

%feature("docstring")  simuPOP::statNumOfAffected::activate " 

Usage:
  x.activate(yes=True)

";

%feature("docstring")  simuPOP::statNumOfAffected::numOfAffected " 

Usage:
  x.numOfAffected()

";

%feature("docstring")  simuPOP::statNumOfAffected::numOfUnaffected " 

Usage:
  x.numOfUnaffected()

";

%feature("docstring")  simuPOP::statNumOfAffected::numOfAffected " 

Usage:
  x.numOfAffected(subPop)

";

%feature("docstring")  simuPOP::statNumOfAffected::numOfUnaffected " 

Usage:
  x.numOfUnaffected(subPop)

";

%feature("docstring")  simuPOP::statNumOfAffected::apply " 

Usage:
  x.apply(&pop)

";

%ignore simuPOP::statNumOfAlleles;

%feature("docstring")  simuPOP::statNumOfAlleles::statNumOfAlleles " 

Usage:
  statNumOfAlleles(&calc, &atLoci=vectori)

";

%feature("docstring")  simuPOP::statNumOfAlleles::apply " 

Usage:
  x.apply(&pop)

";

%ignore simuPOP::statNumOfMale;

%feature("docstring")  simuPOP::statNumOfMale::statNumOfMale " 

Usage:
  statNumOfMale(numOfMale=False)

";

%feature("docstring")  simuPOP::statNumOfMale::activate " 

Usage:
  x.activate(yes=True)

";

%feature("docstring")  simuPOP::statNumOfMale::numOfMale " 

Usage:
  x.numOfMale()

";

%feature("docstring")  simuPOP::statNumOfMale::numOfFemale " 

Usage:
  x.numOfFemale()

";

%feature("docstring")  simuPOP::statNumOfMale::numOfMale " 

Usage:
  x.numOfMale(subPop)

";

%feature("docstring")  simuPOP::statNumOfMale::numOfFemale " 

Usage:
  x.numOfFemale(subPop)

";

%feature("docstring")  simuPOP::statNumOfMale::apply " 

Usage:
  x.apply(&pop)

";

%feature("docstring") simuPOP::Stator "

NOTE: the default output for stator is\"\", i.e., no output i.e., stator will write to shared variables and
unless specified by output=\">\"etc, no output will be generated. this class will also list ALL
statistics and its names?

";

%feature("docstring")  simuPOP::Stator::Stator " 

constructor. default to be always active. default to have NO output
(shared variables will be set.) phase: if we treat Aa!=aA, default
is false, i.e, Aa=aA.

Usage:
  stator(output=\"\", outputExpr=\"\", headers=[], stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
    sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::Stator::~Stator " 

destructor

Usage:
  x.~Stator()

";

%feature("docstring")  simuPOP::Stator::clone " 

this function is very important

Usage:
  x.clone()

";

%ignore simuPOP::statPopSize;

%feature("docstring")  simuPOP::statPopSize::statPopSize " 

Usage:
  statPopSize(popSize=False)

";

%feature("docstring")  simuPOP::statPopSize::activate " 

Usage:
  x.activate()

";

%feature("docstring")  simuPOP::statPopSize::apply " 

Usage:
  x.apply(&pop)

";

%ignore simuPOP::StreamElem;

%feature("docstring")  simuPOP::StreamElem::StreamElem " 

Stream element, can be of different
types///////////////////////////////////////////////////////////.

Usage:
  streamElem(&name, readable, realAppend, useString)

";

%feature("docstring")  simuPOP::StreamElem::StreamElem " 

copy constructor, we need to clear rhs.m_stream to avoid closing
file too early this is techniquely advanced (and dangerous)

Usage:
  streamElem(&rhs)

";

%feature("docstring")  simuPOP::StreamElem::~StreamElem " 

destructor

Usage:
  x.~StreamElem()

";

%feature("docstring")  simuPOP::StreamElem::makeReadable " 

the file was write-only, re-open it as read-write

Usage:
  x.makeReadable()

";

%feature("docstring")  simuPOP::StreamElem::makeAppend " 

change the append status

Usage:
  x.makeAppend(append)

";

%ignore simuPOP::StreamElem::stream();

%ignore simuPOP::StreamElem::type();

%ignore simuPOP::StreamElem::info();

%ignore simuPOP::StreamElem::append();

%ignore simuPOP::StreamProvider;

%feature("docstring")  simuPOP::StreamProvider::StreamProvider " 

Stream
provider///////////////////////////////////////////////////////////.

Usage:
  streamProvider(&output, &outputExpr)

";

%ignore simuPOP::StreamProvider::~StreamProvider();

%ignore simuPOP::StreamProvider::setOutput(const string &output, const string &outputExpr);

%ignore simuPOP::StreamProvider::noOutput();

%ignore simuPOP::StreamProvider::getOstream(PyObject *dict=NULL, bool readable=false);

%feature("docstring")  simuPOP::StreamProvider::closeOstream " 

close ostream and delete ostream pointer.. if it is a ofstream.

Usage:
  x.closeOstream()

";

%feature("docstring") swig_const_info "

";

%feature("docstring") swig_globalvar "

";

%feature("docstring") swig_type_info "

";

%feature("docstring") swig_varlinkobject "

";

%feature("docstring") simuPOP::SystemError "

exception, thrown if system error occurs

";

%feature("docstring")  simuPOP::SystemError::SystemError " 

Usage:
  systemError(msg)

";

%feature("docstring") simuPOP::Tagger "

tagger is a during mating operator that tag individual with various
information. Potential usages are 1. record parenting information to
track pedigree. 2. tag a individual/allele and monitor its spread in
the population etc. 3...

Details:
  Bo Peng

";

%feature("docstring")  simuPOP::Tagger::Tagger " 

constructor. default to be always active but no output.

Usage:
  tagger(begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL,
    sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::Tagger::~Tagger " 

destructor

Usage:
  x.~Tagger()

";

%feature("docstring")  simuPOP::Tagger::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring") simuPOP::TerminateIf "

terminate according to a condition which can be, e.g.
any(alleleNum0) == 0 all(alleleNum1)>0
.5 alleleNum0{2} == 0 etc. More specifically, as of now, function
supported any, all, sum operator supported<><=>=
 == When the condition is true, a shared variable var=\"t
erminate\"will be set to current generation.

";

%feature("docstring")  simuPOP::TerminateIf::TerminateIf " 

Usage:
  terminateIf(condition=\"\", var=\"terminate\", output=\"\",
    outputExpr=\"\", headers=[], stage=PostMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::TerminateIf::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::TerminateIf::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::TerminateIf::apply " 

check all alleles in vector allele if they are fixed.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::TerminateIf::outputHeader " 

before evoluion, header function will be called to output header?
This should be controlled in Simulator::evolve.

Usage:
  x.outputHeader(&pop)

";

%feature("docstring")  simuPOP::TerminateIf::~TerminateIf " 

Usage:
  x.~TerminateIf()

";

%feature("docstring") simuPOP::Terminator "

";

%feature("docstring")  simuPOP::Terminator::Terminator " 

constructor. default to be always active.

Usage:
  terminator(output=\">\", outputExpr=\"\", headers=[],
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL, sep=\"\\\\t\")

";

%feature("docstring")  simuPOP::Terminator::~Terminator " 

destructor

Usage:
  x.~Terminator()

";

%feature("docstring")  simuPOP::Terminator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring") simuPOP::TypeError "

exception, thrown if type mismatch

";

%feature("docstring")  simuPOP::TypeError::TypeError " 

Usage:
  typeError(msg)

";

%feature("docstring") simuPOP::ValueError "

exception, thrown if value of range etc

";

%feature("docstring")  simuPOP::ValueError::ValueError " 

Usage:
  valueError(msg)

";

%ignore simuPOP::haploKey(const vectori &seq);

%feature("docstring")  simuPOP::turnOnDebug " 

set debug area, default to turn all code on

Usage:
  turnOnDebug(code)

";

%feature("docstring")  simuPOP::turnOffDebug " 

turn off debug, default to turn all code off

Usage:
  turnOffDebug(code)

";

%feature("docstring")  simuPOP::listDebugCode " 

show all dbg codes (print to cout)

Usage:
  listDebugCode()

";

%ignore simuPOP::simuPOP_kbhit(void);

%ignore simuPOP::simuPOP_getch(void);

%ignore simuPOP::PyObj_As_Bool(PyObject *obj, bool &val);

%ignore simuPOP::PyObj_As_Int(PyObject *obj, int &val);

%ignore simuPOP::PyObj_As_Double(PyObject *obj, double &val);

%ignore simuPOP::PyObj_As_String(PyObject *obj, string &val);

%ignore simuPOP::PyObj_As_StrDict(PyObject *obj, strDict &val);

%ignore simuPOP::PyObj_As_Array(PyObject *obj, vectorf &val);

%ignore simuPOP::PyObj_As_IntDict(PyObject *obj, intDict &val);

%ignore simuPOP::PyObj_Is_IntNumArray(PyObject *obj);

%ignore simuPOP::PyObj_Is_DoubleNumArray(PyObject *obj);

%ignore simuPOP::PyObj_Is_ALLELENumArray(PyObject *obj);

%ignore simuPOP::Int_Vec_As_NumArray(int dim, int *buf, bool copyOver);

%ignore simuPOP::Double_Vec_As_NumArray(int dim, double *buf, bool copyOver);

%ignore simuPOP::Allele_Vec_As_NumArray(int dim, Allele *buf, bool copyOver);

%ignore simuPOP::NumArray_Size(PyObject *obj);

%ignore simuPOP::NumArray_Data(PyObject *obj);

%ignore simuPOP::mainVars();

%ignore simuPOP::moduleVars();

%feature("docstring")  simuPOP::ExpressionValueAsType " 

Usage:
  ExpressionValueAsType(bool, Bool, false)

";

%feature("docstring")  simuPOP::ExpressionValueAsType " 

Usage:
  ExpressionValueAsType(int, Int, 0)

";

%feature("docstring")  simuPOP::ExpressionValueAsType " 

Usage:
  ExpressionValueAsType(double, Double, 0.0)

";

%feature("docstring")  simuPOP::ExpressionValueAsType " 

Usage:
  ExpressionValueAsType(string, String, "")

";

%feature("docstring")  simuPOP::ExpressionValueAsType " 

Usage:
  ExpressionValueAsType(vectorf, Array, )

";

%feature("docstring")  simuPOP::ExpressionValueAsType " 

Usage:
  ExpressionValueAsType(strDict, StrDict, )

";

%feature("docstring")  simuPOP::ExpressionValueAsType " 

Usage:
  ExpressionValueAsType(intDict, IntDict, )

";

%feature("docstring")  simuPOP::ostreamManager " 

return ostream manager

Usage:
  ostreamManager()

";

%feature("docstring")  simuPOP::rng " 

return the global RNG

Usage:
  rng()

";

%feature("docstring")  simuPOP::setRNG " 

set global rng this is temporary since rng()m
ight not exist in the future

Usage:
  setRNG(r, seed)

";

%feature("docstring")  simuPOP::listAllRNG " 

list all available RNG.

Usage:
  listAllRNG()

";

%feature("docstring")  simuPOP::bitSet " 

show turned on bits

Usage:
  bitSet(&bs)

";

%feature("docstring")  simuPOP::gsl_error_handler " 

Global debug and initialization related functions///////////////////

Usage:
  gsl_error_handler(*reason, *, int, gsl_errno)

";

%feature("docstring")  simuPOP::g_cnull " 

null stream

Usage:
  g_cnull(&g_nullStreamBuf)

";

%feature("docstring")  simuPOP::cnull " 

return null stream

Usage:
  cnull()

";

%feature("docstring")  simuPOP::setLogOutput " 

set default output (\"\"m
eans standard output)

Usage:
  setLogOutput(filename)

";

%feature("docstring")  simuPOP::ver " 

Global debug and initialization related
functions///////////////////////////////////////////////////////////
return version infomation of simuPOP.

Usage:
  ver()

";

%ignore simuPOP::initialize();

%feature("docstring")  simuPOP::vecSeparator " 

Usage:
  vecSeparator()

";

%feature("docstring")  simuPOP::Int_Vec_As_NumArray " 

Usage:
  Int_Vec_As_NumArray(dim, *buf)

";

%feature("docstring")  simuPOP::Double_Vec_As_NumArray " 

Usage:
  Double_Vec_As_NumArray(dim, *buf)

";

%feature("docstring")  simuPOP::ALLELE_Vec_As_NumArray " 

Usage:
  ALLELE_Vec_As_NumArray(dim, *buf)

";

%feature("docstring")  simuPOP::setTypeVarImpl " 

Usage:
  setTypeVarImpl(bool, Bool)

";

%feature("docstring")  simuPOP::setTypeVarImpl " 

Usage:
  setTypeVarImpl(int, Int)

";

%feature("docstring")  simuPOP::setTypeVarImpl " 

Usage:
  setTypeVarImpl(double, Double)

";

%feature("docstring")  simuPOP::setTypeVarImpl " 

Usage:
  setTypeVarImpl(&, String)

";

%feature("docstring")  simuPOP::setTypeVarImpl " 

Usage:
  setTypeVarImpl(&, Array)

";

%feature("docstring")  simuPOP::setTypeVarImpl " 

Usage:
  setTypeVarImpl(&, StrDict)

";

%feature("docstring")  simuPOP::setTypeVarImpl " 

Usage:
  setTypeVarImpl(&, IntDict)

";

%feature("docstring")  simuPOP::setTypeVarImpl " 

Usage:
  setTypeVarImpl(*, Variant)

";

%feature("docstring")  simuPOP::setDoubleNumArrayVar " 

Usage:
  setDoubleNumArrayVar(&name, dim, *buf)

";

%feature("docstring")  simuPOP::setDoubleNumArrayVar " 

Usage:
  setDoubleNumArrayVar(&name, dim, *buf, rep)

";

%feature("docstring")  simuPOP::setIntNumArrayVar " 

Usage:
  setIntNumArrayVar(&name, dim, *buf)

";

%feature("docstring")  simuPOP::setIntNumArrayVar " 

Usage:
  setIntNumArrayVar(&name, dim, *buf, rep)

";

%feature("docstring")  simuPOP::getVarAsTypeImpl " 

Usage:
  getVarAsTypeImpl(bool, Bool)

";

%feature("docstring")  simuPOP::getVarAsTypeImpl " 

Usage:
  getVarAsTypeImpl(int, Int)

";

%feature("docstring")  simuPOP::getVarAsTypeImpl " 

Usage:
  getVarAsTypeImpl(double, Double)

";

%feature("docstring")  simuPOP::getVarAsTypeImpl " 

Usage:
  getVarAsTypeImpl(string, String)

";

%feature("docstring")  simuPOP::getVarAsTypeImpl " 

Usage:
  getVarAsTypeImpl(vectorf, Array)

";

%feature("docstring")  simuPOP::getVarAsTypeImpl " 

Usage:
  getVarAsTypeImpl(strDict, StrDict)

";

%feature("docstring")  simuPOP::getVarAsDoubleNumArray " 

Usage:
  getVarAsDoubleNumArray(&name, *&buf)

";

%feature("docstring")  simuPOP::getVarAsDoubleNumArray " 

Usage:
  getVarAsDoubleNumArray(&name, *&buf, rep)

";

%feature("docstring")  simuPOP::getVarAsIntNumArray " 

Usage:
  getVarAsIntNumArray(&name, *&buf)

";

%feature("docstring")  simuPOP::getVarAsIntNumArray " 

Usage:
  getVarAsIntNumArray(&name, *&buf, rep)

";

%feature("docstring")  simuPOP::listVars " 

Usage:
  listVars(rep=-1)

";

%feature("docstring")  simuPOP::evaluate " 

Usage:
  evaluate(expr, stmts, rep, main)

";

%feature("docstring")  simuPOP::resizeSharedDicts " 

internal only create and destroy numRep()s
haredVars

Usage:
  resizeSharedDicts()

";

%feature("docstring")  simuPOP::my_new_handler " 

new handler

Usage:
  my_new_handler()

";

%feature("docstring")  simuPOP::rep " 

Usage:
  rep()

";

%feature("docstring")  simuPOP::repRef " 

Usage:
  repRef()

";

%feature("docstring")  simuPOP::setRep " 

Usage:
  setRep(val)

";

%feature("docstring")  simuPOP::numRep " 

Usage:
  numRep()

";

%feature("docstring")  simuPOP::setNumRep " 

Usage:
  setNumRep(val)

";

%feature("docstring")  simuPOP::gen " 

Usage:
  gen()

";

%feature("docstring")  simuPOP::setGen " 

Usage:
  setGen(val)

";

%feature("docstring")  simuPOP::incGen " 

Usage:
  incGen()

";

%feature("docstring")  simuPOP::endGen " 

Usage:
  endGen()

";

%feature("docstring")  simuPOP::setEndGen " 

Usage:
  setEndGen(val)

";

%feature("docstring")  simuPOP::grp " 

Usage:
  grp()

";

%feature("docstring")  simuPOP::grps " 

Usage:
  grps()

";

%feature("docstring")  simuPOP::setGrp " 

Usage:
  setGrp(val)

";

%feature("docstring")  simuPOP::integrityCheck " 

Usage:
  integrityCheck(&msg)

";

%feature("docstring")  is_carrayobject " 

Usage:
  is_carrayobject(*op)

";

%feature("docstring")  c_getitem " 

Usage:
  c_getitem(*ap, i)

";

%feature("docstring")  c_setitem " 

Usage:
  c_setitem(*ap, i, *v)

";

%feature("docstring")  b_getitem " 

Usage:
  b_getitem(*ap, i)

";

%feature("docstring")  b_setitem " 

Usage:
  b_setitem(*ap, i, *v)

";

%feature("docstring")  BB_getitem " 

Usage:
  BB_getitem(*ap, i)

";

%feature("docstring")  BB_setitem " 

Usage:
  BB_setitem(*ap, i, *v)

";

%feature("docstring")  h_getitem " 

Usage:
  h_getitem(*ap, i)

";

%feature("docstring")  h_setitem " 

Usage:
  h_setitem(*ap, i, *v)

";

%feature("docstring")  HH_getitem " 

Usage:
  HH_getitem(*ap, i)

";

%feature("docstring")  HH_setitem " 

Usage:
  HH_setitem(*ap, i, *v)

";

%feature("docstring")  i_getitem " 

Usage:
  i_getitem(*ap, i)

";

%feature("docstring")  i_setitem " 

Usage:
  i_setitem(*ap, i, *v)

";

%feature("docstring")  II_getitem " 

Usage:
  II_getitem(*ap, i)

";

%feature("docstring")  II_setitem " 

Usage:
  II_setitem(*ap, i, *v)

";

%feature("docstring")  l_getitem " 

Usage:
  l_getitem(*ap, i)

";

%feature("docstring")  l_setitem " 

Usage:
  l_setitem(*ap, i, *v)

";

%feature("docstring")  LL_getitem " 

Usage:
  LL_getitem(*ap, i)

";

%feature("docstring")  LL_setitem " 

Usage:
  LL_setitem(*ap, i, *v)

";

%feature("docstring")  f_getitem " 

Usage:
  f_getitem(*ap, i)

";

%feature("docstring")  f_setitem " 

Usage:
  f_setitem(*ap, i, *v)

";

%feature("docstring")  d_getitem " 

Usage:
  d_getitem(*ap, i)

";

%feature("docstring")  d_setitem " 

Usage:
  d_setitem(*ap, i, *v)

";

%feature("docstring")  newarrayobject " 

Usage:
  newarrayobject(size, *descr)

";

%feature("docstring")  getarrayitem " 

Usage:
  getarrayitem(*op, i)

";

%feature("docstring")  ins1 " 

Usage:
  ins1(*self, where, *v)

";

%feature("docstring")  array_dealloc " 

Usage:
  array_dealloc(*op)

";

%feature("docstring")  array_richcompare " 

Usage:
  array_richcompare(*v, *w, op)

";

%feature("docstring")  array_length " 

Usage:
  array_length(*a)

";

%feature("docstring")  array_item " 

Usage:
  array_item(*a, i)

";

%feature("docstring")  array_slice " 

Usage:
  array_slice(*a, ilow, ihigh)

";

%feature("docstring")  array_concat " 

Usage:
  array_concat(*a, *bb)

";

%feature("docstring")  array_repeat " 

Usage:
  array_repeat(*a, n)

";

%feature("docstring")  array_ass_slice " 

Usage:
  array_ass_slice(*a, ilow, ihigh, *v)

";

%feature("docstring")  array_ass_item " 

Usage:
  array_ass_item(*a, i, *v)

";

%feature("docstring")  setarrayitem " 

Usage:
  setarrayitem(*a, i, *v)

";

%feature("docstring")  ins " 

Usage:
  ins(*self, where, *v)

";

%feature("docstring")  array_count " 

Usage:
  array_count(*self, *args)

";

%feature("docstring")  array_index " 

Usage:
  array_index(*self, *args)

";

%feature("docstring")  array_remove " 

Usage:
  array_remove(*self, *args)

";

%feature("docstring")  array_pop " 

Usage:
  array_pop(*self, *args)

";

%feature("docstring")  array_extend " 

Usage:
  array_extend(*self, *args)

";

%feature("docstring")  array_insert " 

Usage:
  array_insert(*self, *args)

";

%feature("docstring")  array_buffer_info " 

Usage:
  array_buffer_info(*self, *args)

";

%feature("docstring")  array_append " 

Usage:
  array_append(*self, *args)

";

%feature("docstring")  array_byteswap " 

Usage:
  array_byteswap(*self, *args)

";

%feature("docstring")  array_reverse " 

Usage:
  array_reverse(*self, *args)

";

%feature("docstring")  array_fromfile " 

Usage:
  array_fromfile(*self, *args)

";

%feature("docstring")  array_tofile " 

Usage:
  array_tofile(*self, *args)

";

%feature("docstring")  array_fromlist " 

Usage:
  array_fromlist(*self, *args)

";

%feature("docstring")  array_tolist " 

Usage:
  array_tolist(*self, *args)

";

%feature("docstring")  array_fromstring " 

Usage:
  array_fromstring(*self, *args)

";

%feature("docstring")  array_tostring " 

Usage:
  array_tostring(*self, *args)

";

%feature("docstring")  array_getattr " 

Usage:
  array_getattr(*a, *name)

";

%feature("docstring")  array_print " 

Usage:
  array_print(*a, *fp, flags)

";

%feature("docstring")  array_repr " 

Usage:
  array_repr(*a)

";

%feature("docstring")  array_buffer_getreadbuf " 

Usage:
  array_buffer_getreadbuf(*self, index, **ptr)

";

%feature("docstring")  array_buffer_getwritebuf " 

Usage:
  array_buffer_getwritebuf(*self, index, **ptr)

";

%feature("docstring")  array_buffer_getsegcount " 

Usage:
  array_buffer_getsegcount(*self, *lenp)

";

%feature("docstring")  a_array " 

Usage:
  a_array(*self, *args)

";

%feature("docstring")  initcarray " 

Usage:
  initcarray(void)

";

%feature("docstring")  carray_length " 

Usage:
  carray_length(*a)

";

%feature("docstring")  carray_type " 

Usage:
  carray_type(*a)

";

%feature("docstring")  carray_data " 

Usage:
  carray_data(*a)

";

%feature("docstring")  newcarrayobjectfrommem " 

Usage:
  newcarrayobjectfrommem(type, size, *ptr, copyOver)

";

%feature("docstring")  ownmem " 

Usage:
  ownmem(*obj)

";

%feature("docstring")  match " 

Usage:
  match(indRanges)

";

%feature("docstring")  variables " 

Usage:
  variables(linux)

";

%feature("docstring")  binary " 

Usage:
  binary(port)

";

%feature("docstring")  SWIG_TypeNameComp " 

Usage:
  SWIG_TypeNameComp(*f1, *l1, *f2, *l2)

";

%feature("docstring")  SWIG_TypeEquiv " 

Usage:
  SWIG_TypeEquiv(*nb, *tb)

";

%feature("docstring")  SWIG_TypeRegister " 

Usage:
  SWIG_TypeRegister(*ti)

";

%feature("docstring")  SWIG_TypeCheck " 

Usage:
  SWIG_TypeCheck(*c, *ty)

";

%feature("docstring")  SWIG_TypeCast " 

Usage:
  SWIG_TypeCast(*ty, *ptr)

";

%feature("docstring")  SWIG_TypeDynamicCast " 

Usage:
  SWIG_TypeDynamicCast(*ty, **ptr)

";

%feature("docstring")  SWIG_TypeName " 

Usage:
  SWIG_TypeName(*ty)

";

%feature("docstring")  SWIG_TypePrettyName " 

Usage:
  SWIG_TypePrettyName(*type)

";

%feature("docstring")  SWIG_TypeQuery " 

Usage:
  SWIG_TypeQuery(*name)

";

%feature("docstring")  SWIG_TypeClientData " 

Usage:
  SWIG_TypeClientData(*ti, *clientdata)

";

%feature("docstring")  SWIG_PackData " 

Usage:
  SWIG_PackData(*c, *ptr, sz)

";

%feature("docstring")  SWIG_UnpackData " 

Usage:
  SWIG_UnpackData(*c, *ptr, sz)

";

%feature("docstring")  SWIG_PropagateClientData " 

Usage:
  SWIG_PropagateClientData(*type)

";

%feature("docstring")  swig_varlink_repr " 

Usage:
  swig_varlink_repr(*v)

";

%feature("docstring")  swig_varlink_print " 

Usage:
  swig_varlink_print(*v, *fp, flags)

";

%feature("docstring")  swig_varlin " 

Usage:
  swig_varlink_getattr(*v, *n)

";

%feature("docstring")  swig_varlink_setattr " 

Usage:
  swig_varlink_setattr(*v, *n, *p)

";

%feature("docstring")  SWIG_Python_newvarlink " 

Usage:
  SWIG_Python_newvarlink(void)

";

%feature("docstring")  SWIG_Python_addvarlink " 

Usage:
  SWIG_Python_addvarlink(*p, *name, *get_attr, *p)

";

%feature("docstring")  SWIG_Python_TypeError " 

Usage:
  SWIG_Python_TypeError(*type, *obj)

";

%feature("docstring")  SWIG_Python_NullRef " 

Usage:
  SWIG_Python_NullRef(*type)

";

%feature("docstring")  SWIG_Python_AddErrMesg " 

Usage:
  SWIG_Python_AddErrMesg(*mesg, infront)

";

%feature("docstring")  SWIG_Python_ArgFail " 

Usage:
  SWIG_Python_ArgFail(argnum)

";

%feature("docstring")  SWIG_Python_ConvertPtr " 

Usage:
  SWIG_Python_ConvertPtr(*obj, **ptr, *ty, flags)

";

%feature("docstring")  SWIG_Python_MustGetPtr " 

Usage:
  SWIG_Python_MustGetPtr(*obj, *ty, argnum, flags)

";

%feature("docstring")  SWIG_Python_ConvertPacked " 

Usage:
  SWIG_Python_ConvertPacked(*obj, *ptr, sz, *ty, flags)

";

%feature("docstring")  SWIG_Python_PointerStr " 

Usage:
  SWIG_Python_PointerStr(*buff, *ptr, *name, bsz)

";

%feature("docstring")  SWIG_Python_NewPointerObj " 

Usage:
  SWIG_Python_NewPointerObj(*ptr, *type, own)

";

%feature("docstring")  SWIG_Python_NewPackedObj " 

Usage:
  SWIG_Python_NewPackedObj(*ptr, sz, *type)

";

%feature("docstring")  SWIG_Python_InstallConstants " 

Usage:
  SWIG_Python_InstallConstants(*d, constants[])

";

%feature("docstring")  SWIG_Python_FixMethods " 

Usage:
  SWIG_Python_FixMethods(*methods, *const_table, **types,
    **types_initial)

";

%feature("docstring")  PyModule_AddObject " 

Usage:
  PyModule_AddObject(*m, *name, *o)

";

%feature("docstring")  SWIG_Python_LookupTypePointer " 

Usage:
  SWIG_Python_LookupTypePointer(***type_list_handle)

";

%feature("docstring")  SWIG_exception_ " 

Usage:
  SWIG_exception_(code, *msg)

";

%feature("docstring")  SWIG_FromCharPtr " 

Usage:
  SWIG_FromCharPtr(*cptr)

";

%feature("docstring")  loadPopulation " 

Usage:
  loadPopulation(file, format=\"text\")

";

%feature("docstring")  loadSimulator " 

Usage:
  loadSimulator(file, &mate, format=\"text\")

";

%feature("docstring")  _wrap_swig_set " 

Usage:
  _wrap_swig_set(*_val)

";

%feature("docstring")  _wrap_swig_get " 

Usage:
  _wrap_swig_get(void)

";

%feature("docstring")  _wrap_std_set " 

Usage:
  _wrap_std_set(*_val)

";

%feature("docstring")  _wrap_std_get " 

Usage:
  _wrap_std_get(void)

";

%feature("docstring")  _wrap_simuPOP_set " 

Usage:
  _wrap_simuPOP_set(*_val)

";

%feature("docstring")  _wrap_simuPOP_get " 

Usage:
  _wrap_simuPOP_get(void)

";

%feature("docstring")  _wrap_loadPopulation " 

Usage:
  _wrap_loadPopulation(*self, *args)

";

%feature("docstring")  _wrap_loadSimulator " 

Usage:
  _wrap_loadSimulator(*self, *args)

";

%feature("docstring")  SWIGEXPORT " 

Usage:
  SWIGEXPORT(void)

";

%feature("docstring")  simuPOP::newcarrayobjectfrommem " 

Usage:
  newcarrayobjectfrommem(type, size, *ptr, copyOver)

";

%feature("docstring")  simuPOP::is_carrayobject " 

Usage:
  is_carrayobject(*)

";

%feature("docstring")  simuPOP::carray_length " 

Usage:
  carray_length(*a)

";

%feature("docstring")  simuPOP::carray_type " 

Usage:
  carray_type(*a)

";

%feature("docstring")  simuPOP::carray_data " 

Usage:
  carray_data(*a)

";

%feature("docstring")  simuPOP::initcarray " 

Usage:
  initcarray(void)

";

%feature("docstring")  simuPOP::ownmem " 

Usage:
  ownmem(*a)

";

%feature("docstring")  simuPOP::getModuleDict " 

Usage:
  getModuleDict()

";

%feature("docstring")  simuPOP::getMainDict " 

Usage:
  getMainDict()

";

