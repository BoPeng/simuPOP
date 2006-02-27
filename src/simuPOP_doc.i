

%feature("docstring") simuPOP::affectedSibpairSample "

thrink population accroding to some outside value

";

%feature("docstring")  simuPOP::affectedSibpairSample::affectedSibpairSample " 

draw cases and controls

Usage:
  affectedSibpairSample(size=[], chooseUnaffected=False,
    countOnly=False, &name=\"sample\", &nameExpr=\"\", times=1,
    &saveAs=\"\", &saveAsExpr=\"\", &format=\"auto\",
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  size:  number of affected sibpairs to be sampled. Can be a number or
      an array. If a number is given, it is the total number of
      sibpairs, ignoring population structure. Otherwise, given number
      of sibpairs are sampled from subpopulations. If size is
      unspecified, this operator will return all affected sibpairs.
      
  countOnly:  set variables about number of affected sibpairs, do not
      actually draw the sample
      
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

%feature("docstring")  simuPOP::affectedSibpairSample::~affectedSibpairSample " 

destructor

Usage:
  x.~affectedSibpairSample()

";

%feature("docstring")  simuPOP::affectedSibpairSample::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::affectedSibpairSample::preparesample " 

Usage:
  x.preparesample(&pop)

";

%feature("docstring")  simuPOP::affectedSibpairSample::drawsample " 

Usage:
  x.drawsample(&pop)

Details:
  sibs when m_size.size() <= 1

";

%feature("docstring")  simuPOP::affectedSibpairSample::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") arraydescr "

";

%feature("docstring") arrayobject "

";

%feature("docstring") arrayobject::iterator "

";

%ignore simuPOP::BernulliTrials;

%ignore simuPOP::BernulliTrials::BernulliTrials(RNG &rng);

%ignore simuPOP::BernulliTrials::BernulliTrials(RNG &rng, const vectorf &prob, ULONG trials);

%ignore simuPOP::BernulliTrials::~BernulliTrials();

%ignore simuPOP::BernulliTrials::size() ;

%feature("docstring")  simuPOP::BernulliTrials::prob " 

Usage:
  x.prob()

";

%ignore simuPOP::BernulliTrials::setParameter(const vectorf &prob, ULONG trials);

%ignore simuPOP::BernulliTrials::doTrial();

%ignore simuPOP::BernulliTrials::curTrial();

%feature("docstring")  simuPOP::BernulliTrials::trial " 

get a trial corresponding to m_prob.

Usage:
  x.trial()

";

%feature("docstring")  simuPOP::BernulliTrials::succ " 

return succeed trials for p[index] fail when m_cur is not 0. (i.e.,
has retrieve the table through trial()C
PPONLY

Usage:
  x.succ(index)

";

%ignore simuPOP::BernulliTrials::probabilities();

%ignore simuPOP::BernulliTrials::numTrials();

%feature("docstring") simuPOP::binomialSelection "

Details:
  binomial random selection
  No sex. Choose one individual from last generation.
  1. numOffspring protocol is honored 2. population size changes are
  allowed 3. selection is possible.
  So this works just like a sexless random mating. If ploidy is one,
  this is chromosomal mating.

";

%feature("docstring")  simuPOP::binomialSelection::binomialSelection " 

constructor

Usage:
  binomialSelection(numOffspring=1., *numOffspringFunc=NULL,
    maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
    newSubPopSizeExpr=\"\", *newSubPopSizeFunc=NULL)

";

%feature("docstring")  simuPOP::binomialSelection::~binomialSelection " 

destructor

Usage:
  x.~binomialSelection()

";

%feature("docstring")  simuPOP::binomialSelection::clone " 

 clone() const. The same as mating::clone() const.

Usage:
  x.clone()

See Also:
   mating::clone() const

";

%feature("docstring")  simuPOP::binomialSelection::__repr__ " 

return name of the mating type

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::binomialSelection::submitScratch " 

Usage:
  x.submitScratch(&pop, &scratch)

";

%feature("docstring")  simuPOP::binomialSelection::mate " 

do the mating.

Usage:
  x.mate(&pop, &scratch, &ops, submit)

Arguments:

  pop:  population
      
  scratch:  scratch population, will be used in this mating scheme.
      
  ops:  during mating operators

Details:

  determine if mate()will generate offspring genotype
  use deep copy!!!!!!!

Value:
  return false when mating fails.

";

%feature("docstring") simuPOP::caseControlSample "

thrink population accroding to some outside value

";

%feature("docstring")  simuPOP::caseControlSample::caseControlSample " 

draw cases and controls

Usage:
  caseControlSample(&cases=vectori, &controls=vectori,
    spSample=False, &name=\"sample\", &nameExpr=\"\", times=1,
    &saveAs=\"\", &saveAsExpr=\"\", &format=\"auto\",
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  cases:  number of cases, or an array of number of cases from each
      subpopulation.
      
  controls:  number of controls, or an array of number of controls
      from each subpopulation.
      
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

%feature("docstring")  simuPOP::caseControlSample::~caseControlSample " 

destructor

Usage:
  x.~caseControlSample()

";

%feature("docstring")  simuPOP::caseControlSample::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::caseControlSample::preparesample " 

Usage:
  x.preparesample(&pop)

";

%feature("docstring")  simuPOP::caseControlSample::drawsample " 

Usage:
  x.drawsample(&pop)

";

%feature("docstring")  simuPOP::caseControlSample::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::continueIf "

terminate according to a condition which can be, e.g.
any(alleleNum0) == 0 all(alleleNum1) > 0.5 alleleNum0{2} == 0 etc.
When the condition is true, a shared variable var=\"terminate\" will
be set to current generation.

";

%feature("docstring")  simuPOP::continueIf::continueIf " 

Usage:
  continueIf(condition=\"\", message=\"\", var=\"terminate\",
    output=\"\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::continueIf::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::continueIf::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::continueIf::apply " 

check all alleles in vector allele if they are fixed.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::continueIf::~continueIf " 

Usage:
  x.~continueIf()

";

%feature("docstring") simuPOP::controlledMating "

Details:
  controlled mating

";

%feature("docstring")  simuPOP::controlledMating::controlledMating " 

Controlled mating, control allele frequency at a locus.

Usage:
  controlledMating(&matingScheme, loci, alleles, *freqFunc,
    range=0.01)

Arguments:

  mating:  a mating scheme.
      
  loci:  loci at which allele frequency is controlled. Note that
      controlling the allele frequencies at several loci may take a
      long time.
      
  alleles:  alleles to control at each loci. Should have the same
      length as loci
      
  freqFunc:  frequency boundaries. If the length of the return value
      equals the size of loci, the range for loci will be [value0,
      value0+range], [value1, value1+range] etc. If the length of the
      return value is 2 times size of loci, it will be interpreted as
      [low1, high1, low2, high2 ...]

";

%ignore simuPOP::controlledMating::controlledMating(const controlledMating &rhs);

%feature("docstring")  simuPOP::controlledMating::~controlledMating " 

destructor

Usage:
  x.~controlledMating()

";

%feature("docstring")  simuPOP::controlledMating::clone " 

 clone() const. Generate a copy of itself and return pointer this
is to make sure the object is persistent and will not be freed by
python.

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::controlledMating::isCompatible " 

check if the mating type is compatible with population structure

Usage:
  x.isCompatible(&pop)

Details:
  possible things to check:need certain types of individual (age,
  sex etc)
  need resizeable population...

";

%feature("docstring")  simuPOP::controlledMating::__repr__ " 

return name of the mating type

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::controlledMating::countAlleles " 

Usage:
  x.countAlleles(&pop, &loci, &alleles)

";

%feature("docstring")  simuPOP::controlledMating::mate " 

mate: This is not supposed to be called for base mating class.

Usage:
  x.mate(&pop, &scratch, &ops, submit)

Arguments:

  pop:  population
      
  scratch:  scratch population
      
  ops:  during mating operators

Value:
  return false when mating fail.

";

%feature("docstring") simuPOP::dumper "

dump the content of a population.

";

%feature("docstring")  simuPOP::dumper::dumper " 

dump population

Usage:
  dumper(alleleOnly=False, infoOnly=False, ancestralPops=False,
    dispWidth=1, max=100, &chrom=vectori, &loci=vectori, &subPop=[],
    &indRange=[], output=\">\", outputExpr=\"\", stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  alleleOnly:  only display allele
      
  infoOnly:  only display info
      
  dispWidth:  width of allele display, default to 1
      
  ancestralPops:  whether or not display ancestral populations,
      default to False
      
  chrom:  chromsoome(s) to display
      
  loci:  loci to display
      
  subPop:  only display subPop(s)
      
  indRange:  range(s) of individuals to display
      
  max:  max number of individuals to display, default to 100. This is
      to avoid careless dump of huge populations.
      
  output:  output file, default to standard output.
      
  outputExpr:  and other parameters: refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::dumper::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::dumper::alleleOnly " 

only show alleles (not structure, gene information?

Usage:
  x.alleleOnly()

";

%feature("docstring")  simuPOP::dumper::setAlleleOnly " 

Usage:
  x.setAlleleOnly(alleleOnly)

";

%feature("docstring")  simuPOP::dumper::infoOnly " 

only show info

Usage:
  x.infoOnly()

";

%feature("docstring")  simuPOP::dumper::setInfoOnly " 

Usage:
  x.setInfoOnly(infoOnly)

";

%feature("docstring")  simuPOP::dumper::apply " 

Usage:
  x.apply(&pop)

Details:
  dump population structure
  dump all genotypic info

";

%feature("docstring")  simuPOP::dumper::~dumper " 

Usage:
  x.~dumper()

";

%feature("docstring")  simuPOP::dumper::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::Exception "

exception handler. Exceptions will be passed to Python.

";

%feature("docstring")  simuPOP::Exception::Exception " 

constructor msg error message

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

%ignore simuPOP::GenoStructure::GenoStructure(UINT ploidy, const vectoru &loci, bool sexChrom, const vectorf &lociPos, const vectorstr &alleleNames, const vectorstr &lociNames, UINT maxAllele);

%ignore simuPOP::GenoStructure::GenoStructure(const GenoStructure &rhs);

%feature("docstring")  simuPOP::GenoStructure::~GenoStructure " 

destructor, do nothing.

Usage:
  x.~GenoStructure()

";

%feature("docstring") simuPOP::GenoStruTrait "

genoStruTrait

Details:
  A trait class that maintain a static array of geno structure, and
  provide interfaces around a GenoStructureI
  ndex.

";

%feature("docstring")  simuPOP::GenoStruTrait::GenoStruTrait " 

constructor, but m_genoStruIdx will be set later.

Usage:
  genoStruTrait()

";

%ignore simuPOP::GenoStruTrait::setGenoStructure(UINT ploidy, const vectoru &loci, bool sexChrom, const vectorf &lociPos, const vectorstr &alleleNames, const vectorstr &lociNames, UINT maxAllele);

%feature("docstring")  simuPOP::GenoStruTrait::setGenoStructure " 

set an existing geno structure, simply use it This is NOT efficient!
(but has to be used when, for example, loading a structure from file

Usage:
  x.setGenoStructure(&rhs)

";

%ignore simuPOP::GenoStruTrait::setGenoStruIdx(size_t idx);

%ignore simuPOP::GenoStruTrait::genoStru() ;

%ignore simuPOP::GenoStruTrait::genoStruIdx() ;

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
  x.numLoci()

";

%feature("docstring")  simuPOP::GenoStruTrait::sexChrom " 

whether or not the last chromosome is sex chromosome

Usage:
  x.sexChrom()

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

%feature("docstring")  simuPOP::GenoStruTrait::locusPos " 

locus distance.

Usage:
  x.locusPos()

";

%feature("docstring")  simuPOP::GenoStruTrait::arrLociPos " 

expose loci distance

Usage:
  x.arrLociPos()

";

%feature("docstring")  simuPOP::GenoStruTrait::arrLociPos " 

expose loci distance of a chromosome

Usage:
  x.arrLociPos(chrom)

";

%feature("docstring")  simuPOP::GenoStruTrait::numChrom " 

number of chromosome

Usage:
  x.numChrom()

";

%feature("docstring")  simuPOP::GenoStruTrait::chromBegin " 

chromosome index of chromosome chrom

Usage:
  x.chromBegin()

";

%feature("docstring")  simuPOP::GenoStruTrait::chromEnd " 

chromosome index of chromosome chrom

Usage:
  x.chromEnd()

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
  x.chromLocusPair()

";

%feature("docstring")  simuPOP::GenoStruTrait::alleleName " 

return allele name

Usage:
  x.alleleName()

";

%feature("docstring")  simuPOP::GenoStruTrait::alleleNames " 

allele names

Usage:
  x.alleleNames()

";

%feature("docstring")  simuPOP::GenoStruTrait::locusName " 

return locus name

Usage:
  x.locusName()

";

%feature("docstring")  simuPOP::GenoStruTrait::maxAllele " 

Usage:
  x.maxAllele()

";

%feature("docstring")  simuPOP::GenoStruTrait::setMaxAllele " 

Usage:
  x.setMaxAllele(maxAllele)

";

%feature("docstring")  simuPOP::GenoStruTrait::swap " 

Usage:
  x.swap(&rhs)

";

%feature("docstring") simuPOP::gsmMutator "

stepwise mutation model.

Details:
  Generalized Stepwise mutation model (GSM) assumes that alleles are
  represented by integer values and that a mutation either increases
  or decreases the allele value by a random value.

See Also:
  Kimura & Ohta 1978

";

%feature("docstring")  simuPOP::gsmMutator::gsmMutator " 

Usage:
  gsmMutator(rate=[], atLoci=vectori, maxAllele=0, incProb=0.5, p=0,
    *func=NULL, output=\">\", outputExpr=\"\", stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  rate::  mutation rate
      
  incProb:  probability to increase allele state. Default to 0.5
      
  atLoci:  and other parameters: refer to help(mutator),
      help(baseOperator.__init__)
      
  func:  return number of steps. no parameter

Details:
  The generalized stepwise mutation model (GMM) is developed for
  allozymes. It provides better description for these kinds of
  evolutionary processes.

";

%feature("docstring")  simuPOP::gsmMutator::~gsmMutator " 

Usage:
  x.~gsmMutator()

";

%feature("docstring")  simuPOP::gsmMutator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::gsmMutator::mutate " 

how to mutate a single allele. this is usually the only function
that need to be defined by the subclasses.

Usage:
  x.mutate(allele)

";

%feature("docstring")  simuPOP::gsmMutator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::ifElse "

";

%feature("docstring")  simuPOP::ifElse::ifElse " 

conditional operator

Usage:
  ifElse(&cond, *ifOp=NULL, *elseOp=NULL, output=\">\",
    outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL)

Arguments:

  cond:  expression, will be treated as bool variable.
      
  ifOp:  if operator, be called when expr is true
      
  elseOp:  else operator, be called when expr is false

";

%feature("docstring")  simuPOP::ifElse::~ifElse " 

destructor

Usage:
  x.~ifElse()

";

%ignore simuPOP::ifElse::ifElse(const ifElse &rhs);

%feature("docstring")  simuPOP::ifElse::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::ifElse::applyWithScratch " 

simply output some info providing interface to apply operator before
during or after mating.

Usage:
  x.applyWithScratch(&pop, &scratch, stage)

";

%feature("docstring")  simuPOP::ifElse::applyDuringMating " 

give pop, offspring, pop and mom.

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring")  simuPOP::ifElse::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::ifElse::__repr__ " 

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

%feature("docstring") simuPOP::individual "

Basic individual class.

Details:
  class individual withgenotypic information
  shared genotypic structure info (through a GenoStructurep
  ointer)
  flags about sex, affected status
  an internal info field and a tag field of any type (template)
  other individuals will be derived from this class, adding age info
  etc.
   Notethatindividual DOES NOT manage memory. It will use a
  pointer passed from class population. This causes A LOT of trouble
  and I have not evaluated how much benefic I get.
  individual is a template class taking a tag parameter. tag can be a
  integer, a pair or any object to track individual information. The
  prolem with our template based design + SWIG makes instantiation of
  mutiple tag types difficult. :-(
  operator = uses shallow copy. This is required by sort algorithm
  since otherwise individuals are non-copiable. However, in population
  memory management, it is showtimes required that genotypic
  information within one subPop should go together. This is done by a
  shollow_copied flag for each individual and for all individuals.
  population might have to re-arrange individuals to solve this
  problem.
  output of individual can be adjusted by setOutputDelimeter.
  Usage info: (for population classes developers)for individuals are
  created, you are responsible to set its genotypic pointer and
  genotypic information. This is done by
    setInfo()and  info()can be
  used for any temporarypurpose.
   Tagare usually std::pair(int,int) used by class tagger to track
  informations like pedigree structure.

";

%feature("docstring")  simuPOP::individual::individual " 

default constructor, Tag field need a default constructor

Usage:
  individual()

";

%ignore simuPOP::individual::individual(const individual &ind);

%feature("docstring")  simuPOP::individual::~individual " 

destructor. Do nothing.

Usage:
  x.~individual()

";

%ignore simuPOP::individual::setGenoPtr(GenoIterator pos);

%feature("docstring")  simuPOP::individual::copyFrom " 

Deep copy! Important!

Usage:
  x.copyFrom(&rhs)

";

%ignore simuPOP::individual::genoPtr() ;

%feature("docstring")  simuPOP::individual::arrGenotype " 

return genotype as python Numeric.array object This is the whole
genotype (all)

Usage:
  x.arrGenotype()

";

%feature("docstring")  simuPOP::individual::arrGenotype " 

return genotype as python Numeric.array object This is the p'th copy
of chromosomes

Usage:
  x.arrGenotype(p)

";

%feature("docstring")  simuPOP::individual::arrGenotype " 

return genotype as python Numeric.array object This is the ch
chromosome of the pth copy of chromosome

Usage:
  x.arrGenotype(p, ch)

";

%feature("docstring")  simuPOP::individual::allele " 

get allele from an index

Usage:
  x.allele()

Arguments:

  index:  index from the beginning of genotypic info

";

%feature("docstring")  simuPOP::individual::allele " 

get allele from an index, on the pth set of chromosome

Usage:
  x.allele(index, )

Arguments:

  index:  index from the begining of the p'th set of chromosomes.
      
  p:  on p'th set of chromosomes, default to 0

";

%feature("docstring")  simuPOP::individual::allele " 

Usage:
  x.allele(index, p, )

";

%feature("docstring")  simuPOP::individual::alleleChar " 

Usage:
  x.alleleChar()

";

%feature("docstring")  simuPOP::individual::alleleChar " 

get allele from an index, on the pth set of chromosome

Usage:
  x.alleleChar(index, )

Arguments:

  index:  index from the begining of the p'th set of chromosomes.
      
  p:  on p'th set of chromosomes, p=0 by default

";

%feature("docstring")  simuPOP::individual::alleleChar " 

get allele from an index, on the pth set of chromosome

Usage:
  x.alleleChar(index, p, )

Arguments:

  index:  index from the begining of the p'th set of chromosomes.
      
  p:  on p'th set of chromosomes, p=0 by default

";

%feature("docstring")  simuPOP::individual::setAllele " 

set allele from an index.

Usage:
  x.setAllele(allele, index)

Arguments:

  index:  index from the begining of genotype

";

%feature("docstring")  simuPOP::individual::setAllele " 

set allele from an index.

Usage:
  x.setAllele(allele, index, p)

Arguments:

  allele:  allele to set
      
  index:  index from the begining of genotype
      
  p:  on p'th set of chromosome, p=0 by default

";

%feature("docstring")  simuPOP::individual::setAllele " 

Usage:
  x.setAllele(allele, index, p, ch)

";

%feature("docstring")  simuPOP::individual::tag " 

return tag

Usage:
  x.tag()

";

%feature("docstring")  simuPOP::individual::setTag " 

set tag

Usage:
  x.setTag(tag)

";

%feature("docstring")  simuPOP::individual::setTag " 

set a single number?

Usage:
  x.setTag(tag)

";

%feature("docstring")  simuPOP::individual::sex " 

sex?

Usage:
  x.sex()

";

%feature("docstring")  simuPOP::individual::sexChar " 

return M or F for sex, for display purpose

Usage:
  x.sexChar()

";

%feature("docstring")  simuPOP::individual::setSex " 

set sex

Usage:
  x.setSex(sex)

";

%feature("docstring")  simuPOP::individual::affected " 

affected?

Usage:
  x.affected()

";

%feature("docstring")  simuPOP::individual::unaffected " 

unaffected?

Usage:
  x.unaffected()

";

%feature("docstring")  simuPOP::individual::affectedChar " 

return A or U for affected/Unaffected, for display purpose

Usage:
  x.affectedChar()

";

%feature("docstring")  simuPOP::individual::setAffected " 

set affected status

Usage:
  x.setAffected(affected)

";

%feature("docstring")  simuPOP::individual::info " 

get info

Usage:
  x.info()

";

%feature("docstring")  simuPOP::individual::setInfo " 

set info

Usage:
  x.setInfo(info)

";

%ignore simuPOP::individual::genoBegin() ;

%ignore simuPOP::individual::genoEnd() ;

%ignore simuPOP::individual::genoBegin(UINT p) ;

%ignore simuPOP::individual::genoEnd(UINT p) ;

%ignore simuPOP::individual::genoBegin(UINT p, UINT chrom) ;

%ignore simuPOP::individual::genoEnd(UINT p, UINT chrom) ;

%feature("docstring")  simuPOP::individual::__cmp__ " 

Usage:
  x.__cmp__()

";

%feature("docstring")  simuPOP::individual::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::individual::swap " 

swap individuals

Usage:
  x.swap(&ind, swapContent=True)

Arguments:

  ind:  individual to be swapped in
      
  swapContent:  swapContent or only the pointers.
      
  The guideline is that if we swap individuals across subpopulation,
  we should swap content. Otherwise, swap pointers. (There is no order
  right now within subpopulation so the later case is rare, at best.

Details:
  The default behavior is swapping all info, but not the position of
  genotypic info. If swapContent is false, pointer to genotypic info
  is swapped instead. This will lead to better performance for
  swapping but may affected performance of allele counting.

";

%ignore simuPOP::individual::shallowCopied() ;

%ignore simuPOP::individual::setShallowCopied(bool shallowCopied);

%ignore simuPOP::individual::display(ostream &out, int width=1, const vectori &chrom=vectori(), const vectori &loci=vectori());

%ignore simuPOP::individual::shallowCopiedFlagOn();

%ignore simuPOP::individual::clearShallowCopiedFlag();

%feature("docstring") simuPOP::inheritTagger "

inherite tag from parents. If both parents have tags, use fathers.

";

%feature("docstring")  simuPOP::inheritTagger::inheritTagger " 

constructor. default to be always active.

Usage:
  inheritTagger(mode=TAG_Paternal, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::inheritTagger::~inheritTagger " 

Usage:
  x.~inheritTagger()

";

%feature("docstring")  simuPOP::inheritTagger::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::inheritTagger::applyDuringMating " 

give pop, offspring, pop and mom.

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring")  simuPOP::inheritTagger::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring") simuPOP::initByFreq "

initialize genotype by allele frequency and sex by male frequency

";

%feature("docstring")  simuPOP::initByFreq::initByFreq " 

randomly assign alleles according to allele frequency

Usage:
  initByFreq(&alleleFreq=[], identicalInds=False, &subPop=[],
    indRange=intMatrix, &atLoci=[], atPloidy=-1, maleFreq=0.5,
    &sex=vectori, stage=PreMating, begin=0, end=1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL)

Arguments:

  alleleFreq:  an array of allele frequencies. Must add up to 1; or a
      matrix of allele frequencies, each row corresponse a
      subpopulation.
      
  subPop:  an array of applicable subpopulations. default to all
      
  indRange:  a [begin, end] pair of range of individuals; or an array
      of [begin, end] pairs.
      
  identicalInds:  whether or not make individual genotype identical in
      all subpopulation. If true, this operator will randomly generate
      genotype for an individual and spread it to the whole
      subpopulation.
      
  atLoci:  a vector of loci indices. If empty, apply to all loci
      
  atPloidy:  initialize which copy of chromosomes. Default to all.
      
  maleFreq:  male frequency. Default to 0.5.
      
  sex:  an arry of sex [Male, Female, Male]... for individuals. The
      length of sex will not be checked. If length of sex is shorter
      than number of individuals, sex will be reused from the
      beginning.
      
  stages:  is set to PreMating. Other parameters please see
      help(baseOperator.__init__)

Details:
  This operator randomly assign alleles according to given allele
  frequency. Allele frequencies can differ by subpop. Sex is also
  assigned randomly.

";

%feature("docstring")  simuPOP::initByFreq::~initByFreq " 

Usage:
  x.~initByFreq()

";

%feature("docstring")  simuPOP::initByFreq::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::initByFreq::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::initByFreq::apply " 

Usage:
  x.apply(&pop)

Details:
  initialize m_ranges

";

%feature("docstring") simuPOP::initByValue "

initialize genotype by value and then copy to all individuals

";

%feature("docstring")  simuPOP::initByValue::initByValue " 

Usage:
  initByValue(value=intMatrix, atLoci=[], atPloidy=-1, subPop=[],
    indRange=intMatrix, &proportions=[], maleFreq=0.5, &sex=vectori,
    stage=PreMating, begin=0, end=1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

";

%feature("docstring")  simuPOP::initByValue::~initByValue " 

Usage:
  x.~initByValue()

";

%feature("docstring")  simuPOP::initByValue::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::initByValue::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::initByValue::apply " 

Usage:
  x.apply(&pop)

Details:
  fixme: check length of src?
  atLoci is in effect

";

%feature("docstring") simuPOP::initializer "

initialize alleles at the start of generation.

Details:
  Bo Peng

";

%feature("docstring")  simuPOP::initializer::initializer " 

constructor. default to be always active.

Usage:
  initializer(&subPop=[], indRange=intMatrix, &atLoci=[],
    atPloidy=-1, maleFreq=0.5, &sex=vectori, stage=PreMating, begin=0,
    end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::initializer::~initializer " 

destructor

Usage:
  x.~initializer()

";

%feature("docstring")  simuPOP::initializer::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::initializer::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::initializer::setRanges " 

Usage:
  x.setRanges(&pop)

";

%feature("docstring")  simuPOP::initializer::initSexIter " 

Usage:
  x.initSexIter()

";

%feature("docstring")  simuPOP::initializer::nextSex " 

Usage:
  x.nextSex()

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

%feature("docstring") simuPOP::kamMutator "

K-Allele Model mutator.

Details:
  Under this model, there are K (here refers as maxAllele) possible
  allele states, and any allele has a constant probability
  (rate/(K-1)) of mutation towards any of the K-1 allelic states.

Note:
  the theoretical mutation rate is rates/(K-1) towards any of the K-1
  allelic states. So rates is actually the probability to mutate!

See Also:
  Crow & Kimura 1970

";

%feature("docstring")  simuPOP::kamMutator::kamMutator " 

K-Allele Model mutator.

Usage:
  kamMutator(rate=[], atLoci=vectori, maxAllele=0, output=\">\",
    outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL)

Arguments:

  rate:  mutation rate. It is 'probability to mutate'. The actual
      mutation rate to any of the other K-1 allelic states are
      rates/(K-1)!
      
  atLoci:  and other parameters: refer to help(mutator),
      help(baseOperator.__init__)
      
  maxAllele:  maxAllele that can be mutated to. For binary libraries
      allelic states will be [0, maxAllele]. For others, they are [1,
      maxAllele]

";

%feature("docstring")  simuPOP::kamMutator::~kamMutator " 

Usage:
  x.~kamMutator()

";

%feature("docstring")  simuPOP::kamMutator::mutate " 

mutate to a state other than current state with equal probability

Usage:
  x.mutate(allele)

";

%feature("docstring")  simuPOP::kamMutator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::kamMutator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::maPenetrance "

penetrance according to genotype at one locus

Details:
  multiple allele selector. This selector group alleles to disease and
  wild type and return penetrance to AA,Aa,aa. (A is wildtype).

";

%feature("docstring")  simuPOP::maPenetrance::maPenetrance " 

create a multiple allele selector (penetrance according to diseased
or wildtype alleles)

Usage:
  maPenetrance(loci, &penet, &wildtype, exposePenetrance=False,
    stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  loci:  the loci index.
      
  penetrance:  an array of penetrance of AA,Aa,aa. A is the wild type
      group. In the case of multiple loci, fitness should be in the
      order of BB Bb bb AA 1 2 3 Aa 4 5 6 aa 7 8 9
      
  wildtype:  an array of alleles in the wildtype group. Anything else
      is disease allele., default = [1]
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::maPenetrance::~maPenetrance " 

Usage:
  x.~maPenetrance()

";

%feature("docstring")  simuPOP::maPenetrance::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::maPenetrance::penet " 

currently assuming diploid

Usage:
  x.penet(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::maPenetrance::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::mapPenetrance "

penetrance according to genotype at one locus

Details:
  map selector. Assign penetrance value according to a given
  dictionary.

";

%feature("docstring")  simuPOP::mapPenetrance::mapPenetrance " 

create a map penetrance function (penetrance according to genotype
at one locus

Usage:
  mapPenetrance(loci, &penet, phase=False, exposePenetrance=False,
    stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  loci:  the loci index. The genotype of this locus will be axamed.
      
  penetrance:  a dictionary of penetrance. The genotype must be in the
      form of 'a-b' for single locus.
      
  phase:  if true, a/b and b/a will have different penetrance value.
      Default to false.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::mapPenetrance::~mapPenetrance " 

Usage:
  x.~mapPenetrance()

";

%feature("docstring")  simuPOP::mapPenetrance::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::mapPenetrance::penet " 

currently assuming diploid

Usage:
  x.penet(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::mapPenetrance::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::mapQuanTrait "

quantitative trait according to genotype at one locus

Details:
  map selector. Assign qtrait value according to a given dictionary.

";

%feature("docstring")  simuPOP::mapQuanTrait::mapQuanTrait " 

create a map selector (quantitative trait according to genotype at
one locus

Usage:
  mapQuanTrait(loci, &qtrait, sigma=0, phase=False,
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  loci:  the loci.
      
  qtrait:  a dictionary of qtrait. The genotype must be in the form of
      'a-b'. This is the mean of quantitative trait. The actual trait
      value will be N(mean, sigma^2) For multiple loci, the form is
      'a-b|c-d|e-f' etc.
      
  sigma:  standard deviation of the environmental facotr N(0,sigma^2).
      
  phase:  if true, a/b and b/a will have different qtrait value.
      Default to false.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::mapQuanTrait::~mapQuanTrait " 

Usage:
  x.~mapQuanTrait()

";

%feature("docstring")  simuPOP::mapQuanTrait::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::mapQuanTrait::qtrait " 

currently assuming diploid

Usage:
  x.qtrait(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::mapQuanTrait::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::mapSelector "

selection according to genotype at one locus

Details:
  map selector. Assign fitness value according to a given dictionary.

";

%feature("docstring")  simuPOP::mapSelector::mapSelector " 

create a map selector (selection according to genotype at one locus

Usage:
  mapSelector(loci, &fitness, phase=False, stage=PreMating, begin=0,
    end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  loci:  the loci index. The genotype of this locus will be axamed.
      
  fitness:  a dictionary of fitness. The genotype must be in the form
      of 'a-b' for single locus, and 'a-b|c-d|e-f' for multi-locus..
      
  phase:  if true, a/b and b/a will have different fitness value.
      Default to false.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::mapSelector::~mapSelector " 

Usage:
  x.~mapSelector()

";

%feature("docstring")  simuPOP::mapSelector::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::mapSelector::indFitness " 

currently assuming diploid

Usage:
  x.indFitness(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::mapSelector::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::maQuanTrait "

quantitative trait according to genotype at one locus

Details:
  multiple allele selector. This selector group alleles to disease and
  wild type and return qtrait to AA,Aa,aa. (A is wildtype).

";

%feature("docstring")  simuPOP::maQuanTrait::maQuanTrait " 

create a multiple allele selector (quantitative trait according to
diseased or wildtype alleles)

Usage:
  maQuanTrait(loci, &qtrait, &wildtype, &sigma=[], stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  qtrait:  an array of qtrait of AA,Aa,aa. A is the wild type group.
      
  sigma:  an array of standard deviation for each of the trait
      genotype (AA, Aa, aa)
      
  wildtype:  an array of alleles in the wildtype group. Anything else
      is disease allele. default = [1]
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::maQuanTrait::~maQuanTrait " 

destructor

Usage:
  x.~maQuanTrait()

";

%feature("docstring")  simuPOP::maQuanTrait::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::maQuanTrait::qtrait " 

currently assuming diploid

Usage:
  x.qtrait(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::maQuanTrait::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::maSelector "

selection according to genotype at one locus

Details:
  multiple allele selector. This selector group alleles to disease and
  wild type and return fitness to AA,Aa,aa. (A is wildtype).

";

%feature("docstring")  simuPOP::maSelector::maSelector " 

create a multiple allele selector (selection according to diseased
or wildtype alleles) Note that maSelectoro
nly work for diploid population now.

Usage:
  maSelector(loci, &fitness, &wildtype, stage=PreMating, begin=0,
    end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  locus:  the locus index. The genotype of this locus will be axamed.
      
  loci:  the loci index.
      
  fitness:  For the single locus case, fitness is an array of fitness
      of AA,Aa,aa. A is the wild type group. In the case of multiple
      loci, fitness should be in the order of BB Bb bb AA 1 2 3 Aa 4 5
      6 aa 7 8 9 The length for such table is 3^(#loci).
      
  wildtype:  an array of alleles in the wildtype group. Anything else
      is disease allele. default = [1] NOTE that wildtype at all loci
      are the same.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::maSelector::~maSelector " 

Usage:
  x.~maSelector()

";

%feature("docstring")  simuPOP::maSelector::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::maSelector::indFitness " 

currently assuming diploid

Usage:
  x.indFitness(*ind)

Details:
  get genotype of ind

";

%feature("docstring")  simuPOP::maSelector::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::mating "

Details:
  The mating classes describe various mating scheme --- a required
  parameter of simulator.

";

%feature("docstring")  simuPOP::mating::isCompatible " 

check if the mating type is compatible with population structure

Usage:
  x.isCompatible(&pop)

Details:
  possible things to check:need certain types of individual (age,
  sex etc)
  need resizeable population...

";

%feature("docstring")  simuPOP::mating::mating " 

constructor

Usage:
  mating(numOffspring=1.0, *numOffspringFunc=NULL,
    maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
    newSubPopSizeExpr=\"\", *newSubPopSizeFunc=NULL)

";

%feature("docstring")  simuPOP::mating::mating " 

Usage:
  mating(&rhs)

";

%feature("docstring")  simuPOP::mating::~mating " 

destructor

Usage:
  x.~mating()

";

%feature("docstring")  simuPOP::mating::clone " 

 clone() const. Generate a copy of itself and return a pointer

Usage:
  x.clone()

Details:
  This function is important because Python automatically release an
  object after it is used.
  For example:will fail since mating()i
  s released after the first line being executed.
  With the help of clone() const, the C++
  implementation can avoid this problem by

";

%feature("docstring")  simuPOP::mating::__repr__ " 

return name of the mating type used primarily in logging.

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::mating::submitScratch " 

Usage:
  x.submitScratch(&pop, &scratch)

";

%feature("docstring")  simuPOP::mating::mate " 

mate: This is not supposed to be called for base mating class.

Usage:
  x.mate(&pop, &scratch, &ops, submit)

Arguments:

  pop:  population
      
  scratch:  scratch population
      
  ops:  during mating operators

Value:
  return false when mating fail.

";

%feature("docstring")  simuPOP::mating::numOffspring " 

Usage:
  x.numOffspring(gen)

";

%feature("docstring")  simuPOP::mating::resetNumOffspring " 

Usage:
  x.resetNumOffspring()

";

%feature("docstring")  simuPOP::mating::formOffGenotype " 

whether or not to generate offspring genotype this is true when none
of the during-mating operator can do this.

Usage:
  x.formOffGenotype(&ops)

";

%feature("docstring")  simuPOP::mating::prepareScratchPop " 

dealing with pop/subPop size change, copy of structure etc.

Usage:
  x.prepareScratchPop(&pop, &scratch)

Details:
  force new subPop size?

";

%feature("docstring") simuPOP::mergeSubPops "

merge subpopulations

";

%feature("docstring")  simuPOP::mergeSubPops::mergeSubPops " 

merge subpopulations

Usage:
  mergeSubPops(subPops=[], removeEmptySubPops=False,
    stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  subPops:  subpops to be merged, default to all subpops.

";

%feature("docstring")  simuPOP::mergeSubPops::~mergeSubPops " 

destructor

Usage:
  x.~mergeSubPops()

";

%feature("docstring")  simuPOP::mergeSubPops::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::mergeSubPops::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::mergeSubPops::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::migrator "

Details:
  Bo Peng

";

%feature("docstring")  simuPOP::migrator::migrator " 

create a migrator

Usage:
  migrator(&rate, mode=MigrByProbability, fromSubPop=[],
    toSubPop=[], stage=PreMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL)

Arguments:

  rate:  migration rate, proportion or count. Determined by parameter
      mode. rate should be a m by n matrix. If a number is given, the
      migration rate will be r*ones(m,n).
      
  mode:  one of MigrByProbability (default), MigrByProportion or
      MigrByCounts
      
  fromSubPop:  an array of 'from' subpops, default to all
      subpopulations. If a single subpop is specified, [] can be
      ignored. I.e., [a] is equvalent to a.
      
  toSubPop:  an array of 'to' subpops, default to all subpopulations.
      If a single subpop is specified, [] can be ignored.
      
  stage:  is default to PreMating. For details about other parameters,
      please refer to help(baseOperator.__init__)
      
  rate is a matrix with dimensions determined by fromSubPop and
  toSubPop. By default, rate is a matrix with element (i,j) being the
  migration rate, probability or count from subpop i to subpop j. If
  fromSubPop and/or toSubPop are given, migration only happen between
  these subpopulations. An extreme case is 'point migration'
  rate=[[r]], fromSubPop=a, toSubPop=b
  which migrate from subpop a to b with given rate r.

";

%feature("docstring")  simuPOP::migrator::~migrator " 

destructor

Usage:
  x.~migrator()

";

%feature("docstring")  simuPOP::migrator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::migrator::rate " 

return rate

Usage:
  x.rate()

";

%feature("docstring")  simuPOP::migrator::setRates " 

set migration rate

Usage:
  x.setRates(&rate, mode)

Details:
  format 0-0 0-1 0-2, 1-0 1-1 1-2, 2-0, 2-1, 2-2. for mode 1 or 2,
  00,11,22 will be set automatically. regardless of input.
  set r[i][i]--- may need to extend rate (to add i->i)

";

%feature("docstring")  simuPOP::migrator::apply " 

Usage:
  x.apply(&pop)

Details:
  2nd, or 3rd method
  create a vector and assign indices, then random shuffle and assign
  info
  for all subPop.

";

%feature("docstring")  simuPOP::migrator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::mlPenetrance "

penetrance according to genotype at multiple loci multiplicative
model

Details:
  multiple loci selector. This selector takes several selectors and
  multiply their penetrance values... e.g. mlmpenetrance(
  [mappenetrance(...), mapenetrance(...) ])

";

%feature("docstring")  simuPOP::mlPenetrance::mlPenetrance " 

multiple loci selector using a multiplicative model.

Usage:
  mlPenetrance(peneOps, mode=PEN_Multiplicative,
    exposePenetrance=False, stage=DuringMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  selectors:  a list of selectors.
      
  mode:  one of PEN_Multiplicative, PEN_Additive, PEN_Heterogeneity

";

%feature("docstring")  simuPOP::mlPenetrance::~mlPenetrance " 

Usage:
  x.~mlPenetrance()

";

%feature("docstring")  simuPOP::mlPenetrance::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::mlPenetrance::penet " 

currently assuming diploid

Usage:
  x.penet(*ind)

";

%feature("docstring")  simuPOP::mlPenetrance::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::mlQuanTrait "

quantitative trait according to genotype at multiple loci
multiplicative model

Details:
  multiple loci selector. This selector takes several selectors and
  multiply their qtrait values... e.g. mlmquanTrait(
  [mapquanTrait(...), maquanTrait(...) ])

";

%feature("docstring")  simuPOP::mlQuanTrait::mlQuanTrait " 

multiple loci selector using a multiplicative model.

Usage:
  mlQuanTrait(qtraits, mode=QT_Multiplicative, sigma=0,
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  qtraits:  a list of qtraits.

";

%feature("docstring")  simuPOP::mlQuanTrait::~mlQuanTrait " 

Usage:
  x.~mlQuanTrait()

";

%feature("docstring")  simuPOP::mlQuanTrait::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::mlQuanTrait::qtrait " 

currently assuming diploid

Usage:
  x.qtrait(*ind)

";

%feature("docstring")  simuPOP::mlQuanTrait::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::mlSelector "

selection according to genotype at multiple loci multiplicative
model

Details:
  multiple loci selector. This selector takes several selectors and
  multiply their fitness values... e.g. mlmselector(
  [mapselector(...), maselector(...) ])

";

%feature("docstring")  simuPOP::mlSelector::mlSelector " 

multiple loci selector using a multiplicative model.

Usage:
  mlSelector(selectors, mode=SEL_Multiplicative, stage=PreMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  selectors:  a list of selectors.

";

%feature("docstring")  simuPOP::mlSelector::~mlSelector " 

Usage:
  x.~mlSelector()

";

%feature("docstring")  simuPOP::mlSelector::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::mlSelector::indFitness " 

currently assuming diploid

Usage:
  x.indFitness(*ind)

Details:
  fixme

";

%feature("docstring")  simuPOP::mlSelector::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::mutator "

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

%feature("docstring")  simuPOP::mutator::mutator " 

create a mutator All mutators have the following common parameters.
However, the actual meaning of these parameters may vary according
to different model. Check the manual for details!!!
(help(kamMutator) for example.)

Usage:
  mutator(rate=[], atLoci=vectori, maxAllele=0, output=\">\",
    outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL)

Arguments:

  rate:  single rate for all applicable loci (atLoci). Will be ignored
      if rates is specified; or it can be an array of rates, the same
      length as atLoci.
      
  atLoci:  a vector of loci index. Can be ignored only when single
      rate is specified. Default to all loci.
      
  maxAllele:  max allowable allele. Interpreted by each sub mutaor
      class. Default to pop.maxAllele().

";

%feature("docstring")  simuPOP::mutator::~mutator " 

destructor

Usage:
  x.~mutator()

";

%feature("docstring")  simuPOP::mutator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::mutator::rate " 

return mutation rate

Usage:
  x.rate()

";

%feature("docstring")  simuPOP::mutator::setRate " 

set an array of rates

Usage:
  x.setRate(rate, atLoci=vectori)

";

%feature("docstring")  simuPOP::mutator::maxAllele " 

return max allowable allele number

Usage:
  x.maxAllele()

";

%feature("docstring")  simuPOP::mutator::setMaxAllele " 

Usage:
  x.setMaxAllele(maxAllele)

";

%feature("docstring")  simuPOP::mutator::mutationCount " 

return mutation count

Usage:
  x.mutationCount(locus)

";

%feature("docstring")  simuPOP::mutator::mutationCounts " 

return mutation counts

Usage:
  x.mutationCounts()

";

%feature("docstring")  simuPOP::mutator::mutate " 

how to mutate a single allele. this is usually the only function
that need to be defined by the subclasses.

Usage:
  x.mutate(allele)

";

%feature("docstring")  simuPOP::mutator::apply " 

apply!

Usage:
  x.apply(&pop)

";

%feature("docstring") simuPOP::noMating "

Details:
  No mating. No subpopulation change. During mating operator will be
  applied, but the return values are not checked.

";

%feature("docstring")  simuPOP::noMating::noMating " 

constructor, no new subPopsize parameter

Usage:
  noMating()

";

%feature("docstring")  simuPOP::noMating::~noMating " 

destructor

Usage:
  x.~noMating()

";

%feature("docstring")  simuPOP::noMating::clone " 

 clone() const. The same as mating::clone() const.

Usage:
  x.clone()

See Also:
   mating::clone() const

";

%feature("docstring")  simuPOP::noMating::__repr__ " 

return name of the mating type

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::noMating::submitScratch " 

Usage:
  x.submitScratch(&pop, &scratch)

";

%feature("docstring")  simuPOP::noMating::mate " 

do the mating. --- no mating :-)

Usage:
  x.mate(&pop, &scratch, &ops, submit)

Details:
  All individuals will be passed to during mating operators but no one
  will die (ignore during mating failing signal).

";

%feature("docstring") simuPOP::noneOp "

";

%feature("docstring")  simuPOP::noneOp::noneOp " 

do nothing

Usage:
  noneOp(output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
    end=0, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::noneOp::~noneOp " 

destructor

Usage:
  x.~noneOp()

";

%feature("docstring")  simuPOP::noneOp::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::noneOp::applyWithScratch " 

simply output some info providing interface to apply operator before
during or after mating.

Usage:
  x.applyWithScratch(&pop, &scratch, stage)

";

%feature("docstring")  simuPOP::noneOp::applyDuringMating " 

give pop, offspring, pop and mom.

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring")  simuPOP::noneOp::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::noneOp::__repr__ " 

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

%feature("docstring") simuPOP::Operator "

base class of all classes that manipulate populations.

Details:
  Operators are object that act on populations. They can be applied to
  populations directly using apply()member f
  unction, but most of the time they are managed and applied by a
  simulator.
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
  table-like outputs.
  filename can have the following format:
  1. 'filename' this file will be closed after each use. I.e., if
  several operators output to the same file, only the last one will
  succeed.
  2. '>filename' the same as 'filaname'
  3. '>>filename' The file will be created at the beginning of
  evolution ( simulator::evolve) and close at the end.
  Several operators can output to this file to form a table.
  4. '>>>filename' The same as '>>filename' except that the file will
  not be cleared at the beginning of evolution if it is not empty.
  5. '>' out put to standard output.
  6. '' supress output.
  Most operators are applied to every replicate of a simulator during
  evolution. However, you can apply opertors to one or a group of
  replicates only. For example, you can initialize different
  replicates with different initial values and then start evolution.
  c.f. simulator::setGroup.
  Please refer to help(baseOperator) and help(baseOperator.__init__)
  for detailed information about member functions and parameters.
  Bo Peng

";

%feature("docstring")  simuPOP::Operator::Operator " 

create an operator (this function is not supposed to be called
directly)

Usage:
  operator(output, outputExpr, stage, begin, end, step, at, rep,
    grp)

Arguments:

  output:  a string of output filename. Different operators will have
      different default output (most commonly '>' or '')
      
  outputExpr:  an expression that determines output filename
      dynamically.
      
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
  set in simulator::setGroup()

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

%ignore simuPOP::Operator::atRepr();

%feature("docstring")  simuPOP::Operator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::Operator::noOutput " 

if output=\">\". used internally

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

%feature("docstring") simuPOP::outputer "

outputer is a (special) subclass of Operatort
hat will output files with different format.

Details:
  Bo Peng

";

%feature("docstring")  simuPOP::outputer::outputer " 

constructor. default to be always active.

Usage:
  outputer(output=\">\", outputExpr=\"\", stage=PostMating, begin=0,
    end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::outputer::~outputer " 

destructor

Usage:
  x.~outputer()

";

%feature("docstring")  simuPOP::outputer::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring") simuPOP::outputHelper "

";

%feature("docstring")  simuPOP::outputHelper::outputHelper " 

Usage:
  outputHelper(str=\"\\\\n\", output=\">\", outputExpr=\"\",
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

";

%feature("docstring")  simuPOP::outputHelper::apply " 

simply output some info

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::outputHelper::~outputHelper " 

Usage:
  x.~outputHelper()

";

%feature("docstring")  simuPOP::outputHelper::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::outputHelper::setString " 

set output string.

Usage:
  x.setString(str)

";

%feature("docstring")  simuPOP::outputHelper::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::parentsTagger "

inherite tag from parents. If both parents have tags, use fathers.

";

%feature("docstring")  simuPOP::parentsTagger::parentsTagger " 

constructor. default to be always active. string can be any string
(m_Delimiter will be ignored for this class.) r will be replicate
number g will be generation number.

Usage:
  parentsTagger(begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

";

%feature("docstring")  simuPOP::parentsTagger::~parentsTagger " 

Usage:
  x.~parentsTagger()

";

%feature("docstring")  simuPOP::parentsTagger::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::parentsTagger::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::parentsTagger::applyDuringMating " 

give pop, offspring, pop and mom.

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring") simuPOP::pause "

";

%feature("docstring")  simuPOP::pause::pause " 

stop simulation. press q to exit and any other key to continue

Usage:
  pause(prompt=True, stopOnKeyStroke=False, exposePop=True,
    popName=\"pop\", output=\">\", outputExpr=\"\", stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_LAST, grp=GRP_ALL)

Arguments:

  prompt:  if true (default), print prompt message.
      
  stopOnKeyStroke:  if true, goon if no key was pressed
      
  exposePop:  whether or not expose pop to user namespace, only useful
      when user choose 's' at pause. Default to true.
      
  popName:  by which name the population is exposed? default to 'pop'

";

%feature("docstring")  simuPOP::pause::~pause " 

destructor

Usage:
  x.~pause()

";

%feature("docstring")  simuPOP::pause::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pause::apply " 

simply output some info

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::pause::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::penetrance "

penetrance

Details:
  Please refer to the user's guide for details.

";

%feature("docstring")  simuPOP::penetrance::penetrance " 

constructor. default to be always active. default to post mating

Usage:
  penetrance(exposePenetrance=False, stage=DuringMating, begin=0,
    end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::penetrance::~penetrance " 

destructor

Usage:
  x.~penetrance()

";

%feature("docstring")  simuPOP::penetrance::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::penetrance::penet " 

calculate/return penetrance etc

Usage:
  x.penet(*ind)

";

%feature("docstring")  simuPOP::penetrance::apply " 

set pentrance to all individuals and record penetrance if requested.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::penetrance::applyDuringMating " 

set penetrance to all individual

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring")  simuPOP::penetrance::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::pointMutator "

point mutator

Details:
  mutate specified individuals at specified loci to spcified allele.
  I.e., this is a non-random mutator used to introduce disease etc.

";

%feature("docstring")  simuPOP::pointMutator::pointMutator " 

mutate once

Usage:
  pointMutator(atLoci, toAllele, atPloidy=[], inds=[], output=\">\",
    outputExpr=\"\", stage=PostMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL)

Arguments:

  atLoci:  a vector of loci index.
      
  inds:  mutate 'inds' individuals
      
  toAllele:  mutate to 'toAllele'

";

%feature("docstring")  simuPOP::pointMutator::~pointMutator " 

destructor

Usage:
  x.~pointMutator()

";

%feature("docstring")  simuPOP::pointMutator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pointMutator::apply " 

apply!

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::pointMutator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::pointMutator::mutationCount " 

return mutation count

Usage:
  x.mutationCount(locus)

";

%feature("docstring")  simuPOP::pointMutator::mutationCounts " 

return mutation counts

Usage:
  x.mutationCounts()

";

%feature("docstring") simuPOP::population "

a collection of individuals with subpopulation structure

Details:
  Please refer to user's Guide for details about this object.

";

%feature("docstring")  simuPOP::population::population " 

create a population object with given size and genotypic structure

Usage:
  population(size=0, ploidy=2, &loci=[], sexChrom=False,
    &lociPos=[], &subPop=[], ancestralDepth=0, &alleleNames=[],
    &lociNames=[], maxAllele=MaxAllele)

Arguments:

  size:  population size. Can be ignored if subPop is specified. In
      that case, size is sum of subPop. Default to 0.
      
  ploidy:  number of sets of chromosomes. Default to 2 (diploid).
      
  loci:  an array of numbers of loci on each chromosome. If not
      specified, assume a single locus on one chromosome. Number of
      chromosomes is determined by the size of this array.
      
  lociPos:  an array of loci distance for each locus. You can also use
      a nested array to specify loci distance for each chromosome. (
      [1,2,3,4,5] or [[1,2],[3,4,5]] are both allowed for loci=[2,3])
      The default values are 1, 2, etc. on each chromosome.
      
  subPop:  an array of subpopulation sizes. Default value is [size]
      which means a single subpopulation of the whole population. If
      both size and subPop are given, sum of subPop should agree with
      size.
      
  ancestralDepth:  number of ancestral populations to keep. Default to
      0, meaning only current generation will be available. If -1 is
      given, all ancestral populations will be saved. This may exhaust
      your RAM pretty quickly though.
      
  alleleNames:  an array of allele names. The first element should be
      given to invalid/unknown allele. For example, for a locus with
      alleles A,C,T,G, you can specify alleleName as
      ('_','A','C','T','G'). Note that simuPOPu
      ses 1,2,3,4 internally and these names will only be used for
      display purpose.
      
  lociNames:  an array or a matrix (separated by chromosome) of names
      for each loci. Default to \"locX-X\" where X-X is chromosome-
      loci index starting from 1. This info is rarely used.
      
  maxAllele:  maximum allele number. Default to the max allowed allele
      states of current library (standard or long allele version)

Value:
  no return value. Exceptionwill be thrown is
  wrong parameters are given.

See Also:
   simulator, baseOperator, matings
  chemes

Examples:
  popInit.log does not exist

";

%ignore simuPOP::population::population(const population &rhs);

%feature("docstring")  simuPOP::population::clone " 

Usage:
  x.clone(keepAncestralPops=True)

";

%feature("docstring")  simuPOP::population::swap " 

SWAP population swap the content of two populations.

Usage:
  x.swap(&rhs)

";

%feature("docstring")  simuPOP::population::~population " 

destroy a population

Usage:
  x.~population()

";

%feature("docstring")  simuPOP::population::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::population::__cmp__ " 

Usage:
  x.__cmp__()

";

%feature("docstring")  simuPOP::population::setSubPopStru " 

set population/subpopulation given subpopulation sizes subPopSize an
array of subpopulation sizes the population may or may not change
according to parameter allowPopSizeChange if sum of subPopSize does
not match popSize. allowPopSizeChange if true, popSize can change to
sum of subPopSize. none migration, mating

Usage:
  x.setSubPopStru(&newSubPopSizes, allowPopSizeChange=False)

";

%feature("docstring")  simuPOP::population::numSubPop " 

number of sub populations.

Usage:
  x.numSubPop()

Value:
  number of subpopulations (>=1)

";

%feature("docstring")  simuPOP::population::subPopSize " 

get size of subpopulation subPop

Usage:
  x.subPopSize()

Arguments:

  subPop:  index of subpopulation (start from 0)

Value:
  size of subpopulation subPop

";

%feature("docstring")  simuPOP::population::subPopSizes " 

get size of all subpopulations

Usage:
  x.subPopSizes()

Value:
  an array of size of subpopulations

";
%feature("docstring")  simuPOP::population::popSize " 

get population size

Usage:
  x.popSize()

Value:
  total number of individuals in this population

";

%feature("docstring")  simuPOP::population::absIndIndex " 

absolute index of individual at a subpopulation

Usage:
  x.absIndIndex(ind, )

Arguments:

  index:  index of individual at subpopulation subPop
      
  subPop:  subpopulation index

Value:
  absolute index of individual at subpopulation subPop

See Also:
   subPopIndPair

";

%feature("docstring")  simuPOP::population::subPopIndPair " 

subPop and relative index of an individual

Usage:
  x.subPopIndPair(ind)

";

%feature("docstring")  simuPOP::population::subPopBegin " 

beginning index for subpopulation subPop

Usage:
  x.subPopBegin()

Arguments:

  subPop:  subpopulation index

Details:

  absIndIndex(index, subPop) = subPopBegin(subPop) + index

Value:
  beginning index of this subpopulation

See Also:
   absIndIndex

";

%feature("docstring")  simuPOP::population::subPopEnd " 

ending index for subpopulation subPop

Usage:
  x.subPopEnd()

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
%feature("docstring")  simuPOP::population::ind " 

refernce to individual ind in subpopulation subPop

Usage:
  x.ind(ind, subPop=0)

Arguments:

  ind:  individual index within subPop
      
  subPop:  subpopulation index

Value:
  reference to an individual

";

%ignore simuPOP::population::indBegin();

%ignore simuPOP::population::indEnd();

%ignore simuPOP::population::indBegin(UINT subPop);

%ignore simuPOP::population::indEnd(UINT subPop);

%ignore simuPOP::population::alleleBegin(UINT locus);

%ignore simuPOP::population::alleleEnd(UINT locus);

%ignore simuPOP::population::alleleBegin(UINT locus, UINT subPop);

%ignore simuPOP::population::alleleEnd(UINT locus, UINT subPop);

%ignore simuPOP::population::genoBegin();

%ignore simuPOP::population::genoEnd();

%ignore simuPOP::population::genoBegin(UINT subPop);

%ignore simuPOP::population::genoEnd(UINT subPop);

%ignore simuPOP::population::indGenoBegin(ULONG ind) ;

%ignore simuPOP::population::indGenoEnd(ULONG ind) ;

%ignore simuPOP::population::indGenoBegin(ULONG ind, UINT subPop) ;

%ignore simuPOP::population::indGenoEnd(ULONG ind, UINT subPop) ;

%feature("docstring")  simuPOP::population::arrGenotype " 

get the whole genotype. individuals will be in order before exposing
their genotypes.

Usage:
  x.arrGenotype()

";

%feature("docstring")  simuPOP::population::arrGenotype " 

get the whole genotype. individuals will be in order before exposing
their genotypes.

Usage:
  x.arrGenotype(subPop)

";
%feature("docstring")  simuPOP::population::exposeIndInfo " 

Usage:
  x.exposeIndInfo(name=\"info\")

Details:
  brief return individual info in pop namespace

";

%feature("docstring")  simuPOP::population::exposeAffectedness " 

Usage:
  x.exposeAffectedness(name=\"affected\")

Details:
  brief return individual affected status in pop namespace

";

%feature("docstring")  simuPOP::population::setIndInfo " 

Usage:
  x.setIndInfo(&info)

Arguments:

  info:  an array of info values, should have length of pop size

Details:
  This function set info field of all individuals. Info can be used to
  sort individuals and set subpopulations. Therefore, the following
  code will do a migration:
  setIndInfo([ an_array_of_info_value ])
   setSubPopByIndInfo()

See Also:
   individual::setInfo, individual::info,
   info

";

%feature("docstring")  simuPOP::population::setIndInfoWithSubPopID " 

set individual info with their subpopulation id.

Usage:
  x.setIndInfoWithSubPopID()

Details:
  set individual info by subpop id.

";

%feature("docstring")  simuPOP::population::setSubPopByIndInfo " 

adjust subpopulation according to individual info values

Usage:
  x.setSubPopByIndInfo(info=vectori)

Arguments:

  info:  optional info that will be used to set sub pop

Details:
  assume individual has subpop index in their info value and put them
  into corresponding subpopulations.

  rebuild index

Note:
  individual with negative info will be removed!

See Also:
  setInfo,

";

%feature("docstring")  simuPOP::population::splitSubPop " 

split population

Usage:
  x.splitSubPop(which, sizes, subPopID=[])

Details:
  split subpopulation 'which' into subpopulations with size specified
  in subPops, optionally given subPopID. The subPOP ID of Other subpop
  will be kept. For example, if subpop 1 of 0 1 2 3 is split into
  three parts, the new subpop id will be 0 (1 4 5) 2 3.

Note:
  subpop with negative id will be removed. So, you can shrink one
  subpop by split and set one of the new subpop with negative id.

";

%feature("docstring")  simuPOP::population::splitSubPopByProportion " 

split population

Usage:
  x.splitSubPopByProportion(which, proportions, subPopID=[])

Details:
  split subpopulation 'which' into subpopulations with specified
  'proportions', optionally given subPopID.

Note:
  subpop with negative id will be removed. So, you can shrink one
  subpop by split and set one of the new subpop with negative id.

";

%feature("docstring")  simuPOP::population::removeEmptySubPops " 

Usage:
  x.removeEmptySubPops()

Details:
  remove empty subpops, this will adjust subPOP ID of other subpops
  rebuild index

";

%feature("docstring")  simuPOP::population::removeSubPops " 

Usage:
  x.removeSubPops(&subPops=[], shiftSubPopID=True,
    removeEmptySubPops=False)

Details:
  remove subpop, adjust subpop numbers so that there will be no
  'empty' subpops left

";

%feature("docstring")  simuPOP::population::removeIndividuals " 

Usage:
  x.removeIndividuals(&inds=[], subPop=-1, removeEmptySubPops=False)

Details:
  remove subpop, adjust subpop numbers so that there will be no
  'empty' subpops left

";

%feature("docstring")  simuPOP::population::mergeSubPops " 

merge population

Usage:
  x.mergeSubPops(subPops=[], removeEmptySubPops=False)

Details:
  merge subpopulations, subpop id will be the ID of the first in array
  subPops all subpopulation will take the id of the first one.

";

%feature("docstring")  simuPOP::population::reorderSubPops " 

reorder subpopulations

Usage:
  x.reorderSubPops(&order=[], &rank=[], removeEmptySubPops=False)

Arguments:

  order:  new order of the subpopulations. For examples, 3 2 0 1 means
      subpop3, subpop2, subpop0, subpop1 will be the new layout.
      
  rank:  you can also specify new rank for each subpop. For example,
      3,2,0,1 means the original subpopulations will have new ID
      3,2,0,1. To achive order 3,2,0,1. the rank should be 1 0 2 3.

";

%feature("docstring")  simuPOP::population::newPopByIndInfo " 

Usage:
  x.newPopByIndInfo(keepAncestralPops=True, info=vectori,
    removeEmptySubPops=False)

Details:
  form a new population according to info, info can be given directly

";

%feature("docstring")  simuPOP::population::removeLoci " 

Usage:
  x.removeLoci(&remove=[], &keep=[])

";

%feature("docstring")  simuPOP::population::newPopWithPartialLoci " 

Usage:
  x.newPopWithPartialLoci(&remove=[], &keep=[])

Details:
  get a new population with selected loci

";

%feature("docstring")  simuPOP::population::fitness " 

Usage:
  x.fitness()

";

%feature("docstring")  simuPOP::population::arrFitness " 

Usage:
  x.arrFitness()

";

%feature("docstring")  simuPOP::population::pushAndDiscard " 

Usage:
  x.pushAndDiscard(&rhs, force=False)

";

%feature("docstring")  simuPOP::population::ancestralDepth " 

Usage:
  x.ancestralDepth()

";

%feature("docstring")  simuPOP::population::setAncestralDepth " 

set ancestral depth, can be -1

Usage:
  x.setAncestralDepth(depth)

";

%feature("docstring")  simuPOP::population::ancestralPop " 

Usage:
  x.ancestralPop()

";

%feature("docstring")  simuPOP::population::useAncestralPop " 

Usage:
  x.useAncestralPop(idx)

";

%feature("docstring")  simuPOP::population::equalTo " 

compare two populations

Usage:
  x.equalTo(&rhs)

";

%ignore simuPOP::population::adjustGenoPosition(bool deep=false);

%feature("docstring")  simuPOP::population::savePopulation " 

save population to a file

Usage:
  x.savePopulation(&filename, &format=\"auto\")

Arguments:

  filename:  save to filename
      
  format:  format to save. Can be one of 'text', 'bin', 'xml'
      
  The default format is 'text' but the output is not suppored to be
  read. 'bin' has smaller size and should be used for large
  populations. 'xml' format is most readable and should be used when
  you would like to convert simuPOPpopulatio
  ns to other formats.

See Also:
  global function loadPopulation

";

%ignore simuPOP::population::loadPopulation(const string &filename, const string &format="auto");

%feature("docstring")  simuPOP::population::rep " 

Usage:
  x.rep()

";

%ignore simuPOP::population::setRep(int rep, bool setVar=true);

%feature("docstring")  simuPOP::population::grp " 

Usage:
  x.grp()

";

%ignore simuPOP::population::setGrp(int grp, bool setVar=true);

%feature("docstring")  simuPOP::population::gen " 

Usage:
  x.gen()

";

%ignore simuPOP::population::setGen(ULONG gen, bool setVar=true);

%feature("docstring")  simuPOP::population::vars " 

return variables of this population if subPop is given, return
dictionary for specified subpopulation.

Usage:
  x.vars(subPop=-1)

";

%ignore simuPOP::population::dict(int subPop=-1);

%ignore simuPOP::population::setDict(PyObject *dict);

%feature("docstring")  simuPOP::population::hasVar " 

Usage:
  x.hasVar(&name)

";

%feature("docstring")  simuPOP::population::removeVar " 

CPPNLY.

Usage:
  x.removeVar(&name)

";

%feature("docstring")  simuPOP::population::checkRefCount " 

Usage:
  x.checkRefCount()

";

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

%ignore simuPOP::population::varsAsString() ;

%ignore simuPOP::population::varsFromString(const string &vars);

%feature("docstring")  simuPOP::population::evaluate " 

evaluate python statment/expressions

Usage:
  x.evaluate(&expr=\"\", &stmts=\"\")

Details:
  this function evaluate python expressions and return as string
  representing the result

";

%feature("docstring")  simuPOP::population::execute " 

Usage:
  x.execute(&stmts=\"\")

";

%feature("docstring") simuPOP::pyEval "

evaluate an expression.

";

%feature("docstring")  simuPOP::pyEval::pyEval " 

evaluate expr/statments in local replicate namespace

Usage:
  pyEval(&expr=\"\", &stmts=\"\", &preStmts=\"\", &postStmts=\"\",
    exposePop=False, &name=\"\", output=\">\", outputExpr=\"\",
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  expr:  expression to be evaluated. Its result will be sent to
      output.
      
  stmts:  statement that will be executed before the expression.
      
  preStmts:  statement that will be executed when the operaot is
      constructed.
      
  postStmts:  statement that will be executed when the operator is
      destroyed.
      
  exposePop:  if true, expose current pop as variable \"pop\"
      
  name:  used to let pure python operator to identify themselves
      
  output:  default to \">\" . i.e., output to standard output. for
      usage of other parameters, see help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::pyEval::~pyEval " 

Usage:
  x.~pyEval()

";

%feature("docstring")  simuPOP::pyEval::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pyEval::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::pyEval::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::pyEval::name " 

Usage:
  x.name()

";

%feature("docstring") simuPOP::pyExec "

evaluate an expression.

";

%feature("docstring")  simuPOP::pyExec::pyExec " 

evaluate statments in local replicate namespace, no return value

Usage:
  pyExec(&stmts=\"\", &preStmts=\"\", &postStmts=\"\",
    exposePop=False, &name=\"\", output=\">\", outputExpr=\"\",
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  stmts:  statements (a single or multi-line string) that will be
      evaluated before the expression.
      
  preStmts:  statement that will be executed when the operaot is
      constructed.
      
  postStmts:  statement that will be executed when the operator is
      destroyed.
      
  exposePop:  if true, expose current pop as variable \"pop\"
      
  output:  default to \">\" . i.e., output to standard output. for
      usage of other parameters, see help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::pyExec::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pyExec::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::pyInit "

";

%feature("docstring")  simuPOP::pyInit::pyInit " 

initialize populations using given user function.

Usage:
  pyInit(*func, subPop=[], atLoci=[], atPloidy=-1,
    indRange=intMatrix, maleFreq=0.5, &sex=vectori, stage=PreMating,
    begin=0, end=1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  func:  a python function with parameter (index, ploidy, subpop)
      index is the allele index (0 ~ totNumLoci()-1), ploidy (index to
      copy of chromosomes), subpop (subpop number). The return value
      of this function should be a integer.
      
  atLoci:  a vector of loci indices. If empty, apply to all loci
      
  atPloidy:  initialize which copy of chromosomes. Default to all.
      
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

%feature("docstring")  simuPOP::pyInit::~pyInit " 

Usage:
  x.~pyInit()

";

%ignore simuPOP::pyInit::pyInit(const pyInit &rhs);

%feature("docstring")  simuPOP::pyInit::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pyInit::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::pyInit::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring") simuPOP::pyMating "

Details:
  Hybrid mating scheme.

";

%feature("docstring")  simuPOP::pyMating::pyMating " 

constructor

Usage:
  pyMating(*mateFunc, keepSubPopStru=True)

Arguments:

  mateFunc:  a python function that takes a population as parameter
      and return a list of parental indices. Currently, only diploid
      population is supported.
      
  keepSubPopStru:  if true, mating is strictly between subpop and
      subpopulation structure will be kept. Otherwise, mating is
      considered as in the big population regardless of subpopulation
      strcture. The resulting population does not have subpopulation
      strcture.
      
  Note that these indices should be absolte indices and mating across
  subpops is not allowed.
  In this way, you can organize arbitrary complex mating scheme (but
  also with considerable work.)

Details:
  This operator takes only one parameter: a mate function. During
  mating, this function will be called with pop as parameter. The
  expected return value should be a list (or carray) of indices to
  parents in the order of dad1,mom1,dad2,mom2,...

";

%feature("docstring")  simuPOP::pyMating::~pyMating " 

destructor

Usage:
  x.~pyMating()

";

%feature("docstring")  simuPOP::pyMating::pyMating " 

Usage:
  pyMating(&rhs)

";

%feature("docstring")  simuPOP::pyMating::clone " 

 clone() const. The same as mating::clone() const.

Usage:
  x.clone()

See Also:
   mating::clone() const

";

%feature("docstring")  simuPOP::pyMating::__repr__ " 

return name of the mating type

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::pyMating::mate " 

do the mating with specified mating function.

Usage:
  x.mate(&pop, &scratch, &ops)

Details:
  All individuals will be passed to during mating operators but no one
  will die (ignore during mating failing signal).

";

%feature("docstring") simuPOP::pyMigrator "

migrate using given info vector

Details:
  You can use directmigrator to accomplish any migration: that is to
  say you directly specify subpopulation numbers for each individual
  and this operator will do the rest.

";

%feature("docstring")  simuPOP::pyMigrator::pyMigrator " 

create a directmigrator

Usage:
  pyMigrator(*subPopID=NULL, stage=PreMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  subPopID:  a 1-d array (list, typle, carray). Must has length
      greater or equal to population size.
      
  stage:  is default to PreMating, please refer to
      help(baseOperator.__init__) for details about other parameters.

Details:
  This operator accept a one-dimensional Numeric Python int array.
  (created by Numeric.array ). The contend of the array will be
  considered as subpopulation id.

";

%feature("docstring")  simuPOP::pyMigrator::~pyMigrator " 

destructor

Usage:
  x.~pyMigrator()

";

%ignore simuPOP::pyMigrator::pyMigrator(const pyMigrator &rhs);

%feature("docstring")  simuPOP::pyMigrator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pyMigrator::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::pyMigrator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::pyMutator "

mixed mutation model . has not been implemented.

";

%feature("docstring")  simuPOP::pyMutator::pyMutator " 

Usage:
  pyMutator(rate=[], atLoci=vectori, maxAllele=0, *func=NULL,
    output=\">\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::pyMutator::~pyMutator " 

Usage:
  x.~pyMutator()

";

%ignore simuPOP::pyMutator::pyMutator(const pyMutator &rhs);

%feature("docstring")  simuPOP::pyMutator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pyMutator::mutate " 

how to mutate a single allele. this is usually the only function
that need to be defined by the subclasses.

Usage:
  x.mutate(allele)

";

%feature("docstring")  simuPOP::pyMutator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::pyOperator "

";

%feature("docstring")  simuPOP::pyOperator::pyOperator " 

python operator, using a function that accept a population object

Usage:
  pyOperator(*func, *param=NULL, stage=PostMating,
    formOffGenotype=False, passOffspringOnly=False, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  func:  a python function that accept a population and perform
      whatever operation it wants to.
      
  para:  any python object that will be passed to func after pop
      parameter. Multiple parameter can be passed as a tuple.
      
  formOffGenotype:  if stage=DuringMating, set this parameter to false
      will disallow random mating to set genotype.
      
  passOffspringOnly:  Default to false. If true, p
      yOperatorwill expect a function of form func(off, param),
      instead of func(pop, off, dad, mon, param) when
      passOffspringOnly is false. Since many duringMating pyOperatoro
      nly need access to offspring, this will imporve efficiency.
      
  Note: (FIXME) output to output or outputExpr is not yet supported.
  Ideally, this func will take two parameters with pop and then a
  filehandle to output, however, differentiating output, append etc is
  too troublesome right now.

";

%feature("docstring")  simuPOP::pyOperator::~pyOperator " 

destructor

Usage:
  x.~pyOperator()

";

%ignore simuPOP::pyOperator::pyOperator(const pyOperator &rhs);

%feature("docstring")  simuPOP::pyOperator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pyOperator::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::pyOperator::applyDuringMating " 

give pop, offspring, pop and mom.

Usage:
  x.applyDuringMating(&pop, offspring, *dad=NULL, *mom=NULL)

";

%feature("docstring")  simuPOP::pyOperator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::pyPenetrance "

penetrance using user supplied function

Details:
  Assign penetrance value by calling a user supplied function

";

%feature("docstring")  simuPOP::pyPenetrance::pyPenetrance " 

provide locus and penetrance for 11, 12, 13 (in the form of
dictionary)

Usage:
  pyPenetrance(loci, *func, exposePenetrance=False,
    stage=DuringMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  loci:  susceptibility loci. The genotype at these loci will be
      passed to func.
      
  func:  a Python function that accept genotypes at susceptibility
      loci and return penetrance value.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::pyPenetrance::~pyPenetrance " 

destructor

Usage:
  x.~pyPenetrance()

";

%ignore simuPOP::pyPenetrance::pyPenetrance(const pyPenetrance &rhs);

%feature("docstring")  simuPOP::pyPenetrance::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pyPenetrance::penet " 

currently assuming diploid

Usage:
  x.penet(*ind)

";

%feature("docstring")  simuPOP::pyPenetrance::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::pyQuanTrait "

quantitative trait using user supplied function

Details:
  Assign qtrait value by calling a user supplied function

";

%feature("docstring")  simuPOP::pyQuanTrait::pyQuanTrait " 

provide locus and qtrait for 11, 12, 13 (in the form of dictionary)

Usage:
  pyQuanTrait(loci, *func, stage=PostMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  loci:  susceptibility loci. The genotype at these loci will be
      passed to func.
      
  func:  a Python function that accept genotypes at susceptibility
      loci and return qtrait value.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::pyQuanTrait::~pyQuanTrait " 

Usage:
  x.~pyQuanTrait()

";

%ignore simuPOP::pyQuanTrait::pyQuanTrait(const pyQuanTrait &rhs);

%feature("docstring")  simuPOP::pyQuanTrait::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pyQuanTrait::qtrait " 

currently assuming diploid

Usage:
  x.qtrait(*ind)

";

%feature("docstring")  simuPOP::pyQuanTrait::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::pySample "

thrink population accroding to some outside value

";

%feature("docstring")  simuPOP::pySample::pySample " 

create a python sampler

Usage:
  pySample(*keep, keepAncestralPops, &name=\"sample\",
    &nameExpr=\"\", times=1, &saveAs=\"\", &saveAsExpr=\"\",
    &format=\"auto\", stage=PostMating, begin=0, end=-1, step=1,
    at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  keep:  a carray of the length of population. its values will be
      assigned to info.  and other parameters please see
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::pySample::~pySample " 

destructor

Usage:
  x.~pySample()

";

%ignore simuPOP::pySample::pySample(const pySample &rhs);

%feature("docstring")  simuPOP::pySample::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pySample::drawsample " 

Usage:
  x.drawsample(&pop)

";

%feature("docstring")  simuPOP::pySample::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::pySelector "

selection using user supplied function

Details:
  Assign fitness value by calling a user supplied function

";

%feature("docstring")  simuPOP::pySelector::pySelector " 

provide locus and fitness for 11, 12, 13 (in the form of dictionary)

Usage:
  pySelector(loci, *func, stage=PreMating, begin=0, end=-1, step=1,
    at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  loci:  susceptibility loci. The genotype at these loci will be
      passed to func.
      
  func:  a Python function that accept genotypes at susceptibility
      loci and return fitness value.
      
  output:  and other parameters please refer to
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::pySelector::~pySelector " 

destructor

Usage:
  x.~pySelector()

";

%ignore simuPOP::pySelector::pySelector(const pySelector &rhs);

%feature("docstring")  simuPOP::pySelector::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pySelector::indFitness " 

currently assuming diploid

Usage:
  x.indFitness(*ind)

";

%feature("docstring")  simuPOP::pySelector::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::pySubset "

thrink population accroding to some outside value

";

%feature("docstring")  simuPOP::pySubset::pySubset " 

create a directmigrator

Usage:
  pySubset(*keep=NULL, stage=PostMating, begin=0, end=-1, step=1,
    at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  keep:  a carray of the length of population. its values will be
      assigned to info.  and other parameters please see
      help(baseOperator.__init__)

";

%feature("docstring")  simuPOP::pySubset::~pySubset " 

destructor

Usage:
  x.~pySubset()

";

%ignore simuPOP::pySubset::pySubset(const pySubset &rhs);

%feature("docstring")  simuPOP::pySubset::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::pySubset::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::pySubset::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") PySwigClientData "

";

%feature("docstring") PySwigObject "

";

%feature("docstring") PySwigPacked "

";

%feature("docstring") simuPOP::PythonCoutBuf "

create a stream buf that print to python sys.stdout cout will be
redirected to this buf to really output to python console.

";

%feature("docstring")  simuPOP::PythonCoutBuf::PythonCoutBuf " 

Usage:
  pythonCoutBuf()

";

%feature("docstring") simuPOP::quanTrait "

quantitative trait

Details:
  Genetic quantitative trait is tricky to simulate. In simuPOP,
   I employee an ability (fitness) to mate approach. Namely, the
  probability that an individual will be chosen for mating is
  proportional to its fitness value. More specifically,
  PreMating selectors assign fitness values to each individual.
  
  Sexless mating (e.g. binomialSelection) : i
  ndividuals are chosen at probabilities that are proportional to
  their fitness values. More specifically, if there are N individuals
  with fitness values $f_i, i=1,...,N $,
  individual $i$will have probability $
   \\\\frac{f_i}{\\\\sum_{j=1}^N f_j} $to be chosen to be passed to
  the next generation.
  
  Random mating with sex (e.g. randommating): males and females are
  separated and each are chosen as described above.
  
  Please refer to the user's guide for details.

";

%feature("docstring")  simuPOP::quanTrait::quanTrait " 

constructor. default to be always active.

Usage:
  quanTrait(stage=PostMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::quanTrait::~quanTrait " 

destructor

Usage:
  x.~quanTrait()

";

%feature("docstring")  simuPOP::quanTrait::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::quanTrait::qtrait " 

calculate/return quantitative trait etc

Usage:
  x.qtrait(*ind)

";

%feature("docstring")  simuPOP::quanTrait::apply " 

set qtrait to all individual

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::quanTrait::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::randomMating "

Details:
  basic sexual random mating.
  Within each subpopulation, choose male and female randomly randmly
  get one copy of chromosome from father/mother.
  require: sexed individual; ploidy == 2
  apply during mating operators and put into the next generation.
  if ignoreParentsSex is set, parents will be chosen regardless of
  sex.
  Otherwise, male and female will be collected and be chosen randomly.
  If there is no male or female in a subpopulation, if
  m_UseSameSexIfUniSex is true, an warning will be generated and same
  sex mating (?) will be used otherwise, randomMatingw
  ill return false.
  if there is no during mating operator to copy alleles, a direct copy
  will be used.

";

%feature("docstring")  simuPOP::randomMating::randomMating " 

create a random mating scheme

Usage:
  randomMating(numOffspring=1., *numOffspringFunc=NULL,
    maxNumOffspring=0, mode=MATE_NumOffspring, newSubPopSize=[],
    *newSubPopSizeFunc=NULL, newSubPopSizeExpr=\"\",
    contWhenUniSex=True)

Arguments:

  numOffspring:  
  number:  of offspring or p in some modes
      
  numOffspringFunc:  
  a:  python function that determine number of offspring or p
      depending on mode
      
  maxNumOffspring:  used when mode=MATE_BinomialDistribution
      
  mode:  one of MATE_NumOffspring , MATE_NumOffspringEachFamily,
      MATE_GeometricDistribution, MATE_PoissonDistribution,
      MATE_BinomialDistribution
      
  newSubPopSize:  an array of subpopulation sizes, should have the
      same number of subpopulations as current population
      
  newSubPopSizeExpr:  an expression that will be evaluated as an array
      of subpop sizes
      
  newSubPopSizeFunc:  an function that have parameter gen and oldSize
      (current subpop size).
      
  contWhenUniSex:  continue when there is only one sex in the
      population, default to true

";

%feature("docstring")  simuPOP::randomMating::~randomMating " 

destructor

Usage:
  x.~randomMating()

";

%feature("docstring")  simuPOP::randomMating::clone " 

 clone() const. Generate a copy of itself and return pointer this
is to make sure the object is persistent and will not be freed by
python.

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::randomMating::isCompatible " 

check if the mating type is compatible with population structure

Usage:
  x.isCompatible(&pop)

Details:
  possible things to check:need certain types of individual (age,
  sex etc)
  need resizeable population...

";

%feature("docstring")  simuPOP::randomMating::__repr__ " 

return name of the mating type

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::randomMating::submitScratch " 

Usage:
  x.submitScratch(&pop, &scratch)

";

%feature("docstring")  simuPOP::randomMating::mate " 

do the mating. parameters see mating::mate.

Usage:
  x.mate(&pop, &scratch, &ops, submit)

Details:
  Within each subpopulation, choose male and female randomly randmly
  get one copy of chromosome from father/mother.
  require: sexed individual; ploidy == 2
  apply during mating operators and put into the next generation.
  Otherwise, male and female will be collected and be chosen randomly.
  If there is no male or female in a subpopulation,
  if m_contWhenUniSex is true, an warning will be generated and same
  sex mating (?) will be used
  otherwise, randomMatingwill return false.
  
  determine if any during-mating operator will generate offspring
  genotype
  random mating happens within each subpopulation
  now, all individuals of needToFind sex is collected
  if selection is on
  apply all during mating operators

";

%feature("docstring") simuPOP::randomSample "

thrink population accroding to some outside value

";

%feature("docstring")  simuPOP::randomSample::randomSample " 

draw random sample, regardless of affected status

Usage:
  randomSample(size=[], &name=\"sample\", &nameExpr=\"\", times=1,
    &saveAs=\"\", &saveAsExpr=\"\", &format=\"auto\",
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  size:  size of sample. It can be either a number, representing the
      overall sample size, regardless of population strucutre; or an
      array, representing number of samples drawn from each
      subpopulation.
      
  stage:  and other parameters please see help(baseOperator.__init__)
      
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

Note:
  ancestral populations will not be copied to the samples

";

%feature("docstring")  simuPOP::randomSample::~randomSample " 

destructor

Usage:
  x.~randomSample()

";

%feature("docstring")  simuPOP::randomSample::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::randomSample::preparesample " 

value checking

Usage:
  x.preparesample(&pop)

";

%feature("docstring")  simuPOP::randomSample::drawsample " 

Usage:
  x.drawsample(&pop)

";

%feature("docstring")  simuPOP::randomSample::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::recombinator "

Recombination.

Details:
  only works for diploids (and for females in haplodiploids)
  population.
  
  Free recombination between loci. Loci behave completely
  independently.
  
  otherwise there will be some linkage between loci, user need to
  specify physical recombination rate between adjacent loci (ie
  between locus n and n+1)
  
  The recombination rate must be comprised between 0.0 and 0.5.
  
  A recombination rate of 0.0 means that the loci are completely
  linked and thus behave together as a single linked locus.
  
  A recombination rate of 0.5 is equivalent to free recombination.
  
  All values in between will represent various linkage intensities
  between adjacent pairs of loci. The recombination rate is equivalent
  to 1-linkage and represents the probability that the allele at the
  next locus is randomly drawn.

";

%feature("docstring")  simuPOP::recombinator::recombinator " 

recombine chromosomes from parents

Usage:
  recombinator(intensity=-1, rate=[], afterLoci=[],
    maleIntensity=-1, maleRate=[], maleAfterLoci=[], begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  intensity:  recombination rate per unit of loci distance. I.e., the
      really recombination rate between two loci is determined by
      intensity*loci distance between them.
      
  rate:  recombination rate regardless of loci distance; it can also
      be an array of recombination rates. Must be the same length as
      afterLoci or totNumOfLoci(). If totNumLociThe last item can be
      ignored.
      
  afterLoci:  an array of loci index. If rates is also specified, they
      should have the same length. Default to all loci (but
      meaningless for those loci locate at the end of chromosome.) If
      given, afterLoci should be ordered, and can not include loci at
      the end of a chromosome.  recombination rate for male
      individuals. If given, parameter rate will be considered as
      female rate.  recombination intensity for male individuals. If
      given, parameter intensity will be considered as female
      intensity.  if given, male will recombine at different
      locations. This is rarely used.

Note:
  rate gives recombination rate PER unit, rates gives plain rates.
  there is no recombination between sex chromosomes of male
  individuals. (sexChrom=True).

";

%feature("docstring")  simuPOP::recombinator::~recombinator " 

Usage:
  x.~recombinator()

";

%feature("docstring")  simuPOP::recombinator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::recombinator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::recombinator::prepareRecRates " 

this function takes intensity, rate, afterLoci, ... inputs and
return a bernulli trailer and a recBeforeLoci vector.

Usage:
  x.prepareRecRates(&pop, intensity, rate, afterLoci, sexChrom,
    &recBeforeLoci, &vecP)

Details:
  get loci distance * rate and then recombinant points
  get loci distance * rate and then recombinant points
  initialize recombination counter, This will count recombination
  events after each locus.

";

%feature("docstring")  simuPOP::recombinator::recCount " 

return recombination count

Usage:
  x.recCount(locus)

";

%feature("docstring")  simuPOP::recombinator::recCounts " 

return recombination counts

Usage:
  x.recCounts()

";

%feature("docstring")  simuPOP::recombinator::recombine " 

Usage:
  x.recombine(*parent, offspring, offPloidy, &bt, &recBeforeLoci,
    setSex=False)

";

%feature("docstring")  simuPOP::recombinator::applyDuringMating " 

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

%feature("docstring")  simuPOP::RNG::max " 

Usage:
  x.max()

";

%feature("docstring")  simuPOP::RNG::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::RNG::randGet " 

return min, 1, 2, ... max

Usage:
  x.randGet()

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

return geometric with parameter p (k>=1)

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

%feature("docstring")  simuPOP::RNG::randBinomial " 

binomial distribution

Usage:
  x.randBinomial(n, p)

";

%feature("docstring")  simuPOP::RNG::randPoisson " 

Poisson distribution.

Usage:
  x.randPoisson(p)

";

%feature("docstring") simuPOP::sample "

sample from population and save samples sample operator will
generate a new subpopulation in pop namespace.

";

%feature("docstring")  simuPOP::sample::sample " 

create a sample

Usage:
  sample(&name=\"sample\", &nameExpr=\"\", times=1, &saveAs=\"\",
    &saveAsExpr=\"\", &format=\"auto\", stage=PostMating, begin=0,
    end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

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

%feature("docstring")  simuPOP::sample::~sample " 

destructor

Usage:
  x.~sample()

";

%feature("docstring")  simuPOP::sample::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::sample::preparesample " 

Usage:
  x.preparesample(&)

";

%feature("docstring")  simuPOP::sample::drawsample " 

Usage:
  x.drawsample(&pop)

";

%feature("docstring")  simuPOP::sample::samples " 

return the samples

Usage:
  x.samples(&pop)

";

%feature("docstring")  simuPOP::sample::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::sample::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::savePopulation "

save population to a file

";

%feature("docstring")  simuPOP::savePopulation::savePopulation " 

Usage:
  savePopulation(output=\"\", outputExpr=\"\", format=\"bin\",
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

";

%feature("docstring")  simuPOP::savePopulation::~savePopulation " 

Usage:
  x.~savePopulation()

";

%feature("docstring")  simuPOP::savePopulation::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::savePopulation::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::savePopulation::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::selector "

selection

Details:
  Genetic selection is tricky to simulate. In simuPOP,
   I employee an ability (fitness) to mate approach. Namely, the
  probability that an individual will be chosen for mating is
  proportional to its fitness value. More specifically,
  PreMating selectors assign fitness values to each individual.
  
  Sexless mating (e.g. binomialSelection) : i
  ndividuals are chosen at probabilities that are proportional to
  their fitness values. More specifically, if there are N individuals
  with fitness values $f_i, i=1,...,N $,
  individual $i$will have probability $
   \\\\frac{f_i}{\\\\sum_{j=1}^N f_j} $to be chosen to be passed to
  the next generation.
  
  Random mating with sex (e.g. randommating): males and females are
  separated and each are chosen as described above.
  
  Please refer to the user's guide for details.

";

%feature("docstring")  simuPOP::selector::selector " 

constructor. default to be always active.

Usage:
  selector(stage=PreMating, begin=0, end=-1, step=1, at=[],
    rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::selector::~selector " 

destructor

Usage:
  x.~selector()

";

%feature("docstring")  simuPOP::selector::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::selector::indFitness " 

calculate/return w11 etc

Usage:
  x.indFitness(*ind)

";

%feature("docstring")  simuPOP::selector::apply " 

set fitness to all individual

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::selector::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::setAncestralDepth "

";

%feature("docstring")  simuPOP::setAncestralDepth::setAncestralDepth " 

timer if called, output time passed since last calling time.

Usage:
  setAncestralDepth(depth, output=\">\", outputExpr=\"\",
    stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

";

%feature("docstring")  simuPOP::setAncestralDepth::~setAncestralDepth " 

destructor

Usage:
  x.~setAncestralDepth()

";

%feature("docstring")  simuPOP::setAncestralDepth::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::setAncestralDepth::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::setAncestralDepth::__repr__ " 

Usage:
  x.__repr__()

";

%ignore simuPOP::SharedVariables;

%ignore simuPOP::SharedVariables::SharedVariables();

%ignore simuPOP::SharedVariables::SharedVariables(PyObject *dict, bool ownVars);

%ignore simuPOP::SharedVariables::SharedVariables(const SharedVariables &rhs);

%ignore simuPOP::SharedVariables::swap(SharedVariables &rhs);

%ignore simuPOP::SharedVariables::~SharedVariables();

%feature("docstring")  simuPOP::SharedVariables::clear " 

Usage:
  x.clear()

";

%feature("docstring")  simuPOP::SharedVariables::setDict " 

Usage:
  x.setDict(*dict)

";

%feature("docstring")  simuPOP::SharedVariables::checkRefCount " 

check the ref count of this shared variable recursively. report
anything that has refcount > 1

Usage:
  x.checkRefCount()

";

%feature("docstring")  simuPOP::SharedVariables::setVar " 

setvars C++ ==> Python

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

%feature("docstring")  simuPOP::SharedVariables::dict " 

Usage:
  x.dict()

";

%ignore simuPOP::SharedVariables::asString() ;

%feature("docstring")  simuPOP::SharedVariables::fromString " 

Usage:
  x.fromString(&vars)

";

%feature("docstring") simuPOP::simulator "

simulator manage several replicates of a population, evolve them
using given mating scheme and operators.

Details:
  simulators combine three important components of simuPOP:
   population, mating scheme and operators together. A simulator is
  created with an instance of population, a replicate number and a
  mating scheme. It makes 'rep' replicates of this population and
  control the evolution process of these populations.
  The most important functions of a simulator is of course evolve().
   It accepts arrays of operators as its parameters, among which,
  'preOps' and 'postOps' will be applied to the populations at the
  begining and end of evolution, whereas 'operators'will be applied at
  every generation.
  simulators separates operators into pre-, during- and post- mating
  operators. During evolution, simulator first apply all pre-mating
  operators and then call the mate() function of the given mating
  scheme, which will call during-mating operators during the birth of
  each offsrping. After the mating is finished, post-mating operators
  are applied in the order they apprear in the operator list.
  Since operators can apply to specific replicate or replicates group,
  and might not be active at all time, the isActive(m_curRep,
  m_numRep, m_gen, end, grp()) function of each
  operator is called before it is applied to the populations.
  simulators can evolve a given number of generations (the 'end'
  parameter of evolve), or evolve indefinitely using a certain type of
  operators called terminators. In this case, one or more terminators
  will check the status of evolution and determine if the simulation
  should be stopped. An obvious example of such a terminator is a
  fixation-checker.
  Finally, a simulator can be saved to a file in the format of 'text',
  'bin', or 'xml'. This enables us to stop a simulation and ressume it
  at another time or on another machine. It is also a good idea to
  save a snapshot of a simulation every several generations.

";

%feature("docstring")  simuPOP::simulator::simulator " 

create a simulator

Usage:
  simulator(&pop, &matingScheme, varName=\"simuVars\", rep=1,
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

%feature("docstring")  simuPOP::simulator::~simulator " 

destroy a simulator along with all its populations

Usage:
  x.~simulator()

Note:
  pop = simulator::population() returns temporary reference to an
  internal population. After a simulator evolves another genertion or
  after the simulator is destroyed, this referenced population should
  *not* be used.

";

%feature("docstring")  simuPOP::simulator::pop " 

the 'rep' replicate of this simulator

Usage:
  x.pop(rep)

Arguments:

  rep:  number of replicate.

Value:
  reference to population rep.

Note:
  The returned reference is temporary in the sense that the refered
  population will be invalid after another round of evolution.
  Therefore, the use of this function should be limited to 'immediate
  after retrival'. If you would like to get a persistent population,
  use getPopulation(rep)

";

%feature("docstring")  simuPOP::simulator::getPopulation " 

Usage:
  x.getPopulation(rep)

Arguments:

  rep:  number of replicate.

Details:
  this function returns a copy of population rep

Value:
  reference to a population.

";

%feature("docstring")  simuPOP::simulator::setPopulation " 

Usage:
  x.setPopulation(&pop, rep)

";

%ignore simuPOP::simulator::curRep() ;

%feature("docstring")  simuPOP::simulator::numRep " 

Usage:
  x.numRep()

";

%feature("docstring")  simuPOP::simulator::gen " 

Usage:
  x.gen()

";

%ignore simuPOP::simulator::grp();

%feature("docstring")  simuPOP::simulator::group " 

return group indices

Usage:
  x.group()

";

%feature("docstring")  simuPOP::simulator::setGroup " 

set groups for replicates

Usage:
  x.setGroup(&grp)

";

%feature("docstring")  simuPOP::simulator::setGen " 

set generation number

Usage:
  x.setGen(gen)

Arguments:

  gen:  new generation number

Note:
  this will also set shared variable gen

";

%feature("docstring")  simuPOP::simulator::step " 

evolve one step

Usage:
  x.step(&ops=[], &preOps=[], &postOps=[], steps=1)

See Also:
   simulator::evolve()

";

%feature("docstring")  simuPOP::simulator::evolve " 

evolve till 'end' generation subject to given operators

Usage:
  x.evolve(&ops, &preOps=[], &postOps=[], end=-1, dryrun=False)

Arguments:

  ops:  operators that will be applied at all generations. Of course
      they might not be active at all generations.
      
  preOps:  operators that will be applied before evolution. e
      volve()function will *not* check if they are active.
      
  postOps:  operators that will be applied after the last generation,
      or after a replicate is terminated.
      
  end:  ending generation. Should be -1 (no ending generation) or a
      number greater than current generation number. When end=-1,
      simulator can only be stopped by terminators.

Details:
  
  'operators' will be applied in the order of:
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

%feature("docstring")  simuPOP::simulator::apply " 

apply some ops, geneartion of population does not change No mating
is allowed.

Usage:
  x.apply(ops, dryrun=False)

Arguments:

  ops:  operators that will be applied at all generations. Of course
      they might not be active at all generations.

Details:
  
  pre-mating oeprators are applied before post-mating operators. no
  during-mating operators are allowed.

Value:
  True if evolution finishs successfully.

";

%feature("docstring")  simuPOP::simulator::setStopIfOneRepStop " 

stop if one replicate stops or not

Usage:
  x.setStopIfOneRepStop(on=True)

Arguments:

  on:  turn on or off stopIfOneRepStop
      
  if set, the simulator will stop evolution if one replicate stops.
  This is sometimes useful.

";

%feature("docstring")  simuPOP::simulator::stopIfOneRepStop " 

Usage:
  x.stopIfOneRepStop()

";

%feature("docstring")  simuPOP::simulator::setApplyOpToStoppedReps " 

apply ops even if rep stops

Usage:
  x.setApplyOpToStoppedReps(on=True)

Arguments:

  on:  turn on or off applyOpToStoppedReps flag
      
  if set, the simulator will continue to apply operators to all
  stopped repicates until all replicates are marked stopped. This is
  sometimes useful.

";

%feature("docstring")  simuPOP::simulator::applyOpToStoppedReps " 

Usage:
  x.applyOpToStoppedReps()

";

%feature("docstring")  simuPOP::simulator::vars " 

get simulator namespace, if rep > 0 is given, return replicate rep
namespace

Usage:
  x.vars(rep, subPop=-1)

";

%feature("docstring")  simuPOP::simulator::saveSimulator " 

save simulator in 'text','bin' or 'xml' format

Usage:
  x.saveSimulator(filename, format=\"auto\")

Arguments:

  filename:  save to filename
      
  format:  format to save. Can be one of 'text', 'bin', 'xml'
      
  The default format is 'text' but the output is not suppored to be
  read. 'bin' has smaller size and should be used for large
  populations. 'xml' format is most readable and should be used when
  you would like to convert simuPOPpopulatio
  ns to other formats.

See Also:
  global function loadsimulator

";

%ignore simuPOP::simulator::loadSimulator(string filename, string format="auto");

%feature("docstring")  simuPOP::simulator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::smmMutator "

stepwise mutation model.

Details:
  Stepwise mutation model (SMM) assumes that alleles are represented
  by integer values and that a mutation either increases or decreases
  the allele value by one. For variable number tandem repeats loci
  (VNTR), the allele value is generally taken as the number of tandem
  repeats in the DNA sequence.

See Also:
  Kimura & Ohta 1978

";

%feature("docstring")  simuPOP::smmMutator::smmMutator " 

Usage:
  smmMutator(rate=[], atLoci=vectori, maxAllele=0, incProb=0.5,
    output=\">\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  rate::  mutation rate
      
  incProb:  probability to increase allele state. Default to 1
      
  atLoci:  and other parameters: refer to help(mutator),
      help(baseOperator.__init__)

Details:
  The stepwise mutation model (SMM) is developed for allozymes. It
  provides better description for these kinds of evolutionary
  processes.

";

%feature("docstring")  simuPOP::smmMutator::~smmMutator " 

Usage:
  x.~smmMutator()

";

%feature("docstring")  simuPOP::smmMutator::mutate " 

how to mutate a single allele. this is usually the only function
that need to be defined by the subclasses.

Usage:
  x.mutate(allele)

";

%feature("docstring")  simuPOP::smmMutator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::smmMutator::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::splitSubPop "

split subpopulation

";

%feature("docstring")  simuPOP::splitSubPop::splitSubPop " 

split a subpopulation (or whole population as subpop 0)

Usage:
  splitSubPop(which=0, sizes=[], proportions=[], subPopID=[],
    stage=PreMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

Arguments:

  which:  which subpop to split (if there is no subpop structure, 0 is
      the only subpop)
      
  subPop:  new subpop size, should add up to the size of subpop to be
      splitted
      
  proportions:  proportions of new subpop. (use one of subPop or
      proportions). Should add up to one.
      
  subPopID:  optional. new subpop IDs. If given, should have the same
      length as subPop or proportions. SInce subpop with negative id
      will be removed. You can remove part of a subpop by setting a
      new negative id.

";

%feature("docstring")  simuPOP::splitSubPop::~splitSubPop " 

destructor

Usage:
  x.~splitSubPop()

";

%feature("docstring")  simuPOP::splitSubPop::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::splitSubPop::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::splitSubPop::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::spread "

initialize genotype by value and then copy to all individuals

";

%feature("docstring")  simuPOP::spread::spread " 

Usage:
  spread(ind, subPop=[], stage=PreMating, begin=0, end=1, step=1,
    at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::spread::~spread " 

Usage:
  x.~spread()

";

%feature("docstring")  simuPOP::spread::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::spread::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::spread::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring") simuPOP::stat "

";

%feature("docstring")  simuPOP::stat::stat " 

create an stat

Usage:
  stat(popSize=False, numOfMale=False, numOfAffected=False,
    numOfAlleles=vectori, alleleFreq=vectori, heteroFreq=vectori,
    expHetero=vectori, homoFreq=vectori, genoFreq=vectori,
    haploFreq=intMatrix, LD=intMatrix, Fst=vectori,
    relGroups=intMatrix, relLoci=vectori, relBySubPop=False,
    relMethod=vectori, relMinScored=10, hasPhase=False,
    midValues=False, output=\"\", outputExpr=\"\", stage=PostMating,
    begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

Arguments:

  popSize:  whether or not calculate population sizes. will set
      numSubPop, subPopSize, popSize, subPop[sp]['popSize']
      
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
      
  homoFreq:  an array of loci at which homozygosity number and
      frequcies will be calculated
      
  expHetero:  an array of loci at which expected heterozygosity will
      be calculated.
      
  haploFreq:  a matrix of haplotypes (allele sequence) to count
      
  format: haploFreq = [ [ 0,1,2 ], [1,2] ]
  All haplotypes on loci 012, 12 will be counted. If only one
  haplotype is specified, the outer [] can be ommited. I.e.,
  haploFreq=[0,1] is acceptable.

  LD:  apply LD, LD' and r2, given LD=[ [locus1 locus2], [ locus1
      locus2 allele1 allele2], ,...] If two numbers are given, D, D'
      and r2 overall possible allele pairs will be calculated and
      saved as AvgLD etc. If four numbers are given, D, D' and r2
      using specified alleles are provided. If only one item is
      specified, the outer [] can be ignored. I.e., LD=[locus1 locus2]
      is acceptable.
      
  Fst:  calculate Fst. Fis and Fit will be given as a side product.
      
  format Fst = [ 0, 1 ] Fst calculated at locus 0 and 1. Since any
  allele can be used to calculate Fst, Fst[0] will be the average over
  all alleles as suggested by Weir and Cockerham.

  relGroups:  calculated pairwise relatedness between groups.
      relGroups can be in the form of either [[1,2,3],[4,5],[7,8]]
      (gourps of individuals) or [1,3,4] (use subpopulations).
      
  relMethod:  method used to calculate relatedness. Can be either
      REL_Queller or REL_Lynch.
      
  relLoci:  loci on which relatedness calues are calculated.
      
  hasPhase:  if a/b and b/a are the same genotype. default is false.
      
  midValues:  whether or not post intermediate results. Default to
      false. For example, Fst will need to calculate allele frequency,
      if midValues is set to true, allele frequencies will be posted
      as well. This will help debuging and sometimes derived
      statistics.
      
  others:  there is NO output for this operator. other parameters
      please see help(baseOperator.__init__)

Note:
  previous provide 'search for all allele/genotype' option but this
  has proven to be troublesome. In this version, everything should be
  explicitly specified.

";

%feature("docstring")  simuPOP::stat::~stat " 

Usage:
  x.~stat()

";

%feature("docstring")  simuPOP::stat::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::stat::apply " 

count various statistics. use m_alleles etc to save (potentially)
time to resize all these variables.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::stat::__repr__ " 

Usage:
  x.__repr__()

";

%ignore simuPOP::statAlleleFreq;

%feature("docstring")  simuPOP::statAlleleFreq::statAlleleFreq " 

Usage:
  statAlleleFreq(&atLoci=vectori)

";

%feature("docstring")  simuPOP::statAlleleFreq::addLocus " 

Usage:
  x.addLocus(locus, post=False)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleNumAll " 

Usage:
  x.alleleNumAll()

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleNumVec " 

Usage:
  x.alleleNumVec(loc)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleNum " 

Usage:
  x.alleleNum(allele, loc)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleFreqAll " 

Usage:
  x.alleleFreqAll()

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleFreqVec " 

Usage:
  x.alleleFreqVec(loc)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleFreq " 

Usage:
  x.alleleFreq(allele, loc)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleNumAll " 

Usage:
  x.alleleNumAll(subPop)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleNumVec " 

Usage:
  x.alleleNumVec(loc, subPop)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleNum " 

Usage:
  x.alleleNum(allele, loc, subPop)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleFreqAll " 

Usage:
  x.alleleFreqAll(subPop)

";

%feature("docstring")  simuPOP::statAlleleFreq::alleleFreqVec " 

Usage:
  x.alleleFreqVec(loc, subPop)

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

%ignore simuPOP::statExpHetero;

%feature("docstring")  simuPOP::statExpHetero::statExpHetero " 

Usage:
  statExpHetero(&alleleFreq, &expHetero=vectori)

";

%feature("docstring")  simuPOP::statExpHetero::apply " 

Usage:
  x.apply(&pop)

";

%ignore simuPOP::statFst;

%feature("docstring")  simuPOP::statFst::statFst " 

Usage:
  statFst(&alleleFreq, &heteroFreq, &Fst=vectori, midValues=False)

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
  statGenoFreq(&genoFreq=vectori, phase=True)

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
  x.addHaplotype(&haplo, post=False)

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
  statHeteroFreq(&heteroFreq=vectori, &homoFreq=vectori)

";

%feature("docstring")  simuPOP::statHeteroFreq::addLocus " 

Usage:
  x.addLocus(locus, post=False)

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
  statLD(&alleleFreq, &haploFreq, &LD=intMatrix, midValues=False)

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

%feature("docstring") simuPOP::stator "

NOTE: the default output for stator is \"\", i.e., no output i.e.,
stator will write to shared variables and unless specified by
output=\">\" etc, no output will be generated. this class will also
list ALL statistics and its names?

";

%feature("docstring")  simuPOP::stator::stator " 

constructor. default to be always active. default to have NO output
(shared variables will be set.) phase: if we treat Aa!=aA, default
is false, i.e, Aa=aA.

Usage:
  stator(output=\"\", outputExpr=\"\", stage=PostMating, begin=0,
    end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::stator::~stator " 

destructor

Usage:
  x.~stator()

";

%feature("docstring")  simuPOP::stator::clone " 

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

%feature("docstring") simuPOP::statRelatedness "

the relatedness measure between two individuals/families using
Queller and Goodnight or Lynch's method. or internal relatedness
values

";

%feature("docstring")  simuPOP::statRelatedness::statRelatedness " 

calculate relatedness measures between elements in groups

Usage:
  statRelatedness(&alleleFreq, &groups=intMatrix, useSubPop=False,
    &loci=vectori, method=vectori, minScored=10, midValues=False)

Arguments:

  groups:  can be [ [1,2,3],[4,5,6],[7,8,9]] as three groups of
      individuals; or [ 1 3 4] as three subpopulations. To specify
      between individual relatedness, use [[1],[2],[3]] (the first
      form). If this parameter is ignored, this operator calculate
      relatedness between all subpopulations.
      
  method:  can be REL_Queller, REL_Lynch, REL_IR, REL_D2 or REL_Rel.
      Please refer to the manual for details.

";

%feature("docstring")  simuPOP::statRelatedness::relQueller " 

Usage:
  x.relQueller(ind1, ind2)

";

%feature("docstring")  simuPOP::statRelatedness::relLynch " 

Usage:
  x.relLynch(ind1, ind2)

";

%feature("docstring")  simuPOP::statRelatedness::relIR " 

Usage:
  x.relIR(ind1, locus)

";

%feature("docstring")  simuPOP::statRelatedness::relD2 " 

Usage:
  x.relD2(ind1, locus)

";

%feature("docstring")  simuPOP::statRelatedness::relRel " 

Usage:
  x.relRel(ind1, ind2, locus)

";

%feature("docstring")  simuPOP::statRelatedness::groupRelatedness " 

for group i and locus j otherwise

Usage:
  x.groupRelatedness(&pop, i, j, method)

";

%feature("docstring")  simuPOP::statRelatedness::apply " 

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

%feature("docstring") swig_cast_info "

";

%feature("docstring") swig_const_info "

";

%feature("docstring") swig_module_info "

";

%feature("docstring") swig_type_info "

";

%feature("docstring") simuPOP::SystemError "

exception, thrown if system error occurs

";

%feature("docstring")  simuPOP::SystemError::SystemError " 

Usage:
  systemError(msg)

";

%feature("docstring") simuPOP::tagger "

tagger is a during mating operator that tag individual with various
information. Potential usages are 1. record parenting information to
track pedigree. 2. tag a individual/allele and monitor its spread in
the population etc. 3...

Details:
  Bo Peng

";

%feature("docstring")  simuPOP::tagger::tagger " 

constructor. default to be always active but no output.

Usage:
  tagger(begin=0, end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::tagger::~tagger " 

destructor

Usage:
  x.~tagger()

";

%feature("docstring")  simuPOP::tagger::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring") simuPOP::terminateIf "

terminate according to a condition which can be, e.g.
any(alleleNum0) == 0 all(alleleNum1) > 0.5 alleleNum0{2} == 0 etc.
When the condition is true, a shared variable var=\"terminate\" will
be set to current generation.

";

%feature("docstring")  simuPOP::terminateIf::terminateIf " 

Usage:
  terminateIf(condition=\"\", message=\"\", var=\"terminate\",
    output=\"\", outputExpr=\"\", stage=PostMating, begin=0, end=-1,
    step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::terminateIf::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::terminateIf::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring")  simuPOP::terminateIf::apply " 

check all alleles in vector allele if they are fixed.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::terminateIf::~terminateIf " 

Usage:
  x.~terminateIf()

";

%feature("docstring") simuPOP::terminator "

";

%feature("docstring")  simuPOP::terminator::terminator " 

constructor. default to be always active.

Usage:
  terminator(message=\"\", output=\">\", outputExpr=\"\",
    stage=PostMating, begin=0, end=-1, step=1, at=[], rep=REP_ALL,
    grp=GRP_ALL)

";

%feature("docstring")  simuPOP::terminator::~terminator " 

destructor

Usage:
  x.~terminator()

";

%feature("docstring")  simuPOP::terminator::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::terminator::message " 

Usage:
  x.message()

";

%feature("docstring") simuPOP::ticToc "

";

%feature("docstring")  simuPOP::ticToc::ticToc " 

timer if called, output time passed since last calling time.

Usage:
  ticToc(output=\">\", outputExpr=\"\", stage=PreMating, begin=0,
    end=-1, step=1, at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::ticToc::~ticToc " 

destructor

Usage:
  x.~ticToc()

";

%feature("docstring")  simuPOP::ticToc::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::ticToc::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::ticToc::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::trajectory "

";

%feature("docstring")  simuPOP::trajectory::trajectory " 

Usage:
  trajectory(nTraj=1)

";

%feature("docstring")  simuPOP::trajectory::numTraj " 

Usage:
  x.numTraj()

";

%feature("docstring")  simuPOP::trajectory::maxLen " 

Usage:
  x.maxLen()

";

%feature("docstring")  simuPOP::trajectory::setTraj " 

Usage:
  x.setTraj(&freq, idx=0)

";

%feature("docstring")  simuPOP::trajectory::traj " 

Usage:
  x.traj(idx=0)

";

%feature("docstring") simuPOP::turnOffDebugOp "

";

%feature("docstring")  simuPOP::turnOffDebugOp::turnOffDebugOp " 

turn on debug

Usage:
  turnOffDebugOp(code, stage=PreMating, begin=0, end=-1, step=1,
    at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::turnOffDebugOp::~turnOffDebugOp " 

destructor

Usage:
  x.~turnOffDebugOp()

";

%feature("docstring")  simuPOP::turnOffDebugOp::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::turnOffDebugOp::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::turnOffDebugOp::__repr__ " 

Usage:
  x.__repr__()

";

%feature("docstring") simuPOP::turnOnDebugOp "

";

%feature("docstring")  simuPOP::turnOnDebugOp::turnOnDebugOp " 

turn on debug

Usage:
  turnOnDebugOp(code, stage=PreMating, begin=0, end=-1, step=1,
    at=[], rep=REP_ALL, grp=GRP_ALL)

";

%feature("docstring")  simuPOP::turnOnDebugOp::~turnOnDebugOp " 

destructor

Usage:
  x.~turnOnDebugOp()

";

%feature("docstring")  simuPOP::turnOnDebugOp::clone " 

this function is very important

Usage:
  x.clone()

";

%feature("docstring")  simuPOP::turnOnDebugOp::apply " 

apply to one population, does not check if the oeprator is
activated.

Usage:
  x.apply(&pop)

";

%feature("docstring")  simuPOP::turnOnDebugOp::__repr__ " 

Usage:
  x.__repr__()

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

%feature("docstring") simuPOP::Weightedsampler "

";

%feature("docstring")  simuPOP::Weightedsampler::Weightedsampler " 

Usage:
  weightedsampler(&rng, &weight=[], fast=True)

";

%feature("docstring")  simuPOP::Weightedsampler::~Weightedsampler " 

Usage:
  x.~Weightedsampler()

";

%feature("docstring")  simuPOP::Weightedsampler::set " 

FIXME: consider adopting R's implementation. They may be quicker.

Usage:
  x.set(&weight)

";

%feature("docstring")  simuPOP::Weightedsampler::biSearch " 

Usage:
  x.biSearch(a)

";

%feature("docstring")  simuPOP::Weightedsampler::get " 

Usage:
  x.get()

";

%feature("docstring")  simuPOP::Weightedsampler::get " 

Usage:
  x.get(&res, shift=0)

";

%feature("docstring")  simuPOP::Weightedsampler::get " 

Usage:
  x.get(beg, end, shift=0)

";

%feature("docstring")  simuPOP::Weightedsampler::q " 

Usage:
  x.q()

";

%feature("docstring")  simuPOP::Weightedsampler::a " 

Usage:
  x.a()

";

%feature("docstring")  simuPOP::FreqTrajectoryStoch " 

Usage:
  FreqTrajectoryStoch(freq, N=0, *NtFunc=NULL, s=[], *sFunc=NULL,
    T=100000)

";

%feature("docstring")  simuPOP::FreqTrajectoryMultiStoch " 

Usage:
  FreqTrajectoryMultiStoch(freq=[], N=0, *NtFunc=NULL, s=[],
    *sFunc=NULL, T=100000)

";

%feature("docstring")  simuPOP::FreqTrajectorySelSim " 

Usage:
  FreqTrajectorySelSim(sel, Ne, freq, dom_h, selection)

";

%feature("docstring")  simuPOP::FreqTrajectoryForward " 

Usage:
  FreqTrajectoryForward(lowbound, highbound, disAge, grate, N0,
    seleCo)

";

%feature("docstring")  simuPOP::LoadPopulation " 

Usage:
  LoadPopulation(&file, &format=\"auto\")

";

%feature("docstring")  simuPOP::LoadSimulator " 

Usage:
  LoadSimulator(&file, &mate, format=\"auto\")

";

%ignore simuPOP::haploKey(const vectori &seq);

%feature("docstring")  simuPOP::TurnOnDebug " 

set debug code, default to turn all code on

Usage:
  TurnOnDebug(code)

";

%feature("docstring")  simuPOP::TurnOffDebug " 

turn off debug, default to turn all code off

Usage:
  TurnOffDebug(code)

";

%feature("docstring")  simuPOP::debug " 

test if one code is turned on

Usage:
  debug(code)

";

%feature("docstring")  simuPOP::ListDebugCode " 

show all dbg codes (print to cout)

Usage:
  ListDebugCode()

";

%feature("docstring")  simuPOP::dbgString " 

dbg string for a code

Usage:
  dbgString(code)

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

%ignore simuPOP::PyObj_Is_AlleleNumArray(PyObject *obj);

%ignore simuPOP::Int_Vec_As_NumArray(vectori::iterator begin, vectori::iterator end);

%ignore simuPOP::Double_Vec_As_NumArray(vectorf::iterator begin, vectorf::iterator end);

%ignore simuPOP::Allele_Vec_As_NumArray(GenoIterator begin, GenoIterator end);

%ignore simuPOP::NumArray_Size(PyObject *obj);

%ignore simuPOP::NumArray_Data(PyObject *obj);

%feature("docstring")  simuPOP::save_none " 

Usage:
  save_none(&str)

";

%feature("docstring")  simuPOP::load_none " 

Usage:
  load_none(&str, &offset)

";

%feature("docstring")  simuPOP::save_int " 

Usage:
  save_int(&str, *args)

";

%feature("docstring")  simuPOP::load_int " 

Usage:
  load_int(&str, &offset)

";

%feature("docstring")  simuPOP::save_long " 

Usage:
  save_long(&str, *args)

";

%feature("docstring")  simuPOP::load_long " 

Usage:
  load_long(&str, &offset)

";

%feature("docstring")  simuPOP::save_float " 

Usage:
  save_float(&str, *args)

";

%feature("docstring")  simuPOP::load_float " 

Usage:
  load_float(&str, &offset)

";

%feature("docstring")  simuPOP::save_string " 

Usage:
  save_string(&str, *args)

";

%feature("docstring")  simuPOP::load_string " 

Usage:
  load_string(&str, &offset)

";

%feature("docstring")  simuPOP::saveObj " 

Usage:
  saveObj(&str, *args)

";

%feature("docstring")  simuPOP::loadObj " 

Usage:
  loadObj(&vars, &offset)

";

%feature("docstring")  simuPOP::save_dict " 

Usage:
  save_dict(&str, *args)

";

%feature("docstring")  simuPOP::load_dict " 

Usage:
  load_dict(&vars, &offset)

";

%feature("docstring")  simuPOP::save_list " 

Usage:
  save_list(&str, *args)

";

%feature("docstring")  simuPOP::load_list " 

Usage:
  load_list(&vars, &offset)

";

%ignore simuPOP::mainVars();

%ignore simuPOP::moduleVars();

%ignore simuPOP::pyPopObj(void *p);

%ignore simuPOP::pyIndObj(void *p);

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

%ignore simuPOP::ostreamManager();

%feature("docstring")  simuPOP::rng " 

currently, return a global RNG.

Usage:
  rng()

";

%feature("docstring")  simuPOP::setRNG " 

set random number generator

Usage:
  setRNG(r, seed)

";

%feature("docstring")  simuPOP::listAllRNG " 

list all available RNG.

Usage:
  listAllRNG()

";

%feature("docstring")  simuPOP::bitSet " 

Random number
generator///////////////////////////////////////////////////////////
for output purpose only.

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

%ignore simuPOP::cnull();

%feature("docstring")  simuPOP::setLogOutput " 

set standard output to (default standard Python output)

Usage:
  setLogOutput(filename)

";

%ignore simuPOP::initialize();

%feature("docstring")  simuPOP::simuRev " 

Global debug and initialization related
functions///////////////////////////////////////////////////////////
return svn revision.

Usage:
  simuRev()

";

%feature("docstring")  simuPOP::simuVer " 

return version infomation of simuPOP

Usage:
  simuVer()

";

%feature("docstring")  simuPOP::optimized " 

Usage:
  optimized()

";

%feature("docstring")  simuPOP::supportXML " 

Usage:
  supportXML()

";

%feature("docstring")  simuPOP::alleleType " 

Usage:
  alleleType()

";

%feature("docstring")  simuPOP::compileCompiler " 

Usage:
  compileCompiler()

";

%feature("docstring")  simuPOP::compileDate " 

Usage:
  compileDate()

";

%feature("docstring")  simuPOP::compilePyVersion " 

Usage:
  compilePyVersion()

";

%feature("docstring")  simuPOP::compilePlatForm " 

Usage:
  compilePlatForm()

";

%feature("docstring")  is_carrayobject " 

Usage:
  is_carrayobject(*op)

";

%feature("docstring")  a_getitem " 

Usage:
  a_getitem(*ap, i)

";

%feature("docstring")  a_setitem " 

Usage:
  a_setitem(*ap, i, *v)

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

%feature("docstring")  carray_new " 

Usage:
  carray_new(*type, *args, *kwds)

";

%feature("docstring")  carray_init " 

Usage:
  carray_init(*type, *args, *kwds)

";

%feature("docstring")  newcarrayobject " 

Usage:
  newcarrayobject(*ptr, type, size)

";

%feature("docstring")  newcarrayiterobject " 

Usage:
  newcarrayiterobject(begin, end)

";

%feature("docstring")  getarrayitem " 

Usage:
  getarrayitem(*op, i)

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

%feature("docstring")  array_concat " 

Usage:
  array_concat(*a, *bb)

";

%feature("docstring")  array_repeat " 

Usage:
  array_repeat(*a, n)

";

%feature("docstring")  array_item " 

Usage:
  array_item(*a, i)

";

%feature("docstring")  array_slice " 

Usage:
  array_slice(*a, ilow, ihigh)

";

%feature("docstring")  array_ass_slice " 

Usage:
  array_ass_slice(*a, ilow, ihigh, *v)

";

%feature("docstring")  array_ass_item " 

Usage:
  array_ass_item(*a, i, *v)

";

%feature("docstring")  array_count " 

Usage:
  array_count(*self, *args)

";

%feature("docstring")  array_index " 

Usage:
  array_index(*self, *args)

";

%feature("docstring")  array_tolist " 

Usage:
  array_tolist(*self, *args)

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

%feature("docstring")  initcarray " 

Usage:
  initcarray(void)

";

%feature("docstring")  carray_length " 

Usage:
  carray_length(*a)

";

%feature("docstring")  carray_itemsize " 

Usage:
  carray_itemsize(*a)

";

%feature("docstring")  carray_type " 

Usage:
  carray_type(*a)

";

%feature("docstring")  carray_data " 

Usage:
  carray_data(*a)

";

%feature("docstring")  dvars " 

Usage:
  dvars()

";

%feature("docstring")  variables " 

Usage:
  variables(linux)

";

%feature("docstring")  binary " 

Usage:
  binary(port)

";

%feature("docstring")  help " 

Usage:
  help(macOS)

";

%feature("docstring")  simuPOP::FreqTrajectoryStoch " 

Usage:
  FreqTrajectoryStoch(freq, N=0, *NtFunc=NULL, s=[], *sFunc=NULL,
    T=100000)

";

%feature("docstring")  simuPOP::FreqTrajectoryMultiStoch " 

Usage:
  FreqTrajectoryMultiStoch(freq=[], N=0, *NtFunc=NULL, s=[],
    *sFunc=NULL, T=100000)

";

%feature("docstring")  simuPOP::FreqTrajectorySelSim " 

Usage:
  FreqTrajectorySelSim(sel, Ne, freq, dom_h, selection)

";

%feature("docstring")  simuPOP::FreqTrajectoryForward " 

Usage:
  FreqTrajectoryForward(lowbound, highbound, disAge, grate, N0,
    seleCo)

";

%feature("docstring")  simuPOP::LoadPopulation " 

Usage:
  LoadPopulation(&file, &format=\"auto\")

";

%feature("docstring")  simuPOP::LoadSimulator " 

Usage:
  LoadSimulator(&file, &mate, format=\"auto\")

";

%ignore simuPOP::haploKey(const vectori &seq);

%feature("docstring")  SWIG_TypeNameComp " 

Usage:
  SWIG_TypeNameComp(*f1, *l1, *f2, *l2)

";

%feature("docstring")  SWIG_TypeEquiv " 

Usage:
  SWIG_TypeEquiv(*nb, *tb)

";

%feature("docstring")  SWIG_TypeCompare " 

Usage:
  SWIG_TypeCompare(*nb, *tb)

";

%feature("docstring")  SWIG_TypeCheck " 

Usage:
  SWIG_TypeCheck(*c, *ty)

";

%feature("docstring")  SWIG_TypeCheckStruct " 

Usage:
  SWIG_TypeCheckStruct(*from, *into)

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

%feature("docstring")  SWIG_TypeClientData " 

Usage:
  SWIG_TypeClientData(*ti, *clientdata)

";

%feature("docstring")  SWIG_TypeNewClientData " 

Usage:
  SWIG_TypeNewClientData(*ti, *clientdata)

";

%feature("docstring")  SWIG_MangledTypeQueryModule " 

Usage:
  SWIG_MangledTypeQueryModule(*start, *end, *name)

";

%feature("docstring")  SWIG_TypeQueryModule " 

Usage:
  SWIG_TypeQueryModule(*start, *end, *name)

";

%feature("docstring")  SWIG_PackData " 

Usage:
  SWIG_PackData(*c, *ptr, sz)

";

%feature("docstring")  SWIG_UnpackData " 

Usage:
  SWIG_UnpackData(*c, *ptr, sz)

";

%feature("docstring")  SWIG_PackVoidPtr " 

Usage:
  SWIG_PackVoidPtr(*buff, *ptr, *name, bsz)

";

%feature("docstring")  SWIG_UnpackVoidPtr " 

Usage:
  SWIG_UnpackVoidPtr(*c, **ptr, *name)

";

%feature("docstring")  SWIG_PackDataName " 

Usage:
  SWIG_PackDataName(*buff, *ptr, sz, *name, bsz)

";

%feature("docstring")  SWIG_UnpackDataName " 

Usage:
  SWIG_UnpackDataName(*c, *ptr, sz, *name)

";

%feature("docstring")  PyString_FromFormat " 

Usage:
  PyString_FromFormat(*fmt, ...)

";

%feature("docstring")  SWIG_Python_ErrorType " 

Usage:
  SWIG_Python_ErrorType(code)

";

%feature("docstring")  SWIG_Python_AddErrorMsg " 

Usage:
  SWIG_Python_AddErrorMsg(*mesg)

";

%feature("docstring")  SWIG_Python_SetErrorObj " 

Usage:
  SWIG_Python_SetErrorObj(*errtype, *obj)

";

%feature("docstring")  SWIG_Python_SetErrorMsg " 

Usage:
  SWIG_Python_SetErrorMsg(*errtype, *msg)

";

%feature("docstring")  SWIG_Python_SetConstant " 

Usage:
  SWIG_Python_SetConstant(*d, *name, *obj)

";

%feature("docstring")  SWIG_Python_AppendOutput " 

Usage:
  SWIG_Python_AppendOutput(*result, *obj)

";

%feature("docstring")  SWIG_Python_UnpackTuple " 

Usage:
  SWIG_Python_UnpackTuple(*args, *name, min, max, **objs)

";

%feature("docstring")  SWIG_Py_Void " 

Usage:
  SWIG_Py_Void(void)

";

%feature("docstring")  SWIG_Python_CheckImplicit " 

Usage:
  SWIG_Python_CheckImplicit(*ty)

";

%feature("docstring")  SWIG_Python_ExceptionType " 

Usage:
  SWIG_Python_ExceptionType(*desc)

";

%feature("docstring")  PySwigClientData_New " 

Usage:
  PySwigClientData_New(*obj)

";

%feature("docstring")  PySwigClientData_Del " 

Usage:
  PySwigClientData_Del(*data)

";

%feature("docstring")  PySwigObject_long " 

Usage:
  PySwigObject_long(*v)

";

%feature("docstring")  PySwigObject_format " 

Usage:
  PySwigObject_format(*fmt, *v)

";

%feature("docstring")  PySwigObject_oct " 

Usage:
  PySwigObject_oct(*v)

";

%feature("docstring")  PySwigObject_hex " 

Usage:
  PySwigObject_hex(*v)

";

%feature("docstring")  PySwigObject_repr " 

Usage:
  PySwigObject_repr(*v)

";

%feature("docstring")  PySwigObject_print " 

Usage:
  PySwigObject_print(*v, *fp, flags)

";

%feature("docstring")  PySwigObject_str " 

Usage:
  PySwigObject_str(*v)

";

%feature("docstring")  PySwigObject_compare " 

Usage:
  PySwigObject_compare(*v, *w)

";

%feature("docstring")  _PySwigObject_type " 

Usage:
  _PySwigObject_type(void)

";

%feature("docstring")  PySwigObject_type " 

Usage:
  PySwigObject_type(void)

";

%feature("docstring")  PySwigObject_Check " 

Usage:
  PySwigObject_Check(*op)

";

%feature("docstring")  PySwigObject_New " 

Usage:
  PySwigObject_New(*ptr, *ty, own)

";

%feature("docstring")  PySwigObject_dealloc " 

Usage:
  PySwigObject_dealloc(*v)

";

%feature("docstring")  PySwigObject_append " 

Usage:
  PySwigObject_append(*v, *next)

";

%feature("docstring")  PySwigObject_next " 

Usage:
  PySwigObject_next(*v, args)

";

%feature("docstring")  PySwigObject_disown " 

Usage:
  PySwigObject_disown(*v, args)

";

%feature("docstring")  PySwigObject_acquire " 

Usage:
  PySwigObject_acquire(*v, args)

";

%feature("docstring")  PySwigObject_own " 

Usage:
  PySwigObject_own(*v, *args)

";

%feature("docstring")  PySwigObject_getattr " 

Usage:
  PySwigObject_getattr(*sobj, *name)

";

%feature("docstring")  PySwigPacked_print " 

Usage:
  PySwigPacked_print(*v, *fp, flags)

";

%feature("docstring")  PySwigPacked_repr " 

Usage:
  PySwigPacked_repr(*v)

";

%feature("docstring")  PySwigPacked_str " 

Usage:
  PySwigPacked_str(*v)

";

%feature("docstring")  PySwigPacked_compare " 

Usage:
  PySwigPacked_compare(*v, *w)

";

%feature("docstring")  _PySwigPacked_type " 

Usage:
  _PySwigPacked_type(void)

";

%feature("docstring")  PySwigPacked_type " 

Usage:
  PySwigPacked_type(void)

";

%feature("docstring")  PySwigPacked_Check " 

Usage:
  PySwigPacked_Check(*op)

";

%feature("docstring")  PySwigPacked_dealloc " 

Usage:
  PySwigPacked_dealloc(*v)

";

%feature("docstring")  PySwigPacked_New " 

Usage:
  PySwigPacked_New(*ptr, size, *ty)

";

%feature("docstring")  PySwigPacked_UnpackData " 

Usage:
  PySwigPacked_UnpackData(*obj, *ptr, size)

";

%feature("docstring")  _SWIG_This " 

Usage:
  _SWIG_This(void)

";

%feature("docstring")  SWIG_This " 

Usage:
  SWIG_This(void)

";

%feature("docstring")  SWIG_Python_GetSwigThis " 

Usage:
  SWIG_Python_GetSwigThis(*pyobj)

";

%feature("docstring")  SWIG_Python_AcquirePtr " 

Usage:
  SWIG_Python_AcquirePtr(*obj, own)

";

%feature("docstring")  SWIG_Python_ConvertPtrAndOwn " 

Usage:
  SWIG_Python_ConvertPtrAndOwn(*obj, **ptr, *ty, flags, *own)

";

%feature("docstring")  SWIG_Python_ConvertFunctionPtr " 

Usage:
  SWIG_Python_ConvertFunctionPtr(*obj, **ptr, *ty)

";

%feature("docstring")  SWIG_Python_ConvertPacked " 

Usage:
  SWIG_Python_ConvertPacked(*obj, *ptr, sz, *ty)

";

%feature("docstring")  SWIG_Python_NewShadowInstance " 

Usage:
  SWIG_Python_NewShadowInstance(*data, *swig_this)

";

%feature("docstring")  SWIG_Python_NewPointerObj " 

Usage:
  SWIG_Python_NewPointerObj(*ptr, *type, flags)

";

%feature("docstring")  SWIG_Python_NewPackedObj " 

Usage:
  SWIG_Python_NewPackedObj(*ptr, sz, *type)

";

%feature("docstring")  SWIG_Python_GetModule " 

Usage:
  SWIG_Python_GetModule(void)

";

%feature("docstring")  PyModule_AddObject " 

Usage:
  PyModule_AddObject(*m, *name, *o)

";

%feature("docstring")  SWIG_Python_DestroyModule " 

Usage:
  SWIG_Python_DestroyModule(*vptr)

";

%feature("docstring")  SWIG_Python_SetModule " 

Usage:
  SWIG_Python_SetModule(*swig_module)

";

%feature("docstring")  SWIG_Python_AddErrMesg " 

Usage:
  SWIG_Python_AddErrMesg(*mesg, infront)

";

%feature("docstring")  SWIG_Python_ArgFail " 

Usage:
  SWIG_Python_ArgFail(argnum)

";

%feature("docstring")  PySwigObject_GetDesc " 

Usage:
  PySwigObject_GetDesc(*self)

";

%feature("docstring")  SWIG_Python_TypeError " 

Usage:
  SWIG_Python_TypeError(*type, *obj)

";

%feature("docstring")  SWIG_Python_MustGetPtr " 

Usage:
  SWIG_Python_MustGetPtr(*obj, *ty, argnum, flags)

";

%feature("docstring")  SWIG_TypeQuery " 

Usage:
  SWIG_TypeQuery(*name)

";

%feature("docstring")  SWIG_MangledTypeQuery " 

Usage:
  SWIG_MangledTypeQuery(*name)

";

%feature("docstring")  simuPOP::newcarrayobject " 

Usage:
  newcarrayobject(*buf, type, size)

";

%feature("docstring")  simuPOP::newcarrayiterobject " 

Usage:
  newcarrayiterobject(begin, end)

";

%feature("docstring")  simuPOP::is_carrayobject " 

Usage:
  is_carrayobject(*)

";

%feature("docstring")  simuPOP::carray_length " 

Usage:
  carray_length(*a)

";

%feature("docstring")  simuPOP::carray_itemsize " 

Usage:
  carray_itemsize(*a)

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

%feature("docstring")  simuPOP::TurnOnDebug " 

set debug code, default to turn all code on

Usage:
  TurnOnDebug(code)

";

%feature("docstring")  simuPOP::TurnOffDebug " 

turn off debug, default to turn all code off

Usage:
  TurnOffDebug(code)

";

%feature("docstring")  simuPOP::debug " 

test if one code is turned on

Usage:
  debug(code)

";

%feature("docstring")  simuPOP::ListDebugCode " 

show all dbg codes (print to cout)

Usage:
  ListDebugCode()

";

%feature("docstring")  simuPOP::dbgString " 

dbg string for a code

Usage:
  dbgString(code)

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

%ignore simuPOP::PyObj_Is_AlleleNumArray(PyObject *obj);

%ignore simuPOP::Int_Vec_As_NumArray(vectori::iterator begin, vectori::iterator end);

%ignore simuPOP::Double_Vec_As_NumArray(vectorf::iterator begin, vectorf::iterator end);

%ignore simuPOP::Allele_Vec_As_NumArray(GenoIterator begin, GenoIterator end);

%ignore simuPOP::NumArray_Size(PyObject *obj);

%ignore simuPOP::NumArray_Data(PyObject *obj);

%feature("docstring")  simuPOP::save_none " 

Usage:
  save_none(&str)

";

%feature("docstring")  simuPOP::load_none " 

Usage:
  load_none(&str, &offset)

";

%feature("docstring")  simuPOP::save_int " 

Usage:
  save_int(&str, *args)

";

%feature("docstring")  simuPOP::load_int " 

Usage:
  load_int(&str, &offset)

";

%feature("docstring")  simuPOP::save_long " 

Usage:
  save_long(&str, *args)

";

%feature("docstring")  simuPOP::load_long " 

Usage:
  load_long(&str, &offset)

";

%feature("docstring")  simuPOP::save_float " 

Usage:
  save_float(&str, *args)

";

%feature("docstring")  simuPOP::load_float " 

Usage:
  load_float(&str, &offset)

";

%feature("docstring")  simuPOP::save_string " 

Usage:
  save_string(&str, *args)

";

%feature("docstring")  simuPOP::load_string " 

Usage:
  load_string(&str, &offset)

";

%feature("docstring")  simuPOP::saveObj " 

Usage:
  saveObj(&str, *args)

";

%feature("docstring")  simuPOP::loadObj " 

Usage:
  loadObj(&vars, &offset)

";

%feature("docstring")  simuPOP::save_dict " 

Usage:
  save_dict(&str, *args)

";

%feature("docstring")  simuPOP::load_dict " 

Usage:
  load_dict(&vars, &offset)

";

%feature("docstring")  simuPOP::save_list " 

Usage:
  save_list(&str, *args)

";

%feature("docstring")  simuPOP::load_list " 

Usage:
  load_list(&vars, &offset)

";

%ignore simuPOP::mainVars();

%ignore simuPOP::moduleVars();

%ignore simuPOP::pyPopObj(void *p);

%ignore simuPOP::pyIndObj(void *p);

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

%ignore simuPOP::ostreamManager();

%feature("docstring")  simuPOP::rng " 

currently, return a global RNG.

Usage:
  rng()

";

%feature("docstring")  simuPOP::setRNG " 

set random number generator

Usage:
  setRNG(r, seed)

";

%feature("docstring")  simuPOP::listAllRNG " 

list all available RNG.

Usage:
  listAllRNG()

";

%feature("docstring")  simuPOP::bitSet " 

Random number
generator///////////////////////////////////////////////////////////
for output purpose only.

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

%ignore simuPOP::cnull();

%feature("docstring")  simuPOP::setLogOutput " 

set standard output to (default standard Python output)

Usage:
  setLogOutput(filename)

";

%ignore simuPOP::initialize();

%feature("docstring")  simuPOP::simuRev " 

Global debug and initialization related
functions///////////////////////////////////////////////////////////
return svn revision.

Usage:
  simuRev()

";

%feature("docstring")  simuPOP::simuVer " 

return version infomation of simuPOP

Usage:
  simuVer()

";

%feature("docstring")  simuPOP::optimized " 

Usage:
  optimized()

";

%feature("docstring")  simuPOP::supportXML " 

Usage:
  supportXML()

";

%feature("docstring")  simuPOP::alleleType " 

Usage:
  alleleType()

";

%feature("docstring")  simuPOP::compileCompiler " 

Usage:
  compileCompiler()

";

%feature("docstring")  simuPOP::compileDate " 

Usage:
  compileDate()

";

%feature("docstring")  simuPOP::compilePyVersion " 

Usage:
  compilePyVersion()

";

%feature("docstring")  simuPOP::compilePyVersion " 

Usage:
  compilePyVersion()

";

