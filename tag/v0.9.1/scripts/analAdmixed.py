#!/usr/bin/env python
'''
This script implements a few quantitative trait, penetrance and ascertainment
models that are used to analyze populations generated from simuAdmixture.py.
Because arbitrary models can be applied to such a population, this script
is intended only to
1) record models that are used by examples in Peng 2008
2) serve as examples on how to apply such models.


A multi-locus quantitative trait model
---------------------------------------

Using this method, a list of chromosomes with disease pre-disposing locus
is given (if there are more than one disease locus, repeat the chromosome
number), along with desired allele frequency (and allowed window around this
frequency). From given positions (default to 0.5, from the middle of 
each chromosome), this script will locate loci that satify allele frequency
requirement.

Then, a list of 'percentage of variation that is explained by each disease
locus' is required. An additive quantitative trait model is used which 
yields trait value -a, 0 and a for genotype AA, Aa and aa respectively.
Assume allele frequency p, a locus would contribute 2p(1-p)a^2 variance to
the total trait variance. With given variance v, a is calculated using
sqrt(v/(2p(1-p))).

The rest of the variance is simulated using a normal distribution.


A trait-derived affection status model
---------------------------------------

Using a above quantitative trait model, a cutoff value can be given which
assign affection status to individuals according to their trait value.


With a given penetrance model
------------------------------

It is also possible to use an explicit penetrance model to assign affection
status. Using this approach, a list of markers (by names) is expected and 
a penetrance matrix is given. For a two-locus model, the penetrance matrix
should be given as p_AABB, p_AABb, p_AAbb, p_AaBB, p_AaBb, p_Aabb, p_aaBB,
p_aaBb, p_aabb, which p_XXXX is the probability of affected for genotype
XXXX.


Draw a Case-control sample
------------------------------

If affection status is assigned to each individual, case-control samples
can be drawn from the resulting population. One can simply specify number
of cases and controls, or a list of cases and controls from each population.
For example, [0,0,500], [0,0,500] would ignore the first two populations
and draw 500 cases and 500 controls.


Draw a random sample
---------------------

If a quantitative trait is specified, a random sample can be drawn from
the admixed population. Given number of random individuals can be drawn from
from all populations (a single number is given), or from each population
(a list is given).

The samples are saved in Merlin QTDT format, which can be easily converted
to other formats.

Bo Peng (bpeng@mdanderson.org)

$LastChangedDate: 2007-10-12 14:48:57 -0500 (Fri, 12 Oct 2007) $
$Rev: 1090 $

'''

from simuOpt import *
setOptions(alleleType='binary')
from simuPOP import *
from hapMapUtil import getMarkersFromName, getMarkersFromRange, evolveHapMap

import os, sys, types, exceptions, getopt, math, random
from simuUtil import SaveQTDT

HapMap_pops = ['CEU', 'YRI', 'JPT+CHB']

options = [
    {'arg': 'h',
     'longarg': 'help',
     'default': False, 
     'description': 'Print this usage message.',
     'allowedTypes': [types.NoneType, type(True)],
     'jump': -1                    # if -h is specified, ignore any other parameters.
    },
    {'separator': 'Population to analyze'},
    {'longarg': 'dir=',
     'default': 'simu',
     'useDefault': False,
     'allowedTypes': [types.StringType],
     'label': 'Analyses directory',
     'description': '''Directory where all the samples and results will be stored''',
    },
    {'longarg': 'popName=',
     'default': 'admixed.bin',
     'label': 'Population file',
     'useDefault': False,
     'description': '''Name of the population to be analyzed. It is by default relative to parameter
                dir, but an absolute path can be given.''',
     'allowedTypes': [types.StringType],
    },
    {'separator': 'Quantitative trait model'},
    {'longarg': 'DSLtrait=',
     'default': [],
     'label': 'Disease susceptibility loci',
     'useDefault': True,
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Names of disease susceptibility loci for the quantitative trait model.
                If given, parameters chromWithDSL, freqDSL and freqDev will be ignored. 
                Otherwise, these DSL will be determined by allele frequency automatically.'''
    },
    {'longarg': 'chromWithDSL=',
     'default': [1,2],
     'label': 'Chromosomes with DSL',
     'useDefault': True,
     'allowedTypes': [types.TupleType, types.ListType],
     'description': '''Chromosomes with DSL. The chromosomes are counted from the loaded
                population and indexed from 1. For example, if chromosomes 5, 6, 8 are
                used in the simulation, [1, 3] here would refer to chromosome 5 and 8.'''
    },
    {'longarg': 'freqDSL=',
     'default': 0.1,
     'label': 'MAF of DSL',
     'useDefault': True,
     'allowedTypes': [types.FloatType],
    },
    {'longarg': 'freqDev=',
     'default': 0.01,
     'useDefault': True,
     'label': 'Allowed deviation from dslFreq',
     'allowedTypes': [types.FloatType],
    },
    {'longarg': 'dslVar=',
     'default': [0.05, 0.1],
     'label': 'proportion of Variance to explain',
     'allowedTypes': [types.ListType, types.TupleType],
     'description': '''Proption of variance explained by each DSL,
                The a in the additive formula is determined by
                this variable and allale frequency'''
    },
    {'longarg': 'cutoff=',
     'default': 0,
     'useDefault': True,
     'label': 'Cutoff value to determine affection status',
     'allowedTypes': [types.IntType, types.FloatType],
     'description': 'Cutoff value used to determine affection status'
    },      
    {'separator': 'Penetrance model'},
    {'longarg': 'DSLpene=',
     'label': 'Disease susceptibility loci',
     'default': [],
     'useDefault': True,
     'allowedTypes': [types.TupleType, types.ListType],
     'description': '''A list of markers (by name) that will be used to 
                determine affection status'''
    },
    {'longarg': 'peneFunc=',
     'default': 'None',
     'label': 'Penetrance function',
     'allowedTypes': [types.StringType],
     'description': ''' Penetrance functions to be applied to the final
        population. Two penetrance fucntions are provided, namely recessive
        or additive single-locus model with heterogeneity multi-locus model. 
        You can define another customized penetrance functions by modifying
        this script. ''',
     'validate':    valueOneOf(['recessive', 'additive', 'custom', 'None']),
     'chooseOneOf': [ 'recessive', 'additive', 'custom', 'None']
    },
    {'longarg': 'peneParam=',
     'default': [0.5],
     'label': 'Penetrance parameters',
     'description': '''Penetrance parameter for all DSL. An array of parameter 
        can be given to each DSL. The meaning of this parameter differ by 
        penetrance model. For a recessive model, the penetrance is 0,0,p for 
        genotype AA,Aa,aa (a is disease allele) respectively. For an additive 
        model, the penetrance is 0,p/2,p respectively.''',
     'allowedTypes': [types.ListType, types.TupleType],
     'validate':   valueOneOf([ 
             valueBetween(0,1), valueListOf(valueBetween(0,1))] )
    },
    {'separator': 'Draw samples'},
    {'longarg': 'ccSampleSize=',
     'default': [100,100],
     'label': 'Case-control sample size',
     'description': '''Sample size. If sampleType == 'random' and a single number s
                is given, a random sample of size s will be drawn from the whole
                population. If a list [s1, s2, ...] is given, s1 random individuals will
                be drawn from the first population, s2 random individuals will be
                drawn from the second population, and so on.
                If sampleType == 'case-control' and a two-number list 
                [s1, s2] is given, s1 case and s2 controls will be drawn from the
                whole population. If [list1, list2] is given, list1 is the number
                of cases from each subpopulation, and list2 is the number of controls
                from each subpopulation. The second form obviously can not be
                used for the Hybrid Isolation model because there is no subpopulation
                in the final population. Note that elements in list1 and list2 can
                be zero so you can sample from only one of the subpopulations.''',
      'allowedTypes': [types.IntType, types.LongType, types.TupleType, types.ListType],
      'validate': valueOr(valueGE(0), valueListOf(valueOr(valueGE(0), valueListOf(valueGE(0)))))
    },
    {'longarg': 'ccSampleName=',
     'default': 'case-control',
     'useDefault': True,
     'allowedTypes': [types.StringType],
     'label': 'Name of the case-control sample',
     'description': 'Name of the case-control sample',
    },
    {'longarg': 'randomSampleSize=',
     'default': 100,
     'label': 'Random sample size',
     'description': '''Size of a random sample''',
     'allowedTypes': [types.IntType, types.LongType, types.TupleType, types.ListType],
     'validate': valueOr(valueGE(0), valueListOf(valueOr(valueGE(0), valueListOf(valueGE(0)))))
    },
    {'longarg': 'randomSampleName=',
     'default': 'random',
     'useDefault': True,
     'allowedTypes': [types.StringType],
     'label': 'Name of the random sample',
     'description': 'Name of the random sample',
    },
]


def setQuanTrait(pop, DSLtrait, chromWithDSL, p, sd, vars, cutoff, name):
    '''Set quantitative trait and affection status if a cutoff value
    is given
    
    DSLtrait: names of each disease susceptibility loci. If given,
        ignore parameters chromWithDSL, p, sd.
    
    chromWithDSL: chromosomes with DSL, chromosomes should be
        indexed from 0

    p: target allele frequency of the DSL

    sd: allowed deviation from p

    vars: variance for each marker

    cutoff: set an individual as 'affected' if its trait value is greater
        than this cutoff value.
    '''
    #
    if len(DSLtrait) == 0:
        numDSL = len(chromWithDSL)
        DSL = []
        sign = []
        for ch1 in chromWithDSL:
            ch = ch1 - 1
            DSL.append(pop.chromBegin(ch) + pop.numLoci(ch)/2)
            for i in range(DSL[-1], pop.chromEnd(ch)):
                af = 1 - pop.dvars().alleleFreq[i][0]
                # allele 1 is the minor allele
                if af > p - sd and af < p + sd:
                    DSL[-1] = i
                    sign.append(1)
                    break
                # allele 0 is the minor allele
                elif 1 - af > p - sd and 1 - af < p + sd:
                    DSL[-1] = i
                    sign.append(-1)
                    break
    else:
        try:
            DSL = pop.lociByNames(DSLtrait)
            numDSL = len(DSL)
        except:
            print 'Can not find one of the DSL %s in this population' % \
                ', '.join(DSLtrait)
            sys.exit(1)
        sign = []
        for i in DSL:
            af = 1 - pop.dvars().alleleFreq[i][0]
            # allele 1 is the minor allele
            if af < 0.5:
                sign.append(1)
            # allele 0 is the minor allele
            else:
                sign.append(-1)
    maf = [min(pop.dvars().alleleFreq[x][0], 1-pop.dvars().alleleFreq[x][0]) for x in DSL]
    print 'Using DSL %s with minor allele frequency %s' % (DSL, maf)
    #
    dsl = open(os.path.join(name, 'DSL.lst'), 'w')
    print >> dsl, "name\tchrom\tindex\tlocation\tminor allele frequency\tpercentage of variation"
    for i,d in enumerate(DSL):
        print >> dsl, "%s\t%d\t%d\t%.5f\t%.3f\t%.3f" % (pop.locusName(d),
            pop.chromLocusPair(d)[0] + 1, d - pop.chromBegin(pop.chromLocusPair(d)[0]), pop.locusPos(d),
            pop.dvars().alleleFreq[d][0], vars[i])
    dsl.close()
    # applying quantitative trait
    table = [[-math.sqrt(vars[i]/(2.*maf[i]*(1-maf[i])))*sign[i], 
        0, sign[i]*math.sqrt(vars[i]/(2.*maf[i]*(1-maf[i])))] \
            for i in range(numDSL)]
    def traitFunc(geno):
        return sum([table[x][geno[2*x]+geno[2*x+1]] for x in range(numDSL)]) + \
            random.normalvariate(0, math.sqrt(1-sum(vars)))
    pop.addInfoField('qtrait')
    print "Assigning quantitative trait..."
    PyQuanTrait(pop, loci=DSL, func=traitFunc)
    print "Mean quantitative trait is ", sum(pop.indInfo('qtrait', False))/pop.popSize()
    print "Minimal quantitative trait is ", min(pop.indInfo('qtrait', False))
    print "Maximal quantitative trait is ", max(pop.indInfo('qtrait', False))
    # set affection status
    qidx = pop.infoIdx('qtrait')
    for ind in pop.individuals():
        if ind.info(qidx) > cutoff:
            ind.setAffected(True)
        else:
            ind.setAffected(False)
    print 'Using cutoff value %.2f' % cutoff
    Stat(pop, numOfAffected=True)
    print "There are %d (%.2f percent) affected individuals." % (pop.dvars().numOfAffected, pop.dvars().numOfAffected*100.0/pop.popSize())


# penetrance generator functions. They will return a penetrance function
# with given penetrance parameter
def recessive(para):
    ''' recessive single-locus, heterogeneity multi-locus '''
    def func(geno):
        val = 1
        for i in range(len(geno)/2):
            if geno[i*2] + geno[i*2+1] == 2:
                val *= 1 - para[i]
        return 1-val
    return func
    

def additive(para):
    ''' additive single-locus, heterogeneity multi-locus '''
    def func(geno):
        val = 1
        for i in range(len(geno)/2):
            val *= 1 - (geno[i*2]+geno[i*2+1])*para[i]/2.
        return 1-val
    return func


# if you need some specialized penetrance function, modify this
# function here.
# NOTE:
# 
# 1. geno is the genptype at DSL. For example, if your DSL is [5,10]
#     geno will be something like [0,1,1,1] where 0,1 is the genotype at 
#     locus 5 and 1,1 is the genotype at locus 10.
# 2. in simuComplexDisease.py, 0 is wild type, 1 is disease allele.
def custom(para):
    ''' quantitative trait '''
    def func(geno):
        return 1
    return func


short_desc = '''This program analyze an admixed population based on 
two or more HapMap populations. Please refer to simuAdmixture.py on
how to simulate such a population'''

# determine which script to run.
if __name__ == '__main__':
    # 
    # seed population does not exist, generate it
    allParam = getParam(options, short_desc, __doc__, nCol=2)
    # when user click cancel ...
    if len(allParam) == 0:
        sys.exit(1)
    # -h or --help
    if allParam[0]:    
        print usage(options, __doc__)
        sys.exit(0)
    # 
    (name, popName,
      DSLtrait, chromWithDSL, freqDSL, freqDev, dslVar, cutoff,
      DSLpene, peneFunc, peneParam, 
      ccSampleSize, ccSampleName, randomSampleSize, randomSampleName) \
        = allParam[1:]
    # simulation name?
    if not os.path.isdir(name):
        print 'Directory %s does not exist.'
        sys.exit(1)
    cfgFile = os.path.join(name, name + '_anal.cfg')
    print 'Save configuration to', cfgFile
    # save current configuration
    saveConfig(options, cfgFile, allParam)
    #
    if os.path.isabs(popName):
        popFile = popName
    else:
        popFile = os.path.join(name, popName)
    #
    if not os.path.isfile(popFile):
        print 'Population file %s does not exist.' % popFile
        sys.exit(1)
    #
    print 'Loading population %s for analyses ' % popFile
    pop = LoadPopulation(popFile)
    #
    print 'ch with dsl', chromWithDSL
    if len(chromWithDSL) > 0:
        # assign case/control status and quantitative trait
        setQuanTrait(pop, DSLtrait, chromWithDSL, freqDSL, freqDev, 
            dslVar, cutoff, name)
    #
    if len(DSLpene) > 0:
        if 'recessive' == peneFunc:
            pene_func = recessive(para)
        elif 'additive' == peneFunc:
            pene_func = additive(para)
        elif 'custom' == peneFunc:
            pene_func = custom(para)
        else:
            pene_func = None
        #
        if pene_func is not None:
            PyPenetrance(pop, pop.lociByNames(DSLpene), func=pene_func)
    #
    # draw sample
    def comb(geno):
        return geno[0]+geno[1]
    if ccSampleSize != [0, 0]:
        print 'Drawing a case control sample with size ', ccSampleSize
        (samples,) = CaseControlSample(pop, cases=ccSampleSize[0],
            controls=ccSampleSize[1])
        print "Saving case control sample to %s " % os.path.join(name, 'case_control')
        SaveQTDT(samples, output=os.path.join(name, 'case_control'), affectionCode=['1', '2'], 
                fields=['affection', 'qtrait'], combine=comb, header=True)
    if randomSampleSize != 0:
        print 'Drawing a random sample with size ', randomSampleSize
        (ran,) = RandomSample(pop, size=randomSampleSize)
        # random sample
        print "Saving random sample to %s ..." % os.path.join(name, 'random')
        SaveQTDT(ran, output=os.path.join(name, 'random'), affectionCode=['1', '2'], 
            fields=['qtrait'], combine=comb, header=True) 


