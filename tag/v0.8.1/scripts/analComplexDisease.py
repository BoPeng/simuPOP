#!/usr/bin/env python
#
# Purpose: Analyze population generated by simuComplexDisease.py
#
# Bo Peng (bpeng@rice.edu)
#
# $LastChangedDate: 2005-10-31 17:29:34 -0600 (Mon, 31 Oct 2005) $
# $Rev: 78 $
#
# Known bugs:
#     None
# 

"""
This program analyzes a population generated by simuComplexDisease.py. 
It basically provide a user interface to relevant functions in 
simuUtil.py.

The program is written in Python using simuPOP modules. For more information,
please visit simuPOP website http://simupop.sourceforge.net .
"""

import simuOpt, simuUtil
import os, sys, types, exceptions, os.path, operator, time
#
# declare all options. getParam will use these information to get parameters
# from a tk/wxPython-based dialog, command line, config file or user input
#
# detailed information about these fields is given in the simuPOP reference
# manual.
options = [
    {'arg': 'h',
     'longarg': 'help',
     'useDefault': True,
     'default': False, 
     'description': 'Print this usage message.',
     'allowedTypes': [types.NoneType, type(True)],
     'jump': -1                    # if -h is specified, ignore any other parameters.
    },
    {'longarg': 'markerType=',
     'default': 'SNP',
     'allowedTypes': [types.StringType],
     'label': 'Marker type used',
     'description': '''Marker type used to generated the sample. This is
        important since the file formats are not compatible between 
        binary and standard simuPOP modules''',
     'validate':    simuOpt.valueOneOf([ 'microsatellite', 'SNP']),
     'chooseOneOf': ['microsatellite', 'SNP']
    }, 
    {'longarg': 'dataset=',
     'default': 'simu.txt',
     'allowedTypes': [types.StringType],
     'label': 'Dataset to analyze',
     'description': 'Dataset generated by simuComplexDisease.py. ',
     'validate':    simuOpt.valueValidFile()
    },
    {'longarg': 'peneFunc=',
     'default': 'additive',
     'label': 'Penetrance function',
     'allowedTypes': [types.StringType],
     'description': ''' Penetrance functions to be applied to the final
        population. Two penetrance fucntions are provided, namely recessive
        or additive single-locus model with heterogeneity multi-locus model. 
        You can define another customized penetrance functions by modifying
        this script. ''',
     'validate':    simuOpt.valueOneOf(['recessive', 'additive', 'custom', 'None']),
     'chooseOneOf': [ 'recessive', 'additive', 'custom']
    },
    {'longarg': 'qtraitFunc=',
     'default': 'None',
     'label': 'Quantitative trait model',
     'allowedTypes': [types.StringType],
     'description': '''How to calculate quantitative trait. There is no default
        function so one has to provide a customized one.
        ''',
     'validate':    simuOpt.valueOneOf(['custom', 'None']),
     'chooseOneOf': [ 'custom', 'None']
    },
    {'longarg': 'parameter=',
     'default': [0.5],
     'label': 'Penetrance parameters',
     'description': '''Penetrance parameter for all DSL. An array of parameter 
        can be given to each DSL. The meaning of this parameter differ by 
        penetrance model. For a recessive model, the penetrance is 0,0,p for 
        genotype AA,Aa,aa (a is disease allele) respectively. For an additive 
        model, the penetrance is 0,p/2,p respectively.''',
     'allowedTypes': [types.ListType, types.TupleType],
     'validate':   simuOpt.valueOneOf([ 
             simuOpt.valueBetween(0,1), simuOpt.valueListOf(simuOpt.valueBetween(0,1))] )
    },
    {'longarg': 'sampleSize=',
     'default': 800,
     'label':    'Sample size',
     'allowedTypes':    [types.IntType, types.LongType],
     'description':    '''Size of the samples, that will mean N/4 affected 
        sibpair families (of size 4), N/2 cases and controls. Up to N total number
        of individuals for large pedigrees.''',
     'validate':    simuOpt.valueGT(1)
    },
    {'longarg': 'outputDir=',
     'default': '.',
     'allowedTypes': [types.StringType],
     'label': 'Output directory',
     'description': 'Directory into which the datasets will be saved. ',
     'validate':    simuOpt.valueValidDir()
    },
    {'longarg': 'analyses=',
     'default': [],
     'label': 'Gene mapping methods',
     'allowedTypes': [types.TupleType, types.ListType],
     'description': ''' Gene mapping methods to apply. geneHunter is needed for
        TDT and Linkage methods, and R/RPy are needed for chisq association 
        tests.''',
     'allowedTypes': [types.ListType, types.TupleType],
     'chooseFrom': [ 
        'affectedSibs_GH_TDT', 
        'affectedSibs_GH_linkage', 
        'affectedSibs_merlin_linkage', 
        'caseControl_Association', 
        'qtraitSibs_merlin_reg',
        'qtraitSibs_merlin_vc',
        'largePeds_merlin_reg', 
        'largePeds_merlin_vc',
        ]
    },
    # another two hidden parameter
    {'longarg': 'saveConfig=',
     'default': 'anal.cfg',
     'allowedTypes': [types.StringType, types.NoneType],
     'label': 'Save configuration',
     'description': 'Save current paremeter set to specified file.'
    },
    {'arg': 'v',
     'longarg': 'verbose',
     'useDefault': True,
     'default': False,
     'allowedTypes': [types.NoneType, types.IntType],
     'description': 'Verbose mode.'
    },
]

outputVars = {
    # basic information
    'dataset': 'name of the datafile',
    'logfile': 'Name of the log file generated by simuComplexDisease.py',
    'DSL': 'Location of DSL in the populations. Use them for the whole population',
    'DSLafter': 'Paramter you give to simuComplexDisease.py. They are the one you should use for the samples',
    # evolution related parameters
    'numSubPop': 'Number of subpopulation',
    'recRate': 'recombination rate',
    'initSize': 'Size of initial founder population',
    'endingSize': 'Population size of the last generation',
    'burninGen': 'Length of burnin stage',
    'splitGen': 'Generation at which population is split to subpopulations',
    'mixingGen': 'Generation at which migration is allowed',
    'endingGen': 'total evolution length',
    'migrModel': 'migration model',
    'migrRate': 'migration rate',
    'mutaModel': 'mutation model',
    'mutaRate': 'mutation rate',
    # population statistics
    'Fst': 'F_st: measure of population differentiation',
    'Fis': 'F_is',
    'Fit': 'F_it',
    'AvgFst': 'Average of F_st estimated from all loci',
    'AvgFis': 'Average of F_is estimated from all loci',
    'AvgFit': 'Average of F_it estimated from all loci',
    'AvgHetero': 'Averge heterozygosity',
    'Fprime': '',
    'K': 'Disease prevalence (penetrance dependent)',
    'Ks': 'Sibling recurrence risk',
    'Ls': 'Sibling recurrence risk ratio', 
    'P11': 'Pr( (N,N) | affected ): proportion of NN (normal) genotype among affected individuals',
    'P12': 'Pr( (N,S) | affected ): proportion of NS (normal, susceptible) or SN genotype among affected individuals',
    'P22': 'Pr( (S,S) | affected ): proportion of SS (susceptible) genotype among affected individuals',
    'LD': 'Linkage disequilibrium between DSL and its surrounding markers',
    'LD_prime': "D' measure of Linkage disequilibrium",
    'R2': 'R^2 measure of Linkage disequilibrium',
    'alleleFreq': 'allele frequency of DSL',    
    'mono': 'Number of fixed or lost loci on the last chromosome',
    # result of gene mapping methods
    'affectedSibs_GH_TDT' : '', 
    'affectedSibs_GH_linkage' : '', 
    'affectedSibs_merlin_linkage' : '', 
    'caseControl_Association' : '', 
    'qtraitSibs_merlin_reg' : '',
    'qtraitSibs_merlin_vc' : '',
    'largePeds_merlin_reg' : '', 
    'largePeds_merlin_vc' : '',
}

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


# this customized qtrait function provide quantitative trait
# for example2 in the plos paper.
# The basic form is
#
# Q = X1 + X2 + X3 + Env
# 
def customQtrait(para):
    ''' quantitative trait '''
    def func(geno):
        x = random.normalvariate(0, 0.5)
        for i in range(len(geno)/2):
            x += random.normalvariate(geno[2*i]+geno[2*i+1], 0.5)
        return x
    return func


def getOptions(details=__doc__):
    ''' get options from options structure,
        if this module is imported, instead of ran directly,
        user can specify parameter in some other way.
    '''
    # get all parameters, __doc__ is used for help info
    allParam = simuOpt.getParam(options, 
        '''    This program simulates the evolution of a complex common disease, subject 
     to the impact of mutation, migration, recombination and population size change. 
     Click 'help' for more information about the evolutionary scenario.''', details, nCol=1)
    # when user click cancel ...
    if len(allParam) == 0:
        sys.exit(1)
    # -h or --help
    if allParam[0]:    
        print simuOpt.usage(options, __doc__)
        sys.exit(0)
    # --saveConfig
    simuOpt.saveConfig(options, allParam[-2], allParam)
    # --verbose or -v (these is no beautifying of [floats]
    if allParam[-1]:                 # verbose
        simuOpt.printConfig(options, allParam)
    # return the rest of the parameters
    return allParam[1:-2]

def popStat(pop):
    'Calculate population statistics '
    # K -- populaiton prevalance
    print "Calculating population statistics "
    Stat(pop, numOfAffected=True, alleleFreq=range(pop.totNumLoci()))
    result = {}
    result['K'] = pop.dvars().numOfAffected * 1.0 / pop.popSize()
    # P11 = [ ] = proportion of 11 | affected, 
    # P12 = [ ] = proportion of 12 | affected
    DSL = pop.dvars().DSL
    P11 = [0.]*len(DSL)
    P12 = [0.]*len(DSL)
    P22 = [0.]*len(DSL)
    for ind in range(pop.popSize()):
        if pop.individual(ind).affected():
            for x in range(len(DSL)):
                s1 = pop.individual(ind).allele(DSL[x], 0)
                s2 = pop.individual(ind).allele(DSL[x], 1)
                if s1 == 0 and s2 == 0:
                    P11[x] += 1
                elif s1 == 1 and s2 == 1:
                    P22[x] += 1
                else:
                    P12[x] += 1
    # calculate the number of fixed/lost loci
    result['mono'] = 0
    for x in range(pop.chromBegin(pop.numChrom()-1), pop.chromEnd(pop.numChrom()-1)):
        #print  pop.dvars().alleleNum[x][0]
        if pop.dvars().alleleNum[x][0] == 0 or pop.dvars().alleleNum[x][0] == pop.popSize() * pop.ploidy():
            result['mono'] += 1
    N = pop.dvars().numOfAffected
    # some analysis does not have affected, or all affected
    if N == 0 or N == pop.popSize():
        return result
    result['P11'] = [ x/N for x in P11 ]
    result['P12'] = [ x/N for x in P12 ]
    result['P22'] = [ x/N for x in P22 ]
    result['Fprime'] = [ (P12[x]/2. + P22[x])/N for x in range(len(DSL)) ]
    # Ks = Pr(Xs=1 | Xp=1 ) = Pr(Xs=1, Xp=1) | P(Xp=1)
    Xsp = 0.
    for ind in range(pop.popSize()/2):
        s1 = pop.individual(ind*2).affected()
        s2 = pop.individual(ind*2+1).affected()
        if s1 and s2:
            Xsp += 1
    result['Ks'] = 2*Xsp / pop.dvars().numOfAffected
    # Lambda = Ks/K
    result['Ls'] = result['Ks'] / result['K']
    return result


if __name__ == '__main__':
    allParam = getOptions()
    # unpack options
    (markerType, dataset, peneFunc, qtraitFunc, parameter, N, outputDir, analyses) = allParam
    # load simuPOP libraries
    if markerType == 'microsatellite':
        simuOpt.setOptions(alleleType='short', quiet=True)
    else:
        simuOpt.setOptions(alleleType='binary', quiet=True)
    #
    from simuPOP import *
    from simuUtil import *
    #
    print "Loading population", dataset
    pop = LoadPopulation(dataset)
    nDSL = len(pop.dvars().DSL)
    if len(parameter) == 1:
        para = parameter * nDSL
    elif len(parameter) == nDSL:
        para = parameter
    else:
        print "Length of penetrance/quantitative trait parameter should be one or the number of DSL"
        sys.exit(0)
    #
    res = {}
    res.update({
        'dataset':  dataset,
        'logfile':  dataset[0:-4] + '.log',
        'alleleFreq': [1- pop.dvars().alleleFreq[i][0] for i in pop.dvars().DSL],
        'pene==Func': peneFunc,
    })
    # calculate LD (should be already there, just to make sure)
    Stat(pop, LD=[[x,x+1] for x in pop.dvars().DSL] + [[x,x+5] for x in pop.dvars().DSL])
    # get all the variables from pop
    res.update(pop.vars())

    #res = analyzePopulation(dataset, peneFunc, qtraitFunc, parameter, N, numSample, outputDir, 
    #    analyses)
    # penetrance function.
    # for debugging purpose
    keep_temp = True
    if 'recessive' == peneFunc:
        pene_func = recessive(para)
    elif 'additive' == peneFunc:
        pene_func = additive(para)
    elif 'custom' == peneFunc:
        pene_func = custom(para)
    else:
        pene_func = None
    if 'affectedSibs_GH_TDT' in analyses:
        res['affectedSibs_GH_TDT'] = Sibpair_TDT_gh(pop, N, pene_func, keep_temp=keep_temp)
    if 'affectedSibs_GH_linkage' in analyses:
        res['affectedSibs_GH_linkage'] = Sibpair_LOD_gh(pop, N, pene_func, keep_temp=keep_temp)
    if 'affectedSibs_merlin_linkage' in analyses: 
        res['affectedSibs_merlin_linkage'] = Sibpair_LOD_merlin(pop, N, pene_func, keep_temp=keep_temp)
    if 'caseControl_Association' in analyses: 
        res['caseControl_Association'] = CaseControl_ChiSq(pop, N, pene_func)
    if 'qtraitSibs_merlin_reg' in analyses:
        res['qtraitSibs_merlin_reg'] = QtraitSibs_Reg_merlin(pop, N, customQtrait(para), keep_temp=keep_temp)
    if 'qtraitSibs_merlin_vc' in analyses:
        res['qtraitSibs_merlin_vc'] = QtraitSibs_VC_merlin(pop, N, customQtrait(para), keep_temp=keep_temp)
    if 'largePeds_merlin_reg' in analyses: 
        res['largePeds_merlin_reg'] = LargePeds_Reg_merlin(pop, N, customQtrait(para), keep_temp=keep_temp)
    if 'largePeds_merlin_vc' in analyses:
        res['largePeds_merlin_vc'] = LargePeds_VC_merlin(pop, N, customQtrait(para), keep_temp=keep_temp)
    # calculate population statistics like prevalence
    res.update(popStat(pop))
    print
    print "Writing results to file %s/res.py" % outputDir
    resFile = open(os.path.join(outputDir, 'res.py'), 'w')
    print >> resFile, "# analysis of population %s, at %s" % (dataset, time.asctime())
    print >> resFile
    for key in outputVars.keys():
        if res.has_key(key):
            # description
            print >> resFile, '# %s' % outputVars[key]
            if type(res[key]) == type(''):
                print >> resFile, '%s = "%s"' % (key, res[key])
            else:
                print >> resFile, '%s = %s' % (key, str(res[key]))
            print >> resFile
    resFile.close()
    print 'Done'

