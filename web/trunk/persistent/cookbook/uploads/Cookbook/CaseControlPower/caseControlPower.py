#!/usr/bin/env python
'''
This program is a power calculator for case-control association studies. Compared
to other power calculators, this program is unique in that it assumes that we
know the family history of cases and/or contorls so although only samples are
genotyped and used for statistical analysis, this study design has higher power
than regular case control studies because the probands have higher probability
to have disease alleles at the disease causing locus.

Cases and controls should be specified as a dictionary of cases or controls with
different types. The pedigree types are  presented as a string of affection status
of father, mother, proband, and optional siblings of the proband. The affection
status is specified as 'A' for affected, 'U' for unaffected', and '*' for unknown.
For example, to calculate the power of two mixing pedigree types, you can use 
  cases = "{'**A': 100, '**AA': 200}"
which has 100 regular cases and 200 cases with an affected sibling. Of course, 
regular controls has the third letter as 'U'.

This script can carry out the following functions:

1. Power calculation:
Given the number of cases of each pedigree type and genotype relative risk,
calculate statistical power. You should specify number of cases and controls.

2. Calculate minimal detectable relative risk.
Given number of cases and controls and power, calculate minimal detectable
relative risk. 

3. Calculate needed samples from ratio between ctrols / cases
Given power and relative risk, calculate number of cases and controls from 
a ratio.

4. Calculate needed cased from fixed controls.
Given power and relative risk, calculate number of cases and controls from 
a fixed number of controls.

4. Calculate needed cased from fixed cases.
Given power and relative risk, calculate number of cases and controls from 
a fixed number of cases.

This program allows the users to variate the following passed-in parameters:

K: Disease prevalence
p: Disease allele frequency
x: Marker allele frequency
LD: LD(D')
alpha: Significant level
power: Power
Grr: Genotype relative risk

so analyses on multiple disease models can be performed all at once.

'''

import simuOpt, os, sys, types, time, itertools
simuOpt.setOptions(quiet=True)
from simuPOP.gsl import gsl_cdf_ugaussian_Pinv, gsl_cdf_ugaussian_P
from math import ceil, sqrt

options = [
    {'separator': 'Disease model'},
    {'longarg':'mode=',
     'default':'Additive',
     'label':'Mode of inheritance',
     'chooseOneOf':['Dominant', 'Recessive', 'Additive', 'Multiplicative', 'Additive_AA', 'LogAdditive'],
     'validate':simuOpt.valueOneOf(['Dominant', 'Recessive', 'Additive', 'Multiplicative', 'Additive_AA', 'LogAdditive']),
     'description': '''Disease models with different relationship between risks of AA, AB and BB,
        and different definitions of genotype relative risk. Let A be the disease allele, and
        r_AA, r_AB and r_BB be the penetrance of genotypes AA, AB and BB respectively, six disease
        types are provided for this analysis. Additive_AA and LogAdditive are alternatives of
        types Additive and Multiplicative by using alternative definitions of genotype relative risks.
        | Dominant:  r_AB=r_AA, Grr = r_AA / r_BB
        | Recessive: r_AB=r_BB, Grr = r_AA / r_BB
        | Multiplicative: r_AB=sqrt(r_AA * r_BB), Grr = r_AB / r_BB
        | Additive:  r_AB=(r_AA + r_BB)/2, Grr = r_AB / r_BB
        | Additive_AA: r_AB=(r_AA + r_BB)/2, Grr = r_AA / r_BB
        | LogAdditive: r_AB=sqrt(r_AA * r_BB), Grr = r_AA / r_BB
        ''',
    },
    {'longarg':'K=',
     'default':0.05,
     'label':'Disease prevalence',
     'allowedTypes':[types.ListType, types.TupleType],
     'validate':simuOpt.valueListOf(simuOpt.valueBetween(0., 1)),                   
    },
    {'longarg':'p=',
     'default':0.15,
     'label':'Disease allele frequency',
     'allowedTypes':[types.ListType, types.TupleType],
     'validate':simuOpt.valueListOf(simuOpt.valueBetween(0., 1)),
    },
    {'separator': 'Marker locus information'},
    {'longarg':'x=',
     'default': 0,
     'label': "Marker allele frequency if D' != 1",
     'allowedTypes':[types.ListType, types.TupleType],
     'validate':simuOpt.valueListOf(simuOpt.valueBetween(0., 1)),
     'description': '''Marker allele frequency. If LD is 1, this parameter is ignored
        because it will be assumed to be the same as disease allele frequency.''',
    },
    {'longarg':'genoP=',
     'default':0,
     'allowedTypes':[types.ListType, types.TupleType],
     'validate':simuOpt.valueListOf(simuOpt.valueBetween(0., 1.)),
     'description': '''This option is hidden. If specified, it overrides parameter
        p and specifies the frequency of causal genotypes so that
        | p=sqrt(genoP) for a recessive model (frequency of AA).
        | p=1-sqrt(1-genoP) for other models (frequency of AA + AB)
        ''',
    },
    {'longarg':'genoX=',
     'default':0,
     'allowedTypes':[types.ListType, types.TupleType],
     'validate':simuOpt.valueListOf(simuOpt.valueBetween(0., 1.)),
     'description': '''This option is hidden. If specified, it overrides paramter
        x and specifies the frequency of causal genotypes at the marker locus
        so that
        | p=sqrt(genoX) for a recessive model (frequency of AA).
        | p=1-sqrt(1-genoX) for other models (frequency of AA + AB)
        ''',
    },
    {'longarg':'Dprime=',
     'default': 1,
     'label': "LD (D')",
     'allowedTypes':[types.FloatType, types.IntType],
     'validate':simuOpt.valueBetween(-1., 1.),
     'description': '''Linkage disequilibrium (D') between marker locus and disease locus.
        Only one of parameters Dprime and R2 should be specified.''',
    },
    {'longarg':'R2=',
     'default': 1,
     'label': "LD (R2)",
     'allowedTypes':[types.FloatType, types.IntType],
     'validate':simuOpt.valueBetween(-1., 1.),
     'description': '''Linkage disequilibrium (R2) between marker locus and disease locus.
        Only one of the parameters Dprime and R2 should be specified.''',
    },
    {'separator': 'Type of analysis'},
    {'longarg':'analysis=',
     'default':'Statistical power',
     'label': 'Type of analysis',
     'allowedTypes': [types.StringType],
     'chooseOneOf': ['Statistical power', 'Minimal detectable relative risk', 'Sample size from #ctrl/#case',
        'Sample size from fixed #ctrl', 'Sample size from fixed #case'],
     'description': '''Perform different types of calculations, each require the specification of
        different parameters. More specifically,
        |Statistical power: fix number of cases and controls, alpha, Grr,
            calculate power. Ignore ratio.
        |Minimal detectable relative risk: fix number of cases and controls,
            alpha, power, calculate Grr. Ignore ratio.
        |Sample size from #ctrl/#case: number of cases and controls treats
            as weights, fix ratio, alpha, Grr, power, calculate sample size.
        |Sample size from fixed #ctrl: Fix number of controls, alpha, Grr, power
            and calculate number of cases. Ignore ratio. Numbers of different
            types of cases are treated as weights.
        |Sample size from fixed #cases: Fix number of cases, alpha, Grr, power
            and calculate number of controls. Ignore ratio. Number of different
            types of controls are treated as weights.
        ''',
    },
    {'longarg':'cases=',
     'default': {'**A': 1000},
     'label':'#cases',
     'allowedTypes':[types.DictType, types.NoneType],
     'description': '''Type and number of cases. The count will be considered as 
        weight in 'sample size from #ctrl/#case', and 'sample size from fixed #ctrl'
        analysis'''
    },
    {'longarg':'controls=',
     'default': {'**U': 1000},
     'label': '#controls',
     'allowedTypes':[types.DictType, types.NoneType],
     'description': '''Type and number of controls. The count will be considered 
        as weight in 'Sample size from #ctrl/#case'. '''
    },
    {'longarg':'ratio=',
     'default':1,
     'label':'#controls / #cases for analysis 3',
     'allowedTypes':[types.FloatType, types.IntType],
     'validate':simuOpt.valueGE(0),
     'description': '''Ratio, use for calculate sample size from #ctrl/#case''',
    },
    {'longarg':'alpha=',
     'default':1e-7,
     'label':'Significant level',
     'allowedTypes':[types.ListType, types.TupleType],
     'validate':simuOpt.valueListOf(simuOpt.valueBetween(0., 1)),
    },
    {'longarg':'power=',
     'default': 0.8,
     'label':'Power',
     'allowedTypes':[types.ListType, types.TupleType],
     'validate':simuOpt.valueListOf(simuOpt.valueBetween(0., 1)),
     'description': '''Power of a case-control assocition test if number cases
        and relative risk are specified''',
    },
    {'longarg':'Grr=',
     'default': 1.2,
     'label':'Genotype Relative risk',
     'allowedTypes':[types.ListType, types.TupleType],
     'validate':simuOpt.valueGE(1),
     'description': '''Relative genotype risk for a disease model. The exact definition
        of this quantity depends on the disease model.''',
     },
]


class pedProbabilities:
    '''
    1. Calculate individual probability of being affected conditioned on affection
    status of other relatives in a given pedigree.
    2. Calculate an individual's genotypic probability given affection status of
    a pedigree.
    '''
    def __init__(self, pedType, p, r, x, Dprime=1, R2=1, *args, **kwargs):
        '''
        Parameters description:
        pedType:    pedigree type, which is a sequence of affection status of father
                    mother, proband, and its siblings.
        p:          allele frquency for 'A' (high risk allele)
        r:          penetrances for genotypes 'AA', 'AB', 'BB'
        x:          marker allele frequency
        Dprime:     LD (D') between marker and disease predisposing locus
        R2:         LD (R2) between marker and disease predisposing locus
        '''
        self.pedType = pedType
        self.p = p
        self.r = r
        # probabilities of all possible genotype.
        self.genoProbs = {}
        self._getWholeGenoProbs()
        self.x = x
        self.Dprime = Dprime
        self.R2 = R2
        if self.Dprime == 1.0 and self.R2 == 1.0:
            self.x = self.p
        self._getHaploFreq()

    def _getHaploFreq(self):
        '''
        Calculate hyplotype frequencies given LD, marker allele frequency
        and disease allele frequency

        D = LD / D_max

        ** STEP 2 in documentation. **
        '''
        # calculate haplotype frequencies of AX, AY, BX and BY
        self.h = {}
        if self.Dprime < 1.0:
            D = self.LD * min(self.p * ( 1- self.x), (1 - self.p) * self.x)
        else:
            D = sqrt(self.R2 * self.p * (1 - self.p) * self.x * (1 - self.x))
        self.h['AX'] = D + self.x * self.p
        # check if self.h[0] - xp is positive?
        if self.h['AX'] - self.x * self.p <= 0:
            raise ValueError('h11-xp can NOT be negative, adjust p and x to restrict it to be positive')
        self.h['AY'] = self.p - self.h['AX']
        self.h['BX'] = self.x - self.h['AX']
        self.h['BY'] = (1 - self.p) - self.h['BX']
        # to avoid round error
        for k in self.h.keys():
            if self.h[k] > -1e-9 and self.h[k] <0:
                self.h[k] = 0.

    def PrGeno(self, geno):
        '''
        Calculate the probability of genotype 'AA', 'AB' or 'BB'
        '''
        if geno == 'AA':
            return self.p ** 2
        elif geno == 'AB':
            return 2 * self.p * (1 - self.p)
        elif geno == 'BB':
            return (1 - self.p) ** 2
        else:
            # this should not happen
            assert 0

    def _getWholeGenoProbs(self):
        '''
        Compute probability of all possible genotypic structures given the
        number of offspring in self.pedType.

        ** Part of STEP 3 in the documentation **
        '''
        for key in itertools.product(* tuple([('AA', 'AB', 'BB') for x in self.pedType])):
            if key[:2] == ('AA', 'AA'):
                pOffspring = {'AA':1, 'AB':0, 'BB':0}
            elif key[:2] in [('AA', 'AB'), ('AB', 'AA')]:
                pOffspring = {'AA':0.5, 'AB':0.5, 'BB':0}
            elif key[:2] in [('AA', 'BB'), ('BB', 'AA')]:
                pOffspring = {'AA':0, 'AB':1, 'BB':0}
            elif key[:2] == ('AB', 'AB'):
                pOffspring = {'AA':0.25, 'AB':0.5, 'BB':0.25}
            elif key[:2] in [('AB', 'BB'), ('BB', 'AB')]:
                pOffspring = {'AA':0, 'AB':0.5, 'BB':0.5}
            elif key[:2] == ('BB', 'BB'):
                pOffspring = {'AA':0, 'AB':0, 'BB':1}
            else:
                # should have listed all cases
                assert 0
            self.genoProbs[key] = self.PrGeno(key[0]) * self.PrGeno(key[1])
            for off in key[2:]:
                self.genoProbs[key] *= pOffspring[off]
    
    def PrPedGivenGeno(self, geno):
        '''
          Pr(affection status of all memebers | genotype at proband)
        = Pr(affection status of all members with genotype at prob) / Pr(geno)

        ** Part of STEP 3 in the documentation **
        '''
        TotProb = 0
        for key in itertools.product(* tuple([('AA', 'AB', 'BB') for x in self.pedType])):
            # proband with specified genotype
            if key[2] == geno:
                Prob = self.genoProbs[key]
                for ind,status in enumerate(self.pedType):
                    if status == 'A':
                        Prob *= self.r[key[ind]]
                    elif status == 'U':
                        Prob *= 1 - self.r[key[ind]]
                    # the '*' case does not count
                TotProb += Prob
        assert TotProb >= 0
        return TotProb / self.PrGeno(geno)

    def PrPed(self):
        '''
        Pr(affection status of all members) = sum (aff | geno all) Pr(geno all)

        ** Part of STEP 4 in the documentation **
        '''
        TotProb = 0
        for key in itertools.product(* tuple([('AA', 'AB', 'BB') for x in self.pedType])):
            Prob = self.genoProbs[key]
            for ind,status in enumerate(self.pedType):
                if status == 'A':
                    Prob *= self.r[key[ind]]
                elif status == 'U':
                    Prob *= 1 - self.r[key[ind]]
                # the '*' case does not count
            TotProb += Prob
            assert Prob >= 0
        assert TotProb >= 0
        return TotProb

    def PrGenoGivenPed(self, geno):
        '''
          Pr(genotype at proband | affection status of all members)
        = Pr(ped | geno) * Pr(geno) / Pr(ped)

        ** Part of STEP 4 in the documentation **
        '''
        return self.PrPedGivenGeno(geno) * self.PrGeno(geno) / self.PrPed()

    def PrMarkerGivenPed(self, geno):
        '''
          Pr(genotype at marker | affection status)
        = combination of Pr(genotype at proband | affection status)

        ** STEP 5 in the documentation **
        '''
        if geno == 'XX':
            return self.PrGenoGivenPed('AA')*(self.h['AX']**2)/self.PrGeno('AA') + \
                2 * self.PrGenoGivenPed('AB')*self.h['AX']*self.h['BX']/self.PrGeno('AB') + \
                self.PrGenoGivenPed('BB')*(self.h['BX']**2)/self.PrGeno('BB')
        elif geno == 'XY':
            return 2 * self.PrGenoGivenPed('AA')*self.h['AX']*self.h['AY']/self.PrGeno('AA') + \
                2 * self.PrGenoGivenPed('AB')*(self.h['AX']*self.h['BY'] + self.h['BX']*self.h['AY'])/self.PrGeno('AB') + \
                2 * self.PrGenoGivenPed('BB')*self.h['BX']*self.h['BY']/self.PrGeno('BB')
        elif geno == 'YY':
            return self.PrGenoGivenPed('AA')*(self.h['AY']**2)/self.PrGeno('AA') + \
                2 * self.PrGenoGivenPed('AB')*self.h['AY']*self.h['BY']/self.PrGeno('AB') + \
                self.PrGenoGivenPed('BB')*(self.h['BY']**2)/self.PrGeno('BB')
        else:
            assert 0
      
    def PrAlleleGivenPed(self, allele):
        '''
          Pr(allele | affection status)
        = Pr(XX | aff) + Pr(XY | aff) / 2

        ** STEP 6 in the documentation **
        '''
        if allele == 'X':
            return self.PrMarkerGivenPed('XX') + self.PrMarkerGivenPed('XY') / 2.
        else:
            return self.PrMarkerGivenPed('YY') + self.PrMarkerGivenPed('XY') / 2.


def RelativeRiskToPenetrance(K, p, mode, Grr):
    '''Convert a disease model to Penetrance. The reverse will be needed when calculating
    minimal detectable relative risk
    ** STEP 1 in documentation. **
    '''
    if Grr < 1:
        raise ValueError('Can not calculate power with genotype relative risk < 1')
    q = 1 - p
    r = {}
    if mode == 'Recessive':
        r['BB'] = r['AB'] = K/((p**2)*Grr + 2*p*q + q**2)
        r['AA'] = Grr * r['BB']
    elif mode == 'Additive':
        r['BB'] = K/((p**2)*(2*Grr - 1) + 2*p*q*Grr + q**2)
        r['AB'] = Grr * r['BB']
        r['AA'] = 2. * r['AB'] - r['BB']
    elif mode == 'Dominant':
        r['BB'] = K/((p**2)*Grr + 2*p*q*Grr + q**2)
        r['AA'] = r['AB'] = Grr * r['BB']
    elif mode == 'Multiplicative':
        r['BB'] = K/((p**2)*Grr*Grr + 2*p*q*Grr + q**2)
        r['AB'] = Grr*r['BB']
        r['AA'] = r['AB']*r['AB']/r['BB']
    elif mode == 'Additive_AA':
        r['BB'] = K/((p**2)*Grr + p*q*(Grr+1) + q**2)
        r['AA'] = Grr * r['BB']
        r['AB'] = (r['BB'] + r['AA']) / 2.
    elif mode == 'LogAdditive':
        r['BB'] = K/((p**2)*Grr + 2*p*q*sqrt(Grr) + q**2)
        r['AB'] = sqrt(Grr)*r['BB']
        r['AA'] = r['AB']*r['AB']/r['BB']
    else:
        raise ValueError('Non-supported mode: ' + mode)
    for k in r.keys():
        if r[k] > 1.:
            # print 'Warning: disease model needs to higher than 1 penetrance. Assuming 1.'
            r[k] = 1.
    assert r['AA'] >= 0 and r['AA'] <= 1
    assert r['AB'] >= 0 and r['AB'] <= 1
    assert r['BB'] >= 0 and r['BB'] <= 1
    assert r['AA'] >= r['AB'] and r['AB'] >= r['BB']
    return r
 

class powerCalculator:
    '''
    Compute the statistical power or sample size and detectable Relative Risk
    for Case-Control association studies.
    '''
    def __init__(self, mode, K, p=0, x=0, Dprime=1, R2=1, genoP=0, genoX=0, logger=None):
        '''
        Parameters description:
        mode:    Mode of inheritance (dominant, recessive, additive, and multiplicative)
        K:       Prevalence
        p:       Disease allele frequency at risk at the disease locus
        x:       Marker allele frequency (allele X)
        genoP:   If specifies, frequency of causal genotype (p will be ignored)
        genoX:   If specifies, frequency of causal marker (x will be ignored)
        LD:      Linkage disequilibrium D'
        logger:  If a logging object is given, it will be used to output more detailed information.
        '''
        self.mode = mode
        self.K = K
        self.Dprime = Dprime
        self.R2 = R2
        self.p = p
        self.x = x
        if genoP != 0:
            if self.mode == 'Recessive':
                self.p = sqrt(genoP)
            else:
                self.p = 1 - sqrt(1 - genoP)
        if genoX != 0:
            if self.mode == 'Recessive':
                self.x = sqrt(genoX)
            else:
                self.x = 1 - sqrt(1 - genoX)
        self.r = {}
        self.logger = logger
        if logger:
            print 'P:', self.p, ' X:', self.x


    def expectedFreq(self, samples):
        '''
        Calculate expected frequency of specified pedigree structure and counts.

        samples should be a dictionary of pedigree types with values as weights.
        '''
        allP = 0
        for ped,cnt in samples.iteritems():
            # change aa to AA.
            probs = pedProbabilities(ped.upper(), self.p, self.r, self.x, self.Dprime, self.R2)
            p = probs.PrAlleleGivenPed('X')
            assert p >= 0
            if self.logger:
                self.logger.info('Pedigree %s has Pr(X): %.4f' % (ped, p))
            if cnt == 0:
                print 'Warning: Pedigree %s has zero weight' % ped
            allP += p * cnt
        allP /= sum(samples.values())
        if self.logger:
            self.logger.info('Overall Pr(X): %.4f' % allP)
        return allP


    def getPower(self, Grr, cases, controls, alpha):
        '''
        Compute the statistical power given cases
        cases:     a dictionary of number of cases for each type of pedigree.
        controls:  a dictionary of number of controls for each type of pedigree.
        '''
        # calculate relative risk
        self.r = RelativeRiskToPenetrance(self.K, self.p, self.mode, Grr)
        if self.logger:
            self.logger.info('Penetrance: ' + str(self.r))
        #
        # p1: for controls
        if len(controls) == 0:
            raise ValueError('Number of controls should be given to calculate statistical power')
        if len(cases) == 0:
            raise ValueError('Number of cases should be given to calculate statistical power')
        #
        p1 = self.expectedFreq(controls)
        N1 = sum(controls.values())
        p2 = self.expectedFreq(cases)
        N2 = sum(cases.values())
        if N1 == 0:
            raise ValueError('Number of controls should be given to calculate statistical power')
        if N2 == 0:
            raise ValueError('Number of cases should be given to calculate statistical power')
        #
        ratio = float(N1) / N2
        #
        #p0 = (p2 + ratio * p1) / (1 + ratio)
        #sigmaH0 = (p0 * (1 - p0) * (1 + 1/ratio))**0.5
        #
        sigmaH1 = (p2 * (1 - p2) + p1 * (1 - p1) / ratio)**0.5
        # negative sign is canceled because gsl_cdf(0.05/2) is 0.025, not 0.975.
        v1 = (gsl_cdf_ugaussian_Pinv(alpha/2.)*sigmaH1 + ((2*N2)**0.5)*(p2-p1)) / sigmaH1
        v2 = (gsl_cdf_ugaussian_Pinv(alpha/2.)*sigmaH1 - ((2*N2)**0.5)*(p2-p1)) / sigmaH1
        power = (gsl_cdf_ugaussian_P(v1)) + (gsl_cdf_ugaussian_P(v2))
        if self.logger:
            self.logger.info('p ctrl: %.4f p case: %.4f sigma_H1: %.7f statistic: %.4f power: %.4f' \
                    % (p1, p2, sigmaH1/((2*N2)**0.5), ((2*N2)**0.5)*(p1-p2)/sigmaH1, power))
        return power

    def getRelativeRisk(self, cases, controls, alpha, power):
        '''
        Using the bisection method to calculate and return a Detectable Relative Risk(DDR)
        '''
        # find a possible range of Grr(minRR, maxRR)
        minGrr = 1
        maxGrr = 10
        while maxGrr < 10000:
            if self.getPower(maxGrr, cases, controls, alpha) < power:
                maxGrr *= 2
            else:
                break
        Grr = 1
        while True:
            if maxGrr - minGrr < 1e-7:
                return Grr
            p = self.getPower(Grr, cases, controls, alpha)
            if p > power:
                maxGrr = Grr
                Grr = (minGrr + maxGrr) / 2.
            elif p <= power:
                minGrr = Grr
                Grr = (minGrr + maxGrr) / 2.
            if self.logger:
                self.logger.info('Grr: %.4f (min: %.4f, max: %.4f), Power: %.4f' % (Grr, minGrr, maxGrr, p))
        return 0


    def getSampleSizeFromRatio(self, Grr, cases, controls, ratio, alpha, power):
        '''
        Return required sample size given power, ratio of number of controls
        to number of cases and pedigrees from where cases will be drawn. 
        '''
        # calculate relative risk
        self.r = RelativeRiskToPenetrance(self.K, self.p, self.mode, Grr)
        if self.logger:
            self.logger.info('Penetrance: ' + str(self.r))
        # get expected p from cases and controls
        p1 = self.expectedFreq(controls)
        p2 = self.expectedFreq(cases)
        #p0 = (p2 + ratio*p1) / (1 + ratio)
        #sigmaH0 = (p0 * (1 - p0) * (1 + 1/ratio))**0.5
        sigmaH1 = (p2 * (1 - p2) + p1 * (1 - p1) / ratio)**0.5
        # return number of cases
        N2 = (((gsl_cdf_ugaussian_Pinv(alpha/2)*sigmaH1 + 
                gsl_cdf_ugaussian_Pinv(1-power)*sigmaH1)/(p2-p1))**2)/2
        N1 = ratio * N2
        ret = {}
        caseCnt = sum(cases.values())
        for ped,cnt in cases.iteritems():
            ret[ped] = int(ceil(N2 * cnt / caseCnt))
        ctrlCnt = sum(controls.values())
        for ped,cnt in controls.iteritems():
            ret[ped] = int(ceil(N1 * cnt / ctrlCnt))
        return ret
    
    def getSampleSizeFromControls(self, Grr, cases, controls, alpha, power):
        '''
        Return required sample size given power and number of controls.
        '''
        # calculate relative risk
        #
        if sum(controls.values()) == 0:
            raise ValueError("# of controls must be specified for this analysis")
        caseCnt = sum(cases.values())
        if caseCnt == 0:
            raise ValueError('Weight of controls are zero')
        pp = 0
        N2_lower = 1
        N2_upper = 1
        # find upper bound
        while True:
            tmp_cases = {}
            for key,cnt in cases.iteritems():
                tmp_cases[key] = int(ceil(float(N2_upper) * cnt / caseCnt))
            pp = self.getPower(Grr, tmp_cases, controls, alpha)
            if pp < power:
                N2_lower = N2_upper
                N2_upper *= 2
                if N2_upper >= 1e10:
                    raise ValueError('Specified power is not achievable. Please adjust disease model or number of controls')
            else:
                break
        # bisection method
        N2 = N2_lower
        while True:
            if N2_upper - N2_lower <= 1:
                break
            tmp_cases = {}
            for key,cnt in cases.iteritems():
                tmp_cases[key] = int(ceil(float(N2) * cnt / caseCnt))
            pp = self.getPower(Grr, tmp_cases, controls, alpha)
            if pp > power:
                N2_upper = N2
                N2 = (N2_lower + N2_upper) / 2
            elif pp <= power:
                N2_lower = N2
                N2 = (N2_lower + N2_upper) / 2
            #print 'Cases %s has power %.2f (low: %d, high: %d)' % (tmp_cases, pp, N2_lower, N2_upper)
        #
        ret = {}
        # known controls
        ret.update(controls)
        # cases
        for ped,cnt in cases.iteritems():
            ret[ped] = int(ceil(N2) * cnt / caseCnt)
        return ret

    def getSampleSizeFromCases(self, Grr, cases, controls, alpha, power):
        '''
        Return required sample size given power and number of controls.
        '''
        # calculate relative risk
        #
        if sum(cases.values()) == 0:
            raise ValueError("# of cases must be specified for this analysis")
        ctrlCnt = sum(controls.values())
        if ctrlCnt == 0:
            raise ValueError('Weight of controls are zero')
        pp = 0
        N1_lower = 1
        N1_upper = 1
        # find upper bound
        while True:
            tmp_controls = {}
            for key,cnt in controls.iteritems():
                tmp_controls[key] = int(ceil(float(N1_upper) * cnt / ctrlCnt))
            pp = self.getPower(Grr, cases, tmp_controls, alpha)
            if pp < power:
                N1_lower = N1_upper
                N1_upper *= 2
                if N1_upper >= 1e10:
                    raise ValueError('Specified power is not achievable. Please adjust disease model or number of controls')
            else:
                break
        # bisection method
        N1 = N1_lower
        while True:
            if N1_upper - N1_lower <= 1:
                break
            tmp_controls = {}
            for key,cnt in controls.iteritems():
                tmp_controls[key] = int(ceil(float(N1) * cnt / ctrlCnt))
            pp = self.getPower(Grr, cases, tmp_controls, alpha)
            if pp > power:
                N1_upper = N1
                N1 = (N1_lower + N1_upper) / 2
            elif pp <= power:
                N1_lower = N1
                N1 = (N1_lower + N1_upper) / 2
            #print 'Controls %s has power %.2f (low: %d, high: %d)' % (tmp_cases, pp, N1_lower, N1_upper)
        #
        ret = {}
        # known cases
        ret.update(cases)
        # controls
        for ped,cnt in controls.iteritems():
            ret[ped] = int(ceil(N1) * cnt / ctrlCnt)
        return ret

if __name__ == '__main__':
    # get all parameters
    pars = simuOpt.Params(options, 'A power calculator for case control association studies with know family histories',
        __doc__)
    if not pars.getParam():
        sys.exit(0)
    if pars.analysis == 'Statistical power':
        # calculate power
        print "K\tp (or genoype p)\tx (or genotype x)\tD' (or R2)\talpha\tGrr\tPower"
        for K, p, x, genoP, genoX, alpha, Grr in itertools.product(
                pars.K, pars.p, pars.x, pars.genoP, pars.genoX, 
                pars.alpha, pars.Grr):
            print '%.3f\t%.3f\t%.3f\t%.2f\t%.8f\t%.2f\t%f' % \
                (K, p if genoP == 0. else genoP, x if genoX == 0. else genoX,
                        pars.Dprime if pars.Dprime == 1 else pars.R2, alpha, Grr, 
                powerCalculator(mode=pars.mode, K=K, p=p, x=x,
                    Dprime=pars.Dprime, R2=pars.R2, genoP=genoP, genoX=genoX).getPower(
                Grr=Grr, cases=pars.cases, controls=pars.controls, alpha=alpha))
    elif pars.analysis == 'Minimal detectable relative risk':
        print "K\tp (or genoype p)\tx (or genotype x)\tD' (or R2)\talpha\tPower\tGrr"
        for K, p, x, genoP, genoX, alpha, power in itertools.product(
                pars.K, pars.p, pars.x, pars.genoP, pars.genoX, 
                pars.alpha, pars.power):
            print '%.3f\t%.3f\t%.3f\t%.2f\t%.8f\t%.2f\t%f' % \
                (K, p if genoP == 0. else genoP, x if genoX == 0. else genoX,
                        pars.Dprime if pars.Dprime == 1 else pars.R2, alpha, power,
                powerCalculator(mode=pars.mode, K=K, p =p, x=x,
                    Dprime=pars.Dprime, R2=pars.R2, genoP=genoP, genoX=genoX).getRelativeRisk(
                    cases=pars.cases, controls=pars.controls, alpha=alpha, power=power))
    elif pars.analysis == 'Sample size from #ctrl/#case':
        print "K\tp (or genoype p)\tx (or genotype x)\tD' (or R2)\talpha\tPower\tGrr\tSampleSize"
        for K, p, x, genoP, genoX, alpha, power, Grr in itertools.product(
                pars.K, pars.p, pars.x, pars.genoP, pars.genoX, 
                pars.alpha, pars.power, pars.Grr):
            print '%.3f\t%.3f\t%.3f\t%.2f\t%.8f\t%.2f\t%f\t%s' % \
                (K, p if genoP == 0. else genoP, x if genoX == 0. else genoX,
                        pars.Dprime if pars.Dprime == 1 else pars.R2, alpha, power, Grr,
                powerCalculator(mode=pars.mode, K=K, p =p, x=x,
                    Dprime=pars.Dprime, R2=pars.R2, genoP=genoP, genoX=genoX).getSampleSizeFromRatio(
                    Grr=Grr, cases=pars.cases, controls=pars.controls, ratio=pars.ratio,
                    alpha=alpha, power=power))
    elif pars.analysis == 'Sample size from fixed #ctrl':
        print "K\tp (or genoype p)\tx (or genotype x)\tD' (or R2)\talpha\tPower\tGrr\tSampleSize"
        for K, p, x, genoP, genoX, alpha, power, Grr in itertools.product(
                pars.K, pars.p, pars.x, pars.genoP, pars.genoX, 
                pars.alpha, pars.power, pars.Grr):
            print '%.3f\t%.3f\t%.3f\t%.2f\t%.8f\t%.2f\t%f\t%s' % \
                (K, p if genoP == 0. else genoP, x if genoX == 0. else genoX,
                        pars.Dprime if pars.Dprime == 1 else pars.R2, alpha, power, Grr,
                powerCalculator(mode=pars.mode, K=K, p =p, x=x,
                    Dprime=pars.Dprime, R2=pars.R2, genoP=genoP, genoX=genoX).getSampleSizeFromControls(
                    Grr=Grr, cases=pars.cases, controls=pars.controls, alpha =alpha, power =power))
    elif pars.analysis == 'Sample size from fixed #case':
        print "K\tp (or genoype p)\tx (or genotype x)\tD' (or R2)\talpha\tPower\tGrr\tSampleSize"
        for K, p, x, genoP, genoX, alpha, power, Grr in itertools.product(
                pars.K, pars.p, pars.x, pars.genoP, pars.genoX, 
                pars.alpha, pars.power, pars.Grr):
            print '%.3f\t%.3f\t%.3f\t%.2f\t%.8f\t%.2f\t%f\t%s' % \
                (K, p if genoP == 0. else genoP, x if genoX == 0 else genoX,
                        pars.Dprime if pars.Dprime == 1 else pars.R2, alpha, power, Grr,
                powerCalculator(mode=pars.mode, K=K, p =p, x=x,
                    Dprime=pars.Dprime, R2=pars.R2, genoP=genoP, genoX=genoX).getSampleSizeFromCases(
                    Grr=Grr, cases=cases, controls=controls, alpha =alpha, power =power))
    else:
        raise ValueError('Wrong analysis type ' + pars.analysis)

