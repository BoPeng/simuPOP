#!/usr/bin/env python
#
# This script is used to perform all analysis and draw figures for Peng et al, 2010.
#
# To use this script, run this script with function name as arguments. For example,
#
# > PowerAnalysis.py Figure1
# > PowerAnalysis.py Table1
#
from caseControlPower import *

from itertools import product
import math

# rpy is used to draw figures 1 and 2.
from rpy import *
# under windows, this option is needed.
if os.name == 'nt':
    r.options(windowsBuffered=False)
    if int(r.R_Version()['svn rev']) >= 48333:
        r.windows_options(buffered=False)


defaultAlpha = 0.05/500000

def Figure1():
    '''Draw figure 1'''
    ratio = 1
    K = 0.1

    all_freq = [0.05, 0.1, 0.25, 0.5]
    all_Grr = [1.2, 1.3, 1.4]
    all_cases = [{'**A': 1000}, {'**AA': 1000}, {'*AAA': 1000}]
    controls = {'**U': 1000}

    power = {}
    sampleSize = {}
    for freq, Grr, cases in product(all_freq, all_Grr, all_cases):
        cal = powerCalculator(mode='Additive', K=K, p=freq, x=freq)
        power[(Grr, freq, cases.keys()[0])] = cal.getPower(Grr=Grr, cases=cases, controls=controls,
            alpha=defaultAlpha)
        sampleSize[(Grr, freq, cases.keys()[0])] = cal.getSampleSizeFromRatio(Grr=Grr, cases=cases,
            controls=controls, ratio=1, power=0.8, alpha=defaultAlpha)[cases.keys()[0]]

    Grr_shift = {all_Grr[0]:0, all_Grr[1]:1, all_Grr[2]:2}
    freq_shift = {all_freq[0]:0, all_freq[1]:0.20, all_freq[2]:0.40, all_freq[3]:0.60}
    ped_shift = {'**A':0, '**AA':0.05, '*AAA':0.10}
    cols = {'**A': 'green', '**AA':'red', '*AAA': 'blue'}

    r.postscript('Fig1.eps', width=8, height=6, horizontal=True)
    #r.pdf('Figure1.pdf', width=8, height=6)
    #r.png('Figure1.png', width=800, height=600)
    r.par(mfrow=[2,1], mar=[3, 4, 1.5, 0], oma=[3,1,1,1])
    r.plot(0, 0, axes=False, main='Statistical power with 1000 cases and 1000 controls', xlab='', ylab='Power (%)',
            xlim=[0, 2.75], ylim=[0, 1.3], type='l', yaxs='i')
    at = []
    labels = []
    for Grr in all_Grr:
        for freq in all_freq:
            for ped in [x.keys()[0] for x in all_cases]:
                at.append(Grr_shift[Grr] + freq_shift[freq] + ped_shift[ped] + 0.05)
                labels.append('%.2f' % freq)

    for Grr in all_Grr:
        for freq in all_freq:
            ypos = [power[(Grr, freq, x.keys()[0])] + 0.03 for x in all_cases]
            if ypos[1] - ypos[0] < 0.06:
                ypos[1] = ypos[0] + 0.06
            if ypos[2] - ypos[1] < 0.06:
                ypos[2] = ypos[1] + 0.06
            for idx,ped in enumerate([x.keys()[0] for x in all_cases]):
                r.rect(Grr_shift[Grr] + freq_shift[freq] + ped_shift[ped],
                            0, Grr_shift[Grr] + freq_shift[freq] + ped_shift[ped] + 0.05,
                            power[(Grr, freq, ped)],
                            col=cols[ped])
                label = '%.1f' % (power[(Grr, freq, ped)]*100)

                r.text(Grr_shift[Grr] + freq_shift[freq] + ped_shift[ped] + 0.045,
                            ypos[idx] + 0.02, label, cex=0.6, adj=1)

    for Grr in all_Grr:
        r.mtext('Allele frequency', side=1, line=2, at = Grr_shift[Grr] + 0.25)
        r.text('relative risk = %.1f' % Grr, x = Grr_shift[Grr] + 0.25, y= 1.25, cex=0.8)

    r.axis(1, lwd_ticks=0, at=at, labels=labels)
    r.axis(2, ['0', '20', '40', '60', '80', '100'], at=[0, 0.2, 0.4, 0.6, 0.8, 1], las=1)

    # sample size
    r.plot(0, 0, axes=False, main='Required sample size to achieve expected power of 80%', xlab='', ylab='Number of cases',
            type='l', xlim=[0, 2.75], ylim=[0, max(sampleSize.values()) + 5000], yaxs='i')

    for Grr in all_Grr:
        for freq in all_freq:
            ypos = [sampleSize[(Grr, freq, x.keys()[0])] + 700 for x in all_cases]
            if ypos[1] - ypos[2] < 1100:
                ypos[1] = ypos[2] + 1100
            if ypos[0] - ypos[1] < 1100:
                ypos[0] = ypos[1] + 1100
            for idx,ped in enumerate([x.keys()[0] for x in all_cases]):
                r.rect(Grr_shift[Grr] + freq_shift[freq] + ped_shift[ped],
                            0, Grr_shift[Grr] + freq_shift[freq] + ped_shift[ped] + 0.05,
                            sampleSize[(Grr, freq, ped)],
                            col=cols[ped])
                r.text(Grr_shift[Grr] + freq_shift[freq] + ped_shift[ped] + 0.001,
                            ypos[idx],
                            '%d' % sampleSize[(Grr, freq, ped)], cex=0.6, adj=0)
    r.axis(1, lwd_ticks=0, at=at, labels=labels, adj=1)
    r.axis(2, at=[0, 10000, 20000], labels=['0', '10k', '20k'], las=1)

    for Grr in all_Grr:
        r.mtext('Allele frequency', side=1, line=2, at = Grr_shift[Grr] + 0.25)
        r.text('relative risk = %.1f' % Grr, x = Grr_shift[Grr] + 0.25, y= 35000, cex=0.8)
    r.mtext('**a', side=1, line=4, at=0.49, adj=0)
    r.mtext('**aa', side=1, line=4, at=1, adj=0)
    r.mtext('*aaa', side=1, line=4, at=1.6, adj=0)
    r.par(family="HersheySerif")
    r(r"mtext('\\#H2327', side=1, line=4.3, at=0.45, cex=2, col='green', adj=1)")
    r(r"mtext('\\#H2327', side=1, line=4.3, at=0.96, cex=2, col='red', adj=1)")
    r(r"mtext('\\#H2327', side=1, line=4.3, at=1.56, cex=2, col='blue', adj=1)")
    r.dev_off()


def caseRatio(ped1, ped2, K, p, mode, Grr, power=0.8, alpha=defaultAlpha):
    '''
    For a particular configuration, calculate power of N case, N control,
    then use N control + ?? cases to achieve the same power, return N/??
    '''
    cal = powerCalculator(mode=mode, K=K, p=p, x=p)
    N = cal.getSampleSizeFromRatio(Grr=Grr, cases={'**A': 1}, controls={'**U': 1}, ratio=1, alpha=alpha, power=power)['**U']
    size = cal.getSampleSizeFromRatio(Grr=Grr, cases={ped1: 1}, controls={ped2: 1}, ratio=1,
        power=0.8, alpha=alpha)
    return N, N * 1.0 / size[ped1]

def caseRatioFixedControl(ped1, ped2, K, p, mode, Grr, power=0.8, alpha=defaultAlpha):
    '''
    For a particular configuration, calculate power of N case, N control,
    then use N control + ?? cases to achieve the same power, return N/??
    '''
    cal = powerCalculator(mode=mode, K=K, p=p, x=p)
    N = cal.getSampleSizeFromRatio(Grr=Grr, cases={'**A': 1}, controls={ped2: 1}, ratio=1, alpha=alpha, power=power)[ped2]
    size = cal.getSampleSizeFromControls(Grr=Grr, cases={ped1: 1}, controls={ped2: N},
        power=0.8, alpha=alpha)
    return N, size[ped1] * 1.0 / N

from itertools import product

# use a global one ...
L_all = {}
def allValues(allGrr, allK, allp, allMode, allPower, allFam, allAlpha):
    L = {}
    for K, p, mode, Grr, power, fam, alpha in product(
            allK,
            allp,
            allMode,
            allGrr, 
            allPower,
            allFam,
            allAlpha):
            key = (Grr, K, p, mode, power, fam, alpha)
            if not L_all.has_key(key):
                L[key] = caseRatio(fam, '**U', K, p, mode, Grr, power, alpha)[1]
            else:
                L[key] = L_all[key]
    return L

def Figure2():
    # typical values
    allGrr =  (1.4, )
    allK =    (0.02, )
    allp =    (0.1,)
    allMode = ('Additive',)
    allPower = (0.8,)
    allFam = ('**AA',)
    allAlpha = (defaultAlpha,)
    #
    r.postscript('Fig2.eps', width=4, height=15, horizontal=False)
    #r.pdf('Figure2.pdf', width=8, height=6)
    r.par(mfrow=[5,1], mar=[4, 4, 1.5, 1], oma=[2,1,1,1])
    # Lambda vs. # of affected siblings
    myFam = ('**AA', '**AAA', '**AAAA', '**AAAAA')
    myK = (0.1, 0.02)
    myGrr = (1.4, 1.8)
    L = allValues(myGrr, myK, allp, allMode, allPower, myFam, allAlpha)
    r.plot(0, 0, axes=False, main='', xlab='Number of affected siblings', ylab=r('expression(lambda)'),
        xlim=[1, 4], ylim=[0, 12], type='n', yaxs='i')
    r.axis(1, [1, 2, 3, 4], [x.lower() for x in myFam], adj=0.5)
    r.axis(2)
    cols = ['blue', 'red']
    pchs = [2, 3]
    for Grr, K, p in product(myGrr, myK, allp):
        L = allValues([Grr], [K], [p], allMode, allPower, myFam, allAlpha)
        r.lines(range(1, 5),
                [L[(Grr, K, p, allMode[0], allPower[0], fam, allAlpha[0])] for fam in myFam],
                col=cols[myGrr.index(Grr)],
                pch=pchs[myK.index(K)],
                lty=1,
                type='b')
    r.legend(x=1.0, y=12, legend=['Grr=1.4', 'Grr=1.8'], col=cols, lty=1, bty='n')
    r.legend(x=2.2, y=12, legend=['K=0.1', 'K=0.02'], pch=pchs, bty='n')
    r.mtext('a', 3, line=0, adj=1.05)
    #
    # Lambda vs. Grr differet by fam
    myGrr = [1 + x/10. for x in range(1,10)]
    myFam = ('**AA', '**AAA', '**AAAA')
    myK = (0.1, 0.02)
    L = allValues(myGrr, myK, allp, allMode, allPower, myFam, allAlpha)
    r.plot(0, 0, axes=False, main='', xlab='Genotype relative risk', ylab=r('expression(lambda)'),
        xlim=[1, 2], ylim=[0, 12], type='n', yaxs='i')
    r.axis(1, myGrr, [str(x) for x in myGrr])
    r.axis(2)
    cols = ['blue', 'red', 'green']
    pchs = [2, 3]
    for K, fam in product(myK, myFam):
        L = allValues(myGrr, [K], allp, allMode, allPower, [fam], allAlpha)
        r.lines(myGrr,
                [L[(Grr, K, allp[0], allMode[0], allPower[0], fam, allAlpha[0])] for Grr in myGrr],
                col=cols[myFam.index(fam)],
                pch=pchs[myK.index(K)],
                lty=1,
                type='b')
    r.legend(x=1.0, y=12, legend=['**aa', '**aaa', '**aaaa'], col=cols, lty=1, bty='n')
    r.legend(x=1.4, y=12, legend=['K=0.1', 'K=0.02'], pch=pchs, bty='n')
    r.mtext('b', 3, line=0, adj=1.05)
    #
    # Expected Power
    myPower = [x/10. for x in range(3, 10)]
    myGrr = [1.4, 1.8]
    myFam = ('**AA', '**AAA', '**AAAA')
    L = allValues(myGrr, allK, allp, ['Additive'], myPower, myFam, allAlpha)
    r.plot(0, 0, axes=False, main='', xlab='Expected statistical power (%)', ylab=r('expression(lambda)'),
        xlim=[min(myPower), max(myPower)], ylim=[0, 12], type='n', yaxs='i')
    r.axis(1, myPower, [str(int(x*100)) for x in myPower])
    r.axis(2)
    cols = ['blue', 'red', 'green']
    pchs = [2, 3]
    for Grr, fam in product(myGrr, myFam):
        L = allValues([Grr], allK, allp, allMode, myPower, [fam], allAlpha)
        r.lines(myPower,
                [L[(Grr, allK[0], allp[0], allMode[0], power, fam, allAlpha[0])] for power in myPower],
                col=cols[myFam.index(fam)],
                pch=pchs[myGrr.index(Grr)],
                lty=1,
                type='b')
    r.legend(x=0.3, y=12, legend=['**aa', '**aaa', '**aaaa'], col=cols, lty=1, bty='n')
    r.legend(x=0.55, y=12, legend=['Grr=1.4', 'Grr=1.8'], pch=pchs, bty='n')
    r.mtext('c', 3, line=0, adj=1.05)
    #
    # Alpha
    myAlpha = [1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7]
    myGrr = [1.2, 1.4]
    myFam = ('**AA', '**AAA', '**AAAA')
    L = allValues(myGrr, allK, allp, ['Additive'], allPower, myFam, myAlpha)
    r.plot(myAlpha[0], 0, axes=False, main='', xlab='Significance level', ylab=r('expression(lambda)'),
        xlim=[min(myAlpha), max(myAlpha)], ylim=[0, 12], type='n', yaxs='i', log='x')
    r('axis(1, c(1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7), c(expression(10^-4), expression(3*"x"*10^-5),'
        'expression(10^-5), expression(3*"x"*10^-6), expression(10^-6), expression(3*"x"*10^-7), expression(10^-7)))')
    r.axis(2)
    cols = ['blue', 'red', 'green']
    pchs = [2, 3]
    for Grr, fam in product(myGrr, myFam):
        L = allValues([Grr], allK, allp, allMode, allPower, [fam], myAlpha)
        r.lines(myAlpha,
                [L[(Grr, allK[0], allp[0], allMode[0], allPower[0], fam, alpha)] for alpha in myAlpha],
                col=cols[myFam.index(fam)],
                pch=pchs[myGrr.index(Grr)],
                lty=1,
                type='b')
    r.legend(x=1e-7, y=12, legend=['**aa', '**aaa', '**aaaa'], col=cols, lty=1, bty='n')
    r.legend(x=1e-6, y=12, legend=['Grr=1.2', 'Grr=1.4'], pch=pchs, bty='n')
    r.mtext('d', 3, line=0, adj=1.05)
    #
    # allele frequency
    myP = [x/100. for x in range(2, 15, 2)]
    myGrr = [1.4, 1.8]
    myFam = ('**AA', '**AAA', '**AAAA')
    L = allValues(myGrr, allK, allp, ['Additive'], myPower, myFam, allAlpha)
    r.plot(0, 0, axes=False, main='', xlab='Disease allele frequency', ylab=r('expression(lambda)'),
        xlim=[min(myP), max(myP)], ylim=[0, 12], type='n', yaxs='i')
    r.axis(1)
    r.axis(2)
    cols = ['blue', 'red', 'green']
    pchs = [2, 3]
    for Grr, fam in product(myGrr, myFam):
        L = allValues([Grr], allK, myP, allMode, allPower, [fam], allAlpha)
        r.lines(myP,
                [L[(Grr, allK[0], p, allMode[0], allPower[0], fam, allAlpha[0])] for p in myP],
                col=cols[myFam.index(fam)],
                pch=pchs[myGrr.index(Grr)],
                lty=1,
                type='b')
    r.legend(x=0.11, y=12, legend=['**aa', '**aaa', '**aaaa'], col=cols, lty=1, bty='n')
    r.legend(x=0.075, y=12, legend=['Grr=1.4', 'Grr=1.8'], pch=pchs, bty='n')
    r.mtext('e', 3, line=0, adj=1.05)
    # allele frequency
    r.dev_off()


 
def Table1():
    for K, p, mode, Grr in [
            (0.1, 0.1, 'Additive', 1.1),
            (0.1, 0.1, 'Additive', 1.2),
            (0.1, 0.1, 'Additive', 1.4),
            (0.05, 0.02, 'Additive', 1.2),
            (0.05, 0.02, 'Additive', 1.4),
            (0.05, 0.02, 'Additive', 1.8),
            (0.01, 0.02, 'Additive', 1.5),
            (0.01, 0.02, 'Additive', 2.0),
            (0.01, 0.02, 'Additive', 4.0),
            #
            (0.1, 0.1, 'Recessive', 1.1),
            (0.1, 0.1, 'Recessive', 1.2),
            (0.1, 0.1, 'Recessive', 1.4),
            (0.05, 0.02, 'Recessive', 1.2),
            (0.05, 0.02, 'Recessive', 1.4),
            (0.05, 0.02, 'Recessive', 1.8),
            (0.01, 0.02, 'Recessive', 1.5),
            (0.01, 0.02, 'Recessive', 2.0),
            (0.01, 0.02, 'Recessive', 4.0),
            #
            (0.1, 0.1, 'Dominant', 1.1),
            (0.1, 0.1, 'Dominant', 1.2),
            (0.1, 0.1, 'Dominant', 1.4),
            (0.05, 0.02, 'Dominant', 1.2),
            (0.05, 0.02, 'Dominant', 1.4),
            (0.05, 0.02, 'Dominant', 1.8),
            (0.01, 0.02, 'Dominant', 1.5),
            (0.01, 0.02, 'Dominant', 2.0),
            (0.01, 0.02, 'Dominant', 4.0)]:
        N, r1 = caseRatio('UUA', '**U',  K, p, mode, Grr)
        N, r2 = caseRatio('UAA', '**U',  K, p, mode, Grr)
        N, r3 = caseRatio('*AA', '**U',  K, p, mode, Grr)
        N, r4 = caseRatio('**AA', '**U', K, p, mode, Grr)
        N, r5 = caseRatio('AUAA', '**U',  K, p, mode, Grr)
        N, r6 = caseRatio('AUAAA', '**U', K, p, mode, Grr)
        N, r7 = caseRatio('**A', '*UU',  K, p, mode, Grr)
        N, r8 = caseRatio('**A', 'UUU', K, p, mode, Grr)
        N, r9 = caseRatio('**AA', '**UU', K, p, mode, Grr)
        print '%.2f, %.2f, %10s, %.1f, %10d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f' % \
                (K, p, mode, Grr, N, r1, r2, r3, r4, r5, r6, r7, r8, r9)
 
def Table1FixedCtrl():
    for K, p, mode, Grr in [
            (0.1, 0.1, 'Additive', 1.1),
            (0.1, 0.1, 'Additive', 1.2),
            (0.1, 0.1, 'Additive', 1.4),
            (0.05, 0.02, 'Additive', 1.2),
            (0.05, 0.02, 'Additive', 1.4),
            (0.05, 0.02, 'Additive', 1.8),
            (0.01, 0.02, 'Additive', 1.5),
            (0.01, 0.02, 'Additive', 2.0),
            (0.01, 0.02, 'Additive', 4.0),
            #
            (0.1, 0.1, 'Recessive', 1.1),
            (0.1, 0.1, 'Recessive', 1.2),
            (0.1, 0.1, 'Recessive', 1.4),
            (0.05, 0.02, 'Recessive', 1.2),
            (0.05, 0.02, 'Recessive', 1.4),
            (0.05, 0.02, 'Recessive', 1.8),
            (0.01, 0.02, 'Recessive', 1.5),
            (0.01, 0.02, 'Recessive', 2.0),
            (0.01, 0.02, 'Recessive', 4.0),
            #
            (0.1, 0.1, 'Dominant', 1.1),
            (0.1, 0.1, 'Dominant', 1.2),
            (0.1, 0.1, 'Dominant', 1.4),
            (0.05, 0.02, 'Dominant', 1.2),
            (0.05, 0.02, 'Dominant', 1.4),
            (0.05, 0.02, 'Dominant', 1.8),
            (0.01, 0.02, 'Dominant', 1.5),
            (0.01, 0.02, 'Dominant', 2.0),
            (0.01, 0.02, 'Dominant', 4.0)]:
        N, r1 = caseRatioFixedControl('UUA', '**U',  K, p, mode, Grr)
        N, r2 = caseRatioFixedControl('UAA', '**U',  K, p, mode, Grr)
        N, r3 = caseRatioFixedControl('*AA', '**U',  K, p, mode, Grr)
        N, r4 = caseRatioFixedControl('**AA', '**U', K, p, mode, Grr)
        N, r5 = caseRatioFixedControl('AUAA', '**U',  K, p, mode, Grr)
        N, r6 = caseRatioFixedControl('AUAAA', '**U', K, p, mode, Grr)
        N, r7 = caseRatioFixedControl('**A', '*UU',  K, p, mode, Grr)
        N, r8 = caseRatioFixedControl('**A', 'UUU', K, p, mode, Grr)
        N, r9 = caseRatioFixedControl('**AA', '**UU', K, p, mode, Grr)
        print '%.2f, %.2f, %10s, %.1f, %10d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f' % \
                (K, p, mode, Grr, N, r1, r2, r3, r4, r5, r6, r7, r8, r9)


import logging
def example():
    K = 0.01
    p = 0.1
    Grr = 1.5
    mode = 'Additive'

    # using a logger object to get more details.
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('')
    p_U = pedProbabilities('**U', mode=mode, p=p, 
            r=RelativeRiskToPenetrance(K, p, mode, Grr), x=p).PrAlleleGivenPed('X')
    p_A = pedProbabilities('**A', mode=mode, p=p, 
            r=RelativeRiskToPenetrance(K, p, mode, Grr), x=p).PrAlleleGivenPed('X')
    p_AA = pedProbabilities('**AA', mode=mode, p=p, 
            r=RelativeRiskToPenetrance(K, p, mode, Grr), x=p).PrAlleleGivenPed('X')
    pow_A = powerCalculator(mode, K=K, p=p, x=p, logger=logger).getPower(Grr, {'**U':1000}, {'**A':1000}, alpha=defaultAlpha)
    pow_AA = powerCalculator(mode, K=K, p=p, x=p, logger=logger).getPower(Grr, {'**U':1000}, {'**AA':1000}, alpha=defaultAlpha)
    pow_2000 = powerCalculator(mode, K=K, p=p, x=p, logger=logger).getPower(Grr, {'**U':2000}, {'**A':2000}, alpha=defaultAlpha)


def Houlston():
    for ped in ['**A', '**AA', '**AAA', '**AAAA']:
        for p in [0.005, 0.01, 0.02, 0.05]:
            cal = powerCalculator(mode='Dominant', K=0.05, p=p, x=p)
            N = cal.getSampleSizeFromRatio(Grr=2.0, cases={ped: 1}, controls={'**U': 1},
                    ratio=2, alpha=0.0001, power=0.95)
            p1 = cal.expectedFreq({ped:1})
            p2 = cal.expectedFreq({'**U':1})
            print '%d (%.1f%%),' % (N[ped], p1*100),
        print


if __name__ == '__main__':
    cmd = sys.argv[1]
    if not cmd.endswith(')'):
        cmd += '()'
    eval(cmd)


