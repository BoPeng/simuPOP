#!/usr/bin/env python

#
# This script is used to draw some figures that will be used in simuPOP
# documentations. Most of such figures are draw by hand, using Dia and
# exported to PNG format. However, a few figures are easier to be 
# drawn programatically.
#
#
from rpy import *

def drawGenotype():
    '''Draw a figure that demonstrates the memory layout of
       individual genotypes.
    '''
    r.pdf('genotype.pdf', height=4, width=8)
    ch0 = 10
    ch1 = 15
    ch = ch0+ch1
    pad = 10
    r.par(mar=[0, 0, 0, 0])
    r.plot(0, 0, axes=False, type='n', main='', xlab='', ylab='',
        xlim=[-(ch+pad), (ch+pad)], ylim=[-8, 4])
    for y in [0, 1]:
        r.lines([-ch, ch], [y, y], lty=1, lwd=2)
        r.lines([-(ch+pad), (ch+pad)], [y, y], lty=2)
    for i in range(-ch, ch):
        r.lines([i, i], [0, 1], lty=2)
    for i in [-ch, 0, ch]:
        r.lines([i, i], [0, 1], lty=1, lwd=2)
    for i in [-ch1, ch0]:
        r.lines([i, i], [0, 1], lty=1)
    r.text(-(ch+0.5), 0.5, 'previous ind', adj=1)
    r.text(ch+0.5, 0.5, 'next ind', adj=0)
    for x in [-ch, -ch1, 0, ch0]:
        r.text(x, 1.5, '<', adj=0)
    for x in [-ch1, 0, ch0, ch]:
        r.text(x, 1.5, '>', adj=1)
    for x in [-(ch1 + ch0/2), ch0/2]:
        r.text(x, 1.5, 'ch 0')
    for x in [-ch1/2., ch0 + ch1/2.]:
        r.text(x, 1.5, 'ch 1')
    for x in [-ch, -ch1, 0, ch0, ch]:
        r.lines([x, x], [1, 2], lty=1)
    #
    for x in [-ch, 0]:
        r.text(x, 2.5, '<', adj=0)
    for x in [0, ch]:
        r.text(x, 2.5, '>', adj=1)
    for idx,x in enumerate([-ch/2., ch/2.]):
        r.text(x, 2.5, 'ploidy %d' % idx)
    for x in [-ch, 0, ch]:
        r.lines([x, x], [2, 3], lty=1)
    #
    r.rect(-ch1+2, 0, -ch1+3, 1, density=50)
    r.rect(-(ch1+pad), -1.5, -(ch1+pad)+1, -0.5, density=50)
    r.text(-(ch1+pad) + 2, -1, 'allele(%d), allele(%d, 0) or allele(2, 0, 1)' % (ch0+2, ch0+2), adj=0)
    #
    r.rect(0, 0, ch0, 1, density=20)
    r.rect(-(ch1+pad), -2, -(ch1+pad)+1, -3, density=20)
    r.text(-(ch1+pad) + 2, -2.5, 'genotype(1, 0)', adj=0)
    #
    r.text(-(ch+pad), -4.5, 'Single-allele read:  allele(idx), allele(idx, p), allele(idx, p, ch)', adj=0)
    r.text(-(ch+pad), -5.5, 'Single-allele write: setAllele(allele, idx), setAllele(allele, idx, p), setAllele(allele, idx, p, ch)', adj=0)
    r.text(-(ch+pad), -6.5, 'Batch read:  genotype(), genotype(p), genotype(p, ch)', adj=0)
    r.text(-(ch+pad), -7.5, 'Batch write: setGenotype(), setGenotype(p), setGenotype(p, ch)', adj=0)
    r.dev_off()

if __name__ == '__main__':
    drawGenotype()
