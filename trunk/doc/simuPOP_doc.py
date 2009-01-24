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
    r.png('genotype.png', height=400, width=800)
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


def drawChromType():
    r.png('chromType.png', height=300, width=800)
    r.par(mar=[0, 0, 0, 0])
    r.plot(0, 0, axes=False, type='n', main='', xlab='', ylab='',
        ylim=[-6, 6], xlim=[-40, 40])
    pad = 0.5
    for x,y in [(-30, 1), (-30, 3), (-30, -2), (-30, -4), (5, 2), (5, -3)]:
        r.rect(x, y, x+6 - pad, y+0.5)
        r.rect(x+6, y, x+12 - pad, y+0.5)
        r.rect(x+12, y, x+16 - pad, y+0.5)
        r.rect(x+16, y, x+20 - pad, y+0.5)
        r.rect(x+20, y, x+24 - pad, y+0.5)
        r.text(x+3-pad/2, y+1, 'A')
        r.text(x+9-pad/2, y+1, 'X')
        r.text(x+14-pad/2, y+1, 'Y')
        r.text(x+18-pad/2, y+1, 'C')
        r.text(x+22-pad/2, y+1, 'C')
    r.text(-31, 2.5, 'Maternal\nchromosomes', adj=1)
    r.text(-31, -2.5, 'Paternal\nchromosomes', adj=1)
    r.text(5, 4, 'Maternally inherited chromosomes', adj=0)
    r.text(5, -1, 'Paternally inherited chromosomes', adj=0)
    # autosome
    for x1,x2,y1,col in [
        # maternal autosome
        (-30, -28, 3, 'green'), (5, 7, 2, 'green'),
        (-28, -24-pad, 1, 'red'), (7, 11-pad, 2, 'red'),
        # paternal autosome
        (-30, -25.5, -2, 'green'), (5, 9.5, -3, 'green'),
        (-25.5, -24-pad, -4, 'red'), (9.5, 11-pad, -3, 'red'),
        # maternal chromosome X
        (-24, -20, 3, 'green'), (11, 15, 2, 'green'),
        (-20, -18-pad, 1, 'red'), (15, 17-pad, 2, 'red'),
        # paternal chromsome Y
        (-18, -14-pad, -4, 'green'), (17, 21-pad, -3, 'green'),
        # mitochondria 1, 2
        (-14, -10-pad, 3, 'blue'), #(-10, -6-pad, 3, 'green'),
        (21, 25-pad, 2, 'blue'), (25, 29-pad, 2, 'blue'),
        # gray out
        (-18, -14-pad, 3, 'gray'), (-18, -14-pad, 1, 'gray'),
        (-24, -18-pad, -4, 'gray'), (-18, -14-pad, -2, 'gray'),
        (-14, -10-pad, 1, 'gray'), (-10, -6-pad, 1, 'gray'),
        (-14, -10-pad, -4, 'gray'), (-10, -6-pad, -4, 'gray'),
        (21, 25-pad, -3, 'gray'), (25, 29-pad, -3, 'gray'),
        (11, 17-pad, -3, 'gray'), (17, 21-pad, 2, 'gray'),
        # legend
        (-40, -38, -6, 'red'), (-38, -36, -6, 'blue'), (-36, -34, -6, 'green'),
        (-5, -3, -6, 'white'),
        (30, 32, -6, 'gray'), 
        ]: 
        r.rect(x1, y1, x2, y1+0.5, col=col)
    r.text(-33, -5.75, 'Transmitted chromosome regions', adj=0)
    r.text(-2, -5.75, 'Untransmitted chromosome regions', adj=0)
    r.text(33, -5.75, 'Unused', adj=0)
    r.dev_off()

if __name__ == '__main__':
    if sys.argv[1] == 'genotype':
        drawGenotype()
    elif sys.argv[1] == 'chromType':
        drawChromType()
