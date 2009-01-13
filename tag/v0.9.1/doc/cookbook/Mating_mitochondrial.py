#!/usr/bin/env python

'''
File: Mating_mitochondrial.py
Author: Bo Peng (bpeng@mdanderson.org)

Purpose:
  This script demonstrates how to use mitochondrial chromosomes.

$Date: 2008-12-14 02:19:13 -0500 (Sun, 14 Dec 2008) $
$Revision: 2114 $
$HeadURL: https://simupop.svn.sourceforge.net/svnroot/simupop/trunk/doc/cookbook/Mating_mitochondrial.py $
'''

from simuPOP import *

# This example demonstrates how to simulate mitochondrial chromosomes. Such
# chromosomes are inherited maternally. If there are multiple copies of
# mitochondrial chromosomes, they are transmitted randomly.
#
# simuPOP does not define Mitochondrial as an internally recognized chromosome
# type so you have to define a number of 'Customized' chromososomes. Such
# chromosomes are not handled by 'normal' genotype transmitters such as
# mendelianGenoTransmitter. In this case, a mitochondrialGenoTransmitter is
# used, which, by default, treats all customized chromosomes as mitochondrial
# chromosomes and transmits them properly. A parameter chroms can be used to
# specify which chromosomes are mitochondrial and allow, possibly, another
# during-mating operator to handle another type of chromosome.
#
# This example passes a recombinator and a mitochondrialGenoTransmitter in
# the evolve() function of a simulator. The recombinator will replace the
# default genotype transmitter, namely mendelianGenoTransmitter, to transmit
# non-customized chromosomes with recombination. mitochondrialGenoTransmitter
# then handles the rest of the chromosomes. Note that it is also possible to
# pass these operators to the ops parameter of a mating scheme, e.g.,
#
#   randomMating(ops=[recombinator(), mitochondrialGenoTransmitter()])
#
# This is required if more than one mating schemes are used and you do not
# want to use these operators in all mating schemes.
#

def simuMitochondrial(N, numMito=3, gen=10):
    '''
    '''
    pop = population(N, loci=[5]*(3 + numMito),
        # one autosome, two sex chromosomes, and numMito mitochondrial chromosomes
        chromTypes=[Autosome, ChromosomeX, ChromosomeY] + [Customized]*numMito,
        # record indexes of parents for verification purpose
        ancGen=1, infoFields=['father_idx', 'mother_idx'])

    simu = simulator(pop, randomMating())

    simu.evolve(
        preOps=[
            # initialize alleles 0, 1, 2, 3 with different frequencies
            initByFreq([0.4] + [0.2]*3),
        ],
        ops=[
            recombinator(rate=0.1),
            mitochondrialGenoTransmitter(),
            parentsTagger(),
            dumper(structure=False),
        ],
        gen = gen
    )
    return simu.extract(0)


if __name__ == '__main__':
    simuMitochondrial(10, 3, 2)

# A possible output:
##
# Subpopulation 0 (unnamed):
#    0: FU 20001 33000 _____ 23300 23300 00131 | 02203 20300 _____ 00000 00000 00000 |  3 8
#    1: FU 21330 22311 _____ 00223 00223 00030 | 00030 01100 _____ 00000 00000 00000 |  2 7
#    2: MU 03222 11210 _____ 33002 02203 02203 | 00030 _____ 22230 00000 00000 00000 |  2 6
#    3: MU 01130 00013 _____ 20003 00233 00232 | 21131 _____ 23200 00000 00000 00000 |  1 0
#    4: FU 20001 00300 _____ 23300 00131 00131 | 21131 20031 _____ 00000 00000 00000 |  1 8
#    5: MU 03302 00000 _____ 23300 11001 11001 | 00030 _____ 22230 00000 00000 00000 |  2 8
#    6: MU 00310 11003 _____ 03000 03000 03000 | 13331 _____ 23200 00000 00000 00000 |  1 4
#    7: MU 00031 11000 _____ 03000 23030 03000 | 21313 _____ 13031 00000 00000 00000 |  9 4
#    8: FU 03030 21220 _____ 31012 02203 02203 | 30200 01100 _____ 00000 00000 00000 |  2 6
#    9: MU 22200 30003 _____ 23000 23310 23310 | 30230 _____ 22230 00000 00000 00000 |  2 5
# End of individual genotype.
# 
# Genotype of individuals in the present generation:
# Subpopulation 0 (unnamed):
#    0: MU 21330 22311 _____ 00223 00030 00223 | 21131 _____ 23200 00000 00000 00000 |  3 1
#    1: MU 20001 00301 _____ 00131 23300 00131 | 13310 _____ 23200 00000 00000 00000 |  6 4
#    2: MU 21131 00300 _____ 23300 00131 23300 | 00310 _____ 23200 00000 00000 00000 |  6 4
#    3: MU 20203 33000 _____ 23300 00131 23300 | 21313 _____ 13031 00000 00000 00000 |  7 0
#    4: MU 20001 20000 _____ 23300 23300 23300 | 00031 _____ 13031 00000 00000 00000 |  7 0
#    5: MU 20001 20030 _____ 00131 23300 00131 | 03222 _____ 22230 00000 00000 00000 |  2 4
#    6: MU 03030 01100 _____ 31012 31012 31012 | 30200 _____ 22230 00000 00000 00000 |  9 8
#    7: MU 21330 01100 _____ 00030 00030 00223 | 00011 _____ 13031 00000 00000 00000 |  7 1
#    8: FU 02203 20300 _____ 00131 00131 23300 | 13331 11003 _____ 00000 00000 00000 |  6 0
#    9: FU 30200 01100 _____ 02203 31012 02203 | 13330 11003 _____ 00000 00000 00000 |  6 8
# End of individual genotype.
# 
# * Male individuals have chromosome X (ploidy 0), and Y (ploidy 1)
# * Female individuals have chromosome X (ploidy 0), and X (ploidy 1)
# * The second homologous copy of mitochondrial chromosomes is unused.
# * The last two numbers are 'father_idx' and 'mother_idx' which are indexes of
#   each individual in the parental generation.
# 
# Taking an example of the last individual (9) in the present generation
# * Autosome 30200 is the second copy of parent 8 (0030 and 30200).
# * Autosome 13330 is the recombined copy of parent 6 (00310 and 13331)
# * First chromosome X 01100 is the second copy of chromosome X of parent 8 (21220 and 01100)
# * Second chromosome X 11003 is the first copy of chromosome X of parent 6 (11003).
#   No recombination is possible.
# * Three mitochondrials 02202, 31012, 02203 are randomly chosen from
#   mitochondrial chromosomes of parent 8 (31012, 02203, and 02203).

