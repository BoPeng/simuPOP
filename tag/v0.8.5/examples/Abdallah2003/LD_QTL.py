# simulation for
# Abdallah2003
#
#

import simuOpt
simuOpt.setOptions(optimized=True)

from simuPOP import *
from simuUtil import *
from simuRPy import *

# if validPops.bin already exists, load
# it directly. Otherwise, re-simulate
#
# population size
sz = 1000
# SNP or multiallelic marker? maxAllele=1 or 5
ma = 5
try:
    validPops = LoadPopulations('validPops.bin')
    samples = len(validPops)
except:
    # parameter
    # number of valid populations to be sampled.
    samples = 1000

    # 30 marker loci with the a QTL locus (location 0),
    simu = simulator(population(size=sz, ploidy=2, 
        loci=[31], lociPos=[[x/30.*3 for x in range(0,31)]],
        infoFields=['fitness']),
        randomMating(), rep=1)

    valid = 0
    fixed = 0
    validPops = []
    while True:
        if valid >= samples:
            break

        simu.setGen(0)
        simu.evolve(
            preOps = [
                # init SNP
                #initByFreq([.5,.5], atLoci=range(1,31)),
                # init MST
                initByValue(value=[[0]*30,[1]*30,[2]*30,[3]*30,[4]*30],
                    atLoci=range(1, 31), proportions=[.2]*5),
                # init QTL
                initByValue([0], atLoci=[0]),
            ],
            ops = [
                # mutation
                smmMutator(rate=1e-4, atLoci=range(1,31), maxAllele=99),
                # calculate LD
                stat(alleleFreq=range(0,31)),
                # LD=[[0,x] for x in [1,10,25]]),
                #varPlotter("[LD_prime[0][1],LD_prime[0][10],LD_prime[0][25]]",
                #    varDim=3, byRep=1, history=True, update=10, ylim=[0,1]),
                # introduce mutation at 100 gen to 10 individuals
                pointMutator(inds=[0], atLoci=[0], toAllele=1, at = [100]),
                # recombination, with rate 0.001. 
                recombinator(rate=0.001),
                # selective advantage
                mapSelector(locus=0, fitness={'0-0':1, '0-1':1.5, '1-1':2},
                    begin = 100, end=110),
                # check fixation
                terminateIf('alleleFreq[0][0]==1.', begin = 110),
                #pyEval('gen,alleleFreq[0][0]', begin = 100),
                #endl(rep=REP_LAST, begin = 100),
            ],
            end=200
        )
        if simu.gen() < 200:
            print "fixed at gen ", simu.gen(), "(count: " , fixed, ")"
            fixed += 1
        else:     # save valid population
            print "valid (count:", valid, ")"
            valid += 1
            validPops.append(simu.getPopulation(0))

    # save these valid populations for easy access
    SavePopulations(validPops, 'validPops.bin')

# now I have a list of valid population, then what?

for i in range(0,samples):
    # calculate statistics
    s = stat(LD=[ [0,x] for x in range(1,31) ], alleleFreq=range(0,31),
        midValues=False)
    s.apply(validPops[i])

# separate into groups by allelefreq
LDprime = []
R2 = []
var_LD = []
var_R2 = []
count = [0]*5
for i in range(0,5):
    LDprime.append([0]*30 )
    R2.append([0]*30 )    
    var_LD.append([0]*30 )
    var_R2.append([0]*30 )
    
# get average 
for i in range(0,samples):
    fq = 1. - validPops[i].dvars().alleleFreq[0][1]
    if fq < 0.05:
        idx = 0
    elif fq < 0.10:
        idx = 1
    elif fq < 0.15:
        idx =2
    elif fq < 0.20:
        idx = 3
    else:
        idx = 4
    count[idx] += 1
    for j in range(1,31):
        LDprime[idx][j-1] += validPops[i].dvars().LD_prime[0][j]
        R2[idx][j-1] += sqrt(validPops[i].dvars().R2[0][j])
        var_LD[idx][j-1] += validPops[i].dvars().LD_prime[0][j]**2
        var_R2[idx][j-1] += validPops[i].dvars().R2[0][j]
    
        
# get average and variance
for i in range(0,5):
    if count[i] > 0:
        for j in range(0,30):
            LDprime[i][j] /= count[i]
            R2[i][j] /= count[i]
            var_LD[i][j] = (var_LD[i][j] - count[i]*(LDprime[i][j]**2))/(count[i]-1)
            var_R2[i][j] = (var_R2[i][j] - count[i]*(R2[i][j]**2))/(count[i]-1)


# plot average LD'
r.X11()
r.par(mfrow=[2,2])
# mean LD'
r.plot(LDprime[0], type='l', ylab="mean LD'", xlab='dist', ylim=[0,1],
    lty=1, main='mean LD')
for i in range(1,5):
    r.lines(LDprime[i], type='l', lty=i+1)
r.legend(x=22, y=1, lty=range(1,6), legend=['p<0.05','0.05<p<0.10',
    '0.10<p<0.15','0.15<p<0.20','p>0.20'])
# mean R2
r.plot(R2[0], type='l', ylab="mean R", xlab='dist', ylim=[0,1],
    lty=1, main="mean R")
for i in range(1,5):
    r.lines(R2[i], type='l', lty=i+1)

# variance LD
r.plot(var_LD[0], type='l', ylab="var LD'", xlab='dist', ylim=[0,.1],
    lty=1, main="variance LD")
for i in range(1,5):
    r.lines(var_LD[i], type='l', lty=i+1)

# variance of R2
r.plot(var_R2[0], type='l', ylab="var R'", xlab='dist', ylim=[0,.02],
    lty=1, main="variance R")
for i in range(1,5):
    r.lines(var_R2[i], type='l', lty=i+1)



# apply quantitative trait
for i in range(0,samples):
    # calculate statistics
    q = mapQuanTrait(locus=0, qtrait={'0-0':1, '0-1':1.5, '1-1':2},
        sigma=.1)
    q.apply(validPops[i])

# regression model:
#     yi = mu + sum_j bj xij + ei
#
#
# for each sample (validPops[i]),
# there are sz=1000 individuals
# so sz QTLs
# mu: mean QTL
# i : individual
# yi: qtl of individual i
# x_ij: number of allele j of a marker for individual i
#     NOT alleleNum[marker][allele]
#     IS        population(i).allele(M,0) == j
#         + population(i).allele(M,1) == j
#     where
#         M is locus index, 0, 1 are ploidy index, j is allele
def xij(ind, locus, allele):
    return (ind.allele(locus, 0) == allele) + \
        (ind.allele(locus, 1) == allele)


pvalue = [[0]*30]*samples
x = carray('i', [0]*sz*5)
# so, for each sample
for s in range(0, samples):
    print 'sample ', s
    pop = validPops[s]
    # each locus 1-30
    r.assign("y", pop.dvars().qtrait)
    for loc in range(1, 31):
        # each individual
        for i in range(0, sz):
            ind = pop.individual(i)
            # j can be 1,2,3,4,5
            for ale in range(1,6):
                x[i*5+ale-1] = xij(ind, loc, ale)
        # send x to R and form a matrix
        r.assign('xij', x)
        r("""xij = matrix(xij, ncol=5, byrow=TRUE)
                 res=anova(lm(y~xij))$'Pr(>F)'[1]""")
        pvalue[s][loc-1] = r("res")
    # write. and make sure I will get partial result
    # even if the process is disrupted.
    pfile = open('pvalues.log','a')
    pfile.write(str(pvalue[s])+'\n')
    pfile.close()

# find the average of p-values and plot
r.assign('pvalue', pvalue)
r(""" pvalue = do.call('cbind', pvalue)
    mp = apply(pvalue, 2, mean)
    plot(1:30, mp)""")


        
