from simuPOP import *
from math import sqrt
from random import uniform, randint, normalvariate
from sys import argv, exit

#BATCHTESTING SKIP

if len(argv)<1:
    print "Syntax: epiState"
    print "epiState=0 wild-recessive / =1 wild-dominant / =2 random "
    exit(-1)
else:
    epiState = int(argv[1])
    #0 wild-recessive
    #1 wild-dominant
        #2 random
    
MATERADIUS=30
WALKRADIUS=10
POPSIZE=1000
maxX = 100.0
maxY = 100.0
LOCIPERDIM=5
 
pop = Population(POPSIZE, 2,
    loci=[1]*LOCIPERDIM*2,
    infoFields = ['x', 'y']
)

for ind in pop.individuals():
    ind.setInfo(uniform(0,maxX), pop.infoIdx('x'))
    ind.setInfo(uniform(0,maxY), pop.infoIdx('y'))


def geoChooser(pop, sp):
    while (True):
        gender = 2
        while gender != 1:
            g1Idx = randint(0, pop.popSize()-1)
            gender = pop.individual(g1Idx).sex()
        g1 = pop.individual(g1Idx)
        g1x = g1.info('x')
        g1y = g1.info('y')
        g2Cases = []
        for idx in range(0, pop.popSize()-1):
            tempG2 = pop.individual(idx)
            if tempG2.sex()==2:
                tg2x = tempG2.info('x')
                tg2y = tempG2.info('y')
                dist =  sqrt((g1x-tg2x)**2 + (g1y-tg2y)**2)
                if dist<=MATERADIUS: g2Cases.append(idx)
        if len(g2Cases)==0: continue
        yield g1Idx, g2Cases[randint(0,len(g2Cases)-1)]


def dmp(pop):
    print "pSize", pop.popSize()    
    for ind in pop.individuals():
        print ind.info(pop.infoIdx('x')),
        print ind.info(pop.infoIdx('y')),
        numLoci = LOCIPERDIM*2
        numLoci *= 2
        for i in range(numLoci): print ind.allele(i),
        print
        pass
    return True



def getOptimalCoord(maxPos, numMuts):
    return numMuts*maxPos/LOCIPERDIM



def placeIndividual(pop, off, dad, mom):
    cX = (dad.info('x') + mom.info('x'))/2
    cY = (dad.info('y') + mom.info('y'))/2
    #This is really a SQUARE
    x=-1
    y=-1
    xMut, yMut = countMutants(off)
    optX = getOptimalCoord(maxX, xMut)
    optY = getOptimalCoord(maxY, yMut)
    while x<0 or x>maxX: # or x<cX-WALKRADIUS or x>cX+WALKRADIUS:
        #We will assure that the most distant extreme is within 68?%
        d1X = abs((cX-WALKRADIUS)-optX) 
        d2X = abs((cX+WALKRADIUS)-optX)
        if d1X>d2X:
            bigger=cX-WALKRADIUS
            smaller=cX+WALKRADIUS
        else:
            bigger=cX+WALKRADIUS
            smaller=cX-WALKRADIUS
        if smaller<0: smaller=0
        if bigger<0: bigger=maxX
        x = normalvariate(smaller, WALKRADIUS)
    while y<0 or y>maxY: # or y<cY-WALKRADIUS or y>cY+WALKRADIUS:
        #We will assure that the most distant extreme is within 68?%
        d1Y = abs((cY-WALKRADIUS)-optY) 
        d2Y = abs((cY+WALKRADIUS)-optY)
        if d1Y>d2Y:
            bigger=cY-WALKRADIUS
            smaller=cY+WALKRADIUS
        else:
            bigger=cY+WALKRADIUS
            smaller=cY-WALKRADIUS
        if smaller<0: smaller=0
        if bigger<0: bigger=maxY
        y = normalvariate(smaller, WALKRADIUS)

    off.setInfo(x, pop.infoIdx('x'))
    off.setInfo(y, pop.infoIdx('y'))
    return True




def countMutants(ind):
    xMut = 0
    yMut = 0
    for i in range(LOCIPERDIM):
        if epiState==0: #wild-recissive
            if (ind.allele(i) == 1 or ind.allele(i+LOCIPERDIM*2) == 1): xMut+=1
            if (ind.allele(i+LOCIPERDIM) == 1 or ind.allele(i+LOCIPERDIM+LOCIPERDIM*2) == 1): yMut+=1
        elif epiState==1: #wild-dominant
            if (ind.allele(i) == 1 and ind.allele(i+LOCIPERDIM*2) == 1): xMut+=1
            if (ind.allele(i+LOCIPERDIM) == 1 and ind.allele(i+LOCIPERDIM+LOCIPERDIM*2) == 1): yMut+=1
        elif epiState ==2: #random
            xMut += ind.allele(i+randint(0,1)*2*LOCIPERDIM)
            yMut += ind.allele(i+LOCIPERDIM+randint(0,1)*2*LOCIPERDIM)
    return xMut, yMut


def getPenalty(maxPos, pos, numMuts):
    return 1.0 - abs(pos - getOptimalCoord(maxPos, numMuts))/maxPos

def killUnfit(pop):
    indivs = []
    for indPos in range(pop.popSize()):
        ind     = pop.individual(indPos)
        x       = ind.info('x')
        y       = ind.info('y')
        xMut, yMut = countMutants(ind)
        xFactor = getPenalty(maxX, x, xMut) 
        yFactor = getPenalty(maxY, y, yMut) 
        luck    = uniform(0,1)
        if luck>xFactor*yFactor:
            #print x , y, xFactor, yFactor, luck
            indivs.append(indPos)
    pop.removeIndividuals(indivs)
    return True



pop.evolve(
    initOps = [
        InitGenotype(freq=[1,0])
    ],
    preOps = [
        PyOperator(killUnfit),
        KamMutator(k=2, rates=[0.01]*LOCIPERDIM*2, loci=range(LOCIPERDIM*2)),
    ],
    matingScheme = HomoMating(
        PyParentsChooser(geoChooser),
        OffspringGenerator(ops=[
            MendelianGenoTransmitter(),#numOffspring=(UniformDistribution, 2, 4)),
            PyOperator(placeIndividual)]),
        subPopSize=[POPSIZE]
    ),
    postOps = PyOperator(dmp),
    gen = 1000
)
