#Copyright 2009 - Tiago Antao
#License is GNU Public License v3 - GPL v3
from pylab import *
from sys import stdin

#BATCHTESTING SKIP

currImg = -1
l = stdin.readline()
LOCIPERDIM=5
INTERVAL=5

#Colours based on just the first allele (ignores the second)
#Should be OK

def getCombos(n):
    if n > 1:
        moreCombos = getCombos(n-1)
        for combo in moreCombos:
            yield '0_' + combo
            yield '1_' + combo
    else:
        yield '0'
        yield '1'

while l != '':
    l = l.rstrip()
    if l.startswith('pSize'):
        if currImg != -1 and currImg % INTERVAL == 0:
            currImgStr = str(currImg)
            currImgStr = '0'*(4-len(currImgStr)) + currImgStr
            title(currImgStr + "  " + popSize)
            for combo in getCombos(LOCIPERDIM*2):
                print combo
                colors = map (lambda x: float(x), combo.split('_'))
                g=0
                b=0
                for i in range(LOCIPERDIM):
                    g += 1.0*colors[i]/(1.0*LOCIPERDIM)
                    b += 1.0*colors[i+LOCIPERDIM]/(1.0*LOCIPERDIM)
                exec 'plot(x' + combo + ', y' + combo + ", ',', color=(0," +str(g)+ "," +str(b)+"))"
                #plot(x,y,'b,')
            print
            axis([0,100,0,100])
            savefig('img' + currImgStr + '.png')
            clf()
        popSize = l.split(' ')[1]
        for combo in getCombos(LOCIPERDIM*2):
            exec 'x' + combo + '= []'
            exec 'y' + combo + '= []'
        currImg += 1
        l = stdin.readline()
        continue
    if currImg==-1:
        l = stdin.readline()
        continue
    toks = l.split(' ')
    exec 'x'+'_'.join(toks[2:LOCIPERDIM*2+2])+'.append(float(toks[0]))'
    exec 'y'+'_'.join(toks[2:LOCIPERDIM*2+2])+'.append(float(toks[1]))'
    l = stdin.readline()
    
