    # write map file
    print >> mapOut, 'CHROMOSOME MARKER POSITION'
    for marker in loci:
        print >> mapOut, '%d\t%s\t%f' % (pop.chromLocusPair(marker)[0] + 1, 
            pop.locusName(marker), pop.locusPos(marker))
    mapOut.close()
