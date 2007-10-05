    # write ped file
    def sexCode(ind):
        if ind.sex() == Male:
            return 1
        else:
            return 2
    # disease status: in linkage affected is 2, unaffected is 1
    def affectedCode(ind):
        if ind.affected():
            return 'a'
        else:
            return 'u'
    #
    pldy = pop.ploidy()
    def writeInd(ind, famID, id, fa, mo):
        print >> pedOut, '%d %d %d %d %d' % (famID, id, fa, mo, sexCode(ind)),
        if outputAffectation:
            print >> pedOut, affectedCode(ind),
        for f in fields:
            print >> pedOut, '%.3f' % ind.info(f),
        for marker in loci:
            for p in range(pldy):
                print >> pedOut, "%d" % (ind.allele(marker, p) + shift), 
        print >> pedOut

