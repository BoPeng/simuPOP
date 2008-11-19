    # number of pedigrees
    # get unique pedgree id numbers
    from sets import Set
    peds = Set(pop.indInfo('pedindex', False))
    # do not count peds -1
    peds.discard(-1)
    #
    newPedIdx = 1
    #
    for ped in peds:
        id = 1
        # -1 means no parents
        pastmap = {-1:0}
        # go from generation 2, 1, 0 (for example)
        for anc in range(pop.ancestralDepth(), -1, -1):
            newmap = {-1:0}
            pop.useAncestralPop(anc)
            # find all individual in this pedigree
            for i in range(pop.popSize()):
                ind = pop.individual(i)
                if ind.info('pedindex') == ped:
                    dad = int(ind.info('father_idx'))
                    mom = int(ind.info('mother_idx'))
                    if dad == mom and dad != -1:
                        print ("Something wrong with pedigree %d, father and mother " + \
                            "idx are the same: %s") % (ped, dad)
                    writeInd(ind, newPedIdx, id, pastmap.setdefault(dad, 0), \
                        pastmap.setdefault(mom, 0))
                    newmap[i] = id
                    id += 1
            pastmap = newmap
        newPedIdx += 1
    pedOut.close()
