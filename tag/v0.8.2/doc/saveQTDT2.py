    # write dat file
    # 
    if 'affection' in fields:
        outputAffectation = True
        fields.remove('affection')
        print >> datOut, 'A\taffection'
    else:
        outputAffectation = False
    for f in fields:
        print >> datOut, 'T\t%s' % f
    for marker in loci:
        print >> datOut, 'M\t%s' % pop.locusName(marker)
    datOut.close()
