def VC_merlin(file, merlin='merlin'):
    ''' run variance component method 
        file: file.ped, file.dat, file.map and file,mdl are expected.
            file can contain directory name.
    '''
    cmd = 'merlin -d %s.dat -p %s.ped -m %s.map --pair --vc' % (file, file, file)
    resline = re.compile('\s+([\d.+-]+|na)\s+([\d.+-]+|na)%\s+([\d.+-]+|na)\s+([\d.+-]+|na)\s+([\d.+-]+|na)')
    print "Running", cmd
    fout = os.popen(cmd)
    pvalues = []
    for line in fout.readlines():
        try:
            # currently we only record pvalue
            (pos, h2, chisq, lod, pvalue) = resline.match(line).groups()
            try:
                pvalues.append(float(pvalue))
            except:
                # na?
                pvalues.append(-1)
        except AttributeError:
            pass
    fout.close()
    return pvalues

