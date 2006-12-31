#!/usr/bin/env python

import os, sys, shutil
from snapshot import *

def build_src():
    ver = open(last_version_file).read()

    if ver == 'snapshot':
        rev = cmdOutput('svnversion .')
        if ':' in rev:
            rev = rev.split(':')[0]
    else:
        rev = open(last_revision_file).read()
    #
    # replace simuPOP.release file
    (old_ver, old_rev) = writeReleaseFile(ver, rev)
    # build source
    os.system('python setup.py sdist --formats=gztar,zip')
    # write old release file back
    writeReleaseFile(old_ver, old_rev)
    # coppy files
    shutil.copy('dist/simuPOP-%s.tar.gz' % ver, '%s/simuPOP-%s.tar.gz' % (download_directory, ver))
    shutil.copy('dist/simuPOP-%s.zip' % ver, '%s/simuPOP-%s.zip' % (download_directory, ver))

if __name__ == '__main__':
    os.chdir('..')
    build_src()
