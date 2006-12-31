#!/usr/bin/env python

import os, sys, shutil
from snapshot import *

def build_linux():
    ver = open(last_version_file).read()
    #
    if not os.path.isfile('%s/simuPOP-%s.zip' % (download_directory, ver)):
        print 'Source package does not exist. Please run build_src.py first'
        sys.exit(1)
    # 
    # build
    print 'Copying source package to user temp directory...'
    os.system('/bin/rm -rf simuPOP-%s' % ver)
    shutil.copy('%s/simuPOP-%s-src.tar.gz' % (download_directory, ver), 'simuPOP-%s-src.tar.gz' % ver)
    print 'Unpacking ...'
    os.system('tar zxf simuPOP-%s-src.tar.gz' % ver)
    os.chdir('simuPOP-%s' % ver)
    print 'Building ...'
    os.system('python setup.py bdist --formats=gztar,rpm')
    # coppy files
    shutil.copy('dist/simuPOP-%s-1.x86_64.rpm' % ver, '%s/simuPOP-%s-1.x86_64.rpm' % (download_directory, ver))
    shutil.copy('dist/simuPOP-%s-1.src.rpm' % ver, '%s/simuPOP-%s-1.src.rpm' % (download_directory, ver))

if __name__ == '__main__':
    os.chdir(user_tmp_directory)
    build_linux()
