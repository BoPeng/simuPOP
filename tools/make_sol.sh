#!/bin/sh

# make binary libraries for soalris platform
#
# This is part of simuPOP
# 
# Bo Peng (bpeng@rice.edu)
#
# Oct 12, 2004
#
# Note: use local gsl copy to build binary distribution
#
VER=`cat /var/www/html/simuPOP/download/latestversion`
ssh thor.stat.rice.edu '/bin/rm -rf temp &&  mkdir temp'
scp /var/www/html/simuPOP/download/simuPOP-$VER-src.tar.gz thor:temp

UNCOMPRESS="/usr/local/bin/tar zxf simuPOP-$VER-src.tar.gz"
BUILD="setenv CXX /usr/site/gcc-3.3/bin/g++; setenv CC /usr/site/gcc-3.3/bin/gcc; python setup.py bdist_dumb"
ssh -X thor.stat.rice.edu "cd temp && $UNCOMPRESS && cd simuPOP-$VER && $BUILD"
scp thor:temp/simuPOP-$VER/dist/simuPOP-$VER*.tar.gz /var/www/html/simuPOP/download/simuPOP-$VER-sol-py23.tar.gz

