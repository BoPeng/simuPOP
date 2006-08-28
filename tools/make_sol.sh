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
ssh cedric '/bin/rm -rf temp &&  mkdir temp'
scp /var/www/html/simuPOP/download/simuPOP-$VER-src.tar.gz cedric:temp

UNCOMPRESS="gunzip -d simuPOP-$VER-src.tar.gz; tar xf simuPOP-$VER-src.tar"
BUILD="python setup.py bdist_dumb"
ssh -X cedric "cd temp && $UNCOMPRESS && cd simuPOP-$VER && $BUILD"
scp cedric:temp/simuPOP-$VER/dist/simuPOP-$VER*.tar.gz /var/www/html/simuPOP/download/simuPOP-$VER-sol-py23.tar.gz

