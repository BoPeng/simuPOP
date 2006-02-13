#!/bin/sh

# make binary libraries for soalris platform
#
# This is part of simuPOP
# 
# Bo Peng (bpeng@rice.edu)
#
# Oct 12, 2004
#
# Note: use local gsl distribution to build binary
# 

VER=`cat /var/www/html/simuPOP/download/latestversion`
ssh apple-20-147 '/bin/rm -rf temp &&  mkdir temp && /bin/rm -rf simuPOP'
scp /var/www/html/simuPOP/download/simuPOP-$VER-src.tar.gz apple-20-147:temp

UNCOMPRESS="tar zxf simuPOP-$VER-src.tar.gz"
BUILD="perl -pi.bak -e 's/WITH_SERIALIZATION=1/WITH_SERIALIZATION=0/' setup.py; python setup.py bdist --format=gztar"
ssh -X apple-20-147.stat.rice.edu "cd temp && $UNCOMPRESS && cd simuPOP-$VER && $BUILD"
scp apple-20-147:temp/simuPOP-$VER/dist/* /var/www/html/simuPOP/download/

