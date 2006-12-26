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
ssh apple-20-147.stat.rice.edu '/bin/rm -rf temp &&  mkdir temp && /bin/rm -rf simuPOP'
scp /var/www/html/simuPOP/download/simuPOP-$VER-src.tar.gz apple-20-147.stat.rice.edu:temp

UNCOMPRESS="tar zxf simuPOP-$VER-src.tar.gz"
BUILD="python setup.py build_ext --include-dirs=/Users/bpeng/boost/include/boost-1_33_1 --library-dirs=/Users/bpeng/boost/lib bdist_dumb"
ssh -X apple-20-147.stat.rice.edu "cd temp && $UNCOMPRESS && cd simuPOP-$VER && $BUILD"
scp apple-20-147.stat.rice.edu:temp/simuPOP-$VER/dist/* /var/www/html/simuPOP/download/

