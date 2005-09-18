#!/bin/sh

# make binary libraries for soalris platform
#
# This is part of simuPOP
# 
# Bo Peng (bpeng@rice.edu)
#
# Oct 12, 2004
#

VER=`cat /var/www/html/simuPOP/download/latestversion`
cd /home/bpeng/temp
/bin/rm -rf simuPOP-$VER &&  /bin/rm -rf simuPOP
cp /var/www/html/simuPOP/download/simuPOP-$VER-src.tar.gz .

UNCOMPRESS="tar zxf simuPOP-$VER-src.tar.gz"
BUILD="python setup.py bdist --formats=gztar,rpm"
$UNCOMPRESS && cd simuPOP-$VER && rm -rf dist && $BUILD 
cp dist/simuPOP-$VER.linux-i686.tar.gz /var/www/html/simuPOP/download/simuPOP-$VER.linux-i686-py23.tar.gz
cp dist/simuPOP-$VER-1.i386.rpm /var/www/html/simuPOP/download/simuPOP-$VER.linux-i386-py23.rpm
cp dist/simuPOP-$VER-1.src.rpm /var/www/html/simuPOP/download/simuPOP-$VER-src.rpm

