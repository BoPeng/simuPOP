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
cd /home/bpeng/tmp
/bin/rm -rf simuPOP-$VER &&  /bin/rm -rf simuPOP
cp /var/www/html/simuPOP/download/simuPOP-$VER-src.tar.gz .

UNCOMPRESS="tar zxf simuPOP-$VER-src.tar.gz"
# this is python 2.3
BUILD="/usr/bin/python setup.py bdist --formats=gztar,rpm"
$UNCOMPRESS && cd simuPOP-$VER && rm -rf dist && $BUILD 
cp dist/simuPOP-$VER.linux-x86_64.tar.gz /var/www/html/simuPOP/download/simuPOP-$VER.linux-x86_64-py23.tar.gz
cp dist/simuPOP-$VER-1.x86_64.rpm /var/www/html/simuPOP/download/simuPOP-$VER.linux-x86_64-py23.rpm
cp dist/simuPOP-$VER-1.src.rpm /var/www/html/simuPOP/download/simuPOP-$VER-src.rpm

