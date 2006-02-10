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


# set display, since this script will be launched from crontab.
vmrun start /vmware/Mandrake/Mandrake\ Linux.vmx
# wait for the vm to start.
sleep 30

VER=`cat /var/www/html/simuPOP/download/latestversion`
ssh -p 8891 vmmdk "/bin/rm -rf simuPOP-$VER.zip simuPOP-$VER.tar.gz simuPOP-$VER"
scp -P 8891 /var/www/html/simuPOP/download/simuPOP-$VER-src.tar.gz vmmdk:

UNCOMPRESS="tar zxf simuPOP-$VER-src.tar.gz"
BUILD="python setup.py bdist --formats=gztar,rpm"
ssh -X -p 8891 vmmdk "$UNCOMPRESS && cd simuPOP-$VER && $BUILD"
scp -P 8891 vmmdk:simuPOP-$VER/dist/simuPOP-$VER.linux-i686.tar.gz /var/www/html/simuPOP/download/simuPOP-$VER-mdk-py23.tar.gz
scp -P 8891 vmmdk:simuPOP-$VER/dist/simuPOP-$VER*86.rpm /var/www/html/simuPOP/download/simuPOP-$VER-mdk-i586-py23.rpm

vmrun suspend /vmware/Mandrake/Mandrake\ Linux.vmx

