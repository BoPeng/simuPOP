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
vmrun start /vmware/RHEL4/RHEL4.vmx
# wait for the vm to start.
sleep 30

VER=`cat /var/www/html/simuPOP/download/latestversion`
ssh -p 8892 vmrhel4 "/bin/rm -rf simuPOP-$VER.zip simuPOP-$VER.tar.gz simuPOP-$VER"
scp -P 8892 /var/www/html/simuPOP/download/simuPOP-$VER-src.tar.gz vmrhel4:

UNCOMPRESS="tar zxf simuPOP-$VER-src.tar.gz"
BUILD="python setup.py config --include-dirs=/usr/include/linux bdist --formats=gztar,rpm"
ssh -X -p 8892 vmrhel4 "$UNCOMPRESS && cd simuPOP-$VER && $BUILD"
scp -P 8892 vmrhel4:simuPOP-$VER/dist/simuPOP-$VER.linux-i686.tar.gz /var/www/html/simuPOP/download/simuPOP-$VER-rhel4-py23.tar.gz
scp -P 8892 vmrhel4:simuPOP-$VER/dist/simuPOP-$VER*86.rpm /var/www/html/simuPOP/download/simuPOP-$VER-rhel4-i586-py23.rpm

vmrun suspend /vmware/RHEL4/RHEL4.vmx

