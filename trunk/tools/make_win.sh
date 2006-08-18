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


# start vm
vmrun start /vmware/WinXPPro/Windows\ XP\ Professional.vmx
sleep 30

VER=`cat /var/www/html/simuPOP/download/latestversion`
ssh -p 8890 vmwin "/bin/rm -rf simuPOP-$VER.zip simuPOP-$VER.tar.gz simuPOP-$VER"
scp -P 8890 /var/www/html/simuPOP/download/simuPOP-$VER-src.tar.gz vmwin:

UNCOMPRESS="tar zxf simuPOP-$VER-src.tar.gz"
BUILD23="/cygdrive/c/python23/python.exe setup.py build --compiler=mingw32 bdist_wininst"
BUILD24="/cygdrive/c/python24/python.exe setup.py build --compiler=mingw32 bdist_wininst"
ssh -X -p 8890 vmwin "$UNCOMPRESS && cd simuPOP-$VER && $BUILD23 "
scp -P 8890 vmwin:simuPOP-$VER/dist/*.exe /var/www/html/simuPOP/download/
ssh -X -p 8890 vmwin "cd ~; rm -rf simuPOP-$VER && $UNCOMPRESS && cd simuPOP-$VER && $BUILD24"
scp -P 8890 vmwin:simuPOP-$VER/dist/*.exe /var/www/html/simuPOP/download/
vmrun suspend /vmware/winXPPro/Windows\ XP\ Professional.vmx

