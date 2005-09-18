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


scp -P 8890 /var/www/html/simuPOP/download/rpy-0.4.1-R-2.1.0.tar.gz localhost:

UNCOMPRESS="tar zxf rpy-0.4.1-R-2.1.0.tar.gz"
BUILD="/cygdrive/c/python23/python.exe setup.py build --compiler=mingw32 bdist_wininst"
BUILD="/cygdrive/c/python23/python.exe setupall.py build --compiler=mingw32 bdist_wininst"
ssh -X -p 8890 vmwin "$UNCOMPRESS && cd rpy-0.4.1-R-2.1.0 && $BUILD"
scp -P 8890 vmwin:simuPOP-$VER/dist/*.exe /var/www/html/simuPOP/download/

