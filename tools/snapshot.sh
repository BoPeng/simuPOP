#!/bin/sh

# make src and binary distribution on all platforms 
# (currently mac, linux and solaris) as daily snapshot
#
# This is part of simuPOP
# 
# Bo Peng (bpeng@rice.edu)
#

export PATH=$PATH:/usr/site/icc/bin:.:/home/bpeng/bin:/usr/site/icc/bin:/usr/local/bin:/usr/bin:/bin:/usr/site/icc/bin:.:/home/bpeng/bin:/usr/site/icc/bin:/usr/local/bin:/usr/bin:/bin:/usr/site/icc/bin:.:/home/bpeng/bin:/usr/site/icc/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin:/usr/games:/usr/X11R6/bin:/usr/X11R6/bin:/usr/bin/X11:/sbin:/usr/sbin:/home/bpeng/research/randfam/randfam:/usr/java/j2sdk1.4.2_01/jre/bin:/usr/X11R6/bin:/usr/bin/X11:/sbin:/usr/sbin:/home/bpeng/research/randfam/randfam:/usr/java/j2sdk1.4.2_01/jre/bin:/usr/X11R6/bin:/usr/bin/X11:/sbin:/usr/sbin:/home/bpeng/research/randfam/randfam:/usr/java/j2sdk1.4.2_01/jre/bin

# first, make distribution:
cd ..
find . -name '*.pyc' -exec rm -f {} \;
find . -name 'core.*' -exec rm -f {} \;
find . -name '*.pyo' -exec rm -f {} \;
#find . -name '*.so' -exec rm -f {} \;
find . -name '*.bak' -exec rm -f {} \;
find . -name '.#*' -exec rm -f {} \;
find . -name '*~' -exec rm -f {} \;
find . -name '#*#' -exec rm -f {} \;

# snapshot version.
if test "x$1" == 'x'; then
  SIMUPOP_VER=snapshot
else
  SIMUPOP_VER=$1
fi
echo "Building " $SIMUPOP_VER
SIMUPOP_REV=`svnversion .`

echo SIMUPOP_VER = '"'$SIMUPOP_VER'"' > simuPOP.release
echo SIMUPOP_REV = '"'$SIMUPOP_REV'"' >> simuPOP.release

# export, to be used in Doxyfile
# but nobody else will sue it later
export SIMUPOP_VER
export SIMUPOP_REV

# make docstring
doxygen Doxyfile
python tools/doxy2swig.py /var/www/html/simuPOP_doc/xml/index.xml tmp.i
perl tools/processDocString.pl tmp.i > src/simuPOP/download.i
rm -f tmp.i
cd doc
postGuide.sh 
cd ../tools

# remove previous build
if test -f /var/www/html/simuPOP/download/simuPOP-$SIMUPOP_VER-src.tar.gz ; then
  rm -f /var/www/html/simuPOP/download/simuPOP-$SIMUPOP_VER-src.tar.gz
fi

# build source distribution
make_src.sh

# distribute source
# this file will be used by other batch files
echo $SIMUPOP_VER > /var/www/html/simuPOP/download/latestversion

# check files
if test ! -f /var/www/html/simuPOP/download/simuPOP-$SIMUPOP_VER-src.tar.gz ; then
  echo "Can not make src distribution."
  exit
fi

# make binary distributions
make_linux.sh

make_win.sh

make_sol.sh

make_mdk.sh

make_mac.sh

# copy files to thor:
scp /var/www/html/simuPOP/download/simuPOP-$SIMUPOP_VER* thor:public_html/simuPOP
scp /var/www/html/simuPOP_doc/*.pdf thor:public_html/simuPOP

