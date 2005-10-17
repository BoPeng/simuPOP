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
SIMUPOP_VER=snapshot
SIMUPOP_REV=`svnversion .`

# export, to be used in Doxyfile
export SIMUPOP_VER
export SIMUPOP_REV

echo "Get revision number and setup src/utility.cpp"
# perl -pi.bak -e "s/unknown-version/$SIMUPOP_VER/" src/simupop_cfg.h
perl -pi.bak -e "s/^#define SIMUPOP_VER.*$/#define SIMUPOP_VER \"$SIMUPOP_VER\"/" src/utility.cpp
perl -pi.bak -e "s/^#define SIMUPOP_REV.*$/#define SIMUPOP_REV \"$SIMUPOP_REV\"/" src/utility.cpp
perl -pi.bak -e "s/^SIMUPOP_VER=.*$/SIMUPOP_VER=\"$SIMUPOP_VER\"/" setup.py

# make docstring
doxygen Doxyfile
python tools/doxy2swig.py /var/www/html/simuPOP_doc/xml/index.xml tmp.i
perl tools/processDocString.pl tmp.i > src/simuPOP/download.i
rm -f tmp.i
cd doc
postGuide.sh 
cd ../tools

# remove previous build
rm -f /var/www/html/simuPOP/download/simuPOP-$SIMUPOP_VER-src.tar.gz

# build source distribution
make_src.sh

# distribute source

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

# this is to avoid that svn will update utility.cpp every time when a revision 
# number has changed.
perl -pi.bak -e "s/^#define SIMUPOP_VER.*$/#define SIMUPOP_VER \"snapshot\"/" ../src/utility.cpp
perl -pi.bak -e "s/^#define SIMUPOP_REV.*$/#define SIMUPOP_REV \"0\"/" ../src/utility.cpp
perl -pi.bak -e "s/^SIMUPOP_VER=.*$/SIMUPOP_VER=\"snapshot\"/" ../setup.py
