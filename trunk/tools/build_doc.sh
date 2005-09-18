#!/bin/csh

# make src and binary distribution on all platforms 
# (currently mac, linux and solaris)
#
# This is part of simuPOP
# 
# Bo Peng (bpeng@rice.edu)
#

# first, make distribution:

# make docstring

# 1. build xml
cd ..
doxygen Doxyfile
cd tools
 
# 2. get simuPOP_doc.i
python doxy2swig.py /var/www/html/simuPOP_doc/doc/xml/index.xml tmp.i
perl processDocString.pl tmp.i > ../src/simuPOP_doc.i

# 3. post manual
cd ../doc
postGuide.sh 
 
  
# 4. compile library (regernerate .py file)
cd ../src
rm -f simuPOP_wrap.cpp 
make -f Makefile.old simuPOP_wrap.cpp
