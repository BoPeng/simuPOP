#!/bin/sh

cd ..
VER=`cat /var/www/html/simuPOP/download/latestversion`
rm -rf dist

# generate wrap files for other libraries
#SWIG='swig  -shadow -python -keyword -w-503,-312,-511,-362,-383,-384,-389,-315,-525 -nodefault -c++'
SWIG='swig -O -templatereduce -shadow -python -c++ -keyword -nodefaultctor -w-503,-312,-511,-362,-383,-384,-389,-315,-525'
echo "Generating SWIG wrap files"
swig -python -external-runtime src/swigpyrun.h
$SWIG   -o src/simuPOP_std_wrap.cpp src/simuPOP_std.i
$SWIG   -DOPTIMIZED -o src/simuPOP_op_wrap.cpp src/simuPOP_op.i
$SWIG   -DLONGALLELE -o src/simuPOP_la_wrap.cpp src/simuPOP_la.i
$SWIG   -DLONGALLELE -DOPTIMIZED -o src/simuPOP_laop_wrap.cpp src/simuPOP_laop.i

python setup.py sdist
cp dist/simuPOP-$VER.tar.gz /var/www/html/simuPOP/download/simuPOP-$VER-src.tar.gz
cp dist/simuPOP-$VER.zip /var/www/html/simuPOP/download/simuPOP-$VER-src.zip
