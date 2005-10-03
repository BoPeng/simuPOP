#!/bin/sh

cd ..
VER=`cat /var/www/html/simuPOP/download/latestversion`
rm -rf dist

# generate wrap files for other libraries
SWIG='swig  -shadow -python -keyword -w-503,-312,-511,-362,-383,-384,-389,-315,-525 -nodefault -c++'
echo "Generating SWIG wrap files"
$SWIG   -o src/simuPOP_std_wrap.cpp src/simuPOP_std.i
$SWIG   -DOPTIMIZED -o src/simuPOP_op_wrap.cpp src/simuPOP_op.i
$SWIG   -DLONGALLELE -o src/simuPOP_la_wrap.cpp src/simuPOP_la.i
$SWIG   -DLONGALLELE -DOPTIMIZED -o src/simuPOP_laop_wrap.cpp src/simuPOP_laop.i
perl -pi.bak -e 's/static PyMethodDef SwigMethods\[\] = \{/static PyMethodDef SwigMethods[] = {\n  {"carray", a_array, METH_VARARGS, a_array_doc}, \/* added by Bo Peng *\//' src/simuPOP_std_wrap.cpp
perl -pi.bak -e 's/static PyMethodDef SwigMethods\[\] = \{/static PyMethodDef SwigMethods[] = {\n  {"carray", a_array, METH_VARARGS, a_array_doc}, \/* added by Bo Peng *\//' src/simuPOP_op_wrap.cpp
perl -pi.bak -e 's/static PyMethodDef SwigMethods\[\] = \{/static PyMethodDef SwigMethods[] = {\n  {"carray", a_array, METH_VARARGS, a_array_doc}, \/* added by Bo Peng *\//' src/simuPOP_la_wrap.cpp
perl -pi.bak -e 's/static PyMethodDef SwigMethods\[\] = \{/static PyMethodDef SwigMethods[] = {\n  {"carray", a_array, METH_VARARGS, a_array_doc}, \/* added by Bo Peng *\//' src/simuPOP_laop_wrap.cpp

python setup.py sdist
cp dist/simuPOP-$VER.tar.gz /var/www/html/simuPOP/download/simuPOP-$VER-src.tar.gz
cp dist/simuPOP-$VER.zip /var/www/html/simuPOP/download/simuPOP-$VER-src.zip
