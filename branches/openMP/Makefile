#
# Makefile for processing simuPOP
#


rebuild:
	bcpp src/*.h src/*.cpp
	doxygen
	tools/doxy2swig.py
	scons install -j4 std
	tools/doxy2swig.py
	scons install -j4
	cd doc; python runSampleCode.py; make dist_svn

rebuild3:
	python3 setup.py install

clean:
	@rm -rf build
	@rm -f src/*wrap* src/simuPOP_*.py
	@find . -name '*~' -exec rm -f {} \;
	@find . -name '*.bak' -exec rm -f {} \;
	@find . -name '*.py?' -exec rm -f {} \;
	@rm -f MANIFEST doc/*.aux doc/*.log doc/*.idx

OBJ = src/simuPOP_std_wrap.o src/main.o\
    src/utility.o\
    src/genoStru.o\
    src/individual.o\
    src/population.o\
    src/simulator.o\
    src/mating.o\
    src/operator.o\
    src/initializer.o\
    src/migrator.o\
    src/outputer.o\
    src/selector.o\
    src/penetrance.o\
    src/qtrait.o\
    src/stator.o\
    src/mutator.o\
    src/transmitter.o\
    src/tagger.o\
    src/pedigree.o\
    src/virtualSubPop.o\
    src/sandbox.o

build/std_common/libcommon.a:
	scons build std

src/simuPOP_std_wrap.cpp:
	swig -O -templatereduce -shadow -python -c++ -keyword -nodefaultctor -w-503,-312,-511,-362,-383,-384,-389,-315,-509,-525 -Isrc \
		-outdir src -o src/simuPOP_std_wrap.cpp src/simuPOP_std.i 

config.h:
	cp config_linux.h config.h

src/swigpyrun.h:
	swig -python -external-runtime src/swigpyrun.h

test_simuPOP: config.h src/swigpyrun.h $(OBJ) build/std_common/libcommon.a
	g++ -o test_simuPOP $(OBJ) build/std_common/libcommon.a -lz -lstdc++ -L/usr/lib  -lpython2.6

%.o: %.cpp
	g++ -O3 -Wall -DNDEBUG -DSIMUPOP_MODULE=simuPOP_std -D_SECURE_SCL=1 -D_SIMUPOP_VER=1.0.6svn -DSIMUPOP_REV=9999 \
		-I. -Igsl -Igsl/cdf -Igsl/specfunc -I/usr/local/include/python2.7 -Iboost_1_44_0 -c $< -o $@

