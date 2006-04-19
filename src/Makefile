#
# Make file to build simuPOP
#
# This is an alternative way to 'python setup.py build', using this
# can save some compile time.
# 
# I include path for two systems, but for convenience.
PYTHON_INC = -I/usr/include/python2.3/ -I/usr/local/python-2.3.3/include/python2.3
BOOST_INC = -I../
GSL_INC = -I../
INCLUDES= -I. $(PYTHON_INC) $(BOOST_INC) $(GSL_INC)
SERIAL_SRC=serialization
IOSTREAMS_SRC=iostreams

CC = gcc
CXX = g++
LINK = g++ -lz
CFLAGS = -Wall -pthread -fno-strict-aliasing -DNDEBUG -O3 -g  -fPIC 
# for swig 1.3.28 cvs
SWIG = swig -O -templatereduce -shadow -python -keyword -nodefaultctor -c++ -w-503,-312,-511,-362,-383,-384,-389,-315,-525
# for swig 1.3.27
#SWIG = swig -keyword -w-312,-401,-503,-511,-362,-383,-384,-389,-315,-525 -nodefault 

RM = rm -f

DEF_std=
DEF_op=-DOPTIMIZED
DEF_la=-DLONGALLELE 
DEF_laop=-DOPTIMIZED -DLONGALLELE 
DEF_ba=-DBINARYALLELE
DEF_baop=-DOPTIMIZED -DBINARYALLELE

HEADFILES=../config.h simupop_cfg.h migrator.h population.h simulator.h terminator.h \
individual.h mutator.h stator.h  utility.h \
initializer.h operator.h recombinator.h arraymodule.c \
mating.h  outputer.h selector.h  tagger.h swigpyrun.h

SERIAL_OBJS= \
$(SERIAL_SRC)/basic_archive.o \
$(SERIAL_SRC)/basic_xml_archive.o \
$(SERIAL_SRC)/basic_iarchive.o \
$(SERIAL_SRC)/basic_oarchive.o \
$(SERIAL_SRC)/basic_serializer_map.o \
$(SERIAL_SRC)/basic_text_iprimitive.o \
$(SERIAL_SRC)/basic_text_oprimitive.o \
$(SERIAL_SRC)/text_iarchive.o \
$(SERIAL_SRC)/text_oarchive.o \
$(SERIAL_SRC)/binary_iarchive.o \
$(SERIAL_SRC)/binary_oarchive.o \
$(SERIAL_SRC)/xml_grammar.o \
$(SERIAL_SRC)/xml_iarchive.o \
$(SERIAL_SRC)/xml_oarchive.o \
$(SERIAL_SRC)/extended_type_info.o \
$(SERIAL_SRC)/extended_type_info_no_rtti.o \
$(SERIAL_SRC)/extended_type_info_typeid.o \
$(SERIAL_SRC)/void_cast.o \
$(SERIAL_SRC)/polymorphic_iarchive.o \
$(SERIAL_SRC)/polymorphic_oarchive.o \
$(SERIAL_SRC)/stl_port.o \
$(SERIAL_SRC)/basic_pointer_iserializer.o \
$(SERIAL_SRC)/basic_iserializer.o \
$(SERIAL_SRC)/basic_oserializer.o \
$(SERIAL_SRC)/basic_pointer_oserializer.o \
$(SERIAL_SRC)/basic_archive_impl.o

IOSTREAMS_OBJS = \
$(IOSTREAMS_SRC)/file_descriptor.o \
$(IOSTREAMS_SRC)/mapped_file.o \
$(IOSTREAMS_SRC)/zlib.o

GSL_OBJS= ../gsl/sys/infnan.o ../gsl/sys/coerce.o ../gsl/sys/fdiv.o ../gsl/sys/pow_int.o \
../gsl/sys/fcmp.o ../gsl/specfunc/psi.o ../gsl/specfunc/trig.o ../gsl/specfunc/exp.o \
../gsl/specfunc/log.o ../gsl/specfunc/erfc.o ../gsl/specfunc/zeta.o ../gsl/specfunc/elementary.o \
../gsl/specfunc/gamma.o ../gsl/rng/borosh13.o ../gsl/rng/fishman2x.o ../gsl/rng/mt.o ../gsl/rng/rand.o \
../gsl/rng/ranmar.o ../gsl/rng/types.o ../gsl/rng/cmrg.o ../gsl/rng/gfsr4.o ../gsl/rng/r250.o \
../gsl/rng/random.o ../gsl/rng/rng.o ../gsl/rng/uni32.o ../gsl/rng/coveyou.o ../gsl/rng/knuthran2.o \
../gsl/rng/ran0.o ../gsl/rng/randu.o ../gsl/rng/slatec.o ../gsl/rng/uni.o ../gsl/rng/default.o \
../gsl/rng/knuthran.o ../gsl/rng/ran1.o ../gsl/rng/ranf.o ../gsl/rng/taus113.o ../gsl/rng/vax.o \
../gsl/rng/file.o ../gsl/rng/lecuyer21.o ../gsl/rng/ran2.o ../gsl/rng/ranlux.o ../gsl/rng/taus.o \
../gsl/rng/waterman14.o ../gsl/rng/fishman18.o ../gsl/rng/minstd.o ../gsl/rng/ran3.o ../gsl/rng/ranlxd.o \
../gsl/rng/transputer.o ../gsl/rng/zuf.o ../gsl/rng/fishman20.o ../gsl/rng/mrg.o ../gsl/rng/rand48.o \
../gsl/rng/ranlxs.o ../gsl/rng/tt.o ../gsl/randist/nbinomial.o ../gsl/randist/beta.o ../gsl/randist/exponential.o \
../gsl/randist/geometric.o ../gsl/randist/binomial.o ../gsl/randist/poisson.o ../gsl/randist/rdgamma.o \
../gsl/randist/chisq.o ../gsl/randist/gauss.o ../gsl/error.o ../gsl/randist/multinomial.o


OBJS = $(GSL_OBJS) $(SERIAL_OBJS) $(IOSTREAMS_OBJS)

std: _simuPOP_std.so

op: _simuPOP_op.so

la: _simuPOP_la.so

laop: _simuPOP_laop.so

ba: _simuPOP_ba.so

baop: _simuPOP_baop.so

all: std op la laop ba baop

swigpyrun.h:
	$(SWIG) -python -external-runtime swigpyrun.h
	
######################## std
std_OBJS=utility_std.o migrator_std.o population_std.o simulator_std.o terminator_std.o \
individual_std.o mutator_std.o stator_std.o  \
initializer_std.o operator_std.o recombinator_std.o \
mating_std.o  outputer_std.o selector_std.o  tagger_std.o  
_simuPOP_std.so: $(OBJS) $(std_OBJS) simuPOP_std_wrap.o 
	$(LINK) -shared $(OBJS) $(std_OBJS) simuPOP_std_wrap.o -o _simuPOP_std.so
simuPOP_std_wrap.o: simuPOP_std_wrap.cpp $(HEADFILES)
	$(CXX) -g $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_std $(DEF_std) $(CFLAGS) -c $<
simuPOP_std_wrap.cpp: simuPOP_std.i simuPOP_common.i $(HEADFILES)
	$(SWIG) -shadow -c++ -python -o simuPOP_std_wrap.cpp simuPOP_std.i
%_std.o: %_std.cpp %.h simupop_cfg.h utility.h swigpyrun.h
	$(CXX) $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_std $(DEF_std) $(CFLAGS) -c $< -o $@
%_std.cpp: %.cpp
	cp -f $< $@

######################### OP
op_OBJS=utility_op.o migrator_op.o population_op.o simulator_op.o terminator_op.o \
individual_op.o mutator_op.o stator_op.o  \
initializer_op.o operator_op.o recombinator_op.o \
mating_op.o  outputer_op.o selector_op.o  tagger_op.o 
_simuPOP_op.so: $(OBJS) $(op_OBJS) simuPOP_op_wrap.o 
	$(LINK) -shared $(OBJS) $(op_OBJS) simuPOP_op_wrap.o -o _simuPOP_op.so
simuPOP_op_wrap.o: simuPOP_op_wrap.cpp $(HEADFILES)
	$(CXX) -g $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_op $(DEF_op) $(CFLAGS) -c $<
simuPOP_op_wrap.cpp: simuPOP_op.i simuPOP_common.i $(HEADFILES)
	$(SWIG) $(DEF_op) -shadow -c++ -python -o simuPOP_op_wrap.cpp simuPOP_op.i
%_op.o: %_op.cpp %.h simupop_cfg.h utility.h swigpyrun.h
	$(CXX) $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_op $(DEF_op) $(CFLAGS) -c $< -o $@
%_op.cpp: %.cpp
	cp -f $< $@


######################### LA
la_OBJS=utility_la.o migrator_la.o population_la.o simulator_la.o terminator_la.o \
individual_la.o mutator_la.o stator_la.o  \
initializer_la.o operator_la.o recombinator_la.o \
mating_la.o  outputer_la.o selector_la.o  tagger_la.o 
_simuPOP_la.so: $(OBJS) $(la_OBJS) simuPOP_la_wrap.o 
	$(LINK) -shared $(OBJS) $(la_OBJS) simuPOP_la_wrap.o -o _simuPOP_la.so
simuPOP_la_wrap.o: simuPOP_la_wrap.cpp $(HEADFILES)
	$(CXX) -g $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_la $(DEF_la) $(CFLAGS) -c $<
simuPOP_la_wrap.cpp: simuPOP_la.i simuPOP_common.i $(HEADFILES)
	$(SWIG) $(DEF_la) -shadow -c++ -python -o simuPOP_la_wrap.cpp simuPOP_la.i
%_la.o: %_la.cpp %.h simupop_cfg.h utility.h swigpyrun.h
	$(CXX) $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_la $(DEF_la) $(CFLAGS) -c $< -o $@
%_la.cpp: %.cpp
	cp -f $< $@
	
######################### LAOP
laop_OBJS=utility_laop.o migrator_laop.o population_laop.o simulator_laop.o terminator_laop.o \
individual_laop.o mutator_laop.o stator_laop.o  \
initializer_laop.o operator_laop.o recombinator_laop.o \
mating_laop.o  outputer_laop.o selector_laop.o  tagger_laop.o 
_simuPOP_laop.so: $(OBJS) $(laop_OBJS) simuPOP_laop_wrap.o 
	$(LINK) -shared $(OBJS) $(laop_OBJS) simuPOP_laop_wrap.o -o _simuPOP_laop.so
simuPOP_laop_wrap.o: simuPOP_laop_wrap.cpp $(HEADFILES)
	$(CXX) -g $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_laop $(DEF_laop) $(CFLAGS) -c $<
simuPOP_laop_wrap.cpp: simuPOP_laop.i simuPOP_common.i $(HEADFILES)
	$(SWIG) $(DEF_laop) -shadow -c++ -python -o simuPOP_laop_wrap.cpp simuPOP_laop.i
%_laop.o: %_laop.cpp %.h simupop_cfg.h utility.h swigpyrun.h
	$(CXX) $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_laop $(DEF_laop) $(CFLAGS) -c $< -o $@
%_laop.cpp: %.cpp
	cp -f $< $@

######################### BA
ba_OBJS=utility_ba.o migrator_ba.o population_ba.o simulator_ba.o terminator_ba.o \
individual_ba.o mutator_ba.o stator_ba.o  \
initializer_ba.o operator_ba.o recombinator_ba.o \
mating_ba.o  outputer_ba.o selector_ba.o  tagger_ba.o 
_simuPOP_ba.so: $(OBJS) $(ba_OBJS) simuPOP_ba_wrap.o 
	$(LINK) -shared $(OBJS) $(ba_OBJS) simuPOP_ba_wrap.o -o _simuPOP_ba.so
simuPOP_ba_wrap.o: simuPOP_ba_wrap.cpp $(HEADFILES)
	$(CXX) -g $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_ba $(DEF_ba) $(CFLAGS) -c $<
simuPOP_ba_wrap.cpp: simuPOP_ba.i simuPOP_common.i $(HEADFILES)
	$(SWIG) $(DEF_ba) -shadow -c++ -python -o simuPOP_ba_wrap.cpp simuPOP_ba.i
%_ba.o: %_ba.cpp %.h simupop_cfg.h utility.h swigpyrun.h
	$(CXX) $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_ba $(DEF_ba) $(CFLAGS) -c $< -o $@
%_ba.cpp: %.cpp
	cp -f $< $@


######################### BAOP
baop_OBJS=utility_baop.o migrator_baop.o population_baop.o simulator_baop.o terminator_baop.o \
individual_baop.o mutator_baop.o stator_baop.o  \
initializer_baop.o operator_baop.o recombinator_baop.o \
mating_baop.o  outputer_baop.o selector_baop.o  tagger_baop.o 
_simuPOP_baop.so: $(OBJS) $(baop_OBJS) simuPOP_baop_wrap.o 
	$(LINK) -shared $(OBJS) $(baop_OBJS) simuPOP_baop_wrap.o -o _simuPOP_baop.so
simuPOP_baop_wrap.o: simuPOP_baop_wrap.cpp $(HEADFILES)
	$(CXX) -g $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_baop $(DEF_baop) $(CFLAGS) -c $<
simuPOP_baop_wrap.cpp: simuPOP_baop.i simuPOP_common.i $(HEADFILES)
	$(SWIG) $(DEF_baop) -shadow -c++ -python -o simuPOP_baop_wrap.cpp simuPOP_baop.i
%_baop.o: %_baop.cpp %.h simupop_cfg.h utility.h swigpyrun.h
	$(CXX) $(INCLUDES) -DSIMUPOP_MODUE=simuPOP_baop $(DEF_baop) $(CFLAGS) -c $< -o $@
%_baop.cpp: %.cpp
	cp -f $< $@

#########################
%.o: %.c
	$(CC) $(INCLUDES) $(CFLAGS) -c $< -o $@
%.o: %.c
	$(CC) $(INCLUDES) $(CFLAGS) -c $< -o $@
%.o: %.cpp 
	$(CXX) $(INCLUDES) $(CFLAGS) -c $< -o $@
clean:
	$(RM) _simuPOP_*.so *_wrap.cpp core.* swigpyrun.h *.o *_std.cpp *_op.cpp *_la.cpp *_laop.cpp *_ba.cpp *_baop.cpp
