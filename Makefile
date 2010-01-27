# Source user overrides
-include Makefile.user

#
# Note: you should override these in Makefile.user
#
CCUDA?=/opt/cuda/bin/nvcc -arch=sm_13
CXX?=g++
CCUDAFLAGS?=
CXXFLAGS?=-g -O0 -I /opt/cuda/include -I ./src
LDFLAGS?=-L /opt/cuda/lib64

LIBPEYTON=src/astro/BinaryStream.o
OBJECTS= src/swarmlib.o src/swarmlog.o src/swarm.cu_o src/cux/cux.o $(LIBPEYTON)
EXEOBJECTS = src/swarm.o src/swarmdump.o
ALLOBJECTS = $(OBJECTS) $(EXEOBJECTS)
CUDA_DEPS=src/swarmlib.cu src/swarmlog.h

all: bin/swarm bin/swarmdump

test: all
	(cd run && (test -f data.0 || ../bin/easyGen.py) && ../bin/swarm && ../bin/swarmdump)

clean:
	rm -f $(ALLOBJECTS) *.linkinfo src/*/*.o src/swarm.cu bin/swarmdump bin/swarm integrators/*/*.o

tidy: clean
	rm -f *~ src/*~ src/astro/*~ src/cux/*~ integrators/*/*~ DEADJOE
	rm -f run/data.* run/observeTimes.dat run/*~ run/output.bin
	rm -f .gitignore~

#### Integrators
include integrators/*/Makefile.mk
####

bin/swarm: $(OBJECTS) src/swarm.o
	$(CXX) -rdynamic -o $@ $(LDFLAGS) -lcuda -lcudart $(OBJECTS) src/swarm.o
bin/swarmdump: $(OBJECTS) src/swarmdump.o
	$(CXX) -rdynamic -o $@ $(LDFLAGS) -lcuda -lcudart $(OBJECTS)  src/swarmdump.o

src/cux/cux.o: src/cux/cuda_rng.h src/cux/cux.h src/cux/cux_lowlevel.h

src/swarm.o: src/swarm.h
bin/swarmdump.o: src/swarm.h src/swarmio.h

src/swarmlib.o: src/swarm.h src/swarmio.h
src/swarmlog.o: src/swarm.h src/swarmio.h src/swarmlog.h
src/swarm.cu: src/swarm.h $(CUDA_DEPS)
	echo "// AUTO-GENERATED FILE. DO NOT EDIT BY HAND!!!" > $@
	./bin/combine_cu_files.sh src/swarmlib.cu integrators/*/*.cu >> $@

%.cu_o:%.cu
	$(CCUDA) -c $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) $(PRECOMPILE) $< -o $@

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $(PRECOMPILE) $< -o $@

