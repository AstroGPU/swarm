#CCUDA=/opt/cuda/bin/nvcc
#CXX=g++
#CCUDAFLAGS=--device-emulation -DTHROW_IS_ABORT
#CXXFLAGS=-g -O0 -I /opt/cuda/include -I ./src
#LDFLAGS=-L /opt/cuda/lib64 
CCUDA=/usr/local/cuda/bin/nvcc
CXX=g++
CCUDAFLAGS=--device-emulation -DTHROW_IS_ABORT
CXXFLAGS=-g -O0 -I /usr/local/cuda/include -I ./src
LDFLAGS=-L /usr/local/cuda/lib

LIBPEYTON=src/astro/BinaryStream.o
OBJECTS= src/swarmlib.o src/swarm.cu_o $(LIBPEYTON)
EXEOBJECTS = src/swarm.o src/swarmdump.o
ALLOBJECTS = $(OBJECTS) $(EXEOBJECTS)

all: bin/swarm bin/swarmdump

#### Integrators
include integrators/*/Makefile.mk
####

bin/swarm: $(OBJECTS) src/swarm.o
	$(CXX) -rdynamic -o $@ $(LDFLAGS) -lcuda -lcudart $(OBJECTS) src/swarm.o

src/swarm.o: src/swarm.h

test: all
	(cd run && (test -f data.0 || ../bin/easyGen.py) && ../bin/swarm && ../bin/swarmdump)


bin/swarmdump: $(OBJECTS) src/swarmdump.o
	$(CXX) -rdynamic -o $@ $(LDFLAGS) -lcuda -lcudart $(OBJECTS)  src/swarmdump.o
bin/swarmdump.o: src/swarm.h src/swarmio.h


src/swarmlib.o: src/swarm.h src/swarmio.h
src/swarm.cu: src/swarm.h $(CUDA_DEPS)
	echo "// AUTO-GENERATED FILE. DO NOT EDIT BY HAND!!!" > $@
	./bin/combine_cu_files.sh src/swarmlib.cu integrators/*/*.cu >> $@

clean:
	rm -f $(ALLOBJECTS) *.linkinfo src/astro/*.o src/swarm.cu bin/swarmdump bin/swarm integrators/*/*.o

tidy: clean
	rm -f *~ src/*~ src/astro/*~ integrators/*/*~ DEADJOE
	rm -f run/data.* run/observeTimes.dat

%.cu_o:%.cu
	$(CCUDA) -c $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) $(PRECOMPILE) $< -o $@

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $(PRECOMPILE) $< -o $@

