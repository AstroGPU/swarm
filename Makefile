CCUDA=/usr/local/cuda/bin/nvcc
CXX=g++
CCUDAFLAGS=--device-emulation -DTHROW_IS_ABORT
CXXFLAGS=-g -O0 -I /usr/local/cuda/include -I ./src
LDFLAGS=-L /usr/local/cuda/lib

OBJECTS=src/swarm.o src/swarmlib.o src/swarm.cu_o

all: bin/swarm

#### Integrators
include integrators/*/Makefile.mk
####

bench: all
	(cd run && ../bin/easyGen.py && ../bin/swarm)

bin/swarm: $(OBJECTS)
	$(CXX) -rdynamic -o $@ $(LDFLAGS) -lcuda -lcudart $(OBJECTS)

src/swarm.o: src/swarm.h
src/swarmlib.o: src/swarm.h
src/swarm.cu: src/swarm.h $(CUDA_DEPS)
	echo "// AUTO-GENERATED FILE. DO NOT EDIT BY HAND!!!" > $@
	./bin/combine_cu_files.sh src/swarmlib.cu integrators/*/*.cu >> $@

clean:
	rm -f $(OBJECTS) *.linkinfo src/swarm.cu bin/swarm integrators/*/*.o

tidy: clean
	rm -f *~ src/*~ integrators/*/*~ DEADJOE
	rm -f run/data.* run/observeTimes.dat

%.cu_o:%.cu
	$(CCUDA) -c $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) $(PRECOMPILE) $< -o $@

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $(PRECOMPILE) $< -o $@
