CCUDA=/usr/local/cuda/bin/nvcc
CXX=g++
#CCUDAFLAGS=--device-emulation -DTHROW_IS_ABORT
CXXFLAGS=-g -O0 -I /usr/local/cuda/include -I .
LDFLAGS=-L /usr/local/cuda/lib

OBJECTS=swarm.o swarm.cu_o swarmlib.o

all: swarm

#### Integrators
include integrators/*/Makefile.mk
####

swarm: $(OBJECTS)
	$(CXX) -rdynamic -o $@ $(LDFLAGS) -lcuda -lcudart $(OBJECTS)

swarm.o: swarm.h
swarmlib.o: swarm.h
swarm.cu: swarm.h $(CUDA_DEPS)
	echo "// AUTO-GENERATED FILE. DO NOT EDIT BY HAND!!!" > $@
	./combine_cu_files.sh swarmlib.cu integrators/*/*.cu >> $@

clean:
	rm -f $(OBJECTS) *.linkinfo swarm.cu swarm integrators/*/*.o

tidy: clean
	rm -f *~ integrators/*/*~ DEADJOE

%.cu_o:%.cu
	$(CCUDA) -c $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) $(PRECOMPILE) $< -o $@

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $(PRECOMPILE) $< -o $@
