CCUDA=/usr/local/cuda/bin/nvcc
CXX=g++
CCUDAFLAGS=--device-emulation
IFLAGSCU=-I /usr/local/cuda/include
LDFLAGS=-L /usr/local/cuda/lib

SOURCES=swarm.cpp swarmlib.cpp swarm.cu
OBJECTS=swarm.o swarm.cu_o swarmlib.o

all: swarm

swarm: $(OBJECTS)
	$(CXX) -o $@ $(LDFLAGS) -lcuda -lcudart $(OBJECTS) 

swarm.cpp: swarm.h
swarmlib.cpp: swarm.h integrators.h
swarm.cu: swarm.h integrators.h

clean:
	rm -f $(OBJECTS) *.linkinfo swarm

%.cu_o:%.cu
	$(CCUDA) -c $(CCUDAFLAGS) $(IFLAGSCU) $(DEBUG) $(PRECOMPILE) $< -o $@

%.o:%.cpp
	$(CXX) -c $(IFLAGSCU) $(PRECOMPILE) $<
