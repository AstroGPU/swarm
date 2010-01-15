CCUDA=/usr/local/cuda/bin/nvcc
CXX=g++
CCUDAFLAGS=--device-emulation
IFLAGSCU=-I /usr/local/cuda/include
LDFLAGS=-L /usr/local/cuda/lib

SOURCES=swarm.cpp swarmlib.cpp swarm.cu
OBJECTS=swarm.o swarm.cu_o swarmlib.o

all: swarm

swarm: $(OBJECTS)
	$(CXX) -g -o $@ $(LDFLAGS) -lcuda -lcudart $(OBJECTS) 

swarm.o: swarm.cpp swarm.h
swarmlib.o: swarmlib.cpp swarm.h integrators.h
swarm.cu_o: swarm.cu swarm.h integrators.h

clean:
	rm -f $(OBJECTS) *.linkinfo swarm

%.cu_o:%.cu
	$(CCUDA) -g -O0 -c $(CCUDAFLAGS) $(IFLAGSCU) $(DEBUG) $(PRECOMPILE) $< -o $@

%.o:%.cpp
	$(CXX) -g -O0 -c $(IFLAGSCU) $(PRECOMPILE) $<
