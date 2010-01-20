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

OBJECTS= src/swarmlib.o src/swarm.cu_o  # removed swarm.o since not common to different test program
EXEOBJECTS = src/swarm.o src/swarm_test_hermite_cpu.o src/swarm_test_hermite_gpu.o
ALLOBJECTS = $(OBJECTS) $(EXEOBJECTS)

all: bin/swarm bin/swarm_test_hermite_cpu bin/swarm_test_hermite_gpu

#### Integrators
include integrators/*/Makefile.mk
####

bench: all
	(cd run && ../bin/easyGen.py && ../bin/swarm)

bin/swarm: $(OBJECTS) src/swarm.o 
	$(CXX) -rdynamic -o $@ $(LDFLAGS) -lcuda -lcudart $(OBJECTS) src/swarm.o

src/swarm.o: src/swarm.h

test: all
	(cd run && (test -f data.0 || ../bin/easyGen.py) && ../bin/swarm_test_hermite_gpu)


bin/swarm_test_hermite_cpu: $(OBJECTS) src/swarm_test_hermite_cpu.o 
	$(CXX) -rdynamic -o $@ $(LDFLAGS) -lcuda -lcudart $(OBJECTS)  src/swarm_test_hermite_cpu.o
bin/swarm_test_hermite_cpu.o: src/swarm.h

bin/swarm_test_hermite_gpu: $(OBJECTS) src/swarm_test_hermite_gpu.o 
	$(CXX) -rdynamic -o $@ $(LDFLAGS) -lcuda -lcudart $(OBJECTS)  src/swarm_test_hermite_gpu.o
bin/swarm_test_hermite_gpu.o: src/swarm.h


src/swarmlib.o: src/swarm.h
src/swarm.cu: src/swarm.h $(CUDA_DEPS)
	echo "// AUTO-GENERATED FILE. DO NOT EDIT BY HAND!!!" > $@
	./bin/combine_cu_files.sh src/swarmlib.cu integrators/*/*.cu >> $@

clean:
	rm -f $(ALLOBJECTS) *.linkinfo src/swarm.cu bin/swarm integrators/*/*.o bin/swarm_test_hermite_cpu bin/swarm_test_hermite_gpu

tidy: clean
	rm -f *~ src/*~ integrators/*/*~ DEADJOE
	rm -f run/data.* run/observeTimes.dat

%.cu_o:%.cu
	$(CCUDA) -c $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) $(PRECOMPILE) $< -o $@

%.o:%.cpp
	$(CXX) -c $(CXXFLAGS) $(PRECOMPILE) $< -o $@

