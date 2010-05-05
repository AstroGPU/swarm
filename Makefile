#
# swarm build system (mjuric, 2010/01/31)
#
# Requirements:
#	- GNU make
#	- GNU ld
#	- GNU GCC g++
#	- nvcc
#
# - sources local definitions from Makefile.user (if exists)
# - builds bin/swarmlib.so from LIBSWARM_SOURCES and LIBSWARM_CUDA
# - includes integrators/*/Makefile.mk files for integrator implementations
# - builds bin/app, where apps are read from $(APPS) and built from $(app_SOURCES)
# - does the automatic dependency tracking for all .cpp and .cu sources
# - shows simple/pretty output (linux kernel-style), unless VERBOSE=1 is 
#   defined (e.g., `env VERBOSE=1 make')
#
# Assuming everything is set up correctly in Makefile.user, (re)building
# everything should be as easy as running `make'
#

# Source user overrides
-include Makefile.user

###
### libswarm library
###
LIBSWARM_SOURCES=src/astro/BinaryStream.cpp src/astro/MemoryMap.cpp src/astro/Util.cpp src/swarmlib.cpp src/swarmlog.cpp src/swarmio.cpp src/cux/cux.cpp
LIBSWARM_CUDA=src/swarmlib.cu

###
### documentation in asciidoc format
###
DOC_INPUT=docs/README.txt docs/build_system.txt docs/configuration_file.txt docs/for_developers.txt docs/organization.txt docs/performance.txt docs/eventlog.txt docs/snapshotting.txt # asciidoc general documentation
MAN_INPUT=docs/swarm_tutorial_cpu.man docs/swarm_tutorial_gpu.man docs/swarm_tutorial_compare.man docs/swarm_tutorial_benchmark.man docs/swarm_tutorial_montecarlo.man docs/swarm_scatter_demo.man  docs/swarm.man docs/swarmquery.man  docs/swarmquerykeplerian.man docs/swarm_test_energy.man docs/ic_generators.man # asciidoc man page documentation

####
#### Integrator pieces of libswarm
####
include integrators/*/Makefile.mk

###
### Applications
###
APPS+=swarm 
swarm_SOURCES=src/swarm.cpp

APPS+=swarmquery
swarmquery_SOURCES=src/swarmquery.cpp

APPS+=swarmquerykeplerian
swarmquerykeplerian_SOURCES=src/swarmquerykeplerian.cpp

APPS+=swarm_test_energy
swarm_test_energy_SOURCES=src/swarm_test_energy.cpp

APPS+=swarm_tutorial_cpu
swarm_tutorial_cpu_SOURCES=src/tutorials/swarm_tutorial_cpu.cpp 

APPS+=swarm_tutorial_gpu 
swarm_tutorial_gpu_SOURCES=src/tutorials/swarm_tutorial_gpu.cpp 

APPS+=swarm_tutorial_compare 
swarm_tutorial_compare_SOURCES=src/tutorials/swarm_tutorial_compare.cpp 

APPS+=swarm_tutorial_benchmark
swarm_tutorial_benchmark_SOURCES=src/tutorials/swarm_tutorial_benchmark.cpp 

APPS+=swarm_tutorial_montecarlo
swarm_tutorial_montecarlo_SOURCES=src/tutorials/swarm_tutorial_montecarlo.cpp 

# APPS+=swarm_tutorial_montecarlo_rv
# swarm_tutorial_montecarlo_rv_SOURCES=src/tutorials/swarm_tutorial_montecarlo_rv.cpp 

APPS+=swarm_scatter_demo
swarm_scatter_demo_SOURCES=src/scatter_demo/swarm_scatter_demo.cpp

##########################################################
#
#  You shouldn't need to touch anything below this point
#
##########################################################

#
# Defaults: you should override these in Makefile.user
#
CUDAPATH?=/opt/cuda
CCUDA?=$(CUDAPATH)/bin/nvcc -arch=sm_13
CXX?=g++
CCUDAFLAGS?=
CCUDADIAGFLAGS?=-Xcudafe --diag_suppress=subscript_out_of_range -Xcudafe --diag_suppress=partial_override  -Xcudafe --diag_suppress=initialization_not_reachable
DEVEMU?=
CXXFLAGS?=-g -O0 -I $(CUDAPATH)/include -I ./src
LDFLAGS?=-L $(CUDAPATH)/lib64  -L /usr/lib64
INTEGRATORCFG?=integrator.cfg
VERBOSE?=0

# for libswarm
LDFLAGS+=-L./bin
CXXFLAGS+= -I ./src
CXX+=-fPIC

# Link command just adds the required bits to compiler command
LINK=$(CXX) -Wl,-rpath,$(BIN) -rdynamic $(LDFLAGS) -lcuda -lcudart -lswarm -lgsl -lgslcblas -lboost_program_options -lboost_regex

SWARM_SOURCES := $(foreach app,$(APPS),$($(app)_SOURCES))	# Collect all app sources
SWARM_OBJECTS=$(SWARM_SOURCES:.cpp=.o)				# All app objects
EXE=$(addprefix bin/, $(APPS))			# bin/swarm bin/swarm_test_energy
LIBSWARM_OBJECTS=$(LIBSWARM_SOURCES:.cpp=.o)	# libswarm objects
OBJECTS=$(SWARM_OBJECTS) $(LIBSWARM_OBJECTS)	# all objects
SOURCES=$(SWARM_SOURCES) $(LIBSWARM_SOURCES)	# all sources
DOC_OUTPUT=$(DOC_INPUT:.txt=.html)              # html version of asciidoc text documentation
MAN_OUTPUT=$(MAN_INPUT:.man=.html)              # html version of asciidoc man pages

BIN=$(shell pwd)/bin

# User-interface beautification (useless but cute ;-))
# For full output, define VERBOSE=1 in Makefile.user or run
# make with `env VERBOSE=1 make'
#
ifneq ($(VERBOSE),1)
DEPUI=@ echo "[ DEP  ] $@ " &&
CXXUI=@ echo "[ CXX  ] $@ " &&
NVCCUI=@ echo "[ NVCC ] $@ " &&
LINKUI=@ echo "[ LINK ] $@ " &&
ARUI=@ echo "[ AR   ] $@ " &&
GENUI=@ echo "[ GEN  ] $@ " &&
CLEANUI=@ echo "[ CLEAN ]" &&
TIDYUI=@ echo "[ TIDY  ]" &&
endif

all: $(APPS)

#
# rules for libswarm.so
#

src/autogen_dont_edit.cu: $(LIBSWARM_CUDA)
	@ echo "/*! \file autogen_dont_edit.cu" > $@
	@ echo " *  \brief merged cu files for nvcc" >> $@
	@ echo " * " >> $@
	@ echo " * AUTO-GENERATED FILE. DO NOT EDIT BY HAND!!!" >> $@
	@ echo " */" >> $@
	$(GENUI) ./scripts/combine_cu_files.sh $(LIBSWARM_CUDA) >> $@

src/autogen_dont_edit.o: src/autogen_dont_edit.cu_o
	$(GENUI) cp src/autogen_dont_edit.cu_o src/autogen_dont_edit.o

bin/libswarm.so: src/autogen_dont_edit.o $(LIBSWARM_OBJECTS)
	$(NVCCUI) $(CCUDA) -Xcompiler -fPIC $(DEVEMU) $(CCUDADIAGFLAGS) $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) -shared -o $@ $^

#
# Utilities
#

apps:  $(APPS)

clean-doc:
	rm -rf reference; rm -f README.html docs/*.html

clean-test:
	$(CLEANUI) rm -rf test-outputs test-outputs.tar.gz

test: all
	@ ./scripts/run_tests.sh test || $(MAKE) --no-print-directory test-feedback

test-create: all
	@ ./scripts/run_tests.sh create

test-feedback:
	@ $(CCUDA) --version  > test-outputs/nvcc.info 2>&1
	@ $(CXX)   -v         > test-outputs/gcc.info 2>&1
	@ uname -a            > test-outputs/kernel.info 2>&1
	@ dmesg | grep NVIDIA > test-outputs/driver.info 2>&1
	@ rm -f test-outputs.tar.gz && tar czf test-outputs.tar.gz test-outputs
	@ echo
	@ echo Some tests have failed. We\'ve packed the outputs of these tests
	@ echo to test-outputs.tar.gz. Please e-mail it to astrogpu@gmail.com
	@ echo to help us diagnose the problems you\'ve encountered, and 
	@ echo improve the future releases.
	@ echo "                                     Swarm-NG Developers"

.PHONY: test test-ref

benchmark-quick: bin/swarm_tutorial_benchmark
	cd run; rm -f benchmark.out ; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -p 1 -s 3840 | tee benchmark.out

benchmark: bin/swarm_tutorial_benchmark
	rm -f run/benchmark.out ;
	@ echo "# Benchmarking num systems (Please be patient)" |tee run/benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -s   960 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -s  1920 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -s  3840 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -s  7680 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -s 15360 2>> benchmark.out
	@ echo "# Benchmarking num bodys per system (Please be patient)" |tee -a run/benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -s  3840 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 4 -s  3840 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 5 -s  3840 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 6 -s  3840 2>> benchmark.out
	@ echo "# Benchmarking blocksize (Please be patient)" |tee -a run/benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -p 1 -b  16 -s  7680 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -p 1 -b  32 -s  7680 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -p 1 -b  64 -s  7680 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -p 1 -b  96 -s  7680 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -p 1 -b 128 -s  7680 2>> benchmark.out
	@ echo "# Benchmarking preision (Please be patient)" |tee -a run/benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -p 1 -s  3840 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -p 2 -s  3840 2>> benchmark.out
	cd run; ../bin/swarm_tutorial_benchmark -t 1. -n 3 -p 3 -s  3840 2>> benchmark.out
	@ echo "# Thank you for your patience" | tee -a run/benchmark.out

clean: clean-test
	$(CLEANUI) rm -f *.linkinfo $(OBJECTS) $(EXE) $(OBJECTS:.o=.d) bin/libswarm.so src/autogen_dont_edit.* bin/Makefile.d

tidy: clean
	$(TIDYUI) rm -f *~ .*~ src/*~ src/astro/*~ src/cux/*~ integrators/*/*~ docs/*~ DEADJOE run/data.* run/observeTimes.dat run/*~ run/*.bin run/*.idx feedback.* run/benchmark.out
	$(TIDYUI) rmdir bin

info:
	@ echo LIBSWARM_SOURCES=$(LIBSWARM_SOURCES)
	@ echo SWARM_SOURCES=$(SWARM_SOURCES)
	@ echo LIBSWARM_OBJECTS=$(LIBSWARM_OBJECTS)
	@ echo SWARM_OBJECTS=$(SWARM_OBJECTS)
	@ echo SOURCES=$(SOURCES)
	@ echo OBJECTS=$(OBJECTS)
	@ echo APPS=$(APPS)
	@ echo BIN=$(BIN)

doc:    doc-asciidoc doc-doxygen

doc-asciidoc: $(DOC_OUTPUT) $(MAN_OUTPUT)

doc-doxygen: Doxyfile $(SOURCES) 
	doxygen 

feedback:
	@ echo "About to create feedback.tgz which will contain info that may help debug problems with setup and demos"
	./scripts/make_feedback.sh
	@ echo "Done."

#
# Build patterns
#

# Static pattern rule to let make know that 'make swarm' should mean 'make bin/swarm'
$(APPS): %: bin/%

# Executables
bin/Makefile.d: Makefile
	$(GENUI) mkdir -p bin && ./scripts/generate_app_makefiles.sh $(APPS) > $@
-include bin/Makefile.d

# CUDA object files
%.cu_o:%.cu
	$(NVCCUI) $(CCUDA) -Xcompiler -fPIC -c $(DEVEMU) $(CCUDADIAGFLAGS) $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) $< -o $@ 2>&1 | ./scripts/silence_nvcc_warnings.sh

# C++ object files
%.o:%.cpp
	$(CXXUI) $(CXX) -c $(CXXFLAGS) $< -o $@

#
# Dependency generation code, as advised by the GNU make manual:
#   http://www.ipp.mpg.de/~dpc/gmake/make_42.html
# The -MT and -odir pieces are there because we're building files
#   in subdirectories
# The sed piece is to add .d file as a target (so that deps are
#   automatically remade any time any file on which the source
#   depends changes)
#

%.d: %.cpp
	$(DEPUI) $(CXX) -M -MT "$@ $(subst .cpp,.o,$<)" $(CXXFLAGS) $< > $@

# The 1st sed adds cu_d file as a target to be remade if any of the deps change
# The 2nd sed is a workaround for nvcc 2.3 (or gcc?) bug where the file directory is listed as a dependency
# The 3rd sed fixes a problem where multiple slashes (i.e. //) may be present in the target, because $(dir $<) leaves the trailing slash, and nvcc -odir expects there to be none
%.cu_d: %.cu
	$(DEPUI) $(CCUDA) -M -odir $(dir $<) $(CCUDADIAGFLAGS) $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) $< \
		| sed 's,\($$*\)\.o[ :]*,\1.cu_o $@ : ,g' \
		| sed 's,.*/ \\,    \\,g' \
		| sed 's,//,/,g' \
		 > $@

# Rules for generating html version of AsciiDoc documentation
%.html: %.txt
	asciidoc $<

%.html: %.man
	asciidoc $<

# Include all auto-generated dependencies, unless we're clean-ing or tidy-ing
#  (see http://ftp.gnu.org/old-gnu/Manuals/make-3.79.1/html_node/make_88.html )
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),tidy)
-include $(subst .cpp,.d,$(SOURCES))
-include src/autogen_dont_edit.cu_d
endif
endif


