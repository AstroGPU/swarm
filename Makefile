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
LIBSWARM_SOURCES=src/astro/BinaryStream.cpp src/swarmlib.cpp src/swarmlog.cpp src/cux/cux.cpp
LIBSWARM_CUDA=src/swarmlib.cu

####
#### Integrator pieces of libswarm
####
include integrators/*/Makefile.mk

###
### Applications
###
APPS+=swarm
swarm_SOURCES=src/swarm.cpp

APPS+=swarmdump
swarmdump_SOURCES=src/swarmdump.cpp

##########################################################
#
#  You shouldn't need to touch anything below this point
#
##########################################################

#
# Defaults: you should override these in Makefile.user
#
CCUDA?=/opt/cuda/bin/nvcc -arch=sm_13
CXX?=g++
CCUDAFLAGS?=
DEVEMU?=
CXXFLAGS?=-g -O0 -I /opt/cuda/include -I ./src
LDFLAGS?=-L /opt/cuda/lib64
INTEGRATORCFG?=integrator.cfg
VERBOSE?=0

# for libswarm
LDFLAGS+=-L./bin
CXXFLAGS+= -I ./src
CXX+=-fPIC

# Link command just adds the required bits to compiler command
LINK=$(CXX) -Wl,-rpath,$(BIN) -rdynamic $(LDFLAGS) -lcuda -lcudart -lswarm 

SWARM_SOURCES := $(foreach app,$(APPS),$($(app)_SOURCES))	# Collect all app sources
SWARM_OBJECTS=$(SWARM_SOURCES:.cpp=.o)				# All app objects
EXE=$(addprefix bin/, $(APPS))			# bin/swarm bin/swarmdump
LIBSWARM_OBJECTS=$(LIBSWARM_SOURCES:.cpp=.o)	# libswarm objects
OBJECTS=$(SWARM_OBJECTS) $(LIBSWARM_OBJECTS)	# all objects
SOURCES=$(SWARM_SOURCES) $(LIBSWARM_SOURCES)	# all sources

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
	@ echo "// AUTO-GENERATED FILE. DO NOT EDIT BY HAND!!!" > $@
	$(GENUI) ./bin/combine_cu_files.sh $(LIBSWARM_CUDA) >> $@

src/autogen_dont_edit.o: src/autogen_dont_edit.cu_o
	$(GENUI) cp src/autogen_dont_edit.cu_o src/autogen_dont_edit.o

bin/libswarm.so: src/autogen_dont_edit.o $(LIBSWARM_OBJECTS)
	$(NVCCUI) $(CCUDA) -Xcompiler -fPIC $(DEVEMU) $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) -shared -o $@ $^

#
# Utilities
#

test-dataset:
	(cd run && (test -f data.0 || ../bin/easyGen.py))

test: all test-dataset
	(cd run && ../bin/swarm $(INTEGRATORCFG) && ../bin/swarmdump)

clean:
	$(CLEANUI) rm -f *.linkinfo $(OBJECTS) $(EXE) $(OBJECTS:.o=.d) bin/libswarm.so src/autogen_dont_edit.* bin/Makefile.d

tidy: clean
	$(TIDYUI) rm -f *~ .*~ src/*~ src/astro/*~ src/cux/*~ integrators/*/*~ DEADJOE run/data.* run/observeTimes.dat run/*~ run/output.bin run/*.bin

info:
	@ echo LIBSWARM_SOURCES=$(LIBSWARM_SOURCES)
	@ echo SWARM_SOURCES=$(SWARM_SOURCES)
	@ echo LIBSWARM_OBJECTS=$(LIBSWARM_OBJECTS)
	@ echo SWARM_OBJECTS=$(SWARM_OBJECTS)
	@ echo SOURCES=$(SOURCES)
	@ echo OBJECTS=$(OBJECTS)
	@ echo APPS=$(APPS)
	@ echo BIN=$(BIN)

#
# Build patterns
#

# Static pattern rule to let make know that 'make swarm' should mean 'make bin/swarm'
$(APPS): %: bin/%

# Executables
bin/Makefile.d: Makefile
	$(GENUI) ./bin/generate_app_makefiles.sh $(APPS) > $@
-include bin/Makefile.d

# CUDA object files
%.cu_o:%.cu
	$(NVCCUI) $(CCUDA) -Xcompiler -fPIC -c $(DEVEMU) $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) $< -o $@

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
	$(DEPUI) $(CCUDA) -M -odir $(dir $<) $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) $< \
		| sed 's,\($$*\)\.o[ :]*,\1.cu_o $@ : ,g' \
		| sed 's,.*/ \\,    \\,g' \
		| sed 's,//,/,g' \
		 > $@


# Include all auto-generated dependencies, unless we're clean-ing or tidy-ing
#  (see http://ftp.gnu.org/old-gnu/Manuals/make-3.79.1/html_node/make_88.html )
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),tidy)
-include $(subst .cpp,.d,$(SOURCES))
-include src/autogen_dont_edit.cu_d
endif
endif
