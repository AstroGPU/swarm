#
# swarm build system (mjuric, 2010/01/31)
#
# Requirements:
#	- GNU make
#	- GCC g++
#	- nvcc
#
# - sources local definitions from Makefile.user (if exists)
# - builds src/swarmlib.a from LIBSWARM_SOURCES and LIBSWARM_CUDA
# - includes integrators/*/Makefile.mk files for integrator implementations
# - builds bin/app, where app are read from $(APPS) and built from $(app_SOURCES)
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

# Always append this 
CXXFLAGS+= -I ./src

SWARM_SOURCES := $(foreach app,$(APPS),$($(app)_SOURCES))	# Collect all app sources
SWARM_OBJECTS=$(SWARM_SOURCES:.cpp=.o)				# All app objects
EXE=$(addprefix bin/, $(APPS))			# bin/swarm bin/swarmdump
LIBSWARM_OBJECTS=$(LIBSWARM_SOURCES:.cpp=.o)	# libswarm objects
OBJECTS=$(SWARM_OBJECTS) $(LIBSWARM_OBJECTS)	# all objects
SOURCES=$(SWARM_SOURCES) $(LIBSWARM_SOURCES)	# all sources

# User-interface beautification (useless but cute ;-))
# For full output, define VERBOSE=1 in Makefile.user or run
# make with `env VERBOSE=1 make'
#
ifneq ($(VERBOSE),1)
DEPUI=@ echo "[ DEP  ] $< " &&
CXXUI=@ echo "[ CXX  ] $< " &&
NVCCUI=@ echo "[ NVCC ] $< " &&
LINKUI=@ echo "[ LINK ] $< " &&
ARUI=@ echo "[ AR   ] $< " &&
GENUI=@ echo "[ GEN  ] $< " &&
CLEANUI=@ echo "[ CLEAN ]" &&
TIDYUI=@ echo "[ TIDY  ]" &&
endif

all: $(APPS)

#
# rules for libswarm.a
#

src/swarm.cu: $(LIBSWARM_CUDA)
	@ echo "// AUTO-GENERATED FILE. DO NOT EDIT BY HAND!!!" > $@
	$(GENUI) ./bin/combine_cu_files.sh $(LIBSWARM_CUDA) >> $@

src/libswarm.a: $(LIBSWARM_OBJECTS) src/swarm.cu_o
	$(ARUI) ar rcs $@ $(LIBSWARM_OBJECTS)


#
# Utilities
#

test: all
	(cd run && (test -f data.0 || ../bin/easyGen.py) && ../bin/swarm $(INTEGRATORCFG) && ../bin/swarmdump)

clean:
	$(CLEANUI) rm -f *.linkinfo $(OBJECTS) $(EXE) $(OBJECTS:.o=.d) src/swarmlib.a src/swarm.cu

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

#
# Build patterns
#

# Static pattern rule to let make know that 'make swarm' should mean 'make bin/swarm'
$(APPS): %: bin/%

# Executables
bin/%: src/libswarm.a src/%.o $($_OBJECTS)
	$(LINKUI) $(CXX) -rdynamic -o $@ $(LDFLAGS) -lcuda -lcudart $^ src/libswarm.a

# CUDA object files
%.cu_o:%.cu
	$(CXXUI) $(CCUDA) -c $(DEVEMU) $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) $< -o $@

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

%.cu_d: %.cu
	$(DEPUI) $(CCUDA) -M -odir $(dir $<) $(CCUDAFLAGS) $(CXXFLAGS) $(DEBUG) $< | sed 's,$$*.o,& $@,g' > $@


# Include all auto-generated dependencies, unless we're clean-ing or tidy-ing
#  (see http://ftp.gnu.org/old-gnu/Manuals/make-3.79.1/html_node/make_88.html )
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),tidy)
-include $(subst .cpp,.d,$(SOURCES))
-include src/swarm.cu_d
endif
endif
