#
# Makefile for your integrator. It will be included from the main makefile
# residing in the source directory. You must therefore specify paths w.r.t.
# the source directory, not the local one!
#
# RULES:
#    1.) There must be one and only one .cu file.
#    2.) You must append the name of this file, together with any files it
#        may depends, to CUDA_DEPS makefile variable.
#    3.) If defining additional source files, add them to
#        OBJECTS variable, _with full paths!_. Any dependency rules
#        they require must also refer to them using full paths.
#

CUDA_DEPS+=integrators/mvs/mvs.cu integrators/mvs/mvs.h

OBJECTS+=integrators/mvs/mvs.o

integrators/mvs/mvs.o : integrators/mvs/mvs.h
