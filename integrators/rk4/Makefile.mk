#
# Makefile for your integrator. It will be included from the main makefile
# residing in the source directory. You must therefore specify paths w.r.t.
# the source directory, not the local one!
#
# RULES:
#    1.) C++ sources must end in .cpp; CUDA sources must end in .cu
#    2.) Add C++ sources to LIBSWARM_SOURCES
#    3.) Add CUDA sources to LIBSWARM_CUDA
#    4.) Do not add include files (the main makefile will do that)
#

LIBSWARM_CUDA+=integrators/rk4/rk4.cu
LIBSWARM_SOURCES+=integrators/rk4/rk4.cpp
