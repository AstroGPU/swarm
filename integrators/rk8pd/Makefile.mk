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

# rk8pd code can't compile, so don't try yet
# If pass pointer for each array of dydx, then too many bytes passed to kernel
# If use 4-d array, then not supported by cuxDevicePtr and cuxDeviceAutoPtr
# Potential workaround:  make a simple array of cuxDevicePtr<double,3>'s ?
#
#LIBSWARM_CUDA+=integrators/rk8pd/rk8pd.cu
#LIBSWARM_SOURCES+=integrators/rk8pd/rk8pd.cpp
