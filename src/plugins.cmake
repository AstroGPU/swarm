## All the plugins are specified here. To add a new plugin, one should add
## one line with the file name of the main file, the ID a boolean value and a description text
## The line looks like:
### ADD_PLUGIN(path/to/main_file.cu  Plugin_ID  TRUE  "Plug-in Description")
## If you substitute TRUE with FALSE then the plugin will be disabled by default
## However, the user can enable the plugin from the CMake GUI without changing 
## the source code

# CPU plugins
ADD_PLUGIN(plugins/hermite_cpu.cpp Hermite_CPU TRUE "Hermite CPU Integrator")
if(OPENMP_FOUND)
	ADD_PLUGIN(plugins/hermite_omp.cpp Hermite_OMP TRUE "Hermite OpenMP Integrator")
endif()

# GPU Integrators
ADD_PLUGIN(plugins/hermite.cu Hermite TRUE  "Hermite w/ Fixed Time step GPU Integrator")
ADD_PLUGIN(plugins/hermite_adap.cu Hermite_Adap TRUE  "Hermite w/ Adaptive Time step GPU Integrator")
ADD_PLUGIN(plugins/rkck_adaptive.cu RKCK_Adaptive TRUE  "Runge-Kutta Adaptive time step Integrator")
ADD_PLUGIN(plugins/rkck_fixed.cu    RKCK_Fixed    TRUE  "Runge-Kutta Fixed time step Integrator")

# Propagators
ADD_PLUGIN(plugins/mvs.cu MVS TRUE  "Mixed Variable Symplectic Integrator")

# Experimental
ADD_PLUGIN(plugins/euler.cu Euler FALSE "Euler Integrator")
ADD_PLUGIN(plugins/verlet.cu Verlet FALSE "Verlet Integrator")
ADD_PLUGIN(plugins/midpoint.cu Midpoint FALSE "Midpoint Integrator")


## Writer plugins

if(BDB_FOUND)
	INCLUDE_DIRECTORIES(${BDB_INCLUDE_DIR})
	ADD_PLUGIN(swarm/log/bdb_writer.cpp BDB_Writer TRUE "Berkeley DB writer")
endif()

ADD_PLUGIN(swarm/log/binary_writer.cpp Binary_Writer TRUE "Binary file writer")
ADD_PLUGIN(swarm/log/host_array_writer.cpp Host_Array_Writer TRUE "Writer to the host arrays")
