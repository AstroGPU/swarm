## All the plugins are specified here. To add a new plugin, one should add
## one line with the file name of the main file, the ID a boolean value and a description text
## The line looks like:
### ADD_PLUGIN(path/to/main_file.cu  Plugin_ID  TRUE  "Plug-in Description")
## If you substitute TRUE with FALSE then the plugin will be disabled by default
## However, the user can enable the plugin from the CMake GUI without changing 
## the source code

# CPU integrators
ADD_PLUGIN(integrators/hermite_cpu.cpp Hermite_CPU TRUE "Hermite CPU Integrator")

# GPU Integrators
ADD_PLUGIN(integrators/hermite.cu Hermite TRUE  "Hermite w/ Fixed Time step GPU Integrator")
ADD_PLUGIN(integrators/hermite_adap.cu Hermite_Adap TRUE  "Hermite w/ Adaptive Time step GPU Integrator")
ADD_PLUGIN(integrators/rkck.cu RKCK TRUE  "Runge-Kutta Fixed/Adaptive time step Integrator")

# Propagators
ADD_PLUGIN(propagators/mvs.cu MVS TRUE  "Mixed Variable Symplectic Integrator")

# Experimental
ADD_PLUGIN(propagators/euler.cu Euler FALSE "Euler Integrator")
ADD_PLUGIN(propagators/verlet.cu Verlet FALSE "Verlet Integrator")
ADD_PLUGIN(propagators/midpoint.cu Midpoint FALSE "Midpoint Integrator")

