
## @package swarmng.libswarmng_ext
# Swarm as a python module
#
# The reason for lib prefix is because CMake automatically
# adds lib prefix to the name of the target
#


# This is the documentation for libswarmng_ext
# It is included here because we cannot put it in module.cpp


## Initialize the Swarm library, registering the GPU and Logging subsystem.
#
def init(cfg):
  pass

##  Generate an trivial ensemble with planets in circular orbits
#   Following parameters are used from the cfg
#   - nsys : number of system
#   - nbod : number of bodies
def generate_ensemble(cfg):
  pass

## Specialization of std::map to hold all our configuration attributes
#
class Config:
  pass