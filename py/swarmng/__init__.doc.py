## @file __init__.doc.py 
# This is the documentation for libswarmng_ext, but because
# all the symbols from libswarmng_ext are included in the swarmng module
# we document it under swarmng package.
#
# The implementation of all of these is in src/python/module.cpp but
# we have to separate the documentation and put it here.
# 
# To read the documentation generated from this file refer to @ref swarmng
# 

## @package swarmng

# The reason for lib prefix is because CMake automatically
# adds lib prefix to the name of the target
#


## Initialize the Swarm library, registering the GPU and Logging subsystem.
#
#  The initialization procedure must be run before
#  any integration can be done. However, data structure
#  manipulation and loading can be done before swarm.init
#
#  If you are trying to Swarm on a system without GPUs
#  try setting nogpu = 1 in the configuration.
#
def init(cfg): pass

##  Generate an trivial ensemble with planets in circular orbits
#   Following parameters are used from the cfg
#   - nsys : number of system
#   - nbod : number of bodies
#   For more details refer to \ref swarm.generate_ensemble
def generate_ensemble(cfg): pass

## Synchronize all CUDA kernels
#
# This function is rarely used after 
# running an integration to make sure that
# the data structures have been updated.
#
# Note: this is an alias to \c cudaThreadSynchronize from CUDA runtime.
def sync(): pass

## Returns Keplerian coordinates (as a list) from position and
# velocity of a body in cartesian coordinates. 
# 
# The returned list has the following structure:
# @code{.py}
# [ a, e, i, O, w, M ]
# @endcode
#
def keplerian_for_cartesian(x,y,z,vx,vy,vz,GM): pass
## Returns cartesian coordinates (position and velocity) of
# an object for the given Keplerian coordinates.
#
# The return value is a list of following structure:
# @code{.py}
# [ x, y, z, vx, vy, vz ]
# @endcode
# 
def cartesian_for_keplerian(a,e,i,O,w,M): pass


## Specialization of std::map to hold all our configuration attributes
#  
#  to construct an object of this type easily refer to \ref swarmng.config
class Config:
  
  ## Load a Config object from a file name. 
  #
  # File format is similar to INI format.
  @staticmethod
  def load(file_name): pass
  
## A planetary system within an ensemble
#  
#  This class can be treated as a collection
#  of bodies. the bodies can be accessed using
#  brackets or a for loop.
class System:
  ## Time of the system
  time = property
  ## Unique identifier for the system
  id = property
  ## Total kinetic and potential energy of the planetary system
  total_energy = property
  ## State of the system
  #   - 0 means active
  #   - 1 means inactive
  #   - -1 means disabled
  #   - other numbers may have special meanings
  #
  state = property
  ## Collection of attributes of the system of type \ref System.Attributes
  attributes = property
  
  ## Number of bodies in the system
  #  Usage
  #  @code{.py}
  #  >>> len(self)
  #  @endcode
  def __len__(self):pass
  ## Equivalent to state == 0
  def is_active(self):pass
  ## Equivalent to state == 1
  def is_inactive(self):pass
  ## Equivalent to state != -1
  def is_enabled(self):pass
  ## Set state = 0
  def set_active(self):pass
  ## Set state = 1
  def set_inactive(self):pass
  ## Set state = -1
  def set_disabled(self):pass
  
  ## Get the ith body of the system
  #  Usage
  #  @code{.py}
  #  >>> self[i]
  #  @endcode
  def __getitem__(self,i):pass


  ## Attributes of the system, it is just
  # a list of floating point values. You 
  # can use [] and len to extract values from it
  # or iterate over it with a for loop.
  class Attribute:
    ## Get the ith value
    #  Usage
    #  @code{.py}
    #  >>> self[i]
    #  @endcode
    def __getitem__(self,i):pass
    ## Set the ith value
    #  Usage
    #  @code{.py}
    #  >>> self[i] = 3
    #  @endcode
    def __setitem__(self,i,v):pass
    ## Length of the list
    #  Usage
    #  @code{.py}
    #  >>> len(self)
    #  @endcode
    def __len__(self):pass
  
## A body object within a system
class Body:
  ## mass of the body (floating point)
  mass = property

  ##  Position: a list of 3 floating point values: x,y,z
  pos = property

  ##  Velocity: a list of 3 floating point values: vx,vy,vz
  vel = property
  
  ## Distance to (0,0,0)
  def distance_to_origin(self): pass

  ## Magnitude of the velocity
  # basically x*x + y*y + z*z
  def speed(self):pass

  ## Representetive for one component (x,y,z) of the object
  # Contains position and velocity for that component.
  class Components:
    ## Position for the specified component
    pos = property
    ## Velocity for the specified property
    vel = property
    
  ## Attributes of the body, it is just
  # a list of floating point values. You 
  # can use [] and len to extract values from it
  # or iterate over it with a for loop.
  class Attribute:
    ## Get the ith value
    #  Usage
    #  @code{.py}
    #  >>> self[i]
    #  @endcode
    def __getitem__(self,i):pass
    ## Set the ith value
    #  Usage
    #  @code{.py}
    #  >>> self[i] = 3
    #  @endcode
    def __setitem__(self,i,v):pass
    ## Length of the list
    #  Usage
    #  @code{.py}
    #  >>> len(self)
    #  @endcode
    def __len__(self):pass
    
 
## Abstract ensemble data structure.
#
#  The ensemble class has all the accessor methods
#  it is only abstract in terms of storage. Because
#  an ensemble can be stored on GPU memory or system memory.
#
#  To create an ensemble refer to DefaultEnsemble.create
class Ensemble:
  ## Number of systems in the ensemble
  nsys = property
  ## Number of bodies per system
  nbod = property
  ## Get the ith system
  #  Usage
  #  @code{.py}
  #  >>> self[i]
  #  @endcode
  def __getitem__(self,i):pass
  ## Number of systems in the ensemble equal to nsys
  #  Usage
  #  @code{.py}
  #  >>> len(self)
  #  @endcode
  def __len__(self):pass

## The default implementation of ensemble data structor
# that stores data in system memory.
class DefaultEnsemble(Ensemble):
  ## @static
  #  Create an ensemble with specified
  #  number of systems and bodies
  #
  #  Note that an ensemble is not resizable and all
  #  systems have the same number of bodies.
  def create(number_of_bodies,number_of_systems):pass
  ## Save a binary representation of the whole ensemble to a file
  def save_to_bin(self, fileName):pass
  ## Save a textual representation of the whole ensemble to a file
  def save_to_text(self, fileName):pass
  ## @static
  #  Load an ensemble from a text file
  def load_from_text(fileName):pass
  ## @static
  #  Load an ensemble from a binary file
  def load_from_bin(fileName):pass
  
## An ODE integration algorithms
#
# The different implementations of ODE integration methods
# provide the same interface, this class is abstract and
# the create method loads a specific implementation based
# on the configuration that is passed to the method.
class Integrator: 
  ## The ensemble data structure to operate on
  ensemble = property
  ## All of the systems will be integrated to this
  # specified time.
  destination_time = property
  ## @static
  #  Create an integrator object from the configuration specified
  #  the only mandatory items is an identifier for the
  #  integrator, namely the value 'integrator'.
  #
  #  For more information refer to \ref swarm.integrator.create
  def create(cfg):pass
  ## Run the integration
  def integrate(self):pass
  
## GPU accelerated integrator
#
#  WARNING: this class is not thoroughly tested
#  for most intents and purposes, just use the regular 
#  Integrator class.
class GpuIntegrator(Integrator):
  ## GPU ensemble that is used for the integration
  #
  # this seems redundant since ensemble property is already
  # defined in Integrator
  ensemble = property
  
  ## Same as Integrator.create, but returns an instance of GpuIntegrator
  def create(cfg):pass
  ## The default integrate method updates the GPU ensemble
  # every time. The core_integrate just launches the kernel.
  def core_integrate(self):pass
  ## Update the ensemble in system memory from GPU ensemble
  def download_ensemble(self):pass
  ## Update the GPU ensemble from the ensemble in system memory.
  def upload_ensemble(self):pass
  
## Compare two ensembles and find the maximum of energy
# conservation error amongst systems.
def find_max_energy_conservation_error(ref,ens):pass
