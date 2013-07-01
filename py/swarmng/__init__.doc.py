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
#  @arg \c cfg : an object of type @ref swarmng.Config, for inline creation
#  of Config objects look at @ref swarmng.config.
#
#  For more info on what configuration options are available for swarmng.init
#  refer to @ref Configuration page.
#
#  If you are trying to Swarm on a system without GPUs
#  try setting nogpu = 1 in the configuration.
#
def init(cfg): pass

##  Generate an trivial ensemble with planets in circular orbits
# 
#   @arg @c cfg : a @ref Config object with properties for creating the ensemble
#   only `nsys` and `nbod` are mandatory.
#   
#   Returns a @ref DefaultEnsemble populated with generated systems.
#
#   This function is a wrapper for \ref swarm.generate_ensemble. For more details refer to it.
def generate_ensemble(cfg): pass

## Synchronize all CUDA kernels
#
# This function is rarely used after 
# running an integration to make sure that
# the data structures have been updated.
#
# Note: this is an wrapper for \c cudaThreadSynchronize from CUDA runtime.
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
#  To construct an object of this type easily refer to @ref config
#
#  For a complete list of configuration options used in Swarm refer to @ref Configuration
class Config:
  
  ## Load a Config object from a file name. 
  # 
  # @arg @c filename : path to a text file in swarm config format. The file format is similar to INI format.
  @staticmethod
  def load(filename): pass
  
## A planetary system within an ensemble. 
#
#  To obtain an instance of System class you should index into an ensemble.
#  
#  This class can be treated as a collection
#  of bodies. The bodies can be accessed using
#  brackets or a for loop.
class System:
  ## Time of the system (floating point value)
  time = property
  ## Unique integer identifier for the system
  id = property
  ## Total kinetic and potential energy of the planetary system (floating point value)
  total_energy = property
  ## Integer value representing the state of the system
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

  ##  Position: a list of 3 floating point values: [ x,y,z ]
  pos = property

  ##  Velocity: a list of 3 floating point values: [ vx,vy,vz ]
  vel = property
  
  ## Distance to (0,0,0) (floating point value)
  def distance_to_origin(self): pass

  ## Magnitude of the velocity (floating point value)
  def speed(self):pass

  ## Representetive for one component (x,y,z) of the object
  # Contains position and velocity for that component.
  class Components:
    ## Position for the specified component (single floating point value)
    pos = property
    ## Velocity for the specified property (single floating point value)
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
#  To create an ensemble refer to DefaultEnsemble.create or DefaultEnsemble.load_from_bin or DefaultEnsemble.load_from_text
class Ensemble:
  ## Number of systems in the ensemble (integer value)
  nsys = property
  ## Number of bodies per system (integer value)
  nbod = property
  ## Get the ith system as a @ref System object.
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
  #  @arg @c number_of_bodies : number of bodies per system
  #  @arg @c number_of_systems: total number of systems in the ensemble.
  #
  #  Note that an ensemble is not resizable and all
  #  systems have the same number of bodies.
  def create(number_of_bodies,number_of_systems):pass
  ## Save a binary representation of the whole ensemble to a file
  # @arg @c fileName : name of the file to save the contents to
  def save_to_bin(self, fileName):pass
  ## Save a textual representation of the whole ensemble to a file
  # @arg @c fileName : name of the file to save the contents to
  def save_to_text(self, fileName):pass
  ## @static
  #  Load an ensemble from a text file
  # 
  #  returns a @ref DefaultEnsemble
  def load_from_text(fileName):pass
  ## @static
  #  Load an ensemble from a binary file
  #
  #  returns a @ref DefaultEnsemble
  def load_from_bin(fileName):pass
  
## An ODE integration algorithms
#
# The different implementations of ODE integration methods
# provide the same interface, this class is abstract and
# the create method loads a specific implementation based
# on the configuration that is passed to the method.
class Integrator: 
  ## The ensemble data structure to operate on, of type @ref Ensemble
  ensemble = property
  ## All of the systems will be integrated to this
  # specified time. (scalar floating point value)
  destination_time = property
  ## @static
  #  Create an integrator object from the configuration specified
  #  the only mandatory items is an identifier for the
  #  integrator, namely the value 'integrator'.
  #
  #  @arg @c cfg : configuration object of type @ref swarmng.Config for selecting
  #  and configuring the integrator.
  #
  #  For more information refer to \ref swarm.integrator.create
  def create(cfg):pass
  ## Run the integration up to the specified @ref destination_time
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
  # defined in Integrator.
  ensemble = property
  
  ## Same as @ref Integrator.create, but returns an instance of GpuIntegrator
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
#
# @arg ref : Reference ensemble, possibly from the initial conditions, (should be of type @ref swarmng.Ensemble)
# @arg ens : Ensemble to check (should be of type @ref swarmng.Ensemble)
# Returns one number, that is the maximum of energy conservation error.
# Energy conservation error for a system is the difference in total energy compared
# to the reference system normalized by the amount of energy in the reference system.
# 
def find_max_energy_conservation_error(ref,ens):pass
