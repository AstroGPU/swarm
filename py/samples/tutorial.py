# Header
# header 


# @page TutorialPython Beginner Python Tutorial
# First thing is to import swarmng. It is located in <SOURCE DIRECTORY>/py/swarmng.py.
# It takes care of loading the libswarmng_ext shared object. Please run the scripts
# from your build directory so swarmng.py can find your libraries.
import swarmng
import os

# swarm functions cannot take hashes for config. So 
# we have to create the special config object (a C++ std::map)
# using the mkConfig function
# This first line initializes Swarm-NG. 'nogpu' option
# allows us to run this script on a machine without GPU
swarmng.init(swarmng.mkConfig({ 'nogpu' : 1 }))

# Now we can load an ensemble from a text file. The ensemble file is in the swarm source directory:
ensemble_filename = os.path.join(swarmng.SWARMDIR , 'test/integrators/test.3.in.txt')
# We use a simple call that loads the file and returns the data structure
ref = swarmng.DefaultEnsemble.load_from_text( ensemble_filename )

# We can also look into the data structure and view the values
for s in range(ref.nsys):
  for b in range(ref.nbod):
    print("pos: %s\tvel:%s" % (ref[s][b].pos, ref[s][b].vel))

# We can also modify the values
ref[0][0].pos = [ 1.0, 2.0 , 5.0 ]

# We can easily make a copy of the whole data structure in this example
# we need to have two ensembles: reference (ref) and a working copy (ens)
ens = ref.clone()

# Now we have to select our integrator, swarmng.Integrator.create function
# instantiates an integrator for us.
# In this case, we use a CPU based implementation of Hermite PEC2 integrator(hermite_cpu)
# To see what integrators are available, 'run bin/swarm -p' from command line.
integrator_cfg = {
  'integrator'       : 'hermite_cpu',
  'time_step'        :  0.001       ,
  }
integ = swarmng.Integrator.create(swarmng.mkConfig(integrator_cfg))

# we set up the parameters to the integrator
integ.ensemble = ens
integ.destination_time = 1.0

# Now do the integration
integ.integrate()


# A very common task is to check the energy conservation 
# the function swarmng.find_max_energy_conservation_error does this
# task for us, given two ensembles.
max_deltaE = swarmng.find_max_energy_conservation_error( ens, ref)
print("Max energy conservation error %g" % max_deltaE)

