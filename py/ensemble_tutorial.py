#!/usr/bin/env python2

# @page TutorialPythonEnsemble
#
# In this tutorial, we will see how to generate test cases and save/load them from files.
#
# The first step is always to import swarmng module, make sure that it is in PYTHONPATH.
import swarmng
#
#
# For a trivial test case, we create an ensemble of systems with star at the origin and planets orbiting in
# in perfect circular orbits.

# The parameters are as follows:
# * `nsys` : number of systems in the ensemble
# * `nbod` : number of bodies in the ensemble
# * `spacing_factor` : ratio of the radiuses of orbits of two consecutive planets.
# * `planet_mass` : ratio of mass of the planet to the sun.
#
def make_test_case(nsys = 16, nbod = 3 , spacing_factor = 1.4, planet_mass = 0.001, ejection_factor = 1, seed = None):
  random.seed(seed)
  d = swarmng.DefaultEnsemble.create(nbod,nsys)
  for i in range(0,d.nsys):
      s = d[i]
      s.id = i
      s.set_active()
      s.time = 0
      s[0].pos = [ 0, 0, 0 ]
      s[0].vel = [ 0, 0, 0 ]
      s[0].mass = 1
      fill(s.attributes, 0)
      for j in range(0,d.nbod):
        fill(s[j].attributes, 0)

      for j in range(1,d.nbod):
          r = spacing_factor ** (j-1)
          v = sqrt(1/r) * ejection_factor
          phi = random.uniform(0,2*pi)
          s[j].pos = [  r*cos(phi), r*sin(phi), 0 ]
          s[j].vel = [ -v*sin(phi), v*cos(phi), 0 ]
          s[j].mass = planet_mass 
  return d;

if __name__ == '__main__':
  ens = make_test_case()
#
# You can easily save an ensemble object to a text file:
  ens.save_to_text("sample.txt")
# Alternatively, you can save to a binary format. The binary format is not portable because they
# depend not only on the machine but on the flags that are used when compiling the Swarm-NG library.
# The only
# advantage of binary format is fast loading which makes a big difference when the ensemble
# contains thousands of systems. The main purpose of binary files is to save snapshots and resume
# integration later. For collaborative work or sharing data, please use text format instead.
  ens.save_to_bin("sample.bin")