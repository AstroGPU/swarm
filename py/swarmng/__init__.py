## @package swarmng
#  Python interface to the Swarm-NG library
#  This package loads and verifies the libswarmng_ext.so and
#  adds some API that is implemented in Python. 
#
#  The extensive module for opening log files is also included in this package.
#  

import sys
import os
import imp
from ctypes import CDLL, c_char_p
import platform

# First find the source directory of Swarm
SWARMDIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

# We have some guesses as where to find the library
DEFAULT_PLACES = [ 'lib', os.path.join(SWARMDIR,'build/lib'), os.path.join(SWARMDIR,'lib') ]

# Find and load the so library
_, fn, _ = imp.find_module('libswarmng_ext', sys.path + DEFAULT_PLACES)
so_dll = CDLL(fn)

# Get the functions that allow you to read the Python versions
so_hexversion = so_dll.python_hexversion
so_version    = so_dll.python_version
so_version.restype = c_char_p

# Test the Python version
if so_hexversion() != sys.hexversion :
  print("The swarmng library was build for Python version %s, but you are running %s" % (so_version(), platform.python_version()))
  sys.exit(1)
else:
  imp.load_dynamic('libswarmng_ext', fn)

from libswarmng_ext import *




## Convert hashes and keyword arguments to a Swarm-NG 
#    config object
#  @code{.py}
#  >>> swarmng.config({ 'integrator' : 'hermite_cpu' }, { 'nsys' : 16, 'nbod' : 3 } , nogpu = 1)
#  @endcode
#
def config(*l,**kw):
    c = Config()
    def addHash(h):
        for k,v in h.items():
            c[k] = str(v)

    for h in l:
        addHash(h)
    addHash(kw)
    return c

def mkConfig(h):
  return config(h)
