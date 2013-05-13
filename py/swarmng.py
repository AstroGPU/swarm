import sys
import os

SWARMDIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('lib')
sys.path.append(os.path.join(SWARMDIR,'build/lib'))
sys.path.append(os.path.join(SWARMDIR,'lib'))

from numpy import *
from libswarmng_ext import *

def mkConfig(h):
  """ Convert a regular python hash to a Swarm-NG
   config script"""
  c = Config()
  for k in h :
    c[k] = str(h[k])
  return c

