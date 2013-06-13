from math import *
import random
import swarmng
import unittest
from os import system,path

TESTDIR = path.dirname(path.realpath(__file__))

class abstract:
    class IntegrationTest(unittest.TestCase):
        cfg = None
        destination_time = 1.0
        def runTest(self):
            swarmng.init(self.cfg)
            integ = swarmng.Integrator.create( self.cfg )
            self.ref = self.createEnsemble()
            self.ens = self.ref.clone()

            integ.ensemble = self.ens
            integ.destination_time = self.destination_time
            integ.integrate()

            self.examine()
