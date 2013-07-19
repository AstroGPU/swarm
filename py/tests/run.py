#!/usr/bin/env python2
## @file run.py Runner for all the unittests. Execute this file from the interpreter to run all the unit tests.
#
#  for more information on command-line arguments, look at 
#  <a href="http://docs.python.org/2/library/unittest.html#command-line-interface">Python unittest library documentation</a>.
import unittest, sys
from os.path import dirname, realpath
sys.path.append(dirname(dirname(realpath(__file__))))

from trivial import *
from parabolic_collision import *
from bdb_concurrent import BDBConcurrencyTest

if __name__ == '__main__':
	unittest.main()
