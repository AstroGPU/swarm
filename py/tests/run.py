#!/usr/bin/env python2

import unittest, sys
from os.path import dirname, realpath
sys.path.append(dirname(dirname(realpath(__file__))))

from trivial import *
from parabolic_collision import *

if __name__ == '__main__':
	unittest.main()
