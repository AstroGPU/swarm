#!/usr/bin/env python2

import unittest, sys
from os.path import dirname, realpath
sys.path.append(dirname(dirname(realpath(__file__))))

from test_ejection import *

if __name__ == '__main__':
	unittest.main()
