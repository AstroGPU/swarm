#!/usr/bin/env python2
# -*- coding: utf8 -*-
## @file range_type.py Support routines for the defining ranges to use with queries.

## @package swarmng.range_type
#  Classes for defining arbitrary ranges easily in Python
#

import argparse

SINGLE, INTERVAL, UNIVERSAL = range(3)
## Simple data structure to specify a range and test against it
#  This is different from Python range object since it has methods
#  for testing if a number is in range and also supports infinite ranges.
#
#  To create a Range object: use one of the following static methods:
#    - Range.single
#    - Range.interval
#    - Range.universal
#
#  One of the main uses of Range is to test if it contains a number in the range
#  @code{.py}
#  >>> r = Range.interval(5,8)
#  >>> r.contains(7)
#  true
#  >>> r.contains(9)
#  false
#  @endcode{.py}
#
#  A range object can be iterated on just like the Python range object
#  @code{.py}
#  >>> for i in Range.interval(1,10):
#  ...  print(i)
#  1
#  2
#  3
#  4
#  5
#  6
#  7
#  8
#  9
#  10
#  @endcode
#  Note that these ranges contain the upper bound.
#
#
#  \TODO: make the constructor private
class Range:
    ## Create a range that contains only one number
    #  The upper and lower bound are set to \c s.
    @staticmethod
    def single(s):
        x = Range()
        x.type = SINGLE
        x.singleton = s
        return x
      
    ## Create a range from a tuple. 
    # This function can take one tuple as an argument
    # or two arguments.
    # @code{.py}
    # >>> t = (3,8)
    # >>> Range.interval(t)
    # >>> Range.interval(10,20)
    # @endcode
    @staticmethod
    def interval(t, u = None):
        if u != None :
          tpl = (t,u)
        else:
          tpl = t
          
        x = Range()
        x.type = INTERVAL
        x.interval = t
        return x
    ## Create a range that contains everything.
    #  The contain method for this object always returns true.
    #  The upper and lower for this object are undefined.
    @staticmethod
    def universal():
        x = Range()
        x.type = UNIVERSAL
        return x

    def isInterval(self):
        return self.type == INTERVAL
    def isSingle(self):
        return self.type == SINGLE
    def isUniversal(self):
        return self.type == UNIVERSAL


    ## Test that the range contains the number \c i
    def contains(self,i):
        if(self.type == SINGLE):
            return i == self.singleton
        elif(self.type == INTERVAL):
            a,b = self.interval
            return a <= i <= b
        elif(self.type == UNIVERSAL):
            return True

    ## Create an iterator limiting this range to lower and upper.
    #  this is very useful for cases where UNIVERSAL ranges are
    #  allowed, because the UNIVERSAL range will be clamped to
    #  the provided lower and uppper bound.
    def with_limits(self,lower,upper):
        if(self.type == SINGLE):
            if(lower <= self.single <= upper):
                yield self.singleton
        elif(self.type == INTERVAL):
            a,b = self.interval
            for i in range(max(a,lower),min(upper,b)+1):
                yield i
        else:
            for i in range(lower,upper+1):
                yield i



    def __iter__(self):
        if(self.type == SINGLE):
            yield self.singleton
        elif(self.type == INTERVAL):
            a,b = self.interval
            for i in range(a,b+1):
                yield i
        else:
            raise StopIteration

    ## Lower bound for the range.
    #  WARNING it is undefined for Universal range type.
    def lower(self):
        if(self.type == SINGLE):
            return self.singleton
        elif(self.type == INTERVAL):
            return self.interval[0]
        else:
            raise ArgumentError

    ## Upper bound for the range.
    #  WARNING it is undefined for Universal range type.
    def upper(self):
        if(self.type == SINGLE):
            return self.singleton
        elif(self.type == INTERVAL):
            return self.interval[1]
        else:
            raise ArgumentError

    ## Lower/Upper bound as a tuple.
    #  WARNING it is undefined for Universal range type.
    def ulPair(self):
        if(self.type == SINGLE):
            return (self.singleton,self.singleton)
        elif(self.type == INTERVAL):
            return self.interval
        else:
            raise ArgumentError

    def __repr__(self):
        if(self.type == SINGLE):
            return "[ " + repr(self.singleton) + "]"
        elif(self.type == INTERVAL):
            a,b = self.interval
            return "[ " + repr(a) + ".." + repr(b) + "]"
        elif(self.type == UNIVERSAL):
            return "[ .. ]"

class RangeType(object):
    def __init__(self,type=int):
        self.type = type
    def __call__(self, string):
        x = string.split("..")
        if (len(x) == 2):
            return Range.interval((self.type(x[0]),self.type(x[1])))
        elif( len(x) == 1):
            return Range.single(self.type(x[0]))
        else:
            raise argparse.ArgumentTypeError("Invalid range '%s' valid range should be like '1..2'" % string) 
