import argparse

SINGLE, INTERVAL, UNIVERSAL = range(3)
class Range:
    @staticmethod
    def single(s):
        x = Range()
        x.type = SINGLE
        x.singleton = s
        return x
    @staticmethod
    def interval(tuple):
        x = Range()
        x.type = INTERVAL
        x.interval = tuple
        return x
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


    def contains(self,i):
        if(self.type == SINGLE):
            return i == self.singleton
        elif(self.type == INTERVAL):
            a,b = self.interval
            return a <= i <= b
        elif(self.type == UNIVERSAL):
            return True

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

    def lower(self):
        if(self.type == SINGLE):
            return self.singleton
        elif(self.type == INTERVAL):
            return self.interval[0]

    def upper(self):
        if(self.type == SINGLE):
            return self.singleton
        elif(self.type == INTERVAL):
            return self.interval[1]

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
