## @package swarmng.logrecord
#  Contains the support functions for accessing the C++ logrecord data strcture

from bsddb3.db import *
from struct import *
from collections import namedtuple
from keplerian import calc_keplerian_for_cartesian

## Caluclate Keplerian coordinates for a planet with respect to a predefined center
#  
#  \arg \c planet : a tuple of structure ((x,y,z),(vx,vy,vz),mass) mass is not used in the calculations
#  \arg \c center : the supposed center of the orbit that planet is orbiting around
#
#  \TODO: change this function to use the C implementation of calc_keplerian_for_cartesian
def keplerian_for_cartesian(planet,center):
    ((x,y,z),(vx,vy,vz),mass) = planet
    ((cx,cy,cz),(cvx,cvy,cvz),cmass) = center
    pv = ( (x-cx,y-cy,z-cz), (vx-cvx,vy-cvy,vz-cvz) )
    m  = mass + cmass
    return calc_keplerian_for_cartesian(pv,m)

## Compute the physical center of mass (weighted average of properties by mass) of an array of planets
#
#  \arg \c bodies an array of structure [((x,y,z),(vx,vy,vz),mass),...]
#
def center_of_mass(bodies):
    cx = cy = cz = cvx = cvy = cvz = cmass = 0
    for ((x,y,z),(vx,vy,vz),mass) in bodies:
        cx += x; cy += y; cz += z;
        cvx += vx; cvy += vy; cvz += vz;
        cmass += mass

    return ((cx/cmass,cy/cmass,cz/cmass),(cvx/cmass,cvy/cmass,cvz/cmass),cmass)

## Take a list and return a new iterable that yields elements from the original list with indices
def with_index(list,starting_index = 0):
    return zip(range(starting_index,starting_index+len(list)),list)


## Data structure equivalent of the logrecord defined in Swarm-NG C++ library
#  This class does not use any C++ code, intead it uses pack/unpack to directly parse
#  the binary format into Python variables.
#
#  Currently, the data structure stores a snapshot of the
#  state of a planetary system.
#
#  The constructor of this class does not do anything. 
#  Static method from_binary should be used to parse binary strings
#  and create instances of this class
#
#  Once instantiated properly, this class has following properties
#  - bodies : List of bodies (planets and star)
#  - time   : Time of the snapshot in AU
#  - sys    : Serial number of the system
#  - flags  : The current state of system, different codes may be used by
#  any software that writes the logrecord
#
#
#  \TODO: make the constructor private
#
class LogRecord:

    ## Data structure for properties of a body (planet or star) in a system
    Body = namedtuple('Body', ['position', 'velocity', 'mass'])


    ## Return the star of a planetary system
    #  helper function for bodies_in_keplerian
    def star(self):
        return self.bodies[0]

    ## Barycenter, a.k.a center of mass, of the planetary system
    #  helper function for bodies_in_keplerian
    def barycenter(self):
        return center_of_mass(self.bodies)

    ## The origin with the mass of the star
    #  helper function for bodies_in_keplerian
    def origin(self):
        return LogRecord.Body((0,0,0),(0,0,0),self.bodies[0].mass)


    ##    Convert logrecord to keplerian coordinates, this method
    #     converts the bodies one-by-one. The coordinates are calculated with
    #     respect to a center. The possible options for a center  are 
    #     l.star(), l.barycenter() and l.origin()
    #     The results can be used in a for loop. This returns triplets as
    #      - Ordinal number of the planet (first planet has number 1)
    #      - Cartesian coordinates of the planet as an object that has properties
    #      position, velocity and mass
    #      - Keplerian coordinates of the planet as an object with properties
    #      a, e, i, O, w, M
    def bodies_in_keplerian(self,center):
        """
        """
        bs = with_index(self.bodies)
        for i,b in bs[1:]:
            yield i,b,keplerian_for_cartesian(b,center)

    ## Iterator for planets properties. For each planet, a triple of index,carteisan coordinates, keplerian coordinates 
    #  is returned. The keplerian coordinates are computed with
    #  Jacobi semantics, meaning that center for every planet is 
    #  the center of mass of the interior planetory system (star and all interior planets)
    def bodies_in_keplerian_jacobi(self):
        bs = with_index(self.bodies)
        for i,b in bs[1:]:
            center = center_of_mass(self.bodies[0:i])
            yield i,b,keplerian_for_cartesian(b,center)


    ## Parse a binary string representing a C++ logrecord struct
    #  and return a LogRecord object. Only snapshot logrecords are
    #  supported at the moment
    @staticmethod
    def from_binary(s):
        (msgid, length) = unpack('ii',s[0:8]) # msgid, len
        l = LogRecord()
        l.msgid = msgid
        if msgid == 1 :
            (time,sys,flags,nbod) = unpack('diii',s[8:28])
            bodies = []
            for b in range(1,nbod+1):
                (x,y,z,vx,vy,vz,mass,body_id) = \
                        unpack('dddddddi4x',s[(32+(b-1)*64):(32+(b)*64)])
                bodies.append(LogRecord.Body((x,y,z),(vx,vy,vz) ,mass))
            l.time = time
            l.sys = sys
            l.flags = flags
            l.bodies = bodies
            return l
        else:
            return None
          
    def as_map(self):
        return {'time':self.time, 'sys':self.sys,
        'flags':self.flags, 'bodies':self.bodies }

    ## String representation of the planetary system
    #  the representation is shown like a Hash.
    def __repr__(self):
        return self.as_map().__repr__();


