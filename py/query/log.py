from bsddb3.db import *
from struct import *
from collections import namedtuple
from keplerian import calc_keplerian_for_cartesian


def keplerian_for_cartesian(planet,center):
    ((x,y,z),(vx,vy,vz),mass) = planet
    ((cx,cy,cz),(cvx,cvy,cvz),cmass) = center
    pv = ( (x-cx,y-cy,z-cz), (vx-cvx,vy-cvy,vz-cvz) )
    m  = mass + cmass
    return calc_keplerian_for_cartesian(pv,m)

def center_of_mass(bodies):
    cx = cy = cz = cvx = cvy = cvz = cmass = 0
    for ((x,y,z),(vx,vy,vz),mass) in bodies:
        cx += x; cy += y; cz += z;
        cvx += vx; cvy += vy; cvz += vz;
        cmass += mass

    return ((cx/cmass,cy/cmass,cz/cmass),(cvx/cmass,cvy/cmass,cvz/cmass),cmass)

def with_index(list,starting_index = 0):
    return zip(range(starting_index,starting_index+len(list)),list)

class LogRecord:

    Body = namedtuple('Body', ['position', 'velocity', 'mass'])


    def star(self):
        return self.bodies[0]

    def barycenter(self):
        return center_of_mass(self.bodies)

    def origin(self):
        return LogRecord.Body((0,0,0),(0,0,0),self.bodies[0].mass)

    def bodies_in_keplerian(self,center):
        """
        Convert logrecord to keplerian coordinates, this method
        converts the bodies one-by-one. The coordinates are calculated with
        respect to a center. The possible options for a center  are 
        l.star(), l.barycenter() and l.origin()
        The results can be used in a for loop. This returns triplets as
         - Ordinal number of the planet (first planet has number 1)
         - Cartesian coordinates of the planet as an object that has properties
         position, velocity and mass
         - Keplerian coordinates of the planet as an object with properties
         a, e, i, O, w, M
        """
        bs = with_index(self.bodies)
        for i,b in bs[1:]:
            yield i,b,keplerian_for_cartesian(b,center)

    def bodies_in_keplerian_jacobi(self):
        bs = with_index(self.bodies)
        for i,b in bs[1:]:
            center = center_of_mass(self.bodies[0:i])
            yield i,b,keplerian_for_cartesian(b,center)


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

    def __repr__(self):
        return self.as_map().__repr__();


