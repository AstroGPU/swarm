#!/usr/bin/env python2
# -*- coding: utf8 -*-
## @file collision_course.py Testing collision detection module by putting planets on the same orbit in opposite directions
from common import *
from math import sqrt
from random import uniform

def sphericalToCartesian(r,theta, phi = 0):
  return [ r*cos(theta)*cos(phi), r*sin(theta)*cos(phi), r*sin(phi) ]


class CollisionCourseTest(abstract.IntegrationTest):

  cfg = swarmng.config(
    integrator="hermite_cpu",
    time_step_factor=0.017,
    min_time_step=0.001,
    time_step=0.0001,
    destination_time=1,
    deactivate_on_collision=1,
    collision_radius=.05
    )
  
  def createEnsemble(self):
    nsys = 16
    nbod = 4
    orbit_radius = 1.0
    orbit_velocity = sqrt(1/orbit_radius)
    
    ens = swarmng.DefaultEnsemble.create(nbod,nsys)

    for i,s in enumerate(ens):
      # Trivial attributes of the system
      s.id = i
      s.time = 0
      s.set_active()
      
      # select a random phase to put the first body
      phi_base = uniform(0,2*pi)
      
      for j,b in enumerate(s):
        if(j == 0):
          # Set the central body as a large mass stationary object
          b.pos = [0, 0, 0]
          b.vel = [0, 0, 0]
          b.mass = 1.0
        else:
          ## distribute the planets uniformly around the orbit by different phases
          phi = phi_base + (j-1)*2*pi/(ens.nbod-1) 
          b.mass = 0.001
          b.pos = sphericalToCartesian(orbit_radius, phi) 
          b.vel = sphericalToCartesian( -1**j * orbit_velocity, phi+ pi/2)
                
      return ens
    
  def examine(self):
    for s in self.ens:
      self.assertEqual(s.state, -1)
    
    
    
