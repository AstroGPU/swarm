from common import *
from math import sqrt
def norm(*args):
    return sqrt(sum(x * x for x in args))

class ParabolicTest(abstract.IntegrationTest):
    """
    collision happens at exactly t= 1.885638833885 for nbod=3
    """
    cfg = swarmng.config(
            integrator = "irk2",
            nbod       = 3,
            time_step  = .01,
            nogpu      = 0,
            log_writer= "bdb",
            log_output_db = "prab.db",
            deactivate_on_collision= 1
            )
    required_destination_time = 1.885638833885
    destination_time = 21

    def createEnsemble(self):
        nsys = 4
        nbod = 3
        R = 5
        mu = 1
        mass_planet = 0.001


        ens = swarmng.DefaultEnsemble.create(nbod,nsys)
        for s in ens:
            s.time = 0
            s.id = 0
            s.state = 0
            s[0].pos = [ 0, 0, 0 ]
            s[0].vel = [ 0, 0, 0 ]
            s[0].mass = 1
            s[0].attributes[0] = 1
            for j in range(1,ens.nbod):
                x = ((j % 2)*2-1)*10*((j+1)/2)
                y = x * x / 4 / R - R
                vmag = sqrt(2*mu/norm(x,y))
                vdirx = (x/abs(x))/norm(1,x/2/R)
                vdiry = abs(x)/2/R / norm(1,x/2/R)
                s[j].pos = [ x, y, 0 ]
                s[j].vel = [ -vmag*vdirx, -vmag*vdiry, 0 ]
                s[j].mass = mass_planet
                s[j].attributes[0] = 1e-3
        return ens


    def examine(self):
        rs = self.ref[0]
        s = self.ens[0]
        print rs[0].pos, rs[1].pos, rs[2].pos, rs[1].vel, rs[2].vel
        print s[0].pos, s[1].pos, s[2].pos
        print s.time
        max_deltaE = swarmng.find_max_energy_conservation_error(self.ens,self.ref)
        print("Max energy conservation error %g" % max_deltaE)

        #x1,y1,z1 = s[1].pos
        #x2,y2,z2 = s[2].pos
        #self.assertEqual(s.state, -1)
        #self.assertLess(norm(x2+x1,y2+y1,z2+z1),0.011)
        #self.assertAlmostEqual(s.time,self.required_destination_time,2)


