from common import *


##  Testing the monitor stop_on_ejection, some planets are put in orbits
##  very close to one another (1.01 separation) and they are expected to
##  eject out of the system. Since not all of the eject, in the
##  verification stage we check that if the body has ejected, the system
##  should have been deactivated

def fill(list_like,initial_value):
    for i in range(0,len(list_like)):
        list_like[i] = initial_value


def make_test_case(nsys = 16, nbod = 3 , spacing_factor = 1.4, planet_mass = 0.001, ejection_factor = 1, seed = None):
    random.seed(seed)
    d = swarmng.DefaultEnsemble.create(nbod,nsys)
    for i in range(0,d.nsys):
        s = d[i]
        s.id = i
        s.set_active()
        s.time = 0
        s[0].pos = [ 0, 0, 0 ]
        s[0].vel = [ 0, 0, 0 ]
        s[0].mass = 1
        fill(s.attributes, 0)
        for j in range(0,d.nbod):
          fill(s[j].attributes, 0)

        for j in range(1,d.nbod):
            r = spacing_factor ** (j-1)
            v = sqrt(1/r) * ejection_factor
            phi = random.uniform(0,2*pi)
            s[j].pos = [  r*cos(phi), r*sin(phi), 0 ]
            s[j].vel = [ -v*sin(phi), v*cos(phi), 0 ]
            s[j].mass = planet_mass 
    return d;


def norm(l):
	return sqrt(sum(x**2 for x in l))


class BasicIntegration(abstract.IntegrationTest):
    cfg = swarmng.config(
            integrator = 'hermite_cpu',
            time_step  = 1e-4,
            nogpu      = 1
            )
    def createEnsemble(self):
        return make_test_case(nsys = 16, nbod = 3, spacing_factor=1.4)
    def examine(self):
        max_deltaE = swarmng.find_max_energy_conservation_error( self.ens, self.ref)
        self.assertLess(max_deltaE, 1e-13)


class InitialConditions(unittest.TestCase):
    def runTest(self):
        ref = make_test_case(nsys = 16, nbod = 3, spacing_factor=1.4, seed = 14321)
        ref.save_to_text("sample.txt")
        r = system("diff sample.txt '{0}'".format(path.join(TESTDIR,"ref.txt")) )
        self.assertEqual(r, 0)


class EjectionTest(abstract.IntegrationTest):
    RMAX = 10
    cfg = swarmng.config(
            integrator = "hermite_cpu",
            time_step  = 1e-3,
            nogpu      = 1,
            deactivate_on_ejection = 1,
            rmax       = RMAX
            )
    destination_time = 100

    def createEnsemble(self):
        return make_test_case(nsys=20, nbod = 6, spacing_factor=1.01);

    def examine(self):
        for sys in self.ens:
            for b in sys:
                if( b.distance_to_origin() > self.RMAX ):
                    self.assertEqual(sys.state, -1, "a system with an ejected body should be disabled")





