# Use ${SAMPLES} if you need to reference any sample input or configuration file

ENABLE_TESTING()

ADD_TEST(NAME basic COMMAND swarm integrate --defaults )

MACRO(TEST_INTEGRATOR title cfgname)
	ADD_TEST(NAME Verify_${title}
		COMMAND swarm test -c ${SAMPLES}/${cfgname} -i ${SAMPLES}/test.in -o ${SAMPLES}/test.out)
ENDMACRO(TEST_INTEGRATOR title cfgname)


TEST_INTEGRATOR(Hermite_CPU hermite_cpu.cfg)
TEST_INTEGRATOR(Hermite_GPU hermite.cfg)
TEST_INTEGRATOR(Runge_Kutta_Fixed_Time_Step rkckf.cfg)
TEST_INTEGRATOR(Runge_Kutta_Adaptive_Time_Step rkcka.cfg)
TEST_INTEGRATOR(Euler euler.cfg)
