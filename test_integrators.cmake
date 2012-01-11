# Test all integrators as they are listed here
#   1. Test for every number of bodies
#   2. Test for stability of the integrators
#
#  Note:
#  1. There are sample input/output files only up to 9 bodies/system
#  2. Configurations for running each integrator should be provided

# List of integrators (name of configuration file) that are tested versus all number of bodies
SET(TEST_integrator_nbod_list CPU Hermite Hermite_Adaptive Runge_Kutta_Fixed_Time_Step Runge_Kutta_Adaptive_Time_Step CACHE LIST "List of integrators to be tested vs. nbod")

# List of integrators (name of configuration file) that are tested for stability (run for long time)
SET(TEST_integrator_TEST_stability_list Hermite Hermite_Adaptive Runge_Kutta_Fixed_Time_Step Runge_Kutta_Adaptive_Time_Step CACHE LIST "List of integrators to be tested for stability")

# Number of systems generated for the stability test
SET(TEST_stability_nsys 16 CACHE STRING "Number of systems for stability test")
# duration of the stability test
SET(TEST_stability_duration 100 CACHE STRING "Duration of the stability test")

SET(TEST_INTEGRATOR_DIR ${TESTDIR}/integrators)

##################### MACROS START HERE #########################################

MACRO(TEST_INTEGRATOR title nbod)
	ADD_TEST(NAME Verify_${title}_${nbod}_Bodies
		COMMAND swarm test -c ${TEST_INTEGRATOR_DIR}/${title}.cfg -I ${TEST_INTEGRATOR_DIR}/test.${nbod}.in.txt -O ${TEST_INTEGRATOR_DIR}/test.${nbod}.out.txt -v 100 nbod=${nbod}) 
ENDMACRO(TEST_INTEGRATOR)


MACRO(TEST_INTEGRATOR_STABILITY title nbod nsys)
	ADD_TEST(NAME Test_Stability_${title}_${nsys}_Systems_${nbod}_Bodies
		COMMAND swarm integrate -v 100 -c ${TEST_INTEGRATOR_DIR}/${title}.cfg -d ${TEST_stability_duration} -n 10 nbod=${nbod} nsys=${nsys} ) 
ENDMACRO(TEST_INTEGRATOR_STABILITY)

##### Test Integrator X nbod for pre-calculated scenarios #######################################
LIST(LENGTH TEST_integrator_nbod_list integ_list_length)
SET(i 0)
WHILE(i LESS ${integ_list_length})
	LIST(GET TEST_integrator_nbod_list ${i} item)
	SET(nbod 3)
	WHILE(NOT nbod GREATER ${MAX_NBODIES})
		TEST_INTEGRATOR(${item} ${nbod})
		MATH( EXPR nbod "${nbod} + 1")
	ENDWHILE(NOT nbod GREATER ${MAX_NBODIES})
	MATH( EXPR i "${i} + 1" )
ENDWHILE(i LESS ${integ_list_length})

#### Test Integrator X nbod for stability ######################################################
LIST(LENGTH TEST_integrator_TEST_stability_list integ_list_length)
SET(i 0)
WHILE(i LESS ${integ_list_length})
	LIST(GET TEST_integrator_TEST_stability_list ${i} item)
	SET(nbod 3)
	WHILE(NOT nbod GREATER ${MAX_NBODIES})
		TEST_INTEGRATOR_STABILITY(${item} ${nbod} ${TEST_stability_nsys})
		MATH( EXPR nbod "${nbod} + 1")
	ENDWHILE(NOT nbod GREATER ${MAX_NBODIES})
	MATH( EXPR i "${i} + 1" )
ENDWHILE(i LESS ${integ_list_length})

