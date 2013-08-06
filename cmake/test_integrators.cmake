# Test all integrators as they are listed here
#   1. Test for every number of bodies
#   2. Test for stability of the integrators
#
#  Note:
#  1. There are sample input/output files only up to 9 bodies/system
#  2. Configurations for running each integrator should be provided

# Directory where the tests are located
SET(TEST_INTEGRATOR_DIR ${TESTDIR}/integrators)

SET(TEST_integrator_nbod TRUE CACHE BOOL "Test all integrators versus all number of bodies")
SET(TEST_stability TRUE CACHE BOOL "Test stability of all integrators")
SET(TEST_stability_nsys 16 CACHE STRING "Number of systems for stability test")
SET(TEST_stability_duration 100 CACHE STRING "Duration of the stability test")
#SET(TEST_integrator_nbod_list CACHE LIST "List of integrators to be tested vs. nbod")
#SET(TEST_stability_list CACHE LIST "List of integrators (name of configuration file) that are tested for stability (run for long time)")

## Automatically test all plugins that have a cfg file.

FOREACH(id ${SWARM_PLUGINS})
	IF((${PLUGIN_${id}}) AND (EXISTS "${TEST_INTEGRATOR_DIR}/${id}.cfg"))
		LIST(APPEND TEST_integrator_nbod_list ${id})
		LIST(APPEND TEST_stability_list ${id})
	ENDIF()
ENDFOREACH()

MESSAGE("Will Test: ${TEST_integrator_nbod_list}")


##################### MACROS START HERE #########################################

MACRO(TEST_INTEGRATOR title nbod)
	ADD_TEST(NAME Verify_${title}_${nbod}_Bodies
		COMMAND swarm test -c ${TEST_INTEGRATOR_DIR}/${title}.cfg -v 100
        -I ${TEST_INTEGRATOR_DIR}/test.${nbod}.in.txt 
        -O ${TEST_INTEGRATOR_DIR}/test.${nbod}.out.txt 
        nbod=${nbod} destination_time=1 pos_threshold=1E-8 vel_threshold=1E-8) 
ENDMACRO(TEST_INTEGRATOR)


MACRO(TEST_INTEGRATOR_STABILITY title nbod nsys)
	ADD_TEST(NAME Test_Stability_${title}_${nsys}_Systems_${nbod}_Bodies
		COMMAND swarm integrate -v 100 -c ${TEST_INTEGRATOR_DIR}/${title}.cfg -d ${TEST_stability_duration} -n 10 nbod=${nbod} nsys=${nsys} ) 
ENDMACRO(TEST_INTEGRATOR_STABILITY)

##### Test Integrator X nbod for pre-calculated scenarios #######################################
IF(TEST_integrator_nbod)
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
ENDIF()

#### Test Integrator X nbod for stability ######################################################
IF(TEST_stability)
	LIST(LENGTH TEST_stability_list integ_list_length)
	SET(i 0)
	WHILE(i LESS ${integ_list_length})
		LIST(GET TEST_stability_list ${i} item)
		SET(nbod 3)
		WHILE(NOT nbod GREATER ${MAX_NBODIES})
			TEST_INTEGRATOR_STABILITY(${item} ${nbod} ${TEST_stability_nsys})
			MATH( EXPR nbod "${nbod} + 1")
		ENDWHILE(NOT nbod GREATER ${MAX_NBODIES})
		MATH( EXPR i "${i} + 1" )
	ENDWHILE(i LESS ${integ_list_length})
ENDIF()
