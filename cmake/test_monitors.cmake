
SET(TEST_monitors_nbod TRUE CACHE BOOL "Test the monitors against different number of bodies")

SET(TEST_MONITOR_DIR ${TESTDIR}/monitors)

ADD_CUSTOM_TARGET(testcases)

MACRO(TEST_MONITOR title subject nbod)
	SET(DD ${TEST_MONITOR_DIR}/${title})
	ADD_TEST(NAME Test_${title}_on_${subject}_with_${nbod}_Bodies
		COMMAND swarm test -c ${DD}/${subject}.cfg -I ${DD}/${nbod}.in.txt -O ${DD}/${subject}.${nbod}.out.txt -v 100 nbod=${nbod}) 

	# Create a make target for generating the test case data
	ADD_CUSTOM_TARGET(Generate_Testcase_${title}_on_${subject}_with_${nbod}_Bodies
		COMMAND swarm integrate -c ${DD}/${subject}.cfg -I ${DD}/${nbod}.in.txt -O ${DD}/${subject}.${nbod}.out.txt -v 100 nbod=${nbod}) 
	ADD_DEPENDENCIES(Generate_Testcase_${title}_on_${subject}_with_${nbod}_Bodies Generate_Testcase_${title}_with_${nbod}_Bodies)
	ADD_DEPENDENCIES(testcases Generate_Testcase_${title}_on_${subject}_with_${nbod}_Bodies)


ENDMACRO()

MACRO(ADD_MONITOR_TESTCASE title nbod)
	SET(DD ${TEST_MONITOR_DIR}/${title})
	ADD_CUSTOM_TARGET(Generate_Testcase_${title}_with_${nbod}_Bodies
		COMMAND ${DD}/generate ${nbod} ${DD}/${nbod}.in.txt)
ENDMACRO()

IF(TEST_monitors_nbod)
#	ADD_MONITOR_TESTCASE(close_orbits 4)
#	ADD_MONITOR_TESTCASE(collision_course 4)
#	TEST_MONITOR(close_orbits HermiteEjection 4)
#	TEST_MONITOR(close_orbits MVSEjection 4)
#	TEST_MONITOR(close_orbits HermiteCloseEncounter 4)
#	TEST_MONITOR(collision_course HermiteCloseEncounter 4)
#	ADD_MONITOR_TESTCASE(straight_line 3)
#	TEST_MONITOR(straight_line HermiteCloseEncounter 3)

	ADD_MONITOR_TESTCASE(hyperbolic_orbits 3)
	TEST_MONITOR(hyperbolic_orbits HermiteEjection 3)
	
	ADD_MONITOR_TESTCASE(parabolic_collision 3)
	TEST_MONITOR(parabolic_collision HermiteCloseEncounter 3)
ENDIF()

