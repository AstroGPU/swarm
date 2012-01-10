# Use ${SAMPLES} if you need to reference any sample input or configuration file

ENABLE_TESTING()
MACRO(TEST_INTEGRATOR title nbod)
	ADD_TEST(NAME Verify_${title}_${nbod}_Bodies
		COMMAND swarm test -c ${SAMPLES}/${title}.cfg -I ${SAMPLES}/test.${nbod}.in.txt -O ${SAMPLES}/test.${nbod}.out.txt -v 100 nbod=${nbod}) 
ENDMACRO(TEST_INTEGRATOR)



############ ACTUAL TEST Cases Begin Here

ADD_TEST(NAME basic COMMAND swarm integrate --defaults )


TEST_INTEGRATOR(CPU 3)

SET(nbod 3)
WHILE(NOT nbod GREATER ${MAX_NBODIES})
	TEST_INTEGRATOR(Hermite ${nbod})
	TEST_INTEGRATOR(Hermite_Adaptive ${nbod})
	TEST_INTEGRATOR(Runge_Kutta_Fixed_Time_Step ${nbod})
	TEST_INTEGRATOR(Runge_Kutta_Adaptive_Time_Step ${nbod})
	TEST_INTEGRATOR(Euler ${nbod})
	TEST_INTEGRATOR(Midpoint ${nbod})
	TEST_INTEGRATOR(Verlet ${nbod})
	MATH( EXPR nbod "${nbod} + 1")
ENDWHILE(NOT nbod GREATER ${MAX_NBODIES})


## Add your test cases after this line (or include a file containing test cases)

## Including a macro file (should be in cmake format)
# INCLUDE(mytestcases.cmake)

## Add a test case
# ADD_TEST(NAME blah COMMAND swarm .... )
