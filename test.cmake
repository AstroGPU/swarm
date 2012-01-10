# Use ${SAMPLES} if you need to reference any sample input or configuration file

ENABLE_TESTING()
MACRO(TEST_INTEGRATOR title nbod)
	ADD_TEST(NAME Verify_${title}_${nbod}_Bodies
		COMMAND swarm test -c ${SAMPLES}/${title}.cfg -I ${SAMPLES}/test.${nbod}.in.txt -O ${SAMPLES}/test.${nbod}.out.txt -v 100 nbod=${nbod}) 
ENDMACRO(TEST_INTEGRATOR)


MACRO(TEST_INTEGRATOR_STABILITY title nbod nsys)
	ADD_TEST(NAME Test_Stability_${title}_${nsys}_Systems_${nbod}_Bodies
		COMMAND swarm integrate -v 100 -c ${SAMPLES}/${title}.cfg -d 100 -n 10 nbod=${nbod} nsys=${nsys} ) 
ENDMACRO(TEST_INTEGRATOR_STABILITY)


############ ACTUAL TEST Cases Begin Here

ADD_TEST(NAME basic COMMAND swarm integrate --defaults )
TEST_INTEGRATOR(CPU 3)

##### Test Integrator X nbod for pre-calculated scenarios
SET(integrators Hermite Hermite_Adaptive Runge_Kutta_Fixed_Time_Step Runge_Kutta_Adaptive_Time_Step Euler Midpoint Verlet)
LIST(LENGTH integrators integ_list_length)
SET(i 0)
WHILE(i LESS ${integ_list_length})
	LIST(GET integrators ${i} item)
	SET(nbod 3)
	WHILE(NOT nbod GREATER ${MAX_NBODIES})
		TEST_INTEGRATOR(${item} ${nbod})
		MATH( EXPR nbod "${nbod} + 1")
	ENDWHILE(NOT nbod GREATER ${MAX_NBODIES})
	MATH( EXPR i "${i} + 1" )
ENDWHILE(i LESS ${integ_list_length})

#### Test Integrator X nbod for stability
SET(integrators Hermite Hermite_Adaptive Runge_Kutta_Fixed_Time_Step Runge_Kutta_Adaptive_Time_Step )
LIST(LENGTH integrators integ_list_length)
SET(i 0)
WHILE(i LESS ${integ_list_length})
	LIST(GET integrators ${i} item)
	SET(nbod 3)
	WHILE(NOT nbod GREATER ${MAX_NBODIES})
		TEST_INTEGRATOR_STABILITY(${item} ${nbod} 1000)
		MATH( EXPR nbod "${nbod} + 1")
	ENDWHILE(NOT nbod GREATER ${MAX_NBODIES})
	MATH( EXPR i "${i} + 1" )
ENDWHILE(i LESS ${integ_list_length})



## Add your test cases after this line (or include a file containing test cases)

## Including a macro file (should be in cmake format)
# INCLUDE(mytestcases.cmake)

## Add a test case
# ADD_TEST(NAME blah COMMAND swarm .... )

# For each sample
#SET(l 40 20 30 10)
#LIST(LENGTH l list_length)
#SET(i 0)
#WHILE(i LESS ${list_length})
#	LIST(GET l ${i} item)
#	MESSAGE(${item})
#	MATH( EXPR i "${i} + 1" )
#ENDWHILE(i LESS ${list_length})

