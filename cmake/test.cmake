# Use ${SAMPLES} if you need to reference any sample input or configuration file

ENABLE_TESTING()

MACRO(TEST_SCENARIO title config)
ADD_TEST(NAME ${title}_${config} 
	COMMAND swarm test -c ${TESTDIR}/${title}/${config}.cfg -I ${TESTDIR}/${title}/in.txt -O ${TESTDIR}/${title}/out.txt)
ENDMACRO(TEST_SCENARIO)

############ ACTUAL TEST Cases Begin Here

ADD_TEST(NAME Basic_integration_on_CPU COMMAND swarm integrate --defaults --nogpu integrator=hermite_cpu)
ADD_TEST(NAME Basic_integration_on_GPU COMMAND swarm integrate --defaults )

ADD_TEST(NAME "BDB"
    COMMAND "${CMAKE_SOURCE_DIR}/test/bdb/bdb.sh" )


INCLUDE(cmake/test_integrators.cmake)


INCLUDE(cmake/test_monitors.cmake)

ADD_TEST(NAME "Python_Tests"
	COMMAND "${CMAKE_SOURCE_DIR}/py/tests/run.py")

# TEST_SCENARIO makes it easy to create scenarios and add it to the system
# The first argument in the name of the folder and the second argument is the name of the 
# config file (without extension). Just put in.txt and out.txt in that folder alongside with
# configuration files and TEST_SCENARIO integrates the configuration and verifies it against 
# the test case (in.txt/out.txt).
TEST_SCENARIO(sample CPU)

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

