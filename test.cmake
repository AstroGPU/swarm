# Use ${SAMPLES} if you need to reference any sample input or configuration file

ENABLE_TESTING()

MACRO(TEST_SCENARIO title config)
ADD_TEST(NAME ${title}_${config} 
	COMMAND swarm test -c ${TESTDIR}/${title}/${config}.cfg -I ${TESTDIR}/${title}/in.txt -O ${TESTDIR}/${title}/out.txt)
ENDMACRO(TEST_SCENARIO)

############ ACTUAL TEST Cases Begin Here

ADD_TEST(NAME basic COMMAND swarm integrate --defaults )

INCLUDE(test_integrators.cmake)


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

