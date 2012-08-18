# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	add_custom_target(doc
		COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
		)
else(DOXYGEN_FOUND)
	add_custom_target(doc COMMENT "Doxygen not found, cannot make
	documentation")
endif(DOXYGEN_FOUND)


MACRO(tutorial_python pageName )
	ADD_CUSTOM_TARGET(document_${pageName}
	COMMAND mkdir -p ${CMAKE_BINARY_DIR}/docs
	COMMAND cat ${ARGN} | ${CMAKE_SOURCE_DIR}/scripts/python2doxygen  > ${CMAKE_BINARY_DIR}/docs/${pageName}.dox
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	)
	ADD_DEPENDENCIES(doc document_${pageName})
ENDMACRO()

MACRO(tutorial_file pageName )
	ADD_CUSTOM_TARGET(document_${pageName}
	COMMAND mkdir -p ${CMAKE_BINARY_DIR}/docs
	COMMAND cat ${ARGN} | ${CMAKE_SOURCE_DIR}/scripts/program2doxygen  > ${CMAKE_BINARY_DIR}/docs/${pageName}.dox
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	)
	ADD_DEPENDENCIES(doc document_${pageName})
ENDMACRO()
