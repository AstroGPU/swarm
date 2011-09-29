# Use ${SAMPLES} if you need to reference any sample input or configuration file

ADD_CUSTOM_TARGET(benchmark
	COMMAND swarm benchmark nsys 4000 8000 16000 32000 64000 -d 10)
ADD_DEPENDENCIES(benchmark swarm)

ADD_CUSTOM_TARGET(document_tutorials
	COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/program2doxygen < ${CMAKE_CURRENT_SOURCE_DIR}/src/tutorials/swarm_tutorial_simple.cpp > ${CMAKE_CURRENT_BINARY_DIR}/docs/swarm_tutorial_simple.dox
	COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/program2doxygen < ${CMAKE_CURRENT_SOURCE_DIR}/src/tutorials/tutorial_gpu.cpp > ${CMAKE_CURRENT_BINARY_DIR}/docs/tutorial_gpu.dox)

