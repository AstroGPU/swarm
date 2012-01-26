# Use ${SAMPLES} if you need to reference any sample input or configuration file


ADD_CUSTOM_TARGET(document_tutorials
	COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/program2doxygen < ${CMAKE_CURRENT_SOURCE_DIR}/src/tutorials/swarm_tutorial_simple.cpp > ${CMAKE_CURRENT_BINARY_DIR}/docs/swarm_tutorial_simple.dox
	COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/program2doxygen < ${CMAKE_CURRENT_SOURCE_DIR}/src/tutorials/tutorial_gpu.cpp > ${CMAKE_CURRENT_BINARY_DIR}/docs/tutorial_gpu.dox)
add_dependencies(doc document_tutorials)


SET(BENCHMARK_DIR ${TESTDIR}/benchmark)
MACRO(BENCHMARK_INTEGRATOR cfg)

	ADD_CUSTOM_TARGET(benchmark_${cfg})
	ADD_DEPENDENCIES(benchmark benchmark_${cfg})

	MACRO(BENCHMARK_INTEGRATOR_NBOD nbod)
		ADD_CUSTOM_TARGET(benchmark_${cfg}_${nbod}
			COMMAND swarm benchmark --range nsys=4000,8000,16000,32000 nbod=${nbod} -c ${BENCHMARK_DIR}/${cfg}.cfg)
		ADD_DEPENDENCIES(benchmark_${cfg}_${nbod} swarm)
		ADD_DEPENDENCIES(benchmark_${cfg} benchmark_${cfg}_${nbod})	
	ENDMACRO()

	BENCHMARK_INTEGRATOR_NBOD(3)
	BENCHMARK_INTEGRATOR_NBOD(4)
	BENCHMARK_INTEGRATOR_NBOD(5)

ENDMACRO()

ADD_CUSTOM_TARGET(benchmark)

BENCHMARK_INTEGRATOR(Hermite)
BENCHMARK_INTEGRATOR(Hermite_Adaptive)
