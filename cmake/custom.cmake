# Use ${SAMPLES} if you need to reference any sample input or configuration file





SET(BENCHMARK_DIR ${TESTDIR}/benchmark)
MACRO(BENCHMARK_INTEGRATOR cfg)

	ADD_CUSTOM_TARGET(benchmark_${cfg})
	ADD_DEPENDENCIES(benchmark benchmark_${cfg})

	MACRO(BENCHMARK_INTEGRATOR_NBOD nbod)
		ADD_CUSTOM_TARGET(benchmark_${cfg}_${nbod}
			COMMAND swarm benchmark --range nsys=1000,2000,4000,8000,16000 nbod=${nbod} -c ${BENCHMARK_DIR}/${cfg}.cfg)
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

SET(cfg Hermite)
SET(nbod 3)

ADD_CUSTOM_TARGET(benchmark_block_${cfg}
	COMMAND swarm benchmark --range systems_ber_block=1,2,4,8,16,32 nsys=16000 nbod=${nbod} -c ${BENCHMARK_DIR}/${cfg}.cfg)
ADD_DEPENDENCIES(benchmark_block_${cfg} swarm)
ADD_DEPENDENCIES(benchmark benchmark_block_${cfg})	

ADD_CUSTOM_TARGET(pytest COMMAND "${CMAKE_SOURCE_DIR}/py/runtests.py")

