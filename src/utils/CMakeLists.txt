SWARM_ADD_EXECUTABLE(stability_test stability_test.cpp )
SWARM_ADD_EXECUTABLE(swarm_benchmark swarm_benchmark.cpp )
SWARM_ADD_EXECUTABLE(swarm_benchmark_gpu swarm_benchmark_gpu.cpp )

SWARM_ADD_EXECUTABLE(verify_integrator verify_integrator.cpp )

SWARM_ADD_EXECUTABLE(test_energy test_energy.cpp)
SWARM_ADD_EXECUTABLE(swarmquery swarmquery.cpp)
SWARM_ADD_EXECUTABLE(swarmquerykeplerian swarmquerykeplerian.cpp)
SWARM_ADD_EXECUTABLE(snapshot2text snapshot2text.cpp 
	 ../swarm/snapshot.cpp)
SWARM_ADD_EXECUTABLE(text2snapshot text2snapshot.cpp 
	 ../swarm/snapshot.cpp)

