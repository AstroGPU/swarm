###
### Applications
###
APPS+=swarm 
swarm_SOURCES=src/utils/swarm.cpp src/utils/utils.cpp

####APPS+=swarmquery
####swarmquery_SOURCES=src/swarmquery.cpp

####APPS+=swarmquerykeplerian
####swarmquerykeplerian_SOURCES=src/swarmquerykeplerian.cpp

APPS+=test_energy
test_energy_SOURCES=src/utils/test_energy.cpp

####APPS+=swarm_benchmark
####swarm_benchmark_SOURCES=src/tutorials/swarm_benchmark.cpp src/tutorials/utils.cpp

APPS+=stability_test
stability_test_SOURCES=src/utils/stability_test.cpp  src/utils/utils.cpp
