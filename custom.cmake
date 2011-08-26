# Use ${SAMPLES} if you need to reference any sample input or configuration file

ADD_CUSTOM_TARGET(benchmark
	COMMAND swarm benchmark nsys 4000 8000 16000 32000 64000 -d 10)
ADD_DEPENDENCIES(benchmark swarm)
