#!/bin/bash

# Thorough testing of swarm, it consists of running make test
# for different compile-time configuration parameters.
# We still don't know what to test. Maybe the fast math approximations
# can be included in this.
CHUNKSIZE_LIST="1 2 4 8 16"
OUTPUT=Thorough-Testing

mkdir -p $OUTPUT

for C in $CHUNKSIZE_LIST 
do
	cmake . -DSHMEM_CHUNK_SIZE=$C -DENSEMBLE_CHUNK_SIZE=$C
	
	make swarm
	make test
	cp Testing/Temporary/LastTest.log $OUTPUT/Log-$C.log
	cp Testing/Temporary/LastTestsFailed.log $OUTPUT/Log-Failed-$C.log
done
