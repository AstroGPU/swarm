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
	PREFIX="$OUTPUT/$C"
	
	{ cmake . -DSHMEM_CHUNK_SIZE=$C -DENSEMBLE_CHUNK_SIZE=$C ; make swarm; } 2>&1 | tee $PREFIX-Build.log;
	
	make test  2>&1 | tee $PREFIX-Test.log
	cp Testing/Temporary/LastTest.log $PREFIX-Detail.log
	cp Testing/Temporary/LastTestsFailed.log $PREFIX-Failed.log
done
