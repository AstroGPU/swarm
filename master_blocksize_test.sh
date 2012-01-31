#!/bin/bash

for CS in 1 2 4 8 16 32
do
	cmake . -DSHMEM_CHUNK_SIZE=$CS -DENSEMBLE_CHUNK_SIZE=$CS
	grep CHUNK_SIZE CMakeCache.txt
	make clean
	make swarm
	bin/swarm -c blocktest.cfg benchmark --range system_per_block=${CS}..${CS}..64
done
