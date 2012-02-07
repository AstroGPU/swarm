#!/bin/bash

for CS in 4 8 16 32
do
	cmake . -DSHMEM_CHUNK_SIZE=$CS -DENSEMBLE_CHUNK_SIZE=$CS > /dev/null
	grep CHUNK_SIZE CMakeCache.txt
	make clean > /dev/null
	make swarm > /dev/null
	for n in 3 4 5 6 7 8 9
	do
		bin/swarm -c blocktest.cfg benchmark --range system_per_block=${CS}..${CS}..64 nbod=$n
	done
	echo "----------------------------------------------------------------------------------"
done
