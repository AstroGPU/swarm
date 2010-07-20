#!/bin/sh
export CUDA_PROFILE=1
export CUDA_PROFILE_CONFIG="profiler.cfg"
export CUDA_PROFILE_LOG="cudalog.txt"
../scripts/easyGen.py
../bin/swarm integrator.cfg
../bin/swarm_test_energy log.bin > e1.txt
../bin/swarmquery log.bin  > s1.txt
tail -n3048 s1.txt > sd1.txt

export CUDA_PROFILE_LOG="cudalog-ref.txt"
../bin/swarm integrator-ref.cfg
../bin/swarm_test_energy log.bin > e2.txt
../bin/swarmquery log.bin > s2.txt
tail -n3048 s2.txt > sd2.txt

rm data.*
