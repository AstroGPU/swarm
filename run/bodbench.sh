out=bodies-benchmark.txt
nsys=3840
PRE=$1
rm -f $out
export CUDA_PROFILE=1
export CUDA_PROFILE_CONFIG="profiler.cfg" 
for ((i=3;$i<7;i++))
do
	export CUDA_PROFILE_LOG="benchmark.$i.clog"
	../bin/swarm_tutorial_benchmark --cfg ${PRE}integrator.cfg --nocpu -t 1. -n $i -s $nsys -b  64 2>> $out
done
grep "(integration)" $out  > ${PRE}$out.times
grep swarm benchmark.*.clog > ${PRE}benchmark.clog.brief
