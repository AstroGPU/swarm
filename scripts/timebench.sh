nsys=3840
nbod=3
PRE=$1
out=${PRE}times-benchmark.txt
rm -f $out
export CUDA_PROFILE=1
export CUDA_PROFILE_CONFIG="profiler.cfg" 
for ((i=1;$i<65;i++))
#for i in 1 2 4 8 15 16 32 64 
do
	export CUDA_PROFILE_LOG="benchmark.$i.clog"
	../bin/swarm_tutorial_benchmark --cfg ${PRE}integrator.cfg --nocpu -t $i -n $nbod -s $nsys -p 2 -b  64 2>> $out
done
grep "Max" $out  > ${PRE}$out.times
grep swarm benchmark.*.clog > ${PRE}benchmark.clog.brief
rm benchmark.*.clog
