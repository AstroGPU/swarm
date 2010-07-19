out=bodies-benchmark.txt
rm -f $out
for ((i=3;$i<10;i++))
do
	../bin/swarm_tutorial_benchmark --nocpu -t 1. -n $i -s 7680 2>> $out
done
grep "(integration)" $out  > $out.times
