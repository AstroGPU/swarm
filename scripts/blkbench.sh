out=block-benchmark.txt
rm -f $out
for  i in 64 128 256 512
do
	../bin/swarm_tutorial_benchmark --nocpu -t 1. -n 6 -s 3840 -b $i 2>> $out
done
grep "(integration)" $out  > $out.times
