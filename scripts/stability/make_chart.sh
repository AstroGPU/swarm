if [ $# -lt 2 ]
then
	echo "Usage: $0 'filespattern' outputname"
	return 1
fi

pattern=$1
output=$2
for f in $pattern
do
	tail -n +8 $f > $f.isolated
	#| sed -n 's/^[^,]*,\([^,]*\),.*$/\1/p' 
done


nawk -v c=2 -f io.awk $pattern.isolated | sort -n > $output.csv
python average.py $output.csv > datafile.averaged
gnuplot stability.plot
rm datafile.averaged
mv plot.eps $output.eps
ps2pdf $output.eps
