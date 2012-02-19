if [ $# -lt 2 ]
then
	echo    "Usage: $0 <nbodies> <Directory>"
	echo -e "\t nbodies: Number of bodies to make chart for"
	echo -e "\t Directory: Directory containing results"
	exit
fi

#Number of bodies
nb=$1
#Input directory
dir=$2

base=`dirname $0`

for int in $dir/*
do
	f=$int/log_$nb
	[ -f "$f.csv" ] && tail -n +8 $f.csv | python $base/average.py > $f.avg.csv
done

PLOT=$dir/$nb-body.plot
PLOT_OUTPUT=$dir/plot-$nb.eps
PLOT_OUTPUT_PDF=$dir/plot-$nb.pdf

cat > $PLOT <<EOF
#set size 4, 5
set terminal postscript eps enhanced color dashed lw 1 "Helvetica" 14
set output '$PLOT_OUTPUT'
set key left top
set title 'Stability Plot for $nb bodies'
EOF

{
	echo -n "plot "
	M=
	for int in $dir/*
	do
		f=$int/log_$nb.avg.csv
		if [ -f "$f" ]
		then
			[ -n "$M" ] && echo -n ", "
			M=1
			echo -n "\"$f\" using 1:3 title \"`cat $int/title.txt`\" with lines"
		fi
	done
	echo
} >> $PLOT
gnuplot $PLOT
ps2pdf $PLOT_OUTPUT $PLOT_OUTPUT_PDF
