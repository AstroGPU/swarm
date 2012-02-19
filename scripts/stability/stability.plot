#set size 4, 5
set terminal postscript eps enhanced color dashed lw 1 "Helvetica" 14
set output 'plot.eps'
set key left top
set title 'Averaged Plot for MVS'
plot "datafile.averaged" using 1:2 title "3 Body" with lines, \
     "datafile.averaged" using 1:3 title "4 Body" with lines, \
     "datafile.averaged" using 1:4 title "5 Body" with lines, \
     "datafile.averaged" using 1:5 title "6 Body" with lines
