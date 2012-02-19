
AVERAGE_COUNT=10
import os,sys
import fileinput
import Gnuplot

i = 0
for line in fileinput.input():
	columns = map( float , filter( lambda x: len(x) > 0 , line.split(',')))
#	print "see this: ",  reduce(lambda x,y: str(x) + ', ' + str(y), columns)
	if( i == 0 ):
		total = columns
		i = 1
	elif ( i == AVERAGE_COUNT-1):
		i = 0
		average = map( lambda x: x / AVERAGE_COUNT, total )
		print reduce(lambda x,y: str(x) + ', ' + str(y), average)
	else:
		i = i + 1
		total = map( lambda x,y: x + y , columns, total )


