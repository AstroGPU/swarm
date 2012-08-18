#!/bin/bash
SHORTDATE='date +%Y-%m-%d'
OUT=stability_`$SHORTDATE`
NBODY_RANGE="3 4 5 6"
base=`dirname $0`
mkdir -p $OUT

{
for INT in Hermite
#_Adaptive Runge_Kutta_Adaptive_Time_Step MVS
do
	echo "Starting stability test for $INT"
	$base/stability.sh $OUT $INT "$NBODY_RANGE"
	echo $INT | tr _ " " > $OUT/$INT/title.txt
	echo "$INT is done"
done
} | tee $OUT/log.txt
for n in $NBODY_RANGE
do
$base/make_chart.sh $n $OUT
done
