#!/bin/sh


FILE=$1
TR=$2

SRCDIR=`dirname $0`/..
SWARM=bin/swarm

echo Generating report from BDB log
time $SWARM query -f "$FILE" -t $TR > "$FILE.c.txt"
echo Generating report from BDB log using Python
time python $SRCDIR/py/query/swarm-query.py -d "$FILE" -t $TR -b 1..2 > "$FILE.py.txt"

if diff "$FILE.c.txt" "$FILE.py.txt"
then
    echo "Passed"
    return 0
else
    echo "Fail"
    return 1
fi

