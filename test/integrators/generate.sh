#!/bin/bash

DIR=`dirname $0`

for ((i=3;$i<=13;i++))
do
  bin/swarm generate -c $DIR/generate.cfg nbod=$i -O $DIR/test.$i.in.txt
  bin/swarm integrate -c $DIR/Hermite_CPU.cfg -I $DIR/test.$i.in.txt -O $DIR/test.$i.out.txt
done