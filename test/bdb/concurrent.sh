#!/bin/bash

# Testing that concurrent reads from a file works
# 
# swarm command is run to generate a log file and it is scheduled to remove it
# immedaitely after the integration is done.
#
# this test only passes if the query is quick enough to open the file and get the result.
#
#
# 
TESTDIR=`dirname $0`


cat > testing_log.cfg <<EOF
log_output_db=testing_log.db
log_output=testing_log.bin
log_interval=1
destination_time=10
integrator=hermite_cpu_log
time_step=0.01
nogpu=1
EOF

SWARM=bin/swarm

rm -f testing_log.db

$SWARM integrate -n 1 -I $TESTDIR/test.4.in.txt  -c testing_log.cfg  log_writer=bdb && rm testing_log.db testing_log.cfg  &
sleep .1
stat testing_log.db
$SWARM query -f testing_log.db 

wait
