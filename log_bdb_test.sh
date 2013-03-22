#bin/swarm integrate nsys=100 nbod=3 log_writer=bdb log_output_db=testing_bdb log_interval=0.1 -d 1 integrator=hermite_cpu_log time_step=0.01

cat > testing_log.cfg <<EOF
nsys=100
nbod=3
log_output_db=testing_log.db
log_output=testing_log.bin
log_interval=0.1
destination_time=1
integrator=hermite_cpu_log
time_step=0.01
EOF

SWARM=bin/swarm

$SWARM integrate --nogpu -c testing_log.cfg  log_writer=bdb
$SWARM integrate --nogpu -c testing_log.cfg  log_writer=binary
$SWARM query -f testing_log.db > testing_log.db.txt
$SWARM query -f testing_log.bin.raw > testing_log.bin.txt


if diff testing_log.db.txt testing_log.bin.txt
then
    echo "Passed"
else
    echo "Fail"
fi
