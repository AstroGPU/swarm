#bin/swarm integrate nsys=100 nbod=3 log_writer=bdb log_output_db=testing_bdb log_interval=0.1 -d 1 integrator=hermite_cpu_log time_step=0.01

cat > testing_log.cfg <<EOF
log_output_db=testing_log.db
log_output=testing_log.bin
log_interval=0.1
destination_time=1
integrator=hermite_adap_log
time_step=0.01
time_step_factor=170e-4
min_time_step=0.00001
max_time_step=0.01
EOF

SWARM=bin/swarm

rm -f testing_log.db testing_log.bin.raw

$SWARM generate -o test.ens nsys=100 nbod=3
$SWARM integrate -i test.ens -c testing_log.cfg  log_writer=bdb
$SWARM integrate -i test.ens -c testing_log.cfg  log_writer=binary
$SWARM query -f testing_log.db > testing_log.db.txt
$SWARM query -f testing_log.bin.raw > testing_log.bin.txt


diff testing_log.db.txt testing_log.bin.txt
