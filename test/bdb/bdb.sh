#bin/swarm integrate nsys=100 nbod=3 log_writer=bdb log_output_db=testing_bdb log_interval=0.1 -d 1 integrator=hermite_cpu_log time_step=0.01

TESTDIR=`dirname $0`

OUTPUTDIR=Testing

cat > testing_log.cfg <<EOF
log_output_db=$OUTPUTDIR/testing_log.db
log_output=$OUTPUTDIR/testing_log.bin
log_interval=0.1
destination_time=1
integrator=hermite_cpu_log
time_step=0.01
time_step_factor=170e-4
min_time_step=0.00001
max_time_step=0.01
nogpu=1
EOF

SWARM=bin/swarm

rm -f $OUTPUTDIR/testing_log.db $OUTPUTDIR/testing_log.bin.raw $OUTPUTDIR/__db.*

$SWARM integrate -I $TESTDIR/test.4.in.txt  -c testing_log.cfg  log_writer=bdb
$SWARM integrate -I $TESTDIR/test.4.in.txt  -c testing_log.cfg  log_writer=binary
$SWARM query -f $OUTPUTDIR/testing_log.db > $OUTPUTDIR/testing_log.db.txt
$SWARM query -f $OUTPUTDIR/testing_log.bin > $OUTPUTDIR/testing_log.bin.txt


diff $OUTPUTDIR/testing_log.db.txt $TESTDIR/test.4.ref.txt && diff $OUTPUTDIR/testing_log.bin.txt $TESTDIR/test.4.ref.txt
