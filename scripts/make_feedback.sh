#/bin/csh
rm -f feedback.tgz;
make info >& feedback.make.info
make apps >& feedback.make.apps
make doc >& feedback.make.doc
make test >& feedback.make.test
make benchmark >& feedback.make.benchmark
tar czf feedback.tgz feedback.* run/benchmark.out test-outputs
rm -f feedback.make.*
echo "Created feedback.tgz"

