#/bin/csh
rm -f feedback.tgz;
make clean >& feedback.make.clean
make tidy >& feedback.make.tidy
make info >& feedback.make.info
make apps >& feedback.make.apps
make doc-asciidoc >& feedback.make.asciidoc
make doc-doxygen >& feedback.make.doxygen
make test >& feedback.make.test
make benchmark-quick |& tee feedback.make.benchmark
tar czf feedback.tgz feedback.make.* run/benchmark.out test-outputs
rm -f feedback.make.*
echo "Created feedback.tgz"

