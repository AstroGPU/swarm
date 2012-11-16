#/bin/bash
rm -f feedback.tgz
rm -rf Testing
make clean 2>&1 | tee feedback.make.clean
make swarm 2>&1 | tee feedback.make.swarm
make test  2>&1 | tee feedback.make.test
tar -czf feedback.tgz feedback.make.* Testing
rm -f feedback.make.*
echo "Created feedback.tgz"

