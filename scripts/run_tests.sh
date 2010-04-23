#!/bin/bash
#
# Usage: scripts/run_tests.sh <test|create>
#
#
# Sets up the output directories (in test-outputs), the environment, and
# runs the .test (or .create) scripts found in test/*/ subdirectories.  It
# redirects stdout and stderr of each script to a $scripbasename.log file in
# the output directory.
#
# Returns nonzero (logical false) if any test fails, 0 otherwise
#

if [ "$1" != "test" -a "$1" != "create" ]; then
	cat <<EOT
Usage: ./scripts/run_tests.sh <test|create>
EOT
exit -1
fi

#
# Remove any extisting test outputs
#
rm -rf test-outputs && mkdir test-outputs || { echo "Failed to create 'test-outputs' output directory. Aborting" && exit -1; }

FAIL=0

# Loop through all tests in test subdirectory, creating
# a directory in test-outputs for each one, and executing it
# from there.
for exe in $(cd test && ls -B */*.$1 2>/dev/null)
do
	test -x "test/$exe" || continue;

	export PREFIX=`pwd`
	export RUN="$PREFIX/test-outputs/$exe"
	export TEST=`dirname "$PREFIX/test/$exe"`
	export LOG="$RUN.log"

	echo -n "[ TEST ] $exe ... "

	(
		mkdir -p "$RUN" && cd "$RUN" &&
		"$PREFIX/test/$exe" > "$LOG" 2>&1
	) && {
		echo "PASSED";
	} || { 
		echo "FAILED   (output in $LOG)";
		FAIL=1; 
	}
done

exit $FAIL
