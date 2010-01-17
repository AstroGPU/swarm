#!/bin/bash
# Usage: mkcudafile.sh file1.cu file2.cu ... > single_file.cu
#

for i in "$@"
do
	echo "#include \"$i\""
done;
