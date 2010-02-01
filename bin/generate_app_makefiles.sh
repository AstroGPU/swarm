#!/bin/bash
# DO NOT CALL - called when needed by make, from bin/Makefile.d rule
#
# Usage: generate_app_makefiles.sh app1 app2 app3 ...
#

for i in "$@"
do
	echo "${i}_OBJECTS=\$(subst .cpp,.o,\$(${i}_SOURCES))"
	echo "bin/$i: bin/libswarm.so \$(${i}_OBJECTS)"
	echo -e "\t\$(LINKUI) \$(LINK) -o \$@ \$(${i}_OBJECTS)"
	echo
done;
