#!/bin/bash

# Name of the integrator to do the stability test for
# The integrator config file should be in samples directory.
INTEGRATOR=$1

# Range of number of bodies to do stablity test for
# space sparated list of numbers
NBODY_RANGE="3 4 5 6"

INPUT_GENERATION_CONFIG=stability.cfg
INTEGRATOR_CFG=samples/$INTEGRATOR.cfg

DESTINATION_TIME=1e+4
INTERVAL=1e+3

[ ! -f "$INPUT_GENERATION_CONFIG" ] \
	&& echo "Config file $INPUT_GENERATION_CONFIG does not exist" \
	&& return 1
[ ! -f "$INTEGRATOR_CFG" ] \
	&& echo "Integrator config file not found: $INTEGRATOR_CFG" \
	&& return 1

echo "# Generator config: "
cat $INPUT_GENERATION_CONFIG
echo
echo "# Integrator config [$INTEGRATOR_CFG]: "
cat $INTEGRATOR_CFG


SHORTDATE='date +%Y-%m-%d'
DATE='date +%Y-%m-%d,%H:%m:%S'
OUTDIR=${INTEGRATOR}_stability_`$SHORTDATE`
mkdir -p $OUTDIR

for i in $NBODY_RANGE
do
	echo "ITERATION $i " `$DATE` "-----------------------------------------------------------------"
	INPUT_FILE=input_stability_${i}.bin
	OUTPUT_FILE=$OUTDIR/output_${i}.bin
	LOG_FILE=$OUTDIR/log_${i}.csv

	[ ! -f "$INPUT_FILE" ]  && echo "Generating input file from $INPUT_GENERATION_CONFIG" \
		&& bin/swarm generate -c $INPUT_GENERATION_CONFIG -o $INPUT_FILE nbod=$i
	echo "Stability test for $INTEGRATOR with $i bodies from input file $INPUT_FILE"
	echo -e "\tLog Output is written to $LOG_FILE"
	echo -e "\tEnsemble snapshot after integration is written to $OUTPUT_FILE"
	bin/swarm integrate -i $INPUT_FILE -o $OUTPUT_FILE \
		-d $DESTINATION_TIME  -n $INTERVAL \
		-c $INTEGRATOR_CFG > $LOG_FILE
	echo "Integration is done "	`$DATE`
done
