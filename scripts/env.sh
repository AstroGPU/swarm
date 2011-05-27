#!/bin/bash

pushd $(dirname $0) > /dev/null
export SWARMDIR=`pwd`
popd > /dev/null
export PATH=$PATH:$SWARMDIR/bin:$SWARMDIR/scripts
. $SWARMDIR/Makefile.user
