#!/bin/bash

export CUDAPATH=/opt/cuda4.0rc
export CUDAARCH=-gencode 'arch=compute_13,code=sm_13'

#pushd $(dirname $0)/.. > /dev/null
export SWARMDIR=`pwd`
#popd > /dev/null
export PATH=$PATH:$SWARMDIR/bin:$SWARMDIR/scripts

export LD_LIBRARY_PATH=$CUDAPATH/lib64:$CUDAPATH/lib
