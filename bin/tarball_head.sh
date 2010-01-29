#!/bin/sh

mkdir -p tarballs

branch=`git branch --no-color 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/\1/'`
commit=`git rev-parse --short HEAD`
prefix="swarm-$branch-$commit"

git archive --prefix=$prefix/ HEAD | gzip > tarballs/$prefix.tar.gz

echo Tarballed to tarballs/$prefix.tar.gz
