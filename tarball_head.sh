#!/bin/sh

mkdir -p tarballs
git archive --prefix=swarm-`git rev-parse --short HEAD`/ HEAD | gzip > tarballs/swarm-`git rev-parse --short HEAD`.tar.gz
