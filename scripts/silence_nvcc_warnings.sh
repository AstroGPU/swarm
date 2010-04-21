#!/bin/bash
#
# Silence annoying nvcc warnings that I seem to be unable to reliably turn off
#

grep -v "warning: statement is unreachable" | \
grep -v "Warning: Cannot tell what pointer points to, assuming global memory space" | \
grep -v "Loop was not unrolled, unexpected control flow construct" | \
sed '/^$/d'
