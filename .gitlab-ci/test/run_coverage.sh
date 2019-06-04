#!/bin/bash

ARCH=$1
OUTPUT_DIR=$2

. /etc/profile.d/modules.sh
module load gcc/4.9

mkdir $OUTPUT_DIR
lcov --no-checksum --capture --directory $ARCH --output-file coverage.info 
lcov -r coverage.info  ftimings* perf*  --output-file coverage.info
genhtml coverage.info --output-directory $OUTPUT_DIR


  
