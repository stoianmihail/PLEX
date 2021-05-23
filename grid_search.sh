#!/usr/bin/env bash

mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make && cd ..

SPLINE_ERROR="2 4 8 16 32 64 128 256 512 1024"
NUM_BINS="2 4 8 16 32 64 128 256 512 1024"
MAX_ERROR="2 4 8 16 32 64 128 256 512 1024"

DATA_SETS="../SOSD/data/books_200M_uint64"

for DATA_SET in $DATA_SETS; do
  for SPLINE in $SPLINE_ERROR; do
    for BIN in $NUM_BINS; do
      for ERROR in $MAX_ERROR; do
        ./build/bench_end_to_end $DATA_SET ${DATA_SET}_equality_lookups_10M $SPLINE $BIN $ERROR 0 0
        ./build/bench_end_to_end $DATA_SET ${DATA_SET}_equality_lookups_10M $SPLINE $BIN $ERROR 0 1
        ./build/bench_end_to_end $DATA_SET ${DATA_SET}_equality_lookups_10M $SPLINE $BIN $ERROR 1 0
      done;  
    done;
  done;
done;
