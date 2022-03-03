#!/usr/bin/env bash
source activate cellbender

# Move to directory containing raw matrix in mtx format
cd /path/to/workdir
cd CELLBENDER

# Make sure the matrix, barcodes and genes are same as how cellranger outputs them
# Loop over samples (each in its own directory) and run CellBender
for dir in `ls -d $PWD/*`
do

cellbender remove-background \
     --input $dir \
     --output $dir"/CellBender_out.h5" --cuda

done





