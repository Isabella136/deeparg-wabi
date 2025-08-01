#!/bin/bash

calculate_mapping_frequencies="../../../calculate_mapping_frequencies.py"
mapfile -t biosamples < ./real_samples.txt
for biosample in "${biosamples[@]}"; do
    cd $biosample/deeparg_results/arg_alignment_identity_30
    python $calculate_mapping_frequencies
    cd ../arg_alignment_identity_50
    python $calculate_mapping_frequencies
    cd ../arg_alignment_identity_80
    python $calculate_mapping_frequencies
    cd ../../spades/deeparg_results/arg_alignment_identity_30
    python ../$calculate_mapping_frequencies
    cd ../arg_alignment_identity_50
    python ../$calculate_mapping_frequencies
    cd ../arg_alignment_identity_80
    python ../$calculate_mapping_frequencies
    cd ../../../../
done