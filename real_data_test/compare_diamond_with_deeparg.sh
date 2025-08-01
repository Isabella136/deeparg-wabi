#!/bin/bash

check_mapping="../../../../test/check_mapping.py"
summarize="../../summarize_deeparg_ls_mapping_differences.py"
mapfile -t biosamples < ./real_samples.txt
for biosample in "${biosamples[@]}"; do
    cd $biosample/deeparg_results/arg_alignment_identity_30
    python $check_mapping
    cd ../arg_alignment_identity_50
    python $check_mapping
    cd ../arg_alignment_identity_80
    python $check_mapping
    cd ../../spades/deeparg_results/arg_alignment_identity_30
    python ../$check_mapping
    cd ../arg_alignment_identity_50
    python ../$check_mapping
    cd ../arg_alignment_identity_80
    python ../$check_mapping
    cd ../../
    python $summarize
    cd ../../
done