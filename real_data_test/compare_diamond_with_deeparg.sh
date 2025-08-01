#!/bin/bash

check_mapping="../../../../test/check_mapping.py"
mapfile -t biosamples < ./real_samples.txt
for biosample in "${biosamples[@]}"; do
    cd $biosample/deeparg_results/arg_alignment_identity_30
    diamond view -a X.align.daa -o X.align.daa.tsv 
    python $check_mapping
    cd ../arg_alignment_identity_50
    diamond view -a X.align.daa -o X.align.daa.tsv 
    python $check_mapping
    cd ../arg_alignment_identity_80
    diamond view -a X.align.daa -o X.align.daa.tsv 
    python $check_mapping
    cd ../../../
done