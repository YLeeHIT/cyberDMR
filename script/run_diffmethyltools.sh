#!/bin/bash

source ~/miniconda3/bin/activate DiffMethylTools
indir=$1
diffmethyltools_result=$2
test -d ${diffmethyltools_result} || mkdir ${diffmethyltools_result}

cd ${indir}
diffmethyltools_py="~/software/DiffMethylTools/DiffMethylTools.py"
group1="treatment"
group2="control"

start_time=$(date +%s)
python ${diffmethyltools_py} all_analysis \
    --case_data_file $(ls noh_sorted_treatment_sample*tsv) \
    --ctr_data_file $(ls noh_sorted_control_sample*tsv) \
    --case_data_chromosome_column_index 0 \
    --ctr_data_chromosome_column_index 0 \
    --case_data_position_start_column_index 1 \
    --ctr_data_position_start_column_index 1 \
    --case_data_positive_methylation_count_column_index 3 \
    --case_data_negative_methylation_count_column_index 4 \
    --ctr_data_positive_methylation_count_column_index 3 \
    --ctr_data_negative_methylation_count_column_index 4 \
    --case_data_separator $'\t' \
    --ctr_data_separator $'\t' \
    --max_q_value 0.05 \
    --abs_min_diff 0.1 \
    2>"diffmethyltools_error.log"

end_time=$(date +%s)
duration=$((end_time - start_time))

outfile="${indir}/data/generate_DMR_0.csv"
outname="${group1}_vs_${group2}"
awk -F"," -v OFS="\t" 'NR>1{print $1,$2,$3,$4,$6}' ${outfile} > ${diffmethyltools_result}/${outname}.bed
