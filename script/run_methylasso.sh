#!/bin/bash
source ~/miniconda3/bin/activate DM-methylasso
indir=$1
methylasso_result=$2
test -d ${methylasso_result} || mkdir ${methylasso_result}

cd ${indir}
methylasso_r="~/software/methylasso/methylasso-main/MethyLasso.R"
group1="treatment"
group2="control"

start_time=$(date +%s)
time Rscript ${methylasso_r} \
    --n1 $group1 \
    --c1 $(ls noh_sorted_treatment_sample*tsv |tr '\n' ','|sed 's/,$//') \
    --n2 $group2 \
    --c2 $(ls noh_sorted_control_sample*tsv |tr '\n' ','|sed 's/,$//') \
    --cov 4 \
    --meth 5 \
    --q 0.05 \
    -c 5 \
    -d 0.1 \
    -n 5 \
    -o ${methylasso_result}

end_time=$(date +%s)
duration=$((end_time - start_time))

outname="${group1}_vs_${group2}"
awk 'NR>1{print $1"\t"$2"\t"$3"\t"$9"\t"$11}' ${methylasso_result}/${outname}_dmrs.tsv > ${methylasso_result}/${outname}.bed
