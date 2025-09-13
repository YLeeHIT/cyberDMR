#!/bin/bash

source ~/miniconda3/bin/activate HOMEenv
inhome=$1

outdir=$2
test -d ${outdir} || mkdir ${outdir}

home_master="~/software/HOME/HOME-master"
cat <(echo "treatment") <(ls -t ${inhome}/*treatment*tsv) |awk 'NR==1{line=$0; next}{printf "%s\t", line; line=$0} END{print line}' > ${outdir}/treatment.lab
cat <(echo "control") <(ls -t ${inhome}/*control*tsv) |awk 'NR==1{line=$0; next}{printf "%s\t", line; line=$0} END{print line}' > ${outdir}/control.lab
cat ${outdir}/treatment.lab ${outdir}/control.lab > ${home_master}/treatment_control.lab
cp ${home_master}/treatment_control.lab ${outdir}/treatment_control.lab

outname=HOME_loose
cd ${home_master}

start_time=$(date +%s)
time HOME-pairwise -t CG -mc 3 -d 0.1 -sc 0.1 -ml 50 -npp 1 -i treatment_control.lab -o ${outname}
end_time=$(date +%s)
duration=$((end_time - start_time))

home_result=${home_master}/${outname}/HOME_pairwise_DMRs/treatment_VS_control/HOME_DMRs*txt
awk 'NR>1{print "chr"$1"\t"$2"\t"$3"\t"$5"\t"$8}' ${home_result} > ${outdir}/${outname}.bed
