#!/bin/bash

g1=$1
g2=$2
rootdir=$3
threads=$4
indir=${rootdir}/input
outdir=${rootdir}/output
labdir=${rootdir}/group

# Step one: construct input file
echo "### Step one: construct input file ###"
g1_ID=$(cat ${labdir}/${g1})
g2_ID=$(cat ${labdir}/${g2})
test -d ${indir} && echo "input dir exist" || mkdir -p ${indir}
cd ${rootdir}/raw
metilene_input.pl --in1 ${g1_ID} --in2 ${g2_ID} --h1 ${g1} --h2 ${g2} --out ${indir}/${g1}_vs_${g2}.bed

# Step two: find DMR
echo "### Step two: find DMR ###"
test -d ${outdir} && echo "output dir exist" || mkdir -p ${outdir}
infile=${indir}/${g1}_vs_${g2}.bed
outfile=${outdir}/${g1}_vs_${g2}_DMRs.txt
maxdist=300
mincpgs=5
minDMR=0.1
## finding DMRs

metilene -M ${maxdist} -m ${mincpgs} -d ${minDMR} -t ${threads} -f 1 -a ${g1} -b ${g2} ${infile} |sort -V -k1,1 -k2,2n |\
            awk 'BEGIN{print "chr\tstart\tstop\tq-value\tdelta\tnum\tpMWU\tp2D\tmeang1\tmeang2"}{print $0}'> ${outfile}

# Step three: filter
echo "### Step three: filter ###"
ffile=${outdir}/${g1}_vs_${g2}_DMRs.filter.txt

awk 'BEGIN{print "chr\tstart\tstop\tq-value\tdelta\tnum\tpMWU\tp2D\tmeang1\tmeang2"}
    {if($4<0.05)print $0}' ${outfile} > ${ffile}

echo -e "### END ###"

