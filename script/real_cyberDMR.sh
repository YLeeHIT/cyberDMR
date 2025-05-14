#!/bin/bash
# Example:
#   bash ./real_cyberDMR.sh treatment control 1
group1=$1
group2=$2
threads=$3

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root_dir=$(dirname ${script_dir})
data_dir=${root_dir}/data/simulate_data/result/formatted_cyberDMR
echo "${data_dir}"
echo -e "\n##### Doing ${group1} and ${group2} with cyberDMR #####"
outcyber="${data_dir}/cyberDMR_result"
test -d ${outcyber} || mkdir ${outcyber}

cd ${data_dir}
echo "$(pwd)"
ls noh_*${group1}_*tsv noh_*${group2}_*tsv > cyberDMR_result/raw.lab 
awk -v dir=$(pwd) '{split($1,y,"_");split(y[5],x,".");print y[4]"_"NR"\t"y[3]"\t"dir"/"$1}' cyberDMR_result/raw.lab |sort -k2,2r > cyberDMR_result/in_cyber.lab

cd "cyberDMR_result"
py_script=${script_dir}/cyberDMR_main.py
start_time=$(date +%s)
time python ${py_script} --out_dir ${outcyber} --group1 ${group1} --group2 ${group2} --threads ${threads}
end_time=$(date +%s)
duration=$((end_time - start_time))

