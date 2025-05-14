#!/bin/bash
# Example:
#   bash ./simulate_data.sh -o ../data/simulate_data -t 100

# Default parameter
total_dmr=10000
mean_delta=0.25
n_control=10
n_treatment=10
coverage_mean=30
coverage_std=5
output_dir="./output"
chr_name="chr1"
start_pos=100000
length_mean=1000
length_std=100
max_cpgs=100
dmr_per=0.19
dmr_notable_per=0.01
dmr_inconsis_per=0
dmr_sub_per=0
density="mix"
dense_ratio=0.3
seed=42
threads=1   # only use for cyberDMR
simulate_py="./simulated_data.py"
merge_sh="./merage_simulated_samples.sh"
detect_sh="./detect_DMR.sh"
integrate_sh="./integrate_data.sh"


# Help function
print_help() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo " -t, --total_dmr NUM Total DMRs (default: $total_dmr)"
    echo " -d, --mean_delta NUM Mean delta methylation (default: $mean_delta)"
    echo " -c, --n_control NUM Number of control samples (default: $n_control)"
    echo " -e, --n_treatment NUM Number of treatment samples (default: $n_treatment)"
    echo " -m, --coverage_mean NUM Mean coverage (default: $coverage_mean)"
    echo " -s, --coverage_std NUM Coverage standard deviation (default: $coverage_std)"
    echo " -o, --output_dir PATH Output directory (default: $output_dir)"
    echo " -r, --chr_name STR Chromosome name (default: $chr_name)"
    echo " -p, --start_pos NUM Start position (default: $start_pos)"
    echo " -l, --length_mean NUM DMR length mean (default: $length_mean)"
    echo " -z, --length_std NUM DMR length std (default: $length_std)"
    echo " -x, --max_cpgs NUM Max CpGs per DMR (default: $max_cpgs)"
    echo " -q, --dmr_per NUM Proportion of good DMRs (default: $dmr_per)"
    echo " -n, --dmr_notable_per NUM Proportion of notable DMRs (default: $dmr_notable_per)"
    echo " -i, --dmr_inconsis_per NUM Proportion of inconsistent DMRs (default: $dmr_inconsis_per)"
    echo " -u, --dmr_sub_per NUM Proportion of sub DMRs (default: $dmr_sub_per)"
    echo " -y, --density STR Density type: mix/dense/sparse (default: $density)"
    echo " -a, --dense_ratio NUM Dense region ratio (default: $dense_ratio)"
    echo " -S, --seed NUM Random seed (default: $seed)"
    echo " -h, --help Show help"
    exit 0
}

# Parameter
while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--total_dmr) total_dmr="$2"; shift 2 ;;
        -d|--mean_delta) mean_delta="$2"; shift 2 ;;
        -c|--n_control) n_control="$2"; shift 2 ;;
        -e|--n_treatment) n_treatment="$2"; shift 2 ;;
        -m|--coverage_mean) coverage_mean="$2"; shift 2 ;;
        -s|--coverage_std) coverage_std="$2"; shift 2 ;;
        -o|--output_dir) output_dir="$2"; shift 2 ;;
        -r|--chr_name) chr_name="$2"; shift 2 ;;
        -p|--start_pos) start_pos="$2"; shift 2 ;;
        -l|--length_mean) length_mean="$2"; shift 2 ;;
        -z|--length_std) length_std="$2"; shift 2 ;;
        -x|--max_cpgs) max_cpgs="$2"; shift 2 ;;
        -q|--dmr_per) dmr_per="$2"; shift 2 ;;
        -n|--dmr_notable_per) dmr_notable_per="$2"; shift 2 ;;
        -i|--dmr_inconsis_per) dmr_inconsis_per="$2"; shift 2 ;;
        -u|--dmr_sub_per) dmr_sub_per="$2"; shift 2 ;;
        -y|--density) density="$2"; shift 2 ;;
        -a|--dense_ratio) dense_ratio="$2"; shift 2 ;;
        -S|--seed) seed="$2"; shift 2 ;;
        -h|--help) print_help ;;
        *) echo "Unkown parameter: $1"; print_help ;;
    esac
done

# Step One: simulate data
echo "Step One: Simulate data"
python "${simulate_py}" \
    --total_dmr "$total_dmr" \
    --mean_delta "$mean_delta" \
    --n_control "$n_control" \
    --n_treatment "$n_treatment" \
    --coverage_mean "$coverage_mean" \
    --coverage_std "$coverage_std" \
    --output_dir "$output_dir" \
    --chr_name "$chr_name" \
    --start_pos "$start_pos" \
    --length_mean "$length_mean" \
    --length_std "$length_std" \
    --max_cpgs "$max_cpgs" \
    --dmr_per "$dmr_per" \
    --dmr_notable_per "$dmr_notable_per" \
    --dmr_inconsis_per "$dmr_inconsis_per" \
    --dmr_sub_per "$dmr_sub_per" \
    --density "$density" \
    --dense_ratio "$dense_ratio" \
    --seed "$seed"

parafile="${output_dir}/para.log"
awk 'NR>1{len=$3-$2;print len"\t"$4}' ${output_dir}/DMRs.txt | \
    datamash mean 1,2 |awk '{print "\nSimulated result:\nMean(len)\t"$1"\nMean(CpGs)\t"$2}' >> ${parafile}

# Step Two: merge data
echo "Step Two: Merge data and convert format"
bash ${merge_sh} ${output_dir}

# Step Three: detect DMR
echo "Step Three: Detect DMRs"
#bash ${detect_sh} ${output_dir} ${chr_name} ${threads} all

# Step Four: integrate result
echo "Step Four: Integrate result"
#bash ${integrate_sh} ${output_dir}

echo "##### All process has finished #####"

