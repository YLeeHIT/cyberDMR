#!/bin/bash

# ================================================================
# Script: simulate_data.sh
# Description: Simulate DMR regions and run cyberDMR detection pipeline
# ================================================================

# -------------------------------
# Default parameters
# -------------------------------
total_dmr=1000
mean_delta=0.25
n_control=10
n_treatment=10
coverage_mean=30
coverage_std=5
output_dir="$(pwd)/output"
chr_name="chr1"
start_pos=100000
length_mean=1000
length_std=100
max_cpgs=200
dmr_per=0.25
dmr_notable_per=0.02
dmr_inconsis_per=0.03
dmr_sub_per=0.05
density="moderate"
min_gap=10
max_gap=50
sample_missing_ratio=0.1
cpg_missing_ratio=0.1
good_precision=20
bad_precision=20
no_delta_methylation=0.08
attempt_per_slot=200
seed=42
threads=1

#SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${CYBERDMR_HOME:-}" && -d "$CYBERDMR_HOME"   ]]; then
    SCRIPT_DIR="$CYBERDMR_HOME"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "$SLURM_SUBMIT_DIR"   ]]; then
    SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
    _src="${BASH_SOURCE[0]:-${0}}"
    while [ -L "$_src"   ]; do
        _dir="$(cd -P "$(dirname "$_src")" && pwd)"
        _src="$(readlink "$_src")"
        [[ "$_src" != /*   ]] && _src="$_dir/$_src"
    done
    SCRIPT_DIR="$(cd -P "$(dirname "$_src")" && pwd)"
fi

echo "[INFO] SLURM_SUBMIT_DIR = ${SLURM_SUBMIT_DIR:-N/A}"
echo "[INFO] SCRIPT_DIR       = ${SCRIPT_DIR}"

simulate_py="${SCRIPT_DIR}/script/simulated_data.py"
merge_sh="${SCRIPT_DIR}/script/merge_simulated_samples.sh"

# -------------------------------
# Help function
# -------------------------------
print_help() {
    echo "Usage: bash $0 [options]"
    echo "Options:"
    echo " -t, --total_dmr NUM Total number of simulated DMRs (default: $total_dmr)"
    echo " -d, --mean_delta NUM Mean methylation delta (default: $mean_delta)"
    echo " -c, --n_control NUM Number of control samples (default: $n_control)"
    echo " -e, --n_treatment NUM Number of treatment samples (default: $n_treatment)"
    echo " -m, --coverage_mean NUM Mean coverage depth (default: $coverage_mean)"
    echo " -s, --coverage_std NUM Coverage standard deviation (default: $coverage_std)"
    echo " -o, --output_dir PATH Output directory (default: $output_dir)"
    echo " -r, --chr_name STR Chromosome name (default: $chr_name)"
    echo " -p, --start_pos NUM Start position for DMR simulation (default: $start_pos)"
    echo " -l, --length_mean NUM Mean DMR length (default: $length_mean)"
    echo " -z, --length_std NUM Standard deviation of DMR length (default: $length_std)"
    echo " -x, --max_cpgs NUM Max CpGs per DMR (default: $max_cpgs)"
    echo " -q, --dmr_per NUM Proportion of good DMRs (default: $dmr_per)"
    echo " -n, --dmr_notable_per NUM Proportion of notable DMRs (default: $dmr_notable_per)"
    echo " -i, --dmr_inconsis_per NUM Proportion of inconsistent DMRs (default: $dmr_inconsis_per)"
    echo " -u, --dmr_sub_per NUM Proportion of sub DMRs (default: $dmr_sub_per)"
    echo " -y, --density STR Density mode: dense/ moderate / sparse (default: $density)"
    echo " -mn, --min_gap NUM Min distance between two cpgs (default: $min_gap)"
    echo " -mx, --max_gap NUM Max distance between two cpgs (default: $max_gap)"
    echo " -sm, --sample_missing NUM Missing ratio of samples (default: $sample_missing_ratio)"
    echo " -cm, --cpg_missing NUM Missing ratio of cpgs (default: $cpg_missing_ratio)"
    echo " -gp, --good_precisoin NUM Parameter of precision in Beta(alpha, beta) model for good-DMR (default: $good_precision)"
    echo " -bp, --bad_precisoin NUM Parameter of precision in Beta(alpha, beta) model for non-DMR (default: $bad_precision)"
    echo " -nd, --no_delta_methylation Methylation delta for non-DMR (default: $no_delta_methylation)"
    echo " -a, --attempt_per_slot NUM Max attempts for generating a region of a given class before skipping (default: $attempt_per_slot)"
    echo " -S, --seed NUM Random seed (default: $seed)"
    echo " -T, --threads NUM Number of threads for cyberDMR (default: $threads)"
    echo " -h, --help Show this help message and exit"
    echo "Example:"
    echo " bash $0 -o $(pwd)/test -t 100"
    exit 0
}

# -------------------------------
# Parse parameters
# -------------------------------
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
        -mn|--min_gap) min_gap="$2"; shift 2 ;;
        -mx|--mix_gap) max_gap="$2"; shift 2 ;;
        -sm|--sample_missing) sample_missing_ratio="$2"; shift 2 ;;
        -cm|--cpg_missing) cpg_missing_ratio="$2"; shift 2 ;;
        -gp|--good_precision) good_precision="$2"; shift 2 ;;
        -bp|--bad_precision) bad_precision="$2"; shift 2 ;;
        -nd|--no_delta_methylation) no_delta_methylation="$2"; shift 2 ;;
        -a|--attempt_per_slot) attempt_per_slot="$2"; shift 2 ;;
        -S|--seed) seed="$2"; shift 2 ;;
        -h|--help) print_help ;;
        *) echo "[ERROR] Unkown parameter: $1"; print_help ;;
    esac
done

# -------------------------------
# Step 1: Simulate DMRs
# -------------------------------
echo "Step 1: Simulate data"
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
    --min_gap "$min_gap" \
    --max_gap "$max_gap" \
    --sample_missing_max "$sample_missing_ratio" \
    --dmr_missing_max "$cpg_missing_ratio" \
    --good_precision "$good_precision" \
    --no_precision "$bad_precision" \
    --no_delta_methylation "$no_delta_methylation" \
    --max_attempts_per_slot "$attempt_per_slot" \
    --seed "$seed"

# -------------------------------
# Step 2: Merge and convert
# -------------------------------
extract_sh="/home/user/liyang/project/methDmr/scripts/extract_DMR.sh"
echo "Step 2: Merge data and convert format"
bash ${merge_sh} ${output_dir}
DMR_file="${output_dir}/DMRs.txt"
bash ${extract_sh} ${DMR_file} ${output_dir}

# -------------------------------
# Step 3: Detect DMR
# -------------------------------
detect_sh="/home/user/liyang/project/methDmr/scripts/detect_DMR.sh"
echo "Step Three: Detect DMRs"
#fvalue=500000
bash ${detect_sh} ${output_dir} ${chr_name} ${threads} all

# -------------------------------
# Step 4: Integrate result
# -------------------------------
integrate_sh="/home/user/liyang/project/methDmr/scripts/integrate_data.sh"
echo "Step Four: Integrate result"
bash ${integrate_sh} ${output_dir}

echo "[INFO] All process has finished"

