#!/bin/bash

# -----------------------------
# Script to prepare input file and run cyberDMR
# Usage:
# bash cyberDMR.sh <indir> <outdir> <group1> <group2> <threads>
# Example:
# bash cyberDMR.sh /data/bed_files /results lethal normal 8
# -----------------------------

# Function: print help message
print_help() {
    echo "Usage: bash cyberDMR.sh <indir> <outdir> <group1> <group2> <threads>"
    echo ""
    echo "Arguments:"
    echo " indir Path to directory containing input BED files"
    echo " outdir Output directory to save results and in_cyber.lab"
    echo " group1 Group name for first condition (e.g., 'lethal')"
    echo " group2 Group name for second condition (e.g., 'normal')"
    echo " threads Number of threads to use in analysis"
    echo ""
    echo "Description:"
    echo " This script scans the input directory for BED files containing"
    echo " the group names, generates the in_cyber.lab file, and runs cyberDMR."
    echo ""
    echo "Example:"
    echo " bash cyberDMR.sh ./input ./output lethal normal 8"
}

# If user requests help
if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    print_help
    exit 0
fi

# -----------------------------
# Parse input arguments
# -----------------------------
indir=$1
outdir=$2
group1=$3
group2=$4
threads=$5

# -----------------------------
# Check argument validity
# -----------------------------
if [[ -z "$indir" || -z "$outdir" || -z "$group1" || -z "$group2" || -z "$threads" ]]; then
    echo "[ERROR] Missing arguments. Use --help for usage."
    exit 1
fi

# -----------------------------
# Get script directory & Python script path
# -----------------------------
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/script"
cyberDMR_py="${script_dir}/cyberDMR_main.py"

echo "[INFO] script_dir is ${script_dir}"

# Create output directory if not exists
test -d "${outdir}" || mkdir "${outdir}"

# Prepare in_cyber.lab
outlab="${outdir}/in_cyber.lab"
> "$outlab"

cd "${indir}" || { echo "[ERROR] Cannot cd into ${indir}"; exit 1; }

# Process group1 files
i=1
for file in *"${group1}"*bed; do
    echo -e "${group1}_${i}\t${group1}\t$(pwd)/$file" >> "$outlab"
    ((i++))
done

# Process group2 files
j=1
for file in *"${group2}"*bed; do
    echo -e "${group2}_${j}\t${group2}\t$(pwd)/$file" >> "$outlab"
    ((j++))
done

# -----------------------------
# Run cyberDMR
# -----------------------------
cd "${outdir}" || exit 1

echo "[INFO] Running cyberDMR..."
# conda activate DM-cyberDMR
python "${cyberDMR_py}" --out_dir "${outdir}" --threads "${threads}" --group1 "${group1}" --group2 "${group2}"
cat ./chr*txt |sort -k1,1V -k2,2n -k3,3n > cyberDMR_result.txt

echo "[INFO] Finished successfully."
