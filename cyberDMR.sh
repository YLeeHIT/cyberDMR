#!/bin/bash
set -euo pipefail

# -----------------------------
# Function: print help message
# -----------------------------
print_help() {
cat <<'USAGE'
Usage:
    bash cyberDMR.sh --indir/--inlab <PATH/TSV> --out-dir <PATH> --group1 <STR> --group2 <STR> [options]

Required arguments:
    -o,  --out-dir       Output directory.
    -g1, --group1        Group1 label.
    -g2, --group2        Group2 label.
    One of the following must be provided:
        -i,   --in-dir      Input directory containing the TSV files (used to auto-generate cyber.lab)
        -lab, --cyber-lab   Path to an existing cyber.lab file

Optional arguments (defaults match cyberDMR.py):
    -t,    --threads        Number of worker processes (default: 8).
    -chr,  --chroms         Chromosome set specification (default: chr1-chr22,chrX,chrY).
    -d,    --delta          Delta threshold for DMR detection (default: 0.1).
    -bdis, --cpg-distance   Maximum CpG distance for blocking (default: 500).
    -ct,   --cpg-count      Minimum number of CpGs per block (default: 5).
    -cov,  --min-cov        Minimum CpG coverage to keep (default: 5).
    -fdis, --max-dist       Maximum distance of adjacent CpGs (default: 500).
    -lab,  --cyber-lab      Path to the cyber.lab file (optional).
    -q,    --qvalue         BH-corrected p-value threshold (default: 0.05).
    -f,    --Fvalue         F statistic threshold (default: 15).

Other:
    -p <python_bin>         Python interpreter (default: python3).
    -h, --help              Show this help message.

Example:
    bash cyberDMR.sh -i ./data/samples_002 -o ./result -g1 treatment -g2 control -t 16 -d 0.1 -q 0.05
    bash cyberDMR.sh -lab ./my_samples.lab -o ./result -g1 treatment -g2 control
USAGE
}

# -----------------------------
# Default values
# -----------------------------
PYTHON_BIN="python3"
threads=8
chroms="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
delta=0.1
cpg_distance=500
cpg_count=5
min_cov=5
max_dist=500
cyber_lab=""
qvalue=0.05
Fvalue=15
indir=""
outdir=""
group1=""
group2=""

# -----------------------------
# Parse arguments
# -----------------------------
ARGS=()
while [[ $# -gt 0  ]]; do
    case "$1" in
        -i|--in-dir) indir="$2"; shift 2;;
        -o|--out-dir) outdir="$2"; shift 2;;
        -g1|--group1) group1="$2"; shift 2;;
        -g2|--group2) group2="$2"; shift 2;;
        -t|--threads) threads="$2"; shift 2;;
        -chr|--chroms) chroms="$2"; shift 2;;
        -d|--delta) delta="$2"; shift 2;;
        -bdis|--cpg-distance) cpg_distance="$2"; shift 2;;
        -ct|--cpg-count) cpg_count="$2"; shift 2;;
        -cov|--min-cov) min_cov="$2"; shift 2;;
        -fdis|--max-dist) max_dist="$2"; shift 2;;
        -lab|--cyber-lab) cyber_lab="$2"; shift 2;;
        -q|--qvalue) qvalue="$2"; shift 2;;
        -f|--Fvalue) Fvalue="$2"; shift 2;;
        -p) PYTHON_BIN="$2"; shift 2;;
        -h|--help) print_help; exit 0;;
        --) shift; ARGS+=("$@"); break;;
        *) echo "[ERROR] Unknown argument: $1"; print_help; exit 1;;
  esac
done


# -----------------------------
# Validate required arguments
# -----------------------------
if [[ -z "$outdir" || -z "$group1" || -z "$group2"  ]]; then
    echo "[ERROR] --out-dir, --group1, and --group2 are required."
    print_help
    exit 1
fi

if [[ -z "$cyber_lab" && -z "$indir"  ]]; then
    echo "[ERROR] Must provide either --cyber-lab (lab file) OR --indir (input directory)."
    exit 1
fi

[[ -n "$indir"  ]] && indir="$(readlink -f "$indir")"
outdir="$(readlink -f "$outdir")"


# -----------------------------
# Resolve paths
# -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cyberDMR_py="${SCRIPT_DIR}/lib/cyberDMR.py"

echo "[INFO] SCRIPT_DIR = ${SCRIPT_DIR}"
echo "[INFO] cyberDMR.py = ${cyberDMR_py}"

# -----------------------------
# Prepare output dir
# -----------------------------
mkdir -p "${outdir}"


# -----------------------------
# Helper: validate provided lab file
#   Format: sample_id <TAB> group_label <TAB> /abs/or/rel/path/to/tsv
#   Rules:
#     - Column 2 unique == 2 groups
#     - Column 3 paths exist and are readable
# -----------------------------
validate_lab() {
    local lab_file="$1"

    if [[ ! -f "$lab_file"  ]]; then
        echo "[ERROR] Provided lab file not found: $lab_file" >&2
        return 1
    fi

    local n_groups
    n_groups=$(awk -F'\t' 'NF>=3{print $2}' "$lab_file" | sort -u | wc -l | tr -d ' ')
    if [[ "$n_groups" -ne 2  ]]; then
        echo "[ERROR] Lab file must contain exactly 2 unique group labels in column 2; got: $n_groups" >&2
        echo "[HINT] Use: cut -f2 '$lab_file' | sort -u  file lab" >&2
        return 1
    fi

    local base_dir
    base_dir="$(cd "$(dirname "$lab_file")" && pwd)"

    local line=0
    while IFS=$'\t' read -r sid glab pth _rest; do
        ((line++)) || true
        [[ -z "$sid" && -z "$glab" && -z "$pth"  ]] && continue
        if [[ -z "$pth"  ]]; then
            echo "[ERROR] Line $line: empty path (col3) in lab file." >&2
            return 1
        fi
        
        if [[ "$pth" = /*  ]]; then
            file_path="$pth"
        else
            file_path="$base_dir/$pth"
        fi
        if [[ ! -r "$file_path"  ]]; then
            echo "[ERROR] Line $line: file not readable: $file_path" >&2
            return 1
        fi
    done < "$lab_file"

    echo "[INFO] Lab file validation passed: $lab_file"
    return 0
}

# -----------------------------
# If user provided --cyber-lab, validate and use it
# Else, generate ${outdir}/in_cyber.lab by scanning indir
# -----------------------------
if [[ -n "${cyber_lab:-}"  ]]; then
    echo "[INFO] Using user-provided lab: $cyber_lab"
    validate_lab "$cyber_lab" || exit 1
    outlab="$cyber_lab"
else
    if [[ ! -d "$indir"  ]]; then
        echo "[ERROR] Input directory not found: $indir" >&2
        exit 1
    fi

    outlab="${outdir}/in_cyber.lab"
    : > "$outlab"  # truncate

    echo "[INFO] Building lab from indir=$indir with group1=$group1 group2=$group2"
    pushd "$indir" >/dev/null || { echo "[ERROR] Cannot cd into ${indir}"; exit 1;  }

    shopt -s nullglob      

    i=1
    g1_matches=(*"${group1}"*)
    for file in "${g1_matches[@]:-}"; do
        [[ -f "$file"  ]] || continue
        echo -e "${group1}_${i}\t${group1}\t$(pwd)/$file" >> "$outlab"
        ((i++))
    done

    j=1
    g2_matches=(*"${group2}"*)
    for file in "${g2_matches[@]:-}"; do
        [[ -f "$file"  ]] || continue
        echo -e "${group2}_${j}\t${group2}\t$(pwd)/$file" >> "$outlab"
        ((j++))
    done

    shopt -u nullglob nocaseglob
    popd >/dev/null

    g1_n=$(awk -F'\t' -v g="$group1" 'NF>=3 && $2==g{c++} END{print c+0}' "$outlab")
    g2_n=$(awk -F'\t' -v g="$group2" 'NF>=3 && $2==g{c++} END{print c+0}' "$outlab")
    if [[ "$g1_n" -eq 0 || "$g2_n" -eq 0  ]]; then
        echo "[ERROR] No matching TSVs found. group1=$group1 count=$g1_n; group2=$group2 count=$g2_n" >&2
        echo "[HINT] Expected filenames like *${group1}*tsv and *${group2}*tsv under: $indir" >&2
        exit 1
    fi

    #validate_lab "$outlab" || exit 1
fi

# -----------------------------
# Run cyberDMR
# -----------------------------
#echo "[INFO] Running cyberDMR ..."
#set -x
#"${PYTHON_BIN:-python3}" "$cyberDMR_py" \
#    --out-dir "${outdir}" \
#    --group1 "${group1}" \
#    --group2 "${group2}" \
#    --threads "${threads}" \
#    --chroms "${chroms:-chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY}" \
#    --delta "${delta:-0.1}" \
#    --cpg-distance "${cpg_distance:-500}" \
#    --cpg-count "${cpg_count:-5}" \
#    --min-cov "${min_cov:-5}" \
#    --max-dist "${max_dist:-500}" \
#    --qvalue "${qvalue:-0.05}" \
#    --Fvalue "${Fvalue:-15}" \
#    --cyber-lab "${outlab}"
#set +x

# -----------------------------
# Post-processing: merge & sort result (if存在 chr*.txt)
# -----------------------------
#if compgen -G "${outdir}/chr"*".txt" > /dev/null; then
#    echo "[INFO] Merging chr*.txt -> ${outdir}/cyberDMR_result.txt"
#    cat "${outdir}/chr"*".txt" | sort -k1,1V -k2,2n -k3,3n > "${outdir}/cyberDMR_result.txt"
#else
