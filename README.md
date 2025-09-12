# cyberDMR

![Version](https://img.shields.io/badge/version-1.1-blue)
![Language](https://img.shields.io/badge/language-python-blue)
![Language](https://img.shields.io/badge/language-shell-4EAA25)
![Language](https://img.shields.io/badge/language-R-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Platform](https://img.shields.io/badge/platform-linux%20|%20macOS-brightgreen)

**cyberDMR** is a robust and high-sensitivity approach for differentially methylated regions detection.

### Introduction
Differentially methylated regions (DMRs) are key genomic features reflecting changes in DNA methylation status. Accurate identification of DMRs is crucial for investigating tissue-specific regulation, disease mechanisms, and population-level epigenetic variation.

### Features
- Base-level smoothing for low-coverage CpGs
- CpG segmentation based on genomic distance and methylation concordance
- Seed-guided clustering for consistent CpG grouping
- Weighted beta regression with LRT for statistical inference
- Identifiying significant DMRs via BH correction and F-statitics

## 1. Installation
```bash
### Clone the repository
git clone https://github.com/YLeeHIT/cyberDMR.git
cd cyberDMR

### create a new conda environment
conda create -n DM-cyberDMR python=3.12 -y
conda activate DM-cyberDMR

### Install required dependencies
pip install -r requirements.txt
```

## 2. Usage

```bash
bash cyberDMR.sh --in-dir <indir> --out-dir <outdir> --group1 <group1> --group2 <group2> [<optional>]
```

Check all available options with:
```bash
bash cyberDMR.sh --help
```

For detailed parameter descriptions, see **3. Arguments**.
For usage examples, see **8. Demo**


## 3. Arguments

| Parameter               | Required | Description                               | Example                |
|------------------———————|----------|-------------------------------------------|------------------------|
| `-o, --out-dir`         |          | Output directory for storing all results  | `./results/`           |
| `-g1, --group1`         |          | Label of group 1 (e.g., treatment)        | `treatment`            |
| `-g2, --group2`         |          | Label of group 2 (e.g., control)          | `control`              |
| `-i, --in-dir`          |          | Input files (auto-generate `cyber.lab`)   | `./input/`             |
| `-lab, --cyber-lab`     |          | Path to an existing `cyber.lab` file      | `./cyber.lab`          |
| `-t, --threads`         |          | Number of worker processes                | `8`                    |
| `-chr, --chroms`        |          | Chromosome set specification              | `chr1-,chr2,chr3`      |
| `-d, --delta`           |          | Delta threshold for DMR detection         | `0.1`                  |
| `-bdis, --cpg-distance` |          | Maximum CpG distance for blocking         | `500`                  |
| `-ct, --cpg-count`      |          | Minimum number of CpGs per block          | `5`                    |
| `-cov, --min-cov`       |          | Minimum CpG coverage to fill              | `5`                    |
| `-fdis, --max-dist`     |          | Maximum distance of adjacent CpGs         | `500`                  |
| `-q, --qvalue`          |          | BH-corrected p-value threshold            | `0.05`                 |
| `-f, --Fvalue`          |          | F statistic threshold                     | `15`                   |

\* One of `--in-dir` or `--cyber-lab` must be provided.

---

### `--out-dir`
Supports both absolute and relative paths.  
This directory will store all output results, including per-chromosome files and the final merged and sorted file `cyberDMR_result.bed`.  

### `--group1`, `--group2`
Names of the two groups must be provided.  
**The experimental group should come first, followed by the control group**, to ensure consistent statistical comparison.  

### `--in-dir`
Supports both absolute and relative paths.  
Should point to the directory containing input files formatted.  
When this parameter is provided, the program will automatically generate an `in_cyber.lab` file. File names must follow strict naming conventions (see [Input](#input)).  

### `--cyber-lab`
If the user has already prepared a `lab` file that meets the **Input** requirements, it can be provided via this parameter instead of using `--in-dir`.  

### `--threads`
Number of worker processes.  
It is recommended to set this equal to the number of chromosomes for best performance.  

### `--delta`
Minimum methylation difference (Δ).  
DMRs with Δ below this threshold will be filtered out.  

### `--cpg-distance`
Maximum CpG distance for **blocking**.  
This parameter affects the blocking process. Suggested range: `300–1000` (default: `500`).  

### `--cpg-count`
Minimum number of CpGs per DMR block.  
Regions with fewer CpGs will be filtered out.  

### `--min-cov`
Minimum CpG coverage for smoothing:  
- Recommended `5` for WGBS data  
- Recommended `3` for ONT data  
When coverage falls below this threshold, smoothing will be applied.  

### `--max-dist`
Maximum distance between adjacent CpGs for **clustering**.  
This parameter affects the clustering process. Suggested range: `300–1000` (default: `500`).  

### `--qvalue`
Benjamini–Hochberg corrected p-value threshold.  
DMRs with q-values above this cutoff will be filtered out.  

### `--Fvalue`
F-statistic threshold.  
- Strict filtering: `20`  
- Relaxed filtering: `5`  

---

## 4. Input format

Before running `cyberDMR.sh`, you can provide the directory containing all sample files using the `--in-dir` option. In this case, cyberDMR will automatically generate the `in_cyber.lab` file.  
Alternatively, you can supply your own lab file with sample paths and grouping information using the `--lab` option. cyberDMR will also recognize this file and proceed with the analysis.

### Input File Requirements
- Input files should be tab-delimited text (`.tsv` or `.bed`-like format) without a header.  
- Each input file name must include the group label (e.g., `HG002_treatment.tsv`, `HG003_control.tsv`).
- Each file should contain exactly four columns in the following order:

1. **Chromosome** (`string`) – e.g., `chr22`  
2. **CpG position** (`integer`) – genomic coordinate (0-based or 1-based)  
3. **Methylation level** (`float`) – value between `0.0` and `1.0`  
4. **Coverage** (`integer`) – positive integer indicating read depth  

**Example** (`in_cyber.lab`):
```
chr1    107908  1.0     25
chr1    107977  1.0     40
chr1    107988  1.0     20
chr1    108918  0.5301  32
chr1    109368  0.5236  30
chr1    109545  0.675   24
chr1    110009  0.5276  33
chr1    113405  0.2748  32
chr1    113828  0.3616  25
chr1    113945  0.3926  31
```

### Lab File Format Requirements

- This file is used to define the grouping of biological replicates, their phenotypic labels, and the corresponding input files.  
- It must strictly follow the format below (tab-delimited, without a header):

1. **Sample ID** – unique identifier for each biological replicate  
2. **Group label** – e.g., `treatment` or `control` (only **two groups** are supported)  
3. **Absolute file path** – path to the input file (including the group label in the filename)  

**Example** (`in_cyber.lab`):
```
139C    lethal  /absolute/path/to/noh_lethal_139C_auto.bed
1601C   lethal  /absolute/path/to/noh_lethal_1601C_auto.bed
349C    lethal  /absolute/path/to/noh_lethal_349C_auto.bed
379C    lethal  /absolute/path/to/noh_lethal_379C_auto.bed
46C lethal  /absolute/path/to/noh_lethal_46C_auto.bed
514C    lethal  /absolute/path/to/noh_lethal_514C_auto.bed
564C    lethal  /absolute/path/to/noh_lethal_564C_auto.bed
1601N   normal  /absolute/path/to/noh_normal_1601N_auto.bed
448N    normal  /absolute/path/to/noh_normal_448N_auto.bed
508N    normal  /absolute/path/to/noh_normal_508N_auto.bed
564N    normal  /absolute/path/to/noh_normal_564N_auto.bed
```

**Note:** Ensure all paths are absolute (not relative), and that group names match the `--group1` and `--group2` arguments when running `cyberDMR.py`.
Once ready, you can run cyberDMR as follows:


## 5. Simulate DMR regions
To generate simulated DMR regions for benchmarking:

```bash
python simulated_data.py \
    --total_dmr 1000 \
    --mean_delta 0.3 \
    --n_control 5 \
    --n_treatment 5 \
    --coverage_mean 30 \
    --coverage_std 5 \
    --output_dir /mnt/data/sample_outputs \
    --chr_name chr1 \
    --start_pos 10000 \
    --length_mean 1000 \
    --length_std 300 \
    --max_cpgs 50 \
    --dmr_per 0.3 \
    --dmr_notable_per 0.05 \
    --dmr_inconsis_per 0.1 \
    --dmr_sub_per 0.05 \
    --density auto \
    --dense_ratio 0.5 \
    --seed 42
```

To generate input files in formats compatible with **cyberDMR**, **Metilene**, **BSmooth**, and **HOME**, run the provided merging script:

```bash
bash merge_simulates_samples.sh -o ../data/simulate_data
```

### Quick Start (Recommended)

To simplify everything, you can run the pre-configured shell script:

```bash
bash ./simulate_data.sh [options]
```

| Parameter               | Required | Description                                               | Default       |
|-------------------------|----------|-----------------------------------------------------------|---------------|
| `--output_dir`          | ✅       | Output directory to store simulated data and results      | *(no default)*|
| `--total_dmr`           | ❌       | Total number of DMR regions to simulate                   | `10000`       |
| `--mean_delta`          | ❌       | Average methylation difference between groups             | `0.25`        |
| `--n_control`           | ❌       | Number of control samples                                 | `10`          |
| `--n_treatment`         | ❌       | Number of treatment samples                               | `10`          |
| `--coverage_mean`       | ❌       | Mean sequencing coverage                                  | `30`          |
| `--coverage_std`        | ❌       | Standard deviation of coverage                            | `5`           |
| `--chr_name`            | ❌       | Chromosome name to simulate DMRs                          | `chr1`        |
| `--start_pos`           | ❌       | Start position for simulation                             | `100000`      |
| `--length_mean`         | ❌       | Mean DMR region length                                    | `1000`        |
| `--length_std`          | ❌       | Standard deviation of DMR length                          | `100`         |
| `--max_cpgs`            | ❌       | Maximum number of CpGs per DMR                            | `100`         |
| `--dmr_per`             | ❌       | Proportion of good DMRs                                   | `0.19`        |
| `--dmr_notable_per`     | ❌       | Proportion of notable DMRs                                | `0.01`        |
| `--dmr_inconsis_per`    | ❌       | Proportion of inconsistent DMRs                           | `0`           |
| `--dmr_sub_per`         | ❌       | Proportion of sub DMRs                                    | `0`           |
| `--density`             | ❌       | CpG density type: `mix`, `dense`, or `sparse`             | `mix`         |
| `--dense_ratio`         | ❌       | Proportion of dense regions (only applies if `mix`)       | `0.3`         |
| `--seed`                | ❌       | Random seed                                               | `42`          |
| `--threads`             | ❌       | Number of threads used by cyberDMR                        | `1`           |

You can use the provided script to automatically generate the input file (`in_cyber.lab`) and run `cyberDMR`.



View help information:

```bash
bash simulate_data.sh -h
bash cyberDMR.sh -h
```

## 6. Example:

```
bash simulate_data.sh -o ../data/simulate_data -t 100
bash cyberDMR.sh ./data/real_data/chr22 ./data/real_data/chr22/cyberDMR_result lethal normal 8
```


## Release Notes

### Release Notes – cyberDMR v1.0

**Release Date:** 2025-05-13
**Status:** Initial release


### Release Notes – cyberDMR v1.1

**Release Date:** 2025-08-5
**Status:** Initial release

- Fixed the "Maximum Likelihood optimization failed" error in certain edge cases during model fitting.
- Added simulated datasets for multiple scenarios to demonstrate tool behavior under different conditions.
- Expanded and clarified usage instructions.


If you use cyberDMR in your research, please cite:

If you use **cyberDMR** in your research, please cite the following paper:

> **Li, Yang**, *et al.*
> **cyberDMR: a robust and high-sensitivity approach for differentially methylated regions detection**
> *Bioinformatics*, 2025 (under review)
> [GitHub Project](https://github.com/YLeeHIT/cyberDMR)

BibTeX:

```bibtex
@article{li2025cyberdmr,
    title = {cyberDMR: a robust and high-sensitivity approach for differentially methylated regions detection},
    author = {Li, Yang and others},
    journal = {Bioinformatics},
    year = {2025},
    note = {Manuscript under review}
}
```

We appreciate your support!

# Contributors

This package is developed and maintaned by [Lee](https://github.com/YLeeHIT) and [Chen](https://github.com/chong-hun). If you want to contribute, please leave an issue or submit a pull request. Thank you.

# License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
