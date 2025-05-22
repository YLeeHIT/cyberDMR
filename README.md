# cyberDMR

![Version](https://img.shields.io/badge/version-1.0.0-blue)
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

# Running cyberDMR

## Installation
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

## Simulate DMR regions
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

## Run cyberDMR
### Input requirement for cyberDMR

Before running `cyberDMR.py`, you need to prepare a sample information file named `in_cyber.lab`.
This file should be a **tab-separated** text file with **three columns**:
1. **Sample ID**
2. **Group name** (e.g., `normal`, `tumor`, `control`, `treatment`)
3. **Absolute path to the BED-format methylation file**

Example (`in_cyber.lab`):
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
```bash
python cyberDMR.py \
    --out_dir cyberDMR_result \
    --threads 4 \
    --group1 lethal \
    --group2 normal
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

```bash
bash cyberDMR.sh <indir> <outdir> <group1> <group2> <threads>
```

| Parameter | Required | Description | Example |
|-------------|----------|----------------------------------------------|----------------|
| `<indir>` | ✅ | Input folder with `.bed` files | `./input/` |
| `<outdir>` | ✅ | Output folder for results | `./results/` |
| `<group1>` | ✅ | First group label (used in filenames) | `lethal` |
| `<group2>` | ✅ | Second group label | `normal` |
| `<threads>` | ✅ | Number of threads for parallel computation | `8` |

View help information:

```bash
bash simulate_data.sh -h
bash cyberDMR.sh -h
```

Example:

```
bash simulate_data.sh -o ../data/simulate_data -t 100
bash cyberDMR.sh ./data/real_data/chr22 ./data/real_data/chr22/cyberDMR_result lethal normal 8
```

# Release Notes

## Release Notes – cyberDMR v1.0.0

**Release Date:** 2025-05-13
**Status:** Initial release

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
