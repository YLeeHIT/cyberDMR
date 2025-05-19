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
git clone https://github.com/YLeeHIT/cyberDMR.git
cd cyberDMR
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
### input requirement for cyberDMR

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

To build an inlab file, you can refer to the following instructions:
```
cd ./cyberDMR/data/real_data/chr22
ls noh_lethal_*bed noh_normal_*bed > raw.lab
awk -v dir=$(pwd) '{split($1,y,"_");print y[3]"\t"y[2]"\t"dir"/"$1}' raw.lab > in_cyber.lab
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

If you want to merge all chromosome files, you can refer to the following instructions:
```
cat ./chr*txt |sort -k1,1V -k2,2n -k3,3n > final_result.txt 
```

### Quick Start (Recommended)

To simplify everything, you can run the pre-configured shell script:

```bash
bash ./simulate_data.sh -o ../data/simulate_data -t 100
```

View help information:
```bash
bash ./simulate_data.sh -h
```

# Release Notes

## Release Notes – cyberDMR v1.0.0

**Release Date:** 2025-05-13
**Status:** Initial release

If you use cyberDMR in your research, please cite:

```bibtex
@article{
    title={cyberDMR: a robust and high-sensitivity approach for differentially methylated regions detection},
    author={Li, Yang and others},
    journal={Bioinformatics},
    year={2025},
    note={Manuscript under review}
}
```

# Contributors

This package is developed and maintaned by [Lee](https://github.com/YLeeHIT) and [Chen](https://github.com/chong-hun). If you want to contribute, please leave an issue or submit a pull request. Thank you.

# License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
