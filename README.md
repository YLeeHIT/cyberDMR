# cyberDMR

![Version](https://img.shields.io/badge/version-1.0.0-blue)
![Language](https://img.shields.io/badge/language-python-blue)
![Language](https://img.shields.io/badge/language-shell-4EAA25)
![Language](https://img.shields.io/badge/language-R-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Platform](https://img.shields.io/badge/platform-linux%20|%20macOS-brightgreen)

**cyberDMR** is a robust and high-resolution framework for detecting differentially methylated regions (DMRs) from whole-genome bisulfite sequencing (WGBS) data.

### Introduction
Differentially methylated regions (DMRs) are key genomic features reflecting changes in DNA methylation status. Accurate identification of DMRs is crucial for investigating tissue-specific regulation, disease mechanisms, and population-level epigenetic variation. However, existing tools encounter difficulties when applied to increasingly complex WGS datasets.

### Features
- Base-level smoothing for low-coverage CpGs
- Seed-guided clustering for consistent CpG grouping
- Weighted beta regression with LRT for statistical inference
- Supports both simulated and real datasets
- Benchmark-ready and highly scalable

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
After generating the simulated methylation data, you can run **cyberDMR** to detect DMRs between two groups:

```bash
python cyberDMR.py \
    --out_dir /mnt/data/sample_outputs/results \
    --threads 4 \
    --group1 control \
    --group2 treatment
```

### Quick Start (Recommended)

To simplify everything, you can run the pre-configured shell script:

```bash
bash simulate_data.sh -o ../data/simulate_data
```

View help information:
```bash
bash simulate_data.sh -h
```

# Release Notes

## Release Notes – cyberDMR v1.0.0

**Release Date:** 2025-05-13
**Status:** Initial release

If you use cyberDMR in your research, please cite:
@article{li2025cyberdmr,
    title={cyberDMR: a robust and high-sensitivity approach for differentially methylated regions detection},
    author={Li, Yang and others},
    journal={Bioinformatics},
    year={2025},
    note={Manuscript under review}
}

# Contributors

This package is developed and maintaned by [Lee](https://github.com/YLeeHIT) and [Chen](https://github.com/chong-hun). If you want to contribute, please leave an issue or submit a pull request. Thank you.

# License
