# cyberDMR

**cyberDMR** is a robust and high-resolution framework for detecting differentially methylated regions (DMRs) from whole-genome bisulfite sequencing (WGBS) data.

### Features
- Base-level smoothing for low-coverage CpGs
- Seed-guided clustering for consistent CpG grouping
- Weighted beta regression with LRT for statistical inference
- Supports both simulated and real datasets
- Benchmark-ready and highly scalable

### Installation
```bash
git clone https://github.com/YLeeHIT/cyberDMR.git
cd cyberDMR
pip install -r requirements.txt
