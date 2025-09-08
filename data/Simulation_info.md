# Simulated Data Information

---

## 1. DMR density

| Add noise | Simulated ID | Parameters |
|-----------|--------------|------------|
| No  | samples_011 | `bash simulate_data.sh --density "dense" --seed 30 -o samples_011` |
|     | samples_012 | `bash simulate_data.sh --density "sparse" --seed 30 -o samples_012` |
|     | samples_013 | `bash simulate_data.sh --dense_ratio 0.3 --seed 30 -o samples_013` |
|     | samples_043 | `bash simulate_data.sh --dense_ratio 0.6 --seed 30 -o samples_043` |
|     | samples_044 | `bash simulate_data.sh --dense_ratio 0.9 --seed 30 -o samples_044` |
| Yes | samples_014 | `bash simulate_data.sh --density "dense" -u 0.01 -i 0.01 --seed 30 -o samples_014` |
|     | samples_015 | `bash simulate_data.sh --density "sparse" -u 0.01 -i 0.01 --seed 30 -o samples_015` |
|     | samples_016 | `bash simulate_data.sh --dense_ratio 0.3 -u 0.01 -i 0.01 --seed 30 -o samples_016` |
|     | samples_045 | `bash simulate_data.sh --dense_ratio 0.6 -u 0.01 -i 0.01 --seed 30 -o samples_045` |
|     | samples_046 | `bash simulate_data.sh --dense_ratio 0.9 -u 0.01 -i 0.01 --seed 30 -o samples_046` |

---

## 2. DMR length

| Add noise | Simulated ID | Parameters |
|-----------|--------------|------------|
| Yes | samples_041 | `bash simulate_data.sh --length_mean 100 --length_std 10 --seed 42 -o samples_041` |
|     | samples_003 | `bash simulate_data.sh --length_mean 200 --length_std 20 --seed 42 -o samples_003` |
|     | samples_004 | `bash simulate_data.sh --length_mean 500 --length_std 50 --seed 42 -o samples_004` |
|     | samples_005 | `bash simulate_data.sh --length_mean 1000 --length_std 100 --seed 42 -o samples_005` |
|     | samples_006 | `bash simulate_data.sh --length_mean 2000 --length_std 200 --seed 42 -o samples_006` |
| No  | samples_042 | `bash simulate_data.sh --length_mean 100 --length_std 10 -u 0.01 -i 0.01 --seed 42 -o samples_042` |
|     | samples_007 | `bash simulate_data.sh --length_mean 200 --length_std 20 -u 0.01 -i 0.01 --seed 42 -o samples_007` |
|     | samples_008 | `bash simulate_data.sh --length_mean 500 --length_std 50 -u 0.01 -i 0.01 --seed 42 -o samples_008` |
|     | samples_009 | `bash simulate_data.sh --length_mean 1000 --length_std 100 -u 0.01 -i 0.01 --seed 42 -o samples_009` |
|     | samples_010 | `bash simulate_data.sh --length_mean 2000 --length_std 200 -u 0.01 -i 0.01 --seed 42 -o samples_010` |

---

## 3. Methylation differences

| Add noise | Simulated ID | Parameters |
|-----------|--------------|------------|
| No  | samples_047 | `bash simulate_data.sh --mean_delta 0.11 --seed 21 -o samples_047` |
|     | samples_021 | `bash simulate_data.sh --mean_delta 0.15 --seed 21 -o samples_021` |
|     | samples_022 | `bash simulate_data.sh --mean_delta 0.25 --seed 21 -o samples_022` |
|     | samples_023 | `bash simulate_data.sh --mean_delta 0.35 --seed 21 -o samples_023` |
|     | samples_024 | `bash simulate_data.sh --mean_delta 0.45 --seed 21 -o samples_024` |
| Yes | samples_048 | `bash simulate_data.sh --mean_delta 0.11 -u 0.01 -i 0.01 --seed 21 -o samples_048` |
|     | samples_017 | `bash simulate_data.sh --mean_delta 0.15 -u 0.01 -i 0.01 --seed 21 -o samples_017` |
|     | samples_018 | `bash simulate_data.sh --mean_delta 0.25 -u 0.01 -i 0.01 --seed 21 -o samples_018` |
|     | samples_019 | `bash simulate_data.sh --mean_delta 0.35 -u 0.01 -i 0.01 --seed 21 -o samples_019` |
|     | samples_020 | `bash simulate_data.sh --mean_delta 0.45 -u 0.01 -i 0.01 --seed 21 -o samples_020` |

---

## 4. CpG coverage

| Add noise | Simulated ID | Parameters |
|-----------|--------------|------------|
| No  | samples_049 | `bash simulate_data.sh --coverage_mean 5 --coverage_std 3 --seed 57 -o samples_049` |
|     | samples_025 | `bash simulate_data.sh --coverage_mean 10 --coverage_std 5 --seed 57 -o samples_025` |
|     | samples_026 | `bash simulate_data.sh --coverage_mean 25 --coverage_std 5 --seed 57 -o samples_026` |
|     | samples_027 | `bash simulate_data.sh --coverage_mean 40 --coverage_std 5 --seed 57 -o samples_027` |
|     | samples_050 | `bash simulate_data.sh --coverage_mean 80 --coverage_std 5 --seed 57 -o samples_050` |
| Yes | samples_053 | `bash simulate_data.sh --coverage_mean 5 --coverage_std 3 -u 0.01 -i 0.01 --seed 57 -o samples_053` |
|     | samples_028 | `bash simulate_data.sh --coverage_mean 10 --coverage_std 5 -u 0.01 -i 0.01 --seed 57 -o samples_028` |
|     | samples_029 | `bash simulate_data.sh --coverage_mean 25 --coverage_std 5 -u 0.01 -i 0.01 --seed 57 -o samples_029` |
|     | samples_030 | `bash simulate_data.sh --coverage_mean 40 --coverage_std 5 -u 0.01 -i 0.01 --seed 57 -o samples_030` |
|     | samples_052 | `bash simulate_data.sh --coverage_mean 80 --coverage_std 5 -u 0.01 -i 0.01 --seed 57 -o samples_052` |

---

## 5. Sample size

| Add noise | Simulated ID | Parameters |
|-----------|--------------|------------|
| No  | samples_001 | `bash simulate_data.sh --n_control 1 --n_treatment 1 --seed 101 -o samples_001` |
|     | samples_031 | `bash simulate_data.sh --n_control 3 --n_treatment 3 --seed 101 -o samples_031` |
|     | samples_032 | `bash simulate_data.sh --n_control 5 --n_treatment 5 --seed 101 -o samples_032` |
|     | samples_033 | `bash simulate_data.sh --n_control 10 --n_treatment 10 --seed 101 -o samples_033` |
|     | samples_034 | `bash simulate_data.sh --n_control 25 --n_treatment 25 --seed 101 -o samples_034` |
|     | samples_039 | `bash simulate_data.sh --n_control 50 --n_treatment 50 --seed 101 -o samples_039` |
| Yes | samples_002 | `bash simulate_data.sh --n_control 1 --n_treatment 1 -u 0.01 -i 0.01 --seed 101 -o samples_002` |
|     | samples_035 | `bash simulate_data.sh --n_control 3 --n_treatment 3 -u 0.01 -i 0.01 --seed 101 -o samples_035` |
|     | samples_036 | `bash simulate_data.sh --n_control 5 --n_treatment 5 -u 0.01 -i 0.01 --seed 101 -o samples_036` |
|     | samples_037 | `bash simulate_data.sh --n_control 10 --n_treatment 10 -u 0.01 -i 0.01 --seed 101 -o samples_037` |
|     | samples_038 | `bash simulate_data.sh --n_control 25 --n_treatment 25 -u 0.01 -i 0.01 --seed 101 -o samples_038` |
|     | samples_040 | `bash simulate_data.sh --n_control 50 --n_treatment 50 -u 0.01 -i 0.01 --seed 101 -o samples_040` |

---

## 6. Metilene_data

| Simulated ID | Parameters |
|--------------|------------|
| metilene_1551   | `Rscript simulate_DMRs_WGBS.R 15 5 1 background/ metilene_1551` |
| metilene_155087 | `Rscript simulate_DMRs_WGBS.R 15 5 0.87 background/ metilene_155087` |
| metilene_15506  | `Rscript simulate_DMRs_WGBS.R 15 5 0.6 background/ metilene_15506` |
| metilene_4031   | `Rscript simulate_DMRs_WGBS.R 40 3 1 background/ metilene_4031` |
| metilene_403087 | `Rscript simulate_DMRs_WGBS.R 40 3 0.87 background/ metilene_403087` |
| metilene_40306  | `Rscript simulate_DMRs_WGBS.R 40 3 0.6 background/ metilene_40306` |

---
