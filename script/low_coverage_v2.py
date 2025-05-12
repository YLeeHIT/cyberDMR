import os
import pandas as pd
import numpy as np
from numba import njit
from sklearn.impute import KNNImputer
from concurrent.futures import ThreadPoolExecutor
from typing import List, Dict

# Pre-generate 80% of sites as fixed positions (assuming the entire genome consists of 1 million bases)
BASE_CPG_SITES = np.sort(np.random.randint(1, 1_000_000, size=800))  # 80%固定位点

def generate_methylation_data(chr_name="chr1", num_points=1000, coverage_range=(1, 20), meth_range=(0, 1)):
    """
    Generate simulated CpG methylation data, ensuring 80% of sites are fixed and 20% of sites vary randomly
    
    Parameters：
    - chr_name: Chromosome name
    - num_points:  Total number of CpG sites to generate
    - coverage_range: Tuple specifying the range for coverage, as (min_coverage, max_coverage)
    - meth_range: Tuple specifying the range for methylation level, as (min_meth, max_meth)

    Returns：
    - pandas DataFrame:A DataFrame with columns: 'Chr', 'Pos', 'Meth_Level', 'Coverage'
    """
    np.random.seed()  # Ensure each call produces different random data

    # Calculate the number of sites for the 80% and 20% portions
    num_fixed = int(num_points * 0.8)
    num_random = num_points - num_fixed

    # 1. **80% fixed sites**（selected from the pre-generated BASE_CPG_SITES）
    fixed_positions = np.random.choice(BASE_CPG_SITES, num_fixed, replace=False)

    # 2. **20% randomly generated sites**
    random_positions = np.random.randint(1, 1_000_000, size=num_random)

    # 3. **Merge and sort**
    all_positions = np.sort(np.concatenate([fixed_positions, random_positions]))

    # 4. **Generate random methylation levels**
    meth_levels = np.random.uniform(meth_range[0], meth_range[1], size=num_points)

    # 5. **Generate coverage (simulated using a Gamma distribution)**
    shape, scale = 2, 5  # Gamma distribution parameters; a right-skewed distribution to simulate actual sequencing depth
    coverage = np.random.gamma(shape, scale, size=num_points).astype(int)
    coverage = np.clip(coverage, coverage_range[0], coverage_range[1])  # Limit the coverage range

    # 6. **Assemble DataFrame**
    df = pd.DataFrame({
        "Chr": [chr_name] * num_points,
        "Pos": all_positions,
        "Meth_Level": meth_levels,
        "Coverage": coverage
    })
    
    return df

def simulate_multiple_samples(output_folder, group1_count=5, group2_count=5, chr_name="chr1", num_points=1000):
    """
    Generates CpG methylation data for multiple samples, organizes them by group1 and group2, and saves the data to a file

    参数：
    - group1_count: The number of samples in group1
    - group2_count: The number of samples in group2
    - chr_name: The chromosome name
    - num_points: The number of data points to generate per sample

    返回：
    - List[Dict] A list of dictionaries. Each dictionary represents a sample and contains keys: 'sample' (sample identifier), 'group' (group identifier), and 'data' (the corresponding CpG methylation data)
    """
    samples = []
    os.makedirs(output_folder, exist_ok=True)

    sample_index = 1 # Global ID
    inlab_records = []

    for i in range(group1_count):
        sample_name = f"Sample_{sample_index}"
        group = "group1"
        df = generate_methylation_data(chr_name, num_points)
        sample_path = os.path.join(output_folder, f"{sample_name}.csv")
        df.to_csv(sample_path, sep='\t', index=False)
        samples.append({
            'sample': sample_name,
            'group': group,
            'data': df
        })
        sample_index += 1
        inlab_records.append(f"{sample_name}\t{group}\t{sample_path}")

    for i in range(group2_count):
        sample_name = f"Sample_{sample_index}"
        group = "group2"
        df = generate_methylation_data(chr_name, num_points)
        sample_path = os.path.join(output_folder, f"{sample_name}.csv")
        df.to_csv(sample_path, sep='\t', index=False)
        samples.append({
            'sample': sample_name,
            'group': group,
            'data': df
        })
        sample_index += 1
        inlab_records.append(f"{sample_name}\t{group}\t{sample_path}")

    # Write out the inlab.txt file
    inlab_path = os.path.join(output_folder, "inlab.txt")
    with open(inlab_path, 'w') as f:
        f.write("\n".join(inlab_records))

    print(f"[✔] inlab.txt 已写出到: {inlab_path}")

    return samples

# KNN + weighting method
def fill_low_coverage_cpg_knn(df, coverage_threshold=5, n_neighbors=3):
    """
    Performs k-Nearest Neighbors (KNN) imputation combined with weighted fusion on low-coverage CpG sites

    Parameters：
        - df:  Input DataFrame containing the columns 'Chr', 'Pos', 'Meth_Level', and 'Coverage'
        - coverage_threshold: The threshold for defining low coverage
        - n_neighbors: The number of nearest neighbors for KNN

    Returns：
        - pandas.DataFrame: The input DataFrame with an added 'Final_Meth_Level' column, representing the imputed/fused methylation levels
    """
    if df.empty or df["Meth_Level"].dropna().empty:
        # Empty table or the 'Meth_Level' column is entirely empty/null
        df = df.copy()
        df["Final_Meth_Level"] = pd.Series(dtype=float)
        df["Weight"] = pd.Series(dtype=int)
        return df

    df = df.copy()
    print(df)
    # Ensure numerical columns are of float type
    df["Coverage"] = pd.to_numeric(df["Coverage"], errors="coerce")
    df["Meth_Level"] = pd.to_numeric(df["Meth_Level"], errors="coerce")

    # Generate an imputed version
    df_knn = df.copy()
    low_cov_mask = df["Coverage"] < coverage_threshold
    print(low_cov_mask)
    print(df_knn)
    # KNN imputation
    imputer = KNNImputer(n_neighbors=n_neighbors)
    df["Meth_Level_KNN"] = imputer.fit_transform(df_knn[["Meth_Level"]])

    # Merge imputed values and original values
    df["Weight"] = np.clip(df["Coverage"] / coverage_threshold, 0, 1)
    df["Final_Meth_Level"] = (
        df["Weight"] * df["Meth_Level"] + (1 - df["Weight"]) * df["Meth_Level_KNN"]
    )

    # Limit to the 0-1 range
    df["Final_Meth_Level"] = df["Final_Meth_Level"].clip(0, 1)
    print(df)
    return df

# Adopt a distance + weighting method
def fill_low_coverage_cpg(df, coverage_threshold=5, max_distance=500):
    """
    Smoothly imputes low-coverage CpG sites using a combination of a distance-weighted average (from neighbors) and weighted fusion with the site's current value

    Parameters:
    - df: DataFrame: Input DataFrame containing the columns 'Chr', 'Pos', 'Meth_Level', and 'Coverage'
    - coverage_threshold: The threshold used to define low coverage
    - max_distance: The maximum distance to consider for neighbors

    Returns:
    - pandas.DataFrame: The input DataFrame with an added 'Final_Meth_Level' column. This column's values are the result of fusing the original (current) methylation level with a predicted level derived from its neighbors
    """
    df = df.copy()
    df = df.sort_values("Pos").reset_index(drop=True)

    # Type conversion
    df["Coverage"] = pd.to_numeric(df["Coverage"], errors="coerce")
    df["Meth_Level"] = pd.to_numeric(df["Meth_Level"], errors="coerce")
    df["Pos"] = pd.to_numeric(df["Pos"], errors="coerce")
    df["Final_Meth_Level"] = df["Meth_Level"]

    rows_to_drop = []

    for i in range(len(df)):
        coverage = df.at[i, "Coverage"]
        beta_i = df.at[i, "Meth_Level"]

        if pd.isna(beta_i):
            continue  # Truly missing, skip (you can add a dedicated strategy for handling)

        if coverage >= coverage_threshold:
            continue  # For high coverage, retain the original value

        pos_i = df.at[i, "Pos"]

        # Find preceding and succeeding neighbors
        prev = next((j for j in range(i - 1, -1, -1)
                     if pd.notna(df.at[j, "Meth_Level"])), None)
        next_ = next((j for j in range(i + 1, len(df))
                      if pd.notna(df.at[j, "Meth_Level"])), None)

        if prev is None or next_ is None:
            rows_to_drop.append(i)
            continue  # If one side (of neighbors) is missing, do not impute

        # Distance calculation
        d1 = abs(pos_i - df.at[prev, "Pos"])
        d2 = abs(df.at[next_, "Pos"] - pos_i)

        if d1 > max_distance or d2 > max_distance:
            rows_to_drop.append(i)  # If the distance is large, it will not be included in subsequent calculations
            continue  # Exceeds distance limit

        # Neighbor value
        beta1 = df.at[prev, "Meth_Level"]
        beta2 = df.at[next_, "Meth_Level"]
        beta_hat = (beta1 / d1 + beta2 / d2) / (1/d1 + 1/d2)

        # Fuse with the current value
        w = min(1.0, coverage / coverage_threshold)
        beta_final = w * beta_i + (1 - w) * beta_hat

        df.at[i, "Final_Meth_Level"] = np.clip(beta_final, 0, 1)

    # Delete rows that cannot be imputed
    df = df.drop(index=rows_to_drop).reset_index(drop=True)
    return df

# Speed up
def fill_low_coverage_cpg_fast(df, coverage_threshold=5, max_distance=500):
    """
    For CpG data that has no missing values, this function applies GIMMEcpg neighbor weighted imputation (highly optimized version)

    Parameters:
    - df: Input DataFrame containing the columns 'Chr', 'Pos', 'Meth_Level', and 'Coverage'
    - coverage_threshold: The coverage threshold
    - max_distance: The maximum distance to consider for neighbors during imputation

    Returns:
    - The DataFrame with an added 'Final_Meth_Level' column. Rows with low coverage that could not be imputed will have been deleted
    """
    df = df.copy()
    df = df.sort_values("Pos").reset_index(drop=True)

    df["Coverage"] = pd.to_numeric(df["Coverage"], errors="coerce")
    df["Meth_Level"] = pd.to_numeric(df["Meth_Level"], errors="coerce")
    df["Pos"] = pd.to_numeric(df["Pos"], errors="coerce")
    df["Final_Meth_Level"] = df["Meth_Level"]

    pos_arr = df["Pos"].values
    meth_arr = df["Meth_Level"].values
    coverage_arr = df["Coverage"].values

    rows_to_drop = []

    for i in range(len(df)):
        coverage = coverage_arr[i]
        beta_i = meth_arr[i]

        if coverage >= coverage_threshold:
            continue  # High coverage: do not process

        pos_i = pos_arr[i]

        # Find the previous point (i-1)
        prev = i - 1 if i > 0 else -1
        next_ = i + 1 if i < len(df) - 1 else -1

        if prev == -1 or next_ == -1:
            rows_to_drop.append(i)
            continue

        d1 = abs(pos_i - pos_arr[prev])
        d2 = abs(pos_arr[next_] - pos_i)

        if d1 > max_distance or d2 > max_distance:
            rows_to_drop.append(i)
            continue

        # Interpolate
        beta1 = meth_arr[prev]
        beta2 = meth_arr[next_]

        # Avoid division by zero
        if d1 == 0 and d2 == 0:
            rows_to_drop.append(i)
            continue
        elif d1 == 0:
            beta_hat = beta1
        elif d2 == 0:
            beta_hat = beta2
        else:
            beta_hat = (beta1 / d1 + beta2 / d2) / (1 / d1 + 1 / d2)

        # Fuse with the current value
        w = min(1.0, coverage / coverage_threshold)
        beta_final = w * beta_i + (1 - w) * beta_hat
        df.at[i, "Final_Meth_Level"] = np.clip(beta_final, 0, 1)

    df = df.drop(index=rows_to_drop).reset_index(drop=True)
    return df

@njit
def _fill_low_cov_numba(pos_arr, meth_arr, coverage_arr, threshold, max_dist):
    n = len(pos_arr)
    final_arr = meth_arr.copy()
    drop_mask = np.zeros(n, dtype=np.bool_)

    for i in range(n):
        cov = coverage_arr[i]
        beta_i = meth_arr[i]

        if cov >= threshold:
            continue  # High coverage: do not process

        pos_i = pos_arr[i]

        # Two adjacent points
        prev = i - 1 if i > 0 else -1
        next_ = i + 1 if i < n - 1 else -1

        if prev == -1 or next_ == -1:
            drop_mask[i] = True
            continue

        d1 = abs(pos_i - pos_arr[prev])
        d2 = abs(pos_arr[next_] - pos_i)

        if d1 > max_dist or d2 > max_dist:
            drop_mask[i] = True
            continue

        beta1 = meth_arr[prev]
        beta2 = meth_arr[next_]

        # Error handling when distance is zero
        if d1 == 0 and d2 == 0:
            drop_mask[i] = True
            continue
        elif d1 == 0:
            beta_hat = beta1
        elif d2 == 0:
            beta_hat = beta2
        else:
            beta_hat = (beta1 / d1 + beta2 / d2) / (1 / d1 + 1 / d2)

        # Weighted fusion of original values
        w = min(1.0, cov / threshold)
        beta_final = w * beta_i + (1 - w) * beta_hat
        final_arr[i] = min(max(beta_final, 0), 1)

    return final_arr, drop_mask

def fill_low_coverage_cpg_numba(df, coverage_threshold=5, max_distance=500):
    """
    A GIMMEcpg imputation algorithm accelerated with Numba(operates under a 'no missing values' assumption for the input data), suitable for large-scale data

    Parameters:
    - df: Input DataFrame containing 'Pos', 'Meth_Level', and 'Coverage' columns
    - coverage_threshold: Threshold defining low coverage 
    - max_distance: Maximum acceptable distance for neighbors to be considered in imputation

    Returns:
    - pandas.DataFrame: The DataFrame after imputation, with an added 'Final_Meth_Level' column. Low-coverage points that could not be imputed are deleted
    """
    df = df.copy()
    df = df.sort_values("Pos").reset_index(drop=True)

    # Convert to NumPy (Numba-compatible).
    pos_arr = pd.to_numeric(df["Pos"], errors="coerce").values.astype(np.float64)
    meth_arr = pd.to_numeric(df["Meth_Level"], errors="coerce").values.astype(np.float64)
    cov_arr = pd.to_numeric(df["Coverage"], errors="coerce").values.astype(np.float64)

    final_arr, drop_mask = _fill_low_cov_numba(pos_arr, meth_arr, cov_arr,
                                                coverage_threshold, max_distance)

    df["Final_Meth_Level"] = final_arr
    df = df.loc[~drop_mask].reset_index(drop=True)
    return df

# Single-sample serial processing
def process_samples_old(sample_data, coverage_threshold=5, max_distance=500):
    """
    Sequentially processes the data for each sample by performing K-Nearest Neighbors (KNN) smoothing

    Parameters:
    - sample_data: List[Dict],A list of dictionaries, where each dictionary contains keys 'sample', 'group', and 'data'. The 'data' key typically holds the sample's measurement data
    - coverage_threshold: The threshold for defining low coverage
    - n_neighbors: The number of nearest neighbors to use for KNN

    Returns:
    - List[Dict]:A list of dictionaries with the same structure as the original sample_data, but where the value associated with the 'data' key is the DataFrame after processing (KNN smoothing)
    """
    processed_results = []

    for entry in sample_data:
        sample = entry['sample']
        group = entry['group']
        df = entry['data']

        try:
            processed_df = fill_low_coverage_cpg_numba(df, coverage_threshold, max_distance)
            print(f"{sample} belongs to {group} has finished")
            processed_results.append({
                'sample': sample,
                'group': group,
                'data': processed_df
            })
        except Exception as e:
            print(f"[] 样本 {sample} 处理失败: {e}")
            # 如果失败可以返回空 DataFrame
            processed_results.append({
                'sample': sample,
                'group': group,
                'data': pd.DataFrame()
            })

    return processed_results

def process_samples(sample_data, coverage_threshold=5, max_distance=500):
    """
    串行处理每个样本的甲基化数据，使用 GIMMEcpg-style 插补 + Numba 加速版本。

    参数：
    - sample_data: List[Dict]，每个 dict 包含 'sample', 'group', 'data'
    - coverage_threshold: 低覆盖度定义阈值
    - max_distance: 插补时左右邻居最大允许距离

    返回：
    - List[Dict]，结构与 sample_data 相同，但 'data' 为插补后的 DataFrame
    """
    processed_results = []

    for entry in sample_data:
        sample = entry['sample']
        group = entry['group']
        df = entry['data']

        try:
            processed_df = fill_low_coverage_cpg_numba(df, coverage_threshold, max_distance)
            print(f"Success {sample} (group: {group}) Interpolation completed")
            processed_results.append({
                'sample': sample,
                'group': group,
                'data': processed_df
            })
        except Exception as e:
            print(f"Failure {sample} : {e}")
            processed_results.append({
                'sample': sample,
                'group': group,
                'data': pd.DataFrame()
            })

    return processed_results

def process_samples_parallel(sample_data, coverage_threshold=5, max_distance=500, num_threads=4):
    """
    多线程并发处理多个样本数据（适合样本数量多，独立性强的情况）

    参数：
    - sample_data: List[Dict]，包含 sample/group/data
    - coverage_threshold: Coverage 阈值
    - max_distance: 插补时邻居最大距离
    - num_threads: 并发线程数（默认 4）

    返回：
    - List[Dict]，顺序与原始 sample_data 保持一致
    """
    def process_one(entry):
        sample = entry['sample']
        group = entry['group']
        df = entry['data']
        try:
            processed_df = fill_low_coverage_cpg_numba(df, coverage_threshold, max_distance)
            print(f"✅ {sample} 插补完成")
            return {'sample': sample, 'group': group, 'data': processed_df}
        except Exception as e:
            print(f"❌ {sample} 处理失败: {e}")
            return {'sample': sample, 'group': group, 'data': pd.DataFrame()}

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(process_one, sample_data))  # 顺序保留

    return results

def merge_samples_fast_old(processed_samples, sample_order):
    """
    采用更快的方法合并所有样本的数据，并确保 `Pos` 对齐，不同样本没有该 `Pos` 时填充 `NaN`。
    
    参数：
    - processed_samples: 经过 fill_low_coverage_cpg 处理后的样本数据（字典形式）
    - sample_order: 样本名称顺序

    返回：
    - 处理后的 DataFrame 缺失值为 `NaN`
    """
    # **1 提取所有 `Pos` 位置**
    all_positions = set()
    for df in processed_samples.values():
        all_positions.update(df["Pos"])

    # **2 统一 `Pos`**
    all_positions = sorted(all_positions)
    base_df = pd.DataFrame({"Pos": all_positions})

    # **3 用 `merge()` 按 `Pos` 对齐，缺失值填充 `NaN`**
    for sample_id in sample_order:
        sample_df = processed_samples[sample_id][["Pos", "Final_Meth_Level"]].copy()
        sample_df["Final_Meth_Level"] = sample_df["Final_Meth_Level"].round(4)  # **保留 4 位小数**
        sample_df.rename(columns={"Final_Meth_Level": sample_id}, inplace=True)
        base_df = base_df.merge(sample_df, on="Pos", how="left")  # `left join`，确保 `Pos` 统一

    base_df.insert(0, "Chr", processed_samples[sample_order[0]]["Chr"].iloc[0])  # 插入 `Chr` 列

    return base_df

def merge_samples_fast(processed_samples: dict) -> pd.DataFrame:
    """
    合并所有样本的 Final_Meth_Level 按 Pos 对齐 没有的位点填充 NaN。

    参数：
        - processed_samples: 字典，键为样本名，值为含有 "Pos" 和 "Final_Meth_Level" 的 DataFrame

    返回：
        - 合并后的 DataFrame 包含 Chr, Pos, 和每个样本的 Final_Meth_Level
    """
    if not processed_samples:
        raise ValueError("processed_samples 为空")

    all_positions = sorted(
        set(int(pos) for entry in processed_samples for pos in entry['data']["Pos"])
    )
    base_df = pd.DataFrame({"Pos": all_positions})

    # 按样本依次合并
    for entry in processed_samples:
        sample_id = entry['sample']
        df = entry['data'][["Pos", "Final_Meth_Level"]].copy()
        df["Pos"] = df["Pos"].astype(int)
        df["Final_Meth_Level"] = df["Final_Meth_Level"].round(4)
        df.rename(columns={"Final_Meth_Level": sample_id}, inplace=True)
        base_df = base_df.merge(df, on="Pos", how="left")

    # 插入 Chr 列（假设所有样本染色体相同）
    chr_name = processed_samples[0]['data']["Chr"].iloc[0]
    base_df.insert(0, "Chr", chr_name)

    return base_df




if __name__ == "__main__":
    
    output_folder = "/home/ly/shell/deepDMR/data/simulate_sample_data/raw2"
    os.makedirs(output_folder, exist_ok=True)  # 如果文件夹不存在，则创建
    
    num_sample = 10
    chr_name = "chr1"
    num_points = 1000
    samples = simulate_multiple_samples(output_folder=output_folder,
                                        group1_count=5, group2_count=5)


    num_threads = 4
    processed_samples = process_samples(samples, coverage_threshold=5, max_distance=500)

    print(processed_samples[0]['data'][['Chr', 'Pos', 'Meth_Level', 'Coverage', 'Final_Meth_Level']].head())
    merged_df = merge_samples_fast(processed_samples)
    print(merged_df.head())

    merged_df.to_csv(os.path.join(output_folder,"merge.txt"), sep='\t',index=False)


