import os
import pandas as pd
import numpy as np
from numba import njit
from sklearn.impute import KNNImputer
from concurrent.futures import ThreadPoolExecutor
from typing import List, Dict

# 预生成80%固定的位点（假设全基因组有100万个碱基）
BASE_CPG_SITES = np.sort(np.random.randint(1, 1_000_000, size=800))  # 80%固定位点

def generate_methylation_data(chr_name="chr1", num_points=1000, coverage_range=(1, 20), meth_range=(0, 1)):
    """
    生成模拟的 CpG 甲基化数据，保证 80% 位点固定 20% 位点随机变动。
    
    参数：
    - chr_name: 染色体名称
    - num_points: 总共要生成的CpG位点数量
    - coverage_range: (min_coverage, max_coverage) 覆盖度的范围
    - meth_range: (min_meth, max_meth) 甲基化水平的范围

    返回：
    - pandas DataFrame 包含 Chr, Pos, Meth_Level, Coverage
    """
    np.random.seed()  # 确保每次调用产生不同随机数据

    # 计算 80% 和 20% 的位点数量
    num_fixed = int(num_points * 0.8)
    num_random = num_points - num_fixed

    # 1. **80% 固定位点**（从预生成的 BASE_CPG_SITES 里选取）
    fixed_positions = np.random.choice(BASE_CPG_SITES, num_fixed, replace=False)

    # 2. **20% 随机生成的位点**
    random_positions = np.random.randint(1, 1_000_000, size=num_random)

    # 3. **合并并排序**
    all_positions = np.sort(np.concatenate([fixed_positions, random_positions]))

    # 4. **生成随机甲基化水平**
    meth_levels = np.random.uniform(meth_range[0], meth_range[1], size=num_points)

    # 5. **生成覆盖度（Gamma 分布模拟）**
    shape, scale = 2, 5  # Gamma 分布参数，右偏分布模拟实际测序深度
    coverage = np.random.gamma(shape, scale, size=num_points).astype(int)
    coverage = np.clip(coverage, coverage_range[0], coverage_range[1])  # 限制覆盖度范围

    # 6. **组装 DataFrame**
    df = pd.DataFrame({
        "Chr": [chr_name] * num_points,
        "Pos": all_positions,
        "Meth_Level": meth_levels,
        "Coverage": coverage
    })
    
    return df

def simulate_multiple_samples(output_folder, group1_count=5, group2_count=5, chr_name="chr1", num_points=1000):
    """
    生成多个样本的 CpG 甲基化数据，并按 group1/group2 组织，保存为文件。

    参数：
    - group1_count: group1 的样本数
    - group2_count: group2 的样本数
    - chr_name: 染色体名
    - num_points: 每个样本的数据点数

    返回：
    - List[Dict] 格式，每个样本包含 sample, group, data
    """
    samples = []
    #output_folder = "/home/ly/shell/deepDMR/data/simulate_sample_data/raw"
    os.makedirs(output_folder, exist_ok=True)

    sample_index = 1 # 全局编号
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

    # 写出 inlab.txt 文件
    inlab_path = os.path.join(output_folder, "inlab.txt")
    with open(inlab_path, 'w') as f:
        f.write("\n".join(inlab_records))

    print(f"[✔] inlab.txt 已写出到: {inlab_path}")

    return samples

# knn + 加权方式
def fill_low_coverage_cpg_knn(df, coverage_threshold=5, n_neighbors=3):
    """
    对低覆盖 CpG 位点进行 KNN 插值 + 加权融合。

    参数：
        - df: 包含 "Chr", "Pos", "Meth_Level", "Coverage" 的 DataFrame
        - coverage_threshold: 低覆盖度定义的阈值
        - n_neighbors: KNN 的近邻数

    返回：
        - 新增 "Final_Meth_Level" 的 DataFrame
    """
    if df.empty or df["Meth_Level"].dropna().empty:
        # 空表或Meth level 全为空
        df = df.copy()
        df["Final_Meth_Level"] = pd.Series(dtype=float)
        df["Weight"] = pd.Series(dtype=int)
        return df

    df = df.copy()
    print(df)
    # 确保数值列为 float 类型
    df["Coverage"] = pd.to_numeric(df["Coverage"], errors="coerce")
    df["Meth_Level"] = pd.to_numeric(df["Meth_Level"], errors="coerce")

    # 生成插值版本
    df_knn = df.copy()
    low_cov_mask = df["Coverage"] < coverage_threshold
    print(low_cov_mask)
    #df_knn.loc[low_cov_mask, "Meth_Level"] = np.nan # 设置低覆盖为缺失
    print(df_knn)
    # KNN 插补
    imputer = KNNImputer(n_neighbors=n_neighbors)
    df["Meth_Level_KNN"] = imputer.fit_transform(df_knn[["Meth_Level"]])

    # 融合插补值与原始值
    df["Weight"] = np.clip(df["Coverage"] / coverage_threshold, 0, 1)
    df["Final_Meth_Level"] = (
        df["Weight"] * df["Meth_Level"] + (1 - df["Weight"]) * df["Meth_Level_KNN"]
    )

    # 限制在 0~1 范围
    df["Final_Meth_Level"] = df["Final_Meth_Level"].clip(0, 1)
    print(df)
    return df

# 采用距离 + 加权方式
def fill_low_coverage_cpg(df, coverage_threshold=5, max_distance=500):
    """
    使用距离加权平均 + 当前值加权融合，对低覆盖 CpG 位点进行平滑填补。

    参数：
    - df: DataFrame，包含 ['Chr', 'Pos', 'Meth_Level', 'Coverage']
    - coverage_threshold: 低覆盖度的定义阈值
    - max_distance: 限制最大邻居距离

    返回：
    - df: 添加 'Final_Meth_Level' 列，融合当前值与邻居预测值
    """
    df = df.copy()
    df = df.sort_values("Pos").reset_index(drop=True)

    # 类型转换
    df["Coverage"] = pd.to_numeric(df["Coverage"], errors="coerce")
    df["Meth_Level"] = pd.to_numeric(df["Meth_Level"], errors="coerce")
    df["Pos"] = pd.to_numeric(df["Pos"], errors="coerce")
    df["Final_Meth_Level"] = df["Meth_Level"]

    rows_to_drop = []

    for i in range(len(df)):
        coverage = df.at[i, "Coverage"]
        beta_i = df.at[i, "Meth_Level"]

        if pd.isna(beta_i):
            continue  # 真缺失，跳过（你可以加专门策略处理）

        if coverage >= coverage_threshold:
            continue  # 高覆盖保留原始值

        pos_i = df.at[i, "Pos"]

        # 找前后邻居
        prev = next((j for j in range(i - 1, -1, -1)
                     if pd.notna(df.at[j, "Meth_Level"])), None)
        next_ = next((j for j in range(i + 1, len(df))
                      if pd.notna(df.at[j, "Meth_Level"])), None)

        if prev is None or next_ is None:
            rows_to_drop.append(i)
            continue  # 缺一边，不插值

        # 距离计算
        d1 = abs(pos_i - df.at[prev, "Pos"])
        d2 = abs(df.at[next_, "Pos"] - pos_i)

        if d1 > max_distance or d2 > max_distance:
            rows_to_drop.append(i)  # 距离大后续没有可能加入到计算当中
            continue  # 超距离

        # 邻居值
        beta1 = df.at[prev, "Meth_Level"]
        beta2 = df.at[next_, "Meth_Level"]
        beta_hat = (beta1 / d1 + beta2 / d2) / (1/d1 + 1/d2)

        # 融合当前值
        w = min(1.0, coverage / coverage_threshold)
        beta_final = w * beta_i + (1 - w) * beta_hat

        df.at[i, "Final_Meth_Level"] = np.clip(beta_final, 0, 1)

    # 删除无法插补的行
    df = df.drop(index=rows_to_drop).reset_index(drop=True)
    return df

# 提速
def fill_low_coverage_cpg_fast(df, coverage_threshold=5, max_distance=500):
    """
    对无缺失值的 CpG 数据，使用 GIMMEcpg 邻居加权插值（极致优化版）

    参数：
    - df: 包含 ['Chr', 'Pos', 'Meth_Level', 'Coverage']
    - coverage_threshold: Coverage 阈值
    - max_distance: 最大插值范围

    返回：
    - 添加 'Final_Meth_Level' 的 DataFrame，低覆盖无法插值的行将被删除
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
            continue  # 高覆盖不处理

        pos_i = pos_arr[i]

        # 找前一个点（i-1）
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

        # 插值
        beta1 = meth_arr[prev]
        beta2 = meth_arr[next_]

        # 避免除以 0
        if d1 == 0 and d2 == 0:
            rows_to_drop.append(i)
            continue
        elif d1 == 0:
            beta_hat = beta1
        elif d2 == 0:
            beta_hat = beta2
        else:
            beta_hat = (beta1 / d1 + beta2 / d2) / (1 / d1 + 1 / d2)

        # 融合当前值
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
            continue  # 高覆盖不处理

        pos_i = pos_arr[i]

        # 相邻两个点
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

        # 距离为 0 时的容错处理
        if d1 == 0 and d2 == 0:
            drop_mask[i] = True
            continue
        elif d1 == 0:
            beta_hat = beta1
        elif d2 == 0:
            beta_hat = beta2
        else:
            beta_hat = (beta1 / d1 + beta2 / d2) / (1 / d1 + 1 / d2)

        # 权重融合原始值
        w = min(1.0, cov / threshold)
        beta_final = w * beta_i + (1 - w) * beta_hat
        final_arr[i] = min(max(beta_final, 0), 1)

    return final_arr, drop_mask

def fill_low_coverage_cpg_numba(df, coverage_threshold=5, max_distance=500):
    """
    使用 Numba 加速的 GIMMEcpg 插补算法（无缺失假设），适合大规模数据

    参数：
    - df: 包含 'Pos', 'Meth_Level', 'Coverage' 的 DataFrame
    - coverage_threshold: 低覆盖定义阈值
    - max_distance: 可接受的邻居最大距离

    返回：
    - 插补后带 'Final_Meth_Level' 的 DataFrame，删除无法插补的低覆盖点
    """
    df = df.copy()
    df = df.sort_values("Pos").reset_index(drop=True)

    # 转为 numpy（支持 numba）
    pos_arr = pd.to_numeric(df["Pos"], errors="coerce").values.astype(np.float64)
    meth_arr = pd.to_numeric(df["Meth_Level"], errors="coerce").values.astype(np.float64)
    cov_arr = pd.to_numeric(df["Coverage"], errors="coerce").values.astype(np.float64)

    final_arr, drop_mask = _fill_low_cov_numba(pos_arr, meth_arr, cov_arr,
                                                coverage_threshold, max_distance)

    df["Final_Meth_Level"] = final_arr
    df = df.loc[~drop_mask].reset_index(drop=True)
    return df

# 单样本串行处理
def process_samples_old(sample_data, coverage_threshold=5, max_distance=500):
    """
    顺序处理每个样本数据，执行 KNN 平滑处理。

    参数：
    - sample_data: List[Dict]，每个 dict 包含 'sample', 'group', 'data'
    - coverage_threshold: 低覆盖度阈值
    - n_neighbors: KNN 的近邻个数

    返回：
    - List[Dict]，结构与原始 sample_data 相同，但 data 为处理后的 DataFrame
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
        #sample_df = processed_samples[sample_id][["Pos", "Final_Meth_Level"]].rename(columns={"Final_Meth_Level": sample_id})
        #base_df = base_df.merge(sample_df, on="Pos", how="left")  # `left join`，确保 `Pos` 统一
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
    # 示例调用
    #simulated_data = generate_methylation_data(chr_name="chr2", num_points=1000)
    #process_data = fill_low_coverage_cpg(simulated_data)

    #print(simulated_data.head(10))
    #print(process_data.head(10))

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


