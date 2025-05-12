import numpy as np
import pandas as pd
import random
import argparse
import shutil
import os
from scipy.stats import beta
from datetime import datetime
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Simulated DMR Data")
    parser.add_argument(
        "--output_dir",
        type=str,
        default="/home/ly/shell/deepDMR/data/simulate_from_python/samples_16",
        help="outdir (default: /home/ly/shell/deepDMR/data/simulate_from_python/samples_15"
    )
    return parser.parse_args()


# 模拟DMR区域的函数
def simulate_dmr_region_with_input_limit(chr_name, region_start, region_end, delta_methylation, max_cpgs=50, density="auto"):
    """
    模拟一个差异甲基化区域 DMR

    参数说明：
    - chr_name: 染色体名称 如 "chr1"
    - region_start: 区域起始位置 整数
    - region_end: 区域终止位置 必须大于起始位置 + 100
    - delta_methylation: 实验组与对照组的平均甲基化水平差值（可正可负）
    - max_cpgs: 区域中最多CpG数量 超过则提前终止, 默认50

    返回：
    一个字典 包括区域信息、CpG位点、两组甲基化值及平均差值
    """

    # 定义每两个CpG之间的最小和最大间隔（单位：bp）
    if density == "dense":
        min_distance = 2
        max_distance = 50
    elif density == "sparse":
        min_distance = 50
        max_distance = 499
    else:
        min_distance = 2
        max_distance = 499

    # 区域长度校验
    if region_end - region_start < 100:
        #raise ValueError("区域长度必须大于100bp")
        return None

    # 生成CpG位点位置
    cpg_sites = []
    pos = region_start
    while pos < region_end:
        if len(cpg_sites) >= max_cpgs:
            break  # 超过最大数量限制，提前终止

        if len(cpg_sites) > 0:
            max_gap = min(max_distance, region_end - pos)
            if max_gap < min_distance:
                break  # 剩余距离太小，无法添加新的CpG
            gap = random.randint(min_distance, max_gap)
        else:
            gap = 0  # 第一个CpG使用起始位置
        pos += gap
        if pos < region_end:
            cpg_sites.append(pos)

    if len(cpg_sites) == 0:
        #raise ValueError("无法在给定区域内放置CpG位点，请检查起始和终止位置。")
        return None

    # 如果因超过最大数量提前终止，则更新终止位置
    region_end_final = min(region_end, cpg_sites[-1])
    region_length = region_end_final - region_start

    # 判断是否满足三个前置条件
    if len(cpg_sites) < 5 or region_length < 100 or abs(delta_methylation) < 0.1:
        return None

    # 设置甲基化水平的均值（使用 beta 分布模拟）
    mean_control = np.random.uniform(0.2, 0.8)
    delta_methylation_modify = round(np.random.normal(loc=delta_methylation, scale=0.01), 3)
    mean_treatment = mean_control + delta_methylation_modify

    # 防止mean_treatment 超过阈值
    mean_treatment = min(max(mean_treatment, 0.001), 0.999)

    # 设置beta分布的形状参数（固定精度）
    precision = 20
    alpha_control = mean_control * precision
    beta_control = (1 - mean_control) * precision
    alpha_treatment = mean_treatment * precision
    beta_treatment = (1 - mean_treatment) * precision

    # 生成甲基化值序列（四舍五入保留三位小数）
    control_vals = beta.rvs(alpha_control, beta_control, size=len(cpg_sites)).round(3)
    treatment_vals = beta.rvs(alpha_treatment, beta_treatment, size=len(cpg_sites)).round(3)

    # 进行一遍一致性校验
    delta_vals = treatment_vals - control_vals
    if (delta_methylation > 0 and np.any(delta_vals < 0)) or (delta_methylation < 0 and np.any(delta_vals > 0)):
        return None

    # 实际计算平均甲基化差值
    actual_delta = round(np.mean(treatment_vals) - np.mean(control_vals), 3)

    # 返回模拟结果
    return {
        "chr": chr_name,
        "start": region_start,
        "end": region_end_final,
        "CpG_sites": cpg_sites,
        "CpG_count": len(cpg_sites),
        "methylation_control": control_vals.tolist(),
        "methylation_treatment": treatment_vals.tolist(),
        "mean_delta_methylation": actual_delta
    }

# 定义 non-DMR 区域模拟函数
def simulate_nondmr_region(chr_name, region_start, region_end, delta_methylation=None, max_cpgs=50, density="auto"):
    """
    模拟一个 non-DMR 区域（非差异甲基化区域）

    条件：
    - CpG数量 < 5 或 delta_methylation < 0.1,（可以人为设置或随机生成）

    参数：
    - chr_name: 染色体名
    - region_start: 起始位置
    - region_end: 终止位置
    - delta_methylation: 可选，实验与对照组甲基化差值，若不提供则自动生成一个 < 0.1 的差值
    - max_cpgs: 最大 CpG 数量（用于控制生成数量）

    返回：
    模拟的 non-DMR 信息字典
    """
    # 定义每两个CpG之间的最小和最大间隔（单位：bp）
    if density == "dense":
        min_distance = 2
        max_distance = 50
    elif density =="sparse":
        min_distance = 50
        max_distance = 499
    else:
        min_distance = 2
        max_distance = 499

    # 如果未指定差值，随机生成一个在 [-0.05, 0.05] 的差值
    if delta_methylation is None:
        delta_methylation = np.round(np.random.uniform(-0.02, 0.02), 3)

    else:
        delta_methylation = np.round(np.random.uniform(-abs(delta_methylation), abs(delta_methylation)), 3)
    # 决定 CpG 数量：可以让其小于5，也可以大于5但差值小于0.1
    allow_large_cpg = np.random.choice([True, False])

    cpg_sites = []
    pos = region_start
    while pos < region_end:
        if len(cpg_sites) >= max_cpgs:
            break
        if not allow_large_cpg and len(cpg_sites) >= 4:
            break  # 控制 CpG 数小于 5

        if len(cpg_sites) > 0:
            max_gap = min(max_distance, region_end - pos)
            if max_gap < min_distance:
                break
            gap = random.randint(min_distance, max_gap)
        else:
            gap = 0
        pos += gap
        if pos < region_end:
            cpg_sites.append(pos)

    if len(cpg_sites) == 0:
        return None

    region_end_final = min(region_end, cpg_sites[-1])

    mean_control = np.random.uniform(0.2, 0.8)
    delta_methylation_modify = round(np.random.normal(loc=delta_methylation, scale=0.01), 3)
    mean_treatment = mean_control + delta_methylation_modify
    
    # 防止mean_treatment 超过阈值
    mean_treatment = min(max(mean_treatment, 0.001), 0.999)

    precision = 20
    alpha_control = mean_control * precision
    beta_control = (1 - mean_control) * precision
    alpha_treatment = mean_treatment * precision
    beta_treatment = (1 - mean_treatment) * precision

    control_vals = beta.rvs(alpha_control, beta_control, size=len(cpg_sites)).round(3)
    treatment_vals = beta.rvs(alpha_treatment, beta_treatment, size=len(cpg_sites)).round(3)

    actual_delta = round(np.mean(treatment_vals) - np.mean(control_vals), 3)

    # 防止波动模拟出真的DMR
    if actual_delta > 0.1 or actual_delta < -0.1:
        return None

    return {
        "chr": chr_name,
        "start": region_start,
        "end": region_end_final,
        "CpG_sites": cpg_sites,
        "CpG_count": len(cpg_sites),
        "methylation_control": control_vals.tolist(),
        "methylation_treatment": treatment_vals.tolist(),
        "mean_delta_methylation": actual_delta
    }

# 模拟不一致DMR
def simulate_dmr_with_inconsistent_points(chr_name, region_start, region_end, delta_methylation, max_cpgs=50, density="auto"):
    """
    模拟一个带有最大程度扰动的不一致DMR。
    每隔4个CpG插入一个方向相反的点 ，无视原始方向 ，强制打乱方向一致性。

    返回：
    包含不一致点扰动后的DMR结构 （插入位置不判断方向）。
    """

    dmr = simulate_dmr_region_with_input_limit(
        chr_name=chr_name,
        region_start=region_start,
        region_end=region_end,
        delta_methylation=delta_methylation,
        max_cpgs=max_cpgs,
        density=density
    )

    if dmr is None:
        return None

    # 每隔4个点插入一个扰动点（例如 3, 7, 11, ...）
    indices = list(range(3, len(dmr["CpG_sites"]), 4))

    for idx in indices:
        # 交换 control 和 treatment 值，制造“方向反转”
        dmr["methylation_control"][idx], dmr["methylation_treatment"][idx] = dmr["methylation_treatment"][idx], dmr["methylation_control"][idx]

    # 更新整体平均差值
    control_arr = np.array(dmr["methylation_control"])
    treatment_arr = np.array(dmr["methylation_treatment"])
    dmr["mean_delta_methylation"] = round(np.mean(treatment_arr - control_arr), 3)

    return dmr

# 模拟不一致DMR，但是再此基础上可以拆分为子subDMR
def simulate_dmr_with_subdmr_points(chr_name, region_start, region_end, delta_methylation, max_cpgs=50):
    """
    模拟一个带有不一致点的 DMR 区域

    参数与原函数相同，额外行为：
    - 在生成好的 DMR 区域中，随机插入 1-3 个不一致点, control 与 treatment 方向相反

    返回：
    与 DMR 模拟函数相同结构的字典
    """
    
    # 检查是否可以把不一致的拆分为子DMR
    def check_split_subdmrs(dmr, flip_indices):

        n_cpgs = len(dmr["CpG_sites"])
        if n_cpgs == 0:
            return None

        # 排序并构建切割点
        split_points = sorted(set(flip_indices))
        segments = []
        start = 0
        for cut in split_points:
            if start < cut:
                segments.append(range(start, cut))  # 不包含 cut 点
            start = cut + 1
        if start < n_cpgs:
            segments.append(range(start, n_cpgs))

        sub_dmrs = []

        for segment in segments:
            if len(segment) < 5:
                continue  # CpG数量不足

            idx_list = list(segment)
            start_pos = dmr["CpG_sites"][idx_list[0]]
            end_pos = dmr["CpG_sites"][idx_list[-1]]
            if end_pos - start_pos < 50:
                continue  # 长度不足

            control = np.array(dmr["methylation_control"])[idx_list]
            treatment = np.array(dmr["methylation_treatment"])[idx_list]
            delta = np.mean(treatment - control)

            if abs(delta) >= 0.1:
                sub_dmrs.append({
                    "chr": dmr["chr"],
                    "start": start_pos,
                    "end": end_pos,
                    "CpG_count": len(idx_list),
                    "mean_delta_methylation": round(delta, 3),
                    "CpG_sites": [dmr["CpG_sites"][i] for i in idx_list],
                    "methylation_control": control.round(3).tolist(),
                    "methylation_treatment": treatment.round(3).tolist()
                })

        return sub_dmrs if sub_dmrs else None

    # 调用原始 DMR 函数
    dmr = simulate_dmr_region_with_input_limit(
        chr_name=chr_name,
        region_start=region_start,
        region_end=region_end,
        delta_methylation=delta_methylation,
        max_cpgs=max_cpgs
    )

    if dmr is None:
        return None  # 模拟失败直接返回

    # 判断差值方向：是 control 大还是 treatment 大
    direction = np.sign(dmr["mean_delta_methylation"])

    # 如果方向为 0（差值极小），则跳过插入
    if direction == 0:
        return None

    # 决定要插入的不一致点数量（1～3）
    n_flip = min(len(dmr["CpG_sites"]), random.randint(1, 3))

    # 随机选出索引来插入不一致点
    indices = random.sample(range(len(dmr["CpG_sites"])), k=n_flip)

    # 交换 selected 点的 control 和 treatment 值，使其“方向相反”
    for idx in indices:
        c_val = dmr["methylation_control"][idx]
        t_val = dmr["methylation_treatment"][idx]
        if (direction > 0 and c_val > t_val) or (direction < 0 and c_val < t_val):
            # 仅在同向情况下才交换，使之变为“不一致”
            dmr["methylation_control"][idx], dmr["methylation_treatment"][idx] = t_val, c_val

    # 更新整体平均差值（方便校验）
    control_arr = np.array(dmr["methylation_control"])
    treatment_arr = np.array(dmr["methylation_treatment"])
    dmr["mean_delta_methylation"] = round(np.mean(treatment_arr - control_arr), 3)

    # 检查是否形成子DMR
    return check_split_subdmrs(dmr, indices)


# 工具函数：从 mean 和 std 推出 beta 分布的 alpha, beta 参数
def beta_params_from_mean_std(mean, std):
    mean = np.clip(mean, 1e-3, 1 - 1e-3)
    var = std ** 2
    common = mean * (1 - mean) / var - 1
    alpha = mean * common
    beta_ = (1 - mean) * common
    return max(alpha, 1e-3), max(beta_, 1e-3)


# 带缺失控制的样本展开函数
def simulate_group_samples_with_missing(
    dmr,
    n_control: int,
    n_treatment: int,
    coverage_mean: int,
    coverage_std: int,
    group_std: float,
    output_dir: str = "/mnt/data/sample_outputs"
):
    os.makedirs(output_dir, exist_ok=True)

    all_samples = {
        f"control_sample_{i}": [] for i in range(1, n_control + 1)
    }
    all_samples.update({
        f"treatment_sample_{i}": [] for i in range(1, n_treatment + 1)
    })

    dmr_missing_rate = np.random.uniform(0.1, 0.5)
    sample_missing_rate = np.random.choice([0.1, 0.2, 0.3])

    cpg_sites = dmr["CpG_sites"]
    n_cpgs = len(cpg_sites)
    n_missing_cpgs = int(np.floor(dmr_missing_rate * n_cpgs))
    missing_cpg_indices = set(random.sample(range(n_cpgs), n_missing_cpgs))

    for group, meth_values, n_samples in [
        ("control", dmr["methylation_control"], n_control),
        ("treatment", dmr["methylation_treatment"], n_treatment)
    ]:
        for sample_id in range(1, n_samples + 1):
            sample_key = f"{group}_sample_{sample_id}"
            for idx, cpg_pos in enumerate(cpg_sites):
                is_missing_cpg = idx in missing_cpg_indices
                is_missing_sample = is_missing_cpg and (np.random.rand() < sample_missing_rate)

                if is_missing_sample:
                    continue
                
                # 采用beta分布
                coverage = int(np.clip(np.random.normal(coverage_mean, coverage_std), 1, 100))
                mean_meth = meth_values[idx]
                #sample_meth = np.clip(np.random.normal(loc=mean_meth, scale=group_std), 0, 1)
                alpha, beta_ = beta_params_from_mean_std(mean_meth, group_std)
                sample_meth = beta.rvs(alpha, beta_)

                all_samples[sample_key].append({
                    "chr": dmr["chr"],
                    "start": cpg_pos,
                    "end": cpg_pos + 1,
                    "coverage": coverage,
                    "methylation_level": round(sample_meth, 4)
                })

    for sample_key, records in all_samples.items():
        df = pd.DataFrame(records)
        df.to_csv(os.path.join(output_dir, f"{sample_key}.tsv"), sep="\t", index=False)

    return {
        "output_dir": output_dir,
        "dmr_missing_rate": round(dmr_missing_rate, 3),
        "sample_missing_rate": sample_missing_rate
    }


# 更新后的验证函数，增加滑窗子区域均值判断逻辑
def validate_nondmr_from_samples(cpg_sites, all_samples, n_control, n_treatment,
                                  delta_threshold=0.1, max_consecutive=3, window_size=5):
    site_index = {pos: idx for idx, pos in enumerate(cpg_sites)}
    control_meth = defaultdict(list)
    treatment_meth = defaultdict(list)

    for sample_key, records in all_samples.items():
        group = "control" if "control" in sample_key else "treatment"
        for record in records:
            pos = record["start"]
            if pos in site_index:
                if group == "control":
                    control_meth[pos].append(record["methylation_level"])
                else:
                    treatment_meth[pos].append(record["methylation_level"])

    deltas = []
    for pos in cpg_sites:
        c_vals = control_meth.get(pos, [])
        t_vals = treatment_meth.get(pos, [])
        if len(c_vals) >= 1 and len(t_vals) >= 1:
            delta = abs(np.mean(t_vals) - np.mean(c_vals))
        else:
            delta = 0
        deltas.append(delta)

    # 1. 检查连续超标个数
    consecutive = 0
    max_consec = 0
    over_threshold_count = 0
    for delta in deltas:
        if delta > delta_threshold:
            over_threshold_count += 1
            consecutive += 1
            max_consec = max(max_consec, consecutive)
        else:
            consecutive = 0

    too_many_consecutive = max_consec > max_consecutive
    too_many_total = over_threshold_count > len(cpg_sites) / 2

    # 2. 滑动窗口检查子区域均值差值
    found_dmr_window = False
    for i in range(len(deltas) - window_size + 1):
        window = deltas[i:i + window_size]
        window_mean_delta = np.mean(window)
        if window_mean_delta > delta_threshold:
            found_dmr_window = True
            break

    passed = not (too_many_consecutive or too_many_total or found_dmr_window)
    return passed, deltas


def simulate_group_samples_with_validation(
    dmr,
    n_control: int,
    n_treatment: int,
    coverage_mean: int,
    coverage_std: int,
    group_std: float,
    output_dir: str = "/mnt/data/sample_outputs",
    max_retries: int = 10
):
    for attempt in range(max_retries):
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)

        all_samples = {
            f"control_sample_{i}": [] for i in range(1, n_control + 1)
        }
        all_samples.update({
            f"treatment_sample_{i}": [] for i in range(1, n_treatment + 1)
        })

        dmr_missing_rate = np.random.uniform(0.1, 0.5)
        sample_missing_rate = np.random.choice([0.1, 0.2, 0.3])

        cpg_sites = dmr["CpG_sites"]
        n_cpgs = len(cpg_sites)
        n_missing_cpgs = int(np.floor(dmr_missing_rate * n_cpgs))
        missing_cpg_indices = set(random.sample(range(n_cpgs), n_missing_cpgs))

        for group, meth_values, n_samples in [
            ("control", dmr["methylation_control"], n_control),
            ("treatment", dmr["methylation_treatment"], n_treatment)
        ]:
            for sample_id in range(1, n_samples + 1):
                sample_key = f"{group}_sample_{sample_id}"
                for idx, cpg_pos in enumerate(cpg_sites):
                    is_missing_cpg = idx in missing_cpg_indices
                    is_missing_sample = is_missing_cpg and (np.random.rand() < sample_missing_rate)

                    if is_missing_sample:
                        continue

                    coverage = int(np.clip(np.random.normal(coverage_mean, coverage_std), 1, 100))
                    mean_meth = meth_values[idx]
                    alpha, beta_ = beta_params_from_mean_std(mean_meth, group_std)
                    sample_meth = beta.rvs(alpha, beta_)

                    all_samples[sample_key].append({
                        "chr": dmr["chr"],
                        "start": cpg_pos,
                        "end": cpg_pos + 1,
                        "coverage": coverage,
                        "methylation_level": round(sample_meth, 4)
                    })

        passed, deltas = validate_nondmr_from_samples(cpg_sites, all_samples, n_control, n_treatment)

        if passed:
            for sample_key, records in all_samples.items():
                df = pd.DataFrame(records)
                df.to_csv(os.path.join(output_dir, f"{sample_key}.tsv"), sep="\t", index=False)

            return {
                "output_dir": output_dir,
                "dmr_missing_rate": round(dmr_missing_rate, 3),
                "sample_missing_rate": sample_missing_rate,
                "delta_check_passed": True,
                "delta_values": deltas
            }

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    return {
        "output_dir": output_dir,
        "delta_check_passed": False,
        "error": "Failed to simulate valid non-DMR after max retries"
    }   


from datetime import datetime
import os

def write_simulation_parameters(
    output_dir: str,
    total_dmr: int,
    mean_delta: float,
    n_control: int,
    n_treatment: int,
    coverage_mean: int,
    coverage_std: int,
    chr_name: str,
    start_pos: int,
    max_cpgs: int,
    dmr_per: float,
    dmr_notable_per: float,
    dmr_inconsis_per: float,
    dmr_sub_per: float,
    seed: int,
    length_mean: int,
    length_std: int,
    density: str = "auto",
    dense_ratio: float = 0.5,
    min_gap: int = 500,
    max_gap: int = 5000,
):
    log_path = os.path.join(output_dir, "para.log")
    with open(log_path, "w") as f:
        f.write("# Simulation Parameters Log\n")
        f.write(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write(f"Chromosome: {chr_name}\n")
        f.write(f"Start position: {start_pos}\n")
        f.write(f"Total simulated regions: {total_dmr}\n")
        f.write(f" - good-DMR: {int(total_dmr * dmr_per)}\n")
        f.write(f" - inconsistent-DMR: {int(total_dmr * dmr_inconsis_per)}\n")
        f.write(f" - sub-DMR: {int(total_dmr * dmr_sub_per)}\n")
        f.write(f" - notable-DMR: {int(total_dmr * dmr_notable_per)}\n")
        f.write(f" - non-DMR: {total_dmr - int(total_dmr * (dmr_per + dmr_inconsis_per + dmr_notable_per + dmr_sub_per))}\n\n")

        f.write(f"DMR density mode: {density}\n")
        if density == "mix":
            f.write(f" - Dense region proportion (dense_ratio): {dense_ratio:.2f}\n")
        f.write(f"Region length: N(mean={length_mean} bp, std={length_std} bp)\n")
        f.write(f"Max CpGs per region: {max_cpgs}\n")
        f.write(f"Region gap between DMRs: {min_gap} ~ {max_gap} bp\n\n")

        f.write(f"Delta methylation mean: {mean_delta}\n")
        f.write(f"Control group samples: {n_control}\n")
        f.write(f"Treatment group samples: {n_treatment}\n")
        f.write(f"Coverage mean: {coverage_mean}\n")
        f.write(f"Coverage std: {coverage_std}\n")
        f.write(f"Random seed: {seed}\n")


def get_random_density(density, dense_ratio):
    if density == "mix":
        return "dense" if random.random() < dense_ratio else "sparse"
    elif density in {"dense", "sparse", "auto"}:
        return density
    else:
        print(f"[Warning] Unknown density mode '{density}', fallback to 'auto'.")
        return "auto"

# 随机打乱生成顺序 + 修正分类输出 + 分类插入标签
def simulate_mixed_regions_randomized(
    total_dmr: int,
    mean_delta: float,
    n_control: int,
    n_treatment: int,
    coverage_mean: int,
    coverage_std: int,
    output_dir: str = "/mnt/data/sample_outputs",
    chr_name: str = "chr1",
    start_pos: int = 10000,
    length_mean: int= 1000,
    length_std: int = 300,
    max_cpgs: int = 50,
    dmr_per: float = 0.3,
    dmr_notable_per: float = 0.05,
    dmr_inconsis_per: float = 0.1,
    dmr_sub_per: float = 0.05,
    density: str = "auto",
    dense_ratio: float = 0.5,
    seed: int = 42
):
    np.random.seed(seed)
    random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    # 计算各类区域数量
    dmr_count = int(total_dmr * dmr_per)
    dmr_notable_count = int(total_dmr * dmr_notable_per)
    dmr_inconsistent_count = int(total_dmr * dmr_inconsis_per)
    dmr_inconsistent_subdmr_count = int(total_dmr * dmr_sub_per)
    nondmr_count = total_dmr - dmr_count - dmr_inconsistent_count - dmr_notable_count - dmr_inconsistent_subdmr_count

    # 构建任务列表（打标签）
    tasks = (
        ["non-DMR"] * nondmr_count +
        ["good-DMR"] * dmr_count +
        ["inconsistent-DMR"] * dmr_inconsistent_count +
        ["sub-DMR"] * dmr_inconsistent_subdmr_count +
        ["notable-DMR"] * dmr_notable_count
    )
    random.shuffle(tasks)

    regions = []
    pos = start_pos

    def advance_position(last_end):
        return last_end + random.randint(500, 2000)

    region_index = 0
    for task in tasks:
        region_length = int(np.random.normal(loc=length_mean, scale=length_std))
        region_length = max(100, region_length)
        end = pos + region_length
        #end = pos + random.randint(200, 2000)

        # 产生一些付的delta，80%概率保持原值（正数），20%概率原值取反（负数）
        sign = 1 if random.random() < 0.8 else -1

        if task == "non-DMR":
            region = simulate_nondmr_region(chr_name, pos, end, max_cpgs=max_cpgs, density=get_random_density(density, dense_ratio))
            group_std = 0.03
        elif task == "inconsistent-DMR":
            delta = round(np.random.normal(loc=mean_delta, scale=0.05), 3) * sign
            region = simulate_dmr_with_inconsistent_points(chr_name, pos, end, delta, max_cpgs=max_cpgs, density=get_random_density(density, dense_ratio))
            group_std = 0.03
        elif task == "sub-DMR":
            delta = round(np.random.normal(loc=mean_delta, scale=0.05), 3) * sign
            region = simulate_dmr_with_subdmr_points(chr_name, pos, end, delta, max_cpgs=max_cpgs)
            group_std = 0.03
        elif task == "notable-DMR":
            delta = round(np.random.normal(loc=mean_delta, scale=0.05), 3) * sign
            region = simulate_dmr_region_with_input_limit(chr_name, pos, end, delta, max_cpgs=max_cpgs, density=get_random_density(density, dense_ratio))
            group_std = round(1.5*abs(delta), 3)
        else:  # standard DMR
            delta = round(np.random.normal(loc=mean_delta, scale=0.05), 3) * sign
            region = simulate_dmr_region_with_input_limit(chr_name, pos, end, delta, max_cpgs=max_cpgs, density=get_random_density(density, dense_ratio))
            group_std = 0.03

        if region:
            if isinstance(region, list):  # 说明是 sub-DMR 情况
                for sub_region in region:
                    sub_region["category"] = task
                    regions.append(sub_region)
        
                    region_output_dir = os.path.join(output_dir, f"{chr_name}_{sub_region['start']}_{sub_region['end']}_{task}")
                    simulate_group_samples_with_missing(
                        dmr=sub_region,
                        n_control=n_control,
                        n_treatment=n_treatment,
                        coverage_mean=coverage_mean,
                        coverage_std=coverage_std,
                        group_std=group_std,
                        output_dir=region_output_dir
                    )
                    region_index += 1
            
            elif task == "non-DMR":
                region["category"] = task
                regions.append(region)

                region_output_dir = os.path.join(output_dir, f"{chr_name}_{region['start']}_{region['end']}_{task}")

                simulate_group_samples_with_validation(
                    dmr=region,
                    n_control=n_control,
                    n_treatment=n_treatment,
                    coverage_mean=coverage_mean,
                    coverage_std=coverage_std,
                    group_std=group_std,
                    output_dir=region_output_dir
                )

            else:
                region["category"] = task
                regions.append(region)
        
                region_output_dir = os.path.join(output_dir, f"{chr_name}_{region['start']}_{region['end']}_{task}")
                simulate_group_samples_with_missing(
                    dmr=region,
                    n_control=n_control,
                    n_treatment=n_treatment,
                    coverage_mean=coverage_mean,
                    coverage_std=coverage_std,
                    group_std=group_std,
                    output_dir=region_output_dir
                )
                region_index += 1
        
            pos = advance_position(end)
        
     # 整理输出
    rows = []
    for region in regions:
        rows.append({
            "chr": region["chr"],
            "start": region["start"],
            "end": region["end"],
            "CpG_count": region["CpG_count"],
            "mean_delta_methylation": region["mean_delta_methylation"],
            "category": region["category"]
        })

    df = pd.DataFrame(rows)

    write_simulation_parameters(
    output_dir=output_dir,
    total_dmr=total_dmr,
    mean_delta=mean_delta,
    n_control=n_control,
    n_treatment=n_treatment,
    coverage_mean=coverage_mean,
    coverage_std=coverage_std,
    chr_name=chr_name,
    start_pos=start_pos,
    max_cpgs=max_cpgs,
    dmr_per=dmr_per,
    dmr_notable_per=dmr_notable_per,
    dmr_inconsis_per=dmr_inconsis_per,
    dmr_sub_per=dmr_sub_per,
    seed=seed,
    length_mean=length_mean,
    length_std=length_std,
    density=density,  # or "sparse" / "auto"
    dense_ratio=dense_ratio,
    min_gap=500,
    max_gap=5000)
    return df



def main():
    parser = argparse.ArgumentParser(description="Simulate DMR regions")

    parser.add_argument("--total_dmr", type=int, default=1000)
    parser.add_argument("--mean_delta", type=float, default=0.3)
    parser.add_argument("--n_control", type=int, default=5)
    parser.add_argument("--n_treatment", type=int, default=5)
    parser.add_argument("--coverage_mean", type=int, default=30)
    parser.add_argument("--coverage_std", type=int, default=5)
    parser.add_argument("--output_dir", type=str, default="/mnt/data/sample_outputs")
    parser.add_argument("--chr_name", type=str, default="chr1")
    parser.add_argument("--start_pos", type=int, default=10000)
    parser.add_argument("--length_mean", type=int, default=1000)
    parser.add_argument("--length_std", type=int, default=300)
    parser.add_argument("--max_cpgs", type=int, default=50)
    parser.add_argument("--dmr_per", type=float, default=0.3)
    parser.add_argument("--dmr_notable_per", type=float, default=0.05)
    parser.add_argument("--dmr_inconsis_per", type=float, default=0.1)
    parser.add_argument("--dmr_sub_per", type=float, default=0.05)
    parser.add_argument("--density", type=str, choices=["dense", "sparse", "auto", "mix"], default="auto")
    parser.add_argument("--dense_ratio", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=42)

    args = parser.parse_args()

    # 调用模拟函数
    df = simulate_mixed_regions_randomized(
        total_dmr=args.total_dmr,
        mean_delta=args.mean_delta,
        n_control=args.n_control,
        n_treatment=args.n_treatment,
        coverage_mean=args.coverage_mean,
        coverage_std=args.coverage_std,
        output_dir=args.output_dir,
        chr_name=args.chr_name,
        start_pos=args.start_pos,
        length_mean=args.length_mean,
        length_std=args.length_std,
        max_cpgs=args.max_cpgs,
        dmr_per=args.dmr_per,
        dmr_notable_per=args.dmr_notable_per,
        dmr_inconsis_per=args.dmr_inconsis_per,
        dmr_sub_per=args.dmr_sub_per,
        density=args.density,
        dense_ratio=args.dense_ratio,
        seed=args.seed
    )

    print("Simulation completed.")
    df.to_csv(f"{args.output_dir}/DMRs.txt",sep="\t",header=True,index=False)

# 主逻辑入口
if __name__ == "__main__":
    #args = parse_args()
    #output_dir = args.output_dir
    #seed = 42
    #result = simulate_mixed_regions_randomized(
    #    total_dmr=2000,
    #    mean_delta=0.3,
    #    n_control=10,
    #    n_treatment=10,
    #    coverage_mean=30,
    #    coverage_std=5,
    #    length_mean=2000,
    #    length_std=300,
    #    output_dir=output_dir,
    #    density="dense"
    #)
    #result.to_csv(f"{output_dir}/DMRs.txt",sep="\t",header=True,index=False)
    main()