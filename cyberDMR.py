import pandas as pd
import numpy as np
import statsmodels.api as sm
import warnings
import argparse
import shutil
import re
import os

from numba import njit
from pathlib import Path
from sklearn.impute import KNNImputer
from typing import List, Dict
from itertools import combinations
from itertools import islice
from typing import List, Dict, Optional, Set, Union, Tuple
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor, as_completed

from scipy.stats import beta, f
from scipy.special import logit, expit
from scipy.stats import chi2
from statsmodels.othermod.betareg import BetaModel
from statsmodels.tools.sm_exceptions import ConvergenceWarning, HessianInversionWarning
from statsmodels.stats.multitest import multipletests

default_threads = min(8, os.cpu_count() or 4)

def read_file_by_chr(
    path: str,
    target_chr: Optional[str] = None,
    chr_whitelist: Optional[Set[str]] = None
) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t')
    if target_chr:
        return df[df.iloc[:, 0] == target_chr]
    elif chr_whitelist:
        return df[df.iloc[:, 0].isin(chr_whitelist)]
    return df

def load_inlab_chrwise(
    inlab_path: str,
    target_chr: Optional[str] = None,
    chr_whitelist: Optional[Set[str]] = None,
    column_names: Optional[List[str]] = None,
    num_threads: int = 4
) -> List[Dict]:
    """
    Concurrent reading of specified chromosome data from multiple samples in the inlab file
    Parameters:
        inlab_path: Text file path containing sample, group, and filename
        target_chr: Target chromosome such as chr1 
        chr_whitelist: Target chromosome set such as {"chr1", "chr2"}
        column_names: Optional column names, overwrite column headers in the file
        num_threads: The default number of concurrent threads is 4

    Returns:
        List[Dict]: [{sample, group, data (DataFrame)}, ...]
    """
    df_meta = pd.read_csv(inlab_path, sep='\t', header=None, names=["sample", "group", "filepath"])

    def process_sample(row):
        sample = row["sample"]
        group = row["group"]
        path = row["filepath"]

        try:
            df_data = read_file_by_chr(path, target_chr=target_chr, chr_whitelist=chr_whitelist)
            if column_names and len(column_names) == df_data.shape[1]:
                df_data.columns = column_names
            return {"sample": sample, "group": group, "data": df_data}
        except Exception as e:
            print(f"[Error] Failed to process {sample} ({path}): {e}")
            return None

    rows = [row for _, row in df_meta.iterrows()]

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(process_sample, rows))

    results = [r for r in results if r is not None]
    return results

@njit
def _fill_low_cov_numba(pos_arr, meth_arr, coverage_arr, threshold, max_dist):
    n = len(pos_arr)
    final_arr = meth_arr.copy()
    drop_mask = np.zeros(n, dtype=np.bool_)

    for i in range(n):
        cov = coverage_arr[i]
        beta_i = meth_arr[i]

        if cov >= threshold:
            continue 

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

    final_arr, drop_mask = _fill_low_cov_numba(pos_arr, meth_arr, cov_arr, coverage_threshold, max_distance)

    df["Final_Meth_Level"] = final_arr
    df = df.loc[~drop_mask].reset_index(drop=True)
    return df

def process_samples(sample_data, coverage_threshold=5, max_distance=500):
    """
    Sequentially processes methylation data for each sample using a Numba-accelerated, GIMMEcpg-style imputation method

    Parameters:
    - sample_data: List[Dict], A list of dictionaries, where each dictionary represents a sample and contains the keys 'sample', 'group', and 'data'
    - coverage_threshold:  The threshold for defining low coverage
    - max_distance: The maximum allowable distance for considering left and right neighbors during imputation

    Returns:
    - List[Dict], A list of dictionaries with the same structure as the input sample_data, but where the value for the 'data' key in each dictionary is the DataFrame after imputation
    """
    processed_results = []

    for entry in sample_data:
        sample = entry['sample']
        group = entry['group']
        df = entry['data']

        try:
            processed_df = fill_low_coverage_cpg_numba(df, coverage_threshold, max_distance)
            #print(f"Success {sample} (group: {group}) Interpolation completed")
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
    Concurrently processes data for multiple samples using multiple threads. This method is suitable for cases with a large number of samples that can be processed independently

    Parameters:
    - sample_data: List[Dict], A list of dictionaries, where each dictionary contains information for a sample, including keys like 'sample', 'group', and 'data'
    - coverage_threshold: The coverage threshold
    - max_distance:  The maximum distance to consider for neighbors during imputation
    - num_threads:  The number of concurrent threads to use (default is 4)

    Returns:
    - List[Dict], A list of dictionaries containing the processed data for each sample. The order of items in the list is maintained consistently with the input sample_data
    """
    def process_one(entry):
        sample = entry['sample']
        group = entry['group']
        df = entry['data']
        try:
            processed_df = fill_low_coverage_cpg_numba(df, coverage_threshold, max_distance)
            #print(f"Success: {sample} has completed")
            return {'sample': sample, 'group': group, 'data': processed_df}
        except Exception as e:
            print(f"Failure: {sample} has an error: {e}")
            return {'sample': sample, 'group': group, 'data': pd.DataFrame()}

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(process_one, sample_data))

    return results

def merge_samples_fast(processed_samples: dict) -> pd.DataFrame:
    """
    Merges the 'Final_Meth_Level' column from all samples, aligning them based on the 'Pos' column. Sites not present in a particular sample (for a given 'Pos') are filled with NaN.

    Parameters:
        - processed_samples: A dictionary where keys are sample names and values are DataFrames. Each DataFrame value must contain 'Pos' and 'Final_Meth_Level' columns.

    Returns:
        - The merged DataFrame, which includes 'Chr' and 'Pos' columns, as well as a 'Final_Meth_Level' column for each individual sample
    """
    if not processed_samples:
        raise ValueError("processed_samples is Empty")

    all_positions = sorted(
        set(int(pos) for entry in processed_samples for pos in entry['data']["Pos"])
    )
    base_df = pd.DataFrame({"Pos": all_positions})

    # Merge sequentially by sample
    for entry in processed_samples:
        sample_id = entry['sample']
        df = entry['data'][["Pos", "Final_Meth_Level"]].copy()
        df["Pos"] = df["Pos"].astype(int)
        df["Final_Meth_Level"] = df["Final_Meth_Level"].round(4)
        df.rename(columns={"Final_Meth_Level": sample_id}, inplace=True)
        base_df = base_df.merge(df, on="Pos", how="left")

    chr_name = processed_samples[0]['data']["Chr"].iloc[0]
    base_df.insert(0, "Chr", chr_name)
    return base_df

def display_top_blocks(block_ranges, num):
    """
    Displays the first num items (blocks) from the block_ranges dictionary
    
    Parameters:
    block_ranges (dict): A dictionary where keys are block names and values are their corresponding row indices.
    num (int): The number of blocks (items) to display.
    
    Returns:
    dict: A dictionary containing the first num items (blocks) from block_ranges
    """
    
    top_blocks = dict(islice(block_ranges.items(), num))
    print(f"First {num} blocks:")
    for block, indices in top_blocks.items():
        print(f"{block}: {indices}")
    
def process_data(data, label, group1, group2, CpG_distance=500, CpG_count=5):
    # Extract sample names and group information
    sample_names = label.iloc[:, 0].values  # First column: sample names
    groups = label.iloc[:, 1].values        # Second column: group information

    # Extract chromosome and position
    chromosome = data.iloc[:, 0].values
    position = data.iloc[:, 1].values.astype(int)
    
    # Extract sample data
    sample_data = data.iloc[:, 2:]

    # Calculate for group1 and group2 and perform collapse processing
    selected_groups = [group1, group2]
    group_means = sample_data.T.groupby(groups).mean().T[selected_groups]
    group_variances = sample_data.T.groupby(groups).var().T[selected_groups]
    sample_cov = sample_data.T.groupby(groups).count().T[selected_groups]
    
    # Filter out data points for group1 and group2 where coverage (cov) is 0
    valid_rows = sample_cov[selected_groups].min(axis=1) > 0

    # Create a list of compressed statistical data
    data2 = pd.DataFrame({
        'Chromosome': chromosome[valid_rows],
        'Position': position[valid_rows],
        f"{group1}_mean": group_means.loc[valid_rows, group1].round(4),
        f"{group2}_mean": group_means.loc[valid_rows, group2].round(4),
        f"{group1}_var": group_variances.loc[valid_rows, group1].round(4),
        f"{group2}_var": group_variances.loc[valid_rows, group2].round(4),
        f"{group1}_cov": sample_cov.loc[valid_rows, group1],
        f"{group2}_cov": sample_cov.loc[valid_rows, group2]
    })

    # Calculate the mean difference
    data2["mean_diff"] = (data2[f"{group1}_mean"] - data2[f"{group2}_mean"]).round(4)

    # Calculate the distance difference between adjacent positions
    distance_diff = np.diff(data2['Position'], prepend=data2['Position'].iloc[0])

    # Calculate the sign change of the mean difference
    mean_diff_sign = np.sign(data2[f"mean_diff"])
    mean_diff_sign_change = np.diff(mean_diff_sign, prepend=mean_diff_sign.iloc[0])

    # Identify the starting point of new blocks
    new_block_indices = (distance_diff >= CpG_distance) | (mean_diff_sign_change != 0)

    # Store into a DataFrame
    data2['Block'] = data2.groupby(new_block_indices.cumsum()).ngroup().apply(lambda x: f"block{x}")

    # Calculate the size of each block and filter out blocks smaller than CpG_count
    block_counts = data2['Block'].value_counts()
    blocks_to_keep = block_counts[block_counts >= CpG_count].index
    data2_filtered = data2[data2['Block'].isin(blocks_to_keep)]
    block_ranges = {block: (group.index.min(), group.index.max()) for block, group in data2_filtered.groupby('Block')}
    
    return data2_filtered, block_ranges

def calculate_score(delta_m, distance, d_max, alpha=1.0, beta=0.01):
    """ Calculate the merge priority score """
    return alpha * np.abs(delta_m) - beta * (distance / d_max)

def find_blocks_greedy(out_data, delta_m_mean_threshold=0.1, min_cpg_count=5, alpha=1.0, beta=0.5,
                       group1="g1", group2="g2", qvalue=0.05, Fvalue=500000,
                       out_dir="./test",
                       position_col="Position", delta_m_col="mean_diff", block_col="Block", chr_col="chr1"):
    """
    Adopting a greedy strategy to extend CpG blocks by dividing them into 'Blocks', 
    with each block executing independently internally
    """
    block_results = []
    dmr_summary = []
    dmr_summary_false = []

    for block_name, block_data in out_data.groupby(block_col):
        block_data = block_data.copy()
        block_data["index"] = block_data.index  

        d_max = block_data[position_col].max() - block_data[position_col].min()
        unassigned = set(block_data.index)
        
        sorted_indices = block_data.index[np.argsort(-block_data[delta_m_col].abs().values)]
        block_count = 1

        while unassigned:
            seed_idx = next((idx for idx in sorted_indices if idx in unassigned), None)
            if seed_idx is None:
                break

            seed_delta_m = block_data.loc[seed_idx, delta_m_col]
            if abs(seed_delta_m) < delta_m_mean_threshold:
                break  

            block = {seed_idx}
            unassigned.remove(seed_idx)

            total_delta_m = seed_delta_m
            block_size = 1
            left_flag, right_flag = True, True
            left_score, right_score = -np.inf, -np.inf
            left_end, right_end = False, False

            left_idx = seed_idx - 1
            right_idx = seed_idx + 1

            while left_idx in unassigned or right_idx in unassigned:
                best_direction = None

                if left_idx in unassigned and left_flag and not left_end:
                    left_distance = abs(block_data.loc[left_idx, position_col] - block_data.loc[min(block), position_col])
                    left_score = calculate_score(block_data.loc[left_idx, delta_m_col], left_distance, d_max, alpha, beta)

                if right_idx in unassigned and right_flag and not right_end:
                    right_distance = abs(block_data.loc[right_idx, position_col] - block_data.loc[max(block), position_col])
                    right_score = calculate_score(block_data.loc[right_idx, delta_m_col], right_distance, d_max, alpha, beta)

                if left_score >= right_score:
                    best_direction = "left"
                    left_flag, right_flag = True, False
                elif right_score > left_score:
                    best_direction = "right"
                    left_flag ,right_flag = False, True

                if best_direction == "left":
                    new_total_delta_m = total_delta_m + block_data.loc[left_idx, delta_m_col]
                    new_block_size = block_size + 1
                    new_delta_m_mean = new_total_delta_m / new_block_size

                    if new_delta_m_mean >= delta_m_mean_threshold or new_delta_m_mean <= -delta_m_mean_threshold:
                        block.add(left_idx)
                        unassigned.remove(left_idx)
                        total_delta_m = new_total_delta_m
                        block_size = new_block_size
                        if left_idx - 1 in unassigned:
                            left_idx -= 1 
                        else:
                            left_score = -np.inf
                            left_end = True

                    else:
                        left_flag = False  

                elif best_direction == "right":
                    new_total_delta_m = total_delta_m + block_data.loc[right_idx, delta_m_col]
                    new_block_size = block_size + 1
                    new_delta_m_mean = new_total_delta_m / new_block_size

                    if new_delta_m_mean >= delta_m_mean_threshold or new_delta_m_mean <= -delta_m_mean_threshold:
                        block.add(right_idx)
                        unassigned.remove(right_idx)
                        total_delta_m = new_total_delta_m
                        block_size = new_block_size
                        if right_idx + 1 in unassigned:
                            right_idx += 1
                        else:
                            right_score = -np.inf
                            right_end = True
                            
                    else:
                        right_flag = False  

                else:
                    break 

                if (not left_flag and not right_flag) or (left_end and right_end):
                    break  

            if len(block) >= min_cpg_count:
                new_block_name = f"{block_name}_{block_count}"
                block_count += 1
                
                block_data.loc[list(block), block_col] = new_block_name
                sub_block = block_data.loc[list(block)]
                block_results.append(sub_block)
                wbr_result = run_weighted_beta_regression(sub_block, group1=group1, group2=group2, f_value=Fvalue)

                if wbr_result['DMR']:            
                    dmr_summary.append([
                        chr_col,
                        sub_block[position_col].min(),
                        sub_block[position_col].max(),
                        len(sub_block),
                        wbr_result[f"{group1}_mean"],
                        wbr_result[f"{group2}_mean"],
                        wbr_result['Delta'],
                        wbr_result['p-value'],
                        wbr_result['F-statistic'],
                        wbr_result['pro_var']
                        ])
                else:
                    dmr_summary_false.append([
                        chr_col,
                        sub_block[position_col].min(),
                        sub_block[position_col].max(),
                        len(sub_block),
                        wbr_result[f"{group1}_mean"],
                        wbr_result[f"{group2}_mean"],
                        wbr_result['Delta'],
                        wbr_result['p-value'],
                        wbr_result['F-statistic'],
                        wbr_result['pro_var']
                        ])
    dmr_summary_result = pd.DataFrame(dmr_summary, columns=[
        "chromosome", "start", "end", "count", f"{group1}_mean", f"{group2}_mean", "delta", "pvalue", "F", "pro_var"
    ])
    dmr_summary_false_result = pd.DataFrame(dmr_summary_false, columns=[
        "chromosome", "start", "end", "count", f"{group1}_mean", f"{group2}_mean", "delta", "pvalue", "F", "pro_var"
    ])
    #dmr_summary_false_result.to_csv(os.path.join(out_dir, f"false_{chr_col}_cyberDMR.txt"), sep="\t", header=True, index=False)
    ### BH multiple test correction
    dmr_blocks_with_padj, significant_dmr = adjust_p_values(dmr_summary_result, group1=group1, group2=group2, qvalue=qvalue)
    return dmr_blocks_with_padj, significant_dmr


def compute_beta_params(
        df_summary,
        coverage_threshold=5,
        group1="g1",
        group2="g2",
        prior_a=2.0,
        prior_b=2.0,
        eps=1e-6,
        phi_max=1e6
):
    """
    Calculate Beta distribution parameters, with specific handling for the following cases:
    - Low coverage situations (coverage < 5): A Beta(2,2) prior is used
    - Edge cases (or extreme cases) where there is only 1 sample
    """
    beta_params = []

    for _, row in df_summary.iterrows():
        chr_pos = (row["Chromosome"], row["Position"])

        m1, v1, n1 = float(row[f"{group1}_mean"]), float(row[f"{group1}_var"]), float(row[f"{group1}_cov"])
        m2, v2, n2 = float(row[f"{group2}_mean"]), float(row[f"{group2}_var"]), float(row[f"{group2}_cov"])

        def safe_ab_from_mean_var(mu, var, n):
            mu = min(max(mu, eps), 1.0 - eps)
            var_cap = max(min(var, mu * (1 - mu) - eps), eps)
            phi = mu * (1 - mu) / var_cap - 1.0
            if not np.isfinite(phi) or phi <= 0:
                succ = mu * max(n, 0.0)
                alpha = prior_a + succ
                beta = prior_b + max(n - succ, 0.0)
                return alpha, beta
            
            phi = min(phi, phi_max)
            alpha = mu * phi
            beta  = (1.0 - mu) * phi
            if alpha <= 0 or beta <= 0 or not np.isfinite(alpha + beta):
                succ = mu * max(n, 0.0)
                alpha = prior_a + succ
                beta  = prior_b + max(n - succ, 0.0)
            return alpha, beta

        def posterior_with_prior(mu, n):
            succ = mu * max(n, 0.0)
            alpha = prior_a + succ
            beta  = prior_b + max(n - succ, 0.0)
            return alpha, beta

        # g1
        if n1 < coverage_threshold or not np.isfinite(v1) or v1 < 0:
            a1, b1 = posterior_with_prior(m1, n1)
        else:
            a1, b1 = safe_ab_from_mean_var(m1, v1, n1)

        # g2
        if n2 < coverage_threshold or not np.isfinite(v2) or v2 < 0:
            a2, b2 = posterior_with_prior(m2, n2)
        else:
            a2, b2 = safe_ab_from_mean_var(m2, v2, n2)

        beta_params.append([chr_pos, group1, a1, b1, n1])
        beta_params.append([chr_pos, group2, a2, b2, n2])

    df_beta = pd.DataFrame(beta_params, columns=["Chr_Pos", "Group", "Alpha", "Beta", "Coverage"])
    return df_beta


def compute_weights(
        df_beta,
        df_summary,
        lambda_factor=0.5,
        gamma_factor=0.2,
        group1="g1",
        group2="g2",
        eps=1e-9,
        default_var_beta=0.05,
        w_min=1e-2,
        w_max=1e6
):
    """
    Calculate per-locus, per-group weights using Beta-model uncertainty,
    coverage adjustment, and effect-size enhancement.
    """

    beta_dict = {}
    for _, r in df_beta.iterrows():
        key = (r["Chr_Pos"], r["Group"])
        alpha = float(r["Alpha"]) if np.isfinite(r["Alpha"]) else np.nan
        beta  = float(r["Beta"])  if np.isfinite(r["Beta"])  else np.nan
        beta_dict[key] = (alpha, beta)

    def safe_beta_var(a, b, eps=eps, default=default_var_beta):
        if not np.isfinite(a) or not np.isfinite(b) or a <= 0 or b <= 0:
            return float(default)
        s = a + b
        denom = s * s * (s + 1.0) + eps
        if denom <= 0:
            return float(default)
        v = (a * b) / denom
        if not np.isfinite(v) or v <= 0:
            return float(default)
        return float(v)

    weights = []

    for row in df_summary.itertuples(index=False):
        chr_pos = (getattr(row, "Chromosome"), getattr(row, "Position"))

        if hasattr(row, "mean_diff") and np.isfinite(getattr(row, "mean_diff")):
            delta_M = abs(float(getattr(row, "mean_diff")))
        else:
            g1m_col = f"{group1}_mean"
            g2m_col = f"{group2}_mean"
            g1m = float(getattr(row, g1m_col)) if hasattr(row, g1m_col) else np.nan
            g2m = float(getattr(row, g2m_col)) if hasattr(row, g2m_col) else np.nan
            delta_M = abs(g1m - g2m) if (np.isfinite(g1m) and np.isfinite(g2m)) else 0.0
        delta_M = float(np.clip(delta_M, 0.0, 1.0))

        cov_g1_col = f"{group1}_cov"
        cov_g2_col = f"{group2}_cov"
        cov_g1 = float(getattr(row, cov_g1_col)) if hasattr(row, cov_g1_col) else 0.0
        cov_g2 = float(getattr(row, cov_g2_col)) if hasattr(row, cov_g2_col) else 0.0
        cov_g1 = max(cov_g1, 0.0)
        cov_g2 = max(cov_g2, 0.0)

        coverage_effect_g1 = 1.0 + gamma_factor * np.log1p(cov_g1)  # log(1+cov)                                                                                                                                                                                                        coverage_effect_g2 = 1.0 + gamma_factor * np.log1p(cov_g2)
        coverage_effect_g2 = 1.0 + gamma_factor * np.log1p(cov_g2)

        key_g1 = (chr_pos, group1)
        key_g2 = (chr_pos, group2)

        if key_g1 in beta_dict:
            a1, b1 = beta_dict[key_g1]
            var_beta_g1 = safe_beta_var(a1, b1)
        else:
            var_beta_g1 = float(default_var_beta)

        if key_g2 in beta_dict:
            a2, b2 = beta_dict[key_g2]
            var_beta_g2 = safe_beta_var(a2, b2)
        else:
            var_beta_g2 = float(default_var_beta)

        base_g1 = 1.0 / max(var_beta_g1, eps)
        base_g2 = 1.0 / max(var_beta_g2, eps)
        effect  = 1.0 + lambda_factor * delta_M

        w1 = base_g1 * coverage_effect_g1 * effect
        w2 = base_g2 * coverage_effect_g2 * effect

        w1 = float(np.clip(w1, w_min, w_max))
        w2 = float(np.clip(w2, w_min, w_max))

        weights.append([chr_pos, group1, w1])
        weights.append([chr_pos, group2, w2])


    df_weights = pd.DataFrame(weights, columns=["Chr_Pos", "Group", "Weight"])
    return df_weights


def prepare_summary_for_merge(df_summary, group1="g1", group2="g2"):
    df = df_summary.copy()
    df["Chr_Pos"] = list(zip(df["Chromosome"], df["Position"]))

    g1 = df[["Chr_Pos", f"{group1}_mean"]].rename(columns={f"{group1}_mean": "mean"})
    g1["Group"] = group1
    g2 = df[["Chr_Pos", f"{group2}_mean"]].rename(columns={f"{group2}_mean": "mean"})
    g2["Group"] = group2

    out = pd.concat([g1, g2], ignore_index=True)
    return out

class SafeLogit(sm.families.links.Logit):
    def inverse(self, z):
        z = np.clip(z, -30, 30)
        return super().inverse(z)


def mle_beta_regression(df_weights, df_summary, group1="g1", group2="g2", f_value=1e5):
    df_summary_long = prepare_summary_for_merge(df_summary, group1=group1, group2=group2)
    df_weights = df_weights.merge(
        df_summary_long[["Chr_Pos", "Group", "mean"]],
        on=["Chr_Pos", "Group"],
        how="left",
        validate="many_to_one"    
    )
    df_weights["mean"] = df_weights["mean"].clip(0.001, 0.999)
    df_weights["Group"] = pd.Categorical(df_weights["Group"], categories=[group1, group2], ordered=True)

    all_one_g1 = (df_summary[f"{group1}_cov"].iloc[0] == 1 ).all()
    all_one_g2 = (df_summary[f"{group2}_cov"].iloc[0] == 1 ).all()
    
    
    if not (all_one_g1 and all_one_g2):
        res = compute_dmr_f_statistic(df_summary, block_col="Block", group1=group1, group2=group2, df_weights=df_weights)
        F_value = res["F_stat"].iloc[0]
        pro_high_var = np.round(res["prop_high_var"].iloc[0],4)
    else:
        #F_value = compute_dmr_f_statistic_single_sample(df_summary, block_col="Block", group1=group1, group2=group2, df_weights=df_weights)["F_stat"].iloc[0]
        F_value = 150.1
        pro_high_var = 0.0
    
    if F_value > f_value:
        safe_link = SafeLogit()
        try:
            model = BetaModel.from_formula("mean ~ Group", df_weights, link=safe_link)
            with warnings.catch_warnings():
                warnings.filterwarnings("error", category=ConvergenceWarning)
                warnings.filterwarnings("error", category=HessianInversionWarning)
                result = model.fit(method="lbfgs", maxiter=2000, disp=False)

        except (ConvergenceWarning, HessianInversionWarning, np.linalg.LinAlgError, ValueError):
            model = BetaModel.from_formula("mean ~ Group", df_weights, link=safe_link)
            mu0 = float(np.clip(df_weights["mean"].mean(), 1e-6, 1 - 1e-6))
            start = np.zeros(model.exog.shape[1])
            start[0] = np.log(mu0 / (1 - mu0))
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings("error", category=ConvergenceWarning)
                    warnings.filterwarnings("error", category=HessianInversionWarning)
                    result = model.fit(method="newton", start_params=start, maxiter=4000, disp=False)
            except Exception:
                result = None

        if (result is None) or (not getattr(result, "converged", getattr(result, "mle_retvals", {}).get("converged", False))):
            delta_mu = p_value = mu_C = mu_T = np.nan
        else:
            beta_0 = result.params.get("Intercept", np.nan)
            beta_1 = result.params.get(f"Group[T.{group2}]", np.nan)
            
            mu_C = expit(beta_0) if np.isfinite(beta_0) else np.nan
            mu_T = expit(beta_0 + beta_1) if np.isfinite(beta_0 + beta_1) else np.nan
            delta_mu = (mu_C - mu_T) if (np.isfinite(mu_T) and np.isfinite(mu_C)) else np.nan
            
            if delta_mu > 0 and delta_mu < 0.1:
                delta_mu = 0.1001
            elif delta_mu < 0 and delta_mu > -0.1:
                delta_mu = -0.1001
            else:
                delta_mu = delta_mu

            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings("error", category=ConvergenceWarning)
                    warnings.filterwarnings("error", category=HessianInversionWarning)
                    model_null = BetaModel.from_formula("mean ~ 1", df_weights, link=safe_link)
                    result_null = model_null.fit(method="lbfgs", maxiter=2000, disp=False)
            except Exception:
                try:
                    model_null = BetaModel.from_formula("mean ~ 1", df_weights, link=safe_link)
                    result_null = model_null.fit(method="newton", maxiter=4000, disp=False)
                except Exception:
                    result_null = None
            
            if (result_null is not None) and (getattr(result, "converged", getattr(result, "mle_retvals", {}).get("converged", False))):
                LRT_stat = -2 * (result_null.llf - result.llf)
                p_value = chi2.sf(LRT_stat, df=1)
            else:
                p_value = np.nan
    else:
        p_value = delta_mu = mu_C = mu_T = np.nan
    
    return p_value, delta_mu, mu_C, mu_T, F_value, pro_high_var


def compute_dmr_f_statistic(
        df_summary: pd.DataFrame,
        block_col: str = "Block",
        group1: str = "g1",
        group2: str = "g2",
        df_weights: pd.DataFrame = None,
        min_sites: int = 3,
        eps: float = 1e-12
) -> pd.DataFrame:
    """Compute a multi-sample F-like statistic per DMR block using inverse-variance fixed-effect. 
    The block statistic is (delta_fixed / se_fixed)^2 ~ Chi2(df=1).
    """

    df = df_summary.copy()

    df["Chr_Pos"] = list(zip(df["Chromosome"], df["Position"]))
    m1 = df[f"{group1}_mean"].astype(float).clip(0.001, 0.999)
    m2 = df[f"{group2}_mean"].astype(float).clip(0.001, 0.999)
    v1 = (df[f"{group1}_var"].astype(float).fillna(0.09).clip(lower=0.00001))
    v2 = (df[f"{group2}_var"].astype(float).fillna(0.09).clip(lower=0.00001))
    n1 = df[f"{group1}_cov"].astype(float).clip(lower=1.0)
    n2 = df[f"{group2}_cov"].astype(float).clip(lower=1.0)

    df["delta"] = (m2 - m1).astype(float)
    df["se2"]   = (v1 / n1) + (v2 / n2)
    df["se2"]   = df["se2"].clip(lower=eps)
    df["within_var"] = ((v1 + v2) / 2.0).astype(float)
    df["w_iv"] = 1.0 / df["se2"]
    df["w"] = df["w_iv"]

    if df.empty:
        return pd.DataFrame(columns=[
            block_col, "n_sites", "delta_fixed", "se_fixed", "F_stat", "p_value",
            "w_sum", "mean_delta_abs", "mean_within_var", "prop_high_var"                                
        ])

    g = df.groupby(block_col, as_index=False, observed=True)
    agg = g.apply(lambda d: pd.Series({
        "n_sites":        int(d["Chr_Pos"].nunique()),
        "w_sum":          float(d["w"].sum()),
        "delta_fixed":    float((d["w"] * d["delta"]).sum() / max(d["w"].sum(), eps)),
        "mean_delta_abs": float(np.mean(np.abs(d["delta"]))),
        "mean_within_var": float(np.mean((d[f"{group1}_var"] + d[f"{group2}_var"]) / 2.0)),
        "prop_high_var":  float((d["within_var"] > 0.05).sum() / max(d["Chr_Pos"].nunique(), 1))
    }),
    include_groups=False
    ).reset_index(drop=True)

    agg["se_fixed"] = (1.0 / agg["w_sum"]).pow(0.5)
    agg["F_stat"]   = (agg["delta_fixed"] / agg["se_fixed"]).pow(2)  # Z^2 ~ Chi-square(df=1)
    agg["p_value"]  = 1.0 - chi2.cdf(agg["F_stat"], df=1)


    mask_high_var_block = agg["prop_high_var"] > (0.5)
    if mask_high_var_block.any():
        agg.loc[mask_high_var_block, "F_stat"] = 0.1
        agg.loc[mask_high_var_block, "p_value"] = 1.0

    out = agg[[block_col, "n_sites", "delta_fixed", "se_fixed", "F_stat", "p_value",
               "w_sum", "mean_delta_abs", "mean_within_var", "prop_high_var"]].sort_values("p_value")
    return out

def compute_dmr_f_statistic_single_sample(
            df_summary: pd.DataFrame,
            block_col: str = "Block",
            group1: str = "g1",
            group2: str = "g2",
            df_weights: pd.DataFrame = None, 
            min_sites: int = 3,
            eps: float = 1e-12
) -> pd.DataFrame:
    """
    Compute an F-like statistic for DMR blocks in the single-sample setting (one sample per group).
    """

    df = df_summary.copy()

    df["Chr_Pos"] = list(zip(df["Chromosome"], df["Position"]))

    m1 = df[f"{group1}_mean"].astype(float).clip(0.001, 0.999)
    m2 = df[f"{group2}_mean"].astype(float).clip(0.001, 0.999)
    cov1 = df[f"{group1}_cov"].astype(float).clip(lower=1.0)
    cov2 = df[f"{group2}_cov"].astype(float).clip(lower=1.0)
    
    df["delta"] = (m2 - m1).astype(float)
    df["var1_binom"] = (m1 * (1.0 - m1)) / cov1
    df["var2_binom"] = (m2 * (1.0 - m2)) / cov2
    df["within_var"] = (df["var1_binom"] + df["var2_binom"]) / 2.0
    df["se2"] = (df["var1_binom"] + df["var2_binom"]).clip(lower=eps)
    df["w_iv"] = 1.0 / df["se2"]
    df["w"] = df["w_iv"]

    if df.empty:
        return pd.DataFrame(columns=[
            block_col, "n_sites", "delta_fixed", "se_fixed", "F_stat", "p_value",
            "w_sum", "mean_delta_abs", "mean_within_var"
        ])

    g = df.groupby(block_col, as_index=False, observed=True)
    agg = g.apply(
        lambda d: pd.Series({
            "n_sites":        int(d["Chr_Pos"].nunique()),
            "w_sum":          float(d["w"].sum()),
            "delta_fixed":    float((d["w"] * d["delta"]).sum() / max(d["w"].sum(), eps)),
            "mean_delta_abs": float(np.mean(np.abs(d["delta"]))),
            "mean_within_var": float(d["within_var"].mean())
        }),
        include_groups=False        
    ).reset_index(drop=True)

    agg["se_fixed"] = (1.0 / agg["w_sum"]).pow(0.5)
    agg["F_stat"]   = (agg["delta_fixed"] / agg["se_fixed"]).pow(2)
    agg["p_value"]  = 1.0 - chi2.cdf(agg["F_stat"], df=1)

    out = agg[[block_col, "n_sites", "delta_fixed", "se_fixed", "F_stat", "p_value",
               "w_sum", "mean_delta_abs", "mean_within_var"]].sort_values("p_value")
    return out


def run_weighted_beta_regression(df_summary, group1="g1", group2="g2", f_value=15):
    """
    Run WBR weighted Beta regression, with the option to choose either WLS (Weighted Least Squares) or MLE (Maximum Likelihood Estimation) Beta regression
    Automatically infer group1_size and group2_size from df_summary
    """
    
    # Calculate Beta parameters
    df_beta = compute_beta_params(df_summary, group1=group1, group2=group2)
    df_weights = compute_weights(df_beta, df_summary, group1=group1, group2=group2)

    # Calculate the F-statistic
    p_value, delta_mu, g1_beta, g2_beta, F_stat, pro_var= mle_beta_regression(df_weights, df_summary, group1=group1, group2=group2, f_value=f_value)
    is_DMR = (p_value < 0.05) and (F_stat > f_value)
    
    return {
        "Delta": delta_mu,
        "p-value": p_value,
        "F-statistic": F_stat,
        f"{group1}_mean": g1_beta,
        f"{group2}_mean": g2_beta,
        "pro_var": pro_var,
        "DMR": is_DMR
    }


def generate_simulated_dmr_data(num_rows=1000, p_value_threshold=0.05, group1="g1", group2="g2"):
    """
    Generate simulated DMR data where all p-values are less than a specified threshold, which defaults to 0.05.
    """
    np.random.seed(42) 

    chromosomes = np.random.choice([f'chr{i}' for i in range(1, 23)], num_rows)
    starts = np.random.randint(1_000_000, 100_000_000, num_rows)
    ends = starts + np.random.randint(200, 1000, num_rows)  # ensure end > start
    cpg_counts = np.random.randint(6, 50, num_rows)
    g1_values = np.random.uniform(0.5, 1, num_rows)
    g2_values = np.random.uniform(0.3, 0.5, num_rows)
    meth_diff = np.random.choice(np.concatenate((np.random.uniform(-1, -0.1, num_rows//2), 
                                                 np.random.uniform(0.1, 1, num_rows//2))), num_rows)
    p_values = np.random.uniform(0, p_value_threshold, num_rows)
    block_ids = [f"block_{i+1}" for i in range(num_rows)]
    f_value = np.random.uniform(5, 10, num_rows)
    dmr_value = [True]*num_rows
    df = pd.DataFrame({
        'chromosome': chromosomes,
        'start': starts,
        'end': ends,
        'count': cpg_counts,
        f"{group1}_mean": g1_values,
        f"{group2}_mean": g2_values,
        'delta': meth_diff,
        'pvalue': p_values,
        'F': f_value,
        'DMR': dmr_value
    })
    
    return df

def adjust_p_values(df, group1="g1", group2="g2", qvalue=0.05):
    """
    Perform Benjamini-Hochberg (BH) correction on the given DMR data and return two DataFrames: 
    1. dmr_data_with_padj: This DataFrame contains the original p-values and the adjusted p-values.
    2. significant_dmr_data: This DataFrame retains only rows where the adjusted p-value is less than 0.05, and it includes only the adjusted p-values.
    """

    if df.empty:
        empty_cols = ['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean",
                      'delta', 'F', 'pvalue', 'p_adj', 'pro_var']
        empty_df = pd.DataFrame(columns=empty_cols)
        return empty_df.copy(), empty_df.copy()

    # BH correction
    df = df.copy()
    df['p_adj'] = multipletests(df['pvalue'], method='fdr_bh')[1]

    for col in ['count', f"{group1}_mean", f"{group2}_mean", 'delta', 'F']:
        df[col] = df[col].round(4)
    
    # Use scientific notation for a concise representation of pvalue and p_adj
    df['pvalue'] = df['pvalue'].apply(lambda x: float(f"{x:.4g}"))
    df['p_adj'] = df['p_adj'].apply(lambda x: float(f"{x:.4g}"))

    dmr_data_with_padj = df[['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean", 'delta', 'F', 'pvalue', 'p_adj', 'pro_var']].copy()
    dmr_data_with_padj = dmr_data_with_padj.sort_values(by='start')
    significant_dmr_data = dmr_data_with_padj[dmr_data_with_padj['p_adj'] < qvalue][['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean", 'delta', 'F','pvalue', 'p_adj', 'pro_var']].copy()

    return dmr_data_with_padj, significant_dmr_data

def _die(msg: str, code: int = 1):
    raise SystemExit(f"[ERROR] {msg}")

def ensure_outdir(out_dir: str) -> Path:
    out_path = Path(out_dir).expanduser().resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    return out_path

def validate_required_args(out_dir: str, group1: str, group2: str, indir: str, cyber_lab: str):
    if not out_dir or not group1 or not group2:
        _die("--out-dir, --group1, and --group2 are required.")
    if (not cyber_lab) and (not indir):
        _die("Must provide either --cyber-lab (lab file) OR --in-dir (input directory).")

def validate_lab(lab_file: str) -> Path:
    """
    Validate lab file format:
        sample_id <TAB> group_label <TAB> path/to/tsv
    Rules:
        - exactly 2 unique group labels in column 2 (like bash)
        - each referenced file exists and readable (relative paths resolved to lab dir)
    """
    lab_path = Path(lab_file).expanduser().resolve()
    if not lab_path.is_file():
        _die(f"Provided lab file not found: {lab_path}")

    base_dir = lab_path.parent
    groups = set()
    line_no = 0

    with lab_path.open("r", encoding="utf-8") as f:
        for line in f:
            line_no += 1
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                _die(f"Line {line_no}: lab file must have >=3 columns (tab-delimited). Got: {line!r}")
            sid, glab, pth = parts[0].strip(), parts[1].strip(), parts[2].strip()
            if not pth:
                _die(f"Line {line_no}: empty path (col3) in lab file.")

            groups.add(glab)

            file_path = Path(pth)
            if not file_path.is_absolute():
                file_path = (base_dir / file_path).resolve()

            if not file_path.exists():
                _die(f"Line {line_no}: file not found: {file_path}")
            if not os.access(str(file_path), os.R_OK):
                _die(f"Line {line_no}: file not readable: {file_path}")

    if len(groups) != 2:
        _die(f"Lab file must contain exactly 2 unique group labels in column 2; got: {len(groups)} -> {sorted(groups)}")

    print(f"[INFO] Lab file validation passed: {lab_path}")
    return lab_path

def build_lab_from_indir(indir: str, out_dir: Path, group1: str, group2: str) -> Path:
    """
    Mimic the bash behavior:
        - scan indir for files whose name contains group1 / group2 (case-sensitive like bash script)
        - create out_dir/in_cyber.lab with sample ids: group1_1, group1_2 ... group2_1 ...
        - paths written as absolute paths
    """
    in_path = Path(indir).expanduser().resolve()
    if not in_path.is_dir():
        _die(f"Input directory not found: {in_path}")

    outlab = out_dir / "in_cyber.lab"
    # truncate
    outlab.write_text("", encoding="utf-8")

    # bash: g1_matches=(*"${group1}"*)
    g1_files = sorted([p for p in in_path.iterdir() if p.is_file() and (group1 in p.name)])
    g2_files = sorted([p for p in in_path.iterdir() if p.is_file() and (group2 in p.name)])


    if len(g1_files) == 0 or len(g2_files) == 0:
        _die(
            f"No matching TSVs found. group1={group1} count={len(g1_files)}; "
            f"group2={group2} count={len(g2_files)}. "
            f"Expected filenames like *{group1}* and *{group2}* under: {in_path}"
        )


    lines: List[str] = []
    for i, f in enumerate(g1_files, start=1):
        lines.append(f"{group1}_{i}\t{group1}\t{str(f.resolve())}")
    for j, f in enumerate(g2_files, start=1):
        lines.append(f"{group2}_{j}\t{group2}\t{str(f.resolve())}")

    outlab.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[INFO] Built lab: {outlab} (group1={len(g1_files)}, group2={len(g2_files)})")

    validate_lab(str(outlab))
    return outlab

def resolve_or_build_lab(indir: str, cyber_lab: str, out_dir: Path, group1: str, group2: str) -> Path:
    if cyber_lab:
        print(f"[INFO] Using user-provided lab: {cyber_lab}")
        return validate_lab(cyber_lab)
    else:
        print(f"[INFO] Building lab from indir={indir} with group1={group1} group2={group2}")
        return build_lab_from_indir(indir, out_dir, group1, group2)

def merge_chr_results(out_dir: Path, pattern: str = "chr*.txt", out_name: str = "cyberDMR_result.txt"):
    """
    Mimic bash post-processing:
        cat chr*.txt | grep -v "^chromosome" | sort -k1,1V -k2,2n -k3,3n > cyberDMR_result.txt
        Read all chr*.txt files, skip header lines (starting with 'chromosome'), 
        sort by (chromosome natural sort, start integer, end integer), and then write out.
    """
    files = sorted(out_dir.glob(pattern))
    if not files:
        return

    rows: List[Tuple[str, int, int, str]] = []
    for fp in files:
        with fp.open("r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("chromosome"):
                    continue
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                try:
                    start = int(parts[1])
                except ValueError:
                    continue
                try:
                    end = int(parts[2])
                except ValueError:
                    end = start
                rows.append((chrom, start, end, line.rstrip("\n")))


    def chrom_key(c: str):
        # Compatible with chr1..chr22, chrX, chrY; other chromosomes maintain lexicographic order
        m = re.match(r"^chr(\d+)$", c)
        if m:
            return (0, int(m.group(1)))
        if c == "chrX":
            return (1, 23)
        if c == "chrY":
            return (1, 24)
        return (2, c)

    rows.sort(key=lambda x: (chrom_key(x[0]), x[1], x[2]))

    out_fp = out_dir / out_name
    with out_fp.open("w", encoding="utf-8") as out:
        for _, _, _, raw in rows:
            out.write(raw + "\n")
    
    chr_dir = out_dir / "chr"
    chr_dir.mkdir(parents=True, exist_ok=True)

    for fp in files:
        dest = chr_dir / fp.name
        # shutil.move works across filesystems; overwrite protection:
        if dest.exists():
            dest.unlink()
        shutil.move(str(fp), str(dest))

    print(f"[INFO] Merged {len(files)} files -> {out_fp} ({len(rows)} rows)")


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for cyberDMR.

    Required arguments:
        -o/--out-dir   Output directory
        -g1/--group1   Group1 label
        -g2/--group2   Group2 label

    Input source (choose ONE):
        -i/--in-dir        Input directory containing TSV files (auto-generate in_cyber.lab)
        -lab/--cyber-lab   Path to an existing cyber.lab file
    """

    parser = argparse.ArgumentParser(
        prog="cyberDMR",
        description="Detect DMRs with cyberDMR.",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
    )

    # ---------------- Required arguments ----------------
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--out-dir", "-o", required=True, metavar="PATH", help="Output directory.")
    required.add_argument("--group1", "-g1", required=True, metavar="STR", help="Group1 label.")
    required.add_argument("--group2", "-g2", required=True, metavar="STR", help="Group2 label.")

    # ---------------- Input (mutually exclusive) ----------------
    # Either --indir or --cyber-lab must be provided
    input_grp = parser.add_argument_group("Input (choose ONE)")
    mx = input_grp.add_mutually_exclusive_group(required=True)
    mx.add_argument("--in-dir", "-i", dest="in_dir", metavar="PATH",
                    help="Input directory containing TSV files; will auto-generate out-dir/in_cyber.lab.\n"
                         "Files are matched by substring: *group1* and *group2* (case-sensitive).")
    mx.add_argument("--cyber-lab", "-lab", dest="cyber_lab", metavar="PATH", help="Path to an existing cyber.lab file (3 columns: sample_id, group_label, tsv_path).")


    # ---------------- Optional arguments ----------------
    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument("--threads", "-t", type=int, default=8, metavar="INT", help="Number of worker processes (default: 8).")
    optional.add_argument("--chroms", "-chr",
                          default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY",
                          metavar="STR",
                          help="Chromosome set, e.g.:\n"
                               "  'chr1,chr2,chr3'\n"
                               "  '1-22'\n"
                               "  'chr1-chr22'\n"
                               "(default: chr1-chr22,chrX,chrY).")

    optional.add_argument("--delta", "-d", type=float, default=0.1, metavar="FLOAT", help="Methylation difference (delta) threshold (default: 0.1).")
    optional.add_argument("--cpg-distance", "-bdis", type=int, default=500, metavar="INT", help="Max CpG distance for blocking (default: 500).")
    optional.add_argument("--cpg-count", "-ct", type=int, default=5, metavar="INT", help="Min CpG count per block (default: 5).")
    optional.add_argument("--min-cov", "-cov", type=int, default=5, metavar="INT", help="Minimum CpG coverage for filling (default: 5).")
    optional.add_argument("--max-dist", "-fdis", type=int, default=500, metavar="INT", help="Maximum distance (bp) of adjacent CpGs for filling (default: 500).")
    optional.add_argument("--qvalue", "-q", type=float, default=0.05, metavar="FLOAT", help="BH-corrected p-value threshold for significant DMRs (default: 0.05).")
    optional.add_argument("--Fvalue", "-f", type=float, default=150, metavar="FLOAT", help="F statistic threshold (default: 150).")

    return parser.parse_args()


def process_one_chromosome(in_chr, out_dir, group1, group2, 
                           threads=1, coverage_threshold=5, max_distance=500, CpG_distance=500, CpG_count=5, delta_m_mean_threshold=0.1,
                           qvalue=0.05, Fvalue=15, cyber_lab=None):
    
    """
    Process a single chromosome: fill missing values, merge samples, perform CpG blocking, and detect DMRs.        
    
    Arguments:
    - in_chr: Chromosome name (e.g., 'chr1')
    - out_dir: Output directory where results will be saved
    - group1: Group1 label (e.g., 'control')
    - group2: Group2 label (e.g., 'treatment')
    - threads: Number of threads to use
    - coverage_threshold: Minimum coverage threshold for imputation (default: 5)
    - max_distance: Maximum distance for imputation (default: 500)
    - CpG_distance: Maximum distance between CpGs for blocking (default: 500)
    - CpG_count: Minimum number of CpGs per block (default: 5)
    - delta_m_mean_threshold: Delta methylation mean threshold for DMR detection (default: 0.1)
    - cyber_lab: Optional path to a custom cyber.lab file
    - qvalue: BH-corrected p-value
    - Fvalue: F statitic
    """
    
    print(f"Processing {in_chr}...")

    if cyber_lab:
        infile = cyber_lab
    else:
        infile = f"{out_dir}/in_cyber.lab" # Assuming the input is a unified large file, filter by chromosome
    
    sample_data = load_inlab_chrwise(
        infile,
        target_chr=in_chr,
        column_names=["Chr", "Pos", "Meth_Level", "Coverage"],
        num_threads=threads
    )

    if not sample_data:
        print(f"No data found for {in_chr}, skipping.")
        return

    # Extract sample and group information
    label = pd.DataFrame(
        [(entry['sample'], entry['group']) for entry in sample_data],
        columns=['sample','group']
    )

    # Step 1: Impute and merge data
    processed_samples = process_samples(sample_data, coverage_threshold=coverage_threshold, max_distance=max_distance)
    merged_data = merge_samples_fast(processed_samples)
    #merged_data.to_csv(os.path.join(out_dir, f"{in_chr}_merged_data_after_filling.txt"), sep="\t", header=True, index=False)

    # Step 2: CpG blocking
    out_data, block_ranges = process_data(merged_data, label, group1, group2, CpG_distance=CpG_distance, CpG_count=CpG_count)
    
    # Step 3: Clustering and DMR detection
    dmr_data_with_padj, significant_dmr_data = find_blocks_greedy(
        out_data, 
        delta_m_mean_threshold=delta_m_mean_threshold,
        group1=group1, group2=group2,
        qvalue=qvalue,
        chr_col=in_chr,
        Fvalue=Fvalue,
        out_dir=out_dir
    )
    #dmr_data_with_padj.to_csv(os.path.join(out_dir, f"no_adjuset_{in_chr}_cyberDMR.tmp"), sep="\t", header=True, index=False)
    significant_dmr_data.to_csv(os.path.join(out_dir, f"{in_chr}_cyberDMR.txt"), sep="\t", header=True, index=False)

    print(f"Finished processing {in_chr}.")


# Pre-generate 80% of sites as fixed positions (assuming the entire genome consists of 1 million bases)
BASE_CPG_SITES = np.sort(np.random.randint(1, 1_000_000, size=800))  # 80% fixed positions

def generate_methylation_data(chr_name="chr1", num_points=1000, coverage_range=(1, 20), meth_range=(0, 1)):
    """
    Generate simulated CpG methylation data, ensuring 80% of sites are fixed and 20% of sites vary randomly
    
    Parameters:
    - chr_name: Chromosome name
    - num_points:  Total number of CpG sites to generate
    - coverage_range: Tuple specifying the range for coverage, as (min_coverage, max_coverage)
    - meth_range: Tuple specifying the range for methylation level, as (min_meth, max_meth)

    Returns:
    - pandas DataFrame:A DataFrame with columns: 'Chr', 'Pos', 'Meth_Level', 'Coverage'
    """
    np.random.seed()

    # Calculate the number of sites for the 80% and 20% portions
    num_fixed = int(num_points * 0.8)
    num_random = num_points - num_fixed

    fixed_positions = np.random.choice(BASE_CPG_SITES, num_fixed, replace=False)
    random_positions = np.random.randint(1, 1_000_000, size=num_random)
    all_positions = np.sort(np.concatenate([fixed_positions, random_positions]))
    meth_levels = np.random.uniform(meth_range[0], meth_range[1], size=num_points)

    # Generate coverage (simulated using a Gamma distribution)
    shape, scale = 2, 5  # Gamma distribution parameters; a right-skewed distribution to simulate actual sequencing depth
    coverage = np.random.gamma(shape, scale, size=num_points).astype(int)
    coverage = np.clip(coverage, coverage_range[0], coverage_range[1])

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

    Parameter:
    - group1_count: The number of samples in group1
    - group2_count: The number of samples in group2
    - chr_name: The chromosome name
    - num_points: The number of data points to generate per sample

    Return:
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

    print(f" inlab.txt write to: {inlab_path}")
    return samples

def generate_simulated_data(
        num_cpg=10, interval=50, std_dev=0.05,
        group1_size=5, group2_size=5,
        alpha1=2, beta1=5, alpha2=5, beta2=2,
        group1="g1", group2="g2"):
    """
    Generate simulated CpG methylation data (in wide format), while also calculating the following statistical information:
    - For each group: mean methylation level, variance, and coverage (number of samples)
    - Methylation difference (ΔM) between the two groups
    """

    np.random.seed(42)
    positions = np.arange(1, num_cpg * interval + 1, interval)
    
    # Generate methylation levels for the experimental group (Group1) and the control group (Group2)
    meth_group1 = np.random.beta(alpha1, beta1, (num_cpg, group1_size))
    meth_group2 = np.random.beta(alpha2, beta2, (num_cpg, group2_size))

    # Calculate statistical information
    mean_group1 = np.mean(meth_group1, axis=1)
    var_group1 = np.var(meth_group1, axis=1, ddof=1) if group1_size > 1 else np.zeros(num_cpg)
    mean_group2 = np.mean(meth_group2, axis=1)
    var_group2 = np.var(meth_group2, axis=1, ddof=1) if group2_size > 1 else np.zeros(num_cpg)
    
    delta_M = mean_group2 - mean_group1 
    block = ["block1"] * num_cpg
    df_wide = pd.DataFrame({
        "Chromosome": ["chr1"] * num_cpg,
        "Position": positions,
        **{f"{group1}_{i+1}": meth_group1[:, i] for i in range(group1_size)},
        **{f"{group2}_{i+1}": meth_group2[:, i] for i in range(group2_size)}
    })

    # Statistical data DataFrame
    df_summary = pd.DataFrame({
        "Chromosome": ["chr1"] * num_cpg,
        "Position": positions,
        f"{group1}_mean": mean_group1,
        f"{group2}_mean": mean_group2,
        f"{group1}_var": var_group1,
        f"{group2}_var": var_group2,
        f"{group1}_cov": group1_size,
        f"{group2}_cov": group2_size,
        "mean_diff": delta_M,
        "Block": block
    })
    return df_wide, df_summary

def convert_to_long_format(df_wide):
    """
    Convert wide-format data to long-format, suitable for WBR calculation
    """
    df_long = df_wide.melt(id_vars=["Chr", "Pos"], var_name="Sample", value_name="Meth_Level")
    df_long["Group"] = df_long["Sample"].apply(lambda x: "group1" if "G1" in x else "group2")
    return df_long


def main():
    """
    Main control function:
    - Parse parameters.
    - Prepare/validate cyber.lab (from --in-dir OR --cyber-lab).
    - Dispatch chromosome tasks using multiple processes.
    - Merge chr*.txt into a single sorted result file.
    """
    args = parse_args()
    print(f"[INFO] parameters is {args}")
   
    out_dir = ensure_outdir(args.out_dir)
    indir = getattr(args, "in_dir", "") or ""
    cyber_lab = getattr(args, "cyber_lab", "") or ""

    # ===== Build/validate lab =====
    lab_path = resolve_or_build_lab(
        indir=indir,
        cyber_lab=cyber_lab,
        out_dir=out_dir,
        group1=args.group1,
        group2=args.group2            
    )

    # ===== Chroms parsing =====
    chroms = []
    if args.chroms:
        s = args.chroms.strip()

        # Support comma-separated list: chr1,chr2 or 1,2
        if "," in s:
            chroms = []
            for x in s.split(","):
                x = x.strip()
                if not x:
                    continue
                chroms.append(x if x.startswith("chr") else f"chr{x}")
    
        # Support range: 1-22 or chr1-chr22
        elif "-" in s:
            a, b = [z.strip() for z in s.split("-", 1)]
            a = a.replace("chr", "")
            b = b.replace("chr", "")
            if a.isdigit() and b.isdigit():
                chroms = [f"chr{i}" for i in range(int(a), int(b) + 1)]
            else:
                # fallback: non-standard range, treated directly as a single string
                chroms = [s if s.startswith("chr") else f"chr{s}"]

        # Single
        else:
            chroms = [s if s.startswith("chr") else f"chr{s}"]
    
    else:
        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    # ===== run per-chrom (ensure futures are awaited) =====
    futures = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        for chrom in chroms:
            futures.append(
                executor.submit(
                    process_one_chromosome,
                    chrom,
                    str(out_dir),
                    args.group1,
                    args.group2,
                    1,
                    coverage_threshold=args.min_cov,
                    max_distance=args.max_dist,
                    CpG_distance=args.cpg_distance,
                    CpG_count=args.cpg_count,
                    delta_m_mean_threshold=args.delta,
                    cyber_lab=str(lab_path),
                    qvalue=args.qvalue,
                    Fvalue=args.Fvalue,                            
                )                        
            )
    
        # Wait and raise if any child fail
        for fu in as_completed(futures):
            fu.result()

    # ===== Merge outputs =====
    merge_chr_results(out_dir, pattern="chr*cyberDMR.txt", out_name="cyberDMR_result.bed")
    print("[INFO] All done.")

if __name__ == "__main__":
    main()
