import os
import pandas as pd
import numpy as np
from numba import njit
from sklearn.impute import KNNImputer
from concurrent.futures import ThreadPoolExecutor
from typing import List, Dict

# Pre-generate 80% of sites as fixed positions (assuming the entire genome consists of 1 million bases)
BASE_CPG_SITES = np.sort(np.random.randint(1, 1_000_000, size=800))  # 80% fixed positions

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
    np.random.seed()

    # Calculate the number of sites for the 80% and 20% portions
    num_fixed = int(num_points * 0.8)
    num_random = num_points - num_fixed

    # 80% fixed sites (selected from the pre-generated BASE_CPG_SITES)
    fixed_positions = np.random.choice(BASE_CPG_SITES, num_fixed, replace=False)

    # 20% randomly generated sites
    random_positions = np.random.randint(1, 1_000_000, size=num_random)

    # Merge and sort
    all_positions = np.sort(np.concatenate([fixed_positions, random_positions]))

    # Generate random methylation levels
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
    - sample_data: List[Dict]，A list of dictionaries, where each dictionary represents a sample and contains the keys 'sample', 'group', and 'data'
    - coverage_threshold:  The threshold for defining low coverage
    - max_distance: The maximum allowable distance for considering left and right neighbors during imputation

    Returns:
    - List[Dict]，A list of dictionaries with the same structure as the input sample_data, but where the value for the 'data' key in each dictionary is the DataFrame after imputation
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
    Concurrently processes data for multiple samples using multiple threads. This method is suitable for cases with a large number of samples that can be processed independently

    Parameters:
    - sample_data: List[Dict]，A list of dictionaries, where each dictionary contains information for a sample, including keys like 'sample', 'group', and 'data'
    - coverage_threshold: The coverage threshold
    - max_distance:  The maximum distance to consider for neighbors during imputation
    - num_threads:  The number of concurrent threads to use (default is 4)

    Returns:
    - List[Dict]，A list of dictionaries containing the processed data for each sample. The order of items in the list is maintained consistently with the input sample_data
    """
    def process_one(entry):
        sample = entry['sample']
        group = entry['group']
        df = entry['data']
        try:
            processed_df = fill_low_coverage_cpg_numba(df, coverage_threshold, max_distance)
            print(f"Success: {sample} has completed")
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

if __name__ == "__main__":
    output_folder = "./test"
    os.makedirs(output_folder, exist_ok=True)  
    num_sample = 10
    chr_name = "chr1"
    num_points = 1000
    samples = simulate_multiple_samples(output_folder=output_folder, group1_count=5, group2_count=5)
    processed_samples = process_samples(samples, coverage_threshold=5, max_distance=500)
    merged_df = merge_samples_fast(processed_samples)
    merged_df.to_csv(os.path.join(output_folder,"merge.txt"), sep='\t',index=False)
