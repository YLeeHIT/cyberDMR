import numpy as np
import pandas as pd
import random
import argparse
import shutil
import os
from datetime import datetime
from scipy.stats import beta
from datetime import datetime
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Simulated DMR Data")
    parser.add_argument(
        "--output_dir",
        type=str,
        default="./",
        help="outdir (default: ./"
    )
    return parser.parse_args()

# Function to simulate DMR (Differentially Methylated Region) regions
def simulate_dmr_region_with_input_limit(chr_name, region_start, region_end, delta_methylation, max_cpgs=50, density="auto"):
    """
    Simulates a Differentially Methylated Region (DMR)

    Parameters:
    - chr_name: Chromosome name, e.g., "chr1"
    - region_start: The start position of the region (integer)
    - region_end: The end position of the region (must be greater than region_start + 100)
    - delta_methylation: The mean methylation level difference between the experimental and control groups (can be positive or negative)
    - max_cpgs: Maximum number of CpGs within the region. If this limit is exceeded, generation may terminate early (default is 50)
    Returns:
     A dictionary containing information about the simulated DMR, including region details, CpG site positions, methylation values for the two groups, and the mean methylation difference
    """

    # Define the minimum and maximum interval (unit: bp) between every two CpGs
    if density == "dense":
        min_distance = 2
        max_distance = 50
    elif density == "sparse":
        min_distance = 50
        max_distance = 499
    else:
        min_distance = 2
        max_distance = 499

    # Region length validation
    if region_end - region_start < 100:
        return None

    # Generate CpG site positions
    cpg_sites = []
    pos = region_start
    while pos < region_end:
        if len(cpg_sites) >= max_cpgs:
            break  # If the maximum number limit is exceeded, terminate early

        if len(cpg_sites) > 0:
            max_gap = min(max_distance, region_end - pos)
            if max_gap < min_distance:
                break  # If the maximum number limit is exceeded, terminate early
            gap = random.randint(min_distance, max_gap)
        else:
            gap = 0  # The first CpG uses the starting position
        pos += gap
        if pos < region_end:
            cpg_sites.append(pos)

    if len(cpg_sites) == 0:
        return None

    # If terminated early due to exceeding the maximum number, update the end position
    region_end_final = min(region_end, cpg_sites[-1])
    region_length = region_end_final - region_start

    # Check if three preconditions are met
    if len(cpg_sites) < 5 or region_length < 100 or abs(delta_methylation) < 0.1:
        return None

    # Set the mean methylation level (simulated using a beta distribution)
    mean_control = np.random.uniform(0.2, 0.8)
    delta_methylation_modify = round(np.random.normal(loc=delta_methylation, scale=0.01), 3)
    mean_treatment = mean_control + delta_methylation_modify

    # Prevent mean_treatment from exceeding the threshold
    mean_treatment = min(max(mean_treatment, 0.001), 0.999)

    # Set the shape parameters for the beta distribution (fixed precision)
    precision = 20
    alpha_control = mean_control * precision
    beta_control = (1 - mean_control) * precision
    alpha_treatment = mean_treatment * precision
    beta_treatment = (1 - mean_treatment) * precision

    # Generate a sequence of methylation values (rounded to three decimal places)
    control_vals = beta.rvs(alpha_control, beta_control, size=len(cpg_sites)).round(3)
    treatment_vals = beta.rvs(alpha_treatment, beta_treatment, size=len(cpg_sites)).round(3)

    # Perform a consistency check
    delta_vals = treatment_vals - control_vals
    if (delta_methylation > 0 and np.any(delta_vals < 0)) or (delta_methylation < 0 and np.any(delta_vals > 0)):
        return None

    # Actually calculate the mean methylation difference
    actual_delta = round(np.mean(treatment_vals) - np.mean(control_vals), 3)

    # Return the simulation result
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

# Define a simulation function for non-DMR regions
def simulate_nondmr_region(chr_name, region_start, region_end, delta_methylation=None, max_cpgs=50, density="auto"):
    """
    Simulates a non-DMR (Non-Differentially Methylated Region)

    Conditions:
    - CpG count is less than 5 or delta_methylation (the difference in methylation levels) is less than 0.1

    Parameters:
    - chr_name: Chromosome name
    - region_start: Start position of the regio
    - region_end: End position of the region
    - delta_methylation: (optional) The methylation difference between experimental and control groups. If not provided, a value less than 0.1 will be automatically generated to ensure non-DMR status.
    - max_cpgs: Maximum number of CpGs to generate within the region (this can be used to control the CpG count, e.g., to keep it below 5 for non-DMR simulation)

    Returns:
    A dictionary containing information about the simulated non-DMR 
    """

    # Define the minimum and maximum interval (unit: bp) between every two CpGs
    if density == "dense":
        min_distance = 2
        max_distance = 50
    elif density =="sparse":
        min_distance = 50
        max_distance = 499
    else:
        min_distance = 2
        max_distance = 499

    # If the difference value is not specified, randomly generate one in the range [-0.05, 0.05]
    if delta_methylation is None:
        delta_methylation = np.round(np.random.uniform(-0.02, 0.02), 3)

    else:
        delta_methylation = np.round(np.random.uniform(-abs(delta_methylation), abs(delta_methylation)), 3)
    # Determine the CpG count: it can be less than 5, or greater than 5 if the difference value is less than 0.1
    allow_large_cpg = np.random.choice([True, False])

    cpg_sites = []
    pos = region_start
    while pos < region_end:
        if len(cpg_sites) >= max_cpgs:
            break
        if not allow_large_cpg and len(cpg_sites) >= 4:
            break  # Control the CpG count to be less than 5

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
    
    # Prevent mean_treatment from exceeding the threshold.
    mean_treatment = min(max(mean_treatment, 0.001), 0.999)

    precision = 20
    alpha_control = mean_control * precision
    beta_control = (1 - mean_control) * precision
    alpha_treatment = mean_treatment * precision
    beta_treatment = (1 - mean_treatment) * precision

    control_vals = beta.rvs(alpha_control, beta_control, size=len(cpg_sites)).round(3)
    treatment_vals = beta.rvs(alpha_treatment, beta_treatment, size=len(cpg_sites)).round(3)

    actual_delta = round(np.mean(treatment_vals) - np.mean(control_vals), 3)

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

# Simulate inconsistent DMRs
def simulate_dmr_with_inconsistent_points(chr_name, region_start, region_end, delta_methylation, max_cpgs=50, density="auto"):
    """
    Simulates an inconsistent DMR (Differentially Methylated Region) with maximal perturbation
    A point with a reversed methylation difference direction is inserted every 4 CpGs. This is done regardless of the original direction at these specific CpG sites, thereby forcibly disrupting overall directional consistency.

    Returns:
    The DMR structure after perturbation with these inconsistent points. (Note: The original direction at the insertion points is not considered when the 'opposite direction' is enforced.)
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

    # Insert a perturbation point every 4 points (e.g., at indices/positions 3, 7, 11, ...)
    indices = list(range(3, len(dmr["CpG_sites"]), 4))

    for idx in indices:
        # Swap control and treatment values to create a 'direction reversal'
        dmr["methylation_control"][idx], dmr["methylation_treatment"][idx] = dmr["methylation_treatment"][idx], dmr["methylation_control"][idx]

    # Update the overall mean difference
    control_arr = np.array(dmr["methylation_control"])
    treatment_arr = np.array(dmr["methylation_treatment"])
    dmr["mean_delta_methylation"] = round(np.mean(treatment_arr - control_arr), 3)

    return dmr

# Simulate inconsistent DMRs, which can then be further split into sub-DMRs based on this
def simulate_dmr_with_subdmr_points(chr_name, region_start, region_end, delta_methylation, max_cpgs=50):
    """
    Simulates a DMR (Differentially Methylated Region) that includes inconsistent points

    Parameters are the same as those of the base/original DMR simulation function, with the following additional behavior:
    - Within the generated DMR region, 1-3 inconsistent points are randomly inserted. At these points, the direction of the methylation difference between 'control' and 'treatment' is opposite to that of the overall DMR

    Returns:
    A dictionary with the same structure as that returned by the standard/original DMR simulation function
    """
    
    # Check if the inconsistent DMRs can be split into sub-DMRs
    def check_split_subdmrs(dmr, flip_indices):

        n_cpgs = len(dmr["CpG_sites"])
        if n_cpgs == 0:
            return None

        # Sort and construct cutting points (or split points)
        split_points = sorted(set(flip_indices))
        segments = []
        start = 0
        for cut in split_points:
            if start < cut:
                segments.append(range(start, cut))  # Does not include cut points
            start = cut + 1
        if start < n_cpgs:
            segments.append(range(start, n_cpgs))

        sub_dmrs = []

        for segment in segments:
            if len(segment) < 5:
                continue  # Insufficient CpG count

            idx_list = list(segment)
            start_pos = dmr["CpG_sites"][idx_list[0]]
            end_pos = dmr["CpG_sites"][idx_list[-1]]
            if end_pos - start_pos < 50:
                continue  # Insufficient length

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

    # Call the original DMR function
    dmr = simulate_dmr_region_with_input_limit(
        chr_name=chr_name,
        region_start=region_start,
        region_end=region_end,
        delta_methylation=delta_methylation,
        max_cpgs=max_cpgs
    )

    if dmr is None:
        return None  # If simulation fails, return directly

    # Determine the direction of the difference: whether control is greater or treatment is greater
    direction = np.sign(dmr["mean_delta_methylation"])

    # If the direction is 0 (i.e., difference is negligible), skip insertion
    if direction == 0:
        return None

    # Decide the number of inconsistent points to insert (1 to 3)
    n_flip = min(len(dmr["CpG_sites"]), random.randint(1, 3))

    # Randomly select indices to insert inconsistent points
    indices = random.sample(range(len(dmr["CpG_sites"])), k=n_flip)

    # Swap control and treatment values for the selected points, making their 'direction opposite'
    for idx in indices:
        c_val = dmr["methylation_control"][idx]
        t_val = dmr["methylation_treatment"][idx]
        if (direction > 0 and c_val > t_val) or (direction < 0 and c_val < t_val):
            # Swap only if they originally have the same direction, to make them 'inconsistent'
            dmr["methylation_control"][idx], dmr["methylation_treatment"][idx] = t_val, c_val

    # Update the overall mean difference (for easier validation)
    control_arr = np.array(dmr["methylation_control"])
    treatment_arr = np.array(dmr["methylation_treatment"])
    dmr["mean_delta_methylation"] = round(np.mean(treatment_arr - control_arr), 3)

    # Check if sub-DMRs are formed
    return check_split_subdmrs(dmr, indices)

# Utility function: Derive alpha and beta parameters for a beta distribution from mean and std (standard deviation)
def beta_params_from_mean_std(mean, std):
    mean = np.clip(mean, 1e-3, 1 - 1e-3)
    var = std ** 2
    common = mean * (1 - mean) / var - 1
    alpha = mean * common
    beta_ = (1 - mean) * common
    return max(alpha, 1e-3), max(beta_, 1e-3)

# Sample expansion function with missing value control
def simulate_group_samples_with_missing(
    dmr,
    n_control: int,
    n_treatment: int,
    coverage_mean: int,
    coverage_std: int,
    group_std: float,
    output_dir: str = "./"
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
                
                # Use beta distribution
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

    for sample_key, records in all_samples.items():
        df = pd.DataFrame(records)
        df.to_csv(os.path.join(output_dir, f"{sample_key}.tsv"), sep="\t", index=False)

    return {
        "output_dir": output_dir,
        "dmr_missing_rate": round(dmr_missing_rate, 3),
        "sample_missing_rate": sample_missing_rate
    }

# Updated validation function, adding logic to assess mean values of sliding window sub-regions
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

    # Check the number of consecutive out-of-bounds/limit-exceeding occurrences
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

    # Check the mean difference of sub-regions using a sliding window
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
    output_dir: str = "./",
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

# Randomly shuffle generation order + correct classification output + insert labels by category
def simulate_mixed_regions_randomized(
    total_dmr: int,
    mean_delta: float,
    n_control: int,
    n_treatment: int,
    coverage_mean: int,
    coverage_std: int,
    output_dir: str = "./",
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

    # Calculate the number of regions for each type
    dmr_count = int(total_dmr * dmr_per)
    dmr_notable_count = int(total_dmr * dmr_notable_per)
    dmr_inconsistent_count = int(total_dmr * dmr_inconsis_per)
    dmr_inconsistent_subdmr_count = int(total_dmr * dmr_sub_per)
    nondmr_count = total_dmr - dmr_count - dmr_inconsistent_count - dmr_notable_count - dmr_inconsistent_subdmr_count

    # Construct a task list (for labeling/tagging)
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
        
        # Generate some negative deltas: 80% probability to keep the original value (positive), 20% probability to negate the original value (making it negative)
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
            if isinstance(region, list):  # Indicate that it is a sub-DMR scenario
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
        
     # Organize the output
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
    density=density, 
    dense_ratio=dense_ratio,
    min_gap=500,
    max_gap=5000)
    return df

def main():
    parser = argparse.ArgumentParser(description="Simulate DMR regions")

    parser.add_argument("--total_dmr", type=int, default=100)
    parser.add_argument("--mean_delta", type=float, default=0.3)
    parser.add_argument("--n_control", type=int, default=5)
    parser.add_argument("--n_treatment", type=int, default=5)
    parser.add_argument("--coverage_mean", type=int, default=30)
    parser.add_argument("--coverage_std", type=int, default=5)
    parser.add_argument("--output_dir", type=str, default="./")
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

    # Call the simulation function
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

if __name__ == "__main__":
    main()
