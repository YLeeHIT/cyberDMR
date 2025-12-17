import numpy as np
import pandas as pd
import random
import argparse
import os
from math import log
from datetime import datetime
from scipy.stats import beta
from datetime import datetime
from collections import defaultdict
from scipy.stats import truncnorm
from scipy.stats import beta as beta_dist
from typing import Optional, Dict, Any

def _draw_gap(mode):
    """
    Select a CpG spacing (bp). Different patterns are assigned to different intervals and distributions
    Use lognormal distribution: gap = round(lognormal(mu, sigma)) and clip to [min, max]
    """
    # Empirical intervals and distribution parameters for various models
    cfg = {
        "island":     {"min": 3,   "max": 15,  "mu": log(8),   "sigma": 0.45},
        "shore":      {"min": 15,  "max": 50,  "mu": log(30),  "sigma": 0.5},
        "shelf":      {"min": 50,  "max": 100,  "mu": log(100),  "sigma": 0.6},
        "open_sea":   {"min": 100,  "max": 500, "mu": log(150), "sigma": 0.7},
    }
    if mode not in cfg:
        mode = "moderate"
    p = cfg[mode]
    val = int(round(np.random.lognormal(mean=p["mu"], sigma=p["sigma"])))
    return max(p["min"], min(p["max"], val))

def generate_cpg_sites(region_start, region_end, mode="moderate", max_cpgs=50):
    if region_end - region_start < 2:
        return []

    modes = [mode]
    if mode in ("moderate", "dense", "sparse"):
        if mode == "dense":
            modes = ["island", "shore"]
        elif mode == "sparse":
            modes = ["shelf", "open_sea"]
        else:
            modes = ["shore", "shelf"]

    L = region_end - region_start
    splits = np.linspace(region_start, region_end, num=len(modes)+1, dtype=int)

    cpg_sites = []
    for i, m in enumerate(modes):
        s, e = splits[i], splits[i+1]
        pos = s
        first = True
        while pos < e:
            if len(cpg_sites) >= max_cpgs:
                return cpg_sites
            gap = 0 if first else _draw_gap(m)
            first = False
            pos = pos + gap
            if pos < e:
                cpg_sites.append(pos)

    return cpg_sites[:max_cpgs]


def add_flanking_cpgs(
    cpg_sites,
    control_vals,
    treatment_vals,
    region_start,
    region_end_final,
    direction,
    alpha_control,
    beta_control,
    flank_bp=500,
    flank_max_points=10,
    n_opposite=2,
    flank_delta_cap=0.1,
    rng=None
):
    """
    Add a small number of CpG points in the flanking regions (±flank_bp) with small deltas.
    - Total number of new CpGs across both sides: 0..flank_max_points (uniform)
    - |delta| < flank_delta_cap
    - The first 'n_opposite' flanks are forced to have delta opposite to the main 'direction'
    - Remaining flanks have random sign
    """

    if rng is None:
        rng = np.random.default_rng()
            
    # Decide how many flank CpGs to insert (total across both sides)
    n_flank = int(rng.integers(2, flank_max_points + 1))
    if n_flank == 0:
        return cpg_sites, control_vals, treatment_vals, 0
                                
    # Split left/right counts
    n_left = n_flank // 2
    n_right = n_flank - n_left
                                            
    # Sample flank positions
    left_positions = []
    if n_left > 0:
        left_positions = rng.integers(region_start - flank_bp, region_start, size=n_left).tolist()
                                                                    
    right_positions = []
    if n_right > 0:
        right_positions = rng.integers(region_end_final + 1, region_end_final + flank_bp + 1, size=n_right).tolist()

    flank_positions = sorted(left_positions) + sorted(right_positions)

    # Build delta signs for flanks
    # First n_opposite are opposite to 'direction'; rest random
    n_opposite = min(n_opposite, n_flank)
    opp_sign = -np.sign(direction) if direction != 0 else -1
    flank_signs = [opp_sign] * n_opposite
    if n_flank > n_opposite:
        flank_signs += rng.choice([+1, -1], size=(n_flank - n_opposite)).tolist()
    flank_signs = np.array(flank_signs, dtype=int)

    # Small magnitudes: (0, flank_delta_cap], avoid exactly 0
    low_mag = min(0.01, max(1e-4, 0.1 * flank_delta_cap))
    flank_mags = rng.uniform(low_mag, flank_delta_cap-0.01, size=n_flank)
                                                    
    # Sample control values for flanks from the same Beta distribution (stylistic consistency)
    flank_control = beta_dist.rvs(alpha_control, beta_control, size=n_flank, random_state=rng)
    flank_control = np.clip(np.round(flank_control, 3), 0.001, 0.999)
                                                                
    # Apply small deltas with chosen signs
    flank_delta_vals = flank_signs * flank_mags
    flank_treatment = flank_control + flank_delta_vals
    flank_treatment = np.clip(np.round(flank_treatment, 3), 0.001, 0.999)
                                                                                
    # Merge with main arrays and sort by genomic coordinate
    all_positions = np.array(cpg_sites + flank_positions, dtype=int)
    all_control = np.concatenate([control_vals, flank_control])
    all_treatment = np.concatenate([treatment_vals, flank_treatment])
                                                                                                
    order = np.argsort(all_positions)
    new_cpg_sites = all_positions[order].tolist()
    new_control_vals = all_control[order]
    new_treatment_vals = all_treatment[order]

    return new_cpg_sites, new_control_vals, new_treatment_vals, n_flank

def add_trailing_flanking_cpgs(
    cpg_sites,
    control_vals,
    treatment_vals,
    region_start,
    region_end_final,
    direction,
    alpha_control,
    beta_control,
    window_bp=500,
    min_points=5,
    max_points=10,
    n_opposite=2,
    delta_cap=0.099,
    mode = "downstream",
    rng=None,
):
    """
    Add 3–10 CpGs ONLY within the trailing (end) window of the region (last 500 bp).
    - Positions are inside [region_end_final - window_bp + 1, region_end_final].
    - First `n_opposite` added CpGs have Δ opposite to `direction`; the rest random sign.
    - |Δ| for all added CpGs < delta_cap.
    - Ensures no duplicate positions with existing cpg_sites.
    """
    if rng is None:
        rng = np.random.default_rng()

    # Decide how many to add
    n_add = int(rng.integers(min_points, max_points + 1))

    # Trailing window inside the region
    if mode == "inside":
        left = max(region_start, region_end_final - window_bp + 1)
        right = region_end_final
    elif mode == "downstream":
        left = region_end_final + 2
        right = left + window_bp
    else:
        raise ValueError("inject mode must be inside or downstream")
    
    if left > right:
        return cpg_sites, control_vals, treatment_vals, 0

    # Candidate positions in trailing window, exclude existing CpGs
    existing = set(int(p) for p in cpg_sites)
    candidates = np.array([p for p in range(left, right, 2) if p not in existing], dtype=int)
    if candidates.size == 0:
        return cpg_sites, control_vals, treatment_vals, 0

    # If not enough unique positions, cap n_add
    n_add = min(n_add, candidates.size)
    if n_add <= 0:
        return cpg_sites, control_vals, treatment_vals, 0

    # Sample unique positions
    new_pos = rng.choice(candidates, size=n_add, replace=False)
    new_pos.sort()

    # Build signs: first n_opposite opposite to main direction, rest random
    n_op = min(n_opposite, n_add)
    opp_sign = -np.sign(direction) if direction != 0 else -1
    signs = np.empty(n_add, dtype=int)
    signs[:n_op] = opp_sign
    if n_add > n_op:
        signs[n_op:] = rng.choice([+1, -1], size=n_add - n_op)
        signs[-2] = -signs[-1]

    # Small magnitudes in (low, delta_cap)
    delta_cap = float(min(delta_cap, 0.1))
    low_mag = min(0.01, max(1e-4, 0.5 * delta_cap))
    mags = np.clip(rng.uniform(low_mag, delta_cap, size=n_add),1e-4,0.09)

    # Sample control values for the new CpGs (consistent with main block)
    new_ctrl = beta_dist.rvs(alpha_control, beta_control, size=n_add, random_state=rng)
    new_ctrl = np.clip(np.round(new_ctrl, 3), 0.001, 0.999)

    # Apply small deltas
    new_treat = new_ctrl + signs * mags
    new_treat = np.clip(np.round(new_treat, 3), 0.001, 0.999)

    new_delta = new_treat - new_ctrl

    # Merge and sort by genomic position
    all_pos = np.concatenate([np.asarray(cpg_sites, dtype=int), new_pos])
    all_ctrl = np.concatenate([np.asarray(control_vals), new_ctrl])
    all_treat = np.concatenate([np.asarray(treatment_vals), new_treat])

    order = np.argsort(all_pos)
    cpg_sites_new = all_pos[order].tolist()
    control_vals_new = all_ctrl[order]
    treatment_vals_new = all_treat[order]

    return cpg_sites_new, control_vals_new, treatment_vals_new, int(n_add)


def add_trailing_flanking_cpgs_v2(
    cpg_sites,
    control_vals,
    treatment_vals,
    region_start,
    region_end_final,
    direction,
    alpha_control,
    beta_control,
    mean_val,
    window_bp=500,
    min_points=5,
    max_points=10,
    n_opposite=2,
    delta_cap=0.099,
    max_dev = 0.15,
    rand_band = 0.1,
    mode = "downstream",
    rng=None,
):
    """
    Add 3–10 CpGs ONLY within the trailing (end) window of the region (last 500 bp).
    - Positions are inside [region_end_final - window_bp + 1, region_end_final].
    - First `n_opposite` added CpGs have Δ opposite to `direction`; the rest random sign.
    - |Δ| for all added CpGs < delta_cap.
    - Ensures no duplicate positions with existing cpg_sites.
    """
    if rng is None:
        rng = np.random.default_rng()

    # Decide how many to add
    n_add = int(rng.integers(min_points, max_points + 1))

    # Trailing window inside the region
    if mode == "inside":
        left = max(region_start, region_end_final - window_bp + 1)
        right = region_end_final
    elif mode == "downstream":
        left = region_end_final + 2
        right = left + window_bp
    else:
        raise ValueError("inject mode must be inside or downstream")
    
    if left > right:
        return cpg_sites, control_vals, treatment_vals, 0

    # Candidate positions in trailing window, exclude existing CpGs
    existing = set(int(p) for p in cpg_sites)
    candidates = np.array([p for p in range(left, right, 2) if p not in existing], dtype=int)
    if candidates.size == 0:
        return cpg_sites, control_vals, treatment_vals, 0

    # If not enough unique positions, cap n_add
    n_add = min(n_add, candidates.size)
    if n_add <= 0:
        return cpg_sites, control_vals, treatment_vals, 0

    # Sample unique positions
    new_pos = rng.choice(candidates, size=n_add, replace=False)
    new_pos.sort()

    # Build signs: first n_opposite opposite to main direction, rest random
    n_op = min(n_opposite, n_add)
    opp_sign = -np.sign(direction) if direction != 0 else -1
    signs = np.empty(n_add, dtype=int)
    signs[:n_op] = opp_sign
    if n_add > n_op:
        signs[n_op:] = rng.choice([+1, -1], size=n_add - n_op)
        signs[-2] = -signs[-1]

    # Small magnitudes in (low, delta_cap)
    #delta_cap = float(min(delta_cap, 0.1))
    #low_mag = min(0.01, max(1e-4, 0.5 * delta_cap))
    mags = np.clip(rng.uniform(0.5*delta_cap, delta_cap, size=n_add),1e-4,0.09)

    # Sample control values for the new CpGs (consistent with main block)
        
    new_ctrl = beta_dist.rvs(alpha_control, beta_control, size=n_add, random_state=rng)
    new_ctrl = np.clip(np.round(new_ctrl, 3), 0.001, 0.999)

    new_ctrl = np.asarray(new_ctrl, dtype=float)
    n = len(new_ctrl)

    is_outlier = np.abs(new_ctrl - mean_val) > max_dev

    low = max(0.0, mean_val - rand_band)
    high = min(1.0, mean_val + rand_band)

    for i in range(n):
        if not is_outlier[i]:
            continue
        neighbors = []
    
        if i - 1 >= 0 and not is_outlier[i - 1]:
            neighbors.append(new_ctrl[i - 1])
        
        if i + 1 < n and not is_outlier[i + 1]:
            neighbors.append(new_ctrl[i + 1])

        if len(neighbors) == 2:
            new_ctrl[i] = 0.5 * (neighbors[0] + neighbors[1])
        else:
            new_ctrl[i] = np.random.uniform(low, high)

    # Apply small deltas
    new_ctrl = np.clip(np.round(new_ctrl, 3), 0.1, 0.9)    
    new_treat = new_ctrl + signs * mags
    new_treat = np.clip(np.round(new_treat, 3), 0.001, 0.999)

    new_delta = new_treat - new_ctrl

    # Merge and sort by genomic position
    all_pos = np.concatenate([np.asarray(cpg_sites, dtype=int), new_pos])
    all_ctrl = np.concatenate([np.asarray(control_vals), new_ctrl])
    all_treat = np.concatenate([np.asarray(treatment_vals), new_treat])

    order = np.argsort(all_pos)
    cpg_sites_new = all_pos[order].tolist()
    control_vals_new = all_ctrl[order]
    treatment_vals_new = all_treat[order]

    return cpg_sites_new, control_vals_new, treatment_vals_new, int(n_add)


def simulate_dmr_region_with_input_limit(
    chr_name, region_start, region_end, delta_methylation,
    type="good-DMR",
    max_cpgs=50, density="moderate",
    blockiness=False,
    flip_rate_allow=0.05,
    flip_magnitude_max=0.05,
    precision=20,
    min_cpg=5, min_len=50, min_abs_delta=0.1,
    no_delta_methylation=0.08,
    eps=1e-3
):

    """
    Simulate a DMR
    density: "dense" / "sparse" / "moderate"
    """

    if region_end - region_start < 100:
        return None

    cpg_sites = generate_cpg_sites(region_start, region_end, mode=density, max_cpgs=max_cpgs)

    if not cpg_sites:
        return None

    region_end_final = min(region_end, cpg_sites[-1])
    region_length = region_end_final - region_start

    if len(cpg_sites) < min_cpg or region_length < min_len or abs(delta_methylation) < min_abs_delta:
        return None

    mean_control = sample_bimodal_control_mean(
        w_high=0.75,
        m_low=0.10,
        m_high=0.80,
        c_low=18.0,
        c_high=20.0,
        jitter=0.0
        ) 
    
    # Beta shape parameter and sampling
    
    target_sd = 0.1
    var_target = target_sd ** 2
    min_precision = mean_control * (1.0 - mean_control) / var_target - 1.0
    
    effective_precision = max(precision, min_precision)
    
    alpha_control = mean_control * effective_precision
    beta_control  = (1.0 - mean_control) * effective_precision
    
    control_vals = beta.rvs(alpha_control, beta_control, size=len(cpg_sites))
    
    mean_val = float(mean_control)
    max_dev = 0.15
    rand_band = 0.1

    control_vals = np.asarray(control_vals, dtype=float)
    n = len(control_vals)

    is_outlier = np.abs(control_vals - mean_val) > max_dev

    low = max(0.0, mean_val - rand_band)
    high = min(1.0, mean_val + rand_band)

    for i in range(n):
        if not is_outlier[i]:
            continue
        neighbors = []
    
        if i - 1 >= 0 and not is_outlier[i - 1]:
            neighbors.append(control_vals[i - 1])
        
        if i + 1 < n and not is_outlier[i + 1]:
            neighbors.append(control_vals[i + 1])

        if len(neighbors) == 2:
            control_vals[i] = 0.5 * (neighbors[0] + neighbors[1])
        else:
            control_vals[i] = np.random.uniform(low, high)

    #precision = precision
    #alpha_control = mean_control * precision
    #beta_control = (1.0 - mean_control) * precision
    #control_vals = beta.rvs(alpha_control, beta_control, size=len(cpg_sites))
    #control_vals = np.clip(control_vals, 0.001, 0.999)
    
    if mean_control <= 0.4:
        direction = +1
    elif mean_control >= 0.6:
        direction = -1
    else:
        direction = random.choice([+1, -1])

    additive = np.random.normal(loc=delta_methylation, scale=0.03, size=len(cpg_sites))
    additive = np.clip(additive, 0.0001, 0.9999)

    if direction == +1:
        # upregulation：control < treatment
        treatment_vals = np.minimum(control_vals + additive, 0.999)
    else:
        # downregulation:control > treatment
        treatment_vals = np.maximum(control_vals - additive, 0.001)

    control_vals = np.round(control_vals, 3)
    treatment_vals = np.round(treatment_vals ,3)
    #control_vals = np.clip(control_vals, 0.0, 1.0)
    #treatment_vals = np.clip(treatment_vals, 0.0, 1.0)

    # Allow a very small amount of minor reverse movement
    delta_vals = treatment_vals - control_vals
    
    if blockiness:
        if direction == +1:
            reverse_mask = (delta_vals < -eps)
            weak_reverse = reverse_mask & (np.abs(delta_vals) <= flip_magnitude_max)
        else:
            reverse_mask = (delta_vals >  eps)
            weak_reverse = reverse_mask & (np.abs(delta_vals) <= flip_magnitude_max)
    
        reverse_rate       = np.sum(reverse_mask) / len(delta_vals)
        weak_reverse_rate  = np.sum(weak_reverse) / len(delta_vals)
        strong_reverse_any = (reverse_rate - weak_reverse_rate) > 0
    
        if strong_reverse_any or (reverse_rate > flip_rate_allow):
            return None
    else:
        reverse_rate = 0
        weak_reverse_rate = 0        
        if not (np.all(delta_vals <= 0) or np.all(delta_vals >= 0)):
            return None

    # Actual Δ review
    actual_delta = float(np.mean(treatment_vals) - np.mean(control_vals))
    if abs(actual_delta) < min_abs_delta:
        return None

    # Insert boundary points
    if type in {"good-DMR", "sub-DMR"}:
        max_points=10
        min_points=5
        n_opposite=2
        window_bp=500
        flank_cpg_sites, flank_control_vals, flank_treatment_vals, n_flank = add_trailing_flanking_cpgs_v2(
            cpg_sites=cpg_sites,
            control_vals=control_vals,
            treatment_vals=treatment_vals,
            region_start=region_start,
            region_end_final=region_end_final,
            direction=direction,
            alpha_control=alpha_control,
            beta_control=beta_control,
            mean_val=mean_val,
            window_bp=window_bp,
            max_points=max_points,
            min_points=min_points,
            n_opposite=n_opposite,
            delta_cap=min(0.099, no_delta_methylation),
            rng=None
            )
        
    elif type in {"notable-DMR"}:
        max_points=2
        min_points=1
        n_opposite=1
        window_bp=100
        flank_cpg_sites, flank_control_vals, flank_treatment_vals, n_flank = add_trailing_flanking_cpgs_v2(
            cpg_sites=cpg_sites,
            control_vals=control_vals,
            treatment_vals=treatment_vals,
            region_start=region_start,
            region_end_final=region_end_final,
            direction=direction,
            alpha_control=alpha_control,
            beta_control=beta_control,
            mean_val=mean_val,
            window_bp=window_bp,
            max_points=max_points,
            min_points=min_points,
            n_opposite=n_opposite,
            delta_cap=min(0.099, no_delta_methylation),
            rng=None
            )
    else:
        max_points=2
        min_points=1
        n_opposite=1
        window_bp=100
        flank_cpg_sites = flank_control_vals = flank_treatment_vals = np.array([])
        n_flank = 0
    
    return {
        "chr": chr_name,
        "start": region_start,
        "end": region_end_final,
        "CpG_sites": cpg_sites,
        "CpG_count": len(cpg_sites),
        "methylation_control": control_vals.tolist(),
        "methylation_treatment": treatment_vals.tolist(),
        "mean_delta_methylation": round(actual_delta, 3),
        "abs_mean_delta": np.abs(actual_delta),
        "flank_num": n_flank,
        "flank_CpG_sites": flank_cpg_sites,
        "flank_methylation_control": flank_control_vals.tolist(),            
        "flank_methylation_treatment": flank_treatment_vals.tolist(),
        "direction": direction,
        "alpha_control": alpha_control,
        "beta_control": beta_control,
        "max_points": max_points,
        "min_points": min_points,
        "n_opposite": n_opposite,
        "window_bp": window_bp           
    }

def sample_bimodal_control_mean(
    w_high: float = 0.75,
    m_low: float = 0.10,
    m_high: float = 0.80,
    c_low: float = 18.0,
    c_high: float = 20.0,
    jitter: float = 0.0
) -> float:
    """
    Sample mean_control from a bimodal distribution:
    - Mixture of Betas
    - The mode can be precisely controlled to be m through (alpha, beta) = (1 + c*m, 1 + c*(1-m))
    - c Controls the sharpness of the peak (>0 and the larger the value, the "sharper" the peak)    
    """
    
    # Calculate two Beta shape parameters such that their mode values are m_low and m_high respectively
    a_low  = 1.0 + c_low  * m_low
    b_low  = 1.0 + c_low  * (1.0 - m_low)
    a_high = 1.0 + c_high * m_high
    b_high = 1.0 + c_high * (1.0 - m_high)

    # Decide which peak to select based on its weight
    if np.random.rand() < w_high:
        x = beta.rvs(a_high, b_high)
    else:
        x = beta.rvs(a_low, b_low)
    
    # Gaussian perturbation makes peaks more "natural"
    if jitter > 0:
        x = x + np.random.normal(0.0, jitter)
    
    return float(np.clip(x, 0.001, 0.999))


def simulate_nondmr_region(chr_name, region_start, region_end, no_delta_methylation=None, max_cpgs=50, precision=20,density="moderate"):
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

    # If the difference value is not specified, randomly generate one in the range [-0.05, 0.05]
    if no_delta_methylation is None:
        no_delta_methylation = round(truncnorm.rvs((0 - 0.05) / 0.02, (0.1 - 0.05) / 0.02, loc=0.05, scale=0.02),3)

    else:
        no_delta_methylation = np.round(abs(no_delta_methylation), 3)
    
    cpg_sites = generate_cpg_sites(region_start, region_end, mode=density, max_cpgs=max_cpgs)

    if not cpg_sites:
        return None

    region_end_final = min(region_end, cpg_sites[-1])
    
    mean_control = sample_bimodal_control_mean(
        w_high=0.75,   # 0.8 This peak accounts for 75%
        m_low=0.10,    # Low peak mode ~0.1
        m_high=0.80,   # Peak mode ~0.8
        c_low=18.0,    # Sharpness
        c_high=20.0,
        jitter=0.01
    )

    target_sd = 0.1
    var_target = target_sd ** 2
    min_precision = mean_control * (1.0 - mean_control) / var_target - 1.0
    
    effective_precision = max(precision, min_precision)
    
    alpha_control = mean_control * effective_precision
    beta_control  = (1.0 - mean_control) * effective_precision
    
    control_vals = beta.rvs(alpha_control, beta_control, size=len(cpg_sites))
    
    mean_val = float(mean_control)
    max_dev = 0.15
    rand_band = 0.1

    control_vals = np.asarray(control_vals, dtype=float)
    n = len(control_vals)

    is_outlier = np.abs(control_vals - mean_val) > max_dev

    low = max(0.0, mean_val - rand_band)
    high = min(1.0, mean_val + rand_band)

    for i in range(n):
        if not is_outlier[i]:
            continue
        neighbors = []
    
        if i - 1 >= 0 and not is_outlier[i - 1]:
            neighbors.append(control_vals[i - 1])
        
        if i + 1 < n and not is_outlier[i + 1]:
            neighbors.append(control_vals[i + 1])

        if len(neighbors) == 2:
            control_vals[i] = 0.5 * (neighbors[0] + neighbors[1])
        else:
            control_vals[i] = np.random.uniform(low, high)


    #precision = precision
    #alpha_control = mean_control * precision
    #beta_control = (1.0 - mean_control) * precision
    #control_vals = beta.rvs(alpha_control, beta_control, size=len(cpg_sites))
    #control_vals = np.clip(control_vals, 0.001, 0.999)
    
    # Ensure that both upward and downward adjustments exist
    if mean_control <= 0.4:
        direction = +1
    elif mean_control >= 0.6:
        direction = -1
    else:
        direction = random.choice([+1, -1])

    additive = np.random.uniform(0.001, no_delta_methylation,size=len(cpg_sites))
    additive = np.clip(additive, 0.0001, 0.9999)

    if direction == +1:
        # upregulation：control < treatment
        treatment_vals = np.minimum(control_vals + additive, 0.999)
    else:
        # downregulation:control > treatment
        treatment_vals = np.maximum(control_vals - additive, 0.001)

    control_vals = np.round(control_vals, 3)
    treatment_vals = np.round(treatment_vals ,3)
    
    actual_delta = round(float(np.mean(treatment_vals) - np.mean(control_vals)), 3)

    if actual_delta > 0.1 or actual_delta < -0.1:
        return None

    max_points=2
    min_points=1
    n_opposite=1
    window_bp=100

    flank_cpg_sites, flank_control_vals, flank_treatment_vals, n_flank = add_trailing_flanking_cpgs_v2(
        cpg_sites=cpg_sites,
        control_vals=control_vals,
        treatment_vals=treatment_vals,
        region_start=region_start,
        region_end_final=region_end_final,
        direction=direction,
        alpha_control=alpha_control,
        beta_control=beta_control,
        mean_val=mean_val,
        window_bp=window_bp,
        max_points=max_points,
        min_points=min_points,
        n_opposite=n_opposite,
        delta_cap=min(0.099, no_delta_methylation),
        rng=None,
        )

    flank_cpg_delta = flank_treatment_vals - flank_control_vals

    return {
        "chr": chr_name,
        "start": region_start,
        "end": region_end_final,
        "CpG_sites": cpg_sites,
        "CpG_count": len(cpg_sites),
        "methylation_control": control_vals.tolist(),
        "methylation_treatment": treatment_vals.tolist(),
        "mean_delta_methylation": actual_delta,
        "flank_num": n_flank,
        "flank_CpG_sites": flank_cpg_sites,
        "flank_methylation_control": flank_control_vals.tolist(),            
        "flank_methylation_treatment": flank_treatment_vals.tolist()
    }


def inject_reverse_blocks(
    dmr: dict,
    gap_min: int = 2, 
    gap_max: int = 4,
    block_min: int = 1,
    block_max: int = 3,
    eps: float = 1e-9
) -> None:
    """
    Insert "reverse block" on the existing DMR according to the rules:
    - The length of the inverse block is randomly chosen within the range of [block_min, block_max] (inclusive)
    - Insert the next block after every [gap_min, gap_max] points (inclusive)
    - Inverted method: Invert delta;
    """
    ctrl = np.asarray(dmr["methylation_control"], dtype=float)
    trt  = np.asarray(dmr["methylation_treatment"], dtype=float)
    n    = len(ctrl)
    if n == 0:
        return

    i = 0
    while True:
        # Skip a random interval and start a block
        gap = random.randint(gap_min, gap_max)
        start = i + gap
        if start >= n:
            break
        
        blen = random.randint(block_min, block_max)
        end = min(n, start + blen)
        
        tmp = trt[start:end].copy()
        trt[start:end] = ctrl[start:end]
        ctrl[start:end] = tmp        
        i = end

    dmr["methylation_control"]   = np.round(ctrl, 3).tolist()
    dmr["methylation_treatment"] = np.round(trt, 3).tolist()
    dmr_delta = np.asarray(dmr["methylation_treatment"], dtype=float) - np.asarray(dmr["methylation_control"], dtype=float) 
    dmr["mean_delta_methylation"] = round(float(np.mean(trt - ctrl)), 3)


# Simulate inconsistent DMRs
def simulate_dmr_with_inconsistent_points(chr_name, region_start, region_end, delta_methylation, type="inconsistent-DMR", max_cpgs=50, cpg_blocking=False, density="moderate",no_delta_methylation=0.08):
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
        type=type,
        max_cpgs=max_cpgs,
        no_delta_methylation=no_delta_methylation,
        density=density
    )

    if dmr is None:
        return None

    inject_reverse_blocks(
        dmr,
        gap_min=1,gap_max=3,
        block_min=1,block_max=3
        )

    control_arr = np.array(dmr["methylation_control"])
    treatment_arr = np.array(dmr["methylation_treatment"])
    new_delta = treatment_arr - control_arr
    dmr["mean_delta_methylation"] = round(np.mean(new_delta), 3)
    dmr["abs_mean_delta"] = round(np.mean(np.abs(new_delta)),3)
    if new_delta[-1] > 0:
        direction = +1
    else:
        direction = -1
    
    flank_cpg_sites, flank_control_vals, flank_treatment_vals, n_flank = add_trailing_flanking_cpgs(
        cpg_sites=np.array(dmr["CpG_sites"]),
        control_vals=dmr["methylation_control"],
        treatment_vals=dmr["methylation_treatment"],
        region_start=region_start,
        region_end_final=region_end,
        direction=direction,
        alpha_control=dmr["alpha_control"],
        beta_control=dmr["beta_control"],
        window_bp=dmr["window_bp"],
        max_points=dmr["max_points"],
        min_points=dmr["min_points"],
        n_opposite=dmr["n_opposite"],
        delta_cap=min(0.099, abs(delta_methylation)*0.1),
        rng=None
        )

    new_delta2 = flank_treatment_vals - flank_control_vals
    dmr["flank_num"] = n_flank
    dmr["flank_CpG_sites"] = flank_cpg_sites
    dmr["flank_methylation_control"] = flank_control_vals.tolist()
    dmr["flank_methylation_treatment"] = flank_treatment_vals.tolist()
    dmr["direction"] = direction
    return dmr


def shrink_ends_to_subdmr_soft(
    dmr: Dict[str, Any],
    target_max_abs_delta: float = 0.1,
    max_shrink_delta: float = 0.05,
    eps: float = 1e-9
) -> Optional[Dict[str, Any]]:
    """
    Start with "gentle shrinkage" from both ends: Set Δ at the endpoint to sign(meanΔ) * U(0, max_shrink_delta),
    Gradually pull the overall mean Δ towards 0 (but without reversing it) until it is just about to fall below target_max_abs_delta;
    Then, roll back to the previous step and output the intermediate continuous sub-segments as sub-DMRs, ensuring that |mean Δ| > target_max_abs_delta.

    Note: Here, "Δ" refers to per-site delta = treatment - control.
    """

    #sites_base = np.asarray(dmr["CpG_sites"], dtype=int)
    sites = np.asarray(dmr["CpG_sites"], dtype=int)
    flank_num = dmr["flank_num"]
    ctrl  = np.asarray(dmr["methylation_control"], dtype=float)
    trt   = np.asarray(dmr["methylation_treatment"], dtype=float)
    #ctrl  = np.asarray(dmr["flank_methylation_control"], dtype=float)
    #trt   = np.asarray(dmr["flank_methylation_treatment"], dtype=float)

    n = len(sites)
    if n == 0:
        return None

    delta = trt - ctrl
    mean_delta = float(np.mean(delta))

    if abs(mean_delta) <= target_max_abs_delta + eps:
        return None

    modified_stack = []
    left, right = 0, n - 1

    # Gradual shrinkage: Each time, replace the endpoint with a larger |Δ| with a smaller one in the same direction, 
    # Δ (U(0, max_shrink_delta))
    while True:
        delta = trt - ctrl
        mean_delta = float(np.mean(delta))
        if abs(mean_delta) <= target_max_abs_delta + eps:
            break

        if left > right:
            # Extreme case: After processing at both ends, the value still does not fall below the threshold
            break

        # Skip endpoints that are "nearly zero"
        while left <= right and abs(delta[left]) <= eps:
            left += 1
        while left <= right and abs(delta[right]) <= eps:
            right -= 1
        if left > right:
            break

        choose_left = abs(delta[left]) >= abs(delta[right])
        idx = left if choose_left else right

        sign_md = 1.0 if mean_delta >= 0 else -1.0
        tiny = np.random.uniform(0.0, max_shrink_delta)
        target_delta = sign_md * tiny

        trt_old = trt[idx]
        trt[idx] = ctrl[idx] + target_delta
        trt[idx] = float(np.clip(trt[idx], 0.001, 0.999))
        modified_stack.append((idx, trt_old))

        if choose_left:
            left += 1
        else:
            right -= 1

    if not modified_stack:
        return None

    # Take a step back: revert to the last modification and ensure that |mean Δ| > threshold
    idx_last, trt_last = modified_stack.pop()
    trt[idx_last] = float(np.clip(trt_last, 0.001, 0.999))

    # The overall Δ after rollback must be greater than the threshold, otherwise it is impossible to construct a sub-DMR
    delta = trt - ctrl
    mean_delta = float(np.mean(delta))
    if abs(mean_delta) <= target_max_abs_delta + eps:
        return None

    sub = dict(dmr)
    sub["methylation_control"]   = np.round(ctrl, 3).tolist()
    sub["methylation_treatment"] = np.round(trt, 3).tolist()
    sub["mean_delta_methylation"] = round(float(np.mean(trt - ctrl)), 3)
    n = len(sub["methylation_control"])
    sub["flank_methylation_control"][:n] = sub["methylation_control"]
    sub["flank_methylation_treatment"][:n] = sub["methylation_treatment"]
    
    sub["is_subdmr"] = True

    if abs(sub["mean_delta_methylation"]) <= target_max_abs_delta + eps:
        return None

    return sub

def simulate_dmr_with_subdmr_points(chr_name, region_start, region_end, delta_methylation, type="sub-DMR", max_cpgs=50,no_delta_methylation=0.08):
    dmr = simulate_dmr_region_with_input_limit(
        chr_name=chr_name,
        region_start=region_start,
        region_end=region_end,
        delta_methylation=delta_methylation,
        type=type,
        max_cpgs=max_cpgs,
        no_delta_methylation=no_delta_methylation,
        blockiness=False,
    )

    if dmr is None:
        return None  # If simulation fails, return directly
    
    sub = shrink_ends_to_subdmr_soft(dmr, target_max_abs_delta=0.1)
    return sub

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
    cat,
    dmr,
    n_control: int,
    n_treatment: int,
    coverage_mean: int,
    coverage_std: int,
    group_std: float,
    dmr_missing_max: float=0.1,
    sample_missing_max: float=0.1,
    output_dir: str = "./"
):
    os.makedirs(output_dir, exist_ok=True)

    all_samples = {
        f"control_sample_{i}": [] for i in range(1, n_control + 1)
    }
    all_samples.update({
        f"treatment_sample_{i}": [] for i in range(1, n_treatment + 1)
    })

    dmr_missing_rate = np.random.uniform(0, dmr_missing_max) 
    sample_missing_rate = np.random.uniform(0, sample_missing_max)

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
                #alpha, beta_ = beta_params_from_mean_std(mean_meth, group_std)
                #sample_meth = beta.rvs(alpha, beta_)
                raw_meth = np.random.normal(loc=mean_meth, scale=group_std)
                max_dev = 0.3 
                low = max(0.0, mean_meth - max_dev)
                high = min(1.0, mean_meth + max_dev)
                sample_meth = float(np.clip(raw_meth, low, high))


                all_samples[sample_key].append({
                    "chr": dmr["chr"],
                    "start": cpg_pos,
                    "end": cpg_pos + 1,
                    "coverage": coverage,
                    "methylation_level": round(sample_meth, 4),
                    "type": cat,
                })

    for sample_key, records in all_samples.items():
        df = pd.DataFrame(records)
        df.to_csv(os.path.join(output_dir, f"{sample_key}.tsv"), sep="\t", index=False)

    return {
        "output_dir": output_dir,
        "dmr_missing_rate": round(dmr_missing_rate, 3),
        "sample_missing_rate": sample_missing_rate
    }

def simulate_group_samples_with_missing_add_flanking(
    cat,
    dmr,
    n_control: int,
    n_treatment: int,
    coverage_mean: int,
    coverage_std: int,
    group_std: float,
    dmr_missing_max: float=0.1,
    sample_missing_max: float=0.1,
    output_dir: str = "./"
):
    os.makedirs(output_dir, exist_ok=True)

    all_samples = {
        f"control_sample_{i}": [] for i in range(1, n_control + 1)
    }
    all_samples.update({
        f"treatment_sample_{i}": [] for i in range(1, n_treatment + 1)
    })

    dmr_missing_rate = np.random.uniform(0, dmr_missing_max) 
    sample_missing_rate = np.random.uniform(0, sample_missing_max)

    cpg_sites = dmr["CpG_sites"]
    n_cpgs = len(cpg_sites)
    flank_cpg_sites = dmr["flank_CpG_sites"]
    n_flank_cpgs = len(flank_cpg_sites)
    n_missing_cpgs = int(np.floor(dmr_missing_rate * n_cpgs))
    missing_cpg_indices = set(random.sample(range(n_cpgs), n_missing_cpgs))

    for group, meth_values, n_samples in [
        ("control", dmr["flank_methylation_control"], n_control),
        ("treatment", dmr["flank_methylation_treatment"], n_treatment)
    ]:
        for sample_id in range(1, n_samples + 1):
            sample_key = f"{group}_sample_{sample_id}"
            #for idx, cpg_pos in enumerate(cpg_sites):
            for idx, cpg_pos in enumerate(flank_cpg_sites):
                is_missing_cpg = idx in missing_cpg_indices
                is_missing_sample = is_missing_cpg and (np.random.rand() < sample_missing_rate)

                if is_missing_sample:
                    continue
                
                # Use beta distribution
                coverage = int(np.clip(np.random.normal(coverage_mean, coverage_std), 1, 100))
                mean_meth = meth_values[idx]
                #alpha, beta_ = beta_params_from_mean_std(mean_meth, group_std)
                #sample_meth = beta.rvs(alpha, beta_)
                raw_meth = np.random.normal(loc=mean_meth, scale=group_std)
                max_dev = 0.4 
                low = max(0.0, mean_meth - max_dev)
                high = min(1.0, mean_meth + max_dev)
                sample_meth = float(np.clip(raw_meth, low, high))

                if cpg_pos in cpg_sites:
                    all_samples[sample_key].append({
                        "chr": dmr["chr"],
                        "start": cpg_pos,
                        "end": cpg_pos + 1,
                        "coverage": coverage,
                        "methylation_level": round(sample_meth, 4),
                        "type": cat
                    })
                else:
                    all_samples[sample_key].append({
                        "chr": dmr["chr"],
                        "start": cpg_pos,
                        "end": cpg_pos + 1,
                        "coverage": coverage,
                        "methylation_level": round(sample_meth, 4),
                        "type": "flank-DMR"
                    })

    for sample_key, records in all_samples.items():
        df = pd.DataFrame(records)
        df.to_csv(os.path.join(output_dir, f"{sample_key}.tsv"), sep="\t", index=False)

    return {
        "output_dir": output_dir,
        "dmr_missing_rate": round(dmr_missing_rate, 3),
        "sample_missing_rate": sample_missing_rate
    }


def simulate_group_samples_with_missing_add_flanking_v2(
    cat,
    dmr,
    n_control: int,
    n_treatment: int,
    coverage_mean: int,
    coverage_std: int,
    group_std: float,
    dmr_missing_max: float = 0.1,
    sample_missing_max: float = 0.1,
    output_dir: str = "./"
):
    os.makedirs(output_dir, exist_ok=True)

    # 是否存在单样本组：只要有一组样本数为1，就启用“单样本模式”
    single_group_mode = (n_control == 1) or (n_treatment == 1)

    all_samples = {
        f"control_sample_{i}": [] for i in range(1, n_control + 1)                
    }
    all_samples.update({
        f"treatment_sample_{i}": [] for i in range(1, n_treatment + 1)                
    })

    cpg_sites = dmr["CpG_sites"]
    n_cpgs = len(cpg_sites)
    flank_cpg_sites = dmr["flank_CpG_sites"]
    n_flank_cpgs = len(flank_cpg_sites)

    # 如果不是单样本模式，按原逻辑生成缺失率和缺失位点
    if not single_group_mode:
        dmr_missing_rate = np.random.uniform(0, dmr_missing_max)
        sample_missing_rate = np.random.uniform(0, sample_missing_max)

        n_missing_cpgs = int(np.floor(dmr_missing_rate * n_cpgs))
        # 注意：这里仍然是按原来逻辑，从 DMR 内 CpG 的 index 中抽
        missing_cpg_indices = set(random.sample(range(n_cpgs), n_missing_cpgs))
    else:
        # 单样本模式：不启用缺失
        dmr_missing_rate = 0.0
        sample_missing_rate = 0.0
        missing_cpg_indices = set()

    for group, meth_values, n_samples in [
                ("control", dmr["flank_methylation_control"], n_control),
                ("treatment", dmr["flank_methylation_treatment"], n_treatment)
    ]:
        for sample_id in range(1, n_samples + 1):
            sample_key = f"{group}_sample_{sample_id}"

            for idx, cpg_pos in enumerate(flank_cpg_sites):
                # 非单样本模式下才应用缺失逻辑
                if not single_group_mode:
                    is_missing_cpg = idx in missing_cpg_indices
                    is_missing_sample = is_missing_cpg and (np.random.rand() < sample_missing_rate)
                    if is_missing_sample:
                        continue

                # 覆盖度保持原来的随机逻辑
                coverage = int(np.clip(np.random.normal(coverage_mean, coverage_std), 1, 100))

                mean_meth = meth_values[idx]

                if single_group_mode and n_samples == 1:
                    # 单样本组：直接用 mean_meth，不再按方差模拟
                    sample_meth = float(np.clip(mean_meth, 0.0, 1.0))
                else:
                    # 原逻辑：在 mean_meth 周围按 group_std 波动
                    raw_meth = np.random.normal(loc=mean_meth, scale=group_std)
                    max_dev = 0.4
                    low = max(0.0, mean_meth - max_dev)
                    high = min(1.0, mean_meth + max_dev)
                    sample_meth = float(np.clip(raw_meth, low, high))

                record = {
                    "chr": dmr["chr"],
                    "start": cpg_pos,
                    "end": cpg_pos + 1,
                    "coverage": coverage,
                    "methylation_level": round(sample_meth, 4),                
                }

                if cpg_pos in cpg_sites:
                    record["type"] = cat
                else:
                    record["type"] = "flank-DMR"

                all_samples[sample_key].append(record)

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
            return {
                "output_dir": output_dir,
                "dmr_missing_rate": round(dmr_missing_rate, 3),
                "sample_missing_rate": sample_missing_rate,
                "delta_check_passed": True,
                "delta_values": deltas
                }

    return {
        "output_dir": output_dir,
        "delta_check_passed": False,
        "error": "Failed to simulate valid non-DMR after max retries"
        }   

def write_simulation_parameters(
    df: pd.DataFrame,
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
    density: str = "moderate",
    min_gap: int = 10,
    max_gap: int = 50,
    sample_missing: float = 0.1,
    cpg_missing: float = 0.1
):
    log_path = os.path.join(output_dir, "para.log")
    
    total_dmr = len(df)
    counts = df["category"].value_counts().to_dict()
    
    with open(log_path, "w") as f:
        f.write("# Simulation Parameters Log\n")
        f.write(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write(f"Chromosome: {chr_name}\n")
        f.write(f"Start position: {start_pos}\n")
        f.write(f"Total simulated regions: {total_dmr}\n")
        f.write(f" - good-DMR: {counts.get('good-DMR', 0)} ({dmr_per} of all DMRs)\n")
        f.write(f" - inconsistent-DMR: {counts.get('inconsistent-DMR', 0)} ({dmr_inconsis_per} of all DMRs)\n")
        f.write(f" - sub-DMR: {counts.get('sub-DMR', 0)} ({dmr_sub_per} of all DMRs)\n")
        f.write(f" - notable-DMR: {counts.get('notable-DMR', 0)} ({dmr_notable_per} of all DMRs)\n")
        f.write(f" - non-DMR: {counts.get('non-DMR', 0)}\n\n")

        f.write(f"DMR density mode: {density}\n")
        f.write(f"Region length: N(mean={length_mean} bp, std={length_std} bp)\n")
        f.write(f"Max CpGs per region: {max_cpgs}\n")
        f.write(f"CpG distance : {min_gap} ~ {max_gap} bp\n")
        f.write(f"Sample missing: {sample_missing}\n")
        f.write(f"CpG missing: {cpg_missing}\n")

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


FORBID_AFTER = {"sub-DMR", "good-DMR", "inconsistent-DMR"}

def allowed_next_categories(last_cat, remaining: dict):
    base = [k for k, v in remaining.items() if v > 0]

    if last_cat in FORBID_AFTER:
        return [k for k in base if k != "notable-DMR"]

    if last_cat == "notable-DMR":
        return ["non-DMR"]

    return base


def calc_sampling_weights(candidates, remaining, last_cat):
    strict_total = sum(v for k, v in remaining.items() if k != "non-DMR")
    strict_total = max(strict_total, 1)

    weights = []
    for k in candidates:
        if k == "notable-DMR":
            w = max(2, remaining.get("notable-DMR", 0) * 2)
        elif k == "non-DMR":
            w = max(1, remaining.get("non-DMR", 0))
        else:
            w = remaining.get(k, 0)
        weights.append(max(w, 0))
    
    return weights

# Randomly shuffle generation order + correct classification output + insert labels by category
def simulate_mixed_regions_randomized(
    total_dmr: int,
    mean_delta: float,
    n_control: int,
    n_treatment: int,
    coverage_mean: int = 30,
    coverage_std: int = 5,
    output_dir: str = "./out",
    chr_name: str = "chr1",
    start_pos: int = 10000,
    length_mean: int= 1000,
    length_std: int = 300,
    max_cpgs: int = 200,
    dmr_per: float = 0.3,
    dmr_notable_per: float = 0.05,
    dmr_inconsis_per: float = 0.1,
    dmr_sub_per: float = 0.05,
    density: str = "moderate",
    min_gap: int = 10,
    max_gap: int =50,
    cpg_blockiness: bool = False,
    good_flip_rate_allow: float = 0.05,
    good_flip_magnitede_max: float = 0.05,
    good_precision: int = 20,
    no_delta_methylation: float = None,
    no_precision: int = 20,
    dmr_missing_max: float = 0.1,
    sample_missing_max: float = 0.1,
    seed: int = 42,
    max_attempts_per_slot: int =  200
):
    np.random.seed(seed)
    random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    # Target quota
    target = {
        "non-DMR": int(total_dmr - int(total_dmr*dmr_per) - int(total_dmr*dmr_inconsis_per) - int(total_dmr*dmr_notable_per) - int(total_dmr*dmr_sub_per)),
        "good-DMR": int(total_dmr * dmr_per),
        "inconsistent-DMR": int(total_dmr * dmr_inconsis_per),
        "sub-DMR": int(total_dmr * dmr_sub_per),
        "notable-DMR": int(total_dmr * dmr_notable_per),    
    }

    # If the sum does not equal total_dmr due to rounding errors, the difference should be added to non-DMR
    diff = total_dmr - sum(target.values())
    if diff != 0:
        target["non-DMR"] = max(0, target["non-DMR"] + diff)
    
    remaining = target.copy()

    def unfinished_categories():
        return [k for k,v in remaining.items() if k != "non-DMR" and v > 0]

    def advance_position(last_end, min_gap, max_gap):
        return last_end + random.randint(min_gap, max_gap)

    regions = []
    pos = start_pos
    attempts_for_slot = {k: 0 for k in target}
    last_cat = None
    real_end = start_pos

    while unfinished_categories():
        candidates = allowed_next_categories(last_cat, remaining)
            
        strict_cats_left = [k for k, v in remaining.items() if k != "non-DMR" and v > 0]
        only_notable_left = (len(strict_cats_left) == 1 and strict_cats_left[0] == "notable-DMR")

        if not candidates or only_notable_left:
            bump = max(1, 3 * len(strict_cats_left))
            remaining["non-DMR"] = remaining.get("non-DMR", 0) + bump
            
            candidates = allowed_next_categories(last_cat, remaining)

        # Weighted sampling of remaining quota (only sampling from the candidate set)
        weights = calc_sampling_weights(candidates, remaining, last_cat)
        cat = random.choices(candidates, weights=weights, k=1)[0]

        # Exceeding the attempt limit protection
        attempts_for_slot[cat] += 1
        if attempts_for_slot[cat] > max_attempts_per_slot:
            raise RuntimeError(
                f"Category '{cat}' could not be filled to target after {max_attempts_per_slot} attempts. "
                "Consider relaxing validation/parameters or enlarging search space."
            )

        # Interval length extraction
        region_length = int(np.random.normal(loc=length_mean, scale=length_std))
        region_length = max(100, region_length)
        
        if cat in {"good-DMR","sub-DMR"}:
            pos=real_end + random.randint(int(min_gap/2),min_gap)
            end = pos + region_length
        else:
            pos=real_end + random.randint(min_gap, max_gap)
            end = pos + region_length
        
        #end = pos + region_length

        region = None
        group_std = 0.05
        effective_delta = 0.0

        if cat == "non-DMR":
            region = simulate_nondmr_region(
                chr_name, pos, end,
                max_cpgs=max_cpgs,
                density=density,
                precision=no_precision,
                no_delta_methylation=no_delta_methylation
            )
            group_std = 0.05
            effective_delta = 0.0

        elif cat == "inconsistent-DMR":
            delta = float(np.clip(np.random.normal(loc=mean_delta, scale=0.05), 0.1, 0.9))
            region = simulate_dmr_with_inconsistent_points(
                chr_name, pos, end, delta,
                type=cat,
                max_cpgs=max_cpgs,
                density = density,
                no_delta_methylation=no_delta_methylation
            )
            group_std = 0.05
            effective_delta = delta

        elif cat == "sub-DMR":
            delta = float(np.clip(np.random.normal(loc=mean_delta, scale=0.05), 0.1, 0.9))
            region = simulate_dmr_with_subdmr_points(
                chr_name, pos, end, delta,
                type=cat,
                no_delta_methylation=no_delta_methylation,
                max_cpgs=max_cpgs
            )
            group_std = 0.05
            effective_delta = delta 

        elif cat == "notable-DMR":
            delta = float(np.clip(np.random.normal(loc=mean_delta, scale=0.05), 1e-4, 0.9))
            region = simulate_dmr_region_with_input_limit(
                chr_name, pos, end, delta,
                type=cat,
                max_cpgs=max_cpgs,
                blockiness=False,
                precision=good_precision,
                no_delta_methylation=no_delta_methylation,
                density=density                        
            )
            group_std = float(min(0.5, round(1.5 * abs(delta), 3)))
            #group_std = float(min(0.25, round(abs(delta), 3)))
            effective_delta = delta

        else:  # good-DMR
            delta = float(np.clip(np.random.normal(loc=mean_delta, scale=0.05), 0.1, 0.9))
            region = simulate_dmr_region_with_input_limit(
                chr_name, pos, end, delta,
                type=cat,
                max_cpgs=max_cpgs,
                density = density,
                blockiness = cpg_blockiness,
                flip_rate_allow=good_flip_rate_allow,
                flip_magnitude_max=good_flip_magnitede_max,
                precision=good_precision,
                no_delta_methylation=no_delta_methylation,
                min_cpg=5, min_len=50, min_abs_delta=0.1,
                eps=1e-3
            )
            group_std = 0.05
            effective_delta = delta

        accepted_this_round = []

        if region:
            candidate_list = region if isinstance(region, list) else [region]

            # non-DMR needs to be verified first; other categories can directly pass as candidates
            if cat == "non-DMR":
                for r in candidate_list:
                    if remaining[cat] <= 0:
                        break
                    if "mean_delta_methylation" not in r:
                        r["mean_delta_methylation"] = 0.0
                    r["category"] = cat

                    # Perform verification
                    validation = simulate_group_samples_with_validation(
                        dmr=r,
                        n_control=n_control,
                        n_treatment=n_treatment,
                        coverage_mean=coverage_mean,
                        coverage_std=coverage_std,
                        group_std=group_std,
                        output_dir=output_dir
                    )
                    passed = True if validation is None else bool(validation.get("delta_check_passed", True))
                    if passed:
                        accepted_this_round.append(r)
                    # If it fails, just ignore it and wait for the next round of attempt

            else:
                take = min(remaining[cat], len(candidate_list))
                for r in candidate_list[:take]:
                    if "mean_delta_methylation" not in r:
                        r["mean_delta_methylation"] = float(effective_delta)
                    r["category"] = cat
                    accepted_this_round.append(r)

        # For the "admitted" areas, documents are actually written, lists are recorded, and quotas are deducted
        for r in accepted_this_round:
            region_output_dir = os.path.join(output_dir, f"{chr_name}_{r['start']}_{r['end']}_{cat}")
            os.makedirs(region_output_dir, exist_ok=True)

            if cat in {"good-DMR","non-DMR","inconsistent-DMR","notable-DMR","sub-DMR"}:
                simulate_group_samples_with_missing_add_flanking_v2(
                    cat=cat,
                    dmr=r,
                    n_control=n_control,
                    n_treatment=n_treatment,
                    coverage_mean=coverage_mean,
                    coverage_std=coverage_std,
                    group_std=group_std,
                    dmr_missing_max=dmr_missing_max,
                    sample_missing_max=sample_missing_max,
                    output_dir=region_output_dir                                                   
                )
            else:
                simulate_group_samples_with_missing_add_flanking(
                    cat=cat,
                    dmr=r,
                    n_control=n_control,
                    n_treatment=n_treatment,
                    coverage_mean=coverage_mean,
                    coverage_std=coverage_std,
                    group_std=group_std,
                    dmr_missing_max=dmr_missing_max,
                    sample_missing_max=sample_missing_max,
                    output_dir=region_output_dir                                    
                )

            regions.append(r)
            remaining[cat] -= 1
            attempts_for_slot[cat] = 0
        if region is None:
           continue

        real_end = region["flank_CpG_sites"][-1]

        #if last_cat in {"good-DMR", "sub-DMR"}:
       # if cat in {"good-DMR", "sub-DMR"}:
       #     min_gaps, max_gaps = 500+min_gap, 500+max_gap
       # elif cat in {"non-DMR","inconsistent-DMR","notable-DMR"}:
       #     min_gaps, max_gaps = 100+min_gap, 100+max_gap
       # else:
       #     min_gaps, max_gaps = 501, 1000
        last_cat = cat 
        #pos = advance_position(end, min_gap=min_gaps, max_gap=max_gaps)

    # Summary output
    rows = [{
        "chr": r.get("chr", chr_name),
        "start": r["start"],
        "end": r["end"],
        "CpG_count": r.get("CpG_count", np.nan),
        "treatment_methylation": np.round(float(np.mean(r["methylation_treatment"])),3),
        "control_methylation": np.round(float(np.mean(r["methylation_control"])),3),
        "mean_delta_methylation": r.get("mean_delta_methylation", 0.0),
        "category": r["category"],    
    } for r in regions]
    df = pd.DataFrame(rows)

    # Record parameters
    write_simulation_parameters(
        df=df,
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
        min_gap=min_gap,
        max_gap=max_gap    
    )
    
    final_counts = df["category"].value_counts().to_dict()

    # Collect mismatch messages (no exception)
    messages = []
    for k in target:
        got = final_counts.get(k, 0)
        want = target[k]
        if got != want:
            messages.append(f"[Notice] Category '{k}': generated {got} vs target {want}")

    # Report how many extra non-DMR were produced (if any)
    extra_non_dmr = max(0, final_counts.get("non-DMR", 0) - target.get("non-DMR", 0))

    if messages:
        print("\n".join(messages))
        if extra_non_dmr > 0:
            print(f"[Info] To maintain adjacency/spacing constraints, an additional {extra_non_dmr} 'non-DMR' regions were generated beyond the preset quota.")

    return df

def main():
    parser = argparse.ArgumentParser(description="Simulate DMR regions")

    # —— Basic simulation parameters ——
    parser.add_argument("--total_dmr", type=int, default=100, help="total number of simulated DMRs")
    parser.add_argument("--mean_delta", type=float, default=0.3, help="target mean methylation difference (can be +/-)")
    parser.add_argument("--n_control", type=int, default=5, help="number of control samples")
    parser.add_argument("--n_treatment", type=int, default=5, help="number of treatment samples")
    parser.add_argument("--coverage_mean", type=int, default=30, help="mean sequencing coverage per CpG")
    parser.add_argument("--coverage_std", type=int, default=5, help="std of sequencing coverage per CpG")

    # —— Output and interval parameters ——
    parser.add_argument("--output_dir", type=str, default="./out", help="output directory")
    parser.add_argument("--chr_name", type=str, default="chr1", help="chromosome name")
    parser.add_argument("--start_pos", type=int, default=10000, help="start genomic coordinate")
    parser.add_argument("--length_mean", type=int, default=1000, help="mean region length (bp)")
    parser.add_argument("--length_std", type=int, default=300, help="std of region length (bp)")
    parser.add_argument("--max_cpgs", type=int, default=200, help="max CpGs per region (hard cap)")

    # —— Proportion of various DMRs ——
    parser.add_argument("--dmr_per", type=float, default=0.25, help="proportion of good DMRs")
    parser.add_argument("--dmr_notable_per", type=float, default=0.1, help="proportion of notable DMRs")
    parser.add_argument("--dmr_inconsis_per", type=float, default=0.1, help="proportion of inconsistent DMRs")
    parser.add_argument("--dmr_sub_per", type=float, default=0.1, help="proportion of sub DMRs")

    # —— CpG density and block variation ——
    parser.add_argument("--density", type=str, choices=["dense", "sparse", "moderate"], default="moderate", help="CpG density mode")
    parser.add_argument("--cpg_blockiness", action="store_true", help="enable block-wise density changes within a region")
    
    # —— missing CpG and sample control ——
    parser.add_argument("--dmr_missing_max", type=float, default=0.1, help="CpG missing rate in DMR")
    parser.add_argument("--sample_missing_max", type=float, default=0.1, help="sample missing rate within group")

    # —— good-DMR Consistency soft constraints and sampling parameters ——
    parser.add_argument("--good_flip_rate_allow", type=float, default=0.05, help="allowed fraction of reverse-direction CpGs within a good-DMR")
    parser.add_argument("--good_flip_magnitude_max", type=float, default=0.05, help="max |delta| allowed for reverse-direction CpGs in a good-DMR")
    parser.add_argument("--good_precision", type=int, default=20, help="beta precision parameter for good-DMR generation")

    # —— non-DMR Consistency soft constraints and sampling parameters ——
    parser.add_argument("--no_delta_methylation", type=float, default=0.08, help="mean methylation difference for non-DMR")
    parser.add_argument("--no_precision", type=int, default=20, help="beta precision parameter for non-DMR generation")

    # —— Operation and generation control ——
    parser.add_argument("--min_gap", type=int, default=10, help="min inter-region gap when advancing along the chromosome (bp)")
    parser.add_argument("--max_gap", type=int, default=50, help="max inter-region gap when advancing along the chromosome (bp)")
    parser.add_argument("--seed", type=int, default=42, help="random seed")
    parser.add_argument("--max_attempts_per_slot", type=int, default=200, help="max attempts for generating a region of a given class before skipping")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

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
        good_precision=args.good_precision,
        density=args.density,

        dmr_per=args.dmr_per,
        dmr_notable_per=args.dmr_notable_per,
        dmr_inconsis_per=args.dmr_inconsis_per,
        dmr_sub_per=args.dmr_sub_per,

        #cpg_blockiness=args.cpg_blockiness,
        #good_flip_rate_allow=args.good_flip_rate_allow,
        #good_flip_magnitude_max=args.good_flip_magnitude_max,

        no_delta_methylation=args.no_delta_methylation,
        no_precision=args.no_precision,

        dmr_missing_max=args.dmr_missing_max,
        sample_missing_max=args.sample_missing_max,

        min_gap=args.min_gap,
        max_gap=args.max_gap,

        max_attempts_per_slot=args.max_attempts_per_slot,
        seed=args.seed
    )

    print("Simulation completed.")
    df.to_csv(f"{args.output_dir}/DMRs.txt",sep="\t",header=True,index=False)

### 记得删除
def test(output_dir="./test"):
    os.makedirs(output_dir, exist_ok=True)
    
    df = simulate_mixed_regions_randomized(
        total_dmr=200,
        mean_delta=0.3,
        n_control=3,
        n_treatment=3,
        output_dir=output_dir,
        dmr_inconsis_per=0.1,
        dmr_notable_per=0.1,
        dmr_sub_per=0.1,
        no_delta_methylation=0.08,
        cpg_blockiness=False,
        #min_gap=500,
        #max_gap=600,
        density="moderate",
    )
        
    print("Simulation completed.")
    df.to_csv(f"{output_dir}/DMRs.txt",sep="\t",header=True,index=False)

if __name__ == "__main__":
    #test("./test26")
    main()
