import pandas as pd
import numpy as np
import WBR as wbr
import BH_adjust

def calculate_score(delta_m, distance, d_max, alpha=1.0, beta=0.01):
    """ Calculate the merge priority score """
    return alpha * np.abs(delta_m) - beta * (distance / d_max)

def find_blocks_greedy(out_data, delta_m_mean_threshold=0.1, min_cpg_count=5, alpha=1.0, beta=0.5,
                       group1="g1", group2="g2",
                       position_col="Position", delta_m_col="mean_diff", block_col="Block", chr_col="chr1"):
    """
    Adopting a greedy strategy to extend CpG blocks by dividing them into 'Blocks', 
    with each block executing independently internally
    """
    block_results = []
    dmr_summary = []

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
                            left_idx -= 1  # **直接向左扩展**
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
                wbr_result = wbr.run_weighted_beta_regression(sub_block, group1=group1, group2=group2)

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
                        wbr_result['DMR']
                    ])
    dmr_summary_result = pd.DataFrame(dmr_summary, columns=[
        "chromosome", "start", "end", "count", f"{group1}_mean", f"{group2}_mean", "delta", "pvalue", "F", "DMR"
    ])

    ### BH multiple test correction
    dmr_blocks_with_padj, significant_dmr = BH_adjust.adjust_p_values_v2(dmr_summary_result, group1=group1, group2=group2)
    return dmr_blocks_with_padj, significant_dmr

if __name__ == "__main__":
    out_data = pd.DataFrame({
        "Chromosome": ["chr1"] * 25,
        "Position": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500,
                     2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900],
        "north_mean": [0.2, 0.3, 0.25, 0.4, 0.35, 0.38, 0.37, 0.1, 0.15, 0.13, 0.2, 0.22, 0.15, 0.33, 0.44,
                       0.1, 0.12, 0.15, 0.18, 0.2, 0.22, 0.25, 0.27, 0.30, 0.33],
        "xizang_mean": [0.1, 0.2, 0.15, 0.3, 0.25, 0.28, 0.27, 0.3, 0.35, 0.38, 0.40, 0.42, 0.45, 0.48, 0.50,
                        0.05, 0.08, 0.1, 0.12, 0.15, 0.18, 0.20, 0.23, 0.25, 0.28],
        "north_var": [0.01] * 25,
        "xizang_var": [0.01] * 25,
        "north_cov": [10] * 25,
        "xizang_cov": [10] * 25,
        "mean_diff": [0.1, 0.1, 0.2, 0.3, 0.15, 0.4, 0.1,
                      -0.1, -0.2, -0.3, -0.15, -0.01, -0.01, -0.01, -0.01,
                      -0.01, -0.01, -0.03, -0.1, -0.07,
                      0.01, 0.01, 0.2, 0.1, 0.01],
        "Block": ["block0"] * 7 + ["block1"] * 8 + ["block2"] * 5 + ["block3"] * 5
    })

    dmr_blocks_with_padj, significant_dmr = find_blocks_greedy(out_data,chr_col="chr1",group1="north",group2="xizang")
    print(dmr_blocks_with_padj)
