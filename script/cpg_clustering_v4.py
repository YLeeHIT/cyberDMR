import heapq
import numpy as np
import pandas as pd
import cpg_blocking_v2 as blocking

def calculate_score(delta_m, distance, d_max, alpha=1.0, beta=0.1):
    """ 计算扩展方向的优先级得分 """
    return alpha * abs(delta_m) - beta * (distance / d_max)

def find_blocks_greedy(data2_filtered, block_ranges, group1, group2,delta_m_mean_threshold=10, min_cpg_count=5, alpha=1.0, beta=0.5):
    """
    采用贪心策略优化并更新 CpG block，不合并 block
    """
    # 计算最大距离 d_max
    #d_max = data2_filtered["Position"].max() - data2_filtered["Position"].min()

    # 用于存储新的 block 结构（但不合并）
    new_block_ranges = {}

    for block, (start_idx, end_idx) in block_ranges.items():
        # 提取当前 Block 的数据
        block_data = data2_filtered.iloc[start_idx:end_idx + 1].copy()
        print(block_data)
        
        # 计算块内最大位置和最小位置
        d_max = block_data["Position"].max() - block_data["Position"].min()
        
        # 采用最大堆（按 ΔM 绝对值排序）
        heap = [(-abs(delta_m), idx) for idx, delta_m in zip(block_data.index, block_data[f"{group1}_{group2}_mean_diff"])]
        heapq.heapify(heap)

        updated_block = set()
        while heap:
            _, seed_idx = heapq.heappop(heap)
            if seed_idx in updated_block:
                continue
            
            seed_delta_m = block_data.loc[seed_idx, f"{group1}_{group2}_mean_diff"]

            if abs(seed_delta_m) < delta_m_mean_threshold:
                break

            current_block = {seed_idx}
            updated_block.add(seed_idx)

            left_idx, right_idx = seed_idx - 1, seed_idx + 1

            while left_idx >= start_idx or right_idx <= end_idx:
                best_direction = None
                left_score, right_score = -np.inf, -np.inf

                # 限制范围，防止超出 block_ranges
                if left_idx >= start_idx and left_idx not in updated_block:
                    left_delta_m = block_data.loc[left_idx, f"{group1}_{group2}_mean_diff"]
                    left_distance = abs(block_data.loc[left_idx, "Position"] - block_data.loc[start_idx, "Position"])
                    left_score = calculate_score(left_delta_m, left_distance, d_max, alpha, beta)

                if right_idx <= end_idx and right_idx not in updated_block:
                    right_delta_m = block_data.loc[right_idx, f"{group1}_{group2}_mean_diff"]
                    right_distance = abs(block_data.loc[right_idx, "Position"] - block_data.loc[end_idx, "Position"])
                    right_score = calculate_score(right_delta_m, right_distance, d_max, alpha, beta)

                if left_score > right_score:
                    best_direction = "left"
                elif right_score > left_score:
                    best_direction = "right"

                if best_direction == "left":
                    current_block.add(left_idx)
                    updated_block.add(left_idx)
                    left_idx -= 1
                elif best_direction == "right":
                    current_block.add(right_idx)
                    updated_block.add(right_idx)
                    right_idx += 1
                else:
                    break  

            if len(current_block) >= min_cpg_count:
                new_block_ranges[f"block{len(new_block_ranges)}"] = (min(current_block), max(current_block))

    # 更新 `data2_filtered` 的 Block 列
    block_map = {idx: block for block, (start, end) in new_block_ranges.items() for idx in range(start, end + 1)}
    data2_filtered["Block"] = data2_filtered.index.map(block_map)

    return data2_filtered, new_block_ranges

in_data = pd.read_csv("/home/ly/shell/deepDMR/data/mergeData_fromONT_bedtools/sam-han-zang_num-20_head-1w_merge_NA.bed", sep="\t", header=None)
label = pd.read_csv("/home/ly/shell/deepDMR/data/mergeData_fromONT_bedtools/lab.txt", sep="\t", header=None)
CpG_distance = 500  # 设定距离阈值
CpG_count = 5  # 设定最小出现次数

    # 调用函数
out_data, block_ranges = blocking.process_data(in_data, label, "north", "xizang", CpG_distance, CpG_count)

print(out_data.head(10))

# 运行算法
data2_filtered, block_ranges = find_blocks_greedy(out_data, block_ranges, "north", "xizang")

# 打印新的 `Block` 结果
#for block, (start, end) in block_ranges.items():
#    print(f"{block}: Start={start}, End={end}")