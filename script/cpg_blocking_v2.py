import pandas as pd
import numpy as np
from itertools import combinations
from itertools import islice

def display_top_blocks(block_ranges, num):
    """
    显示 block_ranges 字典中前 num 个元素
    
    参数:
    block_ranges (dict): 字典，包含 block 名称和对应的行索引
    num (int): 要展示的 block 数量
    
    返回:
    dict: 包含前 num 个 block 的字典
    """
    # 获取前 num 个元素
    top_blocks = dict(islice(block_ranges.items(), num))

    # 打印前 num 个元素
    print(f"First {num} blocks:")
    for block, indices in top_blocks.items():
        print(f"{block}: {indices}")
    
def process_data(data, label, group1, group2, CpG_distance=500, CpG_count=5):
    # 提取样本名和分组信息
    sample_names = label.iloc[:, 0].values  # 第一列样本名
    groups = label.iloc[:, 1].values        # 第二列分组信息

    # 提取染色体和位置
    chromosome = data.iloc[:, 0].values
    position = data.iloc[:, 1].values.astype(int)
    
    # 提取样本数据
    sample_data = data.iloc[:, 2:]

    # 计算group1和group2 进行collapse处理
    selected_groups = [group1, group2]
    group_means = sample_data.T.groupby(groups).mean().T[selected_groups]
    group_variances = sample_data.T.groupby(groups).var().T[selected_groups]
    sample_cov = sample_data.T.groupby(groups).count().T[selected_groups]
    

    # 过滤掉 group1 和 group2 覆盖度cov为0的点
    #valid_rows = (sample_cov[group1] > 0) & (sample_cov[group2] > 0)
    valid_rows = sample_cov[selected_groups].min(axis=1) > 0

    # 创建压缩后统计数据列表
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

    # 计算均值差异
    data2["mean_diff"] = (data2[f"{group1}_mean"] - data2[f"{group2}_mean"]).round(4)

    # 计算相邻位置的距离差
    distance_diff = np.diff(data2['Position'], prepend=data2['Position'].iloc[0])

    # 计算均值差的符号变化
    mean_diff_sign = np.sign(data2[f"mean_diff"])
    mean_diff_sign_change = np.diff(mean_diff_sign, prepend=mean_diff_sign.iloc[0])

    # 识别新块的起点
    new_block_indices = (distance_diff >= CpG_distance) | (mean_diff_sign_change != 0)

    # 计算 block 编号
    #blocks = np.cumsum(new_block_indices)

    # 存入 DataFrame
    #data2['Block'] = [f"block{i}" for i in blocks]
    data2['Block'] = data2.groupby(new_block_indices.cumsum()).ngroup().apply(lambda x: f"block{x}")


    # 计算每个 block 的大小，并筛选掉小于 CpG_count 的 block
    block_counts = data2['Block'].value_counts()
    blocks_to_keep = block_counts[block_counts >= CpG_count].index
    data2_filtered = data2[data2['Block'].isin(blocks_to_keep)]

    # 存储块的索引，方便后续并行化以及提取
    #block_ranges = {}
    #for block, group in data2_filtered.groupby('Block'):
    #    start_idx = group.index.min()  # 该 block 的起始索引
    #    end_idx = group.index.max()    # 该 block 的结束索引
    #    block_ranges[block] = (start_idx, end_idx)
    
    block_ranges = {block: (group.index.min(), group.index.max()) for block, group in data2_filtered.groupby('Block')}


    return data2_filtered, block_ranges


if __name__ == "__main__":
    # 示例使用
    in_data = pd.read_csv("/home/ly/shell/deepDMR/data/mergeData_fromONT_bedtools/sam-han-zang_num-20_head-1w_merge_NA.bed", sep="\t", header=None)
    label = pd.read_csv("/home/ly/shell/deepDMR/data/mergeData_fromONT_bedtools/lab.txt", sep="\t", header=None)
    CpG_distance = 500  # 设定距离阈值
    CpG_count = 5  # 设定最小出现次数
    group1 = "north"
    group2 = "xizang"

    # 调用函数
    out_data, block_ranges = process_data(in_data, label, group1, group2, CpG_distance, CpG_count)

    # 查看返回的结果
    print(out_data.head(20))
    #print(block_ranges)
    display_top_blocks(block_ranges, 10)

    out_data.to_csv("/home/ly/shell/deepDMR/data/mergeData_fromONT_bedtools/sam-han-zang_num-20_head-1w_merge_stepOne_notrend.out2",sep="\t",header=True,index=False)

