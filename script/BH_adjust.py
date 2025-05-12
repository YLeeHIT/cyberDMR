import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

def generate_simulated_dmr_data(num_rows=1000, p_value_threshold=0.05, group1="g1", group2="g2"):
    """
    生成模拟的DMR数据 其中p值全部小于指定阈值 默认为0.05。
    """
    np.random.seed(42)  # 保证可复现
    
    # 染色体（随机选择1-22号染色体）
    chromosomes = np.random.choice([f'chr{i}' for i in range(1, 23)], num_rows)

    # 生成起始（start）和终止（end）位点，范围在 1M-100M 之间
    starts = np.random.randint(1_000_000, 100_000_000, num_rows)
    ends = starts + np.random.randint(200, 1000, num_rows)  # 确保end > start

    # 生成CpG数量（大于5）
    cpg_counts = np.random.randint(6, 50, num_rows)

    # 生成平均甲基化水平（配平数据）
    g1_values = np.random.uniform(0.5, 1, num_rows)
    g2_values = np.random.uniform(0.3, 0.5, num_rows)
    #setattr(globals(),f"{group1}_mean", g1_values) 
    #setattr(globals(),f"{group2}_mean", g2_values)

    # 生成方差
    #g1_var = [0.01]*num_rows
    #g2_var = [0.01]*num_rows
    #setattr(globals(),f"{group1}_var",g1_var)
    #setattr(globals(),f"{group2}_var",g2_var)

    # 生成覆盖度
    #g1_cov = [10]*num_rows
    #g2_cov = [10]*num_rows
    #setattr(globals(),f"{group1}_cov",g1_cov)
    #setattr(globals(),f"{group2}_cov",g2_cov)


    # 生成甲基化水平的差值（大于0.1或者小于-0.1）
    meth_diff = np.random.choice(np.concatenate((np.random.uniform(-1, -0.1, num_rows//2), 
                                                 np.random.uniform(0.1, 1, num_rows//2))), num_rows)

    # 生成p值（全部小于p_value_threshold）
    p_values = np.random.uniform(0, p_value_threshold, num_rows)

    # 生成block编号
    block_ids = [f"block_{i+1}" for i in range(num_rows)]

    # 生成 F统计量
    f_value = np.random.uniform(5, 10, num_rows)

    # 生成DMR列
    dmr_value = [True]*num_rows
    # 构建DataFrame
    df = pd.DataFrame({
        'chromosome': chromosomes,
        'start': starts,
        'end': ends,
        'count': cpg_counts,
        f"{group1}_mean": g1_values,
        f"{group2}_mean": g2_values,
        'delta': meth_diff,
        'pvalue': p_values,
        #'Block': block_ids
        'F': f_value,
        'DMR': dmr_value
    })
    
    return df



def adjust_p_values_v2(df, group1="g1", group2="g2"):
    """
    对给定的DMR数据进行BH矫正 返回两个DataFrame 
    1. dmr_data_with_padj 包含原始p值和矫正后的p值
    2. significant_dmr_data 仅保留矫正后p值小于0.05的行 并只包含矫正后的p值
    """

    # 空数据处理
    if df.empty:
        empty_cols = ['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean",
                      'delta', 'F', 'pvalue', 'p_adj', 'DMR']
        empty_df = pd.DataFrame(columns=empty_cols)
        return empty_df.copy(), empty_df.copy()

    # BH矫正
    df = df.copy()
    df['p_adj'] = multipletests(df['pvalue'], method='fdr_bh')[1]


    # 保留4位小数
    for col in ['count', f"{group1}_mean", f"{group2}_mean", 'delta', 'F']:
        df[col] = df[col].round(4)
    
    # 对 pvalue 和 p_adj 使用科学记数法简洁表示
    df['pvalue'] = df['pvalue'].apply(lambda x: float(f"{x:.4g}"))
    df['p_adj'] = df['p_adj'].apply(lambda x: float(f"{x:.4g}"))

    # 创建完整的DMR数据，包含p_value和p_adj
    dmr_data_with_padj = df[['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean", 'delta', 'F', 'pvalue', 'p_adj', 'DMR']].copy()

    # 筛选矫正后p值小于0.05的DMR，只保留p_adj列
    significant_dmr_data = dmr_data_with_padj[dmr_data_with_padj['p_adj'] < 0.05][['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean", 'delta', 'F','pvalue', 'p_adj', 'DMR']].copy()

    return dmr_data_with_padj, significant_dmr_data

if __name__ == "__main__":

    # 进行BH矫正并调整输出格式
    dmr_data = generate_simulated_dmr_data(group1="xizang",group2="north")
    print(dmr_data.head(10))
    dmr_data_with_padj, significant_dmr_data = adjust_p_values_v2(dmr_data,group1="xizang",group2="north")

    print(significant_dmr_data)

    dmr_data_with_padj.to_csv("/home/ly/shell/deepDMR/data/simulate_data/BH_simulate_result.txt",sep='\t',index=False)