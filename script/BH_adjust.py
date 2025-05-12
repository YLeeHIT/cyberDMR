import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

def generate_simulated_dmr_data(num_rows=1000, p_value_threshold=0.05, group1="g1", group2="g2"):
    """
    Generate simulated DMR data where all p-values are less than a specified threshold, which defaults to 0.05.
    """
    np.random.seed(42)  # Make it reproducible
    
    # Chromosome (randomly selected from chromosomes 1-22).
    chromosomes = np.random.choice([f'chr{i}' for i in range(1, 23)], num_rows)

    # Generate start (start) and end (end) positions within the range of 1M to 100M.
    starts = np.random.randint(1_000_000, 100_000_000, num_rows)
    ends = starts + np.random.randint(200, 1000, num_rows)  # ensure end > start

    # Generate the number of CpGs (greater than 5).
    cpg_counts = np.random.randint(6, 50, num_rows)

    # 生成平均甲基Generate mean methylation levels (balanced data).
    g1_values = np.random.uniform(0.5, 1, num_rows)
    g2_values = np.random.uniform(0.3, 0.5, num_rows)

    # Generate a methylation level difference (greater than 0.1 or less than -0.1).
    meth_diff = np.random.choice(np.concatenate((np.random.uniform(-1, -0.1, num_rows//2), 
                                                 np.random.uniform(0.1, 1, num_rows//2))), num_rows)

    # Generate p-values (all less than p_value_threshold).
    p_values = np.random.uniform(0, p_value_threshold, num_rows)

    # Generate block ID.
    block_ids = [f"block_{i+1}" for i in range(num_rows)]

    # Generate F-statistic.
    f_value = np.random.uniform(5, 10, num_rows)

    # Generate DMR column.
    dmr_value = [True]*num_rows
    # Construct DataFrame.
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
    Perform Benjamini-Hochberg (BH) correction on the given DMR data and return two DataFrames: 
    1. dmr_data_with_padj：This DataFrame contains the original p-values and the adjusted p-values.
    2. significant_dmr_data：This DataFrame retains only rows where the adjusted p-value is less than 0.05, and it includes only the adjusted p-values.
    """

    # Null data processing
    if df.empty:
        empty_cols = ['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean",
                      'delta', 'F', 'pvalue', 'p_adj', 'DMR']
        empty_df = pd.DataFrame(columns=empty_cols)
        return empty_df.copy(), empty_df.copy()

    # BH correction
    df = df.copy()
    df['p_adj'] = multipletests(df['pvalue'], method='fdr_bh')[1]


    # Keep 4 decimal places
    for col in ['count', f"{group1}_mean", f"{group2}_mean", 'delta', 'F']:
        df[col] = df[col].round(4)
    
    # Use scientific notation for a concise representation of pvalue and p_adj
    df['pvalue'] = df['pvalue'].apply(lambda x: float(f"{x:.4g}"))
    df['p_adj'] = df['p_adj'].apply(lambda x: float(f"{x:.4g}"))

    # Create complete DMR data, including p_value and p_adj
    dmr_data_with_padj = df[['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean", 'delta', 'F', 'pvalue', 'p_adj', 'DMR']].copy()

    # Filter for DMRs where the adjusted p-value is less than 0.05, keeping only the p_adj column
    significant_dmr_data = dmr_data_with_padj[dmr_data_with_padj['p_adj'] < 0.05][['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean", 'delta', 'F','pvalue', 'p_adj', 'DMR']].copy()

    return dmr_data_with_padj, significant_dmr_data

if __name__ == "__main__":

    # Perform BH correction and adjust the output format
    dmr_data = generate_simulated_dmr_data(group1="xizang",group2="north")
    print(dmr_data.head(10))
    dmr_data_with_padj, significant_dmr_data = adjust_p_values_v2(dmr_data,group1="xizang",group2="north")

    print(significant_dmr_data)

    dmr_data_with_padj.to_csv("/home/ly/shell/deepDMR/data/simulate_data/BH_simulate_result.txt",sep='\t',index=False)