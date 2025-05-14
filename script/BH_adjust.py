import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

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

    if df.empty:
        empty_cols = ['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean",
                      'delta', 'F', 'pvalue', 'p_adj', 'DMR']
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

    dmr_data_with_padj = df[['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean", 'delta', 'F', 'pvalue', 'p_adj', 'DMR']].copy()
    significant_dmr_data = dmr_data_with_padj[dmr_data_with_padj['p_adj'] < 0.05][['chromosome', 'start', 'end', 'count', f"{group1}_mean", f"{group2}_mean", 'delta', 'F','pvalue', 'p_adj', 'DMR']].copy()

    return dmr_data_with_padj, significant_dmr_data

if __name__ == "__main__":
    # Perform BH correction and adjust the output format
    dmr_data = generate_simulated_dmr_data(group1="xizang",group2="north")
    dmr_data_with_padj, significant_dmr_data = adjust_p_values_v2(dmr_data,group1="xizang",group2="north")
    print(significant_dmr_data)
