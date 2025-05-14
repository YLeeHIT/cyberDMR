import pandas as pd
import os
from typing import List, Dict, Optional, Set, Union
from concurrent.futures import ThreadPoolExecutor
import low_coverage as filling

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

