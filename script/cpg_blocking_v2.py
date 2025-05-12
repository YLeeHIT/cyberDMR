import pandas as pd
import numpy as np
from itertools import combinations
from itertools import islice

def display_top_blocks(block_ranges, num):
    """
    Displays the first num items (blocks) from the block_ranges dictionary
    
    Parameters:
    block_ranges (dict): A dictionary where keys are block names and values are their corresponding row indices.
    num (int): The number of blocks (items) to display.
    
    Returns:
    dict: A dictionary containing the first num items (blocks) from block_ranges
    """
    # Get the first num elements
    top_blocks = dict(islice(block_ranges.items(), num))

    # Print the first num elements
    print(f"First {num} blocks:")
    for block, indices in top_blocks.items():
        print(f"{block}: {indices}")
    
def process_data(data, label, group1, group2, CpG_distance=500, CpG_count=5):
    # Extract sample names and group information
    sample_names = label.iloc[:, 0].values  # First column: sample names
    groups = label.iloc[:, 1].values        # Second column: group information

    # Extract chromosome and position
    chromosome = data.iloc[:, 0].values
    position = data.iloc[:, 1].values.astype(int)
    
    # Extract sample data
    sample_data = data.iloc[:, 2:]

    # Calculate for group1 and group2 and perform collapse processing
    selected_groups = [group1, group2]
    group_means = sample_data.T.groupby(groups).mean().T[selected_groups]
    group_variances = sample_data.T.groupby(groups).var().T[selected_groups]
    sample_cov = sample_data.T.groupby(groups).count().T[selected_groups]
    

    # Filter out data points for group1 and group2 where coverage (cov) is 0
    valid_rows = sample_cov[selected_groups].min(axis=1) > 0

    # Create a list of compressed statistical data
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

    # Calculate the mean difference
    data2["mean_diff"] = (data2[f"{group1}_mean"] - data2[f"{group2}_mean"]).round(4)

    # Calculate the distance difference between adjacent positions
    distance_diff = np.diff(data2['Position'], prepend=data2['Position'].iloc[0])

    # Calculate the sign change of the mean difference
    mean_diff_sign = np.sign(data2[f"mean_diff"])
    mean_diff_sign_change = np.diff(mean_diff_sign, prepend=mean_diff_sign.iloc[0])

    # Identify the starting point of new blocks
    new_block_indices = (distance_diff >= CpG_distance) | (mean_diff_sign_change != 0)


    # Store into a DataFrame
    data2['Block'] = data2.groupby(new_block_indices.cumsum()).ngroup().apply(lambda x: f"block{x}")


    # Calculate the size of each block and filter out blocks smaller than CpG_count
    block_counts = data2['Block'].value_counts()
    blocks_to_keep = block_counts[block_counts >= CpG_count].index
    data2_filtered = data2[data2['Block'].isin(blocks_to_keep)]

    
    block_ranges = {block: (group.index.min(), group.index.max()) for block, group in data2_filtered.groupby('Block')}


    return data2_filtered, block_ranges


if __name__ == "__main__":
    # Example usage
    in_data = pd.read_csv("/home/ly/shell/deepDMR/data/mergeData_fromONT_bedtools/sam-han-zang_num-20_head-1w_merge_NA.bed", sep="\t", header=None)
    label = pd.read_csv("/home/ly/shell/deepDMR/data/mergeData_fromONT_bedtools/lab.txt", sep="\t", header=None)
    CpG_distance = 500  # Set the distance threshold
    CpG_count = 5  # Set the minimum number of occurrences
    group1 = "north"
    group2 = "xizang"

    # Call the function
    out_data, block_ranges = process_data(in_data, label, group1, group2, CpG_distance, CpG_count)

    # View the returned result
    print(out_data.head(20))
    display_top_blocks(block_ranges, 10)

    out_data.to_csv("/home/ly/shell/deepDMR/data/mergeData_fromONT_bedtools/sam-han-zang_num-20_head-1w_merge_stepOne_notrend.out2",sep="\t",header=True,index=False)

