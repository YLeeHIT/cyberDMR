import pandas as pd
import argparse
import os
from concurrent.futures import ProcessPoolExecutor

# my model
import read_samples as reading
import low_coverage_v2 as filling
import cpg_blocking_v2 as blocking
import cpg_clustering_v5 as clustering

def parse_args():
    """
    Parse command-line arguments, including output directory, number of threads, and names for the experimental and control groups
    """
    parser = argparse.ArgumentParser(description="CyberDMR analysis for all chromosomes")
    parser.add_argument(
        "--out_dir", type=str, required=True,
        help="The directory for output results"
    )
    parser.add_argument(
        "--threads", type=int, default=1,
        help="The number of parallel threads to use"
    )
    parser.add_argument(
        "--group1", type=str, required=True,
        help="Name of group 1 (e.g., control)"
    )
    parser.add_argument(
        "--group2", type=str, required=True,
        help="Name of group 2 (e.g., treatment)"
    )
    return parser.parse_args()

def process_one_chromosome(in_chr, out_dir, group1, group2, threads=1):
    """
    Perform the DMR analysis workflow for a single chromosome
    - in_chr: Chromosome name, e.g., "chr1" 
    - out_dir: Output directory
    - group1/group2: Group names
    - threads: Number of internal parallel threads (mainly for data reading)
    """
    print(f"Processing {in_chr}...")

    infile = f"{out_dir}/in_cyber.lab" # Assuming the input is a unified large file, filter by chromosome
    sample_data = reading.load_inlab_chrwise(
        infile,
        target_chr=in_chr,
        column_names=["Chr", "Pos", "Meth_Level", "Coverage"],
        num_threads=threads
    )

    if not sample_data:
        print(f"No data found for {in_chr}, skipping.")
        return

    # Extract sample and group information
    label = pd.DataFrame(
        [(entry['sample'], entry['group']) for entry in sample_data],
        columns=['sample','group']
    )

    # Step 1: Impute and merge data
    processed_samples = filling.process_samples(
        sample_data, coverage_threshold=5, max_distance=500
    )
    merged_data = filling.merge_samples_fast(processed_samples)

    # Step 2: CpG blocking
    CpG_distance = 500
    CpG_count = 5
    out_data, block_ranges = blocking.process_data(
        merged_data, label, group1, group2, CpG_distance, CpG_count
    )
    
    # Step 3: Clustering and DMR detection
    dmr_data_with_padj, significant_dmr_data = clustering.find_blocks_greedy(
        out_data, delta_m_mean_threshold=0.1,
        group1=group1, group2=group2,
        chr_col=in_chr
    )
    dmr_data_with_padj.to_csv(os.path.join(out_dir, f"{in_chr}_cyberDMR.txt"), sep="\t", header=True, index=False)

    print(f"Finished processing {in_chr}.")

def main():
    """
    Main control function:
    - Parse parameters.
    - Dispatch chromosome tasks using multiple threads
    """
    args = parse_args()

    # Chromosome list, adjustable according to actual needs (or supports a --chroms parameter)
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    # Process different chromosomes using multiple processes
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        for chrom in chroms:
            executor.submit(
                process_one_chromosome,
                chrom, args.out_dir, args.group1, args.group2, 1 # Single-threaded reading within each process
            )


if __name__ == "__main__":
    main()