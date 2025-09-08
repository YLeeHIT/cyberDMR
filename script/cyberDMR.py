import pandas as pd
import argparse
import os
from concurrent.futures import ProcessPoolExecutor

# my model
import read_samples as reading
import low_coverage as filling
import cpg_blocking as blocking
import cpg_clustering as clustering

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for cyberDMR.

    Required arguments:
    --out-dir   Path to the output directory.
    --group1    Group1 label.
    --group2    Group2 label.

    Optional arguments:
    -t,    --threads        Number of worker processes (default: 8).
    -chr,  --chroms         Chromosome set specification (default: '1-22,X,Y').
    -d,    --delta          Delta m-mean threshold for greedy DMR detection (default: 0.1).
    -bdis, --cpg-distance   Maximum CpG distance used in blocking (default: 500).
    -ct,   --cpg-count      Minimum number of CpGs per block (default: 5).
    -cov,  --min-cov        Minimum CpG coverage to keep (default: 5).
    -fdis, --max-dist       Maximum distance of adjacent CpGs (default: 500).
    -lab,  --cyber-lab      Path to the cyber.lab file (optional, default: None).
    -q,    --qvalue         BH-corrected p-value threshold for significant DMRs.
    -f,    --Fvalue         F statistic
    """
    
    parser = argparse.ArgumentParser(prog="cyberDMR", description="Detect DMRs with cyberDMR.")

    # ---------------- Required arguments ----------------
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--out-dir", "-o", required=True, metavar="PATH", help="Output directory.")
    required.add_argument("--group1", "-g1", required=True, metavar="STR", help="Group1 label.")
    required.add_argument("--group2", "-g2", required=True, metavar="STR", help="Group2 label.")

    # ---------------- Optional arguments ----------------
    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument("--threads", "-t", type=int, default=8, metavar="INT", help="Number of worker processes (default: 8).")
    optional.add_argument("--chroms", "-chr", default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY", metavar="STR", help="Chromosome set, e.g. 'chr1-chr22,chrX,chrY' or 'chr22' (default: chr1-chr22,chrX,chrY).")
    optional.add_argument("--delta", "-d", type=float, default=0.1, metavar="FLOAT", help="Methylation difference (delta) threshold (default: 0.1).")
    optional.add_argument("--cpg-distance", "-bdis", type=int, default=500, metavar="INT", help="Max CpG distance for blocking (default: 500).")
    optional.add_argument("--cpg-count", "-ct", type=int, default=5, metavar="INT", help="Min CpG count per block (default: 5).")
    optional.add_argument("--min-cov", "-cov", type=int, default=5, metavar="INT", help="Minimum CpG coverage for filling (default: 5).")
    optional.add_argument("--max-dist", "-fdis", type=int, default=500, metavar="INT", help="Maximum distance (bp) of adjacent CpGs for filling (default: 500).")
    optional.add_argument("--cyber-lab", "-lab",metavar="PATH", help="Path to the cyber.lab file.")
    optional.add_argument("--qvalue", "-q", type=float, default=0.05, metavar="FLOAT", help="BH-corrected p-value threshold for significant DMRs (default: 0.05).")
    optional.add_argument("--Fvalue", "-f", type=float, default=15, metavar="FLOAT", help="F statistic (default: 15).")

    return parser.parse_args()

def process_one_chromosome(in_chr, out_dir, group1, group2, 
                           threads=1, coverage_threshold=5, max_distance=500, CpG_distance=500, CpG_count=5, delta_m_mean_threshold=0.1,
                           qvalue=0.05, Fvalue=15, cyber_lab=None):
    
    """
    Process a single chromosome: fill missing values, merge samples, perform CpG blocking, and detect DMRs.        
    
    Arguments:
    - in_chr: Chromosome name (e.g., 'chr1')
    - out_dir: Output directory where results will be saved
    - group1: Group1 label (e.g., 'control')
    - group2: Group2 label (e.g., 'treatment')
    - threads: Number of threads to use
    - coverage_threshold: Minimum coverage threshold for imputation (default: 5)
    - max_distance: Maximum distance for imputation (default: 500)
    - CpG_distance: Maximum distance between CpGs for blocking (default: 500)
    - CpG_count: Minimum number of CpGs per block (default: 5)
    - delta_m_mean_threshold: Delta methylation mean threshold for DMR detection (default: 0.1)
    - cyber_lab: Optional path to a custom cyber.lab file
    - qvalue: BH-corrected p-value
    - Fvalue: F statitic
    """
    
    print(f"Processing {in_chr}...")

    if cyber_lab:
        infile = cyber_lab
    else:
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
    processed_samples = filling.process_samples(sample_data, coverage_threshold=coverage_threshold, max_distance=max_distance)
    merged_data = filling.merge_samples_fast(processed_samples)
    #merged_data.to_csv(os.path.join(out_dir, f"{in_chr}_merged_data_after_filling.txt"), sep="\t", header=True, index=False)

    # Step 2: CpG blocking
    out_data, block_ranges = blocking.process_data(merged_data, label, group1, group2, CpG_distance=CpG_distance, CpG_count=CpG_count)
    
    # Step 3: Clustering and DMR detection
    dmr_data_with_padj, significant_dmr_data = clustering.find_blocks_greedy(
        out_data, 
        delta_m_mean_threshold=delta_m_mean_threshold,
        group1=group1, group2=group2,
        qvalue=qvalue,
        chr_col=in_chr,
        Fvalue=Fvalue
    )
    #dmr_data_with_padj.to_csv(os.path.join(out_dir, f"{in_chr}_cyberDMR.txt"), sep="\t", header=True, index=False)
    significant_dmr_data.to_csv(os.path.join(out_dir, f"{in_chr}_cyberDMR.txt"), sep="\t", header=True, index=False)

    print(f"Finished processing {in_chr}.")


def test():
    # Test parameters
    chrom = "chr22"
    outdir = "/home/user/liyang/project/2fold/Homo_sapien/exp/pDMR/cyberDMR/result1"
    group1 = "Chinese"
    group2 = "Jew"
    
    # Test additional parameters (using default values or as needed)
    threads = 8
    coverage_threshold = 10  # Example value
    max_distance = 1000  # Example value
    CpG_distance = 500  # Example value
    CpG_count = 7  # Example value
    delta_m_mean_threshold = 0.2  # Example value
    
    # Calling process_one_chromosome with all parameters
    process_one_chromosome(
        chrom, outdir, group1, group2, threads,
        coverage_threshold=coverage_threshold,
        max_distance=max_distance,
        CpG_distance=CpG_distance,
        CpG_count=CpG_count,
        delta_m_mean_threshold=delta_m_mean_threshold
        )


def main():
    """
    Main control function:
    - Parse parameters.
    - Dispatch chromosome tasks using multiple threads
    """
    args = parse_args()
    print(f"parameters is {args}")
    # If a file path is provided, load chroms from the file. 
    # Otherwise, parse chroms directly from the input (either comma-separated or "1-22,X,Y")
    chroms = []
    if args.chroms:  # If chroms parameter is provided
        if "," in args.chroms:
            # Split the chroms string into a list (e.g., "chr1,chr2,chr3")
            chroms = [chrom if chrom.startswith("chr") else f"chr{chrom}" for chrom in args.chroms.split(",")]
        elif "-" in args.chroms:
            # Handle range input, like "1-22"
            start, end = args.chroms.split("-")
            chroms = [f"chr{i}" for i in range(int(start), int(end) + 1)]
        else:
            # Single chromosome input (like "1" or "22"), ensure it's prefixed with "chr"
            chroms = [f"chr{args.chroms}" if not args.chroms.startswith("chr") else args.chroms]
    else:
        # If no chroms parameter is provided, default to chr1 to chr22 and chrX, chrY
        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


    # Process different chromosomes using multiple processes
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        for chrom in chroms:
            executor.submit(
                process_one_chromosome,
                chrom, args.out_dir, args.group1, args.group2, 1,
                coverage_threshold=args.min_cov,
                max_distance=args.max_dist,
                CpG_distance=args.cpg_distance,
                CpG_count=args.cpg_count,
                delta_m_mean_threshold=args.delta,
                cyber_lab=args.cyber_lab,
                qvalue=args.qvalue,
                Fvalue=args.Fvalue
            )

if __name__ == "__main__":
    main()
    #test()
