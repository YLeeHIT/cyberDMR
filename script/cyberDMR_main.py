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
    解析命令行参数，包括输出目录、线程数、实验组和对照组名称。
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
    对单条染色体进行DMR分析流程。
    - in_chr: 染色体名称 如 "chr1" 
    - out_dir: 输出目录
    - group1/group2: 分组名称
    - threads: 内部并行线程数量（主要用于数据读取）
    """
    print(f"Processing {in_chr}...")

    infile = f"{out_dir}/in_cyber.lab" # 假设输入文件为统一的大文件，按染色体筛选
    sample_data = reading.load_inlab_chrwise(
        infile,
        target_chr=in_chr,
        column_names=["Chr", "Pos", "Meth_Level", "Coverage"],
        num_threads=threads
    )

    if not sample_data:
        print(f"No data found for {in_chr}, skipping.")
        return

    # 提取样本与分组信息
    label = pd.DataFrame(
        [(entry['sample'], entry['group']) for entry in sample_data],
        columns=['sample','group']
    )

    # Step 1: 填补、合并数据
    processed_samples = filling.process_samples(
        sample_data, coverage_threshold=5, max_distance=500
    )
    merged_data = filling.merge_samples_fast(processed_samples)

    # Step 2: CpG分块
    CpG_distance = 500
    CpG_count = 5
    out_data, block_ranges = blocking.process_data(
        merged_data, label, group1, group2, CpG_distance, CpG_count
    )
    #out_data.to_csv(os.path.join(out_dir, f"{in_chr}_step2.tsv"), sep="\t", header=True, index=False)

    # Step 3: 聚类和DMR检测
    dmr_data_with_padj, significant_dmr_data = clustering.find_blocks_greedy(
        out_data, delta_m_mean_threshold=0.1,
        group1=group1, group2=group2,
        chr_col=in_chr
    )
    dmr_data_with_padj.to_csv(os.path.join(out_dir, f"{in_chr}_cyberDMR.txt"), sep="\t", header=True, index=False)

    print(f"Finished processing {in_chr}.")

def main():
    """
    主控制函数：
    - 解析参数
    - 多线程分发染色体任务
    """
    args = parse_args()

    # 染色体列表，可根据实际需求调整（或支持 --chroms 参数）
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    # 多进程处理不同染色体
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        for chrom in chroms:
            executor.submit(
                process_one_chromosome,
                chrom, args.out_dir, args.group1, args.group2, 1 # 每个进程内部单线程读取
            )

   # out_dir = "/home/ly/shell/deepDMR/data/real_data_prostate/formatted_cyberDMR/cyberDMR_result"
   # group1 = "lethal"
   # group2 = "normal"
   # threads = 1
   # chroms = ["chr22"]
   # #chrom = "chr10" 
   # #process_one_chromosome(chrom, out_dir, group1, group2, 1)
   # with ProcessPoolExecutor(max_workers=threads) as executor:
   #     for chrom in chroms:
   #         executor.submit(
   #             process_one_chromosome,
   #             chrom, out_dir, group1, group2, 1 # 每个进程内部单线程读取
   #         )

if __name__ == "__main__":
    main()