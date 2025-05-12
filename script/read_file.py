import pandas as pd

chr_whitelist = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

def chr_sort_key(chrom):
    if chrom == "chrX":
        return 23
    elif chrom == "chrY":
        return 24
    else:
        return int(chrom[3:])

#print(sorted(chr_whitelist, key=chr_sort_key))

def read_file_by_chr_list(path, chr_whitelist):
    lines = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            chrom = line.strip().split('\t')[0]
            if chr_whitelist is None or chrom in chr_whitelist:
                lines.append(line)
    return pd.DataFrame(lines)

def read_file_by_chr(path, target_chr):
    rows = []
    with open(path, 'r', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if not parts:
                continue
            chrom = parts[0]
            if chrom == target_chr:
                rows.append(parts)
    return pd.DataFrame(rows)


def load_inlab_chrwise(inlab_path, target_chr, column_names=["Chr","Pos","Meth_Level","Coverage"]):
    df_meta = pd.read_csv(inlab_path, sep='\t', header=None, names=["sample", "group", "filepath"])

    results = []

    for _, row in df_meta.iterrows():
        sample = row['sample']
        group = row['group']
        path = row['filepath']

        df_data  = read_file_by_chr(path, target_chr)

        if column_names and len(column_names) == df_data.shape[1]:
            df_data.columns = column_names   

        results.append({
            'sample': sample,
            'group': group,
            'data': df_data
        })

    return results




sample_list = load_inlab_chrwise("/home/ly/shell/deepDMR/data/simulate_sample_data/raw/in_test.txt", target_chr="chr1")
print(sample_list[0]['sample'])
print(sample_list[0]['group'])
print(sample_list[0]['data'].head())

