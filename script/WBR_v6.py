import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import beta, f
from scipy.special import logit
from scipy.stats import chi2
from statsmodels.othermod.betareg import BetaModel


# **🚀 Step 1: 处理输入数据**
# 生成模拟数据，包含了原始输入矩阵和统计矩阵
def generate_simulated_data(
        num_cpg=10, interval=50, std_dev=0.05,
        group1_size=5, group2_size=5,
        alpha1=2, beta1=5, alpha2=5, beta2=2,
        group1="g1", group2="g2"):
    """
    生成模拟的 CpG 甲基化数据（宽格式），同时计算统计信息：
    - 每个组的甲基化均值、方差、覆盖度（样本数）
    - 两组之间的甲基化差值 ΔM
    """

    np.random.seed(42)
    positions = np.arange(1, num_cpg * interval + 1, interval)
    
    # 生成实验组（Group1）和对照组（Group2）的甲基化水平
    meth_group1 = np.random.beta(alpha1, beta1, (num_cpg, group1_size))
    meth_group2 = np.random.beta(alpha2, beta2, (num_cpg, group2_size))

    # 计算统计信息
    mean_group1 = np.mean(meth_group1, axis=1)
    var_group1 = np.var(meth_group1, axis=1, ddof=1) if group1_size > 1 else np.zeros(num_cpg)
    mean_group2 = np.mean(meth_group2, axis=1)
    var_group2 = np.var(meth_group2, axis=1, ddof=1) if group2_size > 1 else np.zeros(num_cpg)
    
    delta_M = mean_group2 - mean_group1  # 计算两组间甲基化差值
    block = ["block1"] * num_cpg
    # 构造宽格式数据
    df_wide = pd.DataFrame({
        "Chromosome": ["chr1"] * num_cpg,
        "Position": positions,
        **{f"{group1}_{i+1}": meth_group1[:, i] for i in range(group1_size)},
        **{f"{group2}_{i+1}": meth_group2[:, i] for i in range(group2_size)}
    })

    # 统计数据 DataFrame
    df_summary = pd.DataFrame({
        "Chromosome": ["chr1"] * num_cpg,
        "Position": positions,
        f"{group1}_mean": mean_group1,
        f"{group2}_mean": mean_group2,
        f"{group1}_var": var_group1,
        f"{group2}_var": var_group2,
        f"{group1}_cov": group1_size,
        f"{group2}_cov": group2_size,
        "mean_diff": delta_M,
        "Block": block
    })

    return df_wide, df_summary

# 转换为长型数据
def convert_to_long_format(df_wide):
    """
    将宽格式数据转换为长格式，适用于 WBR 计算
    """
    df_long = df_wide.melt(id_vars=["Chr", "Pos"], var_name="Sample", value_name="Meth_Level")
    df_long["Group"] = df_long["Sample"].apply(lambda x: "group1" if "G1" in x else "group2")
    return df_long

# **🚀 Step 2: 计算 Beta 分布参数**
def compute_beta_params(df_summary, coverage_threshold=5, group1="g1", group2="g2"):
    """
    计算 Beta 分布参数，并处理：
    - 低覆盖度情况<5,使用 Beta(2,2) 先验
    - 仅有 1 个样本的极端情况
    """
    beta_params = []

    for _, row in df_summary.iterrows():
        chr_pos = (row["Chromosome"], row["Position"])
        mean_g1, var_g1, cov_g1 = row[f"{group1}_mean"], row[f"{group1}_var"], row[f"{group1}_cov"]
        mean_g2, var_g2, cov_g2 = row[f"{group2}_mean"], row[f"{group2}_var"], row[f"{group2}_cov"]

        # **1. 低覆盖度（n < 5），使用贝叶斯先验**
        if cov_g1 < coverage_threshold:
            alpha_g1 = 2 + mean_g1 * cov_g1
            beta_g1 = 2 + (cov_g1 - mean_g1 * cov_g1)
        else:
            phi_g1 = mean_g1 * (1 - mean_g1) / (var_g1 + 1e-6) - 1 if var_g1 > 0 else 1
            alpha_g1 = mean_g1 * phi_g1
            beta_g1 = (1 - mean_g1) * phi_g1

        if cov_g2 < coverage_threshold:
            alpha_g2 = 2 + mean_g2 * cov_g2
            beta_g2 = 2 + (cov_g2 - mean_g2 * cov_g2)
        else:
            phi_g2 = mean_g2 * (1 - mean_g2) / (var_g2 + 1e-6) - 1 if var_g2 > 0 else 1
            alpha_g2 = mean_g2 * phi_g2
            beta_g2 = (1 - mean_g2) * phi_g2

        beta_params.append([chr_pos, group1, alpha_g1, beta_g1, cov_g1])
        beta_params.append([chr_pos, group2, alpha_g2, beta_g2, cov_g2])

    df_beta = pd.DataFrame(beta_params, columns=["Chr_Pos", "Group", "Alpha", "Beta", "Coverage"])
    return df_beta

# **🚀 Step 3: 计算权重**
def compute_weights(df_beta, df_summary, lambda_factor=0.5, gamma_factor=0.2, group1="g1", group2="g2"):
    """
    计算权重，考虑组内方差、覆盖度（样本数）和甲基化水平差值。
    
    - `lambda_factor`: 控制甲基化水平差值对权重的影响，默认 0.5
    - `gamma_factor`: 控制覆盖度对权重的影响，默认 0.2
    """
    weights = []

    for _, row in df_summary.iterrows():
        chr_pos = (row["Chromosome"], row["Position"])
        delta_M = abs(row["mean_diff"])  # 直接使用预计算的甲基化水平差值
        cov_g1, cov_g2 = row[f"{group1}_cov"], row[f"{group2}_cov"]  # 组1 & 组2的覆盖度

        # 获取 Beta 参数
        sub_beta = df_beta[df_beta["Chr_Pos"] == chr_pos]
        alpha_g1, beta_g1 = sub_beta[sub_beta["Group"] == group1][["Alpha", "Beta"]].values[0]
        alpha_g2, beta_g2 = sub_beta[sub_beta["Group"] == group2][["Alpha", "Beta"]].values[0]

        # 计算 Beta 分布的方差
        var_beta_g1 = (alpha_g1 * beta_g1) / ((alpha_g1 + beta_g1) ** 2 * (alpha_g1 + beta_g1 + 1))
        var_beta_g2 = (alpha_g2 * beta_g2) / ((alpha_g2 + beta_g2) ** 2 * (alpha_g2 + beta_g2 + 1))

        # 计算覆盖度影响
        coverage_effect_g1 = 1 + gamma_factor * np.log(1 + cov_g1)
        coverage_effect_g2 = 1 + gamma_factor * np.log(1 + cov_g2)

        # 计算最终权重
        weight_g1 = (1 / (var_beta_g1 + 1e-6)) * coverage_effect_g1 * (1 + lambda_factor * delta_M)
        weight_g2 = (1 / (var_beta_g2 + 1e-6)) * coverage_effect_g2 * (1 + lambda_factor * delta_M)

        weights.append([chr_pos, group1, weight_g1])
        weights.append([chr_pos, group2, weight_g2])

    df_weights = pd.DataFrame(weights, columns=["Chr_Pos", "Group", "Weight"])
    return df_weights

# 合并统计数组和权重数组
def prepare_summary_for_merge(df_summary, group1="g1", group2="g2"):
    """
    从 `df_summary` 生成仅包含 `Chr_Pos`, `Group`, `Mean` 的 DataFrame。
    - 先按 `Pos` 排序，再按 `Group` 排序（group1 在前）。
    - 生成 `Chr_Pos` 列，与 `df_weights` 格式匹配。
    """

    # ✅ 创建 `Chr_Pos` 列
    df_summary = df_summary.copy()
    df_summary["Chr_Pos"] = df_summary.apply(lambda row: f"({row['Chromosome']}, {row['Position']})", axis=1)

    # ✅ 仅保留需要的列
    df_summary_g1 = df_summary[["Chr_Pos", f"{group1}_mean"]].copy()
    df_summary_g1.rename(columns={f"{group1}_mean": "mean"}, inplace=True)
    df_summary_g1["Group"] = group1

    df_summary_g2 = df_summary[["Chr_Pos", f"{group2}_mean"]].copy()
    df_summary_g2.rename(columns={f"{group2}_mean": "mean"}, inplace=True)
    df_summary_g2["Group"] = group2

    # ✅ 上下拼接数据
    df_summary_grouped = pd.concat([df_summary_g1, df_summary_g2], axis=0, ignore_index=True)

    # ✅ 按 `Pos` 排序，再按 `Group` 排序（group1 在前）
    df_summary_grouped["Position"] = df_summary_grouped["Chr_Pos"].apply(lambda x: int(x.split(", ")[1][:-1]))  # 提取 Pos
    df_summary_grouped = df_summary_grouped.sort_values(by=["Position", "Group"], ascending=[True, True]).drop(columns=["Position"])

    return df_summary_grouped

# **🚀 Step 4: 运行 WBR**

# MLE + LRT
def mle_beta_regression(df_weights, df_summary, group1="g1", group2="g2", f_value=15):
    """
    使用 MLE 进行 Beta 回归，并使用 LRT 计算 p-value。
    """
    
    df_summary = prepare_summary_for_merge(df_summary, group1=group1, group2=group2)
    
    # ✅ 将 `Mean_Group1` 和 `Mean_Group2` 赋值到 `df_weights`
    #df_weights = df_weights.merge(df_summary[["Chr_Pos", "Group", "Mean"]],
    #                              on=["Chr_Pos", "Group"], how="left")
    df_weights["mean"] = df_summary["mean"].values

    # ✅ 确保 Mean 在 (0,1) 之间，避免 logit 计算错误
    df_weights["mean"] = df_weights["mean"].clip(0.001, 0.999)

    # ✅ 运行 Beta 回归
    df_weights["Group"] = df_weights["Group"].astype("category")  # 确保因子变量
    # 显式指定基准组顺序，防止“写反”问题
    df_weights["Group"] = pd.Categorical(df_weights["Group"], categories=[group1, group2])
    #print(df_weights)
    
    F_stat = compute_f_statistic(df_weights, group1=group1, group2=group2)  # ✅ 使用 df_summary
    if F_stat > f_value:
        model = BetaModel.from_formula("mean ~ Group", df_weights, link=sm.families.links.Logit())
        #print(model)
        result = model.fit()

        #print(result.params)
        beta_0 = result.params["Intercept"]
        beta_1 = result.params[f"Group[T.{group2}]"]
        #beta_1 = result.params["Group[T.group2]"]

        # 计算真实甲基化水平差值 Δμ
        mu_C = np.exp(beta_0) / (1 + np.exp(beta_0))
        mu_T = np.exp(beta_0 + beta_1) / (1 + np.exp(beta_0 + beta_1))
        delta_mu = mu_C - mu_T

        # 计算 LRT p-value
        logL_full = result.llf  # 全模型对数似然
        model_null = BetaModel.from_formula("mean ~ 1", df_weights, link=sm.families.links.Logit())
        result_null = model_null.fit()
        logL_null = result_null.llf

        LRT_stat = -2 * (logL_null - logL_full)
        p_value = chi2.sf(LRT_stat, df=1)  # 计算 p-value
    else:
        p_value = delta_mu = mu_C = mu_T = np.nan

    return p_value, delta_mu, mu_C, mu_T, F_stat

# WLS + Wlad 
def wls_regression(df_weights, df_summary, group1="g1", group2="g2"):
    """
    使用 WLS 进行加权最小二乘回归，并使用 Wald 检验计算 p-value。
    """
    
    df_summary = prepare_summary_for_merge(df_summary, group1=group1, group2=group2)
    
    #print(f"summary is :{df_summary}")
    #print(f"weight is :{df_weights}")

    #df_weights = df_weights.merge(df_summary[["Chr_Pos", "Group", "Mean"]],
    #                              on=["Chr_Pos", "Group"], how="left")

    df_weights["mean"] = df_summary["mean"].values


    # ✅ 确保 Mean 在 (0,1) 之间，避免 logit 计算错误
    df_weights["mean"] = df_weights["mean"].clip(0.001, 0.999)
    
    #print(f"weight is :{df_weights}")
    #print(df_weights)

    # ✅ 重新定义 X, y, weights
    df_weights["logit(mean)"] = np.log(df_weights["mean"] / (1 - df_weights["mean"]))
    X = np.where(df_weights["Group"] == group1, 1, 0)
    X = sm.add_constant(X)  # 添加截距项
    y = df_weights["logit(mean)"]
    weights = df_weights["Weight"]

    model = sm.WLS(y, X, weights=weights).fit()
    beta_0 = model.params.iloc[0]
    beta_1 = model.params.iloc[1]

    # 计算真实甲基化水平差值 Δμ
    mu_C = np.exp(beta_0) / (1 + np.exp(beta_0))
    mu_T = np.exp(beta_0 + beta_1) / (1 + np.exp(beta_0 + beta_1))
    delta_mu = mu_C - mu_T

    # 计算 Wald 检验 p-value
    SE_beta1 = model.bse.iloc[1]  # 标准误
    Z_stat = beta_1 / SE_beta1
    p_value = 2 * (1 - chi2.cdf(abs(Z_stat), df=1))

    return p_value, delta_mu, mu_C, mu_T

# F 检验
def compute_f_statistic(df_weights, group1="g1", group2="g2"):
    """
    计算 F 统计量，用于评估组间和组内方差的比值。
    
    参数:
    - df_weights: 包含 `Group` (group1, group2) 和 `Mean` 列的 DataFrame。
    
    返回:
    - F_stat: F 统计量
    """

    # ✅ 按 `Chr_Pos` 和 `Group` 排序，确保匹配正确
    df_weights = df_weights.sort_values(by=["Chr_Pos", "Group"]).reset_index(drop=True)

    # ✅ 计算组间方差 S_between
    mean_group1 = df_weights[df_weights["Group"] == group1]["mean"].mean()
    mean_group2 = df_weights[df_weights["Group"] == group2]["mean"].mean()
    S_between = np.abs(mean_group1 - mean_group2)  # 避免负数

    # ✅ 计算组内方差 S_within
    var_group1 = df_weights[df_weights["Group"] == group1]["mean"].var(ddof=1)
    var_group2 = df_weights[df_weights["Group"] == group2]["mean"].var(ddof=1)
    S_within = (var_group1 + var_group2) / 2  # 取均值，避免某一组方差为 0

    # ✅ 防止 S_within 过小，避免除零错误
    epsilon = 1e-8  # 设置一个小常数，防止分母过小
    S_within = max(S_within, epsilon)

    # ✅ 计算 F 统计量
    F_stat = S_between / S_within

    return F_stat

def compute_f_tight_statistic(df_weights, group1="g1", group2="g2"):
    """
    计算简化版 F 统计量及其对应的 p-value。

    参数:
    - df_weights: 包含 `Group` 和 `mean` 列的 DataFrame
    - group1, group2: 分组名称

    返回:
    - F_stat: F 统计量
    - p_value: F-stat 的右尾概率
    """

    # 排序
    df_weights = df_weights.sort_values(by=["Chr_Pos", "Group"]).reset_index(drop=True)

    # 获取组内均值
    group1_vals = df_weights[df_weights["Group"] == group1]["mean"]
    group2_vals = df_weights[df_weights["Group"] == group2]["mean"]

    mean1 = group1_vals.mean()
    mean2 = group2_vals.mean()

    # 简化组间“方差” = 差值平方（也可以用平方和定义）
    S_between = (mean1 - mean2)**2

    # 组内方差（标准做法）
    var1 = group1_vals.var(ddof=1)
    var2 = group2_vals.var(ddof=1)
    S_within = (var1 + var2) / 2

    # 避免除以0
    epsilon = 1e-8
    S_within = max(S_within, epsilon)

    # F 值
    F_stat = S_between / S_within

    # 自由度
    df1 = 1
    df2 = len(group1_vals) + len(group2_vals) - 2

    # 计算 p-value（右尾）
    p_value = f.sf(F_stat, df1, df2)

    return F_stat, p_value



def run_weighted_beta_regression(df_summary, threshold=5, group1="g1", group2="g2", f_value=15):
    """
    运行 WBR 加权 Beta 回归 ，选择 WLS 或 MLE Beta 回归。
    自动从 `df_summary` 推断 `group1_size` 和 `group2_size`。
    """
    
    # 计算 Beta 参数
    df_beta = compute_beta_params(df_summary, group1=group1, group2=group2)
    #print(df_beta) 
    # 计算权重
    df_weights = compute_weights(df_beta, df_summary, group1=group1, group2=group2)
    #print(df_weights)

    # 获取 block 内 CpG 位点数
    #num_cpg = len(df_summary)

    # 计算 F 统计量
    p_value, delta_mu, g1_beta, g2_beta, F_stat = mle_beta_regression(df_weights, df_summary, group1=group1, group2=group2, f_value=f_value)  # ✅ 使用 df_weights
    #F_stat = compute_f_statistic(df_weights, group1=group1, group2=group2)  # ✅ 使用 df_summary      
    #print(f"Using MLE (Beta Regression) for block with {num_cpg} CpG sites.")
    is_DMR = (p_value < 0.05) and (F_stat > f_value)

    return {
        "Delta": delta_mu,
        "p-value": p_value,
        "F-statistic": F_stat,
        f"{group1}_mean": g1_beta,
        f"{group2}_mean": g2_beta,
        "DMR": is_DMR
    }

   
if __name__ == "__main__":
    # 生成模拟数据（宽格式 + 统计信息）
    df_wide, df_summary = generate_simulated_data()
    print(df_wide.head(10))
    print(df_summary)
    # 运行 WBR 并输出结果
    results = run_weighted_beta_regression(df_summary)
    
    # 打印结果
    print("Results:")
    for key, value in results.items():
        print(f"{key}: {value}")


