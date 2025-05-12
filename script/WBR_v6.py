import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import beta, f
from scipy.special import logit
from scipy.stats import chi2
from statsmodels.othermod.betareg import BetaModel


# **🚀 Step 1: Process input data**
# Generate simulated data, which includes the original input matrix and a statistical matrix
def generate_simulated_data(
        num_cpg=10, interval=50, std_dev=0.05,
        group1_size=5, group2_size=5,
        alpha1=2, beta1=5, alpha2=5, beta2=2,
        group1="g1", group2="g2"):
    """
    Generate simulated CpG methylation data (in wide format), while also calculating the following statistical information:
    - For each group: mean methylation level, variance, and coverage (number of samples)
    - Methylation difference (ΔM) between the two groups
    """

    np.random.seed(42)
    positions = np.arange(1, num_cpg * interval + 1, interval)
    
    # Generate methylation levels for the experimental group (Group1) and the control group (Group2)
    meth_group1 = np.random.beta(alpha1, beta1, (num_cpg, group1_size))
    meth_group2 = np.random.beta(alpha2, beta2, (num_cpg, group2_size))

    # Calculate statistical information
    mean_group1 = np.mean(meth_group1, axis=1)
    var_group1 = np.var(meth_group1, axis=1, ddof=1) if group1_size > 1 else np.zeros(num_cpg)
    mean_group2 = np.mean(meth_group2, axis=1)
    var_group2 = np.var(meth_group2, axis=1, ddof=1) if group2_size > 1 else np.zeros(num_cpg)
    
    delta_M = mean_group2 - mean_group1  # Calculate the methylation difference between the two groups
    block = ["block1"] * num_cpg
    # Construct wide-format data
    df_wide = pd.DataFrame({
        "Chromosome": ["chr1"] * num_cpg,
        "Position": positions,
        **{f"{group1}_{i+1}": meth_group1[:, i] for i in range(group1_size)},
        **{f"{group2}_{i+1}": meth_group2[:, i] for i in range(group2_size)}
    })

    # Statistical data DataFrame
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

# Convert to long-format data
def convert_to_long_format(df_wide):
    """
    Convert wide-format data to long-format, suitable for WBR calculation
    """
    df_long = df_wide.melt(id_vars=["Chr", "Pos"], var_name="Sample", value_name="Meth_Level")
    df_long["Group"] = df_long["Sample"].apply(lambda x: "group1" if "G1" in x else "group2")
    return df_long

# **🚀 Step 2: Calculate Beta distribution parameters**
def compute_beta_params(df_summary, coverage_threshold=5, group1="g1", group2="g2"):
    """
    Calculate Beta distribution parameters, with specific handling for the following cases:
    - Low coverage situations (coverage < 5): A Beta(2,2) prior is used
    - Edge cases (or extreme cases) where there is only 1 sample
    """
    beta_params = []

    for _, row in df_summary.iterrows():
        chr_pos = (row["Chromosome"], row["Position"])
        mean_g1, var_g1, cov_g1 = row[f"{group1}_mean"], row[f"{group1}_var"], row[f"{group1}_cov"]
        mean_g2, var_g2, cov_g2 = row[f"{group2}_mean"], row[f"{group2}_var"], row[f"{group2}_cov"]

        # **1. For low coverage (n < 5), use a Bayesian prior**
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

# **🚀 Step 3: Calculate weights**
def compute_weights(df_beta, df_summary, lambda_factor=0.5, gamma_factor=0.2, group1="g1", group2="g2"):
    """
    Calculate weights, considering within-group variance, coverage (number of samples), and methylation level difference
    
    - `lambda_factor`: Controls the influence of the methylation level difference on the weights (default is 0.5)
    - `gamma_factor`: Controls the influence of coverage (number of samples) on the weights (default is 0.2)
    """
    weights = []

    for _, row in df_summary.iterrows():
        chr_pos = (row["Chromosome"], row["Position"])
        delta_M = abs(row["mean_diff"])  # Directly use the pre-calculated methylation level difference
        cov_g1, cov_g2 = row[f"{group1}_cov"], row[f"{group2}_cov"]  # Coverage of Group 1 & Group 2

        # Get Beta parameters
        sub_beta = df_beta[df_beta["Chr_Pos"] == chr_pos]
        alpha_g1, beta_g1 = sub_beta[sub_beta["Group"] == group1][["Alpha", "Beta"]].values[0]
        alpha_g2, beta_g2 = sub_beta[sub_beta["Group"] == group2][["Alpha", "Beta"]].values[0]

        # Calculate the variance of the Beta distribution
        var_beta_g1 = (alpha_g1 * beta_g1) / ((alpha_g1 + beta_g1) ** 2 * (alpha_g1 + beta_g1 + 1))
        var_beta_g2 = (alpha_g2 * beta_g2) / ((alpha_g2 + beta_g2) ** 2 * (alpha_g2 + beta_g2 + 1))

        # Calculate the coverage effect
        coverage_effect_g1 = 1 + gamma_factor * np.log(1 + cov_g1)
        coverage_effect_g2 = 1 + gamma_factor * np.log(1 + cov_g2)

        # Calculate final weights
        weight_g1 = (1 / (var_beta_g1 + 1e-6)) * coverage_effect_g1 * (1 + lambda_factor * delta_M)
        weight_g2 = (1 / (var_beta_g2 + 1e-6)) * coverage_effect_g2 * (1 + lambda_factor * delta_M)

        weights.append([chr_pos, group1, weight_g1])
        weights.append([chr_pos, group2, weight_g2])

    df_weights = pd.DataFrame(weights, columns=["Chr_Pos", "Group", "Weight"])
    return df_weights

# Merge the statistics array and the weights array
def prepare_summary_for_merge(df_summary, group1="g1", group2="g2"):
    """
    Generate a DataFrame from df_summary that includes only the 'Chr_Pos', 'Group', and 'Mean' columns.
    - First, sort the data by the 'Pos' column, then by the 'Group' column (ensuring that 'group1' comes before other groups in the sort order where applicable)
    - Generate or format the 'Chr_Pos' column to match the format of the 'Chr_Pos' column in the `df_weights` DataFrame
    """

    # ✅ Create a Chr_Pos column
    df_summary = df_summary.copy()
    df_summary["Chr_Pos"] = df_summary.apply(lambda row: f"({row['Chromosome']}, {row['Position']})", axis=1)

    # ✅ Keep only the required columns
    df_summary_g1 = df_summary[["Chr_Pos", f"{group1}_mean"]].copy()
    df_summary_g1.rename(columns={f"{group1}_mean": "mean"}, inplace=True)
    df_summary_g1["Group"] = group1

    df_summary_g2 = df_summary[["Chr_Pos", f"{group2}_mean"]].copy()
    df_summary_g2.rename(columns={f"{group2}_mean": "mean"}, inplace=True)
    df_summary_g2["Group"] = group2

    # ✅ Concatenate data vertically 
    df_summary_grouped = pd.concat([df_summary_g1, df_summary_g2], axis=0, ignore_index=True)

    # ✅ Sort by Pos, then by Group (with group1 first)
    df_summary_grouped["Position"] = df_summary_grouped["Chr_Pos"].apply(lambda x: int(x.split(", ")[1][:-1]))  # 提取 Pos
    df_summary_grouped = df_summary_grouped.sort_values(by=["Position", "Group"], ascending=[True, True]).drop(columns=["Position"])

    return df_summary_grouped

# **🚀 Step 4: Run WBR**

# MLE + LRT
def mle_beta_regression(df_weights, df_summary, group1="g1", group2="g2", f_value=15):
    """
    Use MLE for Beta regression and calculate the p-value using LRT
    """
    
    df_summary = prepare_summary_for_merge(df_summary, group1=group1, group2=group2)
    
    # ✅ Assign Mean_Group1 and Mean_Group2 to df_weights
    df_weights["mean"] = df_summary["mean"].values

    # ✅ Ensure Mean values are strictly between 0 and 1 (i.e., in the open interval (0,1)) to avoid logit calculation errors
    df_weights["mean"] = df_weights["mean"].clip(0.001, 0.999)

    # ✅ Run Beta regression
    df_weights["Group"] = df_weights["Group"].astype("category")  # Ensure variables are treated as factors (factor variables)
    # Explicitly specify the reference group order to prevent issues of 'reversed interpretation' (e.g., misinterpreting the direction of effects)
    df_weights["Group"] = pd.Categorical(df_weights["Group"], categories=[group1, group2])

    
    F_stat = compute_f_statistic(df_weights, group1=group1, group2=group2)  # ✅ Use df_summary
    if F_stat > f_value:
        model = BetaModel.from_formula("mean ~ Group", df_weights, link=sm.families.links.Logit())
        result = model.fit()

        
        beta_0 = result.params["Intercept"]
        beta_1 = result.params[f"Group[T.{group2}]"]
        

        # Calculate the true methylation level difference Δμ
        mu_C = np.exp(beta_0) / (1 + np.exp(beta_0))
        mu_T = np.exp(beta_0 + beta_1) / (1 + np.exp(beta_0 + beta_1))
        delta_mu = mu_C - mu_T

        # Calculate LRT p-value
        logL_full = result.llf  # Full model log-likelihood
        model_null = BetaModel.from_formula("mean ~ 1", df_weights, link=sm.families.links.Logit())
        result_null = model_null.fit()
        logL_null = result_null.llf

        LRT_stat = -2 * (logL_null - logL_full)
        p_value = chi2.sf(LRT_stat, df=1)  # Calculate p-value
    else:
        p_value = delta_mu = mu_C = mu_T = np.nan

    return p_value, delta_mu, mu_C, mu_T, F_stat

# WLS + Wlad 
def wls_regression(df_weights, df_summary, group1="g1", group2="g2"):
    """
    Perform weighted least squares (WLS) regression and use the Wald test to calculate the p-value
    """
    
    df_summary = prepare_summary_for_merge(df_summary, group1=group1, group2=group2)
    

    df_weights["mean"] = df_summary["mean"].values


    # ✅ Ensure Mean values are strictly between 0 and 1 (i.e., in the open interval (0,1)) to avoid logit calculation errors
    df_weights["mean"] = df_weights["mean"].clip(0.001, 0.999)
    

    # ✅ Redefine X, y, weights
    df_weights["logit(mean)"] = np.log(df_weights["mean"] / (1 - df_weights["mean"]))
    X = np.where(df_weights["Group"] == group1, 1, 0)
    X = sm.add_constant(X)  # Add an intercept term
    y = df_weights["logit(mean)"]
    weights = df_weights["Weight"]

    model = sm.WLS(y, X, weights=weights).fit()
    beta_0 = model.params.iloc[0]
    beta_1 = model.params.iloc[1]

    # Calculate the true methylation level difference Δμ
    mu_C = np.exp(beta_0) / (1 + np.exp(beta_0))
    mu_T = np.exp(beta_0 + beta_1) / (1 + np.exp(beta_0 + beta_1))
    delta_mu = mu_C - mu_T

    # Calculate the Wald test p-value
    SE_beta1 = model.bse.iloc[1]  # Standard Error
    Z_stat = beta_1 / SE_beta1
    p_value = 2 * (1 - chi2.cdf(abs(Z_stat), df=1))

    return p_value, delta_mu, mu_C, mu_T

# F-test
def compute_f_statistic(df_weights, group1="g1", group2="g2"):
    """
    Calculates the F-statistic, which is used to assess the ratio of between-group variance to within-group variance
    
    Parameters:
    - df_weights: A DataFrame containing 'Group' (with levels such as 'group1', 'group2') and 'Mean' columns.
    
    Returns:
    - F_stat: The calculated F-statistic
    """

    # ✅ Sort by Chr_Pos and Group to ensure correct matching
    df_weights = df_weights.sort_values(by=["Chr_Pos", "Group"]).reset_index(drop=True)

    # ✅ Calculate between-group variance S_between
    mean_group1 = df_weights[df_weights["Group"] == group1]["mean"].mean()
    mean_group2 = df_weights[df_weights["Group"] == group2]["mean"].mean()
    S_between = np.abs(mean_group1 - mean_group2)  # Avoid negative numbers

    # ✅ Calculate within-group variance S_within
    var_group1 = df_weights[df_weights["Group"] == group1]["mean"].var(ddof=1)
    var_group2 = df_weights[df_weights["Group"] == group2]["mean"].var(ddof=1)
    S_within = (var_group1 + var_group2) / 2  # Take the mean (e.g., of variances), to avoid a scenario where any single group's variance is 0

    # ✅ Prevent S_within (e.g., within-group sum of squares or variance) from being too small to avoid division by zero errors
    epsilon = 1e-8  # Set a small constant (epsilon) to prevent the denominator from being too small
    S_within = max(S_within, epsilon)

    # ✅ Calculate the F-statistic
    F_stat = S_between / S_within

    return F_stat

def compute_f_tight_statistic(df_weights, group1="g1", group2="g2"):
    """
    Calculates a simplified F-statistic and its corresponding p-value

    Parameters:
    - df_weights: A DataFrame containing 'Group' and 'mean' columns
    - group1, group2: The name of the first group

    Returns:
    - F_stat: The calculated simplified F-statistic
    - p_value: The right-tailed probability (p-value) for the F-statistic
    """

    # sort
    df_weights = df_weights.sort_values(by=["Chr_Pos", "Group"]).reset_index(drop=True)

    # Get within-group means
    group1_vals = df_weights[df_weights["Group"] == group1]["mean"]
    group2_vals = df_weights[df_weights["Group"] == group2]["mean"]

    mean1 = group1_vals.mean()
    mean2 = group2_vals.mean()

    # Simplified between-group 'variance' = squared difference (can also be defined using sum of squares)
    S_between = (mean1 - mean2)**2

    # Within-group variance (standard method)
    var1 = group1_vals.var(ddof=1)
    var2 = group2_vals.var(ddof=1)
    S_within = (var1 + var2) / 2

    # Avoid division by zero
    epsilon = 1e-8
    S_within = max(S_within, epsilon)

    # F data
    F_stat = S_between / S_within

    # Degrees of freedom
    df1 = 1
    df2 = len(group1_vals) + len(group2_vals) - 2

    # Calculate p-value(right-tailed)
    p_value = f.sf(F_stat, df1, df2)

    return F_stat, p_value



def run_weighted_beta_regression(df_summary, threshold=5, group1="g1", group2="g2", f_value=15):
    """
    Run WBR weighted Beta regression, with the option to choose either WLS (Weighted Least Squares) or MLE (Maximum Likelihood Estimation) Beta regression
    Automatically infer group1_size and group2_size from df_summary
    """
    
    # Calculate Beta parameters
    df_beta = compute_beta_params(df_summary, group1=group1, group2=group2)
    #print(df_beta) 
    # Calculate weights
    df_weights = compute_weights(df_beta, df_summary, group1=group1, group2=group2)


    # Calculate the F-statistic
    p_value, delta_mu, g1_beta, g2_beta, F_stat = mle_beta_regression(df_weights, df_summary, group1=group1, group2=group2, f_value=f_value)  # ✅ 使用 df_weights
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
    # Generate simulated data (wide format + statistical information)
    df_wide, df_summary = generate_simulated_data()
    print(df_wide.head(10))
    print(df_summary)
    # Run WBR and output the results
    results = run_weighted_beta_regression(df_summary)
    
    # Print result
    print("Results:")
    for key, value in results.items():
        print(f"{key}: {value}")


