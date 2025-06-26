#/src/core/hypo_test.py
# Created on 2025-04-06

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ks_2samp

import seaborn as sns
import matplotlib.pyplot as plt
from src.config import path

def load_energy_file(filepath):

    data = pd.read_csv(filepath, sep='\t', header=None, names=['ID', 'Energy'])
    print(f"Before filtering ... {data.shape}")
    # Remove rows with negative energy values
    data = data[data['Energy'] >= 0]

    print(f"After filtering ... {data.shape}")
    return data

def calculate_z_scores(bound_df, unbound_df):
    # Concatenate both datasets to get combined mean and std
    combined_df = pd.concat([bound_df.assign(Type='Bound'), unbound_df.assign(Type='Unbound')])

    mean_energy = combined_df['Energy'].mean()
    std_energy = combined_df['Energy'].std()

    # Clearly explained Z-score calculation
    combined_df['Z_score'] = (combined_df['Energy'] - mean_energy) / std_energy
    return combined_df

# Statistical testing (Mann-Whitney U-test) clearly shown
def perform_statistical_test(combined_df):
    bound_z = combined_df[combined_df['Type']=='Bound']['Z_score']
    unbound_z = combined_df[combined_df['Type']=='Unbound']['Z_score']

    # Mann-Whitney U test (testing if Bound < Unbound)
    stat, p_value = mannwhitneyu(bound_z, unbound_z, alternative='less')

    print(f"Mann-Whitney U test statistic: {stat:.4f}")
    print(f"P-value (Bound < Unbound): {p_value:.3e}")

    return stat, p_value

# Visualization (clearly intuitive)
def visualize(combined_df):
    sns.violinplot(x='Type', y='Z_score', data=combined_df, palette=['coral', 'skyblue'])
    plt.axhline(0, color='gray', linestyle='--')
    plt.title('Z-score distribution (Bound vs Unbound)')
    plt.ylabel('Z-score (Normalized Free Energy)')
    plt.xlabel('')
    plt.tight_layout()
    plt.show()

# Main analysis function
def main(bound_file, unbound_file):
    print("Loading data...")
    bound_df = load_energy_file(bound_file)
    unbound_df = load_energy_file(unbound_file)

     # Find common IDs between the two dataframes
    common_indices = bound_df.index.intersection(unbound_df.index)

    # Filter both dataframes to keep only rows with common IDs
    bound_df = bound_df.loc[common_indices]
    unbound_df = unbound_df.loc[common_indices]

    print("Calculating Z-scores...")
    combined_df = calculate_z_scores(bound_df, unbound_df)

    print("\nPerforming statistical significance test...")
    perform_statistical_test(combined_df)

    print("\nVisualizing results...")
    visualize(combined_df)

def permutation_test(group1, group2, num_permutations=10000, alternative='less'):
    # Observed difference in median (or mean) energies
    observed_diff = np.median(group1) - np.median(group2)

    combined = np.concatenate([group1, group2])
    count = 0

    perm_diffs = []

    for _ in range(num_permutations):
        np.random.shuffle(combined)
        perm_group1 = combined[:len(group1)]
        perm_group2 = combined[len(group1):]

        perm_diff = np.median(perm_group1) - np.median(perm_group2)
        perm_diffs.append(perm_diff)

        if alternative == 'less' and perm_diff <= observed_diff:
            count += 1
        elif alternative == 'greater' and perm_diff >= observed_diff:
            count += 1
        elif alternative == 'two-sided' and abs(perm_diff) >= abs(observed_diff):
            count += 1

    p_value = count / num_permutations

    # Plot permutation distribution
    plt.hist(perm_diffs, bins=50, color='skyblue', alpha=0.7, label='Permutation diffs')
    plt.axvline(observed_diff, color='red', linestyle='--', label='Observed diff')
    plt.xlabel('Difference in median energies')
    plt.ylabel('Frequency')
    plt.title('Permutation Test')
    plt.legend()
    plt.show()

    return observed_diff, p_value



def ks_test(group1, group2):
    # Perform KS test
    statistic, p_value = ks_2samp(group1, group2, alternative='less')

    # Plot CDF clearly
    plt.figure(figsize=(8,5))
    for data, label, color in zip([group1, group2], ['Bound', 'Unbound'], ['coral', 'skyblue']):
        sorted_data = np.sort(data)
        yvals = np.arange(len(sorted_data))/float(len(sorted_data))
        plt.step(sorted_data, yvals, label=label, color=color)

    plt.xlabel('Free energy values')
    plt.ylabel('Cumulative fraction')
    plt.title('KS Test - Cumulative distributions')
    plt.legend()
    plt.show()

    return statistic, p_value




if __name__ == '__main__':

    # bound_file = path.RESULTS_DIR/ 'bound_free_energy_results_all_MD.txt'
    # unbound_file =  path.RESULTS_DIR/ 'unbound_free_energy_results_all_MD.txt'    

    # bound_file = path.RESULTS_DIR/ 'bound_free_energy_results_all_mixed.txt'
    # unbound_file =  path.RESULTS_DIR/ 'unbound_free_energy_results_all_mixed.txt'

    bound_file = path.RESULTS_DIR/ 'bound_free_energy_results_Olson.txt'
    unbound_file =  path.RESULTS_DIR/ 'unbound_free_energy_results_Olson.txt'


    main(bound_file, unbound_file)





    bound_df = load_energy_file(bound_file)
    unbound_df = load_energy_file(unbound_file)
    print(bound_df.shape)
    print(unbound_df.shape)

    # Find common IDs between the two dataframes
    common_indices = bound_df.index.intersection(unbound_df.index)

       # Filter both dataframes to keep only rows with common IDs
    bound_df = bound_df.loc[common_indices]
    unbound_df = unbound_df.loc[common_indices]


    # observed_diff, p_value = permutation_test(bound_df['Energy'], unbound_df['Energy'], num_permutations=10000, alternative='less')
    # print(f"Observed difference in median energies: {observed_diff:.4f}")
    # print(f"P-value from permutation test: {p_value:.3e}")


    # # Perform KS test
    # statistic, p_value = ks_test(bound_df['Energy'], unbound_df['Energy'])
    # print(f"KS test statistic: {statistic:.4f}")
    # print(f"P-value from KS test: {p_value:.3e}")

    
    # stat, p_value = mannwhitneyu(bound_df['Energy'], unbound_df['Energy'], alternative='less')
    # print(f"Mann-Whitney U test p-value (raw energies): {p_value:.3e}")


    # Plot the distributions of bound and unbound energies
    plt.figure(figsize=(8, 5))
    sns.kdeplot(bound_df['Energy'], label='Bound', color='coral', fill=True, alpha=0.7)
    sns.kdeplot(unbound_df['Energy'], label='Unbound', color='skyblue', fill=True,alpha=0.7)
    plt.title('Retained Nucleosomes vs Free DNA in Sperms - Energy Distribution (Olson)')
    plt.xlabel('Energy')
    plt.ylabel('Density')
    plt.legend()
    plt.tight_layout()
    plt.savefig(path.RESULTS_DIR / 'figs/olson_param_energy_dist.pdf', dpi=300)
    # plt.show()