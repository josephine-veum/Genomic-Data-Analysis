# python script to identify differentially expressed genes
# reads in data from SeqData.txt/ArrayData.txt
# uses Bonferroni correction with t-test (with FDR) with a significance level 0.05

import pandas as pd
import random
from scipy import stats
from statsmodels.stats.multitest import multipletests

# function to remove genes with no expression in either group
# takes in a dataframe, returns a dataframe genes removed
def filter_genes(group):
    for gene in group.columns[2:]:
        if (group[gene].sum() == 0):
            group = group.drop(gene, axis=1)
    return group

# function to split data into two groups
# takes in full dataframe, returns two data frames (group 1 & 2)
def split_data(data):
    group1 = data[data['Group ID'] == 1].copy()
    group2 = data[data['Group ID'] == 2].copy()

    return group1, group2

# function to run t-tests on all genes
# takes in two dataframes (group 1 & group 2)
# returns a dictionary containing the genes and their respective p-values
def t_test(group1, group2):
    p_vals = {}
    for gene in group1.columns[2:]:
        group1_data = group1[gene].to_list()
        group2_data = group2[gene].to_list()

        t_stat, p_val = stats.ttest_ind(group1_data, group2_data) #ttest_ind used for independent samples
        p_vals[gene] = float(p_val)
    return p_vals

# function to apply FDR to p-values
# takes in a dictionary with genes & their p_values and the significance level (alpha)
# returns a dictionary with p-values updated using FDR
def fdr(p_vals, alpha=0.05):
    genes = list(p_vals.keys())
    values = list(p_vals.values())

    _, corrected_vals, _, _ = multipletests(values, alpha=alpha, method='fdr_bh')
    return dict(zip(genes, corrected_vals))

# function to count the number of statistically significant differences in gene expression between groups
# takes in a dictionary of genes and their respective p-values and the significance level (alpha)
# returns dictionary with genes statistically significant differentally expressed genes and their p-values
def significant(p_vals, alpha=0.05):
    significant_vals = {}
    for gene, p_val in p_vals.items():
        if p_val < alpha:
            significant_vals[gene] = p_val
    return significant_vals

# function to print genes and their respective p-values
# takes in a dictionary with gene names and p-values
def print_pvals(p_vals):
    for gene, p_val in p_vals.items():
        print(gene, ": ", p_val)

# function to take n random items from a dictionary
# takes in a dictionary, value n, and a seed for reproductibility
# returns a dicitonary containing the n random items of the inputed dictionary
def dict_slice(my_dict, n, seed=3737):
    random.seed(seed)
    random_keys = random.sample(list(my_dict.keys()), n)
    return {key: my_dict[key] for key in random_keys}
    
# function to apply FDR correction to t-tests from given file
# takes in the filename (SeqData.txt or ArrayData.txt) and number of genes to select
# returns nothing but prints out the number of significant genes and theses genes & their p-values
def apply_FDR(filename, num_genes, alpha=0.05):
    # read in data from SeqData.txt or ArrayData.txt, store in dataframe
    df = pd.read_csv(filename, sep='\t')

    # remove genes with no expression from the dataframe
    df = filter_genes(df)

    # split data into group 1 and group 2
    group1, group2 = split_data(df)

    # apply t-test on the two groups
    # this function returns dictionaries with genes and their p-values
    ttest_pvals = t_test(group1, group2)

    # apply FDR correction on p-values from t-tests using given number of genes
    fdr_pvals = fdr(dict_slice(ttest_pvals, num_genes), alpha=alpha)

    # find significant genes after FDR correction
    fdr_significant = significant(fdr_pvals, alpha=alpha)

    # compute upper bound of FDR correction
    num_significant = len(fdr_significant)
    estimated_nulls = num_genes - num_significant
    upper_bound = (estimated_nulls / num_genes) * alpha

    # print out number of significant genes, what the genes are, and the upper bound of FDR
    print("Significant genes after FDR correction: ", num_significant)
    print_pvals(fdr_significant)
    print("Upper bound of FDR correction: ", upper_bound)

# use FDR correction for t-tests on data from SeqData.txt and ArrayData.txt
# use 20, 50, 100, and 200 for the number of selected genes
if __name__ == "__main__":
    num_genes = [20, 50, 100, 200]
    
    # FDR correction for t-tests from SeqData.txt
    for num in num_genes:
        print("FDR correction for t-tests on SeqData.txt using ", num, " genes")
        apply_FDR("SeqData.txt", num)
        print()

    # FDR correction for t-tests from ArrayData.txt
    for num in num_genes:
        print("FDR correction for t-tests on ArrayData.txt using ", num, " genes")
        apply_FDR("ArrayData.txt", num)
        print()