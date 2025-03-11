# python script to identify differentially expressed genes
# reads in data from SeqData.txt/ArrayData.txt
# uses Bonferroni correction with t-test (with FDR) with a significance level 0.05

import pandas as pd
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

# function to apply bonferroni correction
# takes in a dictionary with genes & their p_values and the significance level (alpha)
# returns a dictionary with p-values updated with bonferroni correction
def bonferroni(p_vals, alpha=0.05):
    genes = list(p_vals.keys())
    values = list(p_vals.values())

    _, corrected_vals, _, _ = multipletests(values, alpha=alpha, method='bonferroni')
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

# function to apply Bonferroni correction to t-tests from given file
# takes in the filename (SeqData.txt or ArrayData.txt)
# returns nothing but prints out the number of significant genes and theses genes & their p-values
def apply_bonferroni(filename):
    # read in data from SeqData.txt or ArrayData.txt, store in dataframe
    df = pd.read_csv(filename, sep='\t')

    # remove genes with no expression from the dataframe
    df = filter_genes(df)

    # split data into group 1 and group 2
    group1, group2 = split_data(df)

    # apply t-test on the two groups
    # this function returns dictionaries with genes and their p-values
    ttest_pvals = t_test(group1, group2)

    # apply Bonferroni correction on p-values from t-tests
    bonferroni_pvals = bonferroni(ttest_pvals)

    # print out number & list of statistically significant differentially expressed genes with alpha=0.05
    bonferroni_significant = significant(bonferroni_pvals)
    print("Significant genes after Bonferroni correction: ", len(bonferroni_significant))
    print_pvals(bonferroni_significant)

# use Bonferroni correction for t-tests on data from SeqData.txt and ArrayData.txt
if __name__ == "__main__":
    print("Bonferroni correction for t-tests on SeqData.txt")
    apply_bonferroni("SeqData.txt")

    print("\nBonferroni correction for t-tests on ArrayData.txt")
    apply_bonferroni("ArrayData.txt")