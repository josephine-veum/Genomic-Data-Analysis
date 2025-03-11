# python script to identify differentially expressed genes
# reads in data from SeqData.txt/ArrayData.txt
# uses t-test and Wilcoxon rank-sum statistics with a significance level of 0.05

import pandas as pd
from scipy import stats
from matplotlib import pyplot as plt

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

# function to run Wilcoxon rank-sum test on all genes
# takes in two dataframes (group 1 & group 2)
# returns a dictionary containing the genes and their respective p-values
def wilcoxon_test(group1, group2):
    p_vals = {}
    for gene in group1.columns[2:]:
        group1_data = group1[gene].to_list()
        group2_data = group2[gene].to_list()

        stat, p_val = stats.mannwhitneyu(group1_data, group2_data) # mannwhitneyu used for independent samples
        p_vals[gene] = float(p_val)
    return p_vals

# function to get the 10 genes with the most significant difference in expression between the two groups
# takes in the dictionary with gene names and their respective p-values
# returns a dictionary containing the genes with the most significant differences in expression
def most_significant(p_vals, num_vals=10):
    small_p_vals = dict(sorted(p_vals.items(), key=lambda item: item[1])[:num_vals])
    return small_p_vals

# function to print genes and their respective p-values
# takes in a dictionary with gene names and p-values
def print_pvals(p_vals):
    for gene, p_val in p_vals.items():
        print(gene, ": ", p_val)

# function to count the number of statistically significant differences in gene expression between groups
# takes in a dictionary of genes and their respective p-values and the significance level (alpha)
# returns dictionary with genes statistically significant differentally expressed genes and their p-values
def significant(p_vals, alpha=0.05):
    significant_vals = {}
    for gene, p_val in p_vals.items():
        if p_val < 0.05:
            significant_vals[gene] = p_val
    return significant_vals

# function to create and save a histogram of p-values and their frequencies
# takes in dictionary with genes & p-values, title for the plot, and filename to save the plot to
def histogram(p_vals, title, filename):
    plt.figure(figsize=(6,4))
    plt.hist(p_vals.values(), bins=50)
    plt.xlabel('p-value')
    plt.ylabel('frequency')
    plt.title(title, fontweight='bold')
    plt.savefig(filename)

# function to run hypothesis tests on all genes
# takes in filename (ArrayData.txt or SeqData.txt), filenames indicating where to save histograms, & histogram titles
# returns nothing but prints out significant genes and their p-values and makes histograms
def hypothesis_test(filename, ttest_save, wilcoxon_save, ttest_title, wilcoxon_title):
    # read in data from SeqData.txt or ArrayData.txt, store in dataframe
    df = pd.read_csv(filename, sep='\t')

    # remove genes with no expression from the dataframe
    df = filter_genes(df)

    # split data into group 1 and group 2
    group1, group2 = split_data(df)

    # apply t-test and Wilcoxon rank-sum test on the two groups
    # these functions return dictionaries with genes and their p-values
    ttest_pvals = t_test(group1, group2)
    wilcoxon_pvals = wilcoxon_test(group1, group2)

    # get the 10 most significant genes & their p-values from each test
    ttest_10 = most_significant(ttest_pvals)
    wilcoxon_10 = most_significant(wilcoxon_pvals)

    # list out the most significant genes & their p-values
    print("Most significant differently expressed genes between groups 1 & 2")
    print("Using t-test")
    print_pvals(ttest_10)
    print("\nUsing Wilcoxon rank-sum test")
    print_pvals(wilcoxon_10)
    
    # calculate number of statistically significant differences in gene expression for each test
    ttest_significant = significant(ttest_pvals)
    wilcoxon_significant = significant(wilcoxon_pvals)

    # print out number of statistically significant differentially expressed genes
    # this number is the same as the number of null hypotheses that were rejected
    print("\nNumber of statistically significant differentially expressed genes (using alpha=0.05)")
    print("Using t-test: ", len(ttest_significant))
    print("Using Wilcoxon rank-sum test: ", len(wilcoxon_significant))

    # for each test, create and save histogram as an image file
    histogram(ttest_pvals, ttest_title, ttest_save)
    histogram(wilcoxon_pvals, wilcoxon_title, wilcoxon_save)

# run hypothesis tests for ArrayData.txt and SeqData.txt
if __name__ == "__main__":
    print("Hypothesis tests for SeqData.txt")
    hypothesis_test("SeqData.txt", "ttest_rna.png", "wilcoxon_rna.png", "p-values for t-test on RNAseq data", "p-values for Wilcoxon rank-sum test on RNAseq data")
    
    print("\nHypothesis tests for ArrayData.txt")
    hypothesis_test("ArrayData.txt", "ttest_array.png", "wilcoxon_array.png", "p-values for t-test on microarray data", "p-values for Wilcoxon rank-sum test on microarray data")