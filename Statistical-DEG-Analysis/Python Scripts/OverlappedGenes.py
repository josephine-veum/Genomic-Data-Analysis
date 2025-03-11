# python script to identify overlapped genes from the significant genes lists
# reads in data from SeqData.txt/ArrayData.txt, runs hypothesis tests, and compares results

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

# function to run hypothesis tests on all genes
# takes in filename (ArrayData.txt or SeqData.txt), filenames indicating where to save histograms, & histogram titles
# returns a dictionary containing the 1000 most significant genes from the t-tests and Wilcoxon rank-sum tests
def hypothesis_test(filename):
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

    # determine and return the 1000 most significant genes & their p-values from each test
    ttest_1000 = most_significant(ttest_pvals, num_vals=1000)
    wilcoxon_1000 = most_significant(wilcoxon_pvals, num_vals=1000)

    return ttest_1000, wilcoxon_1000

# function to plot selected genes versus the number of common genes between RNAseq & microarray significant genes
# takes in dictionary of most significant 1000 genes from each dataset, filename indicating where to save the plot, and title for the plot
def plot_overlaps(rna_1000, array_1000, filename, title):
    num_genes = list(range(10, 1001, 10))
    common_genes = []

    for i in num_genes:
        # take the first i most significant genes from each dataset
        rna_significant = dict(sorted(rna_1000.items(), key=lambda item: item[1])[:i])
        array_significant = dict(sorted(array_1000.items(), key=lambda item: item[1])[:i])

        # determine how many genes overlap & add this value to the common genes list
        overlap = rna_significant.keys() & array_significant.keys()
        common_genes.append(len(overlap))

    plt.figure(figsize=(6,4))
    plt.scatter(num_genes, common_genes)
    plt.xlabel('selected genes')
    plt.ylabel('common genes')
    plt.title(title, fontweight='bold')
    plt.savefig(filename)

# run hypothesis tests for ArrayData.txt and SeqData.txt
if __name__ == "__main__":
    # get the 1000 most significant genes from each dataset using each test
    rna_ttest_1000, rna_wilcoxon_1000 = hypothesis_test("SeqData.txt")
    array_ttest_1000, array_wilcoxon_1000 = hypothesis_test("ArrayData.txt")

    # create plots of the number of selected genes vs. the number of common genes for each test
    plot_overlaps(rna_ttest_1000, array_ttest_1000, "ttest_scatter.png", "Overlapped genes for t-test")
    plot_overlaps(rna_wilcoxon_1000, array_wilcoxon_1000, "wilcoxon_scatter.png", "Overlapped genes for Wilcoxon rank-sum test")