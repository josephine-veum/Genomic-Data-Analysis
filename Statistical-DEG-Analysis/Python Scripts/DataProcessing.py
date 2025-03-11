# python script to preprocess data
# reads the RNAseq & microarray data
# dumps preprocessed data into SeqData.txt and ArrayData.txt

import pandas as pd

# read in data from files to dataframes
RNAseq_df = pd.read_csv('HiSeqV2', sep='\t', header=0)
microarray_df = pd.read_csv('HT_HG-U133A', sep='\t', header=0)
clinical_df = pd.read_csv('ov_tcga_clinical_data.tsv', sep='\t', header=0)

# lists to hold sample IDs for group 1 and group 2
group1 = []
group2 = []

# determine which group a sample should be in
# add the sample ID to the list for that group
for i in range(len(clinical_df['Sample ID'])):
    if clinical_df['Overall Survival (Months)'][i] < 36:
        if clinical_df['Overall Survival Status'][i] == '1:DECEASED':
            group1.append(clinical_df['Sample ID'][i])
    else:
        group2.append(clinical_df['Sample ID'][i])

# get column names for new RNAseq dataframe (IDs & list of genes)
RNAseq_cols = ['Sample ID', 'Group ID']
for gene in RNAseq_df['sample']:
    RNAseq_cols.append(gene)

# get column names for new microarray dataframe (IDs & list of genes)
microarray_cols = ['Sample ID', 'Group ID']
for gene in microarray_df['sample']:
    microarray_cols.append(gene)

# create new dataframes to hold the preprocessed data
RNAseq_df2 = pd.DataFrame(columns=RNAseq_cols)
microarray_df2 = pd.DataFrame(columns=microarray_cols)

# populate new dataframes with relevant data
for sample_id in group1 + group2:
    # add sample ID and group ID (1/2)
    group = 1 if sample_id in group1 else 2
    row = {'Sample ID': sample_id, 'Group ID': group}

    # if the sample is in the RNAseq data, add gene expressions from RNAseq
    if sample_id in RNAseq_df.columns:
        for i in range(len(RNAseq_cols[2:])):
            row.update({RNAseq_cols[i+2]: RNAseq_df[sample_id][i]})
        RNAseq_df2 = pd.concat([RNAseq_df2, pd.DataFrame([row])], ignore_index=True)

    # if the sample is in the microarray data, add gene expressions from microarray
    elif sample_id in microarray_df.columns:
        for i in range(len(microarray_cols[2:])):
            row.update({microarray_cols[i+2]: microarray_df[sample_id][i]})
        microarray_df2 = pd.concat([microarray_df2, pd.DataFrame([row])], ignore_index=True)

# write data to text files
RNAseq_df2.to_csv('SeqData.txt', sep='\t', index=False)
microarray_df2.to_csv('ArrayData.txt', sep='\t', index=False)

# code to help answer part 1 questions
print("Genes in RNAseq gene expression profiles: ", len(RNAseq_cols) - 2) # subtract 2 because sample & group ID
print("Patient samples in RNAseq gene expression profiles: ", len(RNAseq_df2['Sample ID'])) # should be 227
print("Genes in microarray gene expression profiles: ", len(microarray_cols) - 2) # subtract 2 because sample & group ID
print("Patient samples in microarray gene expression profiles: ", len(microarray_df2['Sample ID'])) # should be 229