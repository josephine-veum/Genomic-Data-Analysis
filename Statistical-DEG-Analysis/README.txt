Project Overview:
This project applies statistical methods to identify differentially expressed genes between ovarian 
cancer patients with different survival statuses. It explores gene expression data from two genomic 
platforms: RNA sequencing (Illumina HiSeq 2000) and microarray (Affymetrix HT Human Genome U133a) to 
find markers associated with survival in TCGA ovarian serous cystadenocarcinoma (OV).

Data Sources:
Raw data can be accessed and downloaded through the links provided in the raw_data.txt file.
- RNAseq Data (Illumina HiSeq 2000): Data from UCSC Xena RNAseq
- Microarray Data (Affymetrix HT Human Genome U133a): Data from UCSC Xena Microarray
- Clinical Data: Clinical survival information for TCGA ovarian cancer patients Clinical data

Preprocessed Data:
- The two files correspond to the two genomic technologies used (RNAseq and microarray) and contain the following:
- Sample IDs
- Group IDs (1 or 2)
    - Group 1 contains patients who are dead with overall survivals less than 36 months.
    - Group 2 contains patients whose survivals are over 36 months without considering survival status.
- Gene expression data

Python Scripts:
- DataProcessing.py - loads clinical and gene expression data, processes them, and outputs 
  SeqData.txt (RNAseq data) and ArrayData.txt (microarray data).
- MyHPTest.py - implements the t-test and Wilcoxon rank-sum test to identify differentially expressed genes, 
  lists the top 10 genes, reports the number of selected genes, and plots a histogram of p-values.
- MyBonferroni.py - applies Bonferroni correction for multiple testing and reports the number of selected genes.
- MyFDR.py - applies FDR correction for multiple testing and outputs the upper bound of FDR for selecting 
  20, 50, 100, or 200 genes.
- OverlappedGenes.py - identifies overlapped genes from the RNAseq and microarray lists of significant genes, 
  plots the number of selected genes against the number of common genes for the t-tests and Wilcoxon rank-sum tests.

Graphs & Plots:
- ttest_rna.png - histogram of p-values and their frequencies for the RNAseq data when t-tests are used
  (from MyHPTest.py).
- ttest_array.png - histogram of p-values and their frequencies for the microarray data when t-tests are
  used (from MyHPTest.py).
- wilcoxon_rna.png - histogram of p-values and their frequencies for the RNAseq data when Wilcoxon rank-sum
  tests are used (from MyHPTest.py).
- wilcoxon_array.png - histogram of p-values and their frequencies for the microarray data when Wilcoxon 
  rank-sum tests are used (from MyHPTest.py).
- ttest_scatter.png - scatter plot of number of selected genes against the number of common genes between 
  the RNAseq and microarray data when t-tests are used (from OverlappedGenes.py).
- wilcoxon_scatter.png - scatter plot of number of selected genes agains the number of common genes between 
  the RNAseq and microarray data when Wilcoxon rank-sum tests are used (from OverlappedGenes).  
  
Libraries Used:
- pandas - data manipulation and analysis
- scipy - statistical analysis functions
- multipletests (from statsmodels.stats.multitest) - multiple testing correction methods (Bonferroni, FDR)
- random - random number generation
- matplotlib - plotting graphs and visualizations

How to Run:
- Clone or download this repository.
- Install the required libraries (listed above).
- Download the necessary raw data from the UCSC Xena and cBioPortal links (listed in raw_data.txt).
- Run the following scripts in order:
    - DataProcessing.py: Preprocesses the raw data and generates the SeqData.txt and ArrayData.txt files.
    - MyHPTest.py: Performs differential expression analysis with t-tests and Wilcoxon rank-sum tests, 
      lists top genes, and generates p-value histograms.
    - MyBonferroni.py: Applies Bonferroni correction for multiple testing.
    - MyFDR.py: Applies FDR correction for multiple testing.
    - OverlappedGenes.py: Identifies overlapped genes and generates the scatter plots.

Expected Output:
- Differentially Expressed Genes: A list of the top 10 genes from each test with their associated p-values 
  and histograms of p-values.
- Bonferroni and FDR Results: The number of genes selected after applying corrections.
- Overlap Analysis: Plots showing the overlap of differentially expressed genes between RNAseq and microarray data.