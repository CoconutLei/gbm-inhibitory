ssGSEA.R aims to use GSVA package to analyze expression level of specific gene sets. 
Demo: Use "CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt" as input, the code will ouput the expression level of "EN" and "IN" gene sets across CGGA samples. Expected rum time is less than 10 seconds.

CGGA_survival.R aims to analyze the impact of specific genes or gene sets on  overall survival of glioblastoma patients.
Demo: Use "CGGA.mRNAseq_693_clinical.20200506.txt" as input, the code will ouput the impact of "EN" gene sets on  overall survival of glioblastoma patients. Expected rum time is less than 10 seconds.

findPeak.R aims to mark peak frame of calcium imaging data.
Demo: Use "coactiveB.csv" as input, the code will ouput the peak frame. Expected rum time is less than 10 seconds.

Co-activated.R aims to calculate similarity between cell pairs of calcium imaging data.
Demo: Use "coactiveB.csv" as input, the code will ouput similarity matrix. Expected rum time is less than 10 seconds.

screen.R aims to screen calcium related candidate genes in glioblastoma.
Demo: Use "CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt", "CGGA.mRNAseq_693_clinical.20200506.txt" ,"res_df.csv", "GBM_gene_log_rank_results.csv" ,"maligdf.csv" and "calcium related genes.csv" as input, the code will ouput the candidate genes. Expected rum time is less than 1 minute.
