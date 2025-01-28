# TF-Activity-Omics-Integration
Code of the paper "Integrating gene expression, genomic and phosphoproteomic data to infer transcription factor activity in lung cancer"

**Overview**  
In this analysis, we performed a Master Regulator Analysis (MRA) to assess transcription factor (TF) activity in lung adenocarcinoma samples (1). We adapted the pipeline of the corto R package (2), modifying it to focus on direct target regulons. To identify these regulons, we filtered the target lists of each TF to include only direct targets by leveraging the SEanalysis database (3). The MRA analysis generates a Normalized Enrichment Score (NES) for each TF, indicating whether its direct target genes are significantly up- or down-regulated in tumor tissues compared to normal adjacent tissues. This NES serves as a proxy for the overall change in TF activity across all tumor samples. 
Additionally, we applied LASSO regression to examine the relationship between TF activity and phosphorylation levels in patient-specific phosphoproteomic profiles. The goal is to identify key transcription factors whose activity correlates with changes in protein phosphorylation, which may provide insights into the regulation of critical signaling pathways in the tumor microenvironment.

## Requirements  
Python 3.4.1  
R version 4.4.1  

## Master_regulator_Analysis_using_corto.R  

**Input Data**  
centroids_TFs.txt: list of transcription factors 
RNA-seq_GNAME_VST_TUM.tsv and RNA-seq_GNAME_VST_NAT.tsv: dataset of RNA-seq data normalized by vst R function

**Output Data**  

TUMvsNAT_reg_NAT-regulon_onlyDBCK.rda: regulon of only direct target genes for every TF  
mra_TUMvsNAT_reg_NAT-regulon_onlyDBCK.rda: result of MRA analysis performed on tumor versus normal samples  
mra_TUMvsNAT_top10_vst_regulon_onlyDBCK.png: plot of top 10 TF with highest NES values.  

## LASSO_regression.R  

**Input Data**  
TF_Activity_sum_pos_sub_neg.tsv: transcription factors activity calculated as the sum of the vst expression of positive direct target genes of a TF, subtracting the sum of the vst expression of negative target genes of a TF.  
MSigDB_TF_gene_phosphorylations_pathways_TN_FC.txt: dataframe with the phosphosites and the substrates who partecipate in the same pathway with a TF. The pathway are derived from MSigDB database (4).  
Phospho_TN_FC.txt: A CSV (or appropriate format) containing phosphorylation levels for various phosphosites across multiple patient samples.   

**Output Data**  
Coef_phospho_TN_FC_onlyDBCK_NES_SumPosSubNeg_LassoOnPathwayPhosphorylations_TN_FC.txt: result of LASSO regression with TFs, phsphosites and their contribution to the prediction of the TF activity  
RMSE_Rsquared_phospho_TN_FC_TF_onlyDBCK_NES_SumPosSubNeg_LassoOnPathwayPhosphorylations_TN_FC.txt: result of LASSO regression with Rsquared, Root Mean Square Error (RMSE) and the TFs of the regressions.  
