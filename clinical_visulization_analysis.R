clinical_raw_data <- read.table('../../Public_Data/TCGA_BRCA_clinical_simple/brca_tcga_pan_can_atlas_2018_clinical_data.tsv', header = T, 
                                stringsAsFactors = F, sep = '\t', row.names = 1)


# receptor info is based on 
receptor_raw_data <- read.table('../../Public_Data/TCGA_BRCA_clinical_oncoLnc/TCGA_BRCA_clinical_receptors.txt', 
                                header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
# based on receptor status split into TNBC, Hormone positive (either ER+ or PR+ and HER2-) and other
receptor_info = receptor_raw_data %>% 
  mutate(clinical_types = case_when(er_status_by_ihc == 'Negative' & pr_status_by_ihc == 'Negative' & her2_status_by_ihc == 'Negative' ~ 'TNBC',
                                    (er_status_by_ihc == 'Positive' | pr_status_by_ihc == 'Positive') & her2_status_by_ihc == 'Negative' ~ 'HR_postive',
                                    TRUE ~ 'other'))


MYC_V1_mRNA = input_cbio("../../Public_Data/brca_tcga_subset_tcga/MYC_V1/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt")
G2M_mRNA = input_cbio("../../Public_Data/brca_tcga_subset_tcga/G2M/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt")
E2F_mRNA = input_cbio("../../Public_Data/brca_tcga_subset_tcga/E2F/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt")
EMT_mRNA = input_cbio("../../Public_Data/brca_tcga_subset_tcga/EMT/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt")
Glyco_mRNA = input_cbio("../../Public_Data/brca_tcga_subset_tcga/Glycosis/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt")
TNF_alpha_mRNA = input_cbio("../../Public_Data/brca_tcga_subset_tcga/TNF_alpha/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt")

# split mRNA expression file based on clinical types of TNBC, HR positive and other 
MYC_V1_mRNA_split = df_split(MYC_V1_mRNA, receptor_info, "clinical_types", tosplit = T)
G2M_mRNA_split = df_split(G2M_mRNA, receptor_info, "clinical_types", tosplit = T)
E2F_mRNA_split = df_split(E2F_mRNA, receptor_info, "clinical_types", tosplit = T)
EMT_mRNA_split = df_split(EMT_mRNA, receptor_info, "clinical_types", tosplit = T)
Glyco_mRNA_split = df_split(Glyco_mRNA, receptor_info,"clinical_types", tosplit=T)
TNF_alpha_mRNA_split = df_split(TNF_alpha_mRNA, receptor_info, "clinical_types", tosplit = T)


# T-test between TNBC and HR positive by genes
ttest_by_row(t(MYC_V1_mRNA_split$TNBC), t(MYC_V1_mRNA_split$HR_postive), col_result = c("tValue", "pValue", "mean_of_TNBC", "mean_of_HR_positive"))
ttest_by_row(t(G2M_mRNA_split$TNBC), t(G2M_mRNA_split$HR_postive), col_result = c("tValue", "pValue", "mean_of_TNBC", "mean_of_HR_positive"))
ttest_by_row(t(E2F_mRNA_split$TNBC), t(E2F_mRNA_split$HR_postive), col_result = c("tValue", "pValue", "mean_of_TNBC", "mean_of_HR_positive"))
ttest_by_row(t(EMT_mRNA_split$TNBC), t(EMT_mRNA_split$HR_postive), col_result = c("tValue", "pValue", "mean_of_TNBC", "mean_of_HR_positive"))
ttest_by_row(t(Glyco_mRNA_split$TNBC), t(Glyco_mRNA_split$HR_postive), col_result = c("tValue", "pValue", "mean_of_TNBC", "mean_of_HR_positive"))
ttest_by_row(t(TNF_alpha_mRNA_split$TNBC), t(TNF_alpha_mRNA_split$HR_postive), col_result = c("tValue", "pValue", "mean_of_TNBC", "mean_of_HR_positive"))
# write results to files
# write.csv( ttest_subtypes(MYC_V1_mRNA), "./Processed_data_output/MYC_V1_ttest.csv")
# write.csv( ttest_subtypes(G2M_mRNA), "./Processed_data_output/G2M_ttest.csv")
# write.csv( ttest_subtypes(E2F_mRNA), "./Processed_data_output/E2F_ttest.csv")
# write.csv( ttest_subtypes(EMT_mRNA), "./Processed_data_output/EMT_ttest.csv")
# write.csv( ttest_subtypes(Glyco_mRNA), "./Processed_data_output/Glycosis_ttest.csv")
# write.csv( ttest_subtypes(TNF_alpha_mRNA), "./Processed_data_output/TNF_alpha_ttest.csv" )







# heatmap for mRNA of TCGA BRCA
# MYC_V1 
clinical_subtypes(MYC_V1_mRNA, receptor_info)
# G2M 
clinical_subtypes(G2M_mRNA, receptor_info)
# E2F 
clinical_subtypes(E2F_mRNA, receptor_info)
# EMT 
clinical_subtypes(EMT_mRNA, receptor_info)
# Glycosis 
clinical_subtypes(Glyco_mRNA, receptor_info)
# TNF-alpha 
clinical_subtypes(TNF_alpha_mRNA, receptor_info)
