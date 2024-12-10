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


input_cbio = function(fpath, re_na_col=T) {
  # fpath: file path
  # re_na_col: remove all NA columns
  raw_data = read.csv(fpath, sep='\t', header = T, row.names = 2, check.names = F, stringsAsFactors = F)
  # remove the study id and transpose
  raw_data = t(raw_data[, -1])
  colnames(raw_data) = lapply(colnames(raw_data), FUN=sampletopatient)
  if (re_na_col) {
    return(remove_na_column(raw_data))
  } else {
    return(raw_data)
  } 
}

df_split = function(df, group_info, col2split, tosplit=F, order_info=F) {    
  # df: dataframe to be split  patient id at column and gene symbol at row   
  # col2split: specifiy the column containing the subtype information
  # tosplit:to keep the resulting data integral or split into subset
  # order_info: to rank dataframe
  df_merge = merge(t(df), group_info[, col2split, drop=F], by=0, all.x=T)
  df_merge[,col2split] = replace_na(df_merge[,col2split], "unknown")
  df_merge = df_merge[order(df_merge[, col2split]), ]
  rownames(df_merge) = df_merge[,1]
  df_ready = df_merge[, -c(1, ncol(df_merge))]
  split_info = df_merge[, ncol(df_merge), drop=F]
  if (tosplit) {
    df_splitted = split(df_ready, split_info)
    return(df_splitted)
  } else { 
    return(list(results=df_ready, split_info=split_info))
    }
  }
clinical_subtypes = function(mRNA_data, clinical_info=receptor_info){
  mRNA_TNBC = mRNA_data[, colnames(mRNA_data) %in% rownames(clinical_info[clinical_info$clinical_types == 'TNBC',])]
  mRNA_HR = mRNA_data[, colnames(mRNA_data) %in% rownames(clinical_info[clinical_info$clinical_types == 'HR+',])]
  mRNA_other = mRNA_data[, colnames(mRNA_data) %in% rownames(clinical_info[clinical_info$clinical_types == 'other',])]
  col_map1 = colorRamp2(c(-4, 0, 4), c( "#123456", "white", "#FF0000")   )
  ha_top1 = HeatmapAnnotation(TNBC=rep('a', dim(mRNA_TNBC)[2]),
                              col = list(TNBC= c('a'="orange")), 
                              annotation_name_side = "left", show_legend = F)
  ht1 = Heatmap(mRNA_TNBC, col = col_map1, top_annotation = ha_top1,
                show_column_names = F, show_row_dend = F, show_column_dend = F,
                cluster_columns = F)
  ha2 = HeatmapAnnotation(HR_positive=rep("b", dim(mRNA_HR)[2]),
                          col = list(HR_positive=c("b" = "cyan")), 
                          show_annotation_name = F,
                          show_legend = F)
  ht2 = Heatmap(mRNA_HR, col= col_map1, top_annotation = ha2,
                show_column_names = F, show_row_dend = F, show_column_dend = F, 
                cluster_columns = F, 
                show_heatmap_legend = F)
  ha3 = HeatmapAnnotation(other=rep("c", dim(mRNA_other)[2]),
                          col =list(other=c("c" = "black")), 
                          annotation_name_side = 'right', show_legend = F)
  ht3 = Heatmap(mRNA_other, col= col_map1, top_annotation = ha3,
                show_column_names = F, show_row_dend = F, show_column_dend = F, 
                cluster_columns = F, 
                show_heatmap_legend = F)
  
  return(ht1 + ht2 + ht3)
}


ttest_subtypes <- function(mRNA_data, clinical_info=receptor_info){
  mRNA_TNBC = mRNA_data[, colnames(mRNA_data) %in% rownames(clinical_info[clinical_info$clinical_types == 'TNBC',])]
  mRNA_HR = mRNA_data[, colnames(mRNA_data) %in% rownames(clinical_info[clinical_info$clinical_types == 'HR+',])]
  return(ttest_by_row(mRNA_TNBC, mRNA_HR, col_result=c("tValue", "pValue", "mean of TNBC", "mean of HR+")))
}

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
