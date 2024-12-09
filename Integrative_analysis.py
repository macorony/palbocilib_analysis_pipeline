import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')
from palbociclib_signature_analysis import palbo_signatures
from palbo_RNAseq_analysis import palbo_RNAseq

def build_library(fpath="./input_data/Enrichment/MSigDB_Hallmark_2020.txt"):
    """Build enrichment library for different gene sets"""
    with open(fpath, 'r') as f:
        lines = f.readlines()
    
    library = {}
    for line in lines:
        parts = line.strip().split('\t')
        new_parts = []
        for part in parts:
            if part != '':
                new_parts.append(part)
        library[new_parts[0]] = new_parts[1:]
    return library


def geneset_heatmap(gene_set, signature_data, seq_data, save_csv=False):
    """Analyze gene set in RNA-seq and drug signature datasets"""
    set_in_signature = signature_data[signature_data['Partner'].isin(gene_set)]
    set_in_seq = seq_data[seq_data['GeneName'].isin(gene_set)]
    
    overlap = pd.merge(set_in_signature, set_in_seq, left_on='Partner', right_on='GeneName')
    
    print(f"Genes in CARE data: {len(set_in_signature)}")
    print(f"Genes in sequencing data: {len(set_in_seq)}")
    print(f"Overlapping genes: {len(overlap)}")
    
    # Create heatmap
    plot_data = overlap[['log2FoldChange', 't_value_CDK6', 't_value_RB1']].copy()
    plot_data.index = overlap['Partner']
    
    plt.figure(figsize=(8, 10))
    sns.clustermap(plot_data, cmap='RdYlBu_r', center=0)
    plt.title('Gene Expression Patterns')
    
    # if save_csv:
        # output_name = f"{gene_set.name}_presence.csv"
        # overlap.to_csv(f"./Processed_data_output/{output_name}")
    
    # return overlap

# 
if __name__ == "__main__":
    # Initialize analysis  
    fpath1="./input_data//palbo_signature/CCLE_palbo_CDK6.csv"
    fpath2="./input_data/palbo_signature/CCLE_palbo_RB1.csv" 
    palbo_sig = palbo_signatures(fpath1, fpath2)
    palbo_seq = palbo_RNAseq("./input_data/palbo_RNAseq/palbociclib_gene_list.txt")
    
    geneset_library = build_library
    pathways2plot = ["Myc Targets V1","G2-M Checkpoint","E2F Targets",
                     "Epithelial Mesenchymal Transition", "TNF-alpha Signaling via NF-kB",
                     "Glycolysis"]
    for each in pathways2plot:
        geneset_heatmap(geneset_library[each], palbo_sig.combo, palbo_seq.gene_summary)