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


def set_presence(gene_set, signature_data, seq_data, save_csv=False):
    """Analyze gene set presence in different datasets"""
    set_at_care = signature_data[signature_data['Partner'].isin(gene_set)]
    set_at_seq = seq_data[seq_data['GeneName'].isin(gene_set)]
    
    overlap = pd.merge(set_at_care, set_at_seq, 
                      left_on='Partner', right_on='GeneName')
    
    print(f"Genes in CARE data: {len(set_at_care)}")
    print(f"Genes in sequencing data: {len(set_at_seq)}")
    print(f"Overlapping genes: {len(overlap)}")
    
    # Create heatmap
    plot_data = overlap[['log2FoldChange', 't_value']].copy()
    plot_data.index = overlap['Partner']
    
    plt.figure(figsize=(8, 10))
    sns.heatmap(plot_data, cmap='RdYlBu_r', center=0)
    plt.title('Gene Expression Patterns')
    
    if save_csv:
        output_name = f"{gene_set.name}_presence.csv"
        overlap.to_csv(f"./Processed_data_output/{output_name}")
    
    return overlap

def build_library(fpath="./Enrichment/library/MSigDB_Hallmark_2020.txt"):
    """Build enrichment library"""
    with open(fpath, 'r') as f:
        lines = f.readlines()
    library = {}
    for line in lines:
        parts = line.strip().split('\t')
        library[parts[0]] = parts[1:]
    return library

class ComplexHeatmap:
    """Python implementation of R's ComplexHeatmap"""
    def __init__(self, data, split_cols=None, annotations=None):
        self.data = data
        self.split_cols = split_cols
        self.annotations = annotations
    
    def plot(self):
        fig = plt.figure(figsize=(12, 8))
        
        if self.annotations:
            gs = plt.GridSpec(len(self.annotations) + 1, 1, height_ratios=[1]*len(self.annotations) + [4])
            
            for i, (name, values) in enumerate(self.annotations.items()):
                ax = fig.add_subplot(gs[i])
                if isinstance(values, pd.Series):
                    values.plot(kind='bar', ax=ax)
                else:
                    ax.imshow([values], aspect='auto')
                ax.set_title(name)
        
        ax_main = fig.add_subplot(gs[-1])
        sns.heatmap(self.data, ax=ax_main, cmap='RdBu_r')
        
        if self.split_cols:
            for split in np.unique(self.split_cols):
                plt.axvline(x=(self.split_cols == split).sum(), color='black', linewidth=2)
        
        plt.tight_layout()
        return fig


