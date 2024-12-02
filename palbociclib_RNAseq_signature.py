import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
import plotly.express as px
from statsmodels.stats.multitest import fdrcorrection
import warnings
warnings.filterwarnings('ignore')

def read_and_process_data():
    # CARE data input
    ccle_cdk6 = pd.read_csv("./CARE/CCLE_palbo_CDK6.csv")
    ccle_rb1 = pd.read_csv("./CARE/CCLE_palbo_RB1.csv")
    
    # RNA-Seq data input
    pal_seq = pd.read_csv("./RNA-seq/palbociclib_gene_list.txt", sep='\t')
    cdk4ko_seq = pd.read_csv("./RNA-seq/cdk4ko_seq_gene_list.txt", sep='\t')
    cdk6ko_seq = pd.read_csv("./RNA-seq/cdk6ko_seq_gene_list.txt", sep='\t')
    
    # CRISPR screen data
    crispr_gene_summary = pd.read_csv('./control_normalized_working/palbo_vs_dmso_control_normalized.gene_summary.txt')
    
    # Filter sequencing data by fold change
    pal_seq_select = pal_seq[abs(pal_seq['log2FoldChange']) >= 0.58][['GeneName', 'log2FoldChange', 'padj']]
    cdk4ko_seq_select = cdk4ko_seq[abs(cdk4ko_seq['log2FoldChange']) >= 0.58][['GeneName', 'log2FoldChange', 'padj']]
    cdk6ko_seq_select = cdk6ko_seq[abs(cdk6ko_seq['log2FoldChange']) >= 0.58][['GeneName', 'log2FoldChange', 'padj']]
    
    return ccle_cdk6, ccle_rb1, pal_seq_select, cdk4ko_seq_select, cdk6ko_seq_select, crispr_gene_summary

def mutation_status(data, pattern=r"^.*Mutation$"):
    """Add mutation status column"""
    data['mutation'] = data['Partner'].str.match(pattern)
    return data

def plot_volcano(data, title):
    """Create volcano plot"""
    plt.figure(figsize=(10, 8))
    plt.scatter(data['t_value'], -np.log10(data['p_value']), 
               c=data['mutation'].map({True: 'red', False: 'blue'}),
               alpha=0.6)
    plt.xlabel('t-value')
    plt.ylabel('-log10(p-value)')
    plt.title(title)
    plt.show()

def cross_heatmap(df1, df2, col2draw=None):
    """Generate correlation heatmap between two datasets"""
    if col2draw is None:
        col2draw = ['t_value', 'log2FoldChange']
    
    cross_data = pd.merge(df1, df2, left_on='Partner', right_on='GeneName', how='inner')
    
    # Correlation matrix
    cor_matrix = cross_data[col2draw].corr()
    
    # Plot heatmaps
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    
    # Correlation heatmap
    sns.heatmap(cor_matrix, ax=ax1, cmap='RdBu_r', vmin=-1, vmax=1, 
                annot=True, square=True)
    ax1.set_title('Correlation Matrix')
    
    # Data heatmap
    sns.heatmap(cross_data[col2draw], ax=ax2, cmap='RdBu_r',
                center=0, cbar_kws={'label': 't-score | lfc'})
    ax2.set_title('Data Values')
    
    plt.tight_layout()
    return cross_data

def set_presence(gene_set, care_data, seq_data, save_csv=False):
    """Analyze gene set presence in different datasets"""
    set_at_care = care_data[care_data['Partner'].isin(gene_set)]
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

def main():
    # Read data
    ccle_cdk6, ccle_rb1, pal_seq_select, cdk4ko_seq_select, cdk6ko_seq_select, crispr_gene_summary = read_and_process_data()
    
    # Add mutation status
    ccle_cdk6 = mutation_status(ccle_cdk6)
    ccle_rb1 = mutation_status(ccle_rb1)
    
    # Create volcano plots
    plot_volcano(ccle_cdk6, 'CCLE CDK6 Volcano Plot')
    plot_volcano(ccle_rb1, 'CCLE RB1 Volcano Plot')
    
    # Cross dataset analysis
    cross_data = cross_heatmap(ccle_cdk6, pal_seq_select)
    
    return cross_data

if __name__ == "__main__":
    main()
