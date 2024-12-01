import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy import stats
from adjustText import adjust_text




class CRISPRScreenAnalysis:
    def __init__(self, gene_summary_path, normalized_count_path):
        
        """Initialize CRISPR screen analysis with data files."""
        
        # initialize attributes
        self.gene_summary = pd.read_csv(gene_summary_path, sep='\t')
        self.normalized_count = pd.read_csv(normalized_count_path, sep='\t', index_col=0)
        self.normalized_count = self.normalized_count.iloc[:, 1:]
        # self.process_data()
        
        # initialize storge attributes
        self.combined_data = None
        self.pos_sig_genes = None
        self.neg_sig_genes = None

        self.process_data()

    # def verify_data(self):
   
    def process_data(self):
        """Process and organize CRISRP screen gene summary data."""
        
        # Split positive and negative data
        positive_init = self.gene_summary[self.gene_summary['pos|lfc'] >= 0]
        negative_init = self.gene_summary[self.gene_summary['pos|lfc'] < 0]
        
        # Process positive data
        pos_cols = ['id', 'num', 'pos|score', 'pos|p-value', 'pos|fdr', 
                'pos|rank', 'pos|goodsgrna', 'pos|lfc']
        positive = positive_init[pos_cols]
        # Process negative data
        neg_cols = ['id', 'num', 'neg|score', 'neg|p-value', 'neg|fdr', 
                   'neg|rank', 'neg|goodsgrna', 'neg|lfc']
        negative = negative_init[neg_cols]
        # Rename columns
        rename = ['gene_symbol', 'num', 'score', 'p_value', 'fdr', 
                 'rank', 'goodsgrna', 'lfc']
        positive.columns = rename
        negative.columns = rename
        
        # Combine and process data
        self.combined_data = pd.concat([positive, negative])
        self.combined_data['logFDR'] = -np.log10(self.combined_data['fdr'])
        
        # Filter miRNA entries
        self.combined_data = self.combined_data[
            ~self.combined_data['gene_symbol'].str.contains("hsa-", na=False)
        ]
        print(self.combined_data)

        # Identify significant genes
        self.pos_sig_genes = self.combined_data[
            (self.combined_data['fdr'] <= 0.05) & (self.combined_data['lfc'] > 0)
        ]['gene_symbol']
        self.neg_sig_genes = self.combined_data[
            (self.combined_data['fdr'] <= 0.05) & (self.combined_data['lfc'] < 0)
        ]['gene_symbol']
        
    
    def create_volcano_plot(self, top_n=10, figsize=(12, 8)):
        """Create enhanced volcano plot with labeled top genes."""
        plt.figure(figsize=figsize)
        
        # Plot all points
        sns.scatterplot(data=self.combined_data, x='lfc', y='logFDR', 
                       alpha=0.2, label='Non-significant')
        
        # Plot significant points
        sns.scatterplot(data=self.combined_data[
            self.combined_data['gene_symbol'].isin(self.pos_sig_genes)
        ], x='lfc', y='logFDR', color='red', label='Positive significant')
        
        sns.scatterplot(data=self.combined_data[
            self.combined_data['gene_symbol'].isin(self.neg_sig_genes)
        ], x='lfc', y='logFDR', color='blue', label='Negative significant')
        
        # Add significance threshold line
        plt.axhline(-np.log10(0.05), color='gray', linestyle='--', 
                    linewidth=1, label='Significance threshold')
        
        # Label top genes
        top_genes = self.combined_data.nlargest(top_n, 'logFDR')
        texts = []
        for _, gene in top_genes.iterrows():
            texts.append(plt.text(gene['lfc'], gene['logFDR'], 
                                  gene['gene_symbol']))
        
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray', lw=0.5))
        
        plt.title('CRISPR Screen Volcano Plot', pad=20)
        plt.xlabel('Log Fold Change (LFC)')
        plt.ylabel('-Log10(FDR)')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        return plt.gcf()
        
    def create_correlation_heatmap(self, figsize=(10, 8)):
        """Create correlation heatmap with enhanced visualization."""
        plt.figure(figsize=figsize)
        
        # Calculate correlation matrix
        corr_matrix = self.normalized_count.corr()
        
        # Create mask for upper triangle
        mask = np.triu(np.ones_like(corr_matrix), k=1)
        
        # Create heatmap
        sns.heatmap(corr_matrix, vmin=0.95, vmax=1, annot=True, 
                   cmap='coolwarm', square=True, mask=mask,
                   fmt='.3f', cbar_kws={'label': 'Correlation'})
        
        plt.title('Sample Correlation Heatmap', pad=20)
        plt.tight_layout()
        return plt.gcf()
        
    def perform_pca(self, figsize=(12, 10)):
        """Perform PCA analysis with multiple visualizations."""
        # Perform PCA
        pca = PCA()
        pca_results = pca.fit_transform(self.normalized_count.T)
        
        # Create figure with subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
        
        # Plot PCA scatter
        ax1.scatter(pca_results[:, 0], pca_results[:, 1])
        for i, sample in enumerate(self.normalized_count.columns):
            ax1.annotate(sample, (pca_results[i, 0], pca_results[i, 1]))
        
        ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
        ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
        ax1.set_title('PCA Sample Distribution')
        
        # Plot explained variance
        explained_var = pd.DataFrame({
            'PC': [f'PC{i+1}' for i in range(len(pca.explained_variance_ratio_))],
            'Explained Variance': pca.explained_variance_ratio_
        })
        sns.barplot(data=explained_var, x='PC', y='Explained Variance', ax=ax2)
        ax2.set_title('Explained Variance by Principal Components')
        
        plt.tight_layout()
        return fig
        
    def get_significant_genes(self):
        """Return dictionary of significant genes and their statistics."""
        
        return {
            'positive': self.pos_sig_genes,
            'negative': self.neg_sig_genes,
            'summary': {
                'total_significant': len(self.pos_sig_genes) + len(self.neg_sig_genes),
                'positive_hits': len(self.pos_sig_genes),
                'negative_hits': len(self.neg_sig_genes)
            }
        }

# Usage example
if __name__ == "__main__":
    # Initialize analysis
    gene_summary_fpath = "./input_data/CRISPR_screen/palbo_vs_dmso_control_normalized.gene_summary.txt"
    normalized_counts_fpath = "./input_data/CRISPR_screen/palbo_vs_dmso_control_normalized.normalized.txt"
    palbo_analysis = CRISPRScreenAnalysis(
        gene_summary_fpath, normalized_counts_fpath
        )
    
    # Create visualizations
    volcano_plot = palbo_analysis.create_volcano_plot()
    correlation_heatmap = palbo_analysis.create_correlation_heatmap()
    pca_plots = palbo_analysis.perform_pca()
    
    # Get significant genes
    significant_results = palbo_analysis.get_significant_genes()
    
    # Display results
    print("\nSignificant Gene Summary:")
    print(f"Total significant genes: {significant_results['summary']['total_significant']}")
    print(f"Positive hits: {significant_results['summary']['positive_hits']}")
    print(f"Negative hits: {significant_results['summary']['negative_hits']}")
    
    plt.show()