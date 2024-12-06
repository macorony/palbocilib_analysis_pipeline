import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')


class palbo_RNAseq:
    def __init__(self, fpath):
        self._raw_data = pd.read_table(fpath)
        self.gene_summary = self._raw_data[abs(self._raw_data['log2FoldChange']) >= 0.58][['GeneName', 'log2FoldChange', 'padj']]

    def plot_volcano(self):
        """Create volcano plot"""
        plt.figure(figsize=(10, 8))
        plt.scatter(self.gene_summary['log2FoldChange'], -np.log10(self.gene_summary['padj']), 
                alpha=0.6)
        plt.xlabel('Log Fold Change')
        plt.ylabel('-log10(adjusted P-value)')
        plt.show()