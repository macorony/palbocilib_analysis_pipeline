import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

class palbo_signatures:
    def __init__(self, cdk6_fpath, rb1_fpath):
    # check funtion for correct corresponding data
        self.CDK6 = pd.read_csv(cdk6_fpath)
        self.RB1 = pd.read_csv(rb1_fpath)
        self.combo = pd.merge(self.CDK6, self.RB1, on='Partner', how='inner',
                              suffixes=('_CDK6', '_RB1'))
        self.mutation_status()
    
    def mutation_status(self, pattern=r"^.*Mutation$"):
        """add mutation status"""
        self.CDK6['mutation'] = self.CDK6['Partner'].str.match(pattern)
        self.RB1['mutation'] = self.RB1['Partner'].str.match(pattern)

    def _plot_scatter(self, data, title):
        """Create volcano plot"""
        plt.figure(figsize=(10, 8))
        plt.scatter(data['t_value'], -np.log10(data['p_value']), 
                c=data['mutation'].map({True: 'red', False: 'blue'}),
                alpha=0.6)
        plt.xlabel('t-value')
        plt.ylabel('-log10(p-value)')
        plt.title(title)
        plt.show()
    
    def plot_volcano(self, target):
        if target == "CDK6":
            self._plot_scatter(self.CDK6, "Palbociclib signatures when CDK6 as primary target")
        elif target == "RB1":
            self._plot_scatter(self.RB1, "Palbociclib signatures when RB1 as primary target")
        else:
            raise ValueError("The target must be either CDK6 or RB1")

    def cross_heatmap(self):
        """Generate heatmap between two signature datasets"""        
        # merge two datasets
        cross_data = pd.merge(self.CDK6, self.RB1, on='Partner', how='inner',
                              suffixes=('_CDK6', '_RB1'))
        
        # Correlation matrix
        # cor_matrix = cross_data[col2draw].corr()
        col2draw = cross_data.columns.str.match("t_value_*")
        # data heatmaps
        plt.figure(figsize=(10, 8))
        sns.clustermap(cross_data.loc[:, col2draw], cmap='RdBu_r', row_cluster=True, col_cluster=False,
                    center=0, cbar_kws={'label': 'CDK6 | RB1'})
        plt.title('T score heatmap')
        
        plt.tight_layout()
    

if __name__ == "__main__":
    # Initialize analysis
    CDK6_fpath="./input_data//palbo_signature/CCLE_palbo_CDK6.csv"
    RB1_fpath="./input_data/palbo_signature/CCLE_palbo_RB1.csv" 
    palbo_sig = palbo_signatures(CDK6_fpath, RB1_fpath)

    # Create visualizations
    palbo_sig.plot_volcano()
    palbo_sig.cross_heatmap()
    plt.show()