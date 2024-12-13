import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Optional, Dict, List, Tuple, Union
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class DataPaths:
    """Centralize data paths configuration"""
    BASE_PATH: Path = Path('./input_data')
    CLINICAL_PATH: Path = BASE_PATH / 'clinical_features/brca_tcga_pan_can_atlas_2018_clinical_data.tsv'
    RECEPTOR_PATH: Path = BASE_PATH / 'clinical_features/TCGA_BRCA_clinical_receptors.txt'
    MRNA_BASE: Path = BASE_PATH / 'mRNA_data'

# split the data based molecular and clinical subtype

class TCGA_GeneExpressionData:
    """Handle gene expression data processing and analysis"""
    def __init__(self):
        self.clinical_feature = None
        self.receptor_feature = None
        self.gene_sets = {}
        
    def load_clinical_data(self, paths: DataPaths) -> None:
        """Load and process clinical data"""
        try:
            self.clinical_feature = pd.read_table(paths.CLINICAL_PATH, index_col=0)
            # NaN value as Unidentified 
            self.clinical_feature.fillna("Uni", inplace=True)
            self.receptor_feature = pd.read_table(paths.RECEPTOR_PATH, index_col=0)
            self._process_receptor_data()
        except Exception as e:
            logger.error(f"Error loading clinical data: {e}")
            raise
            
    def _process_receptor_data(self) -> None:
        """Process receptor status into clinical types of TNBC, Hormone Receptor(HR) positive and other"""
        conditions = [
            (self.receptor_feature['er_status_by_ihc'].eq('Negative') & 
             self.receptor_feature['pr_status_by_ihc'].eq('Negative') & 
             self.receptor_feature['her2_status_by_ihc'].eq('Negative')),
            ((self.receptor_feature['er_status_by_ihc'].eq('Positive') | 
              self.receptor_feature['pr_status_by_ihc'].eq('Positive')) & 
             self.receptor_feature['her2_status_by_ihc'].eq('Negative'))
        ]
        choices = ['TNBC', 'HR_positive']
        self.receptor_feature['clinical_types'] = np.select(conditions, choices, default='other')
    
    @staticmethod
    def preprocess_mRNA(fpath: Union[str, Path], remove_na: bool = True) -> pd.DataFrame:
        """Load and process cBioPortal data"""
        try:
            data = pd.read_table(fpath, index_col=1).iloc[:, 1:].T
            data.columns = [col[:12] for col in data.columns]  # sample to patient conversion
            return data.dropna(axis=1, how='all') if remove_na else data
        except Exception as e:
            logger.error(f"Error processing cBioPortal data from {fpath}: {e}")
            raise
            
    def load_gene_sets(self, paths: DataPaths) -> None:
        """Load all gene set expression data"""
        gene_sets = {
            'MYC_V1': 'MYC_V1/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt',
            'G2M': 'G2M/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt',
            'E2F': 'E2F/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt',
            'EMT': 'EMT/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt',
            'Glyco': 'Glycosis/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt',
            'TNF_alpha': 'TNF_alpha/mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt'
        }
        
        for name, file_path in gene_sets.items():
            full_path = paths.MRNA_BASE / file_path
            self.gene_sets[name] = self.preprocess_mRNA(full_path)

    # merge two methods into one 
    def split_by_types(self, data: pd.DataFrame, by_type: Union[str, List[str]] = "clinical") -> Dict[str, pd.DataFrame]:
        """
        Split DataFrame by clinical or molecular subtypes.
        
        Args:
            data: DataFrame to split
            by_type: "clinical" for TNBC/HR classification or "molecular" for BRCA subtypes
            
        Returns:
            Dictionary of split DataFrames by subtype
        """
        TYPE_MAPPINGS = {
            "clinical": {
                "feature": "receptor_feature",
                "column": "clinical_types",
                "categories": ['TNBC', 'HR_positive', 'other']
            },
            "molecular": {
                "feature": "clinical_feature",
                "column": "Subtype",
                "categories": ['BRCA_LumA', 'BRCA_LumB', 'BRCA_Normal', 'BRCA_Basal', 'Uni']
            }
        }
        
        if by_type not in TYPE_MAPPINGS:
            raise ValueError(f"by_type must be one of {list(TYPE_MAPPINGS.keys())}")
        
        type_info = TYPE_MAPPINGS[by_type]
        feature_data = getattr(self, type_info["feature"])

    
        splits = {}
        for category in type_info["categories"]:
            mask = feature_data[type_info["column"]] == category
            patient_ids = feature_data[mask].index
            splits[category] = data.loc[:, data.columns.isin(patient_ids)]
        
        return splits
    
    def run_one_way_anova(self, data: pd.DataFrame, by_type: str = "clinical"):
        """
        Perform one-way ANOVA for each gene across groups
        
        Args:
            data: Gene expression matrix (genes × samples)
            by_type: Type of grouping to use ("clinical" or molecular)

        Returns:
            DataFrame with ANOVA results including F-statistic, p-value, and group means
        """
        if by_type == "clinical":
            groups = self.receptor_feature["clinical_types"]
        elif by_type == "molecular":
            groups = self.clinical_feature["Subtype"]
        else:
            raise ValueError("by_type must be either 'clinical' or 'molecular'")
        
        common_samples = groups.index.intersection(data.columns)
        groups = groups[common_samples]
        data = data[common_samples]
        results = []
        for gene in data.index:
            # Create lists of expression values for each group
            group_data = {group: data.loc[gene, groups == group] for group in groups.unique()}
            # Calculate group means
            group_means = {f"mean_{group}": values.mean() 
                          for group, values in group_data.items()}
            
            # Run ANOVA
            f_stat, p_val = stats.f_oneway(*group_data.values())
            
            # Add to results
            results.append({
                'gene': gene,
                'F_statistic': f_stat,
                'p_value': p_val,
                **group_means  # Include all group means in results
            })
        results_df = pd.DataFrame(results).set_index('gene')
        results_df = results_df.sort_values('p_value')
        # Add multiple testing correction
        results_df['adjusted_p_value'] = stats.false_discovery_control(results_df['p_value'])
        
        return results_df
    
    def analyze_all_genesets_anova(self, by_type: str = "clinical", output_dir: str = "./Processed_data_output") -> Dict[str, pd.DataFrame]:
        """
        Run ANOVA analysis for all gene sets and combine results.
        
        Args:
            by_type: Type of grouping to use ("clinical" or "molecular")
            output_dir: Directory to save results
            
        Returns:
            Dictionary containing ANOVA results for each gene set
        """
        if not self.gene_sets:
            raise ValueError("Gene sets not loaded. Please call load_gene_sets() first.")
        
        # Create output directory if it doesn't exist
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Store results for each gene set
        all_results = {}
        
        # Combined significant results across gene sets
        significant_genes = pd.DataFrame()
        
        for name, data in self.gene_sets.items():
            logger.info(f"Running ANOVA analysis for {name} gene set")
            
            # Run ANOVA
            anova_results = self.run_one_way_anova(data, by_type)
            all_results[name] = anova_results
            
            # Extract significant genes (adjusted p-value < 0.05)
            sig_genes = anova_results[anova_results['adjusted_p_value'] < 0.05].copy()
            sig_genes['gene_set'] = name
            significant_genes = pd.concat([significant_genes, sig_genes])
            
            # Save individual results
            output_file = Path(output_dir) / f"{name}_anova_{by_type}.csv"
            anova_results.to_csv(output_file)
            
        # Save combined significant results
        significant_genes.sort_values('adjusted_p_value').to_csv(
            Path(output_dir) / f"all_significant_genes_anova_{by_type}.csv"
        )
        
        # Generate summary statistics
        summary_stats = pd.DataFrame({
            'gene_set': list(all_results.keys()),
            'total_genes': [len(df) for df in all_results.values()],
            'significant_genes': [len(df[df['adjusted_p_value'] < 0.05]) for df in all_results.values()],
            'min_p_value': [df['p_value'].min() for df in all_results.values()],
            'median_p_value': [df['p_value'].median() for df in all_results.values()]
        })
        
        # Save summary statistics
        summary_stats.to_csv(Path(output_dir) / f"anova_summary_stats_{by_type}.csv", index=False)
        
        return all_results
    def generate_heatmap(self, splits: Dict[str, pd.DataFrame], title: str, 
                        figsize: Tuple[int, int] = (12, 8),
                        cmap: str = "RdBu_r",
                        save_path: Optional[str] = "./output") -> None:
        """
        Generate and save a heatmap visualization of gene expression patterns.
        
        Args:
            splits: Dictionary of split DataFrames by subtype
            title: Title for the heatmap
            figsize: Figure size as (width, height)
            cmap: Color map for the heatmap
            save_path: Directory to save the generated figure
        """
        # Create figure
        plt.figure(figsize=figsize)
        
        # Combine all data and add subtype labels
        combined_data = pd.DataFrame()
        subtype_labels = []
        
        for subtype, data in splits.items():
            combined_data = pd.concat([combined_data, data], axis=1)
            subtype_labels.extend([subtype] * data.shape[1])
            
        # Create color mapping for subtypes
        unique_subtypes = list(splits.keys())
        colors = sns.color_palette("husl", len(unique_subtypes))
        subtype_colors = {subtype: color for subtype, color in zip(unique_subtypes, colors)}
        color_row = [subtype_colors[label] for label in subtype_labels]
        
        # Generate heatmap
        g = sns.clustermap(combined_data,
                          cmap=cmap,
                          col_colors=[color_row],
                          xticklabels=False,
                          yticklabels=True,
                          figsize=figsize,
                          row_cluster=True,
                          col_cluster=False,
                          z_score=0, 
                          vmin=-1.5,
                          vmax=1.5)  # Standardize rows
        
        # Add title
        g.figure.suptitle(title, y=1.02)
        
        # Add legend
        legend_elements = [plt.Rectangle((0,0),1,1, facecolor=color, label=subtype)
                         for subtype, color in subtype_colors.items()]
        g.ax_heatmap.legend(handles=legend_elements,
                           title="Subtypes",
                           bbox_to_anchor=(1.3, 1),
                           loc='upper right')
        
        # Save figure if path provided
        if save_path:
            Path(save_path).mkdir(parents=True, exist_ok=True)
            plt.savefig(f"{save_path}/{title.replace(' ', '_')}.png",
                       bbox_inches='tight',
                       dpi=300)
            
        plt.close()
    
    

def main():
    paths = DataPaths()
    analysis = TCGA_GeneExpressionData()
    
    # Load data
    analysis.load_clinical_data(paths)
    analysis.load_gene_sets(paths)
    
    # Run ANOVA analysis for all gene sets
    clinical_results = analysis.analyze_all_genesets_anova(by_type="clinical")
    molecular_results = analysis.analyze_all_genesets_anova(by_type="molecular")
    
    # Process each gene set individually as before
    for name, data in analysis.gene_sets.items():
        logger.info(f"Processing {name} gene set")
        
        # Split data by clinical types
        splits = analysis.split_by_types(data, "clinical")
        
        # Run t-test analysis
        ttest_results = analysis.run_ttest_analysis(splits)
        
        # Generate visualization
        analysis.generate_heatmap(splits, f"{name} Expression Patterns")
        
        # Save t-test results
        ttest_results.to_csv(f"./Processed_data_output/{name}_ttest.csv")

if __name__ == "__main__":
    main()
