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
    
    def run_one_way_anova(self, data: pd.DataFrame, groups: pd.Series):
        """
        Perform one-way ANOVA for each gene across groups
        
        Args:
            data: Gene expression matrix (genes × samples)
            groups: Series containing group labels
        """
        results = []
        for gene in data.index:
            # Create lists of expression values for each group
            group_values = [data.loc[gene, groups == group] for group in groups.unique()]
            
            # Run ANOVA
            f_stat, p_val = stats.f_oneway(*group_values)
            
            # Add to results
            results.append({
                'gene': gene,
                'F_statistic': f_stat,
                'p_value': p_val
            })
        
        return pd.DataFrame(results).set_index('gene')
   
    
    def run_ttest_analysis(self, data_split: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Perform t-test analysis between TNBC and HR positive samples"""
        tnbc_data = data_split['TNBC'].T
        hr_data = data_split['HR_positive'].T
        
        results = []
        for gene in tnbc_data.columns:
            t_stat, p_val = stats.ttest_ind(tnbc_data[gene], hr_data[gene])
            results.append({
                'gene': gene,
                'tValue': t_stat,
                'pValue': p_val,
                'mean_of_TNBC': tnbc_data[gene].mean(),
                'mean_of_HR_positive': hr_data[gene].mean()
            })
        
        return pd.DataFrame(results).set_index('gene')
    
    def generate_heatmap(self, data_split: Dict[str, pd.DataFrame], title: str) -> None:
        """Generate heatmap visualization"""
        plt.figure(figsize=(15, 10))
        
        # Combine data with clinical types
        combined_data = pd.concat([df for df in data_split.values()], axis=1)
        clinical_types = pd.Series('other', index=combined_data.columns)
        for ctype, df in data_split.items():
            clinical_types[df.columns] = ctype
            
        # Create heatmap
        sns.clustermap(
            combined_data,
            col_colors=pd.get_dummies(clinical_types).reindex(combined_data.columns),
            cmap='RdBu_r',
            center=0,
            figsize=(15, 10),
            xticklabels=False
        )
        plt.title(title)
        plt.show()

def main():
    paths = DataPaths()
    analysis = GeneExpressionData()
    
    # Load data
    analysis.load_clinical_data(paths)
    analysis.load_gene_sets(paths)
    
    # Process each gene set
    for name, data in analysis.gene_sets.items():
        logger.info(f"Processing {name} gene set")
        
        # Split data by clinical types
        splits = analysis.split_by_clinical_types(data)
        
        # Run t-test analysis
        ttest_results = analysis.run_ttest_analysis(splits)
        
        # Generate visualization
        analysis.generate_heatmap(splits, f"{name} Expression Patterns")
        
        # Save results
        ttest_results.to_csv(f"./Processed_data_output/{name}_ttest.csv")

if __name__ == "__main__":
    main()
