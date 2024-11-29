# Palbociclib Analysis Pipeline in Python

## Overview
This repository contains a reproducible pipeline for analyzing CRISPR screen data to identify significant gene hits. The pipeline is implemented in **Python** and includes steps for data preprocessing, statistical analysis, and result visualization.

## Features
- Filters and normalizes gRNA count data.
- Performs fold-change calculations and statistical tests (e.g., t-tests).
- Aggregates gRNA results at the gene level.
- Generates visualizations such as volcano plots and heatmaps.

## Requirements
To run this pipeline, ensure you have the following software and packages installed:

### Software
- **Python** (Version ≥ 3.8)
- Jupyter Notebook (optional but recommended)

### Python Packages
- `pandas`
- `numpy`
- `scipy`
- `matplotlib`
- `seaborn`

You can install these packages using the following command:
```bash
pip install pandas numpy scipy matplotlib seaborn
