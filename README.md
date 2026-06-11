# TinProject — Longitudinal Brain Connectomics in Chronic Tinnitus

Analysis pipeline for tracking large-scale brain network reorganization during real-time fMRI neurofeedback treatment for chronic tinnitus.

## Background

This project investigates how functional brain connectivity changes across longitudinal resting-state and neurofeedback fMRI sessions in patients with chronic tinnitus. The pipeline covers functional connectomics, graph-theoretical network metrics, and brain fingerprinting.

Based on the NeuroTin trial reported in:
> Gninenko et al. (2024). Functional MRI Neurofeedback Outperforms Cognitive Behavioral Therapy for Reducing Tinnitus Distress. *Radiology*, 310(2), e231143.

## Repository structure

```
TinProject/
├── Code/Tessnim/
│   ├── notebooks/          ← numbered analysis notebooks
│   └── utils/              ← shared Python modules and atlas files
├── Results/Tessnim/
│   ├── all_subjects_yeo7_metrics.csv
│   ├── all_subjects_yeo17_metrics.csv
│   ├── figures/            ← exported 
```


## Notebooks

## Analysis workflow

The notebooks follow the progression of the project from data loading to final statistical analysis.

### 1. Data loading and exploratory analyses

- **1. TinLoading.ipynb**  
  Loads Schaefer parcel-wise fMRI time series, computes Fisher z-transformed Pearson functional connectivity (FC) matrices, and aggregates them into Yeo-7 network-level FC matrices.

- **2. TinDiving.ipynb**  
  Exploratory analysis of FC matrices, including partial-correlation (inverse covariance) connectivity and comparison with Pearson correlation.

### 2. Resting-state connectivity analyses

- **3. NF_initial_analysis.ipynb**  
  First implementation of parcel-level segregation and integration metrics with network-level visualizations.

- **4. rsfMRI.ipynb**  
  Computes segregation, integration, and normalized segregation across resting-state runs and sessions.

### 3. Neurofeedback metrics pipeline

- **6. Neurofeedback_ALL.ipynb**  
  Core neurofeedback analysis pipeline. Computes segregation, integration, normalized segregation, and participation coefficient across all subjects and runs.

- **6a. Data building.ipynb**  
  Builds the main analysis dataset, adds negative participation coefficient, removes Run 7, and exports final CSV files.

### 4. Visualization

- **6b. Visualization.ipynb**  
  Session-wise and subject-wise trajectories of network metrics.

- **6g. Visualization violin.ipynb**  
  Violin plots of network metrics across sessions.

- **7a. Yeo17 Visualization.ipynb**  
  Visualization of Yeo-17 network metrics.

### 5. Statistical analyses

- **6c. Statistical testing.ipynb**  
  Spearman trend analyses, FDR correction, and subject clustering.

- **10. Sampling.ipynb**  
  Feedback vs. Transfer comparisons for Yeo-7 networks using KDE visualization and Mann–Whitney tests.

- **11. Sampling yeo 17.ipynb**  
  Same analyses applied to Yeo-17 networks.

- **12. Linear model and trend lines.ipynb**  
  Linear mixed-effects models for longitudinal changes.

- **14. Control Covariates.ipynb**  
  Sensitivity analyses controlling for global FC and condition effects.

### 6. Parcellation and robustness analyses

- **6d. Parcellation influence.ipynb**  
  Evaluates robustness across Schaefer resolutions (100–800 parcels).

- **6e. The Run 7 dilemma.ipynb**  
  Investigates the influence of the anomalous Run 7.

### 7. Yeo-17 data generation

- **7b. Yeo17 800 Data building.ipynb**  
  Computes Yeo-17 metrics on the Schaefer-800 atlas.

- **7c. Yeo17 200 Data building.ipynb**  
  Computes Yeo-17 metrics on the Schaefer-200 atlas.

### 8. Connectome fingerprinting

- **8. Fingerprinting.ipynb**  
  Implements connectome fingerprinting, differential identifiability, PCA reconstruction, and network-specific identification analyses.

### 9. Participation coefficient analyses

- **15. Participation coefficient clean.ipynb**  
  Computes positive and negative participation coefficients and aggregates results at parcel and network levels.
## Utils

| File | Purpose |
|------|---------|
| `metrics.py` | Graph-metric functions (segregation, integration, normalized segregation) |
| `metrics_yeo17.py` | Same metrics recomputed natively with Yeo-17 network labels |
| `plotting.py` | Shared plotting helpers |
| `matched_labels_exact.csv` | Schaefer-200 → Yeo-17 label mapping |
| `matched_labels_exact_800.csv` | Schaefer-800 → Yeo-17 label mapping |

## Methods overview

**Parcellations:** Schaefer atlas (r100–r800); network assignments via Yeo-7 and Yeo-17.

**Functional connectivity:** Pearson correlation (default), partial correlation, mutual information.

**Network metrics:**
- `segregation` — mean within-network FC
- `integration` — mean between-network FC
- `participation_coefficient` — parcel-level cross-network integration

**Statistics:** Spearman correlation for trends, Wilcoxon signed-rank for paired comparisons, FDR correction (Benjamini-Hochberg), effect sizes as rank-biserial correlation.

## Environment

```
Python 3.x
numpy, pandas, scipy, nilearn, matplotlib, seaborn
```
