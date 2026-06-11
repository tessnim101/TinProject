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

## Analysis notebooks

```
1. TinLoading.ipynb                 Data loading and FC matrix construction
2. TinDiving.ipynb                  Exploratory FC analyses (Pearson vs partial correlation)

3. NF_initial_analysis.ipynb        Initial segregation/integration analyses
4. rsfMRI.ipynb                     Resting-state network metrics

6. Neurofeedback_ALL.ipynb          Main neurofeedback analysis pipeline
6a. Data building.ipynb             Final dataset generation (+ negative PC)
6b. Visualization.ipynb             Longitudinal visualizations
6c. Statistical testing.ipynb       Trend analyses and clustering
6d. Parcellation influence.ipynb    Robustness across Schaefer resolutions
6e. The Run 7 dilemma.ipynb         Run 7 quality-control investigation
6g. Visualization violin.ipynb      Distribution visualizations

7a. Yeo17 Visualization.ipynb       Yeo-17 visualizations
7b. Yeo17 800 Data building.ipynb   Yeo-17 metrics (Schaefer-800)
7c. Yeo17 200 Data building.ipynb   Yeo-17 metrics (Schaefer-200)

8. Fingerprinting.ipynb             Connectome fingerprinting analyses

10. Sampling.ipynb                  Feedback vs Transfer (Yeo-7)
11. Sampling yeo 17.ipynb           Feedback vs Transfer (Yeo-17)

12. Linear model and trend lines.ipynb   Linear mixed-effects models
14. Control Covariates.ipynb             Covariate-control analyses

15. Participation coefficient clean.ipynb
    Positive and negative participation coefficient analyses
```

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
