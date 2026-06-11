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

| Notebook | Purpose |
|----------|---------|
| `1. TinLoading.ipynb` | Loads Schaefer parcel-wise resting-state fMRI time series, computes Fisher z-transformed Pearson FC matrices, and aggregates them into Yeo-7 network-to-network connectivity matrices |
| `2. TinDiving.ipynb` | Explores FC matrix properties using partial correlation (inverse covariance), comparing Pearson and partial correlation approaches on Yeo-7 networks |
| `3. NF_initial_analysis.ipynb` | Computes parcel-level segregation and integration across 20 subjects with Fisher z-transform, generates session-averaged visualizations per Yeo-7 network, and exports run-level metrics to CSV |
| `4. rsfMRI.ipynb` | Loads pre-computed FC matrices and computes per-network segregation, integration, and normalized segregation across run/session pairs with group and individual visualizations |
| `6. Neurofeedback_ALL.ipynb` | Loads 20 subjects' neurofeedback FC data, computes segregation/integration/normalized-segregation and participation coefficient per network per run, and exports parcel-level metrics to CSV |
| `6a. Data building.ipynb` | Extends notebook 6 by adding negative participation coefficient (anticorrelated edges), removes Run 7, and exports augmented CSV across 21 subjects and 15 sessions |
| `6b. Visualization.ipynb` | Plots per-subject per-session segregation/integration/normalized-segregation trajectories across runs with group average and per-subject lines, with condition shading (Feedback/Transfer/NoFeedback) |
| `6c. Statistical testing.ipynb` | Tests monotonic session trends via Spearman correlation, applies FDR correction, and clusters subjects by trend pattern using a |ρ| ≥ 0.2 threshold |
| `6d. Parcellation influence.ipynb` | Compares session trends in segregation/integration/normalized-segregation across four Schaefer resolutions (r100–r800) using LME, demonstrating robustness to parcellation choice |
| `6e. The Run 7 dilemma.ipynb` | Investigates the Run 7 anomaly by comparing raw vs. mixed-effects corrected vs. sqrt(n)-weighted data to assess its impact on longitudinal trends |
| `6g. Visualization violin.ipynb` | Generates violin plots for segregation/integration/normalized-segregation per Yeo-7 network across sessions (runs averaged within session) with session means and subject distributions |
| `7a. Yeo17 Visualization.ipynb` | Loads and visualizes segregation/integration/normalized-segregation across 17 Yeo networks for 20 subjects across 15 sessions with per-network subplot grids |
| `7b. Yeo17 800 Data building.ipynb` | Computes Yeo-17 metrics on 800-parcel Schaefer atlas using nilearn atlas labels matched via CSV lookup, exporting segregation/integration/normalized-segregation and parcel-level PC to CSV |
| `7c. Yeo17 200 Data building.ipynb` | Computes Yeo-17 metrics on 200-parcel atlas via `matched_labels_exact.csv`, exporting the same structure as the Yeo-7 CSV but with 17 networks |
| `8. Fingerprinting.ipynb` | Implements Finn et al. connectome fingerprinting: Pearson FC similarity across sessions, differential identifiability (I_diff), PCA reconstruction, and per-network identification accuracy |
| `10. Sampling.ipynb` | Compares Feedback vs. Transfer distributions using bootstrapped KDE overlays and Mann-Whitney U tests (FDR-corrected) on segregation/integration/normalized-segregation per Yeo-7 network |
| `11. Sampling yeo 17.ipynb` | Same as notebook 10 but applied to Yeo-17 networks: bootstrapped KDE + Mann-Whitney U tests on segregation/integration comparing Feedback and Transfer conditions |
| `12. Linear model and trend lines.ipynb` | Fits linear mixed models (session as fixed effect, subject random intercept) on network metrics, with and without global mean FC control, reporting FDR-corrected β and p-values per network |
| `14. Control Covariates.ipynb` | Tests whether session trends remain significant after controlling for global mean FC and condition effects (Transfer/NoFeedback vs. Feedback) using LME |
| `15. Participation coefficient clean.ipynb` | Computes positive and negative participation coefficient per parcel (hub score based on cross-network edge distribution), aggregates to network level, and exports results for Yeo-7 and Yeo-17 |

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
