# Ensemble SDM with Spatial Block Cross-Validation — Prosopis cineraria (Pakistan)

Spatially-validated ensemble species distribution modeling comparing BIOCLIM, GLM, and Random Forest for Khejri (Prosopis cineraria) in Pakistan.

## Key Results
- **4-fold spatial block cross-validation** (longitude-based blocks)
- **GLM**: AUC = 0.774 ± 0.001 (most consistent single algorithm)
- **Ensemble**: AUC = 0.808 ± 0.016 (highest accuracy, lowest variance)
- **Random Forest**: AUC = 0.739 ± 0.033 (high variance across blocks)
- **BIOCLIM**: AUC = 0.738 ± 0.063 (highest variance)

## Files
- `report.pdf` — Full manuscript
- `fig_01_spatial_cv_boxplot.png` — Cross-validation performance
- `fig_01_model_agreement.png` — Model agreement map
- `fig_02_uncertainty.png` — Between-algorithm disagreement
- `fig_03_response_curves.png` — Ecological response curves
- `fig_05_variable_importance.png` — Random Forest variable importance
- `spatial_cv_results.csv` — Performance statistics

## Citation
Naveed, A. (2026). Ensemble Climate-Suitability Modeling for Prosopis cineraria in Pakistan: Spatial Block Cross-Validation and Conservation Priorities. Zenodo. [https://doi.org/10.5281/zenodo.22247414]

## Related Work
- Khejri SDM (BIOCLIM): https://doi.org/10.5281/zenodo.22174696
- Spatial Bias Analysis: https://doi.org/10.5281/zenodo.22140857
- Data Quality Audit: https://doi.org/10.5281/zenodo.22199329
