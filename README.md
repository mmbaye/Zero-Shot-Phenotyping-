# **Foundation AI-Based**  High-Throughput Phenotyping of Sorghum Architectural Traits

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-4.3%2B-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.10%2B-green.svg)](https://python.org)
[![Platform](https://img.shields.io/badge/Platform-Google%20Colab-orange.svg)](https://colab.research.google.com)

> **Paper:** *Foudation AI-Based High-Throughput Phenotyping Reveals Genotypic Variation in Sorghum Panicle and Leaf Architecture Across Eight Varieties*
> Submitted to **Plant Phenomics** (Elsevier/AAAS)

---

## Overview

This repository contains all code and analysis scripts for a pipeline combining the **Segment Anything Model (SAM)** with **consumer-grade LiDAR (iPad Pro M4)** for automated phenotyping of sorghum architectural traits from non-standardized field photographs.

### Key results
- **8 sorghum genotypes** (Senegal field trial)
- **Welch's ANOVA**: F(7, 8.56) = 36.29, p = 1.07 × 10⁻⁵ for panicle length
- **SAM vs LiDAR validation**: r = 0.87, R² = 0.76, p = 0.005
- **LiDAR vs manual validation**: R² = 0.9992, MAE = 0.42 cm
- **Novel quality metric**: Completeness ratio C = A / (π × L/2 × l/2) × 100

---

## Requirements

### Python (Google Colab  GPU T4 recommended)
```
segment-anything
torch torchvision
opencv-python
scikit-learn
matplotlib numpy pandas
```

Install SAM model checkpoint:
```python
# In Colab
!wget https://dl.fbaipublicfiles.com/segment_anything/sam_vit_b_01ec64.pth
```

### R (≥ 4.3)
```r
install.packages(c(
  "tidyverse", "ggplot2", "patchwork", "scales",
  "janitor", "readxl",
  "emmeans", "multcomp", "multcompView",
  "rstatix", "car",
  "FactoMineR", "factoextra",
  "ggrepel"
))
```
---

## Calibration

```python
# In SAM_pipeline_V3.ipynb
TAILLE_CARREAU_CM = 15.   
```
The tile dimension was verified using iPad Pro LiDAR, correcting an initial assumption of 15 cm (correction factor: ×1.632 on all L, l measurements; ×2.663 on areas).
---

## Pipeline overview

```
Field photo
    │
    ▼
Hough calibration (tile → px/cm)
    │
    ▼
Background removal (HSV + CLAHE)
    │
    ▼
SAM segmentation (ViT-B, point prompt)
    │
    ▼
Ellipse fitting (PCA on mask pixels)
    │  → L (length), l (width), A (area), ∠ (angle)
    ▼
Completeness ratio C = A / (π × L/2 × l/2) × 100
    │  → COMPLETE ≥70% | PARTIAL 20–70% | ERROR <20%
    ▼
CSV export
```

---

## Running the analysis

### Step 1 Process images (Colab)
1. Open `python/SAM_pipeline_V3.ipynb` in Google Colab
2. Upload your images
3. Set `TAILLE_CARREAU_CM` to your tile size
4. Run all cells → outputs one CSV per image

### Step 2  Merge outputs (R)
### Step 3 Statistical analysis (R)
---

## Data

The complete dataset is available at Zenodo: **[DOI: to be assigned upon acceptance]**

| File | Description | n |
|------|-------------|---|
| `Genotype_SAM_FINAL.csv` | SAM traits per image | 49 images × 8 genotypes |
| `Panicle_lidar.xlsx` | LiDAR panicle lengths | 22 measurements |
| `yield_data.csv` | Agronomic performance | 48 rows (4 reps × 12 genotypes) |
| `dataset_SAM_LIDAR_Yield_merged.csv` | Merged analysis dataset | 8 genotypes |

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{[author]_2026_sorghum_SAM,
  title   = {Foundation AI-Based High-Throughput Phenotyping Reveals Genotypic Variation
             in Sorghum Panicle and Leaf Architecture Across Eight Varieties},
  author  = {[Authors]},
  journal = {Plant Phenomics},
  year    = {2026},
  doi     = {[to be assigned]}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

## Corresponding author email  
 modou.mbaye@isra.sn 
