# Per-Vessel Myocardial Blood Flow Improvement After CABG Quantified by CT Myocardial Perfusion Imaging

This repository contains the data and analysis scripts associated with the following publication:

**Seresti A**, Castillo E, Fahmy S, Menon K, Ahmed AH, Aung PP, Malik SB, Nguyen PKP, Boyd J, Marsden AL, Khan MO, Nieman K.  
**Per-vessel myocardial blood flow improvement after coronary artery bypass graft surgery quantified by CT myocardial perfusion imaging.**  
*Journal of Cardiovascular Computed Tomography*, 2026.  
https://doi.org/10.1016/j.jcct.2026.05.009

---

## Overview

This study uses dynamic CT myocardial perfusion imaging (CT-MPI) to quantify regional (per-vessel) changes in myocardial blood flow (MBF) before and after coronary artery bypass grafting (CABG). A key methodological contribution is a **median-based normalization** approach to compensate for inter-scan variability in absolute MBF values, enabling reliable pre- vs post-CABG comparisons.

The repository includes:

- Anonymized per-territory MBF data (pre- and post-CABG)
- Python scripts for data processing, normalization, statistical analysis, and figure generation
- Code to reproduce the main results and figures from the paper

---

## Repository Structure


├── data/                  
├── scripts/             
│   ├── dumps/
│   ├── ImageTools/
│   ├── MBFTools/
│   └── SimVascularTools/
├── main.ipynb.         # main notebook for data analysis             
├── figures/               
├── requirements.txt
└── README.md

---

## Requirements

- Python ≥ 3.10
- Core libraries: `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, `statsmodels`

Install dependencies with:

```bash
pip install -r requirements.txt


@article{Seresti2026,
  title   = {Per-vessel myocardial blood flow improvement after coronary artery bypass graft surgery quantified by CT myocardial perfusion imaging},
  author  = {Seresti, Anahita and Castillo, Edgard and Fahmy, Seraphim and Menon, Karthik and Ahmed, Abdrabu Hamid and Aung, Phyo Phyo and Malik, Sachin Basiq and Nguyen, Patricia Kim Phuong and Boyd, Jack and Marsden, Alison L. and Khan, M. Owais and Nieman, Koen},
  journal = {Journal of Cardiovascular Computed Tomography},
  year    = {2026},
  doi     = {10.1016/j.jcct.2026.05.009}
}