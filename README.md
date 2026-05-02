# Sex Differences in Gene Expression Across Human Tissues
### GTEx v8 Transcriptomic Analysis in R — MSc Health Data Science

**Tosan Akpituren** · MSc Health Data Science for Applied Precision Medicine · University of Dundee · 2024  
Supervisor: Dr. Andrew Brown

---

## Overview

This project investigates sex-driven differences in gene expression across three clinically significant human tissues — **aorta artery**, **breast mammary tissue**, and **brain hippocampus** — using publicly available GTEx v8 transcriptomic data. The central question is whether genes influenced by biological sex differences are also the genes contributing to observed disparities in the prevalence, susceptibility, and severity of atherosclerosis, breast cancer, and Alzheimer's disease.

**[View Interactive Portfolio →](https://tosan431.github.io/gtex-sex-differential-expression-r)**

---

## Key Statistics

| Metric | Value |
|---|---|
| Total data points processed | ~21 million |
| Gene variables analysed per tissue | ~19,000 |
| Biological samples | 1,088 (432 aorta · 197 brain · 459 breast) |
| Significant DEGs — Aorta Artery | 801 (Benjamini-Hochberg adjusted) |
| Significant DEGs — Breast Mammary | 5,597 (Bonferroni adjusted) |
| Significant KEGG pathways — Breast | 5 (BH-adjusted, DAVID) |

---

## Tissues & Associated Diseases

| Tissue | Disease of Interest | Samples |
|---|---|---|
| Aorta Artery | Atherosclerosis | 432 |
| Brain Hippocampus | Alzheimer's Disease | 197 |
| Breast Mammary | Breast Cancer | 459 |

---

## Analytical Pipeline

```
1. Data Acquisition     → GTEx v8 TPM expression matrices + phenotype metadata
2. Harmonisation & QC   → Subject ID validation, whitespace removal, gene filtering
3. Differential Expr.   → Gene-wise linear regression (lm) across ~19,000 genes per tissue
4. Multiple Testing     → Benjamini-Hochberg (aorta/brain) · Bonferroni (breast)
5. Annotation & GWAS    → biomaRt Ensembl mapping · Fisher's Exact Test vs GWAS Catalogue
6. Pathway Enrichment   → KEGG via clusterProfiler (aorta) and DAVID (breast)
```

---

## Key Findings

**Aorta Artery**
- 801 significant DEGs identified after Benjamini-Hochberg correction
- clusterProfiler KEGG enrichment identified significant female-upregulated pathways including Ribosome, Coronavirus Disease–COVID-19, and Dopaminergic Synapse
- GWAS orthogonal analysis: odds ratio 1.12 (p = 0.683, non-significant) vs atherosclerosis genes

**Breast Mammary Tissue**
- 5,597 significant DEGs after Bonferroni correction (13,245 under BH — too large for pathway analysis)
- Five statistically significant KEGG pathways (BH-adjusted), all enriched in **female-upregulated genes**:
  1. ECM-Receptor Interaction
  2. Cell Adhesion Molecules
  3. Viral Protein Interaction with Cytokine and Cytokine Receptor
  4. Coronavirus Disease – COVID-19
  5. Basal Cell Carcinoma
- GWAS orthogonal analysis: odds ratio 1.15 (p = 0.084, non-significant) vs breast cancer genes

**Brain Hippocampus**
- Only 13 significant DEGs after BH correction — insufficient for pathway analysis
- Likely reflects smaller sample size (n=197) and male-skewed cohort (73% male)
- Non-parametric Wilcoxon test applied as robustness check; results consistent

---

## Tools & Packages

**Language:** R v4.4.1 (RStudio)

| Package | Purpose |
|---|---|
| `ggplot2` | Volcano plots, Manhattan plots, histograms, boxplots, pathway charts |
| `biomaRt` | Ensembl ID → gene symbol mapping, genomic coordinates |
| `clusterProfiler` | KEGG pathway enrichment analysis |
| `dplyr` / `tidyverse` | Data manipulation and pipeline management |
| `limma` | Differential expression framework |
| `DAVID` (web tool) | KEGG pathway enrichment — breast tissue |
| `GWAS Catalogue` | Orthogonal disease gene validation |

**Statistical methods:** Gene-wise linear regression (`lm`), Wilcoxon rank-sum test, Benjamini-Hochberg FDR correction, Bonferroni correction, Fisher's Exact Test

---

## Repository Structure

```
gtex-sex-differential-expression-r/
│
├── index.html          # Interactive portfolio (live via GitHub Pages)
├── README.md           # This file
```

---

## Data Source

All expression data used in this project are from the **GTEx v8** public resource:

> The GTEx Consortium. *The GTEx Consortium atlas of genetic regulatory effects across human tissues.* Science 369, 1318–1330 (2020).

Data are publicly available at [gtexportal.org](https://gtexportal.org). All samples are pseudonymised; no personally identifiable information is contained in this repository.

---

## Author

**Tosan Akpituren**  
MSc Health Data Science for Applied Precision Medicine  
University of Dundee, 2024

---

*Submitted in partial fulfilment of the requirements for the degree of Master of Science in Health Data Science for Applied Precision Medicine, University of Dundee, August 2024.*
