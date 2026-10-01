# Bovine 3D endometrial gland–stroma model: analysis code and data

Code and processed data for:

> Ribes Martinez E\*, Viana AES\*, Ferronato GA, Nociti RP, Oshiro TSI, De Bem THC, Meirelles FV, Ferraz MAMM.
> *A biomimetic 3D bovine endometrial gland–stroma model responds to ovarian hormones and embryonic signals.*
> (\* equal contribution)

A single script, `run_all.R`, reproduces every statistical result and data figure in the paper from the files in `data/`.

## How to run

1. Install R (≥ 4.3) and, optionally, RStudio.
2. Open this folder in RStudio (or `setwd()` to it).
3. Run `source("run_all.R")`.

The first run installs any missing CRAN and Bioconductor packages. This needs internet access and takes about 10–20 minutes. A full run then takes a few minutes. Results are written to `output/`:

| Output | Contents |
|---|---|
| `output/key_results.txt` | Every number quoted in the manuscript (sample sizes, effect estimates, 95% CIs, p-values) |
| `output/figures/Figure1.*`, `Supplementary_Figure1_*` | EndoECM matrisome composition |
| `output/figures/Figure3.*`, `Figure4.*`, `Figure5.*` | Assembled data figures (PDF, PNG, TIFF 600 dpi) |
| `output/figures/SuppFig2_Reactome_hormone.*`, `SuppFig3_Reactome_embryo.*`, `Supplementary_Figure4_RNAseq_all_libraries.*` | Supplementary figures |
| `output/figures/panels/` | Every panel as a separate file |
| `output/tables/Supplementary_Data_1_*.xlsx`, `Supplementary_Data_2_*.xlsx` | Differential-expression results (DEGs, all genes, sensitivity analysis, normalised counts) |
| `output/tables/*.csv` | Model contrasts, design tables, GO enrichment tables, library QC |
| `output/sessionInfo.txt` | R and package versions used |

## Repository structure

```
run_all.R                  entry point
R/00_setup.R               packages, paths, plotting defaults
R/01_rnaseq.R              DESeq2 differential expression, PCA, MA, volcano, heatmaps (Fig. 4b–e, 5c–f; Suppl. Fig. 4; Suppl. Data 1–2)
R/02_enrichment.R          Reactome (Suppl. Fig. 2–3) and GO biological-process over-representation
R/03_nanoindentation.R     Young's modulus and viscoelasticity (Fig. 3b–e)
R/04_gland_area.R          gland area, hormone and embryo experiments (Fig. 4a, 5a)
R/05_embryo.R              blastocyst hatching (Fig. 5b)
R/06_assemble_figures.R    multi-panel figures
R/07_proteomics.R          EndoECM matrisome re-analysis (Fig. 1, Suppl. Fig. 1)
data/                      input data (see below)
```

## Data

| File | Description | Unit of observation |
|---|---|---|
| `data/rnaseq/counts_gene_level.csv` | Gene-level read counts (featureCounts, uniquely mapped reads in exons; STAR alignment to the Ensembl *Bos taurus* genome) for all nine libraries, including the two excluded low-depth libraries | gene × library |
| `data/rnaseq/sample_metadata.csv` | Library → experiment (Exp1–3), well and condition. B2 = hormone-treated co-culture + embryos; B3 = hormone-treated co-culture; B4 = untreated co-culture | library |
| `data/nanoindentation/young_modulus.csv` | Young's modulus (YM) and effective YM per indentation site; `Batch` = independent preparation, `Replicate`/`Drop` identify the hydrogel; `Hormonal_treatment = yTissue` = native endometrium | indentation site |
| `data/nanoindentation/viscoelasticity.csv` | Storage (E1 = E′) and loss (E2 = E″) moduli at 1, 5, 10 and 20 Hz; consecutive rows of the same hydrogel, module and frequency are successive indentation sites | site × frequency |
| `data/gland_area/gland_area.csv` | Area of individual glands (µm²) per imaging day; `Replicate` = independent co-culture experiment, `Well_group` = hydrogel | gland |
| `data/embryo/hatching.csv` | Hatched / allocated Day 7 blastocysts per group and IVF replicate | group × IVF replicate |
| `data/nanoindentation/native_tissue_viscoelasticity.csv` | E′/E″ of native bovine endometrium (6 cows, follicular and luteal phase; from Ferraz et al., ref. 15), plotted for comparison in Fig. 3c–e | site × frequency |
| `data/proteomics/endoECM_proteomics_maxLFQ.tsv` | DIA-MS MaxLFQ protein intensities of 10 EndoECM hydrogel samples (5 follicular EF*, 5 luteal EL*; ref. 28) | protein × sample |
| `data/proteomics/Bt_Matrisome.csv` | *Bos taurus* matrisome gene list (division, category) | gene |
| `data/figure_assets/Fig3a_timeline.png` | Schematic used as Figure 3a | — |

Raw sequencing reads (FASTQ) are no longer available from the sequencing provider. The gene-level count matrix is the most primary RNA-seq data that remain. Figure 1 and Supplementary Figure 1 are a secondary re-analysis of previously published EndoECM proteomic data (ref. 28; ProteomeXchange/PRIDE PXD053248).

## Analysis decisions

- **RNA-seq.**
  - Genes whose identifier starts with `LOC` or `NON` were removed. Genes with fewer than 10 total counts across the libraries of a contrast were also removed.
  - DESeq2 was run with design `~ condition`, the Wald test and Benjamini–Hochberg adjustment. DEGs are genes with padj < 0.05 and |log2FC| > 1.
  - Libraries with fewer than 10 million gene-assigned counts were excluded: Hydrog.Exp1.B3 (2.7 M) and Hydrog.Exp1.B4 (5.0 M). All retained libraries have 23.8–35.4 M.
  - The Hormone vs No hormone comparison uses B3 vs B4 from Exp2 and Exp3.
  - Embryos were cultured on hormone-stimulated constructs, so the Embryo vs No embryo comparison uses B2 (Exp1–3) vs B3 (Exp2–3).
  - As a sensitivity analysis, both contrasts are repeated with `~ experiment + condition`.
  - PCA uses the 500 most variable genes after the variance-stabilising transformation.
- **Pathway enrichment.** Reactome over-representation (ReactomePA; bovine symbols mapped to human orthologues by symbol) and GO biological process (clusterProfiler, org.Bt.eg.db), run separately for genes with higher and lower expression; universe = all genes tested; BH-adjusted p < 0.05.
- **Nanoindentation.**
  - Indentation sites are sub-samples of hydrogels, and hydrogels are nested within independent preparations.
  - Linear mixed-effects models on log(YM) and log(tan δ) have preparation and hydrogel as random intercepts. Contrasts are reported as ratios with 95% CIs, Holm-adjusted.
  - tan δ = E″/E′ per site, averaged over the four frequencies.
- **Gland area.**
  - Individual glands are sub-samples of hydrogels (wells), and hydrogels are nested within experiments.
  - Linear mixed-effects models on log(area) have experiment and hydrogel as random intercepts.
  - The treatment × day interaction is tested by likelihood-ratio test. Day-specific ratios are Holm-adjusted.
- **Proteomics.** Proteins detected in ≥ 6 of the 10 samples are retained (missing values set to 0) and classified with the *Bos taurus* matrisome list; relative abundance is the share of mean MaxLFQ intensity within the matrisome or within each category.
- **Hatching.**
  - A binomial GLMM uses IVF replicate as a random intercept. If that variance is not estimable, replicate enters as a fixed blocking factor.
  - Results are reported as pairwise odds ratios with Wald 95% CIs, Holm-adjusted.

## License

Code: MIT. Data: CC BY 4.0.
