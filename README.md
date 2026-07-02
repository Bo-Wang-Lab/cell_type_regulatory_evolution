# Flexible use of conserved motif vocabularies constrains genome access in cell type evolution

Code and analysis notebooks for the study described in [Chai et al., bioRxiv 2024.09.03.611027](https://www.biorxiv.org/content/10.1101/2024.09.03.611027v1).

## Overview

This repository contains the Python, R, and bash scripts used to process multiome (scRNA + scATAC) data across flatworms (*Macrostomum lignano*, *Schmidtea mediterranea*, *Schistosoma mansoni*), zebrafish, and mouse, and to train ChromBPNet models of chromatin accessibility per cell type. Notebooks are organised by manuscript figure.

## Repository layout

```
cell_type_regulatory_evolution/
├── environment/                # conda environment files
│   ├── base_python.yml
│   └── base_R.yml
├── notebooks/Figures/          # analysis notebooks, one folder per figure
│   ├── Fig1/                   # single-cell atlases and marker-gene / ArchR analyses per species
│   ├── Fig2/                   # cross-species SAMap comparisons + sub-clustering
│   ├── Fig3/                   # ChromBPNet motif deviations and conservation
│   ├── Fig4/                   # motif co-occurrence, jaccard, scVI, gene phylogenies
│   ├── Fig5/                   # cross-species motif quantification (Smed–Sman)
│   └── Fig6/                   # cross-species motif work in mouse and zebrafish
└── step-by-step-scripts/       # preprocessing pipeline commands and configs
    ├── STARsolo_alignment_parameters       # scRNA read alignment
    ├── chromap_alignment_parameters         # scATAC read alignment
    ├── kallisto_alignment_parameters        # alternative RNA quantification
    ├── peakcalling_consensus_step_by_step.txt  # ENCODE atac-seq-pipeline + iterative overlap
    ├── iterative_overlap_peaks.py           # consensus-peak builder
    └── ChromBPNet_model_step_by_step.txt    # end-to-end ChromBPNet training recipe
```

## Environments

Two conda environments cover most notebooks:

```bash
conda env create -f environment/base_python.yml   # base_python
conda env create -f environment/base_R.yml        # base_R
```

Two figures use dedicated Docker containers rather than the conda envs:

| Figure | Environment | Container |
|--------|-------------|-----------|
| Fig 1 | `base_python` / `base_R` | — |
| **Fig 2** | **SAMap Docker** | [`jessedgibson/samap-test`](https://hub.docker.com/r/jessedgibson/samap-test) |
| Fig 3 | `base_python` / `base_R` | — |
| Fig 4 | `base_python` / `base_R` | — |
| **Fig 5** | **ChromBPNet Docker** | [`kundajelab/chrombpnet`](https://hub.docker.com/r/kundajelab/chrombpnet) |
| **Fig 6** | **ChromBPNet Docker** | [`kundajelab/chrombpnet`](https://hub.docker.com/r/kundajelab/chrombpnet) |

All ChromBPNet model training (see `step-by-step-scripts/ChromBPNet_model_step_by_step.txt`) was performed inside the ChromBPNet Docker container — the same image linked above.

**Pull the containers:**
```bash
docker pull jessedgibson/samap-test
docker pull kundajelab/chrombpnet
```

## Preprocessing pipeline

Before running any of the figure notebooks, raw reads are processed through the following steps (see `step-by-step-scripts/`):

1. **Read alignment.** scRNA with STARsolo (`STARsolo_alignment_parameters`) or alternatively kallisto (`kallisto_alignment_parameters`); scATAC with chromap (`chromap_alignment_parameters`).
2. **Peak calling and consensus set.** ENCODE atac-seq-pipeline per cell type, then merged into an iterative-overlap consensus with `iterative_overlap_peaks.py`. See `peakcalling_consensus_step_by_step.txt` for the end-to-end recipe.
3. **ChromBPNet training.** Per-species / per-cell-type bias and accessibility models, following `ChromBPNet_model_step_by_step.txt` inside the ChromBPNet Docker container.

## Data availability

- **Raw sequencing reads:** SRA BioProject [PRJNA1153627](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1153627).
- **Pre-trained ChromBPNet models:** Figshare [https://doi.org/10.6084/m9.figshare.32861690](https://doi.org/10.6084/m9.figshare.32861690).
- **Processed h5ad files:** available from the authors upon request.
