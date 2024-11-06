# Flexible use of conserved motif vocabularies constrains genome access in cell type evolution

A combination of Python, R, and bash scripts used in the analysis described here - https://www.biorxiv.org/content/10.1101/2024.09.03.611027v1

Complete descriptions of the packages installed in the Python and R environments used in this study (suitable for installing with `conda env create -f`) can be found in the `base_python.yml`, and `base_R.yml` files, respectively. Any analyses performed using [SAM](https://github.com/atarashansky/self-assembling-manifold)/[SAMap](https://github.com/atarashansky/SAMap) were run instead within the Docker container located [here](https://hub.docker.com/r/jessedgibson/samap-test), and analyses using [ChromBPNet](https://github.com/kundajelab/chrombpnet) were performed using the Docker container located [here](https://hub.docker.com/r/kundajelab/chrombpnet)

Notebooks provided can reproduce all data panels in the manuscript when run on the data available through SRA (PRJNA1153627), which will also be available as both raw fastq, processed h5ad files, and pre-trained ChromBPNet model h5 files through Dryad upon publication.
