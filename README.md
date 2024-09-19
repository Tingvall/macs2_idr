#MACS2_IDR
**Nextflow based pipeline for macs2 peak calling and IDR**

## Installation

####  Create conda environment:
Download the macs2_idr_env.yml file and create a conda environment that contain all packages required to run the pieline:
```bash
conda env create -f macs2_idr_env.yml
```
Activate conda environment:
```bash
conda activate macs2_idr_env
```

####  Launch pipeline:
Download the pipeline (including macs2_idr.nf & nextflow.config). 
```bash
nextflow run macs2_idr.nf --help
```
