# MACS2_IDR
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


## Usage

#### Input

#### Running the pipeline

There are two version of the pipeline depending on the mode of peak calling in Macs2 (narrow/broad). Typical commands for running the pipelines are as follows:

Narrow mode
```bash
nextflow run macs2_idr.nf --samples sample_info.txt --outdir peaks --genome_size 2652783500 --macs_q 0.05 --idr_threshold 0.05
```

Broad mode
```bash
nextflow run macs2_idr.nf --samples sample_info.txt --outdir peaks --genome_size 2652783500 --macs_q 0.05 --idr_threshold 0.05
```

#### Arguments
| Argument | Description |
| --- | --- |
| `--samples` | Path to text file specifying inputs. (For example see: [peak_info.txt](example_files/peak_info.txt)).|
| `--genome_size` | Size of the genome. Default 2652783500 (mm10).|
| `--macs_q` | q-value treshold for Macs2 peakcalling. Default: 0.05.|
| `--skip_idr` | Specify if IDR should be skipped. Default: false.|
| `--idr_threshold` | Treshold for IDR. Default: 0.05.|
| `--cores` | Number of cores to use. Deatult: 8|
| `--help` | Display help message.|


## Output
All outputs are placed in the direcory specified by `--outdir`. Depending on the options, a number of different subdirectories will be created within this directory:
- `<outdir>/`
  - `peaks/`: Peak files for each replicate.
  - `idr/`: Output from IDR analysis, including optimal peak set and QC.
  - `bigwigs/`: Bigwig files for each replicates and pooled samples for visualization in genome browser.
