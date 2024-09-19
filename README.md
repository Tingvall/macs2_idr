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
The typical command for running the pipeline is as follows:
```bash
nextflow run macs2_idr.nf --samples sample_info.txt --outdir peaks
```

#### Arguments
| Argument | Description |
| --- | --- |
| `--samples` | Path to text file specifying inputs. (For example see: [peak_into.txt](example_files/peak_info.txt))|
| `--genome_size` | |
| `--macs_q` | |
| `--idr_threshold` | |
| `--skip_idr` | |
| `--cores` | |
| `--help` | |


## Output
