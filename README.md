# CapMux

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)


**CapMux**: a Snakemake pipeline for early demultiplexing of split-pool scRNA-seq data into sample-resolved outputs


Single-cell RNA sequencing methods based on split-pool combinatorial barcoding enable high-throughput profiling, but sample identity is often introduced during early barcoding steps rather than through the final Illumina library index. Consequently, reads from multiple biological samples remain pooled until relatively late stages of data processing, complicating per-sample analysis and selective extraction of samples of interest. Here, I present CapMux, a Snakemake-based pipeline for processing split-pool scRNA-seq data, from raw sequencing files to sample-resolved outputs. CapMux supports workflows starting from either Illumina BCL files or pre-generated FASTQ files and reconstructs sample identity by integrating sub-library index information with experiment-specific barcoding plate layout. The pipeline was developed for the CapSeq method but is configurable for related scRNA-seq combinatorial barcoding designs through specification of barcode positions, linker structure, and experimental layout. CapMux resolves pooled data into outputs for each sample, enabling independent quality-control summaries, mapping statistics, count matrices, and downstream visualizations. Runtime benchmarking indicated that secondary demultiplexing added only a modest computational overhead. Together, these results show that CapMux provides a practical and adaptable framework for recovering sample-level resolution from split-pool scRNA-seq data.



## Features
- Converts index-centric sequencing data into sample-centric libraries
- Supports custom split-pool or combinatorial indexing designs (up to 3 barcode segments and 2 UMIs in Read 1)
- Optional demultiplexing step based on BC1
- Works from Illumina raw output folder or pre-existing FASTQs
- Reproducible execution via Snakemake
- STARsolo mapping


## Installation and usage

Install snakemake via mamba and activate environment
```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
```

Clone the repository into the desired run folder
```bash
git clone https://github.com/Boyoron/capmux.git run
```

Perform the run of the workflow using snakemake
```bash
cd run
snakemake --cores 16 --use-conda
```

Edit `config.yaml` before running to match your experiment.



## Configuration (`config.yaml`)

All pipeline behavior is controlled through `config.yaml`.  
A detailed description with comments available in `config_comments.yaml`.

Below is a short description of each parameter.

---

### Inputs / Outputs

| Parameter | Allowed Values | Description |
|--------|------|---------|
| `run_dir` | string | Path to Illumina folder containing BCL files. Used only when `run_mode = "bcl"`. |
| `out_dir` | string | Output root directory where the pipeline results will be stored. |
| `sample_sheet` | string | Extended sample sheet (`.xlsx`) used by the workflow to describe barcode‑to‑sample mapping. |
| `tso_seq` | string | Template-switch oligo (TSO) sequence used in reverse transcription reaction protocol (used for trimming). |

An example `extended_sample_sheet.xlsx` is provided in the repository:
```
assets/sample_sheets/example_extended_sample_sheet.xlsx
```

---

### Demultiplexing Selection

| Parameter | Allowed Values | Description |
|--------|----------------|---------|
| `demux_by` | `index`, `bc1` | `index`: Illumina demux only. `bc1`: additional BC1 demultiplexing step. |
| `run_mode` | `bcl`, `fastq` | `bcl`: start from Illumina BCL files. `fastq`: start from existing FASTQs. |
| `fastq_dir` | string | Folder containing FASTQ files (used only when `run_mode = fastq`). |
| `demux_index` | `i7`, `i5`, `i7_i5` | Index reads used for Illumina demultiplexing. |
| `use_bases_mask` | string | Custom `bcl2fastq` mask. Leave "" for default behaviour. |


### Barcode / UMI Parsing

Controls how barcode segments and UMIs are extracted from Read 1 (R1).

| Parameter | Allowed Values | Description |
|--------|------|---------|
| `barcode.allow_mismatches` | int | Allowed mismatches when matching BC1 to whitelist. Used only when `demux_by = "bc1"`.|
| `barcode.segments` | dict | 1‑based coordinates on R1 defining barcode segments and UMI positions. |
| `barcode.whitelists` | dict | Paths to whitelist files for each barcode segment. |

**Segment Rules**

- Coordinates are **1‑based**
- Disable a segment using `start: 0, len: 0`
- Supports up to 3 BC segments + 2 UMIs
- index whitelist is optional. Leave "" if unused or all indexes are provided in extended sample sheet (.xlsx)

Example:
```yaml
barcode:
  allow_mismatches: 1
  segments:
    bc1:  {start: 32, len: 8}
    bc2:  {start: 18, len: 10}
    bc3:  {start: 5,  len: 8}
    umi1: {start: 1,  len: 4}
    umi2: {start: 40, len: 4}
  whitelists:
    bc1: "assets/barcodes/bc1_list.txt"
    bc2: "assets/barcodes/bc2_list.txt"
    bc3: "assets/barcodes/bc3_list.txt"
    index: "" 
```

---

### STARsolo Parameters

| Parameter | Allowed Values | Description |
|--------|----------------|---------|
| `star.features` | `Gene`, `GeneFull` | Feature type to quantify. |
| `star.cb_match` | `Exact`, `1MM`, `1MM_multi`, `EditDist_2` | Cell barcode matching strategy. |
| `star.index` | string | Path to STAR genome index directory. |

---

### Resources

| Parameter | Allowed Values | Description |
|--------|------|---------|
| `threads` | int | Default number of CPU threads. |
| `memory` | string | Default RAM allocation (e.g., `64G`). |

---



## Pipeline Output Structure

The results folder of the pipeline is structured as follows:

```
results
├── demuxed/                  # FASTQ files after demultiplexing
├── logs/                     # Logs from all steps
├── mapped/                   # STARsolo mapping outputs
├── merged_fastq/             # bc read assembly from R1 (e.g. demux_by = bc1, BC1-BC2-BC3-UMI) or from R1 and I1 (demux_by = bc1, BC1-BC2-BC3-BC4-UMI)
├── qc/                       # FastQC and MultiQC reports
├── trimmed_fastq/            # FASTQ files trimmed from TSO sequence 
├── bc1_to_sample.csv         # Sample sheet generated from 'bc1' sheet (demux_by = bc1)
├── sample_sheet.csv          # Sample sheet generated to run Illumina bcl2fastq demultiplexing
├── demux_root.txt            # File contains the name of demux root folder
├── sample_ids_bc1.txt        # List of sample IDs after bc1 demultiplexing (demux_by = bc1)
└── sample_ids.txt            # List of sample IDs after i7 (bcl2fastq) demultiplexing
```


## Version

**v1.0.0** – initial release  


## Citation

If you use CapMux, please cite: 

> Denis Baronas et al., **High-throughput single cell omics using semipermeable capsules.** *Science* DOI:https://doi.org/10.1126/science.ady7227 (2025)


## License
[![Licence](https://img.shields.io/github/license/Ileriayo/markdown-badges?style=for-the-badge)](./LICENSE)

## Contact

denis.baronas [at] gmc [dot] vu [dot] lt