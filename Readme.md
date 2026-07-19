# snp_call_nf

Nextflow pipeline for joint variant calling on *Plasmodium* Whole Genome Sequencing (WGS) datasets.

---

## Table of Contents

1. [Overview](#overview)
2. [Software Environment & Installation](#software-environment--installation)
3. [Quick Start: Test the Pipeline](#quick-start-test-the-pipeline)
4. [Running on Your Own Data](#running-on-your-own-data)
5. [Pipeline Options](#pipeline-options)
6. [Input & Output Files](#input--output-files)
7. [Debugging](#debugging)
8. [Workflow Chart](#workflow-chart)
9. [Citations](#citations)

---

## Overview

The pipeline follows [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels) and the [MalariaGEN Pf6 data-generating methods](https://ngs.sanger.ac.uk//production/malaria/pfcommunityproject/Pf6/Pf_6_extended_methods.pdf).

**Pipeline steps:**

1. **Host read removal** — Raw reads are mapped to the human GRCh38 reference genome; aligning reads are discarded.
2. **Parasite alignment** — Remaining reads are mapped to the *P. falciparum* 3D7 reference (PlasmoDB_44).
3. **BAM processing** — `MarkDuplicates` and `BaseRecalibrator` (GATK) produce analysis-ready BAMs.
4. **Per-sample calling** — `HaplotypeCaller` in GVCF mode generates per-isolate variant calls.
5. **Joint calling** — Per-sample GVCFs are combined via `GenotypeGVCFs` into a multi-sample VCF.
6. **Variant filtration** — Hard-filtering or VQSR (`VariantRecalibrator`) retains high-quality variants.

### Species support

The `main` branch targets **_Plasmodium falciparum_**. The pipeline is expected to work with other *Plasmodium* species when the appropriate configuration and reference files are provided (see `nextflow.config`).

A `vivax` branch for **_Plasmodium vivax_** is under development: [vivax branch](https://github.com/bguo068/snp_call_nf/tree/vivax).

---

## Software Environment & Installation

**Supported platforms:** macOS (Intel) and Linux.  
**Estimated install time:** 5–20 minutes.

The Nextflow engine and the pipeline may require different Java versions. Two separate Conda environments are used to avoid conflicts:
- `env/nf.yaml` — Nextflow engine
- `env/snp_call_nf.yaml` — Pipeline tools

### Step 1: Clone the repository

```sh
cd YOUR_WORKING_DIR          # replace with your actual path
git clone https://github.com/bguo068/snp_call_nf.git
cd snp_call_nf
```

### Step 2: Install Conda (if needed)

Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

### Step 3: Create the Conda environments

```sh
# Install Nextflow environment
conda env create -f env/nf.yaml

# Install pipeline environment
conda env create -f env/snp_call_nf.yaml
```

> **Note:** If you already have Nextflow installed, you can skip the `nf` environment and pin the version directly:
> ```sh
> nextflow self-update
> export NXF_VER=26.04.6
> nextflow main.nf [args]
> ```

> **Apptainer/Singularity users:** See [`env/README.md`](env/README.md) for container-based setup instructions.

---

## Quick Start: Test the Pipeline

A tiny test dataset is included in the `test_data/` folder. Use it to verify your installation works.

### Local / HPC (single node)

```sh
conda activate nf
nextflow main.nf
```

To test with VQSR enabled:

```sh
conda activate nf
nextflow main.nf --vqsr true -profile standard,test
```

### SGE cluster

```sh
conda activate nf
nextflow main.nf -profile sge
```

### Slurm cluster

```sh
conda activate nf
nextflow main.nf -profile slurm
```

### Larger test dataset (ENA downloads)

To test with real (but public) ENA data:

```sh
cd ena_data
bash download_ena_data.sh      # may take hours to download
cd ..
conda activate nf
nextflow main.nf --fq_map ena_data/ena_fastq_map.tsv
```

---

## Running on Your Own Data

### 1. Prepare reference files

**Internal users (IGS / Rosalind):** Symlink existing reference files.

| Server    | Command                                                                  |
|-----------|--------------------------------------------------------------------------|
| IGS       | `ln -s /local/projects-t3/toconnor_grp/bing.guo/ref/* ref/`             |
| IGS (alt) | `ln -s /local/projects-t2/CVD/Takala-Harrison/Cambodia_Bing/ref/* ref/` |
| Rosalind  | `ln -s /local/data/Malaria/Projects/Takala-Harrison/Cambodia_Bing/ref/* ref/` |

**External users:** Generate reference files from scratch.

```sh
cd ref
conda activate snp_call_nf
python3 prep_ref_files.py
conda deactivate
cd ..
```

### 2. Prepare your sample sheet

Edit `fastq_map.tsv` (or create a new one) with your sample information.  
See [fastq_map.tsv](fastq_map.tsv) for an example.  
The file is **tab-delimited** with five columns:

| Column   | Type    | Description                                                          |
|----------|---------|----------------------------------------------------------------------|
| Sample   | string  | Unique sample identifier                                             |
| HostId   | integer | Host genome index (0-based); see `params.host` in `nextflow.config`  |
| Run      | string  | Run / replicate identifier (e.g. `r1`)                               |
| MateId   | integer | `0` for single-end; `1` or `2` for paired-end                        |
| Fastq    | string  | Path to the FASTQ file                                               |

### 3. Run the pipeline

```sh
conda activate nf
nextflow main.nf --fq_map your_fastq_map.tsv -profile slurm
```

Replace `-profile slurm` with `-profile sge` or omit it for local execution.

> **Before running on a cluster**, edit `nextflow.config` to match your lab's settings:
> - **SGE:** eg., update `clusterOptions = "-P toconnor-lab -cwd -V"`
> - **Slurm:** eg., update `clusterOptions = "--account=cvd"`

---

## Pipeline Options

### Variant filtration

| Flag          | Default | Description                                        |
|---------------|---------|----------------------------------------------------|
| `--vqsr true` | off     | Enable VQSR-based variant filtration               |

If VQSR fails due to convergence issues, adjust `--vqsr_opts` (see `nextflow.config`).

### Parallelization

Control how the genome is divided for parallel joint calling via `--split`. The available intervals are defined in `nextflow.config` under `params.genome_intervals`.

| Value          | Behavior                                                                 |
|----------------|--------------------------------------------------------------------------|
| `chromosomes`  | **(default)** One interval per chromosome (14 intervals for *P. falciparum*) |
| `intervals`    | Sub-chromosomal intervals — splits larger chromosomes for finer-grained parallelism (44 intervals) |

Example:

```sh
nextflow main.nf --split intervals
```

### Subset / partial runs

Use these to stop the pipeline at a specific stage:

| Flag                     | Stops after…                                                |
|--------------------------|-------------------------------------------------------------|
| `--parasite_reads_only`  | Host read removal — unmapped reads saved to FASTQ           |
| `--coverage_only`        | Coverage computation — before `GATK_APPLY_BQSR`             |
| `--gvcf_only`            | Per-sample GVCF generation — before joint calling            |

### Host read filtering method

| Flag                   | Behavior                                                                                              |
|------------------------|-------------------------------------------------------------------------------------------------------|
| _(default)_            | Reads are aligned to human genome first; aligning reads are removed before parasite mapping.          |
| `--use_concat_genome`  | Reads are aligned to a concatenated human+parasite genome. Reads mapping to parasite chromosomes are retained and realigned; human-mapping reads are discarded.  |

---

## Input & Output Files

### Key inputs

| File                | Purpose                                             |
|---------------------|-----------------------------------------------------|
| `fastq_map.tsv`     | Sample sheet (see [Running on Your Own Data](#1-prepare-your-sample-sheet)) |
| `nextflow.config`   | Pipeline configuration (profiles, params, resources) |
| `main.nf`           | Main Nextflow workflow script                       |

### Key outputs

| Path                              | Description                                              |
|-----------------------------------|----------------------------------------------------------|
| `result/readlen_raw/`             | Raw read lengths per sample/run                           |
| `result/flagstat_raw/`            | Flagstat for host-genome alignments                       |
| `result/flagstat_parasite/`       | Flagstat for parasite-genome alignments                   |
| `result/recalibrated/`            | Analysis-ready BAM files (post host-removal & recalibration) |
| `result/recal_bam_coverage/`      | Read coverage from analysis-ready BAMs                    |
| `result/recal_bam_flagstat/`      | Flagstat from analysis-ready BAMs                         |
| `result/gvcf/`                    | Per-sample GVCF files                                     |
| `result/hardfilt_vcf/`            | Joint-called VCF with hard-filter annotations             |
| `result/vqsrfilt_vcf/`            | Joint-called VCF with VQSR annotations (only when `--vqsr true` is used) |

You can use either `result/hardfilt_vcf/` or `result/vqsrfilt_vcf/` as your final filter annotated call set.

---

## Debugging

Run a stub workflow to quickly validate workflow logic without executing real jobs:

```sh
conda activate nf
NXF_VER=26.04.6 nextflow main.nf -stub --vqsr true -without-conda
```

---

## Workflow Chart

![flowchart](./flowchart.png)

Note: the VQSR/use concatenated genome features are not shown in the workflow chart.

---

## Citations

This pipeline was originally developed for the `posseleff` project. If you find it useful, please cite:

> Guo, B., Borda, V., Laboulaye, R. et al. Strong positive selection biases identity-by-descent-based inferences of recent demography and population structure in *Plasmodium falciparum*. *Nat Commun* **15**, 2499 (2024).  
> <https://doi.org/10.1038/s41467-024-46659-0>
