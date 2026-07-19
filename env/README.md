# Environment Setup

This directory provides three ways to provision the `snp_call_nf` software environment.

## Contents

| File | Purpose |
|------|---------|
| `snp_call_nf.yaml` | Conda environment spec (reference / direct use) |
| `pixi.toml` | Pixi project manifest (used by Docker build) |
| `Dockerfile` | Multi-stage Docker build with pixi |
| `build.sh` | One-command builder for Docker image + Apptainer SIF |

## Included Tools

| Tool | Version |
|------|---------|
| bowtie2 | 2.4.4 |
| GATK4 | 4.2.2.0 |
| bedtools | 2.30.0 |
| samtools | 1.13 |
| bcftools | 1.13 |
| requests | latest |

---

## Option 1 — Conda (native)

Use the YAML file directly with Conda / Mamba / Micromamba:

```bash
# Create the environment
conda env create -f env/snp_call_nf.yaml

# Activate it
conda activate snp_call_nf
```

Or with Mamba (faster):

```bash
mamba env create -f env/snp_call_nf.yaml
mamba activate snp_call_nf
```

---

## Option 2 — Docker

### Prerequisites

- Docker 20.10+

### Build

```bash
# From the project root:
./env/build.sh docker
```

This builds a minimal image (~952 MB) using a two-stage process:
1. **Builder stage** — installs packages via [pixi](https://prefix.dev/pixi), then strips headers, man pages, static libs, and other cruft.
2. **Runtime stage** — copies only the cleaned environment into `debian:bookworm-slim`.

### Run

```bash
# Interactive shell
docker run --rm -it snp_call_nf:latest

# Run a specific tool
docker run --rm snp_call_nf:latest gatk --version
docker run --rm snp_call_nf:latest samtools --version

# Mount data and run a pipeline step
docker run --rm \
    -v /path/to/data:/data \
    snp_call_nf:latest \
    bowtie2 -x /data/ref -U /data/reads.fq -S /data/out.sam
```

### Custom image name / tag

```bash
IMAGE_NAME=my_registry/snp_call_nf IMAGE_TAG=v1.0 ./env/build.sh docker
```

---

## Option 3 — Apptainer (SIF)

NOTE: If you are a IGS server, you can skip the following and use a copy of the SIF image located at

`/local/projects-t3/CVD/public_data/snp_call_nf_RESOURCES/snp_call_nf.sif`

### Prerequisites

- Apptainer >=1.3.6
- An existing Docker image (`snp_call_nf:latest` — see Option 2)

### Build

```bash
# From the project root:
./env/build.sh apptainer
```

This converts the local Docker image into a single-file SIF container (~491 MB).

To build both Docker + SIF in one go:

```bash
./env/build.sh          # same as: ./env/build.sh all
```

### Run

```bash
# Interactive shell
apptainer shell snp_call_nf.sif

# Run a specific tool
apptainer exec snp_call_nf.sif gatk --version
apptainer exec snp_call_nf.sif samtools --version

# Mount data and run a pipeline step
apptainer exec --bind /path/to/data:/data \
    snp_call_nf.sif \
    bowtie2 -x /data/ref -U /data/reads.fq -S /data/out.sam
```

### Using with Nextflow

```sh
# under the project root directory
mkdir -p tmpdir
nextflow main.nf \
    -process.container=$PWD/snp_call_nf.sif \
    -without-conda \
    -with-apptainer \
    -process.containerOptions="-B $(realpath .)/tmpdir -B $(realpath .)/ref"

# some $PWD does not work but $(realpath .) does, likely nextflow resolves paths differently
# 
# Alternatively add a profile, see nextflow.config:
# apptainer {
#     apptainer {
#         enabled    = true
#         autoMounts = true
#         runOptions = "-B ${params.gatk_tmpdir} -B ${projectDir}/ref"
#     }
#     process {
#         container = "${projectDir}/snp_call_nf.sif"
#     }
# }
# 
# then run 

mkdir -p tmpdir
nextflow main.nf -without-conda  --with-apptainer  -profile apptainer,standard 
# or 
# nextflow main.nf -without-conda  --with-apptainer  -profile apptainer,slurm 
```


## Updating Dependencies

1. Edit the version pins in `pixi.toml` (and optionally `snp_call_nf.yaml`).
2. Rebuild:

   ```bash
   ./env/build.sh
   ```

3. If using Conda directly, recreate the environment:

   ```bash
   conda env remove -n snp_call_nf
   conda env create -f env/snp_call_nf.yaml
   ```
