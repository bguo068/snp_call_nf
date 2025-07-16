# snp_call_nf 

The pipeline is designed to perform joint variant calling on large Plasmodium
Whole Genome Sequencing (WGS) datasets. It follows `GATK` best practices and the
[MalariaGEN Pf6 data-generating
methods]((https://ngs.sanger.ac.uk//production/malaria/pfcommunityproject/Pf6/Pf_6_extended_methods.pdf)).
Briefly, raw reads are first mapped to the human GRCh38 reference genome to
remove host reads, with the remaining reads being mapped to the Pf3D7 reference
genome (PlasmoDB_44). The mapped reads are then processed using `GATK`'s
`MarkDuplicates` and `BaseRecalibrator` tools. Following this, analysis-ready
mapped reads for each isolate are used to generate per-sample calls
(`HaplotypeCaller` /GVCF mode). These per-sample calls are combined and run
through a joint-call step (`GenotypeGVCFs`) to obtain unfiltered multi-sample
VCFs. A machine learning-based variant filtration strategy (VQSR via `GATK`'s
`VariantRecalibrator`), or/and hard-filteration strategy, can then be used to
retain high-quality variants.

The `main` branch of the repository is employed for *Plasmodium falciparum*.
However, the pipeline is expected to work with other *Plasmodium* species,
provided the corresponding configuration and reference files are given (See
nextflow.config). A separate `vivax` branch with configuration and reference
files for *Plasmodium vivax* under development and can be checked with the
[link](https://github.com/bguo068/snp_call_nf/tree/vivax).


# Software environtment

The pipeline has been tested on MacOS and Linux system. The software
dependencies (including the version numbers of used software) are defined within the
`env/nf.yaml` Conda recipe (for Nextflow engine) and the `env/snp_call_nf.yaml`
Conda recipe (for the pipeline itself). Installtion instruction is listed below.
The estimated time of instation is about 5-20 minutes.

1. Change to a working folder that is large enough to store the snp call result
files. Git clone the pipeline and change directory to the pipeline folder
```sh
cd YOUR_WORKING_DIR # replace `YOUR_WORKING_DIR` with your real path
git clone https://github.com/bguo068/snp_call_nf.git
cd snp_call_nf
```
2. Install conda from [here](https://docs.conda.io/en/latest/miniconda.html) if you have not
3. Install the `nf` and the `snp_call_nf` conda environments
```sh
# NOTE: The nextflow engine and the pipeline may need different version of java.
# We use two different Conda environments to address the conflict.
# install nextflow
conda env create -f env/nf.yaml
# install snp_call_nf
conda env create -f env/snp_call_nf.yaml
```
# How to run the pipeline?

1. Link the reference files (internal users) or prepare them by yourself
(**external** users)
- Link the ref files on IGS server
```
ln -s /local/projects-t3/toconnor_grp/bing.guo/ref/* ref/
```
or
```
ln -s /local/projects-t2/CVD/Takala-Harrison/Cambodia_Bing/ref/* ref/
```

- Link the ref files on Rosalind, the reference file can be linked by running
```
ln -s /local/data/Malaria/Projects/Takala-Harrison/Cambodia_Bing/ref/* ref/
```
- Prepare ref file by yourself:
```sh
cd ref
conda activate snp_call_nf
python3 prep_ref_files.py
conda deactivate
cd ..
```

2. Run the pipeline
    - Test it on HPC (local): `conda activate nf; nextflow main.nf`
        - This will use a tiny **test dataset** from `test_data` folder
        - To test on a larger **test data**, please run `cd ena_data;
        ena_data/download_ena_data.sh` (it may take hours to download all the
        real data from ENA). Once downloaded, change directory to project folder
        (where `main.nf` is located) by run `cd ..`, and the pipeline can be run
        with `conda activate nf; nextflow main.nf --fq_map
        ena_data/ena_fastq_map.tsv`
    - Test it on SGE server: `conda activate nf; nextflow main.nf -profile sge`
        - This will use a small dataset from `test_data` folder
    - You will need to edit `fastq_map.tsv` file to include the raw
    reads(`fastq.gz` files) of your own samples.

## Optional arguments
3. Split chromosomes to better parallelize joint call:
    - by default, the genome is split by chromosomes
    - you can specify cmd line option `--split intervals` to split the chromosome into more 
    intervals.

4. Enable `vqsr` variant filtering. By default, `vqsr` is not enabled. To enable
this option, you can specify `--vqsr true` to the nextflow command line

There are options to run parts of the pipeline:

- If '--parasite_reads_only' is specified, the pipeline will stop after host
reads are removed and the remaining reads are saved to FQ files.
- If '--coverage_only' is specified, the pipeline will stop after completing
the necessary steps to produce the reads coverage, before the GATK_APPLY_BQSR
process.
- If '--gvcf_only' is specified, the pipeline will stop after generating
per-sample gVCF files.

# Important input and output files

1. Main input file is `./fastq_map.tsv`
    - Five columns delimited by tab: string, interger, string, interger, string
    - `HostId` is the index of host genomes from 0, see `params.host` in nextflow.config file
    - `MateId` can be 0 for single-end sequencing, or 1 and 2 for pair-end sequencing
2. Main configureation file is `./nextflow.config`
    - For SEG users, be sure to edit sge config about `clusterOptions = "-P toconnor-lab -cwd -V"` to reflect your lab specifc sge qsub option
3. Main pipeline script is `./main.nf`
4. Main output files/folders:
    - `result/read_length` folder: report the raw read length for each samples/runs
    - `result/flagstat_host` and  `result/flagstat_parasite`: flagstat of
    aligned reads (aligned with host genome and parasite genome respectively)
    - `result/recalibrated`: analysis ready bam files
    - `result/coverage`: read converage based on the analysis-ready bam files 
    - `result/flagstat`: bam flatstat based on the analysis-ready bam files 
    - `result/gvcf`: single-sample vcf files
    - `result/hardfilt`: multiple-sample (joint-call) vcf files with hard filterating annotations
    - `result/vqsrfilt`: multiple-sample (joint-call) vcf files with vqsr-based filterating annotations.
   You can decide to use one of these, `result/hardfilt` and `result/vqsrfilt`.

# Workflow chart

![flowchar](./flowchart.png)

# Citations
This pipeline was originally developed for the `posseleff` project. 
If you find this pipeline useful, please consider citing our preprint:
> Guo, B., Borda, V., Laboulaye, R., Spring, M. D., Wojnarski, M., Vesely, B. A., Silva, J. C.,
> Waters, N. C., O'Connor, T. D., & Takala-Harrison, S. (2023). Strong Positive Selection Biases
> Identity-By-Descent-Based Inferences of Recent Demography and Population Structure in
> Plasmodium falciparum. bioRxiv : the preprint server for biology, 2023.07.14.549114.
> https://doi.org/10.1101/2023.07.14.549114
