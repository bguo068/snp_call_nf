# How to run the pipeline?

1. Install conda from [here](https://docs.conda.io/en/latest/miniconda.html)
2. Install the `nf` and the `snp_call_nf` conda environments:
```sh
# install nextflow
conda env create -f env/nf.ymal
# install snp_call_nf
conda env create -f env/snp_call_nf.yaml
```
3. Change to a working folder that is large enough to store the snp call result
files. Git clone the pipeline and change directory to the pipeline folder
```sh
cd YOUR_WORKING_DIR # replace `YOUR_WORKING_DIR` with your real path
git clone git@github.com:gbinux/snp_call_nf.git
cd snp_call_nf
```
4. Link the reference files or prepare them by yourself

- Link the ref files on IGS server
```
ln -s /local/projects-t3/toconnor_grp/bing.guo/ref/* ref/
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

5. Run the pipeline
    - Test it on HPC (local): `conda activate snp_call_nf; nextflow main.nf`
    - Test it on SGE server: `conda activate snp_call_nf; nextflow main.nf -profile sge`

6. Split chromosomes to better parallelize joint call:
    - by default, the genome is split by chromosomes
    - you can specify cmd line option `--split intervals` to split the chromosome into more 
    intervals.

# Important files

1. Main input file is `./fastq_map.tsv`
    - Five columns: string, interger, string, interger, string
    - `HostId` is the index of host genomes from 0, see `params.host` in nextflow.config file
    - `MateId` can be 0 for single-end sequencing, or 1 and 2 for pair-end sequencing
2. Main configureation file is `./nextflow.config`
    - Be sure to edit sge config about `clusterOptions = "-P toconnor-lab -cwd -V"` to reflect your lab specifc sge qsub option
3. Main pipeline script is `./main.nf`

# Workflow chart

![flowchar](./flowchart.png)
