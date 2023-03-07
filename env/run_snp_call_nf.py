#! /usr/bin/env python3

import sys
from subprocess import run
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(
    "A wrapper script to faciliate running `snp_call_nf` pipeline"
)
parser.add_argument(
    "--fqmap",
    type=str,
    default=None,
    help="""path to the tsv table of input fastq files. The first row is ignored; 
    for the following lines, each should have 5 fields.
       col1: sample name of the fastq file;
       col2: an index of host species, use 0 if samples are isolated from human;
       col3: the run id, useful when the samples were sequnced in different runs;
       col4: the mate id, should be 1 or 2 for pair-end sequencing, otherwise use 0;
       col5: the path to the fastq.gz file.""",
)
parser.add_argument(
    "--test", action="store_true", help="when set, run pipeline on a small test sets"
)
parser.add_argument(
    "--ena", action="store_true", help="when set, run pipeline on 6 real WGS samples"
)


def check_fastq_map(fn):
    with open(fn, "r") as f:
        for line in f:
            fields = line.rstrip().split("\t")
            for e in fields:
                assert e != ""
            assert len(fields) == 5
            sample = fields[0]
            hostid = fields[1]
            runid = fields[2]
            mateid = fields[3]
            path = fields[4]


def lnk_testdata_fastq_map():
    dir = "/local/data/Malaria/Projects/Takala-Harrison/AFRIMS/snp_call_nf/"
    fn = f"{dir}/test_data"
    if Path("test_data").exists():
        Path("test_data").unlink()
    Path("test_data").symlink_to(fn, target_is_directory=True)
    if Path("fastq_map.tsv").exists():
        Path("fastq_map.tsv").unlink()
    Path("fastq_map.tsv").symlink_to(f"{dir}/fastq_map.tsv")

    return "fastq_map.tsv"


def lnk_ena_fastq_map():
    dir = "/local/data/Malaria/Projects/Takala-Harrison/AFRIMS/snp_call_nf/"
    fn = f"{dir}/ena_data"
        
    if Path("ena_data").exists():
        Path("ena_data").unlink()
    Path("ena_data").symlink_to(fn, target_is_directory=True)
    if Path("fastq_map.tsv").exists():
        Path("fastq_map.tsv").unlink()
    Path("fastq_map.tsv").symlink_to(f"{dir}/ena_data/ena_fastq_map.tsv")

    return "fastq_map.tsv"


def run_nextflow(fn):
    root = "/local/data/Malaria/Projects/Takala-Harrison/AFRIMS"
    cmd = f"""
    eval "$({root}/miniconda3/bin/conda shell.bash hook)"
    conda activate nf
    nextflow {root}/snp_call_nf/main.nf --fq_map {fn} -resume \
            --gatk_tmpdir /local/data/Malaria/Projects/Takala-Harrison/AFRIMS/tmp
    """
    print(cmd)
    assert run(cmd, shell=True).returncode == 0


if __name__ == "__main__":
    args = parser.parse_args()
    fn = args.fqmap
    
    if fn is None and not args.test and not args.ena:
        print(
            "if fqmap not specified, either test or ena flag should be set", file=sys.stderr
        )
        sys.exit(-1)
    elif fn is None and args.test:
        fn = lnk_testdata_fastq_map()
        run_nextflow(fn)
    elif fn is None and args.ena:
        fn = lnk_ena_fastq_map()
        run_nextflow(fn)
    else:
        run_nextflow(fn)
