#! /usr/bin/env python
import ftplib
import sys
from pathlib import Path
from subprocess import run
from typing import ClassVar

import requests


class URLS:
    host_idx = 0
    pf_idx = 1
    concat_idx = 2
    genomes: ClassVar = [
        {
            "prefix": "host/hg38",
            "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
            "local_fn": "host/hg38.fa.gz",
            "fasta_fn": "host/hg38.fasta",
        },
        {
            "prefix": "PlasmoDB-44_Pfalciparum3D7_Genome",
            "url": "https://github.com/bguo068/snpeff_annot/releases/download/v26.09.14-alpha/PlasmoDB-44_Pfalciparum3D7.tgz",
            "local_fn": "PlasmoDB-44_Pfalciparum3D7_Genome.fasta",
            "fasta_fn": "PlasmoDB-44_Pfalciparum3D7_Genome.fasta",
        },
        {
            "prefix": "concat/hg38_pf3d7v44",
            "fasta_fn": "concat/hg38_pf3d7v44.fasta",
        },
    ]
    pf_crosses_v1: ClassVar = {
        "server": "ngs.sanger.ac.uk",
        "remotedir": "/production/malaria/pf-crosses/1.0/",
        "vcf_files": [
            "3d7_hb3.gatk.final.vcf.gz",
            "7g8_gb4.gatk.final.vcf.gz",
            "hb3_dd2.gatk.final.vcf.gz",
        ],
        "orig_files": [
            "pf_crosses_v1/3d7_hb3.gatk.final.vcf.gz",
            "pf_crosses_v1/7g8_gb4.gatk.final.vcf.gz",
            "pf_crosses_v1/hb3_dd2.gatk.final.vcf.gz",
        ],
        "filt_files": [
            "pf_crosses_v1/pass_3d7_hb3.gatk.final.vcf.gz",
            "pf_crosses_v1/pass_7g8_gb4.gatk.final.vcf.gz",
            "pf_crosses_v1/pass_hb3_dd2.gatk.final.vcf.gz",
        ],
        "known_variants_vcf": "pf_crosses_v1/known_variants.vcf",
    }


def download_genome_fasta_files():
    for genome in URLS.genomes:
        prefix = genome["prefix"]
        if "concat" in prefix:
            continue
        # local file
        url: str = genome["url"]
        local_fn: str = genome["local_fn"]
        target_fn: str = genome["fasta_fn"]

        Path(local_fn).parent.mkdir(parents=True, exist_ok=True)

        # skip if already exists
        if Path(target_fn).exists():
            continue

        if "https://github.com/bguo068/snpeff_annot/releases" in url:
            cmd = f"""
              wget -qO- --max-redirect=20 {url}  | tar -xzf - {target_fn}
            """
            run(cmd, shell=True, check=True)
            continue

        with requests.get(url, allow_redirects=True, stream=True) as response:
            response.raise_for_status()
            with open(local_fn, "wb") as f:
                f.writelines(response.iter_content(chunk_size=8192))

        if local_fn.endswith(".gz"):
            res = run(
                f"""gzip -dc {local_fn} > {target_fn}; rm {local_fn}""",
                shell=True,
                check=False,
            )
            assert res.returncode == 0


def concat_pf_human_genomes():
    host_g = URLS.genomes[URLS.host_idx]["fasta_fn"]
    pf_g = URLS.genomes[URLS.pf_idx]["fasta_fn"]
    concat_g = URLS.genomes[URLS.concat_idx]["fasta_fn"]
    if Path(concat_g).exists():
        return
    Path(concat_g).parent.mkdir(parents=True, exist_ok=True)
    cmd = f"cat {host_g} {pf_g} > {concat_g}"
    res = run(cmd, shell=True, capture_output=True, check=False)
    if res.returncode != 0:
        print(res.stderr)
        sys.exit(-1)


def index_fasta_files():
    for genome in URLS.genomes:
        fa_prefix = genome["prefix"]
        fa_fn = genome["fasta_fn"]
        # bt2 indices
        if not Path(fa_prefix + ".1.bt2").exists():
            print(Path(fa_prefix + ".1.bt2"))
            cmd = f"bowtie2-build --threads 20 {fa_prefix}.fasta {fa_prefix}"
            print(f"\t building bt2 indices for {fa_fn}")
            res = run(cmd, shell=True, capture_output=True, check=False)
            if res.returncode != 0:
                print(res.stderr)
                sys.exit(-1)
        # fai index
        if not Path(fa_prefix + ".fasta.fai").exists():
            cmd = f"samtools faidx {fa_prefix}.fasta"
            print(f"\t building fai index for {fa_fn}")
            res = run(cmd, shell=True, capture_output=True, check=False)
            if res.returncode != 0:
                print(res.stderr)
                sys.exit(-1)
        # dict index
        if not Path(fa_prefix + ".dict").exists():
            cmd = f"""
            gatk --java-options "-Xmx10G" CreateSequenceDictionary -R {fa_prefix}.fasta -O {fa_prefix}.dict
            """
            print(f"\t building dict index for {fa_fn}")
            res = run(cmd, shell=True, capture_output=True, check=True)
            if res.returncode != 0:
                print(res.stderr)
                sys.exit(-1)


def download_pf_crosses_v1_vcfs():
    # file names and folders
    remote_filename_lst = URLS.pf_crosses_v1["vcf_files"]
    local_filename_lst = URLS.pf_crosses_v1["orig_files"]
    for fn in local_filename_lst:
        Path(fn).parent.mkdir(parents=True, exist_ok=True)

    # DOWNLOAD VCF from pfcrosses1.0
    ftp_server = ftplib.FTP(URLS.pf_crosses_v1["server"])  # ty: ignore[invalid-argument-type]
    ftp_server.login()
    ftp_server.cwd(URLS.pf_crosses_v1["remotedir"])  # ty: ignore[invalid-argument-type]
    for remote_filename, local_filename in zip(remote_filename_lst, local_filename_lst):
        if not Path(local_filename).exists():
            with open(local_filename, "wb") as file:
                print("download " + remote_filename)
                ftp_server.retrbinary("RETR " + remote_filename, file.write)
    ftp_server.quit()

    # index vcf file with
    for local_fn in local_filename_lst:
        idx_fn = local_fn + ".tbi"
        if not Path(idx_fn).exists():
            cmd = f"gatk IndexFeatureFile -I {local_fn}"
            res = run(cmd, shell=True, capture_output=True, text=True, check=False)
            if (
                res.returncode != 0
                or "error" in res.stderr.lower()
                or "error" in res.stdout.lower()
            ):
                print(res.stderr)
                sys.exit(-1)


def known_variants_from_pf_crosses_v1():
    local_filename_lst = URLS.pf_crosses_v1["orig_files"]
    filt_filename_lst = URLS.pf_crosses_v1["filt_files"]

    # keep pass only
    for origf, passf in zip(local_filename_lst, filt_filename_lst):
        print("filter " + origf)
        _res = run(
            f"""
            bcftools view -f PASS {origf} -Oz -o {passf}
            bcftools index {passf}
            tabix {passf}
            """,
            shell=True,
            capture_output=True,
            text=True,
            check=True,
        )

    # merge vcfs and generate .idx file
    print("merge filtered vcfs and generate.idx file")
    filt_filename_str = " ".join(filt_filename_lst)
    known_variants_vcf = URLS.pf_crosses_v1["known_variants_vcf"]
    cmd = (
        f"""
        bcftools merge --force-samples {filt_filename_str}  -Ov -o {known_variants_vcf}
        gatk IndexFeatureFile -I {known_variants_vcf}
        """,
    )
    res = run(cmd, shell=True, capture_output=True, text=True, check=True)
    if (
        res.returncode != 0
        or "error" in res.stderr.lower()
        or "error" in res.stdout.lower()
    ):
        print(res.stderr)
        sys.exit(-1)


def check_if_all_files_exist_already():
    files = """
    PlasmoDB-44_Pfalciparum3D7_Genome.1.bt2
    PlasmoDB-44_Pfalciparum3D7_Genome.2.bt2
    PlasmoDB-44_Pfalciparum3D7_Genome.3.bt2
    PlasmoDB-44_Pfalciparum3D7_Genome.4.bt2
    PlasmoDB-44_Pfalciparum3D7_Genome.dict
    PlasmoDB-44_Pfalciparum3D7_Genome.fasta
    PlasmoDB-44_Pfalciparum3D7_Genome.fasta.fai
    PlasmoDB-44_Pfalciparum3D7_Genome.rev.1.bt2
    PlasmoDB-44_Pfalciparum3D7_Genome.rev.2.bt2
    concat/hg38_pf3d7v44.1.bt2
    concat/hg38_pf3d7v44.2.bt2
    concat/hg38_pf3d7v44.3.bt2
    concat/hg38_pf3d7v44.4.bt2
    concat/hg38_pf3d7v44.dict
    concat/hg38_pf3d7v44.fasta
    concat/hg38_pf3d7v44.fasta.fai
    concat/hg38_pf3d7v44.rev.1.bt2
    concat/hg38_pf3d7v44.rev.2.bt2
    host/hg38.1.bt2
    host/hg38.2.bt2
    host/hg38.3.bt2
    host/hg38.4.bt2
    host/hg38.dict
    host/hg38.fasta
    host/hg38.fasta.fai
    host/hg38.rev.1.bt2
    host/hg38.rev.2.bt2
    jacob2014_chip_sites.vcf
    jacob2014_chip_sites.vcf.idx
    pf_crosses_v1/3d7_hb3.gatk.final.vcf.gz
    pf_crosses_v1/3d7_hb3.gatk.final.vcf.gz.tbi
    pf_crosses_v1/7g8_gb4.gatk.final.vcf.gz
    pf_crosses_v1/7g8_gb4.gatk.final.vcf.gz.tbi
    pf_crosses_v1/hb3_dd2.gatk.final.vcf.gz
    pf_crosses_v1/hb3_dd2.gatk.final.vcf.gz.tbi
    pf_crosses_v1/known_variants.vcf
    pf_crosses_v1/known_variants.vcf.idx
    pf_crosses_v1/pass_3d7_hb3.gatk.final.vcf.gz
    pf_crosses_v1/pass_3d7_hb3.gatk.final.vcf.gz.csi
    pf_crosses_v1/pass_3d7_hb3.gatk.final.vcf.gz.tbi
    pf_crosses_v1/pass_7g8_gb4.gatk.final.vcf.gz
    pf_crosses_v1/pass_7g8_gb4.gatk.final.vcf.gz.csi
    pf_crosses_v1/pass_7g8_gb4.gatk.final.vcf.gz.tbi
    pf_crosses_v1/pass_hb3_dd2.gatk.final.vcf.gz
    pf_crosses_v1/pass_hb3_dd2.gatk.final.vcf.gz.csi
    pf_crosses_v1/pass_hb3_dd2.gatk.final.vcf.gz.tbi
    pf_crosses_v1/tmp/pass_3d7_hb3.gatk.final.vcf.gz
    pf_crosses_v1/tmp/pass_3d7_hb3.gatk.final.vcf.gz.csi
    pf_crosses_v1/tmp/pass_3d7_hb3.gatk.final.vcf.gz.tbi
    pf_crosses_v1/tmp/pass_7g8_gb4.gatk.final.vcf.gz
    pf_crosses_v1/tmp/pass_7g8_gb4.gatk.final.vcf.gz.csi
    pf_crosses_v1/tmp/pass_7g8_gb4.gatk.final.vcf.gz.tbi
    pf_crosses_v1/tmp/pass_hb3_dd2.gatk.final.vcf.gz
    pf_crosses_v1/tmp/pass_hb3_dd2.gatk.final.vcf.gz.csi
    pf_crosses_v1/tmp/pass_hb3_dd2.gatk.final.vcf.gz.tbi
    """
    for p in files.strip().split():
        if not Path(p).exists():
            return False
    return True


if __name__ == "__main__":
    if check_if_all_files_exist_already():
        print("all ref files already exist. exit with 0")
        sys.exit(0)

    print("downloading fasta files")
    download_genome_fasta_files()
    print("Done obtaining fasta files")

    print("concat host and parasite genome to a pseudo genome")
    concat_pf_human_genomes()
    print("Done concatenating the pseudo genome")

    print("indexing will take a long time ....")
    index_fasta_files()
    print("Done indexing")

    print("download and index pf crosses v1.0 vcf files")
    download_pf_crosses_v1_vcfs()
    print("Done download/index pf crosses 1.0")

    print("generating known variants vcf from pf crosses 1.0")
    known_variants_from_pf_crosses_v1()
    print("Done known variants vcf")
