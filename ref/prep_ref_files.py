#! /usr/bin/env python
import ftplib
import os
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


def tmp_path(fn) -> str:
    """Temporary path next to `fn`, keeping the final extension so tools that
    infer the format from the extension still work: a/b.dict -> a/b.tmp.dict.

    Outputs are written to tmp paths and only renamed to their final names once
    complete, so an interrupted run never leaves a partial file that a later
    run would mistake for a finished one."""
    p = Path(fn)
    return str(p.with_name(f"{p.stem}.tmp{p.suffix}"))


def run_bash(cmd: str):
    # bash with pipefail so failures on the left side of a pipe are not masked
    return run(
        f"set -euo pipefail\n{cmd}",
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
        check=False,
    )


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

        target_tmp = tmp_path(target_fn)

        if "https://github.com/bguo068/snpeff_annot/releases" in url:
            # extract the member to stdout so it lands in the tmp file
            cmd = f"""
              wget -qO- --max-redirect=20 {url}  | tar -xzf - -O {target_fn} > {target_tmp}
            """
            res = run_bash(cmd)
            if res.returncode != 0:
                print(res.stderr)
                sys.exit(-1)
            os.replace(target_tmp, target_fn)
            continue

        local_tmp = tmp_path(local_fn)
        with requests.get(url, allow_redirects=True, stream=True) as response:
            response.raise_for_status()
            with open(local_tmp, "wb") as f:
                f.writelines(response.iter_content(chunk_size=8192))
        os.replace(local_tmp, local_fn)

        if local_fn.endswith(".gz"):
            res = run_bash(f"gzip -dc {local_fn} > {target_tmp}")
            if res.returncode != 0:
                print(res.stderr)
                sys.exit(-1)
            os.replace(target_tmp, target_fn)
            Path(local_fn).unlink()


def concat_pf_human_genomes():
    host_g = URLS.genomes[URLS.host_idx]["fasta_fn"]
    pf_g = URLS.genomes[URLS.pf_idx]["fasta_fn"]
    concat_g = URLS.genomes[URLS.concat_idx]["fasta_fn"]
    if Path(concat_g).exists():
        return
    Path(concat_g).parent.mkdir(parents=True, exist_ok=True)
    concat_tmp = tmp_path(concat_g)
    cmd = f"cat {host_g} {pf_g} > {concat_tmp}"
    res = run(cmd, shell=True, capture_output=True, check=False)
    if res.returncode != 0:
        print(res.stderr)
        sys.exit(-1)
    os.replace(concat_tmp, concat_g)


def index_fasta_files():
    for genome in URLS.genomes:
        fa_prefix = genome["prefix"]
        fa_fn = genome["fasta_fn"]
        # bt2 indices
        if not Path(fa_prefix + ".1.bt2").exists():
            print(Path(fa_prefix + ".1.bt2"))
            tmp_prefix = fa_prefix + ".tmp"
            cmd = f"bowtie2-build --threads 20 {fa_prefix}.fasta {tmp_prefix}"
            print(f"\t building bt2 indices for {fa_fn}")
            res = run(cmd, shell=True, capture_output=True, check=False)
            if res.returncode != 0:
                print(res.stderr)
                sys.exit(-1)
            # rename .1.bt2 last as it is the file checked for existence above
            tmp_files = sorted(
                Path(tmp_prefix).parent.glob(Path(tmp_prefix).name + ".*bt2*"),
                key=lambda p: p.name.endswith(".1.bt2"),
            )
            for p in tmp_files:
                os.replace(p, fa_prefix + p.name[len(Path(tmp_prefix).name) :])
        # fai index
        fai_fn = fa_prefix + ".fasta.fai"
        if not Path(fai_fn).exists():
            fai_tmp = tmp_path(fai_fn)
            cmd = f"samtools faidx {fa_prefix}.fasta --fai-idx {fai_tmp}"
            print(f"\t building fai index for {fa_fn}")
            res = run(cmd, shell=True, capture_output=True, check=False)
            if res.returncode != 0:
                print(res.stderr)
                sys.exit(-1)
            os.replace(fai_tmp, fai_fn)
        # dict index
        dict_fn = fa_prefix + ".dict"
        if not Path(dict_fn).exists():
            dict_tmp = tmp_path(dict_fn)
            Path(dict_tmp).unlink(missing_ok=True)  # gatk refuses to overwrite
            cmd = f"""
            gatk --java-options "-Xmx10G" CreateSequenceDictionary -R {fa_prefix}.fasta -O {dict_tmp}
            """
            print(f"\t building dict index for {fa_fn}")
            res = run(cmd, shell=True, capture_output=True, check=False)
            if res.returncode != 0:
                print(res.stderr)
                sys.exit(-1)
            os.replace(dict_tmp, dict_fn)


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
            local_tmp = tmp_path(local_filename)
            with open(local_tmp, "wb") as file:
                print("download " + remote_filename)
                ftp_server.retrbinary("RETR " + remote_filename, file.write)
            os.replace(local_tmp, local_filename)
    ftp_server.quit()

    # index vcf file with
    for local_fn in local_filename_lst:
        idx_fn = local_fn + ".tbi"
        if not Path(idx_fn).exists():
            idx_tmp = tmp_path(idx_fn)
            cmd = f"gatk IndexFeatureFile -I {local_fn} -O {idx_tmp}"
            res = run(cmd, shell=True, capture_output=True, text=True, check=False)
            if (
                res.returncode != 0
                or "error" in res.stderr.lower()
                or "error" in res.stdout.lower()
            ):
                print(res.stderr)
                sys.exit(-1)
            os.replace(idx_tmp, idx_fn)


def known_variants_from_pf_crosses_v1():
    local_filename_lst = URLS.pf_crosses_v1["orig_files"]
    filt_filename_lst = URLS.pf_crosses_v1["filt_files"]

    # keep pass only
    for origf, passf in zip(local_filename_lst, filt_filename_lst):
        print("filter " + origf)
        pass_tmp = tmp_path(passf)
        res = run_bash(
            f"""
            bcftools view -f PASS {origf} -Oz -o {pass_tmp}
            bcftools index -f {pass_tmp}
            tabix -f {pass_tmp}
            """
        )
        if res.returncode != 0:
            print(res.stderr)
            sys.exit(-1)
        # move indices first so the vcf only appears once everything is done
        os.replace(pass_tmp + ".csi", passf + ".csi")
        os.replace(pass_tmp + ".tbi", passf + ".tbi")
        os.replace(pass_tmp, passf)

    # merge vcfs and generate .idx file
    print("merge filtered vcfs and generate.idx file")
    filt_filename_str = " ".join(filt_filename_lst)
    known_variants_vcf = URLS.pf_crosses_v1["known_variants_vcf"]
    known_variants_tmp = tmp_path(known_variants_vcf)
    cmd = f"""
        bcftools merge --force-samples {filt_filename_str}  -Ov -o {known_variants_tmp}
        gatk IndexFeatureFile -I {known_variants_tmp}
        """
    res = run_bash(cmd)
    if (
        res.returncode != 0
        or "error" in res.stderr.lower()
        or "error" in res.stdout.lower()
    ):
        print(res.stderr)
        sys.exit(-1)
    os.replace(known_variants_tmp + ".idx", known_variants_vcf + ".idx")
    os.replace(known_variants_tmp, known_variants_vcf)


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
