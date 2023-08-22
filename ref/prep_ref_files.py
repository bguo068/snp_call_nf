import requests
import ftplib
from pathlib import Path
from subprocess import run


class URLS:
    genomes = [
        {
            "prefix": "host/hg38",
            "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
            "local_fn": "host/hg38.fa.gz",
            "fasta_fn": "host/hg38.fasta",
        },
        {
            "prefix": "PlasmoDB-44_Pfalciparum3D7_Genome",
            "url": "https://plasmodb.org/common/downloads/release-44/Pfalciparum3D7/fasta/data/PlasmoDB-44_Pfalciparum3D7_Genome.fasta",
            "local_fn": "PlasmoDB-44_Pfalciparum3D7_Genome.fasta",
            "fasta_fn": "PlasmoDB-44_Pfalciparum3D7_Genome.fasta",
        },
    ]
    pf_crosses_v1 = {
        "server": "ngs.sanger.ac.uk",
        "remotedir": "/production/malaria/pf-crosses/1.0/",
        "vcf_files": "3d7_hb3.gatk.final.vcf.gz 7g8_gb4.gatk.final.vcf.gz hb3_dd2.gatk.final.vcf.gz".split(),
        "orig_files": [
            "pf_crosses_v1/3d7_hb3.gatk.final.vcf.gz",
            "pf_crosses_v1/7g8_gb4.gatk.final.vcf.gz",
            "pf_crosses_v1/hb3_dd2.gatk.final.vcf.gz",
        ],
        "filt_files": [
            "pf_crosses_v1/3d7_hb3.gatk.final.pass.vcf.gz",
            "pf_crosses_v1/7g8_gb4.gatk.final.pass.vcf.gz",
            "pf_crosses_v1/hb3_dd2.gatk.final.pass.vcf.gz",
        ],
        "known_variants_vcf": "pf_crosses_v1/known_variants.vcf",
    }
    pv_malgen = {
        "server": "ngs.sanger.ac.uk",
        "remotedir": "/production/malaria/Resource/30/Pv4_vcf",
        "vcf_files": """
            Pv4_PvP01_01_v1.vcf.gz
            Pv4_PvP01_02_v1.vcf.gz
            Pv4_PvP01_03_v1.vcf.gz 
            Pv4_PvP01_04_v1.vcf.gz
            Pv4_PvP01_05_v1.vcf.gz
            Pv4_PvP01_06_v1.vcf.gz
            Pv4_PvP01_07_v1.vcf.gz
            Pv4_PvP01_08_v1.vcf.gz
            Pv4_PvP01_09_v1.vcf.gz
            Pv4_PvP01_10_v1.vcf.gz
            Pv4_PvP01_11_v1.vcf.gz
            Pv4_PvP01_12_v1.vcf.gz
            Pv4_PvP01_13_v1.vcf.gz
            Pv4_PvP01_14_v1.vcf.gz""".split(),
         "orig_files": [
            "pv_malgen/Pv4_PvP01_01_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_02_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_03_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_04_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_05_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_06_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_07_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_08_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_09_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_10_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_11_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_12_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_13_v1.vcf.gz",
            "pv_malgen/Pv4_PvP01_14_v1.vcf.gz",
        ],
        "filt_files": [
            "pv_malgen/Pv4_PvP01_01_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_02_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_03_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_04_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_05_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_06_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_07_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_08_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_09_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_10_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_11_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_12_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_13_v1.pass.vcf.gz",
            "pv_malgen/Pv4_PvP01_14_v1.pass.vcf.gz",
        ],
        "known_variants_vcf": "pv_malgen/known_variants.vcf",
    }


def download_genome_fasta_files():
    for genome in URLS.genomes:
        # local file
        url: str = genome["url"]
        local_fn: str = genome["local_fn"]
        target_fn: str = genome["fasta_fn"]

        Path(local_fn).parent.mkdir(parents=True, exist_ok=True)

        # skip if already exists
        if Path(target_fn).exists():
            continue

        with requests.get(url, allow_redirects=True, stream=True) as response:
            response.raise_for_status()
            with open(local_fn, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

        if local_fn.endswith(".gz"):
            res = run(
                f"""gzip -dc {local_fn} > {target_fn}; rm {local_fn}""", shell=True
            )
            assert res.returncode == 0


def index_fasta_files():
    for genome in URLS.genomes:
        fa_prefix = genome["prefix"]
        fa_fn = genome["fasta_fn"]
        # bt2 indices
        if not Path(fa_prefix + ".1.bt2").exists():
            print(Path(fa_prefix + ".1.bt2"))
            cmd = f"bowtie2-build --threads 20 {fa_prefix}.fasta {fa_prefix}"
            print(f"\t building bt2 indices for {fa_fn}")
            res = run(cmd, shell=True, capture_output=True)
            if res.returncode != 0:
                print(res.stderr)
                exit(-1)
        # fai index
        if not Path(fa_prefix + ".fasta.fai").exists():
            cmd = f"samtools faidx {fa_prefix}.fasta"
            print(f"\t building fai index for {fa_fn}")
            res = run(cmd, shell=True, capture_output=True)
            if res.returncode != 0:
                print(res.stderr)
                exit(-1)
        # dict index
        if not Path(fa_prefix + ".dict").exists():
            cmd = f"""
            gatk --java-options "-Xmx10G" CreateSequenceDictionary -R {fa_prefix}.fasta -O {fa_prefix}.dict
            """
            print(f"\t building dict index for {fa_fn}")
            res = run(cmd, shell=True, capture_output=True)
            if res.returncode != 0:
                print(res.stderr)
                exit(-1)


def download_pf_crosses_v1_vcfs():
    # file names and folders
    remote_filename_lst = URLS.pf_crosses_v1["vcf_files"]
    local_filename_lst = URLS.pf_crosses_v1["orig_files"]
    for fn in local_filename_lst:
        Path(fn).parent.mkdir(parents=True, exist_ok=True)

    # DOWNLOAD VCF from pfcrosses1.0
    ftp_server = ftplib.FTP(URLS.pf_crosses_v1["server"])
    ftp_server.login()
    ftp_server.cwd(URLS.pf_crosses_v1["remotedir"])
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
            res = run(cmd, shell=True, capture_output=True, text=True)
            if (
                res.returncode != 0
                or "error" in res.stderr.lower()
                or "error" in res.stdout.lower()
            ):
                print(res.stderr)
                exit(-1)


def known_variants_from_pf_crosses_v1():
    local_filename_lst = URLS.pf_crosses_v1["orig_files"]
    filt_filename_lst = URLS.pf_crosses_v1["filt_files"]

    # keep pass only
    for origf, passf in zip(local_filename_lst, filt_filename_lst):
        print("filter " + origf)
        res = run(
            f"""
            bcftools view -f PASS {origf} -Oz -o {passf}
            bcftools index {passf}
            """,
            shell=True,
            capture_output=True,
            text=True,
        )
        assert res.returncode == 0

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
    res = run(cmd, shell=True, capture_output=True, text=True)
    if (
        res.returncode != 0
        or "error" in res.stderr.lower()
        or "error" in res.stdout.lower()
    ):
        print(res.stderr)
        exit(-1)

def download_pv_malgen_vcfs():
    # file names and folders
    remote_filename_lst = URLS.pv_malgen["vcf_files"]
    local_filename_lst = URLS.pv_malgen["orig_files"]
    for fn in local_filename_lst:
        Path(fn).parent.mkdir(parents=True, exist_ok=True)

    # download vcfs
    ftp_server = ftplib.FTP(URLS.pv_malgen["server"])
    ftp_server.login()
    ftp_server.cwd(URLS.pv_malgen["remotedir"])
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
            res = run(cmd, shell=True, capture_output=True, text=True)
            if (
                res.returncode != 0
                or "error" in res.stderr.lower()
                or "error" in res.stdout.lower()
            ):
                print(res.stderr)
                exit(-1)
                
def known_variants_from_pv_malgen():
    local_filename_lst = URLS.pv_malgen["orig_files"]
    filt_filename_lst = URLS.pv_malgen["filt_files"]

    # keep pass only
    for origf, passf in zip(local_filename_lst, filt_filename_lst):
        print("filter " + origf)
        res = run(
            f"""
            bcftools view -f PASS {origf} -Oz -o {passf}
            bcftools index {passf}
            """,
            shell=True,
            capture_output=True,
            text=True,
        )
        assert res.returncode == 0

    # merge vcfs and generate .idx file
    print("merge filtered vcfs and generate.idx file")
    filt_filename_str = " ".join(filt_filename_lst)
    known_variants_vcf = URLS.pv_malgen["known_variants_vcf"]
    cmd = (
        f"""
        bcftools concat {filt_filename_str}  -Ov -o {known_variants_vcf}
        gatk IndexFeatureFile -I {known_variants_vcf}
        """,
    )
    res = run(cmd, shell=True, capture_output=True, text=True)
    if (
        res.returncode != 0
        or "error" in res.stderr.lower()
        or "error" in res.stdout.lower()
    ):
        print(res.stderr)
        exit(-1)


if __name__ == "__main__":

    print("downloading fasta files")
    download_genome_fasta_files()
    print("Done obtaining fasta files")

    print("indexing will take a long time ....")
    index_fasta_files()
    print("Done indexing")

    print("download and index pf crosses v1.0 vcf files")
    download_pf_crosses_v1_vcfs()
    print("Done download/index pf crosses 1.0")

    print("generating known variants vcf from pf crosses 1.0")
    known_variants_from_pf_crosses_v1()
    print("Done known variants vcf")
