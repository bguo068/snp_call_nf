from subprocess import run
import ftplib
from pathlib import Path

# file names and folders
Path("orig").mkdir(parents=True, exist_ok=True)
Path("filt").mkdir(parents=True, exist_ok=True)
remote_filename_lst = "3d7_hb3.gatk.final.vcf.gz 7g8_gb4.gatk.final.vcf.gz hb3_dd2.gatk.final.vcf.gz".split()
local_filename_lst = [f"orig/{fn}" for fn in remote_filename_lst]
filt_filename_lst = [
    f"filt/{Path(fn).name.replace('final', 'final.pass')}" for fn in remote_filename_lst
]


# DOWNLOAD VCF from pfcrosses1.0

ftp_server = ftplib.FTP("ngs.sanger.ac.uk")
ftp_server.login()
ftp_server.cwd("/production/malaria/pf-crosses/1.0/")
for remote_filename, local_filename in zip(remote_filename_lst, local_filename_lst):
    with open(local_filename, "wb") as file:
        print("download " + remote_filename)
        ftp_server.retrbinary("RETR " + remote_filename, file.write)
ftp_server.quit()


# FILTER and MERGE VCF files

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

print("merge filtered vcfs")
run(
    "bcftools merge --force-samples {}  -Oz -o known_variants.vcf.gz".format(
        " ".join(filt_filename_lst)
    ),
    shell=True,
)
assert res.returncode == 0
