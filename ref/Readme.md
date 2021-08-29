# Reference Genome Fasta Files and Index Files

1. `.fasta` files:
- hg38: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
- panTro6: https://hgdownload.soe.ucsc.edu/goldenPath/panTro6/bigZips/panTro6.fa.gz
- Download using wget
- Decompress with tar xf
- Rename to xx.fasta

2. fasta index,  bowtie2 index and gatk index files
- bowtie2 index (bt files): `bowtie2-build --threads 8 ref_genome.fasta ref_genome`
- fastq index (fai file): `samtoools faidx hg38.fasta`
- gatk .dict file: 
`gatk --java-options "-Xmx5G" CreateSequenceDictionary -R hg38.fasta -O hg38.dict`

# Known-site vcf files

1. Resource is mentioned in the `"##_COMMENT4": "KNOWN SITES RESOURCES"` block from
[](https:
//github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/JointGenotypingWf.hg38.inputs.json)
 and [](https:
//github.com/gatk-workflows/gatk4-genome-processing-pipeline/blob/master/WholeGenomeGermlineSingleSample.inputs.json)

2. Here 3 are used:
```
Homo_sapiens_assembly38.dbsnp138.vcf
Homo_sapiens_assembly38.known_indels.vcf.gz
Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```
The URLs are:
```
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
```
3. They can be downloaded by 
`gsutil cp gs://my-awesome-bucket/kitten.png Desktop/kitten2.png`

TODO: need to find the resource of the known site vcf for P.f.

# Local copies of these ref files can be found on 

- Rosalind:
- Thanos/IGS


