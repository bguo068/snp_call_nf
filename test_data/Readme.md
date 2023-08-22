# Notes for creating small test data (fastq files)

Assumption: you have already aligned bam files so that we know the coordinates
of the reads

Steps:

1. Filter the bam file by region
    ```sh
    samtools index in.bam -@ 10
    samtools view -b  in.bam "chr_name:start-end" > filt.bam
    ```
2. Subsample reads using `-s` command (adjust the decimal part control file size):
    ```
    samtools view -s 42.10 -b filt.bam > filt_ss.bam
    ```
3. Sort BAM read by name:
    ```
    samtools sort -n filt_ss.bam -o filt_ss_sortn.bam
    ```
4. Convert a filtered BAM file back to paired-end FASTQ files
    ```
    samtools bam2fq -1 s_1.fastq.gz -2 s_2.fastq.gz -0 /dev/null filt_ss_sortn.bam
    ```