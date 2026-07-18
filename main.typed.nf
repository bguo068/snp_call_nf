#! /usr/bin/env nextflow

nextflow.enable.types = true
nextflow.enable.dsl = 2

def print_parameters() -> Void {
    println(
        """======================== PARAMETERS ==========================
fq_map:\t${file(params.fq_map)}
known_variants:\t${params.known_sites}
split:\t${params.split}
filtering:
\thard:\t ${params.hard}
\tvqsr:\t ${params.vqsr}
use_concat_genome: ${params.use_concat_genome}
outdir:\t${params.outdir}
tmpdir:\t${params.gatk_tmpdir}
parasite_reads_only:\t${params.parasite_reads_only}
coverage_only:\t ${params.coverage_only}
gvcf_only:\t${params.gvcf_only}
============================================================="""
    )
    return null
}

record SampleMeta {
    Sample: String
    Run: String
    Host_ref_prefix: Path?
    Concat_ref_prefix: Path?
    Parasite_ref_prefix: Path?
    is_paired: Integer
}

process BOWTIE2_ALIGN_TO_CONCAT_GENOME {
    tag "${meta.Sample}~${meta.Run}"

    input:
    tuple(meta: SampleMeta, fastq: List<Path>)

    output:
    rec = record(
        concat_bam: tuple(meta, file("*to_concat.bam")),
        parasite_bam: tuple(meta, file("*to_parasite.bam")),
        flagstat: tuple(meta, file("*.flagstat")),
        read_length: tuple(meta, file("*_read_length.txt")),
    )

    script:
    def reads = meta.is_paired == 1 ? "-1 ${fastq[0]} -2 ${fastq[1]}" : "-U ${fastq}"
    def read_group = "--rg-id ${meta.Sample} --rg SM:${meta.Sample} --rg PL:Illumina"
    """
    #! /usr/bin/env bash

    set -eEux -o pipefail
    # align to host
    bowtie2 -x ${meta.Concat_ref_prefix}  ${reads} ${read_group} -p ${task.cpus}|\
        samtools view -q 0 -bS|\
        samtools sort --threads ${task.cpus} -o ${meta.Sample}~${meta.Run}~to_concat.bam -
    
    # get read length
    set +o pipefail
    samtools view ${meta.Sample}~${meta.Run}~to_concat.bam|\
         awk '{a+=length(\$10)} {if (NR>100000) {print (a/NR); exit}} END{print (a/NR)}' \
        > ${meta.Sample}~${meta.Run}_read_length.txt
    set -o pipefail

    # get parasite chromosome names
    parasite_chrnames=`grep -E '^>' ${meta.Parasite_ref_prefix}.fasta.sed 's:^>\\s*::'.sed 's:\\s.*\$::'.tr '\\n' ' '`
    echo \$parasite_chrnames

    # subset to alignment to parasite chromosome
    samtools index ${meta.Sample}~${meta.Run}~to_concat.bam
    samtools view -b ${meta.Sample}~${meta.Run}~to_concat.bam  \${parasite_chrnames}  -o ${meta.Sample}~${meta.Run}~to_parasite.bam  

    # get flatstats
    samtools flagstat ${meta.Sample}~${meta.Run}~to_parasite.bam \
        > ~to_parasite.bam.flagstat
    rm -f sorted.bam
    """

    stub:
    """ 
    touch ${meta.Sample}~${meta.Run}~to_concat.bam{,.flagstat}
    touch ${meta.Sample}~${meta.Run}_read_length.txt
    """
}



process BOWTIE2_ALIGN_TO_HOST {
    tag "${meta.Sample}~${meta.Run}"

    input:
    tuple(meta: SampleMeta, fastq: List<Path>)

    output:
    record(
        meta: meta,
        bam: file("*.bam"),
        flagstat: file("*.flagstat"),
        read_length: file("*_read_length.txt"),
    )

    script:
    def reads = meta.is_paired == 1 ? "-1 ${fastq[0]} -2 ${fastq[1]}" : "-U ${fastq}"
    """
    bowtie2 -x ${meta.Host_ref_prefix}  ${reads} -p ${task.cpus}|\
        samtools view -q 0 -bS > ${meta.Sample}~${meta.Run}~to_host.bam
    samtools flagstat ${meta.Sample}~${meta.Run}~to_host.bam \
        > ${meta.Sample}~${meta.Run}~to_host.bam.flagstat
    samtools view ${meta.Sample}~${meta.Run}~to_host.bam|\
        head -100000|awk '{a+=length(\$10)}END{print (a/NR)}' \
        > ${meta.Sample}~${meta.Run}_read_length.txt

    """

    stub:
    """ 
    touch ${meta.Sample}~${meta.Run}~to_host.bam{,.flagstat}
    touch ${meta.Sample}~${meta.Run}_read_length.txt
    """
}

process SAMTOOLS_VIEW_RM_HOST_READS {
    tag "${meta.Sample}~${meta.Run}"

    input:
    tuple(meta: SampleMeta, bam: Path)

    output:
    bam = tuple(meta, file("*.bam"))

    script:
    def flag = meta.is_paired ? "-f 12 -F 256" : "-f 4"
    """
    samtools view -b ${flag} ${bam}|\
        samtools sort -@ ${task.cpus} -n -o ${meta.Sample}~${meta.Run}~unmapped_sorted.bam
    """

    stub:
    """ touch ${meta.Sample}~${meta.Run}~unmapped_sorted.bam """
}

process SAMTOOLS_FASTQ {
    tag "${meta.Sample}~${meta.Run}"

    input:
    tuple(meta: SampleMeta, bam: Path)

    output:
    fq = tuple(meta, files("*.fastq.gz").toSorted())

    script:
    if (meta.is_paired) {
        """
    samtools fastq -@ ${task.cpus} ${bam} \
        -1 f1.fastq.gz -2 f2.fastq.gz -0 /dev/null -s /dev/null -n	
    """
    }
    else {
        """
    samtools fastq -@ ${task.cpus} ${bam} \
        -o /dev/null -s /dev/null -0 f.fastq.gz -n  
    """
    }

    stub:
    if (meta.is_paired) {
        """ touch f1.fastq.gz f2.fastq.gz """
    }
    else {
        """ touch f.fastq.gz """
    }
}

process BOWTIE2_ALIGN_TO_PARASITE {
    tag "${meta.Sample}~${meta.Run}"

    input:
    tuple(meta: SampleMeta, fastq: List<Path>)
    parasite_ref_prefix: String

    output:
    record(
        meta: meta,
        bam: file("*.bam"),
        flagstat: file("*.flagstat"),
    )

    script:
    def reads = meta.is_paired == 1 ? "-1 ${fastq[0]} -2 ${fastq[1]}" : "-U ${fastq}"
    def read_group = "--rg-id ${meta.Sample} --rg SM:${meta.Sample} --rg PL:Illumina"
    """
    bowtie2 -x ${parasite_ref_prefix} ${reads} ${read_group} -p ${task.cpus}|\
        samtools view -q 0 -bS > ${meta.Sample}~${meta.Run}~to_parasite.bam
    samtools flagstat ${meta.Sample}~${meta.Run}~to_parasite.bam \
        > ${meta.Sample}~${meta.Run}~to_parasite.bam.flagstat
    """

    stub:
    """ touch ${meta.Sample}~${meta.Run}~to_parasite.bam{,.flagstat} """
}

process PICARD_MERGE_SORT_BAMS {
    tag "${sample}"

    input:
    tuple(sample: String, bams: Bag<Path>)

    output:
    tuple(sample, file("merged_sorted.bam"))

    script:
    def to_merge = bams.size() > 1
    def in_bams = bams.join(" -I ")
    def merge_opts = "-VALIDATION_STRINGENCY LENIENT -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
    def sort_opts = "-VALIDATION_STRINGENCY LENIENT -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
    """
    if ${to_merge}; then
        gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.toGiga()}G" \
	MergeSamFiles -I ${in_bams} -O merged.bam ${merge_opts}
    else
       ln -s ${bams.join()} merged.bam
    fi

    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.toGiga()}G" \
        SortSam -I merged.bam -O merged_sorted.bam --SORT_ORDER coordinate ${sort_opts}
    rm merged.bam
    """

    stub:
    """ touch merged_sorted.bam """
}

process PICARD_MARK_DUPLICATES {
    tag "${sample}"

    input:
    tuple(sample: String, bam: Path)

    output:
    bam = tuple(sample, file("*dedup.bam"))

    script:
    def opts = "--USE_JDK_DEFLATER true --USE_JDK_INFLATER true"
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.toGiga()}G" \
    MarkDuplicates \
	--REMOVE_DUPLICATES ${opts} \
        -I ${bam} -O ${sample}_dedup.bam --METRICS_FILE dedup_metrics.txt
    """

    stub:
    """ touch ${sample}_dedup.bam """
}

process GATK_BASE_RECALIBRATOR {
    tag "${sample}"

    input:
    tuple(sample: String, bam: Path)
    ref: String
    known_sites: List<String>

    output:
    bam_rectbl = tuple(sample, file(bam.name), file("recal_data.table"))

    script:
    def known_sites_str = known_sites.join(" --known-sites ")
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir}  -Xmx${task.memory.toGiga()}G" \
        BaseRecalibrator -I ${bam} -O recal_data.table \
	-R ${ref} --known-sites ${known_sites_str} 
    """

    stub:
    """ touch recal_data.table """
}

process GATK_APPLY_BQSR {
    tag "${sample}"
    publishDir "${params.outdir}/recalibrated"

    input:
    tuple(sample: String, bam: Path, recal_table: Path)
    ref: String

    output:
    bam = tuple(sample, file("*recalibrated.bam"))

    script:
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.toGiga()}G" \
        ApplyBQSR -I ${bam} -O ${sample}_recalibrated.bam \
        -R ${ref} --bqsr-recal-file ${recal_table}
    """

    stub:
    """ touch ${sample}_recalibrated.bam """
}

process BEDTOOLS_GENOMECOV {
    tag "${sample}"
    publishDir "${params.outdir}/coverage"

    input:
    tuple(sample: String, bam: Path)
    ref: String

    output:
    record(
        bedgraph: file("*recalibrated.coverage.BedGraph.gz"),
        summary: file("*.coverage.summary.txt"),
    )

    script:
    """
    bedtools genomecov -bg -ibam ${bam} -g ${ref}|\
        gzip -v > ${sample}_recalibrated.coverage.BedGraph.gz

    # summarize the coverage
    zcat ${sample}_recalibrated.coverage.BedGraph.gz|\
        awk -v chrom_reg="${params.chrom_reg}" \
            -v genome_size_bp=${params.genome_size_bp} \
            -v sample="${sample}" \
        '  BEGIN {
                # exit in begin block wont affect the end block; need call exit again in end block
                if(chrom_reg == "" || genome_size_bp == "" ||  sample == "") {exit 1;}

                cov5x=0; cov10x=0; cov25x=0; cov50x=0; cov100x=0; area = 0; lines_used = 0;
            }
            NF == 4 && \$1 ~ chrom_reg {
                L = \$3 - \$2;
                d = \$4;
                if (d >=   5) {cov5x   += L; }
                if (d >=  10) {cov10x  += L; }
                if (d >=  25) {cov25x  += L; }
                if (d >=  50) {cov50x  += L; }
                if (d >= 100) {cov100x += L; }
                area += L * d;
                lines_used += 1;
            }
            END {
                # exit in begin block wont affect the end block; need call exit again in end block
                if(chrom_reg == "" || genome_size_bp == "" ||  sample == "") {print "input variables are not set"; exit 1;}
                avg_depth = area / genome_size_bp
                OFS = ","
                # print header
                print "#sample", "cov5x", "cov10x", "cov25x", "cov50x", "cov100x", "area", "genome_size_bp", "avg_depth", "lines_used", "lines_total"
                # print values
                print sample, cov5x, cov10x, cov25x, cov50x, cov100x, area, genome_size_bp, avg_depth, lines_used, NR
            } ' \
        >  ${sample}_recalibrated.coverage.summary.txt

    """

    stub:
    """ 
    touch ${sample}_recalibrated.coverage.BedGraph.gz 
    touch ${sample}_recalibrated.coverage.summary.txt
    """
}

process SAMTOOLS_FLAGSTAT {
    tag "${sample}"
    publishDir "${params.outdir}/flagstat"

    input:
    tuple(sample: String, bam: Path)

    output:
    file("*.flagstat")

    script:
    """
    samtools flagstat ${bam} > ${sample}.flagstat
    """

    stub:
    """ touch ${sample}.flagstat """
}

process GATK_HAPLOTYPE_CALLER {
    tag "${sample}"
    publishDir "${params.outdir}/gvcf"

    input:
    tuple(sample: String, bam: Path)
    ref: String

    output:
    gvcf = tuple(sample, file("*.g.vcf"), file("*.g.vcf.idx"))

    script:
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.toGiga()}G" \
        HaplotypeCaller \
	-I ${bam} -O ${sample}.g.vcf -R ${ref} -ERC GVCF
    """

    stub:
    """ touch ${sample}.g.vcf{,.idx} """
}

process GATK_GENOMICS_DB_IMPORT {
    tag "${interval}"

    input:
    interval: String
    gvcf_map: Path

    output:
    dbdir = file("*", type: 'dir', maxDepth: 1)

    script:
    def dbname = interval.replaceAll(":", "~")
    def mem = Math.round(task.memory.toGiga() * 0.75)
    // the rest of memory for c/c++ library
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${mem}G" \
        GenomicsDBImport \
        --batch-size 100 --reader-threads 5 --consolidate true \
        --sample-name-map ${gvcf_map} \
        --genomicsdb-workspace-path ${dbname} \
        -L ${interval}
    """

    stub:
    def dbname = interval.replaceAll(":", "~")
    """ mkdir ${dbname} """
}

process GATK_GENOTYPE_GVCFS {
    tag "${db.getName()}"

    input:
    db: Path
    ref: String

    output:
    db_vcf = tuple(env('DBNAME'), file("*.vcf"), file("*.idx"))

    script:
    def dbname = db.getName()
    def mem = Math.round(task.memory.toGiga() * 0.75)
    // the rest of memory for TileDB library

    // use environmental variable to pass value to output channel
    """
    DBNAME=${dbname}
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${mem}G" \
        GenotypeGVCFs -V gendb://${dbname} -O ${dbname}.vcf -R ${ref}
    """

    stub:
    def dbname = db.getName()
    """ touch ${dbname}.vcf{,.idx}; DBNAME="${dbname}" """
}

process GATK_SELECT_VARIANTS {
    tag "${dbname}"

    input:
    tuple(dbname: String, vcf: Path, idx: Path)
    ref: String

    output:
    db_vcf = tuple(dbname, file("*.snp.vcf"), file("*.snp.vcf.idx"))

    script:
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.toGiga()}G" \
        SelectVariants --select-type-to-include SNP \
        -R ${ref} -V ${vcf} -O ${dbname}.snp.vcf
    """

    stub:
    """touch ${dbname}.snp.vcf{,.idx} """
}

process GATK_VARIANT_FILTRATION {
    tag "${dbname}"

    input:
    tuple(dbname: String, vcf: Path, idx: Path)
    ref: String

    output:
    tuple(dbname, file("*.snp.hardfilt.vcf"), file("*.snp.hardfilt.vcf.idx"))

    script:
    def filter_str = ''
    params.hard_filters.each { filt -> filter_str += " -filter \"${filt.filter}\" --filter-name \"${filt.name}\" " }
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.toGiga()}G" VariantFiltration \
    ${filter_str} \
	-R ${ref} \
	-V ${vcf} \
	-O ${dbname}.snp.hardfilt.vcf
    """

    stub:
    """ touch ${dbname}.snp.hardfilt.vcf{,.idx}"""
}

def get_vqsr_resources(resources) -> String {
    def str = ''
    resources.each { it ->
        def type_str = ''
        if (it.type == "truth") {
            type_str = "known=false,training=true,truth=true"
        }
        else if (it.type == "training") {
            type_str = "known=false,training=true,truth=false"
        }
        else if (it.type == "known") {
            type_str = "known=true,training=false,truth=false"
        }
        else {
            print("vqsr resources type can only be one of known, training, known")
            System.exit(-10)
        }
        def vcf_file = file(it.vcf)

        str += " --resource:${it.name},${type_str},prior=${it.prior} ${vcf_file} "
    }
    return str
}

process GATK_VARIANT_RECALIBRATOR {
    tag "${dbname}"

    input:
    tuple(dbname: String, vcf: Path, idx: Path)
    resources: Map<String, String>
    opts: String
    mode: String
    ref: String

    output:
    recalvcf = tuple(dbname, file("*.recal.vcf"), file("*.tranches"))

    script:
    def resources_str = get_vqsr_resources(resources)
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.toGiga()}G" VariantRecalibrator \
        -R ${ref} -V ${vcf} ${resources_str} ${opts} -mode ${mode} \
        -O ${dbname}.recal.vcf --tranches-file ${dbname}.all.tranches
    """

    stub:
    """ touch ${dbname}.recal.vcf  ${dbname}.all.tranches """
}
process GATK_APPLY_VQSR {
    // publishDir "${params.outdir}/vqsrfilt"
    tag "${dbname}"

    input:
    tuple(dbname: String, vcf: Path, vcf_idx: Path, recal: Path, tranches: Path)
    mode: String
    ref: String

    output:
    tuple(dbname, file("*.vqsrfilt.*.vcf"))

    script:
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.toGiga()}G" IndexFeatureFile -I ${recal}
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.toGiga()}G" ApplyVQSR \
        -R ${ref} -V ${vcf} --recal-file ${recal} --tranches-file ${tranches} -mode ${mode} \
        --output ${dbname}.vqsrfilt.${mode}.vcf
    """

    stub:
    """ touch ${dbname}.vqsrfilt.${mode}.vcf """
}

workflow {

    // print key parameters
    print_parameters()

    // Get paths
    // Not use String instead Path type so that associated files can be corrected found
    def paths = [:]
    paths.parasite = [
        fasta: file(params.parasite.fasta).toUriString(),
        fasta_prefix: file(params.parasite.fasta_prefix).toUriString(),
    ]
    paths.host = [fasta: [], fasta_prefix: []]
    params.host.fasta.each { it -> paths.host.fasta.add(file(it).toUriString()) }
    params.host.fasta_prefix.each { it -> paths.host.fasta_prefix.add(file(it).toUriString()) }
    paths.concat = [fasta: [], fasta_prefix: []]
    params.concat.fasta.each { it -> paths.concat.fasta.add(file(it).toUriString()) }
    params.concat.fasta_prefix.each { it -> paths.concat.fasta_prefix.add(file(it).toUriString()) }
    paths.known_sites = []
    params.known_sites.each { it -> paths.known_sites.add(file(it).toUriString()) }

    // prepare tmpdir for gatk
    def tmpdir = file(params.gatk_tmpdir)
    if (!tmpdir.exists()) {
        tmpdir.mkdirs()
    }

    if (params.use_concat_genome) {
        // Prepare input chanel
        // csv contains sample, host_id, run, mate_id, fq_path
        input_ch = channel.fromPath(params.fq_map)
            .flatMap { csv -> csv.splitCsv(skip: 1, sep: '\t') }
            .map { row ->
                def fields = row as List<?>
                tuple(
                    tuple(fields[0].toString(), fields[1].toInteger(), fields[2].String()),
                    fields[4].toString(),
                )
            }
            .groupBy()
            .map { key, fq_unsorted ->
                def (sample, host_id, run) = key
                def fq = fq_unsorted.toSorted()
                assert (fq.size() == 1) || (fq.size() == 2) : "number of fq files for each run can only be 1 or 2! Found ${fq.size()} for ${sample}/${run}\n"
                def meta = record(
                    Sample: sample,
                    Run: run,
                    Concat_ref_prefix: paths.concat.fasta_prefix[host_id.toInteger()],
                    Parasite_ref_prefix: paths.parasite.fasta_prefix,
                    is_paired: fq.size() == 2 ? 1 : 0,
                )
                def fq_paths = fq.collect { it -> file(it) }
                return tuple(meta, fq_paths)
            }
        out_BOWTIE2_ALIGN_TO_CONCAT_GENOME = BOWTIE2_ALIGN_TO_CONCAT_GENOME(input_ch)
        out_SAMTOOLS_FASTQ = SAMTOOLS_FASTQ(out_BOWTIE2_ALIGN_TO_CONCAT_GENOME.map { rec -> rec.parasite_bam })
    }
    else {
        // Prepare input chanel
        // csv contains sample, host_id, run, mate_id, fq_path
        input_ch = channel.fromPath(params.fq_map)
            .flatMap { csv -> csv.splitCsv(skip: 1, sep: '\t') }
            .map { row ->
                def fields = row as List<?>
                tuple(
                    tuple(fields[0].toString(), fields[1].toInteger(), fields[2].toString()),
                    fields[4].toString(),
                )
            }
            .groupBy()
            .map { key, fq_unsorted ->
                def (sample, host_id, run) = key
                def fq = fq_unsorted.toSorted()
                assert (fq.size() == 1) || (fq.size() == 2) : "number of fq files for each run can only be 1 or 2! Found ${fq.size()} for ${sample}/${run}\n"
                def meta = record(
                    Sample: sample,
                    Run: run,
                    Host_ref_prefix: paths.host.fasta_prefix[host_id.toInteger()],
                    is_paired: fq.size() == 2 ? 1 : 0,
                )
                def fq_paths = fq.collect { it -> file(it) }
                return tuple(meta, fq_paths)
            }

        // Remove reads mapped to host and split unmapped to fastq files
        out_BOWTIE2_ALIGN_TO_HOST = BOWTIE2_ALIGN_TO_HOST(input_ch)
        out_SAMTOOLS_VIEW_RM_HOST_READS = SAMTOOLS_VIEW_RM_HOST_READS(out_BOWTIE2_ALIGN_TO_HOST.map { rec -> tuple(rec.meta, rec.bam) })
        out_SAMTOOLS_FASTQ = SAMTOOLS_FASTQ(out_SAMTOOLS_VIEW_RM_HOST_READS)
    }

    ch_parasite_reads = out_SAMTOOLS_FASTQ
    if (params.parasite_reads_only) {
        ch_parasite_reads = channel.empty()
    }
    // Align to parasite genome
    out_BOWTIE2_ALIGN_TO_PARASITE = BOWTIE2_ALIGN_TO_PARASITE(ch_parasite_reads, paths.parasite.fasta_prefix)

    // Non-blocking grouped gathering
    // csv contains sample, host_id, run, mate_id, fq_path
    n_run = channel.fromPath(params.fq_map)
        .flatMap { csv -> csv.splitCsv(skip: 1, sep: '\t') }
        .map { row ->
            def fields = row as List<?>
            // sample, run
            tuple(fields[0].toString(), fields[2].toString())
        }
        .unique()
        .groupBy()
        .map { sample, runs ->
            record(sample: sample, num: runs.size())
        }

    merge_input = out_BOWTIE2_ALIGN_TO_PARASITE.map { rec -> record(sample: rec.meta.Sample, bam: rec.bam) }.join(n_run, by: "sample").map { rec -> tuple(rec.sample, rec.num, rec.bam) }.groupBy()

    // For each sample, merge and sort bam files of different runs 
    out_PICARD_MERGE_SORT_BAMS = PICARD_MERGE_SORT_BAMS(merge_input)

    out_PICARD_MARK_DUPLICATES = PICARD_MARK_DUPLICATES(out_PICARD_MERGE_SORT_BAMS)

    // Base recalibration
    out_GATK_BASE_RECALIBRATOR = GATK_BASE_RECALIBRATOR(out_PICARD_MARK_DUPLICATES, paths.parasite.fasta, paths.known_sites)
    out_GATK_APPLY_BQSR = GATK_APPLY_BQSR(out_GATK_BASE_RECALIBRATOR, paths.parasite.fasta)

    // Generate stat files
    BEDTOOLS_GENOMECOV(out_GATK_APPLY_BQSR, paths.parasite.fasta)
    SAMTOOLS_FLAGSTAT(out_GATK_APPLY_BQSR)

    // stop before GATK_HAPLOTYPE_CALLER if only coverage inforation is need
    def ch_bqsr_bam: Channel<Tuple<String, Path>>
    if (params.coverage_only) {
        ch_bqsr_bam = channel.fromList([])
    }
    else {
        ch_bqsr_bam = out_GATK_APPLY_BQSR
    }

    // Generate gvcf
    out_GATK_HAPLOTYPE_CALLER = GATK_HAPLOTYPE_CALLER(ch_bqsr_bam, paths.parasite.fasta)

    if (!params.gvcf_only) {
        // Collect information to make gvcf_map file
        gvcf_map_ch = out_GATK_HAPLOTYPE_CALLER
            .map { sample, gvcf, _idx -> "${sample}\t${gvcf}" }
            .collect()
            .map { lines ->
                def out_file = file("gvcf_map.txt")
                out_file.text = lines.toSorted().join("\n") + "\n"
                return out_file
            }
        // .collectFile(name: "gvcf_map.txt", newLine: true, sort: true).first()


        // Import gvcf files to genomicsdb
        def interval_ch: Channel<String>
        interval_ch = channel.fromList(params.genome_intervals[params.split])
        out_GATK_GENOMICS_DB_IMPORT = GATK_GENOMICS_DB_IMPORT(interval_ch, gvcf_map_ch)

        // Genotype gvcf genomics db
        out_GATK_GENOTYPE_GVCFS = GATK_GENOTYPE_GVCFS(out_GATK_GENOMICS_DB_IMPORT, paths.parasite.fasta)

        // Select SNP only
        out_GATK_SELECT_VARIANTS = GATK_SELECT_VARIANTS(out_GATK_GENOTYPE_GVCFS, paths.parasite.fasta)

        // Hard filtering and Keep 'PASS' variants
        if (params.hard == true) {
            GATK_VARIANT_FILTRATION(out_GATK_SELECT_VARIANTS, paths.parasite.fasta)
        }

        // VQSR variant filtering
        if (params.vqsr == true) {
            out_GATK_VARIANT_RECALIBRATOR = GATK_VARIANT_RECALIBRATOR(
                out_GATK_SELECT_VARIANTS,
                params.vqsr_resources,
                params.vqsr_opts,
                params.vqsr_mode,
                paths.parasite.fasta,
            )
            GATK_APPLY_VQSR(
                out_GATK_SELECT_VARIANTS.map { dbname, vcf, vcf_idx -> record(dbname: dbname, vcf: vcf, vcf_idx: vcf_idx) }.join(
                    out_GATK_VARIANT_RECALIBRATOR.map { dbname, recal, tranches ->
                        record(dbname: dbname, recal: recal, tranches: tranches)
                    },
                    by: "dbname"
                ).map { rec -> tuple(rec.dbname, rec.vcf, rec.vcf_idx, rec.recal, rec.tranches) },
                params.vqsr_mode,
                paths.parasite.fasta,
            )
        }
    }
}
