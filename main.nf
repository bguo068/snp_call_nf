#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

def print_parameters() {
    println("""======================== PARAMETERS ==========================
fq_map:\t${file(params.fq_map)}
known_variants:\t${params.known_sites}
split:\t${params.split}
filtering:
\thard:\t ${params.hard}
\tvqsr:\t ${params.vqsr}
coverage_only:\t ${params.coverage_only}
outdir:\t${params.outdir}
tmpdir:\t${params.gatk_tmpdir}
gvcf_only:\t${params.gvcf_only}
=============================================================""")
}

process BOWTIE2_ALIGN_TO_HOST {
    tag "$meta.Sample~$meta.Run"
    publishDir "$params.outdir/flagstat_host", pattern: "*.flagstat"
    publishDir "$params.outdir/read_length",   pattern: "*_read_length.txt"
    input:
    tuple val(meta), path(fastq)
    output:
    tuple val(meta), path("*.bam"),             emit: bam
    tuple val(meta), path("*.flagstat"),        emit: flagstat
    tuple val(meta), path("*_read_length.txt"), emit: read_length
    shell:
    def reads = meta.is_paired == 1 ? "-1 ${fastq[0]} -2 ${fastq[1]}" : "-U $fastq"
    """
    bowtie2 -x $meta.Host_ref_prefix  $reads | \
        samtools view -q 0 -bS > ${meta.Sample}~${meta.Run}~to_host.bam
    samtools flagstat ${meta.Sample}~${meta.Run}~to_host.bam \
        > ${meta.Sample}~${meta.Run}~to_host.bam.flagstat
    samtools view ${meta.Sample}~${meta.Run}~to_host.bam | \
        head -100000 | awk '{a+=length(\$10)}END{print (a/NR)}' \
        > ${meta.Sample}~${meta.Run}_read_length.txt

    """
    stub:
    """ 
    touch ${meta.Sample}~${meta.Run}~to_host.bam{,.flagstat}
    touch ${meta.Sample}~${meta.Run}_read_length.txt
    """
}

process SAMTOOLS_VIEW_RM_HOST_READS {
    tag "$meta.Sample~$meta.Run"
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path("*.bam")
    shell:
    def flag = meta.is_paired ? "-f 12 -F 256" : "-f 4"
    """
    samtools view -b $flag $bam | \
        samtools sort -n -o ${meta.Sample}~${meta.Run}~unmapped_sorted.bam
    """
    stub:
    """ touch ${meta.Sample}~${meta.Run}~unmapped_sorted.bam """
}

process SAMTOOLS_FASTQ {
    tag "$meta.Sample~$meta.Run"
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path("*.fastq.gz")
    shell:
    if(meta.is_paired)
    """
    samtools fastq -@ $task.cpus $bam \
        -1 f1.fastq.gz -2 f2.fastq.gz -0 /dev/null -s /dev/null -n	
    """
    else
    """
    samtools fastq -@ $task.cpus $bam \
        -o /dev/null -s /dev/null -0 f.fastq.gz -n  
    """
    stub:
    if(meta.is_paired) """ touch f1.fastq.gz f2.fastq.gz """
    else               """ touch f.fastq.gz """
}

process BOWTIE2_ALIGN_TO_PARASITE {
    tag "$meta.Sample~$meta.Run"
    publishDir "$params.outdir/flagstat_parasite", pattern: "*.flagstat"
    input:
    tuple val(meta), path(fastq)
    val (parasite_ref_prefix)
    output:
    tuple val(meta), path("*.bam"),       emit: bam
    tuple val(meta), path("*.flagstat"),  emit: flagstat
    shell:
    def reads = meta.is_paired == 1 ? "-1 ${fastq[0]} -2 ${fastq[1]}" : "-U $fastq"
    def read_group = "--rg-id ${meta.Sample} --rg SM:${meta.Sample} --rg PL:Illumina"
    """
    bowtie2 -x $parasite_ref_prefix $reads $read_group | \
        samtools view -q 0 -bS > ${meta.Sample}~${meta.Run}~to_parasite.bam
    samtools flagstat ${meta.Sample}~${meta.Run}~to_parasite.bam \
        > ${meta.Sample}~${meta.Run}~to_parasite.bam.flagstat
    """
    stub:
    """ touch ${meta.Sample}~${meta.Run}~to_parasite.bam{,.flagstat} """
}

process PICARD_MERGE_SORT_BAMS{
    tag "$sample"
    input:
    tuple val(sample), path(bams)
    output:
    tuple val(sample), path("merged_sorted.bam")
    shell:
    def to_merge =  bams.size() > 1
    def in_bams = bams.join(" -I ")
    def merge_opts = "-VALIDATION_STRINGENCY LENIENT -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
    def sort_opts = "-VALIDATION_STRINGENCY LENIENT -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
    """
    if $to_merge; then
        gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.giga}G" \
	MergeSamFiles -I $in_bams -O merged.bam $merge_opts
    else
       ln -s $bams merged.bam
    fi

    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.giga}G" \
        SortSam -I merged.bam -O merged_sorted.bam --SORT_ORDER coordinate $sort_opts
    rm merged.bam
    """
    stub:
    """ touch merged_sorted.bam """
}

process PICARD_MARK_DUPLICATES {
    tag "$sample"
    input:
    tuple val(sample), path(bam)
    output:
    tuple val(sample), path("*dedup.bam")
    shell:
    def opts = "--USE_JDK_DEFLATER true --USE_JDK_INFLATER true"
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.giga}G" \
    MarkDuplicates \
	--REMOVE_DUPLICATES $opts \
        -I $bam -O ${sample}_dedup.bam --METRICS_FILE dedup_metrics.txt
    """
    stub:
    """ touch ${sample}_dedup.bam """
}

process GATK_BASE_RECALIBRATOR {
    tag "$sample"
    input:
    tuple val(sample), path(bam)
    val(ref)
    val(known_sites)
    output:
    tuple val(sample), path(bam), path("recal_data.table")
    shell:
    def known_sites_str = known_sites.join(" --known_sites ")
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir}  -Xmx${task.memory.giga}G" \
        BaseRecalibrator -I $bam -O recal_data.table \
	-R ${ref} --known-sites ${known_sites_str} 
    """
    stub: 
    """ touch recal_data.table """
}

process GATK_APPLY_BQSR {
    tag "$sample"
    publishDir "$params.outdir/recalibrated"
    input:
    tuple val(sample), path(bam), path(recal_table)
    val(ref)
    output: 
    tuple val(sample), path("*recalibrated.bam")
    shell:
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.giga}G" \
        ApplyBQSR -I $bam -O ${sample}_recalibrated.bam \
        -R ${ref} --bqsr-recal-file $recal_table
    """
    stub: 
    """ touch ${sample}_recalibrated.bam """
}

process BEDTOOLS_GENOMECOV {
    tag "$sample"
    publishDir "$params.outdir/coverage"
    input:
    tuple val(sample), path(bam)
    val(ref)
    output:
    path("*recalibrated.coverage.BedGraph.gz")
    path("*.coverage.summary.txt")
    shell:
    """
    bedtools genomecov -bg -ibam $bam -g ${ref} | \
        gzip -v > ${sample}_recalibrated.coverage.BedGraph.gz

    # summarize the coverage
    zcat ${sample}_recalibrated.coverage.BedGraph.gz | \
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
    tag "$sample"
    publishDir "$params.outdir/flagstat"
    input:
    tuple val(sample), path(bam)
    output:
    path("*.flagstat")
    shell: 
    """
    samtools flagstat $bam > ${sample}.flagstat
    """
    stub:
    """ touch ${sample}.flagstat """
}

process GATK_HAPLOTYPE_CALLER {
    tag "$sample"
    publishDir "$params.outdir/gvcf"
    input: 
    tuple val(sample), path(bam)
    val(ref)
    output:
    tuple val(sample), path("*.g.vcf"), path("*.g.vcf.idx")
    shell:
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.giga}G" \
        HaplotypeCaller \
	-I $bam -O ${sample}.g.vcf -R ${ref} -ERC GVCF
    """
    stub:
    """ touch ${sample}.g.vcf{,.idx} """
}

process GATK_GENOMICS_DB_IMPORT {
    tag "$interval"
    input:
    val(interval)
    path(gvcf_map)
    output:
    path("*", type: 'dir', maxDepth: 1)
    shell:
    def dbname = interval.replaceAll(":", "~")
    def mem = Math.round(task.memory.giga * 0.75) // the rest of memory for c/c++ library
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${mem}G" \
        GenomicsDBImport \
        --batch-size 100 --reader-threads 5 --consolidate true \
        --sample-name-map ${gvcf_map} \
        --genomicsdb-workspace-path $dbname \
        -L $interval
    """
    stub:
    def dbname = interval.replaceAll(":", "~")
    """ mkdir $dbname """
}

process GATK_GENOTYPE_GVCFS {
    tag "${db.getName()}"
    input:
    path(db)
    val(ref)
    output:
    tuple env(DBNAME), path("*.vcf"), path("*.idx")
    shell:
    def dbname = db.getName()
    def mem = Math.round(task.memory.giga * 0.75) // the rest of memory for TileDB library

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
    tag "$dbname"
    input:
    tuple val(dbname), path(vcf), path(idx)
    val(ref)
    output:
    tuple val(dbname), path("*.snp.vcf"), path("*.snp.vcf.idx")
    shell:
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.giga}G" \
        SelectVariants --select-type-to-include SNP \
        -R $ref -V $vcf -O ${dbname}.snp.vcf
    """
    stub:
    """touch ${dbname}.snp.vcf{,.idx} """
}

process GATK_VARIANT_FILTRATION {
    tag "$dbname"
    publishDir "$params.outdir/hardfilt"
    input:
    tuple val(dbname), path(vcf), path(idx)
    val(ref)
    output:
    tuple val(dbname), path("*.snp.hardfilt.vcf"), path("*.snp.hardfilt.vcf.idx")
    shell:
    def filter_str = ''
    params.hard_filters.each {filter_str+= " -filter \"${it.filter}\" --filter-name \"${it.name}\" "}
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.giga}G" VariantFiltration \
    ${filter_str} \
	-R ${ref} \
	-V ${vcf} \
	-O ${dbname}.snp.hardfilt.vcf
    """
    stub:
    """ touch ${dbname}.snp.hardfilt.vcf{,.idx}"""
}

def get_vqsr_resources (resources) {
    def str = ''
    resources.each {it ->
        def type_str = ''
        if(it.type == "truth") {
            type_str = "known=false,training=true,truth=true"
        }
        else if (it.type== "training") {
            type_str = "known=false,training=true,truth=false"
        }
        else if (it.type== "known") {
            type_str = "known=true,training=false,truth=false"
        }
        else {
            print "vqsr resources type can only be one of known, training, known"
            System.exit(-10)
        } 
        def vcf_file = file(it.vcf)
    
        str += " --resource:$it.name,$type_str,prior=$it.prior $vcf_file " 
    }
    return str
}

process GATK_VARIANT_RECALIBRATOR {
    tag "$dbname"
    input:
    tuple val(dbname), path(vcf), path(idx)
    val(resources)
    val(opts)
    val(mode)
    val(ref)
    output:
    tuple val(dbname), path("*.recal.vcf"), path("*.tranches")
    shell:
    def resources_str = get_vqsr_resources(resources)
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.giga}G" VariantRecalibrator \
        -R $ref -V $vcf $resources_str $opts -mode $mode \
        -O ${dbname}.recal.vcf --tranches-file ${dbname}.all.tranches
    """
    stub:
    """ touch ${dbname}.recal.vcf  ${dbname}.all.tranches """
}
process GATK_APPLY_VQSR {
    publishDir "$params.outdir/vqsrfilt"
    tag "$dbname"
    input:
    tuple val(dbname),  path(vcf), path(vcf_idx), path(recal), path(tranches)
    val(mode)
    val(ref)
    output:
    tuple val(dbname), path("*.vqsrfilt.*.vcf")
    shell:
    """
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.giga}G" IndexFeatureFile -I $recal
    gatk --java-options "-Djava.io.tmpdir=${params.gatk_tmpdir} -Xmx${task.memory.giga}G" ApplyVQSR \
        -R $ref -V $vcf --recal-file $recal --tranches-file $tranches -mode $mode \
        --output ${dbname}.vqsrfilt.${mode}.vcf
    """
    stub:
    """ touch ${dbname}.vqsrfilt.${mode}.vcf """
}

workflow {

    // print key parameters
    print_parameters()

    // Get paths
    def paths = [:]
    paths.parasite = [fasta : file(params.parasite.fasta),  \
            fasta_prefix: file(params.parasite.fasta_prefix) ]
    paths.host = [fasta : [], fasta_prefix: [] ]
    params.host.fasta.each {it -> paths.host.fasta.add(file(it))}
    params.host.fasta_prefix.each {it -> paths.host.fasta_prefix.add(file(it))}
    paths.known_sites = []
    params.known_sites.each {it -> paths.known_sites.add(file(it))}

    // prepare tmpdir for gatk
    def tmpdir = file(params.gatk_tmpdir)
    if (!tmpdir.exists()) {
        tmpdir.mkdirs()
    }

    // Prepare input chanel
    input_ch = channel.fromPath(params.fq_map) \
        | splitCsv(skip:1,  sep: '\t') \
        | groupTuple (by: [0, 1, 2], sort: true) \
        | map {
            sample, host_id, run, mate_id, fq ->
            assert (fq.size() == 1) || (fq.size() == 2): \
                "number of fq files for each run can only be 1 or 2! Found ${fq.size()} for ${sample}/${run}\n"
            def meta = [:]
            meta.Sample = sample
            meta.Run = run
            meta.Host_ref_prefix = paths.host.fasta_prefix[host_id.toInteger()]
            meta.is_paired = fq.size() == 2 ? 1:0
            def fq_paths = []
            fq.each(it->fq_paths.add(file(it)))
            return [meta, fq_paths] 
        } 

    // Remove reads mapped to host and split unmapped to fastq files
    input_ch |  BOWTIE2_ALIGN_TO_HOST
    BOWTIE2_ALIGN_TO_HOST.out.bam  | SAMTOOLS_VIEW_RM_HOST_READS | SAMTOOLS_FASTQ

    // Align to parasite genome
    BOWTIE2_ALIGN_TO_PARASITE( SAMTOOLS_FASTQ.out, paths.parasite.fasta_prefix) 
    
    // Non-blocking grouped gathering
    n_run = channel.fromPath(params.fq_map) | splitCsv(skip:1,  sep: '\t') \
        | map {sample, host_id, run, mate_id, fq -> [sample, run]} \
        | unique | groupTuple | map {sample, runs -> [sample, runs.size()]} 

    merge_input = BOWTIE2_ALIGN_TO_PARASITE.out.bam \
        | map {meta, bam->[meta.Sample, bam]} \
        | combine(n_run, by: 0) \
        | map {sample, bam, sz -> [groupKey(sample, sz), bam] } \
        | groupTuple 
    
    // For each sample, merge and sort bam files of different runs 
    merge_input | PICARD_MERGE_SORT_BAMS | PICARD_MARK_DUPLICATES

    // Base recalibration
    GATK_BASE_RECALIBRATOR(PICARD_MARK_DUPLICATES.out, paths.parasite.fasta, paths.known_sites )
    GATK_APPLY_BQSR (GATK_BASE_RECALIBRATOR.out, paths.parasite.fasta)

    // Generate stat files
    BEDTOOLS_GENOMECOV(GATK_APPLY_BQSR.out, paths.parasite.fasta)
    SAMTOOLS_FLAGSTAT(GATK_APPLY_BQSR.out)

    // stop before GATK_HAPLOTYPE_CALLER if only coverage inforation is need
    if (params.coverage_only) {
        ch_bqsr_bam = Channel.empty()
    }
    else {
        ch_bqsr_bam = GATK_APPLY_BQSR.out
    }

    // Generate gvcf
    GATK_HAPLOTYPE_CALLER(ch_bqsr_bam, paths.parasite.fasta)

    if (!params.gvcf_only)
    {
        // Collect information to make gvcf_map file
        gvcf_map_ch = GATK_HAPLOTYPE_CALLER.out \
            | map{sample, gvcf, idx -> "$sample\t$gvcf"} \
            | collectFile(name: "gvcf_map.txt", newLine: true, sort: true) // sort by nautral ordering instead of hash on content
            | first

        // Import gvcf files to genomicsdb
        interval_ch = channel.fromList(params.genome_intervals[params.split])
        GATK_GENOMICS_DB_IMPORT(interval_ch, gvcf_map_ch) 

        // Genotype gvcf genomics db
        GATK_GENOTYPE_GVCFS(GATK_GENOMICS_DB_IMPORT.out, paths.parasite.fasta) 

        // Select SNP only
        GATK_SELECT_VARIANTS(GATK_GENOTYPE_GVCFS.out, paths.parasite.fasta) 

        // Hard filtering and Keep 'PASS' variants
        if (params.hard == true) {
            GATK_VARIANT_FILTRATION(GATK_SELECT_VARIANTS.out, paths.parasite.fasta) 
        }

        // VQSR variant filtering
        if (params.vqsr == true) {
            GATK_VARIANT_RECALIBRATOR( 
                GATK_SELECT_VARIANTS.out, 
                params.vqsr_resources, 
                params.vqsr_opts, 
                params.vqsr_mode, 
                paths.parasite.fasta)
            GATK_APPLY_VQSR(
                GATK_SELECT_VARIANTS.out.combine(GATK_VARIANT_RECALIBRATOR.out, by: 0),
                params.vqsr_mode, 
                paths.parasite.fasta)
        }
    }
}
