nextflow.enable.dsl = 2

process R_MOIMIX_FWS {
    publishDir "result"
    conda "$projectDir/modules/moimix_env.yaml"
    input: 
    path(vcf)
    output:
    path("*.fws.txt")
    script:
    """#! /usr/bin/env Rscript

    dir.create("$projectDir/Rlib", showWarnings = FALSE) 
    .libPaths("$projectDir/Rlib")

    if( !require("moimix") ){
        BiocManager::install("bahlolab/moimix")
    }

    library(SeqArray)
    library("moimix")

    seqVCF2GDS("$vcf", "snps.gds")
    isolates <- seqOpen("snps.gds")
    sample.id <- seqGetData(isolates, "sample.id")
    fws_all <- getFws(isolates)
    df = data.frame(sample.id, fws_all)
    write.table(df, "out.fws.txt", sep="\t", row.names=FALSE, quote=FALSE)

    sessionInfo()
    """
}

