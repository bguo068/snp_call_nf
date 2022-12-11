#! /usr/bin/env bash

# Download the ENA tool jar file
# This is the best tool I found to download files via ENA run accession number.
# 1. It will figure out the url and number of fastq files per run
# 2. It also generate logs to verify which files are downloaded correctly
# 3. It will repspect files already downloaded when re-run.

if [ ! -f ena-file-downloader.jar ]; then
    wget http://ftp.ebi.ac.uk/pub/databases/ena/tools/ena-file-downloader.zip
    unzip ena-file-downloader.zip
    rm ena-file-downloader.zip run.bat run.sh 
fi

# Download fastq files (to download multiple runs, --acesssions=A,B,C,D)
# 
# For testing purposes, we will download the following samples/runs from MalariaGen Pf6
# 1. a sample failed QC:                ERR1099208	        FP0026-C	WAF	    0.45	False
# 2. a sample with multiple runs:       ERR018904,ERR015361	PA0035-C	WAF	    89.9	True
# 3. samples from different regions:    ERR676507	        PE0309-C	EAF	    89.87	True
#                                       ERR404155	        QG0016-C	CAF	    88.65	True
#                                       ERR1172588	        QJ0006-C	WAF	    89.43	True
#                                       ERR036596	        PR0001-CW	SAS	    89.62	True
#                                       ERR015419	        PH0057-C	ESEA	88.91	True 
#                                       ERR019041	        PD0009-01	WSEA	77.01	True
#                                       ERR022856	        PP0010-C	SAM	    88.92	True
#                                       ERR018922	        PP0006-C	Lab	    87.87	False	Suspected_lab_strain
# MalariaGen Pf6 meta: 
#   wget ftp://ngs.sanger.ac.uk/production/malaria/pfcommunityproject/Pf6/Pf_6_samples.txt


java -jar ena-file-downloader.jar --format=READS_FASTQ --protocol=FTP --asperaLocation=null \
    --location=$PWD \
    --accessions=ERR1099208,ERR018904,ERR015361,ERR676507,ERR404155,ERR1172588,ERR036596,ERR015419,ERR019041,ERR022856,ERR018922