#!/bin/bash

# create a new directory for RNA_map
cd /projects/compbio/users/xran2/tool/RNA_map

#-------------------------------------------#
#                                           #
#      Method 1: use conda to install       #
#                                           #
#-------------------------------------------#

# Check current conda env
conda env list

# create a new env called rna_map
conda create --name rna_map
conda activate rna_map
conda install -c bioconda hisat2 subread 
conda install -c bioconda samtools

# Check for installation 
hisat2 --version
featureCounts -v
samtools -v

#-------------------------------------------#
#                                           #
#      Method 2: install local              #
#                                           #
#-------------------------------------------#

# Download HISAT2 
wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
unzip download
#add /projects/compbio/users/xran2/tool/RNA_map/hisat2-2.2.1 to PATH
export PATH=$PATH:/projects/compbio/users/xran2/tool/RNA_map/hisat2-2.2.1
hisat2 --version

# Download featureCounts
wget https://github.com/ShiLab-Bioinformatics/subread/releases/download/2.0.2/subread-2.0.2-Linux-x86_64.tar.gz
tar -xzvf subread-2.0.2-Linux-x86_64.tar.gz
export PATH=$PATH:/projects/compbio/users/xran2/tool/RNA_map/subread-2.0.2-Linux-x86_64/bin
featureCounts -v

# Download samtools
wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
tar -xjvf samtools-1.20.tar.bz2
export PATH=$PATH:/projects/compbio/users/xran2/tool/RNA_map/samtools-1.20
samtools 

#-------------------------------------------#
#                                           #
#      Part 2: Make reference               #
#                                           #
#-------------------------------------------#

# Download Human reference data
# Human reference (GRCh38) 
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xzvf grch38_genome.tar.gz
export PATH=$PATH:/projects/compbio/users/xran2/tool/RNA_map/hisat2-2.2.1 # add hisat2 to PATH
cd grch38
./make_grch38.sh
cd ..

# Download Human reference gtf file
# Version 106
wget https://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz  # Download GTF file version 106
gunzip Homo_sapiens.GRCh38.106.gtf.gz

# # Version 113
# wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz # Download GTF file version 113
# gunzip Homo_sapiens.GRCh38.113.gtf.gz

# Create Gene ID and Gene Name Mapping file
awk '$3 == "gene" {print $10 "\t" $14}' Homo_sapiens.GRCh38.106.gtf | sed 's/"//g' | sed 's/;//g' | sort > gene_id_name_map_sorted.txt
# awk '$3 == "gene" {print $10 "\t" $14}' Homo_sapiens.GRCh38.113.gtf | sed 's/"//g' | sed 's/;//g' | sort > gene_id_name_map_sorted.txt

