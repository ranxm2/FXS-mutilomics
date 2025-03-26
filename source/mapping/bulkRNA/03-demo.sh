#!/bin/bash


# Path to the tools of You Inited
tool_path="/projects/compbio/users/xran2/tool/RNA_map"

# Path to the FASTQ files
fastq_path="/projects/compbio/users/xran2/wen/22q/02-bulkRNA-synaptosome/data"

# Path to the result directory
result_path="/projects/compbio/users/xran2/wen/22q/02-bulkRNA-synaptosome/result"


#-------------------------------------------#
#                                           #
# Method 1: Load Sample information from CSV#
#                                           #
#-------------------------------------------#

SLURM_ARRAY_TASK_ID=1
echo "Processing File index: ${SLURM_ARRAY_TASK_ID}"
SAMPLE_NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ../code/file_list.csv | tr -d ',')
echo "Sample name: $SAMPLE_NAME"


#-------------------------------------------#
#                                           #
# Method 2: Given Sample information        #
#                                           #
#-------------------------------------------#

# SAMPLE_NAME = "Sample_1"

# ADD tool pathway into the PATH
export PATH=$PATH:${tool_path}/hisat2-2.2.1 # 
export PATH=$PATH:${tool_path}/subread-2.0.2-Linux-x86_64/bin 
export PATH=$PATH:${tool_path}/samtools-1.20

# ADD reference genome pathway
HISAT2_INDEX="${tool_path}/grch38/genome"
featureCounts_INDEX="${tool_path}/grch38/Homo_sapiens.GRCh38.106.gtf"
Gene_name_map_Index="${tool_path}/gene_id_name_map_sorted.txt"
BAM_DIR="01-BAM"
Feature_DIR="02-Feature"
Count_DIR="03-Counts"
mkdir -p ${BAM_DIR} ${Feature_DIR} ${Count_DIR}

echo "Processing: ${SAMPLE_NAME}"

echo " " 
echo "-------------------------------------------------"
echo "------- Step 1: HISAT2 Mapping Mutil Mapping ----"
echo "-------------------------------------------------"
echo "Step 1 start at: $(date)"
start_time=$(date +%s)


hisat2 -q --rna-strandness R -x ${HISAT2_INDEX} \
    -1 ${fastq_path}/${SAMPLE_NAME}_R1_001.fastq.gz \
    -2 ${fastq_path}/${SAMPLE_NAME}_R2_001.fastq.gz  | \
    samtools sort -o ${BAM_DIR}/${SAMPLE_NAME}.bam
                                
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 1 finish at: $(date)"
echo "Step 1 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."



echo " " 
echo "--------------------------------------------------------------"
echo "------- Step 2: Use featureCounts to generate gene counts ----"
echo "--------------------------------------------------------------"
echo "Step 2 start at: $(date)"
start_time=$(date +%s)


featureCounts -p -S 2 -a ${featureCounts_INDEX} \
              -o ${Feature_DIR}/${SAMPLE_NAME}_feature.txt ${BAM_DIR}/${SAMPLE_NAME}.bam

                                
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 2 finish at: $(date)"
echo "Step 2 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."




echo " " 
echo "--------------------------------------------------------------"
echo "------- Step 3: Exctract gene counts and sort ---------------"
echo "--------------------------------------------------------------"
echo "Step 3 start at: $(date)"
start_time=$(date +%s)

tail -n +3 ${Feature_DIR}/${SAMPLE_NAME}_feature.txt| cut -f 1,7 | sort > ${Feature_DIR}/${SAMPLE_NAME}_all_gene_counts_sorted.txt
                                
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 3 finish at: $(date)"
echo "Step 3 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."



echo " " 
echo "--------------------------------------------------------------"
echo "------- Step 4: Add gene names to the count file ------------"
echo "--------------------------------------------------------------"
echo "Step 4 start at: $(date)"
start_time=$(date +%s)

join -1 1 -2 1 ${Feature_DIR}/${SAMPLE_NAME}_all_gene_counts_sorted.txt ${Gene_name_map_Index} > ${Feature_DIR}/${SAMPLE_NAME}_all_gene_with_names.txt

                                
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 4 finish at: $(date)"
echo "Step 4 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."




echo " " 
echo "--------------------------------------------------------------"
echo "------- Step 5: Convert to CSV format -----------------------"
echo "--------------------------------------------------------------"
echo "Step 5 start at: $(date)"
start_time=$(date +%s)

echo -e "GeneID,Count,GeneName" > ${Count_DIR}/${SAMPLE_NAME}_final_gene_with_names.csv
awk '{print $1 "," $2 "," $3}' ${Feature_DIR}/${SAMPLE_NAME}_all_gene_with_names.txt >> ${Count_DIR}/${SAMPLE_NAME}_final_gene_with_names.csv
                                
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Step 5 finish at: $(date)"
echo "Step 5 used:  $hours hour(s), $minutes minute(s), and $seconds second(s)."







