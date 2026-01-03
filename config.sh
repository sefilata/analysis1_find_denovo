#!/bin/bash

# Directories
base_dir=./			# CHANGE THIS LINE
code_dir=${base_dir}/code
data_dir=${base_dir}/data
analysis_dir=${base_dir}/analysis
log_dir=${base_dir}/log

# General
sample_karyotype_list=${base_dir}/samples.tsv
skip_sample_list=${base_dir}/skip_samples.tsv
trio_list=${base_dir}/trios.tsv
skip_trio_list=${base_dir}/skip_trios.tsv
reference_genome=/path/to/reference.fa.gz			# CHANGE THIS LINE
repeat_list=/path/to/repeats.bed.gz			# CHANGE THIS LINE
# hifi mapped bam paths are supposed to be constructed as: 
# 	${data_path_prefix}${sample}${data_path_suffix}
data_path_prefix=${data_dir}/hifi/mapped/CHM13/			# CHANGE THIS LINE
data_path_suffix=.CHM13.haplotagged.bam			# CHANGE THIS LINE

# Tool paths
bcftools=/path/to/bcftools			# CHANGE THIS LINE
samtools=/path/to/samtools			# CHANGE THIS LINE
trgt=/path/to/trgt			# CHANGE THIS LINE
trgt_denovo=/path/to/trgt_denovo			# CHANGE THIS LINE

# TRGT
trgt_parallel=32
trgt_bash=${code_dir}/01_trgt.sh
trgt_log=${log_dir}/01_trgt.log
trgt_dir=${analysis_dir}/01_trgt

# TRGT-denovo
trgt_denovo_parallel=32
trgt_denovo_bash=${code_dir}/02_trgtdenovo.sh
trgt_denovo_log=${log_dir}/02_trgtdenovo.log
trgt_denovo_dir=${analysis_dir}/02_trgtdenovo

# Validate and count DNMs
validate_count_cpp1=${code_dir}/exec/03_validate_count.out
validate_count_cpp2=${code_dir}/exec/03-2_validate_count_GC.out
validate_count_bash=${code_dir}/03_validate_count.sh
validate_count_log=${log_dir}/03_validate_count.log
validate_count_dir1=${analysis_dir}/03_validate_count
validate_count_dir2=${analysis_dir}/03_validate_count_GC
