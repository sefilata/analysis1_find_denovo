#!/bin/bash

config=../config.sh
source ${config}

cd ${base_dir}
mkdir -p ${analysis_dir}
mkdir -p ${log_dir}
make -C ${code_dir}

# TRGT
mkdir -p ${trgt_dir}
bash ${trgt_bash} ${config}

# TRGT-denovo
mkdir -p ${trgt_denovo_dir}
bash ${trgt_denovo_bash} ${config}

# Validate and count DNMs (classified by unit length)
mkdir -p ${validate_count_dir1}
bash ${validate_count_bash} ${config} "03_validate_count"

# Validate and count DNMs (classified by GC content)
mkdir -p ${validate_count_dir2}
bash ${validate_count_bash} ${config} "03-2_count_GC"


