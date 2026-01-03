#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=all

source ${1:-../config.sh}

# samples
mapfile -t sample_list < <(awk -F'\t' -v i=1 '{print $i}' "$trio_list")
mapfile -t father_list < <(awk -F'\t' -v i=2 '{print $i}' "$trio_list")
mapfile -t mother_list < <(awk -F'\t' -v i=3 '{print $i}' "$trio_list")
skips_list=()
if [[ -f "$skip_trio_list" ]]; then
	mapfile -t skips_list < <(awk -F'\t' -v i=1 '{print $i}' "$skip_trio_list")
fi
n=${#sample_list[@]}


# count settings
# exec_mode="03_validate_count"
# exec_mode="03-2_count_GC"
exec_mode=${2:-"03_validate_count"}
if [ "$exec_mode" == "03_validate_count" ]; then
	count_exec=${validate_count_cpp1}
	out_dir=${validate_count_dir1}
elif [ "$exec_mode" == "03-2_count_GC" ]; then
	count_exec=${validate_count_cpp2}
	out_dir=${validate_count_dir2}
else
	echo "[ERROR] Unknown exec_mode: ${exec_mode}" 2>&1 | tee -a ${validate_count_log}
	exit 1
fi

# execution
cd ${exec_idf}/
for ((i = 0; i < n; i++)); do
	if [[ " ${skips[*]} " =~ " $i " ]]; then
		continue
	fi

	sample=${sample_list[$i]}
	father=${father_list[$i]}
	mother=${mother_list[$i]}

	mkdir -p ${sample}
	cd ${sample}

	echo "Started processing sample: ${sample}, Father: ${father}, Mother: ${mother} at $(date)" >> ${validate_count_log} 2>&1

	father_vcf="${trgt_dir}/${father}/${father}.sorted.vcf"
	mother_vcf="${trgt_dir}/${mother}/${mother}.sorted.vcf"
	child_vcf="${trgt_dir}/${sample}/${sample}.sorted.vcf"
	trgt_denovo="${trgt_denovo_dir}/${sample}/${sample}_denovo_TRs.tsv"

	# decompress
	if [ ! -f "$father_vcf" ]; then
		gunzip -c ${father_vcf}.gz > ${father_vcf}
	fi
	if [ ! -f "$mother_vcf" ]; then
		gunzip -c ${mother_vcf}.gz > ${mother_vcf}
	fi
	if [ ! -f "$child_vcf" ]; then
		gunzip -c ${child_vcf}.gz > ${child_vcf}
	fi

	/usr/bin/time -v ${count_exec} ${trgt_denovo} ${child_vcf} ${father_vcf} ${mother_vcf} >> ${trgt_denovo_log} 2>&1

	awk -F'\t' 'NR==1 {print; next} {for(i=34;i<=NF;i++) if($i!=0) {print; next}}' validate_count.tsv > filter_all.tsv
	awk -F'\t' '{if(NR==1 || $33==1) print} END {print $0}' validate_count_filter_all.tsv > filter_denovo.tsv
	
	cd ../
done

echo "Completed processing all samples with mode ${exec_mode} at $(date)" >> ${validate_count_log} 2>&1

