#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
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

# execution
cd ${trgt_denovo_dir}
for ((i = 0; i < n; i++)); do
	sample=${sample_list[$i]}
	father=${father_list[$i]}
	mother=${mother_list[$i]}
	if printf '%s\n' "${skips_list[@]}" | grep -qx "${sample}"; then
		echo "Skipping sample: ${sample}" >> ${trgt_denovo_log} 2>&1
		continue
	fi

	mkdir -p ${sample}/
	cd ${sample}/
	echo "Started processing sample: $sample, Father: $father, Mother: $mother at $(date)" >> ${trgt_denovo_log} 2>&1
	time $trgtdenovo_exec trio \
		-@ ${trgt_denovo_parallel} \
		--reference ${reference_genome} \
		--bed ${repeat_list} \
		--father ${trgt_dir}/${father}/${father} \
		--mother ${trgt_dir}/${mother}/${mother} \
		--child ${trgt_dir}/${sample}/${sample} \
		--out ${sample}_denovo_TRs.tsv >> ${trgt_denovo_log} 2>&1
	cd ../
done

echo "Completed processing all samples at $(date)"
