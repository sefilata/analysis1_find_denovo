#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=all

source ${1:-../config.sh}

# samples
mapfile -t sample_list < <(awk -F'\t' -v i=1 '{print $i}' "$sample_karyotype_list")
mapfile -t karyotype_list < <(awk -F'\t' -v i=2 '{print $i}' "$sample_karyotype_list")
skips_list=()
if [[ -f "$skip_sample_list" ]]; then
	mapfile -t skips_list < <(awk -F'\t' -v i=1 '{print $i}' "$skip_sample_list")
fi
n=${#sample_list[@]}

# execution
cd ${trgt_dir}
for ((i = 0; i < n; i++)); do
	sample=${sample_list[$i]}
	karyotype=${karyotype_list[$i]}
	if printf '%s\n' "${skips_list[@]}" | grep -qx "${sample}"; then
		echo "Skipping sample: ${sample}" >> ${trgt_log} 2>&1
		continue
	fi

	mkdir -p ${sample}/
	cd ${sample}/
	echo "Started processing sample: ${sample}, Karyotype: ${karyotype} at $(date)" >> ${trgt_log} 2>&1

	reads=${data_path_prefix}${sample}${data_path_suffix}
	time $trgt genotype \
		--genome ${reference_genome} \
		--threads ${trgt_parallel} \
		--repeats ${repeat_list} \
		--reads ${reads} \
		--output-prefix ${sample} \
		--karyotype ${karyotype} >> ${trgt_log} 2>&1

	bcftools sort -m 3072M -Ob -o ${sample}.sorted.vcf.gz ${sample}.vcf.gz >> ${trgt_log} 2>&1
	bcftools index --threads 4 ${sample}.sorted.vcf.gz >> ${trgt_log} 2>&1
	samtools sort -@ 32 -o ${sample}.spanning.sorted.bam ${sample}.spanning.bam >> ${trgt_log} 2>&1
	samtools index -@ 32 ${sample}.spanning.sorted.bam >> ${trgt_log} 2>&1
	cd ../
done

echo "Completed processing all samples at $(date)"
