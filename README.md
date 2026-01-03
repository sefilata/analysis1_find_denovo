# About this pipeline

This pipeline estimates de novo duplication and contraction mutation rates per repeat unit using family-based data (trios).  
The input data are assumed to be PacBio HiFi reads.

The original analysis was designed to detect de novo mutations in CEPH1463 pedigree data, following the methodology described in Porubsky et al. (2025). This pipeline is not limited to CEPH1463 pedigree data and can be applied to other trio- or family-based high-fidelity long-read sequencing dataset.

---

# Usage

## Requirements

This pipeline is intended to run in a Linux environment and requires the following software:

- Bash
- GCC (with C++17 support)
- bcftools
- samtools
- TRGT
- TRGT-denovo

---

## Inputs

### Sequencing data

Mapped HiFi reads are required and must be placed under the `data/` directory.  
The directory structure must follow the patterns specified by the variables `data_path_prefix` and `data_path_suffix` in `config.sh`.

### Reference files

- A reference genome (FASTA + index)
- A BED file specifying tandem repeat regions

### Sample and pedigree information

- **`samples.tsv`**  
  A tab-separated file without a header:
  - Row 1: sample names  
  - Row 2: karyotypes (`XX` or `XY`)  
  - Row 3 and beyond: optional comments (ignored by the pipeline)

- **`trios.tsv`**  
  A tab-separated file defining family relationships:
  - Row 1: child samples  
  - Row 2: father samples  
  - Row 3: mother samples  
  - Row 4 and beyond: optional comments (ignored by the pipeline)

- **Skip file (optional)**  
  A list of samples to exclude from the analysis.

---

## Running the pipeline

1. Edit the required variables in `config.sh`.
2. Run the pipeline script:

```bash
bash pipeline.sh
```

---

## Outputs

- **`filter_denovo.tsv`**  
  A list of de novo duplication and contraction mutations detected in the target tandem repeats.

- Rows after the `is_denovo` field represent unit-based annotations, such as repeat unit length or GC content.
- The final row, labeled `Statistics`, summarizes the total number of repeat units across all target TRs, which serves as the denominator for calculating per-unit mutation rates.
- `result_only.tsv` contains the same mutation list but excludes filtering-related information.

- **`filter_all.tsv`**  
  A list of all target tandem repeats, together with the counts of repeat units for each unit length or GC content.


# Example
For your information, we used the following environment and data.
```bash
# download sequence data
aws s3 --no-sign-request cp s3://platinum-pedigree-data/ ./data/ --recursive

# download reference data
aws s3 --no-sign-request cp s3://human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz /path/to/reference.fa.gz
samtools faidx /path/to/reference.fa.gz

# download the list of tandem repeats
wget https://zenodo.org/records/13178746/files/chm13v2.0_maskedY_rCRS.platinumTRs-v1.0.trgt.annotations.bed.gz -O /path/to/repeats.bed.gz

# Rust environment for TRGT
conda create -n rust_trgt
conda activete rust_trgt
mamba install -c conda-forge rust -y
mamba install -c conda-forge clangdev -y
echo 'export LIBCLANG_PATH=$CONDA_PREFIX/lib' >> ~/.bashrc
source ~/.bashrc
```

# References
- Porubsky, D., Dashnow, H., Sasani, T.A. et al. Human de novo mutation rates from a four-generation pedigree reference. Nature 643, 427–436 (2025). https://doi.org/10.1038/s41586-025-08922-2
