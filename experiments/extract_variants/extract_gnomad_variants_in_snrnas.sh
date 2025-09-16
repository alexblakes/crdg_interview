#!/usr/bin/env bash
# Extract variants in gnomAD v3 genomes which overlap snRNA genes
set -euo pipefail

DIR_GNOMAD="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/gnomad/v3.1.2/genomes"
export DIR_TMP="/scratch/m40482ab/gnomad_snrna_variants"
export FILE_BED="data/interim/snrnas_merged.bed"
FILE_OUT="data/interim/gnomad_snrna_variants.vcf.gz"
export suffix=".snrna_variants.vcf.gz"
FILE_STATS="data/stats/gnomad_snrna_bcftools_stats.txt"

function extract_snrna_variants() {
    local vcf=$1
    local base=$(basename -s .vcf.bgz $vcf)
    local out="${DIR_TMP}/${base}${suffix}"
    
    bcftools view $vcf -R $FILE_BED -Oz -o $out -W=tbi
}

export -f extract_snrna_variants


# Empty the tmp directory
[[ -d $DIR_TMP ]] && rm -rf $DIR_TMP
mkdir -p $DIR_TMP

# Extract variants in parallel
find $DIR_GNOMAD -type f -name "*.vcf.bgz" \
| sort -V \
| parallel --jobs 80% extract_snrna_variants {}

# Combine variants, get summary statistics
interim_vcf_paths=$(find $DIR_TMP -type f -name "*${suffix}" | sort -V)
bcftools concat -n -Oz -o $FILE_OUT $interim_vcf_paths && tabix $FILE_OUT
bcftools stats $FILE_OUT > $FILE_STATS