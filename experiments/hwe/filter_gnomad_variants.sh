#!/usr/bin/env bash
# Tidy gnomAD variants for HWE / inbreeding coefficient analysis
# Script adapted from experiments/extract_variants/filter_gnomad_variants.sh
set -euo pipefail

source src/utils.sh

FILE_IN="data/interim/gnomad_snrna_variants.vcf.gz"
FILE_OUT="data/interim/gnomad_snrna_variants_tidy_inbreeding_coeff.tsv"
FILE_BED="data/interim/snrnas.bed"
FILE_HEADER="data/interim/header_bed_info.txt" # Created earlier in the pipeline

# Add header line to output file
echo -e "chrom\tpos\tref\talt\tfilter\tac\tan\taf\tnhomalt\tac_nfe\tan_nfe\taf_nfe\tnhomalt_nfe\taf_popmax\tallele_type\tcadd_phred\tcoi\tbed_info" \
> $FILE_OUT

# Filter gnomAD variants, select a subset of INFO fields, and annotate with 
# snRNA data from the BED file
bcftools view -Ou $FILE_IN \
| tee >(log "BCFtools count, input file:\n" "$(bcftools +counts -)") \
| bcftools view -i 'FILTER="PASS" || FILTER="InbreedingCoeff"' -Ou \
| tee >(log "BCFtools count, filter PASS or InbreedingCoeff:\n" "$(bcftools +counts -)") \
| bcftools view -i 'INFO/variant_type="snv"' -Ou \
| tee >(log "BCFtools count, include SNVs only:\n" "$(bcftools +counts -)") \
| bcftools view -i 'INFO/n_alt_alleles=1' -Ou \
| tee >(log "BCFtools count, include n_alt_alleles=1:\n" "$(bcftools +counts -)") \
| bcftools annotate \
    -x ^INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt,INFO/AC_nfe,INFO/AN_nfe,INFO/AF_nfe,INFO/nhomalt_nfe,INFO/AF_popmax,INFO/allele_type,INFO/cadd_phred,INFO/InbreedingCoeff \
    -Ou \
| bcftools annotate \
    --annotations $FILE_BED \
    --columns CHROM,FROM,TO,bed_info \
    --header-lines $FILE_HEADER \
    -Ou \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/nhomalt\t%INFO/AC_nfe\t%INFO/AN_nfe\t%INFO/AF_nfe\t%INFO/nhomalt_nfe\t%INFO/AF_popmax\t%INFO/allele_type\t%INFO/cadd_phred\t%INFO/InbreedingCoeff\t%bed_info\n' \
| tee >(log "Filter value counts:\n" "$(cut -f5 | sort | uniq -c)") \
>> $FILE_OUT

sleep 1