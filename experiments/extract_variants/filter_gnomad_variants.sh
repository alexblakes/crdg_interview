#!/usr/bin/env bash
# Tidy gnomAD variants for downstream analysis
# Notes:
#   The VEP annotations in gnomAD v3 are outdated. Annotate with GENCODE v48 
#   instead.
set -euo pipefail

FILE_IN="data/interim/gnomad_snrna_variants.vcf.gz"
FILE_OUT="data/interim/gnomad_snrna_variants_tidy.tsv"
FILE_BED="data/interim/snrnas.bed"
FILE_HEADER="data/interim/header_bed_info.txt"

# Add header line for new INFO field
echo '##INFO=<ID=bed_info,Number=1,Type=String,Description="Ensembl gene ID, HGNC symbol, length.">' \
> $FILE_HEADER

# Add header line to output file
echo -e "chrom\tpos\tref\talt\tac\tan\taf\tnhomalt\tac_nfe\tan_nfe\taf_nfe\tnhomalt_nfe\taf_popmax\tallele_type\tcadd_phred\tbed_info" \
> $FILE_OUT

# Filter gnomAD variants, select a subset of INFO fields, and annotate with 
# snRNA data from the BED file
bcftools view -i 'FILTER="PASS"' -Ou $FILE_IN \
| bcftools annotate \
    -x ^INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt,INFO/AC_nfe,INFO/AN_nfe,INFO/AF_nfe,INFO/nhomalt_nfe,INFO/AF_popmax,INFO/allele_type,INFO/cadd_phred \
    -Ou \
| bcftools annotate \
    --annotations $FILE_BED \
    --columns CHROM,FROM,TO,bed_info \
    --header-lines $FILE_HEADER \
    -Ou \
| bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/nhomalt\t%INFO/AC_nfe\t%INFO/AN_nfe\t%INFO/AF_nfe\t%INFO/nhomalt_nfe\t%INFO/AF_popmax\t%INFO/allele_type\t%INFO/cadd_phred\t%bed_info\n' \
>> $FILE_OUT