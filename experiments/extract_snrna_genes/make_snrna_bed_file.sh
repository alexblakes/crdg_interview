#!/usr/bin/env bash
# Extract snRNA genes and snRNA pseudogenes from GENCODE.
# Notes:
#   "snRNA_pseudogene" is a valid biotype but not used in GENCODE v48.

set -euo pipefail

source src/utils.sh # Logging

FILE_GENCODE="data/raw/gencode.v48.annotation.gtf.gz"
FILE_OUT="data/interim/snrnas.bed"
script="${0##*/}"

echo "Running ${script}"
echo "Writing to $FILE_OUT"

# grep for snRNA exons
# log whether each snRNA has only one exon
# awk
#   Output bed file
#   Comma-separated name column of Ensembl gene ID, HGNC symbol, length
< $FILE_GENCODE zcat \
| grep -v "^#" \
| grep 'gene_type "snRNA"' \
| tee >(log "Unique features: " $(cut -f 3 | sort -u)) \
| awk -v FS="\t" '$3=="exon"' \
| tr -s ";" " " \
| tr " " "\t" \
| tee >(log "Number of snRNA exons: " $(wc -l)) \
| tee >(log "Unique snRNA gene IDs: " $(cut -f 10 | sort -u | wc -l)) \
| cut -f 1,4,5,10,16 \
| sed 's/\"//g' \
| awk -v OFS="\t" '$2-=1 {print $1, $2, $3, $4"\,"$5"\,"$3-$2}' \
| tee >(log 'Number of pseudogenes (HGNC Symbol ends with "P"): ' "$(grep P$ | wc -l)") \
| tee >(log 'Number of "expressed" snRNA genes (HGNC Symbol does not end with "P"): ' "$(grep -v P$ | wc -l)") \
> $FILE_OUT

sleep 0.1 # Print log messages before exit

