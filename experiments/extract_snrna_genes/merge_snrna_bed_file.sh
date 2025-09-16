#!/usr/bin/env bash
# Merge overlapping intervals in the snRNA bed file
set -euo pipefail
source src/utils.sh

FILE_IN="data/interim/snrnas.bed"
FILE_OUT="data/interim/snrnas_merged.bed"

< $FILE_IN tee >(log "Input lines: " $(wc -l)) \
| bedtools merge -i - \
| tee >(log "Output lines: " $(wc -l)) \
> $FILE_OUT

sleep 0.1 # Wait for logs to print