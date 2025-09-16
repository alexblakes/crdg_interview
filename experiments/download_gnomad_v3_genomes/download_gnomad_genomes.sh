#!/bin/bash
# Script to download gnomAD VCF and TBI for one chromosome.
set -euo pipefail

url=$1
dir_download=$2

file_name=$(basename $url)
file_path="${dir_download}/${file_name}"
file_md5="${file_path}.md5"

wget -O - $url \
| tee $file_path \
| md5sum \
| awk -v fname=$file_name '{print $1"\t"fname}' \
> $file_md5

sleep 1