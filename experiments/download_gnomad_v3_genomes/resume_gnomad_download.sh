#!/bin/bash
# Resume download of a gnomAD VCF after a timeout.
set -euo pipefail

url=$1
dir_download=$2

file_name=$(basename $url)
file_path="${dir_download}/${file_name}"
file_md5="${file_path}.md5"

wget -c -O $file_path $url

md5sum $file_path \
| awk -v fname=$file_name '{print $1"\t"fname}' \
> $file_md5

sleep 1