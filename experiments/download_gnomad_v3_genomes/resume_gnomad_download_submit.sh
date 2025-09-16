#!/bin/bash --login
# Submission script to download gnomAD v3 genomes
# Run on the login node
# Run from crdg_interview directory

set -euo pipefail

SCRIPT=$0
JOB_SCRIPT="./experiments/download_gnomad_v3_genomes/resume_gnomad_download.sh"
DIR_GNOMAD_V3="/mnt/bmh01-rds/Ellingford_gene/public_data_resources/gnomad/v3.1.2/genomes"
now=$(date +'%Y-%m-%dT%H:%M:%S:%3N')
DIR_LOGS="data/logs/cluster/download_gnomad_v3_${now}"

mkdir -p $DIR_LOGS

# Submit a download job for each contig
replace_string="<replace_me>"
template_url_vcf="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.${replace_string}.vcf.bgz"
template_url_tbi="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.${replace_string}.vcf.bgz.tbi"

chroms=(chr{1..22} chr{X,Y})

for chrom in ${chroms[@]}; do
    url_vcf=${template_url_vcf/${replace_string}/${chrom}}
    url_tbi=${template_url_tbi/${replace_string}/${chrom}}

    for url in $url_tbi $url_vcf; do
        sbatch \
            -o "${DIR_LOGS}/${chrom}.%jobid.o" \
            -e "${DIR_LOGS}/${chrom}.%jobid.e" \
            --partition="serial" \
            --time="24:00:00" \
            $JOB_SCRIPT $url $DIR_GNOMAD_V3
    done
done
