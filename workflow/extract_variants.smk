rule all:
    input:
        "data/interim/gnomad_snrna_variants_tidy.tsv",

rule extract_gnomad_variants:
    input:
        "experiments/extract_variants/extract_gnomad_variants_in_snrnas.sh",
        "data/interim/snrnas_merged.bed",
    output:
        "data/interim/gnomad_snrna_variants.vcf.gz",
        "data/interim/gnomad_snrna_variants.vcf.gz.tbi",
        "data/stats/gnomad_snrna_bcftools_stats.txt",
    conda:
        "vep"
    shell:
        "bash experiments/extract_variants/extract_gnomad_variants_in_snrnas.sh"


rule tidy_gnomad_variants:
    input:
        "data/interim/gnomad_snrna_variants.vcf.gz",
        "data/interim/snrnas.bed",
        "experiments/extract_variants/filter_gnomad_variants.sh", 
    output:
        "data/interim/gnomad_snrna_variants_tidy.tsv",
        "data/interim/header_bed_info.txt",
    conda: "vep"
    shell: "bash experiments/extract_variants/filter_gnomad_variants.sh"