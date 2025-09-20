rule all:
    input:
        "data/final/gnomad_snrna_variants_hwe_stats.tsv",


rule filter_gnomad_variants:
    input:
        "data/interim/gnomad_snrna_variants.vcf.gz",
        "experiments/excess_heterozygosity/filter_gnomad_variants.sh",
        "data/interim/snrnas.bed",
        "data/interim/header_bed_info.txt",
    output:
        "data/interim/gnomad_snrna_variants_tidy_inbreeding_coeff.tsv",
    conda:
        "vep"
    shell:
        "bash experiments/excess_heterozygosity/filter_gnomad_variants.sh"


rule get_excess_het_stats:
    input:
        "data/interim/gnomad_snrna_variants_tidy_inbreeding_coeff.tsv",
        "experiments/excess_heterozygosity/get_oe_stats.py",
    output:
        "data/final/gnomad_snrna_variants_hwe_stats.tsv",
    log:
        "data/logs/experiments.excess_heterozygosity.get_oe_stats.log",
    conda:
        "crdg_test"
    shell:
        "python3 -m experiments.excess_heterozygosity.get_oe_stats 2&> {log}"
