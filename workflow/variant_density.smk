rule all:
    input:
        "data/final/snrna_variant_counts.tsv",


rule parse_gnomad_variants:
    input:
        "data/interim/gnomad_snrna_variants_tidy.tsv",
        "experiments/variant_density/parse_gnomad_variants.py",
    output:
        "data/interim/gnomad_snrna_variants_tidy_density.tsv",
    conda:
        "crdg_test"
    log:
        "data/logs/experiments.variant_density.parse_gnomad_variants.log",
    shell:
        "python3 -m experiments.variant_density.parse_gnomad_variants 2&> {log}"


rule get_variant_density:
    input:
        "data/interim/gnomad_snrna_variants_tidy_density.tsv",
    output:
        "data/final/snrna_variant_counts.tsv",
    conda:
        "crdg_test"
    log:
        "data/logs/experiments.variant_density.get_variant_density.log",
    shell:
        "python3 -m experiments.variant_density.get_variant_density 2&> {log}"
