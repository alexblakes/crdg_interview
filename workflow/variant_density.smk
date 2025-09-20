rule all:
    input:
        "data/plots/snrna_variant_density.png",


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
        "experiments/variant_density/get_variant_density.py",
    output:
        "data/final/snrna_variant_counts.tsv",
    conda:
        "crdg_test"
    log:
        "data/logs/experiments.variant_density.get_variant_density.log",
    shell:
        "python3 -m experiments.variant_density.get_variant_density 2&> {log}"


rule plot_variant_density:
    input:
        "data/final/snrna_variant_counts.tsv",
        "experiments/variant_density/plot_variant_density.ipynb",
    output:
        "data/plots/snrna_variant_density.png",
        "data/plots/snrna_variant_density.svg",
    conda:
        "crdg_test"
    log:
        "data/logs/experiments.variant_density.plot_variant_density.log",
    notebook:
        "../experiments/variant_density/plot_variant_density.ipynb"
