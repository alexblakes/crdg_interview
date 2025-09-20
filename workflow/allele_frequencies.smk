rule all:
    input:
        "data/plots/snrna_allele_frequency.png",
        "data/plots/snrna_allele_frequency.svg",


rule filter_snvs:
    input:
        "data/final/gnomad_snrna_variants_hwe_stats.tsv",
        "experiments/allele_frequencies/filter_snvs.py",
    output:
        "data/final/gnomad_snrna_variants_hwe_stats_snvs.tsv",
    log:
        "data/logs/experiments.allele_frequencies.filter_snvs.log",
    conda:
        "crdg_test"
    shell:
        "python3 -m experiments.allele_frequencies.filter_snvs 2&>{log}"


rule get_proportion_singletons:
    input:
        "data/final/gnomad_snrna_variants_hwe_stats_snvs.tsv",
        "experiments/allele_frequencies/get_proportion_singletons.py",
    output:
        "data/stats/snv_proportion_singletons_by_gene_type.tsv",
    log:
        "data/logs/experiments.allele_frequencies.get_proportion_singletons.log",
    conda:
        "crdg_test"
    shell:
        "python3 -m experiments.allele_frequencies.get_proportion_singletons 2&>{log}"


rule plot_allele_frequencies:
    input:
        "data/stats/snv_proportion_singletons_by_gene_type.tsv",
        "data/final/gnomad_snrna_variants_hwe_stats_snvs.tsv",
        "experiments/allele_frequencies/plot_allele_frequencies.ipynb",
    output:
        "data/plots/snrna_allele_frequency.png",
        "data/plots/snrna_allele_frequency.svg",
    conda:
        "crdg_test"
    notebook:
        "../experiments/allele_frequencies/plot_allele_frequencies.ipynb"
