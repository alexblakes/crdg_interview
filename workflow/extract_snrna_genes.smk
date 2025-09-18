rule all:
    input:
        "data/interim/snrnas_merged.bed",


rule make_snrna_bed_file:
    input:
        "data/raw/gencode.v48.annotation.gtf.gz",
        "experiments/extract_snrna_genes/make_snrna_bed_file.sh",
    output:
        "data/interim/snrnas.bed",
    conda:
        "vep"
    shell:
        "bash experiments/extract_snrna_genes/make_snrna_bed_file.sh"


rule merge_bed_intervals:
    input:
        "data/interim/snrnas.bed",
        "experiments/extract_snrna_genes/merge_snrna_bed_file.sh",
    output:
        "data/interim/snrnas_merged.bed",
    conda:
        "vep"
    shell:
        "bash experiments/extract_snrna_genes/merge_snrna_bed_file.sh"
