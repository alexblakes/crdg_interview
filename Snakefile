module extract_snrna_genes:
    snakefile:
        "workflow/extract_snrna_genes.smk"
    config:
        {}


use rule * from extract_snrna_genes as extract_snrna_genes_*


module extract_variants:
    snakefile:
        "workflow/extract_variants.smk"
    config:
        {}


use rule * from extract_variants as extract_variants_*


module variant_density:
    snakefile:
        "workflow/variant_density.smk"
    config:
        {}


use rule * from variant_density as tidy_data_for_plots_*


rule all:
    input:
        rules.extract_snrna_genes_all.input,
        rules.extract_variants_all.input,
        rules.tidy_data_for_plots_all.input,
    default_target: True
