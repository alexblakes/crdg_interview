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


use rule * from variant_density as variant_density_*


module excess_heterozygosity:
    snakefile:
        "workflow/excess_heterozygosity.smk"
    config:
        {}


use rule * from excess_heterozygosity as excess_heterozygosity_*


module allele_frequencies:
    snakefile:
        "workflow/allele_frequencies.smk"
    config:
        {}


use rule * from allele_frequencies as allele_frequencies_*


rule all:
    input:
        rules.extract_snrna_genes_all.input,
        rules.extract_variants_all.input,
        rules.variant_density_all.input,
        rules.excess_heterozygosity_all.input,
        rules.allele_frequencies_all.input,
    default_target: True
