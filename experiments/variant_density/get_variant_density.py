# Get variant density per snRNA gene (using ENSG IDs)

import pandas as pd

import src

FILE_IN = "data/interim/gnomad_snrna_variants_tidier.tsv"
FILE_OUT = "data/final/snrna_variant_counts.tsv"


def main():
    return (
        pd.read_csv(FILE_IN, sep="\t")
        .groupby(["ensg", "allele_type"])
        .agg(
            symbol=("symbol", "first"),
            gene_type=("gene_type", "first"),
            length=("length", "first"),
            n_variants=("chrom", "count"),
        )
        .assign(variants_per_nt=lambda x: x["n_variants"] / x["length"])
        .check.info()
        .check.head()
        .check.write(FILE_OUT, index=True)
    )

if __name__ == "__main__":
    main()