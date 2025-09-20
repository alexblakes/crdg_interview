# Filter for gnomAD SNVs

import pandas as pd
import src

FILE_IN = "data/final/gnomad_snrna_variants_hwe_stats.tsv"
FILE_OUT = "data/final/gnomad_snrna_variants_hwe_stats_snvs.tsv"


def main():
    return (
        pd.read_csv(FILE_IN, sep="\t", usecols=["af", "ac", "allele_type", "gene_type"])
        .check.nrows(check_name="Input variant count")
        .loc[lambda x: x["allele_type"] == "SNV"]
        .drop(columns=["allele_type"])
        .check.nrows(check_name="SNV count")
        .check.value_counts("gene_type", check_name="SNVs per gene type:")
        .check.write(FILE_OUT, index=False)
    )


if __name__ == "__main__":
    main()
