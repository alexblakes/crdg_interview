import pandas as pd
import pandas_checks as pdc
import numpy as np
import src

FILE_IN = "data/interim/gnomad_snrna_variants_tidy.tsv"
FILE_OUT = "data/interim/gnomad_snrna_variants_tidy_density.tsv"

# **NOTE** Several genes have the same HGNC symbol, but different ENSG IDs.
#   Grouping should be done by ENSG ID. Specifically U1, U2, U4, U6, U7.


def parse_gnomad_variants(df):
    return (
        df.pipe(
            lambda x: x.assign(
                **pd.DataFrame(
                    x["bed_info"].str.split(",").to_list(),
                    columns=["ensg", "symbol", "length"],
                    index=x.index,
                )
            ).drop(columns=["bed_info"])
        )
        .astype({"length": "int"})
        .assign(
            allele_type=lambda x: np.where(x["allele_type"] == "snv", "SNV", "Indel"),
            gene_type=lambda x: np.where(
                x["symbol"].str.endswith("P"), "Pseudogene", "snRNA"
            ),
        )
        .check.value_counts("allele_type", check_name="Allele type value counts:")
        .check.function(
            lambda x: x.drop_duplicates("ensg")["gene_type"].value_counts(),
            check_name="Gene type value counts:",
        )
        .check.ndups(
            subset=["chrom", "pos", "ref", "alt"], check_name="Duplicate variants"
        )
        .check.nrows(check_name="Number of variants")
        .check.nunique("ensg", check_name="Number of genes")
        .check.assert_greater_than(
            0, subset="ac", pass_message="Pass assert: all AC > 0", verbose=True
        )
        .check.function(
            lambda x: x[["ensg", "symbol"]]
            .drop_duplicates()
            .loc[lambda x: x.duplicated("symbol", keep=False)]
            .sort_values("symbol")
            .value_counts("symbol"),
            check_name="Duplicate gene symbols, different ENSG IDs:",
        )
        .check.function(
            lambda x: x[["chrom", "ensg"]]
            .drop_duplicates()
            .loc[lambda x: x.duplicated("ensg", keep=False)],
            check_name="ENSG IDs on multiple chromosomes:",
        )
    )


def main():
    return (
        pd.read_csv(FILE_IN, sep="\t", na_values=".")
        .pipe(parse_gnomad_variants)
        .check.write(FILE_OUT, index=False)
    )


if __name__ == "__main__":
    main()
