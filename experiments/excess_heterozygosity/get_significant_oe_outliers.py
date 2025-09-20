# Get Bonferroni-significant OE outliers

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import src
from src import _constants as C

FILE_IN = "data/final/gnomad_snrna_variants_hwe_stats.tsv"
FILE_OUT = "data/final/oe_outliers_bonferroni_significant.tsv"


def main():
    return (
        pd.read_csv(FILE_IN, sep="\t")
        .check.nrows(check_name="Input rows")
        .loc[lambda x: x["bfr_sig"]]
        .check.nrows(check_name="Bonferroni significant HWE outliers")
        .check.function(
            lambda x: (x["oe"] > 1).sum(), check_name="Bonferroni significant with OE>1"
        )
        .check.function(
            lambda x: x.loc[x["oe"] > 1],
            check_name="Bonferroni significant variants with O/E>1",
        )
        .assign(
            rank=lambda x: x["oe"].rank(method="first", ascending=False).astype(int)
        )
        .loc[:, ["oe", "rank", "gene_type"]]
        .check.write(FILE_OUT, index=False)
    )


if __name__ == "__main__":
    main()
