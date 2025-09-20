# Get proportion of singletons for gnomAD SNVs in snRNA genes and pseudogenes

import pandas as pd
from scipy import stats
from statsmodels.stats import proportion

import src
from src import pandas_utils as pdu

FILE_IN = "data/final/gnomad_snrna_variants_hwe_stats_snvs.tsv"
FILE_OUT = "data/stats/snv_proportion_singletons_by_gene_type.tsv"


def get_proportion_ci(row):
    count = int(row["n_singletons"])
    n = int(row["n_variants"])
    ci_lo, ci_hi = stats.binomtest(count, n).proportion_ci(
        confidence_level=0.95, method="exact"
    )
    return ci_lo, ci_hi


def main():
    return (
        pd.read_csv(FILE_IN, sep="\t")
        .groupby("gene_type")
        .agg(n_singletons=("ac", lambda x: (x == 1).sum()), n_variants=("ac", "size"))
        .astype({"n_variants": int, "n_singletons": int})
        .assign(proportion_singletons=lambda x: x["n_singletons"] / x["n_variants"])
        .pipe(
            pdu.assign_with_per_row_fn, get_proportion_ci, new_cols=["ci_lo", "ci_hi"]
        )
        .assign(
            err_lo=lambda x: x["proportion_singletons"] - x["ci_lo"],
            err_hi=lambda x: x["ci_hi"] - x["proportion_singletons"],
        )
        .check.function(lambda x: x, check_name="Proportion of singletons:")
        .check.write(FILE_OUT, index=True)
    )


if __name__ == "__main__":
    main()
