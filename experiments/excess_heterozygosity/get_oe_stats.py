# Get number of observed and expected heterozygous genotypes, plus Fisher's
# exact statistics

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats import contingency_tables

from experiments.variant_density import parse_gnomad_variants
import src
from src import pandas_utils as pdu

FILE_IN = "data/interim/gnomad_snrna_variants_tidy_inbreeding_coeff.tsv"
FILE_OUT = "data/final/gnomad_snrna_variants_hwe_stats.tsv"


def per_row_chi2(row):
    n_genotypes = row["an"] / 2
    obs = row["het_obs"]
    obs_remainder = n_genotypes - obs
    exp = row["het_exp"]
    exp_remainder = n_genotypes - exp

    try:
        odds_ratio, pval = stats.fisher_exact(
            [[obs, obs_remainder], [exp, exp_remainder]]
        )
    except ValueError:
        odds_ratio, pval = np.nan, np.nan

    return odds_ratio, pval


def per_row_fisher_exact(row):
    n_genotypes = row["an"] / 2
    obs = row["het_obs"]
    obs_remainder = n_genotypes - obs
    exp = row["het_exp"]
    exp_remainder = n_genotypes - exp

    table = [[obs, obs_remainder], [exp, exp_remainder]]
    table2x2 = contingency_tables.Table2x2(table)

    try:
        odds_ratio = table2x2.oddsratio
        pval = stats.fisher_exact(table).pvalue
    except ValueError:
        odds_ratio, pval = np.nan, np.nan

    return odds_ratio, pval


def main():
    return (
        pd.read_csv(FILE_IN, sep="\t")
        .check.nrows(check_name="Input variants")
        .loc[lambda x: ~x["chrom"].isin(["chrX", "chrY"])]
        .check.nrows(check_name="Variants after dropping sex chromosomes")
        .check.value_counts("filter", check_name="Variant counts by filter:")
        .pipe(parse_gnomad_variants.parse_gnomad_variants)
        .loc[
            :,
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                "ac",
                "an",
                "af",
                "nhomalt",
                "allele_type",
                "coi",
                "symbol",
                "ensg",
                "gene_type",
            ],
        ]
        .assign(
            oe=lambda x: 1 - x["coi"],
            het_obs=lambda x: x["ac"] - 2 * x["nhomalt"],
            het_exp=lambda x: x["het_obs"] / x["oe"],
        )
        .check.function(
            lambda x: x["het_exp"].dropna().shape[0],
            check_name='Variants with valid "het_exp" values: ',
        )
        .pipe(
            pdu.assign_with_per_row_fn,
            per_row_fisher_exact,
            new_cols=["odds_ratio", "p_val"],
        )
        .dropna(subset="p_val")
        .check.nrows(check_name="Variants with valid chi squared P values:")
        .check.function(
            lambda x: f"Bonferroni significance threshold = {0.05 / len(x):.2e}"
        )
        .sort_values("oe", ascending=False)  # Important for ranking
        .assign(
            bfr_sig=lambda x: np.where(x["p_val"] < 0.05 / len(x), True, False),
            rank_oe=lambda x: x["oe"].rank(ascending=False, method="first").astype(int),
            rank_p=lambda x: x["p_val"]
            .rank(ascending=True, method="first")
            .astype(int),
        )
        .check.function(
            lambda x: x["bfr_sig"].sum(),
            check_name="Number of Bonferroni-significant variants:",
        )
        .check.write(FILE_OUT, index=False)
    )


if __name__ == "__main__":
    main()
