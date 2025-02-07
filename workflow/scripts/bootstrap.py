from itertools import combinations
import math

import polars as pl
from statsmodels.stats.nonparametric import rank_compare_2indep

pl.set_random_seed(snakemake.params.seed)

data = pl.read_parquet(snakemake.input.data)

var_values = data.select(snakemake.params.vars).unique(maintain_order=True)

value_col = snakemake.params.value

comparisons = list(
    combinations(
        var_values.select(pl.struct(snakemake.params.vars).alias("case")).get_column(
            "case"
        ),
        2,
    )
)


def bootstrap(data):
    for _ in range(snakemake.params.n_bootstraps):
        sample = data.group_by(*snakemake.params.vars).map_groups(
            lambda group: group.sample(n=group.height, with_replacement=True)
        )
        yield sample


def sample_effect_sizes(sample):
    means = sample.group_by(snakemake.params.vars).mean()

    def comp_effect_size(comp):
        a = means.filter([pl.col(var) == comp[0][var] for var in snakemake.params.vars])
        b = means.filter([pl.col(var) == comp[1][var] for var in snakemake.params.vars])
        a_val = a.select(pl.col(value_col)).item()
        b_val = b.select(pl.col(value_col)).item()
        var_cols = pl.col("*").exclude("index", value_col)
        return {
            "group_a": a.select(var_cols).row(0),
            "group_b": b.select(var_cols).row(0),
            "log2_fold_change": math.log2(b_val / a_val),
        }

    effect_sizes = pl.from_dicts(comp_effect_size(comp) for comp in comparisons)
    return effect_sizes


class ComparisonSummary:
    def __init__(self, comparison):
        (self.group_a, self.group_b), self.comparison = comparison

    def hist(self):
        fc = self.comparison.get_column("log2_fold_change")
        fc_max = fc.max()
        n = 10
        if fc.min() >= fc.max():
            fc_max += 0.1
            n = 2
        fc.hist(bins=pl.linear_space(start=fc.min(), end=fc_max, num_samples=n))
        return (
            self.comparison.get_column("log2_fold_change")
            .hist(bin_count=30)
            .with_columns(
                pl.lit(self.group_a).alias("group_a"),
                pl.lit(self.group_b).alias("group_b"),
                ((pl.col("breakpoint") + pl.col("breakpoint").shift()) / 2.0).alias(
                    "center"
                ),
            )
            .filter(~pl.col("center").is_null() & ~pl.col("center").is_infinite())
        )

    def confidence_interval(self):
        fc = self.comparison.get_column("log2_fold_change")
        return pl.DataFrame(
            {
                "lower": fc.quantile(0.05),
                "upper": fc.quantile(0.95),
                "group_a": [self.group_a],
                "group_b": [self.group_b],
            }
        )

def brunner_munzel_pvalue(comparison):
    group_a = comparison["group_a"]
    group_b = comparison["group_b"]
    x = data.filter(pl.col("case") == group_a).get_column(value_col)
    y = data.filter(pl.col("case") == group_b).get_column(value_col)
    return rank_compare_2indep(x, y, use_t=False).pvalue

effect_sizes = sample_effect_sizes(data)

bootstraps = pl.concat([sample_effect_sizes(sample) for sample in bootstrap(data)])

bootstrap_summaries = [
    ComparisonSummary(comparison)
    for comparison in bootstraps.group_by("group_a", "group_b", maintain_order=True)
]
bootstrap_hists = pl.concat([summary.hist() for summary in bootstrap_summaries])
bootstrap_cis = pl.concat(
    [summary.confidence_interval() for summary in bootstrap_summaries]
)

data = data.with_columns(pl.concat_list(snakemake.params.vars).alias("case"))

bootstrap_cis = bootstrap_cis.with_columns(pl.struct("group_a", "group_b").map_elements(brunner_munzel_pvalue, return_dtype=float).alias("brunner_munzel_pvalue"))

bootstrap_hists.write_parquet(snakemake.output.hists)
bootstrap_cis.write_parquet(snakemake.output.cis)
