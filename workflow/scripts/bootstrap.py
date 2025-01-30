from itertools import combinations
import math

import polars as pl


data = pl.read_csv(snakemake.input[0], separator="\t")

def bootstrap(data):
    for _ in range(snakemake.params.n_bootstraps):
        sample = data.group_by(*snakemake.params.vars).map_groups(lambda group: group.sample(n=group.height, with_replacement=True))
        yield sample


def sample_effect_sizes(sample):
    means = sample.group_by(snakemake.params.vars).mean().with_row_index()
    assert len(means.columns) == 4
    idx = means.get_column("index").to_list()
    comps = list(combinations(idx, 2))

    def comp_effect_size(comp):
        a = means.filter(pl.col("index") == comp[0])
        b = means.filter(pl.col("index") == comp[1])
        a_val = a.select(pl.col("value")).item()
        b_val = b.select(pl.col("value")).item()
        var_cols = pl.col("*").exclude("index", "value")
        return {"group_a": a.select(var_cols).row(0), "group_b": b.select(var_cols).row(0),  "log2_fold_change": math.log2(b_val / a_val)}

    effect_sizes = pl.from_dicts(comp_effect_size(comp) for comp in comps)
    return effect_sizes


# def accumulated_effect_size(sample):
#     sizes = sample_effect_sizes(sample)
#     breakpoint()
#     return sizes.select(pl.col("log2_fold_change").sum()).unique().item()


class ComparisonSummary:
    def __init__(self, comparison):
        (self.group_a, self.group_b), self.comparison = comparison

    def hist(self):
        fc = self.comparison.get_column("log2_fold_change")
        fc.hist(bins=pl.linear_space(start=fc.min(), end=fc.max(), num_samples=10))
        return self.comparison.get_column("log2_fold_change").hist(bin_count=30).with_columns(
            pl.lit(self.group_a).alias("group_a"),
            pl.lit(self.group_b).alias("group_b"),
            ((pl.col("breakpoint") + pl.col("breakpoint").shift()) / 2.0).alias("center"),
        ).filter(
            ~pl.col("center").is_null() & ~pl.col("center").is_infinite()
        )
    
    def confidence_interval(self):
        fc = self.comparison.get_column("log2_fold_change")
        return pl.DataFrame({"lower": fc.quantile(0.05), "upper": fc.quantile(0.95), "group_a": [self.group_a], "group_b": [self.group_b]})

effect_sizes = sample_effect_sizes(data)

bootstraps = pl.concat([sample_effect_sizes(sample) for sample in bootstrap(data)])

bootstrap_summaries = [ComparisonSummary(comparison) for comparison in bootstraps.group_by("group_a", "group_b")]
bootstrap_hists = pl.concat([summary.hist() for summary in bootstrap_summaries])
bootstrap_cis = pl.concat([summary.confidence_interval() for summary in bootstrap_summaries])

bootstrap_hists.write_parquet(snakemake.output.hists)
bootstrap_cis.write_parquet(snakemake.output.cis)