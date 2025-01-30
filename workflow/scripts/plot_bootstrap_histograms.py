import polars as pl
import altair as alt

hists = pl.read_parquet(snakemake.input[0]).with_columns(
    [pl.col(group).list.get(i).alias(f"{group}_{varname}") for i, varname in enumerate(snakemake.params.vars) for group in ["group_a", "group_b"]],
).filter(pl.col("count") > 0)

alt.Chart(hists).mark_tick(tooltip=True).encode(
    alt.X("center").axis(title="log2 fold change"),
    alt.Size("count"),
    alt.Row("group_a_diet", type="nominal", title=None),
    alt.Y("group_b_diet", type="nominal", title=None),
    alt.Color("group_a_genotype", title=None),
    alt.Column("group_b_genotype", title=None),
).resolve_scale(x="independent", y="independent").properties(
    width=500,
).save(snakemake.output[0])