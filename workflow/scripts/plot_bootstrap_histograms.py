import math
import polars as pl
import altair as alt

hists = (
    pl.read_parquet(snakemake.input[0]).with_row_index()
    .with_columns(pl.col("group_a").list.join(":"), pl.col("group_b").list.join(":"))
    .with_columns(
        pl.concat_list("group_a", "group_b").list.join(" vs ").alias("comparison")
    )
)

log2fcs = hists.get_column("center")

def fold_change_ticks():
    fc_lower = 2 ** log2fcs.min()
    fc_upper = 2 ** log2fcs.max()


    if fc_lower < 1.0:
        yield from map(math.log2, reversed(range(1, math.ceil(fc_lower) + 1)))
    else:

    if log2fc == 0:
        return 0
    if lower:
        if log2fc < 0:
            return math.log2fc(1 / math.ceil(2 ** -log2fc))
        else:
            return math.log2fc(math.floor(2 ** log2fc))
    else:
        if log2fc < 0:
            return math.log2fc(1 / math.floor(2 ** -log2fc))
        else:
            return math.log2fc(math.ceil(2 ** log2fc))



ticks = list(range(math.log2(math.floor(2 ** log2fc.min())), math.log2(math.ceil(2 ** log2fc.max()))))

label_charts = [
    alt.Chart(hists.select("comparison", group).unique(maintain_order=True)).mark_text(align="right").encode(
        alt.Y("comparison", type="nominal", sort=None).axis(None).title("comparison"),
        alt.Text(group),
        alt.Color(group, sort=None).legend(title=None),
    ).properties(
        view=alt.ViewConfig(strokeWidth=0),
    )
    for group in ["group_a", "group_b"]
]

hist_chart = alt.Chart(hists).mark_tick(tooltip=True).encode(
    alt.X("center").axis(title="log2 fold change", ticks=ticks),
    alt.Size("count").legend(None),
    alt.Y("comparison", type="nominal", sort=None).axis(labels=False, title=None),
)

alt.hconcat(*(label_charts + [hist_chart]), spacing=0).save(
    snakemake.output[0]
)

# hists = (
#     pl.read_parquet(snakemake.input[0])
#     .with_columns(
#         [
#             pl.col(group).list.get(i).alias(f"{group}_{varname}")
#             for i, varname in enumerate(snakemake.params.vars)
#             for group in ["group_a", "group_b"]
#         ],
#     )
#     .filter(pl.col("count") > 0)
#     .with_columns(
#         pl.concat_list("group_a_genotype", "group_b_genotype")
#         .list.join("/")
#         .alias("genotype")
#     )
# )
# breakpoint()

# alt.Chart(hists).mark_tick(tooltip=True).encode(
#     alt.X("center").axis(title="log2 fold change"),
#     alt.Size("count"),
#     alt.Row("group_a_diet", type="nominal", title=None),
#     alt.Y("group_b_diet", type="nominal", title=None),
#     alt.Color("genotype", title=None),
#     alt.YOffset("genotype"),
# ).resolve_scale(x="shared", y="independent").properties(
#     width=500,
# ).save(
#     snakemake.output[0]
# )
