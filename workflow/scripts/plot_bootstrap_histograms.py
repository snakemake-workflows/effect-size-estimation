import sys
sys.stderr = open(snakemake.log[0], "w")

import math
import polars as pl
import altair as alt

hists = (
    pl.read_parquet(snakemake.input[0])
    .with_row_index()
    .with_columns(pl.col("group_a").list.join(":"), pl.col("group_b").list.join(":"))
    .with_columns(
        pl.concat_list("group_a", "group_b").list.join(" vs ").alias("comparison")
    )
)

log2fcs = hists.get_column("center")


# we calculate the ticks explicitly such that they have integer values when translating
# them to normal fold changes in the label expression below.
def fold_change_ticks():
    fc_lower = 2 ** log2fcs.min()
    fc_upper = 2 ** log2fcs.max()

    if fc_lower < 1.0:
        if fc_upper < 1.0:
            yield from map(
                lambda v: math.log2(1 / v),
                reversed(range(math.ceil(1 / fc_upper), math.ceil(1 / fc_lower) + 1)),
            )
        else:
            yield from map(
                lambda v: math.log2(1 / v),
                reversed(range(1, math.ceil(1 / fc_lower) + 1)),
            )
            yield from map(math.log2, range(2, math.ceil(fc_upper) + 1))
    else:
        yield from map(math.log2, range(math.floor(fc_lower), math.ceil(fc_upper) + 1))


ticks = list(fold_change_ticks())

label_charts = [
    alt.Chart(hists.select("comparison", group).unique(maintain_order=True))
    .mark_text(align="right")
    .encode(
        alt.Y("comparison", type="nominal", sort=None).axis(None).title("comparison"),
        alt.Text(group),
        alt.Color(group, sort=None).legend(title=None),
    )
    .properties(
        view=alt.ViewConfig(strokeWidth=0),
    )
    for group in ["group_a", "group_b"]
]

hist_chart = (
    alt.Chart(hists)
    .mark_tick(tooltip=True)
    .encode(
        alt.X("center").axis(
            title="fold change",
            values=ticks,
            labelOverlap="greedy",
            labelExpr="datum.value < 0 ? '⏷' + format(pow(2, -datum.value), ',.0f')"
            " : (datum.value == 0.0 ? '=' : "
            "('⏶' + format(pow(2, datum.value), ',.0f')))",
        ).axis(grid=True),
        alt.Size("count").legend(None),
        alt.Y("comparison", type="nominal", sort=None).axis(labels=False, title=None),
    )
)

alt.hconcat(*(label_charts + [hist_chart]), spacing=0).interactive().save(snakemake.output[0])
