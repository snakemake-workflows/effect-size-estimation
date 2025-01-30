import math
import polars as pl
import altair as alt

EPSILON = 0.1

assert len(snakemake.params.vars) == 2, "only two variables are supported for now"

data = pl.read_csv(snakemake.input.data, separator="\t")
var_values = set(data.get_column(snakemake.params.vars[1]))
var_indexes = {value: i for i, value in enumerate(sorted(var_values))}
data = data.with_columns(pl.col(snakemake.params.vars[1]).replace_strict(
    var_indexes
).map_elements(lambda val: val).alias("index"))


def get_conservative_log2_fold_change(ci):
    lower, upper = ci["lower"], ci["upper"]
    if lower <= 0.0 <= upper:
        return 0.0
    else:
        return lower if abs(lower) < abs(upper) else upper


cis = pl.read_parquet(snakemake.input.cis).with_columns(
    pl.struct(["lower", "upper"]).map_elements(get_conservative_log2_fold_change).alias("conservative_log2_fold_change")
).with_columns(
    [
        pl.col(group).list.get(i).alias(f"{group}_{varname}")
        for group in ["group_a", "group_b"]
        for i, varname in enumerate(snakemake.params.vars)
    ],
)

cis = cis.with_columns(
    *[
        pl.Series([
            var_indexes[value]
            for value in cis.get_column(f"{group}_{snakemake.params.vars[1]}")
        ]).alias(f"{group}_index")
        for group in ["group_a", "group_b"]
    ]
).with_row_index()

# generate data frame with two rows for eahc group_a, group_b pair, one with
# the group_a values and the corresponding snakemake.params.vars[0] value and the corresponding index,
# and one with group_b values and the corresponding snakemake.params.vars[0] value and the corresponding index

dist_chart = alt.Chart(data.with_columns()).mark_circle().encode(
    alt.X(snakemake.params.vars[0]),
    alt.Y("value").axis(grid=False),
    alt.XOffset("index"),
    alt.Color(snakemake.params.vars[1]),
).properties(width=200, height=200, view=alt.ViewConfig(stroke=None),)

fc_chart = alt.layer(*([
    alt.Chart(
        pl.concat([
            comparison.select(
                [pl.col(f"{group}_{var}").alias(var) for var in snakemake.params.vars] +
                [
                    pl.col(f"{group}_index").alias("index"),
                    pl.col("conservative_log2_fold_change"),
                    pl.concat_list("group_a", "group_b").list.join(",").alias("comparison")
                ],
            )
            for group in ["group_a", "group_b"]
        ])
    ).mark_line(tooltip=True).encode(
        alt.X(snakemake.params.vars[0]).axis(None),
        alt.XOffset("index"),
        alt.Y("comparison").axis(None),
        alt.Color("conservative_log2_fold_change").scale(scheme="viridis"),
    ).properties(
        width=150,
        height=100,
        view=alt.ViewConfig(stroke=None),
    )
    for _, comparison in cis.group_by("group_a", "group_b")
]))

def fmt_fold_change(log2_fold_change: float) -> str:
    if log2_fold_change > 0:
        direction = "⏶"
        fold_change = 2 ** log2_fold_change
    else:
        direction = "⏷"
        fold_change = 1 / (2 ** log2_fold_change)
    fold_change = round(fold_change)
    if fold_change == 1:
        return "="
    else:
        return f"{direction}{fold_change:.0f}"


comp_chart = alt.Chart(cis.with_columns(pl.col("conservative_log2_fold_change").map_elements(fmt_fold_change).alias("fold change"))).mark_text().encode(
    alt.X(f"group_a_{snakemake.params.vars[0]}").axis(None),
    alt.XOffset("group_a_index"),
    alt.Y(f"group_b_{snakemake.params.vars[0]}").axis(title=None),
    alt.YOffset("group_b_index"),
    alt.Text("fold change"),
    alt.Color(f"group_b_{snakemake.params.vars[1]}").legend(None),
).properties(width=200, height=160, view=alt.ViewConfig(stroke=None),)

(comp_chart & dist_chart).resolve_scale(y="independent", x="shared", color="independent").save(snakemake.output[0])
# comp_chart.save(snakemake.output[0])