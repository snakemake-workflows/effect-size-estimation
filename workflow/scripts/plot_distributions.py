import sys
sys.stderr = open(snakemake.log[0], "w")

from collections import defaultdict
import math
import polars as pl
import altair as alt

pl.Config.set_tbl_rows(100)

EPSILON = 0.1
VAR_SEP = ", "
FIRST_VAR_SEP = ": "

vars = snakemake.params.vars
mode = snakemake.wildcards.mode

assert (
    snakemake.params.min_fold_change >= 1.0
), "min_fold_change must be greater than 1.0"
min_conservative_log2_fold_change = math.log2(snakemake.params.min_fold_change)

# Stably sort the data such that the first variable is grouped together
# This is necessary because we draw a rule to label x-axis with first variable
# values.
data = pl.read_parquet(snakemake.input.data).with_row_count("idx").with_columns(
    pl.first("idx").over(vars[0]).alias("group_first_idx")
).sort(["group_first_idx", "idx"]).drop(["idx", "group_first_idx"])

if len(vars) > 2:
    combined_var = VAR_SEP.join(vars[1:])
    # combine vars[1:] into a single variable
    data = data.with_columns(
        pl.concat_list(vars[1:]).list.join(VAR_SEP).alias(combined_var),
    )
    vars = [vars[0], combined_var]

color_col = "case" if mode == "all" else vars[1]

var_values = (
    data.get_column(vars[1]).unique(maintain_order=True).to_list()
)
var_indexes = {value: i for i, value in enumerate(var_values)}
data = data.with_columns(
    pl.col(vars[1]).replace_strict(var_indexes).alias("index"),
    pl.concat_list(vars).list.join(FIRST_VAR_SEP).alias("case"),
)

color_order = data.get_column(color_col).unique(maintain_order=True).to_list()

def get_conservative_log2_fold_change(ci) -> float:
    lower, upper = ci["lower"], ci["upper"]
    if lower <= 0.0 <= upper:
        return 0.0
    else:
        return lower if abs(lower) < abs(upper) else upper


cis = (
    pl.read_parquet(snakemake.input.cis)
    .with_columns(
        pl.struct(["lower", "upper"])
        .map_elements(get_conservative_log2_fold_change, return_dtype=float)
        .alias("conservative_log2_fold_change")
    )
    .with_columns(
        [
            pl.col(f"group_{group}").list.slice(
                # the first element can be taken as is, the second and
                # potential rest is joined with ','
                i, i + 1 if i == 0 else None
            ).list.join(VAR_SEP).alias(f"{varname}_{group}")
            for group in ["a", "b"]
            for i, varname in enumerate(vars)
        ],
    )
)

cis = cis.with_columns(
    *[
        pl.Series(
            [
                var_indexes[value]
                for value in cis.get_column(f"{vars[1]}_{group}")
            ]
        ).alias(f"index_{group}")
        for group in ["a", "b"]
    ]
)


def fmt_fold_change(effect) -> str:
    log2_fold_change = effect["conservative_log2_fold_change"]
    pvalue = effect["brunner_munzel_adjusted_pvalue"]
    if log2_fold_change > 0:
        direction = "⏶"
        fold_change = 2**log2_fold_change
    else:
        direction = "⏷"
        fold_change = 1 / (2**log2_fold_change)
    if log2_fold_change == 0:
        return "="
    else:
        if fold_change < 1.05:
            return "≈"
        else:
            if mode == "all":
                pvalmark = ""
            else:
                pvalmark = "*"
                if pvalue <= 0.0001:
                    pvalmark *= 4
                elif pvalue <= 0.001:
                    pvalmark *= 3
                elif pvalue <= 0.01:
                    pvalmark *= 2
                elif pvalue > 0.05:
                    pvalmark = ""
            return f"{direction}{fold_change:.1f} {pvalmark}"


cis = cis.with_columns(
    pl.struct("conservative_log2_fold_change", "brunner_munzel_adjusted_pvalue")
    .map_elements(fmt_fold_change, return_dtype=str)
    .alias("fold change"),
    pl.concat_list([f"{var}_a" for var in vars])
    .list.join(FIRST_VAR_SEP)
    .alias("case_a"),
    pl.concat_list([f"{var}_b" for var in vars])
    .list.join(FIRST_VAR_SEP)
    .alias("case_b"),
)


# generate data frame with two rows for each group_a, group_b pair, one with
# the group_a values and the corresponding vars[0] value and the corresponding index,
# and one with group_b values and the corresponding vars[0] value and the corresponding index
color_spec = alt.Color(color_col, type="nominal").scale(
    domain=color_order, range=snakemake.params.color_scheme
)
if snakemake.wildcards.legend == "yes":
    color_spec = color_spec.legend(title=None)
else:
    color_spec = color_spec.legend(None)

dist_chart = (
    alt.Chart(data)
    .mark_circle(tooltip=True)
    .encode(
        alt.X("case", type="nominal", sort=None).axis(None),
        alt.Y(snakemake.params.value, type="quantitative")
        .scale(zero=False)
        .axis(grid=False, title=None),
        color_spec,
    )
    .properties(
        width=230,
        height=120,
        view=alt.ViewConfig(stroke=None),
        title=snakemake.params.value,
    )
)

if mode == "selected":
    # add an underline for each variable in vars[0]
    # alt.X should be the value of the first value of case of the respective group in the data frame
    # alt.X2 should be the value of the last value of case of the respective group in the data frame
    underline_data = data.group_by(vars[0], maintain_order=True).agg(
        [
            pl.col(vars[0]).first().alias("label"),
            pl.col("case").first().alias("x"),
            pl.col("case").last().alias("x2"),
        ]
    )

    dist_chart += alt.Chart(underline_data).mark_rule(strokeWidth=0.5).encode(
        alt.X("x", type="nominal", sort=None).axis(None),
        alt.X2("x2"),
        alt.Y(value=-2),
    ) + alt.Chart(underline_data).mark_text(dy=-4, align="left", fontSize=8).encode(
        alt.X("x", type="nominal", sort=None).axis(None),
        alt.Text("label"),
        alt.Y(value=-2),
    )


def get_all_effect_chart():
    return (
        alt.Chart(cis)
        .mark_text(fontSize=10, align="center")
        .encode(
            alt.X("case_a", type="nominal", sort=None).axis(None),
            alt.Y("case_b", type="nominal", sort=None).axis(None).scale(reverse=True),
            alt.Text("fold change"),
            alt.Color("case_b", type="nominal").scale(domain=color_order),
        )
        .properties(
            width=230,
            height=120,
            view=alt.ViewConfig(stroke=None),
        )
    )


def get_selected_effect_chart():
    comparisons = pl.read_csv(snakemake.input.comparisons, separator="\t")

    if len(snakemake.params.vars) > 2:
        # combine vars[1:] into a single variable
        comparisons = comparisons.with_columns(
            [
                pl.concat_list([f"{var}_{group}" for var in snakemake.params.vars[1:]]).list.join(VAR_SEP).alias(
                    f"{vars[1]}_{group}"
                )
                for group in ["a", "b"]
            ]
        )
        print(comparisons)

    def swap_colname(col):
        return col[:-1] + ("b" if col.endswith("_a") else "a")

    # using how="diagonal" to concat the dataframes based on column names,
    # not their order, see https://github.com/pola-rs/polars/issues/5789
    comparisons = pl.concat(
        [
            comparisons,
            comparisons.rename({col: swap_colname(col) for col in comparisons.columns}),
        ],
        how="diagonal",
    )
    joincols = [
        f"{var}_{group}"
        for var in vars
        for group in ["a", "b"]
    ]
    print(comparisons.select(joincols))
    print(cis.select(joincols))

    selected_cis = (
        cis.join(
            comparisons,
            how="semi",
            on=joincols,
        )
        .filter(
            pl.col("fold change") != "=",
            pl.col("fold change") != "≈",
            pl.col("conservative_log2_fold_change").abs()
            >= min_conservative_log2_fold_change,
        )
        .with_row_index()
    )
    print(selected_cis)

    case_idx = {
        case: i
        for i, case in enumerate(data.get_column("case").unique(maintain_order=True))
    }

    placements = defaultdict(list)
    placements[0].append(0)
    for i in range(1, selected_cis.height):
        case_a = selected_cis.item(i, "case_a")
        placed = False
        for placement, rows in placements.items():
            max_case = max(case_idx[selected_cis.item(j, "case_b")] for j in rows)
            if case_idx[case_a] > max_case:
                placements[placement].append(i)
                placed = True
                break
        if not placed:
            placements[len(placements)].append(i)

    row_placements = {
        row: placement for placement, rows in placements.items() for row in rows
    }

    selected_cis = selected_cis.with_columns(
        pl.Series([row_placements[row] for row in range(selected_cis.height)]).alias(
            "placement"
        )
    )
    if selected_cis.is_empty():
        return None

    base = (
        alt.Chart(selected_cis)
        .encode(
            alt.X("case_a", type="nominal", sort=None).axis(None),
            alt.Y("placement", type="ordinal").axis(None),
        )
        .properties(view=alt.ViewConfig(stroke=None))
    )

    return base.mark_rule().encode(alt.X2("case_b")) + base.mark_text(
        dy=-6, align="left", fontSize=8
    ).encode(alt.Text("fold change"))


if mode == "all":
    chart = dist_chart & get_all_effect_chart()
else:
    effects = get_selected_effect_chart()
    if effects is not None:
        chart = dist_chart & effects
    else:
        chart = dist_chart

chart.configure_concat(spacing=0).resolve_scale(
    y="independent", x="shared", color="shared"
).save(snakemake.output[0])
# comp_chart.save(snakemake.output[0])
