import sys
sys.stderr = open(snakemake.log[0], "w")

import math
import polars as pl
from scipy.stats import ttest_ind


def validate_prediction(test: str, ci_path: str, data_path: str):
    ci = pl.read_parquet(ci_path)

    test_config = snakemake.config["tests"][test]

    expected_log2fc = test_config["log2fc"]

    lower, upper, _, _ = ci.row(0)

    data = pl.read_csv(data_path, separator="\t")

    vars = snakemake.config["datasets"][test]["variables"]

    means = data.group_by(vars).mean()
    assert means.height == 2
    naive_log2fc = math.log2(means.row(1)[-1] / means.row(0)[-1])

    group_data = [group["value"].to_list() for _, group in data.group_by(vars)]

    ttest = ttest_ind(*group_data)

    return pl.DataFrame(
        {
            "test": [test],
            "seed": test_config["seed"],
            "expected_log2fc": [float(expected_log2fc)],
            "naive_log2fc": naive_log2fc,
            "ttest_pval": ttest.pvalue,
            "lower": [float(lower)],
            "upper": [float(upper)],
            "approx_within_ci": [
                math.floor(lower * 10) / 10
                <= expected_log2fc
                <= math.ceil(upper * 10) / 10
            ],
        }
    )


predictions = pl.concat(
    [
        validate_prediction(test, ci_path, data_path)
        for test, ci_path, data_path in zip(
            snakemake.params.tests, snakemake.input.cis, snakemake.input.raw
        )
    ]
)

invalid = predictions.filter(~pl.col("approx_within_ci"))

predictions.write_csv(snakemake.output[0], separator="\t")

if invalid.height > 0:
    raise ValueError(
        f"Some predictions are not within the confidence interval:\n{invalid}"
    )
