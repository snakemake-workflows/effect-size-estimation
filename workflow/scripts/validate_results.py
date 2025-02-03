import polars as pl

data = pl.read_parquet(snakemake.input[0])

expected_log2fc = snakemake.params.log2fc

lower, upper, _, _ = data.row(0)

if expected_log2fc < lower or expected_log2fc > upper:
    raise ValueError("Log2 fold change not within confidence interval")