import sys
sys.stderr = open(snakemake.log[0], "w")

import polars as pl

# Read the parquet file
df = pl.read_parquet(snakemake.input[0])

# Convert list columns to strings
for col in df.columns:
    if df[col].dtype == pl.List:
        df = df.with_columns(pl.col(col).list.join(",").alias(col))

# Write to TSV
df.write_csv(snakemake.output[0], separator="\t")