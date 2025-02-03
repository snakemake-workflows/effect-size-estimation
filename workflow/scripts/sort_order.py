import polars as pl

data = pl.read_csv(snakemake.input[0], separator="\t")

for var, values in reversed(snakemake.params.order.items()):
    cols = data.columns
    if var not in cols:
        raise ValueError(
            f"Variable {var} specified in config file was not found in {snakemake.input[0]}."
        )
    data = (
        data.join(pl.DataFrame({var: values, "index": range(len(values))}), on=var)
        .sort("index")
        .select(pl.col("*").exclude("index"))
    )

data.write_parquet(snakemake.output[0])
