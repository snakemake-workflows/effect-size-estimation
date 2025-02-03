import polars as pl

data = pl.read_csv(snakemake.input[0], separator="\t")

for var, values in reversed(snakemake.params.order.items()):
    cols = data.columns
    if var not in cols:
        raise ValueError(
            f"Variable {var} specified in config file was not found in {snakemake.input[0]}."
        )
    joined = data.join(pl.DataFrame({var: values, "index": range(len(values))}), how="left", on=var)
    if joined.filter(pl.col("index").is_null()).height > 0:
        raise ValueError(
            "Dataset variable values are missing in the order defined in the config file."
        )
    data = (
        joined.sort("index")
        .select(pl.col("*").exclude("index"))
    )

data.write_parquet(snakemake.output[0])
