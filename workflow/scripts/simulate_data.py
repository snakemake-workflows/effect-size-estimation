import polars as pl
import random
from scipy.stats import norm

random.seed(snakemake.params.seed)

var_names = [f"var{i}" for i in range(snakemake.params.n_vars)]
var_values = [val for val in ["a", "b"] for _ in range(snakemake.params.n_points)]
values = [float(value) for factor in [1, 2 ** snakemake.params.log2fc] for value in norm.rvs(loc=5.0 * factor, size=snakemake.params.n_points)]

data = pl.DataFrame({var: var_values for var in var_names} | {"value": values})
data.write_csv(snakemake.output[0], separator="\t")
