def get_value_column(wildcards):
    dataset = wildcards.dataset
    definition = lookup(f"datasets/{dataset}/value", within=config, default=None)
    if definition is None:
        n_vars = len(lookup(f"datasets/{dataset}/variables", within=config))
        data_path = lookup(f"datasets/{dataset}/data", within=config)
        with open(data_path) as f:
            header = f.readline().strip().split("\t")
            col_idx = lookup(f"datasets/{dataset}/value_column_index", within=config, default=None)
            if col_idx is not None:
                return header[col_idx]
            elif len(header) == n_vars + 1:
                return header[-1]
            else:
                raise WorkflowError(
                    "No value column specified in config but could not infer it from "
                    "data file {data_path}. Automatic inference only works if the data "
                    "file contains exactly #variables + 1 columns."
                )
    else:
        return definition