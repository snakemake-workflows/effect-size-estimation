def get_value_column(wildcards):
    dataset = wildcards.dataset
    definition = lookup(f"datasets/{dataset}/value", within=config, default=None)
    if definition is None:
        n_vars = len(lookup(f"datasets/{dataset}/variables", within=config))
        data_path = lookup(f"datasets/{dataset}/path", within=config)
        with open(data_path) as f:
            header = f.readline().strip().split("\t")
            if len(header) == n_vars + 1:
                return header[-1]
            else:
                raise WorkflowError(
                    "No value column specified in config but could not infer it from "
                    "data file {data_path}. Automatic inference only works if the data "
                    "file contains exactly #variables + 1 columns."
                )