rule sort_data:
    input:
        local(lookup("datasets/{dataset}/data", within=config)),
    output:
        "results/data/{dataset}.sorted.parquet",
    log:
        "logs/sort_data/{dataset}.log",
    params:
        order=lookup("datasets/{dataset}/order", within=config),
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/sort_order.py"


rule bootstrap:
    input:
        data="results/data/{dataset}.sorted.parquet",
    output:
        hists="results/bootstrap/histograms/{dataset}.parquet",
        cis="results/bootstrap/confidence_intervals/{dataset}.parquet",
    log:
        "logs/bootstrap/{dataset}.log",
    params:
        vars=lookup("datasets/{dataset}/variables", within=config),
        n_bootstraps=config["n_bootstraps"],
        seed=config["seed"],
        value=get_value_column,
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/bootstrap.py"
