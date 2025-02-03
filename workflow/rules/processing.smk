rule sort_data:
    input:
        lookup("datasets/{dataset}/path", within=config),
    output:
        "results/data/{dataset}.sorted.parquet",
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
    params:
        vars=lookup("datasets/{dataset}/vars", within=config),
        n_bootstraps=config["n_bootstraps"],
        seed=config["seed"],
        value=lookup("datasets/{dataset}/value", within=config),
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/bootstrap.py"
