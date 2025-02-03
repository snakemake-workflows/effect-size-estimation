rule plot_effect_bootstrap_histograms:
    input:
        "results/bootstrap/histograms/{dataset}.parquet"
    output:
        "results/plots/{dataset}/bootstrap_histograms.html"
    params:
        vars=lookup("datasets/{dataset}/vars", within=config),
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/plot_bootstrap_histograms.py"


rule plot_dists_and_effects:
    input:
        cis="results/bootstrap/confidence_intervals/{dataset}.parquet",
        data="results/data/{dataset}.sorted.parquet",
    output:
        "results/plots/{dataset}/distributions.html"
    params:
        vars=lookup("datasets/{dataset}/vars", within=config),
        value=lookup("datasets/{dataset}/value", within=config),
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/plot_distributions.py"