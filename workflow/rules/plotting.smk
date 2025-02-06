rule plot_effect_bootstrap_histograms:
    input:
        "results/bootstrap/histograms/{dataset}.parquet",
    output:
        report(
            "results/plots/{dataset}/bootstrap_histograms.html",
            category="{dataset}",
            labels={
                "plot": "Posterior fold change distributions",
                "effects": "all",
            },
            caption="../report/effect_histograms.rst",
        ),
    params:
        vars=lookup("datasets/{dataset}/variables", within=config),
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/plot_bootstrap_histograms.py"


rule plot_dists_and_effects:
    input:
        cis="results/bootstrap/confidence_intervals/{dataset}.parquet",
        data="results/data/{dataset}.sorted.parquet",
        comparisons=local(lookup("datasets/{dataset}/comparisons", within=config, default=[])),
    output:
        report(
            "results/plots/{dataset}/distributions_{mode,all|selected}_effects.html",
            category="{dataset}",
            labels={
                "plot": "Data distributions and conservative fold changes",
                "effects": "{mode}",
            },
            caption="../report/distributions.rst",
        ),
    params:
        vars=lookup("datasets/{dataset}/variables", within=config),
        value=get_value_column,
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/plot_distributions.py"
