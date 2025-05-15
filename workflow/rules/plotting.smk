rule plot_effect_bootstrap_histograms:
    input:
        "results/bootstrap/histograms/{dataset}.parquet",
    output:
        report(
            "results/plots/{dataset}/bootstrap_histograms.html",
            category="Posterior fold change distribution plots",
            labels={
                "dataset": "{dataset}",
            },
            caption="../report/effect_histograms.rst",
        ),
    log:
        "logs/plot_effect_bootstrap_histograms/{dataset}.log",
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
            "results/plots/{dataset}/distributions_{mode,all|selected}_legend_{legend,yes|no}_effects.html",
            category="Distribution and fold change plots ({mode} comparisons)",
            labels={
                "dataset": "{dataset}",
                "legend": "{legend}",
            },
            caption="../report/distributions.rst",
        ),
    log:
        "logs/plot_dists_and_effects/{dataset}_{mode}_{legend}.log",
    params:
        vars=lookup("datasets/{dataset}/variables", within=config),
        value=get_value_column,
        color_scheme=lookup(dpath="color_scheme", within=config),
        min_fold_change=lookup(dpath="min_fold_change", within=config),
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/plot_distributions.py"
