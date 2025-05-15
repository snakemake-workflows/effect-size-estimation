rule parquet_to_tsv:
    input:
        "{prefix}.parquet",
    output:
        "{prefix}.tsv"
    log:
        "logs/parquet_to_tsv/{prefix}.log"
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/parquet_to_tsv.py"


rule datavzrd_report:
    input:
        config=workflow.source_path("../resources/datavzrd_config.yaml"),
        cis="results/bootstrap/confidence_intervals/{dataset}.tsv",
        data="results/data/{dataset}.sorted.tsv",
    output:
        report(
            directory("results/datavzrd/{dataset}"),
            htmlindex="index.html",
            category="Tabular data",
            labels={
                "dataset": "{dataset}",
            },
            caption="../report/datavzrd.rst",
        )
    log:
        "logs/datavzrd_report/{dataset}.log"
    params:
        value_column=get_value_column,
        variables=lookup("datasets/{dataset}/variables", within=config),
        effect_plot_selected="results/plots/{dataset}/distributions_selected_legend_yes_effects.html",
        effect_plot_all="results/plots/{dataset}/distributions_all_legend_yes_effects.html",
        hist_plot="results/plots/{dataset}/bootstrap_histograms.html",
    wrapper:
        "v5.8.2/utils/datavzrd"
