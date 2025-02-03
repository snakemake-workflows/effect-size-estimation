tests = lookup("tests", within=config, default={})

for test_name, test in tests.items():
    config["datasets"][test_name] = {
        "path": f"resources/simulation/n_vars:2.n_points:{test['n_points']}_seed:{test['seed']}_log2fc:{test['log2fc']}.tsv",
        "variables": ["var0", "var1"],
        "value": "value",
        "order": {
            "var0": ["a", "b"],
            "var1": ["a", "b"],
        },
    }


rule all_tests:
    input:
        "results/tests.tsv"

rule simulate_data:
    output:
        "resources/simulation/n_vars:{n_vars}.n_points:{n_points}_seed:{seed}_log2fc:{log2fc}.tsv",
    params:
        n_vars=evaluate("int({n_vars})"),
        n_points=evaluate("int({n_points})"),
        seed=evaluate("int({seed})"),
        log2fc=evaluate("float({log2fc})"),
    conda:
        "../envs/simulation.yaml"
    script:
        "../scripts/simulate_data.py"


rule validate_results:
    input:
        cis=collect("results/bootstrap/confidence_intervals/{test}.parquet", test=tests),
        raw=[lookup(f"datasets/{test}/path", within=config) for test in tests],
    log:
        "results/tests.tsv",
    params:
        tests=tests,
    conda:
        "../envs/simulation.yaml"
    script:
        "../scripts/validate_results.py"