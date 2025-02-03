tests = lookup("tests", within=config, default={})

for test_name, test in tests.items():
    config["datasets"][test_name] = {
        "path": f"resources/simulation/n_vars:2.n_points:{test['n_points']}_seed:{test['seed']}_log2fc:{test['log2fc']}.tsv",
        "vars": ["var0", "var1"],
        "value": "value",
        "order": {
            "var0": ["a", "b"],
            "var1": ["a", "b"],
        },
    }


rule all_tests:
    input:
        collect("results/tests/{test}.passed", test=tests),


rule simulate_data:
    output:
        "resources/simulation/n_vars:{n_vars}.n_points:{n_points}_seed:{seed}_log2fc:{log2fc}.tsv",
    params:
        n_vars=evaluate("int({n_vars})"),
        n_points=evaluate("int({n_points})"),
        seed=evaluate("int({seed})"),
        log2fc=evaluate("float({log2fc})"),
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/simulate_data.py"


rule validate_results:
    input:
        "results/bootstrap/confidence_intervals/{test}.parquet",
    output:
        touch("results/tests/{test}.passed"),
    params:
        log2fc=lookup("tests/{test}/log2fc", within=config),
    script:
        "../scripts/validate_results.py"
