To configure the workflow, edit the config directory.
It contains the `config/config.yaml` which has as central element the `datasets` section.
In there, modify or add new subsections, each representing one dataset to process, with the intended dataset name as the key.

* For each dataset define the path to a data table in `.tsv` format as can be seen in the example.
* Define the columns that contain the variables to compare against using the `variables` key.
* Define the column that contains the measurement data that shall be compared in terms of fold change and significance using the `value` key. The `value` key is optional and will be inferred automatically as the only non-variable column if left out.
* Define the order in which variable values are compared and occur in the generated visualizations via the `order` key.

Outside of the `datasets` section you can set various technical parameters as instructed inside of the template `config/config.yaml`.