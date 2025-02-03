Measured {{ snakemake.params.value }} across {{ snakemake.params.vars|join("/") }}.
The upper panel shows the data distributions.
Below each distribution, we depict a conservative, uncertainty-aware estimate of the fold change against each other distribution right of it.
The estimate can be interpreted as the strongest fold change such that the true fold change is stronger with a probability of :math:`\geq 95%`.
Depicted fold changes are annotated as
up-regulation: ``⏶``, down-regulation: ``⏷``, no change: ``=``, and almost no change (up/down fold change :math:`<1.05`): ≈.

The estimate is obtained as follows:
We compute 1000 bootstraps of each distribution, and then calculate an empirical posterior distribution of the fold change between each pair of groups.
If the :math:`95%` credible interval of the empirical posterior encloses :math:`1.0` the conservative fold change is :math:`1.0` (i.e. indicating "no change").
Otherwise, the conservative fold change is either the :math:`5%` or the :math:`95%` quantile, whichever is closer to :math:`1.0` (calculated in log2 space such that the distances are symmetric).