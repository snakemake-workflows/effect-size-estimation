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

The reliability of the conservative fold change estimate depends on how well the measured data points reflect the true (and unknown) distribution of the data.
In particular, the estimate is conservative in the sense that it is designed to avoid false positives, i.e. to avoid reporting a too high fold change.
However, this also means that the estimate can underestimate the true fold change.
If, e.g. by coincidence, the data points of two compared distributions are not representative, the estimate can be also too high without being able to properly control for that.
The more data points one can consider, the smaller is the chance that this happens.
Unlike, classical p-value based approaches though, the estimate provided here does not suffer from several well known issues like the winners curse, dichotomization into significant/non-significant (and the naturally occurring irreproducibility at the boundaries of the used signficance threshold), and most importantly the potential violation of assumptions about the underlying distributions (e.g. normality, homogeneous variance, etc.).
See `Halsey et al. (2015) <https://doi.org/10.1038/nmeth.3288>`__ for a comprehensive discussion of these issues.

{% if snakemake.wildcards.mode == "selected" %}
In addition, we performed a Brunner-Munzel test in order to obtain a p-value for each shown comparison.
The p-value tests the null hypothesis that for a pair of random values :math:`X` and :math:`Y` from the compared distributions the probability that :math:`X>Y` is equal to the probability that :math:`Y>X`.
We represent p-values :math:`\leq 0.0001` with ``****``, p-values :math:`\leq 0.001` with ``***``, p-values :math:`\leq 0.01` with ``**``, and p-values :math:`\leq 0.05` with ``*``.
The Brunner-Munzel test was chosen because it does not make assumptions about the underlying distributions (like e.g. equal variance).
{% endif %}