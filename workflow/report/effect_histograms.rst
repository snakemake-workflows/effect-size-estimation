Histograms of empirical posterior distributions of log2 fold changes between groups of measured {{ snakemake.params.value }} across {{ snakemake.params.vars|join("/") }}.
Distributions were obtained by computing 1000 bootstraps of each group of measured values.

Horizontal axis labels are converted to fold changes by taking :math:`2^f` for any log2 fold change :math:`f`.
Direction is annotated as ``⏶`` for up-regulation, ``⏷`` for down-regulation, and ``=`` for no change.