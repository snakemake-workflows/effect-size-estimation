This Snakemake workflow quantifies effect sizes between pairs of conditions while accounting for uncertainty. The analysis consists of two independent components: bootstrapping-based estimation of the empirical posterior distribution of the log2 fold change and statistical significance testing using the non-parametric Brunner-Munzel test. The workflow is entirely non-parametric and does not make assumptions about the distributions of the underlying data.

Workflow Steps:

* Bootstrap-Based Fold Change Estimation:
  * For each pairwise condition comparison, resample the measured values 1000 times with replacement while preserving the original sample size.
  * Compute the log2 fold change for each resampling iteration.
  * Aggregate these values to construct the empirical posterior distribution.
  * Determine a conservative fold change estimate by selecting the log2 fold change closest to 0, ensuring that 95% of the posterior distribution values are more extreme.

* Statistical Significance Testing:
  * Independently of the bootstrapping procedure, apply the Brunner-Munzel test to compute p-values for the null hypothesis that the probability of one condition yielding a greater value than the other is equal for both conditions.
  * As a non-parametric test, this approach remains valid even when the data do not follow a normal distribution or have unequal variances.

* Result Visualization & Reporting:
  * Generate primary plots displaying the conservative fold change estimates and p-values (using star notation).
  * Produce additional plots showing the empirical posterior distribution of the fold changes.
  * Create comprehensive datavzrd reports showing all pairwise comparisons, including effect sizes and confidence intervals, and underlying data.
  * All results are listed along with the respectively used code, software, and parameters.
  * For plots, legends can be activated/deactivated via the toggle button at the top left of the sidebar.