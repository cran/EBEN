
<!-- README.md is generated from the source: README.Rmd -->

# Empirical Bayesian Elastic Net (EBEN) for Generalized Linear Models <img src="man/figures/logo.png" width="100" align="right" />

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/EBEN)](https://cran.r-project.org/package=EBEN)[![](https://cranlogs.r-pkg.org/badges/EBEN)](https://CRAN.R-project.org/package=EBEN)

We provide extremely efficient procedures for fitting the empirical
Bayesian methods with lasso and elastic net hierarchical priors for
linear regression (gaussian), and logistic regression (binomial) models.
`EBEN` is a sister package to `EBglmnet` (available in CRAN). Both
packages share key features include:

- sparse variable selection and effect estimation via generalized linear
  regression models;
- high dimensionality with p\>\>n; and
- significance test (with output of `p-value`) for nonzero effects; and
- closed-form solution for Bayesian variance estimation in an iterative
  cooridinate descent algorithm estimating the Bayesian means.

The implementation enables extremely efficient computation comparable
with that of `glmnet` package.

### When you need `EBEN`

While `EBglmnet` offers generic functions for a broad range of use
cases, `EBEN` takes care of the following special cases:

- two-way interaction terms (`epistasis`) are included with
  `epis = TRUE`: for input independent parameter `X` with n x p
  dimension, the functions will evaluate p(p-1)/2 additional parameters;
- group Empirical Bayesian Lasso are avaiable with `group = TRUE`: the
  penalty parameter for the group of p(p-1)/2 parameters are weighted
  with group size in comparing with the group origin p variables.

### Further readings

Details may be found in Huang A. and Liu D ([2016](#ref-package)), Huang
A., Xu S., and Cai X. ([2015](#ref-EBEN)), Huang A.
([2014](#ref-dissertation)), Huang A., Xu S., and Cai X.
([2013](#ref-EBlasso)), and Cai X., Huang A., and Xu S.,
([2011](#ref-QTL)).

### Version notes

Version 5.1 is a major release with several new features, including:

- group Empirical Bayesian Lasso (EBlasso) and built-in two-way
  interaction support moved to `EBEN` package.
- BLAS/Lapack routines are updated according to R-API change.

## References

<div id="refs" class="references">

<div id="ref-package">

<p>
Huang A., Liu D., (2016) <br> EBglmnet: a comprehensive R package for
sparse generalized linear regression models <br> Bioinformatics, Volume
37, Issue 11, Pages 1627–1629
</p>

</div>

<div id="ref-EBEN">

<p>
Huang A., Xu S., and Cai X. (2015). <br> Empirical Bayesian elastic net
for multiple quantitative trait locus mapping.</a><br>
<em>Heredity</em>, Vol. 114(1), 107-115.
</p>

</div>

<div id="ref-dissertation">

<p>
Huang A. (2014) <br> Sparse Model Learning for Inferring Genotype and
Phenotype Associations. <br> Ph.D Dissertation, University of Miami,
Coral Gables, FL, USA.
</p>

</div>

<div id="ref-EBlasso">

<p>
Huang A., Xu S., and Cai X. (2013). <br> Empirical Bayesian
LASSO-logistic regression for multiple binary trait locus mapping.
</a><br> <em>BMC Genetics</em>, 14(1),5.
</p>

</div>

<div id="ref-QTL">

<p>
Cai X., Huang A., and Xu S., (2011). <br> Fast empirical Bayesian LASSO
for multiple quantitative trait locus mapping. </a><br> <em>BMC
Bioinformatics</em>, 12(1),211.
</p>

</div>

</div>
