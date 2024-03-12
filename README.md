# Cluster-Robust Heteroskedasticity-Consistent Standard Errors

This repository explains and illustrates a variety of cluster-robust heteroskedasticity-consistent standard error estimators. It was created by <a rel="creator" href="https://twitter.com/DuzgunDilsiz">Düzgün Dilsiz</a> (University of Basel) for teaching purposes. Note that there is also [a website version](https://d-dilsiz.github.io/cluster-robust/).

Heteroskedasticity-Consistent (HC) standard errors (SE) allow different units having different variance of the error term, rather than assuming all units having the same variance of the error term. In addition to this, cluster-robust HC SE allow for correlation within a cluster. In the following the estimated asymptotic variance for the panel data fixed effects estimator is presented for the standard `small sample size correction` (`Stata` and `R`), and the `HC0-HC3` estimators, as computed by the function `vcovHC` (`plm` package) in `R`.

STATA SSS: $$\widehat{Avar}\left[\widehat \beta_{FE}\right]$ = $\frac{N}{N-1}$ $\frac{NT-1}{NT-K-1}$ $\left(\ddot X' \ddot X \right)^{-1}$ $\Sigma_{i=1}^{N}{\ddot X_{i}'\widehat{\ddot u_{i}} \widehat{\ddot u_{i}'}\ddot X_{i}}$ $\left(\ddot X' \ddot X\right)^{-1}$$

R SSS: $\widehat{Avar}\left[\widehat \beta_{FE}\right]$ = $\frac{N}{N-1}$ $\frac{NT-1}{NT-K}$ $\left(\ddot X' \ddot X \right)^{-1}$ $\Sigma_{i=1}^{N}{\ddot X_{i}'\widehat{\ddot u_{i}} \widehat{\ddot u_{i}'}\ddot X_{i}}$ $\left(\ddot X' \ddot X\right)^{-1}$

R HC0: $\widehat{Avar}\left[\widehat \beta_{FE}\right]$ = $\left(\ddot X' \ddot X \right)^{-1}$ $\Sigma_{i=1}^{N}{\ddot X_{i}'\widehat{\ddot u_{i}} \widehat{\ddot u_{i}'}\ddot X_{i}}$ $\left(\ddot X' \ddot X\right)^{-1}$

R HC1: $\widehat{Avar}\left[\widehat \beta_{FE}\right]$ = $\frac{NT}{NT-K}$ $\left(\ddot X' \ddot X \right)^{-1}$ $\Sigma_{i=1}^{N}{\ddot X_{i}'\widehat{\ddot u_{i}} \widehat{\ddot u_{i}'}\ddot X_{i}}$ $\left(\ddot X' \ddot X\right)^{-1}$

R HC2: $\widehat{Avar}\left[\widehat \beta_{FE}\right]$ = $\left(\ddot X' \ddot X \right)^{-1}$ $\Sigma_{i=1}^{N}{\ddot X_{i}'\tilde{\ddot u_{i}} \tilde{\ddot u_{i}}'\ddot X_{i}}$ $\left(\ddot X' \ddot X\right)^{-1}$

R HC3: $\widehat{Avar}\left[\widehat \beta_{FE}\right]$ = $\left(\ddot X' \ddot X \right)^{-1}$ $\Sigma_{i=1}^{N}{\ddot X_{i}'\tilde{\ddot u_{i}} \tilde{\ddot u_{i}}'\ddot X_{i}}$ $\left(\ddot X' \ddot X\right)^{-1}$

where $\tilde{\ddot u_{it}}=\widehat{\ddot u_{i}}(1-\ddot h_{it})^{-\delta_{m}}$ where m=0.5 for `HC2` and $m=1$ for `HC3`, and $\ddot h_{it}$ indicating the diagonal elements of $\ddot H_{jj}$ of the hat matrix $\ddot H = \ddot X\left(\ddot X' \ddot X \right)^{-1}\ddot X'$.

In the code, you find a small example with crime data from the United States where all the estimators are computed and it is shown that the results are identical to when using the provided functions directly.

## File structure

This repository consists of 2 folders:
* [code](https://github.com/d-dilsiz/cluster-robust/tree/main/code): provides the `R`-code for cluster-robust `HC0-HC3` standard error estimators with examples
* [pictures](https://github.com/d-dilsiz/cluster-robust/tree/main/pictures): resources for readme.md

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
