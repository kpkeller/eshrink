# eshrink
[![CRAN Status](http://www.r-pkg.org/badges/version/eshrink)](https://cran.r-project.org/package=eshrink)


Shrinkage estimators for estimating regression parameters

This R package provides functions for estimating the penalization parameter
for shrinkage estimators using the approach of Keller and Rice (2017). 
The package currently contains functionality for ridge regressiona and the LASSO.
The penalty parameter is selected by minimizing bias and/or variance in data
generated from the posterior predictive distribution.

## References
Keller JP and Rice KM. (2017) Selecting Shrinkage Parameters for Effect Estimation:
the Multi-Ethnic Study of Atherosclerosis. *American Journal of Epidemiology*. 
https://doi.org/10.1093/aje/kwx225
